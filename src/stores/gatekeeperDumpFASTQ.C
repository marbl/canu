
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2012, J. Craig Venter Institute.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id$";

#include "AS_global.H"
#include "gkStore.H"
#include "AS_UTL_fileIO.H"
#include "AS_UTL_fasta.H"


//  Write sequence in multiple formats.  This used to write to four fastq files, the .1, .2, .paired and .unmated.
//  It's left around for future expansion to .fastq and .bax.h5.
//
class libOutput {
public:
  libOutput(char const *outPrefix, char const *libName) {
    strcpy(_p, outPrefix);
    strcpy(_n, libName);
  };

  ~libOutput() {
    if (_FASTA)
      fclose(_FASTA);
    if (_FASTQ)
      fclose(_FASTQ);
  };

  FILE  *getFASTQ(void) {

    if (_FASTQ)
      return(_FASTQ);

    return(openFASTQ());
  };

  FILE  *openFASTQ(void) {
    char  N[FILENAME_MAX];

    sprintf(N, "%s.%s.fastq", _p, _n);

    errno = 0;
    _FASTQ = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Failed to open FASTQ output file '%s': %s\n", N, strerror(errno)), exit(1);

    return(_FASTQ);
  };


  FILE  *getFASTA(void) {

    if (_FASTA)
      return(_FASTA);

    return(openFASTA());
  };

  FILE  *openFASTA(void) {
    char  N[FILENAME_MAX];

    sprintf(N, "%s.%s.fasta", _p, _n);

    errno = 0;
    _FASTA = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Failed to open FASTA output file '%s': %s\n", N, strerror(errno)), exit(1);

    return(_FASTA);
  };

private:
  char   _p[FILENAME_MAX];
  char   _n[FILENAME_MAX];

  FILE  *_FASTA;
  FILE  *_FASTQ;
};    




int
main(int argc, char **argv) {
  char            *gkpStoreName      = NULL;
  char            *outPrefix         = NULL;

  uint32           libToDump         = 0;

  uint32           bgnID             = 0;
  uint32           endID             = AS_MAX_READS;

  bool             dumpAllBases      = true;
  bool             dumpAllReads      = true;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-l") == 0) {
      libToDump = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endID  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];

    } else {
      err++;
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (gkpStoreName == NULL)
    err++;
  if (outPrefix == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s [...] -o fastq-prefix -g gkpStore\n", argv[0]);
    fprintf(stderr, "  -g gkpStore\n");
    fprintf(stderr, "  -o fastq-prefix     write files fastq-prefix.1.fastq, fastq-prefix.2.fastq, fastq-prefix.paired.fastq, fastq-prefix.unmated.fastq\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  -l libToDump        output only fragments in library number libToDump (NOT IMPLEMENTED)\n");
    fprintf(stderr, "  -b id               output starting at fragment id\n");
    fprintf(stderr, "  -e id               output stopping after fragment id\n");
    fprintf(stderr, "  \n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-g) supplied.\n");
    if (outPrefix == NULL)
      fprintf(stderr, "ERROR: no output prefix (-o) supplied.\n");
    exit(1);
  }

  gkStore    *gkpStore  = new gkStore(gkpStoreName);
  uint32      numReads  = gkpStore->gkStore_getNumReads();
  uint32      numLibs   = gkpStore->gkStore_getNumLibraries();

  if (numReads < endID)
    endID = numReads;

  fprintf(stderr, "Dumping reads from %u to %u (inclusive).\n", bgnID, endID - 1);

  libOutput   **out     = new libOutput * [numLibs];

  out[0] = NULL;  //  There isn't a zeroth library.

  for (uint32 i=1; i<=numLibs; i++)
    out[i] = new libOutput(outPrefix, gkpStore->gkStore_getLibrary(i)->gkLibrary_libraryName());


  for (uint32 rid=bgnID; rid<=endID; rid++) {
    gkRead      *read     = gkpStore->gkStore_getRead(rid);
    gkReadData  *readData = new gkReadData;

    uint32   lclr   = read->gkRead_clearRegionBegin();
    uint32   rclr   = read->gkRead_clearRegionEnd();

    uint32   ldump  = read->gkRead_clearRegionBegin();
    uint32   rdump  = read->gkRead_clearRegionEnd();

    if (dumpAllBases) {
      ldump = 0;
      rdump = read->gkRead_sequenceLength();
    }

    uint32  libID = read->gkRead_libraryID();

    if ((dumpAllReads == false) && (read->gkRead_isDeleted() == true))
      //  Fragment is deleted, don't dump.
      continue;

    if ((libToDump != 0) && (libID == libToDump))
      //  Fragment isn't marked for dumping, don't dump.
      continue;

    if ((dumpAllBases == false) && (lclr >= rclr))
      //  Fragment has null or invalid clear range, don't dump.
      continue;

    gkpStore->gkStore_loadReadData(read, readData);

    char *seq = readData->gkReadData_getSequence()  + ldump;
    char *qlt = readData->gkReadData_getQualities() + ldump;

    seq[rdump - ldump] = 0;
    qlt[rdump - ldump] = 0;

    //  Soft mask not-clear bases

    if (dumpAllBases == true) {
      for (uint32 i=0; i<lclr; i++)
        seq[i] += (seq[i] >= 'A') ? 'a' - 'A' : 0;

      for (uint32 i=lclr; i<rclr; i++)
        seq[i] += (seq[i] >= 'A') ? 0 : 'A' - 'a';

      for (uint32 i=rclr; seq[i]; i++)
        seq[i] += (seq[i] >= 'A') ? 'a' - 'A' : 0;
    }

    AS_UTL_writeFastQ(out[libID]->getFASTQ(), seq, (rdump - ldump), qlt, (rdump - ldump),
                      "@"F_U32" clr="F_U32","F_U32"\n",
                      rid, lclr, rclr);

    AS_UTL_writeFastQ(stdout, seq, (rdump - ldump), qlt, (rdump - ldump),
                      "@"F_U32" clr="F_U32","F_U32"\n",
                      rid, lclr, rclr);
  }

  delete gkpStore;

  exit(0);
}

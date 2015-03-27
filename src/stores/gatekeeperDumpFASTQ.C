
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

#include "clearRangeFile.H"

//  Write sequence in multiple formats.  This used to write to four fastq files, the .1, .2, .paired and .unmated.
//  It's left around for future expansion to .fastq and .bax.h5.
//
class libOutput {
public:
  libOutput(char const *outPrefix, char const *libName = NULL) {
    strcpy(_p, outPrefix);

    if (libName)
      strcpy(_n, libName);
    else
      _n[0] = 0;

    _FASTA = NULL;
    _FASTQ = NULL;
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

    if (_n[0])
      sprintf(N, "%s.%s.fastq", _p, _n);
    else
      sprintf(N, "%s.fastq", _p);

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

    if (_n[0])
      sprintf(N, "%s.%s.fasta", _p, _n);
    else
      sprintf(N, "%s.fasta", _p);

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

  char            *clrName           = NULL;

  uint32           libToDump         = 0;

  uint32           bgnID             = 1;
  uint32           endID             = AS_MAX_READS;

  bool             dumpAllReads      = false;
  bool             dumpAllBases      = false;
  bool             dumpOnlyDeleted   = false;

  bool             dumpFASTQ         = true;
  bool             dumpFASTA         = false;

  bool             withLibName       = true;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];


    } else if (strcmp(argv[arg], "-c") == 0) {
      clrName = argv[++arg];


    } else if (strcmp(argv[arg], "-l") == 0) {
      libToDump = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endID  = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-r") == 0) {
      bgnID  = atoi(argv[++arg]);
      endID  = bgnID;


    } else if (strcmp(argv[arg], "-allreads") == 0) {
      dumpAllReads    = true;

    } else if (strcmp(argv[arg], "-allbases") == 0) {
      dumpAllBases    = true;

    } else if (strcmp(argv[arg], "-onlydeleted") == 0) {
      dumpOnlyDeleted = true;
      dumpAllReads    = true;  //  Otherwise we won't report the deleted reads!


    } else if (strcmp(argv[arg], "-fastq") == 0) {
      dumpFASTQ       = true;
      dumpFASTA       = false;

    } else if (strcmp(argv[arg], "-fasta") == 0) {
      dumpFASTQ       = false;
      dumpFASTA       = true;

    } else if (strcmp(argv[arg], "-nolibname") == 0) {
      withLibName     = false;


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
    fprintf(stderr, "  -G gkpStore\n");
    fprintf(stderr, "  -o fastq-prefix     write files fastq-prefix.(libname).fastq, ...\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -l libToDump        output only read in library number libToDump (NOT IMPLEMENTED)\n");
    fprintf(stderr, "  -b id               output starting at read 'id'\n");
    fprintf(stderr, "  -e id               output stopping after read 'id'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c clearFile        clear range file from OBT modules\n");
    fprintf(stderr, "  -allreads           if a clear range file, lower case mask the deleted reads\n");
    fprintf(stderr, "  -allbases           if a clear range file, lower case mask the non-clear bases\n");
    fprintf(stderr, "  -onlydeleted        if a clear range file, only output deleted reads (the entire read)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -r id               output only the single read 'id'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -fastq              output is FASTQ format (with extension .fastq, default)\n");
    fprintf(stderr, "  -fasta              output is FASTA format (with extension .fasta)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nolibname          don't include the library name in the output file name\n");
    fprintf(stderr, "\n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-G) supplied.\n");
    if (outPrefix == NULL)
      fprintf(stderr, "ERROR: no output prefix (-o) supplied.\n");

    exit(1);
  }

  gkStore        *gkpStore  = new gkStore(gkpStoreName);
  uint32          numReads  = gkpStore->gkStore_getNumReads();
  uint32          numLibs   = gkpStore->gkStore_getNumLibraries();

  clearRangeFile *clrRange  = (clrName == NULL) ? NULL : new clearRangeFile(clrName, gkpStore);

  if (bgnID < 1)
    bgnID = 1;

  if (numReads < endID)
    endID = numReads;

  if (endID < bgnID)
    fprintf(stderr, "No reads to dump; reversed ranges make no sense: bgn="F_U32" end="F_U32"??\n", bgnID, endID);




  fprintf(stderr, "Dumping reads from %u to %u (inclusive).\n", bgnID, endID);

  libOutput   **out = new libOutput * [numLibs + 1];

  out[0] = NULL;  //  There isn't a zeroth library.

  for (uint32 i=1; i<=numLibs; i++)
    out[i] = new libOutput(outPrefix, (withLibName == false) ? NULL : gkpStore->gkStore_getLibrary(i)->gkLibrary_libraryName());

  gkReadData   *readData = new gkReadData;

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    gkRead      *read   = gkpStore->gkStore_getRead(rid);

    uint32       libID  = read->gkRead_libraryID();

    uint32       lclr   = 0;
    uint32       rclr   = read->gkRead_sequenceLength();
    bool         ignore = false;

    //  If a clear range file is supplied, grab the clear range.  If it hasn't been set, default
    //  back to the entire read.

    if (clrRange) {
      lclr   = clrRange->bgn(rid);
      rclr   = clrRange->end(rid);
      ignore = clrRange->isDeleted(rid);
    }

    //  Abort if we're not dumping anything from this read
    //   - not in a library we care about
    //   - deleted, and not dumping all reads
    //   - not deleted, but only reporting deleted reads

    if (((libToDump != 0) && (libID == libToDump)) ||
        ((dumpAllReads == false) && (ignore == true)) ||
        ((dumpOnlyDeleted == true) && (ignore == false)))
      continue;

    //  And if we're told to ignore the read, and here, then the read was deleted and we're printing
    //  all reads.  Reset the clear range to the whole read, the clear range is invalid.

    if (ignore) {
      rclr = read->gkRead_sequenceLength();
      lclr = 0;
    }

    //  Grab the sequence and quality.

    gkpStore->gkStore_loadReadData(read, readData);

    char   *seq = readData->gkReadData_getSequence();
    char   *qlt = readData->gkReadData_getQualities();
    uint32  len = rclr - lclr;

    //  Soft mask not-clear bases

    if (dumpAllBases == true) {
      for (uint32 i=0; i<lclr; i++)
        seq[i] += (seq[i] >= 'A') ? 'a' - 'A' : 0;

      for (uint32 i=lclr; i<rclr; i++)
        seq[i] += (seq[i] >= 'A') ? 0 : 'A' - 'a';

      for (uint32 i=rclr; seq[i]; i++)
        seq[i] += (seq[i] >= 'A') ? 'a' - 'A' : 0;

      rclr = read->gkRead_sequenceLength();
      lclr = 0;
    }

    //  Chop off the ends we're not printing.

    seq += lclr;

    seq[len] = 0;
    qlt[len] = 0;

    //  And print the read.
    if (dumpFASTQ)
      AS_UTL_writeFastQ(out[libID]->getFASTQ(), seq, len, qlt, len,
                        "@"F_U32" clr="F_U32","F_U32"\n",
                        rid, lclr, rclr);

    if (dumpFASTA)
      AS_UTL_writeFastA(out[libID]->getFASTA(), seq, len, 0,
                        ">"F_U32" clr="F_U32","F_U32"\n",
                        rid, lclr, rclr);
  }

  delete   readData;

  for (uint32 i=1; i<=numLibs; i++)
    delete out[i];

  delete [] out;

  delete    gkpStore;

  exit(0);
}

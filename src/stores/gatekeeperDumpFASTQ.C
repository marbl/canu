
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_GKP/gkpStoreDumpFASTQ.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-FEB-06 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-21 to 2015-JUN-03
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-23
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "AS_UTL_decodeRange.H"
#include "AS_UTL_fileIO.H"
#include "AS_UTL_fasta.H"

#include "clearRangeFile.H"

//  Write sequence in multiple formats.  This used to write to four fastq files, the .1, .2, .paired and .unmated.
//  It's left around for future expansion to .fastq and .bax.h5.
//
class libOutput {
public:
  libOutput(char const *outPrefix, char const *outSuffix, char const *libName = NULL) {
    strcpy(_p, outPrefix);

    if (outSuffix[0])
      sprintf(_s, ".%s", outSuffix);
    else
      _s[0] = 0;

    if (libName)
      strcpy(_n, libName);
    else
      _n[0] = 0;

    _WRITER = NULL;
    _FASTA  = NULL;
    _FASTQ  = NULL;
  };

  ~libOutput() {
    if (_WRITER)
      delete _WRITER;
  };

  FILE  *getFASTQ(void) {

    if (_FASTQ)
      return(_FASTQ);

    return(openFASTQ());
  };

  FILE  *openFASTQ(void) {
    char  N[FILENAME_MAX];

    if (_n[0])
      sprintf(N, "%s.%s.fastq%s", _p, _n, _s);
    else
      sprintf(N, "%s.fastq%s", _p, _s);

    if ((_p[0] == '-') && (_p[1] == 0)) {
      sprintf(N, "(stdout)");
      _FASTQ = stdout;
    }

    else {
      _WRITER = new compressedFileWriter(N);
      _FASTQ  = _WRITER->file();
    }

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
      sprintf(N, "%s.%s.fasta%s", _p, _n, _s);
    else
      sprintf(N, "%s.fasta%s", _p, _s);

    if ((_p[0] == '-') && (_p[1] == 0)) {
      sprintf(N, "(stdout)");
      _FASTA = stdout;
    }

    else {
      _WRITER = new compressedFileWriter(N);
      _FASTA  = _WRITER->file();
    }

    return(_FASTA);
  };

private:
  char   _p[FILENAME_MAX];
  char   _s[FILENAME_MAX];
  char   _n[FILENAME_MAX];

  compressedFileWriter  *_WRITER;
  FILE                  *_FASTA;
  FILE                  *_FASTQ;
};




char *
scanPrefix(char *prefix) {
  int32 len = strlen(prefix);

  if ((len > 3) && (strcasecmp(prefix + len - 3, ".gz") == 0)) {
    prefix[len-3] = 0;
    return(prefix + len - 2);
  }

  if ((len > 4) && (strcasecmp(prefix + len - 4, ".bz2") == 0)) {
    prefix[len-4] = 0;
    return(prefix + len - 3);
  }

  if ((len > 3) && (strcasecmp(prefix + len - 3, ".xz") == 0)) {
    prefix[len-3] = 0;
    return(prefix + len - 2);
  }

  return(prefix + len);
}




int
main(int argc, char **argv) {
  char            *gkpStoreName      = NULL;
  char            *outPrefix         = NULL;
  char            *outSuffix         = NULL;

  char            *clrName           = NULL;

  uint32           libToDump         = 0;

  uint32           bgnID             = 1;
  uint32           endID             = UINT32_MAX;

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
      outSuffix = scanPrefix(outPrefix);


    } else if (strcmp(argv[arg], "-c") == 0) {
      clrName = argv[++arg];


    } else if (strcmp(argv[arg], "-l") == 0) {
      libToDump = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-b") == 0) {   //  DEPRECATED!
      bgnID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {   //  DEPRECATED!
      endID  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-r") == 0) {
      AS_UTL_decodeRange(argv[++arg], bgnID, endID);

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
    fprintf(stderr, "                      if fastq-prefix is '-', all sequences output to stdout\n");
    fprintf(stderr, "                      if fastq-prefix ends in .gz, .bz2 or .xz, output is compressed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -l libToDump        output only read in library number libToDump (NOT IMPLEMENTED)\n");
    fprintf(stderr, "  -r id[-id]          output only the single read 'id', or the specified range of ids\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c clearFile        clear range file from OBT modules\n");
    fprintf(stderr, "  -allreads           if a clear range file, lower case mask the deleted reads\n");
    fprintf(stderr, "  -allbases           if a clear range file, lower case mask the non-clear bases\n");
    fprintf(stderr, "  -onlydeleted        if a clear range file, only output deleted reads (the entire read)\n");
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

  gkStore        *gkpStore  = gkStore::gkStore_open(gkpStoreName);
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

  //  Allocate outputs.  If withLibName == false, all reads will artificially be in lib zero, the
  //  other files won't ever be created.  Otherwise, the zeroth file won't ever be created.

  out[0] = new libOutput(outPrefix, outSuffix, NULL);

  for (uint32 i=1; i<=numLibs; i++)
    out[i] = new libOutput(outPrefix, outSuffix, gkpStore->gkStore_getLibrary(i)->gkLibrary_libraryName());

  //  Grab a new readData, and iterate through reads to dump.

  gkReadData   *readData = new gkReadData;

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    gkRead      *read   = gkpStore->gkStore_getRead(rid);

    uint32       libID  = (withLibName == false) ? 0 : read->gkRead_libraryID();

    uint32       flen   = read->gkRead_sequenceLength();
    uint32       lclr   = 0;
    uint32       rclr   = flen;
    bool         ignore = false;

    //fprintf(stderr, "READ %u claims id %u length %u in lib %u\n", rid, read->gkRead_readID(), read->gkRead_sequenceLength(), libID);

    //  If a clear range file is supplied, grab the clear range.  If it hasn't been set, the default
    //  is the entire read.

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
      lclr = 0;
      rclr = read->gkRead_sequenceLength();
    }

    //  Grab the sequence and quality.

    gkpStore->gkStore_loadReadData(read, readData);

    char   *seq  = readData->gkReadData_getSequence();
    char   *qlt  = readData->gkReadData_getQualities();
    uint32  clen = rclr - lclr;

    //  Soft mask not-clear bases

    if (dumpAllBases == true) {
      for (uint32 i=0; i<lclr; i++)
        seq[i] += (seq[i] >= 'A') ? 'a' - 'A' : 0;

      for (uint32 i=lclr; i<rclr; i++)
        seq[i] += (seq[i] >= 'A') ? 0 : 'A' - 'a';

      for (uint32 i=rclr; flen; i++)
        seq[i] += (seq[i] >= 'A') ? 'a' - 'A' : 0;

      lclr = 0;
      rclr = flen;
    }

    //  Chop off the ends we're not printing.

    seq += lclr;
    qlt += lclr;

    seq[clen] = 0;
    qlt[clen] = 0;

    //  Print the read.

    if (dumpFASTA)
      AS_UTL_writeFastA(out[libID]->getFASTA(), seq, clen, 100,
                        ">"F_U32" clr="F_U32","F_U32"\n",
                        rid, lclr, rclr);

    if (dumpFASTQ)
      AS_UTL_writeFastQ(out[libID]->getFASTQ(), seq, clen, qlt, clen,
                        "@"F_U32" clr="F_U32","F_U32"\n",
                        rid, lclr, rclr);
  }

  delete clrRange;

  delete readData;

  for (uint32 i=0; i<=numLibs; i++)
    delete out[i];
  delete [] out;

  gkpStore->gkStore_close();

  exit(0);
}


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
 *    src/stores/gatekeeperDumpFASTQ.C
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
#include "sqStore.H"
#include "files.H"
#include "strings.H"
#include "sequence.H"

#include "clearRangeFile.H"

//  Write sequence in multiple formats.  This used to write to four fastq files, the .1, .2, .paired and .unmated.
//  It's left around for future expansion to .fastq and .bax.h5.
//
class libOutput {
public:
  libOutput(char const *outPrefix, char const *outSuffix, char const *libName = NULL) {
    strcpy(_p, outPrefix);

    if (outSuffix[0])
      snprintf(_s, FILENAME_MAX, ".%s", outSuffix);
    else
      _s[0] = 0;

    if (libName)
      strncpy(_n, libName, FILENAME_MAX-1);
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
      snprintf(N, FILENAME_MAX, "%s.%s.fastq%s", _p, _n, _s);
    else
      snprintf(N, FILENAME_MAX, "%s.fastq%s", _p, _s);

    if ((_p[0] == '-') && (_p[1] == 0)) {
      snprintf(N, FILENAME_MAX, "(stdout)");
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
      snprintf(N, FILENAME_MAX, "%s.%s.fasta%s", _p, _n, _s);
    else
      snprintf(N, FILENAME_MAX, "%s.fasta%s", _p, _s);

    if ((_p[0] == '-') && (_p[1] == 0)) {
      snprintf(N, FILENAME_MAX, "(stdout)");
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
  char            *seqStoreName      = NULL;
  uint32           seqStorePart      = UINT32_MAX;

  char            *outPrefix         = NULL;
  char            *outSuffix         = NULL;

  char            *clrName           = NULL;

  uint32           libToDump         = 0;

  uint32           bgnID             = 1;
  uint32           endID             = UINT32_MAX;

  bool             dumpRaw           = false;
  bool             dumpCorrected     = false;
  bool             dumpTrimmed       = false;

  bool             dumpAllReads      = false;
  bool             dumpAllBases      = false;
  bool             dumpOnlyDeleted   = false;

  bool             dumpFASTQ         = true;
  bool             dumpFASTA         = false;

  bool             withLibName       = true;
  bool             withReadName      = true;

  bool             asReverse         = false;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStoreName = argv[++arg];

      if ((arg+1 < argc) && (argv[arg+1][0] != '-'))
        seqStorePart = atoi(argv[++arg]);

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
      decodeRange(argv[++arg], bgnID, endID);

    } else if (strcmp(argv[arg], "-raw") == 0) {
      dumpRaw = true;

    } else if (strcmp(argv[arg], "-corrected") == 0) {
      dumpCorrected = true;

    } else if (strcmp(argv[arg], "-trimmed") == 0) {
      dumpTrimmed = true;

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

    } else if (strcmp(argv[arg], "-noreadname") == 0) {
      withReadName    = false;

    } else if (strcmp(argv[arg], "-reverse") == 0) {
      asReverse       = true;


    } else {
      err++;
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (seqStoreName == NULL)
    err++;
  if (outPrefix == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s [...] -o fastq-prefix -g seqStore\n", argv[0]);
    fprintf(stderr, "  -S seqStore\n");
    fprintf(stderr, "  -o fastq-prefix     write files fastq-prefix.(libname).fastq, ...\n");
    fprintf(stderr, "                      if fastq-prefix is '-', all sequences output to stdout\n");
    fprintf(stderr, "                      if fastq-prefix ends in .gz, .bz2 or .xz, output is compressed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -fastq              output is FASTQ format (with extension .fastq, default)\n");
    fprintf(stderr, "  -fasta              output is FASTA format (with extension .fasta)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nolibname          don't include the library name in the output file name\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -noreadname         don't include the read name in the sequence header.  header will be:\n");
    fprintf(stderr, "                        '>original-name id=<seqID> clr=<bgn>,<end>   with names\n");
    fprintf(stderr, "                        '>read<seqID> clr=<bgn>,<end>                without names\n");
    fprintf(stderr, "\n");
    fprintf(stderr, " -reverse             Dump the reverse-complement of the read.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -l libToDump        output only read in library number libToDump\n");
    fprintf(stderr, "  -r id[-id]          output only the single read 'id', or the specified range of ids\n");
    fprintf(stderr, "\n");
    fprintf(stderr, " The default is to dump the latest version of each read.  You can force it to dump:\n");
    fprintf(stderr, "  -raw                Dump only raw reads.\n");
    fprintf(stderr, "  -corrected          Dump only corrected reads.\n");
    fprintf(stderr, "  -trimmed            Dump only trimmed reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, " An Overlap Based Trimming clear range file can be supplied (BUT NOT WIDELY TESTED):\n");
    fprintf(stderr, "  -c clearFile        clear range file from OBT modules\n");
    fprintf(stderr, "  -allreads           if a clear range file, lower case mask the deleted reads\n");
    fprintf(stderr, "  -allbases           if a clear range file, lower case mask the non-clear bases\n");
    fprintf(stderr, "  -onlydeleted        if a clear range file, only output deleted reads (the entire read)\n");

    if (seqStoreName == NULL)
      fprintf(stderr, "ERROR: no seqStore (-S) supplied.\n");
    if (outPrefix == NULL)
      fprintf(stderr, "ERROR: no output prefix (-o) supplied.\n");

    exit(1);
  }

  sqStore        *seqStore  = sqStore::sqStore_open(seqStoreName, sqStore_readOnly, seqStorePart);
  uint32          numReads  = seqStore->sqStore_getNumReads();
  uint32          numLibs   = seqStore->sqStore_getNumLibraries();

  clearRangeFile *clrRange  = (clrName == NULL) ? NULL : new clearRangeFile(clrName, seqStore);

  if (bgnID < 1)
    bgnID = 1;

  if (numReads < endID)
    endID = numReads;

  if (endID < bgnID)
    fprintf(stderr, "No reads to dump; reversed ranges make no sense: bgn=" F_U32 " end=" F_U32 "??\n", bgnID, endID), exit(1);

  if ((dumpRaw == true) && (seqStore->sqStore_getNumRawReads() == 0)) {
    fprintf(stderr, "WARNING:\n");
    fprintf(stderr, "WARNING:  No raw reads in store.\n");
    fprintf(stderr, "WARNING:\n");
  }

  if ((dumpCorrected == true) && (seqStore->sqStore_getNumCorrectedReads() == 0)) {
    fprintf(stderr, "WARNING:\n");
    fprintf(stderr, "WARNING:  No corrected reads in store.\n");
    fprintf(stderr, "WARNING:\n");
  }

  if ((dumpTrimmed == true) && (seqStore->sqStore_getNumTrimmedReads() == 0)) {
    fprintf(stderr, "WARNING:\n");
    fprintf(stderr, "WARNING:  No trimmed reads in store.\n");
    fprintf(stderr, "WARNING:\n");
  }


  if (seqStorePart == UINT32_MAX)
    fprintf(stderr, "Dumping reads from %u to %u (inclusive).\n", bgnID, endID);
  else
    fprintf(stderr, "Dumping reads from %u to %u (inclusive) from partition %u.\n", bgnID, endID, seqStorePart);


  libOutput   **out = new libOutput * [numLibs + 1];

  //  Allocate outputs.  If withLibName == false, all reads will artificially be in lib zero, the
  //  other files won't ever be created.  Otherwise, the zeroth file won't ever be created.

  out[0] = new libOutput(outPrefix, outSuffix, NULL);

  for (uint32 i=1; i<=numLibs; i++)
    out[i] = new libOutput(outPrefix, outSuffix, seqStore->sqStore_getLibrary(i)->sqLibrary_libraryName());

  //  Grab a new readData, and iterate through reads to dump.

  sqReadData   *readData  = new sqReadData;

  char         *qltString = new char [AS_MAX_READLEN + 1];

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    sqRead      *read   = seqStore->sqStore_getRead(rid);

    if ((read == NULL) ||
        (seqStore->sqStore_readInPartition(rid) == false))
      continue;

    uint32       libID  = (withLibName == false) ? 0 : read->sqRead_libraryID();

    uint32       flen   = read->sqRead_sequenceLength();

    if (dumpRaw == true)
      flen = read->sqRead_sequenceLength(sqRead_raw);

    if (dumpCorrected == true)
      flen = read->sqRead_sequenceLength(sqRead_corrected);

    if (dumpTrimmed == true)
      flen = read->sqRead_sequenceLength(sqRead_trimmed);

    uint32       lclr   = 0;
    uint32       rclr   = flen;
    bool         ignore = false;

    //  If a clear range file is supplied, grab the clear range.  If it hasn't been set, the default
    //  is the entire read.

    if (clrRange) {
      lclr   = clrRange->bgn(rid);
      rclr   = clrRange->end(rid);
      ignore = clrRange->isDeleted(rid);
    }

    //  Abort if we're not dumping anything from this read

    if (((libToDump != 0) && (libID == libToDump)) ||            //   - not in a library we care about
        ((dumpAllReads == false) && (ignore == true)) ||         //   - deleted, and not dumping all reads
        ((dumpOnlyDeleted == true) && (ignore == false)))        //   - not deleted, but only reporting deleted reads
      continue;

    //  If the read length is zero, then the read has been removed from this set.

    if ((dumpAllReads == false) && (flen == 0))
      continue;

    //  And if we're told to ignore the read, and here, then the read was deleted and we're printing
    //  all reads.  Reset the clear range to the whole read, the clear range is invalid.

    if (ignore) {
      lclr = 0;
      rclr = flen;
    }

    uint32  clen = rclr - lclr;

    //  Grab the _latest_ sequence and quality.

    seqStore->sqStore_loadReadData(read, readData);

    char   *name = readData->sqReadData_getName();

    char   *seq  = readData->sqReadData_getSequence();
    uint8  *qlt8 = readData->sqReadData_getQualities();
    char   *qlt  = qltString;

    //  Grab the specified sequence and quality, if specified.

    if (dumpRaw == true) {
      seq  = readData->sqReadData_getRawSequence();
      qlt8 = readData->sqReadData_getRawQualities();
    }

    if (dumpCorrected == true) {
      seq  = readData->sqReadData_getCorrectedSequence();
      qlt8 = readData->sqReadData_getCorrectedQualities();
    }

    if (dumpTrimmed == true) {
      seq  = readData->sqReadData_getTrimmedSequence();
      qlt8 = readData->sqReadData_getTrimmedQualities();
    }

    //  Soft mask not-clear bases.

    if (dumpAllBases == true) {
      for (uint32 i=0; i<lclr; i++)
        seq[i] += (seq[i] >= 'A') ? 'a' - 'A' : 0;

      for (uint32 i=lclr; i<rclr; i++)
        seq[i] += (seq[i] >= 'A') ? 0 : 'A' - 'a';

      for (uint32 i=rclr; i<flen; i++)
        seq[i] += (seq[i] >= 'A') ? 'a' - 'A' : 0;

      lclr = 0;
      rclr = flen;
    }

    //  Create the QV string.

    for (uint32 i=0; i<flen; i++)
      qlt[i] = '!' + qlt8[i];

    //  Chop off the ends we're not printing.

    seq += lclr;
    qlt += lclr;

    seq[clen] = 0;
    qlt[clen] = 0;

    //  And maybe reverse complement it.

    if (asReverse)
      reverseComplement(seq, qlt, clen);

    //  Print the read.

    if (dumpFASTA)  //  Dear GCC:  I'm NOT ambiguous
      if ((withReadName == true) && (name != NULL))
        AS_UTL_writeFastA(out[libID]->getFASTA(), seq, clen, 0,
                        ">%s id=" F_U32 " clr=" F_U32 "," F_U32 "\n",
                        name, rid, lclr, rclr);
      else
        AS_UTL_writeFastA(out[libID]->getFASTA(), seq, clen, 0,
                          ">read" F_U32 " clr=" F_U32 "," F_U32 "\n",
                          rid, lclr, rclr);

    if (dumpFASTQ)  //  Dear GCC:  I'm NOT ambiguous
      if ((withReadName == true) && (name != NULL))
        AS_UTL_writeFastQ(out[libID]->getFASTQ(), seq, clen, qlt, clen,
                          "@%s id=" F_U32 " clr=" F_U32 "," F_U32 "\n",
                          name,
                          rid, lclr, rclr);
      else
        AS_UTL_writeFastQ(out[libID]->getFASTQ(), seq, clen, qlt, clen,
                          "@read" F_U32 " clr=" F_U32 "," F_U32 "\n",
                          rid, lclr, rclr);
  }

  delete clrRange;

  delete readData;

  delete [] qltString;

  for (uint32 i=0; i<=numLibs; i++)
    delete out[i];
  delete [] out;

  seqStore->sqStore_close();

  exit(0);
}

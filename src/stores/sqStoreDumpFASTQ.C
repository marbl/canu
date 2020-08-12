
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"
#include "sqStore.H"
#include "files.H"
#include "strings.H"
#include "sequence.H"

#include "clearRangeFile.H"


//  Write sequence in multiple formats.  This used to write to four fastq
//  files, the .1, .2, .paired and .unmated.  It's left around for future
//  expansion to .fastq and .bax.h5.
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

  char            *outPrefix         = NULL;
  char            *outSuffix         = NULL;

  uint32           libToDump         = 0;

  uint32           bgnID             = 1;
  uint32           endID             = UINT32_MAX;

  sqRead_which     readType          = sqRead_unset;

#if 0
  char            *clrName           = NULL;
  bool             dumpAllReads      = false;
  bool             dumpAllBases      = false;
  bool             dumpOnlyDeleted   = false;
#endif

  bool             dumpFASTQ         = true;

  bool             withLibName       = true;
  bool             withReadName      = true;

  bool             asReverse         = false;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];
      outSuffix = scanPrefix(outPrefix);


    } else if (strcmp(argv[arg], "-l") == 0) {
      libToDump = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], bgnID, endID);


    } else if (strcmp(argv[arg], "-raw") == 0) {
      readType &= ~sqRead_corrected;
      readType |=  sqRead_raw;

    } else if (strcmp(argv[arg], "-corrected") == 0) {
      readType &= ~sqRead_raw;
      readType |=  sqRead_corrected;

    } else if (strcmp(argv[arg], "-trimmed") == 0) {
      readType |=  sqRead_trimmed;

    } else if (strcmp(argv[arg], "-compressed") == 0) {
      readType |=  sqRead_compressed;

    } else if (strcmp(argv[arg], "-normal") == 0) {
      readType |=  sqRead_normal;

#if 0
    } else if (strcmp(argv[arg], "-c") == 0) {
      clrName = argv[++arg];

    } else if (strcmp(argv[arg], "-allreads") == 0) {
      dumpAllReads    = true;

    } else if (strcmp(argv[arg], "-allbases") == 0) {
      dumpAllBases    = true;

    } else if (strcmp(argv[arg], "-onlydeleted") == 0) {
      dumpOnlyDeleted = true;
      dumpAllReads    = true;  //  Otherwise we won't report the deleted reads!
#endif


    } else if (strcmp(argv[arg], "-fastq") == 0) {
      dumpFASTQ       = true;

    } else if (strcmp(argv[arg], "-fasta") == 0) {
      dumpFASTQ       = false;


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
    fprintf(stderr, "usage: %s [...] -o out-prefix -g seqStore\n", argv[0]);
    fprintf(stderr, "  -S seqStore\n");
    fprintf(stderr, "  -o out-prefix       write files out-prefix.(libname).fastq, ...\n");
    fprintf(stderr, "                      if out-prefix is '-', all sequences output to stdout\n");
    fprintf(stderr, "                      if out-prefix ends in .gz, .bz2 or .xz, output is compressed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -fastq              output is FASTQ format (with extension .fastq, default)\n");
    fprintf(stderr, "                      (note that QVs are not stored, and are invalid)\n");
    fprintf(stderr, "\n");
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
    fprintf(stderr, "  -l id               output only read in library number 'id'\n");
    fprintf(stderr, "  -r id[-id]          output only the single read 'id', or the specified range of ids\n");
    fprintf(stderr, "\n");
    fprintf(stderr, " The default is to dump the latest version of each read.  You can force it to dump:\n");
    fprintf(stderr, "  -raw                Dump raw reads.\n");
    fprintf(stderr, "  -corrected          Dump corrected reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -trimmed            Dump the trimmed version of the raw/corrected read.\n");
    fprintf(stderr, "  -compressed         Dump the compressed version of the raw/corrected read.\n");
    fprintf(stderr, "  -normal             Dump the uncompressed version of the raw/corrected read.\n");
    fprintf(stderr, "                        (for stores that are by default compressing homopolymers)\n");
    fprintf(stderr, "\n");
#if 0
    fprintf(stderr, " An Overlap Based Trimming clear range file can be supplied (BUT NOT WIDELY TESTED):\n");
    fprintf(stderr, "  -c clearFile        clear range file from OBT modules\n");
    fprintf(stderr, "  -allreads           if a clear range file, lower case mask the deleted reads\n");
    fprintf(stderr, "  -allbases           if a clear range file, lower case mask the non-clear bases\n");
    fprintf(stderr, "  -onlydeleted        if a clear range file, only output deleted reads (the entire read)\n");
#endif

    if (seqStoreName == NULL)
      fprintf(stderr, "ERROR: no seqStore (-S) supplied.\n");
    if (outPrefix == NULL)
      fprintf(stderr, "ERROR: no output prefix (-o) supplied.\n");

    exit(1);
  }

  sqRead_setDefaultVersion(readType);

  sqStore        *seqStore  = new sqStore(seqStoreName, sqStore_readOnly);
  uint32          numReads  = seqStore->sqStore_lastReadID();
  uint32          numLibs   = seqStore->sqStore_lastLibraryID();

  fprintf(stderr, "Opened seqStore '%s' for '%s' reads.\n", seqStoreName, sqRead_getDefaultVersion());

#if 0
  clearRangeFile *clrRange  = (clrName == NULL) ? NULL : new clearRangeFile(clrName, seqStore);
#endif

  if (bgnID < 1)
    bgnID = 1;

  if (numReads < endID)
    endID = numReads;

  if (endID < bgnID)
    fprintf(stderr, "No reads to dump; reversed ranges make no sense: bgn=" F_U32 " end=" F_U32 "??\n", bgnID, endID), exit(1);


  fprintf(stderr, "Dumping %s reads from %u to %u (inclusive).\n",
          toString(sqRead_defaultVersion), bgnID, endID);


  //  Allocate outputs.  If withLibName == false, all reads will artificially be in lib zero, the
  //  other files won't ever be created.  Otherwise, the zeroth file won't ever be created.

  libOutput   **out = new libOutput * [numLibs + 1];

  out[0] = new libOutput(outPrefix, outSuffix, NULL);

  for (uint32 i=1; i<=numLibs; i++)
    out[i] = new libOutput(outPrefix, outSuffix, seqStore->sqStore_getLibrary(i)->sqLibrary_libraryName());



  sqRead       *read      = new sqRead();
  char         *readName  = new char [1024];
  char         *seq       = new char [AS_MAX_READLEN + 1];
  char         *qlt       = new char [AS_MAX_READLEN + 1];

  for (uint32 rid=bgnID; rid<=endID; rid++) {
    uint32       libID  = seqStore->sqStore_getLibraryIDForRead(rid);

    //  Skip the read if it isn't in our library.
    if ((libToDump != 0) && (libID != libToDump))
      continue;

    //  Skip the read if it is ignored or not valid.
    if ((seqStore->sqStore_isValidRead(rid)   == false) ||
        (seqStore->sqStore_isIgnoredRead(rid) == true))
      continue;

    //  Dump the read.  The store does all trimming and compressing, we just
    //  need to (maybe) reverse-complement it, and print it.

    seqStore->sqStore_getRead(rid, read);       //  Load the sequence data.

    uint32   seqLen = seqStore->sqStore_getReadLength(rid);
    char    *S      = read->sqRead_sequence();

    for (uint32 i=0; i<seqLen; i++) {           //  Create a QV string.
      seq[i] = S[i];
      qlt[i] = '!';
    }
    seq[seqLen] = 0;
    qlt[seqLen] = 0;

    if (asReverse)                              //  Reverse complement?
      reverseComplement(seq, qlt, seqLen);

    //  Print the read.

    uint32  outid = (withLibName == false) ? 0 : libID;

    if (withReadName)
      snprintf(readName, 1024, "%s id=" F_U32, read->sqRead_name(), rid);
    else
      snprintf(readName, 1024, "read" F_U32, rid);

    if (dumpFASTQ)
      AS_UTL_writeFastQ(out[outid]->getFASTQ(), seq, seqLen, qlt, seqLen, "@%s\n", readName);
    else
      AS_UTL_writeFastA(out[outid]->getFASTA(), seq, seqLen, 0,           ">%s\n", readName);
  }

  delete    read;
  delete [] readName;
  delete [] qlt;

  for (uint32 i=0; i<=numLibs; i++)
    delete out[i];
  delete [] out;

  delete seqStore;

  exit(0);
}

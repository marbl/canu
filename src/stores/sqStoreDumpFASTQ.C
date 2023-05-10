
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


class params {
public:
  params() {
  }
  ~params() {
    delete    seqStore;
    delete    read;
    delete [] readName;
    delete [] seq;
    delete [] qlt;
  }

  char            *seqStoreName  = nullptr;
  sqStore         *seqStore      = nullptr;
  uint32           numReads      = 0;
  uint32           numLibs       = 0;

  sqRead          *read          = new sqRead();
  char            *readName      = new char [1024];
  char            *seq           = new char [AS_MAX_READLEN + 1];
  char            *qlt           = new char [AS_MAX_READLEN + 1];

  char            *outPrefix     = nullptr;
  char            *outSuffix     = nullptr;

  uint32           libToDump     = 0;

  uint32           bgnID         = 1;
  uint32           endID         = uint32max;
  std::set<uint32> setIDs;

  sqRead_which     readType      = sqRead_unset;

  bool             dumpFASTQ     = true;

  bool             withLibName   = true;
  bool             withReadName  = true;

  bool             asReverse     = false;
};




//  Write sequence in multiple formats.  This used to write to four fastq
//  files, the .1, .2, .paired and .unmated.  It's left around for future
//  expansion to .fastq and .bax.h5.
//
class libOutput {
public:
  libOutput(char const *outPrefix, char const *outSuffix, char const *libName = nullptr) {
    strcpy(_p, outPrefix);

    if (outSuffix[0])
      snprintf(_s, FILENAME_MAX, ".%s", outSuffix);
    else
      _s[0] = 0;

    if (libName)
      strncpy(_n, libName, FILENAME_MAX-1);
    else
      _n[0] = 0;

    _WRITER = nullptr;
    _FASTA  = nullptr;
    _FASTQ  = nullptr;
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



//  If arg has commas, or if bgnID or endID have been set, store the IDs to dump in setIDs.
void
decodeRange(char const *arg, uint32 &bgnID, uint32 &endID, std::set<uint32> &setIDs) {
  char const *comma = strchr(arg, ',');

  if ((bgnID != 1) || (endID != uint32max)) {
    for (uint32 ii=bgnID; ii<endID; ii++)
      setIDs.insert(ii);

    bgnID = 1;
    endID = uint32max;
  }

  if ((comma != nullptr) ||
      (setIDs.size() > 0))
    decodeRange(arg, setIDs);
  else
    decodeRange(arg, bgnID, endID);
}



void
dumpRead(params &p, libOutput **out, uint32 rid) {
  uint32       libID  = p.seqStore->sqStore_getLibraryIDForRead(rid);

  //  Skip the read if it isn't in our library.
  if ((p.libToDump != 0) && (libID != p.libToDump))
    return;

  //  Skip the read if it is ignored or not valid.
  if ((p.seqStore->sqStore_isValidRead(rid)   == false) ||
      (p.seqStore->sqStore_isIgnoredRead(rid) == true))
    return;

  //  Skip the read if it's out of range.
  if (p.numReads < rid)
    return;

  //  Dump the read.  The store does all trimming and compressing, we just
  //  need to (maybe) reverse-complement it, and print it.

  p.seqStore->sqStore_getRead(rid, p.read);     //  Load the sequence data.

  uint32   seqLen = p.seqStore->sqStore_getReadLength(rid);
  char    *S      = p.read->sqRead_sequence();

  for (uint32 i=0; i<seqLen; i++) {             //  Create a QV string.
    p.seq[i] = S[i];
    p.qlt[i] = '!';
  }
  p.seq[seqLen] = 0;
  p.qlt[seqLen] = 0;

  if (p.asReverse)                              //  Reverse complement?
    reverseComplement(p.seq, p.qlt, seqLen);

  //  Print the read.

  uint32  outid = (p.withLibName == false) ? 0 : libID;

  if (p.withReadName)
    snprintf(p.readName, 1024, "%s id=" F_U32, p.read->sqRead_name(), rid);
  else
    snprintf(p.readName, 1024, "read" F_U32, rid);

  if (p.dumpFASTQ)
    outputFASTQ(out[outid]->getFASTQ(), p.seq, p.qlt, seqLen,    "%s", p.readName);
  else
    outputFASTA(out[outid]->getFASTA(), p.seq,        seqLen, 0, "%s", p.readName);
}



int
main(int argc, char **argv) {
  params           p;

  argc = AS_configure(argc, argv, 1);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      p.seqStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      p.outPrefix = argv[++arg];
      p.outSuffix = scanPrefix(p.outPrefix);

    } else if (strcmp(argv[arg], "-l") == 0) {
      p.libToDump = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], p.bgnID, p.endID, p.setIDs);


    } else if (strcmp(argv[arg], "-raw") == 0) {
      p.readType &= ~sqRead_corrected;
      p.readType |=  sqRead_raw;
    } else if (strcmp(argv[arg], "-corrected") == 0) {
      p.readType &= ~sqRead_raw;
      p.readType |=  sqRead_corrected;
    } else if (strcmp(argv[arg], "-trimmed") == 0) {
      p.readType |=  sqRead_trimmed;
    } else if (strcmp(argv[arg], "-compressed") == 0) {
      p.readType |=  sqRead_compressed;
    } else if (strcmp(argv[arg], "-normal") == 0) {
      p.readType |=  sqRead_normal;


    } else if (strcmp(argv[arg], "-fastq") == 0) {
      p.dumpFASTQ       = true;
    } else if (strcmp(argv[arg], "-fasta") == 0) {
      p.dumpFASTQ       = false;


    } else if (strcmp(argv[arg], "-nolibname") == 0) {
      p.withLibName     = false;

    } else if (strcmp(argv[arg], "-noreadname") == 0) {
      p.withReadName    = false;

    } else if (strcmp(argv[arg], "-reverse") == 0) {
      p.asReverse       = true;


    } else {
      err++;
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (p.seqStoreName == nullptr)
    err++;
  if (p.outPrefix == nullptr)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -S seqStore -o out-prefix [...]\n", argv[0]);
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

    if (p.seqStoreName == nullptr)
      fprintf(stderr, "ERROR: no seqStore (-S) supplied.\n");
    if (p.outPrefix == nullptr)
      fprintf(stderr, "ERROR: no output prefix (-o) supplied.\n");

    exit(1);
  }

  //  Open the store.

  sqRead_setDefaultVersion(p.readType);

  p.seqStore  = new sqStore(p.seqStoreName, sqStore_readOnly);
  p.numReads  = p.seqStore->sqStore_lastReadID();
  p.numLibs   = p.seqStore->sqStore_lastLibraryID();

  fprintf(stderr, "Opened seqStore '%s' for '%s' reads.\n", p.seqStoreName, sqRead_getDefaultVersion());

  //  Check ranges.

  if (p.setIDs.size() > 0) {
    uint32  nInv = 0;
    uint32  minID = uint32max;
    uint32  maxID = 0;

    for (auto id=p.setIDs.begin(); id != p.setIDs.end(); id++) {
      minID = std::min(*id, minID);
      maxID = std::max(*id, maxID);

      if ((*id < 1) ||
          (p.numReads < *id)) {
        fprintf(stderr, "ERROR: id %u is not a valid read (range %u-%u).\n", *id, 1, p.numReads);
        nInv++;
      }

      if (nInv > 0)
        fprintf(stderr, "Invalid read IDs detected.  No reads dumped.\n"), exit(1);
    }

    fprintf(stderr, "Dumping %lu %s reads between %u and %u (inclusive).\n",
            p.setIDs.size(), toString(sqRead_defaultVersion), minID, maxID);
  }
  else {
    if (p.bgnID < 1)           p.bgnID = 1;
    if (p.numReads < p.endID)  p.endID = p.numReads;

    if (p.endID < p.bgnID)
      fprintf(stderr, "No reads to dump; reversed ranges make no sense: bgn=" F_U32 " end=" F_U32 "??\n", p.bgnID, p.endID), exit(1);

    fprintf(stderr, "Dumping %s reads from %u to %u (inclusive).\n",
            toString(sqRead_defaultVersion), p.bgnID, p.endID);
  }

  //  Allocate outputs.  If withLibName == false, all reads will artificially be in lib zero, the
  //  other files won't ever be created.  Otherwise, the zeroth file won't ever be created.

  libOutput   **out = new libOutput * [p.numLibs + 1];

  out[0] = new libOutput(p.outPrefix, p.outSuffix, nullptr);

  for (uint32 i=1; i<=p.numLibs; i++)
    out[i] = new libOutput(p.outPrefix, p.outSuffix, p.seqStore->sqStore_getLibrary(i)->sqLibrary_libraryName());

  //  Dump!

  if (p.setIDs.size() > 0)
    for (auto id=p.setIDs.begin(); id != p.setIDs.end(); id++)
      dumpRead(p, out, *id);

  else
    for (uint32 rid=p.bgnID; rid<=p.endID; rid++)
      dumpRead(p, out, rid);

  //  Cleanup.

  for (uint32 i=0; i<=p.numLibs; i++)
    delete out[i];
  delete [] out;

  exit(0);
}

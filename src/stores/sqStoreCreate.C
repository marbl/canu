
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

#include "system.H"
#include "files.H"
#include "strings.H"
#include "math.H"

#include "sqStore.H"

#include <vector>
#include <algorithm>


//  A list of the letters that we accept in sequences.
uint32  validSeq[256] = {0};


//  A list of files that should be loaded into one library.
class seqLib {
public:
  seqLib(char *name, sqLibrary_tech tech, sqRead_which stat) {
    _name = name;
    _tech = tech;
    _stat = stat;
  };
  ~seqLib() {
  }

  void             addFile(char *name) {
    _files.push_back(name);
  };

public:
  char                *_name;
  sqLibrary_tech       _tech;
  sqRead_which         _stat;

  std::vector<char *>  _files;
};



class loadStats {
public:
  loadStats() {
    nINVALID = nSHORT = nLONG = nLOADED = nWARNINGS = 0;
    bINVALID = bSHORT = bLONG = bLOADED = 0;
  };
  ~loadStats() {
  };

  void         import(loadStats &that) {
    nINVALID  += that.nINVALID;   bINVALID  += that.bINVALID;
    nSHORT    += that.nSHORT;     bSHORT    += that.bSHORT;
    nLONG     += that.nLONG;      bLONG     += that.bLONG;
    nLOADED   += that.nLOADED;    bLOADED   += that.bLOADED;
    nWARNINGS += that.nWARNINGS;
  };

#define PERC(x,t)  (t > 0) ? (100.0 * x / t) : (0.0)

  void         displayTableHeader(FILE *F) {
    fprintf(stderr, "\n");
    fprintf(stderr, "               reads               bases\n");
    fprintf(stderr, "---------- --------- ------ ------------ ------\n");
  };

  void         displayTable(FILE *F, char *fileName = NULL) {
    uint32 nTotal = nINVALID + nSHORT + nLONG + nLOADED;
    uint64 bTotal = bINVALID + bSHORT + bLONG + bLOADED;

    if (fileName)
      fprintf(F, "%-10s %9" F_U32P " %5.1f%% %12" F_U64P " %5.1f%%  %s\n",
              "Loaded",
              nLOADED, PERC(nLOADED, nTotal),
              bLOADED, PERC(bLOADED, bTotal),
              fileName);
    else
      fprintf(F, "%-10s %9" F_U32P " %5.1f%% %12" F_U64P " %5.1f%%\n",
              "Loaded",
              nLOADED, PERC(nLOADED, nTotal),
              bLOADED, PERC(bLOADED, bTotal));


    if (nSHORT > 0)
      fprintf(F, "%-10s %9" F_U32P " %5.1f%% %12" F_U64P " %5.1f%%\n",
              "Short",
              nSHORT, PERC(nSHORT, nTotal),
              bSHORT, PERC(bSHORT, bTotal));

    if (nLONG > 0)
      fprintf(F, "%-10s %9" F_U32P " %5.1f%% %12" F_U64P " %5.1f%%\n",
              "Long",
              nLONG, PERC(nLONG, nTotal),
              bLONG, PERC(bLONG, bTotal));

    if (nINVALID > 0)
      fprintf(F, "%-10s %9" F_U32P " %5.1f%% %12" F_U64P " %5.1f%%\n",
              "Invalid",
              nINVALID, PERC(nINVALID, nTotal),
              bINVALID, PERC(bINVALID, bTotal));

    if (nWARNINGS > 0)
      fprintf(F, "%-10s %9" F_U32P "\n",
              "Warnings",
              nWARNINGS);

    fprintf(F, "\n");
  };

  uint32       nINVALID, nSHORT, nLONG, nLOADED, nWARNINGS;
  uint64       bINVALID, bSHORT, bLONG, bLOADED;
};





uint64
trimBgn(dnaSeq &sq, uint64 bgn, uint64 end) {
  char const *bases = sq.bases();

  while ((bgn <= end) && ((bases[bgn] == 'N') ||
                          (bases[bgn] == 'n')))
    bgn++;

  return(bgn);
}



uint64
trimEnd(dnaSeq &sq, uint64 bgn, uint64 end) {
  char const *bases = sq.bases();

  end--;  //  So end is now the last base in the sequence.

  while ((bgn <= end) && ((bases[end] == 'N') ||
                          (bases[end] == 'n')))
    end--;

  end++;  //  So now end is the C-style end.

  return(end);
}



uint32
checkInvalid(dnaSeq &sq, uint64 bgn, uint64 end) {
  uint32       invalid = 0;
  char const  *bases   = sq.bases();

  //  The conversion from U to T that was here is now in sqReadDataWriter
  //  setRawBases().
  for (uint32 ii=bgn; ii<end; ii++)         //  We cast to uint8 because char is signed so values >128 get interpreted as negative
    if (validSeq[uint8(bases[ii])] == 0)    //  and access invalid memory locations in the array.
      invalid++;

  return(invalid);
}



void
loadReads(sqStore          *seqStore,
          sqLibrary        *seqLibrary,
          sqRead_which      readStat,
          uint32            minReadLength,
          bool              homopolyCompress,
          FILE             *nameMap,
          FILE             *errorLog,
          char             *fileName,
          loadStats        &stats) {

  //fprintf(stderr, "  %s:\n", fileName);

  loadStats    filestats;

  dnaSeqFile  *SF = new dnaSeqFile(fileName);
  dnaSeq       sq;

  while (SF->loadSequence(sq) == true) {

    //  Check for and log parsing errors.

    if (sq.wasError() == true) {
      fprintf(errorLog, "error reading sequence at/before '%s' in file '%s'.\n",
              sq.ident(), fileName);
      filestats.nWARNINGS += 1;
    }

    if (sq.wasReSync() == true) {
      fprintf(errorLog, "lost sync reading before sequence '%s' in file '%s'.\n",
              sq.ident(), fileName);
      filestats.nWARNINGS += 1;
    }

    //  Trim Ns from the ends of the sequence.
    uint64  bgn = trimBgn(sq, 0,   sq.length());
    uint64  end = trimEnd(sq, bgn, sq.length());

    if ((bgn > 0) && (end < sq.length()))
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - trimmed " F_U64 " non-ACGT bases from the 5' and " F_U64 " non-ACGT bases from the 3' end.\n",
              sq.ident(), sq.length(), fileName, bgn, sq.length() - end);

    else if (bgn > 0)
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - trimmed " F_U64 " non-ACGT bases from the 5' end.\n",
              sq.ident(), sq.length(), fileName, bgn);

    else if (end < sq.length())
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - trimmed " F_U64 " non-ACGT bases from the 3' end.\n",
              sq.ident(), sq.length(), fileName, sq.length() - end);

    //  Check for invalid bases.

    uint32 invalid = checkInvalid(sq, bgn, end);

    if (invalid > 0) {
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - contains %u invalid letters, skipping.\n",
              sq.ident(), sq.length(), fileName, invalid);

      filestats.nINVALID += 1;
      filestats.bINVALID += sq.length();

      continue;
    }

    //  Create a writer for the read data and load bases.

    sqReadDataWriter *rdw  = seqStore->sqStore_createEmptyRead(seqLibrary, sq.ident());

    if (readStat & sqRead_raw)
      rdw->sqReadDataWriter_setRawBases(sq.bases() + bgn, end - bgn);
    else
      rdw->sqReadDataWriter_setCorrectedBases(sq.bases() + bgn, end - bgn);

    //  Get the (homopolymer compressed) length of the sequence we just loaded.

    uint32 rLen = (readStat & sqRead_raw) ? rdw->sqReadDataWriter_getRawLength(homopolyCompress)
                                          : rdw->sqReadDataWriter_getCorrectedLength(homopolyCompress);

    //  Drop any sequences that are short...
    if      (rLen < minReadLength) {
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - too short, skipping.\n",
              sq.ident(), sq.length(), fileName);

      filestats.nSHORT += 1;
      filestats.bSHORT += sq.length();
    }

    //  ...or too long...
    else if (rLen > AS_MAX_READLEN - 2) {
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - too long, skipping.\n",
              sq.ident(), sq.length(), fileName);

      filestats.nLONG += 1;
      filestats.bLONG += sq.length();
    }

    //  ...and add the ones that are Just Right!
    //
    //  With the read added to the store, set trim points.  Presently,
    //  trimming only occurs on corrected reads, but later we need to allow
    //  trimmed raw reads.
    //
    //  Finally, update the nameMap and save some silly statistics.
    else {
      seqStore->sqStore_addRead(rdw);

      if (readStat & sqRead_trimmed) {
        uint32      rid  = seqStore->sqStore_lastReadID();
        sqReadSeq  *nseq = seqStore->sqStore_getReadSeq(rid, sqRead_corrected);
        sqReadSeq  *cseq = seqStore->sqStore_getReadSeq(rid, sqRead_corrected | sqRead_compressed);

        nseq->sqReadSeq_setAllClear();
        cseq->sqReadSeq_setAllClear();
      }

      fprintf(nameMap, F_U32"\t%s%s%s\n",
              seqStore->sqStore_lastReadID(),
              sq.ident(),
              (sq.flags()[0] == 0) ? "" : " ",
              sq.flags());

      filestats.nLOADED += 1;
      filestats.bLOADED += end - bgn;
    }

    //  All done with this read.  Delete the writer and continue.
    delete rdw;
  }

  delete SF;

  //  Write status to the screen
  filestats.displayTable(stderr, fileName);

  //  Add the just loaded numbers to the global numbers
  stats.import(filestats);
};





bool
createStore(const char            *seqStoreName,
            std::vector<seqLib>   &libraries,
            uint32                 minReadLength,
            bool                   homopolyCompress) {

  sqStore     *seqStore     = new sqStore(seqStoreName, sqStore_create);   //  sqStore_extend MIGHT work
  sqRead      *seqRead      = NULL;
  sqLibrary   *seqLibrary   = NULL;

  uint32       inLineLen    = 1024;
  char         inLine[1024] = { 0 };

  FILE        *errorLog = merylutil::openOutputFile(seqStoreName, '/', "errorLog");
  FILE        *nameMap  = merylutil::openOutputFile(seqStoreName, '/', "readNames.txt");

  loadStats    stats;

  if (homopolyCompress) {
    FILE *H = merylutil::openOutputFile(seqStoreName, '/', "homopolymerCompression");
    merylutil::closeFile(H);
  }

  for (uint32 ll=0; ll<libraries.size(); ll++) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Creating library '%s' for %s %s reads.\n",
            libraries[ll]._name,
            toString(libraries[ll]._tech),
            toString(libraries[ll]._stat));

    if ((libraries[ll]._tech == sqTechType_pacbio_hifi) &&
        (libraries[ll]._stat  & sqRead_raw)) {
      fprintf(stderr, "ERROR: HiFi reads must be loaded as 'corrected'.\n");
      exit(1);
    }

    stats.displayTableHeader(stderr);

    seqLibrary = seqStore->sqStore_addEmptyLibrary(libraries[ll]._name, libraries[ll]._tech);


    for (uint32 ff=0; ff<libraries[ll]._files.size(); ff++) {
      char *file = libraries[ll]._files[ff];

      if (fileExists(file) == false) {
        fprintf(stderr, "ERROR:  sequence file '%s' not found.\n", file);
      }

      else {
        loadReads(seqStore,
                  seqLibrary,
                  libraries[ll]._stat,
                  minReadLength,
                  homopolyCompress,
                  nameMap,
                  errorLog,
                  file,
                  stats);
      }
    }
  }


  delete seqStore;

  merylutil::closeFile(nameMap,  seqStoreName, '/', "readNames.txt");
  merylutil::closeFile(errorLog, seqStoreName, '/', "errorLog");

  fprintf(stderr, "\n");
  fprintf(stderr, "All reads processed.\n");

  stats.displayTableHeader(stderr);
  stats.displayTable(stderr);

  return(true);
}



struct rl_t {
  uint32   readID;
  uint32   length;
  double   score;
};



bool
deleteShortReads(const char *seqStoreName,
                 uint64      genomeSize,
                 double      desiredCoverage,
                 double      lengthBias,
                 uint32      randomSeed) {
  mtRandom     mtctx;
  uint64       desiredBases = (uint64)floor(genomeSize * desiredCoverage);

  if (desiredBases == 0)
    return(true);

  if (randomSeed != 0)
    mtctx.mtSetSeed(randomSeed);

  //
  //  Open the store for modification.  Request whatever read is the default,
  //  in the uncompressed format.
  //

  sqStore     *seqStore     = new sqStore(seqStoreName, sqStore_extend);

  sqRead_setDefaultVersionExclude(sqRead_compressed);

  uint32       nReads       = seqStore->sqStore_lastReadID();
  rl_t        *readLen      = new rl_t [nReads + 1];

  uint32       readsInpt = 0, readsKept = 0,  readsRmvd = 0;
  uint64       basesInpt = 0, basesKept = 0,  basesRmvd = 0;

  //
  //  Initialize our list of read scores, and report a summary.
  //

  readLen[0].readID = UINT32_MAX;    //  Initialize the non-existent
  readLen[0].length = UINT32_MAX;    //  zeroth read.
  readLen[0].score  = 0;

  for (uint32 ii=1; ii<nReads+1; ii++) {
    uint32  len = seqStore->sqStore_getReadLength(ii);

    readLen[ii].readID = ii;
    readLen[ii].length = len;
    readLen[ii].score  = mtctx.mtRandomRealOpen() * pow(len, lengthBias);

    basesInpt += len;
    readsInpt += 1;
  }

  //
  //  Sample reads if we have more than the desired coverage.
  //

  if (desiredBases < basesInpt) {
    fprintf(stderr, "\n");
    fprintf(stderr, "EXCESSIVE COVERAGE DETECTED.  Sampling reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "For genome size of %12" F_U64P " bases,\n", genomeSize);
    fprintf(stderr, "            retain %12" F_U64P " bases (%.2fX coverage).\n", desiredBases, desiredCoverage);
    fprintf(stderr, "\n");
    fprintf(stderr, "Found    %9" F_U32P " reads with %12" F_U64P " bases (%.2fX coverage).\n",  readsInpt, basesInpt, (double)basesInpt / genomeSize);

    //
    //  Sort the array by score, and keep only the longest reads
    //  up to coverage * genomeSize.
    //

    auto  byScore = [](const rl_t &a, const rl_t &b)   { return(a.score > b.score); };

    std::sort(readLen, readLen + nReads+1, byScore);

    FILE *LOG = merylutil::openOutputFile(seqStoreName, '/', "filteredReads");

    fprintf(LOG, "readID  ordinal   length        score\n");
    fprintf(LOG, "------- ------- -------- ------------\n");

    for (uint32 ii=0; ii<nReads; ii++) {        //  Start at ii=0.  The zeroth read has
      fprintf(LOG, "%-7u %7u %8u %12.4e%s\n",   //  min score, and is now at position
              readLen[ii].readID, ii,           //  [nReads], while the longest read is
              readLen[ii].length,               //  at position [0].
              readLen[ii].score,
              (basesKept < desiredBases) ? "" : " REMOVED");

      if (basesKept < desiredBases) {
        basesKept += readLen[ii].length;
        readsKept += 1;
      }

      else {
        seqStore->sqStore_setIgnored(readLen[ii].readID, true, true);

        basesRmvd += readLen[ii].length;
        readsRmvd += 1;
      }
    }

    fprintf(stderr, "Dropped  %9" F_U32P " reads with %12" F_U64P " bases (%.2fX coverage).\n",  readsRmvd, basesRmvd, (double)basesRmvd / genomeSize);
    fprintf(stderr, "Retained %9" F_U32P " reads with %12" F_U64P " bases (%.2fX coverage).\n",  readsKept, basesKept, (double)basesKept / genomeSize);

    merylutil::closeFile(LOG);
  }

  delete [] readLen;
  delete    seqStore;

  return(true);
}



int
addFiles(char **argv, int arg, int argc, std::vector<char const *> &err, seqLib &lib) {
  char  *opt = argv[arg++];   //  e.g., -pacbio-raw
  char  *nam = argv[arg++];   //  e.g., LIB_NAME

  //  First word must be a library name and not a file.

  if (fileExists(nam) == true) {
    char *s = new char [1024];
    snprintf(s, 1024, "Option %s expects library name; '%s' is a file.\n", opt, nam);
    err.push_back(s);
  }

  //  All later words, until the end of the options or the next option, must be files.

  while ((arg < argc) &&
         (argv[arg][0] != '-')) {

    if (fileExists(argv[arg]) == false) {
      char *s = new char [1024];
      snprintf(s, 1024, "Option %s: file '%s' not found.\n", opt, argv[arg]);
      err.push_back(s);
    }

    lib.addFile(argv[arg++]);
  }

  //  Return the arg that leaves us right on the last file.

  return(arg - 1);
}



int
main(int argc, char **argv) {
  char const          *seqStoreName      = NULL;

  uint32               minReadLength     = 0;
  uint64               genomeSize        = 0;
  double               desiredCoverage   = 0;
  double               lengthBias        = 1.0;
  uint32               randomSeed        = 0;

  bool                 homopolyCompress  = false;

  std::vector<seqLib>  libraries;

  sqRead_which         readStatus        = sqRead_raw;

  //  Initialize the global.

  validSeq['a'] = validSeq['c'] = validSeq['g'] = validSeq['t'] = validSeq['u'] = validSeq['n'] = 1;
  validSeq['A'] = validSeq['C'] = validSeq['G'] = validSeq['T'] = validSeq['U'] = validSeq['N'] = 1;

  //  Parse options.

  argc = AS_configure(argc, argv, 1);

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if        (strcmp(argv[arg], "-o") == 0) {
      seqStoreName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-minlength") == 0) {
      minReadLength = atoi(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-genomesize") == 0) {
      genomeSize = atof(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-coverage") == 0) {
      desiredCoverage = atof(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-bias") == 0) {
      lengthBias = atof(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-seed") == 0) {
      randomSeed = atof(argv[++arg]);
    }


    else if (strcmp(argv[arg], "-homopolycompress") == 0) {
      homopolyCompress = true;
    }

    else if (strcmp(argv[arg], "-raw") == 0) {
      readStatus &= ~sqRead_corrected;
      readStatus |=  sqRead_raw;
    }

    else if (strcmp(argv[arg], "-corrected") == 0) {
      readStatus &= ~sqRead_raw;
      readStatus |=  sqRead_corrected;
    }

    else if (strcmp(argv[arg], "-untrimmed") == 0) {
      readStatus &= ~sqRead_trimmed;
    }

    else if (strcmp(argv[arg], "-trimmed") == 0) {
      readStatus |=  sqRead_trimmed;
    }


    else if (strcmp(argv[arg], "-pacbio") == 0) {
      seqLib  lib(argv[arg+1], sqTechType_pacbio, readStatus);
      arg = addFiles(argv, arg, argc, err, lib);
      libraries.push_back(lib);
    }

    else if (strcmp(argv[arg], "-nanopore") == 0) {
      seqLib  lib(argv[arg+1], sqTechType_nanopore, readStatus);
      arg = addFiles(argv, arg, argc, err, lib);
      libraries.push_back(lib);
    }

    else if (strcmp(argv[arg], "-pacbio-hifi") == 0) {
      sqRead_which   rs = readStatus;

      readStatus &= ~sqRead_raw;         //  HiFi MUST be loaded as corrected.
      readStatus |=  sqRead_corrected;   //

      seqLib  lib(argv[arg+1], sqTechType_pacbio_hifi, readStatus);
      arg = addFiles(argv, arg, argc, err, lib);
      libraries.push_back(lib);

      readStatus = rs;
    }


    else if (fileExists(argv[arg]) == true) {
      char *s = new char [1024];
      snprintf(s, 1024, "File with no library supplied: '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  if (seqStoreName == NULL)
    err.push_back("ERROR: no seqStore (-o) supplied.\n");

  if (libraries.size() == 0)
    err.push_back("ERROR: no input libraries (-pacbio-raw, etc) supplied.\n");

  if ((desiredCoverage > 0) && (genomeSize == 0))
    err.push_back("ERROR: no genome size (-genomesize) set, needed for coverage filtering (-coverage) to work.\n");

  if (err.size() > 0) {
    int32 al = strlen(argv[0]);

    fprintf(stderr, "usage: %s -o S.seqStore\n", argv[0]);
    fprintf(stderr, "       %*s [options] \\\n", al, "");
    fprintf(stderr, "       %*s [[processing-options] technology-option libName reads ...] ...\n", al, "");
    fprintf(stderr, "  -o S.seqStore          create output S.seqStore and load reads into it\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -minlength L           discard reads shorter than L (regardless of coverage)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -homopolycompress      set up for accessing homopolymer compressed reads\n");
    fprintf(stderr, "                         by default; also compute coverage and filter lengths\n");
    fprintf(stderr, "                         using the compressed read sequence.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "COVERAGE FILTERING\n");
    fprintf(stderr, "  When more than C coverage in reads is supplied, random reads are removed\n");
    fprintf(stderr, "  until coverage is C.  Bias B will remove shorter (B > 0) or longer (B < 0)\n");
    fprintf(stderr, "  reads preferentially.  B=0 will remove random reads.  Default is B=1.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -genomesize G          expected genome size (needed to compute coverage)\n");
    fprintf(stderr, "  -coverage C            desired coverage in long reads\n");
    fprintf(stderr, "  -bias B                remove shorter (B > 0) or longer (B < 0) reads.\n");
    fprintf(stderr, "  -seed S                seed the pseudo random number generator with S\n");
    fprintf(stderr, "                           1 <= S <= 4294967295\n");
    fprintf(stderr, "                           S = 0 will use a seed derived from the time and process id\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "READ SPECIFICATION\n");
    fprintf(stderr, "  Reads are supplied as a collection of libraries.  Each library should contain\n");
    fprintf(stderr, "  all the reads from one sequencing experiment (e.g., sample collection, sample\n");
    fprintf(stderr, "  preperation, sequencing run).  A library is created when any of the 'read\n");
    fprintf(stderr, "  technology' options is encountered, and will use whatever 'processing state'\n");
    fprintf(stderr, "  have been already supplied.  The first word after a 'read technology' option\n");
    fprintf(stderr, "  must be the name of the library.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Note that -pacbio-hifi will force -corrected status.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Example:  '-raw -pacbio LIBRARY_1 file.fasta'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -raw                   set the 'processing state' of the reads\n");
    fprintf(stderr, "  -corrected               next on the command line.\n");
    fprintf(stderr, "  -untrimmed\n");
    fprintf(stderr, "  -trimmed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nanopore              set the 'read technology' of the reads\n");
    fprintf(stderr, "  -pacbio                  next on the command line.\n");
    fprintf(stderr, "  -pacbio-hifi\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  createStore(seqStoreName, libraries, minReadLength, homopolyCompress);

  deleteShortReads(seqStoreName, genomeSize, desiredCoverage, lengthBias, randomSeed);

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");
  exit(0);
}

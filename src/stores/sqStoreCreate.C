
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

#include "mt19937ar.H"

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
  char              *_name;
  sqLibrary_tech     _tech;
  sqRead_which       _stat;

  vector<char *>     _files;
};



class loadStats {
public:
  loadStats() {
    nINVALID = nSHORT = nLONG = nLOADED=0;
    bINVALID = bSHORT = bLONG = bLOADED=0;
  };
  ~loadStats() {
  };

  void         import(loadStats &that) {
    nINVALID += that.nINVALID;   bINVALID += that.bINVALID;
    nSHORT   += that.nSHORT;     bSHORT   += that.bSHORT;
    nLONG    += that.nLONG;      bLONG    += that.bLONG;
    nLOADED  += that.nLOADED;    bLOADED  += that.bLOADED;
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

    fprintf(F, "\n");
  };

  uint32       nINVALID, nSHORT, nLONG, nLOADED;
  uint64       bINVALID, bSHORT, bLONG, bLOADED;
};





uint64
trimBgn(dnaSeq &sq, uint64 bgn, uint64 end) {
  char *bases = sq.bases();

  while ((bgn <= end) && ((bases[bgn] == 'N') ||
                          (bases[bgn] == 'n'))) {
    //bases[bgn] = 0;
    bgn++;
  }

  return(bgn);
}



uint64
trimEnd(dnaSeq &sq, uint64 bgn, uint64 end) {
  char *bases = sq.bases();

  end--;  //  So end is now the last base in the sequence.

  while ((bgn <= end) && ((bases[end] == 'N') ||
                          (bases[end] == 'n'))) {
    //bases[end] = 0;
    end--;
  }

  end++;  //  So now end is the C-style end.

  return(end);
}



uint32
checkInvalid(dnaSeq &sq, uint64 bgn, uint64 end) {
  uint32   invalid = 0;
  char    *bases   = sq.bases();

  for (uint32 ii=bgn; ii<end; ii++) {
    // special case Us
    if (bases[ii] == 'U' || bases[ii] == 'u')
       bases[ii] = 'T';
    if (validSeq[bases[ii]] == 0)
      invalid++;
  }

  return(invalid);
}



void
loadReads(sqStore          *seqStore,
          sqLibrary        *seqLibrary,
          sqRead_which      readStat,
          uint32            minReadLength,
          FILE             *nameMap,
          FILE             *errorLog,
          char             *fileName,
          loadStats        &stats) {

  //fprintf(stderr, "  %s:\n", fileName);

  loadStats    filestats;

  dnaSeqFile  *SF = new dnaSeqFile(fileName);
  dnaSeq       sq;

  while (SF->loadSequence(sq) == true) {

    //  Trim Ns from the ends of the sequence.
    uint64  bgn = trimBgn(sq, 0,   sq.length());
    uint64  end = trimEnd(sq, bgn, sq.length());

    if ((bgn > 0) && (end < sq.length()))
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - trimmed " F_U64 " non-ACGT bases from the 5' and " F_U64 " non-ACGT bases from the 3' end.\n",
              sq.name(), sq.length(), fileName, bgn, sq.length() - end);

    else if (bgn > 0)
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - trimmed " F_U64 " non-ACGT bases from the 5' end.\n",
              sq.name(), sq.length(), fileName, bgn);

    else if (end < sq.length())
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - trimmed " F_U64 " non-ACGT bases from the 3' end.\n",
              sq.name(), sq.length(), fileName, sq.length() - end);


    //  Check for invalid bases.
    uint32 invalid = checkInvalid(sq, bgn, end);

    if (invalid > 0) {
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - contains %u invalid letters, skipping.\n",
              sq.name(), sq.length(), fileName, invalid);

      filestats.nINVALID += 1;
      filestats.bINVALID += sq.length();

      continue;
    }


    //  Drop any sequences that are short.  This just sets length to zero, which we then skip.
    if (end - bgn < minReadLength) {
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - too short, skipping.\n",
              sq.name(), sq.length(), fileName);

      filestats.nSHORT += 1;
      filestats.bSHORT += sq.length();

      continue;
    }


    //  Warn if this sequence is too long.
    if (end - bgn > AS_MAX_READLEN - 2) {
      fprintf(errorLog, "read '%s' of length " F_U64 " in file '%s' - too long, skipping.\n",
              sq.name(), sq.length(), fileName);

      filestats.nLONG += 1;
      filestats.bLONG += sq.length();

      continue;
    }

    //  Create a writer for the read data and load bases.

    sqReadDataWriter *rdw = seqStore->sqStore_addEmptyRead(seqLibrary, sq.name());

    if (readStat & sqRead_raw) {
      rdw->sqReadDataWriter_setRawBases(sq.bases() + bgn, end - bgn);
    } else {
      rdw->sqReadDataWriter_setCorrectedBases(sq.bases() + bgn, end - bgn);
    }

    seqStore->sqStore_addRead(rdw);

    delete rdw;

    //  Now that the read is added to the store, we can set trim points.
    //  Presently, trimming only occurs on corrected reads, but later we
    //  need to allow trimmed raw reads.

    if (readStat & sqRead_trimmed) {
      uint32      rid  = seqStore->sqStore_lastReadID();
      sqReadSeq  *nseq = seqStore->sqStore_getReadSeq(rid, sqRead_corrected);
      sqReadSeq  *cseq = seqStore->sqStore_getReadSeq(rid, sqRead_corrected | sqRead_compressed);

      nseq->sqReadSeq_setAllClear();
      cseq->sqReadSeq_setAllClear();
    }

    //  And also update our nameMap.

    fprintf(nameMap, F_U32"\t%s\n", seqStore->sqStore_lastReadID(), sq.name());

    //  Save some silly statistics.

    filestats.nLOADED += 1;
    filestats.bLOADED += end - bgn;
  }

  delete SF;

  //  Write status to the screen
  filestats.displayTable(stderr, fileName);

  //  Add the just loaded numbers to the global numbers
  stats.import(filestats);
};





bool
createStore(const char       *seqStoreName,
            vector<seqLib>   &libraries,
            uint32            minReadLength) {

  sqStore     *seqStore     = new sqStore(seqStoreName, sqStore_create);   //  sqStore_extend MIGHT work
  sqRead      *seqRead      = NULL;
  sqLibrary   *seqLibrary   = NULL;

  uint32       inLineLen    = 1024;
  char         inLine[1024] = { 0 };

  FILE        *errorLog = AS_UTL_openOutputFile(seqStoreName, '/', "errorLog");
  FILE        *nameMap  = AS_UTL_openOutputFile(seqStoreName, '/', "readNames.txt");

  loadStats    stats;

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
                  nameMap,
                  errorLog,
                  file,
                  stats);
      }
    }
  }


  delete seqStore;

  AS_UTL_closeFile(nameMap,  seqStoreName, '/', "readNames.txt");
  AS_UTL_closeFile(errorLog, seqStoreName, '/', "errorLog");

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

bool  byScore(const rl_t &a, const rl_t &b)   { return(a.score > b.score); }


bool
deleteShortReads(const char *seqStoreName,
                 uint64      genomeSize,
                 double      desiredCoverage,
                 double      lengthBias) {
  mtRandom     mtctx;
  uint64       desiredBases = (uint64)floor(genomeSize * desiredCoverage);

  if (desiredBases == 0)
    return(true);

  //
  //  Open the store for modification.
  //

  sqStore     *seqStore     = new sqStore(seqStoreName, sqStore_extend);

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

    sort(readLen, readLen + nReads+1, byScore);

    fprintf(stdout, "readID    length        score\n");
    fprintf(stdout, "------- -------- ------------\n");

    for (uint32 ii=0; ii<nReads; ii++) {         //  Start at ii=0.  The zeroth read has
      fprintf(stdout, "%-7u %8u %12.4f%s\n",     //  min score, and is now at position
              readLen[ii].readID,                //  [nReads], while the longest read is
              readLen[ii].length,                //  at position [0].
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
  }

  delete [] readLen;
  delete    seqStore;

  return(true);
}



int
addFiles(char **argv, int arg, int argc, vector<char const *> &err, seqLib &lib) {
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
  char const      *seqStoreName      = NULL;

  uint32           minReadLength     = 0;
  uint64           genomeSize        = 0;
  double           desiredCoverage   = 0;
  double           lengthBias        = 1.0;

  vector<seqLib>   libraries;

  sqRead_which     readStatus        = sqRead_raw;

  //  Initialize the global.

  validSeq['a'] = validSeq['c'] = validSeq['g'] = validSeq['t'] = validSeq['n'] = 1;
  validSeq['A'] = validSeq['C'] = validSeq['G'] = validSeq['T'] = validSeq['N'] = 1;

  //  Parse options.

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
  while (arg < argc) {
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

    arg++;
  }

  if (seqStoreName == NULL)
    err.push_back("ERROR: no seqStore (-o) supplied.\n");

  if (libraries.size() == 0)
    err.push_back("ERROR: no input libraries (-pacbio-raw, etc) supplied.\n");

  if ((desiredCoverage > 0) && (genomeSize == 0))
    err.push_back("ERROR: no genome size (-genomesize) set, needed for coverage filtering (-coverage) to work.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -o seqStore [-minlength L] [-genomesize G -coverage C] [-pacbio-raw NAME file1 ...]\n", argv[0]);
    fprintf(stderr, "  -o seqStore            load raw reads into new seqStore\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  -minlength L           discard reads shorter than L\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  -genomesize G          expected genome size, for keeping only the longest reads\n");
    fprintf(stderr, "  -coverage C            desired coverage in long reads\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  Reads are supplied as a collection of libraries.  Each library should\n");
    fprintf(stderr, "  contain all the reads from one sequencing experiment (e.g., sample collection,\n");
    fprintf(stderr, "  sample preperation, sequencing run).\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  The library is specified as a sequencing technology, a read processing status,\n");
    fprintf(stderr, "  and a unique library name.\n");
    fprintf(stderr, "    -technology-status LIBRARY_NAME seqFile1 [seqFile2] [...] \n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  Valid combinations of technology and status are:\n");
    fprintf(stderr, "    -pacbio-raw\n");
    fprintf(stderr, "    -pacbio-corrected\n");
    fprintf(stderr, "    -pacbio-trimmed\n");
    fprintf(stderr, "    -pacbio-raw\n");
    fprintf(stderr, "    -nanopore-raw\n");
    fprintf(stderr, "    -nanopore-corrected\n");
    fprintf(stderr, "    -nanopore-trimmed\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  createStore(seqStoreName, libraries, minReadLength);

  deleteShortReads(seqStoreName, genomeSize, desiredCoverage, lengthBias);

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");
  exit(0);
}

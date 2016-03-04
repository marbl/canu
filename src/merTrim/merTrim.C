
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
 *    src/AS_MER/merTrim.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-FEB-22 to 2014-APR-11
 *      are Copyright 2010-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2010-JUL-12
 *      are Copyright 2010 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
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
#include "AS_UTL_reverseComplement.H"
#include "AS_PER_gkpStore.H"
#include "AS_PER_encodeSequenceQuality.H"  //  QUALITY_MAX, QV conversion
#include "AS_OVS_overlapStore.H"

#include <algorithm>

#include "AS_MER_gkpStore_to_FastABase.H"

#include "bio++.H"
#include "sweatShop.H"
#include "existDB.H"
#include "positionDB.H"
#include "libmeryl.H"
#include "logMsg.H"

#include "merTrimResult.H"

//  There is a serious bug in storing the kmer count in the existDB that merTrim 1.40 exposes.  kmer
//  must be at least r1950 to fix the bug.
#if !defined(EXISTDB_H_VERSION) || (EXISTDB_H_VERSION < 1960)
#error kmer needs to be updated to at least r1960
#error note that the kmer svn url changed in mid December 2012, the old url does not have r1960
#endif

uint32  VERBOSE = 0;

#define ALLGOOD 1
#define ALLCRAP 2
#define ATTEMPTCORRECTION 3

//  Doesn't work, gets different results.
#undef USE_MERSTREAM_REBUILD

#undef TEST_TESTBASE

char *createAdapterString(bool adapIllumina, bool adap454);


class mertrimGlobalData {
public:
  mertrimGlobalData() {
    gkpPath                   = 0L;
    fqInputPath               = 0L;
    fqOutputPath              = 0L;

    fqVerifyPath              = 0L;

    merSize                   = 22;

    merCountsFile             = 0L;
    merCountsCache            = false;
    adapCountsFile            = 0L;
    adapIllumina              = false;
    adap454                   = false;

    compression               = 0;  //  Does not work!
    numThreads                = 4;

    beVerbose                 = false;
    forceCorrection           = false;
    correctMismatch           = true;
    correctIndel              = true;

    actualCoverage            = 0;          //  Estimate of coverage
    minCorrectFraction        = 1.0 / 3.0;  //  Base can be corrected if less than 1/3 coverage
    minCorrect                = 0;          //
    minVerifiedFraction       = 1.0 / 4.0;  //  Base can be corrected to something with only 1/4 coverage
    minVerified               = 0;          //

    endTrimDefault            = true;
    endTrimNum                = 0;

#if 0
    //  These never worked well.
    endTrimWinScale[0]        = 0.50;
    endTrimErrAllow[0]        = 2;

    endTrimWinScale[1]        = 0.25;
    endTrimErrAllow[1]        = 0;
#endif

    endTrimQV                 = '2';

    doTrimming                = true;

    discardZeroCoverage       = false;
    discardImperfectCoverage  = false;

    trimImperfectCoverage     = true;

    gkRead                    = NULL;

    fqInput                   = NULL;
    fqOutput                  = NULL;
    fqVerify                  = NULL;
    fqLog                     = NULL;

    genomicDB                 = NULL;
    adapterDB                 = NULL;

    resPath                   = NULL;
    resFile                   = NULL;

    gktBgn                    = 0;
    gktEnd                    = 0;
    gktCur                    = 0;
  };

  ~mertrimGlobalData() {
    delete gkRead;

    delete fqInput;
    delete fqOutput;
    delete fqVerify;
    delete fqLog;

    delete genomicDB;
    delete adapterDB;

    if (resFile != NULL)
      fclose(resFile);
  };

  void              initializeGatekeeper(void) {

    if (gkpPath == NULL)
      return;

    fprintf(stderr, "opening gkStore '%s'\n", gkpPath);
    gkRead  = gkStore::gkStore_open(gkpPath, FALSE, FALSE);

    if (gktBgn == 0) {
      gktBgn = 1;
      gktEnd = gkRead->gkStore_getNumFragments();
    }

    gktCur = gktBgn;

    if (gktBgn > gktEnd)
      fprintf(stderr, "ERROR: invalid range:  -b ("F_U32") >= -e ("F_U32").\n",
              gktBgn, gktEnd), exit(1);
    if (gktEnd > gkRead->gkStore_getNumFragments())
      fprintf(stderr, "ERROR: invalid range:  -e ("F_U32") > num frags ("F_U32").\n",
              gktEnd, gkRead->gkStore_getNumFragments()), exit(1);

    errno = 0;
    resFile = fopen(resPath, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", resPath, strerror(errno)), exit(1);
  };

  void              initializeFASTQ(void) {

    if (fqInputPath == NULL)
      return;

    if (fqOutputPath == NULL)
      return;

    fqInput  = new compressedFileReader(fqInputPath);
    fqOutput = new compressedFileWriter(fqOutputPath);

    if (fqVerifyPath)
      fqVerify = new compressedFileReader(fqVerifyPath);

    char fqName[FILENAME_MAX];

    sprintf(fqName, "%s.log", fqOutputPath);

    fqLog    = new compressedFileWriter(fqName);
  };

  void              initialize(void) {
    initializeGatekeeper();
    initializeFASTQ();

    if (actualCoverage == 0) {
      merylStreamReader  *MF = new merylStreamReader(merCountsFile);

      uint32  i  = 0;
      uint32  iX = 0;

      //fprintf(stderr, "distinct: "F_U64"\n", MF->numberOfDistinctMers());
      //fprintf(stderr, "unique:   "F_U64"\n", MF->numberOfUniqueMers());
      //fprintf(stderr, "total:    "F_U64"\n", MF->numberOfTotalMers());

      //fprintf(stderr, "Xcoverage zero 1 0 "F_U64"\n", MF->histogram(1));

      for (i=2; (i < MF->histogramLength()) && (MF->histogram(i-1) > MF->histogram(i)); i++)
        //fprintf(stderr, "Xcoverage drop "F_U32" "F_U64" "F_U64"\n", i, MF->histogram(i-1), MF->histogram(i));
        ;

      iX = i - 1;

      for (; i < MF->histogramLength(); i++) {
        if (MF->histogram(iX) < MF->histogram(i)) {
          //fprintf(stderr, "Xcoverage incr "F_U32" "F_U64" "F_U64"\n", i, MF->histogram(iX), MF->histogram(i));
          iX = i;
        } else {
          //fprintf(stderr, "Xcoverage drop "F_U32" "F_U64" "F_U64"\n", i, MF->histogram(iX), MF->histogram(i));
        }
      }

      fprintf(stderr, "Guessed X coverage is "F_U32"\n", iX);

      delete MF;

      actualCoverage = iX;
    }

    if (minCorrectFraction > 0)
      minCorrect = (uint32)floor(minCorrectFraction * actualCoverage);

    if (minVerifiedFraction > 0)
      minVerified = (uint32)floor(minVerifiedFraction * actualCoverage);

    fprintf(stderr, "Use minCorrect="F_U32" minVerified="F_U32"\n", minCorrect, minVerified);

    if (minCorrect < minVerified) {
      fprintf(stderr, "WARNING!\n");
      fprintf(stderr, "WARNING!  minVerified (-verified) should be less than minCorrect (-correct).\n");
      fprintf(stderr, "WARNING!\n");
    }

    if (adapCountsFile) {
      fprintf(stderr, "loading adapter mer database.\n");
      adapterDB = new existDB(adapCountsFile, merSize, existDBcounts, 0, UINT32_MAX);
      adapterDB->printState(stderr);

    } else if (adapIllumina || adap454) {
      fprintf(stderr, "creating adapter mer database.\n");

      char *adapter = createAdapterString(adapIllumina, adap454);

      adapterDB = new existDB(adapter, merSize, existDBcanonical | existDBcounts);
      //adapterDB->printState(stderr);

      delete [] adapter;

    } else {
      fprintf(stderr, "not searching for adapter.\n");
    }

    char  cacheName[FILENAME_MAX];
    sprintf(cacheName, "%s.merTrimDB", merCountsFile);

    if (AS_UTL_fileExists(cacheName, FALSE, FALSE)) {
      fprintf(stderr, "loading genome mer database from cache '%s'.\n", cacheName);
      genomicDB = new existDB(cacheName);

    } else if (merCountsFile) {
      fprintf(stderr, "loading genome mer database from meryl '%s'.\n", merCountsFile);
      genomicDB = new existDB(merCountsFile, merSize, existDBcounts, MIN(minCorrect, minVerified), UINT32_MAX);

      if (merCountsCache) {
        fprintf(stderr, "saving genome mer database to cache '%s'.\n", cacheName);
        genomicDB->saveState(cacheName);
      }
    }
  };

public:

  //  Command line parameters
  //
  char         *gkpPath;
  char         *fqInputPath;
  char         *fqOutputPath;

  char         *fqVerifyPath;

  uint32        merSize;

  char         *merCountsFile;
  bool          merCountsCache;
  char         *adapCountsFile;

  bool          adapIllumina;
  bool          adap454;

  uint32        compression;
  uint32        numThreads;

  bool          beVerbose;
  bool          forceCorrection;
  bool          correctMismatch;
  bool          correctIndel;

  char         *resPath;
  FILE         *resFile;

  bool          endTrimDefault;
  uint32        endTrimNum;
  double        endTrimWinScale[16];
  uint32        endTrimErrAllow[16];

  char          endTrimQV;

  bool          doTrimming;

  bool          discardZeroCoverage;
  bool          discardImperfectCoverage;

  bool          trimImperfectCoverage;

  //  Global data
  //
  gkStore      *gkRead;

  uint32        actualCoverage;
  double        minCorrectFraction;
  uint32        minCorrect;
  double        minVerifiedFraction;
  uint32        minVerified;

  compressedFileReader  *fqInput;
  compressedFileWriter  *fqOutput;
  compressedFileReader  *fqVerify;
  compressedFileWriter  *fqLog;

  existDB      *genomicDB;
  existDB      *adapterDB;

  //  Input State
  //
  uint32        gktBgn;
  uint32        gktCur;
  uint32        gktEnd;
};




class mertrimThreadData {
public:
  mertrimThreadData(mertrimGlobalData *g) {
    kb     = new kMerBuilder(g->merSize, g->compression, 0L);
  };
  ~mertrimThreadData() {
    delete kb;
  };

public:
  kMerBuilder  *kb;
};




class mertrimComputation {
public:
  mertrimComputation() : log(false, 131072) {
    readName   = NULL;

    origSeq    = NULL;
    origQlt    = NULL;
    corrSeq    = NULL;
    corrQlt    = NULL;
    seqMap     = NULL;

    verifyName = NULL;
    verifySeq  = NULL;
    verifyErr   = NULL;

    rMS        = NULL;

    disconnect = NULL;
    coverage   = NULL;
    adapter    = NULL;
    corrected  = NULL;

    eDB        = NULL;
  }
  ~mertrimComputation() {
    delete [] readName;
    delete [] origSeq;
    delete [] origQlt;
    delete [] corrSeq;
    delete [] corrQlt;

    delete [] verifyName;
    delete [] verifySeq;
    delete [] verifyErr;

    delete [] seqMap;

    delete    rMS;

    delete [] disconnect;
    delete [] coverage;
    delete [] adapter;
    delete [] corrected;
  }


  void   initializeGatekeeper(mertrimGlobalData *g_) {
    g  = g_;

    readIID  = fr.gkFragment_getReadIID();
    seqLen   = fr.gkFragment_getSequenceLength();
    allocLen = seqLen + seqLen;

    readName = NULL;

    origSeq  = new char   [allocLen];
    origQlt  = new char   [allocLen];
    corrSeq  = new char   [allocLen];
    corrQlt  = new char   [allocLen];

    verifyName = NULL;
    verifySeq  = NULL;
    verifyErr  = NULL;

    seqMap   = new uint32 [allocLen];

    nMersExpected = 0;
    nMersTested   = 0;
    nMersFound    = 0;
    nMersCorrect  = 0;

    rMS = NULL;

    disconnect = NULL;
    coverage   = NULL;
    adapter    = NULL;
    corrected  = NULL;

    eDB        = NULL;

    strcpy(origSeq, fr.gkFragment_getSequence());
    strcpy(origQlt, fr.gkFragment_getQuality());
    strcpy(corrSeq, fr.gkFragment_getSequence());
    strcpy(corrQlt, fr.gkFragment_getQuality());

    //  Replace Ns with a random low-quality base.  This is necessary, since the mer routines
    //  will not make a mer for N, and we never see it to correct it.

    char  letters[4] = { 'A', 'C', 'G', 'T' };

    //  Not really replacing N with random ACGT, but good enough for us.
    for (uint32 i=0; i<seqLen; i++)
      if (corrSeq[i] == 'N') {
        corrSeq[i] = letters[i & 0x03];
        corrQlt[i] = 0;
      }

    clrBgn = 0;
    clrEnd = seqLen;

    garbageInInput        = false;
    gapInConfirmedKmers   = false;
    hasNoConfirmedKmers   = false;
    imperfectKmerCoverage = false;

    containsAdapter       = false;
    containsAdapterCount  = false;
    containsAdapterFixed  = false;
    containsAdapterBgn    = false;
    containsAdapterEnd    = false;

    suspectedChimer       = false;
    suspectedChimerBgn    = 0;
    suspectedChimerEnd    = 0;

    for (uint32 i=0; i<allocLen; i++)
      seqMap[i] = i;
  };


  bool    initializeFASTQ(mertrimGlobalData *g_) {
    g  = g_;

    readIID  = g->gktCur++;
    seqLen   = 0;
    allocLen = AS_READ_MAX_NORMAL_LEN + AS_READ_MAX_NORMAL_LEN + 1;  //  Used for seq/qlt storage only
    allocLen = 65536;

    readName   = new char   [1024];

    origSeq    = new char   [allocLen];
    origQlt    = new char   [allocLen];
    corrSeq    = new char   [allocLen];
    corrQlt    = new char   [allocLen];

    verifyName = NULL;  //  probably redundant
    verifySeq  = NULL;
    verifyErr  = NULL;

    seqMap     = new uint32 [allocLen];

    nMersExpected = 0;
    nMersTested   = 0;
    nMersFound    = 0;
    nMersCorrect  = 0;

    rMS = NULL;

    disconnect = NULL;
    coverage   = NULL;
    adapter    = NULL;
    corrected  = NULL;

    eDB        = NULL;

    //  Load the answer, if supplied (uses the real read storage space as temporary)

    if (g->fqVerify) {
      verifyName = new char   [1024];
      verifySeq  = new char   [allocLen];
      verifyErr  = new char   [allocLen];

      memset(verifyName, 1024,     0);
      memset(verifySeq,  allocLen, 0);  //  Needed, for verified reads shorter than real reads
      memset(verifyErr,  allocLen, '-');

      fgets(verifyName, 1024,     g->fqVerify->file());
      fgets(verifySeq,  allocLen, g->fqVerify->file());
      fgets(origQlt,    allocLen, g->fqVerify->file());  //  qv name line, ignored
      fgets(origQlt,    allocLen, g->fqVerify->file());

      chomp(verifyName);
      chomp(verifySeq);
    }

    //  Load a read to correct

    fgets(readName, 1024,     g->fqInput->file());
    fgets(origSeq,  allocLen, g->fqInput->file());
    fgets(origQlt,  allocLen, g->fqInput->file());  //  qv name line, ignored
    fgets(origQlt,  allocLen, g->fqInput->file());

    if (feof(g->fqInput->file()))
      return(false);

    chomp(readName);
    chomp(origSeq);
    chomp(origQlt);

    //  Adjust QVs to CA encoding, bases to upper case, and non-acgt to acgt.

#warning ASSUMING SANGER QV ENCODING

    for (uint32 i=0; origQlt[i]; i++) {
      if ('a' <= origSeq[i])
        origSeq[i] += 'A' - 'a';

      if (origQlt[i] < '!')
        fprintf(stderr, "ERROR: invalid QV '%c' (%d) in read '%s': '%s'\n",
                origQlt[i], origQlt[i], readName, origQlt);
      //  Our Sanger reads (dumped as fastq from gkpStore) have QV's higher than this.
      //if ('J' < origQlt[i])
      //  fprintf(stderr, "ERROR: invalid QV '%c' (%d) in read '%s': '%s'\n",
      //          origQlt[i], origQlt[i], readName, origQlt);

      origQlt[i] -= '!';
    }

    //  Copy to the corrected sequence

    strcpy(corrSeq, origSeq);
    strcpy(corrQlt, origQlt);

    //  Replace Ns with a random low-quality base.  This is necessary, since the mer routines
    //  will not make a mer for N, and we never see it to correct it.

    uint32  numReplace = 0;
    char    letters[4] = { 'A', 'C', 'G', 'T' };

    //  Not really replacing N with random ACGT, but good enough for us.
    for (uint32 i=0; corrSeq[i]; i++)
      if ((corrSeq[i] != 'A') &&
          (corrSeq[i] != 'C') &&
          (corrSeq[i] != 'G') &&
          (corrSeq[i] != 'T')) {
        numReplace++;
        corrSeq[i] = letters[i & 0x03];
        corrQlt[i] = 0;
      }

    seqLen   = strlen(origSeq);
    allocLen = seqLen + seqLen;

    clrBgn = 0;
    clrEnd = seqLen;

    garbageInInput        = false;
    gapInConfirmedKmers   = false;
    hasNoConfirmedKmers   = false;
    imperfectKmerCoverage = false;

    containsAdapter       = false;
    containsAdapterCount  = false;
    containsAdapterFixed  = false;
    containsAdapterBgn    = false;
    containsAdapterEnd    = false;

    suspectedChimer       = false;
    suspectedChimerBgn    = 0;
    suspectedChimerEnd    = 0;

    for (uint32 i=0; i<allocLen; i++)
      seqMap[i] = i;

    if (numReplace >= AS_OVL_ERROR_RATE * seqLen) {
      garbageInInput = true;
      clrBgn = 0;
      clrEnd = 0;
    }

    return(true);
  };


  uint32     evaluate(void);

  void       reverse(void);
  void       analyze(void);

  uint32     testBases(char *bases, uint32 basesLen);
  uint32     testBaseChange(uint32 pos, char replacement);
  uint32     testBaseIndel(uint32 pos, char replacement);

  bool       correctMismatch(uint32 pos, uint32 mNum, uint32 mExtra, bool isReversed);
  bool       correctIndel(uint32 pos, uint32 mNum, uint32 mExtra, bool isReversed);

  void       searchAdapter(bool isReversed);
  void       scoreAdapter(void);

  void       attemptCorrection(bool isReversed);

  void       attemptTrimming(bool doTrimming, char endTrimQV);
  void       attemptTrimming5End(uint32 *errorPos, uint32 endWindow, uint32 endAllowed);
  void       attemptTrimming3End(uint32 *errorPos, uint32 endWindow, uint32 endAllowed);

  uint32     getClrBgn(void) { return(seqMap[clrBgn]); };
  uint32     getClrEnd(void) { return(seqMap[clrEnd]); };
  uint32     getSeqLen(void) { return(seqMap[seqLen]); };

  void       dump(char *label);

  //  Public for the writer.
  gkFragment           fr;

  mertrimGlobalData   *g;
  mertrimThreadData   *t;

  AS_IID     readIID;
  uint32     seqLen;
  uint32     allocLen;

  char      *readName;

  char      *origSeq;
  char      *origQlt;
  char      *corrSeq;
  char      *corrQlt;

  char      *verifyName;
  char      *verifySeq;
  char      *verifyErr;

  uint32    *seqMap;

  uint32     nMersExpected;
  uint32     nMersTested;
  uint32     nMersFound;
  uint32     nMersCorrect;

  merStream *rMS;  //  kmers in the read, for searching against genomic kmers

  uint32     clrBgn;
  uint32     clrEnd;

  bool       garbageInInput;
  bool       gapInConfirmedKmers;
  bool       hasNoConfirmedKmers;
  bool       imperfectKmerCoverage;

  bool       containsAdapter;       //  Read hits adapter kmers
  uint32     containsAdapterCount;  //  Number of uncorrected adapter kmers hit
  uint32     containsAdapterFixed;  //  Number of corrected adapter kmers hit
  uint32     containsAdapterBgn;    //  Location of adapter
  uint32     containsAdapterEnd;    //  Location of adapter

  bool       suspectedChimer;
  uint32     suspectedChimerBgn;
  uint32     suspectedChimerEnd;

  uint32    *disconnect;  //  per base - a hole before this base
  uint32    *coverage;    //  per base - mer coverage
  uint32    *adapter;     //  per base - mer coverage in adapter kmers
  uint32    *corrected;   //  per base - type of correction here

  existDB   *eDB;

  uint32     nHole;  //  Number of spaces (between bases) with no mer coverage
  uint32     nCorr;  //  Number of bases corrected
  uint32     nFail;  //  Number of bases uncorrected because no answer found
  uint32     nConf;  //  Number of bases uncorrected because multiple answers found

  char       merstring[256];

  logMsg     log;
};




//  Scan the sequence, counting the number of kmers verified.  If we find all of them, we're done.
//
uint32
mertrimComputation::evaluate(void) {

  if (VERBOSE > 1)
    log.add("\nPROCESS read %d name %s\n", readIID, readName);

  if (garbageInInput == true)
    return(ALLCRAP);

  if (rMS == NULL)
    rMS = new merStream(t->kb, new seqStream(corrSeq, seqLen), false, true);

  rMS->rewind();

  nMersExpected = clrEnd - clrBgn - g->merSize + 1;
  nMersTested   = 0;
  nMersFound    = 0;
  nMersCorrect  = 0;

  while ((rMS->nextMer()) &&
         (rMS->thePositionInSequence() + g->merSize - 1 < clrEnd)) {
    if (rMS->thePositionInSequence() < clrBgn)
      //  Mer before the clear range begins.
      continue;

    if (clrEnd <= rMS->thePositionInSequence() + g->merSize - 1)
      //  Mer after the clear range ends
      continue;

    nMersTested++;

    //log.add("pos %d count %d\n",
    //        rMS->thePositionInSequence() + g->merSize - 1,
    //        eDB->count(rMS->theCMer()));

    if (eDB->count(rMS->theCMer()) >= g->minCorrect)
      //  We don't need to correct this kmer.
      nMersCorrect++;

    if (eDB->count(rMS->theCMer()) >= g->minVerified)
      //  We trust this mer.
      nMersFound++;
  }

  if (VERBOSE > 0)
    log.add("INITIAL read %u %s len %u has %u mers, %u correct and %u trusted.\n",
            readIID, readName, seqLen, nMersTested, nMersCorrect, nMersFound);

  if (nMersCorrect == nMersExpected)
    //  All mers correct, read is 100% verified!
    return(ALLGOOD);

  if (nMersFound == 0) {
    //  No trusted kMers found, read is 100% garbage (or 100% unique).
    hasNoConfirmedKmers = true;
    return(ALLCRAP);
  }

  //  Attempt correction.
  return(ATTEMPTCORRECTION);
}



void
mertrimComputation::reverse(void) {
  uint32  c = 0;
  uint32 *s = NULL;
  uint32 *S = NULL;

  reverseComplement(corrSeq, corrQlt, seqLen);

  uint32  cb = seqLen - clrBgn;
  uint32  ce = seqLen - clrEnd;

  clrBgn = ce;
  clrEnd = cb;

  delete rMS;
  rMS = new merStream(t->kb, new seqStream(corrSeq, seqLen), false, true);

  if (corrected) {
    s = corrected;
    S = corrected + seqLen - 1;

    while (s < S) {
      if (*s == 'X')  *s = 0;
      if (*S == 'X')  *S = 0;

      c    = *s;
      *s++ =  *S;
      *S-- =  c;
    }
  }

  if (adapter) {
    s = adapter;
    S = adapter + seqLen - 1;

    while (s < S) {
      c    = *s;
      *s++ =  *S;
      *S-- =  c;
    }
  }

  s = seqMap;
  S = seqMap + seqLen - 1;

  while (s < S) {
    c    = *s;
    *s++ =  *S;
    *S-- =  c;
  }
}


void
mertrimComputation::analyze(void) {

  if (rMS == NULL)
    return;

  rMS->rewind();

  if (coverage == NULL)
    coverage = new uint32 [allocLen];
  if (disconnect == NULL)
    disconnect = new uint32 [allocLen];

  memset(coverage,   0, sizeof(uint32) * (allocLen));
  memset(disconnect, 0, sizeof(uint32) * (allocLen));

  while (rMS->nextMer()) {
    uint32  posBgn = rMS->thePositionInSequence();
    uint32  posEnd = rMS->thePositionInSequence() + g->merSize;

    assert(posEnd <= seqLen);

    if (eDB->count(rMS->theCMer()) < g->minVerified)
      //  This mer is too weak for us.  SKip it.
      continue;

    //  If we aren't the first mer, then there should be coverage for our first base.  If not,
    //  we have found a correctable error, an uncorrectable error, or a chimeric read.
    if ((posBgn > 0) &&
        (coverage[posBgn-1] > 0) && (coverage[posBgn] == 0))
      disconnect[posBgn-1] = disconnect[posBgn] = 'D';

    //  Add coverage for the good mer.
    for (uint32 add=posBgn; add<posEnd; add++)
      coverage[add]++;

  }  //  Over all mers

  rMS->rewind();

  if (VERBOSE > 1)
    dump("ANALYZE");
}





bool
mertrimComputation::correctMismatch(uint32 pos, uint32 mNum, uint32 mExtra, bool isReversed) {
  uint32 nA = (corrSeq[pos] != 'A') ? testBaseChange(pos, 'A') : 0;
  uint32 nC = (corrSeq[pos] != 'C') ? testBaseChange(pos, 'C') : 0;
  uint32 nG = (corrSeq[pos] != 'G') ? testBaseChange(pos, 'G') : 0;
  uint32 nT = (corrSeq[pos] != 'T') ? testBaseChange(pos, 'T') : 0;
  uint32 rB = 0;  //  Base to change to
  uint32 rV = 0;  //  Count of that kmer evidence
  uint32 nR = 0;

  if (VERBOSE > 2) {
    if (nA > mNum + mExtra)  log.add("testA at %d -- %d req=%d\n", pos, nA, mNum + mExtra);
    if (nC > mNum + mExtra)  log.add("testC at %d -- %d req=%d\n", pos, nC, mNum + mExtra);
    if (nG > mNum + mExtra)  log.add("testG at %d -- %d req=%d\n", pos, nG, mNum + mExtra);
    if (nT > mNum + mExtra)  log.add("testT at %d -- %d req=%d\n", pos, nT, mNum + mExtra);
  }  //  VERBOSE

  //  If we found a single perfectly correct choice, ignore all the other solutions.

  if (nA == g->merSize)  nR++;
  if (nC == g->merSize)  nR++;
  if (nG == g->merSize)  nR++;
  if (nT == g->merSize)  nR++;

  if (nR == 1) {
    if (nA != g->merSize)  nA = 0;
    if (nC != g->merSize)  nC = 0;
    if (nG != g->merSize)  nG = 0;
    if (nT != g->merSize)  nT = 0;
  }

  //  Count the number of viable solutions.

  nR = 0;

  if (nA > mNum + mExtra)  { nR++;  rB = 'A';  rV = nA; }
  if (nC > mNum + mExtra)  { nR++;  rB = 'C';  rV = nC; }
  if (nG > mNum + mExtra)  { nR++;  rB = 'G';  rV = nG; }
  if (nT > mNum + mExtra)  { nR++;  rB = 'T';  rV = nT; }

  if (nR == 0)
    //  Nothing viable, keep the base as is.
    return(false);

  //  Something to change to.  By definition, this is a stronger mer than the original.

  if (nR > 1) {
    //  Multiple solutions.  Pick the most common.  If we don't do this, the confirmed kmer
    //  coverage drops at this location, and we trim the read.  This would result in
    //  zero coverage at every variation.

    uint32  mm = MAX(MAX(nA, nC), MAX(nG, nT));

    if (nA == mm)   { rB = 'A';  rV = nA; }
    if (nC == mm)   { rB = 'C';  rV = nC; }
    if (nG == mm)   { rB = 'G';  rV = nG; }
    if (nT == mm)   { rB = 'T';  rV = nT; }

    //corrected[pos] = 'X';
    //return(false);
  }

  //  One solution!  Correct it.

  if (VERBOSE > 0) {
    if (nR > 1)
      log.add("Correct read %d at position %d from %c (%u) to %c (%u) (QV %d) (%s) (multiple choices nA=%d nC=%d nG=%d nT=%d)\n",
              readIID,
              (isReversed == false) ? pos : seqLen - pos,
              corrSeq[pos], mNum,
              rB,           rV,
              corrQlt[pos],
              (isReversed == false) ? "fwd" : "rev",
              nA, nC, nG, nT);
    else
      log.add("Correct read %d at position %d from %c (%u) to %c (%u) (QV %d) (%s)\n",
              readIID,
              (isReversed == false) ? pos : seqLen - pos,
              corrSeq[pos], mNum,
              rB,           rV,
              corrQlt[pos],
              (isReversed == false) ? "fwd" : "rev");
  }  //  VERBOSE

  corrSeq[pos] = rB;

  corrected[pos] = 'C';

  //  Rebuild the merStream to use the corrected sequence, then move to the same position in
  //  the sequence.  This is done because we cannot simply change the string -- we need to
  //  change the state of the kMerBuilder associated with the merStream, and we can't do that.

#ifdef USE_MERSTREAM_REBUILD
  rMS->rebuild();
#else
  pos = rMS->thePositionInSequence();

  delete rMS;
  rMS = new merStream(t->kb, new seqStream(corrSeq, seqLen), false, true);
  rMS->nextMer();

  while (pos != rMS->thePositionInSequence())
    rMS->nextMer();
#endif

  return(true);
}



bool
mertrimComputation::correctIndel(uint32 pos, uint32 mNum, uint32 mExtra, bool isReversed) {
  uint32 nD = testBaseIndel(pos, '-');
  uint32 nA = testBaseIndel(pos, 'A');
  uint32 nC = testBaseIndel(pos, 'C');
  uint32 nG = testBaseIndel(pos, 'G');
  uint32 nT = testBaseIndel(pos, 'T');
  char   rB = 0;
  uint32 rV = 0;
  uint32 nR = 0;

  if (nD > mNum + mExtra)  { nR++;  rB = '-';  rV = nD;  }
  if (nA > mNum + mExtra)  { nR++;  rB = 'A';  rV = nA;  }
  if (nC > mNum + mExtra)  { nR++;  rB = 'C';  rV = nC;  }
  if (nG > mNum + mExtra)  { nR++;  rB = 'G';  rV = nG;  }
  if (nT > mNum + mExtra)  { nR++;  rB = 'T';  rV = nT;  }

  if (VERBOSE > 2) {
    if (nD > mNum + mExtra)  log.add("test-- %d -- %d req=%d\n", pos, nD, mNum + mExtra);
    if (nA > mNum + mExtra)  log.add("test+A %d -- %d req=%d\n", pos, nA, mNum + mExtra);
    if (nC > mNum + mExtra)  log.add("test+C %d -- %d req=%d\n", pos, nC, mNum + mExtra);
    if (nG > mNum + mExtra)  log.add("test+G %d -- %d req=%d\n", pos, nG, mNum + mExtra);
    if (nT > mNum + mExtra)  log.add("test+T %d -- %d req=%d\n", pos, nT, mNum + mExtra);
  }  //  VERBOSE

  if (nR == 0)
    //  No solutions.
    return(false);

  if (nR > 1) {
    //  Multiple solutions.
    corrected[pos] = 'X';
    return(false);
  }

  //  One solution.  Make a correction.  Either a deletion or an insert.

  if (nD > mNum + mExtra) {
    if (VERBOSE > 0) {
      log.add("Correct read %d at position %d from %c (%u) to DELETE (%u) (QV %d) (%s)\n",
              readIID,
              (isReversed == false) ? pos : seqLen - pos,
              corrSeq[pos], mNum,
              rV,
              corrQlt[pos],
              (isReversed == false) ? "fwd" : "rev");

    }  //  VERBOSE
    for (uint32 i=pos; i<seqLen; i++) {
      corrSeq[i] = corrSeq[i+1];
      corrQlt[i] = corrQlt[i+1];
      if (adapter)
        adapter[i] = adapter[i+1];
      seqMap[i]  = seqMap[i+1];
    }

    seqLen--;
    clrEnd--;

    corrected[pos] = 'D';

  } else {
    if (VERBOSE > 0) {
      log.add("Correct read %d at position %d from . (%u) to INSERT %c (%u) (%s)\n",
              readIID,
              (isReversed == false) ? pos : seqLen - pos,
              mNum,
              rB, rV,
              (isReversed == false) ? "fwd" : "rev");

    }  //  VERBOSE
    for (uint32 i=seqLen+1; i>pos; i--) {
      corrSeq[i] = corrSeq[i-1];
      corrQlt[i] = corrQlt[i-1];
      if (adapter)
        adapter[i] = adapter[i-1];
      seqMap[i]  = seqMap[i-1];
    }

    corrSeq[pos] = rB;
    corrQlt[pos] = '5';
    if (adapter)
      adapter[pos] = 0;
    seqMap[pos]  = seqMap[pos-1];

    seqLen++;
    clrEnd++;

    corrected[pos] = 'I';
  }

  //  Rebuild the merstream.  When we call analyze() the stream gets rewound.  If we do not
  //  restore our position, it is possible to get stuck in an infinite loop, inserting and
  //  deleting the same base over and over.  This happens because we don't explicitly require
  //  that the mer we are at be found, just that we find enough mers to make the change.
  //
  //  So, on the next pass through, we'd encounter the same mer we didn't find before, attempt
  //  to change it again, and possibly delete the base we inserted.
  //
  //  test+A 57 -- 6 req=2
  //  Correct read 6 at position 57 INSERT A
  //
  //  test-- 58 -- 3 req=2
  //  Correct read 6 at position 58 from A to DELETE (QV 5)
  //
  //  The first time through, we insert an A (with 6 mers agreeing).  The second time through,
  //  since we didn't fix the mer we were at, our choice is to delete the base we inserted (3
  //  mers tell us to do so).

#ifdef USE_MERSTREAM_REBUILD
  rMS->rebuild();
#else
  pos = rMS->thePositionInSequence();

  delete rMS;
  rMS = new merStream(t->kb, new seqStream(corrSeq, seqLen), false, true);

  analyze();

  rMS->nextMer();
  while (pos != rMS->thePositionInSequence())
    rMS->nextMer();
#endif

  return(true);
}




void
mertrimComputation::searchAdapter(bool isReversed) {

  if (corrected == NULL) {
    corrected = new uint32 [allocLen];
    memset(corrected, 0, sizeof(uint32) * (allocLen));
  }

  if (rMS == NULL)
    rMS = new merStream(t->kb, new seqStream(corrSeq, seqLen), false, true);

  rMS->rewind();

  //  A combination of evaluate() and attemptCorrection().  Find any adapter kmers in the
  //  read, do any corrections, then mark (in array 'adapter' the location of those adapter
  //  bases

  while (rMS->nextMer()) {
    uint32  pos   = rMS->thePositionInSequence() + g->merSize - 1;
    uint32  count = eDB->count(rMS->theCMer());

    if (count >= 1) {
      //  Mer exists, no need to correct.
      containsAdapter = true;
      continue;
    }

    uint32 mNum = testBaseChange(pos, corrSeq[pos]);

    //  Test if we can repair the sequence with a single base change.
    if (g->correctMismatch)
      if (correctMismatch(pos, mNum, 1, isReversed)) {
        containsAdapter = true;
        containsAdapterFixed++;
      }
  }
}



void
mertrimComputation::scoreAdapter(void) {

  if (containsAdapter == false)
    return;

  assert(adapter == NULL);

  adapter = new uint32 [allocLen];
  memset(adapter, 0, sizeof(uint32) * (allocLen));

  assert(clrBgn == 0);
  assert(clrEnd == seqLen);

  containsAdapterBgn = seqLen;
  containsAdapterEnd = 0;

  rMS->rewind();

  while (rMS->nextMer()) {
    uint32  bgn   = rMS->thePositionInSequence();
    uint32  end   = bgn + g->merSize - 1;
    uint32  count = eDB->count(rMS->theCMer());

    if (count == 0)
      continue;

    containsAdapterCount++;

    containsAdapterBgn = MIN(containsAdapterBgn, bgn);
    containsAdapterEnd = MAX(containsAdapterEnd, end + 1);

    if (VERBOSE > 1)
      log.add("ADAPTER at "F_U32","F_U32" ["F_U32","F_U32"]\n",
              bgn, end, containsAdapterBgn, containsAdapterEnd);

    for (uint32 a=bgn; a<=end; a++)
      adapter[a]++;
  }

  if (VERBOSE)
    dump("ADAPTERSEARCH");
}



void
mertrimComputation::attemptCorrection(bool isReversed) {

  if (corrected == NULL) {
    corrected = new uint32 [allocLen];
    memset(corrected, 0, sizeof(uint32) * (allocLen));
  }

  assert(coverage);
  assert(disconnect);
  assert(corrected);

  assert(rMS != NULL);

  rMS->rewind();

  while (rMS->nextMer()) {
    uint32  pos   = rMS->thePositionInSequence() + g->merSize - 1;
    uint32  count = eDB->count(rMS->theCMer());

    //log.add("MER at %d is %s has count %d %s\n",
    //        pos,
    //        rMS->theFMer().merToString(merstring),
    //        (count >= g->minCorrect) ? "CORRECT" : "ERROR",
    //        count);

    if (count >= g->minCorrect)
      //  Mer is correct, no need to correct it!
      continue;

    //  State the minimum number of mers we'd accept as evidence any change we make is correct.  The
    //  penalty for OVER correcting (too low a threshold) is potentially severe -- we could
    //  insert/delete over and over and over eventually blowing up.  The penalty for UNDER
    //  correcting, however, is that we trim a read too aggressively.
    //
    //  Being strictly greater than before works for mismatches and deletions.
    //
    //  For insertions, especially insertions in single nucleotide runs, this doesn't work so well.
    //  There are two cases, an insertion before a run, and an insertion after a run.  Below, X
    //  represents ACGT, and the A's are the run.
    //    XXXXXAAAAA: inserting an A after the X's, returns the same sequence, and the count of
    //                good mers doesn't change.
    //    AAAAAXXXXX: inserting an A before the X's changes the sequence, returning AAAAAAXXXX,
    //                and it is possible (likely) that this new sequence will have a higher count.
    //
    //  We therefore require that we find at least TWO more good mers before accepting a change.
    //
    //  The drawback of this is that we cannot correct two adjacent errors.  The first error
    //  (the one we're currently working on) is corrected and adds one to the coverage count,
    //  but then we hit that second error and do not find any more mers.
    //
    //  A solution would be to retry any base we cannot correct and allow a positive change of
    //  one mer to accept the change.  (in other words, change +1 below to +0).

    uint32 mNum = testBaseChange(pos, corrSeq[pos]);

    //  Test if we can repair the sequence with a single base change.
    if (g->correctMismatch)
      if (correctMismatch(pos, mNum, 1, isReversed))
        continue;

    if ((g->correctIndel) &&
        (g->merSize     < pos) &&
        (pos            < seqLen - g->merSize))
      if (correctIndel(pos, mNum, 3, isReversed))
        continue;
  }

  if (VERBOSE > 1) {
    dump("POSTCORRECT");
  }  //  VERBOSE
}



uint32
mertrimComputation::testBases(char *bases, uint32 basesLen) {
  uint32  offset       = 0;
  uint32  numConfirmed = 0;

  //
  //  UNTESTED with KMER_WORDS != 1
  //

  kMer F(g->merSize);
  kMer R(g->merSize);

  for (uint32 i=1; i<g->merSize && offset<basesLen; i++, offset++) {
    F += letterToBits[bases[offset]];
    R -= letterToBits[complementSymbol[bases[offset]]];
  }

  for (uint32 i=0; i<g->merSize && offset<basesLen; i++, offset++) {
    F += letterToBits[bases[offset]];
    R -= letterToBits[complementSymbol[bases[offset]]];

    F.mask(true);
    R.mask(false);

    if (F < R) {
      if (eDB->count(F) >= g->minVerified)
        numConfirmed++;
    } else {
      if (eDB->count(R) >= g->minVerified)
        numConfirmed++;
    }
  }

  return(numConfirmed);
}



//  Attempt to change the base at pos to make the kmers spanning it agree.
//  Returns the number of kmers validated, and the letter to change to.
//
uint32
mertrimComputation::testBaseChange(uint32 pos, char replacement) {
  uint32   numConfirmed = 0;
  char     originalBase = corrSeq[pos];
  uint32   offset       = pos + 1 - g->merSize;

  corrSeq[pos] = replacement;

  numConfirmed = testBases(corrSeq + offset, MIN(seqLen - offset, 2 * g->merSize - 1));

#ifdef TEST_TESTBASE
  {
    uint32 oldConfirmed = 0;

    merStream  *localms = new merStream(new kMerBuilder(merSize, compression, 0L),
                                        new seqStream(corrSeq + offset, seqLen - offset),
                                        true,
                                        true);

    //  Test
    for (uint32 i=0; i<g->merSize && localms->nextMer(); i++)
      if (eDB->count(localms->theCMer()) >= g->minVerified)
        oldConfirmed++;

    delete localms;

    assert(oldConfirmed == numConfirmed);
  }
#endif

  corrSeq[pos] = originalBase;

  //if (numConfirmed > 0)
  //  log.add("testBaseChange() pos=%d replacement=%c confirmed=%d\n",
  //          pos, replacement, numConfirmed);

  return(numConfirmed);
}




uint32
mertrimComputation::testBaseIndel(uint32 pos, char replacement) {
  uint32   numConfirmed = 0;
  char     testStr[128] = {0};
  uint32   len          = 0;
  uint32   offset       = pos + 1 - g->merSize;
  uint32   limit        = g->merSize * 2 - 1;

  assert(2 * g->merSize < 120);  //  Overly pessimistic

  //  Copy the first merSize bases.

  while (len < g->merSize - 1)
    testStr[len++] = corrSeq[offset++];

  //  Copy the second merSize bases, but overwrite the last base in the first copy (if we're testing
  //  a de;etion) or insert a replacement base.

  if (replacement == '-') {
    offset++;
  } else {
    testStr[len++] = replacement;
  }

  //  Copy the rest of the bases.

  while ((len < limit) && (corrSeq[offset]))
    testStr[len++] = corrSeq[offset++];

  numConfirmed = testBases(testStr, len);

#ifdef TEST_TESTBASE
  {
    uint32 oldConfirmed = 0;

    merStream  *localms = new merStream(new kMerBuilder(merSize, compression, 0L),
                                        new seqStream(testStr, len),
                                        true,
                                        true);

    //  Test
    for (uint32 i=0; i<g->merSize && localms->nextMer(); i++)
      if (existDB->count(localms->theCMer()) >= g->minVerified)
        oldConfirmed++;

    delete localms;

    assert(oldConfirmed == numConfirmed);
  }
#endif

  //if (numConfirmed > 0)
  //  log.add("testBaseIndel() pos=%d replacement=%c confirmed=%d\n",
  //          pos, replacement, numConfirmed);

  return(numConfirmed);
}




void
mertrimComputation::attemptTrimming5End(uint32 *errorPos, uint32 endWindow, uint32 errAllow) {

  for (bool doTrim = true; doTrim; ) {
    uint32  endFound   = 0;
    uint32  endTrimPos = 0;

    uint32   bgn = clrBgn;
    uint32   end = clrBgn + endWindow;

    if (end > clrEnd)
      end = clrEnd;

    for (uint32 i=bgn; i<end; i++) {
      if (errorPos[i]) {
        endFound++;
        endTrimPos = i + 1;
      }
    }

    if (VERBOSE > 1)
      log.add("BGNTRIM found=%u pos=%u from %u to %u\n",
              endFound, endTrimPos,
              clrBgn, clrBgn + endWindow);

    if (endFound > errAllow)
      clrBgn = endTrimPos;

    else
      doTrim = false;
  }
}


void
mertrimComputation::attemptTrimming3End(uint32 *errorPos, uint32 endWindow, uint32 errAllow) {

  for (bool doTrim = true; doTrim; ) {
    uint32  endFound   = 0;
    uint32  endTrimPos = clrEnd;

    uint32   bgn = clrBgn;
    uint32   end = clrEnd;

    if (clrBgn + endWindow <= clrEnd)
      bgn = clrEnd - endWindow;

    assert(bgn >= clrBgn);
    assert(bgn <= clrEnd);

    for (int32 i=bgn; i<end; i++) {
      if (errorPos[i]) {
        endFound++;
        if (i < endTrimPos)
          endTrimPos = i;
      }
    }

    if (VERBOSE > 1)
      log.add("ENDTRIM found=%u pos=%u from %u to %u\n",
              endFound, endTrimPos,
              clrEnd - endWindow, clrEnd);

    if (endFound > errAllow)
      clrEnd = endTrimPos;

    else
      doTrim = false;
  }
}


void
mertrimComputation::attemptTrimming(bool doTrimming, char endTrimQV) {

  //  Just bail if the read is all junk.  Nothing to do here.
  //
  if ((garbageInInput) ||
      (hasNoConfirmedKmers)) {
    clrBgn = 0;
    clrEnd = 0;
    return;
  }

  if ((clrBgn == 0) &&
      (clrEnd == 0))
    return;

  assert(coverage != NULL);

  //  Lop off the ends with no confirmed kmers or a QV less than 3.  The Illumina '2' qv ('B' in
  //  Illumina 1.5+ encodings) is defined as "rest of read is < QV 15; do not use".
  //
  if (doTrimming) {
    while ((clrBgn < clrEnd) && ((coverage[clrBgn] == 0) || (corrQlt[clrBgn] <= endTrimQV)))
      clrBgn++;

    while ((clrEnd > clrBgn) && ((coverage[clrEnd-1] == 0) || (corrQlt[clrEnd-1] <= endTrimQV)))
      clrEnd--;
  }

  //log.add("TRIM: %d,%d (lop-off-ends)\n", clrBgn, clrEnd);
  //dump("TRIM");

  //  Deal with adapter.  We'll pick the biggest end that is adapter free for our sequence.
  //  If both ends have adapter, so be it; the read gets trashed.

  if (containsAdapter) {
    containsAdapterBgn = seqLen;
    containsAdapterEnd = 0;

    for (uint32 i=0; i<seqLen; i++) {
      if ((adapter[i] > 0) && (i < containsAdapterBgn))
        containsAdapterBgn = i;
      if ((adapter[i] > 0) && (containsAdapterEnd < i))
        containsAdapterEnd = i + 1;
    }

    if (containsAdapterBgn < clrBgn)
      containsAdapterBgn = clrBgn;
    if (clrEnd < containsAdapterEnd)
      containsAdapterEnd = clrEnd;

    uint32  bgnClear = containsAdapterBgn - clrBgn;
    uint32  endClear = clrEnd - containsAdapterEnd;

    if (bgnClear > endClear)
      //  Start is bigger
      clrEnd = containsAdapterBgn;
    else
      //  End is bigger
      clrBgn = containsAdapterEnd;

    if (clrBgn >= clrEnd) {
      clrBgn = 0;
      clrEnd = 0;
    }

    //log.add("TRIM: %d,%d (adapter)\n", clrBgn, clrEnd);
    //dump("TRIM");
  }

  if (doTrimming == false)
    return;

  //  Lop off the ends with no confirmed kmers (again)
  //
  while ((clrBgn < clrEnd) && (coverage[clrBgn] == 0))
    clrBgn++;

  while ((clrEnd > clrBgn) && (coverage[clrEnd-1] == 0))
    clrEnd--;

  //log.add("TRIM: %d,%d (lop-off-ends-2)\n", clrBgn, clrEnd);
  //dump("TRIM");

  //  True if there is an error at position i.
  //
  uint32  *errorPos = new uint32 [seqLen];

  for (uint32 i=0; i<seqLen; i++)
    errorPos[i] = ((corrected[i] == 'C') || (corrected[i] == 'I') || (corrected[i] == 'D'));


  //  If there are more than 'errAllow' corrections in the last 'winScale' (x kmerSize) bases of the
  //  read, trim all the errors off.

  for (uint32 i=0; i<g->endTrimNum; i++) {
    attemptTrimming5End(errorPos, g->merSize * g->endTrimWinScale[i], g->endTrimErrAllow[i]);
    attemptTrimming3End(errorPos, g->merSize * g->endTrimWinScale[i], g->endTrimErrAllow[i]);
  }

  delete [] errorPos;

  //  If there are zero coverage areas in the interior, trash the whole read.

  if (g->discardZeroCoverage == true) {
    for (uint32 i=clrBgn; i<clrEnd; i++)
      if (coverage[i] == 0)
        gapInConfirmedKmers = true;
  }

  //  If the coverage isn't perfect (ramp up, constant, ramp down), trash the whole read.

  if ((g->discardImperfectCoverage == true) &&
      (clrBgn < clrEnd) &&
      (clrEnd > 0)) {
    uint32  bgn = clrBgn;
    uint32  end = clrEnd - 1;

    while ((bgn + 1 < end) && (coverage[bgn] < coverage[bgn+1]))
      bgn++;

    while ((bgn + 1 < end) && (coverage[end-1] > coverage[end]))
      end--;

    uint32  bgnc = coverage[bgn];
    uint32  endc = coverage[end];

    if (VERBOSE)
      log.add("IMPERFECT: bgn=%u %u  end=%u %u\n",
              bgn, bgnc,
              end, endc);

    if (bgnc != endc)
      imperfectKmerCoverage = true;

    else
      for (uint32 i=bgn; i<=end; i++)
        if (coverage[i] != bgnc)
          imperfectKmerCoverage = true;
  }

  //  If the coverage isn't perfect (see above), trim out the ends that make it imperfect.
  //  This leaves interior imperfect regions, but there should be enough to get
  //  an overlap on either end.

  if ((g->trimImperfectCoverage == true) &&
      (clrBgn < clrEnd) &&
      (clrEnd > 0)) {
    uint32  bgn = clrBgn;
    uint32  end = clrEnd - 1;

    //  Search from the ends in until we find the highest coverage.

    while ((bgn <= end) && (coverage[bgn] < g->merSize - 4))
      bgn++;

    while ((bgn <= end) && (coverage[end] < g->merSize - 4))
      end--;

    end++;

    //  Check that everything between is nice and high.  Why the -2?  This lets
    //  us tolerate one or two 'missing' kmers.

    uint32  mbgn = bgn;  //  Maximal bgn/end region
    uint32  mend = bgn;

    uint32  tbgn = bgn;  //  Currently testing region
    uint32  tend = bgn;

    while (tend < end) {
      if (coverage[tend] < g->merSize - 4) {
        if ((mend - mbgn) < (tend - tbgn)) {
          mbgn = tbgn;
          mend = tend;
        }

        tbgn = tend + 1;  //  Next region starts on the next position.
      }

      tend++;
    }

    //  If we found a sub region, reset to it (after remembering to check the last sub region).

    if (mbgn < mend) {
      if ((mend - mbgn) < (tend - tbgn)) {
        mbgn = tbgn;
        mend = tend;
      }

      if (VERBOSE)
        log.add("Reset clr from %d,%d to %d,%d\n", bgn, end, mbgn, mend);

      bgn = mbgn;
      end = mend;
    }

    //  Then search towards the ends, as long as the coverage is strictly decreasing.

    if (bgn < end) {
      while ((bgn > clrBgn) && (coverage[bgn-1] < coverage[bgn]) && (coverage[bgn-1] > 0))
        bgn--;

      while ((end < clrEnd) && (coverage[end] < coverage[end-1]) && (coverage[end] > 0))
        end++;

    } else {
      imperfectKmerCoverage = true;
    }

    clrBgn = bgn;
    clrEnd = end;
  }

  //  Make sense of the clear ranges.  If they're invalid, reset to 0,0.

  if (clrBgn >= clrEnd) {
    clrBgn = 0;
    clrEnd = 0;
  }

  if (VERBOSE > 1) {
    log.add("TRIM: %d,%d (post)\n", clrBgn, clrEnd);
    dump("TRIM");
  }  //  VERBOSE
}



void
mertrimComputation::dump(char *label) {
  char    *logLine = new char [4 * seqLen];
  uint32   logPos = 0;
  uint32   bogus  = (clrEnd == 0) ? 0 : UINT32_MAX;

  log.add("%s read %d %s len %d (trim %d-%d)\n", label, readIID, readName, seqLen, clrBgn, clrEnd);

  logPos = 0;
  for (uint32 i=0; origSeq[i]; i++) {
    if (i == clrBgn) { logLine[logPos++] = '-'; logLine[logPos++] = '['; }
    if (i == bogus)  { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }  //  read deleted
    logLine[logPos++] = origSeq[i];
    if (i+1 == clrEnd) { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
  }
  strcpy(logLine + logPos, " (ORI)\n");
  log.add(logLine);

  logPos = 0;
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn) { logLine[logPos++] = '-'; logLine[logPos++] = '['; }
    if (i == bogus)  { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
    logLine[logPos++] = corrSeq[i];
    if (i+1 == clrEnd) { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
  }
  strcpy(logLine + logPos, " (SEQ)\n");
  log.add(logLine);

  if (corrSeq && verifySeq) {
    uint32 i=0;

    for (; (i<seqLen) && (corrSeq[i] != 0) && (verifySeq[i] != 0); i++)
      verifyErr[i] = (corrSeq[i] == verifySeq[i]) ? '.' : '!';

    for (; i<seqLen; i++) {
      verifyErr[i] = '-';
      verifySeq[i] = '-';
    }

    logPos = 0;
    for (uint32 i=0; i<seqLen; i++) {
      if (i == clrBgn) { logLine[logPos++] = '-'; logLine[logPos++] = '['; }
      if (i == bogus)  { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
      logLine[logPos++] = (verifyErr) ? verifyErr[i] : ' ';
      if (i+1 == clrEnd) { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
    }
    strcpy(logLine + logPos, " (VAL)\n");
    log.add(logLine);

    logPos = 0;
    for (uint32 i=0; i<seqLen; i++) {
      if (i == clrBgn) { logLine[logPos++] = '-'; logLine[logPos++] = '['; }
      if (i == bogus)  { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
      logLine[logPos++] = (verifySeq && verifySeq[0]) ? verifySeq[i] : '-';
      if (i+1 == clrEnd) { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
    }
    strcpy(logLine + logPos, " (VAL)\n");
    log.add(logLine);
  }

  logPos = 0;
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn) { logLine[logPos++] = '-'; logLine[logPos++] = '['; }
    if (i == bogus)  { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
    logLine[logPos++] = corrQlt[i];
    if (i+1 == clrEnd) { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
  }
  strcpy(logLine + logPos, " (QLT)\n");
  log.add(logLine);

  logPos = 0;
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn) { logLine[logPos++] = '-'; logLine[logPos++] = '['; }
    if (i == bogus)  { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
    logLine[logPos++] = (coverage) ? coverage[i] + 'A' : 'A';
    if (i+1 == clrEnd) { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
  }
  strcpy(logLine + logPos, " (COVERAGE)\n");
  log.add(logLine);

  logPos = 0;
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn) { logLine[logPos++] = '-'; logLine[logPos++] = '['; }
    if (i == bogus)  { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
    logLine[logPos++] = (corrected && corrected[i]) ? corrected[i] : ' ';
    if (i+1 == clrEnd) { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
  }
  strcpy(logLine + logPos, " (CORRECTIONS)\n");
  log.add(logLine);

  logPos = 0;
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn) { logLine[logPos++] = '-'; logLine[logPos++] = '['; }
    if (i == bogus)  { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
    logLine[logPos++] = (disconnect && disconnect[i]) ? disconnect[i] : ' ';
    if (i+1 == clrEnd) { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
  }
  strcpy(logLine + logPos, " (DISCONNECTION)\n");
  log.add(logLine);

  logPos = 0;
  for (uint32 i=0; i<seqLen; i++) {
    if (i == clrBgn) { logLine[logPos++] = '-'; logLine[logPos++] = '['; }
    if (i == bogus)  { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
    logLine[logPos++] = (adapter && adapter[i]) ? adapter[i] + 'A' : ' ';
    if (i+1 == clrEnd) { logLine[logPos++] = ']'; logLine[logPos++] = '-'; }
  }
  strcpy(logLine + logPos, " (ADAPTER)\n");
  log.add(logLine);

  delete [] logLine;
}





void
mertrimWorker(void *G, void *T, void *S) {
  mertrimGlobalData    *g = (mertrimGlobalData  *)G;
  mertrimThreadData    *t = (mertrimThreadData  *)T;
  mertrimComputation   *s = (mertrimComputation *)S;

  s->t = t;

  s->eDB = g->genomicDB;

  uint32  eval = s->evaluate();

  //  Attempt correction if there are kmers to correct from.

  if (eval == ATTEMPTCORRECTION) {
    s->analyze();
    s->attemptCorrection(false);
    s->reverse();

    s->analyze();
    s->attemptCorrection(true);
    s->reverse();
  }

  //  Search for linker/adapter.  This needs to be after correction, since indel screws
  //  up the clear ranges we set in scoreAdapter().

  if ((eval != ALLCRAP) &&
      (g->adapterDB != NULL)) {
    s->eDB = g->adapterDB;

    s->searchAdapter(false);
    s->reverse();

    s->searchAdapter(true);
    s->reverse();

    s->scoreAdapter();

    s->eDB = g->genomicDB;
  }

  //  Attempt trimming if the read wasn't perfect

  if (eval != ALLGOOD) {
    s->analyze();
    s->attemptTrimming(g->doTrimming, g->endTrimQV);
  }

  if (VERBOSE)
    s->dump("FINAL");
}


mertrimComputation *
mertrimReaderGatekeeper(mertrimGlobalData *g) {
  mertrimComputation   *s = NULL;

  while ((g->gktCur <= g->gktEnd) &&
         (s == NULL)) {
    s = new mertrimComputation();

    g->gkRead->gkStore_getFragment(g->gktCur, &s->fr, GKFRAGMENT_QLT);
    g->gktCur++;

    if ((g->forceCorrection) ||
        (g->gkRead->gkStore_getLibrary(s->fr.gkFragment_getLibraryIID())->doTrim_initialMerBased)) {
      s->initializeGatekeeper(g);
    } else {
      delete s;
      s = NULL;
    }
  }

  return(s);
}


mertrimComputation *
mertrimReaderFASTQ(mertrimGlobalData *g) {
  mertrimComputation   *s = new mertrimComputation();

  if (s->initializeFASTQ(g) == false) {
    delete s;
    s = NULL;
  }

  return(s);
}


void *
mertrimReader(void *G) {
  mertrimGlobalData    *g = (mertrimGlobalData  *)G;
  mertrimComputation   *s = NULL;

  if (g->gkRead)
    s = mertrimReaderGatekeeper(g);

  if (g->fqInput)
    s = mertrimReaderFASTQ(g);

  return(s);
}



void
mertrimWriterGatekeeper(mertrimGlobalData *g, mertrimComputation *s) {
  mertrimResult         res;

  res.readIID = s->fr.gkFragment_getReadIID();

  if ((s->getClrEnd() <= s->getClrBgn()) ||
      (s->getClrEnd() - s->getClrBgn() < AS_READ_MIN_LEN))
    res.deleted = true;
  else
    res.deleted = false;

  res.clrBgn  = s->getClrBgn();
  res.clrEnd  = s->getClrEnd();

  res.chimer  = s->suspectedChimer;
  res.chmBgn  = s->suspectedChimerBgn;
  res.chmEnd  = s->suspectedChimerEnd;

  res.writeResult(g->resFile);
}


void
mertrimWriterFASTQ(mertrimGlobalData *g, mertrimComputation *s) {

  //  Note that getClrBgn/getClrEnd return positions in the ORIGINAL read, not the correct read.
  //  DO NOT USE HERE!

  uint32   seqOffset = 0;
  char     label[256];

  label[0] = 0;

  if (s->garbageInInput == true) {
    strcat(label, "DEL-GARBAGE");

    seqOffset = 0;

    s->corrSeq[0] = 0;
    s->corrQlt[0] = 0;

    s->clrBgn = 0;
    s->clrEnd = 0;

    goto outputFastq;
  }

  if (s->hasNoConfirmedKmers == true) {
    strcat(label, "DEL-NO-KMER");

    seqOffset = 0;

    s->corrSeq[0] = 0;
    s->corrQlt[0] = 0;

    s->clrBgn = 0;
    s->clrEnd = 0;

    goto outputFastq;
  }

  if ((s->containsAdapter == true) &&
      (s->suspectedChimer == false)) {
    strcat(label, "ADAPTERTRIM");

    seqOffset = s->clrBgn;

    s->corrSeq[s->clrEnd] = 0;
    s->corrQlt[s->clrEnd] = 0;

    goto outputFastq;
  }

  if ((s->containsAdapter == false) &&
      (s->suspectedChimer == true)) {
    if (s->suspectedChimerBgn - s->clrBgn >= AS_READ_MIN_LEN) {
      strcpy(label, "MP");  //  Junction read, longest portion would make an MP pair

      seqOffset = s->clrBgn;

      s->corrSeq[s->suspectedChimerBgn] = 0;
      s->corrQlt[s->suspectedChimerBgn] = 0;

    } else if (s->clrEnd - s->suspectedChimerEnd >= AS_READ_MIN_LEN) {
      strcpy(label, "PE");  //  Junction read, longest portion would make a PE pair

      seqOffset = s->suspectedChimerEnd;

      s->corrSeq[s->clrEnd] = 0;
      s->corrQlt[s->clrEnd] = 0;

    } else {
      strcpy(label, "MP-PE-SHORT");//  Junction read, neither side is long enough to save, so save nothing.

      seqOffset = 0;

      s->corrSeq[0] = 0;
      s->corrQlt[0] = 0;

      s->clrBgn = 0;
      s->clrEnd = 0;
    }

    goto outputFastq;
  }

  if (s->gapInConfirmedKmers == true) {
    strcat(label, "DEL-ZERO-COV");

    seqOffset = 0;

    s->corrSeq[0] = 0;
    s->corrQlt[0] = 0;

    s->clrBgn = 0;
    s->clrEnd = 0;

    goto outputFastq;
  }

  if (s->imperfectKmerCoverage == true) {
    strcat(label, "DEL-INPERFECT-COV");

    seqOffset = 0;

    s->corrSeq[0] = 0;
    s->corrQlt[0] = 0;

    s->clrBgn = 0;
    s->clrEnd = 0;

    goto outputFastq;
  }

  if ((s->suspectedChimer       == false) &&
      (s->containsAdapter       == false)) {
    strcat(label, "CLEAN");

    seqOffset = s->clrBgn;

    s->corrSeq[s->clrEnd] = 0;
    s->corrQlt[s->clrEnd] = 0;

    goto outputFastq;
  }

  assert(s->suspectedChimer);
  assert(s->containsAdapter);

  strcpy(label, "JUNCTION+ADAPTER");

  seqOffset = 0;

  s->corrSeq[0] = 0;
  s->corrQlt[0] = 0;



 outputFastq:

  //  If there is a verify sequence, verify.

  uint32  nVerifyErrors = 0;
#warning unfinished verify
  if (s->verifySeq) {
  }

  fprintf(g->fqLog->file(), F_U32"\t"F_U32"\tchimer\t%c\t"F_U32"\t"F_U32"\tadapter\t%c\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t%s\t%s\n",
          s->clrBgn,
          s->clrEnd,
          s->suspectedChimer ? 't' : 'f',
          s->suspectedChimerBgn,
          s->suspectedChimerEnd,
          s->containsAdapter ? 't' : 'f',
          s->containsAdapterCount,
          s->containsAdapterFixed,
          s->containsAdapterBgn,
          s->containsAdapterEnd,
          label,
          s->readName);

  //  Convert from CA QV to Sanger QV

  for (uint32 i=0; s->corrQlt[i]; i++)
    s->corrQlt[i] += '!';

  //  If no sequence, write a single N, with low QV (! == lowest).

  if ((s->corrSeq[seqOffset] == 0) ||
      (s->corrQlt[seqOffset] == 0))
    fprintf(g->fqOutput->file(), "%s type=%s\nN\n+\n!\n",
            s->readName,
            label);
  else
    fprintf(g->fqOutput->file(), "%s type=%s\n%s\n+\n%s\n",
            s->readName,
            label,
            s->corrSeq + seqOffset,
            s->corrQlt + seqOffset);

  if (VERBOSE)
    s->log.add("RESULT read %d len %d (trim %d-%d) %s\n", s->readIID, s->seqLen, s->clrBgn, s->clrEnd, label);
  s->log.fwrite(stdout);
}


void
mertrimWriter(void *G, void *S) {
  mertrimGlobalData    *g = (mertrimGlobalData  *)G;
  mertrimComputation   *s = (mertrimComputation *)S;

  //  The get*() functions return positions in the original uncorrected sequence.  They
  //  map from positions in the corrected sequence (which has inserts and deletes) back
  //  to the original sequence.
  //
  assert(s->getClrBgn() <= s->getClrEnd());
  assert(s->getClrEnd() <= s->getSeqLen());

  if (g->resFile)
    mertrimWriterGatekeeper(g, s);

  if (g->fqOutput)
    mertrimWriterFASTQ(g, s);

  delete s;
}










int
main(int argc, char **argv) {
  mertrimGlobalData  *g        = new mertrimGlobalData;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      g->gkpPath = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      g->fqInputPath = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      g->fqVerifyPath = argv[++arg];

    } else if (strcmp(argv[arg], "-m") == 0) {
      g->merSize = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-mc") == 0) {
      g->merCountsFile = argv[++arg];

    } else if (strcmp(argv[arg], "-enablecache") == 0) {
      g->merCountsCache = true;

    } else if (strcmp(argv[arg], "-mC") == 0) {
      g->adapCountsFile = argv[++arg];

    } else if (strcmp(argv[arg], "-mCillumina") == 0) {
      g->adapIllumina = 1;

    } else if (strcmp(argv[arg], "-mC454") == 0) {
      g->adap454 = 2;

    } else if (strcmp(argv[arg], "-t") == 0) {
      g->numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-b") == 0) {
      g->gktBgn = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      g->gktEnd = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-v") == 0) {
      g->beVerbose = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      VERBOSE++;

    } else if (strcmp(argv[arg], "-coverage") == 0) {
      g->actualCoverage = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-correct") == 0) {
      g->minCorrectFraction = atof(argv[++arg]);
      if (g->minCorrectFraction >= 1) {
        g->minCorrect         = (uint32)g->minCorrectFraction;
        g->minCorrectFraction = 0;
      }

    } else if (strcmp(argv[arg], "-evidence") == 0) {
      g->minVerifiedFraction = atof(argv[++arg]);
      if (g->minVerifiedFraction >= 1) {
        g->minVerified         = (uint32)g->minVerifiedFraction;
        g->minVerifiedFraction = 0;
      }

    } else if (strcmp(argv[arg], "-endtrim") == 0) {
      if (g->endTrimDefault == true)
        g->endTrimNum = 0;
      g->endTrimWinScale[g->endTrimNum] = atof(argv[++arg]);
      g->endTrimErrAllow[g->endTrimNum] = atoi(argv[++arg]);
      g->endTrimNum++;

    } else if (strcmp(argv[arg], "-notrimming") == 0) {
      g->doTrimming = false;

    } else if (strcmp(argv[arg], "-discardzero") == 0) {
      g->discardZeroCoverage = true;

    } else if (strcmp(argv[arg], "-discardimperfect") == 0) {
      g->discardImperfectCoverage = true;

    } else if (strcmp(argv[arg], "-notrimimperfect") == 0) {
      g->trimImperfectCoverage = false;

    } else if (strcmp(argv[arg], "-endtrimqv") == 0) {
      g->endTrimQV = argv[++arg][0];

    } else if (strcmp(argv[arg], "-f") == 0) {
      g->forceCorrection = true;

    } else if (strcmp(argv[arg], "-NM") == 0) {
      g->correctMismatch = false;

    } else if (strcmp(argv[arg], "-NI") == 0) {
      g->correctIndel = false;

    } else if (strcmp(argv[arg], "-o") == 0) {
      g->fqOutputPath = argv[++arg];
      g->resPath      = argv[arg];

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((g->gkpPath == 0L) && (g->fqInputPath == 0L))
    err++;

  if (err) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -F reads.fastq       input reads\n");
    fprintf(stderr, "  -o reads.fastq       output reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -T reads.fasta       truth reads for validation\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -m ms                mer size\n");
    fprintf(stderr, "  -mc counts           kmer database (in 'counts.mcdat' and 'counts.mcidx')\n");
    fprintf(stderr, "  -enablecache         dump the final kmer data to 'counts.merTrimDB'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -coverage C\n");
    fprintf(stderr, "  -correct n           mers with count below n can be changed\n");
    fprintf(stderr, "                         (that is, count >= n are correct mers)\n");
    fprintf(stderr, "  -evidence n          mers with count at least n will be used for changes\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -mC adapter.fasta    screen for these adapter sequences\n");
    fprintf(stderr, "  -mCillumina          screen for common Illumina adapter sequences\n");
    fprintf(stderr, "  -mC454               screen for common 454 adapter and linker sequences\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -endtrim             (undocumented)\n");
    fprintf(stderr, "  -notrimming          do only correction, no trimming\n");
    fprintf(stderr, "  -discardzero         trash the whole read if coverage drops to zero in the middle\n");
    fprintf(stderr, "  -discardimperfect    trash the whole read if coverage isn't perfect\n");
    fprintf(stderr, "  -notrimimperfect     do NOT trim off ends that make the coverage imperfect\n");
    fprintf(stderr, "  -endtrimqv Q         trim ends of reads if they are below qv Q (Sanger encoded; default '2')\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -NM                  do NOT correct mismatch errors\n");
    fprintf(stderr, "  -NI                  do NOT correct indel errors\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t T                 use T CPU cores\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -v                   report progress to stderr\n");
    fprintf(stderr, "  -V                   report trimming evidence to stdout (more -V -> more reports)\n");
    fprintf(stderr, "\n");

    exit(1);
  }

  gkpStoreFile::registerFile();

  g->initialize();

  gkFragment   fr;

#if 0
  //  DEBUG, non-threaded version.
  speedCounter SC(" Trimming: %11.0f reads -- %7.5f reads/second\r", 1.0, 0x1fff, true);

  mertrimThreadData *t = new mertrimThreadData(g);

  g->tBgn = 140222;
  g->tCur = 140222;
  g->tEnd = 140222;

  mertrimComputation *s = (mertrimComputation *)mertrimReader(g);
  while (s) {
    mertrimWorker(g, t, s);
    mertrimWriter(g, s);
    SC.tick();
    s = (mertrimComputation *)mertrimReader(g);
  }

  delete t;
#else
  //  PRODUCTION, threaded version
  sweatShop *ss = new sweatShop(mertrimReader, mertrimWorker, mertrimWriter);

  ss->setLoaderQueueSize(16384);
  ss->setWriterQueueSize(1024);

  ss->setNumberOfWorkers(g->numThreads);

  for (uint32 w=0; w<g->numThreads; w++)
    ss->setThreadData(w, new mertrimThreadData(g));  //  these leak

  ss->run(g, g->beVerbose);  //  true == verbose
#endif

  delete g;

  fprintf(stderr, "\nSuccess!  Bye.\n");

  return(0);
}

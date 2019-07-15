
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2019-APR-22
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */


class trReadData {
public:
  trReadData() {
    clrBgn         = 0;
    clrEnd         = 0;

    rawLength      = 0;
    trimmedLength  = 0;
  };
  ~trReadData() {
  };

  void      setClear(int32 bgn, int32 end) {
    clrBgn        = bgn;
    clrEnd        = end;
    trimmedLength = end - bgn;
  };

public:
  int32     clrBgn;
  int32     clrEnd;

  int32     rawLength;
  int32     trimmedLength;
};



class trGlobalData {
public:
  trGlobalData() {
    numThreads          = 1;

    maxErate            = 0.30;
    memLimit            = UINT64_MAX;

    bgnID               = 0;
    curID               = 0;
    endID               = UINT32_MAX;

    minReadLength       = 1000;

    minOverlapLength    = 1000;
    minOverlapCoverage  = 2;

    verboseTrim         = 0;
    verboseAlign        = 0;

    readData            = NULL;

    seqStoreName        = NULL;
    seqStore            = NULL;
    seqCache            = NULL;

    ovlStoreName        = NULL;
    ovlStore            = NULL;

    outFileName         = NULL;
    outFile             = NULL;

    clearRangesFileName = NULL;
  };

  void    initialize(sqStore_mode mode=sqStore_readOnly) {

    //  Open the store, assuming the reads are trimmed already.

    fprintf(stderr, "Opening seqStore '%s' in mode '%s'\n", seqStoreName, toString(mode));

    sqRead_setDefaultVersion(sqRead_raw);
    seqStore  = sqStore::sqStore_open(seqStoreName, mode);

    //  Load all the reads.  Regardless of trim status, we ALWAYS want
    //  to load raw reads, because we ALWAYS need to adjust overlaps
    //  from raw reads to trimmed reads.

    fprintf(stderr, "Loading all reads.\n");

    seqCache  = new sqCache(seqStore, sqRead_raw, memLimit);
    seqCache->sqCache_loadReads();

    //  Open overlaps.

    if (ovlStoreName) {
      fprintf(stderr, "Reading overlaps from store '%s'.\n", ovlStoreName);

      ovlStore = new ovStore(ovlStoreName, seqStore);
    }

    if (outFileName) {
      fprintf(stderr, "Writing overlaps to file '%s'.\n", outFileName);

      outFile = new ovFile(seqStore, outFileName, ovFileFullWrite);
    }

    //  Check parameters.

    if (bgnID == 0)                                bgnID = 1;
    if (endID  > seqStore->sqStore_getNumReads())  endID = seqStore->sqStore_getNumReads();

    curID = bgnID;

    //  Trimming data.  Initialize clear range and lengths to the full read.
    //  Data is loaded in main loop.

    readData = new trReadData [seqStore->sqStore_getNumReads() + 1];

    for (uint32 ii=0; ii<seqStore->sqStore_getNumReads() + 1; ii++) {
      sqRead  *read = seqStore->sqStore_getRead(ii);

      readData[ii].clrBgn        = read->sqRead_clearBgn();
      readData[ii].clrEnd        = read->sqRead_clearEnd();
      readData[ii].rawLength     = read->sqRead_sequenceLength(sqRead_raw);
      readData[ii].trimmedLength = read->sqRead_sequenceLength(sqRead_trimmed);
    }
  };

  void    resetOverlapIteration(void) {
    ovlStore->setRange(curID = bgnID, endID);
  };

  ~trGlobalData() {
    seqStore->sqStore_close();

    delete [] readData;

    delete seqCache;
    delete ovlStore;

    delete outFile;
  };

  //  Parameters

  uint32             numThreads;

  double             maxErate;
  uint64             memLimit;

  uint32             bgnID;  //  INCLUSIVE range of reads to process.
  uint32             curID;  //    (currently loading id)
  uint32             endID;

  uint32             minReadLength;

  uint32             minOverlapLength;
  uint32             minOverlapCoverage;

  uint32             verboseTrim;
  uint32             verboseAlign;

  //  Statistics

  //  Trimming

  trReadData        *readData;      //  Data for internal use.

  //  Inputs

  char              *seqStoreName;
  sqStore           *seqStore;
  sqCache           *seqCache;

  char              *ovlStoreName;
  ovStore           *ovlStore;

  //  Outputs

  char              *outFileName;
  ovFile            *outFile;

  char              *clearRangesFileName;
};



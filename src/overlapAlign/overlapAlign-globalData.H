
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

    isGarbage      = true;
    isContained    = false;
    isSuspicious   = false;
  };
  ~trReadData() {
  };


public:
  int32     clrBgn;
  int32     clrEnd;

  int32     rawLength;
  int32     trimmedLength;

  bool      isGarbage;
  bool      isContained;
  bool      isSuspicious;
};



class trGlobalData {
public:
  trGlobalData() {
    numThreads          = 1;

    maxErate            = 0.12;
    memLimit            = UINT64_MAX;

    bgnID               = 0;
    curID               = 0;
    endID               = UINT32_MAX;

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
    ovlFile             = NULL;

    outStoreName        = NULL;
    outStore            = NULL;
    outFile             = NULL;

    trimFileName        = NULL;
    trimFile            = NULL;

    correctionFileName  = NULL;
    correctionFile      = NULL;
  };

  void    initialize(void) {

    //sqRead_setDefaultVersion(sqRead_raw);

    //  Open the store, assuming the reads are trimmed already.

    fprintf(stderr, "Opening seqStore '%s'\n", seqStoreName);
    seqStore  = sqStore::sqStore_open(seqStoreName, sqStore_readOnly);

    //  But if no trimmed reads, reopen the store to do trimming.

#if 1
    if (seqStore->sqStore_getNumTrimmedReads() == 0) {
      seqStore->sqStore_close();

      //sqRead_setDefaultVersion(sqRead_raw);

      fprintf(stderr, "Opening seqStore '%s' (for trimming)\n", seqStoreName);
      seqStore  = sqStore::sqStore_open(seqStoreName, sqStore_extend);
    }
#endif

    //  Load all the reads.  Regardless of trim status, we ALWAYS want
    //  to load raw reads, because we ALWAYS need to adjust overlaps
    //  from raw reads to trimmed reads.

    fprintf(stderr, "Loading all reads.\n");
    seqCache  = new sqCache(seqStore, sqRead_raw, memLimit);
    seqCache->sqCache_loadReads();

    //  Open overlaps.

    if (directoryExists(ovlStoreName)) {
      fprintf(stderr, "Reading overlaps from store '%s' and writing to '%s'\n",
              ovlStoreName, outStoreName);

      ovlStore =                  new ovStore      (ovlStoreName, seqStore);
      outStore = (outStoreName) ? new ovStoreWriter(outStoreName, seqStore) : NULL;
    }

    else {
      fprintf(stderr, "Reading overlaps from file '%s' and writing to '%s'\n",
              ovlStoreName, outStoreName);

      ovlFile =                  new ovFile(seqStore, ovlStoreName, ovFileFull);
      outFile = (outStoreName) ? new ovFile(seqStore, outStoreName, ovFileFullWrite) : NULL;
    }

    //  Parameters

    if (bgnID == 0)                                bgnID = 1;
    if (endID  > seqStore->sqStore_getNumReads())  endID = seqStore->sqStore_getNumReads();

    curID = bgnID;

    //  Trimming data.  Initialize clear range and lengths to the full read.

    readData = new trReadData [seqStore->sqStore_getNumReads() + 1];

    for (uint32 ii=0; ii<seqStore->sqStore_getNumReads() + 1; ii++) {
      readData[ii].clrEnd        = seqCache->sqCache_getLength(ii);
      readData[ii].rawLength     = seqCache->sqCache_getLength(ii);
      readData[ii].trimmedLength = seqCache->sqCache_getLength(ii);
    }

    //  Outputs
  };

  void    resetOverlapIteration(void) {
    ovlStore->setRange(curID = bgnID, endID);
  };

  ~trGlobalData() {
    seqStore->sqStore_close();

    delete [] readData;

    delete seqCache;
    delete ovlStore;
    delete ovlFile;

    delete outStore;
    delete outFile;

    AS_UTL_closeFile(trimFile);
    AS_UTL_closeFile(correctionFile);
  };

  //  Parameters

  uint32             numThreads;

  double             maxErate;
  uint64             memLimit;

  uint32             bgnID;  //  INCLUSIVE range of reads to process.
  uint32             curID;  //    (currently loading id)
  uint32             endID;

  uint32             minOverlapLength;
  uint32             minOverlapCoverage;

  uint32             verboseTrim;
  uint32             verboseAlign;

  //  Statistics

  //  Trimming

  trReadData        *readData;

  //  Inputs

  char              *seqStoreName;
  sqStore           *seqStore;
  sqCache           *seqCache;

  char              *ovlStoreName;
  ovStore           *ovlStore;
  ovFile            *ovlFile;

  //  Outputs

  char              *outStoreName;
  ovStoreWriter     *outStore;
  ovFile            *outFile;

  char              *trimFileName;
  FILE              *trimFile;

  char              *correctionFileName;
  FILE              *correctionFile;
};



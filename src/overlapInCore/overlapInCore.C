
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

#include "overlapInCore.H"
#include "strings.H"


uint64  hashFunctions::MASK = uint64zero;
uint32  hashFunctions::HSF1 = 0;
uint32  hashFunctions::HSF2 = 0;
uint32  hashFunctions::SV1  = 0;
uint32  hashFunctions::SV2  = 0;
uint32  hashFunctions::SV3  = 0;


int
main(int argc, char **argv) {
  uint32    numThreads         = 1;

  double    maxErate           = 0.06;
  double    maxAlignErate      = 0.06;
  uint32    minOlapLen         = 0;

  uint32    bgnHashID          = 1;
  uint32    endHashID          = uint32max;

  uint32    bgnRefID           = 1;
  uint32    endRefID           = uint32max;

  uint32    kmerSize           = 22;
  char     *frequentMersPath   = nullptr;

  uint64    filterByKmerCount  = 0;

  uint32    minOverlapLen      = 0;
  uint64    overlapLimit       = uint64max;

  bool      uniqOverlapPerPair = true;
  bool      partialOverlaps    = false;
  bool      checkHopeless      = true;

  uint32    hashBits           = 22;
  uint64    hashDataLenMax     = 100000000;
  double    hashMaxLoad        = 0.6;

  char     *seqStorePath       = nullptr;
  char     *overlapOutPath     = nullptr;
  char     *statsOutPath       = nullptr;

  argc = AS_configure(argc, argv);

  std::vector<char const *>  err;
  for (int arg=1; arg < argc; arg++) {
    if        (strcmp(argv[arg], "-partial") == 0) {
      partialOverlaps = true;

    } else if (strcmp(argv[arg], "-h") == 0) {
      decodeRange(argv[++arg], bgnHashID, endHashID);

    } else if (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], bgnRefID, endRefID);

    } else if (strcmp(argv[arg], "-k") == 0) {
      arg++;

      if (fileExists(argv[arg]) == true)
        frequentMersPath = argv[arg];
      else
        kmerSize = strtouint32(argv[arg]);

    } else if (strcmp(argv[arg], "-l") == 0) {
      overlapLimit = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-m") == 0) {
      uniqOverlapPerPair = false;
    } else if (strcmp(argv[arg], "-u") == 0) {
      uniqOverlapPerPair = true;

    } else if (strcmp(argv[arg], "--hashbits") == 0) {
      hashBits = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "--hashdatalen") == 0) {
      hashDataLenMax = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "--hashload") == 0) {
      hashMaxLoad = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      overlapOutPath = argv[++arg];

    } else if (strcmp(argv[arg], "-s") == 0) {
      statsOutPath = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      numThreads = setNumThreads(argv[++arg]);


    } else if (strcmp(argv[arg], "--minlength") == 0) {
      minOverlapLen = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "--minkmers") == 0) {
      filterByKmerCount = 1;

    } else if (strcmp(argv[arg], "--maxerate") == 0) {
      maxErate = strtodouble(argv[++arg]);
    } else if (strcmp(argv[arg], "--alignnoise") == 0) {
      maxAlignErate = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-z") == 0) {
      checkHopeless = false;

    } else if (seqStorePath == nullptr) {
      seqStorePath = argv[arg];

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  if (filterByKmerCount == 1)      filterByKmerCount = int(floor(exp(-1.0 * (double)kmerSize * maxErate) * (minOverlapLen - kmerSize + 1)));

  if (kmerSize == 0)               err.push_back("Kmer length (-k) must be supplied.\n");
  if (overlapOutPath == nullptr)   err.push_back("No overlap output (-o) supplied.\n");
  if (seqStorePath == nullptr)     err.push_back("No seqStore input (last word in the command line) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [options] <seqStorePath>\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  (( VERY OUT OF DATE ))\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-b <fn>     in contig mode, specify the output file\n");
    fprintf(stderr, "-c          contig mode.  Use 2 frag stores.  First is\n");
    fprintf(stderr, "            for reads; second is for contigs\n");
    fprintf(stderr, "-partial    do partial overlaps\n");
    fprintf(stderr, "-h <range>  to specify fragments to put in hash table\n");
    fprintf(stderr, "            Implies LSF mode (no changes to frag store)\n");
    fprintf(stderr, "-I          designate a file of frag iids to limit olaps to\n");
    fprintf(stderr, "            (Contig mode only)\n");
    fprintf(stderr, "-k          if one or two digits, the length of a kmer, otherwise\n");
    fprintf(stderr, "            the filename containing a list of kmers to ignore in\n");
    fprintf(stderr, "            the hash table\n");
    fprintf(stderr, "-l          specify the maximum number of overlaps per\n");
    fprintf(stderr, "            fragment-end per batch of fragments.\n");
    fprintf(stderr, "-m          allow multiple overlaps per oriented fragment pair\n");
    fprintf(stderr, "-M          specify memory size.  Valid values are '8GB', '4GB',\n");
    fprintf(stderr, "            '2GB', '1GB', '256MB'.  (Not for Contig mode)\n");
    fprintf(stderr, "-o          specify output file name\n");
    fprintf(stderr, "-P          write protoIO output (if not -partial)\n");
    fprintf(stderr, "-r <range>  specify old fragments to overlap\n");
    fprintf(stderr, "-t <n>      use <n> parallel threads\n");
    fprintf(stderr, "-u          allow only 1 overlap per oriented fragment pair\n");
    fprintf(stderr, "-z          skip the hopeless check (also skipped at > 0.06)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--maxerate <n>     only output overlaps with fraction <n> or less error (e.g., 0.06 == 6%%)\n");
    fprintf(stderr, "--alignnoise <n>\n");
    fprintf(stderr, "--minlength <n>    only output overlaps of <n> or more bases\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--hashbits n       Use n bits for the hash mask.\n");
    fprintf(stderr, "--hashdatalen n    Load at most n bytes into the hash table at one time.\n");
    fprintf(stderr, "--hashload f       Load to at most 0.0 < f < 1.0 capacity (default 0.7).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "--readsperbatch n  Force batch size to n.\n");
    fprintf(stderr, "--readsperthread n Force each thread to process n reads.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  //
  //  Open inputs.  We need this early so we can get the number of reads to
  //  check read ranges.
  //

  sqStore        *readStore = new sqStore(seqStorePath);
  sqCache        *readCache = new sqCache(readStore);

  if (bgnHashID == 0)                                 bgnHashID = 1;
  if (bgnHashID  > readStore->sqStore_lastReadID())   bgnHashID = readStore->sqStore_lastReadID();

  if (endHashID == 0)                                 endHashID = 1;
  if (endHashID  > readStore->sqStore_lastReadID())   endHashID = readStore->sqStore_lastReadID();

  if (bgnRefID  == 0)                                 bgnRefID  = 1;
  if (bgnRefID   > readStore->sqStore_lastReadID())   bgnRefID  = readStore->sqStore_lastReadID();

  if (endRefID  == 0)                                 endRefID  = 1;
  if (endRefID   > readStore->sqStore_lastReadID())   endRefID  = readStore->sqStore_lastReadID();

  if (endHashID < bgnHashID) {
    fprintf(stderr, "ERROR: Invalid hash id range bgn=%u > end=%u.\n", bgnHashID, endHashID);
    return(1);
  }

  if (endRefID < bgnRefID) {
    fprintf(stderr, "ERROR: Invalid reference id range bgn=%u > end=%u.\n", bgnRefID, endRefID);
    return(1);
  }

  //if (endHashID < bgnRefID) {
  //  fprintf(stderr, "ERROR: overlaps only computed for reference ID <= hash ID.\n");
  //  fprintf(stderr, "ERROR: supplied ranges will result in no overlaps being computed.\n");
  //}

  //
  //  Update some of the parameters.
  //

  if (maxErate > 0.06)
    checkHopeless = false;



  //  
  //  Construct all the objects we need, and set parameters in them.
  //

  omp_set_num_threads(numThreads);

  fprintf(stderr, "Opening output '%s'.\n", overlapOutPath);
  ovFile     *OF = new ovFile(readStore, overlapOutPath, ovFileFullWrite);

  //fprintf(stderr, "Creating %u work areas.\n", numThreads);
  //workArea   *WA = new workArea [numThreads];

  //#pragma omp parallel for    //  Does this still need threads?
  //for (uint32 i=0;  i<G.Num_PThreads;  i++)
  //  WA[i].initialize(i, readStore, readCache);

  fprintf(stderr, "Loading reference reads %u-%u inclusive.\n", bgnRefID, endRefID);
  readCache->sqCache_loadReads(bgnRefID, endRefID, true);

  //
  //  Loop over the hash reads.
  //   - load as many as possible into the hash table.  the last read loaded
  //     is returned into 'endLoaded'.  at the end of the loop we use this to
  //     start the next loop.
  //   - search for overlaps
  //   - repeat until all hash reads are processed
  //

  overlapHitStats  statistics;

  for (uint32 bgnLoaded=bgnHashID; bgnLoaded < endHashID; ) {
    fprintf(stderr, "Creating hash table.\n");
    hashTable  *HT = new hashTable(kmerSize, hashBits, hashMaxLoad, minOverlapLen);

    fprintf(stderr, "Loading hash table.\n");
    uint32 endLoaded = HT->loadTable(readStore, bgnLoaded, endHashID, frequentMersPath);

    //  Set up batches for searching.  Each thread is told to process
    //  'perThread' reads at a time.  When they finish that batch, they'll
    //  grab another batch to work on.

    uint32  numBlocks = numThreads * 16;
    uint32  perThread = (endRefID - bgnRefID) / numBlocks + 1;
    uint32 *thrBgn    = new uint32 [numBlocks * 16];
    uint32 *thrEnd    = new uint32 [numBlocks * 16];
    uint32  thrLen    = 0;

    //  Initialize each thread, reset the current position.  curRefID and endRefID are updated, this
    //  cannot be done in the parallel loop!

    for (uint32 bb=bgnRefID; bb <= endRefID; bb += perThread) {
      thrBgn[thrLen] = bb;
      thrEnd[thrLen] = bb + perThread - 1;

      fprintf(stderr, "BLOCK %2u computes ref reads %u - %u\n", thrLen, thrBgn[thrLen], thrEnd[thrLen]);

      thrLen++;
    }

#pragma omp parallel for
    for (uint32 tt=0; tt<thrLen; tt++) {
      uint32     basesLen = 0;
      uint32     basesMax = 0;
      char      *bases    = nullptr;
      workArea  *WA       = new workArea(readStore, readCache, partialOverlaps, maxErate, maxAlignErate);

      fprintf(stderr, "Batch %03u processes reads %9u - %9u.\n", tt, thrBgn[tt], thrEnd[tt]);

      for (uint32 fi=thrBgn[tt]; fi <= thrEnd[tt]; fi++) {
        uint32  readLen = readCache->sqCache_getLength(fi);

        if (readLen < minOverlapLen)
          continue;

        readCache->sqCache_getSequence(fi, bases, basesLen, basesMax);

        for (uint32 i=0; i<readLen; i++)
          bases[i] = tolower(bases[i]);

        WA->findOverlaps   (bases, readLen, fi, true,  HT);
        WA->processOverlaps(bases, readLen, fi, true);

        reverseComplementSequence(bases, readLen);

        WA->findOverlaps   (bases, readLen, fi, false, HT);
        WA->processOverlaps(bases, readLen, fi, false);
      }

#pragma omp critical
      statistics += WA->outputOverlaps(tt, thrBgn[tt], thrEnd[tt], OF);

#if 0
      //fprintf(stderr, "Batch %03u writes    reads " F_U32 "-" F_U32 " (" F_U64 " overlaps " F_U64 "/" F_U64 "/" F_U64 " kmer hits with/without overlap/skipped)\n",
      //        tt, thrBgn[tt], thrEnd[tt],
      //        WA->overlapsLen,
      //        WA->statistics.hitsWithOverlap, WA->statistics.hitsWithoutOverlap, WA->statistics.hitsSkipped);

      {
        for (int zz=0; zz<WA->overlapsLen; zz++)
          OF->writeOverlap(WA->overlaps + zz);
        WA->overlapsLen = 0;

        statistics += WA->statistics;
      }
#endif

      delete [] bases;
      delete    WA;
    }

    //  All done with this hash table.

    delete HT;

    //  Move the hash range to the next batch.

    bgnLoaded = endLoaded + 1;
  }

  //  Close outputs and inputs.

  delete OF;
  delete readCache;
  delete readStore;

  //
  ////////////////////////////////////////
  //

  FILE *stats = stderr;

  if (statsOutPath != NULL) {
    errno = 0;
    stats = fopen(statsOutPath, "w");
    if (errno) {
      fprintf(stderr, "WARNING: failed to open '%s' for writing: %s\n", statsOutPath, strerror(errno));
      stats = stderr;
    }
  }

  fprintf(stats, " Kmer hits without olaps = " F_S64 "\n", statistics.hitsWithoutOverlap);
  fprintf(stats, "    Kmer hits with olaps = " F_S64 "\n", statistics.hitsWithOverlap);
  fprintf(stats, "  Multiple overlaps/pair = " F_S64 "\n", statistics.multipleOverlap);
  fprintf(stats, " Total overlaps produced = " F_S64 "\n", statistics.totalOverlaps);
  fprintf(stats, "      Contained overlaps = " F_S64 "\n", statistics.containedOverlaps);
  fprintf(stats, "       Dovetail overlaps = " F_S64 "\n", statistics.dovetailOverlaps);
  //fprintf(stats, "Rejected by short window = " F_S64 "\n", statistics.Bad_Short_Window_Ct);
  //fprintf(stats, " Rejected by long window = " F_S64 "\n", statistics.Bad_Long_Window_Ct);

  AS_UTL_closeFile(stats, statsOutPath);

  fprintf(stderr, "Bye.\n");

  return(0);
}

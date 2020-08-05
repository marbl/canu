
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
#include "sqCache.H"
#include "ovStore.H"
#include "tgStore.H"

#include "intervalList.H"

#include "sequence.H"

#include "falconConsensus.H"

#include <set>

using namespace std;

//  Define this to recreate the falconConsensus object for each read.
//  This allows the precise memory size needed to process each read to be reported.
//  Performance degradation is severe.
//#define CHECK_MEMORY


//  Duplicated in generateCorrectionLayouts.C
void
loadReadList(char const *readListName, uint32 iidMin, uint32 iidMax, set<uint32> &readList) {
  char  L[1024];

  if (readListName == NULL)
    return;

  fprintf(stderr, "-- Loading list of reads to process from '%s'.\n", readListName);

  //  To give the list some size - if the read list has no elements between iidMin and iidMax,
  //  nothing is inserted below, and then we _think_ no read list was supplied, and try to
  //  process every read, when in fact we should be processing no reads.

  readList.insert(0);

  FILE *R = AS_UTL_openInputFile(readListName);

  for (fgets(L, 1024, R);
       feof(R) == false;
       fgets(L, 1024, R)) {
    uint32  id = strtouint32(L);

    if ((iidMin <= id) &&
        (id     <= iidMax))
      readList.insert(id);
  }

  AS_UTL_closeFile(R, readListName);
}




sqRead *
loadReadData(uint32                     readID,
             sqStore                   *seqStore,
             map<uint32, sqRead *>     &reads) {

  if (reads.count(readID) == 0) {
    reads[readID] = new sqRead;
    seqStore->sqStore_getRead(readID, reads[readID]);
  }

  assert(reads[readID] != NULL);

  return(reads[readID]);
}



void
generateFalconConsensus(falconConsensus           *fc,
                        tgTig                     *layout,
                        sqCache                   *seqCache,
                        map<uint32, sqRead *>     &reads,
                        bool                       trimToAlign,
                        uint32                     minOlapLength) {

  //  What rolls down stairs
  //  alone or in pairs,
  //  rolls over your neighbor's dog?
  //  What's great for a snack,
  //  And fits on your back?
  //  It's log, log, log!

  fprintf(stdout, "%8u %7u %8u", layout->tigID(), layout->length(), layout->numberOfChildren());

  //  Parse the layout and push all the sequences onto our seqs vector.  The first 'evidence'
  //  sequence is the read we're trying to correct.

  falconInput   *evidence = new falconInput [layout->numberOfChildren() + 1];

  uint32         seqLen   = 0;
  uint32         seqMax   = 1048576;
  char          *seq      = new char [seqMax];

  evidence[0].addInput(layout->tigID(),
                       seqCache->sqCache_getSequence(layout->tigID(), seq, seqLen, seqMax),
                       seqCache->sqCache_getLength(layout->tigID()),
                       0,
                       seqCache->sqCache_getLength(layout->tigID()));

  for (uint32 cc=0; cc<layout->numberOfChildren(); cc++) {
    tgPosition  *child = layout->getChild(cc);

    //  Grab a copy of the sequence.

    seqCache->sqCache_getSequence(child->ident(), seq, seqLen, seqMax);

    //  Now screw up the sequence by reverse-complementing and trimming it.

    if (child->isReverse())
      reverseComplementSequence(seq, seqLen);

    uint32  b = 0;
    uint32  e = seqLen;

    if (trimToAlign) {
      b += child->askip();
      e -= child->bskip();
    }

    seq[e] = 0;

    //  Save the read if it is larger than the minimum overlap length.  Anything smaller than this will have zero chance of aligning.

    if (minOlapLength <= e - b)
      evidence[cc+1].addInput(child->ident(), seq + b, e - b, child->min(), child->max());
  }

  delete [] seq;

  //  Loaded all reads, build consensus.

  falconData  *fd = fc->generateConsensus(evidence, layout->numberOfChildren() + 1);

  //  Find the largest stretch of uppercase sequence.  Lowercase sequence denotes MSA coverage was below minOutputCoverage.

  uint32  bgn = 0;
  uint32  end = 0;
  uint32  nrg = 0;

  for (uint32 in=0, bb=0, ee=0; ee<fd->len; ee++) {
    bool   isLower = (('a' <= fd->seq[ee]) && (fd->seq[ee] <= 'z'));
    bool   isLast  = (ee == fd->len - 1);

    if ((in == true) && (isLower || isLast)) {     //  Report the regions we could be saving.
      fprintf(stdout, " %6u-%-6u", bb, ee + isLast);
      nrg++;
    }

    if (isLower) {                                 //  If lowercase, declare that we're not in a
      in = 0;                                      //  good region any more.
    }

    else if (in == 0) {                            //  Otherwise, if not in a region (so the first
      bb = ee;                                     //  uppercase), remember the coordinate and
      in = 1;                                      //  switch to being 'in' a region.
    }

    if ((in == 1) && (ee + 1 - bb > end - bgn)) {  //  If 'in' a good region, remember the longest.
      bgn = bb;                                    //  'ee + 1': if the next letter is lower case
      end = ee + 1;                                //  our bgn,end interval will be set in this iteration
    }
  }

  if (nrg == 0)
    fprintf(stdout, " %6u-%-6u", 0, 0);

  uint32 len = 0;
  uint64 mem = 0;

  fc->analyzeLength(layout, len, mem);

  fprintf(stdout, "(%6u) memory act %10lu est %10lu act/est %.2f", len, fc->getRSS(), mem, fc->getRSS() * 100.0 / mem);
  fprintf(stdout, "\n");

  //  Note where in the full corrected read the output corrected read came from.

  layout->_trimBgn   = (end == 0) ? (0) : (fd->pos[bgn]);            //  Space based (probably).
  layout->_trimEnd   = (end == 0) ? (0) : (fd->pos[end - 1] + 1);    //  Space based.

  //  Update the layout with consensus sequence, positions, et cetera.
  //  If the whole string is lowercase (grrrr!) then bgn == end == 0.

  resizeArrayPair(layout->_bases, layout->_quals, layout->_basesLen, layout->_basesMax, end - bgn + 1, resizeArray_doNothing);

  for (uint32 ii=bgn; ii<end; ii++) {
    layout->_bases[ii-bgn] = fd->seq[ii];
    layout->_quals[ii-bgn] = fd->eqv[ii];
  }

  layout->_layoutLen = end - bgn;
  layout->_basesLen  = end - bgn;

  layout->_bases[layout->_basesLen] = 0;
  layout->_quals[layout->_basesLen] = 0;

  //  One could dump bases and quals here, if so desired.

  ;

  //  Clean up.  Remvoe all the reads[] we've loaded.

  for (map<uint32, sqRead     *>::iterator it=reads.begin(); it != reads.end(); ++it)
    delete it->second;

  reads.clear();

  delete    fd;
  delete [] evidence;
}




int
main(int argc, char **argv) {
  char const       *seqName   = 0L;
  char const       *corName   = 0L;
  uint32            corVers   = 1;

  char const       *exportName = NULL;
  char const       *importName = NULL;

  uint32            errorRate = AS_OVS_encodeEvalue(0.015);

  char const       *outputPrefix = NULL;
  bool              outputCNS    = false;
  bool              outputFASTQ  = false;
  bool              outputLog    = false;

  uint64            memoryLimit = 0;
  uint64            memPerRead  = 0;
  uint32            batchLimit  = 0;
  uint32            readLimit   = 0;

  uint32            idMin = 1;
  uint32            idMax = UINT32_MAX;
  char const       *readListName = NULL;
  set<uint32>       readList;

  uint32            numThreads         = omp_get_max_threads();

  uint32            minOutputCoverage  = 4;
  uint32            minOutputLength    = 1000;
  double            minOlapIdentity    = 0.5;
  double            minOlapLength      = 500;

  bool              trimToAlign        = true;
  bool              restrictToOverlap  = true;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {   //  INPUTS
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      corName = argv[++arg];


    } else if (strcmp(argv[arg], "-p") == 0) {   //  OUTPUTS
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-cns") == 0) {
      outputCNS = true;

    } else if (strcmp(argv[arg], "-fastq") == 0) {
      outputFASTQ = true;

    } else if (strcmp(argv[arg], "-log") == 0) {
      outputLog = true;

    } else if (strcmp(argv[arg], "-partition") == 0) {
      memoryLimit = (uint64)(strtodouble(argv[++arg]) * 1024 * 1024 * 1024);
      memPerRead  = (uint64)(strtodouble(argv[++arg]) * 1024 * 1024 * 1024);
      batchLimit  = strtouint32(argv[++arg]);
      readLimit   = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {   //  COMPUTE RESOURCES
      numThreads = strtouint32(argv[++arg]);


    } else if (strcmp(argv[arg], "-f") == 0) {   //  ALGORITHM OPTIONS
      restrictToOverlap = false;


    } else if (strcmp(argv[arg], "-R") == 0) {   //  READ SELECTION
      readListName = argv[++arg];

    } else if (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], idMin, idMax);


    } else if (strcmp(argv[arg], "-cc") == 0) {   //  CONSENSUS
      minOutputCoverage = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-cl") == 0) {
      minOutputLength = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-oi") == 0) {
      minOlapIdentity = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-ol") == 0) {
      minOlapLength = strtodouble(argv[++arg]);


    } else if (strcmp(argv[arg], "-export") == 0) {   //  DEBUGGING
      exportName = argv[++arg];

    } else if (strcmp(argv[arg], "-import") == 0) {
      importName = argv[++arg];


    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if ((seqName == NULL) && (importName == NULL))
    err.push_back("ERROR: no seqStore input (-S) supplied.\n");

  if ((corName == NULL) && (importName == NULL))
    err.push_back("ERROR: no corStore input (-C) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S seqStore -O ovlStore ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS (all mandatory)\n");
    fprintf(stderr, "  -S seqStore        mandatory path to seqStore\n");
    fprintf(stderr, "  -C corStore        mandatory path to corStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUTS:\n");
    fprintf(stderr, "  -p prefix          output filename prefix\n");
    fprintf(stderr, "  -cns               enable primary output (to 'prefix.cns')\n");
    fprintf(stderr, "  -fastq             enable fastq output (to 'prefix.fastq')\n");
    fprintf(stderr, "  -log               enable (debug) logging output (to 'prefix.log')\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "RESOURCE PARAMETERS:\n");
    fprintf(stderr, "  -t numThreads      number of compute threads to use (default: all)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "ALGORITHM PARAMETERS:\n");
    fprintf(stderr, "  -f                 align evidence to the full read, ignore overlap position\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "READ SELECTION:\n");
    fprintf(stderr, "  -R readsToCorrect  only process reads listed in file 'readsToCorrect'\n");
    fprintf(stderr, "  -r bgn[-end]       only process reads from ID 'bgn' to 'end' (inclusive)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CONSENSUS PARAMETERS:\n");
    fprintf(stderr, "  -cc coverage       output:   minimum consensus coverage needed call a corrected base\n");
    fprintf(stderr, "  -cl length         output:   minimum length of corrected region to output as a corrected read\n");
    fprintf(stderr, "  -oi identity       evidence: minimum identity of an aligned evidence read overlap\n");
    fprintf(stderr, "  -ol length         evidence: minimum length   of an aligned evidence read overlap\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "PARTITIONING SUPPORT:\n");
    fprintf(stderr, "  -partition M m R   configure jobs to fit in M GB memory with not more than R reads per batch,\n");
    fprintf(stderr, "                     allowing m GB memory for processing.  write output to 'prefix.batches'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DEBUGGING SUPPORT:\n");
    fprintf(stderr, "  -export name       write the data used for the computation to file 'name'\n");
    fprintf(stderr, "  -import name       compute using the data in file 'name'\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  omp_set_num_threads(numThreads);

  //  Probably not needed, as sqCache explicitly loads only sqRead_raw, but
  //  setting the default version guarantees that we access only 'raw' reads.

  sqRead_setDefaultVersion(sqRead_raw);

  //  Open inputs.

  sqStore *seqStore = NULL;
  sqCache *seqCache = NULL;
  tgStore *corStore = NULL;

  if (seqName) {
    fprintf(stderr, "-- Opening seqStore '%s'.\n", seqName);
    seqStore = new sqStore(seqName);
    seqCache = new sqCache(seqStore, sqRead_raw);
  }

  if (corName) {
    fprintf(stderr, "-- Opening corStore '%s' version %u.\n", corName, corVers);
    corStore = new tgStore(corName, corVers);
  }

  if ((seqStore) &&
      (seqStore->sqStore_lastReadID() < idMax))        //  Limit the range of processing to the
    idMax = seqStore->sqStore_lastReadID();            //  number of reads in the store.

  loadReadList(readListName, idMin, idMax, readList);   //  Further limit to a set of good reads.

  //  Open any import or export files.

  writeBuffer *exportFile = NULL;
  readBuffer  *importFile = NULL;

  if (exportName) {
    fprintf(stderr, "-- Opening export file '%s'.\n", exportName);
    exportFile = new writeBuffer(exportName, "w");
  }

  if (importName) {
    fprintf(stderr, "-- Opening import file '%s'.\n", importName);
    importFile  = new readBuffer(importName);
  }

  //  Open logging and summary files

  FILE *logFile = NULL;
  FILE *cnsFile = NULL;
  FILE *seqFile = NULL;
  FILE *batFile = NULL;

  cnsFile = AS_UTL_openOutputFile(outputPrefix, '.', "cns",     outputCNS);
  seqFile = AS_UTL_openOutputFile(outputPrefix, '.', "fastq",   outputFASTQ);
  logFile = AS_UTL_openOutputFile(outputPrefix, '.', "log",     outputLog);
  batFile = AS_UTL_openOutputFile(outputPrefix, '.', "batches", memoryLimit > 0);

  //  Initialize processing.
  //
  //  One might be tempted to NOT create a seqCache or a falconConsensus object if we're only
  //  partitioning, but that would be wrong, because partitioning uses these objects to determine
  //  the base amount of memory needed.

  falconConsensus           *fc = new falconConsensus(minOutputCoverage, minOutputLength, minOlapIdentity, minOlapLength, restrictToOverlap);
  map<uint32, sqRead *>      reads;

  if (memoryLimit == 0) {
    fprintf(stdout, "    read    read evidence     corrected\n");
    fprintf(stdout, "      ID  length    reads       regions\n");
    fprintf(stdout, "-------- ------- -------- ------------- ...\n");
  }

  //
  //  If input from a package file, load and process data until there isn't any more.
  //

  if (importFile) {
    tgTig                     *layout = new tgTig();

    FILE  *importedLayouts = AS_UTL_openOutputFile(importName, '.', "layout", (importName != NULL));
    FILE  *importedReads   = AS_UTL_openOutputFile(importName, '.', "fasta",  (importName != NULL));

    while (layout->importData(importFile, reads, NULL, NULL) == true) {
      generateFalconConsensus(fc,
                              layout,
                              seqCache,
                              reads,
                              trimToAlign,
                              minOlapLength);

      if (cnsFile)
        layout->saveToStream(cnsFile);

      if (seqFile)
        layout->dumpFASTQ(seqFile);

      delete layout;
      layout = new tgTig();    //  Next loop needs an existing empty layout.
    }

    AS_UTL_closeFile(importedReads);
    AS_UTL_closeFile(importedLayouts);

    //  We properly should clean up the data loaded.

    delete layout;
  }

  //
  //  Otherwise, if we're just dumping data, just dump the data without processing.
  //

  else if (exportFile) {
    for (uint32 ii=idMin; ii<=idMax; ii++) {
      if ((readList.size() > 0) &&      //  Skip reads not on the read list,
          (readList.count(ii) == 0))    //  if there actually is a read list.
        continue;

      tgTig *layout = corStore->loadTig(ii);

      if (layout) {
        fprintf(stdout, "%8u %7u %8u", layout->tigID(), layout->length(), layout->numberOfChildren());

        layout->exportData(exportFile, seqStore, true);
        corStore->unloadTig(layout->tigID());

        fprintf(stdout, "        DUMPED\n");
      }
    }
  }

  //
  //  If a memory limit, set up partitions.
  //
  //  This computes spands of reads such that the memory needed to load
  //  all the overlapping reads is less than some limit.
  //

  else if (memoryLimit > 0) {
    uint32   lastID   = seqStore->sqStore_lastReadID();
    uint32  *readLens = new uint32 [lastID + 1];
    uint32  *readRefs = new uint32 [lastID + 1];

    //  Load read lengths, convert to an approximate size they'll use when loaded, and initialize references to zero.
    //
    //  It's not ideal, since we use lots of insider knowledge.
    //    12          - chunk header, chunk length, possibly length of data
    //    Length / 4  - 2-bit encoded bases
    //    4           - padding on chunk
    //    cacheEntry  - storage internal to the cache.

    for (uint32 ii=1; ii <= lastID; ii++) {
      readLens[ii] = 12 + seqStore->sqStore_getReadLength(ii, sqRead_raw) / 4 + 4 + sizeof(sqCacheEntry);   //  Round up, and 3 extra uint32.
      readRefs[ii] = 0;
    }

    //  The user is requesting batchLimit batches with at least readLimit reads per batch.
    //  Further, each batch can use no more than memoryLimit GB, assuming memPerRead GB to actually compute the corrected read.

    uint32  readsPerBatch = lastID / batchLimit + 1;

    if (readList.size() > 0)
      readsPerBatch = readList.size() / batchLimit + 1;

    if (readsPerBatch < readLimit)
      readsPerBatch = readLimit;

    //  Analyze each layout, remembering how much memory is needed.

    uint64   memUsedBase = getBytesAllocated();  //  For seqCache, falconConsensus and misc gunk.
    uint64   memUsed     = memUsedBase;
    uint32   nReads      = 0;
    uint32   batchNum    = 1;
    uint32   bgnID       = idMin;

    if (memUsedBase + memPerRead > memoryLimit) {
      fprintf(stderr, "\n");
      fprintf(stderr, "ERROR:  Need at least M=%6.3f GB (with m=%6.3f GB) to compute corrections.\n",
              (memUsedBase + memPerRead) / 1024.0 / 1024.0 / 1024.0,
              memPerRead  / 1024.0 / 1024.0 / 1024.0);
      exit(1);
    }

    //fprintf(stderr, "readsPerBatch %u\n", readsPerBatch);
    //fprintf(stderr, "memoryLimit   %f GB\n", memoryLimit / 1024.0 / 1024.0 / 1024.0);

    fprintf(batFile, "batch     bgnID     endID  nReads  memory (base memory %.3f GB)\n", memUsedBase / 1024.0 / 1024.0 / 1024.0);
    fprintf(batFile, "----- --------- --------- ------- -------\n");

    for (uint32 ii=idMin; ii<=idMax; ii++) {
      if ((readList.size() > 0) &&      //  Skip reads not on the read list,
          (readList.count(ii) == 0))    //  if there actually is a read list.
        continue;

      tgTig *layout = corStore->loadTig(ii);

      if (layout == NULL)
        continue;

      //  Compute how much memory this tig needs needs to store it's reads.
      //  This is an overestimate as it includes singleton reads.  Correctly
      //  accounting for not loading singleton reads will be tricky because
      //  removing earlier tigs could turn reads to singletons.

      uint64   memAdded = readLens[ii];

      readRefs[ii]++;

      for (uint32 cc=0; cc<layout->numberOfChildren(); cc++) {
        tgPosition  *child = layout->getChild(cc);
        uint32       rdID  = child->ident();

        if (readRefs[rdID] == 0)
          memAdded += readLens[rdID];

        readRefs[rdID]++;
      }

      corStore->unloadTig(layout->tigID());

      //  If we're over the limit, report the range and reset.

      if ((memUsed + memAdded > memoryLimit) ||
          (nReads + 1 > readsPerBatch)) {
        fprintf(batFile, "%5u %9u %9u %7u %7.3f\n", batchNum, bgnID, ii-1, nReads, memUsed / 1024.0 / 1024.0 / 1024.0);
        batchNum += 1;
        bgnID     = ii;
        memUsed   = memUsedBase;
        nReads    = 0;

        for (uint32 ii=0; ii <= lastID; ii++)
          readRefs[ii] = 0;
      }

      memUsed += memAdded;
      nReads  += 1;
    }

    //  And one final report for the last block.

    fprintf(batFile, "%5u %9u %9u %7u %7.3f\n", batchNum, bgnID, idMax, nReads, memUsed / 1024.0 / 1024.0 / 1024.0);

    delete [] readRefs;
    delete [] readLens;
  }

  //
  //  Otherwise, load and process from a store, the usual processing loop.
  //

  else {

    //  First, scan all tigs we're going to process and count the number
    //  of times we need each read.  The sqCache can then figure out what
    //  reads to cache, and what reads to load on demand.

    map<uint32,uint32>   readsToLoad;

    for (uint32 ii=idMin; ii<=idMax; ii++) {
      if ((readList.size() > 0) &&      //  Skip reads not on the read list,
          (readList.count(ii) == 0))    //  if there actually is a read list.
        continue;

      tgTig *layout = corStore->loadTig(ii);

      if (layout) {
        readsToLoad[ii]++;

        for (uint32 cc=0; cc<layout->numberOfChildren(); cc++)
          readsToLoad[layout->getChild(cc)->ident()]++;
      }
    }

    seqCache->sqCache_loadReads(readsToLoad);

    //  Now, with all (most) of the read sequences loaded, process.

#ifdef CHECK_MEMORY
    delete fc;
    fc = NULL;
#endif

    for (uint32 ii=idMin; ii<=idMax; ii++) {
      if ((readList.size() > 0) &&      //  Skip reads not on the read list,
          (readList.count(ii) == 0))    //  if there actually is a read list.
        continue;

      tgTig *layout = corStore->loadTig(ii);

      if (layout) {
#ifdef CHECK_MEMORY
        fc = new falconConsensus(minOutputCoverage, minOutputLength, minOlapIdentity, minOlapLength, restrictToOverlap);
#endif

        generateFalconConsensus(fc,
                                layout,
                                seqCache,
                                reads,
                                trimToAlign,
                                minOlapLength);

#ifdef CHECK_MEMORY
        delete fc;
        fc = NULL;
#endif

        if (cnsFile)
          layout->saveToStream(cnsFile);

        if (seqFile)
          layout->dumpFASTQ(seqFile);

        corStore->unloadTig(layout->tigID());
      }
    }
  }

  //  Close files and clean up.

  AS_UTL_closeFile(logFile);
  AS_UTL_closeFile(cnsFile);
  AS_UTL_closeFile(seqFile);
  AS_UTL_closeFile(batFile);

  delete    exportFile;
  delete    importFile;

  delete    fc;
  delete    corStore;

  delete    seqCache;

  delete seqStore;

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  return(0);
}

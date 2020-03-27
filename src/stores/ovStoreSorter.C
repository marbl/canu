
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
#include "ovStore.H"
#include "ovStoreConfig.H"

#include <algorithm>
using namespace std;


void
checkSentinel(const char *ovlName, uint32 sliceNum, ovStoreConfig *config, bool forceRun) {
  char   N[FILENAME_MAX+1];

  //  Check if the user is a moron.

  if ((sliceNum == 0) ||
      (sliceNum > config->numSlices())) {
    fprintf(stderr, "No slice " F_U32 " exists; only slices 1-" F_U32 " exist.\n", sliceNum, config->numSlices());
    exit(1);
  }

  //  Check if done.

  snprintf(N, FILENAME_MAX, "%s/%04u.started", ovlName, sliceNum);

  if ((forceRun == false) &&
      (fileExists(N) == true)) {
    fprintf(stderr, "Job (appears to be) in progress; sentinel file '%s' exists.\n", N);
    exit(1);
  }

  //  Not done and not running, so create a sentinel to say we're running.

  AS_UTL_createEmptyFile(N);
}



void
removeSentinel(const char *ovlName, uint32 sliceNum) {
  char   N[FILENAME_MAX+1];

  snprintf(N, FILENAME_MAX, "%s/%04u.started", ovlName, sliceNum);

  AS_UTL_unlink(N);
}



void
checkMemory(const char *ovlName, uint32 sliceNum, uint64 totOvl, uint64 maxMemory) {

  if (ovOverlapSortSize * totOvl > maxMemory) {
    fprintf(stderr, "ERROR:  Overlaps need %.2f GB memory, but process limited (via -M) to " F_U64 " GB.\n",
            ovOverlapSortSize * totOvl / 1024.0 / 1024.0 / 1024.0, maxMemory >> 30);
    removeSentinel(ovlName, sliceNum);
    exit(1);
  }

  fprintf(stderr, "\n");

  double  memUsed    = totOvl * ovOverlapSortSize / 1024.0 / 1024.0 / 1024.0;
  uint64  memAllowed = maxMemory >> 30;

  if (maxMemory == UINT64_MAX)
    fprintf(stderr, "Loading %10" F_U64P " overlaps using %.2f GB memory.\n", totOvl, memUsed);
  else
    fprintf(stderr, "Loading %10" F_U64P " overlaps using %.2f GB of allowed (-M) " F_U64 " GB memory.\n", totOvl, memUsed, memAllowed);
}







int
main(int argc, char **argv) {
  char const     *ovlName      = NULL;
  char const     *seqName      = NULL;
  char const     *cfgName      = NULL;
  uint32          sliceNum     = UINT32_MAX;

  uint64          maxMemory    = UINT64_MAX;

  bool            deleteIntermediateEarly = false;
  bool            deleteIntermediateLate  = false;
  bool            forceRun = false;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      cfgName = argv[++arg];

    } else if (strcmp(argv[arg], "-s") == 0) {
      sliceNum  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-M") == 0) {
      maxMemory  = (uint64)ceil(atof(argv[++arg]) * 1024.0 * 1024.0 * 1024.0);

    } else if (strcmp(argv[arg], "-deleteearly") == 0) {
      deleteIntermediateEarly = true;

    } else if (strcmp(argv[arg], "-deletelate") == 0) {
      deleteIntermediateLate  = true;

    } else if (strcmp(argv[arg], "-f") == 0) {
      forceRun = true;

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (ovlName == NULL)
    err.push_back("ERROR: No overlap store (-O) supplied.\n");

  if (sliceNum == UINT32_MAX)
    err.push_back("ERROR: no slice number (-F) supplied.\n");

  if (maxMemory < OVSTORE_MEMORY_OVERHEAD + ovOverlapSortSize)
    fprintf(stderr, "ERROR: Memory (-M) must be at least 0.25 GB to account for overhead.\n");  //  , OVSTORE_MEMORY_OVERHEAD / 1024.0 / 1024.0 / 1024.0

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -O asm.ovlStore -S asm.seqStore -C ovStoreConfig -s slice [opts]\n", argv[0]);
    fprintf(stderr, "  -O asm.ovlStore       path to overlap store to create\n");
    fprintf(stderr, "  -S asm.seqStore       path to sequence store\n");
    fprintf(stderr, "  -C config             path to ovStoreConfig configuration file\n");
    fprintf(stderr, "  -s slice              slice to process (1 ... N)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -M m             maximum memory to use, in gigabytes\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -deleteearly     remove intermediates as soon as possible (unsafe)\n");
    fprintf(stderr, "  -deletelate      remove intermediates when outputs exist (safe)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -f               force a recompute, even if the output exists or appears in progress\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  //  Load the config.

  ovStoreConfig  *config = new ovStoreConfig(cfgName);

  //  Check if the sentinel exists (and if the user is a moron).

  checkSentinel(ovlName, sliceNum, config, forceRun);

  //  Not done.  Let's go!

  sqStore             *seq    = new sqStore(seqName);
  ovStoreSliceWriter  *writer = new ovStoreSliceWriter(ovlName, seq, sliceNum, config->numSlices(), config->numBuckets());

  //  Get the number of overlaps in each bucket slice.

  fprintf(stderr, "\n");
  fprintf(stderr, "Finding overlaps.\n");

  uint64 *bucketSizes = new uint64 [config->numBuckets() + 1];     //  The number of overlaps in bucket i
  uint64  totOvl      = writer->loadBucketSizes(bucketSizes);

  //  Fail if we don't have enough memory to process.

  checkMemory(ovlName, sliceNum, totOvl, maxMemory);

  //  Allocatge space for overlaps, and load them.

  ovOverlap *ovls    = new ovOverlap [totOvl];
  uint64     ovlsLen = 0;

  for (uint32 bb=0; bb<=config->numBuckets(); bb++)
    writer->loadOverlapsFromBucket(bb, bucketSizes[bb], ovls, ovlsLen);

  //  Check that we found all the overlaps we were expecting.

  if (ovlsLen != totOvl) {
    fprintf(stderr, "ERROR: read " F_U64 " overlaps, expected " F_U64 "\n", ovlsLen, totOvl);
    exit(1);
  }

  //  Clean up space if told to.

  if (deleteIntermediateEarly)
    writer->removeOverlapSlice();

  //  Sort the overlaps!  Finally!  The parallel STL sort is NOT inplace, and blows up our memory.

  fprintf(stderr, "\n");
  fprintf(stderr, "Sorting.\n");

  sort(ovls, ovls + ovlsLen);

  //  Output to the store.

  fprintf(stderr, "\n");   //  Sorting has no output, so this would generate a distracting extra newline
  fprintf(stderr, "Writing sorted overlaps.\n");

  writer->writeOverlaps(ovls, ovlsLen);

  //  Clean up.  Delete inputs, remove the sentinel, release memory, etc.

  delete [] ovls;
  delete [] bucketSizes;

  delete seq;

  if (deleteIntermediateLate) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Removing bucketized overlaps.\n");
    fprintf(stderr, "\n");

    writer->removeOverlapSlice();
  }

  delete writer;

  removeSentinel(ovlName, sliceNum);

  //  Success!

  fprintf(stderr, "Success!\n");

  return(0);
}

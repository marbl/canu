#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "meryl.H"
#include "merstreamfile.H"
#include "libmeryl.H"
#include "britime.H"


void
adjustHeap(u64bit *M, s64bit i, s64bit n) {
  u64bit   m = M[i];
  s64bit    j = (i << 1) + 1;  //  let j be the left child

  while (j < n) {
    if (j<n-1 && M[j] < M[j+1])
      j++;                   //  j is the larger child

    if (m >= M[j])           //  a position for M[i] has been found
      break;

    M[(j-1)/2] = M[j];       //  Move larger child up a level

    j = (j << 1) + 1;
  }

  M[(j-1)/2] = m;
}




void
runSegment(merylArgs *args, u64bit segment) {
  merStreamFileReader *R = 0L;
  merylStreamWriter   *W = 0L;
  speedCounter        *C = 0L;
  u32bit              *bucketSizes = 0L;
  u64bit              *bucketPointers = 0L;
  u64bit              *merData = 0L;

  //
  //  We can do all allocations up front:
  //    mer data storate (the buckets themselves, plus 64 for slop)
  //    bucket pointers (plus an extra bucket at the end and a little for slop)
  //    bucket size counting space, last because we toss it out quickly
  //
  if (args->beVerbose)
    fprintf(stderr, " Allocating "u64bitFMT"MB for mer storage ("u32bitFMT" bits wide).\n",
            (args->mersPerBatch * args->merDataWidth + 64) >> 23, args->merDataWidth);
  merData = new u64bit [ (args->mersPerBatch * args->merDataWidth + 64) >> 6 ];
  bzero(merData, sizeof(u64bit) * ((args->mersPerBatch * args->merDataWidth + 64) >> 6));

  if (args->beVerbose)
    fprintf(stderr, " Allocating "u64bitFMT"MB for bucket pointer table ("u32bitFMT" bits wide).\n",
            (args->numBuckets * args->bucketPointerWidth + 128) >> 23, args->bucketPointerWidth);
  bucketPointers = new u64bit [(args->numBuckets * args->bucketPointerWidth + 128) >> 6];
  bzero(bucketPointers, sizeof(u64bit) * ((args->numBuckets * args->bucketPointerWidth + 128) >> 6));

  if (args->beVerbose)
    fprintf(stderr, " Allocating "u64bitFMT"MB for counting the size of each bucket.\n", args->numBuckets >> 18);
  bucketSizes = new u32bit [ args->numBuckets ];
  for (u64bit i=args->numBuckets; i--; )
    bucketSizes[i] = u32bitZERO;

  if (args->beVerbose)
    fprintf(stderr, " Counting mers in buckets.\n");

  //  Position the mer stream at the start of this segments' mers.
  //  The last segment goes until the stream runs out of mers,
  //  everybody else does args->mersPerBatch mers.

  C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);
  R = new merStreamFileReader(args->outputFile);
  R->seekToMer(args->mersPerBatch * segment);
  R->setIterationLimit(args->mersPerBatch);

  if (args->doForward) {
    while (R->nextMer()) {

#if 0
      if (R->thePosition() > 157571920)
        fprintf(stderr, u64bitFMT"x"u64bitFMT" -- "u64bitHEX" -- %s\n", R->theSequenceNumber(), R->thePosition(), R->theFMer(), R->theFMerString());
      if ((R->theSequenceNumber() == 1) && (R->thePosition() == 157550003))
        exit(0);
#endif

      bucketSizes[ args->hash(R->theFMer()) ]++;
      C->tick();
    }
  }

  if (args->doReverse) {
    while (R->nextMer()) {
      bucketSizes[ args->hash(R->theRMer()) ]++;
      C->tick();
    }
  }

  if (args->doCanonical) {
    while (R->nextMer()) {
      if (R->theFMer() <= R->theRMer())
        bucketSizes[ args->hash(R->theFMer()) ]++;
      else
        bucketSizes[ args->hash(R->theRMer()) ]++;
      C->tick();
    }
  }

  delete C;
  delete R;

  if (args->beVerbose)
    fprintf(stderr, "\n");



  //
  //  Create the hash index using the counts.  The hash points
  //  to the end of the bucket; when we add a word, we move the
  //  hash bucket pointer down one.
  //
  //  When done, we can deallocate the counting table.
  //
  if (args->beVerbose)
    fprintf(stderr, " Initializing the bucket pointers.\n");

  {
    u64bit i=0;
    u64bit j=0;
    u64bit c=0;

    while (i < args->numBuckets) {
      c += bucketSizes[i++];
      setDecodedValue(bucketPointers, j, args->bucketPointerWidth, c);
      j += args->bucketPointerWidth;
    }

    //  Add the location of the end of the table.  This is not
    //  modified when adding words, but is used to determine
    //  the size of the last bucket.
    //
    setDecodedValue(bucketPointers, j, args->bucketPointerWidth, c);
  }


  //  All done with the counting table, get rid of it.
  //
  delete [] bucketSizes;






  if (args->beVerbose)
    fprintf(stderr, " Filling mers into list.\n");

  C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);
  R = new merStreamFileReader(args->outputFile);
  R->seekToMer(args->mersPerBatch * segment);
  R->setIterationLimit(args->mersPerBatch);

  while (R->nextMer()) {
    u64bit  m = R->theFMer();

    if ((args->doReverse) || (args->doCanonical && (m > R->theRMer())))
      m = R->theRMer();

    setDecodedValue(merData,
                    preDecrementDecodedValue(bucketPointers,
                                             args->hash(m) * args->bucketPointerWidth,
                                             args->bucketPointerWidth) * args->merDataWidth,
                    args->merDataWidth,
                    m);

    C->tick();
  }

  delete C;
  delete R;

  if (args->beVerbose)
    fprintf(stderr, "\n");





  if (args->beVerbose)
    fprintf(stderr, " 6) Writing output.\n");

  char *batchOutputFile = new char [strlen(args->outputFile) + 33];
  sprintf(batchOutputFile, "%s.batch"u64bitFMT, args->outputFile, segment);

  C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);
  W = new merylStreamWriter((args->segmentLimit == 1) ? args->outputFile : batchOutputFile,
                            args->merSize,
                            args->numBuckets_log2,
                            0,
                            0,
                            0);

  //  Sort each bucked into sortedList, then output the mers
  //
  u64bit  *sortedList    = 0L;
  u32bit   sortedListMax = 0;
  u32bit   sortedListLen = 0;

  for (u64bit bucket=0, bucketPos=0; bucket < args->numBuckets; bucket++) {
    u64bit st  = getDecodedValue(bucketPointers, bucketPos, args->bucketPointerWidth);
    bucketPos += args->bucketPointerWidth;
    u64bit ed  = getDecodedValue(bucketPointers, bucketPos, args->bucketPointerWidth);

#ifdef SANITY_CHECKS
    if (ed < st) {
      fprintf(stderr, "ERROR: Bucket "u64bitFMT" ends before it starts!\n", bucket);
      fprintf(stderr, "ERROR: start="u64bitFMT"\n", st);
      fprintf(stderr, "ERROR: end  ="u64bitFMT"\n", ed);
    }
#endif

    //  Nothing here?  Keep going.
    if (ed == st)
      continue;

    sortedListLen = (u32bit)(ed - st);

    //  Allocate more space, if we need to.
    //
    if (sortedListLen > sortedListMax) {
      delete [] sortedList;
      sortedList    = new u64bit [2 * sortedListLen + 1];
      sortedListMax = 2 * sortedListLen;
    }

    //  Unpack the mers into the sorting array
    //
    for (u64bit i=st, J=st*args->merDataWidth; i<ed; i++, J += args->merDataWidth)
      sortedList[i-st] = (bucket << args->merDataWidth) | getDecodedValue(merData, J, args->merDataWidth);

    //  Sort if there is more than one item
    //
    if (sortedListLen > 1) {
      for (s64bit t=(sortedListLen-2)/2; t>=0; t--)
        adjustHeap(sortedList, t, sortedListLen);

      for (s64bit t=sortedListLen-1; t>0; t--) {
        u64bit           tv = sortedList[t];
        sortedList[t]      = sortedList[0];
        sortedList[0]      = tv;

        adjustHeap(sortedList, 0, t);
      }
    }

    //  Dump the list of mers to the file.
    //
    for (u32bit t=0; t<sortedListLen; t++) {
      C->tick();
      W->addMer(sortedList[t]);
    }
  }

  delete C;
  delete W;

  delete [] batchOutputFile;

  delete [] merData;
  delete [] bucketPointers;

  if (args->beVerbose)
    fprintf(stderr, "\n");
}


















void
build(merylArgs *args) {
  bool  fatalError = false;

  if (args->inputFile == 0L)
    fprintf(stderr, "ERROR - no input file specified.\n"), fatalError = true;

  if (args->outputFile == 0L)
    fprintf(stderr, "ERROR - no output file specified.\n"), fatalError = true;

  if ((args->doForward == false) &&
      (args->doReverse == false) &&
      (args->doCanonical == false))
    fprintf(stderr, "ERROR - need to specify at least one of -f, -r, -C\n"), fatalError = true;

  if ((args->doForward && args->doReverse) ||
      (args->doForward && args->doCanonical) ||
      (args->doReverse && args->doCanonical))
    fprintf(stderr, "ERROR - only one of -f, -r and -C may be specified!\n"), fatalError = true;

  if (args->lowCount > args->highCount)
    fprintf(stderr, "ERROR - lowCount > highCount??\n"), fatalError = true;

  if (fatalError)
    exit(1);


  //  Everybody needs to dump the mers to a merStreamFile.  The only
  //  one that doesn't is "sequential, no memory limit", but
  //  optimizing for just that case makes the code too complex.
  //
  //  If the stream already exists, skip this step.
  //
  if (merStreamFileExists(args->outputFile)) {
    fprintf(stderr, "Using existing merStreamFile!\n");
    merStreamFileReader *R = new merStreamFileReader(args->outputFile);
    args->numMersActual = R->numberOfMers();
    delete R;
  } else {
    fprintf(stderr, "Building new merStreamFile!\n");
    merStreamFileBuilder   *B = new merStreamFileBuilder(args->merSize,
                                                         args->inputFile,
                                                         args->outputFile);

    args->numMersActual = B->build(args->beVerbose);

    delete B;
  }


  //  If there is a memory limit, ignore the total number of mers and
  //  pick a value that fits in memory.
  //
  //  Otherwise, if there is a segment limit, split the total number
  //  of mers into n pieces.  Remember, there cannot be both a
  //  memoryLimit and a segmentLimit.
  //
  //  Otherwise, we must be doing it all in one fell swoop.
  //
  if (args->memoryLimit) {
    args->mersPerBatch = estimateNumMersInMemorySize(args->merSize, args->memoryLimit, args->beVerbose);
    args->segmentLimit = (u64bit)ceil((double)args->numMersActual / (double)args->mersPerBatch);
  } else if (args->segmentLimit) {
    args->mersPerBatch = (u64bit)ceil((double)args->numMersActual / (double)args->segmentLimit);
  } else {
    args->mersPerBatch = args->numMersActual;
    args->segmentLimit = 1;
  }


  //  Choose the optimal number of buckets to reduce memory usage.
  //  Yes, this is already done in estimateNumMersInMemorySize() (but
  //  not saved) and we need to do it for the other cases anyway.
  //
  //  We use the number of mers per batch + 1 because we need to store
  //  the first position after the last mer.  That is, if there are
  //  two mers, we will store that the first mer is at position 0, the
  //  second mer is at position 1, and the end of the second mer is at
  //  position 2.
  //
  args->bucketPointerWidth = logBaseTwo_64(args->mersPerBatch + 1);
  args->numBuckets_log2    = optimalNumberOfBuckets(args->merSize, args->mersPerBatch);
  args->numBuckets         = (u64bitONE << args->numBuckets_log2);
  args->merDataWidth       = args->merSize * 2 - args->numBuckets_log2;
  args->bucketPointerMask  = u64bitMASK(args->numBuckets_log2);


#if 0
  if (args->beVerbose) {
    fprintf(stderr, "numMersActual      = "u64bitFMT"\n", args->numMersActual);
    fprintf(stderr, "numMersEstimated   = "u64bitFMT"\n", args->numMersEstimated);
    fprintf(stderr, "numBuckets         = "u64bitFMT"\n", args->numBuckets);
    fprintf(stderr, "numBuckets_log2    = "u32bitFMT"\n", args->numBuckets_log2);
    fprintf(stderr, "mersPerBatch       = "u64bitFMT"\n", args->mersPerBatch);
    fprintf(stderr, "merDataWidth       = "u32bitFMT"\n", args->merDataWidth);
    fprintf(stderr, "merDataMask        = "u64bitFMT"\n", args->merDataMask);
    fprintf(stderr, "bucketPointerWidth = "u32bitFMT"\n", args->bucketPointerWidth);
    fprintf(stderr, "bucketPointerMask  = "u64bitHEX"\n", args->bucketPointerMask);
    fprintf(stderr, "memoryLimit        = "u64bitFMT"\n", args->memoryLimit);
    fprintf(stderr, "segmentLimit       = "u64bitFMT"\n", args->segmentLimit);
    fprintf(stderr, "numThreads         = "u32bitFMT"\n", args->numThreads);
    fprintf(stderr, "batchNumber        = "u32bitFMT"\n", args->batchNumber);
  }
#endif


  //  Three choices:
  //
  //    threaded -- start threads, launch pieces in each thread.  This
  //    thread waits for completion and then merges the results.
  //
  //    batched -- write info file and exit.  Compute and merge is done
  //    on separate invocations.
  //
  //    segmented -- write info file, then do each piece sequentially.
  //    After all pieces finished, do a merge.
  //
  //    

  for (u64bit s=0; s<args->segmentLimit; s++) {

    //  Don't do it if the segment exists already
    //
    char *filename = new char [strlen(args->outputFile) + 17];
    sprintf(filename, "%s.batch"u64bitFMT".mcdat", args->outputFile, s);

    if (fileExists(filename)) {
      fprintf(stderr, "Found result for batch "u64bitFMT" in %s.\n", s, filename);
    } else {
      if ((args->beVerbose) && (args->segmentLimit > 1))
        fprintf(stderr, "Computing segment "u64bitFMT".\n", s);
      runSegment(args, s);
    }

    delete [] filename;
  }
  

  //  If there is more than one segment, merge them to get the output.
  //
  //  We do this by contructing a meryl command line and recursively
  //  (effectively) calling meryl.
  //
  //  The command line is
  //
  //  ./meryl -M add [-v] -s batch1 -s batch2 ... -s batchN -o outputFile
  //
  if (args->segmentLimit > 1) {
    int     argc = 0;
    char  **argv = new char* [7 + 2 * args->segmentLimit];

    argv[argc++] = "meryl-build-merge";
    argv[argc++] = "-M";
    argv[argc++] = "add";

    if (args->beVerbose)
      argv[argc++] = "-v";

    for (u32bit i=0; i<args->segmentLimit; i++) {
      argv[argc++] = "-s";
      argv[argc] = new char [strlen(args->outputFile) + 33];
      sprintf(argv[argc], "%s.batch"u32bitFMT, args->outputFile, i);
      argc++;
    }

    argv[argc++] = "-o";
    argv[argc++] = args->outputFile;

    merylArgs *addArgs = new merylArgs(argc, argv);
    multipleOperations(addArgs);

    //  Cleanup the memory leak.

    //  Remove temporary files
    //
#if 0
    char *filename = new char [strlen(args->outputFile) + 17];
    for (u32bit i=0; i<args->segmentLimit; i++) {
      sprintf(filename, "%s.batch"u32bitFMT, args->outputFile, i);
      unlink(filename);
    }
    delete [] filename;
#endif
  }
}

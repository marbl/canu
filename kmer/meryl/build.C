#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libbri.H"
#include "meryl.H"
#include "outputMer.H"
#include "britime.H"

#include "merstreamfile.H"


void
runSegment(merylArgs *args, u64bit segment);


void
build(merylArgs *args) {
  u64bit         *chck = 0L;
  u64bit         *hash = 0L;

  bool  fatalError = false;

  if (args->inputFile == 0L)
    fprintf(stderr, "ERROR - no input file specified.\n"), fatalError = true;

  if (args->outputFile == 0L)
    fprintf(stderr, "ERROR - no output file specified.\n"), fatalError = true;

  //  these should never happen, unles main() is broken.
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

#if 0
  //  30jan03
  //  It would appear that we need at least 2 bits in the table.
  //
  if ((args->merSize << 1) < args->tblSize + 2) {
    fprintf(stderr, "WARNING:  "u32bitFMT" bit table is too big for "u32bitFMT" bit mers.\n", args->tblSize, (args->merSize << 1));
    fprintf(stderr, "WARNING:  adjusting table size (-t) to "u32bitFMT"\n", (args->merSize << 1) - 2);
    args->tblSize = (args->merSize << 1) - 2;
  }
#endif


  //  Everybody needs to dump the mers to a merStreamFile.  The only
  //  one that doesn't is "sequential, no memory limit", but
  //  optimizing for just that case makes the code too complex.
  //
  merStreamFileBuilder   *B = new merStreamFileBuilder(args->merSize,
                                                       args->inputFile,
                                                       args->outputFile);

  //  Build the merstreamfile, saving the number of mers in that file.
  //
  args->numMersActual = B->build(args->beVerbose);

  delete B;

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
  args->numBuckets         = u64bitONE << args->numBuckets_log2;
  args->merDataWidth       = args->merSize * 2 - args->numBuckets_log2;

  //  Compute some useful values for use elsewhere
  //
  args->merDataWidth      = args->merSize * 2 - args->numBuckets_log2; 
  args->bucketPointerMask = u64bitMASK(args->numBuckets_log2);


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
    runSegment(args, s);
  }
  

  //  Merge the segments together


  delete [] chck;
  delete [] hash;
}








void
runSegment(merylArgs *args, u64bit segment) {

  merStreamFileReader *R = 0L;
  speedCounter        *C = 0L;
  u32bit              *bucketSizes = 0L;
  u64bit              *bucketPointers = 0L;
  u64bit              *merData = 0L;


  //
  //  We can do all allocations up front.
  //



  //  The buckets themselves
  //
  if (args->beVerbose)
    fprintf(stderr, " Allocating "u64bitFMT"MB for mer storage ("u32bitFMT" bits wide).\n",
            (args->mersPerBatch * args->merDataWidth + 64) >> 23, args->merDataWidth);
  merData = new u64bit [ (args->mersPerBatch * args->merDataWidth + 64) >> 6 ];


  //  Bucket pointers
  //
  if (args->beVerbose)
    fprintf(stderr, " Allocating "u64bitFMT"MB for bucket pointer table ("u32bitFMT" bits wide).\n",
            (args->numBuckets * args->bucketPointerWidth + 128) >> 23, args->bucketPointerWidth);
  bucketPointers = new u64bit [(args->numBuckets * args->bucketPointerWidth + 128) >> 6];


  //  Bucket size counting space
  //
  if (args->beVerbose)
    fprintf(stderr, " Allocating "u64bitFMT"MB for counting the size of each bucket.\n", args->numBuckets >> 18);
  bucketSizes = new u32bit [ args->numBuckets ];



  if (args->beVerbose)
    fprintf(stderr, " 1) Counting mers in buckets.\n");

  //  XXX:  Need to pick out a subset of the merstreamfile to use

  for (u64bit i=args->numBuckets; i--; )
    bucketSizes[i] = 0;

  R = new merStreamFileReader(args->outputFile);
  C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);

  if (args->doForward) {
    while (R->nextMer()) {
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
    fprintf(stderr, "\n    Found "u64bitFMT" mers.\n", args->numMersActual);





  //
  //  Create the hash index using the counts.  The hash points
  //  to the end of the bucket; when we add a word, we move the
  //  hash bucket pointer down one.
  //
  //  When done, we can deallocate the counting table.
  //
  if (args->beVerbose)
    fprintf(stderr, " 3) Initializing the bucket pointers.\n");

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


  //  All done with the counting table, get rid of it.
  //
  delete [] bucketSizes;






  if (args->beVerbose)
    fprintf(stderr, " 5) Filling mers into list.\n");

  R = new merStreamFileReader(args->outputFile);
  C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);

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

  if (args->beVerbose)
    fprintf(stderr, "\n");

  delete C;
  delete R;
}

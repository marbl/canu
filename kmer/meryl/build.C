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
  args->numBuckets         = (u64bitONE << args->numBuckets_log2);
  args->merDataWidth       = args->merSize * 2 - args->numBuckets_log2;
  args->bucketPointerMask  = u64bitMASK(args->numBuckets_log2);


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

  //  Bucket pointers - plus 128 for the extra bucket at the end, and a little slop
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
      fprintf(stderr, u64bitHEX" -- %s\n", R->theFMer(), R->theFMerString());
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





















  if (args->beVerbose)
    fprintf(stderr, " 6) Writing output.\n");


  //
  //  XXX:  The writer need to know how many buckets there are!
  //


  merylStreamWriter *W = new merylStreamWriter(args->outputFile,
                                               args->merSize,
                                               args->numBuckets_log2,
                                               0,
                                               0,
                                               0);


  C = new speedCounter(" %7.2f Mbuckets -- %5.2f Mbuckets/second\r", 1000000.0, 0x1fffff, args->beVerbose);


  //
  //  Sort each bucked into sortedList, then output
  //


  u64bit  *sortedList    = 0L;
  u32bit   sortedListMax = 0;
  u32bit   sortedListLen = 0;


  fprintf(stderr, "sorting/writing "u64bitFMT" buckets\n", args->numBuckets);


  //  For each bucket, sort it.  The output is done in the sort.
  //
  for (u64bit bucket=0, bucketPos=0; bucket < args->numBuckets; bucket++) {
    C->tick();

    u64bit st  = getDecodedValue(bucketPointers, bucketPos, args->bucketPointerWidth);
    bucketPos += args->bucketPointerWidth;
    u64bit ed  = getDecodedValue(bucketPointers, bucketPos, args->bucketPointerWidth);

    if (ed < st) {
      fprintf(stderr, "ERROR: Bucket "u64bitFMT" ends before it starts!\n", bucket);
      fprintf(stderr, "ERROR: start="u64bitFMT"\n", st);
      fprintf(stderr, "ERROR: end  ="u64bitFMT"\n", ed);
    }

    fprintf(stderr, "sorting/writing bucket "u64bitFMT" of size "u64bitFMT"\n", bucket, ed-st);

    //  Nothing here?  Keep going.
    if (ed == st)
      continue;

    sortedListLen = (u32bit)(ed - st);

    //  Allocate more space, if we need to.
    //
    if (sortedListLen > sortedListMax) {
      delete [] sortedList;
      sortedList    = new u64bit [sortedListLen + 1];
      sortedListMax = sortedListLen;
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
    for (u32bit t=0; t<sortedListLen; t++)
      W->addMer(sortedList[t]);
  }

  delete C;
  delete W;

  delete [] merData;
  delete [] bucketPointers;

  if (args->beVerbose)
    fprintf(stderr, "\n");
}

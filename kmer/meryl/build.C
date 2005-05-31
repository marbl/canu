#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bio++.H"
#include "meryl.H"
#include "libmeryl.H"

#ifdef ENABLE_THREADS
void runThreaded(merylArgs *args);
#endif

//  Define this to remove temporary files.  You probably want it
//  defined.
#define REMOVE_TEMPORARY_FILES


//  Threaded build times on G5/2GHz, including merStreamFile building,
//  of 167256385 mers.
//
//  1 314.100u  5.900s 5:21.84 99.4%  0+0k 0+ 31io 0pf+0w (2.45, 1.06, 2.67)
//  2 376.620u  9.970s 6:28.84 99.4%  0+0k 0+ 55io 0pf+0w (2.67, 1.07, 2.66, 2.38)
//  4 363.760u 10.670s 6:18.72 98.8%  0+0k 0+ 57io 0pf+0w (3.15, 1.07, 2.48, 2.15)
//  7 347.130u 12.100s 6:31.48 91.7%  0+0k 0+ 95io 0pf+0w (3.5x, 1.05, 2.30, 1.94)
//  8 333.330u 10.410s 5:46.75 99.1%  0+0k 0+ 75io 0pf+0w (3.78, 1.18, 2.65)
// 32 335.240u 12.670s 5:50.42 99.2%  0+0k 0+104io 0pf+0w (5.20, 1.33, 2.41, 1.61)
// 64 338.400u 14.170s 5:54.96 99.3%  0+0k 0+ 96io 0pf+0w (6.6x, 1.50, 2.42, 1.28)
//
//  Notice that writing the output (includes the sort) is constant
//  speed.  Buckets are small enough to fit in cache.
//
//  The slowdown in merge is probably implementation related.
//
//  Times might get N seconds better (N segments) if we stop destroying
//  the merStreamFile for each segment.  Easy enough if only one thread,
//  trouble if more than one.
//
//  md5's of the dumps of all those are the same.  The file sizes are all
//  different (except 7 and 8):
//
//  330677936  t0.mcdat  6116056 t0.mcidx
//  348964824  t2.mcdat  3669912 t2.mcidx
//  367251712  t4.mcdat  2202664 t4.mcidx
//  385538600  t7.mcdat  1281176 t7.mcidx
//  385538600  t8.mcdat  1281176 t8.mcidx
//  422112376  t32.mcdat  418704 t32.mcidx
//  440399264  t64.mcdat  235592 t64.mcidx


//  You probably want this to be the same as KMER_WORDS, but in rare
//  cases, it can be less.
//
#define SORTED_LIST_WIDTH  KMER_WORDS

//  to make the sorted list be wider, we also need to store wide
//  things in the bitpackedarray buckets.  probably easy (do multiple
//  adds of data, each at most 64 bits) but not braindead.

#if SORTED_LIST_WIDTH == 1

typedef u64bit  sortedList_t;

#else

class sortedList_t {
public:
  u64bit    _w[SORTED_LIST_WIDTH];

  bool operator<(sortedList_t &that) {
    for (u32bit i=SORTED_LIST_WIDTH; i--; ) {
      if (_w[i] < that._w[i])  return(true);
      if (_w[i] > that._w[i])  return(false);
    }
    return(false);
  };

  bool operator>=(sortedList_t &that) {
    for (u32bit i=SORTED_LIST_WIDTH; i--; ) {
      if (_w[i] > that._w[i])  return(true);
      if (_w[i] < that._w[i])  return(false);
    }
    return(true);
  };

  sortedList_t &operator=(sortedList_t &that) {
    for (u32bit i=SORTED_LIST_WIDTH; i--; )
      _w[i] = that._w[i];
    return(*this);
  };
};

#endif




void
adjustHeap(sortedList_t *M, s64bit i, s64bit n) {
  sortedList_t   m = M[i];
  s64bit         j = (i << 1) + 1;  //  let j be the left child

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
prepareBatch(merylArgs *args) {
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

  if (args->segmentLimit && args->memoryLimit)
    fprintf(stderr, "ERROR: Only one of -memory and -segments can be specified.\n"), fatalError=true;

  if (fatalError)
    exit(1);

#ifdef ENABLE_THREADS
  if (args->numThreads > 0) {
    //  If we were given no segment or memory limit, but threads, we
    //  really want to create n segments.
    //
    if ((args->segmentLimit == 0) && (args->memoryLimit == 0)) {
      args->segmentLimit = args->numThreads;
    }

    //  If we are given a memory limit and threads, we want to use that much memory
    //  total, not per thread.
    //
    if ((args->memoryLimit > 0) && (args->numThreads > 0)) {
      args->segmentLimit = 0;
      args->memoryLimit /= args->numThreads;
    }
  }
#endif

  //  Everybody needs to dump the mers to a merStreamFile.  The only
  //  one that doesn't is "sequential, no memory limit", but
  //  optimizing for just that case makes the code too complex.
  //
  //  If the stream already exists, skip this step.
  //
  if (merStreamFileExists(args->outputFile)) {
    fprintf(stderr, "Using existing merStreamFile!\n");
    merStreamFileReader *R = new merStreamFileReader(args->outputFile, args->merSize);
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
    fprintf(stderr, "Have a memory limit: mersPerBatch="u64bitFMT" segmentLimit="u64bitFMT"\n", args->mersPerBatch, args->segmentLimit);
  } else if (args->segmentLimit) {
    args->mersPerBatch = (u64bit)ceil((double)args->numMersActual / (double)args->segmentLimit);
    fprintf(stderr, "Have a segment limit: mersPerBatch="u64bitFMT" segmentLimit="u64bitFMT"\n", args->mersPerBatch, args->segmentLimit);
  } else {
    args->mersPerBatch = args->numMersActual;
    args->segmentLimit = 1;
    fprintf(stderr, "Have NO LIMITS!: mersPerBatch="u64bitFMT" segmentLimit="u64bitFMT"\n", args->mersPerBatch, args->segmentLimit);
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
  args->bucketPointerWidth = logBaseTwo64(args->mersPerBatch + 1);
  args->numBuckets_log2    = optimalNumberOfBuckets(args->merSize, args->mersPerBatch);
  args->numBuckets         = (u64bitONE << args->numBuckets_log2);
  args->merDataWidth       = args->merSize * 2 - args->numBuckets_log2;
  //args->bucketPointerMask  = u64bitMASK(args->numBuckets_log2);


  if (args->merDataWidth > SORTED_LIST_WIDTH * 64) {
    fprintf(stderr, "  numMersActual      = "u64bitFMT"\n", args->numMersActual);
    fprintf(stderr, "  mersPerBatch       = "u64bitFMT"\n", args->mersPerBatch);
    fprintf(stderr, "  numBuckets         = "u64bitFMT" ("u32bitFMT" bits)\n", args->numBuckets, args->numBuckets_log2);
    fprintf(stderr, "  bucketPointerWidth = "u32bitFMT"\n", args->bucketPointerWidth);
    fprintf(stderr, "  merDataWidth       = "u32bitFMT"\n", args->merDataWidth);
    fprintf(stderr, "Sorry!  Your data width can't fit into SORTED_LIST_WIDTH=%d 64-bit words.\n",
            SORTED_LIST_WIDTH);
    exit(1);
  }

  if (args->beVerbose) {
    if (args->memoryLimit)
      fprintf(stderr, "Computing "u64bitFMT" segments using "u64bitFMT"MB memory each.\n",
              args->segmentLimit, args->memoryLimit);
    else
      fprintf(stderr, "Computing "u64bitFMT" segments using AS MUCH MEMORY AS NEEDED.\n",
              args->segmentLimit);
    fprintf(stderr, "  numMersActual      = "u64bitFMT"\n", args->numMersActual);
    fprintf(stderr, "  mersPerBatch       = "u64bitFMT"\n", args->mersPerBatch);
    fprintf(stderr, "  numBuckets         = "u64bitFMT" ("u32bitFMT" bits)\n", args->numBuckets, args->numBuckets_log2);
    fprintf(stderr, "  bucketPointerWidth = "u32bitFMT"\n", args->bucketPointerWidth);
    fprintf(stderr, "  merDataWidth       = "u32bitFMT"\n", args->merDataWidth);
  }
}




void
runSegment(merylArgs *args, u64bit segment) {
  merStreamFileReader *R = 0L;
  merylStreamWriter   *W = 0L;
  speedCounter        *C = 0L;
  u32bit              *bucketSizes = 0L;
  u64bit              *bucketPointers = 0L;
  u64bit              *merDataArray[SORTED_LIST_WIDTH] = { 0L };

  //  If this segment exists already, skip it.
  //
  //  XXX:  This should be a command line option.
  //  XXX:  This should check that the files are complete meryl files.
  //
  char *filename = new char [strlen(args->outputFile) + 17];
  sprintf(filename, "%s.batch"u64bitFMT".mcdat", args->outputFile, segment);

  if (fileExists(filename)) {
    fprintf(stderr, "Found result for batch "u64bitFMT" in %s.\n", segment, filename);
    delete [] filename;
    return;
  }

  if ((args->beVerbose) && (args->segmentLimit > 1))
    fprintf(stderr, "Computing segment "u64bitFMT" of "u64bitFMT".\n", segment+1, args->segmentLimit);

  delete [] filename;



  //
  //  We can do all allocations up front:
  //    mer data storate (the buckets themselves, plus 64 for slop)
  //    bucket pointers (plus an extra bucket at the end and a little for slop)
  //    bucket size counting space, last because we toss it out quickly
  //
  if (args->beVerbose)
    fprintf(stderr, " Allocating "u64bitFMT"MB for mer storage ("u32bitFMT" bits wide).\n",
            (args->mersPerBatch * args->merDataWidth + 64) >> 23, args->merDataWidth);

#if SORTED_LIST_WIDTH == 1
  merDataArray[0] = new u64bit [ (args->mersPerBatch * args->merDataWidth + 64) >> 6 ];
#else
  for (u64bit mword=0, width=args->merDataWidth; width > 0; ) {
    if (width >= 64) {
      merDataArray[mword] = new u64bit [ args->mersPerBatch + 1 ];
      width -= 64;
      mword++;
    } else {
      merDataArray[mword] = new u64bit [ (args->mersPerBatch * width + 64) >> 6 ];
      width  = 0;
    }
  }
#endif
  //bzero(merData, sizeof(u64bit) * ((args->mersPerBatch * args->merDataWidth + 64) >> 6));



  if (args->beVerbose)
    fprintf(stderr, " Allocating "u64bitFMT"MB for bucket pointer table ("u32bitFMT" bits wide).\n",
            (args->numBuckets * args->bucketPointerWidth + 128) >> 23, args->bucketPointerWidth);
  bucketPointers = new u64bit [(args->numBuckets * args->bucketPointerWidth + 128) >> 6];
  //bzero(bucketPointers, sizeof(u64bit) * ((args->numBuckets * args->bucketPointerWidth + 128) >> 6));

  if (args->beVerbose)
    fprintf(stderr, " Allocating "u64bitFMT"MB for counting the size of each bucket.\n", args->numBuckets >> 18);
  bucketSizes = new u32bit [ args->numBuckets ];
  for (u64bit i=args->numBuckets; i--; )
    bucketSizes[i] = u32bitZERO;

  //  Position the mer stream at the start of this segments' mers.
  //  The last segment goes until the stream runs out of mers,
  //  everybody else does args->mersPerBatch mers.

  C = new speedCounter(" Counting mers in buckets: %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);
  R = new merStreamFileReader(args->outputFile, args->merSize);
  R->seekToMer(args->mersPerBatch * segment);
  R->setIterationLimit(args->mersPerBatch);

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


  //  Create the hash index using the counts.  The hash points
  //  to the end of the bucket; when we add a word, we move the
  //  hash bucket pointer down one.
  //
  //  When done, we can deallocate the counting table.
  //
  if (args->beVerbose)
    fprintf(stderr, " Creating bucket pointers.\n");

  {
    u64bit mi=0;
    u64bit mj=0;
    u64bit mc=0;

    while (mi < args->numBuckets) {
      mc += bucketSizes[mi++];
      setDecodedValue(bucketPointers, mj, args->bucketPointerWidth, mc);
      mj += args->bucketPointerWidth;
    }

    //  Add the location of the end of the table.  This is not
    //  modified when adding words, but is used to determine
    //  the size of the last bucket.
    //
    setDecodedValue(bucketPointers, mj, args->bucketPointerWidth, mc);
  }


  //  All done with the counting table, get rid of it.
  //
  if (args->beVerbose)
    fprintf(stderr, " Releasing "u64bitFMT"MB from counting the size of each bucket.\n", args->numBuckets >> 18);
  delete [] bucketSizes;


  C = new speedCounter(" Filling mers into list:   %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);
  R = new merStreamFileReader(args->outputFile, args->merSize);
  R->seekToMer(args->mersPerBatch * segment);
  R->setIterationLimit(args->mersPerBatch);

  while (R->nextMer()) {

    kMer const &m =  ((args->doReverse) || (args->doCanonical && (R->theFMer() > R->theRMer()))) ?
      R->theRMer()
      :
      R->theFMer();

    //  additionally, we can use the real bitPackedArray here, since everything
    //  is fixed size.

    u64bit  element = preDecrementDecodedValue(bucketPointers,
                                               args->hash(m) * args->bucketPointerWidth,
                                               args->bucketPointerWidth);

#if SORTED_LIST_WIDTH == 1
    setDecodedValue(merDataArray[0],
                    element * args->merDataWidth,
                    args->merDataWidth,
                    m.endOfMer(args->merDataWidth));
#else
    for (u64bit mword=0, width=args->merDataWidth; width>0; ) {
      if (width >= 64) {
        merDataArray[mword][element] = m.getWord(mword);
        width -= 64;
        mword++;
      } else {
        setDecodedValue(merDataArray[mword],
                        element * width,
                        width,
                        m.getWord(mword) & u64bitMASK(width));
        width = 0;
      }
    }
#endif

    C->tick();
  }

  delete C;
  delete R;


  char *batchOutputFile = new char [strlen(args->outputFile) + 33];
  sprintf(batchOutputFile, "%s.batch"u64bitFMT, args->outputFile, segment);

  C = new speedCounter(" Writing output:           %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);
  W = new merylStreamWriter((args->segmentLimit == 1) ? args->outputFile : batchOutputFile,
                            args->merSize,
                            args->numBuckets_log2);

  //  Sort each bucked into sortedList, then output the mers
  //
  sortedList_t  *sortedList    = 0L;
  u32bit         sortedListMax = 0;
  u32bit         sortedListLen = 0;

  for (u64bit bucket=0, bucketPos=0; bucket < args->numBuckets; bucket++) {
    u64bit st  = getDecodedValue(bucketPointers, bucketPos, args->bucketPointerWidth);
    bucketPos += args->bucketPointerWidth;
    u64bit ed  = getDecodedValue(bucketPointers, bucketPos, args->bucketPointerWidth);

    //fprintf(stderr, "bucket="u64bitFMT" st="u64bitFMT" ed="u64bitFMT"\n", bucket, st, ed);

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
      sortedList    = new sortedList_t [2 * sortedListLen + 1];
      sortedListMax = 2 * sortedListLen;
    }

    //  XXX:  We shouldn't need this, but alpha does.
    //
    //  If we don't use all the words in SORTED_LIST_WIDTH, we
    //  definitely leave the higher ones unset, and that will case the
    //  sort to be random.
    //
    bzero(sortedList, sizeof(sortedList_t) * sortedListLen);

    //  Unpack the mers into the sorting array
    //
#if SORTED_LIST_WIDTH == 1
    for (u64bit i=st, J=st*args->merDataWidth; i<ed; i++, J += args->merDataWidth)
      sortedList[i-st] = getDecodedValue(merDataArray[0],
                                         J,
                                         args->merDataWidth);
#else
    for (u64bit i=st; i<ed; i++) {
      for (u64bit mword=0, width=args->merDataWidth; width>0; ) {
        if (width >= 64) {
          sortedList[i-st]._w[mword] = merDataArray[mword][i];
          width -= 64;
          mword++;
        } else {
          sortedList[i-st]._w[mword] = getDecodedValue(merDataArray[mword],
                                                       i * width,
                                                       width);
          width = 0;
        }
      }
    }
#endif

    //  Sort if there is more than one item
    //
    if (sortedListLen > 1) {
      for (s64bit t=(sortedListLen-2)/2; t>=0; t--)
        adjustHeap(sortedList, t, sortedListLen);

      for (s64bit t=sortedListLen-1; t>0; t--) {
        sortedList_t    tv = sortedList[t];
        sortedList[t]      = sortedList[0];
        sortedList[0]      = tv;

        adjustHeap(sortedList, 0, t);
      }
    }

    //  Dump the list of mers to the file.
    //
    kMer  mer(args->merSize);

    for (u32bit t=0; t<sortedListLen; t++) {
      C->tick();

#if 0
      fprintf(stderr, u32bitFMTW(2)": "u64bitHEX" "u64bitHEX" "u64bitHEX"\n",
              t, sortedList[t]._w[0], sortedList[t]._w[1], sortedList[t]._w[2]);
#endif

      //  Build the complete mer
      //
#if SORTED_LIST_WIDTH == 1
      //mer.setBits(0, args->merDataWidth, sortedList[t]);
      mer.setWord(0, sortedList[t]);
      mer.setBits(args->merDataWidth, args->numBuckets_log2, bucket);
#else
      for (u64bit mword=0; mword < SORTED_LIST_WIDTH; mword++)
        mer.setWord(mword, sortedList[t]._w[mword]);
      mer.setBits(args->merDataWidth, args->numBuckets_log2, bucket);
#endif

      //  Add it
      W->addMer(mer);
    }
  }

  delete [] sortedList;

  delete C;
  delete W;

  delete [] batchOutputFile;

  for (u32bit x=0; x<SORTED_LIST_WIDTH; x++)
    delete [] merDataArray[x];

  delete [] bucketPointers;

  if (args->beVerbose)
    fprintf(stderr, "Segment "u64bitFMT" finished.\n", segment);
}


















void
build(merylArgs *args) {

  if (!args->countBatch && !args->mergeBatch)
    prepareBatch(args);


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

  bool  doMerge = false;

#ifdef ENABLE_THREADS
  if (args->numThreads > 1) {

    //  Run, using threads.  There is a lot of baloney needed, so it's
    //  all in a separate function.
    //
    runThreaded(args);
    doMerge = true;
  } else
#endif
  if (args->configBatch) {

    //  Write out our configuration and exit if we are -configbatch
    //
    args->writeConfig();

    fprintf(stdout, "Batch prepared.  Please run:\n");
    for (u64bit s=0; s<args->segmentLimit; s++)
      fprintf(stdout, "%s -countbatch "u64bitFMT" -o %s\n", args->execName, s, args->outputFile);
    fprintf(stdout, "%s -mergebatch -o %s\n", args->execName, args->outputFile);
  } else   if (args->countBatch) {

    //  Read back the configuration, run the segment and exit if we
    //  are -countbatch
    //
    merylArgs *savedArgs = new merylArgs(args->outputFile);
    savedArgs->beVerbose = args->beVerbose;
    runSegment(savedArgs, args->batchNumber);
    delete savedArgs;
  } else if (args->mergeBatch) {

    //  Check that all the files exist if we are -mergebatch and
    //  continue with execution
    //
    //  MEMORY LEAK!  We should delete this at the end of the
    //  function, but it's a pain, and who cares?
    //
    merylArgs *savedArgs = new merylArgs(args->outputFile);
    savedArgs->beVerbose = args->beVerbose;

    args = savedArgs;

    doMerge = true;
  } else {

    //  No special options given, do all the work here and now
    //
    for (u64bit s=0; s<args->segmentLimit; s++)
      runSegment(args, s);

    doMerge = true;
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
  if ((doMerge) && (args->segmentLimit > 1)) {

    if (args->beVerbose)
      fprintf(stderr, "Merge results.\n");

    int     argc = 0;
    char  **argv = new char* [7 + 2 * args->segmentLimit];
    bool   *arga = new bool  [7 + 2 * args->segmentLimit];

    arga[argc] = false;  argv[argc++] = "meryl-build-merge";
    arga[argc] = false;  argv[argc++] = "-M";
    arga[argc] = false;  argv[argc++] = "add";

    if (args->beVerbose) {
      arga[argc] = false;
      argv[argc++] = "-v";
    }

    for (u32bit i=0; i<args->segmentLimit; i++) {
      arga[argc] = false;
      argv[argc++] = "-s";
      arga[argc] = true;
      argv[argc] = new char [strlen(args->outputFile) + 33];
      sprintf(argv[argc], "%s.batch"u32bitFMT, args->outputFile, i);
      argc++;
    }

    arga[argc] = false;  argv[argc++] = "-o";
    arga[argc] = false;  argv[argc++] = args->outputFile;

    merylArgs *addArgs = new merylArgs(argc, argv);
    multipleOperations(addArgs);

    //  Cleanup the memory leak.
    //
    delete addArgs;
    for (int i=0; i<argc; i++)
      if (arga[i])
        delete [] argv[i];
    delete [] argv;
    delete [] arga;

    //  Remove temporary files
    //
    char *filename = new char [strlen(args->outputFile) + 17];

#ifdef REMOVE_TEMPORARY_FILES
    for (u32bit i=0; i<args->segmentLimit; i++) {
      sprintf(filename, "%s.batch"u32bitFMT".mcidx", args->outputFile, i);
      unlink(filename);
      sprintf(filename, "%s.batch"u32bitFMT".mcdat", args->outputFile, i);
      unlink(filename);
    }
#endif

    delete [] filename;
  }

  //  If we just merged, delete the merstream file
  //
#ifdef REMOVE_TEMPORARY_FILES
  if (doMerge) {
    char *filename = new char [strlen(args->outputFile) + 17];

    sprintf(filename, "%s.merStream", args->outputFile);
    unlink(filename);

    delete [] filename;
  }
#else
    fprintf(stderr, "NOT REMOVING TEMPORARY FILES!\n");
#endif
}

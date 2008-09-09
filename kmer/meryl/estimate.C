#include <stdio.h>
#include <stdlib.h>

#include "bio++.H"
#include "seqStream.H"
#include "merStream.H"
#include "libmeryl.H"
#include "meryl.H"

//  Takes a memory limit in MB, returns the number of mers that we can
//  fit in that memory size, assuming optimalNumberOfBuckets() below
//  uses the same algorithm.
//
u64bit
estimateNumMersInMemorySize(u32bit merSize,
                            u32bit mem,
                            bool   positionsEnabled,
                            bool   beVerbose) {
  u64bit maxN    = 0;
  u64bit bestT   = 0;


  //  For each possible number of buckets, try all poissible pointer
  //  widths.  First we compute the number of mers that fit in a
  //  bucket pointer table of size 2^t storing N bits in the mer data
  //  table, then we check that the number of mers in the mer data
  //  table agrees with the width of the pointer table.


  //  This is the memory size we are trying to fill, in bits.
  //
  u64bit memLimt = ((u64bit)mem) << 23;

  //  Positions consume space too, but only if enabled.
  //
  u64bit posPerMer = 0;
  if (positionsEnabled)
    posPerMer = 32;

  //  Limit the number of entries in the bucket pointer table to
  //  50 bits -- thus, the prefix of each mer is at most 25.
  //
  u32bit  tMax = 2*merSize - 2;
  if (tMax > 50)
    tMax = 50;

  for (u64bit t=2; t < tMax; t++) {

    //  We need to try all possibilities of N, the width of the
    //  bucket pointer table === log2(numMers).
    //
    //  Increased to 40 bits, so we're valid up to 1 trillion mers.
    //
    for (u64bit N=1; N<40; N++) {
      u64bit Nmin = u64bitONE << (N - 1);
      u64bit Nmax = u64bitONE << (N);

      //  The size in bits of the bucket pointer table.
      //
      u64bit bucketsize = (u64bitONE << t) * N;

      //  If our bucket pointer table size hasn't already blown our
      //  memory limit, compute the number of mers that we can stuff
      //  into the list.
      //
      if (memLimt > bucketsize) {

        //  The number of mers we can then fit into the mer data table
        //  is easy to compute.
        //
        //  Even though we allocate merDataArray, bucketPointers,
        //  bucketSizes, we don't use merDataArray until after we
        //  release bucketSizes, and so we only estimate the maximum
        //  in core (not allocated) size.
        //
        u64bit n = (memLimt - bucketsize) / (2*merSize - t + posPerMer);

        //  We can stop now if our computed number of mers is outside the range that
        //  the bucket pointer table can address.
        //
        if ((Nmin <= n) && (n <= Nmax)) {

          //fprintf(stderr, "prefixSize="u64bitFMTW(2)" numMers="u64bitFMTW(9)" memory=%.3fMB\n",
          //        t, n,
          //        (((u64bitONE << t) * logBaseTwo64(n) + n * (2*merSize - t + posPerMer)) >> 3) / 1048576.0);

          //  Remember the settings with the highest number of mers, if
          //  more than zero mers.
          //  
          if ((n > 0) &&
              (maxN < n)) {
            maxN  = n;
            bestT = t;
          }

        }
      }
    }
  }

  if (beVerbose)
    fprintf(stdout, "Can fit "u64bitFMT" mers into table with prefix of "u64bitFMT" bits, using %8.3fMB (%8.3fMB for positions)\n",
            maxN,
            bestT,
            (((u64bitONE << bestT) * logBaseTwo64(maxN) + maxN * (2*merSize - bestT + posPerMer)) >> 3) / 1048576.0,
            ((maxN * posPerMer) >> 3) / 1048576.0);

  return(maxN);
}



u32bit
optimalNumberOfBuckets(u32bit merSize,
                       u64bit numMers,
                       bool   positionsEnabled) {
  u64bit opth   = ~u64bitZERO;
  u64bit opts   = ~u64bitZERO;
  u64bit h      = 0;
  u64bit s      = 0;
  u64bit hwidth = logBaseTwo64(numMers);

  //  Positions consume space too, but only if enabled.  Probably
  //  doesn't matter here.
  //
  u64bit posPerMer = 0;
  if (positionsEnabled)
    posPerMer = 32;

  //  Find the table size (in bits, h) that minimizes memory usage
  //  for the given merSize and numMers
  //
  //  We have two tables:
  //    the bucket pointers num buckets * pointer width   == 2 << h * hwidth
  //    the mer data:       num mers * (mersize - hwidth)
  //
  u64bit hmax = 64 - logBaseTwo64(hwidth + numMers * (2 * merSize - h));
  for (h=2; h<=hmax && h<2*merSize; h++) {
    s = (u64bitONE << h) * hwidth + numMers * (2 * merSize - h + posPerMer);

    //fprintf(stderr, "optimalNumberOfBuckets()-- h="u64bitFMT" s="u64bitFMT"\n", h, s);

    if (s < opts) {
      opth = h;
      opts = s;
    }
  }

  return((u32bit)opth);
}



void
estimate(merylArgs *args) {

  if (args->inputFile) {
    merStream          M(new kMerBuilder(args->merSize, args->merComp),
                         new seqStream(args->inputFile),
                         true, true);
    speedCounter       C(" %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);

    if (args->beVerbose)
      fprintf(stderr, "Counting mers in '%s'\n", args->inputFile);

    args->numMersEstimated = 0;

    while (M.nextMer()) {
      C.tick();
      args->numMersEstimated++;
    }

    C.finish();
  }

  u32bit opth = optimalNumberOfBuckets(args->merSize, args->numMersEstimated, args->positionsEnabled);
  u64bit memu = ((u64bitONE << opth)    * logBaseTwo64(args->numMersEstimated+1) +
                 args->numMersEstimated * (2 * args->merSize - opth));

  fprintf(stderr, u64bitFMT" "u32bitFMT"-mers can be computed using "u64bitFMT"MB memory.\n",
          args->numMersEstimated, args->merSize, memu >> 23);
}

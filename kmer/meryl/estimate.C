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
uint64
estimateNumMersInMemorySize(uint32 merSize,
                            uint32 mem,
                            bool   positionsEnabled,
                            bool   beVerbose) {
  uint64 maxN    = 0;
  uint64 bestT   = 0;


  //  For each possible number of buckets, try all poissible pointer
  //  widths.  First we compute the number of mers that fit in a
  //  bucket pointer table of size 2^t storing N bits in the mer data
  //  table, then we check that the number of mers in the mer data
  //  table agrees with the width of the pointer table.


  //  This is the memory size we are trying to fill, in bits.
  //
  uint64 memLimt = ((uint64)mem) << 23;

  //  Positions consume space too, but only if enabled.
  //
  uint64 posPerMer = 0;
  if (positionsEnabled)
    posPerMer = 32;

  //  Limit the number of entries in the bucket pointer table to
  //  50 bits -- thus, the prefix of each mer is at most 25.
  //
  uint32  tMax = 2*merSize - 2;
  if (tMax > 50)
    tMax = 50;

  for (uint64 t=2; t < tMax; t++) {

    //  We need to try all possibilities of N, the width of the
    //  bucket pointer table === log2(numMers).
    //
    //  Increased to 40 bits, so we're valid up to 1 trillion mers.
    //
    for (uint64 N=1; N<40; N++) {
      uint64 Nmin = uint64ONE << (N - 1);
      uint64 Nmax = uint64ONE << (N);

      //  The size in bits of the bucket pointer table.
      //
      uint64 bucketsize = (uint64ONE << t) * N;

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
        uint64 n = (memLimt - bucketsize) / (2*merSize - t + posPerMer);

        //  We can stop now if our computed number of mers is outside the range that
        //  the bucket pointer table can address.
        //
        if ((Nmin <= n) && (n <= Nmax)) {

          //fprintf(stderr, "prefixSize="uint64FMTW(2)" numMers="uint64FMTW(9)" memory=%.3fMB\n",
          //        t, n,
          //        (((uint64ONE << t) * logBaseTwo64(n) + n * (2*merSize - t + posPerMer)) >> 3) / 1048576.0);

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
    fprintf(stdout, "Can fit "uint64FMT" mers into table with prefix of "uint64FMT" bits, using %8.3fMB (%8.3fMB for positions)\n",
            maxN,
            bestT,
            (((uint64ONE << bestT) * logBaseTwo64(maxN) + maxN * (2*merSize - bestT + posPerMer)) >> 3) / 1048576.0,
            ((maxN * posPerMer) >> 3) / 1048576.0);

  return(maxN);
}



uint32
optimalNumberOfBuckets(uint32 merSize,
                       uint64 numMers,
                       bool   positionsEnabled) {
  uint64 opth   = ~uint64ZERO;
  uint64 opts   = ~uint64ZERO;
  uint64 h      = 0;
  uint64 s      = 0;
  uint64 hwidth = logBaseTwo64(numMers);

  //  Positions consume space too, but only if enabled.  Probably
  //  doesn't matter here.
  //
  uint64 posPerMer = 0;
  if (positionsEnabled)
    posPerMer = 32;

  //  Find the table size (in bits, h) that minimizes memory usage
  //  for the given merSize and numMers
  //
  //  We have two tables:
  //    the bucket pointers num buckets * pointer width   == 2 << h * hwidth
  //    the mer data:       num mers * (mersize - hwidth)
  //
  uint64 hmax = 64 - logBaseTwo64(hwidth + numMers * (2 * merSize - h));
  for (h=2; h<=hmax && h<2*merSize; h++) {
    s = (uint64ONE << h) * hwidth + numMers * (2 * merSize - h + posPerMer);

    //fprintf(stderr, "optimalNumberOfBuckets()-- h="uint64FMT" s="uint64FMT"\n", h, s);

    if (s < opts) {
      opth = h;
      opts = s;
    }
  }

  return((uint32)opth);
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

  uint32 opth = optimalNumberOfBuckets(args->merSize, args->numMersEstimated, args->positionsEnabled);
  uint64 memu = ((uint64ONE << opth)    * logBaseTwo64(args->numMersEstimated+1) +
                 args->numMersEstimated * (2 * args->merSize - opth));

  fprintf(stderr, uint64FMT" "uint32FMT"-mers can be computed using "uint64FMT"MB memory.\n",
          args->numMersEstimated, args->merSize, memu >> 23);
}

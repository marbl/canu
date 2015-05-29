#include "meryl.H"
#include "libleaff/seqCache.H"
#include "libleaff/seqStore.H"
#include "libleaff/merStream.H"

//  Takes a memory limit in MB, returns the number of mers that we can fit in that memory size,
//  assuming optimalNumberOfBuckets() below uses the same algorithm.
//
//  For each possible number of buckets, try all poissible pointer widths.  First we compute the
//  number of mers that fit in a bucket pointer table of size 2^t storing N bits in the mer data
//  table, then we check that the number of mers in the mer data table agrees with the width of the
//  pointer table.
//
uint64
estimateNumMersInMemorySize(uint32 merSize,
                            uint64 mem,
                            bool   positionsEnabled,
                            bool   beVerbose) {
  uint64 maxN    = 0;
  uint64 bestT   = 0;

  uint64 memLimt    = mem * 8;                                //  Memory limit, in bits.
  uint64 posPerMer  = (positionsEnabled == false) ? 0 : 32;   //  Positions consume space, if enabled.
  uint64 tMax       = (merSize > 25) ? 50 : 2 * merSize - 2;  //  Max width of bucket pointer table.

  //  t - prefix stored in the bucket pointer table; number of entries in the table
  //  N - width of a bucket pointer

  for (uint64 t=2; t < tMax; t++) {
    for (uint64 N=1; N<40; N++) {
      uint64 Nmin = uint64ONE << (N - 1);  //  Minimum number of mers we want to fit in the table
      uint64 Nmax = uint64ONE << (N);      //  Maximum number of mers that can fit in the table

      uint64 bucketsize = (uint64ONE << t) * N;  //  Size, in bits, of the pointer table

      uint64 n = (memLimt - bucketsize) / (2*merSize - t + posPerMer);  //  Number of mers we can fit into mer data table.

      if ((memLimt >  bucketsize) &&  //  pointer table small enough to fit in memory
          (n       >  0)          &&  //  at least some space to store mers
          (n       <= Nmax)       &&  //  enough space for the mers in the data table
          (Nmin    <= n)          &&  //  ...but not more than enough space
          (maxN    <  n)) {           //  this value of t fits more mers that any other seen so far
        maxN  = n;
        bestT = t;
      }
    }
  }

  if (beVerbose)
    fprintf(stdout, "Can fit "F_U64" mers into table with prefix of "F_U64" bits, using %.3fMB (%.3fMB for positions)\n",
            maxN,
            bestT,
            (((uint64ONE << bestT) * logBaseTwo64(maxN) + maxN * (2*merSize - bestT + posPerMer)) >> 3) / 1048576.0,
            ((maxN * posPerMer) >> 3) / 1048576.0);

  return(maxN);
}



uint64
estimateMemory(uint32 merSize,
               uint64 numMers,
               bool   positionsEnabled) {

  uint64 posPerMer = (positionsEnabled == false) ? 0 : 32;
  uint64 tMax      = (merSize > 25) ? 50 : 2 * merSize - 2;

  uint64  tMin   = tMax;
  uint64  memMin = UINT64_MAX;

  for (uint64 t=2; t < tMax; t++) {
    uint64  N       = logBaseTwo64(numMers);  //  Width of the bucket pointer table
    uint64  memUsed = ((uint64ONE << t) * logBaseTwo64(numMers) + numMers * (2 * merSize - t + posPerMer)) >> 3;

    if (memUsed < memMin) {
      tMin   = t;
      memMin = memUsed;
    }

    //fprintf(stderr, "t=%2lu N=%2lu memUsed=%16lu -- tMin=%2lu memMin=%16lu\n",
    //        t, N, memUsed, tMin, memMin);
  }

  return(memMin >> 20);
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

    //fprintf(stderr, "optimalNumberOfBuckets()-- h="F_U64" s="F_U64"\n", h, s);

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
  uint64 memu = ((uint64ONE << opth) * logBaseTwo64(args->numMersEstimated+1) + args->numMersEstimated * (2 * args->merSize - opth));

  fprintf(stderr, F_U64" "F_U32"-mers can be computed using "F_U64"MB memory.\n",
          args->numMersEstimated, args->merSize, memu >> 23);
}

#include <stdio.h>
#include <stdlib.h>

#include "bio++.H"
#include "meryl.H"


u64bit
estimateNumMersInMemorySize(u32bit merSize,
                            u32bit mem,
                            bool   beVerbose) {
  u64bit memLimt = ((u64bit)mem) << 23;
  u64bit maxN    = 0;
  u64bit bestT   = 0;

  merSize *= 2;

  //  For each possible number of buckets, try all poissible pointer
  //  widths.  First we compute the number of mers that fit in a
  //  bucket pointer table of size 2^t storing N bits in the mer data
  //  table, then we check that the number of mers in the mer data
  //  table agrees with the width of the pointer table.
  //
  for (u64bit t=2; t < merSize - 2; t++) {
    for (u64bit N=1; N<33; N++) {
      u64bit Nmin = u64bitONE << (N - 1);
      u64bit Nmax = u64bitONE << (N);

      //  The size in bits of the bucket pointer table
      //
      u64bit bucketsize = (u64bitONE << t) * N;

      //  If our bucket pointer table size hasn't already blown our
      //  memory limit, compute the number of mers that we can stuff
      //  into the list.
      //
      if (memLimt > bucketsize) {
        u64bit n = (memLimt - bucketsize) / (merSize - t);

        //  Remember the settings with the highest number of mers, if
        //  more than zero mers but not too few or too many for the
        //  bucket pointer width.
        //  
        if ((n > 0) && (Nmin <= n) && (n <= Nmax) && (maxN < n)) {
          maxN  = n;
          bestT = t;
        }

        //fprintf(stderr, "n="u64bitFMT" t="u64bitFMT"\n", n, t);
      }
    }
  }

  if (beVerbose)
    fprintf(stdout, "Can fit "u64bitFMT" mers into table with t="u64bitFMT" using "u64bitFMT"MB\n",
            maxN,
            bestT,
            ((u64bitONE << bestT) * logBaseTwo32(maxN) + maxN * (merSize - bestT)) >> 23);

  return(maxN);
}



u32bit
optimalNumberOfBuckets(u32bit merSize,
                       u64bit numMers) {
  u64bit opth   = ~u64bitZERO;
  u64bit opts   = ~u64bitZERO;
  u64bit h      = 0;
  u64bit s      = 0;
  u64bit hwidth = logBaseTwo64(numMers);

  //  Find the table size (in bits, h) that minimizes memory usage
  //  for the given merSize and numMers
  //
  //  We have two tables:
  //    the bucket pointers num buckets * pointer width   == 2 << h * hwidth
  //    the mer data:       num mers * (mersize - hwidth)
  //
  u64bit hmax = 64 - logBaseTwo64(hwidth + numMers * (2 * merSize - h));
  for (h=2; h<=hmax && h<2*merSize; h++) {
    s = (u64bitONE << h) * hwidth + numMers * (2 * merSize - h);

    if (s < opts) {
      opth = h;
      opts = s;
    }
  }

  return(opth);
}



void
estimate(merylArgs *args) {

  if (args->inputFile) {
    FastAstream        F(args->inputFile);
    merStream          M(args->merSize, &F);
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

  u32bit opth = optimalNumberOfBuckets(args->merSize, args->numMersEstimated);
  u64bit memu = ((u64bitONE << opth)    * logBaseTwo64(args->numMersEstimated+1) +
                 args->numMersEstimated * (2 * args->merSize - opth));

  fprintf(stderr, u64bitFMT" "u32bitFMT"-mers can be computed using "u64bitFMT"MB memory.\n",
          args->numMersEstimated, args->merSize, memu >> 23);
}

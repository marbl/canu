#include <stdio.h>
#include <stdlib.h>
#include "meryl.H"
#include "libbri.H"


u64bit
estimateNumMersInMemorySize(u32bit merSize,
                            u32bit mem,
                            bool   beVerbose) {
  u64bit memLimt = ((u64bit)mem) << 23;

  u64bit maxN  = 0;
  u64bit bestT = 0;

  //  For each possible number of buckets
  for (u64bit t=2; t < merSize - 2; t++) {

    //  Try all poissible number of mers
    for (u64bit N=1; N<33; N++) {

      //  The minimum and maximum number of mers that could
      //  be stored with these parameters.
      //
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

        if ((n > 0) && (Nmin <= n) && (n <= Nmax) && (maxN < n)) {
          maxN  = n;
          bestT = t;
        }
      }
    }
  }

  if (beVerbose)
    fprintf(stdout, "Can fit "u64bitFMT" mers into table with t="u64bitFMT" using "u64bitFMT"KB\n",
            maxN,
            bestT,
            ((u64bitONE << bestT) * logBaseTwo(maxN) + maxN * (merSize - bestT)) >> 3);

  return(maxN);
}



u32bit
optimalNumberOfBuckets(u32bit merSize,
                       u64bit numMers) {
  u64bit opth   = ~u64bitZERO;
  u64bit opts   = ~u64bitZERO;
  u64bit h      = 0;
  u64bit s      = 0;
  u64bit hwidth = logBaseTwo_64(numMers);

  //  Find the table size (in bits, h) that minimizes memory usage
  //  for the given merSize and numMers
  //
  //  We have two tables:
  //    the bucket pointers num buckets * pointer width   == 2 << h * hwidth
  //    the mer data:       num mers * (mersize - hwidth)
  //    
  for (h=2; h<=32 && h<2*merSize; h++) {
    s = (u64bitONE << h) * h + numMers * (2 * merSize - hwidth);

    //fprintf(stderr, "numBuckets_log2 = "u64bitFMT" -- memory = "u64bitFMT"\n", h, s);

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
    merStream          M(args->merSize, args->inputFile);
    speedCounter       C(" %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);

    if (args->beVerbose)
      fprintf(stderr, "Counting mers in '%s'\n", args->inputFile);

    args->numMersEstimated = 0;

    while (M.nextMer()) {
      C.tick();
      args->numMersEstimated++;
    }

    if (args->beVerbose)
      fprintf(stderr, "\n");
  }


  //  How many bits do we need in the hash table to store all the mers?
  //
  u32bit hBits = 1;
  while ((args->numMersEstimated+1) > (u64bitONE << hBits))
    hBits++;

  u64bit h, hSize, c, cSize;


  //  Find the optimial memory settings, so we can mark it in the output
  //
  u32bit opth = optimalNumberOfBuckets(args->merSize, args->numMersEstimated);


  fprintf(stderr, "For "u64bitFMT" mers ("u32bitFMT" bits/hash), the optimal memory usage is achieved\n", args->numMersEstimated, hBits);
  fprintf(stderr, "by picking the 'h' with the smallest 't'.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " h |  h MB |  c |  c MB |  t MB\n");
  fprintf(stderr, "---+-------+----+-------+------\n");

  for (h=16; h<=32 && h<2*args->merSize; h++) {
    c     = 2 * args->merSize - h;

    hSize = (u64bitONE << h) * hBits;
    cSize = args->numMersEstimated * c;

    hSize >>= 23;
    cSize >>= 23;

#ifdef TRUE64BIT
    fprintf(stderr, "%2lu | %5lu | %2lu | %5lu | %5lu%s\n",
            h, hSize,
            c, cSize,
            hSize + cSize,
            (h==opth) ? " <- smallest memory" : "");
#else
    fprintf(stderr, "%2llu | %5llu | %2llu | %5llu | %5llu%s\n",
            h, hSize,
            c, cSize,
            hSize + cSize,
            (h==opth) ? " <- smallest memory" : "");
#endif
  }
}

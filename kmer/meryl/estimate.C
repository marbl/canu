#include <stdio.h>
#include <stdlib.h>
#include "meryl.H"
#include "libbri.H"


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
            ((u64bitONE << bestT) * logBaseTwo(maxN) + maxN * (merSize - bestT)) >> 23);

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
  for (h=2; h<=64 && h<2*merSize; h++) {
    s = (u64bitONE << h) * hwidth + numMers * (2 * merSize - h);

    //fprintf(stderr, "optimalNumberfOfBuckets()-- numBuckets_log2: "u64bitFMT"  memory: "u64bitFMT"MB\n", h, s >> 23);

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

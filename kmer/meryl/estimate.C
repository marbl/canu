#include <stdio.h>
#include <stdlib.h>
#include "libbri.H"
#include "britime.H"


u32bit
partition(u32bit merSize,
          u32bit mem,
          bool   beVerbose) {
  u64bit memLimt = ((u64bit)mem) << 23;

  u32bit maxN  = 0;
  u32bit bestT = 0;

  for (u64bit t=2; t < merSize - 2; t++) {
    for (u64bit N=1; N<33; N++) {
      u64bit bucketsize = (u64bitONE << t) * N;
      u64bit Nmin = u64bitONE << (N - 1);
      u64bit Nmax = u64bitONE << (N);

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
    fprintf(stdout, "Can fit %10lu mers into table with t=%8lu using %10luKB\n",
            maxN,
            bestT,
            ((u64bitONE << bestT) * logBaseTwo(maxN) + maxN * (merSize - bestT)) >> 3);

  return(bestT);
}


void
estimate(char   *inputFile,
         u32bit  merSize,
         u64bit  numMers,
         bool    beVerbose) {

  if (inputFile) {
    merStream          M(merSize, inputFile);
    speedCounter       C(" %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);

    if (beVerbose)
      fprintf(stderr, "Counting mers in '%s'\n", inputFile);

    numMers = 0;

    while (M.nextMer()) {
      C.tick();
      numMers++;
    }

    if (beVerbose)
      fprintf(stderr, "\n");
  }


  //  How many bits do we need in the hash table to store all the mers?
  //
  u32bit hBits = 1;
  while ((numMers+1) > (u64bitONE << hBits))
    hBits++;

  u64bit h, hSize, c, cSize;


  //  Find the optimial memory settings, so we can mark it in the output
  //
  u32bit opth = ~u32bitZERO;
  u32bit opts = ~u32bitZERO;

  for (h=16; h<=32 && h<2*merSize; h++) {
    c     = 2 * merSize - h;

    hSize = (u64bitONE << h) * hBits;
    cSize = numMers * c;

    hSize >>= 23;
    cSize >>= 23;

    if (hSize + cSize < opts) {
      opth = h;
      opts = hSize + cSize;
    }
  }

  fprintf(stderr, "For "u64bitFMT" mers ("u32bitFMT" bits/hash), the optimal memory usage is achieved\n", numMers, hBits);
  fprintf(stderr, "by picking the 'h' with the smallest 't'.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " h |  h MB |  c |  c MB |  t MB\n");
  fprintf(stderr, "---+-------+----+-------+------\n");

  for (h=16; h<=32 && h<2*merSize; h++) {
    c     = 2 * merSize - h;

    hSize = (u64bitONE << h) * hBits;
    cSize = numMers * c;

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

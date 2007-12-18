#include "util++.H"
#include "bio++.H"


//  These should add up to 2 * mersize - hash size
//
//  And they should all be odd - seems to reduce collisions, and
//  collisions that occur are obviously different sequence.  If any
//  one value is even, then we get collisions that are 2 mismatches
//  different.
//
u32bit  shiftValues[6] = { 25, 33, 27, 29, 27, 21 };


//  Hash a kmer down to n bits.
u64bit
hash(kMer k, u32bit n) {
  u64bit  h = u64bitZERO;

  h  ^= k.getWord(0);
  k >>= shiftValues[0];

  h  ^= ~k.getWord(0);
  k >>= shiftValues[1];

  h  ^= k.getWord(0);
  k >>= shiftValues[2];

  h  ^= ~k.getWord(0);
  k >>= shiftValues[3];

  h  ^= k.getWord(0);
  k >>= shiftValues[4];

  h  ^= ~k.getWord(0);
  k >>= shiftValues[5];

  h  ^= k.getWord(0);

  h &= u64bitMASK(n);

  return(h);
}



//  Same as hash(), but prints out debug info if
//  the hash is all 1's.
//
u64bit hashdebug(kMer k, u32bit n) {

  u64bit h = hash(k, n);

  if (h == u64bitNUMBER(0x0000007fffffffff)) {
    u32bit  shift = n / 2;
    u32bit  iters = 4 * k.getMerSize() / n;

    u64bit  d = k.getWord(0);
    fprintf(stderr, "d["u32bitFMT"] "u64bitHEX" from "u64bitHEX"\n", iters, d, k.getWord(0));

    while (iters--) {
      k >>= shift;
      d ^= k.getWord(0);
      fprintf(stderr, "d["u32bitFMT"] "u64bitHEX" from "u64bitHEX"\n", iters, d, k.getWord(0));
    }
  }
    
  return(h);
}




int
main(int argc, char **argv) {

  FastAstream  *FS = new FastAstream("/home/work/SEQUENCE/B35LC/B35LC.fasta");
  merStream    *MS = new merStream(100, FS);
  speedCounter *SC = new speedCounter(" %8f Mbp (%8.5f Mbp/sec)\r", 1000000, 1000000, true);

  char          merstring[1024];

  while (MS->nextMer()) {
    SC->tick();

    u64bit  h = hash(MS->theFMer(), 39);

    if ((h >> 24) == u64bitNUMBER(0x0000005da3))
      fprintf(stdout, u64bitHEX" %s\n", h, MS->theFMer().merToString(merstring));
  }

  delete MS;
  delete FS;

  return(0);
}

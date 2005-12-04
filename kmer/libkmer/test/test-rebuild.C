#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"

//  Tests a positionDB when using an existDB for masking.
//
//  existDB can be either include or exclude
//  positionDB can use include, exclude or threshold
//

int
main(int argc, char **argv) {

  if (argc != 3) {
    fprintf(stderr, "usage: %s mersize tblsize\n", argv[0]);
    exit(1);
  }

  u32bit  merSize = strtou32bit(argv[1], 0L);
  u32bit  tblSize = strtou32bit(argv[2], 0L);
  u64bit  maxMers = u64bitONE << 25;

#define TRY_ALL
#ifdef TRY_ALL
  for (merSize=8; merSize<33; merSize++) {
    for (tblSize = merSize/2; tblSize < ((7*merSize/4 > 27) ? 27 : 7*merSize/4); tblSize++) {
#endif

      fprintf(stderr, "Testing "u64bitFMT" Mmers at merSize "u32bitFMT" tblSize "u32bitFMT".\n",
              maxMers, merSize, tblSize);

      merStream   *T  = new merStream(merSize, "acgcgactcgagctacgagcgatcacgacgactacgagca", 40);
      positionDB  *P  = new positionDB(T, merSize, 0, tblSize, 0L, 0L, 0, false);

      u64bit p = 0;
      u64bit f = 0;

      speedCounter  *C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, true);

      mt_s   *mts = mtInit(3492);
      u64bit  msk = u64bitMASK(2*merSize);
      u64bit  cnt = maxMers;

      while (cnt--) {
        if (P->checkREBUILD(mtRandom64(mts) & msk) == false) {
          f++;
          fprintf(stderr, "PASS: "u64bitFMT"\n", p);
          fprintf(stderr, "FAIL: "u64bitFMT"\n", f);
          exit(1);
        } else {
          p++;
        }

        C->tick();
      }

      free(mts);

      delete C;
      delete P;
      delete T;
#ifdef TRY_ALL
    }
  }
#endif

  exit(0);
}

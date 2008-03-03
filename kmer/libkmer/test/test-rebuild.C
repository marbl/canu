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

  u64bit  maxMers = u64bitONE << 25;

  for (u32bit merSize=8; merSize<33; merSize++) {
    fprintf(stderr, "Testing "u64bitFMT" Mmers at merSize "u32bitFMT".\n", maxMers, merSize);

    kMerBuilder *K  = new kMerBuilder(merSize);
    merStream   *T  = new merStream(K, "acgcgactcgagctacgagcgatcacgacgactacgagca", 40);
    positionDB  *P  = new positionDB(T, merSize, 0, 0L, 0L, 0L, 0, 0, false, true);

    u64bit p = 0;
    u64bit f = 0;

    mt_s   *mts = mtInit(3492);
    u64bit  msk = u64bitMASK(2*merSize);
    u64bit  cnt = maxMers;

    while (cnt--) {
      if (P->checkREBUILD(mtRandom64(mts) & msk) == false) {
        f++;
      } else {
        p++;
      }
    }

    if (f) {
      fprintf(stderr, "PASS: "u64bitFMT"  FAIL: "u64bitFMT"\n", p, f);
      exit(1);
    }

    free(mts);

    delete P;
    delete T;
  }

  exit(0);
}

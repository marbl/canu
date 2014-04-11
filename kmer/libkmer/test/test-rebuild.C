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

  uint64  maxMers = uint64ONE << 25;

  for (uint32 merSize=8; merSize<33; merSize++) {
    fprintf(stderr, "Testing "uint64FMT" Mmers at merSize "uint32FMT".\n", maxMers, merSize);

    kMerBuilder *K  = new kMerBuilder(merSize);
    merStream   *T  = new merStream(K, "acgcgactcgagctacgagcgatcacgacgactacgagca", 40);
    positionDB  *P  = new positionDB(T, merSize, 0, 0L, 0L, 0L, 0, 0, false, true);

    uint64 p = 0;
    uint64 f = 0;

    mt_s   *mts = mtInit(3492);
    uint64  msk = uint64MASK(2*merSize);
    uint64  cnt = maxMers;

    while (cnt--) {
      if (P->checkREBUILD(mtRandom64(mts) & msk) == false) {
        f++;
      } else {
        p++;
      }
    }

    if (f) {
      fprintf(stderr, "PASS: "uint64FMT"  FAIL: "uint64FMT"\n", p, f);
      exit(1);
    }

    free(mts);

    delete P;
    delete T;
  }

  exit(0);
}

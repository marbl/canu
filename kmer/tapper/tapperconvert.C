#include "tapperTag.H"
#include "tapperHit.H"
#include "tapperGlobalData.H"
#include "tapperThreadData.H"
#include "tapperComputation.H"

int
main(int argc, char **argv) {
  char *resultName = 0L;

  bool  dumpIndex = false;
  bool  dumpFrag  = false;
  bool  dumpSing  = false;
  bool  dumpMate  = false;
  bool  dumpTang  = false;

  bool  allIndex  = false;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-dumpindex", 6) == 0) {
      dumpIndex = true;

    } else if (strncmp(argv[arg], "-dumpfragments", 6) == 0) {
      dumpFrag = true;

    } else if (strncmp(argv[arg], "-dumpsingleton", 6) == 0) {
      dumpSing = true;

    } else if (strncmp(argv[arg], "-dumpmated", 6) == 0) {
      dumpMate = true;

    } else if (strncmp(argv[arg], "-dumptangled", 6) == 0) {
      dumpTang = true;

    } else if (strncmp(argv[arg], "-allindex", 6) == 0) {
      allIndex = true;

    } else if (resultName == 0L) {
      resultName = argv[arg];

    } else {
        err++;
    }

    arg++;
  }
  if ((err) || (resultName == 0L)) {
    fprintf(stderr, "usage: %s [-dumpindex] [-dumpfragments] [-dumpsingletons] [-dumpmated] [-dumptangled] prefix\n", argv[0]);
    fprintf(stderr, "       -allIndex -- also dump index for unmapped fragments\n");
    exit(1);
  }

  tapperAlignmentFile   *AF = new tapperAlignmentFile(resultName, 'r');
  tapperAlignment        AL;

  while (AF->read(&AL)) {
    if ((dumpIndex) &&
        ((allIndex) ||
         ((dumpFrag) && (AL.idx._numFragment           > 0)) ||
         ((dumpFrag) && (AL.idx._numFragmentDiscarded  > 0)) ||
         ((dumpSing) && (AL.idx._numSingleton          > 0)) ||
         ((dumpMate) && (AL.idx._numMated              > 0)) ||
         ((dumpTang) && (AL.idx._numTangled            > 0))))
      AL.idx.print(stdout);

    if (dumpFrag)
      for (u32bit i=0; i<AL.idx._numFragment; i++)
        AL.frag[i].print(stdout, &AL.idx);

    if (dumpSing)
      for (u32bit i=0; i<AL.idx._numSingleton; i++)
        AL.sing[i].print(stdout, &AL.idx);

    if (dumpMate)
      for (u32bit i=0; i<AL.idx._numMated; i++)
        AL.mate[i].print(stdout, &AL.idx);

    if (dumpTang)
      for (u32bit i=0; i<AL.idx._numTangled; i++)
        AL.tang[i].print(stdout, &AL.idx);
  }

  delete AF;

  exit(0);
}

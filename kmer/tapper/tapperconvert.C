#include "tapperTag.H"
#include "tapperResult.H"
#include "tapperAlignment.H"
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

  tapperResultFile   *inp = new tapperResultFile(resultName, 'r');
  tapperResult       *res = new tapperResult;

  while (inp->read(res)) {
    if ((dumpIndex) &&
        ((allIndex) ||
         ((dumpFrag) && (res->idx._numFrag           > 0)) ||
         ((dumpFrag) && (res->idx._numFragDiscarded  > 0)) ||
         ((dumpSing) && (res->idx._numFragSingleton  > 0)) ||
         ((dumpMate) && (res->idx._numMated          > 0)) ||
         ((dumpTang) && (res->idx._numTangled        > 0))))
      res->idx.print(stdout);

    if (dumpFrag)
      for (u32bit i=0; i<res->idx._numFrag; i++)
        res->frag[i].print(stdout, &res->idx);

    if (dumpSing)
      for (u32bit i=0; i<res->idx._numFragSingleton; i++)
        res->sing[i].print(stdout, &res->idx);

    if (dumpMate)
      for (u32bit i=0; i<res->idx._numMated; i++)
        res->mate[i].print(stdout, &res->idx);

    if (dumpTang)
      for (u32bit i=0; i<res->idx._numTangled; i++)
        res->tang[i].print(stdout, &res->idx);
  }

  delete inp;
  delete res;

  exit(0);
}

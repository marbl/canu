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

  char                   fileName[FILENAME_MAX];
  tapperResult           res;
  tapperResultFragment   frag;
  tapperResultSingleton  sing;
  tapperResultMated      mate;
  tapperResultTangled    tang;

  sprintf(fileName, "%s.tapperMappedIndex", resultName);
  recordFile *IDX = new recordFile(fileName, 0, sizeof(tapperResult), 'r');

  sprintf(fileName, "%s.tapperMappedFragment", resultName);
  recordFile *FRAG = new recordFile(fileName, 0, sizeof(tapperResultFragment), 'r');

  sprintf(fileName, "%s.tapperMappedSingleton", resultName);
  recordFile *SING = new recordFile(fileName, 0, sizeof(tapperResultSingleton), 'r');

  sprintf(fileName, "%s.tapperMappedMated", resultName);
  recordFile *MATE = new recordFile(fileName, 0, sizeof(tapperResultMated), 'r');

  sprintf(fileName, "%s.tapperMappedTangled", resultName);
  recordFile *TANG = new recordFile(fileName, 0, sizeof(tapperResultTangled), 'r');


  while (IDX->getRecord(&res) == 1) {
    if ((dumpIndex) &&
        ((allIndex) ||
         ((dumpFrag) && (res._numFragment           > 0)) ||
         ((dumpFrag) && (res._numFragmentDiscarded  > 0)) ||
         ((dumpSing) && (res._numSingleton          > 0)) ||
         ((dumpMate) && (res._numMated              > 0)) ||
         ((dumpTang) && (res._numTangled            > 0))))
      res.print(stdout);

    if (dumpFrag) {
      for (u32bit i=0; i<res._numFragment; i++) {
        if (FRAG->getRecord(&frag) != 1)
          fprintf(stderr, "Failed to read FRAG record.\n"), exit(1);
        frag.print(stdout, &res);
      }
    }

    if (dumpSing) {
      for (u32bit i=0; i<res._numSingleton; i++) {
        if (SING->getRecord(&sing) != 1)
          fprintf(stderr, "Failed to read SING record.\n"), exit(1);
        sing.print(stdout, &res);
      }
    }

    if (dumpMate) {
      for (u32bit i=0; i<res._numMated; i++) {
        if (MATE->getRecord(&mate) != 1)
          fprintf(stderr, "Failed to read MATE record.\n"), exit(1);
        mate.print(stdout, &res);
      }
    }

    if (dumpTang) {
      for (u32bit i=0; i<res._numTangled; i++) {
        if (TANG->getRecord(&tang) != 1)
          fprintf(stderr, "Failed to read TANG record.\n"), exit(1);
        tang.print(stdout, &res);
      }
    }
  }

  delete IDX;
  delete FRAG;
  delete SING;
  delete MATE;
  delete TANG;

  exit(0);
}

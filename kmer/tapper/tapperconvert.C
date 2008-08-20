#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"

#include "tapperTag.H"
#include "tapperHit.H"
#include "tapperGlobalData.H"
#include "tapperThreadData.H"
#include "tapperComputation.H"

int
main(int argc, char **argv) {
  char *resultName = 0L;

  recordFile  *IDX = 0L;
  recordFile  *DAT = 0L;

  bool  dumpFrag = false;
  bool  dumpMate = false;
  bool  dumpSing = false;
  bool  dumpTang = false;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-dumpfragments", 6) == 0) {
      resultName = argv[++arg];
      dumpFrag = true;

    } else if (strncmp(argv[arg], "-dumpmated", 6) == 0) {
      resultName = argv[++arg];
      dumpMate = true;

    } else if (strncmp(argv[arg], "-dumpsingleton", 6) == 0) {
      resultName = argv[++arg];
      dumpSing = true;

    } else if (strncmp(argv[arg], "-dumptangled", 6) == 0) {
      resultName = argv[++arg];
      dumpTang = true;

    } else {
      err++;
    }

    arg++;
  }
  if ((err) || (resultName == 0L)) {
    fprintf(stderr, "usage: %s -dumpXXX prefix\n", argv[0]);
    fprintf(stderr, "       XXX == fragments, mated, singleton or tangled\n");
    exit(1);
  }


  if (dumpFrag) {
    char                   fileName[FILENAME_MAX];
    tapperResult           res;
    tapperResultFragment   frag;
    u16bit                 id[4];

    sprintf(fileName, "%s.tapperMappedIndex", resultName);
    IDX = new recordFile(fileName, 0, sizeof(tapperResult), 'r');

    sprintf(fileName, "%s.tapperMappedFragment", resultName);
    DAT = new recordFile(fileName, 0, sizeof(tapperResultFragment), 'r');

    while (IDX->getRecord(&res) == 1) {
      for (u32bit i=0; i<res._numFragment; i++) {
        if (DAT->getRecord(&frag) != 1)
          fprintf(stderr, "Failed to read DAT record.\n"), exit(1);

        if (frag._tag._tag1)
          decodeTagID(res._tag1id, id);
        else
          decodeTagID(res._tag2id, id);

        fprintf(stdout, "F "u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT" "u32bitFMT" "u32bitFMT" %c "u32bitFMT" "u32bitFMT"/"u32bitFMT"/"u32bitFMT"\n",
                id[0], id[1], id[2], id[3],
                frag._seq,
                frag._pos,
                frag._tag._rev ? 'r' : 'f',
                frag._tag._rank,
                frag._tag._basesMismatch,
                frag._tag._colorMismatch,
                frag._tag._colorInconsistent);
      }
    }
  }


  if (dumpMate) {
    char                   fileName[FILENAME_MAX];
    tapperResult           res;
    tapperResultMated      mate;
    u16bit                 id1[4];
    u16bit                 id2[4];

    sprintf(fileName, "%s.tapperMappedIndex", resultName);
    IDX = new recordFile(fileName, 0, sizeof(tapperResult), 'r');

    sprintf(fileName, "%s.tapperMappedMated", resultName);
    DAT = new recordFile(fileName, 0, sizeof(tapperResultMated), 'r');

    while (IDX->getRecord(&res)) {
      for (u32bit i=0; i<res._numMated; i++) {
        if (DAT->getRecord(&mate) != 1)
          fprintf(stderr, "Failed to read DAT record.\n"), exit(1);

        decodeTagID(res._tag1id, id1);
        decodeTagID(res._tag2id, id2);

        fprintf(stdout, "M "u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT" "u32bitFMT" "u32bitFMT" %c "u32bitFMT" "u32bitFMT"/"u32bitFMT"/"u32bitFMT" "u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT" "u32bitFMT" "u32bitFMT" %c "u32bitFMT" "u32bitFMT"/"u32bitFMT"/"u32bitFMT"\n",
                id1[0], id1[1], id1[2], id1[3],
                mate._seq,
                mate._pos1,
                mate._tag1._rev ? 'r' : 'f',
                mate._tag1._rank,
                mate._tag1._basesMismatch,
                mate._tag1._colorMismatch,
                mate._tag1._colorInconsistent,
                id2[0], id2[1], id2[2], id2[3],
                mate._seq,
                mate._pos2,
                mate._tag2._rev ? 'r' : 'f',
                mate._tag2._rank,
                mate._tag2._basesMismatch,
                mate._tag2._colorMismatch,
                mate._tag2._colorInconsistent);
      }
    }
  }


  if (dumpSing) {
    char                   fileName[FILENAME_MAX];
    tapperResult           res;
    tapperResultSingleton  sing;
    u16bit                 id[4];

    sprintf(fileName, "%s.tapperMappedIndex", resultName);
    IDX = new recordFile(fileName, 0, sizeof(tapperResult), 'r');

    sprintf(fileName, "%s.tapperMappedSingleton", resultName);
    DAT = new recordFile(fileName, 0, sizeof(tapperResultSingleton), 'r');

    while (IDX->getRecord(&res)) {
      for (u32bit i=0; i<res._numSingleton; i++) {
        if (DAT->getRecord(&sing) != 1)
          fprintf(stderr, "Failed to read DAT record.\n"), exit(1);

        if (sing._tag._tag1)
          decodeTagID(res._tag1id, id);
        else
          decodeTagID(res._tag2id, id);

        fprintf(stdout, "F "u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT" "u32bitFMT" "u32bitFMT" %c "u32bitFMT" "u32bitFMT"/"u32bitFMT"/"u32bitFMT"\n",
                id[0], id[1], id[2], id[3],
                sing._seq,
                sing._pos,
                sing._tag._rev ? 'r' : 'f',
                sing._tag._rank,
                sing._tag._basesMismatch,
                sing._tag._colorMismatch,
                sing._tag._colorInconsistent);
      }
    }
  }


  if (dumpTang) {
    char                   fileName[FILENAME_MAX];
    tapperResult           res;
    tapperResultTangled    tang;
    u16bit                 id1[4];
    u16bit                 id2[4];

    sprintf(fileName, "%s.tapperMappedIndex", resultName);
    IDX = new recordFile(fileName, 0, sizeof(tapperResult), 'r');

    sprintf(fileName, "%s.tapperMappedTangled", resultName);
    DAT = new recordFile(fileName, 0, sizeof(tapperResultTangled), 'r');

    while (IDX->getRecord(&res)) {
      for (u32bit i=0; i<res._numTangled; i++) {
        if (DAT->getRecord(&tang) != 1)
          fprintf(stderr, "Failed to read DAT record.\n"), exit(1);

        decodeTagID(res._tag1id, id1);
        decodeTagID(res._tag2id, id2);

        fprintf(stdout, "T "u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT" "u32bitFMT" "u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
                id1[0], id1[1], id1[2], id1[3],
                tang._tag1count,
                id2[0], id2[1], id2[2], id2[3],
                tang._tag2count,
                tang._seq,
                tang._bgn,
                tang._end);
      }
    }
  }


  exit(0);
}

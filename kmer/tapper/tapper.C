#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"

//  Tests a positionDB, looking for hits with three errors.

void
encodeToColor(char *seq, u32bit len) {
}

int
main(int argc, char **argv) {
  char   *genName   = 0L;
  char   *qryName   = 0L;
  u32bit  merSize   = 25;
  u32bit  maxError  = 3;
  bool    beVerbose = false;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-genomic") == 0) {
      genName = argv[++arg];
    } else if (strcmp(argv[arg], "-queries") == 0) {
      qryName = argv[++arg];

    } else if (strcmp(argv[arg], "-mersize") == 0) {
      merSize = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-maxerror") == 0) {
      maxError = strtou32bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-genomic") == 0) {
    } else if (strcmp(argv[arg], "-genomic") == 0) {

    } else if (strcmp(argv[arg], "-verbose") == 0) {
      beVerbose = 1;
    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err > 0) || (genName == 0L) || (qryName == 0L)) {
    fprintf(stderr, "usage: %s -genomic g.fasta -queries q.fasta\n", argv[0]);
    exit(1);
  }

  seqStream  *SS = new seqStream(genName, true);
  merStream  *MS = new merStream(new kMerBuilder(merSize), SS);
  positionDB *PS = new positionDB(MS, merSize, 0, 0L, 0L, 0, beVerbose, true);
  delete MS;  MS=0L;
  delete SS;  SS=0L;

  seqStream  *SQ = new seqStream(qryName, true);
  merStream  *MQ = new merStream(new kMerBuilder(merSize), SQ);

  u64bit *posn    = 0L;
  u64bit  posnMax = 0;
  u64bit  posnLen = 0;

  char    themer[256];

  u64bit  numMers = 0;

  double  startTime = getTime();

  while (MQ->nextMer()) {
    if (PS->getMismatch(MQ->theFMer(), 4, posn, posnMax, posnLen)) {
      fprintf(stderr, u64bitFMT"] found mer %s "u64bitFMT" times; "u64bitFMT" (%f/sec).\n",
              numMers,
              MQ->theFMer().merToString(themer),
              posnLen, posn[0],
              numMers / (0.00001 + getTime() - startTime));
    }

    numMers++;
  }
  
  delete MQ;
  delete SQ;

  delete PS;
  delete MS;
  delete SS;

  exit(0);
}

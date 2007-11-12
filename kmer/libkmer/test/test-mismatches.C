#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"

//  Tests a positionDB, looking for hits with three errors.

//#define MERSIZE 24   //  mersize in bases
//#define TBLSIZE 22   //  table size in bits

#define MERSIZE 27   //  mersize in bases

int
main(int argc, char **argv) {

  if (argc != 3) {
    fprintf(stderr, "usage: %s seq.fasta queries.fasta\n", argv[0]);
    exit(1);
  }

  char *seqName = argv[1];
  char *qryName = argv[2];

  seqStream  *SS = new seqStream(seqName, true);
  merStream  *MS = new merStream(new kMerBuilder(MERSIZE), SS);
  positionDB *PS = new positionDB(MS, MERSIZE, 0, 0L, 0L, 0, true, true);
  delete MS;  MS=0L;
  delete SS;  SS=0L;

  seqStream  *SQ = new seqStream(qryName, true);
  merStream  *MQ = new merStream(new kMerBuilder(MERSIZE), SQ);

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

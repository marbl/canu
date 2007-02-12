#include "trim.H"

//  Read a fragStore, writes multifasta of the trimmed sequences.

//  XXX:  Should use a fragStream

int
main(int argc, char **argv) {
  char   *gkpStore   = 0L;
  u32bit  seqMax = 10240;
  u32bit  qltLen = 0;
  char   *seq    = new char   [seqMax];
  char   *qlt    = new char   [seqMax];
  bool    allbases = false;
  bool    allfrags = false;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-frg", 2) == 0) {
      gkpStore = argv[++arg];
    } else if (strncmp(argv[arg], "-allbases", 5) == 0) {
      allbases = true;
    } else if (strncmp(argv[arg], "-allfrags", 5) == 0) {
      allfrags = true;
    }
    arg++;
  }

  if (!gkpStore) {
    fprintf(stderr, "usage: %s [-allbases] [-allfrags] -frg some.gkpStore\n", argv[0]);
    fprintf(stderr, "  -allbases      Print all the sequence, not just the clear range.\n");
    fprintf(stderr, "  -allfrags      Print all the fragments, including deleted ones.\n");
    exit(1);
  }

  //  Open the store
  //
  GateKeeperStore   *gkp = openGateKeeperStore(gkpStore, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStore);
    exit(1);
  }

  u32bit   firstElem = getFirstElemFragStore(gkp);
  u32bit   lastElem  = getLastElemFragStore(gkp) + 1;

  ReadStruct       *rd = new_ReadStruct();
  unsigned int      deleted = 0;
  unsigned int      clrBeg  = 0;
  unsigned int      clrEnd  = 0;

  for (u32bit elem=firstElem; elem<lastElem; elem++) {
    getFrag(gkp, elem, rd, FRAG_S_SEQ);

    getIsDeleted_ReadStruct(rd, &deleted);

    if (allfrags || !deleted) {
      Fragment_ID      uid = 0;

      getClearRegion_ReadStruct(rd, &clrBeg, &clrEnd, READSTRUCT_OVL);
      getAccID_ReadStruct(rd, &uid);

      if (getSequence_ReadStruct(rd, seq, qlt, seqMax)) {
        fprintf(stderr, "getSequence_ReadStruct() failed.\n");
        exit(1);
      }

      if (allbases) {
        fprintf(stdout, ">"u64bitFMT","u32bitFMT" beg="u32bitFMT" end="u32bitFMT" trimmed=0 deleted=%d\n%s\n",
                (u64bit)uid, elem, clrBeg, clrEnd, deleted,
                seq);
      } else {
        seq[clrEnd] = 0;
        fprintf(stdout, ">"u64bitFMT","u32bitFMT" beg="u32bitFMT" end="u32bitFMT" trimmed=1 deleted=%d\n%s\n",
                (u64bit)uid, elem, clrBeg, clrEnd, deleted,
                seq + clrBeg);
      }
    }
  }

  delete [] seq;
  delete [] qlt;

  closeGateKeeperStore(gkp);
}

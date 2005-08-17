#include "trim.H"

//  Read a fragStore, writes multifasta of the trimmed sequences.

//  XXX:  Should use a fragStream

int
main(int argc, char **argv) {
  char   *frgStore   = 0L;
  u32bit  seqMax = 10240;
  u32bit  qltLen = 0;
  char   *seq    = new char   [seqMax];
  char   *qlt    = new char   [seqMax];

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-frg", 2) == 0) {
      frgStore = argv[++arg];
    }
    arg++;
  }

  if (!frgStore) {
    fprintf(stderr, "usage: %s -frg some.frgStore\n", argv[0]);
    exit(1);
  }

  //  Open the store
  //
  FragStoreHandle   fs = openFragStore(frgStore, "rw+");
  if (fs == NULLSTOREHANDLE) {
    fprintf(stderr, "Failed to open %s\n", frgStore);
    exit(1);
  }

  u32bit   firstElem = getFirstElemFragStore(fs);
  u32bit   lastElem  = getLastElemFragStore(fs) + 1;

  ReadStructp       rd = new_ReadStruct();
  unsigned int      deleted = 0;
  unsigned int      clrBeg  = 0;
  unsigned int      clrEnd  = 0;

  for (u32bit elem=firstElem; elem<lastElem; elem++) {
    getFragStore(fs, elem, FRAG_S_ALL, rd);

    getIsDeleted_ReadStruct(rd, &deleted);

    if (!deleted) {
      getClearRegion_ReadStruct(rd, &clrBeg, &clrEnd, READSTRUCT_OVL);

      if (getSequence_ReadStruct(rd, seq, qlt, seqMax)) {
        fprintf(stderr, "getSequence_ReadStruct() failed.\n");
        exit(1);
      }

      seq[clrEnd] = 0;
      fprintf(stdout, ">"u32bitFMT" beg="u32bitFMT" end="u32bitFMT"\n%s\n",
              elem, clrBeg, clrEnd,
              seq + clrBeg);
    }
  }

  delete [] seq;
  delete [] qlt;

  closeFragStore(fs);
}

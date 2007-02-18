#include "trim.H"

//  Read a fragStore, reports duplicate fragments, removes all but one
//  copy of the duplicate.

void
usage(char *name) {
  fprintf(stderr, "usage: %s -frg some.gkpStore\n", name);
  fprintf(stderr, "\n");
  fprintf(stderr, "  -frg       Operate on this gkpStore\n");
  fprintf(stderr, "  -log       Report the iid, original trim and new quality trim\n");
}


//  12 bytes per frag, 25,000,000 frags will need about 300MB to run.
//
struct fragHash {
  u64bit    hash;
  u32bit    iid;
  u32bit    mate;
  u32bit    dupcount:8;
  u32bit    killme:1;
};


int
fragHashCompare(const void *a, const void *b) {
  fragHash const *A = (fragHash const *)a;
  fragHash const *B = (fragHash const *)b;

  if (A->hash < B->hash) return(-1);
  if (A->hash > B->hash) return( 1);
  return(0);
}


int
main(int argc, char **argv) {
  char   *gkpStore   = 0L;
  FILE   *logFile    = 0L;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-frg", 2) == 0) {
      gkpStore = argv[++arg];
    } else if (strncmp(argv[arg], "-log", 2) == 0) {
      errno=0;
      logFile = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open %s for writing the log: %s\n", argv[arg], strerror(errno)), exit(1);
    } else {
      fprintf(stderr, "Invalid option: '%s'\n", argv[arg]);
      usage(argv[0]);
      exit(1);
    }
    arg++;
  }

  if (!gkpStore) {
    usage(argv[0]);
    exit(1);
  }

  //  Open the store
  //
  GateKeeperStore  *gkp = openGateKeeperStore(gkpStore, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStore);
    exit(1);
  }

  u32bit   firstElem = getFirstElemFragStore(gkp);
  u32bit   lastElem  = getLastElemFragStore(gkp) + 1;

  fragRecord       *fr1 = new_fragRecord();
  fragRecord       *fr2 = new_fragRecord();

  ////////////////////////////////////////

  fragHash   *fh = new fragHash [lastElem - firstElem + 1];
  u32bit      seqMax = 10240;
  char       *seq1   = NULL;
  char       *qlt1   = NULL;
  char       *seq2   = NULL;
  char       *qlt2   = NULL;

  ////////////////////////////////////////

  fprintf(stderr, "Read "u32bitFMT" fragments to build hashes.\n", lastElem - firstElem + 1);

  for (u32bit elem=firstElem; elem<lastElem; elem++) {
    getFrag(gkp, elem, fr1, FRAG_S_SEQ);
    seq1 = getFragRecordSequence(fr1);

    u32bit seqLen   = getFragRecordSequenceLength(fr1);
    u64bit hash     = 0;
    u32bit map[256] = { 0 };

    for (u32bit s=0; s<256; s++)
      map[s] = 1;
    map['A'] = map['a'] = 2;
    map['C'] = map['c'] = 3;
    map['G'] = map['g'] = 4;
    map['T'] = map['t'] = 5;

    for (u64bit s=0; s<seqLen; s++) {
      hash  ^= s * (u64bit)(map[seq1[s]]) * (u64bit)(qlt1[s] - '0');
      hash   = (hash << 5) | (hash >> 59);
    }

    fh[elem-firstElem].hash     = hash;
    fh[elem-firstElem].iid      = elem;
    fh[elem-firstElem].mate     = 0;
    fh[elem-firstElem].dupcount = 0;
    fh[elem-firstElem].killme   = 0;
  }

  ////////////////////////////////////////

  fprintf(stderr, "Sort hashes.\n");

  qsort(fh, lastElem-firstElem+1, sizeof(fragHash), fragHashCompare);

  ////////////////////////////////////////

  fprintf(stderr, "Examine hashes to find collisiosn.\n");

  u32bit   reallyDup;
  u32bit   hashCollisions = 0;
  u32bit   realCollisions = 0;
  u32bit   maxDup         = 0;

  for (u32bit elem=1; elem<lastElem; elem++) {
    if (fh[elem-1].hash == fh[elem].hash) {
      hashCollisions++;
    }
  }

  fprintf(stderr, "Found "u32bitFMT" hash collisions, examining.\n", hashCollisions);

  for (u32bit elem=1; elem<lastElem; elem++) {
    if (fh[elem-1].hash == fh[elem].hash) {

      //  Grab those two fragments, compare sequence and quality directly

      getFrag(gkp, fh[elem-1].iid, fr1, FRAG_S_SEQ);
      seq1 = getFragRecordSequence(fr1);
      qlt1 = getFragRecordQuality(fr1);

      getFrag(gkp, fh[elem].iid, fr2, FRAG_S_SEQ);
      seq2 = getFragRecordSequence(fr2);
      qlt2 = getFragRecordQuality(fr2);

      if ((strcmp(seq1, seq2) == 0) && (strcmp(qlt1, qlt2) == 0)) {
        realCollisions++;

        fh[elem-1].dupcount++;
        fh[elem].dupcount++;

        if (maxDup < fh[elem-1].dupcount)  maxDup = fh[elem-1].dupcount;
        if (maxDup < fh[elem].dupcount)  maxDup = fh[elem].dupcount;

        //fprintf(stderr, "Dup "u32bitFMT" <-> "u32bitFMT" ("u32bitFMT" hash collisions, "u32bitFMT" real collisions.)\n",
        //        fh[elem-1].iid, fh[elem].iid, hashCollisions, realCollisions);

        uint64 uid1=0, uid2=0;
        uid1 = getFragRecordUID(fr1);
        uid2 = getFragRecordUID(fr2);

        fprintf(stdout, u64bitFMT","u64bitFMT"\n", uid1, uid2);
      }
    }
  }

  fprintf(stderr, "Found "u32bitFMT" real collisions (maximum duplication "u32bitFMT").\n",
          realCollisions, maxDup);

  ////////////////////////////////////////

  fprintf(stderr, "Examine collisions to remove duplicates.\n");

  for (u32bit elem=1; elem<lastElem; elem++) {
    if (fh[elem-1].hash == fh[elem].hash) {

      if ((fh[elem-1].dupcount == 1) && (fh[elem].dupcount == 1)) {
        //  Remove the one without a mate, or pick anyone.

        
      }

    }
  }


  closeGateKeeperStore(gkp);
}



#include "util++.H"
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
  uint64    hash;
  uint32    iid;
  uint32    mate;
  uint32    dupcount:8;
  uint32    killme:1;
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

  uint32   firstElem = getFirstElemFragStore(gkp);
  uint32   lastElem  = getLastElemFragStore(gkp) + 1;

  fragRecord       *fr1 = new_fragRecord();
  fragRecord       *fr2 = new_fragRecord();

  ////////////////////////////////////////

  fragHash   *fh = new fragHash [lastElem - firstElem + 1];
  uint32      seqMax = 10240;
  char       *seq1   = NULL;
  char       *qlt1   = NULL;
  char       *seq2   = NULL;
  char       *qlt2   = NULL;

  ////////////////////////////////////////

  fprintf(stderr, "Read "F_U32" fragments to build hashes.\n", lastElem - firstElem + 1);

  for (uint32 elem=firstElem; elem<lastElem; elem++) {
    getFrag(gkp, elem, fr1, FRAG_S_SEQ);
    seq1 = getFragRecordSequence(fr1);

    uint32 seqLen   = getFragRecordSequenceLength(fr1);
    uint64 hash     = 0;
    uint32 map[256] = { 0 };

    for (uint32 s=0; s<256; s++)
      map[s] = 1;
    map['A'] = map['a'] = 2;
    map['C'] = map['c'] = 3;
    map['G'] = map['g'] = 4;
    map['T'] = map['t'] = 5;

    for (uint64 s=0; s<seqLen; s++) {
      hash  ^= s * (uint64)(map[seq1[s]]) * (uint64)(qlt1[s] - '0');
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

  uint32   reallyDup;
  uint32   hashCollisions = 0;
  uint32   realCollisions = 0;
  uint32   maxDup         = 0;

  for (uint32 elem=1; elem<lastElem; elem++) {
    if (fh[elem-1].hash == fh[elem].hash) {
      hashCollisions++;
    }
  }

  fprintf(stderr, "Found "F_U32" hash collisions, examining.\n", hashCollisions);

  for (uint32 elem=1; elem<lastElem; elem++) {
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

        //fprintf(stderr, "Dup "F_U32" <-> "F_U32" ("F_U32" hash collisions, "F_U32" real collisions.)\n",
        //        fh[elem-1].iid, fh[elem].iid, hashCollisions, realCollisions);

        uint64 uid1=0, uid2=0;
        uid1 = getFragRecordUID(fr1);
        uid2 = getFragRecordUID(fr2);

        fprintf(stdout, F_U64","F_U64"\n", uid1, uid2);
      }
    }
  }

  fprintf(stderr, "Found "F_U32" real collisions (maximum duplication "F_U32").\n",
          realCollisions, maxDup);

  ////////////////////////////////////////

  fprintf(stderr, "Examine collisions to remove duplicates.\n");

  for (uint32 elem=1; elem<lastElem; elem++) {
    if (fh[elem-1].hash == fh[elem].hash) {

      if ((fh[elem-1].dupcount == 1) && (fh[elem].dupcount == 1)) {
        //  Remove the one without a mate, or pick anyone.

        
      }

    }
  }


  closeGateKeeperStore(gkp);
}



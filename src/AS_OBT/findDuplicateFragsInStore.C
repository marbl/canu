#include "trim.H"

//  Read a fragStore, reports duplicate fragments, removes all but one
//  copy of the duplicate.

void
usage(char *name) {
  fprintf(stderr, "usage: %s -frg some.frgStore\n", name);
  fprintf(stderr, "\n");
  fprintf(stderr, "  -frg       Operate on this frgStore\n");
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
  char   *frgStore   = 0L;
  FILE   *logFile    = 0L;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-frg", 2) == 0) {
      frgStore = argv[++arg];
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

  if (!frgStore) {
    usage(argv[0]);
    exit(1);
  }

  //  Open the store
  //
  FragStoreHandle   fs = openFragStore(frgStore, "r");
  if (fs == NULLSTOREHANDLE) {
    fprintf(stderr, "Failed to open %s\n", frgStore);
    exit(1);
  }

  u32bit   firstElem = getFirstElemFragStore(fs);
  u32bit   lastElem  = getLastElemFragStore(fs) + 1;

  ReadStructp       rd1 = new_ReadStruct();
  ReadStructp       rd2 = new_ReadStruct();

  ////////////////////////////////////////

  fragHash   *fh = new fragHash [lastElem - firstElem + 1];
  u32bit      seqMax = 10240;
  char       *seq1   = new char   [seqMax];
  char       *qlt1   = new char   [seqMax];
  char       *seq2   = new char   [seqMax];
  char       *qlt2   = new char   [seqMax];

  ////////////////////////////////////////

  fprintf(stderr, "Read "u32bitFMT" fragments to build hashes.\n", lastElem - firstElem + 1);

  for (u32bit elem=firstElem; elem<lastElem; elem++) {
    getFragStore(fs, elem, FRAG_S_ALL, rd1);
    if (getSequence_ReadStruct(rd1, seq1, qlt1, seqMax)) {
      fprintf(stderr, "getSequence_ReadStruct() failed.\n");
      exit(1);
    }

    u32bit seqLen   = strlen(seq1);
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

      getFragStore(fs, fh[elem-1].iid, FRAG_S_ALL, rd1);
      if (getSequence_ReadStruct(rd1, seq1, qlt1, seqMax)) {
        fprintf(stderr, "getSequence_ReadStruct() failed.\n");
        exit(1);
      }

      getFragStore(fs, fh[elem].iid, FRAG_S_ALL, rd2);
      if (getSequence_ReadStruct(rd2, seq2, qlt2, seqMax)) {
        fprintf(stderr, "getSequence_ReadStruct() failed.\n");
        exit(1);
      }

      if ((strcmp(seq1, seq2) == 0) && (strcmp(qlt1, qlt2) == 0)) {
        realCollisions++;

        fh[elem-1].dupcount++;
        fh[elem].dupcount++;

        if (maxDup < fh[elem-1].dupcount)  maxDup = fh[elem-1].dupcount;
        if (maxDup < fh[elem].dupcount)  maxDup = fh[elem].dupcount;

        //fprintf(stderr, "Dup "u32bitFMT" <-> "u32bitFMT" ("u32bitFMT" hash collisions, "u32bitFMT" real collisions.)\n",
        //        fh[elem-1].iid, fh[elem].iid, hashCollisions, realCollisions);

        uint64 uid1=0, uid2=0;
        getAccID_ReadStruct(rd1, &uid1);
        getAccID_ReadStruct(rd2, &uid2);

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


  closeFragStore(fs);
}



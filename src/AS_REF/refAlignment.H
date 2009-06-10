#ifndef REFALIGNMENT_H
#define REFALIGNMENT_H



static
AS_UID
lookupUID(char *uidstr) {
  AS_UID uid = AS_UID_lookup(uidstr, 0L);

  if (AS_UID_compare(uid, AS_UID_undefined()) == 0)
    uid = AS_UID_load(uidstr);

  return(uid);
}


//  Indexed by IID
struct fragmentData {
  AS_UID   fragUID;
  AS_IID   mateIID;
  AS_IID   libraryIID;

  uint32   mappingDataIndex;
};


struct mappingData {
  AS_UID frgUID;  //  redundant
  AS_IID frgIID;
  uint32 frgBgn;
  uint32 frgEnd;

  uint32 fwd;

  AS_UID refUID;
  uint32 refBgn;
  uint32 refEnd;
};



static
int
mappingDataSortByRefPosition(const void *a, const void *b) {
  const mappingData *A = (const mappingData *)a;
  const mappingData *B = (const mappingData *)b;

  int x = AS_UID_compare(A->refUID, B->refUID);

  if (x)
    return(x);

  if (A->refBgn < B->refBgn)
    return(-1);
  if (A->refBgn > B->refBgn)
    return(1);

  if (A->refEnd < B->refEnd)
    return(-1);
  if (A->refEnd > B->refEnd)
    return(1);

  return(0);
}


static
void
mappingDataPrint(mappingData *md) {
  fprintf(stderr, "%s,%d "F_U32"-"F_U32" %c %s "F_U32"-"F_U32"\n",
          AS_UID_toString(md->frgUID), md->frgIID, md->frgBgn, md->frgEnd,
          (md->fwd) ? 'f' : 'r',
          AS_UID_toString(md->refUID), md->refBgn, md->refEnd);
}


static
void
mappingDataParseLine(mappingData *md, char *line) {
  splitToWords   S(line);

  //  posmap: 5 words -- readID refID refBgn refEnd ori
  //  refmap: 7 words -- readID readBgn readEnd ori refId refBgn refEnd

  if (S.numWords() == 5) {
    md->frgUID = lookupUID(S[0]);
    md->frgIID = gkpStore->gkStore_getUIDtoIID(md->frgUID, NULL);
    md->frgBgn = 0;
    md->frgEnd = 0;

    md->fwd = (S[4][0] == 'f') ? 1 : 0;

    md->refUID = lookupUID(S[1]);
    md->refBgn = atoi(S[2]);
    md->refEnd = atoi(S[3]);

  } else if (S.numWords() == 7) {
    md->frgUID = lookupUID(S[0]);
    md->frgIID = gkpStore->gkStore_getUIDtoIID(md->frgUID, NULL);
    md->frgBgn = atoi(S[1]);
    md->frgEnd = atoi(S[2]);

    md->fwd = (S[3][0] == 'f') ? 1 : 0;

    md->refUID = lookupUID(S[4]);
    md->refBgn = atoi(S[5]);
    md->refEnd = atoi(S[6]);

  } else {
    fprintf(stderr, "Unknown line format for line '%s'\n", line);
    exit(1);
  }
}




#endif



#include "findErrors.H"



//  Open and read fragments with IIDs from  Lo_Frag_IID  to
//  Hi_Frag_IID (INCLUSIVE) from  gkpStore_Path  and store them in
//  global  Frag .

//  This shares lots of code with Extract_Needed_Frags.

void
Read_Frags(feParameters   *G,
           gkStore        *gkpStore) {

  //  The original converted to lowercase, and made non-acgt be 'a'.

  char  filter[256];

  for (uint32 i=0; i<256; i++)
    filter[i] = 'a';

  filter['A'] = filter['a'] = 'a';
  filter['C'] = filter['c'] = 'c';
  filter['G'] = filter['g'] = 'g';
  filter['T'] = filter['t'] = 't';

  //  Count the number of bases, so we can do two gigantic allocations for
  //  bases and votes.

  uint64  basesLength = 0;
  uint64  votesLength = 0;
  uint64  readsLoaded = 0;

  fprintf(stderr, "Read_Frags()-- from "F_U32" through "F_U32"\n",
          G->bgnID, G->endID);

  for (uint32 curID=G->bgnID; curID<=G->endID; curID++) {
    gkRead *read = gkpStore->gkStore_getRead(curID);

    basesLength += read->gkRead_sequenceLength() + 1;
    votesLength += read->gkRead_sequenceLength();
  }

  G->readBases = new char          [basesLength];
  G->readVotes = new Vote_Tally_t  [votesLength];             //  NO constructor, MUST INIT
  G->readsLen  = G->endID - G->bgnID + 1;
  G->reads     = new Frag_Info_t   [G->readsLen];             //  Has constructor, no need to init

  memset(G->readBases, 0, sizeof(char)         * basesLength);
  memset(G->readVotes, 0, sizeof(Vote_Tally_t) * votesLength);

  basesLength = 0;
  votesLength = 0;

  gkReadData  readData;

  for (uint32 curID=G->bgnID; curID<=G->endID; curID++) {
    gkRead *read       = gkpStore->gkStore_getRead(curID);

    if (read->gkRead_isDeleted() == true)
      continue;

    gkpStore->gkStore_loadReadData(read, &readData);

    uint32  readLength = read->gkRead_sequenceLength();
    char   *readBases  = readData.gkReadData_getSequence();

    G->reads[curID - G->bgnID].sequence = G->readBases + basesLength;
    G->reads[curID - G->bgnID].vote     = G->readVotes + votesLength;

    basesLength += readLength + 1;
    votesLength += readLength;
    readsLoaded += 1;

    for (uint32 bb=0; bb<readLength; bb++)
      G->reads[curID - G->bgnID].sequence[bb] = filter[readBases[bb]];

    G->reads[curID - G->bgnID].sequence[readLength] = 0;  //  All good reads end.

    G->reads[curID - G->bgnID].clear_len    = readLength;
    G->reads[curID - G->bgnID].shredded     = false;

    G->reads[curID - G->bgnID].left_degree  = 0;
    G->reads[curID - G->bgnID].right_degree = 0;
  }

  fprintf(stderr, "Read_Frags()-- from "F_U32" through "F_U32" -- loaded "F_U64" bases in "F_U64" reads.\n",
          G->bgnID, G->endID-1, basesLength, readsLoaded);
}

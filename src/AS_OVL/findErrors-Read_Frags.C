

#include "findErrors.H"



//  Open and read fragments with IIDs from  Lo_Frag_IID  to
//  Hi_Frag_IID  from  gkpStore_Path  and store them in
//  global  Frag .

//  This shares lots of code with Extract_Needed_Frags.

void
Read_Frags(feParameters   &G,
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

  for (uint32 curID=G.bgnID; curID<G.endID; curID++) {
    gkRead *read = gkpStore->gkStore_getRead(curID);

    basesLength += read->gkRead_clearRegionLength() + 1;
    votesLength += read->gkRead_clearRegionLength();
  }

  G.readBases = new char          [basesLength];
  G.readVotes = new Vote_Tally_t  [votesLength];
  G.reads     = new Frag_Info_t   [G.endID - G.bgnID + 1];
  G.readsLen  = 0;

  basesLength = 0;
  votesLength = 0;

  gkReadData  readData;

  for (uint32 curID=G.bgnID; curID<G.endID; curID++) {
    gkRead *read       = gkpStore->gkStore_getRead(curID);

    if (read->gkRead_isDeleted() == true)
      continue;

    gkpStore->gkStore_loadReadData(read, &readData);

    uint32  readLength = read->gkRead_clearRegionLength();
    char   *readBases  = readData.gkReadData_getSequence();

    G.reads[G.readsLen].sequence = G.readBases + basesLength;
    G.reads[G.readsLen].vote     = G.readVotes + votesLength;

    basesLength += readLength + 1;
    votesLength += readLength;

    for (uint32 bb=0; bb<readLength; bb++)
      G.reads[G.readsLen].sequence[bb] = filter[readBases[bb]];

    G.reads[G.readsLen].sequence[readLength] = 0;  //  All good reads end.

    G.reads[G.readsLen].left_degree  = 0;
    G.reads[G.readsLen].right_degree = 0;

    G.readsLen++;
  }
}

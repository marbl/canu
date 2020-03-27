
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "findErrors.H"



//  Open and read fragments with IIDs from  Lo_Frag_IID  to
//  Hi_Frag_IID (INCLUSIVE) from  seqStore_Path  and store them in
//  global  Frag .

//  This shares lots of code with Extract_Needed_Frags.

void
Read_Frags(feParameters   *G,
           sqStore        *seqStore) {

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


  for (uint32 curID=G->bgnID; curID<=G->endID; curID++) {
    basesLength += seqStore->sqStore_getReadLength(curID) + 1;
    votesLength += seqStore->sqStore_getReadLength(curID);
  }

  uint64  totAlloc = (sizeof(char)         * basesLength +
                      sizeof(Vote_Tally_t) * votesLength +
                      sizeof(Frag_Info_t)  * G->readsLen);

  fprintf(stderr, "Read_Frags()-- Loading target reads " F_U32 " through " F_U32 " with " F_U64 " bases.\n", G->bgnID, G->endID, basesLength);

  G->readBases = new char          [basesLength];
  G->readVotes = new Vote_Tally_t  [votesLength];             //  Has constructor, no need to init
  G->readsLen  = G->endID - G->bgnID + 1;
  G->reads     = new Frag_Info_t   [G->readsLen];             //  Has constructor, no need to init

  memset(G->readBases, 0, sizeof(char)         * basesLength);

  basesLength = 0;
  votesLength = 0;

  sqRead  *read = new sqRead;

  for (uint32 curID=G->bgnID; curID<=G->endID; curID++) {
    seqStore->sqStore_getRead(curID, read);

    uint32  readLength = read->sqRead_length();

    G->reads[curID - G->bgnID].sequence = G->readBases + basesLength;
    G->reads[curID - G->bgnID].vote     = G->readVotes + votesLength;

    basesLength += readLength + 1;
    votesLength += readLength;
    readsLoaded += 1;

    if (readLength > 0) {
      char   *readBases  = read->sqRead_sequence();

      for (uint32 bb=0; bb<readLength; bb++)
        G->reads[curID - G->bgnID].sequence[bb] = filter[readBases[bb]];
    }

    G->reads[curID - G->bgnID].sequence[readLength] = 0;  //  All good reads end.

    G->reads[curID - G->bgnID].clear_len    = readLength;
    G->reads[curID - G->bgnID].shredded     = false;

    G->reads[curID - G->bgnID].left_degree  = 0;
    G->reads[curID - G->bgnID].right_degree = 0;
  }

  delete read;

  fprintf(stderr, "Read_Frags()-- %.3f GB for bases, votes and info.\n", totAlloc / 1024.0 / 1024.0 / 1024.0);
  fprintf(stderr, "\n");
}

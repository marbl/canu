



#include "correctOverlaps.H"


//  Load overlaps with aIID from G->bgnID to G->endID.
//  Overlaps can be unsorted.

void
Read_Olaps(coParameters *G, gkStore *gkpStore) {

  ovStore *ovs = new ovStore(G->ovlStorePath, gkpStore);

  ovs->setRange(G->bgnID, G->endID);
  
  uint64 numolaps  = ovs->numOverlapsInRange();
  uint64 numNormal = 0;
  uint64 numInnie  = 0;

  fprintf(stderr, "Read_Olaps()--  Loading "F_U64" overlaps from '%s' for reads "F_U32" to "F_U32"\n",
          numolaps, G->ovlStorePath, G->bgnID, G->endID);

  G->olaps    = new Olap_Info_t [numolaps];
  G->olapsLen = 0;

  ovOverlap  olap;

  while (ovs->readOverlap(&olap)) {
    G->olaps[G->olapsLen].a_iid  =  olap.a_iid;
    G->olaps[G->olapsLen].b_iid  =  olap.b_iid;
    G->olaps[G->olapsLen].a_hang =  olap.a_hang();
    G->olaps[G->olapsLen].b_hang =  olap.b_hang();
    //G->olaps[G->olapsLen].orient = (olap.flipped()) ? INNIE : NORMAL;
    G->olaps[G->olapsLen].innie  = (olap.flipped() == true);
    G->olaps[G->olapsLen].normal = (olap.flipped() == false);

    G->olaps[G->olapsLen].order  = G->olapsLen;
    G->olaps[G->olapsLen].evalue = olap.evalue();

    numNormal += (G->olaps[G->olapsLen].normal == true);
    numInnie  += (G->olaps[G->olapsLen].innie  == true);

    G->olapsLen++;
  }

  delete ovs;

  fprintf(stderr, "Read_Olaps()--  Loaded "F_U64" overlaps -- "F_U64" normal and "F_U64" innie.\n",
          G->olapsLen, numNormal, numInnie);
}



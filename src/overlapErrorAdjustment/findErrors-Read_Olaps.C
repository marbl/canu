



#include "findErrors.H"


//  Load overlaps with aIID from G->bgnID to G->endID.
//  Overlaps can be unsorted.

void
Read_Olaps(feParameters *G, gkStore *gkpStore) {
  ovStore *ovs = new ovStore(G->ovlStorePath, gkpStore);

  ovs->setRange(G->bgnID, G->endID);

  uint64 numolaps = ovs->numOverlapsInRange();

  fprintf(stderr, "Read_Olaps()-- loading "F_U64" overlaps.\n",
          numolaps);

  G->olaps    = new Olap_Info_t [numolaps];
  G->olapsLen = 0;

  ovOverlap  olap(gkpStore);

  while (ovs->readOverlap(&olap)) {
    G->olaps[G->olapsLen].a_iid  =  olap.a_iid;
    G->olaps[G->olapsLen].b_iid  =  olap.b_iid;
    G->olaps[G->olapsLen].a_hang =  olap.a_hang();
    G->olaps[G->olapsLen].b_hang =  olap.b_hang();
    G->olaps[G->olapsLen].innie  = (olap.flipped() == true);
    G->olaps[G->olapsLen].normal = (olap.flipped() == false);

    //  These are violated if the innie/normal members are signed!
    assert(G->olaps[G->olapsLen].innie != G->olaps[G->olapsLen].normal);
    assert((G->olaps[G->olapsLen].innie == false) ||
           (G->olaps[G->olapsLen].innie == true));
    assert((G->olaps[G->olapsLen].normal == false) ||
           (G->olaps[G->olapsLen].normal == true));

    G->olapsLen++;
  }

  delete ovs;
}



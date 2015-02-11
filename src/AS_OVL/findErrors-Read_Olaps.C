



#include "findErrors.H"


//  Load overlaps with aIID from G.bgnID to G.endID.
//  Overlaps can be unsorted.

void
Read_Olaps(feParameters &G) {
  ovStore *ovs = new ovStore(G.ovlStorePath);

  ovs->setRange(G.bgnID, G.endID);
  
  uint64 numolaps = ovs->numOverlapsInRange();

  G.olaps    = new Olap_Info_t [numolaps];
  G.olapsLen = 0;

  ovsOverlap  olap;

  while (ovs->readOverlap(&olap)) {
    G.olaps[G.olapsLen].a_iid  =  olap.a_iid;
    G.olaps[G.olapsLen].b_iid  =  olap.b_iid;
    G.olaps[G.olapsLen].a_hang =  olap.a_hang();
    G.olaps[G.olapsLen].b_hang =  olap.b_hang();
    //G.olaps[G.olapsLen].orient = (olap.flipped()) ? INNIE : NORMAL;
    G.olaps[G.olapsLen].innie  = (olap.flipped() == true);
    G.olaps[G.olapsLen].normal = (olap.flipped() == false);

    G.olapsLen++;
  }

  delete ovs;
}



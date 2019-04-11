
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-JUN-16 to 2015-JUN-25
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-OCT-17
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "findErrors.H"


//  Load overlaps with aIID from G->bgnID to G->endID.
//  Overlaps can be unsorted.

void
Read_Olaps(feParameters *G, sqStore *seqStore) {
  ovStore *ovs = new ovStore(G->ovlStorePath, seqStore);

  ovs->setRange(G->bgnID, G->endID);

  uint64 numolaps = ovs->numOverlapsInRange();

  fprintf(stderr, "Read_Olaps()-- Loading " F_U64 " overlaps.\n", numolaps);

  G->olaps    = new Olap_Info_t [numolaps];
  G->olapsLen = 0;

  ovOverlap  olap;

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

  fprintf(stderr, "Read_Olaps()-- %.3f GB for overlaps..\n", sizeof(Olap_Info_t) * numolaps / 1024.0 / 1024.0 / 1024.0);
  fprintf(stderr, "\n");
}



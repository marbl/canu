
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
 *    Brian P. Walenz beginning on 2016-MAY-02
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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

  fprintf(stderr, "--Allocate "F_U64" MB for overlaps.\n",
          (sizeof(Olap_Info_t) * numolaps) >> 20);

  G->olaps    = new Olap_Info_t [numolaps];
  G->olapsLen = 0;

  ovOverlap  olap(gkpStore);

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



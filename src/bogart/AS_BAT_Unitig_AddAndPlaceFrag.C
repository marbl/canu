
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
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_Unitig_AddAndPlaceFrag.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010,2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"



//  Given two edges, place fragment node.ident into this unitig using the thickest edge to decide on
//  the placement.  At least one of the edges must be from the node to a fragment in the target
//  unitig.
//
//  Returns true if placement was successful.
//
bool
Unitig::addAndPlaceFrag(int32 fid, BestEdgeOverlap *bestedge5, BestEdgeOverlap *bestedge3, bool report) {
  int32        bidx5 = -1,   bidx3 = -1;
  int32        blen5 =  0,   blen3 =  0;
  ufNode       frag;

  frag.ident             = fid;
  frag.contained         = 0;
  frag.parent            = 0;
  frag.ahang             = 0;
  frag.bhang             = 0;
  frag.position.bgn      = 0;
  frag.position.end      = 0;
  frag.containment_depth = 0;

  //  The length of the overlap depends only on the length of the a frag and the hangs.  We don't
  //  actually care about the real length (except for logging), only which is thicker.

  if ((bestedge5) && (bestedge5->fragId() == 0))
    bestedge5 = NULL;

  if ((bestedge3) && (bestedge3->fragId() == 0))
    bestedge3 = NULL;

  if ((bestedge5) && (fragIn(bestedge5->fragId()) == id())) {
    bidx5 = pathPosition(bestedge5->fragId());
    blen5 = FI->fragmentLength(fid) + ((bestedge5->ahang() < 0) ? bestedge5->bhang() : -bestedge5->ahang());
#ifdef DEBUG_PLACEMENT
    writeLog("addAndPlaceFrag()-- bestedge5:  %d,%d,%d,%d len %d\n",
            bestedge5->fragId(), bestedge5->frag3p, bestedge5->ahang(), bestedge5->bhang(), blen5);
#endif
    assert(bestedge5->fragId() == ufpath[bidx5].ident);
  }

  if ((bestedge3) && (fragIn(bestedge3->fragId()) == id())) {
    bidx3 = pathPosition(bestedge3->fragId());;
    blen3 = FI->fragmentLength(fid) + ((bestedge3->ahang() < 0) ? bestedge3->bhang() : -bestedge3->ahang());
#ifdef DEBUG_PLACEMENT
    writeLog("addAndPlaceFrag()-- bestedge3:  %d,%d,%d,%d len %d\n",
            bestedge3->fragId(), bestedge3->frag3p, bestedge3->ahang(), bestedge3->bhang(), blen3);
#endif
    assert(bestedge3->fragId() == ufpath[bidx3].ident);
  }

  //  Use the longest that exists -- an alternative would be to take the average position, but that
  //  could get messy if the placements are different.  Picking one or the other has a better chance
  //  of working, though it'll fail if the fragment is chimeric or spans something it shouldn't,
  //  etc.

  if ((blen5 == 0) && (blen3 == 0)) {
    writeLog("Unitig::addAndPlaceFrag()-- WARNING:  Failed to place frag %d into unitig %d; no edges to the unitig.\n",
            fid, id());
    return(false);
  }

  if (blen5 < blen3)
    bestedge5 = NULL;
  else
    bestedge3 = NULL;

  //  Compute the placement -- a little scary, as we stuff both placements into the same frag, but
  //  we guarantee only one placement is computed.

  if (placeFrag(frag, bidx5, bestedge5,
                frag, bidx3, bestedge3) == false)
    return(false);

  //  If we just computed a placement before the start of the unitig, we need to shift the unitig to
  //  make space.

  int32 frgBgn = MIN(frag.position.bgn, frag.position.end);

  if (frgBgn < 0) {
    frgBgn = -frgBgn;

    frag.position.bgn += frgBgn;
    frag.position.end += frgBgn;

    _length += frgBgn;

    for (uint32 fi=0; fi<ufpath.size(); fi++) {
      ufNode *tfrg = &ufpath[fi];

      tfrg->position.bgn += frgBgn;
      tfrg->position.end += frgBgn;
    }
  }

  //  Finally, add the fragment.

  addFrag(frag, 0, report);

  return(true);
}


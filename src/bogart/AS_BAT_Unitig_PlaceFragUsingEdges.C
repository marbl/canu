
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
 *    src/AS_BAT/AS_BAT_Unitig_PlaceFragUsingEdges.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2015-MAR-06
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"


#undef  DEBUG_PLACE_FRAG



ufNode
placeFrag_computePlacement(uint32           fragId,
                           bool             frag3p,
                           ufNode          &parent,
                           BestEdgeOverlap *edge) {

  //  We have an 'edge' from 'fragId' end 'frag3p' back to 'parent'.
  //  Use that to compute the placement of 'frag'.

  bool   pFwd = (parent.position.bgn < parent.position.end) ? true : false;
  int32  pMin = (parent.position.bgn < parent.position.end) ? parent.position.bgn : parent.position.end;
  int32  pMax = (parent.position.bgn < parent.position.end) ? parent.position.end : parent.position.bgn;

  assert(pMin < pMax);

  //  Scale the hangs based on the placed versus actual length of the parent read.

  //double  intraScale = (double)(pMax - pMin) / FI->fragmentLength(parent.ident);  //  Within the parent read overlap
  //double  interScale = 1.0;                                                       //  Outside the parent read overlap

  //  We have two cases, adding a read 'to the left' or 'to the right' of the parent read.
  //  'To the right' looks like:
  //
  //        -a    --/-------    frag (edge from frag to parent)
  //        -------v---   -b    parent
  //
  //  'To the left' looks like:
  //
  //      ------\---  +b        frag (edge from frag to parent)
  //      +a  ---v-------       parent
  //
  //  The construction of our edges is such that they tell, via simple addition, the placement
  //  of the next read.  So, in both cases, to flip the edge to be from parent to frag, all we
  //  need to do is negate both hangs.
  //
  //
  //  HOWEVER, if the read we're placing is reverse (edge is from the 3p end), we still need to swap
  //  the hangs when the edge is reversed, because the ahang/bhang labels are assuming the source
  //  read is forward oriented.
  //
  //      -------\-->  1035   This edge has ahang=831 bhang=1035.
  //      831  <==v========
  //
  //      -------^-->  1035   This edge has ahang=1035 and bhang=831
  //      831  <==\========   (rotate the picture 180 degrees to get the 'standard' overlap form)
  //

  //assert(pMin > edge->ahang());
  //assert(pMax > edge->bhang());

  int32  fMin = pMin + ((frag3p == false) ? -edge->ahang() : edge->bhang());   //  * intraScale
  int32  fMax = pMax + ((frag3p == false) ? -edge->bhang() : edge->ahang());   //  * interScale

  assert(fMin < fMax);

  //  We don't know the true length of the overlap, and our hang-based math tends to shrink reads.
  //  Reset the end coordinate using the actual length of the read.

  fMax = fMin + FI->fragmentLength(fragId);


  //  Orientation is a bit more complicated, with eight cases (drawing pictures helps).
  //
  //    edge from frag3p=true  to forward parent 3p -> reverse
  //    edge from frag3p=false to reverse parent 3p -> reverse
  //    edge from frag3p=false to forward parent 5p -> reverse
  //    edge from frag3p=true  to reverse parent 5p -> reverse
  //
  //    edge from frag3p=true  to reverse parent 3p -> forward
  //    edge from frag3p=false to forward parent 3p -> forward
  //    edge from frag3p=false to reverse parent 5p -> forward
  //    edge from frag3p=true  to forward parent 5p -> forward
  //
  bool fFwd = (((frag3p == true)  && (pFwd == true)  && (edge->frag3p() == true))  ||
               ((frag3p == false) && (pFwd == false) && (edge->frag3p() == true))  ||
               ((frag3p == false) && (pFwd == true)  && (edge->frag3p() == false)) ||
               ((frag3p == true)  && (pFwd == false) && (edge->frag3p() == false))) ? false : true;

  ufNode   frag;

  frag.ident        = fragId;
  frag.contained    = 0;
  frag.parent       = edge->fragId();       //  == parent->ident
  frag.ahang        = 0;                    //  Not used in bogart, set on output
  frag.bhang        = 0;                    //  Not used in bogart, set on output
  frag.position.bgn = (fFwd) ? fMin : fMax;
  frag.position.end = (fFwd) ? fMax : fMin;

#ifdef DEBUG_PLACE_FRAG
  writeLog("placeFrag()-- parent %7d pos %7d,%7d -- edge to %7d %c' hangs %7d %7d -- frag %7d %c' -- placed %7d-%7d oriented %s %7d-%7d\n",
           parent.ident, parent.position.bgn, parent.position.end,
           edge->fragId(), (edge->frag3p()) ? '3' : '5', edge->ahang(), edge->bhang(),
           fragId, (frag3p) ? '3' : '5',
           fMin, fMax, (fFwd) ? "rev" : "fwd", frag.position.bgn, frag.position.end);
#endif

  return(frag);
}





//  Place a read into this tig using an edge from the read to some read in this tig.
//    frag   - output placement
//    fragId - read we're trying to place
//    frag3p - edge is from the 3p end of the read?
//    edge   - edge from read to something in this tig
//
bool
Unitig::placeFrag(ufNode          &frag,
                  uint32           fragId,
                  bool             frag3p,
                  BestEdgeOverlap *edge) {

  assert(fragId > 0);
  assert(fragId <= FI->numFragments());

  frag.ident             = fragId;
  frag.contained         = 0;
  frag.parent            = 0;
  frag.ahang             = 0;
  frag.bhang             = 0;
  frag.position.bgn      = 0;
  frag.position.end      = 0;

  if (edge == NULL)
    //  No best edge?  Hard to place without one.
    return(false);

  if (edge->fragId() == 0)
    //  Empty best edge?  Still hard to place.
    return(false);

  if (fragIn(edge->fragId()) != id())
    //  Edge not pointing to a read in this tig?
    return(false);

  //  Grab the index of the parent read.

  uint32 bidx = pathPosition(edge->fragId());
  assert(edge->fragId() == ufpath[bidx].ident);

  //  Now, just compute the placement and return success!

  frag = placeFrag_computePlacement(fragId, frag3p, ufpath[bidx], edge);

  return(true);
}

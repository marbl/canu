
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

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"


#undef  DEBUG_PLACE_FRAG


ufNode
placeFrag_contained(uint32           fragId,
                    ufNode          &parent,
                    BestEdgeOverlap *edge) {

  bool   pFwd  = (parent.position.bgn < parent.position.end) ? true : false;
  int32  pMin  = (parent.position.bgn < parent.position.end) ? parent.position.bgn : parent.position.end;
  int32  pMax  = (parent.position.bgn < parent.position.end) ? parent.position.end : parent.position.bgn;

  assert(pMin < pMax);

  //  Reverse the overlap.  frag3p here means the overlap is flipped.
  int32  ahang = (edge->frag3p() == false) ? -edge->ahang() : edge->bhang();
  int32  bhang = (edge->frag3p() == false) ? -edge->bhang() : edge->ahang();

  //  Depending on the parent orientation...
  //
  //     pMin         pMax        pMin         pMax
  //     ---------------->        <----------------
  //     ahang ----- bhang        bhang ----- ahang
  //     > 0         < 0          < 0         > 0

  int32 fMin = (pFwd == true) ? pMin + ahang : pMin - bhang;
  int32 fMax = (pFwd == true) ? pMax + bhang : pMax - ahang;

  //int32  fMin = pMin + ((frag3p == false) ? -edge->ahang() : edge->bhang());   //  * intraScale
  //int32  fMax = pMax + ((frag3p == false) ? -edge->bhang() : edge->ahang());   //  * interScale

  assert(fMin < fMax);

  //  We don't know the true length of the overlap, and our hang-based math tends to shrink reads.
  //  Reset the end coordinate using the actual length of the read.

  fMax = fMin + FI->fragmentLength(fragId);

  //  Orientation is straightforward, based on the orient of the parent, and the flipped flag.

  bool fFwd = (((pFwd == true)  && (edge->frag3p() == false)) ||  //  parent is fwd, olap is not flipped
               ((pFwd == false) && (edge->frag3p() == true)));    //  parent is rev, olap is     flipped

  ufNode   frag;

  frag.ident        = fragId;
  frag.contained    = 0;
  frag.parent       = edge->fragId();       //  == parent->ident
  frag.ahang        = 0;                    //  Not used in bogart, set on output
  frag.bhang        = 0;                    //  Not used in bogart, set on output
  frag.position.bgn = (fFwd) ? fMin : fMax;
  frag.position.end = (fFwd) ? fMax : fMin;

#ifdef DEBUG_PLACE_FRAG
  writeLog("placeCont()-- parent %7d pos %7d,%7d -- edge to %7d %c' hangs %7d %7d -- frag %7d C' -- placed %7d-%7d oriented %s %7d-%7d\n",
           parent.ident, parent.position.bgn, parent.position.end,
           edge->fragId(), (edge->frag3p()) ? '3' : '5', edge->ahang(), edge->bhang(),
           fragId,
           fMin, fMax, (fFwd) ? "rev" : "fwd", frag.position.bgn, frag.position.end);
#endif

  return(frag);
}





ufNode
placeFrag_dovetail(uint32           fragId,
                   bool             frag3p,
                   ufNode          &parent,
                   BestEdgeOverlap *edge) {

  //  We have an 'edge' from 'fragId' end 'frag3p' back to 'parent'.
  //  Use that to compute the placement of 'frag'.

  bool   pFwd  = (parent.position.bgn < parent.position.end) ? true : false;
  int32  pMin  = (parent.position.bgn < parent.position.end) ? parent.position.bgn : parent.position.end;
  int32  pMax  = (parent.position.bgn < parent.position.end) ? parent.position.end : parent.position.bgn;

  assert(pMin < pMax);

  //  Scale the hangs based on the placed versus actual length of the parent read.

  //double  intraScale = (double)(pMax - pMin) / FI->fragmentLength(parent.ident);  //  Within the parent read overlap
  //double  interScale = 1.0;                                                       //  Outside the parent read overlap

  //  We're given an edge from the read-to-place back to the parent.  Reverse the edge so it points
  //  from the parent to the read-to-place.
  //
  //  The canonical edge is from a forward parent to the child.
  //
  //      -P----\-->      +b
  //      +a  ---v--------C-
  //
  //  To reverse the edge:
  //
  //    If child is forward, swapping the order of the reads results in a canonical overlap.  The
  //    hangs become negative.
  //
  //      -P----\-->      +b    ---->   -a  ---/--------C>
  //      +a  ---v--------C>    ---->   -P----v-->      -b
  //
  //    If child is reverse, swapping the order of the reads results in a backwards canonical
  //    overlap, and we need to flip end-to-end also.  The hangs are swapped.
  //
  //      -P----\-->      +b    ---->   -C--------\-->  +a
  //      +a  <--v--------C-    ---->   +b      <--v----P-
  //
  int32  ahang = (frag3p == false) ? -edge->ahang() : edge->bhang();
  int32  bhang = (frag3p == false) ? -edge->bhang() : edge->ahang();

  //  The read is placed 'to the right' of the parent if
  //    pFwd == true  and edge points to 3' end
  //    pFwd == false and edge points to 5' end
  //
  bool  toRight = (pFwd == edge->frag3p());

  //  If placing 'to the right', we add hangs.  Else, subtract the swapped hangs.

  int32  fMin = 0;
  int32  fMax = 0;

  if (toRight) {
    fMin = pMin + ahang;
    fMax = pMax + bhang;
  } else {
    fMin = pMin - bhang;
    fMax = pMax - ahang;
  }

  //int32  fMin = pMin + ((frag3p == false) ? -edge->ahang() : edge->bhang());   //  * intraScale
  //int32  fMax = pMax + ((frag3p == false) ? -edge->bhang() : edge->ahang());   //  * interScale

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
  writeLog("placeDove()-- parent %7d pos %7d,%7d -- edge to %7d %c' hangs %7d %7d -- frag %7d %c' -- placed %7d-%7d oriented %s %7d-%7d\n",
           parent.ident, parent.position.bgn, parent.position.end,
           edge->fragId(), (edge->frag3p()) ? '3' : '5', edge->ahang(), edge->bhang(),
           fragId, (frag3p) ? '3' : '5',
           fMin, fMax, (fFwd) ? "rev" : "fwd", frag.position.bgn, frag.position.end);
#endif

  return(frag);
}





//  Place a read into this tig using an edge from the read to some read in this tig.
//
bool
Unitig::placeFrag(ufNode          &frag,      //  output placement
                  uint32           fragId,    //  id of read we want to place
                  bool             frag3p,    //  end of read 'edge' is from, meaningless if contained
                  BestEdgeOverlap *edge) {    //  edge to read in this tig

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

  if (((edge->ahang() >= 0) && (edge->bhang() <= 0)) ||
      ((edge->ahang() <= 0) && (edge->bhang() >= 0)))
    frag = placeFrag_contained(fragId, ufpath[bidx], edge);
  else
    frag = placeFrag_dovetail(fragId, frag3p, ufpath[bidx], edge);

  return(true);
}

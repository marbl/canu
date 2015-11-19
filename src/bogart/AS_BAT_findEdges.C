
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
 *    src/AS_BAT/AS_BAT_findEdges.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-NOV-14 to 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

//  Given two fragments that share at least one edge, this will find that edge and construct a new
//  edge to make it mutual.
//
//  For example, if there is a best edge from aFrg 3' to bFrg 5', this will return that edge in a3,
//  and also create the symmetric edge in b5.
//
static
bool
findEdges(ufNode *aFrg, BestEdgeOverlap &a5, BestEdgeOverlap &a3,
                       ufNode *bFrg, BestEdgeOverlap &b5, BestEdgeOverlap &b3) {

  if (OG->isContained(aFrg->ident) ||
      OG->isContained(bFrg->ident))
    return(false);

  //  Grab what edges we have.

  a5 = *OG->getBestEdgeOverlap(aFrg->ident, false);
  a3 = *OG->getBestEdgeOverlap(aFrg->ident, true);
  b5 = *OG->getBestEdgeOverlap(bFrg->ident, false);
  b3 = *OG->getBestEdgeOverlap(bFrg->ident, true);

  //  Erase things that aren't correct

  if (a5.fragId() != bFrg->ident)  a5 = BestEdgeOverlap();
  if (a3.fragId() != bFrg->ident)  a3 = BestEdgeOverlap();
  if (b5.fragId() != aFrg->ident)  b5 = BestEdgeOverlap();
  if (b3.fragId() != aFrg->ident)  b3 = BestEdgeOverlap();

  //  If we have no edges left, there are no edges!

  if ((b5.fragId() != aFrg->ident) && (b3.fragId() != aFrg->ident) &&
      (a5.fragId() != bFrg->ident) && (a3.fragId() != bFrg->ident))
    return(false);

  //  If we found TWO edges for any single fragment....that's madness!  That means the fragment
  //  had best dovetail overlaps to the same other fragment off of both ends.  We'll complain
  //  and return failure.  Ideally, data like this will be cleaned up by OBT, or filtered from
  //  our input.
  //
  if (a5.fragId() == a3.fragId()) {
    writeLog("findEdges()-- frag %d has multiple edges to frag %d - a5 %d/%d' a3 %d/%d'\n",
            aFrg->ident, a5.fragId(),
            a5.fragId(), a5.frag3p() ? 3 : 5,
            a5.fragId(), a5.frag3p() ? 3 : 5);
  }

  if (b5.fragId() == b3.fragId()) {
    writeLog("findEdges()-- frag %d has multiple edges to frag %d - b5 %d/%d' b3 %d/%d'\n",
            bFrg->ident, b5.fragId(),
            b5.fragId(), b5.frag3p() ? 3 : 5,
            b5.fragId(), b5.frag3p() ? 3 : 5);
  }

  if (((a5.fragId() != 0) && (a5.fragId() == a3.fragId())) ||
      ((b5.fragId() != 0) && (b5.fragId() == b3.fragId()))) {
    a5 = BestEdgeOverlap();
    a3 = BestEdgeOverlap();
    b5 = BestEdgeOverlap();
    b3 = BestEdgeOverlap();
    return(false);
  }

  //  Now, populate the other edges using whatever we have.  Best case is that we have two edges
  //  (because we're done).

  assert(((a5.fragId() == bFrg->ident) +
          (a3.fragId() == bFrg->ident) +
          (b5.fragId() == aFrg->ident) +
          (b3.fragId() == aFrg->ident)) <= 2);

  if (((a5.fragId() == bFrg->ident) || (a3.fragId() == bFrg->ident)) &&
      ((b5.fragId() == aFrg->ident) || (b3.fragId() == aFrg->ident)))
    return(true);

  //  Otherwise, we have exactly one edge, and the other one needs to be created.

  assert(((a5.fragId() == bFrg->ident) +
          (a3.fragId() == bFrg->ident) +
          (b5.fragId() == aFrg->ident) +
          (b3.fragId() == aFrg->ident)) == 1);

  if        (a5.fragId() == bFrg->ident) {
    //assert(a5.fragId() == 0);
    assert(a3.fragId() == 0);
    assert(b5.fragId() == 0);
    assert(b3.fragId() == 0);

    //  Edge off of A's 5' end ('false' below)...
    //  ...to B's 3' end (so ANTI or NORMAL -- negate the hangs)
    //  ...to B's 5' end (so INNIE or OUTTIE -- swap the hangs)
    if (a5.frag3p())
      b3.set(aFrg->ident, false, -a5.ahang(), -a5.bhang());
    else
      b5.set(aFrg->ident, false, a5.bhang(), a5.ahang());

  } else if (a3.fragId() == bFrg->ident) {
    assert(a5.fragId() == 0);
    //assert(a3.fragId() == 0);
    assert(b5.fragId() == 0);
    assert(b3.fragId() == 0);

    //  Edge off of A's 3' end ('true' below)...
    //  ...to B's 3' end (so INNIE or OUTTIE -- swap the hangs)
    //  ...to B's 5' end (so ANTI or NORMAL -- negate the hangs)
    if (a3.frag3p())
      b3.set(aFrg->ident, true, a3.bhang(), a3.ahang());
    else
      b5.set(aFrg->ident, true, -a3.ahang(), -a3.bhang());

  } else if (b5.fragId() == aFrg->ident) {
    assert(a5.fragId() == 0);
    assert(a3.fragId() == 0);
    //assert(b5.fragId() == 0);
    assert(b3.fragId() == 0);

    if (b5.frag3p())
      a3.set(bFrg->ident, false, -b5.ahang(), -b5.bhang());
    else
      a5.set(bFrg->ident, false, b5.bhang(), b5.ahang());


  } else if (b3.fragId() == aFrg->ident) {
    assert(a5.fragId() == 0);
    assert(a3.fragId() == 0);
    assert(b5.fragId() == 0);
    //assert(b3.fragId() == 0);

    if (b3.frag3p())
      a3.set(bFrg->ident, true, b3.bhang(), b3.ahang());
    else
      a5.set(bFrg->ident, true, -b3.ahang(), -b3.bhang());

  } else {
    fprintf(stderr, "findEdges()-- Logically impossible!\n");
    assert(0);
  }

  //  And now we should have exactly two edges.

  assert(((a5.fragId() == bFrg->ident) +
          (a3.fragId() == bFrg->ident) +
          (b5.fragId() == aFrg->ident) +
          (b3.fragId() == aFrg->ident)) == 2);

  return(true);
}



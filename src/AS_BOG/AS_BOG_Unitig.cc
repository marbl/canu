
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static const char *rcsid = "$Id: AS_BOG_Unitig.cc,v 1.24 2010-04-20 16:06:53 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_Unitig.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#undef max


extern FragmentInfo     *debugfi;
extern BestOverlapGraph *bog;

// various class static methods and variables
static std::map<uint32,int>* containPartialOrder;

uint32  Unitig::nextId        = 1;
uint32* Unitig::_inUnitig     = NULL;
uint32* Unitig::_pathPosition = NULL;


Unitig::Unitig(bool report){
  _localArrivalRate = -1;
  _covStat          = FLT_MAX;
  _length           = -1;
  _avgRho           = -1;
  dovetail_path_ptr = new DoveTailPath;
  _id               = nextId++;

  if (report)
    fprintf(stderr, "Creating Unitig %d\n", _id);
}


Unitig::~Unitig(void){
  delete dovetail_path_ptr;
}


#warning WHAT REALLLY HAPPENS IF NO BACKBONE NODE, OR NO PREVIOUS BACKBONE NODE

DoveTailNode Unitig::getLastBackboneNode(void) {
  for (int32 fi=dovetail_path_ptr->size()-1; fi >= 0; fi--) {
    DoveTailNode  &node = (*dovetail_path_ptr)[fi];

    if (node.contained)
      continue;

    return(node);
  }

  fprintf(stderr, "Unitig::getLastBackboneNode()--  WARNING:  unitig %d has no backbone nodes, all contained!\n", id());
  DoveTailNode last;
  memset(&last, 0, sizeof(DoveTailNode));
  return(last);
}


DoveTailNode Unitig::getLastBackboneNode(uint32 &prevID) {
  DoveTailNode last;

  memset(&last, 0, sizeof(DoveTailNode));

  prevID = 0;

  for (int32 fi=dovetail_path_ptr->size()-1; (fi >= 0) && (prevID == 0); fi--) {
    DoveTailNode  *node = &(*dovetail_path_ptr)[fi];

    if (node->contained)
      continue;

    if (last.ident == 0)
      //  Save the last dovetail node, but keep looking....
      last = *node;
    else
      //  ...for the next to last ID.
      prevID = node->ident;
  }

  return(last);
}




//  Given an implicit fragment, and at least one best edge to a fragment in this unitig, compute the
//  position of the fragment in this unitig.  If both edges are given, both will independently
//  compute a placement, which might disagree.
//
//  If a placement is not found for an edge, the corresponding bidx value is set to -1.  Otherwise,
//  it is set to the position in the fragment list of the fragment in this unitig (from above).
//
//  Returns true if any placement is found, false otherwise.
//  The bidx value is set to -1 if no placement is found for that end.
//
bool
Unitig::placeFrag(DoveTailNode &frag5, int32 &bidx5, BestEdgeOverlap *bestedge5,
                  DoveTailNode &frag3, int32 &bidx3, BestEdgeOverlap *bestedge3) {
  bool  verbose = false;

#ifdef WITHIMP
  frag5.type         = AS_READ;
#endif
  //frag5.ident
  frag5.contained    = 0;
  frag5.parent       = 0;
  frag5.ahang        = 0;
  frag5.bhang        = 0;
  frag5.position.bgn = 0;
  frag5.position.end = 0;
#ifdef WITHIMP
  frag5.delta_length = 0;
  frag5.delta        = NULL;
#endif

#ifdef WITHIMP
  frag3.type         = AS_READ;
#endif
  //frag3.ident
  frag3.contained    = 0;
  frag3.parent       = 0;
  frag3.ahang        = 0;
  frag3.bhang        = 0;
  frag3.position.bgn = 0;
  frag3.position.end = 0;
#ifdef WITHIMP
  frag3.delta_length = 0;
  frag3.delta        = NULL;
#endif

  assert(frag3.ident > 0);
  assert(frag5.ident > 0);

  assert(frag3.ident <= debugfi->numFragments());
  assert(frag5.ident <= debugfi->numFragments());

  bidx5              = -1;
  bidx3              = -1;

  //  If we have an incoming edge, AND the fragment for that edge is in this unitig, look up its
  //  index.  Otherwise, discard the edge to prevent placement.
  //
  if ((bestedge5) && (fragIn(bestedge5->frag_b_id) == id())) {
    bidx5 = pathPosition(bestedge5->frag_b_id);
    assert(bestedge5->frag_b_id == (*dovetail_path_ptr)[bidx5].ident);
  } else {
    bestedge5 = NULL;
    bidx5     = -1;
  }

  if ((bestedge3) && (fragIn(bestedge3->frag_b_id) == id())) {
    bidx3 = pathPosition(bestedge3->frag_b_id);;
    assert(bestedge3->frag_b_id == (*dovetail_path_ptr)[bidx3].ident);
  } else {
    bestedge3 = NULL;
    bidx3     = -1;

  }
  //  Now, just compute the placement based on edges that exist.

  if ((bestedge5) && (bidx5 != -1)) {
    DoveTailNode *parent = &(*dovetail_path_ptr)[bidx5];

    assert(parent->ident == bestedge5->frag_b_id);

    //  Overlap is stored using 'node' as the A frag, and we negate the hangs to make them relative
    //  to the 'parent'.  (This is opposite from how containment edges are saved.)  A special case
    //  exists when we overlap to the 5' end of the other fragment; we need to flip the overlap to
    //  ensure the (new) A frag is forward.

    int ahang = -bestedge5->ahang;
    int bhang = -bestedge5->bhang;

    if (bestedge5->bend == FIVE_PRIME) {
      ahang = bestedge5->bhang;
      bhang = bestedge5->ahang;
    }

    int  bgn, end;

    //  Place the new fragment using the overlap.  We don't worry about the orientation of the new
    //  fragment, only the location.  Orientation of the parent fragment matters (1) to know which
    //  coordinate is the lower, and (2) to decide if the overlap needs to be flipped (again).

    if (parent->position.bgn < parent->position.end) {
      bgn = parent->position.bgn + ahang;
      end = parent->position.end + bhang;
    } else {
      bgn = parent->position.end - bhang;
      end = parent->position.bgn - ahang;
    }

    assert(bgn < end);

    //  Since we don't know the true length of the overlap, if we use just the hangs to place a
    //  fragment, we typically shrink fragments well below their actual length.  In one case, we
    //  shrank a container enough that the containee was placed in the unitig backwards.
    //
    //  We now revert back to placing the end based on the actual length, but will
    //  adjust to maintain a dovetail relationship.
    //
    //  See comments on other instances of this warning.
    //

#warning not knowing the overlap length really hurts.
    end = bgn + debugfi->fragmentLength(frag5.ident);
    if (end <= MAX(parent->position.bgn, parent->position.end))
      end = MAX(parent->position.bgn, parent->position.end) + 1;

    //  The new frag is reverse if:
    //    the old frag is forward and we hit its 5' end, or
    //    the old frag is reverse and we hit its 3' end.
    //
    //  The new frag is forward if:
    //    the old frag is forward and we hit its 3' end, or
    //    the old frag is reverse and we hit its 5' end.
    //
    bool flip = (((parent->position.bgn < parent->position.end) && (bestedge5->bend == FIVE_PRIME)) ||
                 ((parent->position.end < parent->position.bgn) && (bestedge5->bend == THREE_PRIME)));

    if (verbose)
      fprintf(stderr, "bestedge5:  parent iid %d pos %d,%d   b_iid %d ovl %d,%d,%d  pos %d,%d  flip %d\n",
              parent->ident, parent->position.bgn, parent->position.end,
              bestedge5->frag_b_id, bestedge5->bend, bestedge5->ahang, bestedge5->bhang, bgn, end, flip);

    frag5.contained    = 0;
    frag5.parent       = bestedge5->frag_b_id;
    frag5.ahang        = ahang;
    frag5.bhang        = bhang;
    frag5.position.bgn = (flip) ? end : bgn;
    frag5.position.end = (flip) ? bgn : end;
  }


  if ((bestedge3) && (bidx3 != -1)) {
    DoveTailNode *parent = &(*dovetail_path_ptr)[bidx3];

    assert(parent->ident == bestedge3->frag_b_id);

    int ahang = -bestedge3->ahang;
    int bhang = -bestedge3->bhang;

    if (bestedge3->bend == THREE_PRIME) {
      ahang = bestedge3->bhang;
      bhang = bestedge3->ahang;
    }

    int  bgn, end;

    if (parent->position.bgn < parent->position.end) {
      bgn = parent->position.bgn + ahang;
      end = parent->position.end + bhang;
    } else {
      bgn = parent->position.end - bhang;
      end = parent->position.bgn - ahang;
    }

    assert(bgn < end);

#warning not knowing the overlap length really hurts.
    end = bgn + debugfi->fragmentLength(frag3.ident);
    if (end <= MAX(parent->position.bgn, parent->position.end))
      end = MAX(parent->position.bgn, parent->position.end) + 1;

    //  The new frag is reverse if:
    //    the old frag is forward and we hit its 3' end, or
    //    the old frag is reverse and we hit its 5' end.
    //
    //  The new frag is forward if:
    //    the old frag is forward and we hit its 5' end, or
    //    the old frag is reverse and we hit its 3' end.
    //
    bool flip = (((parent->position.bgn < parent->position.end) && (bestedge3->bend == THREE_PRIME)) ||
                 ((parent->position.end < parent->position.bgn) && (bestedge3->bend == FIVE_PRIME)));

    if (verbose)
      fprintf(stderr, "bestedge3:  parent iid %d %d,%d   b_iid %d ovl %d,%d,%d  pos %d,%d flip %d\n",
              parent->ident, parent->position.bgn, parent->position.end,
              bestedge3->frag_b_id, bestedge3->bend, bestedge3->ahang, bestedge3->bhang, bgn, end, flip);

    frag3.contained    = 0;
    frag3.parent       = bestedge3->frag_b_id;
    frag3.ahang        = ahang;
    frag3.bhang        = bhang;
    frag3.position.bgn = (flip) ? end : bgn;
    frag3.position.end = (flip) ? bgn : end;
  }

  return((bidx5 != -1) || (bidx3 != -1));
}




void
Unitig::addFrag(DoveTailNode node, int offset, bool report) {

  node.position.bgn += offset;
  node.position.end += offset;

  assert(node.ident > 0);

  // keep track of the unitig a frag is in
  _inUnitig[node.ident]     = _id;
  _pathPosition[node.ident] = dovetail_path_ptr->size();

  // keep track of max position in unitig
  int32 frgEnd = MAX(node.position.bgn, node.position.end);
  if (frgEnd > _length)
    _length = frgEnd;

  dovetail_path_ptr->push_back(node);

  if ((report) || (node.position.bgn < 0) || (node.position.end < 0)) {
    int32 len = debugfi->fragmentLength(node.ident);
    int32 pos = (node.position.end > node.position.bgn) ? (node.position.end - node.position.bgn) : (node.position.bgn - node.position.end);

    if (node.contained)
      fprintf(stderr, "Added frag %d (len %d) to unitig %d at %d,%d (idx %d) (lendiff %d) (contained in %d)\n",
              node.ident, len, _id, node.position.bgn, node.position.end,
              dovetail_path_ptr->size() - 1,
              pos - len,
              node.contained);
    else
      fprintf(stderr, "Added frag %d (len %d) to unitig %d at %d,%d (idx %d) (lendiff %d)\n",
              node.ident, len, _id, node.position.bgn, node.position.end,
              dovetail_path_ptr->size() - 1,
              pos - len);
  }

  assert(node.position.bgn >= 0);
  assert(node.position.end >= 0);
}


//  This will add a contained fragment to a unitig, adjusting the position as needed.  It is only
//  needed when moving a contained read from unitig A to unitig B.  It is NOT needed when rebuilding
//  a unitig.
//
bool
Unitig::addContainedFrag(int32 fid, BestContainment *bestcont, bool report) {
  DoveTailNode  frag;
  DoveTailNode *parent = NULL;

#ifdef WITHIMP
  frag.type         = AS_READ;
#endif
  frag.ident        = fid;
  frag.contained    = bestcont->container;
  frag.parent       = bestcont->container;
  frag.ahang        = 0;
  frag.bhang        = 0;
  frag.position.bgn = 0;
  frag.position.end = 0;
#ifdef WITHIMP
  frag.delta_length = 0;
  frag.delta        = NULL;
#endif

  parent = &(*dovetail_path_ptr)[pathPosition(frag.contained)];

#if 0
  //  This block is useful for debugging (maybe).  It is usually triggered only during popBubbles(),
  //  when we try to place a contained fragment into a fragment that has not been moved into the new
  //  unitig yet.  It might be useful if pathPosition ever gets messed up.
  //
  if ((parent == NULL) || (parent->ident != frag.contained)) {
    fprintf(stderr, "WARNING:  Didn't find the correct parent frag (%d) for contained frag %d.\n",
            frag.contained, frag.ident);
    fprintf(stderr, "          Found frag %d instead.\n", (parent == NULL) ? -1 : parent->ident);

    parent = NULL;

    for (int fi=0; fi<dovetail_path_ptr->size(); fi++) {
      DoveTailNode *ix = &(*dovetail_path_ptr)[fi];

      fprintf(stderr, "          path[%4d,%4d] is frag %d %s\n",
              fi, pathPosition(ix->ident),
              ix->ident,
              (ix->ident == frag.contained) ? " CORRECT PARENT!" : "");

      if (ix->ident == frag.contained)
        parent = ix;
    }
  }
#endif

  if ((parent == NULL) || (parent->ident != frag.contained)) {
    fprintf(stderr, "Unitig::addContainedFrag()-- WARNING:  Failed to place frag %d into unitig %d; parent not here.\n",
            fid, id());
    return(false);
  }
  
  //  Adjust orientation.  See comments in AS_BOG_Unitig.cc::placeContains().
  //
  //  NOTE!  Code is duplicated there.
  //
  //  NOTE!  The hangs are from the (parent) container to the (child) containee.  This is opposite
  //  as to how dovetail edges are stored.

  assert(bestcont->a_hang >= 0);
  assert(bestcont->b_hang <= 0);

  if (parent->position.bgn < parent->position.end) {
    //  Container is forward.
    frag.ahang = bestcont->a_hang;
    frag.bhang = bestcont->b_hang;

    if (bestcont->sameOrientation) {
      //  ...and so is containee.
      frag.position.bgn = parent->position.bgn + frag.ahang;
      frag.position.end = parent->position.end + frag.bhang;
    } else {
      //  ...but containee is reverse.
      frag.position.bgn = parent->position.end + frag.bhang;
      frag.position.end = parent->position.bgn + frag.ahang;
    }

  } else {
    //  Container is reverse.
    frag.ahang = -bestcont->b_hang;
    frag.bhang = -bestcont->a_hang;

    if (bestcont->sameOrientation) {
      //  ...and so is containee.
      frag.position.bgn = parent->position.bgn + frag.bhang;
      frag.position.end = parent->position.end + frag.ahang;
    } else {
      //  ...but containee is forward.
      frag.position.bgn = parent->position.end + frag.ahang;
      frag.position.end = parent->position.bgn + frag.bhang;
    }
  }


  //  Containments are particularily painful.  A beautiful example: a fragment of length 253bp is
  //  contained in a fragment of length 251bp (both hangs are zero).  In this case, the
  //  "ahang+length" method fails, placing the contained fragment outside the container (and if
  //  backwards oriented, _BEFORE_ the contained fragment).  The "ahang,bhang" method works here,
  //  but fails on other instances, shrinking deep containments to nothing.
  //
  //  We can use either method first, then adjust using the other method.
  //
  //  We'll use 'ahang,bhang' first (mostly because it was already done, and we need to compute
  //  those values anyway) then reset the end based on the length, limited to maintain a
  //  containment relationship.
  //
#warning not knowing the overlap length really hurts.
  if (frag.position.bgn < frag.position.end) {
    frag.position.end = frag.position.bgn + bog->fragmentLength(frag.ident);
    if (frag.position.end > MAX(parent->position.bgn, parent->position.end))
      frag.position.end = MAX(parent->position.bgn, parent->position.end);
  } else {
    frag.position.bgn = frag.position.end + bog->fragmentLength(frag.ident);
    if (frag.position.bgn > MAX(parent->position.bgn, parent->position.end))
      frag.position.bgn = MAX(parent->position.bgn, parent->position.end);
  }

  //  So we can sort properly, we now abuse the delta_length field, and save the containment depth
  //  of each fragment.

  frag.delta_length = parent->delta_length + 1;


  addFrag(frag, 0, report);

#if 0
  //  Bump that new fragment up to be in the correct spot -- we can't
  //  use the sort() method on Unitig, since we lost the
  //  containPartialOrder.
  //
  int             i = dovetail_path_ptr->size() - 1;
  DoveTailNode   *f = &dovetail_path_ptr->front();

  //  Only needed if the frag we just added (i) begins before the second to last frag.

  if (MIN(f[i].position.bgn, f[i].position.end) < MIN(f[i-1].position.bgn, f[i-1].position.end)) {
    DoveTailNode    containee    = f[i];
    int             containeeMin = MIN(containee.position.bgn, containee.position.end);

    while ((i > 0) &&
           (containee.contained != f[i-1].ident) &&
           (containeeMin < MIN(f[i-1].position.bgn, f[i-1].position.end))) {
      f[i] = f[i-1];
      i--;
    }

    f[i] = containee;
  }
#endif

  return(true);
}


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
  DoveTailNode frag;

#ifdef WITHIMP
  frag.type         = AS_READ;
#endif
  frag.ident        = fid;
  frag.contained    = 0;
  frag.parent       = 0;
  frag.ahang        = 0;
  frag.bhang        = 0;
  frag.position.bgn = 0;
  frag.position.end = 0;
#ifdef WITHIMP
  frag.delta_length = 0;
  frag.delta        = NULL;
#endif

  //  The length of the overlap depends only on the length of the a frag and the hangs.  We don't
  //  actually care about the real length (except for logging), only which is thicker.

  if ((bestedge5) && (fragIn(bestedge5->frag_b_id) == id())) {
    bidx5 = pathPosition(bestedge5->frag_b_id);
    blen5 = debugfi->fragmentLength(fid) + ((bestedge5->ahang < 0) ? bestedge5->bhang : -bestedge5->ahang);
    //fprintf(stderr, "bestedge5:  %d,%d,%d,%d len %d\n", bestedge5->frag_b_id, bestedge5->bend, bestedge5->ahang, bestedge5->bhang, blen5);
    assert(bestedge5->frag_b_id == (*dovetail_path_ptr)[bidx5].ident);
  }

  if ((bestedge3) && (fragIn(bestedge3->frag_b_id) == id())) {
    bidx3 = pathPosition(bestedge3->frag_b_id);;
    blen3 = debugfi->fragmentLength(fid) + ((bestedge3->ahang < 0) ? bestedge3->bhang : -bestedge3->ahang);
    //fprintf(stderr, "bestedge3:  %d,%d,%d,%d len %d\n", bestedge3->frag_b_id, bestedge3->bend, bestedge3->ahang, bestedge3->bhang, blen3);
    assert(bestedge3->frag_b_id == (*dovetail_path_ptr)[bidx3].ident);
  }

  //  Use the longest that exists -- an alternative would be to take the average position, but that
  //  could get messy if the placements are different.  Picking one or the other has a better chance
  //  of working, though it'll fail if the fragment is chimeric or spans something it shouldn't,
  //  etc.

  if ((blen5 == 0) && (blen3 == 0)) {
    fprintf(stderr, "Unitig::addAndPlaceFrag()-- WARNING:  Failed to place frag %d into unitig %d; no edges to the unitig.\n",
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

    for (int fi=0; fi<dovetail_path_ptr->size(); fi++) {
      DoveTailNode *tfrg = &(*dovetail_path_ptr)[fi];

      tfrg->position.bgn += frgBgn;
      tfrg->position.end += frgBgn;
    }
  }

  //  Finally, add the fragment.

  addFrag(frag, 0, report);

  return(true);
}


float Unitig::getAvgRho(FragmentInfo *fi){

  if(dovetail_path_ptr->size() == 1)
    _avgRho = 1;
  if(_avgRho!=-1)
    return(_avgRho);


  // We will compute the average rho.
  //
  // Since rho is the length(unitig) - length(last fragment),
  //   and the length(last fragment) is ambiguous depending on which
  //   direction we are walking the unitig from.  We will take the average
  //   of the rhos through both directions.

  DoveTailPath::const_iterator dtp_iter;
  if (dovetail_path_ptr == NULL || dovetail_path_ptr->empty()) {
    //fprintf(stderr,"NULL dovetailpath\n");
    return 1;
  }

  // Get first fragment's length
  dtp_iter=dovetail_path_ptr->begin();
  int ident1 = dtp_iter->ident;
  long first_frag_len = fi->fragmentLength(dtp_iter->ident);
  assert(first_frag_len > 0);

  // Get last fragment's length
  dtp_iter=dovetail_path_ptr->end();
  dtp_iter--;
  int ident2 = dtp_iter->ident;
  long last_frag_len = fi->fragmentLength(dtp_iter->ident);
  assert(last_frag_len > 0);

  // Get average of first and last fragment lengths
  double avg_frag_len = (last_frag_len + first_frag_len)/2.0;

  // Compute average rho
  long unitig_length=getLength();
  _avgRho = unitig_length - avg_frag_len;

  if (_avgRho <= 0 ) {
    //fprintf(stderr, "Negative Rho ident1 "F_IID" ident2 "F_IID" unitig_length %d first_frag_len %d last_frag_len %d avg_frag_len %f\n",
    //        ident1, ident2, unitig_length, first_frag_len, last_frag_len, avg_frag_len);
    _avgRho = 1;
  }
  return(_avgRho);
}


float Unitig::_globalArrivalRate = -1;

void Unitig::setGlobalArrivalRate(float global_arrival_rate){
  _globalArrivalRate=global_arrival_rate;
}
void Unitig::setLocalArrivalRate(float local_arrival_rate){

  if ( local_arrival_rate < std::numeric_limits<float>::epsilon())
    _localArrivalRate = 0;
  else
    _localArrivalRate = local_arrival_rate;
}
float Unitig::getLocalArrivalRate(FragmentInfo *fi){
  if (_localArrivalRate != -1 )
    return _localArrivalRate;
  setLocalArrivalRate((getNumFrags() - 1) / getAvgRho(fi));
  return _localArrivalRate;
}


float Unitig::getCovStat(FragmentInfo *fi){

  const float ln2=0.69314718055994530941723212145818;

  // Note that we are using numFrags in this calculation.
  //   If the fragments in the unitig are not randomly sampled
  //   from the genome, this calculation will be wrong.
  //   Ie. if fragments being assembled are from a locally
  //   sequenced batch, the region may look artificially repetitive.
  //   The value should really be "number of randomly sampled
  //   fragments in the unitig".

  if(_globalArrivalRate==-1)
    fprintf(stderr, "You have not set the _globalArrivalRate variable.\n");

  if(_covStat == FLT_MAX){
    if(_globalArrivalRate > 0.0){
      _covStat = (getAvgRho(fi) * _globalArrivalRate) - (ln2 * (getNumFrags() -1));
    }else{
      _covStat = 0.0;
    }
  }

  return(_covStat);
}


void Unitig::reverseComplement(bool doSort) {
  DoveTailIter iter   = dovetail_path_ptr->begin();

  //  If there are contained fragments, we need to sort by position to place them correctly after
  //  their containers.  If there are no contained fragments, sorting can break the initial unitig
  //  building.  When two frags start at position zero, we'll exchange the order.  Initial unitig
  //  building depends on having the first fragment added become the last fragment in the unitig
  //  after reversing.

  for(; iter != dovetail_path_ptr->end(); iter++) {
    iter->position.bgn = getLength() - iter->position.bgn;
    iter->position.end = getLength() - iter->position.end;

    //if (iter->contained != 0)
    //  doSort = true;

    assert(iter->position.bgn >= 0);
    assert(iter->position.end >= 0);
  }

  //  We've updated the positions of everything.  Now, sort or reverse the list, and rebuild the
  //  pathPosition map.

  if (doSort) {
    assert(0);
    sort();
  } else {
    std::reverse(dovetail_path_ptr->begin(), dovetail_path_ptr->end());

    for (int fi=0; fi<dovetail_path_ptr->size(); fi++)
      _pathPosition[(*dovetail_path_ptr)[fi].ident] = fi;
  }
}



int
DoveTailNodeCmp(const void *a, const void *b){
  DoveTailNode *impa = (DoveTailNode *)a;
  DoveTailNode *impb = (DoveTailNode *)b;

  int32 abgn = (impa->position.bgn < impa->position.end) ? impa->position.bgn : impa->position.end;
  int32 aend = (impa->position.bgn < impa->position.end) ? impa->position.end : impa->position.bgn;

  int32 bbgn = (impb->position.bgn < impb->position.end) ? impb->position.bgn : impb->position.end;
  int32 bend = (impb->position.bgn < impb->position.end) ? impb->position.end : impb->position.bgn;

  //  NEWSORT does not work.  When bubbles are popped, we add non-contained fragments to
  //  a unitig, but just stick them at the end of the list.  NEWSORT would then maintain
  //  this ordering, which is an error.
  //
#undef NEWSORT

#ifdef NEWSORT
  bool  aIsCont = bog->isContained(impa->ident);
  bool  bIsCont = bog->isContained(impb->ident);

  if ((aIsCont == false) && (bIsCont == false))
    //  Both dovetail nodes, keep same order
    return((int)impa->delta_length - (int)impb->delta_length);
#endif

  if (abgn != bbgn)
    //  Return negative for the one that starts first.
    return(abgn - bbgn);

  if (aend != bend)
    //  Return negative for the one that ends last.
    return(bend - aend);

#ifdef NEWSORT
  if (bIsCont == true)
    //  b is contained in a, so it comes after a.
    return(-1);

 if (aIsCont == true)
    //  a is contained in b, so it comes after b.
    return(1);
#endif

  //  Both contained, fallback on depth added, negative for earliest added
  return((int)impa->delta_length - (int)impb->delta_length);
}


void
Unitig::sort(void) {

#ifdef NEWSORT
  for (int fi=0; fi<dovetail_path_ptr->size(); fi++) {
    DoveTailNode *f = &((*dovetail_path_ptr)[fi]);

    if (bog->isContained(f->ident) == false)
      f->delta_length = fi;
  }
#endif

  qsort( &(dovetail_path_ptr->front()), getNumFrags(), sizeof(DoveTailNode), &DoveTailNodeCmp );

  for (int fi=0; fi<dovetail_path_ptr->size(); fi++)
    _pathPosition[(*dovetail_path_ptr)[fi].ident] = fi;
}

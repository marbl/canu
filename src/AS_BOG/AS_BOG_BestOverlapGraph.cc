
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

#include<iostream>
#include<vector>
#include<limits>
#include<cmath>

#include "AS_BOG_BestOverlapGraph.hh"

#undef max

//  The overlap, pi, exists between A and B:
//
//  A -------------->
//         |||||||||
//  B      ---------------->
//
//  AEnd(pi) is 3'
//  BEnd(pi) is 5'
//

fragment_end_type BestOverlapGraph::AEnd(const OVSoverlap& olap) {
  if (olap.dat.ovl.a_hang < 0 && olap.dat.ovl.b_hang < 0)
    return FIVE_PRIME;
  if (olap.dat.ovl.a_hang > 0 && olap.dat.ovl.b_hang > 0)
    return THREE_PRIME;

  assert(0); // no contained
}

fragment_end_type BestOverlapGraph::BEnd(const OVSoverlap& olap) {
  if (olap.dat.ovl.a_hang < 0 && olap.dat.ovl.b_hang < 0)
    if ( olap.dat.ovl.flipped )
      return FIVE_PRIME;
    else
      return THREE_PRIME;

  if (olap.dat.ovl.a_hang > 0 && olap.dat.ovl.b_hang > 0)
    if ( olap.dat.ovl.flipped )
      return THREE_PRIME;
    else
      return FIVE_PRIME;

  assert(0); // no contained
}


// Create BestOverlapGraph as an array of size max fragments.
//     Assuming that our iuids start at index value of 1.
//
// AS_UTG_ERROR_RATE is fraction error, same as AS_CNS_ERROR_RATE.
//
BestOverlapGraph::BestOverlapGraph(FragmentInfo        *fi,
                                   OverlapStore        *ovlStore,
                                   double               AS_UTG_ERROR_RATE) {
  curFrag = 0;

  _fi = fi;

  _best_overlaps = new BestFragmentOverlap [fi->numFragments() + 1];

  memset(_best_overlaps, 0, sizeof(BestFragmentOverlap) * (fi->numFragments() + 1));

  assert(AS_UTG_ERROR_RATE >= 0.0);
  assert(AS_UTG_ERROR_RATE <= AS_MAX_ERROR_RATE);

  assert(AS_CNS_ERROR_RATE >= 0.0);
  assert(AS_CNS_ERROR_RATE <= AS_MAX_ERROR_RATE);

  fprintf(stderr, "BestOverlapGraph()-- UTG erate %.4f%%, CNS erate %.4f%%\n",
          100.0 * AS_UTG_ERROR_RATE, 100.0 * AS_CNS_ERROR_RATE);

  mismatchCutoff  = AS_OVS_encodeQuality( AS_UTG_ERROR_RATE );
  consensusCutoff = AS_OVS_encodeQuality( AS_CNS_ERROR_RATE );

  // Go through the overlap stream in two passes:
  // first pass finds the containments
  // second pass builds the overlap graph, excluding contained frags

  OVSoverlap olap;

  AS_OVS_resetRangeOverlapStore(ovlStore);
  while  (AS_OVS_readOverlapFromStore(ovlStore, &olap, AS_OVS_TYPE_OVL))
    scoreContainment( olap );

  AS_OVS_resetRangeOverlapStore(ovlStore);
  while  (AS_OVS_readOverlapFromStore(ovlStore, &olap, AS_OVS_TYPE_OVL))
    scoreEdge( olap );

  updateInDegree();

  //  Diagnostic.  Dump the best edges, count the number of contained
  //  reads, etc.
  {
    FILE *BC = fopen("best.contains", "w");
    FILE *BE = fopen("best.edges", "w");

    if ((BC) && (BE)) {
      fprintf(BC, "#fragId\tlibId\tbestCont\n");
      fprintf(BE, "#fragId\tlibId\tbest5\tbest3\n");

      for (int id=1; id<_fi->numFragments() + 1; id++) {
        BestContainment *bestcont  = getBestContainer(id);
        BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, FIVE_PRIME);
        BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, THREE_PRIME);

        if (bestcont)
          fprintf(BC, "%d\t%d\t%d\n", id, _fi->libraryIID(id), bestcont->container);
        else
          fprintf(BE, "%d\t%d\t%d\t%d\n", id, _fi->libraryIID(id), bestedge5->frag_b_id, bestedge3->frag_b_id);
      }

      fclose(BC);
      fclose(BE);
    }
  }
}

BestOverlapGraph::~BestOverlapGraph(){
  delete[] _best_overlaps;
}


//  Given a fragment IUID and which end, returns pointer to
//  BestOverlap node.  
BestEdgeOverlap *BestOverlapGraph::getBestEdgeOverlap(iuid frag_id, fragment_end_type which_end){
  if(which_end == FIVE_PRIME)
    return(&_best_overlaps[frag_id].five_prime);
  else if(which_end == THREE_PRIME){
    return(&_best_overlaps[frag_id].three_prime);
  }
}

BestEdgeOverlap *BestOverlapGraph::getBestEdgeOverlap(FragmentEnd* end) {
  return getBestEdgeOverlap(end->fragId(),end->fragEnd());
}

void BestOverlapGraph::followOverlap(FragmentEnd* end) {
  BestEdgeOverlap* edge = getBestEdgeOverlap(end);
  *end = FragmentEnd(edge->frag_b_id, (edge->bend == FIVE_PRIME) ? THREE_PRIME : FIVE_PRIME);
}

//  Given an overlap, determines which record (iuid and end) and sets the newScore.
void BestOverlapGraph::setBestEdgeOverlap(const OVSoverlap& olap, float newScore) {

  if (AEnd(olap) == THREE_PRIME) {
    _best_overlaps[ olap.a_iid ].three_prime.frag_b_id = olap.b_iid;
    _best_overlaps[ olap.a_iid ].three_prime.score     = newScore;
    _best_overlaps[ olap.a_iid ].three_prime.bend      = BEnd(olap);
    _best_overlaps[ olap.a_iid ].three_prime.ahang     = olap.dat.ovl.a_hang;
    _best_overlaps[ olap.a_iid ].three_prime.bhang     = olap.dat.ovl.b_hang;

  }
  if (AEnd(olap) == FIVE_PRIME) {
    _best_overlaps[ olap.a_iid ].five_prime.frag_b_id = olap.b_iid;
    _best_overlaps[ olap.a_iid ].five_prime.score     = newScore;
    _best_overlaps[ olap.a_iid ].five_prime.bend      = BEnd(olap);
    _best_overlaps[ olap.a_iid ].five_prime.ahang     = olap.dat.ovl.a_hang;
    _best_overlaps[ olap.a_iid ].five_prime.bhang     = olap.dat.ovl.b_hang;
  }
}

void BestOverlapGraph::addContainEdge( iuid contain, iuid otherRead ) {
  _best_containments[ contain ].overlapsAreSorted = false;
  _best_containments[ contain ].overlaps.push_back( otherRead );
}

bool BestOverlapGraph::containHaveEdgeTo( iuid contain, iuid otherRead ) {

  if (_best_containments.find( contain ) == _best_containments.end())
    return(false);

  if (_best_containments[ contain ].overlapsAreSorted == false) {
    std::sort(_best_containments[ contain ].overlaps.begin(),
              _best_containments[ contain ].overlaps.end());
    _best_containments[ contain ].overlapsAreSorted = true;
  }

  return(std::binary_search(_best_containments[ contain ].overlaps.begin(),
                            _best_containments[ contain ].overlaps.end(),
                            otherRead));
}

void BestOverlapGraph::setBestContainer(const OVSoverlap& olap, float newScr) {
  BestContainment newBest;

  newBest.container         = olap.a_iid;
  newBest.contain_score     = newScr;
  newBest.a_hang            = olap.dat.ovl.a_hang;
  newBest.b_hang            = olap.dat.ovl.b_hang;
  newBest.sameOrientation   = olap.dat.ovl.flipped ? false : true;
  newBest.isPlaced          = false;
  newBest.overlapsAreSorted = false;

  //fprintf(stderr, "bestContainer of frag %d is frag %d score = %f.\n", olap.b_iid, olap.a_iid, newScr);

  _best_containments[ olap.b_iid ] = newBest;
}


// Transitively removes redundant containments, so all containees in a container, refer
// to the same container.  Algorithm will go through each element in the list of contained
// fragments, and then for each element follow each container's container.
void BestOverlapGraph::removeTransitiveContainment() {

  // Loop through each containee that has been stored in _best_containments
  for(std::map<iuid,BestContainment>::const_iterator it = _best_containments.begin();
      it != _best_containments.end(); it++) {
    iuid id = it->first;                     // Gets the iuid of the containee under analysis
    BestContainment bst = it->second;        // Gets the BestContainment record of the containee under analysis

    // Get container information based on containee id/containment record
    bool sameOrient = bst.sameOrientation;
    bool useAhang = true;

    // For this std::map:
    //   iuid is the containee iuid
    //   BestContainment.container is the best container iuid
    //       This finds the iterator if the containerOf(id) is contained, else returns
    //       _best_containments.end() if not found.
    //
    // When this find returns, the i2 will point to id's container's container
    std::map<iuid,BestContainment>::iterator i2 =
      _best_containments.find( bst.container);

    // Keep track of which container's we've looked at for each transitive path starting from
    //   the containee under analysis.
    std::map<iuid,BestContainment> found;
    found[bst.container] = bst;

    // Loop while the current container is a containee of another container 
    while ( i2 != _best_containments.end() ) {
      BestContainment nb = i2->second;
      std::cout << id <<" "<<bst.container<<" "<< nb.container<< std::endl;

      // Delete containee under analysis from _best_containments if its container is concontained
      //   by itself.  ie if id's container is contained by id.  This eliminates ciruclar containment.
      if ( nb.container == id ) {
        _best_containments.erase( id );
        std::cout << "Erase self" << std::endl;
        break;
      }

      // Look through the list of containers we've already walked past to find
      //   circular containment greater than one degree away.
      std::map<iuid,BestContainment>::iterator seen= found.find( nb.container );

      // If i2's container has already been traversed
      if ( seen != found.end() ) { 

        // Set id's container to the larger container.
        _best_containments[ id ] = seen->second;
        std::cout << "Circled " << seen->second.container<< std::endl;

        // Remove the container of id's new larger container
        _best_containments.erase( seen->second.container );
        break;
      }

      // Set id's container to id's container's container, in this case
      //   one transitive step away.
      _best_containments[ id ] = nb;

      // Keep track of what containers have been used.
      found[ nb.container ] = nb;

      // Reset the orientation of id's containment.
      if (sameOrient) {
        if (nb.sameOrientation) {
          useAhang = true;
          sameOrient =_best_containments[id].sameOrientation=true;
        } else {
          sameOrient =_best_containments[id].sameOrientation=false;
          useAhang = false;
        }
      } else {
        if (nb.sameOrientation) {
          sameOrient =_best_containments[id].sameOrientation=false;
          useAhang = true;
        } else {
          useAhang = false;
          sameOrient =_best_containments[id].sameOrientation=true;
        }
      }
      found[nb.container].sameOrientation = sameOrient;
      // make hang relative to new container
      _best_containments[ id ].a_hang += useAhang ? bst.a_hang : abs(bst.b_hang); 

      // Iterator to the i2's container
      i2 = _best_containments.find( nb.container);
    }
  }
}


// Frag Store related data structures

void BestOverlapGraph::updateInDegree() {
  if (curFrag == 0)
    return;
  if ( ! isContained(curFrag) ) {
    // Update B's in degree on A's 3' End
    iuid bid = _best_overlaps[ curFrag ].three_prime.frag_b_id;
    switch(_best_overlaps[ curFrag ].three_prime.bend){
      case THREE_PRIME:
        _best_overlaps[ bid ].three_prime.in_degree++; break;
      case FIVE_PRIME:
        _best_overlaps[ bid ].five_prime.in_degree++; break;
      default: assert(0);
    }

    // Update B's in degree on A's 5' End
    bid = _best_overlaps[ curFrag ].five_prime.frag_b_id;
    switch(_best_overlaps[ curFrag ].five_prime.bend){
      case THREE_PRIME:
        _best_overlaps[ bid ].three_prime.in_degree++; break;
      case FIVE_PRIME:
        _best_overlaps[ bid ].five_prime.in_degree++; break;
      default: assert(0);
    }
  }
}
bool BestOverlapGraph::checkForNextFrag(const OVSoverlap& olap) {
  // Update the in_degrees whenever the incoming overlap's A's Fragment IUID changes.
  //   Returns true, if the olap's A fragment ID has changed since the previous call,
  //   else false.
  //
  // This method is only called by processOverlap, should make it private.
  //
  // Since the overlaps are coming in in order of A's iuid, 
  //   its safe for us to assume that once the incoming overlap A iuid changes
  //   that we are completely done processing overlaps from A's
  //   fragment.   This implies that we know that A's best overlap
  //   at this point will not change, so it's safe to update it's overlapping
  //   partner's (B's), in degree for both ends of A.
        
  // In this code, the olap.a_iid is considered the "next" frag, so 
  //   curFrag is the IID of the fragment prior to receiving this olap.
  if (curFrag != olap.a_iid) {

    // update in degree if not contained
    updateInDegree();

    // Set up the curFrag to refer to the incoming (next) fragment IUID	
    curFrag = olap.a_iid;

    // Initialize the overlap score for A
    _best_overlaps[ olap.a_iid ].three_prime.score = 0;
    _best_overlaps[ olap.a_iid ].five_prime.score = 0;
    bestLength = 0;

    // Means that A's fragment ID has changed since the previous
    return true;
  }

  // Means that A's fragment ID is the same
  return false;
}
void BestOverlapGraph::scoreContainment(const OVSoverlap& olap) {

  // in the case of no hang, make the lower frag the container
  if ((olap.dat.ovl.a_hang == 0) &&
      (olap.dat.ovl.b_hang == 0) &&
      (olap.a_iid > olap.b_iid))
    return;

  //  We only care if A contains B.

  if ((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0)) {
    BestContainment *best = getBestContainer(olap.b_iid);

    float newScr = scoreOverlap(olap);

    if (newScr <= 0)
      return;

    if ((NULL == best) ||
        (newScr > best->contain_score)) {
      setBestContainer(olap, newScr);
    }
  }
}

void BestOverlapGraph::scoreEdge(const OVSoverlap& olap) {
  // This function builds the BestOverlapGraph by considering the specified
  //   overlap record.  It's important to remember that the overlaps
  //   must be passed in in sorted order of fragment A's iuid.  Also,
  //   it is also very important that overlaps are doubly listed, ie.
  //   if an overlap exists, an olap must be passed in where A overlaps with
  //   B and an olap must later be passed in where B overlaps A.
  //
  // It calls the virtual score function to score the overlap,
  //   determines whether the overlap is containment or dovetailing,
  //   then stores the overlap in the BestOverlapGraph member variables.

  //  Store real edges from contained frags to help with unhappy
  //  mate splitting
  //
  //  If A was previously found to be a containee, and
  //
  //  a < 0           b < 10
  //        ----------..
  //  ......------------
  //
  //  a > 0           b > -10
  //  ......-----------
  //        ---------.......
  //
  //  That looks wrong.
  //
  //  From Eli: These are contained, but close either way.  We're
  //  storing the non-containment edges for this fragment, plus a
  //  few containment edges that are "close" to being dovetails.  "I
  //  think there are cases when a change in the alignemtn
  //  (consensus) will change which one is contained and screw up
  //  the order, so having this 10 base fudge factor helps things
  //  work out."
  //
  if (isContained(olap.a_iid)) {
    if ((olap.dat.ovl.a_hang < 0 && olap.dat.ovl.b_hang < 10) ||
        (olap.dat.ovl.a_hang > 0 && olap.dat.ovl.b_hang > -10) )
      addContainEdge(olap.a_iid, olap.b_iid);
    return;
  }

  //  Skip contained fragments.
  if (isContained(olap.a_iid) || isContained(olap.b_iid))
    return;

  //  Skip containment overlaps.  Can this happen?  Yup.  How?
  //  The overlap could be above our allowed error.
  //
  if (((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0)) ||
      ((olap.dat.ovl.a_hang <= 0) && (olap.dat.ovl.b_hang >= 0)))
    return;

  float newScr = scoreOverlap(olap);

  //  If the score is 0, the overlap doesn't pass the scoring
  //  criteria at all so don't store the overlap whether or not
  //  it's dovetailing or containment.
  if ( newScr <= 0 )
    return;

  //  Dove tailing overlap

  //  Update the in degree for what A overlaps with if the current
  //  A fragment is different than the previous one.
  checkForNextFrag(olap);

  BestEdgeOverlap *best = getBestEdgeOverlap( olap.a_iid, AEnd(olap));
  short olapLen = olapLength(olap);

  // Store the overlap if:
  //   1.)  The score is better than what is already in the graph
  //   2.)  If the scores are identical, the one with the longer length
  //
  // Since the order of how the overlaps are read in from the overlap
  //   store are by A's increasing iuid, by default, if the score
  //   and length are the same, the iuid of the lower value will be
  //   kept.

  if ((newScr > best->score) ||
      ((newScr == best->score) && (olapLen > bestLength))) {
    setBestEdgeOverlap( olap, newScr );
    bestLength = olapLen;
  }
}

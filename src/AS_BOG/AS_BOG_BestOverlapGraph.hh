
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

#ifndef INCLUDE_AS_BOG_BESTOVERLAPGRAPH
#define INCLUDE_AS_BOG_BESTOVERLAPGRAPH

#include "AS_BOG_Datatypes.hh"

#include <set>

struct BestOverlapGraph {
  BestOverlapGraph(FragmentInfo *fi, OverlapStore *ovlStore, double erate);
  ~BestOverlapGraph();

  // Accessor Functions
  BestEdgeOverlap *getBestEdgeOverlap(iuid frag_id, fragment_end_type which_end);
  BestEdgeOverlap *getBestEdgeOverlap(FragmentEnd*);

  // given a FragmentEnd sets it to the next FragmentEnd after following the
  // best edge
  void followOverlap(FragmentEnd*);

  bool isContained(const iuid fragid) {
    return(_best_contains[fragid].isContained);
  };

  // Given a containee, returns pointer to BestContainment record
  BestContainment *getBestContainer(const iuid fragid) {
    return((isContained(fragid)) ? &_best_contains[fragid] : NULL);
  };

  bool containHaveEdgeTo( iuid, iuid);

  // Graph building methods
  fragment_end_type AEnd(const OVSoverlap& olap);
  fragment_end_type BEnd(const OVSoverlap& olap);
  void processOverlap(const OVSoverlap& olap);

  float scoreOverlap(const OVSoverlap& olap) {

#if 0
    // Computes the score for a Error Rate BOG based on overlap
    // corrected error rate.
    //
    // Error rate is normalized so that the higher the error
    // rate, the lower the score.
    //
    return(100.0 - AS_OVS_decodeQuality(olap.dat.ovl.corr_erate) * 100.0);
#endif

#if 0
    // Computes the score for a Longest Edge BOG based on
    // overlap length only.
    //
    return(olapLength(olap));
#endif

#if 0
    // The log for this is:
    //   Add alternate scoring schema that wasn't tested much for posterity's sake.
    //
    if (olap.dat.ovl.orig_erate > consensusCutoff)
      return 0;
    if (olap.dat.ovl.corr_erate > mismatchCutoff)
      return 0;
    return(olapLength(olap) / sqrt(1 + olap.dat.ovl.corr_erate));
#endif

#if 0
    // Computes the score for a Longest Edge BOG based on
    // overlap length but after applying an an error rate
    // cutoff.
    //
    if (olap.dat.ovl.orig_erate > consensusCutoff)
      return 0;
    if (olap.dat.ovl.corr_erate > mismatchCutoff)
      return 0;
    return olapLength(olap);
#endif

#if 1
    //  BPW's new score
    if (olap.dat.ovl.orig_erate > consensusCutoff)
      return 0;
    if (olap.dat.ovl.corr_erate > mismatchCutoff)
      return 0;

    int a_hang = olap.dat.ovl.a_hang;
    int b_hang = olap.dat.ovl.b_hang;

    int alen = _fi->fragmentLength(olap.a_iid);
    int blen = _fi->fragmentLength(olap.b_iid);

    double qlt = 1.0 - AS_OVS_decodeQuality(olap.dat.ovl.corr_erate) - AS_OVS_decodeQuality(olap.dat.ovl.orig_erate) / 10000;

    //  Containments - the length of the overlaps are all the same.
    //  We return the quality.
    //
    if ((a_hang >= 0) && (b_hang <= 0))
      return(qlt);

    if ((a_hang <= 0) && (b_hang >= 0))
      return(qlt);

    //  Dovetails - the length of the overlap is the score, but we
    //  bias towards lower error.
    //
    if ((a_hang < 0) && (b_hang < 0))
      return(alen + b_hang + 1.0 - qlt);

    if ((a_hang > 0) && (b_hang > 0))
      return(alen - a_hang + 1.0 - qlt);

    //  Fail if we're here.
    assert(0);
#endif
  };

  uint16  olapLength(iuid a_iid, iuid b_iid, short a_hang, short b_hang) {
    int alen = _fi->fragmentLength(a_iid);
    int blen = _fi->fragmentLength(b_iid);

    if (a_hang < 0) {
      if (b_hang < 0 )
        return alen + b_hang;
      else
        return blen + a_hang - b_hang; // spur or containment
    } else {
      if (b_hang < 0 )
        return alen + b_hang - a_hang; // spur or containment
      else
        return alen - a_hang;
    }
  };

  uint16 olapLength(const OVSoverlap& olap) {
    return(olapLength(olap.a_iid, olap.b_iid, olap.dat.ovl.a_hang, olap.dat.ovl.b_hang));
  };

  bool checkForNextFrag(const OVSoverlap& olap);
  void scoreContainment(const OVSoverlap& olap);
  void scoreEdge(const OVSoverlap& olap);

private:
  BestFragmentOverlap *_best_overlaps;
  BestContainment     *_best_contains;
  FragmentInfo        *_fi;

public:
  uint64 mismatchCutoff;
  uint64 consensusCutoff;
}; //BestOverlapGraph



#endif //INCLUDE_AS_BOG_BESTOVERLAPGRAPH


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

static const char *rcsid_INCLUDE_AS_BOG_BESTOVERLAPGRAPH = "$Id: AS_BOG_BestOverlapGraph.hh,v 1.55 2010-01-26 02:23:55 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"

struct BestOverlapGraph {
  BestOverlapGraph(FragmentInfo *fi, OverlapStore *ovlStoreUniq, OverlapStore *ovlStoreRept, double erate);
  ~BestOverlapGraph();

  //  Given a fragment UINT32 and which end, returns pointer to
  //  BestOverlap node.
  BestEdgeOverlap *getBestEdgeOverlap(uint32 frag_id, uint32 which_end) {
    if(which_end == FIVE_PRIME)    return(&_best_overlaps[frag_id].five_prime);
    if(which_end == THREE_PRIME)   return(&_best_overlaps[frag_id].three_prime);
    return(NULL);
  };

  // given a FragmentEnd sets it to the next FragmentEnd after following the
  // best edge
  FragmentEnd   followOverlap(FragmentEnd end) {
    if (end.fragId() == 0)
      return(FragmentEnd());

    BestEdgeOverlap *edge = getBestEdgeOverlap(end.fragId(), end.fragEnd());

    return(FragmentEnd(edge->frag_b_id, (edge->bend == FIVE_PRIME) ? THREE_PRIME : FIVE_PRIME));
  };

  bool isContained(const uint32 fragid) {
    return(_best_contains[fragid].isContained);
  };

  // Given a containee, returns pointer to BestContainment record
  BestContainment *getBestContainer(const uint32 fragid) {
    return((isContained(fragid)) ? &_best_contains[fragid] : NULL);
  };

  bool containHaveEdgeTo(uint32 contain, uint32 otherRead) {
    BestContainment  *c = &_best_contains[contain];
    bool              r = false;

    if (c->isContained) {
      if (c->olapsLen < 16) {
        for (int i=0; i<c->olapsLen; i++)
          if (c->olaps[i] == otherRead) {
            r = true;
            break;
          }
      } else {
        if (c->olapsSorted == false) {
          std::sort(c->olaps, c->olaps + c->olapsLen);
          c->olapsSorted = true;
        }
        r = std::binary_search(c->olaps, c->olaps + c->olapsLen, otherRead);
      }
    }

    return(r);
  };

  // Graph building methods
  uint32  AEnd(const OVSoverlap& olap);
  uint32  BEnd(const OVSoverlap& olap);
  void    processOverlap(const OVSoverlap& olap);

  uint64  scoreOverlap(const OVSoverlap& olap) {

    //  BPW's newer new score.  For
    //  the most part, we use the length of the overlap, but we also
    //  want to break ties with the higher quality overlap.
    //
    //  The high 20 bits are the length of the overlap.
    //  The next 12 are the corrected error rate.
    //  The last 12 are the original error rate.
    //
    //  (Well, 12 == AS_OVS_ERRBITS)

    if (olap.dat.ovl.orig_erate > consensusCutoff)
      return 0;
    if (olap.dat.ovl.corr_erate > mismatchCutoff)
      return 0;

    uint64  leng = 0;
    uint64  corr = (AS_OVS_MAX_ERATE - olap.dat.ovl.corr_erate);
    uint64  orig = (AS_OVS_MAX_ERATE - olap.dat.ovl.orig_erate);

    //  Shift AFTER assigning to a 64-bit value to avoid overflows.
    corr <<= AS_OVS_ERRBITS;

    int a_hang = olap.dat.ovl.a_hang;
    int b_hang = olap.dat.ovl.b_hang;

    //  Containments - the length of the overlaps are all the same.
    //  We return the quality.  We possibly do not need this test, as
    //  the leng is only set for dovetail overlaps, but lets play
    //  safe.
    //
    if (((a_hang >= 0) && (b_hang <= 0)) ||
        ((a_hang <= 0) && (b_hang >= 0)))
      return(corr | orig);

    //  Dovetails - the length of the overlap is the score, but we
    //  bias towards lower error.
    //
    if ((a_hang < 0) && (b_hang < 0))
      leng = _fi->fragmentLength(olap.a_iid) + b_hang;

    if ((a_hang > 0) && (b_hang > 0))
      leng = _fi->fragmentLength(olap.a_iid) - a_hang;

    //  Shift AFTER assigning to a 64-bit value to avoid overflows.
    leng <<= (2 * AS_OVS_ERRBITS);

    return(leng | corr | orig);
  };

  uint32  olapLength(uint32 a_iid, uint32 b_iid, short a_hang, short b_hang) {
    uint32 alen = _fi->fragmentLength(a_iid);
    uint32 blen = _fi->fragmentLength(b_iid);

    if (a_hang < 0) {
      if (b_hang < 0 )
        return alen + b_hang;
      else
        return alen + a_hang - b_hang; // spur or containment
    } else {
      if (b_hang < 0 )
        return alen - a_hang + b_hang; // spur or containment
      else
        return alen - a_hang;
    }
  };

  uint32 olapLength(const OVSoverlap& olap) {
    return(olapLength(olap.a_iid, olap.b_iid, olap.dat.ovl.a_hang, olap.dat.ovl.b_hang));
  };

  uint32 fragmentLength(uint32 id) {
    return(_fi->fragmentLength(id));
  };

  bool checkForNextFrag(const OVSoverlap& olap);
  void scoreContainment(const OVSoverlap& olap);
  void scoreEdge(const OVSoverlap& olap);

private:
  BestFragmentOverlap *_best_overlaps;
  BestContainment     *_best_contains;
  FragmentInfo        *_fi;

  uint64              *_best_overlaps_5p_score;
  uint64              *_best_overlaps_3p_score;
  uint64              *_best_contains_score;

public:
  uint64 mismatchCutoff;
  uint64 consensusCutoff;
}; //BestOverlapGraph



#endif //INCLUDE_AS_BOG_BESTOVERLAPGRAPH


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

static const char *rcsid_INCLUDE_AS_BOG_BESTOVERLAPGRAPH = "$Id$";

#include "AS_BOG_Datatypes.H"


struct BestOverlapGraph {
  BestOverlapGraph(OverlapStore *ovlStoreUniq, OverlapStore *ovlStoreRept, double erate, double elimit, const char *prefix);
  ~BestOverlapGraph();

  //  Given a fragment UINT32 and which end, returns pointer to
  //  BestOverlap node.
  BestEdgeOverlap *getBestEdgeOverlap(uint32 frag_id, bool threePrime) {
    if (threePrime)
      return(_best3 + frag_id);
    return(_best5 + frag_id);
  };

  // given a FragmentEnd sets it to the next FragmentEnd after following the
  // best edge
  FragmentEnd   followOverlap(FragmentEnd end) {
    if (end.fragId() == 0)
      return(FragmentEnd());

    BestEdgeOverlap *edge = getBestEdgeOverlap(end.fragId(), end.frag3p());

    return(FragmentEnd(edge->fragId(), !edge->frag3p()));
  };

  bool isContained(const uint32 fragid) {
    return(_bestC[fragid].isContained);
  };

  // Given a containee, returns pointer to BestContainment record
  BestContainment *getBestContainer(const uint32 fragid) {
    return((isContained(fragid)) ? &_bestC[fragid] : NULL);
  };

  bool containHaveEdgeTo(uint32 contain, uint32 otherRead) {
    BestContainment  *c = &_bestC[contain];
    bool              r = false;

    if ((c->olapsLen == 0) || (c->olaps == NULL))
      return(r);

    if (c->isContained == false)
      return(r);

    if (c->olapsLen < 16) {
      for (uint32 i=0; i<c->olapsLen; i++)
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

    return(r);
  };

  // Graph building methods
  void    processOverlap(const OVSoverlap& olap);

  bool    isOverlapBadQuality(const OVSoverlap& olap) {

    //  The overlap is ALWAYS bad if the original error rate is above what we initially required
    //  overlaps to be at.  We shouldn't have even seen this overlap.  This is a bug in the
    //  overlapper.
    //
    if (olap.dat.ovl.orig_erate > consensusCutoff)
      return(true);

    if ((FI->fragmentLength(olap.a_iid) == 0) ||
        (FI->fragmentLength(olap.b_iid) == 0))
      //  The overlap is bad if it involves deleted fragments.  Shouldn't happen in a normal
      //  assembly, but sometimes us users want to delete fragments after overlaps are generated.
      return(true);

    //  The overlap is GOOD (false == not bad) if the corrected error rate is below the requested
    //  erate.
    //
    if (olap.dat.ovl.corr_erate <= mismatchCutoff) {
      if (logFileFlags & LOG_OVERLAP_QUALITY)
        fprintf(logFile, "OVERLAP GOOD:     %d %d %c  hangs "F_S64" "F_S64" err %.3f %.3f\n",
                olap.a_iid, olap.b_iid,
                olap.dat.ovl.flipped ? 'A' : 'N',
                (int64)olap.dat.ovl.a_hang,
                (int64)olap.dat.ovl.b_hang,
                AS_OVS_decodeQuality(olap.dat.ovl.orig_erate),
                AS_OVS_decodeQuality(olap.dat.ovl.corr_erate));
      return(false);
    }

    //  If we didn't allow fixed-number-of-errors, the overlap is now bad.  Just a slight
    //  optimization.
    //
    if (mismatchLimit <= 0)
      return(true);

    //  There are a few cases where the orig_erate is _better_ than the corr_erate.  That is, the
    //  orig erate is 0% but we 'correct' it to something more than 0%.  Regardless, we probably
    //  want to be using the corrected erate here.

    double olen = olapLength(olap.a_iid, olap.b_iid, olap.dat.ovl.a_hang, olap.dat.ovl.b_hang);
    double nerr = olen * AS_OVS_decodeQuality(olap.dat.ovl.corr_erate);

    assert(nerr >= 0);

    if (nerr <= mismatchLimit) {
      if (logFileFlags & LOG_OVERLAP_QUALITY)
        fprintf(logFile, "OVERLAP SAVED:    %d %d %c  hangs "F_S64" "F_S64" err %.3f %.3f olen %f nerr %f\n",
                olap.a_iid, olap.b_iid,
                olap.dat.ovl.flipped ? 'A' : 'N',
                (int64)olap.dat.ovl.a_hang,
                (int64)olap.dat.ovl.b_hang,
                AS_OVS_decodeQuality(olap.dat.ovl.orig_erate),
                AS_OVS_decodeQuality(olap.dat.ovl.corr_erate),
                olen, nerr);
      return(false);
    }

    if (logFileFlags & LOG_OVERLAP_QUALITY)
      fprintf(logFile, "OVERLAP REJECTED: %d %d %c  hangs "F_S64" "F_S64" err %.3f %.3f olen %f nerr %f\n",
              olap.a_iid, olap.b_iid,
              olap.dat.ovl.flipped ? 'A' : 'N',
              (int64)olap.dat.ovl.a_hang,
              (int64)olap.dat.ovl.b_hang,
              AS_OVS_decodeQuality(olap.dat.ovl.orig_erate),
              AS_OVS_decodeQuality(olap.dat.ovl.corr_erate),
              olen, nerr);
    return(true);
  };

  uint64  scoreOverlap(const OVSoverlap& olap) {

    //  BPW's newer new score.  For the most part, we use the length of the overlap, but we also
    //  want to break ties with the higher quality overlap.
    //
    //  The high 20 bits are the length of the overlap.
    //  The next 12 are the corrected error rate.
    //  The last 12 are the original error rate.
    //
    //  (Well, 12 == AS_OVS_ERRBITS)

    if (isOverlapBadQuality(olap))
      return(0);

    uint64  leng = 0;
    uint64  corr = (AS_OVS_MAX_ERATE - olap.dat.ovl.corr_erate);
    uint64  orig = (AS_OVS_MAX_ERATE - olap.dat.ovl.orig_erate);

    //  Shift AFTER assigning to a 64-bit value to avoid overflows.
    corr <<= AS_OVS_ERRBITS;

    //  Containments - the length of the overlaps are all the same.  We return the quality.
    //
    if (((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0)) ||
        ((olap.dat.ovl.a_hang <= 0) && (olap.dat.ovl.b_hang >= 0)))
      return(corr | orig);

    //  Dovetails - the length of the overlap is the score, but we bias towards lower error.
    //  (again, shift AFTER assigning to avoid overflows)
    //
    leng   = olapLength(olap.a_iid, olap.b_iid, olap.dat.ovl.a_hang, olap.dat.ovl.b_hang);
    leng <<= (2 * AS_OVS_ERRBITS);

    return(leng | corr | orig);
  };


  uint32  olapLength(uint32 a_iid, uint32 b_iid, int32 a_hang, int32 b_hang) {
    int32  alen = FI->fragmentLength(a_iid);
    int32  blen = FI->fragmentLength(b_iid);
    int32  aovl = 0;
    int32  bovl = 0;

    assert(alen > 0);
    assert(blen > 0);

    if (a_hang < 0) {
      //  b_hang < 0      ?     ----------  :     ----
      //                  ?  ----------     :  ----------
      //
      aovl = (b_hang < 0) ? (alen + b_hang) : (alen);
      bovl = (b_hang < 0) ? (blen + a_hang) : (blen + a_hang - b_hang);
    } else {
      //  b_hang < 0      ?  ----------              :  ----------
      //                  ?     ----                 :     ----------
      //
      aovl = (b_hang < 0) ? (alen - a_hang + b_hang) : (alen - a_hang);
      bovl = (b_hang < 0) ? (blen)                   : (blen - b_hang);
    }

    if ((aovl <= 0) || (bovl <= 0) || (aovl > alen) || (bovl > blen)) {
      fprintf(stderr, "WARNING: bogus overlap found for A="F_U32" B="F_U32"\n", a_iid, b_iid);
      fprintf(stderr, "WARNING:                     A len="F_S32" hang="F_S32" ovl="F_S32"\n", alen, a_hang, aovl);
      fprintf(stderr, "WARNING:                     B len="F_S32" hang="F_S32" ovl="F_S32"\n", blen, b_hang, bovl);
    }

    if (aovl < 0)     aovl = 0;
    if (bovl < 0)     bovl = 0;

    if (aovl > alen)  aovl = alen;
    if (bovl > blen)  bovl = blen;

    assert(aovl > 0);
    assert(bovl > 0);
    assert(aovl <= alen);
    assert(bovl <= blen);

    //  AVE does not work.      return((uint32)((aovl, bovl)/2));
    //  MAX does not work.      return((uint32)MAX(aovl, bovl));

    return(aovl);
  };


private:
  void scoreContainment(const OVSoverlap& olap);
  void scoreEdge(const OVSoverlap& olap);

private:
  void save(const char *prefix, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT);
  bool load(const char *prefix, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT);

private:
  BestEdgeOverlap     *_best5;
  BestEdgeOverlap     *_best3;
  BestContainment     *_bestC;

  uint64              *_best5score;
  uint64              *_best3score;
  uint64              *_bestCscore;

public:
  uint64 mismatchCutoff;
  uint64 consensusCutoff;
  double mismatchLimit;
}; //BestOverlapGraph



#endif //INCLUDE_AS_BOG_BESTOVERLAPGRAPH

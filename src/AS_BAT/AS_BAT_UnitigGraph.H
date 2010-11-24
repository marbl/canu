
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

#ifndef INCLUDE_AS_BAT_UNITIGGRAPH
#define INCLUDE_AS_BAT_UNITIGGRAPH

static const char *rcsid_INCLUDE_AS_BAT_UNITIGGRAPH = "$Id: AS_BAT_UnitigGraph.H,v 1.1 2010-11-24 01:03:31 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"

#include "MultiAlignStore.h"

class MateLocation;

const uint32  noUnitig = 0xffffffff;


struct UnitigBreakPoint {
  FragmentEnd fragEnd;          // frag id and which end to break on
  SeqInterval fragPos;          // coordinates in unitig (used to get fwd/rev)

  int         fragsBefore;      // number of frags before the break fragment
  int         fragsAfter;       // number of frags after the break fragment

  FragmentEnd inEnd;            // frag id of the incoming unitig

  int         inSize;           // the size of the incoming unitig
  int         inFrags;          // the number of fragments in incoming unitig

  UnitigBreakPoint() {
    fragEnd      = FragmentEnd();
    fragPos.bgn  = 0;
    fragPos.end  = 0;

    fragsBefore  = 0;
    fragsAfter   = 0;

    inEnd        = FragmentEnd();
    inSize       = 0;
    inFrags      = 0;
  };

  UnitigBreakPoint(uint32 isectFrg, bool isect3p, SeqInterval isectPos,
                   uint32 frgsBefore, uint32 frgsAfter,
                   uint32 invaderFrg, bool invader3p,
                   uint32 invaderTigLen,
                   uint32 invaderTigFrags) {
    fragEnd      = FragmentEnd(isectFrg, isect3p);
    fragPos      = isectPos;

    fragsBefore  = frgsBefore;
    fragsAfter   = frgsAfter;

    inEnd        = FragmentEnd(invaderFrg, invader3p);
    inSize       = invaderTigLen;
    inFrags      = invaderTigFrags;
  };
};

typedef std::vector<UnitigBreakPoint> UnitigBreakPoints;




class overlapPlacement {
public:
  overlapPlacement() {
    frgID                 = 0;

    tigID                 = 0;
    position.bgn          = 0;
    position.end          = 0;

    nReverse              = 0;
    nForward              = 0;

    fCoverage             = 0.0;

    errors                = 0.0;
    aligned               = 0;
    covered.bgn           = 0;
    covered.end           = 0;
  };
  ~overlapPlacement() {
  };

  bool operator<(const overlapPlacement &that) const {
    if (tigID != that.tigID)
      return(tigID < that.tigID);
    return(position < that.position);
  };

public:
  uint32            frgID;      //  ID of the fragment this position is for.

  uint32            tigID;      //  Unitig ID of this placement
  SeqInterval       position;   //  Unitig position of this placement

  uint32            nReverse;   //  Num overlaps aligned reverse
  uint32            nForward;   //  Num overlaps aligned forward

  double            fCoverage;  //  Coverage of the fragment

  double            errors;     //  number of errors in alignments
  uint32            aligned;    //  number of bases in alignments
  SeqInterval       covered;    //  Position of the overlap on the read
};


class joinEntry;


struct UnitigGraph{   
  // This will store the entire set of unitigs that are generated
  // It's just a unitig container.
  UnitigGraph();
  ~UnitigGraph();

  // Call this on a chunk graph pointer to build a unitig graph
  void build(OverlapStore *ovlStoreUniq,
             OverlapStore *ovlStoreRept,
             bool enableIntersectionBreaking,
             bool enableJoining,
             bool enableBubblePopping,
             int32 badMateBreakThreshold,
             char *output_prefix);
  void setParentAndHang();

  void unitigToMA(MultiAlignT *ma, uint32 iumiid, Unitig *utg);

  void writeIUMtoFile(char *fileprefix, char *tigStorePath, uint32 fragment_count_target, bool isFinal=true);
  void writeOVLtoFile(char *fileprefix);
  void writeCGAtoFile(char *fileprefix, float globalARate);

  float getGlobalArrivalRate(long total_random_frags_in_genome=0, long genome_size=0);

  void breakUnitigs(char *output_prefix, bool enableIntersectionBreaking);

  bool findEdgeToUnitig(ufNode *aFrg, BestEdgeOverlap &a5, BestEdgeOverlap &a3,
                        Unitig *tUtg);
  bool findEdges(ufNode *aFrg, BestEdgeOverlap &a5, BestEdgeOverlap &a3,
                 ufNode *bFrg, BestEdgeOverlap &b5, BestEdgeOverlap &b3);

  void joinUnitigs_examineFirst(Unitig *fr, vector<joinEntry> &joins);
  void joinUnitigs_examineLast(Unitig *fr, vector<joinEntry> &joins);
  bool joinUnitigs_confirmJoin(joinEntry *join);
  void joinUnitigs_merge(joinEntry *join);
  void joinUnitigs_append(joinEntry *join);
  void joinUnitigs(bool enableJoining);

  void placeContainsUsingBestOverlaps(void);
  void placeContainsUsingAllOverlaps(OverlapStore *ovlStoreUniq,
                                     OverlapStore *ovlStoreRept,
                                     bool   withMatesToNonContained,
                                     bool   withMatesToUnambiguousContain);

  void placeZombies(void);

  //
  bool validateBubbleFragmentsWithOverlaps(Unitig *bubble,
                                           ufNode &lFrg,
                                           ufNode &rFrg,
                                           Unitig *larger,
                                           OverlapStore *ovlStoreUniq,
                                           OverlapStore *ovlStoreRept);
  bool validateBubbleWithEdges(Unitig *bubble,
                               ufNode &lFrg, BestEdgeOverlap *lEnd,
                               ufNode &rFrg, BestEdgeOverlap *rEnd,
                               Unitig *larger,
                               OverlapStore *ovlStoreUniq,
                               OverlapStore *ovlStoreRept);

  void popIntersectionBubbles(OverlapStore *ovlStoreUniq,
                              OverlapStore *ovlStoreRept);
  void popMateBubbles(OverlapStore *ovlStoreUniq,
                      OverlapStore *ovlStoreRept);
  void popOverlapBubbles(OverlapStore *ovlStoreUniq,
                         OverlapStore *ovlStoreRept);
  //

  bool placeFragUsingOverlaps(uint32 fid,
                              OverlapStore *ovlStoreUniq,
                              OverlapStore *ovlStoreRept,
                              vector<overlapPlacement> &placements);

  void filterBreakPoints(Unitig *,
                         UnitigBreakPoints &);

  UnitigBreakPoint selectSmall(const Unitig *tig,
                               const UnitigBreakPoints &smalls,
                               const UnitigBreakPoint &big,
                               int   &lastBPCoord,
                               int   &lastBPFragNum);

  UnitigVector*  breakUnitigAt(Unitig *, UnitigBreakPoints &);

  void           checkUnitigMembership(void);
  void           reportOverlapsUsed(const char *prefix, const char *name);
  void           reportUnitigs(const char *prefix, const char *name);

  UnitigVector   unitigs;

public:
  //  This is the MateChecker interface
  UnitigBreakPoints* computeMateCoverage(Unitig *utg, MateLocation *mates, int badMateBreakThreshold);

  void evaluateMates(const char *output_prefix, const char *name);

  void moveContains(void);
  void splitDiscontinuousUnitigs(void);


private:
  static const int MIN_BREAK_FRAGS = 1;
  static const int MIN_BREAK_LENGTH = 500;

  void populateUnitig(int32               fragID);

  void populateUnitig(Unitig             *unitig,
                      BestEdgeOverlap    *nextedge);
};

#endif

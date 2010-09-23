
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

#ifndef INCLUDE_AS_BOG_UNITIGGRAPH
#define INCLUDE_AS_BOG_UNITIGGRAPH

static const char *rcsid_INCLUDE_AS_BOG_UNITIGGRAPH = "$Id: AS_BOG_UnitigGraph.hh,v 1.73 2010-09-23 06:13:11 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_Unitig.hh"

#include "MultiAlignStore.h"

typedef std::vector<uint32>                 FragmentList;
typedef std::map<uint32, FragmentList>      FragmentEdgeList;

typedef std::map<uint32, SeqInterval>       FragmentPositionMap;

inline bool isReverse( SeqInterval pos ) {
  return(pos.bgn > pos.end);
}

struct UnitigBreakPoint {
  FragmentEnd fragEnd;          // frag id and which end to break on
  SeqInterval fragPos;          // coordinates in unitig (used to get fwd/rev)

  //  Number of fragments before and after the fragment we break on.
  //  "Before" always includes the frag we break on and any contains
  //  in it, regardless of end.  Please fix that.
  int         fragsBefore;
  int         fragsAfter;

  FragmentEnd inEnd;            // frag id of the incoming unitig

  int         inSize;           // the size of the incoming unitig
  int         inFrags;          // the number of fragments in incoming unitig

  UnitigBreakPoint(uint32 id=0, uint32 end=FIVE_PRIME) {
    fragEnd      = FragmentEnd(id, end);
    fragPos.bgn  = 0;
    fragPos.end  = 0;

    fragsBefore  = 0;
    fragsAfter   = 0;

    inEnd        = FragmentEnd();
    inSize       = 0;
    inFrags      = 0;
  }
};

//  Unfortunately, this does need to be a list.  Search for pop_front.
//
typedef std::list<UnitigBreakPoint> UnitigBreakPoints;



struct UnitigGraph{   
  // This will store the entire set of unitigs that are generated
  // It's just a unitig container.
  UnitigGraph(FragmentInfo *fi, BestOverlapGraph *);
  ~UnitigGraph();

  // Call this on a chunk graph pointer to build a unitig graph
  void build(ChunkGraph *cg_ptr,
             OverlapStore *ovlStoreUniq,
             OverlapStore *ovlStoreRept,
             bool enableIntersectionBreaking,
             bool enableJoining,
             bool enableBubblePopping,
             char *output_prefix);
  void setParentAndHang(ChunkGraph *cg_ptr);

  void unitigToMA(MultiAlignT *ma, uint32 iumiid, Unitig *utg);

  void writeIUMtoFile(char *fileprefix, char *tigStorePath, int fragment_count_target, bool isFinal=true);
  void writeOVLtoFile(char *fileprefix);

  float getGlobalArrivalRate(long total_random_frags_in_genome=0, long genome_size=0);

  void breakUnitigs(ContainerMap &cMap, char *output_prefix, bool enableIntersectionBreaking);
  void joinUnitigs(bool enableJoining);
  void placeContains(void);
  void placeZombies(void);
  void popIntersectionBubbles(OverlapStore *ovlStoreUniq,
                              OverlapStore *ovlStoreRept);
  void popMateBubbles(OverlapStore *ovlStoreUniq,
                      OverlapStore *ovlStoreRept);

  void filterBreakPoints(ContainerMap &cMap,
                         Unitig *,
                         UnitigBreakPoints &);

  UnitigBreakPoint selectSmall(ContainerMap &cMap,
                               const Unitig *tig,
                               const UnitigBreakPoints &smalls,
                               const UnitigBreakPoint &big,
                               int   &lastBPCoord,
                               int   &lastBPFragNum);

  UnitigVector*  breakUnitigAt(ContainerMap &cMap, Unitig *, UnitigBreakPoints &);

  void           checkUnitigMembership(void);
  void           reportOverlapsUsed(const char *filename);
  void           reportUnitigs(const char *filename);

  // Unitigs are the dove tails and their contained fragments
  UnitigVector *unitigs;

  BestOverlapGraph *bog_ptr;

private:
  static const int MIN_BREAK_FRAGS = 1;
  static const int MIN_BREAK_LENGTH = 500;

  void populateUnitig(int32               fragID);

  void populateUnitig(Unitig             *unitig,
                      BestEdgeOverlap    *nextedge);

  FragmentInfo     *_fi;

  //  This is a map from 'invaded fragment' to a list of 'invading fragments';
  //    unitigIntersect[a] = b means that b is invading into a.
  //
  //  unitigIntersect is a FragmentEdgeList, a map<uint32, FragmentList>.
  //                       FragmentList,     a vector<uint32>
  //
  //  selfIntersect is a list of the invading fragments that are known to be intersecting their own
  //  unitig.  When short unitigs are merged into larger unitigs, those short unitigs likely have
  //  saved a unitigIntersect to the exact spot where we merged them in, and we need to NOT split on
  //  those false intersections.  On the other hand, if an invading fragment is listed here, we know
  //  it is a true self-intersection, and allow it.
  //
  FragmentEdgeList      unitigIntersect;
  map<uint32,bool>      selfIntersect;
};

#endif

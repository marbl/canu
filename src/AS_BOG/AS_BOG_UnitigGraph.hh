
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

#include <set>
#include <iostream>
#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_Unitig.hh"

typedef std::vector<iuid>                 FragmentList;
typedef std::map<iuid, FragmentList>      FragmentEdgeList;

typedef std::map<iuid, SeqInterval>       FragmentPositionMap;

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

  UnitigBreakPoint(iuid id=0, fragment_end_type end=FIVE_PRIME) {
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




struct BestEdgeCounts{
  int oneWayBest;
  int dovetail;
  int neither;
  int contained;

  BestEdgeCounts() : oneWayBest(0), dovetail(0), neither(0), contained(0) {}
  BestEdgeCounts operator+=(BestEdgeCounts other){
    oneWayBest += other.oneWayBest;
    dovetail   += other.dovetail;
    neither    += other.neither;
    contained  += other.contained;
  } 
};



struct UnitigGraph{
  // This will store the entire set of unitigs that are generated
  // It's just a unitig container.
  UnitigGraph(FragmentInfo *fi, BestOverlapGraph *);
  ~UnitigGraph();

  // Call this on a chunk graph pointer to build a unitig graph
  void build(ChunkGraph *cg_ptr);

  void writeIUMtoFile(char *filename, int fragment_count_target);
#if 0
  void readIUMsFromFile(const char *filename, iuid maxIID);
#endif

  float getGlobalArrivalRate(long total_random_frags_in_genome=0, long genome_size=0);

  void breakUnitigs();
  void printUnitigBreaks();

  void filterBreakPoints( Unitig *, UnitigBreakPoints &);

  UnitigBreakPoint selectSmall(const Unitig *tig,
                               const UnitigBreakPoints &smalls,
                               const UnitigBreakPoint &big,
                               int   &lastBPCoord,
                               int   &lastBPFragNum);

  UnitigVector* breakUnitigAt( Unitig *, UnitigBreakPoints &);

  // Counts status of best edges internal to a unitig
  BestEdgeCounts countInternalBestEdges( const Unitig *);
  BestEdgeCounts countInternalBestEdges( ); // all unitigs

  void           checkUnitigMembership(void);

  // Unitigs are the dove tails and their contained fragments
  UnitigVector *unitigs;

  BestOverlapGraph *bog_ptr;

private:
  // Given a fragment, it will follow it's overlaps until 
  //   the end, and add them to the unitig
  void populateUnitig(
                      Unitig* unitig, 
                      iuid src_frag_id, 
                      fragment_end_type whichEnd,
                      ChunkGraph *cg_ptr,
                      int offset,
                      bool verbose);

  // Inverts the containment map to key by container, instead of containee
  void _build_container_map(BestContainmentMap*);

  // Build containee list
  ContainerMap *_extract_containees(DoveTailPath *dtp_ptr, 
                                    ContainerMap *cntnrmap_ptr);

  // Compute the global arrival rate based on the unitig rho's.
  float _compute_global_arrival_rate(void);

  FragmentInfo     *_fi;

  FragmentEdgeList  unitigIntersect;
  ContainerMap     *cntnrmap_ptr;
};
		
#endif

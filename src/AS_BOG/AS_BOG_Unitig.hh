
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

#ifndef INCLUDE_AS_BOG_UNITIG
#define INCLUDE_AS_BOG_UNITIG

#include "AS_BOG_Datatypes.hh"
//#include "AS_BOG_ChunkGraph.hh"

#include <vector>
#include <map>
#include <set>


typedef IntMultiPos                       DoveTailNode;
typedef std::vector<DoveTailNode>         DoveTailPath;


typedef DoveTailPath::iterator            DoveTailIter;
typedef DoveTailPath::const_iterator      DoveTailConstIter;


typedef std::vector<iuid>                 ContaineeList;
typedef std::map<iuid, ContaineeList>     ContainerMap;


struct BestOverlapGraph;


struct Unitig{
  Unitig(bool report=false);
  ~Unitig(void);

  // Sort frags by position on the unitig
  void sort();

  // Compute unitig based on given dovetails and containments
  void recomputeFragmentPositions(ContainerMap &, BestOverlapGraph*);
  void computeFragmentPositions(FragmentInfo*, BestOverlapGraph*);

  void shiftCoordinates(int);
  void reverseComplement();
  void reverseComplement(int offset, BestOverlapGraph *);

  // Accessor methods
  float getAvgRho(FragmentInfo *fi);
  static void setGlobalArrivalRate(float global_arrival_rate);
  void setLocalArrivalRate(float local_arrival_rate);
  float getLocalArrivalRate(FragmentInfo *fi);
  float getCovStat(FragmentInfo *fi);
  long getLength(void);
  long getNumFrags(void);
  long getNumRandomFrags(void); // For now, same as numFrags, but should be randomly sampled frag count
  DoveTailNode getLastBackboneNode(iuid&);

  iuid         id(void) { return(_id); };

  void addContainedFrag(DoveTailNode, BestContainment *bestcont, bool report=false);
  void addFrag(DoveTailNode, int offset=0, bool report=false);

  static iuid fragIn(iuid fragId) {
    if (_inUnitig == NULL)
      return 0;
    return _inUnitig[fragId];
  };

  static void resetFragUnitigMap(iuid numFrags) {
    if (_inUnitig == NULL)
      _inUnitig = new iuid[numFrags+1];
    memset(_inUnitig, 0, (numFrags+1) * sizeof(iuid));
  };

  // Public Member Variables
  DoveTailPath *dovetail_path_ptr;

private:
  void placeContains(const ContainerMap &,
                     BestOverlapGraph *,
                     const iuid,
                     const SeqInterval,
                     const int level);

  // Do not access these private variables directly, they may
  // not be computed yet, use accessors!
  //
  float  _avgRho;
  float  _covStat;
  long   _length;
  float  _localArrivalRate;
  iuid   _id;

  static iuid   nextId;
  static float _globalArrivalRate;
  static iuid *_inUnitig;
};


typedef std::vector<Unitig*> UnitigVector;
typedef UnitigVector::iterator UnitigsIter;
typedef UnitigVector::const_iterator UnitigsConstIter;






#endif  //  INCLUDE_AS_BOG_UNITIG

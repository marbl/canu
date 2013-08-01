
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

static const char *rcsid_INCLUDE_AS_BOG_UNITIG = "$Id$";

#include "AS_BOG_Datatypes.H"

//  Derived from IntMultiPos, but removes some of the data (48b in IntMultiPos, 32b in struct
//  ufNode).  The minimum size (bit fields, assuming maximum limits, not using the contained
//  field) seems to be 24b, and is more effort than it is worth (just removing 'contained' would be
//  a chore).
//
//  ufNode is, of course, 'unitig fragment node'.
//
struct ufNode {
  uint32           ident;
  uint32           contained;
  uint32           parent;     //  IID of the fragment we align to

  int32           ahang;       //  If parent defined, these are relative
  int32           bhang;       //  that fragment

  SeqInterval     position;

  uint32           containment_depth;
};

typedef std::vector<ufNode>         ufPath;



//  The ContainerMap stores a list of IDs contained in the index fragment:
//    cMap[index] = list of IDs contained in index
//
//  It is used during splitting to get the total number of fragments in a unitig, the number of
//  fragments before/after the split point.
//
typedef std::vector<uint32>                 ContaineeList;
typedef std::map<uint32, ContaineeList>     ContainerMap;


struct BestOverlapGraph;


struct Unitig{
  Unitig(bool report=false);
  ~Unitig(void);

  void sort(void);
  void reverseComplement(bool doSort=true);

  // getNumRandomFrags() is a placeholder, random frags should not
  // contain guides, or other fragments that are not randomly sampled
  // across the whole genome.

  int32  getLength(void)          { return(_length);       };
  uint32 getNumFrags(void)        { return(ufpath.size()); };
  uint32 getNumRandomFrags(void)  { return(getNumFrags()); };

  ufNode getLastBackboneNode(void);
  ufNode getLastBackboneNode(uint32 &);

  uint32       id(void) { return(_id); };

  bool placeFrag(ufNode &place5, int32 &fidx5, BestEdgeOverlap *bestedge5,
                 ufNode &place3, int32 &fidx3, BestEdgeOverlap *bestedge3);

  bool placeFrag(ufNode &frag, BestContainment *bestcont);

  void addFrag(ufNode node, int offset=0, bool report=false);
  bool addContainedFrag(int32 fid, BestContainment *bestcont, bool report=false);
  bool addAndPlaceFrag(int32 fid, BestEdgeOverlap *bestedge5, BestEdgeOverlap *bestedge3, bool report=false);


  static uint32 fragIn(uint32 fragId) {
    if ((_inUnitig == NULL) || (fragId == 0))
      return 0;
    return _inUnitig[fragId];
  };

  static uint32 pathPosition(uint32 fragId) {
    if ((_pathPosition == NULL) || (fragId == 0))
      return ~0;
    return _pathPosition[fragId];
  };

  static void resetFragUnitigMap(uint32 numFrags) {
    if (_inUnitig == NULL)
      _inUnitig = new uint32[numFrags+1];
    memset(_inUnitig, 0, (numFrags+1) * sizeof(uint32));

    if (_pathPosition == NULL)
      _pathPosition = new uint32[numFrags+1];
    memset(_pathPosition, 0, (numFrags+1) * sizeof(uint32));
  };

  // Public Member Variables
  ufPath   ufpath;

private:
  int32    _length;
  uint32   _id;

  static uint32  _nextId;
  static uint32 *_inUnitig;      //  Maps a fragment iid to a unitig id.
  static uint32 *_pathPosition;  //  Maps a fragment iid to an index in the dovetail path
};


typedef std::vector<Unitig*> UnitigVector;
typedef UnitigVector::iterator UnitigsIter;
typedef UnitigVector::const_iterator UnitigsConstIter;


#endif  //  INCLUDE_AS_BOG_UNITIG

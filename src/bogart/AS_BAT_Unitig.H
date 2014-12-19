
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

#ifndef INCLUDE_AS_BAT_UNITIG
#define INCLUDE_AS_BAT_UNITIG

static const char *rcsid_INCLUDE_AS_BAT_UNITIG = "$Id$";

#include "AS_BAT_Datatypes.H"

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

  int32            ahang;       //  If parent defined, these are relative
  int32            bhang;       //  that fragment

  SeqInterval      position;

  uint32           containment_depth;
};



struct Unitig {
private:
  Unitig() {
    _length = 0;
    _id     = 0;
  };
    
public:
  ~Unitig(void) {
  };

  friend class UnitigVector;

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

  void placeFrag_computePlacement(ufNode          &frag,
                                  int32           &bidx,
                                  BestEdgeOverlap *bestedge,
                                  bool             bestIs3);

  bool placeFrag(ufNode &place5, int32 &fidx5, BestEdgeOverlap *bestedge5,
                 ufNode &place3, int32 &fidx3, BestEdgeOverlap *bestedge3);

  bool placeFrag(ufNode &frag, BestContainment *bestcont);

  void addFrag(ufNode node, int offset=0, bool report=false);
  bool addContainedFrag(int32 fid, BestContainment *bestcont, bool report=false);
  bool addAndPlaceFrag(int32 fid, BestEdgeOverlap *bestedge5, BestEdgeOverlap *bestedge3, bool report=false);

  void bubbleSortLastFrag(void);

  static void removeFrag(int32 fid) {
    _inUnitig[fid] = 0;
    _pathPosition[fid] = ~0;
  };

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
  vector<ufNode>         ufpath;

private:
  int32    _length;
  uint32   _id;

  static uint32 *_inUnitig;      //  Maps a fragment iid to a unitig id.
  static uint32 *_pathPosition;  //  Maps a fragment iid to an index in the dovetail path
};



class UnitigVector {
public:
  UnitigVector() {
    _blockSize    = 1048576;
    _numBlocks    = 1;
    _maxBlocks    = 1024;
    _blocks       = new Unitig ** [_maxBlocks];
    _blocks[0]    = new Unitig  * [_blockSize];
    _blocks[0][0] = NULL;  //  No first unitig.
    _blockNext    = 1;
    _totalUnitigs = 1;
  };
  ~UnitigVector() {
  };

  Unitig *newUnitig(bool verbose) {
    Unitig *u = new Unitig();

#pragma omp critical
    {
      u->_id = _totalUnitigs++;

      if (verbose)
        writeLog("Creating Unitig %d\n", u->_id);

      if (_blockNext >= _blockSize) {
        assert(_numBlocks < _maxBlocks);

        _blocks[_numBlocks] = new Unitig * [_blockSize];

        memset(_blocks[_numBlocks], 0, sizeof(Unitig **) * _blockSize);

        _numBlocks++;
        _blockNext = 0;
      }

      _blocks[_numBlocks-1][_blockNext++] = u;

      //  The rest are just sanity checks.

      assert((u->id() / _blockSize) == (_numBlocks - 1));
      assert((u->id() % _blockSize) == (_blockNext - 1));

      assert(operator[](u->id()) == u);
    }

    return(u);
  };

  size_t  size(void) {
    return(_totalUnitigs);
  };

  Unitig *&operator[](uint32 i) {
    uint32  idx = i / _blockSize;
    uint32  pos = i % _blockSize;

#ifdef CHECK_UNITIG_ARRAY_INDEXING
    if (((i    >= _totalUnitigs)) ||
        ((idx  >= _numBlocks)) ||
        (((pos >= _blockNext) && (idx >= _numBlocks - 1)))) {
      fprintf(stderr, "UnitigVector::operator[]()--  i="F_U32" with totalUnitigs="F_U64"\n", i, _totalUnitigs);
      fprintf(stderr, "UnitigVector::operator[]()--  blockSize="F_U64"\n", _blockSize);
      fprintf(stderr, "UnitigVector::operator[]()--  idx="F_U32" numBlocks="F_U64"\n", idx, _numBlocks);
      fprintf(stderr, "UnitigVector::operator[]()--  pos="F_U32" blockNext="F_U64"\n", pos, _blockNext);
    }
    assert(i    < _totalUnitigs);
    assert((idx < _numBlocks));
    assert((pos < _blockNext) || (idx < _numBlocks - 1));
#endif

    return(_blocks[idx][pos]);
  };

private:
  uint64              _blockSize;

  uint64              _numBlocks;
  uint64              _maxBlocks;
  Unitig           ***_blocks;
  uint64              _blockNext;

  uint64              _totalUnitigs;
};


#endif  //  INCLUDE_AS_BAT_UNITIG

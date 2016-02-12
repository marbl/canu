
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_Unitig.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_UNITIG
#define INCLUDE_AS_BAT_UNITIG

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
    _length        = 0;
    _id            = 0;
    _tigID         = 0;
    _isUnassembled = false;
    _isBubble      = false;
    _isRepeat      = false;
    _isCircular    = false;
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

  uint32 id(void)                 { return(_id); };      //  ID internal to bogart
  uint32 tigID(void)              { return(_tigID); };   //  ID in the output store

  int32  getLength(void)          { return(_length);       };
  uint32 getNumFrags(void)        { return(ufpath.size()); };
  uint32 getNumRandomFrags(void)  { return(getNumFrags()); };

  ufNode getLastBackboneNode(void);
  ufNode getLastBackboneNode(uint32 &);

  void   placeFrag_computePlacement(ufNode          &frag,
                                    int32           &bidx,
                                    BestEdgeOverlap *bestedge,
                                    bool             bestIs3);

  bool   placeFrag(ufNode &place5, int32 &fidx5, BestEdgeOverlap *bestedge5,
                   ufNode &place3, int32 &fidx3, BestEdgeOverlap *bestedge3);

  bool   placeFrag(ufNode &frag, BestContainment *bestcont);

  void   addFrag(ufNode node, int offset=0, bool report=false);
  bool   addContainedFrag(int32 fid, BestContainment *bestcont, bool report=false);
  bool   addAndPlaceFrag(int32 fid, BestEdgeOverlap *bestedge5, BestEdgeOverlap *bestedge3, bool report=false);

  void   bubbleSortLastFrag(void);


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
  int32          _length;
  uint32         _id;
public:
  uint32         _tigID;

private:
  static uint32 *_inUnitig;      //  Maps a fragment iid to a unitig id.
  static uint32 *_pathPosition;  //  Maps a fragment iid to an index in the dovetail path

public:
  //  Classification.  The output is in three files: 'unassembled', 'bubbles', 'contigs' (defined as
  //  not unassembled and not bubble).

  uint32         _isUnassembled;  //  Is a single read or a pseudo singleton
  uint32         _isBubble;       //  Annotation: from a failed bubble pop
  uint32         _isRepeat;       //  Annotation: from an identified repeat region
  uint32         _isCircular;     //  Annotation: has overlap to self
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

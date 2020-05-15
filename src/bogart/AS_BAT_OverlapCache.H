
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef INCLUDE_AS_BAT_OVERLAPCACHE
#define INCLUDE_AS_BAT_OVERLAPCACHE

#include "runtime.H"
#include "ovStore.H"
#include "sqStore.H"

//  CA8 used to re-encode the error rate into a smaller-precision number.  This was
//  confusing and broken (it tried to use a log-based encoding to give more precision
//  to the smaller values).  CA3g gives up and uses all 12 bits of precision.

//  If not enough space for the minimum number of error bits, bump up to a 64-bit word for overlap
//  storage.

//  For storing overlaps in memory.  12 bytes per overlap.
class BAToverlap {
public:
  BAToverlap() {
    evalue    = 0;
    a_hang    = 0;
    b_hang    = 0;
    flipped   = false;

    filtered  = false;
    symmetric = false;

    a_iid     = 0;
    b_iid     = 0;
  };
  ~BAToverlap() {};

  //  Nasty bit of code duplication.

  bool
  isDovetail(void) const {
    return(((a_hang < 0) && (b_hang < 0)) ||
           ((a_hang > 0) && (b_hang > 0)));
  };

  bool
  AEndIs5prime(void) const {                   //  -------->
    return((a_hang < 0) && (b_hang < 0));      //        -------
  };

  bool
  AEndIs3prime(void) const {                   //     -------->
    return((a_hang > 0) && (b_hang > 0));      //  -------
  };

  bool
  AisContainer(void) const {                   //  -------->
    return((a_hang >= 0) && (b_hang <= 0));    //    ----
  };

  bool
  AisContained(void) const {                   //    --->
    return((a_hang <= 0) && (b_hang >= 0));    //  ---------
  };

  bool
  BEndIs3prime(void) const {
    assert(AisContainer() == false);  //  Function is not defined
    assert(AisContained() == false);  //    for containments.
    return((AEndIs5prime() && (flipped == false)) ||   // <===     ------>
           (AEndIs3prime() && (flipped == true)));     //        ---->
  };

  bool
  BEndIs5prime(void) const {
    assert(AisContainer() == false);  //  Function is not defined
    assert(AisContained() == false);  //    for containments.
    return((AEndIs5prime() && (flipped == true)) ||   //          ------>
           (AEndIs3prime() && (flipped == false)));   // <===          ---->
  };

  double
  erate(void) const {
    return(AS_OVS_decodeEvalue(evalue));
  }

  void
  convert(ovOverlap &olap) {
    olap.clear();

    olap.a_iid = a_iid;
    olap.b_iid = b_iid;

    olap.evalue(evalue);

    olap.flipped(flipped);
    olap.a_hang(a_hang);
    olap.b_hang(b_hang);
  }

#if AS_MAX_READLEN_BITS < 24
  uint64      evalue    : AS_MAX_EVALUE_BITS;     //  12
  int64       a_hang    : AS_MAX_READLEN_BITS+1;  //  21+1
  int64       b_hang    : AS_MAX_READLEN_BITS+1;  //  21+1
  uint64      flipped   : 1;                      //   1

  uint64      filtered  : 1;                      //   1
  uint64      symmetric : 1;                      //   1    - twin overlap exists

  uint32      a_iid;
  uint32      b_iid;

#if (AS_MAX_EVALUE_BITS + (AS_MAX_READLEN_BITS + 1) + (AS_MAX_READLEN_BITS + 1) + 1 + 1 + 1 > 64)
#error not enough bits to store overlaps.  decrease AS_MAX_EVALUE_BITS or AS_MAX_READLEN_BITS.
#endif

#else
  int32       a_hang;
  int32       b_hang;

  uint32      evalue    : AS_MAX_EVALUE_BITS;     //  12
  uint32      flipped   : 1;                      //   1
  uint32      filtered  : 1;                      //   1
  uint32      symmetric : 1;                      //   1    - twin overlap exists

  uint32      a_iid;
  uint32      b_iid;
#endif

};



inline
bool
BAToverlap_sortByEvalue(BAToverlap const &a, BAToverlap const &b) {
  return(a.evalue > b.evalue);
}



class OverlapStorage {
public:
  OverlapStorage(uint64 nOvl) {
    _osAllocLen = 1024 * 1024 * 1024 / sizeof(BAToverlap);  //  1GB worth of overlaps
    _osLen      = 0;                            //  osMax is cheap and we overallocate it.
    _osPos      = 0;                            //  If allocLen is small, we can end up with
    _osMax      = 2 * nOvl / _osAllocLen + 2;   //  more blocks than expected, when overlaps
    _os         = new BAToverlap * [_osMax];    //  don't fit in the remaining space.

    memset(_os, 0, sizeof(BAToverlap *) * _osMax);

    _os[0]      = new BAToverlap [_osAllocLen];   //  Alloc first block, keeps getOverlapStorage() simple
  };

  OverlapStorage(OverlapStorage *original) {
    _osAllocLen = original->_osAllocLen;
    _osLen      = 0;
    _osPos      = 0;
    _osMax      = original->_osMax;
    _os         = NULL;
  };

  ~OverlapStorage() {
    if (_os == NULL)
      return;

    for (uint32 ii=0; ii<_osMax; ii++)
      delete [] _os[ii];
    delete [] _os;
  }


  OverlapStorage *reset(void) {
    _osLen = 0;
    _osPos = 0;

    return(this);
  };


  BAToverlap   *get(void) {
    if (_os == NULL)
      return(NULL);
    return(_os[_osLen] + _osPos);
  };


  BAToverlap   *get(uint32 nOlaps) {
    if (_osPos + nOlaps > _osAllocLen) {           //  If we don't fit in the current allocation,
      _osPos = 0;                                  //  move to the next one.
      _osLen++;
    }

    _osPos += nOlaps;                              //  Reserve space for these overlaps.

    assert(_osLen < _osMax);

    if (_os == NULL)                               //  If we're not allowed to allocate,
      return(NULL);                                //  return nothing.

    if (_os[_osLen] == NULL)                       //  Otherwise, make sure we have space and return
      _os[_osLen] = new BAToverlap [_osAllocLen];  //  that space.

    return(_os[_osLen] + _osPos - nOlaps);
  };


  void          advance(OverlapStorage *that) {
    if (((that->_osLen <  _osLen)) ||                            //  That segment before mine, or
        ((that->_osLen == _osLen) && (that->_osPos <= _osPos)))  //  that segment equal and position before mine
      return;                                                    //  So no need to modify

    _osLen = that->_osLen;
    _osPos = that->_osPos;
  };


private:
  uint32                  _osAllocLen;   //  Size of each allocation
  uint32                  _osLen;        //  Current allocation being used
  uint32                  _osPos;        //  Position in current allocation; next free overlap
  uint32                  _osMax;        //  Number of allocations we can make
  BAToverlap            **_os;           //  Allocations
};



class OverlapCache {
public:
  OverlapCache(const char *ovlStorePath,
               const char *prefix,
               double maxErate,
               uint32 minOverlap,
               uint64 maxMemory,
               uint64 genomeSize);
  ~OverlapCache();

  bool         compareOverlaps(const BAToverlap &a, const BAToverlap &b) const; // we can almost do templated but the fields are functions in one and just members in the other

private:
  bool         compareOverlaps(const ovOverlap &a,  const ovOverlap &b) const;

  uint32       filterOverlaps(uint32 aid, uint32 maxOVSerate, uint32 minOverlap, uint32 no);
  uint32       filterDuplicates(uint32 &no);

  void         computeOverlapLimit(ovStore *ovlStore, uint64 genomeSize);
  void         loadOverlaps(ovStore *ovlStore);
  void         symmetrizeOverlaps(void);

public:
  BAToverlap  *getOverlaps(uint32 readIID, uint32 &numOverlaps) {
    numOverlaps = _overlapLen[readIID];
    return(_overlaps[readIID]);
  }

private:
  const char             *_prefix;

  uint64                  _memLimit;       //  Expected max size of bogart
  uint64                  _memReserved;    //  Memory to reserve for processing
  uint64                  _memAvail;       //  Memory available for storing overlaps
  uint64                  _memStore;       //  Memory used to support overlaps
  uint64                  _memOlaps;       //  Memory used to store overlaps

  uint32                 *_overlapLen;
  uint32                 *_overlapMax;
  BAToverlap            **_overlaps;

  //  Instead of allocating space for overlaps per read (which has some visible but unknown size
  //  cost with each allocation), or in a single massive allocation (which we can't resize), we
  //  allocate overlaps in large blocks then set pointers into each block where overlaps for each
  //  read start.  This is managed by OverlapStorage.

  OverlapStorage         *_overlapStorage;

  uint32                  _maxEvalue;  //  Don't load overlaps with high error
  uint32                  _minOverlap; //  Don't load overlaps that are short

  uint32                  _minPer;     //  Minimum number of overlaps to retain for a single read
  uint32                  _maxPer;     //  Maximum number of overlaps to load for a single read

  uint64                 *_minSco;     //  The minimum score accepted for each read

  bool                    _checkSymmetry;

  uint32                  _ovsMax;     //  For loading overlaps
  ovOverlap              *_ovs;        //
  uint64                 *_ovsSco;     //  For scoring overlaps during the load
  uint64                 *_ovsTmp;     //  For picking out a score threshold

  uint64                  _genomeSize;
};



extern OverlapCache     *OC;

#endif  //  INCLUDE_AS_BAT_OVERLAPCACHE

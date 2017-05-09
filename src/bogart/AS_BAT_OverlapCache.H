
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
 *    src/AS_BAT/AS_BAT_OverlapCache.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2011-FEB-15 to 2013-AUG-28
 *      are Copyright 2011-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2011-OCT-13
 *      are Copyright 2011 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz from 2014-OCT-21 to 2015-JUN-16
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_OVERLAPCACHE
#define INCLUDE_AS_BAT_OVERLAPCACHE

#include "AS_global.H"
#include "ovStore.H"
#include "gkStore.H"
#include "memoryMappedFile.H"

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



class OverlapCacheThreadData {
public:
  OverlapCacheThreadData() {
    _batMax  = 1 * 1024 * 1024;  //  At 8B each, this is 8MB
    _bat     = new BAToverlap [_batMax];
  };

  ~OverlapCacheThreadData() {
    delete [] _bat;
  };

  uint32                  _batMax;   //  For returning overlaps
  BAToverlap             *_bat;      //
};


class OverlapCache {
public:
  OverlapCache(gkStore *gkp,
               ovStore *ovlStoreUniq,
               ovStore *ovlStoreRept,
               const char *prefix,
               double maxErate,
               uint32 minOverlap,
               uint64 maxMemory,
               uint64 genomeSize,
               bool dosave);
  ~OverlapCache();

private:
  uint32       findHighestOverlapCount(void);
  void         allocateLoadingSpace(void);

  uint32       filterOverlaps(uint32 maxOVSerate, uint32 minOverlap, uint32 no);
  uint32       filterDuplicates(uint32 &no);

  void         computeOverlapLimit(void);
  void         loadOverlaps(bool doSave);
  void         symmetrizeOverlaps(void);

public:
  BAToverlap  *getOverlaps(uint32 readIID, uint32 &numOverlaps) {
    numOverlaps = _overlapLen[readIID];
    return(_overlaps[readIID]);
  }

private:
  bool         load(void);
  void         save(void);

private:
  const char             *_prefix;

  uint64                  _memLimit;
  uint64                  _memUsed;

  BAToverlap            **_overlaps;
  uint32                 *_overlapLen;
  uint32                 *_overlapMax;

  uint32                  _maxEvalue;  //  Don't load overlaps with high error
  uint32                  _minOverlap; //  Don't load overlaps that are short

  uint32                  _minPer;     //  Minimum number of overlaps to retain for a single read
  uint32                  _maxPer;     //  Maximum number of overlaps to load for a single read

  bool                    _checkSymmetry;

  uint32                  _ovsMax;     //  For loading overlaps
  ovOverlap              *_ovs;        //
  uint64                 *_ovsSco;     //  For scoring overlaps during the load
  uint64                 *_ovsTmp;     //  For picking out a score threshold

  uint64                  _threadMax;
  OverlapCacheThreadData *_thread;

  uint64                  _genomeSize;

  gkStore                *_gkp;
  ovStore                *_ovlStoreUniq;  //  Pointers to input stores
  ovStore                *_ovlStoreRept;
};



extern OverlapCache     *OC;

#endif  //  INCLUDE_AS_BAT_OVERLAPCACHE

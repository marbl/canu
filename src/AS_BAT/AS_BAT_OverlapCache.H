
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

#ifndef INCLUDE_AS_BAT_OVERLAPCACHE
#define INCLUDE_AS_BAT_OVERLAPCACHE

static const char *rcsid_INCLUDE_AS_BAT_OVERLAPCACHE = "$Id$";

#include "memoryMappedFile.H"

//  AS_BAT_ERRBITS stores a reduced-precision overlap error rate using the same method
//  as used by OverlapStore.  It should be at least 7 bits, and never needs to be more
//  than AS_OVS_ERRBITS.
//
//  AS_OVS_HNGBITS is one more than the read length.  We have no control over it here.

#define AS_BAT_ERRBITS_MIN  7
#define AS_BAT_ERRBITS_MAX  AS_OVS_ERRBITS

#if (AS_OVS_HNGBITS_MIN + AS_OVS_HNGBITS + AS_OVS_HNGBITS + 1 > 64)
#error No space for BAToverlapInt.  Decrease read length AS_READ_MAX_NORMAL_LEN_BITS.
#endif

//  If not enough space for the minimum number of error bits, bump up to a 64-bit word for overlap
//  storage.

#if (AS_BAT_ERRBITS_MIN + AS_OVS_HNGBITS + AS_OVS_HNGBITS + 1 > 32)

#define AS_BAT_ERRBITS    MIN((64 - AS_OVS_HNGBITS - AS_OVS_HNGBITS - 1), AS_BAT_ERRBITS_MAX)
#define AS_BAT_MAX_ERATE  ((1 << AS_BAT_ERRBITS) - 1)

//  For storing overlaps in memory.  12 bytes per overlap.
struct BAToverlapInt {
  uint64      error     :AS_BAT_ERRBITS;  //   7 by default
  int64       a_hang    :AS_OVS_HNGBITS;  //  12 by default
  int64       b_hang    :AS_OVS_HNGBITS;  //  12 by default
  uint64      flipped   :1;
  AS_IID      b_iid;
};

#else

//  Otherwise, there is enough space for the minimum number of bits, so grab as many as we can up to
//  the maximum we need.

#define AS_BAT_ERRBITS    MIN((32 - AS_OVS_HNGBITS - AS_OVS_HNGBITS - 1), AS_BAT_ERRBITS_MAX)
#define AS_BAT_MAX_ERATE  ((1 << AS_BAT_ERRBITS) - 1)

//  For storing overlaps in memory.  8 bytes per overlap.
struct BAToverlapInt {
  uint32      error     :AS_BAT_ERRBITS;  //   7 by default
  int32       a_hang    :AS_OVS_HNGBITS;  //  12 by default
  int32       b_hang    :AS_OVS_HNGBITS;  //  12 by default
  uint32      flipped   :1;
  AS_IID      b_iid;
};

#endif


//  For working with overlaps, 32 bytes per overlap.  This data is copied
//  from the overlap storage (from a BAToverlapInt) with the erate expanded,
//  and a_iid added.
struct BAToverlap {
  int32       a_hang;
  int32       b_hang;

  uint32      flipped;

  uint32      errorRaw;  //  Encoded fraction error, 0   -- AS_BAT_MAX_ERATE
  double      error;     //  Decoded fraction error, 0.0 -- 1.0

  AS_IID      a_iid;
  AS_IID      b_iid;
};


inline
bool
BAToverlap_sortByErate(BAToverlap const &a, BAToverlap const &b) {
  return(a.error > b.error);
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
  OverlapCache(OverlapStore *ovlStoreUniq,
               OverlapStore *ovlStoreRept,
               const char *prefix,
               double maxErate,
               double elimit,
               uint64 maxMemory,
               uint32 maxOverlaps,
               bool onlysave,
               bool dosave);
  ~OverlapCache();

  void         computeOverlapLimit(void);
  void         computeErateMaps(double erate, double elimit);

  uint32       filterOverlaps(uint32 maxOVSerate, uint32 no);

  void         loadOverlaps(double erate, double elimit, const char *prefix, bool onlySave, bool doSave);

  BAToverlap  *getOverlaps(uint32 fragIID, uint32 &numOverlaps);

  double       decodeError(uint32 error) {
    return(_BATerate[error]);
  };

  double       findError(uint32 aIID, uint32 bIID);

private:
  bool         load(const char *prefix, double erate, double elimit, uint64 memlimit, uint32 maxOverlaps);
  void         save(const char *prefix, double erate, double elimit, uint64 memlimit, uint32 maxOverlaps);

private:
  uint64                  _memLimit;
  uint64                  _memUsed;

  uint32                  _storMax;  //  Size of the heap
  uint32                  _storLen;  //  Position we are at in this heap
  BAToverlapInt          *_stor;     //  Pointer to the heap (probably the last thing on _heaps)

  vector<BAToverlapInt*>  _heaps;

  memoryMappedFile       *_cacheMMF;

  BAToverlapInt         **_cachePtr; //  Mapping of frag iid to overlaps stored in the heap
  uint32                 *_cacheLen; //  Number of overlaps per frag iid

  uint32                  _maxPer;   //  Maximum number of overlaps to load for a single fragment

  uint32                  _ovsMax;   //  For loading overlaps
  OVSoverlap             *_ovs;      //  
  uint64                 *_ovsSco;   //

  uint64                  _threadMax;
  OverlapCacheThreadData *_thread;

  uint32                 *_OVSerate; //  Map from OVS erate to BAT erate
  double                 *_BATerate; //  Map from BAT erate to fraction error

  OverlapStore           *_ovlStoreUniq;  //  Pointers to input stores
  OverlapStore           *_ovlStoreRept;
};

#endif  //  INCLUDE_AS_BAT_OVERLAPCACHE

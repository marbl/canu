
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

//  CA8 used to re-encode the error rate into a smaller-precision number.  This was
//  confusing and broken (it tried to use a log-based encoding to give more precision
//  to the smaller values).  CA3g gives up and uses all 12 bits of precision.

//  If not enough space for the minimum number of error bits, bump up to a 64-bit word for overlap
//  storage.

//  For storing overlaps in memory.  12 bytes per overlap.
struct BAToverlapInt {
  uint64      evalue    :AS_MAX_EVALUE_BITS;     //  12 by default (same as AS_MAX_EVALUE_BITS)
  int64       a_hang    :AS_MAX_READLEN_BITS+1;  //  21+1 by default
  int64       b_hang    :AS_MAX_READLEN_BITS+1;  //  21+1 by default
  uint64      flipped   :1;
  uint32      b_iid;
};

#if (AS_MAX_EVALUE_BITS + AS_MAX_READLEN_BITS + 1 + AS_MAX_READLEN_BITS + 1 + 1 > 64)
#error not enough bits to store overlaps.  decrease AS_MAX_EVALUE_BITS or AS_MAX_READLEN_BITS.
#endif


//  For working with overlaps, 32 bytes per overlap.  This data is copied
//  from the overlap storage (from a BAToverlapInt) with the erate expanded,
//  and a_iid added.
struct BAToverlap {
  int32       a_hang;
  int32       b_hang;

  uint32      flipped;

  uint32      evalue;    //  Encoded fraction error
  double      erate;     //  Decoded fraction error

  uint32      a_iid;
  uint32      b_iid;
};


inline
bool
BAToverlap_sortByErate(BAToverlap const &a, BAToverlap const &b) {
  return(a.erate > b.erate);
}


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
  OverlapCache(ovStore *ovlStoreUniq,
               ovStore *ovlStoreRept,
               const char *prefix,
               double maxErate,
               uint32 minOverlap,
               uint64 maxMemory,
               uint32 maxOverlaps,
               bool onlysave,
               bool dosave);
  ~OverlapCache();

  void         computeOverlapLimit(void);

  uint32       filterOverlaps(uint32 maxOVSerate, uint32 minOverlap, uint32 no);

  void         loadOverlaps(double erate, uint32 minOverlap, const char *prefix, bool onlySave, bool doSave);

  BAToverlap  *getOverlaps(uint32 fragIID, double maxErate, uint32 &numOverlaps);

  void         removeWeakOverlaps(uint32 *minEvalue5p,
                                  uint32 *minEvalue3p);

  double       findErate(uint32 aIID, uint32 bIID);

private:
  bool         load(const char *prefix, double erate, uint64 memlimit, uint32 maxOverlaps);
  void         save(const char *prefix, double erate, uint64 memlimit, uint32 maxOverlaps);

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
  ovsOverlap             *_ovs;      //  
  uint64                 *_ovsSco;   //  For scoring overlaps during the load
  uint64                 *_ovsTmp;   //  For picking out a score threshold

  uint64                  _threadMax;
  OverlapCacheThreadData *_thread;

  ovStore                *_ovlStoreUniq;  //  Pointers to input stores
  ovStore                *_ovlStoreRept;
};

#endif  //  INCLUDE_AS_BAT_OVERLAPCACHE


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

static const char *rcsid_INCLUDE_AS_BAT_OVERLAPCACHE = "$Id: AS_BAT_OverlapCache.H,v 1.3 2011-08-03 07:43:34 brianwalenz Exp $";


#if (AS_OVS_HNGBITS + AS_OVS_HNGBITS + 1 > 32)

#define AS_BAT_ERRBITS    (64 - AS_OVS_HNGBITS - AS_OVS_HNGBITS - 1)
#define AS_BAT_MAX_ERATE  ((1 << AS_BAT_ERRBITS) - 1)

//  For storing overlaps in memory.  8 bytes per overlap.
struct BAToverlapInt {
  uint64      error     :AS_BAT_ERRBITS;  //   7 by default
  int64       a_hang    :AS_OVS_HNGBITS;  //  12 by default
  int64       b_hang    :AS_OVS_HNGBITS;  //  12 by default
  uint64      flipped   :1;
  AS_IID      b_iid;
};

#else

#define AS_BAT_ERRBITS    (32 - AS_OVS_HNGBITS - AS_OVS_HNGBITS - 1)
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

  uint32      errorRaw;
  double      error;

  AS_IID      a_iid;
  AS_IID      b_iid;
};


class OverlapCache {
public:
  OverlapCache(OverlapStore *ovlStoreUniq, OverlapStore *ovlStoreRept, double maxErate, double elimit);
  ~OverlapCache();

  void         computeErateMaps(double erate, double elimit);
  void         loadOverlaps(double erate, double elimit);

  BAToverlap  *getOverlaps(uint32 fragIID, uint32 &numOverlaps);

private:
  uint32                  _storMax;  //  Size of the heap
  uint32                  _storLen;  //  Position we are at in this heap
  BAToverlapInt          *_stor;     //  Pointer to the heap (probably the last thin on _heaps)

  vector<BAToverlapInt*>  _heaps;

  BAToverlapInt         **_cachePtr; //  Mapping of frag iid to overlaps stored in the heap
  uint32                 *_cacheLen; //  Number of overlaps per frag iid

  uint32                  _ovsMax;   //  For loading overlaps
  OVSoverlap             *_ovs;      //  

  uint32                  _batMax;   //  For returning overlaps
  BAToverlap             *_bat;      //  

  uint32                 *_OVSerate; //  Map from OVS erate to BAT erate
  double                 *_BATerate; //  Map from BAT erate to fraction error

  OverlapStore           *_ovlStoreUniq;  //  Pointers to input stores
  OverlapStore           *_ovlStoreRept;
};

#endif  //  INCLUDE_AS_BAT_OVERLAPCACHE

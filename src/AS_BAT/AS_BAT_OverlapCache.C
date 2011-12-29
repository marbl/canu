
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

static const char *rcsid = "$Id: AS_BAT_OverlapCache.C,v 1.11 2011-12-29 09:26:03 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_OverlapCache.H"

#include <sys/types.h>
#include <sys/sysctl.h>

OverlapCache::OverlapCache(OverlapStore *ovlStoreUniq,
                           OverlapStore *ovlStoreRept,
                           double erate,
                           double elimit,
                           uint64 memlimit,
                           uint32 maxOverlaps) {

  _memLimit = memlimit;
  _memUsed = 0;

  //#if HAVE_SYSCTL && defined HW_PHYSMEM
#if 1
  uint64  physMemory = 0;

  int     mib[2] = { CTL_HW, HW_PHYSMEM };
  size_t  len    = sizeof(uint64);

  errno = 0;

  if (sysctl(mib, 2, &physMemory, &len, NULL, 0) != 0)
    //  failed to get memory size, so what?
    fprintf(stderr, "sysctl() failed to return CTL_HW, HW_PHYSMEM: %s\n", strerror(errno)), exit(1);

  if (len != sizeof(uint64))
    //  wasn't enough space, so what?
    fprintf(stderr, "sysctl() failed to return CTL_HW, HW_PHYSMEM: %s\n", strerror(errno)), exit(1);

#else
  uint64  physPages  = sysconf(_SC_PHYS_PAGES);
  uint64  pageSize   = sysconf(_SC_PAGESIZE);
  uint64  physMemory = physPages * pageSize;

  fprintf(stderr, "PHYS_PAGES = "F_U64"\n", physPages);
  fprintf(stderr, "PAGE_SIZE  = "F_U64"\n", pageSize);
  fprintf(stderr, "MEMORY     = "F_U64"\n", physMemory);
#endif

  if (_memLimit == UINT64_MAX)
    _memLimit = physMemory;

  //  Decide on the default block size.  We want to use large blocks (to reduce the number of
  //  allocations, and load on the allocator) but not so large that we can't fit nicely.
  //
  //    1gb blocks @ 64 -> 64gb
  //  128mb blocks @ 64 ->  8gb
  //
  //  below 8gb we'll use 128mb blocks
  //  from  8gb to 64gb, we'll use _memLimit/64
  //  from 64gb on, we'll use 1gb block

  if      (_memLimit <= (uint64)8 * 1024 * 1024 * 1024)
    _storMax = 128 * 1024 * 1024 / sizeof(BAToverlapInt);

  else if (_memLimit <= (uint64)64 * 1024 * 1024 * 1024)
    _storMax = _memLimit / 64 / sizeof(BAToverlapInt);

  else
    _storMax = (uint64)1024 * 1024 * 1024 / sizeof(BAToverlapInt);

  _storLen  = 0;
  _stor     = NULL;

  _cachePtr = new BAToverlapInt * [FI->numFragments() + 1];
  _cacheLen = new uint32          [FI->numFragments() + 1];

  _memUsed += (FI->numFragments() + 1) * (sizeof(BAToverlapInt *) + sizeof(uint32));

  memset(_cachePtr, 0, sizeof(BAToverlapInt *) * (FI->numFragments() + 1));
  memset(_cacheLen, 0, sizeof(uint32)          * (FI->numFragments() + 1));

  _maxPer  = maxOverlaps;

  _ovsMax  = 1 * 1024 * 1024;  //  At 16B each, this is 16MB
  _ovs     = new OVSoverlap [_ovsMax];
  _ovsSco  = new uint64     [_ovsMax];

  _batMax  = 1 * 1024 * 1024;  //  At 8B each, this is 8MB
  _bat     = new BAToverlap [_batMax];

  _memUsed += _ovsMax * sizeof(OVSoverlap);
  _memUsed += _ovsMax * sizeof(uint64);
  _memUsed += _batMax * sizeof(BAToverlap);

  _OVSerate     = NULL;
  _BATerate     = NULL;

  _ovlStoreUniq = ovlStoreUniq;
  _ovlStoreRept = ovlStoreRept;

  assert(_ovlStoreUniq != NULL);
  assert(_ovlStoreRept == NULL);

  if (_memUsed > _memLimit)
    fprintf(stderr, "OverlapCache()-- ERROR: not enough memory to load ANY overlaps.\n"), exit(1);

  computeOverlapLimit();
  computeErateMaps(erate, elimit);
  loadOverlaps(erate, elimit);
}


OverlapCache::~OverlapCache() {

  delete [] _BATerate;
  delete [] _OVSerate;

  delete [] _bat;
  delete [] _ovs;

  delete [] _cacheLen;
  delete [] _cachePtr;

  for (uint32 i=0; i<_heaps.size(); i++)
    delete [] _heaps[i];
}




//  Decide on limits per fragment.
//
//  From the memory limit, we can compute the average allowed per fragment.  If this is higher than
//  the expected coverage, we'll not fill memory completely as the fragments in unique sequence will
//  have fewer than this number of overlaps.
//
//  We'd like to iterate this, but the unused space computation assumes all fragments are assigned
//  the same amount of memory.  On the next iteration, this isn't true any more.  The benefit is
//  (hopefully) small, and the algorithm is unknown.
//
//  This isn't perfect.  It estimates based on whatever is in the store, not only those overlaps
//  below the error threshold.  Result is that memory usage is far below what it should be.  Easy to
//  fix if we assume all fragments have the same properties (same library, same length, same error
//  rate) but not so easy in reality.  We need big architecture changes to make it easy (grouping
//  reads by library, collecting statistics from the overlaps, etc).
//
//  It also doesn't distinguish between 5' and 3' overlaps - it is possible for all the long
//  overlaps to be off of one end.
//
void
OverlapCache::computeOverlapLimit(void) {


  if (_maxPer < UINT32_MAX) {
    //  -N supplied on the command line, use that instead.
    fprintf(stderr, "OverlapCache()-- _maxPer     = "F_U32" overlaps/frag (from command line)\n", _maxPer);
    return;
  }

  AS_OVS_resetRangeOverlapStore(_ovlStoreUniq);

  //  AS_OVS_numOverlapsPerFrag returns an array that starts at firstIIDrequested.  This is usually
  //  1, unless the first fragment has no overlaps.  In that case, firstIIDrequested will be the
  //  first fragment with overlaps.  This is a terrible interface.

  uint32  frstFrag = _ovlStoreUniq->firstIIDrequested;
  uint32  lastFrag = _ovlStoreUniq->lastIIDrequested;
  uint32  totlFrag = lastFrag - frstFrag + 1;

  fprintf(stderr, "OverlapCache()-- Loading number of overlaps per fragment ["F_U32" to "F_U32"]\n", frstFrag, lastFrag);
  uint32 *numPer    = AS_OVS_numOverlapsPerFrag(_ovlStoreUniq);
  uint32  numPerMax = 0;

  for (uint32 i=0; i<totlFrag; i++)
    if (numPerMax < numPer[i])
      numPerMax = numPer[i];

  _maxPer = (_memLimit - _memUsed) / (FI->numFragments() * sizeof(BAToverlapInt));

  fprintf(stderr, "OverlapCache()--  Initial guess at _maxPer="F_U32" (max of "F_U32") from (memLimit="F_U64" - memUsed="F_U64") / (numFrags="F_U32" * sizeof(OVL)="F_SIZE_T")\n",
          _maxPer, numPerMax, _memLimit, _memUsed, FI->numFragments(), sizeof(BAToverlapInt));

  if (_maxPer < 10)
    fprintf(stderr, "OverlapCache()-- ERROR: not enough memory to load overlaps (_maxPer="F_U32").\n", _maxPer), exit(1);

  uint64  totalLoad  = 0;  //  Total overlaps we would load at this threshold

  uint32  numBelow   = 0;  //  Number below the threshold
  uint64  numBelowS  = 0;  //  Amount of space wasted beacuse of this
  uint32  numEqual   = 0;
  uint32  numAbove   = 0;  //  Number of fragments above the threshold

  uint32  lastMax    = 0;

  uint32  adjust     = 1;

  while (adjust > 0) {
    totalLoad = 0;
    numBelow  = 0;
    numBelowS = 0;
    numEqual  = 0;
    numAbove  = 0;

    for (uint32 i=0; i<totlFrag; i++) {
      if (numPer[i] < _maxPer) {
        numBelow++;
        numBelowS  += _maxPer - MAX(lastMax, numPer[i]);
        totalLoad  += numPer[i];

      } else if (numPer[i] == _maxPer) {
        numEqual++;
        totalLoad  += _maxPer;

      } else {
        numAbove++;
        totalLoad  += _maxPer;
      }
    }

    fprintf(stderr, "OverlapCache()-- _maxPer="F_U32" (numBelow="F_U32" numEqual="F_U32" numAbove="F_U32" totalLoad="F_U64" -- "F_U64" + "F_U64" = "F_U64" < "F_U64"\n",
            _maxPer, numBelow, numEqual, numAbove,
            totalLoad, _memUsed, totalLoad + _memUsed,
            totalLoad * sizeof(BAToverlapInt), _memLimit);

    adjust = 0;

    if (numAbove == 0) {
      //  All done, nothing to do here.

    } else if (_memUsed + totalLoad * sizeof(BAToverlapInt) < _memLimit) {
      //  This limit worked, let's try moving it a little higher.

      lastMax  = _maxPer;
      adjust   = (numAbove == 0) ? (0) : (numBelowS / numAbove);
      _maxPer += adjust;

      if (_maxPer > numPerMax)
        _maxPer = numPerMax;

    } else {
      //  Whoops!  Too high!  Revert to the last and recompute statistics.

      _maxPer = lastMax;

      totalLoad = 0;
      numBelow  = 0;
      numBelowS = 0;
      numEqual  = 0;
      numAbove  = 0;

      for (uint32 i=0; i<totlFrag; i++) {
        if (numPer[i] < _maxPer) {
          numBelow++;
          numBelowS  += _maxPer - numPer[i];
          totalLoad  += numPer[i];

        } else if (numPer[i] == _maxPer) {
          numEqual++;
          totalLoad  += _maxPer;

        } else {
          numAbove++;
          totalLoad  += _maxPer;
        }
      }

      fprintf(stderr, "OverlapCache()-- _maxPer="F_U32" (final)\n", _maxPer);

    }
  }

  //  Report

  fprintf(stderr, "\n");
  fprintf(stderr, "OverlapCache()-- blockSize   = "F_U32" ("F_SIZE_T"MB)\n", _storMax, (_storMax * sizeof(BAToverlapInt)) >> 20);
  fprintf(stderr, "\n");
  fprintf(stderr, "OverlapCache()-- _maxPer     = "F_U32" overlaps/frag\n", _maxPer);
  fprintf(stderr, "OverlapCache()-- numBelow    = "F_U32"\n", numBelow);
  fprintf(stderr, "OverlapCache()-- numEqual    = "F_U32"\n", numEqual);
  fprintf(stderr, "OverlapCache()-- numAbove    = "F_U32"\n", numAbove);
  fprintf(stderr, "OverlapCache()-- totalLoad   = "F_U64"\n", totalLoad);
  fprintf(stderr, "\n");
  fprintf(stderr, "OverlapCache()-- memoryLimit = "F_U64"MB\n", _memLimit >> 20);
  fprintf(stderr, "OverlapCache()-- totalMemory = "F_U64"MB for misc\n", _memUsed >> 20);
  fprintf(stderr, "OverlapCache()-- totalMemory = "F_U64"MB for overlaps\n", (totalLoad * sizeof(BAToverlapInt)) >> 20);
  fprintf(stderr, "OverlapCache()-- totalMemory = "F_U64"MB total\n", (_memUsed + totalLoad * sizeof(BAToverlapInt)) >> 20);
  fprintf(stderr, "\n");

  safe_free(numPer);
}


//  Compute erate the scaling maps.
//
//  OVSerate converts a 12-bit OVS error rate to a log-scaled 7-bit error rate.  We can set this
//  in any order.
//
//  BATerate converts the log-scaled 7-bit error rate to a floating point 'fraction error'.
//  This needs to be set high-to-low -- otherwise BATerate[0] = 0.0001.  This isn't critical;
//  as far as I remember, this is only used for relative scoring of overlap-based placement.
//
//  THERE ARE HOLES in this encoding.  By about _OVSerate[16] the holes go away -- value 42.
//
//    _OVSerate[0] = 0
//    _OVSerate[1] = 0
//    _OVSerate[2] = 10
//    _OVSerate[3] = 16
//    _OVSerate[4] = 21
//
//  The holes shouldn't hurt (and are valgrind clean) since we use OVSerate to
//  convert to .error, and we create BATerate for all valid OVSerate values.
//
void
OverlapCache::computeErateMaps(double erate, double elimit) {
  _OVSerate = new uint32 [1 << AS_OVS_ERRBITS];
  _BATerate = new double [1 << AS_BAT_ERRBITS];

  _memUsed += (1 << AS_OVS_ERRBITS) * sizeof(uint32);
  _memUsed += (1 << AS_BAT_ERRBITS) * sizeof(double);

  for (uint32 i=0; i<1 << AS_OVS_ERRBITS; i++)
    _OVSerate[i] = (uint32)floor((1 << AS_BAT_ERRBITS) * log(i) / log(1 << AS_OVS_ERRBITS));

  assert(_OVSerate[1 << AS_OVS_ERRBITS - 1] < 1 << AS_BAT_ERRBITS);

  for (uint32 i=1 << AS_OVS_ERRBITS; i--; )
    _BATerate[_OVSerate[i]] = AS_OVS_decodeQuality(i);
}




uint32
OverlapCache::filterOverlaps(uint32 maxOVSerate, uint32 no) {
  uint32 ns = 0;

  //  If there are fewer overlaps than the limit, accept all of them...below the error threshold.

  if (no <= _maxPer) {
    for (uint32 ii=0; ii<no; ii++)
      if ((_ovs[ii].dat.ovl.corr_erate <= maxOVSerate) &&
          (FI->fragmentLength(_ovs[ii].a_iid) != 0) &&
          (FI->fragmentLength(_ovs[ii].b_iid) != 0)) {
        _ovsSco[ii] = 1;
        ns++;
      } else {
        _ovsSco[ii] = 0;
      }

    return(ns);
  }

  memset(_ovsSco, 0, sizeof(uint64) * no);

  //  A simple filter based on quality; keep if good.
#if 0
  for (uint32 ii=0; ii<no; ii++)
    if (_ovs[ii].dat.ovl.corr_erate <= maxOVSerate)
      _ovsSco[ii] = 1;
#endif

  //  A simple filter based on length & quality of overlap.

  uint32  maxLen = 0;
  uint32  minLen = UINT32_MAX;

  uint64  ERR_MASK = ((uint64)1 << AS_OVS_ERRBITS) - 1;

  uint32  SALT_BITS = (64 - AS_READ_MAX_NORMAL_LEN_BITS - AS_OVS_ERRBITS);
  uint64  SALT_MASK = (((uint64)1 << SALT_BITS) - 1);

  //fprintf(stderr, "SALT_BITS %d SALT_MASK 0x%08x\n", SALT_BITS, SALT_MASK);

  for (uint32 ii=0; ii<no; ii++) {
    _ovsSco[ii]   = FI->overlapLength(_ovs[ii].a_iid, _ovs[ii].b_iid, _ovs[ii].dat.ovl.a_hang, _ovs[ii].dat.ovl.b_hang);
    _ovsSco[ii] <<= AS_OVS_ERRBITS;
    _ovsSco[ii]  |= (~_ovs[ii].dat.ovl.corr_erate) & ERR_MASK;
    _ovsSco[ii] <<= SALT_BITS;
    _ovsSco[ii]  |= ii & SALT_MASK;

    if ((_ovs[ii].dat.ovl.corr_erate > maxOVSerate) ||
        (FI->fragmentLength(_ovs[ii].a_iid) == 0) ||
        (FI->fragmentLength(_ovs[ii].b_iid) == 0))
      _ovsSco[ii] = 0;
  }

  //  Sort by longest overlap, then lowest error.
  //  
  sort(_ovsSco, _ovsSco + no);

  uint64  cutoff = _ovsSco[no - _maxPer];

  for (uint32 ii=0; ii<no; ii++) {
    _ovsSco[ii]   = FI->overlapLength(_ovs[ii].a_iid, _ovs[ii].b_iid, _ovs[ii].dat.ovl.a_hang, _ovs[ii].dat.ovl.b_hang);
    _ovsSco[ii] <<= AS_OVS_ERRBITS;
    _ovsSco[ii]  |= (~_ovs[ii].dat.ovl.corr_erate) & ERR_MASK;
    _ovsSco[ii] <<= SALT_BITS;
    _ovsSco[ii]  |= ii & SALT_MASK;

    if ((_ovs[ii].dat.ovl.corr_erate > maxOVSerate) ||
        (FI->fragmentLength(_ovs[ii].a_iid) == 0) ||
        (FI->fragmentLength(_ovs[ii].b_iid) == 0))
      _ovsSco[ii] = 0;

    if (_ovsSco[ii] < cutoff)
      _ovsSco[ii] = 0;
  }

  //  Count how many overlaps we saved.
  for (uint32 ii=0; ii<no; ii++)
    if (_ovsSco[ii] > 0)
      ns++;

  if (ns > _maxPer)
    fprintf(stderr, "WARNING: fragment "F_U32" loaded "F_U32" overlas (it has "F_U32" in total); over the limit of "F_U32"\n",
            _ovs[0].a_iid, ns, no, _maxPer);

  return(ns);
}




void
OverlapCache::loadOverlaps(double erate, double elimit) {
  uint64   numTotal    = 0;
  uint64   numLoaded   = 0;
  uint32   numFrags    = 0;
  uint32   numOvl      = 0;
  uint32   maxOVSerate = AS_OVS_encodeQuality(erate);

  assert(_ovlStoreUniq != NULL);
  assert(_ovlStoreRept == NULL);

  AS_OVS_resetRangeOverlapStore(_ovlStoreUniq);

  fprintf(logFile, "OverlapCache()-- Loading overlap information\n");

  //  Could probably easily extend to multiple stores.  Needs to interleave the two store
  //  loads, can't do one after the other as we require all overlaps for a single fragment
  //  be in contiguous memory.
  
  while (1) {

    //  Ask the store how many overlaps exist for this fragment.
    numOvl = AS_OVS_readOverlapsFromStore(_ovlStoreUniq, NULL, 0, AS_OVS_TYPE_ANY);

    numTotal += numOvl;

    if (numOvl == 0)
      //  No overlaps?  We're at the end of the store.
      break;

    //  Resize temporary storage space to hold all these overlaps.
    while (_ovsMax <= numOvl) {
      _memUsed -= (_ovsMax) * sizeof(OVSoverlap);
      _memUsed -= (_ovsMax) * sizeof(uint64);
      _ovsMax *= 2;
      delete [] _ovs;
      delete [] _ovsSco;
      _ovs    = new OVSoverlap [_ovsMax];
      _ovsSco = new uint64     [_ovsMax];
      _memUsed += (_ovsMax) * sizeof(OVSoverlap);
      _memUsed += (_ovsMax) * sizeof(uint64);
    }

    //  Actually load the overlaps.
    uint32  no = AS_OVS_readOverlapsFromStore(_ovlStoreUniq, _ovs, _ovsMax, AS_OVS_TYPE_ANY);
    uint32  ns = filterOverlaps(maxOVSerate, no);

    //  Resize the permament storage space for overlaps.
    if ((_storLen + ns > _storMax) ||
        (_stor == NULL)) {
      _storLen = 0;
      _stor    = new BAToverlapInt [_storMax];
      _heaps.push_back(_stor);

      _memUsed += _storMax * sizeof(BAToverlapInt);
    }

    //  Save a pointer to the start of the overlaps for this fragment, and the number of overlaps
    //  that exist.
    _cachePtr[_ovs[0].a_iid] = _stor + _storLen;
    _cacheLen[_ovs[0].a_iid] = ns;

    numLoaded += ns;

    uint32 storEnd = _storLen + ns;

    //  Finally, append the overlaps to the storage.
    for (uint32 ii=0; ii<no; ii++) {
      if (_ovsSco[ii] == 0)
        continue;

      _stor[_storLen].error   = _OVSerate[_ovs[ii].dat.ovl.corr_erate];
      _stor[_storLen].a_hang  = _ovs[ii].dat.ovl.a_hang;
      _stor[_storLen].b_hang  = _ovs[ii].dat.ovl.b_hang;
      _stor[_storLen].flipped = _ovs[ii].dat.ovl.flipped;
      _stor[_storLen].b_iid   = _ovs[ii].b_iid;

      _storLen++;
    }

    assert(storEnd == _storLen);

    if ((numFrags++ % 1000000) == 0)
      fprintf(logFile, "OverlapCache()-- Loading overlap information fragments:%d total:%12"F_U64P" loaded:%12"F_U64P"\n", _ovs[0].a_iid, numTotal, numLoaded);
  }

  fprintf(logFile, "OverlapCache()-- Loading overlap information total:%12"F_U64P" loaded:%12"F_U64P"\n", numTotal, numLoaded);

  delete [] _ovs;
  _ovs = NULL;
}




BAToverlap *
OverlapCache::getOverlaps(uint32 fragIID, uint32 &numOverlaps) {

  numOverlaps = _cacheLen[fragIID];

  while (_batMax <= numOverlaps) {
    _batMax *= 2;
    delete [] _bat;
    _bat = new BAToverlap [_batMax];
  }

  BAToverlapInt *ptr = _cachePtr[fragIID];

  for (uint32 pos=0; pos < numOverlaps; pos++) {
    _bat[pos].a_hang   = ptr[pos].a_hang;
    _bat[pos].b_hang   = ptr[pos].b_hang;

    _bat[pos].flipped  = ptr[pos].flipped;

    _bat[pos].errorRaw = ptr[pos].error;
    _bat[pos].error    = _BATerate[ptr[pos].error];

    _bat[pos].a_iid    = fragIID;
    _bat[pos].b_iid    = ptr[pos].b_iid;
  }

  return(_bat);
}

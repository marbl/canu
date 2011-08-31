
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

static const char *rcsid = "$Id: AS_BAT_OverlapCache.C,v 1.4 2011-08-31 17:39:47 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_OverlapCache.H"


OverlapCache::OverlapCache(OverlapStore *ovlStoreUniq, OverlapStore *ovlStoreRept, double erate, double elimit) {
  _storMax = 128 * 1024 * 1024;  //  At 8B each, this is 1GB.  ON SMALL ASSEMBLIES THIS IS NEVER FULLY USED!
  _storLen = 0;
  _stor    = new BAToverlapInt [_storMax];

  _heaps.push_back(_stor);

  _cachePtr = new BAToverlapInt * [FI->numFragments() + 1];
  _cacheLen = new uint32          [FI->numFragments() + 1];

  memset(_cachePtr, 0, sizeof(BAToverlapInt *) * (FI->numFragments() + 1));
  memset(_cacheLen, 0, sizeof(uint32)          * (FI->numFragments() + 1));

  _ovsMax  = 1 * 1024 * 1024;  //  At 16B each, this is 16MB
  _ovs     = new OVSoverlap [_ovsMax];

  _batMax  = 1 * 1024 * 1024;  //  At 8B each, this is 8MB
  _bat     = new BAToverlap [_batMax];

  _OVSerate     = NULL;
  _BATerate     = NULL;

  _ovlStoreUniq = ovlStoreUniq;
  _ovlStoreRept = ovlStoreRept;

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

  for (uint32 i=0; i<1 << AS_OVS_ERRBITS; i++)
    _OVSerate[i] = (uint32)floor((1 << AS_BAT_ERRBITS) * log(i) / log(1 << AS_OVS_ERRBITS));

  assert(_OVSerate[1 << AS_OVS_ERRBITS - 1] < 1 << AS_BAT_ERRBITS);

  for (uint32 i=1 << AS_OVS_ERRBITS; i--; )
    _BATerate[_OVSerate[i]] = AS_OVS_decodeQuality(i);
}



void
OverlapCache::loadOverlaps(double erate, double elimit) {
  uint64   numTotal    = 0;
  uint64   numLoaded   = 0;
  uint32   numFrags    = 0;
  uint32   numOvl      = 0;
  uint32   maxOVSErate = AS_OVS_encodeQuality(erate);

  assert(_ovlStoreRept == NULL);
  assert(_ovlStoreUniq != NULL);

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
      _ovsMax *= 2;
      delete [] _ovs;
      _ovs = new OVSoverlap [_ovsMax];
    }

    //  Actually load the overlaps.
    uint32  no = AS_OVS_readOverlapsFromStore(_ovlStoreUniq, _ovs, _ovsMax, AS_OVS_TYPE_ANY);
    uint32  ns = 0;

    //  Count how many overlaps we would keep.
    for (uint32 ii=0; ii<no; ii++)
      if ((_ovs[ii].dat.ovl.corr_erate <= maxOVSErate) &&
          (FI->fragmentLength(_ovs[ii].a_iid) > 0) &&
          (FI->fragmentLength(_ovs[ii].b_iid) > 0))
        ns++;

    //  Resize the permament storage space for overlaps.
    if (_storLen + ns > _storMax) {
      _storLen = 0;
      _stor    = new BAToverlapInt [_storMax];
      _heaps.push_back(_stor);
    }

    //  Save a pointer to the start of the overlaps for this fragment, and the number of overlaps
    //  that exist.
    _cachePtr[_ovs[0].a_iid] = _stor + _storLen;
    _cacheLen[_ovs[0].a_iid] = ns;

    numLoaded += ns;

    //  Finally, append the overlaps to the storage.
    for (uint32 ii=0; ii<no; ii++) {
      if ((_ovs[ii].dat.ovl.corr_erate > maxOVSErate) ||
          (FI->fragmentLength(_ovs[ii].a_iid) == 0) ||
          (FI->fragmentLength(_ovs[ii].b_iid) == 0))
        continue;

      _stor[_storLen].error   = _OVSerate[_ovs[ii].dat.ovl.corr_erate];
      _stor[_storLen].a_hang  = _ovs[ii].dat.ovl.a_hang;
      _stor[_storLen].b_hang  = _ovs[ii].dat.ovl.b_hang;
      _stor[_storLen].flipped = _ovs[ii].dat.ovl.flipped;
      _stor[_storLen].b_iid   = _ovs[ii].b_iid;

      _storLen++;
    }

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

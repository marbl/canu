
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
 *    src/AS_BAT/AS_BAT_OverlapCache.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2011-FEB-15 to 2013-OCT-14
 *      are Copyright 2011-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2012-JAN-11
 *      are Copyright 2012 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz from 2014-AUG-06 to 2015-JUN-25
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-APR-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"  //  sizeof(BestEdgeOverlap)
#include "AS_BAT_Unitig.H"            //  sizeof(ufNode)
#include "AS_BAT_Logging.H"

#include "memoryMappedFile.H"

#include <sys/types.h>

uint64  ovlCacheMagic = 0x65686361436c766fLLU;  //0102030405060708LLU;


#undef TEST_LINEAR_SEARCH


#define  ERR_MASK   (((uint64)1 << AS_MAX_EVALUE_BITS) - 1)

#define  SALT_BITS  (64 - AS_MAX_READLEN_BITS - AS_MAX_EVALUE_BITS)
#define  SALT_MASK  (((uint64)1 << SALT_BITS) - 1)


OverlapCache::OverlapCache(gkStore *gkp,
                           ovStore *ovlStoreUniq,
                           ovStore *ovlStoreRept,
                           const char *prefix,
                           double maxErate,
                           uint32 minOverlap,
                           uint64 memlimit,
                           uint64 genomeSize,
                           bool doSave) {

  _prefix = prefix;

  writeStatus("\n");

  if (memlimit == UINT64_MAX) {
    _memLimit = getPhysicalMemorySize();
    writeStatus("OverlapCache()-- limited to " F_U64 "MB memory (total physical memory).\n", _memLimit >> 20);
  }

  else if (memlimit > 0) {
    _memLimit = memlimit;
    writeStatus("OverlapCache()-- limited to " F_U64 "MB memory (user supplied).\n", _memLimit >> 20);
  }

  else {
    _memLimit = UINT64_MAX;
    writeStatus("OverlapCache()-- using unlimited memory (-M 0).\n");
  }

  //  Need to initialize thread data before we can account for their size.
  _threadMax = omp_get_max_threads();
  _thread    = new OverlapCacheThreadData [_threadMax];

  //  And this too.
  _ovsMax  = 1 * 1024 * 1024;  //  At 16B each, this is 16MB

  //  Account for memory used by read data, best overlaps, and tigs.
  //  The chunk graph is temporary, and should be less than the size of the tigs.

  uint64 memFI = RI->memoryUsage();
  uint64 memBE = RI->numReads() * sizeof(BestEdgeOverlap);
  uint64 memUL = RI->numReads() * sizeof(ufNode);             //  For read positions in tigs
  uint64 memUT = RI->numReads() * sizeof(uint32) / 16;        //  For tigs (assumes 32 read / unitig)
  uint64 memID = RI->numReads() * sizeof(uint32) * 2;         //  For maps of read id to unitig id
  uint64 memEP = RI->numReads() * Unitig::epValueSize() * 2;  //  For error profile

  uint64 memC1 = (RI->numReads() + 1) * (sizeof(BAToverlap *) + sizeof(uint32));
  uint64 memC2 = _ovsMax * (sizeof(ovOverlap) + sizeof(uint64) + sizeof(uint64));
  uint64 memC3 = _threadMax * _thread[0]._batMax * sizeof(BAToverlap);
  uint64 memC4 = (RI->numReads() + 1) * sizeof(uint32);

  uint64 memOS = (_memLimit == getPhysicalMemorySize()) ? (0.1 * getPhysicalMemorySize()) : 0.0;

  uint64 memTT = memFI + memBE + memUL + memUT + memID + memC1 + memC2 + memC3 + memC4 + memOS;

  writeStatus("OverlapCache()-- %7" F_U64P "MB for read data.\n",                      memFI >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for best edges.\n",                     memBE >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for unitig layouts.\n",                 memUL >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for tigs.\n",                           memUT >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for id maps.\n",                        memID >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for error profiles.\n",                 memEP >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for overlap cache pointers.\n",         memC1 >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for overlap cache initial bucket.\n",   memC2 >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for overlap cache thread data.\n",      memC3 >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for number of overlaps per read.\n",    memC4 >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for other processes.\n",                memOS >> 20);
  writeStatus("OverlapCache()-- ---------\n");
  writeStatus("OverlapCache()-- %7" F_U64P "MB for data structures (sum of above).\n", memTT >> 20);

  if (_memLimit <= memTT) {
    int64 defecit = (int64)memTT - (int64)_memLimit;

    writeStatus("OverlapCache()-- %7" F_S64P "MB available for overlaps.\n", defecit);
    writeStatus("OverlapCache()-- Out of memory before loading overlaps; increase -M.\n");
    exit(1);
  }

  _memLimit -= memTT;
  _memUsed   = 0;

  writeStatus("OverlapCache()-- %7" F_U64P "MB available for overlaps.\n",             _memLimit >> 20);
  writeStatus("\n");

  _overlaps   = new BAToverlap * [RI->numReads() + 1];
  _overlapLen = new uint32       [RI->numReads() + 1];
  _overlapMax = new uint32       [RI->numReads() + 1];

  memset(_overlaps,   0, sizeof(BAToverlap *) * (RI->numReads() + 1));
  memset(_overlapLen, 0, sizeof(uint32)       * (RI->numReads() + 1));
  memset(_overlapMax, 0, sizeof(uint32)       * (RI->numReads() + 1));

  _maxEvalue     = AS_OVS_encodeEvalue(maxErate);
  _minOverlap    = minOverlap;

  _minPer  = 0;
  _maxPer  = 0;

  _checkSymmetry = false;

  _ovs     = ovOverlap::allocateOverlaps(NULL, _ovsMax);  //  So can't call bgn or end.
  _ovsSco  = new uint64     [_ovsMax];
  _ovsTmp  = new uint64     [_ovsMax];

  _genomeSize    = genomeSize;

  _gkp          = gkp;
  _ovlStoreUniq = ovlStoreUniq;
  _ovlStoreRept = ovlStoreRept;

  assert(_ovlStoreUniq != NULL);
  assert(_ovlStoreRept == NULL);

  if (_memUsed > _memLimit)
    writeStatus("OverlapCache()-- ERROR: not enough memory to load ANY overlaps.\n"), exit(1);

  computeOverlapLimit();
  loadOverlaps(doSave);
  symmetrizeOverlaps();

  delete [] _ovs;       _ovs    = NULL;
  delete [] _ovsSco;    _ovsSco = NULL;
  delete [] _ovsTmp;    _ovsTmp = NULL;
}


OverlapCache::~OverlapCache() {

  for (uint32 rr=0; rr<RI->numReads(); rr++)
    delete [] _overlaps[rr];

  delete [] _overlaps;
  delete [] _overlapLen;
  delete [] _overlapMax;

  delete [] _ovs;

  delete [] _thread;
}



//  Decide on limits per read.
//
//  From the memory limit, we can compute the average allowed per read.  If this is higher than
//  the expected coverage, we'll not fill memory completely as the reads in unique sequence will
//  have fewer than this number of overlaps.
//
//  We'd like to iterate this, but the unused space computation assumes all reads are assigned
//  the same amount of memory.  On the next iteration, this isn't true any more.  The benefit is
//  (hopefully) small, and the algorithm is unknown.
//
//  This isn't perfect.  It estimates based on whatever is in the store, not only those overlaps
//  below the error threshold.  Result is that memory usage is far below what it should be.  Easy to
//  fix if we assume all reads have the same properties (same library, same length, same error
//  rate) but not so easy in reality.  We need big architecture changes to make it easy (grouping
//  reads by library, collecting statistics from the overlaps, etc).
//
//  It also doesn't distinguish between 5' and 3' overlaps - it is possible for all the long
//  overlaps to be off of one end.
//

void
OverlapCache::computeOverlapLimit(void) {

  _ovlStoreUniq->resetRange();

  //  AS_OVS_numOverlapsPerFrag returns an array that starts at firstIIDrequested.  This is usually
  //  1, unless the first read has no overlaps.  In that case, firstIIDrequested will be the
  //  first read with overlaps.  This is a terrible interface.

  writeStatus("OverlapCache()-- Loading number of overlaps per read.\n");

  uint32  frstRead  = 0;
  uint32  lastRead  = 0;
  uint32 *numPer    = _ovlStoreUniq->numOverlapsPerFrag(frstRead, lastRead);
  uint32  totlRead  = lastRead - frstRead + 1;
  uint32  numPerMax = findHighestOverlapCount();

  uint64  memAvail  = (_memLimit - _memUsed);

  //  Set the minimum number of overlaps per read to 2-3x coverage.

  _minPer = 2 * 3 * RI->numBases() / _genomeSize;

  writeStatus("OverlapCache()--  Retain at least " F_U32 " overlaps/read, based on %.2fx coverage.\n",
              _minPer, (double)RI->numBases() / _genomeSize);

  //  Set the maximum number of overlaps per read to a guess of what it will take to fill up memory.

  _maxPer = memAvail / (RI->numReads() * sizeof(BAToverlap));

  writeStatus("OverlapCache()--  Initial guess at " F_U32 " overlaps/read (maximum " F_U32 " overlaps/read).\n",
              _maxPer, numPerMax);

  if (_maxPer < 10)
    writeStatus("OverlapCache()-- ERROR: not enough memory to load overlaps!.\n"), exit(1);

  uint64  totalLoad  = 0;  //  Total overlaps we would load at this threshold
  uint64  totalOlaps = _ovlStoreUniq->numOverlapsInRange();

  uint32  numBelow   = 0;  //  Number of reads below the threshold
  uint32  numEqual   = 0;
  uint32  numAbove   = 0;  //  Number of reads above the threshold

  uint32  lastMax    = 0;

  uint32  adjust     = 1;

  while (adjust > 0) {
    totalLoad = 0;
    numBelow  = 0;
    numEqual  = 0;
    numAbove  = 0;

    for (uint32 i=0; i<totlRead; i++) {
      if (numPer[i] < _maxPer) {
        numBelow++;
        totalLoad  += numPer[i];

      } else if (numPer[i] == _maxPer) {
        numEqual++;
        totalLoad  += _maxPer;

      } else {
        numAbove++;
        totalLoad  += _maxPer;
      }
    }

    writeStatus("OverlapCache()-- %7" F_U32P " overlaps/read - load all for %7" F_U32P " reads, some for %7" F_U32P " reads - %12" F_U64P " overlaps to load - %4" F_U64P "MB\n",
                _maxPer,
                numBelow + numEqual,
                numAbove,
                totalLoad,
                totalLoad * sizeof(BAToverlap) >> 20);


    //  All done, nothing to do here.
    if ((numAbove == 0) && (totalLoad * sizeof(BAToverlap) < memAvail)) {
      adjust = 0;
    }

    //  This limit worked, let's try moving it a little higher.
    else if (totalLoad * sizeof(BAToverlap) < memAvail) {
      lastMax  = _maxPer;

      adjust   = (memAvail - totalLoad * sizeof(BAToverlap)) / numAbove / sizeof(BAToverlap);
      _maxPer += adjust;

      if (_maxPer > numPerMax)
        _maxPer = numPerMax;
    }

    //  Whoops!  Too high!  Revert to the last and recompute statistics.
    else {
      adjust    = 0;
      _maxPer   = lastMax;

      totalLoad = 0;
      numBelow  = 0;
      numEqual  = 0;
      numAbove  = 0;

      for (uint32 i=0; i<totlRead; i++) {
        if (numPer[i] < _maxPer) {
          numBelow++;
          totalLoad  += numPer[i];

        } else if (numPer[i] == _maxPer) {
          numEqual++;
          totalLoad  += _maxPer;

        } else {
          numAbove++;
          totalLoad  += _maxPer;
        }
      }

      writeStatus("OverlapCache()-- _maxPer=%7" F_U32P " (overestimated, revert to last good and stop)\n", _maxPer);
    }
  }

  //  Report

  writeStatus("\n");
  writeStatus("OverlapCache()-- minPer           = " F_U32 " overlaps/reads\n", _minPer);
  writeStatus("OverlapCache()-- maxPer           = " F_U32 " overlaps/reads\n", _maxPer);
  writeStatus("OverlapCache()-- numBelow         = " F_U32 " reads (all overlaps loaded)\n", numBelow);
  writeStatus("OverlapCache()-- numEqual         = " F_U32 " reads (all overlaps loaded)\n", numEqual);
  writeStatus("OverlapCache()-- numAbove         = " F_U32 " reads (some overlaps loaded)\n", numAbove);
  writeStatus("OverlapCache()-- totalLoad        = " F_U64 " overlaps (%6.2f%%)\n", totalLoad, (totalOlaps > 0) ? (100.0 * totalLoad / totalOlaps) : 0.0);
  writeStatus("\n");
  writeStatus("OverlapCache()-- availForOverlaps = " F_U64 "MB\n", memAvail >> 20);
  writeStatus("OverlapCache()-- totalMemory      = " F_U64 "MB for organization\n", _memUsed >> 20);
  writeStatus("OverlapCache()-- totalMemory      = " F_U64 "MB for overlaps\n", (totalLoad * sizeof(BAToverlap)) >> 20);
  writeStatus("OverlapCache()-- totalMemory      = " F_U64 "MB used\n", (_memUsed + totalLoad * sizeof(BAToverlap)) >> 20);
  writeStatus("\n");

  _checkSymmetry = (numAbove > 0) ? true : false;

  delete [] numPer;
}



uint32
OverlapCache::findHighestOverlapCount(void) {
  uint32  fRead    = 0;
  uint32  lRead    = 0;
  uint32 *numPer   = _ovlStoreUniq->numOverlapsPerFrag(fRead, lRead);
  uint32  totlRead = lRead - fRead + 1;

  uint32  numPerMax = 0;

  for (uint32 i=0; i<totlRead; i++)
    if (numPerMax < numPer[i])
      numPerMax = numPer[i];

  delete [] numPer;

  return(numPerMax);
}



void
OverlapCache::allocateLoadingSpace(void) {

  _ovsMax = findHighestOverlapCount();

  _ovs    = ovOverlap::allocateOverlaps(NULL, _ovsMax);  //  So can't call bgn or end.
  _ovsSco = new uint64     [_ovsMax];
  _ovsTmp = new uint64     [_ovsMax];

  _memUsed += (_ovsMax) * sizeof(ovOverlap);
  _memUsed += (_ovsMax) * sizeof(uint64);
  _memUsed += (_ovsMax) * sizeof(uint64);
}



uint32
OverlapCache::filterDuplicates(uint32 &no) {
  uint32   nFiltered = 0;

  for (uint32 ii=0, jj=1; jj<no; ii++, jj++) {
    if (_ovs[ii].b_iid != _ovs[jj].b_iid)
      continue;

    //  Found duplicate B IDs.  Drop one of them.

    nFiltered++;

    //  Drop the shorter overlap, or the one with the higher erate.

    uint32  iilen = RI->overlapLength(_ovs[ii].a_iid, _ovs[ii].b_iid, _ovs[ii].a_hang(), _ovs[ii].b_hang());
    uint32  jjlen = RI->overlapLength(_ovs[jj].a_iid, _ovs[jj].b_iid, _ovs[jj].a_hang(), _ovs[jj].b_hang());

    if (iilen == jjlen) {
      if (_ovs[ii].evalue() < _ovs[jj].evalue())
        jjlen = 0;
      else
        iilen = 0;
    }

    if (iilen < jjlen)
      _ovs[ii].a_iid = _ovs[ii].b_iid = 0;
    else
      _ovs[jj].a_iid = _ovs[jj].b_iid = 0;
  }

  //  If nothing was filtered, return.

  if (nFiltered == 0)
    return(0);

  //  Squeeze out the filtered overlaps.  We used to just copy the last element over any deleted
  //  ones, leaving the list unsorted, but we're now (Nov 2016) binary searching on it, so can't do
  //  that.

  //  Needs to have it's own log.  Lots of stuff here.
  //writeLog("OverlapCache()-- read %u filtered %u overlaps to the same read pair\n", _ovs[0].a_iid, nFiltered);

  for (uint32 ii=0, jj=0; jj<no; ) {
    if (_ovs[jj].a_iid == 0) {
      jj++;
      continue;
    }

    if (ii != jj)
      _ovs[ii] = _ovs[jj];

    ii++;
    jj++;
  }

  no -= nFiltered;

  //  Check for errors.

  bool  errors = false;

  for (uint32 jj=0; jj<no; jj++)
    if ((_ovs[jj].a_iid == 0) || (_ovs[jj].b_iid == 0))
      errors = true;

  if (errors == false)
    return(nFiltered);

  writeLog("ERROR: filtered overlap found in saved list for read %u.  Filtered %u overlaps.\n", _ovs[0].a_iid, nFiltered);

  for (uint32 jj=0; jj<no + nFiltered; jj++)
    writeLog("OVERLAP  %8d %8d  hangs %5d %5d  erate %.4f\n",
             _ovs[jj].a_iid, _ovs[jj].b_iid, _ovs[jj].a_hang(), _ovs[jj].b_hang(), _ovs[jj].erate());

  flushLog();

  assert(errors == false);
  return(0);
}



uint32
OverlapCache::filterOverlaps(uint32 maxEvalue, uint32 minOverlap, uint32 no) {
  uint32 ns = 0;

  for (uint32 ii=0; ii<no; ii++) {
    _ovsSco[ii] = 0;                                //  Overlaps 'continue'd below will be filtered, even if 'no filtering' is needed.

    if ((RI->readLength(_ovs[ii].a_iid) == 0) ||    //  At least one read in the overlap is deleted
        (RI->readLength(_ovs[ii].b_iid) == 0))
      continue;

    if (_ovs[ii].evalue() > maxEvalue)              //  Too noisy to care
      continue;

    uint32  olen = RI->overlapLength(_ovs[ii].a_iid, _ovs[ii].b_iid, _ovs[ii].a_hang(), _ovs[ii].b_hang());

    if (olen < minOverlap)                          //  Too short to care
      continue;

    //  Just right!

    _ovsSco[ii]   = olen;
    _ovsSco[ii] <<= AS_MAX_EVALUE_BITS;
    _ovsSco[ii]  |= (~_ovs[ii].evalue()) & ERR_MASK;
    _ovsSco[ii] <<= SALT_BITS;
    _ovsSco[ii]  |= ii & SALT_MASK;

    ns++;
  }

  if (ns <= _maxPer)                                //  Fewer overlaps than the limit, no filtering needed.
    return(ns);

  //  Otherwise, filter out the short and low quality overlaps and count how many we saved.

  memcpy(_ovsTmp, _ovsSco, sizeof(uint64) * no);

  sort(_ovsTmp, _ovsTmp + no);

  uint64  minScore = _ovsTmp[no - _maxPer];

  ns = 0;

  for (uint32 ii=0; ii<no; ii++)
    if (_ovsSco[ii] < minScore)
      _ovsSco[ii] = 0;
    else
      ns++;

  assert(ns <= _maxPer);

  return(ns);
}



void
OverlapCache::loadOverlaps(bool doSave) {

  if (load() == true)
    return;

  assert(_ovlStoreUniq != NULL);
  assert(_ovlStoreRept == NULL);

  _ovlStoreUniq->resetRange();

  uint64   numTotal     = 0;
  uint64   numLoaded    = 0;
  uint64   numDups      = 0;
  uint32   numReads     = 0;
  uint64   numStore     = _ovlStoreUniq->numOverlapsInRange();

  if (numStore == 0)
    writeStatus("ERROR: No overlaps in overlap store?\n"), exit(1);

  //  Could probably easily extend to multiple stores.  Needs to interleave the two store
  //  loads, can't do one after the other as we require all overlaps for a single read
  //  be in contiguous memory.

  while (1) {
    uint32  numOvl = _ovlStoreUniq->numberOfOverlaps();   //  Query how many overlaps for the next read.

    if (numOvl == 0)    //  If no overlaps, we're at the end of the store.
      break;

    assert(numOvl <= _ovsMax);

    //  Actually load the overlaps, then detect and remove overlaps between the same pair, then
    //  filter short and low quality overlaps.

    uint32  no = _ovlStoreUniq->readOverlaps(_ovs, _ovsMax);     //  no == total overlaps == numOvl
    uint32  nd = filterDuplicates(no);                           //  nd == duplicated overlaps (no is decreased by this amount)
    uint32  ns = filterOverlaps(_maxEvalue, _minOverlap, no);    //  ns == acceptable overlaps

    //  Allocate space for the overlaps.  Allocate a multiple of 8k, assumed to be the page size.
    //
    //  If we're loading all overlaps (ns == no) we don't need to overallocate.  Otherwise, we're
    //  loading only some of them and might have to make a twin later.
    //
    //  Once allocated copy the good overlaps.

    if (ns > 0) {
      uint32  id = _ovs[0].a_iid;

      _overlapMax[id] = (ns == no) ? (ns) : ((((sizeof(BAToverlap) * ns / 8192) + 1) * 8192) / sizeof(BAToverlap));
      _overlapLen[id] = ns;
      _overlaps[id]   = new BAToverlap [ _overlapMax[id] ];

      _memUsed += _overlapMax[id] * sizeof(BAToverlap);

      uint32  oo=0;

      for (uint32 ii=0; ii<no; ii++) {
        if (_ovsSco[ii] == 0)
          continue;

        _overlaps[id][oo].evalue    = _ovs[ii].evalue();
        _overlaps[id][oo].a_hang    = _ovs[ii].a_hang();
        _overlaps[id][oo].b_hang    = _ovs[ii].b_hang();
        _overlaps[id][oo].flipped   = _ovs[ii].flipped();
        _overlaps[id][oo].filtered  = false;
        _overlaps[id][oo].symmetric = false;
        _overlaps[id][oo].a_iid     = _ovs[ii].a_iid;
        _overlaps[id][oo].b_iid     = _ovs[ii].b_iid;

        assert(_overlaps[id][oo].a_iid != 0);
        assert(_overlaps[id][oo].b_iid != 0);

        oo++;
      }

      assert(oo == _overlapLen[id]);
    }

    //  Keep track of what we loaded and didn't.

    numTotal  += no + nd;   //  Because no was decremented by nd in filterDuplicates()
    numLoaded += ns;
    numDups   += nd;

    if ((numReads++ % 100000) == 0)
      writeStatus("OverlapCache()-- Loading: overlaps processed %12" F_U64P " (%06.2f%%) loaded %12" F_U64P " (%06.2f%%) droppeddupe %12" F_U64P " (%06.2f%%)\n",
                  numTotal,  100.0 * numTotal  / numStore,
                  numLoaded, 100.0 * numLoaded / numStore,
                  numDups,   100.0 * numDups   / numStore);
  }

  writeStatus("OverlapCache()-- Loading: overlaps processed %12" F_U64P " (%06.2f%%) loaded %12" F_U64P " (%06.2f%%) droppeddupe %12" F_U64P " (%06.2f%%)\n",
              numTotal,  100.0 * numTotal  / numStore,
              numLoaded, 100.0 * numLoaded / numStore,
              numDups,   100.0 * numDups   / numStore);

  if (doSave == true)
    save();
}



bool
searchForOverlap(BAToverlap *ovl, uint32 ovlLen, uint32 bID) {

#ifdef TEST_LINEAR_SEARCH
  bool linearSearchFound = false;

  for (uint32 ss=0; ss<ovlLen; ss++)
    if (ovl[ss].b_iid == bID) {
      linearSearchFound = true;
      break;
    }
#endif

  //  If not, these are repeats and we should binary search everything.
  //  There will be no short lists where we could exhaustively search.

  int32  F = 0;
  int32  L = ovlLen - 1;
  int32  M = 0;

  while (F <= L) {
    M = (F + L) / 2;

    if (ovl[M].b_iid == bID) {
      ovl[M].symmetric = true;
#ifdef TEST_LINEAR_SEARCH
      assert(linearSearchFound == true);
#endif
      return(true);
    }

    if (ovl[M].b_iid < bID)
      F = M+1;
    else
      L = M-1;
  }

#ifdef TEST_LINEAR_SEARCH
  assert(linearSearchFound == false);
#endif

  return(false);
}




void
OverlapCache::symmetrizeOverlaps(void) {

  if (_checkSymmetry == false)
    return;

  uint32   *nonsymPerRead = new uint32 [RI->numReads() + 1];  //  Overlap in this read is missing it's twin

  //  For each overlap, see if the twin overlap exists.  It is tempting to skip searching if the
  //  b-read has loaded all overlaps (the overlap we're searching for must exist) but we can't.
  //  We must still mark the oevrlap as being symmetric.

  writeStatus("OverlapCache()-- Symmetrizing overlaps -- finding missing twins.\n");

#pragma omp parallel for schedule(dynamic, RI->numReads() / 1000)
  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {
    nonsymPerRead[rr] = 0;

    if ((rr % 100) == 0)
      fprintf(stderr, " %6.3f%%\r", 100.0 * rr / RI->numReads());

    for (uint32 oo=0; oo<_overlapLen[rr]; oo++) {
      uint32  rb = _overlaps[rr][oo].b_iid;

      if (_overlaps[rr][oo].symmetric == true)   //  If already marked, we're done.
        continue;

      //  Search for the twin overlap, and if found, we're done.  The twin is marked as symmetric in the function.

      if (searchForOverlap(_overlaps[rb], _overlapLen[rb], rr)) {
        _overlaps[rr][oo].symmetric = true;
        continue;
      }

      //  Didn't find a twin.  Count how many overlaps we need to create duplicates of.

      nonsymPerRead[rr]++;
    }
  }

  uint64  nOverlaps = 0;
  uint64  nOnly     = 0;
  uint64  nCritical   = 0;

  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {
    nOverlaps += _overlapLen[rr];
    nOnly     += nonsymPerRead[rr];

    if (_overlapLen[rr] <= _minPer)
      nCritical += nonsymPerRead[rr];
  }

  writeStatus("OverlapCache()--                       -- found %llu missing twins in %llu overlaps, %llu are strong.\n", nOnly, nOverlaps, nCritical);

  //  Score all the overlaps (again) and drop the lower quality ones.  We need to drop half of the
  //  non-twin overlaps, but also want to retain some minimum number.

  //  But, there are a bunch of overlaps that fall below our score threshold that are symmetric.  We
  //  need to keep these, only because figuring out which ones are 'saved' above will be a total
  //  pain in the ass.

  double  fractionToDrop = 0.6;

  uint64  nDropped = 0;

#warning this should be parallelized
  writeStatus("OverlapCache()-- Symmetrizing overlaps -- dropping weak non-twin overlaps.\n");

  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {
    if (_overlapLen[rr] <= _minPer)
      continue;

    if ((rr % 100) == 0)
      fprintf(stderr, " %6.3f%%\r", 100.0 * rr / RI->numReads());

    for (uint32 oo=0; oo<_overlapLen[rr]; oo++) {
      _ovsSco[oo]   = RI->overlapLength( _overlaps[rr][oo].a_iid, _overlaps[rr][oo].b_iid, _overlaps[rr][oo].a_hang, _overlaps[rr][oo].b_hang);
      _ovsSco[oo] <<= AS_MAX_EVALUE_BITS;
      _ovsSco[oo]  |= (~_ovs[oo].evalue()) & ERR_MASK;
      _ovsSco[oo] <<= SALT_BITS;
      _ovsSco[oo]  |= oo & SALT_MASK;

      _ovsTmp[oo] = _ovsSco[oo];
    }

    sort(_ovsTmp, _ovsTmp + _overlapLen[rr]);

    uint32  minIdx   = (uint32)floor(nonsymPerRead[rr] * fractionToDrop);

    if (minIdx < _minPer)
      minIdx = _minPer;

    uint64  minScore = _ovsTmp[minIdx];

    for (uint32 oo=0; oo<_overlapLen[rr]; oo++) {
      if ((_ovsSco[oo] < minScore) && (_overlaps[rr][oo].symmetric == false)) {
        nDropped++;
        _overlapLen[rr]--;
        _overlaps[rr][oo] = _overlaps[rr][_overlapLen[rr]];
        _ovsSco      [oo] = _ovsSco      [_overlapLen[rr]];
        oo--;
      }
    }

    for (uint32 oo=0; oo<_overlapLen[rr]; oo++)
      if (_overlaps[rr][oo].symmetric == false)
        assert(minScore <= _ovsSco[oo]);
  }

  delete [] nonsymPerRead;
  nonsymPerRead = NULL;

  writeStatus("OverlapCache()--                       -- dropped %llu overlaps.\n", nDropped);

  //  Finally, run through all the saved overlaps and count how many we need to add to each read.

  uint32   *toAddPerRead  = new uint32 [RI->numReads() + 1];  //  Overlap needs to be added to this read

  for (uint32 rr=0; rr<RI->numReads()+1; rr++)
    toAddPerRead[rr] = 0;

  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {
    for (uint32 oo=0; oo<_overlapLen[rr]; oo++)
      if (_overlaps[rr][oo].symmetric == false)
        toAddPerRead[_overlaps[rr][oo].b_iid]++;
  }

  uint64  nToAdd = 0;

  for (uint32 rr=0; rr<RI->numReads()+1; rr++)
    nToAdd += toAddPerRead[rr];

  writeStatus("OverlapCache()-- Symmetrizing overlaps -- adding %llu missing twin overlaps.\n", nToAdd);

  //  Expand or shrink space for the overlaps.

  for (uint32 rr=0; rr<RI->numReads()+1; rr++)
    if (_overlapLen[rr] + toAddPerRead[rr] > _overlapMax[rr])
      resizeArray(_overlaps[rr], _overlapLen[rr], _overlapMax[rr], _overlapLen[rr] + toAddPerRead[rr] + 2048);

  //  Copy non-twin overlaps to their twin.

  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {
    if ((rr % 100) == 0)
      fprintf(stderr, " %6.3f%%\r", 100.0 * rr / RI->numReads());

    for (uint32 oo=0; oo<_overlapLen[rr]; oo++) {
      if (_overlaps[rr][oo].symmetric == true)
        continue;

      uint32  rb = _overlaps[rr][oo].b_iid;
      uint32  nn = _overlapLen[rb]++;

      _overlaps[rb][nn].evalue    =  _overlaps[rr][oo].evalue;
      _overlaps[rb][nn].a_hang    = (_overlaps[rr][oo].flipped) ? (_overlaps[rr][oo].b_hang) : (-_overlaps[rr][oo].a_hang);
      _overlaps[rb][nn].b_hang    = (_overlaps[rr][oo].flipped) ? (_overlaps[rr][oo].a_hang) : (-_overlaps[rr][oo].b_hang);
      _overlaps[rb][nn].flipped   =  _overlaps[rr][oo].flipped;

      _overlaps[rb][nn].filtered  =  _overlaps[rr][oo].filtered;
      _overlaps[rb][nn].symmetric =  _overlaps[rr][oo].symmetric = true;

      _overlaps[rb][nn].a_iid     =  _overlaps[rr][oo].b_iid;
      _overlaps[rb][nn].b_iid     =  _overlaps[rr][oo].a_iid;

      assert(_overlapLen[rb] <= _overlapMax[rb]);

      assert(toAddPerRead[rb] > 0);
      toAddPerRead[rb]--;
    }
  }

  for (uint32 rr=0; rr<RI->numReads()+1; rr++)
    assert(toAddPerRead[rr] == 0);

  delete [] toAddPerRead;
  toAddPerRead = NULL;

  for (uint32 rr=0; rr<RI->numReads()+1; rr++)
    if (_overlaps[rr] != NULL) {
      assert(_overlaps[rr][0                ].a_iid == rr);
      assert(_overlaps[rr][_overlapLen[rr]-1].a_iid == rr);
    }

  //  Probably should sort again.  Not sure if anything depends on this.

  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {
  }

  writeStatus("OverlapCache()-- Symmetrizing overlaps -- finished.\n");
}



bool
OverlapCache::load(void) {
  char     name[FILENAME_MAX];
  FILE    *file;
  size_t   numRead;

  snprintf(name, FILENAME_MAX, "%s.ovlCache", _prefix);
  if (AS_UTL_fileExists(name, FALSE, FALSE) == false)
    return(false);

  writeStatus("OverlapCache()-- Loading graph from '%s'.\n", name);

  errno = 0;

  file = fopen(name, "r");
  if (errno)
    writeStatus("OverlapCache()-- Failed to open '%s' for reading: %s\n", name, strerror(errno)), exit(1);

  uint64   magic      = ovlCacheMagic;
  uint32   ovserrbits = AS_MAX_EVALUE_BITS;
  uint32   ovshngbits = AS_MAX_READLEN_BITS + 1;

  AS_UTL_safeRead(file, &magic,      "overlapCache_magic",      sizeof(uint64), 1);
  AS_UTL_safeRead(file, &ovserrbits, "overlapCache_ovserrbits", sizeof(uint32), 1);
  AS_UTL_safeRead(file, &ovshngbits, "overlapCache_ovshngbits", sizeof(uint32), 1);

  if (magic != ovlCacheMagic)
    writeStatus("OverlapCache()-- ERROR:  File '%s' isn't a bogart ovlCache.\n", name), exit(1);

  AS_UTL_safeRead(file, &_memLimit, "overlapCache_memLimit", sizeof(uint64), 1);
  AS_UTL_safeRead(file, &_memUsed,  "overlapCache_memUsed",  sizeof(uint64), 1);
  AS_UTL_safeRead(file, &_maxPer,   "overlapCache_maxPer",   sizeof(uint32), 1);

  _threadMax = omp_get_max_threads();
  _thread    = new OverlapCacheThreadData [_threadMax];

  _overlaps   = new BAToverlap * [RI->numReads() + 1];
  _overlapLen = new uint32       [RI->numReads() + 1];
  _overlapMax = new uint32       [RI->numReads() + 1];

  AS_UTL_safeRead(file, _overlapLen, "overlapCache_len", sizeof(uint32), RI->numReads() + 1);
  AS_UTL_safeRead(file, _overlapMax, "overlapCache_max", sizeof(uint32), RI->numReads() + 1);

  for (uint32 rr=0; rr<RI->numReads() + 1; rr++) {
    if (_overlapLen[rr] == 0)
      continue;

    _overlaps[rr] = new BAToverlap [ _overlapMax[rr] ];
    memset(_overlaps[rr], 0xff, sizeof(BAToverlap) * _overlapMax[rr]);

    AS_UTL_safeRead(file, _overlaps[rr], "overlapCache_ovl", sizeof(BAToverlap), _overlapLen[rr]);

    assert(_overlaps[rr][0].a_iid == rr);
  }

  fclose(file);

  return(true);
}



void
OverlapCache::save(void) {
  char  name[FILENAME_MAX];
  FILE *file;

  snprintf(name, FILENAME_MAX, "%s.ovlCache", _prefix);

  writeStatus("OverlapCache()-- Saving graph to '%s'.\n", name);

  errno = 0;

  file = fopen(name, "w");
  if (errno)
    writeStatus("OverlapCache()-- Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);

  uint64   magic      = ovlCacheMagic;
  uint32   ovserrbits = AS_MAX_EVALUE_BITS;
  uint32   ovshngbits = AS_MAX_READLEN_BITS + 1;

  AS_UTL_safeWrite(file, &magic,       "overlapCache_magic",      sizeof(uint64), 1);
  AS_UTL_safeWrite(file, &ovserrbits,  "overlapCache_ovserrbits", sizeof(uint32), 1);
  AS_UTL_safeWrite(file, &ovshngbits,  "overlapCache_ovshngbits", sizeof(uint32), 1);

  AS_UTL_safeWrite(file, &_memLimit,   "overlapCache_memLimit",   sizeof(uint64), 1);
  AS_UTL_safeWrite(file, &_memUsed,    "overlapCache_memUsed",    sizeof(uint64), 1);
  AS_UTL_safeWrite(file, &_maxPer,     "overlapCache_maxPer",     sizeof(uint32), 1);

  AS_UTL_safeWrite(file,  _overlapLen, "overlapCache_len",        sizeof(uint32), RI->numReads() + 1);
  AS_UTL_safeWrite(file,  _overlapMax, "overlapCache_max",        sizeof(uint32), RI->numReads() + 1);

  for (uint32 rr=0; rr<RI->numReads() + 1; rr++)
    AS_UTL_safeWrite(file,  _overlaps[rr],   "overlapCache_ovl", sizeof(BAToverlap), _overlapLen[rr]);

  fclose(file);
}


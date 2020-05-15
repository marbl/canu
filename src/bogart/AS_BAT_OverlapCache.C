
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

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"  //  sizeof(BestEdgeOverlap)
#include "AS_BAT_Unitig.H"            //  sizeof(ufNode)
#include "AS_BAT_Logging.H"

#include "system.H"

#include <sys/types.h>

uint64  ovlCacheMagic = 0x65686361436c766fLLU;  //0102030405060708LLU;


#undef TEST_LINEAR_SEARCH


#define  ERR_MASK   (((uint64)1 << AS_MAX_EVALUE_BITS) - 1)

#define  SALT_BITS  (64 - AS_MAX_READLEN_BITS - AS_MAX_EVALUE_BITS)
#define  SALT_MASK  (((uint64)1 << SALT_BITS) - 1)


OverlapCache::OverlapCache(const char *ovlStorePath,
                           const char *prefix,
                           double maxErate,
                           uint32 minOverlap,
                           uint64 memlimit,
                           uint64 genomeSize) {

  _prefix = prefix;

  writeStatus("\n");

  if (memlimit == UINT64_MAX) {
    _memLimit = getPhysicalMemorySize();
    writeStatus("OverlapCache()-- limited to " F_U64 "MB memory (total physical memory).\n", _memLimit >> 20);
    writeStatus("\n");
  }

  else if (memlimit > 0) {
    _memLimit = memlimit;
    writeStatus("OverlapCache()-- limited to " F_U64 "MB memory (user supplied).\n", _memLimit >> 20);
    writeStatus("\n");
  }

  else {
    _memLimit = UINT64_MAX;
    writeStatus("OverlapCache()-- using unlimited memory (-M 0).\n");
    writeStatus("\n");
  }

  //  Account for memory used by read data, best overlaps, and tigs.
  //  The chunk graph is temporary, and should be less than the size of the tigs.
  //  Likewise, the buffers used for loading and scoring overlaps aren't accounted for.
  //
  //  NOTES:
  //
  //  memFI - read length,
  //
  //  memUT - worst case, we have one unitig per read.  also, maps of read-to-unitig and read-to-vector-position.
  //
  //  memEP - each read adds two epValue points, the open and close points, and two uint32 pointers
  //  to the data.
  //
  //  memEO - overlaps for computing error profiles.  this is definitely a hack, but I can't think of
  //  any reasonable estimates.  just reserve 25% of memory, which then dominates our accounting.
  //
  //  memOS - make sure we're this much below using all the memory - allows for other stuff to run,
  //  and a little buffer in case we're too big.

  uint64 memFI = RI->memoryUsage();
  uint64 memBE = RI->numReads() * sizeof(BestEdgeOverlap) * 2;
  uint64 memUT = RI->numReads() * sizeof(Unitig) + RI->numReads() * sizeof(uint32) * 2;
  uint64 memUL = RI->numReads() * sizeof(ufNode);

  uint64 memEP = RI->numReads() * sizeof(uint32) * 2 + RI->numReads() * Unitig::epValueSize() * 2;
  uint64 memEO = (_memLimit == UINT64_MAX) ? (0.0) : (0.25 * _memLimit);

  uint64 memOS = (_memLimit < 0.9 * getPhysicalMemorySize()) ? (0.0) : (0.1 * getPhysicalMemorySize());

  uint64 memST = ((RI->numReads() + 1) * (sizeof(BAToverlap *) + sizeof(uint32)) +   //  Cache pointers
                  (RI->numReads() + 1) * sizeof(uint32) +                            //  Num olaps stored per read
                  (RI->numReads() + 1) * sizeof(uint32));                            //  Num olaps allocated per read


  _memReserved = memFI + memBE + memUL + memUT + memEP + memEO + memST + memOS;
  _memStore    = memST;
  _memAvail    = (_memReserved + _memStore < _memLimit) ? (_memLimit - _memReserved - _memStore) : 0;
  _memOlaps    = 0;

  writeStatus("OverlapCache()-- %7" F_U64P "MB for read data.\n",                      memFI >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for best edges.\n",                     memBE >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for tigs.\n",                           memUT >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for tigs - read layouts.\n",            memUL >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for tigs - error profiles.\n",          memEP >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for tigs - error profile overlaps.\n",  memEO >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for other processes.\n",                memOS >> 20);
  writeStatus("OverlapCache()-- ---------\n");
  writeStatus("OverlapCache()-- %7" F_U64P "MB for data structures (sum of above).\n", _memReserved >> 20);
  writeStatus("OverlapCache()-- ---------\n");
  writeStatus("OverlapCache()-- %7" F_U64P "MB for overlap store structure.\n",        _memStore >> 20);
  writeStatus("OverlapCache()-- %7" F_U64P "MB for overlap data.\n",                   _memAvail >> 20);
  writeStatus("OverlapCache()-- ---------\n");
  writeStatus("OverlapCache()-- %7" F_U64P "MB allowed.\n",                            _memLimit >> 20);
  writeStatus("OverlapCache()--\n");

  if (_memAvail == 0) {
    writeStatus("OverlapCache()-- Out of memory before loading overlaps; increase -M.\n");
    exit(1);
  }

  _maxEvalue     = AS_OVS_encodeEvalue(maxErate);
  _minOverlap    = minOverlap;

  //  Allocate space to load overlaps.  With a NULL seqStore we can't call the bgn or end methods.

  _ovsMax  = 0;
  _ovs     = NULL;
  _ovsSco  = NULL;
  _ovsTmp  = NULL;

  //  Allocate pointers to overlaps.

  _overlapLen = new uint32       [RI->numReads() + 1];
  _overlapMax = new uint32       [RI->numReads() + 1];
  _overlaps   = new BAToverlap * [RI->numReads() + 1];

  memset(_overlapLen, 0, sizeof(uint32)       * (RI->numReads() + 1));
  memset(_overlapMax, 0, sizeof(uint32)       * (RI->numReads() + 1));
  memset(_overlaps,   0, sizeof(BAToverlap *) * (RI->numReads() + 1));

  //  Open the overlap store.

  ovStore *ovlStore = new ovStore(ovlStorePath, NULL);

  //  Load overlaps!

  computeOverlapLimit(ovlStore, genomeSize);
  loadOverlaps(ovlStore);

  delete [] _ovs;       _ovs      = NULL;   //  There is a small cost with these arrays that we'd
  delete [] _ovsSco;    _ovsSco   = NULL;   //  like to not have, and a big cost with ovlStore (in that
  delete [] _ovsTmp;    _ovsTmp   = NULL;   //  it loaded updated erates into memory), so release
  delete     ovlStore;   ovlStore = NULL;   //  these before symmetrizing overlaps.

  symmetrizeOverlaps();

  delete [] _minSco;    _minSco   = NULL;
}


OverlapCache::~OverlapCache() {

  delete [] _overlaps;
  delete [] _overlapLen;
  delete [] _overlapMax;

  delete    _overlapStorage;
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
OverlapCache::computeOverlapLimit(ovStore *ovlStore, uint64 genomeSize) {
  uint32  frstRead  = 0;
  uint32  lastRead  = 0;
  uint32 *numPer    = ovlStore->numOverlapsPerRead();

  //  Set the minimum number of overlaps per read to twice coverage.  Then set the maximum number of
  //  overlaps per read to a guess of what it will take to fill up memory.

  _minPer = 2 * RI->numBases() / genomeSize;
  _maxPer = _memAvail / (RI->numReads() * sizeof(BAToverlap));

  writeStatus("OverlapCache()-- Retain at least " F_U32 " overlaps/read, based on %.2fx coverage.\n", _minPer, (double)RI->numBases() / genomeSize);
  writeStatus("OverlapCache()-- Initial guess at " F_U32 " overlaps/read.\n", _maxPer);
  writeStatus("OverlapCache()--\n");

  uint64  totalOlaps = ovlStore->numOverlapsInRange();

  assert(totalOlaps > 0);

  uint64  olapLoad   = 0;  //  Total overlaps we would load at this threshold
  uint64  olapMem    = 0;

  uint32  numBelow   = 0;  //  Number of reads below the threshold
  uint32  numEqual   = 0;
  uint32  numAbove   = 0;  //  Number of reads above the threshold

  writeStatus("OverlapCache()-- Adjusting for sparse overlaps.\n");
  writeStatus("OverlapCache()--\n");
  writeStatus("OverlapCache()--               reads loading olaps          olaps               memory\n");
  writeStatus("OverlapCache()--   olaps/read       all      some          loaded                 free\n");
  writeStatus("OverlapCache()--   ----------   -------   -------     ----------- -------     --------\n");

  while (true) {
    olapLoad = 0;
    numBelow = 0;
    numEqual = 0;
    numAbove = 0;

    for (uint32 i=1; i<=RI->numReads(); i++) {
      if (numPer[i] < _maxPer) {
        numBelow += 1;
        olapLoad += numPer[i];

      } else if (numPer[i] == _maxPer) {
        numEqual += 1;
        olapLoad += _maxPer;

      } else {
        numAbove += 1;
        olapLoad += _maxPer;
      }
    }

    olapMem = olapLoad * sizeof(BAToverlap);

    //  If we're too high, decrease the threshold and compute again.  We shouldn't ever be too high.

    if (_memAvail < olapMem) {
      _maxPer--;
      continue;
    }

    //  Log what we will be loading.

    writeStatus("OverlapCache()--      %7" F_U32P "   %7" F_U32P "   %7" F_U32P "    %12" F_U32P " %6.2f%%    %7" F_U32P " MB\n",
                _maxPer,
                numBelow + numEqual,
                numAbove,
                olapLoad,
                100.0 * olapLoad / totalOlaps,
                (_memAvail - olapMem) >> 20);

    //  If there are no more overlaps to load, we're done.

    if (numAbove == 0)
      break;

    //  Otherwise, there is still (potentially) space left for more overlaps.  Estimate how much
    //  higher we could push the threshold: compute how many more overlaps we could load before
    //  exceeding the memory limit, then assume we'd load that many overlaps for each of the
    //  numAbove reads.

    int64  olapFree  = (_memAvail - olapMem) / sizeof(BAToverlap);
    int64  increase  = olapFree / numAbove;

    if (increase == 0)
      break;

    _maxPer += increase;
  }

  //  We used to (pre 6 Jul 2017) do the symmetry check only if we didn't load all overlaps.
  //  However, symmetry can also break if we use an error rate cutoff because - for reasons not
  //  explored - the error rate on symmetric overlaps differs.  So, just enable this always.
  //
  //  On a moderate coverage human nanopore assembly, it does:
  //
  //    OverlapCache()-- Symmetrizing overlaps -- finding missing twins.
  //    OverlapCache()--                       -- found 8609 missing twins in 51413413 overlaps, 8002 are strong.
  //    OverlapCache()-- Symmetrizing overlaps -- dropping weak non-twin overlaps.
  //    OverlapCache()--                       -- dropped 454 overlaps.
  //    OverlapCache()-- Symmetrizing overlaps -- adding 8155 missing twin overlaps.

  _checkSymmetry = (numAbove > 0) ? true : false;
  _checkSymmetry = true;

  if (_maxPer < _minPer)
    writeStatus("OverlapCache()-- Not enough memory to load the minimum number of overlaps; increase -M.\n"), exit(1);

  delete [] numPer;
}

//return true if o1 is worse than o2
bool OverlapCache::compareOverlaps(const ovOverlap &o1, const ovOverlap &o2) const {
   assert(o1.a_iid == o2.a_iid);
   assert(o1.b_iid == o2.b_iid);
   auto as_tuple = [&](const ovOverlap &o) {
      return std::make_tuple(1 - o.erate(), RI->overlapLength(o.a_iid, o.b_iid, o.a_hang(), o.b_hang()), !o.flipped());
   };
   return as_tuple(o2) > as_tuple(o1);
}


//return true if o1 is worse than o2
bool OverlapCache::compareOverlaps(const BAToverlap &o1, const BAToverlap &o2) const {
   assert(o1.a_iid == o2.a_iid);
   auto as_tuple = [&](const BAToverlap &o) {
      return std::make_tuple(1 - o.erate(), RI->overlapLength(o.a_iid, o.b_iid, o.a_hang, o.b_hang), !o.flipped);
   };
   return as_tuple(o2) > as_tuple(o1);
}

uint32
OverlapCache::filterDuplicates(uint32 &no) {
  uint32   nFiltered = 0;

  for (uint32 ii=0, jj=1; jj<no; ii++, jj++) {
    if (_ovs[ii].b_iid != _ovs[jj].b_iid)
      continue;

    //  Found duplicate B IDs.  Drop one of them.

    nFiltered++;

    //  Drop the weaker overlap.  If a tie, drop the flipped one.

    uint32 to_drop = compareOverlaps(_ovs[ii], _ovs[jj]) ? ii : jj;
    uint32 to_save = (to_drop == ii ? jj : ii);
#if 0
    writeLog("OverlapCache::filterDuplicates()-- Dropping overlap A: %9" F_U64P " B: %9" F_U64P " - score %8.2f - %6.4f%% - %6" F_S32P " %6" F_S32P " - %s\n",
             _ovs[to_drop].a_iid, _ovs[to_drop].b_iid, to_drops, _ovs[to_drop].erate(), _ovs[to_drop].a_hang(), _ovs[to_drop].b_hang(), _ovs[to_drop].flipped() ? "flipped" : "");
    writeLog("OverlapCache::filterDuplicates()-- Saving   overlap A: %9" F_U64P " B: %9" F_U64P " - score %8.2f - %6.4f%% - %6" F_S32P " %6" F_S32P " - %s\n",
             _ovs[to_save].a_iid, _ovs[to_save].b_iid, to_saves, _ovs[to_save].erate(), _ovs[to_save].a_hang(), _ovs[to_save].b_hang(), _ovs[to_save].flipped() ? "flipped" : "");
#endif

    _ovs[to_drop].a_iid = 0;
    _ovs[to_drop].b_iid = 0;
  }

  //  If nothing was filtered, return.

  if (nFiltered == 0)
    return(0);

  //  Squeeze out the filtered overlaps.  Preserve order so we can binary search later.

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



inline
uint64
ovlSco(uint64 olen, uint64 evalue, uint64 ii) {
  uint64  os;

  os   = olen;
  os <<= AS_MAX_EVALUE_BITS;
  os  |= (~evalue) & ERR_MASK;
  os <<= SALT_BITS;
  os  |= ii & SALT_MASK;

  return(os);
}

inline
uint64
ovlScoToLength(uint64 score) {
  return(score >> (AS_MAX_EVALUE_BITS + SALT_BITS));
}



uint32
OverlapCache::filterOverlaps(uint32 aid, uint32 maxEvalue, uint32 minOverlap, uint32 no) {
  uint32 ns        = 0;
  bool   beVerbose = false;

 //beVerbose = (_ovs[0].a_iid == 3514657);

  for (uint32 ii=0; ii<no; ii++) {
    _ovsSco[ii] = 0;                                //  Overlaps 'continue'd below will be filtered, even if 'no filtering' is needed.
    _ovsTmp[ii] = 0;

    if ((RI->readLength(_ovs[ii].a_iid) == 0) ||    //  At least one read in the overlap is deleted
        (RI->readLength(_ovs[ii].b_iid) == 0)) {
      if (beVerbose)
        writeLog("olap %d involves deleted reads - %u %s - %u %s\n",
                ii,
                _ovs[ii].a_iid, (RI->readLength(_ovs[ii].a_iid) == 0) ? "deleted" : "active",
                _ovs[ii].b_iid, (RI->readLength(_ovs[ii].b_iid) == 0) ? "deleted" : "active");
      continue;
    }

    if (_ovs[ii].evalue() > maxEvalue) {            //  Too noisy to care
      if (beVerbose)
        writeLog("olap %d too noisy evalue %f > maxEvalue %f\n",
                ii, AS_OVS_decodeEvalue(_ovs[ii].evalue()), AS_OVS_decodeEvalue(maxEvalue));
      continue;
    }

    uint32  olen = RI->overlapLength(_ovs[ii].a_iid, _ovs[ii].b_iid, _ovs[ii].a_hang(), _ovs[ii].b_hang());

    //  If too short, drop it.

    if (olen < minOverlap) {
      if (beVerbose)
        writeLog("olap %d too short olen %u minOverlap %u\n",
                ii, olen, minOverlap);
      continue;
    }

    //  Just right!

    _ovsTmp[ii] = _ovsSco[ii] = ovlSco(olen, _ovs[ii].evalue(), ii);

    ns++;
  }

  _minSco[aid] = 0;                          //  Minimum score of overlap we'll keep for this read.

  if (ns <= _maxPer)                         //  Fewer overlaps than the limit, no filtering needed.
    return(ns);

  sort(_ovsTmp, _ovsTmp + no);               //  Sort the scores so we can pick a minScore that
  _minSco[aid] = _ovsTmp[no - _maxPer];      //  results in the correct number of overlaps.

  ns = 0;

  for (uint32 ii=0; ii<no; ii++)
    if (_ovsSco[ii] < _minSco[aid])          //  Score too low, flag it as junk.
      _ovsSco[ii] = 0;                       //  We could also do this when copying overlaps to
    else                                     //  storage, except we need to know how many overlaps
      ns++;                                  //  to copy so we can allocate storage.

  assert(ns <= _maxPer);

  return(ns);
}



void
OverlapCache::loadOverlaps(ovStore *ovlStore) {

  writeStatus("OverlapCache()--\n");
  writeStatus("OverlapCache()-- Loading overlaps.\n");
  writeStatus("OverlapCache()--\n");
  writeStatus("OverlapCache()--          read from store           saved in cache\n");
  writeStatus("OverlapCache()--   ------------ ---------   ------------ ---------\n");

  uint64   numTotal     = 0;
  uint64   numLoaded    = 0;
  uint64   numDups      = 0;
  uint32   numReads     = 0;
  uint64   numStore     = ovlStore->numOverlapsInRange();

  assert(numStore > 0);

  _overlapStorage = new OverlapStorage(ovlStore->numOverlapsInRange());

  //  Scan the overlaps, finding the maximum number of overlaps for a single read.  This lets
  //  us pre-allocate space and simplifies the loading process.

  assert(_ovsMax == 0);
  assert(_ovs    == NULL);

  _ovsMax = 0;

  for (uint32 rr=0; rr<RI->numReads()+1; rr++)
    _ovsMax = max(_ovsMax, ovlStore->numOverlaps(rr));

  _minSco  = new uint64    [RI->numReads()+1];

  _ovs     = new ovOverlap [_ovsMax];
  _ovsSco  = new uint64    [_ovsMax];
  _ovsTmp  = new uint64    [_ovsMax];

  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {

    //  Actually load the overlaps, then detect and remove overlaps between
    //  the same pair, then filter short and low quality overlaps.

    uint32  no = ovlStore->loadOverlapsForRead(rr, _ovs, _ovsMax);   //  no == total overlaps == numOvl
    uint32  nd = filterDuplicates(no);                               //  nd == duplicated overlaps (no is decreased by this amount)
    uint32  ns = filterOverlaps(rr, _maxEvalue, _minOverlap, no);    //  ns == acceptable overlaps

    //  If we still have overlaps, get official space to store them and copy
    //  to storage,

    if (ns > 0) {
      uint32  id = _ovs[0].a_iid;

      _overlapMax[id] = ns;
      _overlapLen[id] = ns;
      _overlaps[id]   = _overlapStorage->get(ns);                //  Get space for overlaps.

      _memOlaps += _overlapMax[id] * sizeof(BAToverlap);

      uint32  oo=0;

      for (uint32 ii=0; ii<no; ii++) {
        if (_ovsSco[ii] == 0)                                    //  Skip if it was filtered.
          continue;

        _overlaps[id][oo].evalue    = _ovs[ii].evalue();         //  Or copy to our storage.
        _overlaps[id][oo].a_hang    = _ovs[ii].a_hang();
        _overlaps[id][oo].b_hang    = _ovs[ii].b_hang();
        _overlaps[id][oo].flipped   = _ovs[ii].flipped();
        _overlaps[id][oo].filtered  = false;
        _overlaps[id][oo].symmetric = false;
        _overlaps[id][oo].a_iid     = _ovs[ii].a_iid;
        _overlaps[id][oo].b_iid     = _ovs[ii].b_iid;

        assert(_overlaps[id][oo].a_iid != 0);   //  Guard against some kind of weird error that
        assert(_overlaps[id][oo].b_iid != 0);   //  I can no longer remember.

        oo++;
      }

      assert(oo == _overlapLen[id]);    //  Ensure we got all the overlaps we were supposed to get.
    }

    //  Keep track of what we loaded and didn't.

    numTotal  += no + nd;   //  Because no was decremented by nd in filterDuplicates()
    numLoaded += ns;
    numDups   += nd;

    if ((numReads++ % 100000) == 99999)
      writeStatus("OverlapCache()--   %12" F_U64P " (%06.2f%%)   %12" F_U64P " (%06.2f%%)\n",
                  numTotal,  100.0 * numTotal  / numStore,
                  numLoaded, 100.0 * numLoaded / numStore);
  }

  writeStatus("OverlapCache()--   ------------ ---------   ------------ ---------\n");
  writeStatus("OverlapCache()--   %12" F_U64P " (%06.2f%%)   %12" F_U64P " (%06.2f%%)\n",
              numTotal,  100.0 * numTotal  / numStore,
              numLoaded, 100.0 * numLoaded / numStore);

  writeStatus("OverlapCache()--\n");
  writeStatus("OverlapCache()-- Ignored %lu duplicate overlaps.\n", numDups);
}



//  Binary search a list of overlaps for one matching bID and flipped.
uint32
searchForOverlap(BAToverlap *ovl, uint32 ovlLen, uint32 bID, bool flipped) {
  int32  F = 0;
  int32  L = ovlLen - 1;
  int32  M = 0;

#ifdef TEST_LINEAR_SEARCH
  bool linearSearchFound = false;

  for (uint32 ss=0; ss<ovlLen; ss++)
    if ((ovl[ss].b_iid   == bID) &&
        (ovl[ss].flipped == flipped)) {
      linearSearchFound = true;
      break;
    }
#endif

  while (F <= L) {
    M = (F + L) / 2;

    if ((ovl[M].b_iid   == bID) &&
        (ovl[M].flipped == flipped)) {
#ifdef TEST_LINEAR_SEARCH
      assert(linearSearchFound == true);
#endif
      return(M);
    }

    if (((ovl[M].b_iid  < bID)) ||
        ((ovl[M].b_iid == bID) && (ovl[M].flipped < flipped)))
      F = M+1;
    else
      L = M-1;
  }

#ifdef TEST_LINEAR_SEARCH
  assert(linearSearchFound == false);
#endif

  return(UINT32_MAX);
}




void
OverlapCache::symmetrizeOverlaps(void) {
  uint32  fiLimit    = RI->numReads() + 1;
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 1000 * numThreads) ? numThreads : fiLimit / 999;

  if (_checkSymmetry == false)
    return;

  uint32  *nNonSymPerRead = new uint32 [fiLimit];
  uint32  *nFiltPerRead   = new uint32 [fiLimit];
  uint32  *nMissPerRead   = new uint32 [fiLimit];

  for (uint32 rr=0; rr < fiLimit; rr++) {
    nNonSymPerRead[rr] = 0;
    nFiltPerRead[rr]   = 0;
    nMissPerRead[rr]   = 0;
  }

  writeStatus("OverlapCache()--\n");
  writeStatus("OverlapCache()-- Symmetrizing overlaps.\n");
  writeStatus("OverlapCache()--   Finding missing twins.\n");

  FILE *NSE = AS_UTL_openOutputFile(_prefix, '.', "non-symmetric-error-rates", false);   //  These turn out to be VERY big
  FILE *NTW = AS_UTL_openOutputFile(_prefix, '.', "non-symmetric-overlaps",    false);

  if (NSE) {
    fprintf(NSE, "     aID      bID  a error b error\n");
    fprintf(NSE, "-------- --------  ------- -------\n");
  }

  //  For each overlap, see if the twin overlap exists.  It is tempting to skip searching if the
  //  b-read has loaded all overlaps (the overlap we're searching for must exist) but we can't.
  //  We must still mark the overlap as being symmetric.

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 ra=0; ra < fiLimit; ra++) {
    for (uint32 oa=0; oa<_overlapLen[ra]; oa++) {
      BAToverlap  *ova = &_overlaps[ra][oa];
      uint32       rb  =  _overlaps[ra][oa].b_iid;

      //  If already marked, we're done.

      if (ova->symmetric == true)
        continue;

      //  Search for the twin overlap.

      uint32       ob  = searchForOverlap(_overlaps[rb], _overlapLen[rb], ra, ova->flipped);

      //  If the twin was found, mark both as symmetric and make sure the
      //  error rate is symmetric too.

      if (ob < UINT32_MAX) {
        BAToverlap  *ovb = &_overlaps[rb][ob];

        ova->symmetric = true;   //  I have a twin!
        ovb->symmetric = true;   //  My twin has a twin, me!

        if (ova->evalue != ovb->evalue) {
          if (NSE)
            fprintf(NSE, "%8u %8u  %7.3f %7.3f\n",
                     ra, rb,
                     ova->erate() * 100.0,
                     ovb->erate() * 100.0);

          uint64 ev = min(ova->evalue, ovb->evalue);

          ova->evalue = ev;
          ovb->evalue = ev;

          nNonSymPerRead[ra]++;
        }

        continue;
      }

      //  Otherwise, the twin was not found.

      //  If this overlap is weak according to the other read (that is, if
      //  the other read dropped it for whatever reason) flag that we don't
      //  want to keep it here either.

      uint32  olen = RI->overlapLength(ova->a_iid, ova->b_iid, ova->a_hang, ova->b_hang);
      uint64  osco = ovlSco(olen, ova->evalue, UINT64_MAX);

      if (osco < _minSco[rb]) {
        if (NTW)
          fprintf(NTW, "NO TWIN for %6u -> %6u - length %lu < min %lu - WEAK\n",
                  ra, ova->b_iid, ovlScoToLength(osco), ovlScoToLength(_minSco[rb]));

        nFiltPerRead[ra]++;
        ova->filtered = true;
        continue;
      }

      //  Not weak.  We want to duplicate this in the other read.

      if (NTW)
        fprintf(NTW, "NO TWIN for %6u -> %6u - length %lu >= min %lu\n",
                ra, ova->b_iid, ovlScoToLength(osco), ovlScoToLength(_minSco[rb]));

#pragma omp critical (nMissPerRead)
      nMissPerRead[rb]++;
    }
  }

  AS_UTL_closeFile(NTW);
  AS_UTL_closeFile(NSE);

  uint64   nOverlaps  = 0;
  uint64   nNonSymErr = 0;
  uint64   nWeak      = 0;
  uint64   nMissing   = 0;

  for (uint32 rr=0; rr < fiLimit; rr++) {
    nOverlaps  += _overlapLen[rr];
    nNonSymErr += nNonSymPerRead[rr];
    nWeak      += nFiltPerRead[rr];
    nMissing   += nMissPerRead[rr];
  }

  writeStatus("OverlapCache()--   In %llu overlaps:\n", nOverlaps);
  writeStatus("OverlapCache()--     Found %llu overlaps with non-symmetric error rates.\n", nNonSymErr);
  writeStatus("OverlapCache()--     Found %llu overlaps with missing twins.\n", nMissing);
  writeStatus("OverlapCache()--     Dropped %llu weak missing-twin overlaps.\n", nWeak);

  //
  //  Expand or shrink space for the overlaps.
  //

  //  Allocate new temporary pointers for each read.

  //  The new storage must start after the old storage.  And if it starts after the old storage ends,
  //  we can copy easier.  If not, we just grab some empty overlaps to make space.

  //  A complication occurs at the end of a single segment.  If there isn't enough space in the
  //  current segment for the overlaps, we skip ahead to the next segment without accounting for the
  //  overlaps we skip.  It's possible for the new size to fit into this unused space, which would
  //  then put the old overlaps physically after the new ones.
  //
  //  [ olaps1+unused   | olaps2+unused   | olaps3+unused    |        ] [ olaps4+unused   | ..... ]
  //  [ olaps1+new      | olaps2+new      | olaps3+new | olaps4+new | ] [ olaps5+new  | ....      ]
  //
  //  So, we need to compare not overlap counts, but raw positions in the OverlapStorage object.

  BAToverlap      **nPtr = new BAToverlap * [fiLimit];            //  Pointers to new overlap storage locations
  OverlapStorage   *oldS = new OverlapStorage(_overlapStorage);   //  Recreates the existing layout without allocating anything
  OverlapStorage   *newS = _overlapStorage->reset();              //  Resets pointers for the new layout, using existing space

  for (uint32 rr=1; rr < fiLimit; rr++)
    nPtr[rr] = nullptr;

  //  Compute new pointers for overlap data that begin after the old data.

  for (uint32 rr=1; rr < fiLimit; rr++) {
    uint32  nOvl = _overlapLen[rr] + nMissPerRead[rr] - nFiltPerRead[rr];

    assert(_overlapLen[rr] + nMissPerRead[rr] >= nFiltPerRead[rr]);

    oldS->get(_overlapMax[rr]);   //  Move old storage ahead to the start of the next read.
    newS->advance(oldS);          //  Advance newS so that it is at or after where oldS is.
    nPtr[rr] = newS->get(nOvl);   //  Grab space for the number of overlaps we'll end up with on this read.

    _overlapMax[rr] = nOvl;       //  Update the maximum number of overlaps on this read.
  }

  //  Copy overlap data from the old space to the new space, backwards.
  //  THIS CANNOT BE THREADED.

  writeStatus("OverlapCache()--   Shifting overlaps.\n");

  FILE *NTD = AS_UTL_openOutputFile(_prefix, '.', "non-symmetric-weak-dropped", false);

  for (uint32 rr=fiLimit; rr-- > 0; ) {
    if (_overlapLen[rr] == 0)   //  Skip the asserts!
      continue;

    assert(_overlaps[rr][0                ].a_iid == rr);  //  Ensure that we didn't overwrite an existing
    assert(_overlaps[rr][_overlapLen[rr]-1].a_iid == rr);  //  overlap with overlaps for the reads after it.

    uint32 oo = 0;  //  Position in original overlap list
    uint32 nn = 0;  //  Position in new overlap list

    for (; oo<_overlapLen[rr]; oo++)             //  Over all the original overlaps,
      if (_overlaps[rr][oo].filtered == false)   //  Copy to new storage if they're
        nPtr[rr][nn++] = _overlaps[rr][oo];      //  not filtered.
      else
        if (NTD)
          fprintf(NTD, "DROP overlap a %u b %u\n", _overlaps[rr][oo].a_iid, _overlaps[rr][oo].b_iid);

    assert(nn == _overlapLen[rr] - nFiltPerRead[rr]);

    _overlapLen[rr] = _overlapLen[rr] - nFiltPerRead[rr];
  }

  AS_UTL_closeFile(NTD);

  //  Swap pointers to the pointers and remove the old pointer storage.

  delete [] _overlaps;
  _overlaps = nPtr;

  delete  oldS;

  //  Copy non-twin overlaps to their twin.

  writeStatus("OverlapCache()--   Adding missing twins.\n");

  FILE *NTA = AS_UTL_openOutputFile(_prefix, '.', "non-symmetric-added", false);

  //  This has several concurrency issues.
  //  1)  loop test on overlapLen[ra] can change if we increment overlapLen[rb].
  //      we could access overlaps[ra][oo] before it is copied to in another thread
  //  2)  overlapLen[rb]++
  //  3)  nMissPerRead[rb]

  //#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 ra=0; ra < fiLimit; ra++) {
    for (uint32 oo=0; oo<_overlapLen[ra]; oo++) {
      if (_overlaps[ra][oo].symmetric == true)
        continue;

      uint32  rb = _overlaps[ra][oo].b_iid;
      uint32  nn = _overlapLen[rb]++;

      assert(nn < _overlapMax[rb]);

      if (NTA)
        fprintf(NTA, "add missing twin from read %u -> read %u at pos %u out of %u\n", ra, rb, nn, _overlapMax[rb] - 1);

      _overlaps[rb][nn].evalue    =  _overlaps[ra][oo].evalue;
      _overlaps[rb][nn].a_hang    = (_overlaps[ra][oo].flipped) ? (_overlaps[ra][oo].b_hang) : (-_overlaps[ra][oo].a_hang);
      _overlaps[rb][nn].b_hang    = (_overlaps[ra][oo].flipped) ? (_overlaps[ra][oo].a_hang) : (-_overlaps[ra][oo].b_hang);
      _overlaps[rb][nn].flipped   =  _overlaps[ra][oo].flipped;

      _overlaps[rb][nn].filtered  =  _overlaps[ra][oo].filtered;
      _overlaps[rb][nn].symmetric =  _overlaps[ra][oo].symmetric = true;

      _overlaps[rb][nn].a_iid     =  _overlaps[ra][oo].b_iid;
      _overlaps[rb][nn].b_iid     =  _overlaps[ra][oo].a_iid;

      assert(ra == _overlaps[ra][oo].a_iid);
      assert(rb == _overlaps[ra][oo].b_iid);

      assert(nMissPerRead[rb] > 0);

      nMissPerRead[rb]--;
    }
  }

  AS_UTL_closeFile(NTA);

  //  Check that everything worked.

  for (uint32 rr=0; rr < fiLimit; rr++) {
    assert(nMissPerRead[rr] == 0);

    if (_overlapLen[rr] == 0)
      continue;

    assert(_overlapLen[rr] == _overlapMax[rr]);

    assert(_overlaps[rr][0                ].a_iid == rr);
    assert(_overlaps[rr][_overlapLen[rr]-1].a_iid == rr);
  }

  //  Cleanup.

  delete [] nNonSymPerRead;
  delete [] nFiltPerRead;
  delete [] nMissPerRead;

  //  Probably should sort again.  Not sure if anything depends on this.

  writeStatus("OverlapCache()--   Sorting overlaps.\n");

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 rr=0; rr < fiLimit; rr++)
    if (_overlapLen[rr] > 0)
      sort(_overlaps[rr], _overlaps[rr] + _overlapLen[rr], [](BAToverlap const &a, BAToverlap const &b) {
                                                             return(((a.b_iid == b.b_iid) && (a.flipped < b.flipped)) || (a.b_iid < b.b_iid)); } );

  //  Check that all overlaps are present.

  writeStatus("OverlapCache()--   Checking overlap symmetry.\n");

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 ra=0; ra < fiLimit; ra++) {
    for (uint32 oa=0; oa<_overlapLen[ra]; oa++) {
      BAToverlap  *ova = &_overlaps[ra][oa];
      uint32       rb  =  _overlaps[ra][oa].b_iid;

      uint32       ob  = searchForOverlap(_overlaps[rb], _overlapLen[rb], ra, ova->flipped);

      if (ob == UINT32_MAX) {
        for(uint32 ii=0; ii<_overlapLen[ra]; ii++)
          fprintf(stderr, "olapA %u -> %u flip %c%s\n", _overlaps[ra][ii].a_iid, _overlaps[ra][ii].b_iid, _overlaps[ra][ii].flipped ? 'Y' : 'N', (ii == ra) ? " **" : "");
        for(uint32 ii=0; ii<_overlapLen[rb]; ii++)
          fprintf(stderr, "olapB %u -> %u flip %c\n",   _overlaps[rb][ii].a_iid, _overlaps[rb][ii].b_iid, _overlaps[rb][ii].flipped ? 'Y' : 'N');
      }
      assert(ob != UINT32_MAX);
    }
  }

  writeStatus("OverlapCache()--   Finished.\n");
}



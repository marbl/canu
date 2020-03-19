
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

  if (_maxPer < _minPer)
    writeStatus("OverlapCache()-- Not enough memory to load the minimum number of overlaps; increase -M.\n"), exit(1);

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

  delete [] numPer;
}



uint32
OverlapCache::filterDuplicates(uint32 &no) {
  uint32   nFiltered = 0;

  for (uint32 ii=0, jj=1, dd=0; jj<no; ii++, jj++) {
    if (_ovs[ii].b_iid != _ovs[jj].b_iid)
      continue;

    //  Found duplicate B IDs.  Drop one of them.

    nFiltered++;

    //  Drop the weaker overlap.  If a tie, drop the flipped one.

    double iiSco = RI->overlapLength(_ovs[ii].a_iid, _ovs[ii].b_iid, _ovs[ii].a_hang(), _ovs[ii].b_hang()) * _ovs[ii].erate();
    double jjSco = RI->overlapLength(_ovs[jj].a_iid, _ovs[jj].b_iid, _ovs[jj].a_hang(), _ovs[jj].b_hang()) * _ovs[jj].erate();

    if (iiSco == jjSco) {             //  Hey gcc!  See how nice I was by putting brackets
      if (_ovs[ii].flipped())         //  around this so you don't get confused by the
        iiSco = 0;                    //  non-ambiguous ambiguous else clause?
      else                            //
        jjSco = 0;                    //  You're welcome.
    }

    if (iiSco < jjSco)
      dd = ii;
    else
      dd = jj;

#if 0
    writeLog("OverlapCache::filterDuplicates()-- Dropping overlap A: %9" F_U64P " B: %9" F_U64P " - %6.4f%% - %6" F_S32P " %6" F_S32P " - %s\n",
             _ovs[dd].a_iid,
             _ovs[dd].b_iid,
             _ovs[dd].a_hang(),
             _ovs[dd].b_hang(),
             _ovs[dd].erate(),
             _ovs[dd].flipped() ? "flipped" : "");
#endif

    _ovs[dd].a_iid = 0;
    _ovs[dd].b_iid = 0;
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
  uint32 ns        = 0;
  bool   beVerbose = false;

 //beVerbose = (_ovs[0].a_iid == 3514657);

  for (uint32 ii=0; ii<no; ii++) {
    _ovsSco[ii] = 0;                                //  Overlaps 'continue'd below will be filtered, even if 'no filtering' is needed.

    if ((RI->readLength(_ovs[ii].a_iid) == 0) ||    //  At least one read in the overlap is deleted
        (RI->readLength(_ovs[ii].b_iid) == 0)) {
      if (beVerbose)
        fprintf(stderr, "olap %d involves deleted reads - %u %s - %u %s\n",
                ii,
                _ovs[ii].a_iid, (RI->readLength(_ovs[ii].a_iid) == 0) ? "deleted" : "active",
                _ovs[ii].b_iid, (RI->readLength(_ovs[ii].b_iid) == 0) ? "deleted" : "active");
      continue;
    }

    if (_ovs[ii].evalue() > maxEvalue) {            //  Too noisy to care
      if (beVerbose)
        fprintf(stderr, "olap %d too noisy evalue %f > maxEvalue %f\n",
                ii, AS_OVS_decodeEvalue(_ovs[ii].evalue()), AS_OVS_decodeEvalue(maxEvalue));
      continue;
    }

    uint32  olen = RI->overlapLength(_ovs[ii].a_iid, _ovs[ii].b_iid, _ovs[ii].a_hang(), _ovs[ii].b_hang());

    if (olen < minOverlap) {                        //  Too short to care
      if (beVerbose)
        fprintf(stderr, "olap %d too short olen %u minOverlap %u\n",
                ii, olen, minOverlap);
      continue;
    }

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

  _ovs     = new ovOverlap [_ovsMax];
  _ovsSco  = new uint64    [_ovsMax];
  _ovsTmp  = new uint64    [_ovsMax];


  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {

    //  Actually load the overlaps, then detect and remove overlaps between the same pair, then
    //  filter short and low quality overlaps.

    uint32  no = ovlStore->loadOverlapsForRead(rr, _ovs, _ovsMax);   //  no == total overlaps == numOvl
    uint32  nd = filterDuplicates(no);                               //  nd == duplicated overlaps (no is decreased by this amount)
    uint32  ns = filterOverlaps(_maxEvalue, _minOverlap, no);        //  ns == acceptable overlaps

    //if (_ovs[0].a_iid == 3514657)
    //  fprintf(stderr, "Loaded %u overlaps - no %u nd %u ns %u\n", numOvl, no, nd, ns);

    //  Allocate space for the overlaps.  Allocate a multiple of 8k, assumed to be the page size.
    //
    //  If we're loading all overlaps (ns == no) we don't need to overallocate.  Otherwise, we're
    //  loading only some of them and might have to make a twin later.
    //
    //  Once allocated copy the good overlaps.

    if (ns > 0) {
      uint32  id = _ovs[0].a_iid;

      _overlapMax[id] = ns;
      _overlapLen[id] = ns;
      _overlaps[id]   = _overlapStorage->get(_overlapMax[id]);

      _memOlaps += _overlapMax[id] * sizeof(BAToverlap);

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
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  if (_checkSymmetry == false)
    return;

  uint64  nNonSymErr = 0;
  uint64  nOverlaps  = 0;
  uint64  nOnly      = 0;
  uint64  nCritical  = 0;

  uint32   *nonsymPerRead = new uint32 [RI->numReads() + 1];  //  Overlap in this read is missing it's twin

  //  For each overlap, see if the twin overlap exists.  It is tempting to skip searching if the
  //  b-read has loaded all overlaps (the overlap we're searching for must exist) but we can't.
  //  We must still mark the overlap as being symmetric.

  writeStatus("OverlapCache()--\n");
  writeStatus("OverlapCache()-- Symmetrizing overlaps.\n");
  writeStatus("OverlapCache()--   Finding missing twins.\n");

  FILE *NSE = AS_UTL_openOutputFile(_prefix, '.', "non-symmetric-error-rates");
  FILE *NTW = AS_UTL_openOutputFile(_prefix, '.', "non-symmetric-overlaps");

  if (NSE) {
    fprintf(NSE, "     aID      bID  a error b error\n");
    fprintf(NSE, "-------- --------  ------- -------\n");
  }

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 ra=0; ra<RI->numReads()+1; ra++) {
    nonsymPerRead[ra] = 0;

    for (uint32 oa=0; oa<_overlapLen[ra]; oa++) {
      uint32  rb = _overlaps[ra][oa].b_iid;

      if (_overlaps[ra][oa].symmetric == true)   //  If already marked, we're done.
        continue;

      //  Search for the twin overlap, and if found, we're done.  The twin is marked as symmetric in the function.

      uint32 ob = searchForOverlap(_overlaps[rb], _overlapLen[rb], ra, _overlaps[ra][oa].flipped);

      if (ob < UINT32_MAX) {
        _overlaps[ra][oa].symmetric = true;   //  I have a twin!
        _overlaps[rb][ob].symmetric = true;   //  My twin has a twin, me!

        if (_overlaps[ra][oa].evalue != _overlaps[rb][ob].evalue) {
          if (NSE)
            fprintf(NSE, "%8u %8u  %7.3f %7.3f\n",
                     ra, rb,
                     _overlaps[ra][oa].erate() * 100.0,
                     _overlaps[rb][ob].erate() * 100.0);

          uint64 ev = min(_overlaps[ra][oa].evalue, _overlaps[rb][ob].evalue);

          _overlaps[ra][oa].evalue = ev;
          _overlaps[rb][ob].evalue = ev;

          nNonSymErr++;
        }

        continue;
      }

      //  Didn't find a twin.  Count how many overlaps we need to create duplicates of.

      if (NTW)
        fprintf(NTW, "NO TWIN for %6u vs %6u\n",
                ra, _overlaps[ra][oa].b_iid);

      nonsymPerRead[ra]++;
    }
  }

  AS_UTL_closeFile(NTW);
  AS_UTL_closeFile(NSE);

  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {
    nOverlaps += _overlapLen[rr];
    nOnly     += nonsymPerRead[rr];

    if (_overlapLen[rr] <= _minPer)
      nCritical += nonsymPerRead[rr];
  }

  writeStatus("OverlapCache()--   Found %llu overlaps with non-symmetric error rates.\n", nNonSymErr);
  writeStatus("OverlapCache()--   Found %llu missing twins in %llu overlaps, %llu are strong.\n", nOnly, nOverlaps, nCritical);

  //  Score all the overlaps (again) and drop the lower quality ones.  We need to drop half of the
  //  non-twin overlaps, but also want to retain some minimum number.

  //  But, there are a bunch of overlaps that fall below our score threshold that are symmetric.  We
  //  need to keep these, only because figuring out which ones are 'saved' above will be a total
  //  pain in the ass.

  //  Allocate some scratch space for each thread

  uint64  **ovsScoScratch   = new uint64 * [numThreads];
  uint64  **ovsTmpScratch   = new uint64 * [numThreads];
  uint64   *nDroppedScratch = new uint64   [numThreads];

  for (uint32 tt=0; tt<numThreads; tt++) {
    ovsScoScratch[tt]   = new uint64 [_ovsMax];
    ovsTmpScratch[tt]   = new uint64 [_ovsMax];
    nDroppedScratch[tt] = 0;
  }

  writeStatus("OverlapCache()--   Dropping weak non-twin overlaps; allocated " F_U64 " MB scratch space.\n",
              ((2 * sizeof(uint64 *) + sizeof(uint64)) * numThreads) >> 20);

  //  As advertised, score all the overlaps and drop the weak ones.

  double  fractionToDrop = 0.6;

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {

    if (_overlapLen[rr] <= _minPer)  //  If already too few overlaps, leave them all as is.
      continue;

    uint64 *ovsSco   = ovsScoScratch[omp_get_thread_num()];
    uint64 *ovsTmp   = ovsTmpScratch[omp_get_thread_num()];
    uint64 &nDropped = nDroppedScratch[omp_get_thread_num()];

    for (uint32 oo=0; oo<_overlapLen[rr]; oo++) {
      ovsSco[oo]   = RI->overlapLength( _overlaps[rr][oo].a_iid, _overlaps[rr][oo].b_iid, _overlaps[rr][oo].a_hang, _overlaps[rr][oo].b_hang);
      ovsSco[oo] <<= AS_MAX_EVALUE_BITS;
      ovsSco[oo]  |= (~_overlaps[rr][oo].evalue) & ERR_MASK;
      ovsSco[oo] <<= SALT_BITS;
      ovsSco[oo]  |= oo & SALT_MASK;

      ovsTmp[oo] = ovsSco[oo];
    }

    sort(ovsTmp, ovsTmp + _overlapLen[rr]);

    uint32  minIdx   = (uint32)floor(nonsymPerRead[rr] * fractionToDrop);

    if (minIdx < _minPer)
      minIdx = _minPer;

    uint64  minScore = ovsTmp[minIdx];

    for (uint32 oo=0; oo<_overlapLen[rr]; oo++) {
      if ((ovsSco[oo] < minScore) && (_overlaps[rr][oo].symmetric == false)) {
        nDropped++;
        _overlapLen[rr]--;
        _overlaps[rr][oo] = _overlaps[rr][_overlapLen[rr]];
        ovsSco       [oo] = ovsSco       [_overlapLen[rr]];
        oo--;
      }
    }

    for (uint32 oo=0; oo<_overlapLen[rr]; oo++)
      if (_overlaps[rr][oo].symmetric == false)
        assert(minScore <= ovsSco[oo]);
  }

  //  Are we sane?

  for (uint32 rr=RI->numReads()+1; rr-- > 0; )
    if (_overlapLen[rr] > 0) {
      assert(_overlaps[rr][0                ].a_iid == rr);
      assert(_overlaps[rr][_overlapLen[rr]-1].a_iid == rr);
    }

  //  Cleanup and log results.

  uint64  nDropped = 0;

  for (uint32 ii=0; ii<numThreads; ii++)
    nDropped += nDroppedScratch[ii];


  delete [] nonsymPerRead;
  nonsymPerRead = NULL;

  for (uint32 tt=0; tt<numThreads; tt++) {
    delete [] ovsScoScratch[tt];
    delete [] ovsTmpScratch[tt];
  }

  delete [] ovsScoScratch;
  delete [] ovsTmpScratch;
  delete [] nDroppedScratch;

  writeStatus("OverlapCache()--   Dropped %llu overlaps; scratch space released.\n", nDropped);

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

  writeStatus("OverlapCache()--   Adding %llu missing twin overlaps.\n", nToAdd);

  //
  //  Expand or shrink space for the overlaps.
  //

  //  Allocate new temporary pointers for each read.

  BAToverlap  **nPtr = new BAToverlap * [RI->numReads()+1];

  memset(nPtr, 0, sizeof(BAToverlap *) * (RI->numReads()+1));

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

  OverlapStorage   *oldS = new OverlapStorage(_overlapStorage);  //  Recreates the existing layout without allocating anything
  OverlapStorage   *newS = _overlapStorage;                      //  Resets pointers for the new layout, using existing space

  newS->reset();

  for (uint32 rr=1; rr<RI->numReads()+1; rr++) {
    nPtr[rr] = newS->get(_overlapLen[rr] + toAddPerRead[rr]);     //  Grab the pointer to the new space

    oldS->get(_overlapMax[rr]);                                   //  Move old storages ahead

    newS->advance(oldS);                                          //  Ensure newS is not before where oldS is.

    _overlapMax[rr] = _overlapLen[rr] + toAddPerRead[rr];
  }

  //  With new pointers in hand, copy overlap data - backwards - to the new locations.
  //  (Remeber that the reads are 1..numReads(), not 0..numReads()-1)

  for (uint32 rr=RI->numReads()+1; rr-- > 0; ) {
    if (_overlapLen[rr] == 0)
      continue;

    assert(_overlaps[rr][0                ].a_iid == rr);
    assert(_overlaps[rr][_overlapLen[rr]-1].a_iid == rr);

    for (uint32 oo=_overlapLen[rr]; oo-- > 0; )
      nPtr[rr][oo] = _overlaps[rr][oo];

    assert(_overlaps[rr][0                ].a_iid == rr);
    assert(_overlaps[rr][_overlapLen[rr]-1].a_iid == rr);
  }

  //  Swap pointers to the pointers and cleanup.

  delete [] _overlaps;
  _overlaps = nPtr;

  delete  oldS;
  //      newS is the original _overlapStorage, which we could delete, we'd just lose all the overlaps.

  //  Copy non-twin overlaps to their twin.
  //
  //  This cannot (easily) be parallelized.  We're iterating over overlaps in read rr, but inserting
  //  overlaps into read rb.

  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {
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

  //  Check that everything worked.

  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {
    assert(toAddPerRead[rr] == 0);

    if (_overlapLen[rr] == 0)
      continue;

    assert(_overlaps[rr][0                ].a_iid == rr);
    assert(_overlaps[rr][_overlapLen[rr]-1].a_iid == rr);
  }

  //  Cleanup.

  delete [] toAddPerRead;
  toAddPerRead = NULL;

  //  Probably should sort again.  Not sure if anything depends on this.

  for (uint32 rr=0; rr<RI->numReads()+1; rr++) {
  }

  writeStatus("OverlapCache()--   Finished.\n");
}



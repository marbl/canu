
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

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"  //  sizeof(BestEdgeOverlap)
#include "AS_BAT_Unitig.H"            //  sizeof(ufNode)
#include "AS_BAT_Logging.H"

#include "memoryMappedFile.H"

#include <sys/types.h>

uint64  ovlCacheMagic = 0x65686361436c766fLLU;  //0102030405060708LLU;

#ifndef __CYGWIN__
   #ifndef _WIN32
      #include <sys/sysctl.h>
   #endif
#endif

#ifdef HW_PHYSMEM

uint64
getMemorySize(void) {
  uint64  physMemory = 0;

  int     mib[2] = { CTL_HW, HW_PHYSMEM };
  size_t  len    = sizeof(uint64);

  errno = 0;

  if (sysctl(mib, 2, &physMemory, &len, NULL, 0) != 0)
    //  failed to get memory size, so what?
    fprintf(stderr, "sysctl() failed to return CTL_HW, HW_PHYSMEM: %s\n", strerror(errno)), exit(1);

  if (len != sizeof(uint64)) {
#ifdef HW_MEMSIZE
    mib[1] = HW_MEMSIZE;
    len = sizeof(uint64);
    if (sysctl(mib, 2, &physMemory, &len, NULL, 0) != 0 || len != sizeof(uint64))
#endif
       //  wasn't enough space, so what?
       fprintf(stderr, "sysctl() failed to return CTL_HW, HW_PHYSMEM: %s\n", strerror(errno)), exit(1);
  }

  return(physMemory);
}

#else

uint64
getMemorySize(void) {
  uint64  physPages  = sysconf(_SC_PHYS_PAGES);
  uint64  pageSize   = sysconf(_SC_PAGESIZE);
  uint64  physMemory = physPages * pageSize;

  fprintf(stderr, "PHYS_PAGES = "F_U64"\n", physPages);
  fprintf(stderr, "PAGE_SIZE  = "F_U64"\n", pageSize);
  fprintf(stderr, "MEMORY     = "F_U64"\n", physMemory);

  return(physMemory);
}

#endif



OverlapCache::OverlapCache(ovStore *ovlStoreUniq,
                           ovStore *ovlStoreRept,
                           const char *prefix,
                           double erate,
                           uint32 minOverlap,
                           uint64 memlimit,
                           uint32 maxOverlaps,
                           bool onlySave,
                           bool doSave) {

  _memLimit      = 0;
  _memUsed       = 0;

  _storMax       = 0;
  _storLen       = 0;
  _stor          = NULL;

  _heaps.clear();

  _cacheMMF      = NULL;

  _cachePtr      = NULL;
  _cacheLen      = NULL;

  _maxPer        = 0;

  _ovsMax        = 0;
  _ovs           = NULL;
  _ovsSco        = NULL;
  _ovsTmp        = NULL;

  _threadMax     = 0;
  _thread        = NULL;

  _ovlStoreUniq  = NULL;
  _ovlStoreRept  = NULL;

  if (load(prefix, erate) == true)
    return;

  fprintf(stderr, "\n");

  if (memlimit == UINT64_MAX) {
    _memLimit = getMemorySize();
    fprintf(stderr, "OverlapCache()-- limited to "F_U64"MB memory (total physical memory).\n", _memLimit >> 20);
  } else if (memlimit > 0) {
    _memLimit = memlimit;
    fprintf(stderr, "OverlapCache()-- limited to "F_U64"MB memory (user supplied).\n", _memLimit >> 20);
  } else {
    fprintf(stderr, "OverlapCache()-- using unlimited memory (-M 0).\n");
    _memLimit = UINT64_MAX;
  }

  //  Need to initialize thread data before we can account for their size.
  _threadMax = omp_get_max_threads();
  _thread    = new OverlapCacheThreadData [_threadMax];

  //  And this too.
  _ovsMax  = 1 * 1024 * 1024;  //  At 16B each, this is 16MB

  //  Account for memory used by fragment data, best overlaps, and unitigs.
  //  The chunk graph is temporary, and should be less than the size of the unitigs.

  uint64 memFI = FI->memoryUsage();
  uint64 memBE = FI->numFragments() * sizeof(BestEdgeOverlap);
  uint64 memUL = FI->numFragments() * sizeof(ufNode);             //  For fragment positions in unitigs
  uint64 memUT = FI->numFragments() * sizeof(uint32) / 16;        //  For unitigs (assumes 32 frag / unitig)
  uint64 memID = FI->numFragments() * sizeof(uint32) * 2;         //  For maps of fragment id to unitig id
  uint64 memEP = FI->numFragments() * Unitig::epValueSize() * 2;  //  For error profile

  uint64 memC1 = (FI->numFragments() + 1) * (sizeof(BAToverlapInt *) + sizeof(uint32));
  uint64 memC2 = _ovsMax * (sizeof(ovOverlap) + sizeof(uint64) + sizeof(uint64));
  uint64 memC3 = _threadMax * _thread[0]._batMax * sizeof(BAToverlap);
  uint64 memC4 = (FI->numFragments() + 1) * sizeof(uint32);

  uint64 memOS = (_memLimit == getMemorySize()) ? (0.1 * getMemorySize()) : 0.0;

  uint64 memTT = memFI + memBE + memUL + memUT + memID + memC1 + memC2 + memC3 + memC4 + memOS;

  if (onlySave) {
    fprintf(stderr, "OverlapCache()-- Only saving overlaps, not computing unitigs.\n");
    memBE = 0;
    memUL = 0;
    memUT = 0;
    memID = 0;
    memTT = memFI + memBE + memUL + memUT + memID + memOS + memC1 + memC2 + memC3 + memC4;
  }

  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for fragment data.\n",                  memFI >> 20);
  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for best edges.\n",                     memBE >> 20);
  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for unitig layouts.\n",                 memUL >> 20);
  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for unitigs.\n",                        memUT >> 20);
  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for id maps.\n",                        memID >> 20);
  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for error profiles.\n",                 memEP >> 20);
  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for overlap cache pointers.\n",         memC1 >> 20);
  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for overlap cache initial bucket.\n",   memC2 >> 20);
  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for overlap cache thread data.\n",      memC3 >> 20);
  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for number of overlaps per read.\n",    memC4 >> 20);
  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for other processes.\n",                memOS >> 20);
  fprintf(stderr, "OverlapCache()-- ---------\n");
  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB for data structures (sum of above).\n", memTT >> 20);

  if (_memLimit <= memTT) {
    int64 defecit = (int64)memTT - (int64)_memLimit;

    fprintf(stderr, "OverlapCache()-- %7"F_S64P"MB available for overlaps.\n", defecit);
    fprintf(stderr, "OverlapCache()--  Out of memory before loading overlaps; increase -M.\n");
    exit(1);
  }

  _memLimit -= memTT;
  _memUsed   = 0;

  fprintf(stderr, "OverlapCache()-- %7"F_U64P"MB available for overlaps.\n",             _memLimit >> 20);
  fprintf(stderr, "\n");

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

  _cacheMMF = NULL;

  _cachePtr = new BAToverlapInt * [FI->numFragments() + 1];
  _cacheLen = new uint32          [FI->numFragments() + 1];

  memset(_cachePtr, 0, sizeof(BAToverlapInt *) * (FI->numFragments() + 1));
  memset(_cacheLen, 0, sizeof(uint32)          * (FI->numFragments() + 1));

  _maxPer  = maxOverlaps;

  _ovs     = ovOverlap::allocateOverlaps(NULL, _ovsMax);  //  So can't call bgn or end.
  _ovsSco  = new uint64     [_ovsMax];
  _ovsTmp  = new uint64     [_ovsMax];

  _ovlStoreUniq = ovlStoreUniq;
  _ovlStoreRept = ovlStoreRept;

  assert(_ovlStoreUniq != NULL);
  assert(_ovlStoreRept == NULL);

  if (_memUsed > _memLimit)
    fprintf(stderr, "OverlapCache()-- ERROR: not enough memory to load ANY overlaps.\n"), exit(1);

  computeOverlapLimit();
  loadOverlaps(erate, minOverlap, prefix, onlySave, doSave);

  delete [] _ovs;       _ovs    = NULL;
  delete [] _ovsSco;    _ovsSco = NULL;
  delete [] _ovsTmp;    _ovsTmp = NULL;

  if (doSave == true)
    save(prefix, erate);

  if ((doSave == true) && (onlySave == true))
    fprintf(stderr, "Exiting; only requested to build the overlap graph.\n"), exit(0);
}


OverlapCache::~OverlapCache() {

  if (_cacheMMF) {
    _stor = NULL;
    delete _cacheMMF;
  }

  delete [] _ovs;

  delete [] _thread;

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

  _ovlStoreUniq->resetRange();

  //  AS_OVS_numOverlapsPerFrag returns an array that starts at firstIIDrequested.  This is usually
  //  1, unless the first fragment has no overlaps.  In that case, firstIIDrequested will be the
  //  first fragment with overlaps.  This is a terrible interface.

  fprintf(stderr, "OverlapCache()-- Loading number of overlaps per fragment.\n");

  uint32  frstFrag  = 0;
  uint32  lastFrag  = 0;
  uint32 *numPer    = _ovlStoreUniq->numOverlapsPerFrag(frstFrag, lastFrag);
  uint32  totlFrag  = lastFrag - frstFrag + 1;
  uint32  numPerMax = 0;

  for (uint32 i=0; i<totlFrag; i++)
    if (numPerMax < numPer[i])
      numPerMax = numPer[i];

  _maxPer = (_memLimit - _memUsed) / (FI->numFragments() * sizeof(BAToverlapInt));

  fprintf(stderr, "OverlapCache()--  Initial guess at _maxPer="F_U32" (max of "F_U32") from (memLimit="F_U64" - memUsed="F_U64") / (numFrags="F_U32" * sizeof(OVL)="F_SIZE_T")\n",
          _maxPer, numPerMax, _memLimit, _memUsed, FI->numFragments(), sizeof(BAToverlapInt));

  if (_maxPer < 10)
    fprintf(stderr, "OverlapCache()-- ERROR: not enough memory to load overlaps (_maxPer="F_U32" < 10).\n", _maxPer), exit(1);

  uint64  totalLoad  = 0;  //  Total overlaps we would load at this threshold

  uint32  numBelow   = 0;  //  Number below the threshold
  //uint64  numBelowS  = 0;  //  Amount of space wasted beacuse of this
  uint32  numEqual   = 0;
  uint32  numAbove   = 0;  //  Number of fragments above the threshold

  uint32  lastMax    = 0;

  uint32  adjust     = 1;

  while (adjust > 0) {
    totalLoad = 0;
    numBelow  = 0;
    //numBelowS = 0;
    numEqual  = 0;
    numAbove  = 0;

    for (uint32 i=0; i<totlFrag; i++) {
      if (numPer[i] < _maxPer) {
        numBelow++;
        //numBelowS  += _maxPer - MAX(lastMax, numPer[i]);  //  Number of extra overlaps we could still load; the unused space for this read
        totalLoad  += numPer[i];

      } else if (numPer[i] == _maxPer) {
        numEqual++;
        totalLoad  += _maxPer;

      } else {
        numAbove++;
        totalLoad  += _maxPer;
      }
    }

    fprintf(stderr, "OverlapCache()-- _maxPer=%7"F_U32P" (numBelow="F_U32" numEqual="F_U32" numAbove="F_U32" totalLoad="F_U64" -- "F_U64" + "F_U64" = "F_U64" <? "F_U64"\n",
            _maxPer, numBelow, numEqual, numAbove,
            totalLoad, _memUsed, totalLoad + _memUsed,
            totalLoad * sizeof(BAToverlapInt), _memLimit);


    if ((numAbove == 0) && (_memUsed + totalLoad * sizeof(BAToverlapInt) < _memLimit)) {
      //  All done, nothing to do here.
      adjust = 0;

    } else if (_memUsed + totalLoad * sizeof(BAToverlapInt) < _memLimit) {
      //  This limit worked, let's try moving it a little higher.

      lastMax  = _maxPer;

      adjust   = (_memLimit - _memUsed - totalLoad * sizeof(BAToverlapInt)) / numAbove / sizeof(BAToverlapInt);
      _maxPer += adjust;

      fprintf(stderr, "OverlapCache()--                 ("F_U64" MB free, adjust by "F_U32")\n",
              (_memLimit - _memUsed - totalLoad * sizeof(BAToverlapInt)) >> 20,
              adjust);

      if (_maxPer > numPerMax)
        _maxPer = numPerMax;

    } else {
      //  Whoops!  Too high!  Revert to the last and recompute statistics.

      adjust    = 0;
      _maxPer   = lastMax;

      totalLoad = 0;
      numBelow  = 0;
      numEqual  = 0;
      numAbove  = 0;

      for (uint32 i=0; i<totlFrag; i++) {
        if (numPer[i] < _maxPer) {
          numBelow++;
          //numBelowS  += _maxPer - numPer[i];
          totalLoad  += numPer[i];

        } else if (numPer[i] == _maxPer) {
          numEqual++;
          totalLoad  += _maxPer;

        } else {
          numAbove++;
          totalLoad  += _maxPer;
        }
      }

      fprintf(stderr, "OverlapCache()-- _maxPer=%7"F_U32P" (overestimated, revert to last good and stop)\n", _maxPer);
    }
  }

  //  Report

  fprintf(stderr, "\n");
  fprintf(stderr, "OverlapCache()-- blockSize        = "F_U32" ("F_SIZE_T"MB)\n", _storMax, (_storMax * sizeof(BAToverlapInt)) >> 20);
  fprintf(stderr, "\n");
  fprintf(stderr, "OverlapCache()-- _maxPer          = "F_U32" overlaps/reads\n", _maxPer);
  fprintf(stderr, "OverlapCache()-- numBelow         = "F_U32" reads (all overlaps loaded)\n", numBelow);
  fprintf(stderr, "OverlapCache()-- numEqual         = "F_U32" reads (all overlaps loaded)\n", numEqual);
  fprintf(stderr, "OverlapCache()-- numAbove         = "F_U32" reads (some overlaps loaded)\n", numAbove);
  fprintf(stderr, "OverlapCache()-- totalLoad        = "F_U64" overlaps (%6.2f%%)\n", totalLoad, 100.0 * totalLoad / _ovlStoreUniq->numOverlapsInRange());
  fprintf(stderr, "\n");
  fprintf(stderr, "OverlapCache()-- availForOverlaps = "F_U64"MB\n", _memLimit >> 20);
  fprintf(stderr, "OverlapCache()-- totalMemory      = "F_U64"MB for organization\n", _memUsed >> 20);
  fprintf(stderr, "OverlapCache()-- totalMemory      = "F_U64"MB for overlaps\n", (totalLoad * sizeof(BAToverlapInt)) >> 20);
  fprintf(stderr, "OverlapCache()-- totalMemory      = "F_U64"MB used\n", (_memUsed + totalLoad * sizeof(BAToverlapInt)) >> 20);
  fprintf(stderr, "\n");

  delete [] numPer;
}





uint32
OverlapCache::filterOverlaps(uint32 maxEvalue, uint32 minOverlap, uint32 no) {
  uint32 ns = 0;

  //  Score the overlaps.

  uint64  ERR_MASK = ((uint64)1 << AS_MAX_EVALUE_BITS) - 1;

  uint32  SALT_BITS = (64 - AS_MAX_READLEN_BITS - AS_MAX_EVALUE_BITS);
  uint64  SALT_MASK = (((uint64)1 << SALT_BITS) - 1);

  memset(_ovsSco, 0, sizeof(uint64) * no);

  for (uint32 ii=0; ii<no; ii++) {
    if ((FI->fragmentLength(_ovs[ii].a_iid) == 0) ||
        (FI->fragmentLength(_ovs[ii].b_iid) == 0))
      //  At least one read deleted in the overlap
      continue;

    if (_ovs[ii].evalue() > maxEvalue)
      //  Too noisy.
      continue;

    uint32  olen = FI->overlapLength(_ovs[ii].a_iid, _ovs[ii].b_iid, _ovs[ii].a_hang(), _ovs[ii].b_hang());

    if (olen < minOverlap)
      //  Too short.
      continue;

    //  Just right!

    _ovsSco[ii]   = olen;
    _ovsSco[ii] <<= AS_MAX_EVALUE_BITS;
    _ovsSco[ii]  |= (~_ovs[ii].evalue()) & ERR_MASK;
    _ovsSco[ii] <<= SALT_BITS;
    _ovsSco[ii]  |= ii & SALT_MASK;
    ns++;
  }

  //  If fewer than the limit, keep them all.  Should we reset ovsSco to be 1?  Do we really need ovsTmp?

  memcpy(_ovsTmp, _ovsSco, sizeof(uint64) * no);

  if (ns <= _maxPer)
    return(ns);

  //  Otherwise, filter out the short and low quality.

  sort(_ovsTmp, _ovsTmp + no);

  uint64  cutoff = _ovsTmp[no - _maxPer];

  for (uint32 ii=0; ii<no; ii++)
    if (_ovsSco[ii] < cutoff)
      _ovsSco[ii] = 0;

  //  Count how many overlaps we saved.

  ns = 0;

  for (uint32 ii=0; ii<no; ii++)
    if (_ovsSco[ii] > 0)
      ns++;

  if (ns > _maxPer)
    fprintf(stderr, "WARNING: fragment "F_U32" loaded "F_U32" overlas (it has "F_U32" in total); over the limit of "F_U32"\n",
            _ovs[0].a_iid, ns, no, _maxPer);

  return(ns);
}




void
OverlapCache::loadOverlaps(double erate, uint32 minOverlap, const char *prefix, bool onlySave, bool doSave) {
  uint64   numTotal     = 0;
  uint64   numLoaded    = 0;
  uint32   numFrags     = 0;
  uint32   numOvl       = 0;
  uint32   maxEvalue    = AS_OVS_encodeEvalue(erate);

  FILE    *ovlDat = NULL;

  if (doSave == true) {
    char     name[FILENAME_MAX];

    sprintf(name, "%s.ovlCacheDat", prefix);

    fprintf(stderr, "OverlapCache()-- Saving overlaps to '%s'.\n", name);

    errno = 0;

    ovlDat = fopen(name, "w");
    if (errno)
      fprintf(stderr, "OverlapCache()-- Failed to open '%s' for write: %s\n", name, strerror(errno)), exit(1);
  }

  assert(_ovlStoreUniq != NULL);
  assert(_ovlStoreRept == NULL);

  _ovlStoreUniq->resetRange();

  uint64 numStore = _ovlStoreUniq->numOverlapsInRange();

  writeLog("OverlapCache()-- Loading overlap information\n");

  //  Could probably easily extend to multiple stores.  Needs to interleave the two store
  //  loads, can't do one after the other as we require all overlaps for a single fragment
  //  be in contiguous memory.

  while (1) {

    //  Ask the store how many overlaps exist for this fragment.
    numOvl = _ovlStoreUniq->numberOfOverlaps();

    numTotal += numOvl;

    if (numOvl == 0)
      //  No overlaps?  We're at the end of the store.
      break;

    //  Resize temporary storage space to hold all these overlaps.
    while (_ovsMax <= numOvl) {
      _memUsed -= (_ovsMax) * sizeof(ovOverlap);
      _memUsed -= (_ovsMax) * sizeof(uint64);
      _ovsMax *= 2;
      delete [] _ovs;
      delete [] _ovsSco;
      delete [] _ovsTmp;
      _ovs    = ovOverlap::allocateOverlaps(NULL, _ovsMax);  //  So can't call bgn or end.
      _ovsSco = new uint64     [_ovsMax];
      _ovsTmp = new uint64     [_ovsMax];
      _memUsed += (_ovsMax) * sizeof(ovOverlap);
      _memUsed += (_ovsMax) * sizeof(uint64);
      _memUsed += (_ovsMax) * sizeof(uint64);
    }

    //  Actually load the overlaps.
    uint32  no = _ovlStoreUniq->readOverlaps(_ovs, _ovsMax);
    uint32  ns = filterOverlaps(maxEvalue, minOverlap, no);

    //  Resize the permament storage space for overlaps.
    if ((_storLen + ns > _storMax) ||
        (_stor == NULL)) {

      if ((ovlDat) && (_storLen > 0))
        AS_UTL_safeWrite(ovlDat, _stor, "_stor", sizeof(BAToverlapInt), _storLen);
      if (onlySave)
        delete [] _stor;

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

      _stor[_storLen].evalue  = _ovs[ii].evalue();
      _stor[_storLen].a_hang  = _ovs[ii].a_hang();
      _stor[_storLen].b_hang  = _ovs[ii].b_hang();
      _stor[_storLen].flipped = _ovs[ii].flipped();
      _stor[_storLen].b_iid   = _ovs[ii].b_iid;

      _storLen++;
    }

    assert(storEnd == _storLen);

    if ((numFrags++ % 1000000) == 0)
      writeLog("OverlapCache()-- Loading overlap information: overlaps processed %12"F_U64P" (%06.2f%%) loaded %12"F_U64P" (%06.2f%%) (at read iid %d)\n",
               numTotal,  100.0 * numTotal  / numStore,
               numLoaded, 100.0 * numLoaded / numStore,
               _ovs[0].a_iid);
  }

  if ((ovlDat) && (_storLen > 0))
    AS_UTL_safeWrite(ovlDat, _stor, "_stor", sizeof(BAToverlapInt), _storLen);
  if (onlySave)
    delete [] _stor;

  if (ovlDat)
    fclose(ovlDat);

  writeLog("OverlapCache()-- Loading overlap information: overlaps processed %12"F_U64P" (%06.2f%%) loaded %12"F_U64P" (%06.2f%%)\n",
           numTotal,  100.0 * numTotal  / numStore,
           numLoaded, 100.0 * numLoaded / numStore);
}




BAToverlap *
OverlapCache::getOverlaps(uint32 fragIID, double maxErate, uint32 &numOverlaps) {
  uint32 tid = omp_get_thread_num();

  while (_thread[tid]._batMax <= _cacheLen[fragIID]) {
    _thread[tid]._batMax *= 2;
    delete [] _thread[tid]._bat;
    _thread[tid]._bat = new BAToverlap [_thread[tid]._batMax];
  }

  BAToverlapInt *ptr       = _cachePtr[fragIID];
  uint32         maxEvalue = AS_OVS_encodeEvalue(maxErate);

  numOverlaps = 0;

  for (uint32 pos=0; pos < _cacheLen[fragIID]; pos++) {
    if (ptr[pos].evalue > maxEvalue)
      continue;

    _thread[tid]._bat[numOverlaps].a_hang   = ptr[pos].a_hang;
    _thread[tid]._bat[numOverlaps].b_hang   = ptr[pos].b_hang;

    _thread[tid]._bat[numOverlaps].flipped  = ptr[pos].flipped;

    _thread[tid]._bat[numOverlaps].evalue   = ptr[pos].evalue;
    _thread[tid]._bat[numOverlaps].erate    = AS_OVS_decodeEvalue(ptr[pos].evalue);

    _thread[tid]._bat[numOverlaps].a_iid    = fragIID;
    _thread[tid]._bat[numOverlaps].b_iid    = ptr[pos].b_iid;

    numOverlaps++;
  }

  return(_thread[tid]._bat);
}




void
OverlapCache::removeWeakOverlaps(uint32 *minEvalue5p,
                                 uint32 *minEvalue3p) {

  uint32  fiLimit    = FI->numFragments();

  uint64  saved      = 0;
  uint64  ignored    = 0;
  uint64  removed    = 0;

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32         numOverlaps = _cacheLen[fi];
    BAToverlapInt *ptr         = _cachePtr[fi];

    for (uint32 pos=0; pos < numOverlaps; pos++) {
      uint32  aiid   = fi;
      uint32  biid   = ptr[pos].b_iid;
      uint32  evalue = ptr[pos].evalue;

      //  Ignore contained overlaps.

      if (((ptr[pos].a_hang <= 0) && (ptr[pos].b_hang >= 0)) ||
          ((ptr[pos].a_hang >= 0) && (ptr[pos].b_hang <= 0))) {
        ignored++;
        continue;
      }

      //  Decide which end we need to be looking at.

      uint32  ta = 0;
      uint32  tb = 0;

      if (ptr[pos].a_hang > 0)
        ta = minEvalue3p[aiid];
      else
        ta = minEvalue5p[aiid];

      if (ptr[pos].flipped == false) {
        if (ptr[pos].b_hang > 0)
          tb = minEvalue5p[biid];
        else
          tb = minEvalue3p[biid];

      } else {
        if (ptr[pos].b_hang > 0)
          tb = minEvalue3p[biid];
        else
          tb = minEvalue5p[biid];
      }

      //  If the erate is more than the threshold, 'remove' the overlap by maxing out the erate.

      if ((evalue > ta) ||
          (evalue > tb)) {
        //fprintf(stdout, "OverlapCache::removeWeakOverlaps()--  remove %7d %7d at %.3f\n", aiid, biid, AS_OVS_decodeEvalue(evalue));
        removed++;
        ptr[pos].evalue = AS_MAX_EVALUE;
      } else {
        saved++;
      }
    }
  }

  writeLog("OverlapCache::removeWeakOverlaps()--  removed  "F_U64" weak overlaps.\n", removed);
  writeLog("OverlapCache::removeWeakOverlaps()--  ignored  "F_U64" contained overlaps.\n", ignored);
  writeLog("OverlapCache::removeWeakOverlaps()--  retained "F_U64" strong overlaps.\n", saved);
}




#if 0
double
OverlapCache::findErate(uint32 aIID, uint32 bIID) {

  for (uint32 pos=0; pos < _cacheLen[aIID]; pos++)
    if (_cachePtr[aIID][pos].b_iid == bIID)
      return(AS_OVS_decodeEvalue(_cachePtr[aIID][pos].evalue));

  for (uint32 pos=0; pos < _cacheLen[bIID]; pos++)
    if (_cachePtr[bIID][pos].b_iid == aIID)
      return(AS_OVS_decodeEvalue(_cachePtr[bIID][pos].evalue));

  return(1.0);
}
#endif





bool
OverlapCache::load(const char *prefix, double erate) {
  char     name[FILENAME_MAX];
  FILE    *file;
  size_t   numRead;

  sprintf(name, "%s.ovlCache", prefix);
  if (AS_UTL_fileExists(name, FALSE, FALSE) == false)
    return(false);

  fprintf(stderr, "OverlapCache()-- Loading graph from '%s'.\n", name);

  errno = 0;

  file = fopen(name, "r");
  if (errno)
    fprintf(stderr, "OverlapCache()-- Failed to open '%s' for reading: %s\n", name, strerror(errno)), exit(1);

  uint64   magic      = ovlCacheMagic;
  uint32   ovserrbits = AS_MAX_EVALUE_BITS;
  uint32   ovshngbits = AS_MAX_READLEN_BITS + 1;

  AS_UTL_safeRead(file, &magic,      "overlapCache_magic",      sizeof(uint64), 1);
  AS_UTL_safeRead(file, &ovserrbits, "overlapCache_ovserrbits", sizeof(uint32), 1);
  AS_UTL_safeRead(file, &ovshngbits, "overlapCache_ovshngbits", sizeof(uint32), 1);

  if (magic != ovlCacheMagic)
    fprintf(stderr, "OverlapCache()-- ERROR:  File '%s' isn't a bogart ovlCache.\n", name), exit(1);

  AS_UTL_safeRead(file, &_memLimit, "overlapCache_memLimit", sizeof(uint64), 1);
  AS_UTL_safeRead(file, &_memUsed, "overlapCache_memUsed", sizeof(uint64), 1);

  uint32 unused;  //  Former _batMax, left in for compatibility with old caches.

  AS_UTL_safeRead(file, &_maxPer, "overlapCache_maxPer", sizeof(uint32), 1);
  AS_UTL_safeRead(file, &unused, "overlapCache_batMax", sizeof(uint32), 1);

  _threadMax = omp_get_max_threads();
  _thread    = new OverlapCacheThreadData [_threadMax];

  _cachePtr = new BAToverlapInt * [FI->numFragments() + 1];
  _cacheLen = new uint32          [FI->numFragments() + 1];

  numRead = AS_UTL_safeRead(file,  _cacheLen, "overlapCache_cacheLen", sizeof(uint32), FI->numFragments() + 1);

  if (numRead != FI->numFragments() + 1)
    fprintf(stderr, "OverlapCache()-- Short read loading graph '%s'.  Fail.\n", name), exit(1);

  _ovlStoreUniq = NULL;
  _ovlStoreRept = NULL;

  fclose(file);

  //  Memory map the overlaps

  sprintf(name, "%s.ovlCacheDat", prefix);

  _cacheMMF = new memoryMappedFile(name);

  _stor     = (BAToverlapInt *)_cacheMMF->get(0);

  //  Update pointers into the overlaps

  _cachePtr[0] = _stor;
  for (uint32 fi=1; fi<FI->numFragments() + 1; fi++)
    _cachePtr[fi] = _cachePtr[fi-1] + _cacheLen[fi-1];

  bool    doCleaning = false;
  uint64  nOvl = 0;

  for (uint32 fi=1; fi<FI->numFragments() + 1; fi++) {
    nOvl += _cacheLen[fi];

    if ((FI->fragmentLength(fi) == 0) &&
        (_cacheLen[fi] > 0))
      doCleaning = true;

    if (_cacheLen[fi] == 0)
      _cachePtr[fi] = NULL;
  }

  //  For each fragment, remove any overlaps to deleted fragments.

  writeLog("OverlapCache()-- Loaded "F_U64" overlaps.\n", nOvl);

  if (doCleaning) {
    uint64   nDel = 0;
    uint64   nMod = 0;
    uint64   nOvl = 0;

    writeLog("OverlapCache()-- Freshly deleted fragments detected.  Cleaning overlaps.\n");

    char  N[FILENAME_MAX];

    sprintf(N, "%s.overlapsRemoved.log", prefix);

    errno = 0;
    FILE *F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "OverlapCache()--  Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

    for (uint32 fi=1; fi<FI->numFragments() + 1; fi++) {
      if ((FI->fragmentLength(fi) == 0) &&
          (_cacheLen[fi] > 0)) {
        nDel++;
        fprintf(F, "Removing "F_U32" overlaps from deleted deleted fragment "F_U32"\n", _cacheLen[fi], fi);
        _cachePtr[fi] = NULL;
        _cacheLen[fi] = 0;
      }

      uint32  on = 0;

      for (uint32 oi=0; oi<_cacheLen[fi]; oi++) {
        uint32  iid = _cachePtr[fi][oi].b_iid;
        bool    del = (FI->fragmentLength(iid) == 0);

        if ((del == false) &&
            (on < oi))
          _cachePtr[fi][on] = _cachePtr[fi][oi];

        if (del == false)
          on++;
      }

      if (_cacheLen[fi] != on) {
        nMod++;
        nOvl += _cacheLen[fi] - on;
        fprintf(F, "Removing "F_U32" overlaps from living fragment "F_U32"\n", _cacheLen[fi] - on, fi);
        memset(_cachePtr[fi] + on, 0xff, (_cacheLen[fi] - on) * (sizeof(BAToverlapInt)));
      }

      _cacheLen[fi] = on;
    }

    fclose(F);

    fprintf(stderr, "OverlapCache()-- Removed all overlaps from "F_U64" deleted fragments.  Removed "F_U64" overlaps from "F_U64" alive fragments.\n",
            nDel, nOvl, nMod);
  }

  return(true);
}


void
OverlapCache::save(const char *prefix, double erate) {
  char  name[FILENAME_MAX];
  FILE *file;

  sprintf(name, "%s.ovlCache", prefix);

  fprintf(stderr, "OverlapCache()-- Saving graph to '%s'.\n", name);

  errno = 0;

  file = fopen(name, "w");
  if (errno)
    fprintf(stderr, "OverlapCache()-- Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);

  uint64   magic      = ovlCacheMagic;
  uint32   ovserrbits = AS_MAX_EVALUE_BITS;
  uint32   ovshngbits = AS_MAX_READLEN_BITS + 1;

  AS_UTL_safeWrite(file, &magic,      "overlapCache_magic",      sizeof(uint64), 1);
  AS_UTL_safeWrite(file, &ovserrbits, "overlapCache_ovserrbits", sizeof(uint32), 1);
  AS_UTL_safeWrite(file, &ovshngbits, "overlapCache_ovshngbits", sizeof(uint32), 1);

  AS_UTL_safeWrite(file, &_memLimit, "overlapCache_memLimit", sizeof(uint64), 1);
  AS_UTL_safeWrite(file, &_memUsed, "overlapCache_memUsed", sizeof(uint64), 1);

  AS_UTL_safeWrite(file, &_maxPer, "overlapCache_maxPer", sizeof(uint32), 1);
  AS_UTL_safeWrite(file, &_maxPer, "overlapCache_batMax", sizeof(uint32), 1);  //  COMPATIBILITY, REMOVE

  AS_UTL_safeWrite(file,  _cacheLen, "overlapCache_cacheLen", sizeof(uint32), FI->numFragments() + 1);

  fclose(file);
}

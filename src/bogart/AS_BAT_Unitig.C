
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

#include "runtime.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#undef  SHOW_PROFILE_CONSTRUCTION
#undef  SHOW_PROFILE_CONSTRUCTION_DETAILS

void
Unitig::reverseComplement(bool doSort) {

  //  If there are contained reads, we need to sort by position to place them correctly after
  //  their containers.  If there are no contained reads, sorting can break the initial unitig
  //  building.  When two reads start at position zero, we'll exchange the order.  Initial unitig
  //  building depends on having the first read added become the last read in the unitig
  //  after reversing.

  for (uint32 fi=0; fi<ufpath.size(); fi++) {
    ufNode  *frg = &ufpath[fi];

    frg->position.bgn = getLength() - frg->position.bgn;
    frg->position.end = getLength() - frg->position.end;

    //if (frg->contained != 0)
    //  doSort = true;

    assert(frg->position.bgn >= 0);
    assert(frg->position.end >= 0);
  }

  //  We've updated the positions of everything.  Now, sort or reverse the list, and rebuild the
  //  ufpathIdx map.

  if (doSort) {
    sort();
  } else {
    std::reverse(ufpath.begin(), ufpath.end());

    for (uint32 fi=0; fi<ufpath.size(); fi++)
      _vector->registerRead(ufpath[fi].ident, _id, fi);
  }
}



void
Unitig::cleanUp(void) {

  if (ufpath.size() > 1)
    sort();

  int32   minPos = ufpath[0].position.min();

  if (minPos != 0)
    for (uint32 fi=0; fi<ufpath.size(); fi++) {
      ufpath[fi].position.bgn -= minPos;
      ufpath[fi].position.end -= minPos;
    }

  _length = 0;

  for (uint32 fi=0; fi<ufpath.size(); fi++) {          //  Could use position.max(), but since
    _length = max(_length, ufpath[fi].position.bgn);   //  it too calls max(), there's no win
    _length = max(_length, ufpath[fi].position.end);
  }
}



void
Unitig::computeArrivalRate(const char *UNUSED(prefix),
                           const char *UNUSED(label),
                           vector<int32> *hist) {

  sort();

  for (uint32 fi=0; fi<ufpath.size(); fi++) {
    ufNode     *rdA    = &ufpath[fi];
    bool        rdAfwd = (rdA->position.bgn < rdA->position.end);
    int32       rdAlo  = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
    int32       rdAhi  = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

    for (uint32 fj=1; fj<6; fj++) {
      if (fi + fj < ufpath.size()) {
        ufNode     *rdB    = &ufpath[fi+fj];
        bool        rdBfwd = (rdB->position.bgn < rdB->position.end);
        int32       rdBlo  = (rdBfwd) ? rdB->position.bgn : rdB->position.end;
        int32       rdBhi  = (rdBfwd) ? rdB->position.end : rdB->position.bgn;

        uint32  dist = rdBlo - rdAlo;

        hist[fj].push_back(dist);
      }
    }
  }
}



class epOlapDat {
public:
  epOlapDat() {
    pos   = 0;
    open  = false;
    erate = 0.0;
  };

  epOlapDat(uint32 p, bool o, float e) {
    pos    = p;
    open   = o;
    erate  = e;
  };

  bool operator<(const epOlapDat &that)     const { return(pos < that.pos); };

  uint32  pos   : 31;
  bool    open  :  1;
  float   erate;
};



void
Unitig::computeErrorProfile(const char *UNUSED(prefix), const char *UNUSED(label)) {

#ifdef SHOW_PROFILE_CONSTRUCTION
  writeLog("errorProfile()-- Find error profile for tig " F_U32 " of length " F_U32 " with " F_SIZE_T " reads.\n",
          id(), getLength(), ufpath.size());
#endif

  errorProfile.clear();
  errorProfileIndex.clear();

  //  Count the number of overlaps we need to save.  We do this, instead of growing the array,
  //  because occasionally these are big, and having two around at the same time can blow our
  //  memory.  (Arabidopsis p5 has a tig with 160,246,250 olaps == 1gb memory)

#if 0
  //  A (much) fancier version would merge the overlap detection and errorProfile compute together.
  //  Keep lists of epOlapDat for each read end (some cleverness could probably get rid of the map,
  //  if we just use the index of the read).  Before we process a new read, all data for positions
  //  before this reads start position can be processed and freed.

  map<uint32, uint32>    baseToIndex;

  uint32                *olapsMax = new uint32    [ufpath.size() * 2];
  uint32                *olapsLen = new uint32    [ufpath.size() * 2];
  epOlapDat            **olaps    = new epOlapDat [ufpath.size() * 2];
#endif

  uint64      olapsMax = 0;
  uint64      olapsLen = 0;
  epOlapDat  *olaps    = NULL;

  for (uint32 fi=0; fi<ufpath.size(); fi++) {
    ufNode     *rdA    = &ufpath[fi];
    int32       rdAlo  = rdA->position.min();
    int32       rdAhi  = rdA->position.max();

    uint32      ovlLen =  0;
    BAToverlap *ovl    =  OC->getOverlaps(rdA->ident, ovlLen);

    for (uint32 oi=0; oi<ovlLen; oi++) {
      if (id() != _vector->inUnitig(ovl[oi].b_iid))          //  Reads in different tigs?
        continue;                                            //  Don't care about this overlap.

      ufNode  *rdB    = &ufpath[ _vector->ufpathIdx(ovl[oi].b_iid) ];

      if (rdA->ident < rdB->ident)                           //  Only want to see one overlap
        continue;                                            //  for each pair.

      int32    rdBlo  = rdB->position.min();
      int32    rdBhi  = rdB->position.max();

      if ((rdAhi <= rdBlo) || (rdBhi <= rdAlo))              //  Reads in same tig but not overlapping?
        continue;                                            //  Don't care about this overlap.

      olapsMax += 2;
    }
  }

  // Scan overlaps to find those that we care about, and save their endpoints.

  olaps = new epOlapDat [olapsMax];

  for (uint32 fi=0; fi<ufpath.size(); fi++) {
    ufNode     *rdA    = &ufpath[fi];
    int32       rdAlo  = rdA->position.min();
    int32       rdAhi  = rdA->position.max();

    uint32      ovlLen =  0;
    BAToverlap *ovl    =  OC->getOverlaps(rdA->ident, ovlLen);

    for (uint32 oi=0; oi<ovlLen; oi++) {
      if (id() != _vector->inUnitig(ovl[oi].b_iid))          //  Reads in different tigs?
        continue;                                            //  Don't care about this overlap.

      ufNode  *rdB    = &ufpath[ _vector->ufpathIdx(ovl[oi].b_iid) ];

      if (rdA->ident < rdB->ident)                           //  Only want to see one overlap
        continue;                                            //  for each pair.

      int32    rdBlo  = rdB->position.min();
      int32    rdBhi  = rdB->position.max();

      if ((rdAhi <= rdBlo) || (rdBhi <= rdAlo))              //  Reads in same tig but not overlapping?
        continue;                                            //  Don't care about this overlap.

      uint32 bgn = max(rdAlo, rdBlo);
      uint32 end = min(rdAhi, rdBhi);

#ifdef SHOW_PROFILE_CONSTRUCTION_DETAILS
      writeLog("errorProfile()-- olap %5u read %7u read %7u at %9u-%9u\n",
               oi, rdA->ident, rdB->ident, bgn, end);
#endif

      olaps[olapsLen++] = epOlapDat(bgn, true,  ovl[oi].erate());  //  Save an open event,
      olaps[olapsLen++] = epOlapDat(end, false, ovl[oi].erate());  //  and a close event.
      assert(olapsLen <= olapsMax);
    }
  }

  //  Warn if too few or too many overlaps.

  if (olapsLen == 0) {
    writeLog("WARNING:  tig %u length %u nReads %u has no overlaps.\n", id(), getLength(), ufpath.size());
    for (uint32 fi=0; fi<ufpath.size(); fi++)
      writeLog("WARNING:    read %7u %7u-%-7u\n",
               ufpath[fi].ident,
               ufpath[fi].position.bgn,
               ufpath[fi].position.end);
  }

  if (olapsLen > (uint64)4 * 1024 * 1024 * 1024) {
    writeLog("WARNING:  tig %u length %u nReads %u has " F_U64 " overlaps.\n", id(), getLength(), ufpath.size(), olapsLen);
    for (uint32 fi=0; fi<ufpath.size(); fi++)
      writeLog("WARNING:    read %7u %7u-%-7u\n",
               ufpath[fi].ident,
               ufpath[fi].position.bgn,
               ufpath[fi].position.end);
  }

  //  Sort.

  std::sort(olaps, olaps + olapsLen);

  //  Convert coordinates into intervals.  Conceptually, squish out the duplicate numbers, then
  //  create an interval for every adjacent pair.  We need to add intervals for the first and last
  //  region.  And one more, for convenience, to hold the final 'close' values on intervals that
  //  extend to the end of the unitig.

  if (olapsLen == 0)                                     //  No olaps, so add an interval
    errorProfile.push_back(epValue(0, getLength()));     //  covering the whole tig

  if ((olapsLen > 0) && (olaps[0].pos != 0))             //  Olaps, but missing the first
    errorProfile.push_back(epValue(0, olaps[0].pos));    //  interval, so add it.


  stdDev<float>  curDev;

  for (uint64 bb=0, ee=0; ee<olapsLen; ee++) {
    if (olaps[bb].pos != olaps[ee].pos) {                //  A different position.
      errorProfile.push_back(epValue(olaps[bb].pos,      //  Save the current stats in a new profile entry.
                                     olaps[ee].pos,
                                     curDev.mean(),
                                     curDev.stddev()));
      bb = ee;
    }

    if (olaps[ee].open == true)                          //  Add the new overlap to our running
      curDev.insert(olaps[ee].erate);                    //  std.dev calculation.
    else
      curDev.remove(olaps[ee].erate);

    if ((ee == olapsLen - 1) &&
        (olaps[bb].pos != olaps[ee].pos)) {              //  If the last olap,
      errorProfile.push_back(epValue(olaps[bb].pos,      //  make the final profile entry
                                     olaps[ee].pos,
                                     curDev.mean(),
                                     curDev.stddev()));
    }
  }

  if ((olapsLen > 0) && (olaps[olapsLen-1].pos != getLength()))            //  Olaps, but missing the last
    errorProfile.push_back(epValue(olaps[olapsLen-1].pos, getLength()));   //  interval, so add it.

  errorProfile.push_back(epValue(getLength(), getLength()+1));   //  And one more to make life easier.

#ifdef SHOW_PROFILE_CONSTRUCTION
  writeLog("errorProfile()-- tig %u generated " F_SIZE_T " profile regions from " F_U64 " overlaps.\n", id(), errorProfile.size(), olapsLen);
#endif

  delete [] olaps;

  //  Adjust regions that have no overlaps (mean == 0) to be the average of the adjacent regions.
  //  There are always at least two elements in the profile list: one that starts at coordinate 0,
  //  and the terminating one at coordinate (len, len+1).

  for (uint64 bi=0; bi<errorProfile.size(); bi++) {
    if (errorProfile[bi].mean != 0)
      continue;

    //  Set any initial zero coverage area to the next one.
    if      (bi == 0) {
      errorProfile[bi].mean   = errorProfile[bi+1].mean;
      errorProfile[bi].stddev = errorProfile[bi+1].stddev;
    }

    //  Set intermediate ones to the average.
    else if (bi < errorProfile.size() - 2) {
      //writeLog("errorProfile()-- tig %u no overlap coverage %u-%u\n", id(), errorProfile[bi].bgn, errorProfile[bi].end);

      errorProfile[bi].mean   = (errorProfile[bi-1].mean   + errorProfile[bi+1].mean)   / 2;
      errorProfile[bi].stddev = (errorProfile[bi-1].stddev + errorProfile[bi+1].stddev) / 2;
    }

    //  Set the last two - the last real one and the terminator - to the previous one.
    else {
      errorProfile[bi].mean   = errorProfile[bi-1].mean;
      errorProfile[bi].stddev = errorProfile[bi-1].stddev;
    }
  }


  //  Build an index.
  //    bi - base we are indexing.
  //    pi - profile
  //
  for (uint32 bi=0, pi=0; bi<getLength(); bi += 1000) {
    while ((pi < errorProfile.size()) && (errorProfile[pi].end <= bi))
      pi++;

    if (pi < errorProfile.size()) {
      assert(errorProfile[pi].bgn <= bi);
      assert(bi <  errorProfile[pi].end);

      errorProfileIndex.push_back(pi);
    }
  }

  //writeLog("errorProfile()-- tig %u generated " F_SIZE_T " profile regions with " F_U64 " overlap pieces.\n",
  //         id(), errorProfile.size(), nPieces);
}


//  For the range bgn..end, returns the amount of sequence (as a fraction)
//  that has an estimated max overlap error rate above the 'erate' threshold.
//
//  For bgn..end the range of an overlap with some 'erate', then a low
//  return value would indicate that the average overlap error rate in this
//  region is lower than the supplied 'erate' - that this overlap is too noisy
//  to be placed here.  Likewise, if the return value is 1.0, then the
//  overlap 'erate' is within the same range as the other overlaps in the tig.
//
double
Unitig::overlapConsistentWithTig(double deviations,
                                 uint32 bgn, uint32 end,
                                 double erate) {
  int32  nBelow = 0;
  int32  nAbove = 0;

  assert(bgn <  end);
  assert(bgn <  getLength());
  assert(end <= getLength());

  //  If this is a singleton tig - we should only be here when finding graph edges to repeats -
  //  we've got nothing to go on, so default to 'consistent'.

  if (errorProfile.size() == 0)
    return(1.0);

  //  Coarse search to find the first index that is after our region.

#undef BINARY_SEARCH

#ifdef BINARY_SEARCH

  uint32  min = 0;
  uint32  max = errorProfile.size();
  uint32  pb  = min + (max - min) / 2;

  while ((bgn < errorProfile[pb].bgn) ||
         (errorProfile[pb].end <= bgn)) {

    if (bgn < errorProfile[pb].bgn)
      max = pb;

    if (errorProfile[pb].end <= bgn)
      min = pb;

    assert(min < max);

    pb = min + (max - min) / 2;
  }

#else

  uint32  pbi = bgn / 1000;

  if (errorProfileIndex.size() <= pbi)
    fprintf(stderr, "errorProfileIndex.size() = " F_SIZE_T " but pbi = " F_U32 "\n", errorProfileIndex.size(),pbi);
  assert(pbi < errorProfileIndex.size());

  while ((0 < pbi) && (errorProfile[errorProfileIndex[pbi]].bgn > bgn)) {
    fprintf(stderr, "BAD ESTIMATE for bgn=%u end=%u\n", bgn, end);
    pbi--;
  }

  while ((pbi < errorProfileIndex.size()) && (errorProfile[errorProfileIndex[pbi]].end <= bgn))
    pbi++;

  if (pbi == errorProfileIndex.size()) {
    //fprintf(stderr, "Fell off loop for bgn=%u end=%u last ep bgn=%u end=%u\n",
    //        bgn, end, errorProfile.back().bgn, errorProfile.back().end);
    pbi--;
  }

  //  The region pb points to will contain bgn.

  uint32 pb = errorProfileIndex[pbi];

  //fprintf(stderr, "For bgn=%u end=%u - stopped at pbi=%u errorProfile[%u] = %u-%u (1)\n",
  //        bgn, end, pbi, pb, errorProfile[pb].bgn, errorProfile[pb].end);

  //  Fine tune search to find the exact first region.

  while ((0 < pb) && (bgn < errorProfile[pb].bgn))
    pb--;
  while ((pb < errorProfile.size()) && (errorProfile[pb].end <= bgn))
    pb++;

#endif

  if ((errorProfile[pb].bgn > bgn) ||
      (bgn >=  errorProfile[pb].end))
    fprintf(stderr, "For bgn=%u end=%u - stopped at errorProfile[%u] = %u-%u BOOM\n",
            bgn, end, pb, errorProfile[pb].bgn, errorProfile[pb].end);
  assert(errorProfile[pb].bgn <= bgn);
  assert(bgn <  errorProfile[pb].end);

  //  Sum the number of bases above the supplied erate.

  uint32 pe = pb;

  while ((pe < errorProfile.size()) && (errorProfile[pe].bgn < end)) {
    if (erate <= errorProfile[pe].max(deviations))
      nAbove += errorProfile[pe].end - errorProfile[pe].bgn;
    else
      nBelow += errorProfile[pe].end - errorProfile[pe].bgn;

    pe++;
  }

  //  Adjust for the bits we overcounted in the first and last regions.

  if (pe > 0)   //  Argh.  If this read is fully in the first region (where there
    pe--;       //  is only 1x coverage) then pe==0.


  uint32  bb = bgn - errorProfile[pb].bgn;
  uint32  be = errorProfile[pe].end - end;

  assert(bgn >= errorProfile[pb].bgn);
  assert(errorProfile[pe].end >= end);

  if (erate <= errorProfile[pb].max(deviations))
    nAbove -= bb;
  else
    nBelow -= bb;

  if (erate <= errorProfile[pe].max(deviations))
    nAbove -= be;
  else
    nBelow -= be;

  assert(nAbove >= 0);
  assert(nBelow >= 0);

  return((double)nAbove / (nBelow + nAbove));
}






void
Unitig::reportErrorProfile(const char *prefix, const char *label) {
  char  N[FILENAME_MAX];
  FILE *F;

  if (logFileFlagSet(LOG_ERROR_PROFILES) == false)
    return;

  snprintf(N, FILENAME_MAX, "%s.%s.%08u.profile", prefix, label, id());

  F = AS_UTL_openOutputFile(N);

  for (uint32 ii=0; ii<errorProfile.size(); ii++)
    fprintf(F, "%u %u %.5f +- %.5f\n",
            errorProfile[ii].bgn,  errorProfile[ii].end,
            errorProfile[ii].mean, errorProfile[ii].stddev);

  AS_UTL_closeFile(F, N);

  //  Reporting the index isn't generally useful, only for debugging.

#if 0
  snprintf(N, FILENAME_MAX, "%s.%s.%08u.profile.index", prefix, label, id());

  F = AS_UTL_openOutputFile(N);

  for (uint32 ii=0; ii<errorProfileIndex.size(); ii++) {
    uint32  xx = errorProfileIndex[ii];

    fprintf(F, "index[%u] = %u -- errorProfile[] = %u-%u  %.6f +- %.6f\n",
            ii,
            xx,
            errorProfile[xx].bgn,
            errorProfile[xx].end,
            errorProfile[xx].mean,
            errorProfile[xx].stddev);
  }
  AS_UTL_closeFile(F, N);
#endif
}


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
 *    src/AS_BAT/AS_BAT_Unitig.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
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



class optPos {
public:
  optPos() {
  };
  ~optPos() {
  };

  void    set(ufNode &n) {
    ident = n.ident;
    min   = n.position.min();
    max   = n.position.max();
    fwd   = n.position.isForward();
  };

  uint32  ident;
  double  min;
  double  max;
  bool    fwd;
};



void
Unitig::optimize_initPlace(uint32        ii,
                           optPos       *op,
                           optPos       *np,
                           bool          firstPass,
                           set<uint32>  &failed,
                           bool          beVerbose) {
  uint32       fid     = ufpath[ii].ident;

  if ((firstPass == false) && (failed.count(fid) == 0))  //  If the second pass and not
    return;                                              //  failed, do nothing.

  uint32       ovlLen  = 0;
  BAToverlap  *ovl     = OC->getOverlaps(fid, ovlLen);

  //  Unlike optimize below, we don't have any previous value to use with the read length, so
  //  we don't use read length here.

  double   nmin = 0;
  double   nmax = 0;
  int32    cnt  = 0;

  //  Then process all overlaps.

  for (uint32 oo=0; oo<ovlLen; oo++) {
    uint32  uu = inUnitig (ovl[oo].b_iid);
    uint32  jj = ufpathIdx(ovl[oo].b_iid);

    if (uu != id())   //  Skip if the overlap is to a different tig.
      continue;       //  (the ufpathIdx() call is valid, but using it isn't)

    if (isOverlapping(ufpath[ii].position, ufpath[jj].position) == false)  //  Skip if the reads
      continue;                                                            //  don't overlap

    if ((firstPass) && (jj > ii))  //  We're setting initial positions, so overlaps to reads after
      continue;                    //  us aren't correct, unless we're in the 2nd pass

    //if (beVerbose)
    //  writeLog(" olap %3u - %7u %7u - SUCCESS %8d-%-8d to %8d-%-8d\n",
    //           oo, fid, ovl[oo].b_iid,
    //           ufpath[ii].position.bgn, ufpath[ii].position.end,
    //           ufpath[jj].position.bgn, ufpath[jj].position.end);

    //  Reads overlap.  Compute the position of the read using
    //  the overlap and the other read.

    nmin += (op[ii].fwd) ? (op[jj].min - ovl[oo].a_hang) : (op[jj].min + ovl[oo].b_hang);
    nmax += (op[ii].fwd) ? (op[jj].max - ovl[oo].b_hang) : (op[jj].max + ovl[oo].a_hang);
    cnt  += 1;
  }  //  over all overlaps

  //  If no overlaps found, flag this read for a second pass.  If in the second pass,
  //  not much we can do.

  if ((firstPass == true) && (cnt == 0)) {
    writeLog("Failed to find overlaps for read %u in tig %u at %d-%d (first pass)\n",
             fid, id(), ufpath[ii].position.bgn, ufpath[ii].position.end);
    failed.insert(fid);
    return;
  }

  if ((firstPass == false) && (cnt == 0)) {
    writeLog("Failed to find overlaps for read %u in tig %u at %d-%d (second pass)\n",
             fid, id(), ufpath[ii].position.bgn, ufpath[ii].position.end);
    flushLog();
  }

  assert(cnt > 0);

  op[ii].min = nmin / cnt;
  op[ii].max = nmax / cnt;

  np[ii].min = 0;
  np[ii].max = 0;

  //  The initialization above does very little to enforce read lengths, and the optimization
  //  doesn't put enough weight in the read length to make it stable.  We simply force
  //  the correct read length here.

  op[ii].max = op[ii].min + RI->readLength(ufpath[ii].ident);

  if (beVerbose)
    writeLog("INIT tig %u read %u %9.2f-%-9.2f endDiff %9.2f%s\n",
             id(), op[ii].ident,
             op[ii].min,
             op[ii].max,
             op[ii].min + RI->readLength(ufpath[ii].ident),
             (firstPass == true) ? "" : " SECONDPASS");
}



void
Unitig::optimize(const char *prefix, const char *label) {
  bool   beVerbose = false;

  if (ufpath.size() == 1)
    return;

  if (beVerbose) {
    writeLog("\n");
    writeLog("optimize()-- tig %u\n", id());
    writeLog("\n");
  }

  optPos *pp = NULL;
  optPos *op = new optPos [ufpath.size()];
  optPos *np = new optPos [ufpath.size()];

  memset(op, 0, sizeof(optPos) * ufpath.size());
  memset(np, 0, sizeof(optPos) * ufpath.size());

  for (uint32 ii=0; ii<ufpath.size(); ii++) {
    op[ii].set(ufpath[ii]);
    np[ii].set(ufpath[ii]);
  }

  //
  //  Initialize - one round using only reads/overlaps before us
  //

  op[0].min = np[0].min = 0;                                //  Seed first read at zero
  op[0].max = np[0].max = RI->readLength(ufpath[0].ident);  //  and ending at length.

  if (beVerbose)
    writeLog("INIT tig %u read %u %9.2f-%-9.2f endDiff %9.2f\n",
             id(), op[0].ident,
             op[0].min,
             op[0].max,
             op[0].min + RI->readLength(ufpath[0].ident));

  set<uint32>   failed;

  for (uint32 ii=1; ii<ufpath.size(); ii++)
    optimize_initPlace(ii, op, np, true,  failed, beVerbose);

  for (uint32 ii=1; ii<ufpath.size(); ii++)
    optimize_initPlace(ii, op, np, false, failed, beVerbose);

  //
  //  Optimize
  //

  for (uint32 iter=0; iter<25; iter++) {
    uint64  nOlapsTotal = 0;
    uint64  nOlapsUsed  = 0;

    for (uint32 ii=0; ii<ufpath.size(); ii++) {
      uint32       fid     = op[ii].ident;

      uint32       readLen = RI->readLength(fid);

      uint32       ovlLen  = 0;
      BAToverlap  *ovl     = OC->getOverlaps(fid, ovlLen);

      double       nmin = 0.0;
      double       nmax = 0.0;
      uint32       cnt  = 0;

      if (beVerbose) {
        writeLog("optimize()-- tig %8u read %8u previous  - %9.2f-%-9.2f\n", id(), fid, op[ii].min,           op[ii].max);
        writeLog("optimize()-- tig %8u read %8u length    - %9.2f-%-9.2f\n", id(), fid, op[ii].max - readLen, op[ii].min + readLen);
      }

      //  Process all overlaps.

      for (uint32 oo=0; oo<ovlLen; oo++) {
        uint32  uu = inUnitig (ovl[oo].b_iid);
        uint32  jj = ufpathIdx(ovl[oo].b_iid);

        if (uu != id())   //  Skip if the overlap is to a different tig.
          continue;       //  (the ufpathIdx() call is valid, but using it isn't)

        if (isOverlapping(ufpath[ii].position, ufpath[jj].position) == false)  //  Skip if the reads
          continue;                                                            //  don't overlap

        //  Reads overlap.  Compute the position of the read using
        //  the overlap and the other read.

        double tmin = (op[ii].fwd) ? (op[jj].min - ovl[oo].a_hang) : (op[jj].min + ovl[oo].b_hang);
        double tmax = (op[ii].fwd) ? (op[jj].max - ovl[oo].b_hang) : (op[jj].max + ovl[oo].a_hang);

        if (beVerbose)
          writeLog("optimize()-- tig %8u read %8u olap %4u - %9.2f-%-9.2f\n", id(), fid, oo, tmin, tmax);

        nmin += tmin;
        nmax += tmax;
        cnt  += 1;
      }  //  over all overlaps

      nOlapsTotal  += ovlLen;   //  Save some stats.
      nOlapsUsed   += cnt;

      //  Add in some evidence for the bases in the read.  We want higher weight than the overlaps,
      //  but not enough to swamp the hangs.

      nmin   += cnt/4 * (op[ii].max - readLen);
      nmax   += cnt/4 * (op[ii].min + readLen);
      cnt    += cnt/4;

      //  Find the average and save.

      np[ii].min = nmin / cnt;
      np[ii].max = nmax / cnt;

      double dmin = 2 * (op[ii].min - np[ii].min) / (op[ii].min + np[ii].min);
      double dmax = 2 * (op[ii].max - np[ii].max) / (op[ii].max + np[ii].max);

      if (beVerbose)
        writeLog("optimize()-- tig %8u read %8u           - %9.2f-%-9.2f length %9.2f/%-6u posChange %+6.4f %+6.4f iter %2u\n",
                 id(), fid,
                 np[ii].min, np[ii].max,
                 np[ii].max - np[ii].min, readLen,
                 dmin, dmax,
                 iter);
    }  //  Over all reads in the tig

    //  Offset back to zero.

    int32  z = np[0].min;

    for (uint32 ii=0; ii<ufpath.size(); ii++) {
      np[ii].min -= z;
      np[ii].max -= z;
    }

    //  Decide if we've converged.  We used to compute percent difference in coordinates, but that is
    //  biased by the position of the read.  Just use percent difference from read length.

    uint32  nConverged = 0;
    uint32  nChanged   = 0;

    for (uint32 ii=0; ii<ufpath.size(); ii++) {
      double  minp = 2 * (op[ii].min - np[ii].min) / (RI->readLength(ufpath[ii].ident));
      double  maxp = 2 * (op[ii].max - np[ii].max) / (RI->readLength(ufpath[ii].ident));

      if (minp < 0)  minp = -minp;
      if (maxp < 0)  maxp = -maxp;

      if ((minp < 0.005) && (maxp < 0.005))
        nConverged++;
      else
        nChanged++;
    }

    //  All reads processed, swap op and np to compute the next iteration.

    pp = op;
    op = np;
    np = pp;

    writeLog("optimize()-- tig %8u iter %2u converged %6u changed %6u  olaps %9lu/%9lu\n",
             id(), iter, nConverged, nChanged, nOlapsUsed, nOlapsTotal);

    if (nChanged == 0)
      break;
  }

  //  Update the tig with new positions.  op[] is the result of the last iteration.

  for (uint32 ii=0; ii<ufpath.size(); ii++) {
    if (op[ii].fwd) {
      //  Logging isn't quite correct - np needs to be offset first
      if (beVerbose)
        writeLog("read %6u from %8d,%-8d to %8d,%-8d ->\n", 
                 ufpath[ii].ident,
                 ufpath[ii].position.bgn,
                 ufpath[ii].position.end,
                 (int32)op[ii].min,
                 (int32)op[ii].max);

      ufpath[ii].position.bgn = (int32)op[ii].min;
      ufpath[ii].position.end = (int32)op[ii].max;
    } else {
      if (beVerbose)
        writeLog("read %6u from %8d,%-8d to %8d,%-8d <-\n", 
                 ufpath[ii].ident,
                 ufpath[ii].position.bgn,
                 ufpath[ii].position.end,
                 (int32)op[ii].max,
                 (int32)op[ii].min);

      ufpath[ii].position.bgn = (int32)op[ii].max;
      ufpath[ii].position.end = (int32)op[ii].min;
    }
  }

  delete [] op;
  delete [] np;

  cleanUp();
}



class epOlapDat {
public:
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



void
Unitig::computeErrorProfile(const char *UNUSED(prefix), const char *UNUSED(label)) {

#ifdef SHOW_PROFILE_CONSTRUCTION
  writeLog("errorProfile()-- Find error profile for tig " F_U32 " of length " F_U32 " with " F_SIZE_T " reads.\n",
          id(), getLength(), ufpath.size());
#endif

  errorProfile.clear();
  errorProfileIndex.clear();

  vector<epOlapDat>  olaps;

  // Scan overlaps to find those that we care about, and save their endpoints.

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
      writeLog("errorProfile()-- olap[%u] %u %u begin %u end %u\n", oi, rdA->ident, rdB->ident, bgn, end);
#endif

      olaps.push_back(epOlapDat(bgn, true,  ovl[oi].erate()));  //  Save an open event,
      olaps.push_back(epOlapDat(end, false, ovl[oi].erate()));  //  and a close event.
    }
  }

  //  Warn if no overlaps.

  if (olaps.size() == 0) {
    writeLog("WARNING:  tig %u length %u nReads %u has no overlaps.\n", id(), getLength(), ufpath.size());
    for (uint32 fi=0; fi<ufpath.size(); fi++)
      writeLog("WARNING:    read %7u %7u-%-7u\n",
               ufpath[fi].ident,
               ufpath[fi].position.bgn,
               ufpath[fi].position.end);
  }

  //  Sort.

  std::sort(olaps.begin(), olaps.end());

  //  Convert coordinates into intervals.  Conceptually, squish out the duplicate numbers, then
  //  create an interval for every adjacent pair.  We need to add intervals for the first and last
  //  region.  And one more, for convenience, to hold the final 'close' values on intervals that
  //  extend to the end of the unitig.

  if (olaps.size() == 0)
    errorProfile.push_back(epValue(0, getLength()));

  if ((olaps.size() > 0) && (olaps[0].pos != 0))
    errorProfile.push_back(epValue(0, olaps[0].pos));

  for (uint32 bb=0, ii=1; ii<olaps.size(); ii++) {
    if (olaps[bb].pos == olaps[ii].pos)
      continue;

    errorProfile.push_back(epValue(olaps[bb].pos, olaps[ii].pos));

#ifdef SHOW_PROFILE_CONSTRUCTION_DETAILS
    writeLog("errorProfile()-- tig %u make region [%u-%u] @ %u-%u\n", id(), bb, ii, olaps[bb].pos, olaps[ii].pos);
#endif

    bb = ii;
  }

  if ((olaps.size() > 0) && (olaps[olaps.size()-1].pos != getLength()))
    errorProfile.push_back(epValue(olaps[olaps.size()-1].pos, getLength()));

  errorProfile.push_back(epValue(getLength(), getLength()+1));


#ifdef SHOW_PROFILE_CONSTRUCTION
  writeLog("errorProfile()-- tig %u generated " F_SIZE_T " profile regions from " F_SIZE_T " overlaps.\n", id(), errorProfile.size(), olaps.size());
#endif

  //  Walk both lists, adding positive erates and removing negative erates.

  stdDev<float>  curDev;

  for (uint32 oo=0, ee=0; oo<olaps.size(); oo++) {
    if (olaps[oo].pos != errorProfile[ee].bgn)  //  Move to the next profile if the pos is different.
      ee++;                                     //  By construction, this single step should be all we need.

#ifdef SHOW_PROFILE_CONSTRUCTION_DETAILS
    writeLog("errorProfile()-- olap[%u] @ %u  ep[%u] @ %u  %s %f  %f +- %f size %u\n",
             oo, olaps[oo].pos,
             ee, errorProfile[ee].bgn,
             olaps[oo].open ? "I" : "R",
             olaps[oo].erate,
             curDev.mean(), curDev.variance(), curDev.size());

    if ((olaps[oo].open == false) && (curDev.size() == 0)) {
      for (uint32 fi=0; fi<ufpath.size(); fi++) {
        ufNode  *frg = &ufpath[fi];
        writeLog("read %6u %6u-%6u\n", frg->ident, frg->position.bgn, frg->position.end);
      }

      writeLog("errorProfile()-- remove from empty set?\n");
      flushLog();
    }
#endif

    assert(olaps[oo].pos == errorProfile[ee].bgn);
    assert(oo < olaps.size());
    assert(ee < errorProfile.size());

    if (olaps[oo].open == true)
      curDev.insert(olaps[oo].erate);
    else
      curDev.remove(olaps[oo].erate);

    errorProfile[ee].dev = curDev;
  }

  //  Finalize the values.

  for (uint32 bi=0; bi<errorProfile.size(); bi++)
    errorProfile[bi].dev.finalize();

  //  Adjust regions that have no overlaps (mean == 0) to be the average of the adjacent regions.
  //  There are always at least two elements in the profile list: one that starts at coordinate 0,
  //  and the terminating one at coordinate (len, len+1).

  for (uint32 bi=0; bi<errorProfile.size(); bi++) {
    if (errorProfile[bi].dev.mean() != 0)
      continue;

    //  Set any initial zero coverage area to the next one.
    if      (bi == 0) {
      errorProfile[bi].dev = errorProfile[bi+1].dev;
    }

    //  Set intermediate ones to the average.
    else if (bi < errorProfile.size() - 2) {
      //writeLog("errorProfile()-- tig %u no overlap coverage %u-%u\n", id(), errorProfile[bi].bgn, errorProfile[bi].end);

      errorProfile[bi].dev = stdDev<float>((errorProfile[bi-1].dev.mean()   + errorProfile[bi+1].dev.mean()) / 2,
                                           (errorProfile[bi-1].dev.stddev() + errorProfile[bi+1].dev.stddev()) / 2,
                                           1);
    }

    //  Set the last two - the last real one and the terminator - to the previous one.
    else {
      errorProfile[bi].dev = errorProfile[bi-1].dev;
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

  F = fopen(N, "w");

  if (F) {
    for (uint32 ii=0; ii<errorProfile.size(); ii++)
      fprintf(F, "%u %u %f +- %f (%u overlaps)\n",
              errorProfile[ii].bgn,        errorProfile[ii].end,
              errorProfile[ii].dev.mean(), errorProfile[ii].dev.stddev(),
              errorProfile[ii].dev.size());
    fclose(F);
  }

  //  Reporting the index isn't generally useful, only for debugging.

#if 0
  snprintf(N, FILENAME_MAX, "%s.%s.%08u.profile.index", prefix, label, id());

  F = fopen(N, "w");

  if (F) {
    for (uint32 ii=0; ii<errorProfileIndex.size(); ii++) {
      uint32  xx = errorProfileIndex[ii];

      fprintf(F, "index[%u] = %u -- errorProfile[] = %u-%u  %.6f +- %.6f (%u values)\n",
              ii,
              xx,
              errorProfile[xx].bgn,
              errorProfile[xx].end,
              errorProfile[xx].dev.mean(),
              errorProfile[xx].dev.stddev(),
              errorProfile[xx].dev.size());
    }
    fclose(F);
  }
#endif
}

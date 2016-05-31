
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
#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_BestOverlapGraph.H"

static std::map<uint32,int>* containPartialOrder;

uint32* Unitig::_inUnitig     = NULL;
uint32* Unitig::_pathPosition = NULL;


#undef  SHOW_PROFILE_CONSTRUCTION
#undef  SHOW_PROFILE_CONSTRUCTION_DETAILS

void
Unitig::reverseComplement(bool doSort) {

  //  If there are contained fragments, we need to sort by position to place them correctly after
  //  their containers.  If there are no contained fragments, sorting can break the initial unitig
  //  building.  When two frags start at position zero, we'll exchange the order.  Initial unitig
  //  building depends on having the first fragment added become the last fragment in the unitig
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
  //  pathPosition map.

  if (doSort) {
    sort();
  } else {
    std::reverse(ufpath.begin(), ufpath.end());

    for (uint32 fi=0; fi<ufpath.size(); fi++)
      _pathPosition[ufpath[fi].ident] = fi;
  }
}




class epOlapDat {
public:
  epOlapDat(uint32 p, bool o, double e) {
    pos    = p;
    open   = o;
    erate  = e;
  };

  bool operator<(const epOlapDat &that)     const { return(pos < that.pos); };

  uint32  pos;
  bool    open;
  double  erate;
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






#if 1
void
Unitig::computeErrorProfileApproximate(const char *UNUSED(prefix), const char *UNUSED(label)) {
}
#endif



void
Unitig::computeErrorProfile(const char *UNUSED(prefix), const char *UNUSED(label)) {

#ifdef SHOW_PROFILE_CONSTRUCTION
  writeLog("Find error profile for tig "F_U32" of length "F_U32" with "F_SIZE_T" reads.\n",
          id(), getLength(), ufpath.size());
#endif

  errorProfile.clear();
  errorProfileIndex.clear();

  vector<epOlapDat>  olaps;



  //  Pick a set of reads to use.  We need full coverage in overlaps.




  // Scan overlaps to find those that we care about, and save their endpoints.

  for (uint32 fi=0; fi<ufpath.size(); fi++) {
    ufNode     *rdA    = &ufpath[fi];
    bool        rdAfwd = (rdA->position.bgn < rdA->position.end);
    int32       rdAlo  = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
    int32       rdAhi  = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

    uint32      ovlLen =  0;
    BAToverlap *ovl    =  OC->getOverlaps(rdA->ident, AS_MAX_ERATE, ovlLen);

    uint32      nDiffTig  = 0;
    uint32      nDiffPos  = 0;
    uint32      nIsect    = 0;

    for (uint32 oi=0; oi<ovlLen; oi++) {

      //  Reads in different tigs?  Don't care about this overlap.

      if (id() != Unitig::fragIn(ovl[oi].b_iid)) {
        nDiffTig++;
        continue;
      }

      //  Reads in same tig but not overlapping?  Don't care about this overlap.

      ufNode  *rdB    = &ufpath[ Unitig::pathPosition(ovl[oi].b_iid) ];
      bool     rdBfwd = (rdB->position.bgn < rdB->position.end);
      int32    rdBlo  = (rdBfwd) ? rdB->position.bgn : rdB->position.end;
      int32    rdBhi  = (rdBfwd) ? rdB->position.end : rdB->position.bgn;

      if ((rdAhi < rdBlo) || (rdBhi < rdAlo)) {
        nDiffPos++;
#ifdef SHOW_PROFILE_CONSTRUCTION_DETAILS
        writeLog("diffPos rdA %u=%u %u-%u rdB %u=%u %u-%u\n",
                 ovl[oi].a_iid, rdA->ident, rdAlo, rdAhi,
                 ovl[oi].b_iid, rdB->ident, rdBlo, rdBhi);
#endif
        continue;
      }

      //  Now figure out what region is covered by the overlap.

      int32    tiglo = 0;
      int32    tighi = FI->fragmentLength(rdA->ident);

      if (ovl[oi].a_hang > 0)
        tiglo += ovl[oi].a_hang;  //  Postiive hang!

      if (ovl[oi].b_hang < 0)
        tighi += ovl[oi].b_hang;  //  Negative hang!

      assert(0     <= tiglo);
      assert(0     <= tighi);
      assert(tiglo <= tighi);
      assert(tiglo <= FI->fragmentLength(rdA->ident));
      assert(tighi <= FI->fragmentLength(rdA->ident));

      //  Offset and adjust to tig coordinates

      //  Beacuse the read is placed with a lot of fudging in the positions, we need
      //  to scale the coordinates we compute here.
      double      sc  = (rdAhi - rdAlo) / (double)FI->fragmentLength(rdA->ident);

      uint32      bgn = (uint32)floor(rdAlo + sc * tiglo);
      uint32      end = (uint32)floor(rdAlo + sc * tighi);

      nIsect++;

      olaps.push_back(epOlapDat(bgn, true,  ovl[oi].erate));
      olaps.push_back(epOlapDat(end, false, ovl[oi].erate));
    }

#ifdef SHOW_PROFILE_CONSTRUCTION_DETAILS
    writeLog("tig %u read %u with %u overlaps - diffTig %u diffPos %u intersect %u\n",
             id(), rdA->ident, ovlLen, nDiffTig, nDiffPos, nIsect);
#endif
  }

#ifdef SHOW_PROFILE_CONSTRUCTION
  writeLog("tig %u generated "F_SIZE_T" olaps.\n", id(), olaps.size());
#endif

  //  Sort.

  std::sort(olaps.begin(), olaps.end());

  //  Convert coordinates into intervals.  Conceptually, squish out the duplicate numbers, then
  //  create an interval for every adjacent pair.  We need to add intervals for the first and last
  //  region.  And one more, for convenience, to hold the final 'close' values on intervals that
  //  extend to the end of the unitig.

  if (olaps[0].pos != 0)
    errorProfile.push_back(epValue(0, olaps[0].pos));

  for (uint32 bb=0, ii=1; ii<olaps.size(); ii++) {
    if (olaps[bb].pos == olaps[ii].pos)
      continue;

    errorProfile.push_back(epValue(olaps[bb].pos, olaps[ii].pos));

#ifdef SHOW_PROFILE_CONSTRUCTION_DETAILS
    writeLog("tig %u make region bb=%u ii=%i - %u %u\n", id(), bb, ii, olaps[bb].pos, olaps[ii].pos);
#endif

    bb = ii;
  }

  if (olaps[olaps.size()-1].pos != getLength())
    errorProfile.push_back(epValue(olaps[olaps.size()-1].pos, getLength()));

  errorProfile.push_back(epValue(getLength(), getLength()+1));


#ifdef SHOW_PROFILE_CONSTRUCTION
  writeLog("tig %u generated "F_SIZE_T" profile regions.\n", id(), errorProfile.size());
#endif

  //  Walk both lists, adding positive erates and removing negative erates.

  stdDev<double>  curDev;

  for (uint32 oo=0, ee=0; oo<olaps.size(); oo++) {
    if (olaps[oo].pos != errorProfile[ee].bgn)  //  Move to the next profile if the pos is different.
      ee++;                                     //  By construction, this single step should be all we need.

#ifdef SHOW_PROFILE_CONSTRUCTION_DETAILS
    writeLog("oo=%u bgn=%u -- ee=%u bgn=%u -- olaps.size "F_SIZE_T" errorProfile.size "F_SIZE_T" -- insert %d erate %f\n",
             oo, olaps[oo].pos,
             ee, errorProfile[ee].bgn,
             olaps.size(), errorProfile.size(),
             olaps[oo].open, olaps[oo].erate);
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

  //  Finalize the values.

  for (uint32 bi=0; bi<errorProfile.size(); bi++)
    errorProfile[bi].dev.finalize();

  //writeLog("tig %u generated "F_SIZE_T" profile regions with "F_U64" overlap pieces.\n",
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
    fprintf(stderr, "errorProfileIndex.size() = "F_SIZE_T"\n", errorProfileIndex.size());
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

  sprintf(N, "%s.%s.%08u.profile", prefix, label, id());

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
  sprintf(N, "%s.%s.%08u.profile.index", prefix, label, id());

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

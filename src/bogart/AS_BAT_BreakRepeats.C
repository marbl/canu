
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
 *    src/AS_BAT/AS_BAT_MergeSplitJoin.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2011-FEB-15 to 2014-MAY-03
 *      are Copyright 2011-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2015-AUG-05
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_Breaking.H"

#include "AS_BAT_RepeatJunctionEvidence.H"

#include "intervalList.H"


uint32 SPURIOUS_COVERAGE_THRESHOLD  = 6;   //  Need to have more than this coverage in non-unitig reads aligned to call it a repeat area
uint32 ISECT_NEEDED_TO_BREAK        = 15;  //  Need to have at least  this number of reads confirming a repeat junction
uint32 REGION_END_WEIGHT            = 15;  //  Pretend there are this many intersections at the end points of each repeat region

#undef  LOG_BUBBLE_TESTS
#undef  LOG_BUBBLE_FAILURE
#define LOG_BUBBLE_SUCCESS

omp_lock_t  markRepeat_breakUnitigs_Lock;





void
markRepeats_buildOverlapList(Unitig *target, double erateRepeat, set<uint32> &ovlFrags) {

  ovlFrags.clear();

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode     *frg   = &target->ufpath[fi];
    uint32      ovlLen = 0;
    BAToverlap *ovl    = OC->getOverlaps(frg->ident, erateRepeat, ovlLen);

    for (uint32 i=0; i<ovlLen; i++) {
      if (Unitig::fragIn(ovl[i].b_iid) != target->id())
        ovlFrags.insert(ovl[i].b_iid);
    }
  }
}




void
markRepeats_computeUnitigErrorRate(UnitigVector &unitigs,
                                   double        erateRepeat,
                                   Unitig       *target,
                                   double       &meanError,
                                   double       &stddevError) {

  vector<overlapPlacement>  op;
  vector<double>            error;

  meanError   = 0;
  stddevError = 0;

#undef DUMPERROR
#ifdef DUMPERROR
  char  N[FILENAME_MAX];
  sprintf(N, "error.%08d.dat", target->id());
  FILE *F = fopen(N, "w");
#endif

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode     *frg   = &target->ufpath[fi];
    uint32      bgn   = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
    uint32      end   = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

    placeFragUsingOverlaps(unitigs, erateRepeat, target, frg->ident, op);

    if (op.size() == 0)
      //  Huh?  Couldn't be placed in my own unitig?
      continue;
    assert(op.size() > 0);

    double  minError = 1.0;
    double  corError = 1.0;

    for (uint32 pl=0; pl<op.size(); pl++) {
      assert(op[pl].tigID == target->id());

      double e = op[pl].errors / op[pl].aligned;

      minError = MIN(minError, e);

      if (op[pl].position.bgn < op[pl].position.end) {
        if ((op[pl].position.end < bgn) ||
            (end < op[pl].position.bgn))
          continue;
      } else {
        if ((op[pl].position.bgn < bgn) ||
            (end < op[pl].position.end))
          continue;
      }

      corError = MIN(corError, e);
    }

#ifdef DUMPERROR
    fprintf(F, "%f\t%f\n", corError, minError);
#endif

    if (corError < 1.0)
      error.push_back(corError);
    else
      error.push_back(minError);
  }

#ifdef DUMPERROR
  fclose(F);
#endif

  for (uint32 i=0; i<error.size(); i++)
    meanError += error[i];

  meanError /= error.size();

  for (uint32 i=0; i<error.size(); i++)
    stddevError += (error[i] - meanError) * (error[i] - meanError);

  stddevError = sqrt(stddevError / error.size());

  writeLog("markRepeats_computeUnitigErrorRate()--  tig %d error %f +- %f\n",
           target->id(), meanError, stddevError);
}



void
markRepeats_placeAndProcessOverlaps(UnitigVector                     &unitigs,
                                    double                            erateRepeat,
                                    Unitig                           *target,
                                    double                           meanError,
                                    double                           stddevError,
                                    set<uint32>                      &ovlFrags,
                                    intervalList<int32>              &aligned,
                                    vector<repeatJunctionEvidence>   &evidence) {

  aligned.clear();
  evidence.clear();

  for (set<uint32>::iterator it=ovlFrags.begin(); it!=ovlFrags.end(); it++) {
    uint32  iid = *it;

    vector<overlapPlacement>  op;

    placeFragUsingOverlaps(unitigs, erateRepeat, target, iid, op);

    //  placeFragUsingOverlaps() returns the expected placement for this fragment in 'position', and
    //  the amount of the fragment covered by evidence in 'covered'.
    //
    //  Below we'll try to decipher this into two intervals, covered by evidence and not covered by
    //  evidence.

    for (uint32 pl=0; pl<op.size(); pl++) {
      assert(op[pl].tigID == target->id());

      double erate        = op[pl].errors / op[pl].aligned;
      bool   erateTooHigh = (meanError + 3 * stddevError < erate);

#if 0
      writeLog("markRepeats()-- op[%3d] tig %d frag %d fCoverage %f position %d %d verified %d %d erate %f%s\n",
               pl, op[pl].tigID, op[pl].frgID, op[pl].fCoverage,
               op[pl].position.bgn, op[pl].position.end,
               op[pl].verified.bgn, op[pl].verified.end,
               erate,
               erateTooHigh ? " - TOO HIGH" : "");
#endif

      if (erateTooHigh)
        continue;

      //  Save the aligned portion.  This will be used in conjunction with the partially aligned
      //  fragments to pick out where to split.
      //
      assert(op[pl].verified.bgn >= 0);
      assert(op[pl].verified.bgn <= target->getLength());  //  is bgn unsigned and underflowed?
      assert(op[pl].verified.end <= target->getLength());

      if (op[pl].verified.bgn < op[pl].verified.end)
        aligned.add(op[pl].verified.bgn, op[pl].verified.end - op[pl].verified.bgn);
      else
        aligned.add(op[pl].verified.end, op[pl].verified.bgn - op[pl].verified.end);

      if (op[pl].fCoverage > 0.99)
        //  No worries, fully placed.
        continue;

      if ((op[pl].position.bgn < 0) ||
          (op[pl].position.bgn > target->getLength()) ||
          (op[pl].position.end < 0) ||
          (op[pl].position.end > target->getLength()))
        //  Placed outside the range of the unitig
        continue;

      //  Otherwise, placed in the target unitig, at less than perfect coverage.  Compute
      //  the unitig coordinates that are covered by actual overlaps.

      repeatJunctionEvidence   ev(target, op[pl]);

      if (ev.tigFrag == FragmentEnd())
        //  Didn't pass muster; weak overhangs likely.
        continue;

      evidence.push_back(ev);
    }
  }
}


#if 0
uint32
markRepeats_computeUnitigCoverage(Unitig *tig) {
  intervalList<int32>   coverage;

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode  frg         = tig->ufpath[fi];

    if (frg.position.bgn < frg.position.end)
      coverage.add(frg.position.bgn, frg.position.end - frg.position.bgn);
    else
      coverage.add(frg.position.end, frg.position.bgn - frg.position.end);
  }

  coverage.merge();

  intervalDepth  depth(coverage);

  uint64   minCov = 0;           //  Probably 1 or 2
  uint64   medCov = 0;           //
  uint64   aveCov = 0;           //  If repeats present, might be high
  uint64   modCov = 0;           //  Probably not useful, due to read bias
  uint64   maxCov = UINT32_MAX;  //

  for (uint32 dd=0; dd<depth.numberOfIntervals(); dd++)
    maxCov = MAX(maxCov, depth.de(dd));
  maxCov++;

  uint32   *histogram = new uint32 [maxCov];

  memset(histogram, 0, sizeof(uint32) * maxCov);

  for (uint32 dd=0; dd<depth.numberOfIntervals(); dd++)
    histogram[depth.de(dd)] += depth.hi(dd) - depth.lo(dd);

  for (uint32 dd=0; dd<maxCov; dd++) {
    if (depth.de(dd) == 0)
      continue;

    if (minCov < depth.de(dd))
      minCov = depth.de(dd);

    if (depth.de(dd) < maxCov)
      maxCov = depth.de(dd);

    //  blah.
  }

  return(medCov);
}
#endif



//  Decide on a spurious coverage level using coverage in this unitig, global coverage, and
//  statistics from the potential repeats.
//
//  This is the end consumer of the 'aligned' data.  It is used only to populate 'regions', and
//  'regions' only cares about bgn,end coords, no underlying data.
//
void
markRepeats_filterIntervalsSpannedByFragment(Unitig                    *target,
                                             intervalList<int32>       &aligned,
                                             vector<repeatRegion>      &regions,
                                             uint32                     minOverlap) {
  uint32   tiglen  = target->getLength();

  uint32   spuriousNoiseThreshold = SPURIOUS_COVERAGE_THRESHOLD;
  uint32   filteredBases          = 0;
  uint32   filteredCovered        = 0;

  intervalList<int32>   depth(aligned);

  aligned.merge();  //  Just for a stupid log message

  writeLog("markRepeats()--  filtering low coverage spurious with t=%u in %u repeat regions (%u depth regions)\n",
           spuriousNoiseThreshold, aligned.numberOfIntervals(), depth.numberOfIntervals());

  regions.clear();
  aligned.clear();

  //
  //  Remove low coverage areas by making a new map for the high coverage areas.
  //

  for (uint32 dd=0; dd<depth.numberOfIntervals(); dd++) {
    if (depth.depth(dd) == 0)
      continue;

    if (depth.depth(dd) <= spuriousNoiseThreshold) {
      filteredBases += depth.hi(dd) - depth.lo(dd);
      continue;
    }

    aligned.add(depth.lo(dd), depth.hi(dd) - depth.lo(dd));
  }

  aligned.merge();

  writeLog("markRepeats()--  filtered %u bases, now with %u repeat regions\n",
           filteredBases, aligned.numberOfIntervals());

  //
  //  Adjust region boundaries so they land on the first read end that makes sense.
  //
  //  For the begin, decide if we need to expand or contract the region.  We will expand if the
  //  start of the read before the region is not anchored.  We will contract otherwise.
  //
  //  For the end, it is more complicated because the end points are not sorted.
  //

  for (uint32 i=0; i<aligned.numberOfIntervals(); i++) {
    uint32  intbgn  = aligned.lo(i);
    uint32  intend  = aligned.hi(i);

    for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
      ufNode     *frg    = &target->ufpath[fi];
      uint32      frgbgn = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
      uint32      frgend = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

      if (frgbgn + minOverlap / 2 < intbgn)
        //  Anchored.
        continue;

      if (frgbgn == intbgn) {
        //  Perfect!  Don't change a thing (unless we already expanded to get an unanchored read).
        break;
      }

      if ((frgbgn <= intbgn) && (intbgn <= frgbgn + minOverlap / 2)) {
        //  Not anchored, expand the region to this location
#ifdef VERBOSE_REGION_FITTING
        writeLog("markRepeats()--  region["F_U32"].bgn expands from "F_U32" to "F_U32" at frag "F_U32"\n", i, aligned.lo(i), frgbgn, frg->ident);
#endif
        aligned.lo(i) = frgbgn;
        break;
      }

      if (intbgn <= frgbgn) {
        //  First read begin that is inside the repeat region
#ifdef VERBOSE_REGION_FITTING
        writeLog("markRepeats()--  region["F_U32"].bgn contracts from "F_U32" to "F_U32" at frag "F_U32"\n", i, aligned.lo(i), frgbgn, frg->ident);
#endif
        aligned.lo(i) = frgbgn;
        break;
      }
    }

    uint32  newexp = 0, newexpid = UINT32_MAX;
    uint32  newcnt = 0, newcntid = UINT32_MAX;

    for (uint32 fi=target->ufpath.size(); fi-- > 0; ) {
      ufNode     *frg    = &target->ufpath[fi];
      uint32      frgbgn = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
      uint32      frgend = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

      if (intend + minOverlap / 2 < frgend)
        //  Anchored.
        continue;

      if (frgend == intend) {
        //  Perfect!
        //newexpid = UINT32_MAX;
        newcntid = UINT32_MAX;
        break;
      }

      if ((intend < frgend) && (frgend <= intend + minOverlap / 2) && (newexp < frgend)) {
        //  Not anchored, expand the region to this location (this will pick the largest expansion)
        newexp   = frgend;
        newexpid = frg->ident;
        continue;
      }

      if ((frgend <= intend) && (newcnt < frgend)) {
        //  Pick the largest read end that is within the repeat region
        newcnt   = frgend;
        newcntid = frg->ident;
        continue;
      }

      if (frgbgn + AS_MAX_READLEN < intend)
        //  All done, no more intersections possible
        break;
    }

    //  Expand the region if one was found, otherwise contract if one was found.

    if      (newexpid != UINT32_MAX) {
#ifdef VERBOSE_REGION_FITTING
      writeLog("markRepeats()--  region["F_U32"].end expands from "F_U32" to "F_U32" at frag "F_U32"\n", i, aligned.hi(i), newexp, newexpid);
#endif
      aligned.hi(i) = newexp;
    }

    else if (newcntid != UINT32_MAX) {
#ifdef VERBOSE_REGION_FITTING
      writeLog("markRepeats()--  region["F_U32"].end contracts from "F_U32" to "F_U32" at frag "F_U32"\n", i, aligned.hi(i), newcnt, newcntid);
#endif
      aligned.hi(i) = newcnt;
    }
  }


  {
    uint32 nc = 0;

    for (uint32 i=0; i<aligned.numberOfIntervals(); i++)
      if (aligned.hi(i) < aligned.lo(i))
        nc++;

    writeLog("markRepeats()--  filtered "F_U32" repeat regions after picking read endpoints, now with "F_U32" repeat regions.\n",
             nc, aligned.numberOfIntervals() - nc);
  }

  //
  //  Discard the interval if there is a fragment that contains it with enough overhang to
  //  unambiguously place it in a unitig.
  //

  for (uint32 i=0; i<aligned.numberOfIntervals(); i++) {
    ufNode  *containedIn  = NULL;

    //  If the region is backwards, then the region is contained in a read.  The easiest case to
    //  argue is:
    //
    //           ------------
    //                  -------RRR------
    //                              --------------
    //
    //  The read with the repeat region is anchored on both sides, so it is contracted to the next
    //  begin (for the start) and the previous end for the end)
    //
    if (aligned.hi(i) < aligned.lo(i)) {
      //writeLog("markRepeats()--  repeat alignment "F_U64","F_U64" DISCARD - CONTRACTED TO NULL\n",
      //         aligned.lo(i), aligned.hi(i));
      continue;
    }

    //  Ensure that reads not near the end of the unitig have enough non-repeat sequence to anchor the read in this location.
    //  This is done by increasing the size of the repeat.

    uint32   rptbgn  = aligned.lo(i);
    uint32   rptend  = aligned.hi(i);

    uint32   unique  = minOverlap / 2;

    bool     bgnFull = true;
    bool     endFull = true;

    if (unique <= rptbgn) {
      bgnFull = false;
      rptbgn -= unique;
    }

    if (rptend + unique <= tiglen) {
      endFull = false;
      rptend += unique;
    }

    //  Search for a covering fragment.

    for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
      ufNode     *frg    = &target->ufpath[fi];
      uint32      frgbgn = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
      uint32      frgend = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

      if (frgend < rptbgn)
        //  Fragment is before the region, keep searching for a spanning fragment.
        continue;

      if (rptend < frgbgn)
        //  Fragment is after the region, we're finished.
        break;

      if ((frgbgn <= rptbgn) &&
          (rptend <= frgend)) {
        //  Fragment contains the region with acceptable overhangs into the non-repeat area.
        containedIn = frg;
        break;
      }
    }

    if (containedIn != NULL) {
      filteredCovered++;
      //writeLog("markRepeats()--  repeat alignment %s"F_U64","F_U64"%s DISCARD - CONTAINED IN FRAGMENT "F_U32" "F_S32","F_S32"\n",
      //         (bgnFull ? "(end) " : ""), aligned.lo(i), aligned.hi(i), (endFull ? " (end)" : ""),
      //         containedIn->ident, containedIn->position.bgn, containedIn->position.end);
      continue;
    }

    if ((regions.size() == 0) ||
        (regions.back().end + 100 < aligned.lo(i))) {
      regions.push_back(repeatRegion(aligned.lo(i), aligned.hi(i)));
      //writeLog("markRepeats()--  repeat alignment "F_U32","F_U32"\n",
      //        regions.back().bgn, regions.back().end);

    } else {
      regions.back().bgn       = MIN(regions.back().bgn, aligned.lo(i));
      regions.back().end       = MAX(regions.back().end, aligned.hi(i));
      //writeLog("markRepeats()--  repeat alignment "F_U32","F_U32" (merged from "F_U64","F_U64")\n",
      //        regions.back().bgn, regions.back().end,
      //        aligned.lo(i), aligned.hi(i));
    }
  }

  writeLog("markRepeats()--  filtered %u repeat regions contained in a read, now with %u repeat regions\n",
           filteredCovered, regions.size());
}




//  Unitig fragment is completely within the repeat interval, or is close enough to the edge that
//  maybe it couldn't be placed uniquely.
//
void
markRepeats_findFragsInRegions(Unitig                    *target,
                               vector<repeatRegion>      &regions,
                               set<uint32>               &rptFrags,
                               set<uint32>               &UNUSED(ejtFrags),
                               uint32                     minOverlap) {

  for (uint32 i=0; i<regions.size(); i++) {
    for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
      ufNode     *frg   = &target->ufpath[fi];
      uint32      bgn   = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
      uint32      end   = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

      if (regions[i].bgn == bgn)
        regions[i].rujBgn = repeatUniqueBreakPoint(regions[i].bgn,
                                                   FragmentEnd(frg->ident, end < bgn),
                                                   false);

      if (regions[i].end == end)
        regions[i].rujEnd = repeatUniqueBreakPoint(regions[i].end,
                                                   FragmentEnd(frg->ident, bgn < end),
                                                   true);

      //  If the read has at least minOverlap/2 bases outside the repeat, assume it
      //  is placed correctly.
      if ((bgn + minOverlap/2 < regions[i].bgn) ||
          (regions[i].end + minOverlap/2 < end))
        continue;

      //  Otherwise, the read is 'contained' in a repeat region.  Remember it for later processing.
      rptFrags.insert(frg->ident);

      //  Read is unanchored in a repeat region, toss it out, but place it with the mate.
      //ejtFrags.insert(frg->ident);
    }
  }
}






//  Find any junctions in a region, and append them to list of fragments to split on.
//  This takes multiple lines of evidence (pointing to the same fragment end) and
//  combines them into a list of break points.
void
markRepeats_filterJunctions(Unitig                          *target,
                            vector<repeatRegion>            &regions,
                            vector<repeatJunctionEvidence>  &evidence,
                            vector<repeatUniqueBreakPoint>  &breakpoints) {

  map<FragmentEnd,uint32>                   brkFrags;
  map<FragmentEnd,repeatUniqueBreakPoint>   brkJunct;

  //  Add breakpoints for each of the end points of the repeat region.

  for (uint32 rr=0; rr<regions.size(); rr++) {
    FragmentEnd &rujBgnFrg = regions[rr].rujBgn.breakFrag;
    FragmentEnd &rujEndFrg = regions[rr].rujEnd.breakFrag;

    if (rujBgnFrg.fragId() > 0) {
      brkFrags[rujBgnFrg] = REGION_END_WEIGHT;
      brkJunct[rujBgnFrg] = regions[rr].rujBgn;
    }

    if (rujEndFrg.fragId() > 0) {
      brkFrags[rujEndFrg] = REGION_END_WEIGHT;
      brkJunct[rujEndFrg] = regions[rr].rujEnd;
    }
  }

  //

  sort(evidence.begin(), evidence.end());

  for (uint32 ai=0, bi=0; ai<evidence.size(); ai++) {
    repeatUniqueBreakPoint  ruj;

    if (evidence[ai].is3 == false) {
      assert(evidence[ai].uncovered5bgn < evidence[ai].uncovered5end);
      ruj = repeatUniqueBreakPoint(evidence[ai].point, evidence[ai].tigFrag, false);
    } else {
      assert(evidence[ai].uncovered3bgn < evidence[ai].uncovered3end);
      ruj = repeatUniqueBreakPoint(evidence[ai].point, evidence[ai].tigFrag, true);
    }

    //  Try to associate this junction with one of the repeat regions.  If there is no region,
    //  this is NOT a junction we care about.  The region must have been contained in a fragment.

    while ((bi < regions.size()) &&
           (regions[bi].end < ruj.point))
      //  Advance the region until it ends after the point.
      bi++;

    //  If this point is in the region, the region bgn will be lower (or equal) than the point.  We
    //  already ensured that the region end is after the point.

    if ((bi >= regions.size()) ||
        (ruj.point < regions[bi].bgn))
      continue;

    assert(regions[bi].bgn <= ruj.point);
    assert(ruj.point       <= regions[bi].end);

    //  A new valid break point.

    //  NOTE!  ruj's seem to be different.  We used to save the 5th ruj, and switching to saving the
    //  last showed differences in position (ruj.point) of the break.  The point comes directly from
    //  the evidence[] above, so no surprise.

    brkFrags[evidence[ai].tigFrag]++;
    brkJunct[evidence[ai].tigFrag] = ruj;
  }

  for (map<FragmentEnd,uint32>::iterator it=brkFrags.begin(); it != brkFrags.end(); it++) {
    uint32                  cnt = brkFrags[it->first];
    repeatUniqueBreakPoint  ruj = brkJunct[it->first];

    if (cnt < ISECT_NEEDED_TO_BREAK)
      continue;

    breakpoints.push_back(ruj);
  }

  sort(breakpoints.begin(), breakpoints.end());

  writeLog("markRepeats()--  unitig %d has "F_SIZE_T" interesting junctions:\n",
           target->id(), breakpoints.size());

  for (uint32 ji=0; ji<breakpoints.size(); ji++)
    writeLog("markRepeats()--  junction["F_U32"] at "F_U32"/%c' position "F_U32" repeat %s count "F_U32"\n",
             ji,
             breakpoints[ji].breakFrag.fragId(), breakpoints[ji].breakFrag.frag3p() ? '3' : '5',
             breakpoints[ji].point,
             breakpoints[ji].rptLeft ? "<-" : "->",
             brkFrags[breakpoints[ji].breakFrag]);
}




void
markRepeats_breakUnitigs(UnitigVector                    &unitigs,
                         double                           erateRepeat,
                         Unitig                          *target,
                         vector<overlapPlacement>        &UNUSED(places),
                         vector<repeatUniqueBreakPoint>  &breakpoints,
                         set<uint32>                     &jctFrags,
                         set<uint32>                     &rptFrags,
                         set<uint32>                     &ejtFrags) {

  jctFrags.clear();

  if (breakpoints.size() == 0)
    return;

  uint32   *breakID  = new uint32 [target->ufpath.size()];

  uint32    nextBreakPoint = 0;
  uint32    curr = 1;
  uint32    next = 2;

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode     *frg   = &target->ufpath[fi];
    uint32      bgn   = (frg->position.bgn < frg->position.end) ? frg->position.bgn : frg->position.end;
    uint32      end   = (frg->position.bgn < frg->position.end) ? frg->position.end : frg->position.bgn;

    //  If out of breakpoints, put all the remaining reads into the current tig.
    if (nextBreakPoint >= breakpoints.size()) {
      breakID[fi]  = curr;

    } else if (breakpoints[nextBreakPoint].rptLeft == false) {
      //  Repeat to the right.  If the fragment starts at or after the junction, place this and
      //  future fragments into a new (repeat) unitig.

      if (breakpoints[nextBreakPoint].point <= bgn) {
        nextBreakPoint++;
        curr++;
        next++;
        jctFrags.insert(frg->ident);
      }

      //  If the fragment ends after the next junction, this fragment goes to the next
      //  unitig.  Otherwise, to the current one.
      //
      if ((nextBreakPoint < breakpoints.size()) &&
          (breakpoints[nextBreakPoint].point < end) &&
          (breakpoints[nextBreakPoint].rptLeft == true)) {
        breakID[fi] = next;
      } else {
        breakID[fi] = curr;
      }

    } else {
      //  Repeat to the left.  If the fragment ends before the junction, move to the current unitig

      if (end < breakpoints[nextBreakPoint].point) {
        breakID[fi] = curr;
        jctFrags.insert(frg->ident);
      } else {
        breakID[fi] = next;
      }

      //  Once we pass the junction, update pointers.  We are out of the repeat interval now.
      if (breakpoints[nextBreakPoint].point < bgn) {
        nextBreakPoint++;
        curr++;
        next++;
      }
    }
  }

  //  Append new unitigs.

  vector<Unitig *>      newTigs;
  Unitig              **uidToUnitig = new Unitig * [next + 1];
  uint32               *uidToOffset = new uint32   [next + 1];

  memset(uidToUnitig, 0, sizeof(Unitig *) * (next + 1));
  memset(uidToOffset, 0, sizeof(uint32)   * (next + 1));

  for (uint32 fi=0; fi<target->ufpath.size(); fi++) {
    ufNode  &frg = target->ufpath[fi];
    uint32   bid = breakID[fi];

    if (ejtFrags.count(frg.ident) > 0) {
      writeLog("markRepeats()-- EJECT unanchored frag %u from unitig %u\n",
               frg.ident, target->id());
      target->removeFrag(frg.ident);
      continue;
    }

    if (uidToUnitig[bid] == NULL) {
      uidToUnitig[bid] = unitigs.newUnitig(false);  //  Add a new unitig to the unitigs list
      //uidToUnitig[bid]->_isRepeat = true;

      uidToOffset[bid] = -MIN(frg.position.bgn, frg.position.end);

      newTigs.push_back(uidToUnitig[bid]);  //  For reporting below.
    }

    uidToUnitig[bid]->addFrag(frg, uidToOffset[bid], false);
  }

  delete [] breakID;
  delete [] uidToUnitig;
  delete [] uidToOffset;


  if (newTigs.size() > 0) {
    writeLog("markRepeats()-- SPLIT unitig %d of length %u with %ld fragments into "F_SIZE_T" unitigs:\n",
             target->id(), target->getLength(), target->ufpath.size(),
             newTigs.size());

    for (uint32 ti=0; ti<newTigs.size(); ti++) {
      Unitig *tig   = newTigs[ti];
      uint32  nRept = 0;
      uint32  nUniq = 0;

      for (uint32 fi=0; fi < tig->ufpath.size(); fi++) {
        ufNode  &frg = tig->ufpath[fi];

        if (rptFrags.count(frg.ident) == 1)
          nRept++;
        else
          nUniq++;
      }

      if (nRept > nUniq)
        tig->_isRepeat = true;

      writeLog("markRepeats()--   unitig %u of length %u with %ld fragments (%u %.4f repeat and %u %.4f non-repeat).\n",
               tig->id(),
               tig->getLength(),
               tig->ufpath.size(),
               nRept, (double)nRept / (nRept + nUniq),
               nUniq, (double)nUniq / (nRept + nUniq));
    }

    writeLog("markRepeats()-- DELETE unitig %d\n", target->id());
    unitigs[target->id()] = NULL;
    delete target;
  }

  //  Run back over the ejected frags, and place them at their best location.

  for (set<uint32>::iterator it=ejtFrags.begin(); it!=ejtFrags.end(); it++) {
    writeLog("markRepeats()-- EJECT frag "F_U32"\n", *it);
    placeFragInBestLocation(unitigs, erateRepeat, *it);
  }

  writeLog("markRepeats()-- FINISHED.\n");
}



void
markRepeats_shatterRepeats(UnitigVector   &unitigs,
                           set<uint32>    &jctFrags,
                           set<uint32>    &rptFrags) {

  //  For each junction read (defined to be the first/last read in a repeat unitig), shatter the unitig.

  for (set<uint32>::iterator it=jctFrags.begin(); it!=jctFrags.end(); it++) {
    uint32   iid = *it;
    uint32   ti  = Unitig::fragIn(iid);
    Unitig  *rpt = unitigs[ti];

    if ((ti == 0) || (rpt == NULL))
      //  Already shattered?
      continue;

    writeLog("markRepeats()--  shatter unitig %u with "F_SIZE_T" fragments from repeat frag %u\n",
             ti, rpt->ufpath.size(), iid);

    for (uint32 fi=0; fi<rpt->ufpath.size(); fi++)
      rpt->removeFrag(rpt->ufpath[fi].ident);

    unitigs[ti] = NULL;
    delete rpt;
  }

  //  For each repeat read (defined to be a read contained nearly entirely in a repeat region), count
  //  the number that are still in a unitig.

  for (set<uint32>::iterator it=rptFrags.begin(); it!=rptFrags.end(); it++) {
    uint32   iid = *it;
    uint32   ti  = Unitig::fragIn(iid);
    Unitig  *rpt = unitigs[ti];

    if ((ti == 0) || (rpt == NULL))
      //  Already shattered?
      continue;

    writeLog("markRepeats()--  frag "F_U32" covered by repeats, but still in unitig %d\n",
             iid, ti);
  }
}




//  Build a list of all the fragments that have overlaps to some fragment in this unitig.
//  Exclude from the list fragments that are already in this unitig.  We expect these fragments
//  to have multiple overlaps to the unitig, and we want to examine each one only once.
//
//  Use placeFragUsingOverlaps() to place each of these fragments.  Ideally, we'd restrict to just
//  this unitig, but for now we need to filter the results.  Three results:
//
//    o Fragment is fully contained in this unitig.  Why isn't it a bubble?  It could be contained
//      completely in a repeat, and the repeat is in two different unitigs.
//
//    o Fragment is partially aligned.  This could be indicating a repeat or chimera that we should
//      be splitting.  We save the location of the break, and direction of the unaligned piece.
//
//  After all fragments are 'placed' the list of breaks is examined.
//
//    o A chimer will induce about as many breaks as the local depth.  Breaks will be on both
//      sides of the chimeric point.
//
//    o A repeat will have many breaks on one side, and few to none on the other side.
//
//    o A spur will look like a repeat, but at the end of a unitig, with few to no fragments
//      following.
//
void
markRepeats(UnitigVector &unitigs,
            double erateRepeat,
            Unitig *target,
            uint32 minOverlap,
            bool shatterRepeats) {

  set<uint32>                     ovlFrags;
  set<uint32>                     rptFrags;  //  Frag IIDs of fragments covered by repeat alignments
  set<uint32>                     jctFrags;  //  Frag IIDs of the first/last fragment in a repeat unitig
  set<uint32>                     ejtFrags;  //  Frag IIDs of frags we should eject instead of split

  double                          meanError = 0;
  double                          stddevError = 0;

  intervalList<int32>             aligned;

  vector<overlapPlacement>        places;

  vector<repeatRegion>            regions;
  vector<repeatJunctionEvidence>  evidence;

  vector<repeatUniqueBreakPoint>  breakpoints;

  //  Build a list of all the fragments that have overlaps to this unitig.
  markRepeats_buildOverlapList(target, erateRepeat, ovlFrags);

  //  Decide what a decent alignment should look like.
  markRepeats_computeUnitigErrorRate(unitigs, erateRepeat, target, meanError, stddevError);

  //  For each overlapping fragment, place it and process.
  markRepeats_placeAndProcessOverlaps(unitigs, erateRepeat, target, meanError, stddevError, ovlFrags, aligned, evidence);

  //  Convert 'aligned' into regions, throwing out weak ones and those contained in a fragment.
  markRepeats_filterIntervalsSpannedByFragment(target, aligned, regions, minOverlap);

  markRepeats_findFragsInRegions(target, regions, rptFrags, ejtFrags, minOverlap);

  //  Discard junctions that are not in a remaining region.
  markRepeats_filterJunctions(target, regions, evidence, breakpoints);

  //  Split at whatever junctions remain.

  //  You'd think declaring this a critical region would work, but it resulted in deadlock on
  //    Linux 2.6.32-279.22.1.el6.x86_64
  //    g++ (GCC) 4.7.1
  //#pragma omp critical

  omp_set_lock(&markRepeat_breakUnitigs_Lock);
  markRepeats_breakUnitigs(unitigs, erateRepeat, target, places, breakpoints, jctFrags, rptFrags, ejtFrags);
  omp_unset_lock(&markRepeat_breakUnitigs_Lock);

  //  For each repeat unitig, shatter into fragments (not singleton unitigs) so we can later re-BOG.
  if (shatterRepeats)
    markRepeats_shatterRepeats(unitigs, jctFrags, rptFrags);
}





void
breakRepeats(UnitigVector &unitigs,
             double UNUSED(erateGraph), double erateBubble, double UNUSED(erateMerge), double erateRepeat,
             const char *prefix,
             uint32 minOverlap,
             bool shatterRepeats,
             uint64 genomeSize) {

  //logFileFlags |= LOG_PLACE_FRAG;
  //logFileFlags &= ~LOG_PLACE_FRAG;

  //  Since we create new unitigs for any of the splits, we need to remember
  //  where to stop.  We don't want to re-examine any of the split unitigs.
  //  Reevaluating seems to just trim off a few fragments at the end of the unitig.

  uint32  tiLimit = unitigs.size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  writeLog("repeatDetect()-- working on "F_U32" unitigs, with "F_U32" threads.\n", tiLimit, numThreads);

  omp_init_lock(&markRepeat_breakUnitigs_Lock);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig        *target = unitigs[ti];

    if ((target == NULL) ||
        (target->ufpath.size() < 15) ||
        (target->getLength() < 300))
      continue;

    //writeLog("repeatDetect()-- WORKING on unitig %d/"F_SIZE_T" of length %u with %ld fragments.\n",
    //         target->id(), unitigs.size(), target->getLength(), target->ufpath.size());

    markRepeats(unitigs, erateRepeat, target, minOverlap, shatterRepeats);
  }

  omp_destroy_lock(&markRepeat_breakUnitigs_Lock);

  logFileFlags &= ~LOG_PLACE_FRAG;
}

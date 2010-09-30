
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

static const char *rcsid = "$Id: AS_BOG_Breaking.cc,v 1.5 2010-09-30 11:32:48 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"




UnitigBreakPoint UnitigGraph::selectSmall(ContainerMap &cMap,
                                          const Unitig *tig,
                                          const UnitigBreakPoints &smalls,
                                          const UnitigBreakPoint& big,
                                          int   &lastBPCoord,
                                          int   &lastBPFragNum) {
  UnitigBreakPoint selection;

  double difference = 0.0;
  int  rFrgs   = big.fragsBefore;
  bool bRev    = isReverse(big.fragPos);
  int bContain = 0;

  int right    = (big.fragEnd.fragEnd() == FIVE_PRIME) ? big.fragPos.bgn : big.fragPos.end;

  if (cMap.find(big.fragEnd.fragId()) != cMap.end())
    bContain = cMap[big.fragEnd.fragId()].size();

  if (((bRev == true)  && (big.fragEnd.fragEnd() == THREE_PRIME)) ||
      ((bRev == false) && (big.fragEnd.fragEnd() == FIVE_PRIME)))
    rFrgs -= 1 + bContain;

  for(UnitigBreakPoints::const_iterator sIter = smalls.begin(); sIter != smalls.end(); sIter++) {
    UnitigBreakPoint small = *sIter;

    if (small.fragEnd.fragId() == big.fragEnd.fragId())
      continue;

    bool rev     = isReverse(small.fragPos);
    int sContain = 0;
    int lFrgs    = small.fragsBefore;
    uint32 sid     = small.fragEnd.fragId();

    if (cMap.find(sid) != cMap.end())
      sContain = cMap[sid].size();

    // left side of the frag in the unitig, don't count it
    if (((rev == true)  && (small.fragEnd.fragEnd() == THREE_PRIME)) ||
        ((rev == false) && (small.fragEnd.fragEnd() == FIVE_PRIME)))
      lFrgs -= 1 + sContain;

    if (rFrgs - lFrgs == 1)
      continue;

    // use middle of frag instead of end, to get some overlap
    double bp    = (small.fragPos.end + small.fragPos.bgn) / 2.0;

    double lRate = (lFrgs - lastBPFragNum) / (bp - lastBPCoord);
    double rRate = (rFrgs - lFrgs) / (right - bp);
    double ratio = (lRate > rRate) ? lRate / rRate : rRate / lRate;
    double diff  = fabs( lRate - rRate);

    if ((ratio > 1.8) && (diff > difference)) {
      //fprintf(logFile, "Break frg %7d b %4d l %4d pos b %5d e %5.0f lRate %.4f\n", sid, lastBPFragNum, lFrgs, lastBPCoord, bp, lRate );
      //fprintf(logFile, "     diff %4d r %4d pos %5d rRate %.4f ratio %.2f to frag %7d\n", rFrgs - lFrgs, rFrgs, right, rRate, ratio, big.fragEnd.fragId());
      //fprintf(logFile,"     select frg %d for break on arrival rate diff %.4f\n", sid, diff);
      difference = diff;
      selection = small;
    }
  }

  lastBPCoord   = right;
  lastBPFragNum = rFrgs;

  return selection;
}



void UnitigGraph::filterBreakPoints(ContainerMap &cMap,
                                    Unitig *tig,
                                    UnitigBreakPoints &breaks) {

  int lastBPCoord = 0;
  int lastBPFragNum = 0;

  // Check for sentinal at end, not added by MateChecker
  UnitigBreakPoint fakeEnd = breaks.back();
  if (fakeEnd.inSize == std::numeric_limits<int>::max())
    breaks.pop_back();

  UnitigBreakPoints smallBPs;
  UnitigBreakPoints newBPs;

  for(UnitigBreakPoints::iterator iter = breaks.begin(); iter != breaks.end(); iter++) {
    UnitigBreakPoint nextBP = *iter;

    if ((nextBP.inFrags > UnitigGraph::MIN_BREAK_FRAGS) && (nextBP.inSize > UnitigGraph::MIN_BREAK_LENGTH)) {

      // big one, compare against smalls -- if we haven't
      // seen this break point already.
      //
      if (newBPs.empty() || (nextBP.fragEnd != newBPs.back().fragEnd)) {

        if (smallBPs.empty() == false) {
          //  smalls exist, select one.
          UnitigBreakPoint small = selectSmall(cMap, tig, smallBPs, nextBP, lastBPCoord, lastBPFragNum);
          if (small.fragEnd.fragId() > 0)
            newBPs.push_back(small);
          smallBPs.clear();

        } else {
          //  No smalls.  Update state to move past the current big breakpoint.
          lastBPCoord   = (nextBP.fragEnd.fragEnd() == FIVE_PRIME) ? nextBP.fragPos.bgn : nextBP.fragPos.end;
          lastBPFragNum = nextBP.fragsBefore;

          int bContain = 0;

          if (cMap.find(nextBP.fragEnd.fragId()) != cMap.end())
            bContain = cMap[nextBP.fragEnd.fragId()].size();

          bool bRev = isReverse(nextBP.fragPos);

          if (((bRev == true)  && (nextBP.fragEnd.fragEnd() == THREE_PRIME)) ||
              ((bRev == false) && (nextBP.fragEnd.fragEnd() == FIVE_PRIME)))
            lastBPFragNum -= 1 + bContain;
        }

        //  Haven't seen this break, so add it.
        newBPs.push_back( nextBP );
      }
    } else {
      //  Not a big breaker.  Save the small breaker if we've not seen it yet.
#warning THIS IS VALGRIND UNCLEAN
      if (newBPs.empty() ||
          smallBPs.empty() ||
          ((nextBP.fragEnd != newBPs.back().fragEnd) &&
           (nextBP.fragEnd != smallBPs.back().fragEnd)))
        smallBPs.push_back(nextBP);
    }
  }

  //  If we've got small ones saved, select one.

  if (smallBPs.empty() == false) {
    UnitigBreakPoint small = selectSmall(cMap, tig, smallBPs, fakeEnd, lastBPCoord, lastBPFragNum);
    if (small.fragEnd.fragId() > 0)
      newBPs.push_back(small);
    smallBPs.clear();
  }

  breaks = newBPs;
}

// Doesn't handle contained frags yet
UnitigVector* UnitigGraph::breakUnitigAt(Unitig *tig,
                                         UnitigBreakPoints &breaks) {

  if (breaks.empty())
    return NULL;

  //  This is to catch a bug in....something before us.
  //  Sometimes the breakpoint list contains nonsense that tells
  //  us to split a unitig before the first fragment, or after
  //  the last.
  //
  //  That code is pretty dense (and big) so instead we'll make
  //  one pass through the list of breaks and see if we'd do any
  //  splitting.
  //
  //  It's a huge code duplication of the main loop.
  //
  {
    UnitigBreakPoints breakstmp = breaks;

    UnitigBreakPoint  breakPoint = breakstmp.front();
    breakstmp.pop_front();

    UnitigBreakPoint nextBP;

    int   newUnitigsConstructed = 0;
    int   newUnitigExists = 0;

    if (!breakstmp.empty())
      nextBP = breakstmp.front();

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode  frg = tig->ufpath[fi];

      bool bothEnds = false;
      while ( !breakstmp.empty() && nextBP.fragEnd.fragId() == breakPoint.fragEnd.fragId() ) {
        if (nextBP.fragEnd.fragEnd() != breakPoint.fragEnd.fragEnd())
          bothEnds = true;
        breakstmp.pop_front();
        if (!breakstmp.empty())
          nextBP = breakstmp.front();
      }

      bool reverse = isReverse(frg.position);

      if (breakPoint.fragEnd.fragId() == frg.ident) {

        if (bothEnds) {
          newUnitigsConstructed++;
          newUnitigExists = 0;
        }

        else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && !reverse ||
                 breakPoint.fragEnd.fragEnd() == THREE_PRIME && reverse) {
          newUnitigsConstructed++;
          newUnitigExists = 1;
        }

        else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && reverse ||
                 breakPoint.fragEnd.fragEnd() == THREE_PRIME && !reverse ) {
          if (newUnitigExists == 0)
            newUnitigsConstructed++;
          newUnitigExists = 0;
        } else {
          assert(0);
        }

        if (breakstmp.empty()) {
          breakPoint.fragEnd = FragmentEnd();
        } else {
          breakPoint = breakstmp.front();
          breakstmp.pop_front();
          if (!breakstmp.empty())
            nextBP = breakstmp.front();
        }
      } else {
        if (newUnitigExists == 0) {
          newUnitigsConstructed++;
          newUnitigExists = 1;
        }
      }
    }

    if (newUnitigsConstructed < 2) {
      //fprintf(logFile, "SPLITTING BUG DETECTED!  Adjusting for it.\n");
      return NULL;
    }

  }  //  end of explicit split test




  UnitigVector *splits = new UnitigVector();
  Unitig       *newTig = NULL;
  int           offset = 0;

  //  Emit a log of unitig creation and fragment moves
  //
  UnitigBreakPoint breakPoint = breaks.front();
  breaks.pop_front();

  UnitigBreakPoint nextBP;

  if (!breaks.empty())
    nextBP = breaks.front();

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode  frg = tig->ufpath[fi];

    // reduce multiple breaks at the same fragment end down to one
    bool bothEnds = false;
    while ( !breaks.empty() && nextBP.fragEnd.fragId() == breakPoint.fragEnd.fragId() ) {
      if (nextBP.fragEnd.fragEnd() != breakPoint.fragEnd.fragEnd())
        bothEnds = true;
      breaks.pop_front();
      if (!breaks.empty())
        nextBP = breaks.front();
    }

    bool reverse = isReverse(frg.position);

    if (breakPoint.fragEnd.fragId() == frg.ident) {

      //  There is method to this madness.  The three if
      //  blocks below have the same two basic code blocks
      //  in various orders.
      //
      //  o One code block makes a new unitig.
      //
      //  o The other adds a fragment to it, possibly
      //    remembering the beginning offset.
      //
      //  Finally, we delay making a new unitig until we are
      //  absolutely sure we're going to put fragments into
      //  it.

      //  Notes on fixing the break/join bug -- the bug is when we have a big breaker
      //  come in and chop off a small guy at the end.  We should be joining the two
      //  large unitigs.
      //
      //  To fix this, BPW was going to store enough stuff in the breakPoint so that
      //  the case can be detected here.  We need to know how much is before/after
      //  (depends on which end of the unitig we're at) our incoming point, (e.g.,
      //  fragsBefore and fragsAfter), and the fragments involved (fragEnd and inEnd).
      //
      //  When detected, push fragEnd and inEnd onto a list of joins.  After all breaks
      //  are finished, use the list of joins to get the unitigs to join and, well,
      //  join them.

      if (bothEnds) {
        //
        //  Break at both ends, create a singleton for this fragment.
        //
        if (logFileFlagSet(LOG_INTERSECTION_BREAKING) && newTig)
          fprintf(logFile, "Done with newTig unitig %d has %d fragments.\n", newTig->id(), newTig->ufpath.size());

        if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
          fprintf(logFile, "Break tig %d at both ends of frag %d\n", tig->id(), breakPoint.fragEnd.fragId());

        newTig = new Unitig(logFileFlagSet(LOG_INTERSECTION_BREAKING));  //  always make a new unitig, we put a frag in it now
        splits->push_back(newTig);

        if (newTig->ufpath.empty())
          offset = reverse ? -frg.position.end : -frg.position.bgn;
        newTig->addFrag(frg, offset, logFileFlagSet(LOG_INTERSECTION_BREAKING));

        if (logFileFlagSet(LOG_INTERSECTION_BREAKING) && newTig)
          fprintf(logFile, "Done with newTig unitig %d has %d fragments.\n", newTig->id(), newTig->ufpath.size());

        newTig = NULL;  //  delay until we need to make it
      }

      else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && !reverse ||
               breakPoint.fragEnd.fragEnd() == THREE_PRIME && reverse) {
        //
        //  Break at left end of frg, frg starts new tig
        //
        if (logFileFlagSet(LOG_INTERSECTION_BREAKING) && newTig)
          fprintf(logFile, "Done with newTig unitig %d has %d fragments.\n", newTig->id(), newTig->ufpath.size());

        if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
          fprintf(logFile,"Break tig %d before frag %d\n", tig->id(), breakPoint.fragEnd.fragId());

        newTig = new Unitig(logFileFlagSet(LOG_INTERSECTION_BREAKING));  //  always make a new unitig, we put a frag in it now
        splits->push_back(newTig);

        if (newTig->ufpath.empty())
          offset = reverse ? -frg.position.end : -frg.position.bgn;
        newTig->addFrag(frg, offset, logFileFlagSet(LOG_INTERSECTION_BREAKING));
      }

      else if (breakPoint.fragEnd.fragEnd() ==  FIVE_PRIME && reverse ||
               breakPoint.fragEnd.fragEnd() == THREE_PRIME && !reverse) {
        //
        //  Break at right end of frg, frg goes in existing tig, then make new unitig
        //

        if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
          fprintf(logFile,"Break tig %d up to frag %d\n", tig->id(), breakPoint.fragEnd.fragId());

        if (newTig == NULL) {  //  delayed creation?
          newTig = new Unitig(logFileFlagSet(LOG_INTERSECTION_BREAKING));
          splits->push_back(newTig);
        }

        if (newTig->ufpath.empty())
          offset = reverse ? -frg.position.end : -frg.position.bgn;

        newTig->addFrag(frg, offset, logFileFlagSet(LOG_INTERSECTION_BREAKING));

        if (logFileFlagSet(LOG_INTERSECTION_BREAKING) && newTig)
          fprintf(logFile, "Done with newTig unitig %d has %d fragments.\n", newTig->id(), newTig->ufpath.size());

        newTig = NULL;

        if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
          fprintf(logFile,"Break tig %d after %d\n", tig->id(), breakPoint.fragEnd.fragId());
      } else {
        // logically impossible!
        assert(0);
      }


      // Done breaking, continue adding remaining frags to newTig
      if (breaks.empty()) {
        breakPoint.fragEnd = FragmentEnd();
      } else {
        breakPoint = breaks.front();
        breaks.pop_front();
        if (!breaks.empty())
          nextBP = breaks.front();
      }
    } else {
      //
      //  frag is not a unitig break point fragment, add to current new tig
      //

      if (newTig == NULL) {  //  delayed creation?
        if (logFileFlagSet(LOG_INTERSECTION_BREAKING))
          fprintf(logFile,"Break tig %d up until some later fragment\n", tig->id());
        newTig = new Unitig(logFileFlagSet(LOG_INTERSECTION_BREAKING));
        splits->push_back(newTig);
      }
      if (newTig->ufpath.empty())
        offset = reverse ? -frg.position.end : -frg.position.bgn;
      newTig->addFrag(frg, offset, logFileFlagSet(LOG_INTERSECTION_BREAKING));
    }
  }

  return splits;
}


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

static const char *rcsid = "$Id: AS_BOG_MateChecker.cc,v 1.105 2011-12-29 09:26:03 brianwalenz Exp $";

#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_MateLocation.hh"



//  True if interval a contains interval b.
//
bool contains(SeqInterval a, SeqInterval b) {
  int aMin,aMax,bMin,bMax;
  if (isReverse(a)) { aMin = a.end; aMax = a.bgn; }
  else              { aMin = a.bgn; aMax = a.end; }
  if (isReverse(b)) { bMin = b.end; bMax = b.bgn; }
  else              { bMin = b.bgn; bMax = b.end; }
  if (aMin <= bMin && aMax >= bMax)
    return true;
  else
    return false;
}

//  Returns the intersection of intervals a and b.
//
SeqInterval intersection(SeqInterval a, SeqInterval b) {
  SeqInterval retVal = NULL_SEQ_LOC;
  int aMin,aMax,bMin,bMax;
  if (isReverse(a)) { aMin = a.end; aMax = a.bgn; }
  else              { aMin = a.bgn; aMax = a.end; }
  if (isReverse(b)) { bMin = b.end; bMax = b.bgn; }
  else              { bMin = b.bgn; bMax = b.end; }

  if (aMax < bMin || bMax < aMin)
    return retVal;

  // so now aMax > bMin && bMax > aMin, thus intersection
  retVal.bgn = MAX(aMin, bMin);
  retVal.end = MIN(aMax, bMax);
  return retVal;
}


static
vector<SeqInterval> *
findPeakBad(int32 *badGraph,
            int32 tigLen,
            int32 badMateBreakThreshold) {
  vector<SeqInterval> *peakBads = new vector<SeqInterval>;
  SeqInterval          peak     = {0, 0};
  int32                peakBad  = 0;

#if 0
  int32                numBad   = 0;

  //  NOT TESTED

  for (int32 i=0; i<tigLen; i++) {
    if (badGraph[i] <= badMateBreakThreshold)
      numBad++;
  }

  //  If we found too many spots with bad, assume this is a large repeat and do not split anything.
  if (numBad * 4 > tigLen) {
    fprintf(stderr, "found %d bad spots out of tigLen %d; ignoring all bad regions\n",
            numBad, tigLen);
    return(peakBads);
  }
#endif

  for (int32 i=0; i<tigLen; i++) {
    if (badGraph[i] <= badMateBreakThreshold) {
      //  We are below the bad threshold, start a new bad region, or extend an existing one.

      if (badGraph[i] < peakBad) {
        //  Reset the bad region, we found one that is worse.
        peakBad  = badGraph[i];
        peak.bgn = i;
        peak.end = i;
      }

      if (badGraph[i] <= peakBad)
        //  Extend the bad region into this base.
        peak.end = i;

    } else {
      //  Else, we are above the bad threshold, save any existing bad region and reset.
      if (peakBad < 0) {
        peakBads->push_back(peak);

        peakBad  = 0;
        peak.bgn = 0;
        peak.end = 0;
      }
    }
  }

  //  If there is still a bad region on the stack, save it too.

  if (peakBad < 0)
    peakBads->push_back(peak);

  return(peakBads);
}




// hold over from testing if we should use 5' or 3' for range generation, now must use 3'
UnitigBreakPoints* UnitigGraph::computeMateCoverage(Unitig* tig,
                                                    int badMateBreakThreshold) {
  int tigLen = tig->getLength();

  MateLocation         positions(tig);
  vector<SeqInterval> *fwdBads = findPeakBad(positions.badFwdGraph, tigLen, badMateBreakThreshold);
  vector<SeqInterval> *revBads = findPeakBad(positions.badRevGraph, tigLen, badMateBreakThreshold);

  UnitigBreakPoints* breaks = new UnitigBreakPoints();

  if ((fwdBads->size() == 0) &&
      (revBads->size() == 0)) {
    delete fwdBads;
    delete revBads;
    //fprintf(logFile, "unitig %d no bad peaks\n", tig->id());
    return(breaks);
  }

  if (logFileFlagSet(LOG_MATE_SPLIT_ANALYSIS)) {
    fprintf(logFile, "unitig %d with "F_SIZE_T" fwd and "F_SIZE_T" rev bads\n",
            tig->id(), fwdBads->size(), revBads->size());
    fprintf(logFile, "fwd:");
    for (uint32 i=0; i<fwdBads->size(); i++)
      fprintf(logFile, " %d,%d", (*fwdBads)[i].bgn, (*fwdBads)[i].end);
    fprintf(logFile, "\n");
    fprintf(logFile, "rev:");
    for (uint32 i=0; i<revBads->size(); i++)
      fprintf(logFile, " %d,%d", (*revBads)[i].bgn, (*revBads)[i].end);
    fprintf(logFile, "\n");
  }

  if (logFileFlagSet(LOG_MATE_SPLIT_COVERAGE_PLOT)) {
    if ((tig->getLength() > 150) &&
        (tig->ufpath.size() > 3)) {
      char  filename[FILENAME_MAX] = {0};
      sprintf(filename, "coverageplot/utg%09u.badCoverage", tig->id());

      fprintf(logFile, "%s -- fwdBads "F_SIZE_T" revBads "F_SIZE_T"\n", filename, fwdBads->size(), revBads->size());

      if (AS_UTL_fileExists("coverageplot", TRUE, TRUE) == 0)
        AS_UTL_mkdir("coverageplot");

      FILE *F = fopen(filename, "w");

      for (uint32 i=0; i<tigLen; i++)
        fprintf(F, "%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                i,
                positions.goodGraph[i],
                positions.badFwdGraph[i],
                positions.badRevGraph[i],
                positions.badExternalFwd[i],
                positions.badExternalRev[i],
                positions.badCompressed[i],
                positions.badStretched[i],
                positions.badNormal[i],
                positions.badAnti[i],
                positions.badOuttie[i]);

      fclose(F);
    }
  }

  bool combine = false;
  int32 currBackbonePredecessorEnd = 0;
  int32 currBackboneEnd = 0;
  int32 lastBreakBBEnd = 0;

  uint32 frgidx = 0;

  ufNode       backbone = tig->getLastBackboneNode();
  int32        backBgn  = isReverse(backbone.position) ? backbone.position.end : backbone.position.bgn ;

  vector<SeqInterval>::const_iterator fwdIter = fwdBads->begin();
  vector<SeqInterval>::const_iterator revIter = revBads->begin();

  // Go through the peak bad ranges looking for reads to break on
  while(fwdIter != fwdBads->end() || revIter != revBads->end()) {
    bool isFwdBad = false;
    SeqInterval bad;
    if (revIter == revBads->end() ||
        fwdIter != fwdBads->end() &&  *fwdIter < *revIter) {
      // forward bad group, break at 1st frag
      isFwdBad = true;
      bad = *fwdIter;
      fwdIter++;
      if (lastBreakBBEnd >= bad.bgn) {
        // Skip, instead of combine trying to detect in combine case
        if (logFileFlagSet(LOG_MATE_SPLIT_ANALYSIS))
          fprintf(logFile,"Skip fwd bad range %d %d due to backbone %d\n",
                  bad.bgn, bad.end, lastBreakBBEnd);
        continue;
      }
    } else {                     // reverse bad group, break at last frag
      bad = *revIter;
      if (lastBreakBBEnd >= bad.bgn) {
        // Skip, instead of combine trying to detect in combine case
        if (logFileFlagSet(LOG_MATE_SPLIT_ANALYSIS))
          fprintf(logFile,"Skip rev bad range %d %d due to backbone %d\n",
                  bad.bgn, bad.end, lastBreakBBEnd);
        revIter++;
        continue;
      }
      if (fwdIter != fwdBads->end()) {
        if (fwdIter->bgn < bad.end && bad.end - fwdIter->bgn > 500) {
          // if fwd and reverse bad overlap
          // and end of reverse is far away, do fwd 1st
          isFwdBad = true;
          bad = *fwdIter;
          fwdIter++;
        } else {
          // check for containment relations and skip them
          if (fwdIter->bgn >= bad.bgn && fwdIter->end <= bad.end) {
            if (logFileFlagSet(LOG_MATE_SPLIT_ANALYSIS))
              fprintf(logFile,"Skip fwd bad range %d %d due to contained in rev %d %d\n",
                      fwdIter->bgn, fwdIter->end, bad.bgn, bad.end);

            fwdIter++;
            if (fwdIter == fwdBads->end()) {
              continue;
            }
          } else if (bad.bgn >= fwdIter->bgn && bad.end <= fwdIter->end) {
            if (logFileFlagSet(LOG_MATE_SPLIT_ANALYSIS))
              fprintf(logFile,"Skip rev bad range %d %d due to contained in fwd %d %d\n",
                      bad.bgn, bad.end, fwdIter->bgn, fwdIter->end);
      	    revIter++;
            continue;
          }

          if (fwdIter->bgn < bad.end &&
              fwdIter->end > bad.end &&
              bad.end - fwdIter->end < 200) {
            if (logFileFlagSet(LOG_MATE_SPLIT_ANALYSIS))
              fprintf(logFile,"Combine bad ranges %d - %d with %d - %d\n",
                      bad.bgn, bad.end, fwdIter->bgn, fwdIter->end);
            if (bad.bgn == 0) { // ignore reverse at start of tig
              bad.bgn = fwdIter->bgn;
              bad.end = fwdIter->end;
            } else {
              bad.bgn = bad.end;
              bad.end = fwdIter->bgn;
            }
            fwdIter++;
            combine = true;
          }
          revIter++;
        }
      } else {
        revIter++;
      }
    }

    if (logFileFlagSet(LOG_MATE_SPLIT_ANALYSIS))
      fprintf(logFile,"Bad peak from %d to %d\n",bad.bgn,bad.end);

    for (; frgidx < tig->ufpath.size(); frgidx++) {
      ufNode frag = tig->ufpath[frgidx];
      SeqInterval loc = frag.position;

      if (isReverse(loc)) {
        loc.bgn = frag.position.end;
        loc.end = frag.position.bgn;
      }

      if (logFileFlagSet(LOG_MATE_SPLIT_ANALYSIS))
        fprintf(logFile, "unitig %d frag %d %d,%d bad %d,%d\n",
                tig->id(), frag.ident, loc.bgn, loc.end, bad.bgn, bad.end);

      // keep track of current and previous uncontained contig end
      // so that we can split apart contained reads that don't overlap each other
      if (!OG->isContained(frag.ident)) {
        currBackbonePredecessorEnd = currBackboneEnd;
        currBackboneEnd = MAX(loc.bgn, loc.end);
      }

      bool breakNow = false;
      MateLocationEntry mloc = positions.getById(frag.ident);

      //  If we do go past the split point, whoops, break now and hope it all works out later.
      if ((loc.bgn > bad.end+1) && (loc.end > bad.end+1)) {
        fprintf(logFile, "SPLIT ERROR: unitig %d frag %d %d,%d missed the split point %d,%d.\n",
                tig->id(), frag.ident, loc.bgn, loc.end, bad.bgn, bad.end);
        breakNow = true;
      }
      // Don't want to go past range and break in wrong place
      assert(loc.bgn <= bad.end+1 || loc.end <= bad.end+1);


      if (mloc.mleFrgID1 != 0 && mloc.isGrumpy) { // only break on bad mates
        if (isFwdBad && bad.bgn <= loc.end) {
          breakNow = true;
        } else if (!isFwdBad && (loc.bgn >= /* bad.bgn*/ bad.end) ||
                   (combine && loc.end >  bad.bgn) ||
                   (combine && loc.end == bad.end)) {
          breakNow = true;
        } else if (bad.bgn > backBgn) {
          // fun special case, keep contained frags at end of tig in container
          // instead of in their own new tig where they might not overlap
          breakNow = true;
        }
      }

      bool noContainerToBreak = false;
      bool incrementToNextFragment = true;
      if (breakNow) {
        if (OG->isContained(frag.ident)) {
          // try to find a subsequent uncontained
          while ((frgidx < tig->ufpath.size()) &&
                 (OG->isContained(tig->ufpath[frgidx].ident) == true))
            frgidx++;
          // we couldn't find a downstream frag take upstream one then
          if (frgidx >= tig->ufpath.size()) {
            do {
              frgidx--;
            } while ((frgidx > 0) &&
                     (OG->isContained(tig->ufpath[frgidx].ident) == true));
            if ((frgidx == 0) && (OG->isContained(tig->ufpath[0].ident) == true))
              noContainerToBreak = true;
            else
              currBackboneEnd = currBackbonePredecessorEnd;
          }
          //  This is either a bug, or a bug fix.  The original version would decrement frgidx above
          //  (which, at that time was an iterator) but NOT update the copy of 'frag' or 'loc' (as
          //  we are doing below).  The next block below would then test with the old value of loc.
#warning bug fix or bug
          frag = tig->ufpath[frgidx];
          loc = frag.position;
        }
        if (noContainerToBreak == false) {
          combine = false;
          lastBreakBBEnd = currBackboneEnd;
          if (logFileFlagSet(LOG_MATE_SPLIT_ANALYSIS))
            fprintf(logFile,"Frg to break in peak bad range is %d fwd %d pos (%d,%d) backbone %d\n",
                    frag.ident, isFwdBad, loc.bgn, loc.end, currBackboneEnd);
          uint32 frag3p = true;
          // If reverse mate is 1st and overlaps its mate break at 5'
          if (mloc.mleUtgID2 == tig->id() && isReverse(loc) &&
              !isReverse(mloc.mlePos2) && loc.bgn >= mloc.mlePos2.bgn)
            frag3p = false;

          // either adjust the break position to be before the container so the container travels together with containees
          if (frag3p == false && !isReverse(frag.position)||
              frag3p == true  &&  isReverse(frag.position)) {
            // do nothing we are breaking before the current fragment which is our container
            incrementToNextFragment = false;
          } else {
            // move one back we are breaking after the current fragment
            //  UGH, this is gross nasty code
            if (frgidx > 0)
              frgidx--;
            else
              fprintf(logFile, "DECREMENT ZERO frgidx!\n");
          }
          frag = tig->ufpath[frgidx];
          loc = frag.position;

          UnitigBreakPoint bp(frag.ident, frag3p);
          bp.fragPos = frag.position;
          bp.inSize = 100000;
          bp.inFrags = 10;
          if (logFileFlagSet(LOG_MATE_SPLIT_COVERAGE_PLOT))
            fprintf(logFile, "BREAK unitig %d at fragment %d position %d,%d from MATES #1.\n",
                    tig->id(), frag.ident, bp.fragPos.bgn, bp.fragPos.end);
          breaks->push_back(bp);
        }
      }

      if (lastBreakBBEnd != 0 && lastBreakBBEnd > MAX(loc.bgn,loc.end)) {
        uint32  np = frgidx + 1;

        if (np < tig->ufpath.size()) {
          if (contains(loc, tig->ufpath[np].position)) {
            // Contains the next one, so skip it
          } else {
            SeqInterval overlap = intersection(loc, tig->ufpath[np].position);
            int diff = abs(overlap.end - overlap.bgn);

            //  No overlap between this and the next
            //  frag, or the overlap is tiny, or this
            //  frag is contained, but not contained
            //  in the next frag; Break after this
            //  frag.
            //
            if ((NULL_SEQ_LOC == overlap) ||
                (diff < AS_OVERLAP_MIN_LEN) ||
                (OG->isContained(frag.ident) && !OG->containHaveEdgeTo(frag.ident, tig->ufpath[np].ident))) {

              uint32 frag3p = true;
              if (isReverse(loc))
                frag3p = false;

              UnitigBreakPoint bp(frag.ident, frag3p);
              bp.fragPos = loc;
              bp.inSize = 100001;
              bp.inFrags = 11;
              if (logFileFlagSet(LOG_MATE_SPLIT_COVERAGE_PLOT))
                fprintf(logFile, "BREAK unitig %d at fragment %d position %d,%d from MATES #2.\n",
                        tig->id(), frag.ident, bp.fragPos.bgn, bp.fragPos.end);
              breaks->push_back(bp);
              if (logFileFlagSet(LOG_MATE_SPLIT_ANALYSIS))
                fprintf(logFile,"Might make frg %d singleton, end %d size %u pos %d,%d\n",
                        frag.ident, frag3p, (uint32)breaks->size(), loc.bgn, loc.end);
            }
          }
        }
      }
      if (breakNow) { // Move to next breakpoint
        if (incrementToNextFragment) {
          frgidx++;  // make sure to advance past curr frg
          frag = tig->ufpath[frgidx];
        }
        break;
      }
    }
  }

  delete fwdBads;
  delete revBads;

  return breaks;
}

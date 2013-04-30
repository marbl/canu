
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

static const char *rcsid = "$Id: AS_BAT_BestOverlapGraph.C,v 1.13 2012-11-20 13:20:15 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Unitig.H"

#include "AS_UTL_intervalList.H"



BestOverlapGraph::BestOverlapGraph(double               AS_UTG_ERROR_RATE,
                                   double               AS_UTG_ERROR_LIMIT,
                                   const char          *prefix,
                                   bool                 removeSuspicious,
                                   bool                 removeSpurs) {

  setLogFile(prefix, "bestOverlapGraph");

  writeLog("BestOverlapGraph-- allocating best edges ("F_SIZE_T"MB) and containments ("F_SIZE_T"MB)\n",
           ((2 * sizeof(BestEdgeOverlap) * (FI->numFragments() + 1)) >> 20),
           ((1 * sizeof(BestContainment) * (FI->numFragments() + 1)) >> 20));

  assert(AS_UTG_ERROR_RATE >= 0.0);
  assert(AS_UTG_ERROR_RATE <= AS_MAX_ERROR_RATE);

  assert(AS_CNS_ERROR_RATE >= 0.0);
  assert(AS_CNS_ERROR_RATE <= AS_MAX_ERROR_RATE);

  _bestA = new BestOverlaps [FI->numFragments() + 1];
  _scorA = new BestScores   [FI->numFragments() + 1];

  memset(_bestA, 0, sizeof(BestOverlaps) * (FI->numFragments() + 1));
  memset(_scorA, 0, sizeof(BestScores)   * (FI->numFragments() + 1));

  _restrict = NULL;

  mismatchCutoff  = AS_UTG_ERROR_RATE;
  mismatchLimit   = AS_UTG_ERROR_LIMIT;

  //  Initialize parallelism.

  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  //  PASS 0:  Find suspicious fragments.  For any found, mark as suspicious and don't allow
  //  these to be best overlaps.

  if (removeSuspicious) {
    writeLog("BestOverlapGraph()-- removing suspicious reads from graph, with %d threads.\n", numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
    for (AS_IID fi=1; fi <= fiLimit; fi++) {
      uint32      no  = 0;
      BAToverlap *ovl = OC->getOverlaps(fi, no);

      bool          verified = false;
      intervalList  IL;

      uint32        fLen = FI->fragmentLength(fi);

      for (uint32 ii=0; (ii<no) && (verified == false); ii++) {
        if (isOverlapBadQuality(ovl[ii]))
          //  Yuck.  Don't want to use this crud.
          continue;

        if      ((ovl[ii].a_hang <= 0) && (ovl[ii].b_hang <= 0))
          //  Left side dovetail
          IL.add(0, fLen + ovl[ii].b_hang);

        else if ((ovl[ii].a_hang >= 0) && (ovl[ii].b_hang >= 0))
          //  Right side dovetail
          IL.add(ovl[ii].a_hang, fLen - ovl[ii].a_hang);

        else if ((ovl[ii].a_hang >= 0) && (ovl[ii].b_hang <= 0))
          //  I contain the other
          IL.add(ovl[ii].a_hang, fLen - ovl[ii].a_hang - ovl[ii].b_hang);

        else if ((ovl[ii].a_hang <= 0) && (ovl[ii].b_hang >= 0))
          //  I am contained and thus now perfectly good!
          verified = true;

        else
          //  Huh?  Coding error.
          assert(0);
      }

      if (verified == false) {
        IL.merge();
        verified = (IL.numberOfIntervals() == 1);
      }

      if (verified == false) {
        if (no > 0)
          writeLog("BestOverlapGraph()-- frag "F_U32" is suspicious ("F_U32" overlaps).\n", fi, no);

#pragma omp critical (suspInsert)
        _suspicious.insert(fi);
      }
    }
  }

  //  PASS 1:  Find containments.

  writeLog("BestOverlapGraph()-- analyzing %d fragments for best contains, with %d threads.\n", fiLimit, numThreads);

  if (false && AS_UTL_fileExists("best.contains", false, false)) {
    writeLog("BestOverlapGraph()-- loading best containes from cache.\n");
    assert(0);  //  Not quite done.

  } else {
#pragma omp parallel for schedule(dynamic, blockSize)
    for (AS_IID fi=1; fi <= fiLimit; fi++) {
      uint32      no  = 0;
      BAToverlap *ovl = OC->getOverlaps(fi, no);

      for (uint32 ii=0; ii<no; ii++)
        scoreContainment(ovl[ii]);
    }
  }

  //  PASS 2:  Find dovetails.

  if (false && AS_UTL_fileExists("best.edges", false, false)) {
    writeLog("BestOverlapGraph()-- loading best edges from cache.\n");
    assert(0);  //  Not quite done.

  } else {
    writeLog("BestOverlapGraph()-- analyzing %d fragments for best edges, with %d threads.\n", fiLimit, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
    for (AS_IID fi=1; fi <= fiLimit; fi++) {
      uint32      no  = 0;
      BAToverlap *ovl = OC->getOverlaps(fi, no);

      for (uint32 ii=0; ii<no; ii++)
        scoreEdge(ovl[ii]);
    }
  }

  //  Remove spurs

  if (removeSpurs) {
    writeLog("BestOverlapGraph()-- detecting spur fragments.\n");

    char   *isSpur = new char [fiLimit + 1];

    memset(isSpur, 0, sizeof(char) * (fiLimit + 1));

    for (AS_IID fi=1; fi <= fiLimit; fi++) {
      bool   spur5 = (getBestEdgeOverlap(fi, false)->fragId() == 0);
      bool   spur3 = (getBestEdgeOverlap(fi, true)->fragId()  == 0);

      if ((spur5 == true) || (spur3 == true)) {
        if ((spur5 == false) || (spur3 == false))
          writeLog("BestOverlapGraph()-- frag "F_U32" is a %s spur.\n", fi, (spur5) ? "5'" : "3'");
        isSpur[fi] = true;
      }
    }

    //  Remove best edges, so we can rebuild

    memset(_bestA, 0, sizeof(BestOverlaps) * (fiLimit + 1));
    memset(_scorA, 0, sizeof(BestScores)   * (fiLimit + 1));

    //  Rebuild best edges, ignoring edges to spurs.  We build edges out of spurs, but don't allow edges into them.
    //  This should prevent them from being incorporated into a promiscuous unitig, but still let them be popped
    //  as bubbles (but they shouldn't because they're spurs).

    //  PASS 3:  Find containments.

    writeLog("BestOverlapGraph()-- analyzing %d fragments for best contains, with %d threads.\n", fiLimit, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
    for (AS_IID fi=1; fi <= fiLimit; fi++) {
      uint32      no  = 0;
      BAToverlap *ovl = OC->getOverlaps(fi, no);

      for (uint32 ii=0; ii<no; ii++)
        scoreContainment(ovl[ii]);
    }

    //  PASS 4:  Find dovetails.

    writeLog("BestOverlapGraph()-- analyzing %d fragments for best edges, with %d threads.\n", fiLimit, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
    for (AS_IID fi=1; fi <= fiLimit; fi++) {
      uint32      no  = 0;
      BAToverlap *ovl = OC->getOverlaps(fi, no);

      for (uint32 ii=0; ii<no; ii++)
        if (isSpur[ovl[ii].b_iid] == false)
          scoreEdge(ovl[ii]);
    }

    delete [] isSpur;
  }

  //  Remove temporary scoring data

  delete [] _scorA;
  _scorA = NULL;

  //  Remove dovetail overlaps for contained fragments.

  writeLog("BestOverlapGraph()-- removing best edges for contained fragments, with %d threads.\n", numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (AS_IID fi=1; fi <= fiLimit; fi++) {
    if (isContained(fi) == true) {
      getBestEdgeOverlap(fi, false)->set(0, 0, 0, 0);
      getBestEdgeOverlap(fi, true) ->set(0, 0, 0, 0);
    }
  }

  writeLog("BestOverlapGraph()-- dumping best edges/contains/singletons.\n");

  reportBestEdges();

  setLogFile(prefix, NULL);
}




BestOverlapGraph::BestOverlapGraph(double               AS_UTG_ERROR_RATE,
                                   double               AS_UTG_ERROR_LIMIT,
                                   set<AS_IID>         *restrict) {

  assert(AS_UTG_ERROR_RATE >= 0.0);
  assert(AS_UTG_ERROR_RATE <= AS_MAX_ERROR_RATE);

  assert(AS_CNS_ERROR_RATE >= 0.0);
  assert(AS_CNS_ERROR_RATE <= AS_MAX_ERROR_RATE);

  mismatchCutoff  = AS_UTG_ERROR_RATE;
  mismatchLimit   = AS_UTG_ERROR_LIMIT;

  _bestA = NULL;
  _scorA = NULL;

  _bestM.clear();
  _scorM.clear();

  _restrict = restrict;

  //  PASS 0:  Load the map (necessary?)

#if 0
  for (set<AS_IID>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    AS_IID      fi  = *it;

    _bestM[fi].insert();
    _scorM[fi].insert();
  }
#endif

  //  PASS 1:  Find containments.

  for (set<AS_IID>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    AS_IID      fi  = *it;
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, no);

    for (uint32 ii=0; ii<no; ii++)
      scoreContainment(ovl[ii]);
  }

  //  PASS 2:  Find dovetails.

  for (set<AS_IID>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    AS_IID      fi  = *it;
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, no);

    for (uint32 ii=0; ii<no; ii++)
      scoreEdge(ovl[ii]);
  }

  //  Remove temporary scoring data

  _scorM.clear();

  //  Remove dovetail overlaps for contained fragments.

  for (set<AS_IID>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    AS_IID      fi  = *it;

    if (isContained(fi) == true) {
      getBestEdgeOverlap(fi, false)->set(0, 0, 0, 0);
      getBestEdgeOverlap(fi, true) ->set(0, 0, 0, 0);
    }
  }

  //  Remove spurs

#if 0
  for (set<AS_IID>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    AS_IID      fi  = *it;

    if ((getBestEdgeOverlap(fi, false)->fragId() == 0) ||
        (getBestEdgeOverlap(fi, true)->fragId()  == 0)) {
      getBestEdgeOverlap(fi, false)->set(0, 0, 0, 0);
      getBestEdgeOverlap(fi, true) ->set(0, 0, 0, 0);
    }
  }
#endif
}
















void
BestOverlapGraph::reportBestEdges(void) {
  FILE *BC = fopen("best.contains", "w");
  FILE *BE = fopen("best.edges", "w");
  FILE *BS = fopen("best.singletons", "w");

  if ((BC) && (BE) && (BS)) {
    fprintf(BC, "#fragId\tlibId\tmated\tbestCont\n");
    fprintf(BE, "#fragId\tlibId\tbest5\tbest3\n");
    fprintf(BS, "#fragId\tlibId\tmated\n");

    for (uint32 id=1; id<FI->numFragments() + 1; id++) {
      BestContainment *bestcont  = getBestContainer(id);
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      if (bestcont->isContained)
        fprintf(BC, "%u\t%u\t%c\t%u\n", id, FI->libraryIID(id), (FI->mateIID(id) > 0) ? 'm' : 'f', bestcont->container);
      else if ((bestedge5->fragId() > 0) || (bestedge3->fragId() > 0))
        fprintf(BE, "%u\t%u\t%u\t%c'\t%u\t%c'\n", id, FI->libraryIID(id),
                bestedge5->fragId(), bestedge5->frag3p() ? '3' : '5',
                bestedge3->fragId(), bestedge3->frag3p() ? '3' : '5');
      else
        fprintf(BS, "%u\t%u\t%c\n", id, FI->libraryIID(id), (FI->mateIID(id) > 0) ? 'm' : 'f');

    }

    fclose(BC);
    fclose(BE);
    fclose(BS);
  }
}






void
BestOverlapGraph::scoreContainment(const BAToverlap& olap) {

  if (isOverlapBadQuality(olap))
    //  Yuck.  Don't want to use this crud.
    return;

  if (isOverlapRestricted(olap))
    //  Whoops, don't want this overlap for this BOG
    return;

  if ((olap.a_hang == 0) &&
      (olap.b_hang == 0) &&
      (olap.a_iid > olap.b_iid))
    //  Exact!  Each contains the other.  Make the lower IID the container.
    return;

  if ((olap.a_hang > 0) ||
      (olap.b_hang < 0))
    //  We only save if A is the contained fragment.
    return;

  uint64           newScr = scoreOverlap(olap);

  assert(newScr > 0);

  //  The previous version (1.5) saved if A contained B.  This was breaking the overlap filtering,
  //  because long A fragments containing short B fragments would have those containment overlaps
  //  filtered out.  Version 1.6 reversed what is saved here so that the containment overlap is
  //  associated with the A fragment (the containee).
  //
  //  The hangs will transform the container coordinates into the containee cordinates.

  if (newScr > bestCscore(olap.a_iid)) {
    BestContainment   *c = getBestContainer(olap.a_iid);

    c->container         = olap.b_iid;
    c->isContained       = true;
    c->sameOrientation   = olap.flipped ? false : true;
    c->a_hang            = olap.flipped ? olap.b_hang : -olap.a_hang;
    c->b_hang            = olap.flipped ? olap.a_hang : -olap.b_hang;

    //writeLog("set best for %d from "F_U64" to "F_U64"\n", olap.a_iid, bestCscore(olap.a_iid), newScr);

    bestCscore(olap.a_iid) = newScr;
  }
}



void
BestOverlapGraph::scoreEdge(const BAToverlap& olap) {
  bool   enableLog = false;  //  useful for reporting this stuff only for specific reads

  //if ((olap.a_iid == 97202) || (olap.a_iid == 30701))
  //  enableLog = true;

  if (isOverlapBadQuality(olap)) {
    //  Yuck.  Don't want to use this crud.
    if ((enableLog == true) && ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY)))
      writeLog("OVERLAP BADQ:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- bad quality\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.error);
    return;
  }

  if (isOverlapRestricted(olap)) {
    //  Whoops, don't want this overlap for this BOG
    if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
      writeLog("OVERLAP REST:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- restricted\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.error);
    return;
  }

  if (isSuspicious(olap.b_iid)) {
    //  Whoops, don't want this overlap for this BOG
    if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
      writeLog("OVERLAP SUSP:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- suspicious\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.error);
    return;
  }

  if (((olap.a_hang >= 0) && (olap.b_hang <= 0)) ||
      ((olap.a_hang <= 0) && (olap.b_hang >= 0))) {
    //  Skip containment overlaps.
    if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
      writeLog("OVERLAP CONT:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- container read\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.error);
    return;
  }

  if ((isContained(olap.a_iid) == true) ||
      (isContained(olap.b_iid) == true)) {
    //  Skip contained fragments.
    if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
      writeLog("OVERLAP CONT:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- contained read\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.error);
    return;
  }

  uint64           newScr = scoreOverlap(olap);
  bool             a3p    = AS_BAT_overlapAEndIs3prime(olap);
  BestEdgeOverlap *best   = getBestEdgeOverlap(olap.a_iid, a3p);
  uint64          &score  = (a3p) ? (best3score(olap.a_iid)) : (best5score(olap.a_iid));

  assert(newScr > 0);

  if (newScr <= score)
    return;

  best->set(olap);

  score = newScr;

  if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
    writeLog("OVERLAP GOOD:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- NOW BEST\n",
             olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.error);
}




bool
BestOverlapGraph::isOverlapBadQuality(const BAToverlap& olap) {
  bool   enableLog = false;  //  useful for reporting this stuff only for specific reads

  //if ((olap.a_iid == 97202) || (olap.a_iid == 30701))
  //  enableLog = true;

  if ((FI->fragmentLength(olap.a_iid) == 0) ||
      (FI->fragmentLength(olap.b_iid) == 0))
    //  The overlap is bad if it involves deleted fragments.  Shouldn't happen in a normal
    //  assembly, but sometimes us users want to delete fragments after overlaps are generated.
    return(true);

  //  The overlap is GOOD (false == not bad) if the corrected error rate is below the requested
  //  erate.
  //
  if (olap.error <= mismatchCutoff) {
    if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
      writeLog("OVERLAP GOOD:     %d %d %c  hangs "F_S32" "F_S32" err %.3f\n",
               olap.a_iid, olap.b_iid,
               olap.flipped ? 'A' : 'N',
               olap.a_hang,
               olap.b_hang,
               olap.error);
    return(false);
  }

  //  If we didn't allow fixed-number-of-errors, the overlap is now bad.  Just a slight
  //  optimization.
  //
  if (mismatchLimit <= 0)
    return(true);

  //  There are a few cases where the orig_erate is _better_ than the corr_erate.  That is, the
  //  orig erate is 0% but we 'correct' it to something more than 0%.  Regardless, we probably
  //  want to be using the corrected erate here.

  double olen = FI->overlapLength(olap.a_iid, olap.b_iid, olap.a_hang, olap.b_hang);
  double nerr = olen * olap.error;

  assert(nerr >= 0);

  if (nerr <= mismatchLimit) {
    if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
      writeLog("OVERLAP SAVED:    %d %d %c  hangs "F_S32" "F_S32" err %.3f olen %f nerr %f\n",
               olap.a_iid, olap.b_iid,
               olap.flipped ? 'A' : 'N',
               olap.a_hang,
               olap.b_hang,
               olap.error,
               olen, nerr);
    return(false);
  }

  if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
    writeLog("OVERLAP REJECTED: %d %d %c  hangs "F_S32" "F_S32" err %.3f olen %f nerr %f\n",
             olap.a_iid, olap.b_iid,
             olap.flipped ? 'A' : 'N',
             olap.a_hang,
             olap.b_hang,
             olap.error,
             olen, nerr);
  return(true);
}


//  If no restrictions are known, this overlap is useful if both fragments are not in a unitig
//  already.  Otherwise, we are restricted to just a specific set of fragments (usually a single
//  unitig and all the mated reads).  The overlap is useful if both fragments are in the set.
//
bool
BestOverlapGraph::isOverlapRestricted(const BAToverlap &olap) {

  if (_restrict == NULL) {
    if ((Unitig::fragIn(olap.a_iid) == 0) &&
        (Unitig::fragIn(olap.b_iid) == 0))
      return(false);
    else
      return(true);
  }

  if ((_restrict->count(olap.a_iid) != 0) &&
      (_restrict->count(olap.b_iid) != 0))
    return(false);
  else
    return(true);
}


uint64
BestOverlapGraph::scoreOverlap(const BAToverlap& olap) {

  //  BPW's newer new score.  For the most part, we use the length of the overlap, but we also
  //  want to break ties with the higher quality overlap:
  //    The high bits are the length of the overlap.
  //    The next are the corrected error rate.
  //    The last are the original error rate.
  //
  uint64  leng = 0;
  uint64  corr = AS_BAT_MAX_ERATE - olap.errorRaw;
  uint64  orig = AS_BAT_MAX_ERATE - 0;

  //  Shift AFTER assigning to a 64-bit value to avoid overflows.
  corr <<= AS_OVS_ERRBITS;

  //  Containments - the length of the overlaps are all the same.  We return the quality.
  //
  if (((olap.a_hang >= 0) && (olap.b_hang <= 0)) ||
      ((olap.a_hang <= 0) && (olap.b_hang >= 0)))
    return(corr | orig);

  //  Dovetails - the length of the overlap is the score, but we bias towards lower error.
  //  (again, shift AFTER assigning to avoid overflows)
  //
  leng   = FI->overlapLength(olap.a_iid, olap.b_iid, olap.a_hang, olap.b_hang);
  leng <<= (2 * AS_OVS_ERRBITS);

  return(leng | corr | orig);
}

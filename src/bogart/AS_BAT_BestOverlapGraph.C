
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
 *    src/AS_BAT/AS_BAT_BestOverlapGraph.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2014-JAN-29
 *      are Copyright 2010-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2015-AUG-14
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Unitig.H"

#include "intervalList.H"
#include "stddev.H"



void
getOverlapForBestEdge(uint32 fi, BestEdgeOverlap *edge, BAToverlap &olap) {
  uint32       no    = 0;
  BAToverlap  *ovl   = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

  for (uint32 ii=0; ii<no; ii++)
    if (ovl[ii].b_iid == edge->fragId())
      olap = ovl[ii];
}



void
getOverlapForBestEdges(uint32 fi,
                       BestEdgeOverlap *e5, BAToverlap &ovl5,
                       BestEdgeOverlap *e3, BAToverlap &ovl3) {
  uint32       no    = 0;
  BAToverlap  *ovl   = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

  for (uint32 ii=0; ii<no; ii++) {
    if (ovl[ii].b_iid == e5->fragId())
      ovl5 = ovl[ii];

    if (ovl[ii].b_iid == e3->fragId())
      ovl3 = ovl[ii];
  }
}



void
BestOverlapGraph::removeSuspicious(void) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  writeLog("BestOverlapGraph()-- removing suspicious reads from graph, with %d threads.\n", numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32               no  = 0;
    BAToverlap          *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    bool                 verified = false;
    intervalList<int32>  IL;

    uint32               fLen = FI->fragmentLength(fi);

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



void
BestOverlapGraph::removeHighErrorBestEdges(void) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  stdDev<double>  edgeStats;

  BAToverlap   b5ovl = { 0, 0, 0, 0, 0, 0, 0};
  BAToverlap   b3ovl = { 0, 0, 0, 0, 0, 0, 0};

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap *b5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap *b3 = getBestEdgeOverlap(fi, true);

    getOverlapForBestEdges(fi,
                           b5, b5ovl,
                           b3, b3ovl);

    if (b5->fragId() != 0)
      edgeStats.update(b5ovl.erate);

    if (b3->fragId() != 0)
      edgeStats.update(b3ovl.erate);
  }

  writeLog("removeHighErrorBestEdges()-- with %u points - mean %f stddev %f\n",
           edgeStats.size(), edgeStats.mean(), edgeStats.stddev());

  _mean   = edgeStats.mean();
  _stddev = edgeStats.stddev();

  uint32  retained = 0;
  uint32  removed  = 0;

  //FILE *F = fopen("test.best.edges.erates", "w");

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap *b5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap *b3 = getBestEdgeOverlap(fi, true);

    b5ovl.erate = 1.0;
    b3ovl.erate = 1.0;

    getOverlapForBestEdges(fi,
                           b5, b5ovl,
                           b3, b3ovl);

    //fprintf(F, "%u 5' %f 3' %f\n", fi, b5ovl.erate, b3ovl.erate);

    if (b5->fragId() != 0) {
      if (b5ovl.erate > _mean + 5 * _stddev)
        removed++;
      else
        retained++;
    }

    if (b3->fragId() != 0) {
      if (b3ovl.erate > _mean + 5 * _stddev)
        removed++;
      else
        retained++;
    }
  }

  //fclose(F);

  fprintf(stderr, "removeHighErrorBestEdges()-- %u edges are above threshold, %u edges are acceptable.\n",
          removed, retained);
}



void
BestOverlapGraph::removeLopsidedEdges(void) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  writeLog("BestOverlapGraph()-- removing suspicious edges from graph, with %d threads.\n", numThreads);

  uint32   nSuspicious = 0;
  uint32   nContained  = 0;
  uint32   nSpur       = 0;
  uint32   nMutual     = 0;
  uint32   nAccepted   = 0;
  uint32   nRejected   = 0;

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap *this5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap *this3 = getBestEdgeOverlap(fi, true);

    //  Ignore spurs and contains...and previously detected suspicious reads.  The suspicious reads
    //  do not have best edges back to them, and it's possible to find reads B where best edge A->B
    //  exists, yet no best edge from B exists.

    if (isSuspicious(fi) == true) {
      nSuspicious++;
      continue;
    }

    if (isContained(fi) == true) {
      nContained++;
      continue;
    }

    if ((this5->fragId() == 0) ||
        (this3->fragId() == 0)) {
      nSpur++;
      continue;
    }

    //  Find the overlap for this5 and this3.

    BAToverlap   this5ovl = {0,0,0,0,0,0,0};
    BAToverlap   this3ovl = {0,0,0,0,0,0,0};

    getOverlapForBestEdges(fi,
                           this5, this5ovl,
                           this3, this3ovl);

    //  Compute the length of each overlap.

    int32   this5ovlLen = FI->overlapLength(this5ovl.a_iid, this5ovl.b_iid, this5ovl.a_hang, this5ovl.b_hang);
    int32   this3ovlLen = FI->overlapLength(this3ovl.a_iid, this3ovl.b_iid, this3ovl.a_hang, this3ovl.b_hang);

    //  Find the edges for our best overlaps.

    BestEdgeOverlap *that5 = getBestEdgeOverlap(this5->fragId(), this5->frag3p());
    BestEdgeOverlap *that3 = getBestEdgeOverlap(this3->fragId(), this3->frag3p());

    //  If both point back to us, we're done.

    if ((that5->fragId() == fi) && (that5->frag3p() == false) &&
        (that3->fragId() == fi) && (that3->frag3p() == true)) {
      nMutual++;
      continue;
    }

    //  Both are supposed to exist, but possibly we could be removing some best edge.

    //if ((that5->fragId() == 0) ||
    //    (that3->fragId() == 0)) {
    //  continue;
    //}

    //  Something doesn't agree.  Find those overlaps...

    BAToverlap   that5ovl = {0,0,0,0,0,0,0};
    BAToverlap   that3ovl = {0,0,0,0,0,0,0};

    getOverlapForBestEdge(this5->fragId(), that5, that5ovl);
    getOverlapForBestEdge(this3->fragId(), that3, that3ovl);

    fprintf(stderr, "fi %u this5 %u that5 %u\n", fi, this5->fragId(), that5->fragId());
    fprintf(stderr, "fi %u this3 %u that3 %u\n", fi, this3->fragId(), that3->fragId());

    assert(that5ovl.a_iid != 0);
    assert(that3ovl.a_iid != 0);

    //  ...and their lengths...

    int32  that5ovlLen = FI->overlapLength(that5ovl.a_iid, that5ovl.b_iid, that5ovl.a_hang, that5ovl.b_hang);
    int32  that3ovlLen = FI->overlapLength(that3ovl.a_iid, that3ovl.b_iid, that3ovl.a_hang, that3ovl.b_hang);

    //  ...and compare.

    double  percDiff5 = 200.0 * abs(this5ovlLen - that5ovlLen) / (this5ovlLen + that5ovlLen);
    double  percDiff3 = 200.0 * abs(this3ovlLen - that3ovlLen) / (this3ovlLen + that3ovlLen);

    if ((percDiff5 <= 5) &&
        (percDiff3 <= 5)) {
#if 0
      writeLog("fi %8u -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%% -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%% -- ACCEPTED\n",
               fi,
               this5->fragId(), this5->frag3p() ? '3' : '5', this5ovlLen, that5->fragId(), that5->frag3p() ? '3' : '5', that5ovlLen, percDiff5,
               this3->fragId(), this3->frag3p() ? '3' : '5', this3ovlLen, that3->fragId(), that3->frag3p() ? '3' : '5', that3ovlLen, percDiff3);
#endif
      nAccepted++;

    } else {
#if 0
      writeLog("fi %8u -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%% -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%%\n",
               fi,
               this5->fragId(), this5->frag3p() ? '3' : '5', this5ovlLen, that5->fragId(), that5->frag3p() ? '3' : '5', that5ovlLen, percDiff5,
               this3->fragId(), this3->frag3p() ? '3' : '5', this3ovlLen, that3->fragId(), that3->frag3p() ? '3' : '5', that3ovlLen, percDiff3);
#endif
      nRejected++;

#pragma omp critical (suspInsert)
      _suspicious.insert(fi);
    }
  }

  writeLog("BestOverlapGraph()--  suspicious %u  contained %u  spur %u  mutual-best %u  accepted %u  rejected %u\n",
           nSuspicious, nContained, nSpur, nMutual, nAccepted, nRejected);
}



char *
BestOverlapGraph::removeSpurs(void) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  writeLog("BestOverlapGraph()-- detecting spur fragments.\n");

  char   *isSpur = new char [fiLimit + 1];

  memset(isSpur, 0, sizeof(char) * (fiLimit + 1));

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    bool   spur5 = (getBestEdgeOverlap(fi, false)->fragId() == 0);
    bool   spur3 = (getBestEdgeOverlap(fi, true)->fragId()  == 0);

    if (isContained(fi)) {
      //writeLog("BestOverlapGraph()-- frag "F_U32" is contained - %d %d.\n", fi, spur5, spur3);
      continue;
    }

    if ((spur5 == false) && (spur3 == false))
      //  Edges off of both ends.  Not a spur.
      continue;

    if ((spur5 == true)  && (spur3 == true))
      //  No edges off either end.  Not a spur, just garbage.  EXCEPT that this could also
      //  be a contained read that has no 5'/3' best edges assigned to it.
      //writeLog("BestOverlapGraph()-- frag "F_U32" is singleton, no best edges and not contained.\n", fi);
      continue;

    //  Exactly one end is missing a best edge.  Bad!

    writeLog("BestOverlapGraph()-- frag "F_U32" is a %s spur.\n", fi, (spur5) ? "5'" : "3'");
    isSpur[fi] = true;
  }

  return(isSpur);
}



void
BestOverlapGraph::findEdges(char *isSpur) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  memset(_bestA, 0, sizeof(BestOverlaps) * (fiLimit + 1));
  memset(_scorA, 0, sizeof(BestScores)   * (fiLimit + 1));

  writeLog("BestOverlapGraph()-- analyzing %d fragments for best contains, with %d threads.\n", fiLimit, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    for (uint32 ii=0; ii<no; ii++)
      scoreContainment(ovl[ii]);
  }

  writeLog("BestOverlapGraph()-- analyzing %d fragments for best edges, with %d threads.\n", fiLimit, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    //  Build edges out of spurs, but don't allow edges into them.  This should prevent them from
    //  being incorporated into a promiscuous unitig, but still let them be popped as bubbles (but
    //  they shouldn't because they're spurs).

    for (uint32 ii=0; ii<no; ii++)
      if ((isSpur == NULL) || (isSpur[ovl[ii].b_iid] == false))
        scoreEdge(ovl[ii]);
  }

  delete [] isSpur;
}



void
BestOverlapGraph::removeContainedDovetails(void) {
  uint32  fiLimit    = FI->numFragments();

  writeLog("BestOverlapGraph()-- removing best edges for contained fragments.\n");

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    if (isContained(fi) == true) {
      getBestEdgeOverlap(fi, false)->set(0, 0, 0, 0);
      getBestEdgeOverlap(fi, true) ->set(0, 0, 0, 0);
    }
  }
}



BestOverlapGraph::BestOverlapGraph(double               erate,
                                   const char          *prefix) {

  setLogFile(prefix, "bestOverlapGraph");

  writeLog("BestOverlapGraph-- allocating best edges ("F_SIZE_T"MB) and containments ("F_SIZE_T"MB)\n",
           ((2 * sizeof(BestEdgeOverlap) * (FI->numFragments() + 1)) >> 20),
           ((1 * sizeof(BestContainment) * (FI->numFragments() + 1)) >> 20));

  //  Initialize data.

  _bestA           = new BestOverlaps [FI->numFragments() + 1];
  _scorA           = new BestScores   [FI->numFragments() + 1];

  _mean            = erate;
  _stddev          = 0.0;

  _suspicious.clear();

  _bestM.clear();
  _scorM.clear();

  _restrict        = NULL;
  _restrictEnabled = false;

  _erate           = erate;

  memset(_bestA, 0, sizeof(BestOverlaps) * (FI->numFragments() + 1));
  memset(_scorA, 0, sizeof(BestScores)   * (FI->numFragments() + 1));

  //  Compute best overlaps.

  removeSuspicious();
  findEdges();
  reportBestEdges("TEST1");

  removeHighErrorBestEdges();
  findEdges();
  reportBestEdges("TEST2");

  removeLopsidedEdges();
  findEdges();
  reportBestEdges("TEST3");

  findEdges(removeSpurs());
  reportBestEdges("TEST4");

  delete [] _scorA;  //  Done with the scoring data.
  _scorA = NULL;

  removeContainedDovetails();

  reportBestEdges(prefix);

  setLogFile(prefix, NULL);
}









void
BestOverlapGraph::rebuildBestContainsWithoutSingletons(UnitigVector  &unitigs,
                                                       double         erate,
                                                       const char    *prefix) {

  _erate = erate;

  uint32     fiLimit = FI->numFragments();

  assert(_restrict        == NULL);
  assert(_restrictEnabled == false);

  //  Save the current best containments for a nice log, then clear

  assert(_bestA != NULL);

  BestContainment  *bestCold = new BestContainment [fiLimit + 1];

  for (uint32 fi=0; fi<=fiLimit; fi++) {
    bestCold[fi] = _bestA[fi]._bestC;

    //  Clearing this destroys unitigs??

    if (bestCold[fi].isContained == false) {
      assert(_bestA[fi]._bestC.container       == 0);
      assert(_bestA[fi]._bestC.sameOrientation == false);
      assert(_bestA[fi]._bestC.a_hang          == 0);
      assert(_bestA[fi]._bestC.b_hang          == 0);
    }

    _bestA[fi]._bestC.container       = 0;
    _bestA[fi]._bestC.sameOrientation = false;
    _bestA[fi]._bestC.a_hang          = 0;
    _bestA[fi]._bestC.b_hang          = 0;
  }

  //  Allocate space for new scores

  assert(_scorA == NULL);

  _scorA = new BestScores [fiLimit + 1];

  memset(_scorA, 0, sizeof(BestScores) * (fiLimit + 1));

  //  Rebuild contains ignoring singleton containers

  for (uint32 fi=1; fi<=fiLimit; fi++) {
    uint32      no   = 0;

    if (bestCold[fi].isContained == false)
      continue;

    BAToverlap *ovl  = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    for (uint32 ii=0; ii<no; ii++) {
      uint32     autg = Unitig::fragIn(ovl[ii].a_iid);
      uint32     butg = Unitig::fragIn(ovl[ii].b_iid);

      assert(autg == 0);  //  Contained cannot be placed yet.

      if ((butg != 0) &&
          (unitigs[butg]->ufpath.size() == 1))
        //  Skip; container is a in a unitig, and that unitig is a singleton.
        continue;

      scoreContainment(ovl[ii]);
    }
  }

  delete [] _scorA;
  _scorA = NULL;

  //  Remove best edges for contains (shouldn't be any; we didn't make new ones since the last time we removed)

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    if (isContained(fi) == true) {
      getBestEdgeOverlap(fi, false)->set(0, 0, 0, 0);
      getBestEdgeOverlap(fi, true) ->set(0, 0, 0, 0);
    }
  }

  //  Log changes

  for (uint32 fi=0; fi<=fiLimit; fi++) {
    if ((bestCold[fi].container       != _bestA[fi]._bestC.container) ||
        (bestCold[fi].sameOrientation != _bestA[fi]._bestC.sameOrientation) ||
        (bestCold[fi].a_hang          != _bestA[fi]._bestC.a_hang) ||
        (bestCold[fi].b_hang          != _bestA[fi]._bestC.b_hang))
      writeLog("frag %u changed container from %c %u/%c/%d/%d to %c %u/%c/%d/%d\n",
               fi,
               bestCold[fi].isContained ? 'T' : 'F',
               bestCold[fi].container,
               bestCold[fi].sameOrientation ? 'N' : 'A',
               bestCold[fi].a_hang,
               bestCold[fi].b_hang,
               _bestA[fi]._bestC.isContained ? 'T' : 'F',
               _bestA[fi]._bestC.container,
               _bestA[fi]._bestC.sameOrientation ? 'N' : 'A',
               _bestA[fi]._bestC.a_hang,
               _bestA[fi]._bestC.b_hang);
  }

  delete [] bestCold;
}





BestOverlapGraph::BestOverlapGraph(double               erate,
                                   set<uint32>         *restrict) {

  _erate = erate;

  _bestA = NULL;
  _scorA = NULL;

  _bestM.clear();
  _scorM.clear();

  assert(restrict != NULL);

  _restrict        = restrict;
  _restrictEnabled = true;

  //  PASS 0:  Load the map (necessary?)

#if 0
  for (set<uint32>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    uint32      fi  = *it;

    _bestM[fi].insert();
    _scorM[fi].insert();
  }
#endif

  //  PASS 1:  Find containments.

  for (set<uint32>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    uint32      fi  = *it;
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    for (uint32 ii=0; ii<no; ii++)
      scoreContainment(ovl[ii]);
  }

  //  PASS 2:  Find dovetails.

  for (set<uint32>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    uint32      fi  = *it;
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    for (uint32 ii=0; ii<no; ii++)
      scoreEdge(ovl[ii]);
  }

  //  Remove temporary scoring data

  _scorM.clear();

  //  Remove dovetail overlaps for contained fragments.

  for (set<uint32>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    uint32      fi  = *it;

    if (isContained(fi) == true) {
      getBestEdgeOverlap(fi, false)->set(0, 0, 0, 0);
      getBestEdgeOverlap(fi, true) ->set(0, 0, 0, 0);
    }
  }

  //  Remove spurs

#if 0
  for (set<uint32>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    uint32      fi  = *it;

    if ((getBestEdgeOverlap(fi, false)->fragId() == 0) ||
        (getBestEdgeOverlap(fi, true)->fragId()  == 0)) {
      getBestEdgeOverlap(fi, false)->set(0, 0, 0, 0);
      getBestEdgeOverlap(fi, true) ->set(0, 0, 0, 0);
    }
  }
#endif

  _restrict        = NULL;
  _restrictEnabled = false;
}




void
BestOverlapGraph::reportBestEdges(const char *prefix) {
  char  N[FILENAME_MAX];

  sprintf(N, "%s.best.contains",   prefix);   FILE *BC = fopen(N, "w");
  sprintf(N, "%s.best.edges",      prefix);   FILE *BE = fopen(N, "w");
  sprintf(N, "%s.best.singletons", prefix);   FILE *BS = fopen(N, "w");

  if ((BC) && (BE) && (BS)) {
    fprintf(BC, "#fragId\tlibId\tbestCont\teRate\n");
    fprintf(BE, "#fragId\tlibId\tbest5iid\tbest5end\tbest3iid\tbest3end\teRate5\teRate3\n");
    fprintf(BS, "#fragId\tlibId\n");

    for (uint32 id=1; id<FI->numFragments() + 1; id++) {
      BestContainment *bestcont  = getBestContainer(id);
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      if (bestcont->isContained) {
        double  erate = OC->findErate(id, bestcont->container);

        fprintf(BC, "%u\t%u\t%u\t%6.4f\n", id, FI->libraryIID(id), bestcont->container, erate);
      }

      else if ((bestedge5->fragId() > 0) || (bestedge3->fragId() > 0)) {
        double  erate5 = OC->findErate(id, bestedge5->fragId());
        double  erate3 = OC->findErate(id, bestedge3->fragId());

        fprintf(BE, "%u\t%u\t%u\t%c'\t%u\t%c'\t%6.4f\t%6.4f\n", id, FI->libraryIID(id),
                bestedge5->fragId(), bestedge5->frag3p() ? '3' : '5',
                bestedge3->fragId(), bestedge3->frag3p() ? '3' : '5',
                erate5, erate3);
      }

      else {
        fprintf(BS, "%u\t%u\n", id, FI->libraryIID(id));
      }
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
      writeLog("scoreEdge()-- OVERLAP BADQ:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- bad quality\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
    return;
  }

  if (isOverlapRestricted(olap)) {
    //  Whoops, don't want this overlap for this BOG
    if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
      writeLog("scoreEdge()-- OVERLAP REST:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- restricted\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
    return;
  }

  if (isSuspicious(olap.b_iid)) {
    //  Whoops, don't want this overlap for this BOG
    if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
      writeLog("scoreEdge()-- OVERLAP SUSP:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- suspicious\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
    return;
  }

  if (((olap.a_hang >= 0) && (olap.b_hang <= 0)) ||
      ((olap.a_hang <= 0) && (olap.b_hang >= 0))) {
    //  Skip containment overlaps.
    if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
      writeLog("scoreEdge()-- OVERLAP CONT:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- container read\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
    return;
  }

  if (isContained(olap.b_iid) == true) {
    //  Skip overlaps to contained reads (allow scoring of best edges from contained reads).
    if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
      writeLog("scoreEdge()-- OVERLAP CONT:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- contained read\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
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
             olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
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

  //  The overlap is GOOD (false == not bad) if the error rate is below the allowed erate.
  //  Initially, this is just the erate passed in.  After the first rount of finding edges,
  //  it is reset to the mean and stddev of selected best edges.
  //
  if (olap.erate <= _mean + 5 * _stddev) {
    if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
      writeLog("isOverlapBadQuality()-- OVERLAP GOOD:     %d %d %c  hangs "F_S32" "F_S32" err %.3f\n",
               olap.a_iid, olap.b_iid,
               olap.flipped ? 'A' : 'N',
               olap.a_hang,
               olap.b_hang,
               olap.erate);

    return(false);
  }

  //  CA8.3 would further compute the (expected) number of errors in the alignment, and compare
  //  against another limit.  This was to allow very short overlaps where one error would push the
  //  error rate above a few percent.  canu doesn't do short overlaps.

  if ((enableLog == true) && (logFileFlags & LOG_OVERLAP_QUALITY))
    writeLog("isOverlapBadQuality()-- OVERLAP REJECTED: %d %d %c  hangs "F_S32" "F_S32" err %.3f\n",
             olap.a_iid, olap.b_iid,
             olap.flipped ? 'A' : 'N',
             olap.a_hang,
             olap.b_hang,
             olap.erate);

  return(true);
}



//  If no restrictions are known, this overlap is useful if both fragments are not in a unitig
//  already.  Otherwise, we are restricted to just a specific set of fragments (usually a single
//  unitig and all the mated reads).  The overlap is useful if both fragments are in the set.
//
bool
BestOverlapGraph::isOverlapRestricted(const BAToverlap &olap) {

  if (_restrictEnabled == false)
    return(false);

  assert(_restrict != NULL);

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
  uint64  corr = AS_MAX_EVALUE - olap.evalue;
  uint64  orig = AS_MAX_EVALUE - 0;

  //  Shift AFTER assigning to a 64-bit value to avoid overflows.
  corr <<= AS_MAX_EVALUE_BITS;

  //  Containments - the length of the overlaps are all the same.  We return the quality.
  //
  if (((olap.a_hang >= 0) && (olap.b_hang <= 0)) ||
      ((olap.a_hang <= 0) && (olap.b_hang >= 0)))
    return(corr | orig);

  //  Dovetails - the length of the overlap is the score, but we bias towards lower error.
  //  (again, shift AFTER assigning to avoid overflows)

  assert((olap.a_hang < 0) == (olap.b_hang < 0));

  //  Compute the length of the overlap, as either the official overlap length that currently
  //  takes into account both reads, or as the number of aligned bases on the A read.

#if 0
  leng   = FI->overlapLength(olap.a_iid, olap.b_iid, olap.a_hang, olap.b_hang);
#endif

  if (olap.a_hang > 0)
    leng = FI->fragmentLength(olap.a_iid) - olap.a_hang;
  else
    leng = FI->fragmentLength(olap.a_iid) + olap.b_hang;

  //  Convert the length into an expected number of matches.

#if 0
  assert(olap.erate <= 1.0);
  leng  -= leng * olap.erate;
#endif

  //  And finally shift it to the correct place in the word.

  leng <<= (2 * AS_MAX_EVALUE_BITS);

  return(leng | corr | orig);
}

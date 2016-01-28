
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Unitig.H"

#include "intervalList.H"


//  HACK
uint32  examineOnly = UINT32_MAX;


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
BestOverlapGraph::examineOnlyTopN(void) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  writeLog("BestOverlapGraph()-- analyzing %d fragments for best edges, with %d threads.\n", fiLimit, numThreads);
  writeLog("BestOverlapGraph()-- scoring highest quality %d overlaps.\n", examineOnly);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    uint32      n5  = 0;
    uint32      n3  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    sort(ovl, ovl + no, BAToverlap_sortByErate);

    for (uint32 ii=0; ii<no; ii++) {
      if (((ovl[ii].a_hang >= 0) && (ovl[ii].b_hang <= 0)) ||
          ((ovl[ii].a_hang <= 0) && (ovl[ii].b_hang >= 0)))
        //  Don't do contains here!
        continue;

      //  Process the 5' overlaps.
      if ((n5 < examineOnly) &&
          (ovl[ii].a_hang < 0)) {
        assert(ovl[ii].b_hang < 0);
        n5++;
        scoreEdge(ovl[ii]);
      }

      //  Process the 3' overlaps.
      if ((n3 < examineOnly) &&
          (ovl[ii].a_hang > 0)) {
        assert(ovl[ii].b_hang > 0);
        n3++;
        scoreEdge(ovl[ii]);
      }
    }
  }
}


void
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
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    for (uint32 ii=0; ii<no; ii++)
      if (isSpur[ovl[ii].b_iid] == false)
        scoreContainment(ovl[ii]);
  }

  //  PASS 4:  Find dovetails.

  writeLog("BestOverlapGraph()-- analyzing %d fragments for best edges, with %d threads.\n", fiLimit, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    for (uint32 ii=0; ii<no; ii++)
      if (isSpur[ovl[ii].b_iid] == false)
        scoreEdge(ovl[ii]);
  }

  delete [] isSpur;
}


void
BestOverlapGraph::removeFalseBest(void) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  writeLog("BestOverlapGraph()-- detecting false best overlaps.\n");

  //  WARNING!  This code was quite confused about erate and evalue (back when they were called
  //  error and errorValue or somethihng confusing like that).  Use with caution.

  uint32    *histo5 = new uint32 [AS_MAX_EVALUE + 1];
  uint32    *histo3 = new uint32 [AS_MAX_EVALUE + 1];

  memset(histo5, 0, sizeof(uint32) * (AS_MAX_EVALUE + 1));
  memset(histo3, 0, sizeof(uint32) * (AS_MAX_EVALUE + 1));

  uint32    *erate5 = new uint32 [fiLimit + 1];
  uint32    *erate3 = new uint32 [fiLimit + 1];

  memset(erate5,     0, sizeof(uint32) * (fiLimit + 1));
  memset(erate3,     0, sizeof(uint32) * (fiLimit + 1));

  char   *altBest = new char [fiLimit + 1];
  char   *isBad   = new char [fiLimit + 1];

  memset(altBest, 0, sizeof(char) * (fiLimit + 1));
  memset(isBad,   0, sizeof(char) * (fiLimit + 1));

  //  Compute a histogram of the current best edges, and save the erate of the best for each read.

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32            olapsLen = 0;
    BAToverlap       *olaps    = OC->getOverlaps(fi, AS_MAX_EVALUE, olapsLen);

    BestEdgeOverlap  *ovl5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap  *ovl3 = getBestEdgeOverlap(fi, true);

    for (uint32 oo=0; oo<olapsLen; oo++) {
      assert(fi == olaps[oo].a_iid);

      if (ovl5->fragId() == olaps[oo].b_iid) {
        histo5[olaps[oo].evalue]++;
        erate5[fi] = olaps[oo].erate;
      }

      if (ovl3->fragId() == olaps[oo].b_iid) {
        histo3[olaps[oo].evalue]++;
        erate3[fi] = olaps[oo].erate;
      }
    }
  }

  //  Compute a nice threshold.  Find the mean and stddev of the best edge error rates.

  double  m5 = 0, s5 = 100;
  double  m3 = 0, s3 = 100;

  for (uint32 xx=0; xx<10; xx++) {
    double  mean5   =   0, mean3   =   0;
    uint64  count5  =   0, count3  =   0;
    double  stddev5 = 100, stddev3 = 100;

    for (uint32 er=0; er <= AS_MAX_EVALUE; er++) {
      double ER = AS_OVS_decodeEvalue(er) * 100;

      if (ER <= 0.0)
        continue;

      if ((m5 - 3 * s5 <= ER) &
          (ER          <= m5 + 3 * s5)) {
        mean5  += histo5[er] * ER;
        count5 += histo5[er];
      }

      if ((m3 - 3 * s3 <= ER) &
          (ER          <= m3 + 3 * s3)) {
        mean3  += histo3[er] * ER;
        count3 += histo3[er];
      }
    }

    mean5 /= count5;
    mean3 /= count3;

    for (uint32 er=0; er <= AS_MAX_EVALUE; er++) {
      double ER =  AS_OVS_decodeEvalue(er) * 100;

      if (ER <= 0.0)
        continue;

      if ((m5 - 3 * s5 <= ER) &
          (ER          <= m5 + 3 * s5)) {
        stddev5 += histo5[er] * (ER - mean5) * (ER - mean5);
      }

      if ((m3 - 3 * s3 <= ER) &
          (ER          <= m3 + 3 * s3)) {
        stddev3 += histo3[er] * (ER - mean3) * (ER - mean3);
      }
    }

    stddev5 = sqrt(stddev5 / (count5 - 1));
    stddev3 = sqrt(stddev3 / (count3 - 1));

    m5 = mean5;    m3 = mean3;
    s5 = stddev5;  s3 = stddev3;

    fprintf(stderr, "mean  %.4f +- %.4f --- %.4f +- %.4f     ",
            mean5, stddev5,
            mean3, stddev3);

    fprintf(stderr, "set xtics ( %.4f, %.4f,    %.4f, %.4f)\n",
            mean5 - 3 * stddev5, mean5 + 3 * stddev5,
            mean3 - 3 * stddev3, mean3 + 3 * stddev3);
  }  //  xx 10 times to stabilize


  //  Output the evalue histogram

  {
    char  EN[FILENAME_MAX];

    sprintf(EN, "best.edges.erate.histogram");

    errno = 0;
    FILE *EH = fopen(EN, "w");
    if (errno)
      fprintf(stderr, "BestOverlapGraph()-- failed to open '%s' for writing: %s\n", EN, strerror(errno)), exit(1);

    for (uint32 er=0; er <= AS_MAX_EVALUE; er++) {
      double  ER = AS_OVS_decodeEvalue(er);

      if (ER <= 0.0)
        continue;

      fprintf(EH, "%.4f\t%u\t%u\t%f\t%f\n",
              ER,
              histo5[er],  //  HUH?!?!  This used to output decodedEvalue() of the histogram...nonsense!
              histo3[er],
              fabs(m5 - ER) / s5,   //  Grubb's test
              fabs(m3 - ER) / s3);
    }

    fclose(EH);
  }

  //  For any read with best edge above THRESHOLD error, see if there is an alternate best
  //  that is within the target error range.

  double  erate5thresh = m5 + 2 * s5;  //  Discard best if it is worse than 2 s.d. from mean.
  double  erate3thresh = m3 + 2 * s3;

  for (uint32 fi=1; fi <= fiLimit; fi++) {

    if (erate5[fi] > erate5thresh) {
      fprintf(stderr, "RECOMPUTE frag %u 5'\n", fi);
      isBad[fi] = true;
    }

    if (erate3[fi] > erate3thresh) {
      fprintf(stderr, "RECOMPUTE frag %u 3'\n", fi);
      isBad[fi] = true;
    }
  }

  fprintf(stderr, "thresholds %f %f\n", erate5thresh, erate3thresh);

  double erateCthresh = MIN(erate5thresh, erate3thresh);

#if 0
  //  Remove best edges, so we can rebuild

  memset(_bestA, 0, sizeof(BestOverlaps) * (fiLimit + 1));
  memset(_scorA, 0, sizeof(BestScores)   * (fiLimit + 1));

  //  Rebuild best edges, ignoring edges to spurs.  We build edges out of spurs, but don't allow edges into them.
  //  This should prevent them from being incorporated into a promiscuous unitig, but still let them be popped
  //  as bubbles (but they shouldn't because they're spurs).

  //  PASS 3:  Find containments.

  writeLog("BestOverlapGraph()-- analyzing %d fragments for best contains, with %d threads.\n", fiLimit, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    for (uint32 ii=0; ii<no; ii++)
      assert(ovl[ii].a_iid == fi);

    for (uint32 ii=0; ii<no; ii++)
      if (ovl[ii].error < erateCthresh)
        scoreContainment(ovl[ii]);
  }

  //  PASS 4:  Find dovetails.

  writeLog("BestOverlapGraph()-- analyzing %d fragments for best edges, with %d threads.\n", fiLimit, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    if (isBad[ovl[fi].a_iid] == true)
      continue;

    for (uint32 ii=0; ii<no; ii++)
      assert(ovl[ii].a_iid == fi);

    for (uint32 ii=0; ii<no; ii++)
      if (isBad[ovl[ii].b_iid] == false)
        scoreEdge(ovl[ii]);
  }
#endif

  delete [] histo5;
  delete [] histo3;

  delete [] erate5;
  delete [] erate3;
}










void
BestOverlapGraph::removeWeak(double threshold) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  writeLog("BestOverlapGraph()-- detecting weak overlaps.\n");

  //  For each read, mark an overlap as bad if it falls in the lower
  //  X% of overlaps sorted by identity.

  uint32  *minEvalue5p = new uint32 [fiLimit + 1];
  uint32  *minEvalue3p = new uint32 [fiLimit + 1];

  memset(minEvalue5p, 0, sizeof(uint32) * (fiLimit + 1));
  memset(minEvalue3p, 0, sizeof(uint32) * (fiLimit + 1));

  uint32   evaluesMax  = 1048576;
  uint32   evalues5len = 0;
  uint32  *evalues5    = new uint32 [evaluesMax];
  uint32   evalues3len = 0;
  uint32  *evalues3    = new uint32 [evaluesMax];

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32            olapsLen = 0;
    BAToverlap       *olaps    = OC->getOverlaps(fi, AS_MAX_EVALUE, olapsLen);

    uint64            ovl5sum = 0;
    uint32            ovl5cnt = 0;

    uint64            ovl3sum = 0;
    uint32            ovl3cnt = 0;

    evalues5len = 0;
    evalues5[0] = 0;

    evalues3len = 0;
    evalues3[0] = 0;

    //  Find the error rate histogram for each end.

    for (uint32 oo=0; oo<olapsLen; oo++) {
      assert(fi == olaps[oo].a_iid);

      if ((AS_BAT_overlapAEndIs5prime(olaps[oo])) && (evalues5len < evaluesMax))
        evalues5[evalues5len++] = olaps[oo].evalue;

      if ((AS_BAT_overlapAEndIs3prime(olaps[oo])) && (evalues5len < evaluesMax))
        evalues3[evalues3len++] = olaps[oo].evalue;
    }

    //  Sort by increasing error rate.

    sort(evalues5, evalues5 + evalues5len);
    sort(evalues3, evalues3 + evalues3len);

    //  Pick a min erate for each end.

    minEvalue5p[fi] = evalues5[(int32)(evalues5len - evalues5len * threshold)];
    minEvalue3p[fi] = evalues3[(int32)(evalues3len - evalues3len * threshold)];

    if ((fi % 1000) == 0) {
      fprintf(stderr, "len %d %d t %f vals %d %f %d %f\n", evalues5len, evalues3len, threshold,
              minEvalue5p[fi], AS_OVS_decodeEvalue(minEvalue5p[fi]),
              minEvalue3p[fi], AS_OVS_decodeEvalue(minEvalue3p[fi]));
    }
  }

  delete [] evalues5;
  delete [] evalues3;

  //  Throw this at the OverlapCache, so it can remove overlaps.

  OC->removeWeakOverlaps(minEvalue5p, minEvalue3p);

  delete [] minEvalue5p;
  delete [] minEvalue3p;
}



BestOverlapGraph::BestOverlapGraph(double               erate,
                                   const char          *prefix,
                                   double               doRemoveWeakThreshold,
                                   bool                 doRemoveSuspicious,
                                   bool                 doRemoveSpurs) {

  bool doExamineOnlyTopN = false;
  bool doRemoveFalseBest = false;

  setLogFile(prefix, "bestOverlapGraph");

  writeLog("BestOverlapGraph-- allocating best edges ("F_SIZE_T"MB) and containments ("F_SIZE_T"MB)\n",
           ((2 * sizeof(BestEdgeOverlap) * (FI->numFragments() + 1)) >> 20),
           ((1 * sizeof(BestContainment) * (FI->numFragments() + 1)) >> 20));

  _bestA = new BestOverlaps [FI->numFragments() + 1];
  _scorA = new BestScores   [FI->numFragments() + 1];

  memset(_bestA, 0, sizeof(BestOverlaps) * (FI->numFragments() + 1));
  memset(_scorA, 0, sizeof(BestScores)   * (FI->numFragments() + 1));

  _restrict        = NULL;
  _restrictEnabled = false;

  _erate = erate;

  //  Initialize parallelism.

  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  //  PASS 0:  Find suspicious fragments.  For any found, mark as suspicious and don't allow
  //  these to be best overlaps.

  if (doRemoveWeakThreshold > 0.0)
    removeWeak(doRemoveWeakThreshold);

  if (doRemoveSuspicious)
    removeSuspicious();

  //  PASS 1:  Find containments.

  writeLog("BestOverlapGraph()-- analyzing %d fragments for best contains, with %d threads.\n", fiLimit, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    for (uint32 ii=0; ii<no; ii++)
      scoreContainment(ovl[ii]);
  }

  //  PASS 2:  Find dovetails.

  writeLog("BestOverlapGraph()-- analyzing %d fragments for best edges, with %d threads.\n", fiLimit, numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

    for (uint32 ii=0; ii<no; ii++)
      scoreEdge(ovl[ii]);
  }

  //  Now, several optional refinements.

  if (doExamineOnlyTopN)
    examineOnlyTopN();

  if (doRemoveSpurs)
    removeSpurs();

  if (doRemoveFalseBest)
    removeFalseBest();

  //  Done with the scoring data.

  delete [] _scorA;
  _scorA = NULL;

  //  Finally, remove dovetail overlaps for contained fragments.

  writeLog("BestOverlapGraph()-- removing best edges for contained fragments, with %d threads.\n", numThreads);

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    if (isContained(fi) == true) {
      getBestEdgeOverlap(fi, false)->set(0, 0, 0, 0);
      getBestEdgeOverlap(fi, true) ->set(0, 0, 0, 0);
    }
  }

  writeLog("BestOverlapGraph()-- dumping best edges/contains/singletons.\n");

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

#if 0
    writeLog("set best for %d from "F_U32" score="F_U64" to "F_U32" score="F_U64"\n",
             olap.a_iid,
             c->container, bestCscore(olap.a_iid),
             olap.b_iid,   newScr);
#endif

    c->container         = olap.b_iid;
    c->isContained       = true;
    c->sameOrientation   = olap.flipped ? false : true;
    c->a_hang            = olap.flipped ? olap.b_hang : -olap.a_hang;
    c->b_hang            = olap.flipped ? olap.a_hang : -olap.b_hang;

    bestCscore(olap.a_iid) = newScr;
#if 0
  } else {
    writeLog("NOT best for %d WITH "F_U32" score="F_U64"\n",
             olap.a_iid,
             olap.b_iid,   newScr);
#endif
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

  if ((isContained(olap.a_iid) == true) ||
      (isContained(olap.b_iid) == true)) {
    //  Skip contained fragments.
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

  //  The overlap is GOOD (false == not bad) if the corrected error rate is below the requested
  //  erate.
  //
  if (olap.erate <= _erate) {
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

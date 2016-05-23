
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
 *    Sergey Koren beginning on 2016-MAR-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "intervalList.H"
#include "stddev.H"



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
#pragma omp critical (suspInsert)
      _suspicious.insert(fi);
    }
  }

  writeLog("BestOverlapGraph()-- marked "F_U64" reads as suspicious.\n", _suspicious.size());
}



void
BestOverlapGraph::removeHighErrorBestEdges(void) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  stdDev<double>  edgeStats;

  //  Find the overlap for every best edge.

  double  *erates    = new double [fiLimit + 1 + fiLimit + 1];
  double  *absdev    = new double [fiLimit + 1 + fiLimit + 1];
  uint32   eratesLen = 0;

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap *b5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap *b3 = getBestEdgeOverlap(fi, true);

    if (b5->fragId() != 0)   edgeStats.insert(erates[eratesLen++] = b5->erate());
    if (b3->fragId() != 0)   edgeStats.insert(erates[eratesLen++] = b3->erate());
  }

  _mean   = edgeStats.mean();
  _stddev = edgeStats.stddev();

  writeLog("removeHighErrorBestEdges()-- with %u points - mean %f stddev %f -- would use overlaps below %f fraction error\n",
           edgeStats.size(), _mean, _stddev, _mean + _deviationGraph * _stddev);

  //  Find the median and absolute deviations.

  sort(erates, erates+eratesLen);

  _median = erates[ eratesLen / 2 ];

  for (uint32 ii=0; ii<eratesLen/2; ii++)
    absdev[ii] = _median - erates[ii];

  for (uint32 ii=eratesLen/2; ii<eratesLen; ii++)
    absdev[ii] = erates[ii] - _median;

  sort(absdev, absdev+eratesLen);

  assert(absdev[0] >= 0.0);

  _mad    = absdev[eratesLen/2];

  delete [] absdev;
  delete [] erates;

  writeLog("removeHighErrorBestEdges()-- with %u points - median %f mad %f - would use overlaps below %f fraction error\n",
           edgeStats.size(), _median, _mad, _median + _deviationGraph * 1.4826 * _mad);

  //  The real filtering is done on the next pass through findEdges().  Here, we just report statistics.

  uint32  noedge   = 0;
  uint32  removed  = 0;
  uint32  retained = 0;

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap *b5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap *b3 = getBestEdgeOverlap(fi, true);

    if      (b5->fragId() == 0)
      noedge++;
    else if (b5->erate() > _mean + _deviationGraph * _stddev)
      removed++;
    else
      retained++;

    if      (b3->fragId() == 0)
      noedge++;
    else if (b3->erate() > _mean + _deviationGraph * _stddev)
      removed++;
    else
      retained++;
  }

  writeLog("removeHighErrorBestEdges()-- %u ends have no best edge; %u ends are suspiciously high error; %u ends are acceptable.\n",
           noedge, removed, retained);
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
#pragma omp atomic
      nSuspicious++;
      continue;
    }

    if (isContained(fi) == true) {
#pragma omp atomic
      nContained++;
      continue;
    }

    if ((this5->fragId() == 0) ||
        (this3->fragId() == 0)) {
#pragma omp atomic
      nSpur++;
      continue;
    }

    //  Find the overlap for this5 and this3.

    int32   this5ovlLen = FI->overlapLength(fi, this5->fragId(), this5->ahang(), this5->bhang());
    int32   this3ovlLen = FI->overlapLength(fi, this3->fragId(), this3->ahang(), this3->bhang());

    //  Find the edges for our best overlaps.

    BestEdgeOverlap *that5 = getBestEdgeOverlap(this5->fragId(), this5->frag3p());
    BestEdgeOverlap *that3 = getBestEdgeOverlap(this3->fragId(), this3->frag3p());

    //  If both point back to us, we're done.

    if ((that5->fragId() == fi) && (that5->frag3p() == false) &&
        (that3->fragId() == fi) && (that3->frag3p() == true)) {
#pragma omp atomic
      nMutual++;
      continue;
    }

    //  If there is an overlap to something with no overlaps out of it, that's
    //  a little suspicious.

    if ((that5->fragId() == 0) ||
        (that3->fragId() == 0)) {
      writeLog("WARNING: read %u has overlap to spur - 3' to read %u back to %u - 5' to read %u back to %u\n",
               fi,
               this5->fragId(), that5->fragId(),
               this3->fragId(), that3->fragId());
#pragma omp critical (suspInsert)
      _suspicious.insert(fi);
      continue;
    }

    //  Something doesn't agree.  Find those overlaps...

    int32  that5ovlLen = FI->overlapLength(this5->fragId(), that5->fragId(), that5->ahang(), that5->bhang());
    int32  that3ovlLen = FI->overlapLength(this3->fragId(), that3->fragId(), that3->ahang(), that3->bhang());

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



void
BestOverlapGraph::removeSpurs(void) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  writeLog("BestOverlapGraph()-- detecting spur fragments.\n");

  _spur.clear();

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    bool   spur5 = (getBestEdgeOverlap(fi, false)->fragId() == 0);
    bool   spur3 = (getBestEdgeOverlap(fi, true)->fragId()  == 0);

    if (isContained(fi))
      //  Contained, not a spur.
      continue;

    if ((spur5 == false) && (spur3 == false))
      //  Edges off of both ends.  Not a spur.
      continue;

    if ((spur5 == true)  && (spur3 == true))
      //  No edges off either end.  Not a spur, just garbage.
      continue;

    //  Exactly one end is missing a best edge.  Bad!

    writeLog("BestOverlapGraph()-- frag "F_U32" is a %s spur.\n", fi, (spur5) ? "5'" : "3'");

    _spur.insert(fi);
  }
}



void
BestOverlapGraph::findEdges(void) {
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
      if (_spur.count(ovl[ii].b_iid) == 0)
        scoreEdge(ovl[ii]);
  }
}



void
BestOverlapGraph::removeContainedDovetails(void) {
  uint32  fiLimit    = FI->numFragments();

  writeLog("BestOverlapGraph()-- removing best edges for contained fragments.\n");

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    if (isContained(fi) == true) {
      getBestEdgeOverlap(fi, false)->clear();
      getBestEdgeOverlap(fi, true) ->clear();
    }
  }
}



BestOverlapGraph::BestOverlapGraph(double        erateGraph,
                                   double        deviationGraph,
                                   const char   *prefix) {

  setLogFile(prefix, "bestOverlapGraph");

  writeLog("BestOverlapGraph-- allocating best edges ("F_SIZE_T"MB)\n",
           ((2 * sizeof(BestEdgeOverlap) * (FI->numFragments() + 1)) >> 20));

  _bestA           = new BestOverlaps [FI->numFragments() + 1];  //  Cleared in findEdges()
  _scorA           = new BestScores   [FI->numFragments() + 1];

  _mean            = erateGraph;
  _stddev          = 0.0;

  _median          = erateGraph;
  _mad             = 0.0;

  _suspicious.clear();

  _bestM.clear();
  _scorM.clear();

  _restrict        = NULL;
  _restrictEnabled = false;

  _erateGraph      = erateGraph;
  _deviationGraph  = deviationGraph;

  //  Mark reads as suspicious if they are not fully covered by overlaps.

  removeSuspicious();
  findEdges();

  if (logFileFlagSet(LOG_ALL_BEST_EDGES))
    reportBestEdges(prefix, "best.0.initial");

  //  Analyze the current best edges to set a cutoff on overlap quality used for graph building.

  removeHighErrorBestEdges();
  findEdges();

  if (logFileFlagSet(LOG_ALL_BEST_EDGES))
    reportBestEdges(prefix, "best.1.filtered");

  //  Mark reads as suspicious if the length of the best edge out is very different than the length
  //  of the best edge that should be back to us.  E.g., if readA has best edge to readB (of length
  //  lenAB), but readB has best edge to readC (of length lenBC), and lenAB is much shorter than
  //  lenBC, then something is wrong with readA.
  //
  //  This must come before removeSpurs().

  removeLopsidedEdges();
  findEdges();

  if (logFileFlagSet(LOG_ALL_BEST_EDGES))
    reportBestEdges(prefix, "best.2.cleaned");

  //  Mark reads as spurs, so we don't find best edges to them.

  removeSpurs();
  findEdges();

  reportBestEdges(prefix, logFileFlagSet(LOG_ALL_BEST_EDGES) ? "best.3.final" : "best");

  //  One more pass, to find any ambiguous best edges.

  //  Cleanup the contained reads.  Why?

  removeContainedDovetails();

  //  Done with scoring data.

  delete [] _scorA;
  _scorA = NULL;

  _spur.clear();

  setLogFile(prefix, NULL);
}




BestOverlapGraph::BestOverlapGraph(double         erateGraph,
                                   double         deviationGraph,
                                   set<uint32>   *restrict) {

  _erateGraph     = erateGraph;
  _deviationGraph = deviationGraph;

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
      getBestEdgeOverlap(fi, false)->clear();
      getBestEdgeOverlap(fi, true) ->clear();
    }
  }

  //  Remove spurs

#if 0
  for (set<uint32>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    uint32      fi  = *it;

    if ((getBestEdgeOverlap(fi, false)->fragId() == 0) ||
        (getBestEdgeOverlap(fi, true)->fragId()  == 0)) {
      getBestEdgeOverlap(fi, false)->clear();
      getBestEdgeOverlap(fi, true) ->clear();
    }
  }
#endif

  _restrict        = NULL;
  _restrictEnabled = false;
}




void
BestOverlapGraph::reportBestEdges(const char *prefix, const char *label) {
  char  N[FILENAME_MAX];
  FILE             *BCH = NULL;
  FILE *BE = NULL, *BEH = NULL, *BEG;
  FILE *BS = NULL;
  FILE *SS = NULL;

  sprintf(N, "%s.%s.edges",               prefix, label);   BE = fopen(N, "w");
  sprintf(N, "%s.%s.singletons",          prefix, label);   BS = fopen(N, "w");
  sprintf(N, "%s.%s.edges.suspicious",    prefix, label);   SS = fopen(N, "w");

  sprintf(N, "%s.%s.contains.histogram",  prefix, label);   BCH = fopen(N, "w");
  sprintf(N, "%s.%s.edges.histogram",     prefix, label);   BEH = fopen(N, "w");

  sprintf(N, "%s.%s.edges.gfa",           prefix, label);   BEG = fopen(N, "w");

  if ((BE) && (BS)) {
    fprintf(BE, "#fragId\tlibId\tbest5iid\tbest5end\tbest3iid\tbest3end\teRate5\teRate3\tbest5len\tbest3len\n");
    fprintf(BS, "#fragId\tlibId\n");

    for (uint32 id=1; id<FI->numFragments() + 1; id++) {
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      if ((bestedge5->fragId() == 0) && (bestedge3->fragId() == 0) && (isContained(id) == false)) {
        fprintf(BS, "%u\t%u\n", id, FI->libraryIID(id));
      }

      else if (_suspicious.count(id) > 0) {
        fprintf(SS, "%u\t%u\t%u\t%c'\t%u\t%c'\t%6.4f\t%6.4f\t%u\t%u%s\n", id, FI->libraryIID(id),
          bestedge5->fragId(), bestedge5->frag3p() ? '3' : '5',
                bestedge3->fragId(), bestedge3->frag3p() ? '3' : '5',
                AS_OVS_decodeEvalue(bestedge5->evalue()),
                AS_OVS_decodeEvalue(bestedge3->evalue()),
                (bestedge5->fragId() == 0 ? 0 : FI->overlapLength(id, bestedge5->fragId(), bestedge5->ahang(), bestedge5->bhang())),
                (bestedge3->fragId() == 0 ? 0 : FI->overlapLength(id, bestedge3->fragId(), bestedge3->ahang(), bestedge3->bhang())),
                isContained(id) ? "\tcontained" : "");
      }

      else {
        fprintf(BE, "%u\t%u\t%u\t%c'\t%u\t%c'\t%6.4f\t%6.4f\t%u\t%u%s\n", id, FI->libraryIID(id),
                bestedge5->fragId(), bestedge5->frag3p() ? '3' : '5',
                bestedge3->fragId(), bestedge3->frag3p() ? '3' : '5',
                AS_OVS_decodeEvalue(bestedge5->evalue()),
                AS_OVS_decodeEvalue(bestedge3->evalue()),
                (bestedge5->fragId() == 0 ? 0 : FI->overlapLength(id, bestedge5->fragId(), bestedge5->ahang(), bestedge5->bhang())),
                (bestedge3->fragId() == 0 ? 0 : FI->overlapLength(id, bestedge3->fragId(), bestedge3->ahang(), bestedge3->bhang())),
                isContained(id) ? "\tcontained" : "");
      }
    }

    fclose(BE);
    fclose(BS);
    fclose(SS);
  }


  if (BEG) {
    fprintf(BEG, "H\tVN:Z:bogart/edges\n");

    //  First, write the sequences used.

    for (uint32 id=1; id<FI->numFragments() + 1; id++) {
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      if ((bestedge5->fragId() == 0) && (bestedge3->fragId() == 0) && (isContained(id) == false)) {
        //  Do nothing, a singleton.
      }

      else if (isContained(id) == true) {
        //  Do nothing, a contained read.
      }

      else if (_suspicious.count(id) > 0) {
        //  Do nothing, a suspicious read.
      }

      else {
        //  Report the read, it has best edges - including contained reads.
        fprintf(BEG, "S\tread%08u\t*\tLN:i:%u\n", id, FI->fragmentLength(id));
      }
    }

    //  Now, report edges.  GFA wants edges in exactly this format:
    //
    //       -------------
    //             -------------
    //
    //  with read orientation given by +/-.  Conveniently, this is what we've saved (for the edges).

    for (uint32 id=1; id<FI->numFragments() + 1; id++) {
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      if ((bestedge5->fragId() == 0) && (bestedge3->fragId() == 0) && (isContained(id) == false)) {
        //  Do nothing, a singleton.
      }

      else if (isContained(id) == true) {
        //  Do nothing, a contained read.
      }

      else if (_suspicious.count(id) > 0) {
        //  Do nothing, a suspicious read.
      }

      else {
        if (bestedge5->fragId() != 0) {
          int32  ahang   = bestedge5->ahang();
          int32  bhang   = bestedge5->bhang();
          int32  olaplen = FI->overlapLength(id, bestedge5->fragId(), bestedge5->ahang(), bestedge5->bhang());

          assert((ahang <= 0) && (bhang <= 0));  //  ALL 5' edges should be this.

          fprintf(BEG, "L\tread%08u\t-\tread%08u\t%c\t%uM\n",
                  id,
                  bestedge5->fragId(), bestedge5->frag3p() ? '-' : '+',
                  olaplen);
        }

        if (bestedge3->fragId() != 0) {
          int32  ahang   = bestedge3->ahang();
          int32  bhang   = bestedge3->bhang();
          int32  olaplen = FI->overlapLength(id, bestedge3->fragId(), bestedge3->ahang(), bestedge3->bhang());

          assert((ahang >= 0) && (bhang >= 0));  //  ALL 3' edges should be this.

          fprintf(BEG, "L\tread%08u\t+\tread%08u\t%c\t%uM\n",
                  id,
                  bestedge3->fragId(), bestedge3->frag3p() ? '-' : '+',
                  FI->overlapLength(id, bestedge3->fragId(), bestedge3->ahang(), bestedge3->bhang()));
        }
      }
    }

    fclose(BEG);
  }


  if ((BCH) && (BEH)) {
    double *bc = new double [FI->numFragments() + 1 + FI->numFragments() + 1];
    double *be = new double [FI->numFragments() + 1 + FI->numFragments() + 1];

    uint32  bcl = 0;
    uint32  bel = 0;

    for (uint32 id=1; id<FI->numFragments() + 1; id++) {
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      if (isContained(id)) {
        //bc[bcl++] = bestcont->erate();
#warning what is the error rate of the 'best contained' overlap?
        bc[bcl++] = bestedge5->erate();
        bc[bcl++] = bestedge3->erate();
      }
      else {
        if (bestedge5->fragId() > 0)
          be[bel++] = bestedge5->erate();

        if (bestedge3->fragId() > 0)
          be[bel++] = bestedge3->erate();
      }
    }

    sort(bc, bc+bcl);
    sort(be, be+bel);

    for (uint32 ii=0; ii<bcl; ii++)
      fprintf(BCH, "%f\n", bc[ii]);

    for (uint32 ii=0; ii<bel; ii++)
      fprintf(BEH, "%f\n", be[ii]);

    fclose(BCH);
    fclose(BEH);

    delete [] bc;
    delete [] be;
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

  setContained(olap.a_iid);
}



void
BestOverlapGraph::scoreEdge(const BAToverlap& olap) {
  bool   enableLog = false;  //  useful for reporting this stuff only for specific reads

  //if ((olap.a_iid == 97202) || (olap.a_iid == 30701))
  //  enableLog = true;

  if (isOverlapBadQuality(olap)) {
    //  Yuck.  Don't want to use this crud.
    if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
      writeLog("scoreEdge()-- OVERLAP BADQ:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- bad quality\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
    return;
  }

  if (isOverlapRestricted(olap)) {
    //  Whoops, don't want this overlap for this BOG
    if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
      writeLog("scoreEdge()-- OVERLAP RESTRICT: %d %d %c  hangs "F_S32" "F_S32" err %.3f -- restricted\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
    return;
  }

  if (isSuspicious(olap.b_iid)) {
    //  Whoops, don't want this overlap for this BOG
    if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
      writeLog("scoreEdge()-- OVERLAP SUSP:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- suspicious\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
    return;
  }

  if (((olap.a_hang >= 0) && (olap.b_hang <= 0)) ||
      ((olap.a_hang <= 0) && (olap.b_hang >= 0))) {
    //  Skip containment overlaps.
    if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
      writeLog("scoreEdge()-- OVERLAP CONT:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- container read\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
    return;
  }

  if (isContained(olap.b_iid) == true) {
    //  Skip overlaps to contained reads (allow scoring of best edges from contained reads).
    if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
      writeLog("scoreEdge()-- OVERLAP CONT:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- contained read\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
    return;
  }

  uint64           newScr = scoreOverlap(olap);
  bool             a3p    = olap.AEndIs3prime();
  BestEdgeOverlap *best   = getBestEdgeOverlap(olap.a_iid, a3p);
  uint64          &score  = (a3p) ? (best3score(olap.a_iid)) : (best5score(olap.a_iid));

  assert(newScr > 0);

  if (newScr <= score) {
    if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
      writeLog("scoreEdge()-- OVERLAP GOOD:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- no better than best\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate);
    return;
  }

  best->set(olap);

  score = newScr;

  if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
    writeLog("scoreEdge()-- OVERLAP BEST:     %d %d %c  hangs "F_S32" "F_S32" err %.3f -- NOW BEST\n",
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

  double Tstddev = _mean   + _deviationGraph * _stddev;
  double Tmad    = _median + _deviationGraph * 1.4826 * _mad;

  if (olap.erate <= Tmad || (Tmad == 0 && olap.erate <= Tstddev)) {
    if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
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

  if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
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
  uint64  leng = 0;
  uint64  rate = AS_MAX_EVALUE - olap.evalue;

  assert(olap.evalue <= AS_MAX_EVALUE);
  assert(rate <= AS_MAX_EVALUE);

  //  Containments - the length of the overlaps are all the same.  We return the quality.
  //
  if (((olap.a_hang >= 0) && (olap.b_hang <= 0)) ||
      ((olap.a_hang <= 0) && (olap.b_hang >= 0)))
    return(rate);

  //  Dovetails - the length of the overlap is the score, but we bias towards lower error.
  //  (again, shift AFTER assigning to avoid overflows)

  assert((olap.a_hang < 0) == (olap.b_hang < 0));

  //  Compute the length of the overlap, as either the official overlap length that currently
  //  takes into account both reads, or as the number of aligned bases on the A read.

#if 0
  leng = FI->overlapLength(olap.a_iid, olap.b_iid, olap.a_hang, olap.b_hang);
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

  leng <<= AS_MAX_EVALUE_BITS;

  return(leng | rate);
}

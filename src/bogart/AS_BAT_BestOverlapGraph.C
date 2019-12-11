
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

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "intervalList.H"
#include "stddev.H"



void
BestOverlapGraph::removeReadsWithCoverageGap(const char *UNUSED(prefix)) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  //  Search for reads that have an internal region with no coverage.
  //  If found, flag these as _coverageGap.

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32               no   = 0;
    BAToverlap          *ovl  = OC->getOverlaps(fi, no);

    uint32               fLen = RI->readLength(fi);

    bool                 verified = false;
    intervalList<int32>  IL;

    if (isIgnored(fi) == true)
      continue;

    for (uint32 ii=0; (ii<no) && (verified == false); ii++) {
      if (isOverlapBadQuality(ovl[ii]) == false) {
        if      ((ovl[ii].a_hang <= 0) && (ovl[ii].b_hang <= 0))   //  Left side dovetail
          IL.add(0,
                 fLen - -ovl[ii].b_hang);

        else if ((ovl[ii].a_hang >= 0) && (ovl[ii].b_hang >= 0))   //  Right side dovetail
          IL.add(ovl[ii].a_hang,
                 fLen - ovl[ii].a_hang);

        else if ((ovl[ii].a_hang >= 0) && (ovl[ii].b_hang <= 0))   //  I contain the other
          IL.add(ovl[ii].a_hang,
                 fLen - ovl[ii].a_hang - -ovl[ii].b_hang);

        else if ((ovl[ii].a_hang <= 0) && (ovl[ii].b_hang >= 0))   //  I am contained and thus now perfectly good!
          verified = true;

        else                                                       //  Huh?  Coding error.
          assert(0);
      }
    }

    if ((verified == true) ||            //  If verified, we're done.  It's good.
        (IL.numberOfIntervals() == 0))   //  If no intervals, it's a singleton.  Do nothing.
      continue;

    IL.merge();                          //  Merge the overlapping intervals.

    if (IL.numberOfIntervals() == 1)     //  One interval, so it's good.
      continue;

    setCoverageGap(fi);                  //  Bad regions detected!  Possible chimeric read.
  }
}



void
BestOverlapGraph::findErrorRateThreshold(void) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  //  Find and remember the erate for every best edge.  No way around this,
  //  since we need to find the median value.
  //
  //  If there are no best edges, find the overlap with the most matches and
  //  use that.  This shouldn't happen anymore.

  vector<double>  erates;

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap *b5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap *b3 = getBestEdgeOverlap(fi, true);

    if (isIgnored(fi) == true)
      continue;

    if (b5->readId() != 0)   erates.push_back(b5->erate());
    if (b3->readId() != 0)   erates.push_back(b3->erate());

    if ((b5->readId() == 0) &&
        (b3->readId() == 0)) {
      uint32      no    = 0;
      BAToverlap *ovl   = OC->getOverlaps(fi, no);
      uint32      bestM = 0;
      double      bestE = 0.0;

      for (uint32 oo=0; oo<no; oo++) {
        double  matches = (1 - ovl[oo].erate()) * RI->overlapLength(ovl[oo].a_iid, ovl[oo].b_iid, ovl[oo].a_hang, ovl[oo].b_hang);

        if (bestM < matches) {
          bestM = matches;
          bestE = ovl[oo].erate();
        }
      }

      if (no > 0)
        erates.push_back(bestE);
    }
  }

  sort(erates.begin(), erates.end());

  //  Find mean and stddev with an online calculation.

  stdDev<double>  edgeStats;

  for (uint32 ii=0; ii<erates.size(); ii++)
    edgeStats.insert(erates[ii]);

  _mean   = edgeStats.mean();
  _stddev = edgeStats.stddev();

  //  Find the median and absolute deviations.

  computeMedianAbsoluteDeviation(erates, _median, _mad, true);
  
  //  Compute an error limit based on the median or absolute deviation.
  //    Use median+mad is the median isn't zero.
  //    Use mean+stddev otherwise.
  //
  //    But either way, let the user lower this via the command line.

  double   Tinput  = _errorLimit;
  double   Tmean   = _mean   + _deviationGraph          * _stddev;
  double   Tmad    = _median + _deviationGraph * 1.4826 * _mad;
  double   Tpicked = (_median > 1e-10) ? Tmad : Tmean;

  if (Tinput < Tpicked)
    _errorLimit = Tinput;
  else
    _errorLimit = Tpicked;

  //  The real filtering is done on the next pass through findEdges().  Here, we're just collecting statistics.

  uint32  nFiltered[3] = {0};

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap *b5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap *b3 = getBestEdgeOverlap(fi, true);

    if ((b5->readId() == 0) &&    //  Read has no best edges, probably
        (b3->readId() == 0))      //  not corrected.
      continue;

    uint32  nf = ((b5->erate() > _errorLimit) +
                  (b3->erate() > _errorLimit));

    nFiltered[nf]++;
  }

  writeLog("\n");
  writeLog("ERROR RATES\n");
  writeLog("-----------\n");
  writeLog("\n");
  writeLog("                                            --------threshold-------\n");
  writeLog("%-12u"     "             fraction error     fraction         percent\n", edgeStats.size());
  writeLog("samples                          (1e-5)        error           error\n");
  writeLog("                 ----------------------       ------        --------\n");
  writeLog("command line                             ->   %6.2f       %8.4f%%%s\n",                                                1e5 * Tinput, 100.0 * Tinput, (_errorLimit == Tinput) ? "  (enabled)" : "");
  writeLog("mean + std.dev   %6.2f +- %3.0f * %6.2f  ->   %6.2f       %8.4f%%%s\n", 1e5 * _mean,   _deviationGraph, 1e5 * _stddev, 1e5 * Tmean,  100.0 * Tmean,  (_errorLimit == Tmean)  ? "  (enabled)" : "");
  writeLog("median + mad     %6.2f +- %3.0f * %6.2f  ->   %6.2f       %8.4f%%%s\n", 1e5 * _median, _deviationGraph, 1e5 * _mad,    1e5 * Tmad,   100.0 * Tmad,   (_errorLimit == Tmad)   ? "  (enabled)" : "");
  writeLog("\n");
  writeLog("\n");
  writeLog("BEST EDGE FILTERING\n");
  writeLog("-------------------\n");
  writeLog("\n");
  writeLog("Reads with both best edges below threshold: %9u\n", nFiltered[0]);
  writeLog("Reads with one  best edge  above threshold: %9u\n", nFiltered[1]);
  writeLog("Reads with both best edges above threshold: %9u\n", nFiltered[2]);
  writeLog("\n");
}



void
BestOverlapGraph::removeLopsidedEdges(const char *UNUSED(prefix)) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  //  Compare the size of our edges (this5 and this3) against the edges back
  //  from the reads those point to.
  //
  //  ------------------------->
  //                         |
  //                       this3
  //                         v
  //            <--back3-- -------------------->
  //
  //  If back5 points back to us, our 3' end is good.
  //  If back5 is of a comparable size to this3, our 3' end is good.

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {

    if ((isIgnored(fi)   == true) ||
        (isContained(fi) == true))
      continue;

    BestEdgeOverlap *this5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap *this3 = getBestEdgeOverlap(fi, true);

    //  Find the length of the overlap on either end.

    int32   this5ovlLen = RI->overlapLength(fi, this5->readId(), this5->ahang(), this5->bhang());
    int32   this3ovlLen = RI->overlapLength(fi, this3->readId(), this3->ahang(), this3->bhang());

    //  Find the edges from our best edges.  Note that 'back5' is NOT the
    //  edge out of the 5' end of that read; it is the edge that should point
    //  back to our 5' end.

    BestEdgeOverlap *back5 = getBestEdgeOverlap(this5->readId(), this5->read3p());
    BestEdgeOverlap *back3 = getBestEdgeOverlap(this3->readId(), this3->read3p());

    //  If both point back to us, we're done.  These must be symmetric, else
    //  overlapper is bonkers.

    if ((back5->readId() == fi) && (back5->read3p() == false) &&
        (back3->readId() == fi) && (back3->read3p() == true))
      continue;

    //  Complain loudly if we have a best overlap to a spur.

    if ((this5->readId() != 0) && (back5->readId() == 0))
      writeLog("WARNING: read %u 5' has overlap to spur read %u %c'!\n",
               fi, this5->readId(), this5->read3p() ? '3' : '5');

    if ((this3->readId() != 0) && (back3->readId() == 0))
      writeLog("WARNING: read %u 3' has overlap to spur read %u %c'!\n",
               fi, this3->readId(), this3->read3p() ? '3' : '5');

    //  Compute the length of those best overlaps...

    int32  back5ovlLen = RI->overlapLength(this5->readId(), back5->readId(), back5->ahang(), back5->bhang());
    int32  back3ovlLen = RI->overlapLength(this3->readId(), back3->readId(), back3->ahang(), back3->bhang());

    //  ...and compare.

    double  percDiff5 = 200.0 * abs(this5ovlLen - back5ovlLen) / (this5ovlLen + back5ovlLen);
    double  percDiff3 = 200.0 * abs(this3ovlLen - back3ovlLen) / (this3ovlLen + back3ovlLen);

    setLopsided5(fi, (percDiff5 > 5.0));
    setLopsided3(fi, (percDiff3 > 5.0));

#if 0
    if (isLopsided(fi) == false)
      writeLog("fi %8u -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%% -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%% -- ACCEPTED\n",
               fi,
               this5->readId(), this5->read3p() ? '3' : '5', this5ovlLen, back5->readId(), back5->read3p() ? '3' : '5', back5ovlLen, percDiff5,
               this3->readId(), this3->read3p() ? '3' : '5', this3ovlLen, back3->readId(), back3->read3p() ? '3' : '5', back3ovlLen, percDiff3);
#endif

    if (isLopsided(fi) == true)
      writeLog("lopsidedBest %8u -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%% -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%%\n",
               fi,
               this5->readId(), this5->read3p() ? '3' : '5', this5ovlLen, back5->readId(), back5->read3p() ? '3' : '5', back5ovlLen, percDiff5,
               this3->readId(), this3->read3p() ? '3' : '5', this3ovlLen, back3->readId(), back3->read3p() ? '3' : '5', back3ovlLen, percDiff3);
  }
}





uint32
BestOverlapGraph::spurDistance(BestEdgeOverlap *edge, uint32 limit, uint32 distance) {
  uint32  inid = edge->readId();
  bool    in3p = edge->read3p();
  bool    ot3p = (in3p == false) ? true : false;

  assert(isIgnored(inid) == false);   //  Edge to ignored read, ERROR!

  if (isCoverageGap(inid) == true)    //  Edge to a probably chimeric read?
    return(distance);                 //  Call it a spur.

  if (inid == 0)                      //  If no edge, we've terminated early,
    return(distance);                 //  so we've followed a path to a spur.

  if (distance == limit)              //  If we've hit the limit, not a spur,
    return(UINT32_MAX);               //  return infinity.

  //  Otherwise, follow the path.

  return(spurDistance(getBestEdgeOverlap(inid, ot3p), limit, distance+1));
}



//
//  Find reads that terminate at a spur after some short traversal.
//
//    A read is marked SPUR if there is no edge out of it on either end.
//
//    A read is marked SPURPATH if the path out of the read on that end
//    terminates (after some distance) at a spur read.
//
//  Note that SPUR reads have at least one SPURPATH end flagged.
//
void
BestOverlapGraph::removeSpannedSpurs(const char *prefix, uint32 spurDepth) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  set<uint32>  spurpath5;   //  spurpath is true if the edge out of this end
  set<uint32>  spurpath3;   //           leads to a dead-end spur.
  set<uint32>  spur;        //  spur     is true if this read is a spur.

  //  Compute the distance to a dead end.
  //    If zero, flag the read as a spur.
  //    If small, flag that read end as leading to a spur end.

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    if ((isIgnored(fi)   == true) ||   //  Ignored read, ignore.
        (isContained(fi) == true) ||   //  Contained read, ignore.
        (RI->isValid(fi) == false))    //  Unused read, ignore.
      continue;

    uint32  dist5 = spurDistance(getBestEdgeOverlap(fi, false), spurDepth);
    uint32  dist3 = spurDistance(getBestEdgeOverlap(fi,  true), spurDepth);

#if 0
    if (dist5 == 0)           {  writeLog("read %u 5' is a terminal spur\n", fi);  spur.insert(fi);       }
    if (dist3 == 0)           {  writeLog("read %u 3' is a termainl spur\n", fi);  spur.insert(fi);       }
    if (dist5  < spurDepth)   {  writeLog("read %u 5' is a spur path\n", fi);      spurpath5.insert(fi);  }
    if (dist3  < spurDepth)   {  writeLog("read %u 3' is a spur path\n", fi);      spurpath3.insert(fi);  }
#endif

    if (dist5 == 0)           spur.insert(fi);
    if (dist3 == 0)           spur.insert(fi);
    if (dist5  < spurDepth)   spurpath5.insert(fi);
    if (dist3  < spurDepth)   spurpath3.insert(fi);
  }

  //
  //  Do several iterations of spur path mark removal.
  //

  writeStatus("BestOverlapGraph()--   After initial scan, found:\n");
  writeStatus("BestOverlapGraph()--     %5u spur reads.\n", spur.size());
  writeStatus("BestOverlapGraph()--     %5u 5' spur paths.\n", spurpath5.size());
  writeStatus("BestOverlapGraph()--     %5u 3' spur paths.\n", spurpath3.size());

  uint32  n5pChanged = 1;
  uint32  n3pChanged = 1;

  for (uint32 iter=1; ((n5pChanged + n3pChanged > 0) && (iter <= 2 * spurDepth)); iter++) {
    //writeLog("\n");
    //writeLog("SPUR path removal iteration %u\n", iter);
    //writeLog("\n");

    n5pChanged = 0;
    n3pChanged = 0;

    //
    //  For any read that can get out on a non-spur-path end, mark the
    //  incoming edge as not a spur path:
    //
    //    If the 5' path out of me is not a spur path, then the
    //    edge INTO my 3' end is not a spur path.
    //

    for (uint32 fi=1; fi <= fiLimit; fi++) {
      if ((isIgnored(fi)     == true) ||   //  Ignored read, ignore.
          (isContained(fi)   == true) ||   //  Contained read, ignore.
          (isCoverageGap(fi) == true) ||   //  Suspected chimeric read, ignore.
          (RI->isValid(fi)   == false))    //  Unused read, ignore.
        continue;

      BestEdgeOverlap  *edge5 = getBestEdgeOverlap(fi, false);
      BestEdgeOverlap  *edge3 = getBestEdgeOverlap(fi,  true);

      bool              sp5  = (spurpath5.count(fi) > 0);
      bool              sp3  = (spurpath3.count(fi) > 0);

      bool              sp53 = (spurpath5.count(edge3->readId()) > 0);
      bool              sp33 = (spurpath3.count(edge3->readId()) > 0);
      bool              sp55 = (spurpath5.count(edge5->readId()) > 0);
      bool              sp35 = (spurpath3.count(edge5->readId()) > 0);

      //  Logging, only if enabled, and only if the spur path is actually removed.

#if 0
      if ((sp5 == false) && (edge3->read3p() == false) && (sp53 == true))   writeLog("SPUR path from read %u 5' removed\n", edge3->readId());
      if ((sp5 == false) && (edge3->read3p() ==  true) && (sp33 == true))   writeLog("SPUR path from read %u 3' removed\n", edge3->readId());
      if ((sp3 == false) && (edge5->read3p() == false) && (sp55 == true))   writeLog("SPUR path from read %u 5' removed\n", edge5->readId());
      if ((sp3 == false) && (edge5->read3p() ==  true) && (sp35 == true))   writeLog("SPUR path from read %u 3' removed\n", edge5->readId());
#endif

      //  Remove spur-path marks.

      if ((sp5 == false) && (edge3->read3p() == false) && (sp53 == true))   spurpath5.erase(edge3->readId());
      if ((sp5 == false) && (edge3->read3p() ==  true) && (sp33 == true))   spurpath3.erase(edge3->readId());
      if ((sp3 == false) && (edge5->read3p() == false) && (sp55 == true))   spurpath5.erase(edge5->readId());
      if ((sp3 == false) && (edge5->read3p() ==  true) && (sp35 == true))   spurpath3.erase(edge5->readId());
    }

    //
    //  Compute new edges that don't point to a spur or spurpath, unless
    //  those are the only edges it has.
    //
    //  1) Over all overlaps for this read,
    //     ignore crappy edges, and edges to
    //     contains, etc.                        --------->
    //                                                  |
    //  2) Score the edge if the exit from the          v         (score this edge if
    //     read it points to is not a spur path.      -------->    this end isn't a spur-path) 
    //
    //     Note that spur reads are also spur-path reads.
    //
    //
    //  3) After all overlaps are processed, if
    //     a read still has no best edge, restore
    //     the previous best edge -- this will be
    //     an edge to a spur or spur-path read, but
    //     it's the ONLY path we have.
    //

    memset(_best5score, 0, sizeof(uint64) * (fiLimit + 1));     //  Clear all edge scores.
    memset(_best3score, 0, sizeof(uint64) * (fiLimit + 1));     //  Clear all edge scores.

    for (uint32 fi=1; fi <= fiLimit; fi++) {
      if ((isIgnored(fi)   == true) ||   //  Ignored read, ignore.
          (isContained(fi) == true) ||   //  Contained read, ignore.
          (RI->isValid(fi) == false))    //  Unused read, ignore.
        continue;

      uint32      no  = 0;
      BAToverlap *ovl = OC->getOverlaps(fi, no);

      BestEdgeOverlap  orig5 = *getBestEdgeOverlap(fi, false);    //  Remember the previous best edges.
      BestEdgeOverlap  orig3 = *getBestEdgeOverlap(fi,  true);

      getBestEdgeOverlap(fi, false)->clear();                     //  Then clear them.
      getBestEdgeOverlap(fi,  true)->clear();

      for (uint32 ii=0; ii<no; ii++) {                            //  Over all overlaps for this read,
        if ((ovl[ii].isDovetail()         == false) ||            //    Ignore non-dovetail and crappy overlaps.
            (isOverlapBadQuality(ovl[ii]) == true))               //    They can't form best edges.
          continue;

        if ((isIgnored(ovl[ii].b_iid)     == true) ||             //    Ignore overlaps to ignored reads.
            (isContained(ovl[ii].b_iid)   == true) ||             //    Ignore overlaps to contained and chimeric reads.
            (isCoverageGap(ovl[ii].b_iid) == true) ||             //    We don't want best edges to these reads.
            (RI->isValid(fi)              == false))              //    Ignore overlaps to reads that don't exist.
          continue;

        bool   Aend5 = ovl[ii].AEndIs5prime();
        bool   Aend3 = ovl[ii].AEndIs3prime();

        bool   Bend5 = ovl[ii].BEndIs5prime();
        bool   Bend3 = ovl[ii].BEndIs3prime();

        bool   sp5c  = (spurpath5.count(ovl[ii].b_iid) > 0);
        bool   sp3c  = (spurpath3.count(ovl[ii].b_iid) > 0);

        //  Log the edges we are skipping, if enabled.  This isn't as useful
        //  as you'd think, since it catches EVERY edge in the path to a
        //  spur, not just the best.

#if 0
        if ((Aend5 == true) && (Bend5 ==  true) && (sp3c == true))   writeLog("edge from %u 5' to %u 5' ignored; outgoing 3' edge leads to a spur\n", ovl[ii].a_iid, ovl[ii].b_iid);
        if ((Aend5 == true) && (Bend5 == false) && (sp5c == true))   writeLog("edge from %u 5' to %u 3' ignored; outgoing 5' edge leads to a spur\n", ovl[ii].a_iid, ovl[ii].b_iid);

        if ((Aend3 == true) && (Bend5 ==  true) && (sp3c == true))   writeLog("edge from %u 3' to %u 5' ignored; outgoing 3' edge leads to a spur\n", ovl[ii].a_iid, ovl[ii].b_iid);
        if ((Aend3 == true) && (Bend5 == false) && (sp5c == true))   writeLog("edge from %u 3' to %u 3' ignored; outgoing 5' edge leads to a spur\n", ovl[ii].a_iid, ovl[ii].b_iid);
#endif

        //  Score the edge.

        if ((Aend5 == true) && (Bend5 ==  true) && (sp3c == false))   scoreEdge(ovl[ii]);
        if ((Aend5 == true) && (Bend5 == false) && (sp5c == false))   scoreEdge(ovl[ii]);

        if ((Aend3 == true) && (Bend5 ==  true) && (sp3c == false))   scoreEdge(ovl[ii]);
        if ((Aend3 == true) && (Bend5 == false) && (sp5c == false))   scoreEdge(ovl[ii]);
      }

      //  All edges scored.  If no best edge found, log and reset to whatever
      //  previous edge was there.

#if 0
      if (getBestEdgeOverlap(fi, false)->readId() == 0)   writeLog("restore edge out of %u 5' to previous best %u %c'\n", fi, orig5.readId(), orig5.read3p() ? '3' : '5');
      if (getBestEdgeOverlap(fi,  true)->readId() == 0)   writeLog("restore edge out of %u 3' to previous best %u %c'\n", fi, orig3.readId(), orig3.read3p() ? '3' : '5');
#endif

      if (getBestEdgeOverlap(fi, false)->readId() == 0)   *getBestEdgeOverlap(fi, false) = orig5;
      if (getBestEdgeOverlap(fi,  true)->readId() == 0)   *getBestEdgeOverlap(fi,  true) = orig3;

      //  Count the number of edges we switched from following a spur path to a normal path.
      //  These are just the edges that are different than their original version.

      if (*getBestEdgeOverlap(fi, false) != orig5)  n5pChanged++;
      if (*getBestEdgeOverlap(fi,  true) != orig3)  n3pChanged++;
    }

    //
    //  Done with an interation.  Report a bit of status.
    //

    writeStatus("BestOverlapGraph()--   After iteration %u, found:\n", iter);
    writeStatus("BestOverlapGraph()--     %5u spur reads.\n", spur.size());
    writeStatus("BestOverlapGraph()--     %5u 5' spur paths;  %5u 5' edges changed to avoid a spur path.\n", spurpath5.size(), n5pChanged);
    writeStatus("BestOverlapGraph()--     %5u 3' spur paths;  %5u 5' edges changed to avoid a spur path.\n", spurpath3.size(), n3pChanged);
  }

  //
  //  Set the spur flag for any spur reads, and log the result.
  //

  FILE   *F = AS_UTL_openOutputFile(prefix, '.', "best.spurs");

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    bool s5 = (spurpath5.count(fi) > 0);
    bool s3 = (spurpath3.count(fi) > 0);

    if ((s5 == false) && (s3 == false))  //  No spur mark, so not a spur.
      continue;

    //  Set the spur flag for this read.

    setSpur(fi);

    //  Log it.  It's 'terminal' if there is no edge out of that end.

    bool t5 = (getBestEdgeOverlap(fi, false)->readId() == 0);
    bool t3 = (getBestEdgeOverlap(fi,  true)->readId() == 0);

    if (s5 == true)
      fprintf(F, "%u 5' is a %s spur end\n", fi, (t5) ? "terminal" : "non-terminal");

    if (s3 == true)
      fprintf(F, "%u 3' is a %s spur end\n", fi, (t3) ? "terminal" : "non-terminal");
  }

  fclose(F);

  //
  //  Remove edges from spur reads to good reads.
  //  This prevents spurs from breaking contigs.
  //

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    bool ss = (spur.count(fi) > 0);
    bool s5 = (spurpath5.count(fi) > 0);
    bool s3 = (spurpath3.count(fi) > 0);

    if ((ss == false) &&   //  Read fi isn't a spur, and both
        (s5 == false) &&   //  ends aren't leading to a spur;
        (s3 == false))     //  keep all edges intact.
      continue;

    //  Read fi is either a spur or a spur-path read.

    //  Remove edges from this read to any non-spur read, unless that
    //  non-spur read itself has an edge back to us.  If it does, then, by
    //  construction, this is the only place that read can go, and so our
    //  edge to it is valid.

    //  If our 5' edge points to a non-spur read, delete the edge unless it
    //  points back to us.
    {
      BestEdgeOverlap  *edge5 = getBestEdgeOverlap(fi, false);
      uint32            read5 = edge5->readId();
      bool              spur5 = ((spur.count(read5) + spurpath5.count(read5) + spurpath3.count(read5)) > 0);

      if (spur5 == false) {
        BestEdgeOverlap  *backedge = getBestEdgeOverlap(read5, edge5->read3p());

        if ((backedge->readId() != fi) ||
            (backedge->read3p() != false)) {
          //if (read5 != 0)
          //  writeLog("DELETE edge from spur %u 5' to non-spur %u %c'\n", fi, read5, edge5->read3p() ? '3' : '5');
          edge5->clear();
        }
      }
    }

    //  If our 3' edge points to a non-spur read, delete the edge unless it
    //  points back to us.
    {
      BestEdgeOverlap  *edge3 = getBestEdgeOverlap(fi, true);
      uint32            read3 = edge3->readId();
      bool              spur3 = ((spur.count(read3) + spurpath5.count(read3) + spurpath3.count(read3)) > 0);

      if (spur3 == false) {
        BestEdgeOverlap  *backedge = getBestEdgeOverlap(read3, edge3->read3p());

        if ((backedge->readId() != fi) ||
            (backedge->read3p() != true)) {
          //if (read3 != 0)
          //  writeLog("DELETE edge from spur %u 3' to non-spur %u %c'\n", fi, read3, edge3->read3p() ? '3' : '5');
          edge3->clear();
        }
      }
    }
  }
}



void
BestOverlapGraph::findEdges(void) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  //  Reset our scores.

  memset(_best5score, 0, sizeof(uint64) * (fiLimit + 1));
  memset(_best3score, 0, sizeof(uint64) * (fiLimit + 1));

  //  Remove containment flags so we can recompute them.

  for (uint32 fi=1; fi <= fiLimit; fi++)
    setContained(fi, false);

  //  One pass through all reads to flag any that are in a containment
  //  relationship.

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, no);

    if (isIgnored(fi) == true)
      continue;

    for (uint32 ii=0; ii<no; ii++) {
      if (isOverlapBadQuality(ovl[ii]))      //  Ignore crappy overlaps.
        continue;

      if ((ovl[ii].a_hang == 0) &&           //  If an exact overlap, make
          (ovl[ii].b_hang == 0) &&           //  the lower ID be contained.
          (ovl[ii].a_iid > ovl[ii].b_iid))   //  (Ignore if exact and this
        continue;                            //   ID is larger.)

      if ((ovl[ii].a_hang > 0) ||            //  Ignore if A is not
          (ovl[ii].b_hang < 0))              //  contained in B.
        continue;

      setContained(ovl[ii].a_iid);
    }
  }

  //  A second pass to score and find the best edge for each end.

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, no);

    _best5score[fi] = 0;                     //  Reset scores.
    _best3score[fi] = 0;                     //

    _reads[fi]._best5.clear();               //  Clear existing
    _reads[fi]._best3.clear();               //  best edges.

    if ((isIgnored(fi)   == true) ||         //  Ignore ignored reads.
        (isContained(fi) == true))           //  Ignore contained reads.
      continue;

    for (uint32 ii=0; ii<no; ii++)           //  Compute scores for all overlaps
      scoreEdge(ovl[ii]);                    //  and remember the best.
  }
}



void
BestOverlapGraph::removeContainedDovetails(void) {
  uint32  fiLimit = RI->numReads();
  uint32  fix5    = 0;
  uint32  fix3    = 0;

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    if (isContained(fi) == true) {
      if (getBestEdgeOverlap(fi, false)->readId() != 0)
        fix5++;
      if (getBestEdgeOverlap(fi, true)->readId() != 0)
        fix3++;

      getBestEdgeOverlap(fi, false)->clear();
      getBestEdgeOverlap(fi, true) ->clear();
    }
  }

  fprintf(stderr, "Cleared %u 5' and %u 3' best edges on contained reads.\n", fix5, fix3);
}



BestOverlapGraph::BestOverlapGraph(double            erateGraph,
                                   double            deviationGraph,
                                   const char       *prefix,
                                   bool              filterSuspicious,
                                   bool              filterHighError,
                                   bool              filterLopsided,
                                   bool              filterSpur,
                                   uint32            spurDepth,
                                   BestOverlapGraph *BOG) {

  writeStatus("\n");
  writeStatus("BestOverlapGraph()-- Computing Best Overlap Graph.\n");
  writeStatus("BestOverlapGraph()-- Allocating best edges (" F_SIZE_T "MB).\n",
           ((2 * sizeof(BestEdgeOverlap) * (RI->numReads() + 1)) >> 20));

  _reads               = new BestEdgeRead [RI->numReads() + 1];

  _best5score          = new uint64 [RI->numReads() + 1];   //  Cleared in findEdges().
  _best3score          = new uint64 [RI->numReads() + 1];

  _mean                = erateGraph;
  _stddev              = 0.0;

  _median              = erateGraph;
  _mad                 = 0.0;

  _errorLimit          = erateGraph;

  _erateGraph          = erateGraph;
  _deviationGraph      = deviationGraph;

  //  If there is a BOG supplied, copy the reads from there to here.
  //
  //  This sets the status of orphans and bubbles.  Then we'll further flag
  //  those as ignored.
  //
  //  Edges will be recomputed on the first findEdges() call below.

  if (BOG) {
    memcpy(_reads, BOG->_reads, sizeof(BestEdgeRead) * (RI->numReads() + 1));

    for (uint32 fi=1; fi <= RI->numReads(); fi++) {
      if ((isOrphan(fi) == true) ||
          (isBubble(fi) == true)) {
        writeLog("IGNORE read %u %s %s\n", fi, isOrphan(fi) ? "orphan" : "", isBubble(fi) ? "bubble" : "");
        setIgnored(fi);
      }
    }
  }

  //
  //  Find initial edges, report an initial graph, then analyze those edges
  //  to set a cutoff on overlap quality used for graph building.  Once we
  //  have the magic error rate limit, recompute best edges.
  //
  //  This should be done before removing suspicious reads, so we can skip
  //  high-error overlaps (that would otherwise mask a problematic read).
  //  On one hifi set, there was no difference between doing it as here,
  //  and removing suspicious reads first.
  //

  findEdges();
  reportBestEdges(prefix, "initial");

  if (filterHighError) {
    writeStatus("BestOverlapGraph()-- Filtering high error edges.\n");

    findErrorRateThreshold();
    findEdges();

    writeStatus("BestOverlapGraph()--   Ignore overlaps with more than %.6f%% error.\n", 100.0 * _errorLimit);
  }

  else {
    writeStatus("BestOverlapGraph()-- NOT filtering high error edges.\n");
  }

  //
  //  Mark reads as suspicious if they are not fully covered by (good) overlaps.
  //

  if (filterSuspicious) {
    writeStatus("BestOverlapGraph()-- Filtering reads with a gap in overlap coverage.\n");

    removeReadsWithCoverageGap(prefix);
    findEdges();

    writeStatus("BestOverlapGraph()--   %u reads removed.\n", numCoverageGap());
  }

  else {
    writeStatus("BestOverlapGraph()-- NOT filtering reads with suspicious overlap patterns.\n");
  }

  //
  //  Report initial statistics.
  //

  reportEdgeStatistics(prefix, "INITIAL");

  //
  //  Mark reads as suspicious if the length of the best edge out is very
  //  different than the length of the best edge that should be back to us.
  //  E.g., if readA has best edge to readB (of length lenAB), but readB has
  //  best edge to readC (of length lenBC), and lenAB is much shorter than
  //  lenBC, then something is wrong with readA.
  //

  if (filterLopsided) {
    writeStatus("BestOverlapGraph()-- Filtering reads with lopsided best edges.\n");

    removeLopsidedEdges(prefix);
    findEdges();

    writeStatus("BestOverlapGraph()--   %u reads have lopsided edges.\n", numLopsided());
  }

  else {
    writeStatus("BestOverlapGraph()-- NOT filtering reads with lopsided best edges.\n");
  }

  //
  //  Now that we know the max difference allowed in overlaps, mark reads as
  //  spurs if they're spurs.  Then don't find best edges to them.
  //
  //  Spanned spurs searches the graph for dead ends, then backs up to a
  //  previous branch.  It (supposedly) does this such that the spur path is
  //  preserved unless a branch is found.
  //
  //  Do NOT call findEdges() after this.  It will RESET all the work done by
  //  removeSpannedSpurs().
  //

  if (filterSpur) {
    writeStatus("BestOverlapGraph()-- Filtering spur reads.\n");

    removeSpannedSpurs(prefix, spurDepth);
  }

  else {
    writeStatus("BestOverlapGraph()-- NOT filtering spur reads.\n");
  }

  //
  //  Do some final boring cleanup:
  //    remove best edges associated with contained reads.
  //

  removeContainedDovetails();

  //
  //  Report filtering and final statistics.
  //

  writeLog("\n");
  writeLog("EDGE FILTERING\n");
  writeLog("-------- ------------------------------------------\n");
  writeLog("%8u reads are ignored\n",  numIgnored());
  writeLog("%8u reads have a gap in overlap coverage\n", numCoverageGap());
  writeLog("%8u reads have lopsided best edges\n", numLopsided());

  reportBestEdges(prefix, "best");

  reportEdgeStatistics(prefix, "FINAL");

  //  Done with scoring data.

  delete [] _best5score;    _best5score = NULL;
  delete [] _best3score;    _best3score = NULL;

  setLogFile(prefix, NULL);
}



void
BestOverlapGraph::reportEdgeStatistics(const char *prefix, const char *label) {
  uint32  fiLimit      = RI->numReads();
  uint32  numThreads   = omp_get_max_threads();
  uint32  blockSize    = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  uint32  nContained   = 0;
  uint32  nSingleton   = 0;
  uint32  nSpur        = 0;
  uint32  nSpur1Mutual = 0;
  uint32  nBoth        = 0;
  uint32  nBoth1Mutual = 0;
  uint32  nBoth2Mutual = 0;

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap *this5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap *this3 = getBestEdgeOverlap(fi, true);

    //  Count contained reads

    if (isContained(fi)) {
      nContained++;
      continue;
    }

    //  Count singleton reads

    if ((this5->readId() == 0) && (this3->readId() == 0)) {
      nSingleton++;
      continue;
    }

    //  Compute mutual bestedness

    bool  mutual5 = false;
    bool  mutual3 = false;

    if (this5->readId() != 0) {
      BestEdgeOverlap *that5 = getBestEdgeOverlap(this5->readId(), this5->read3p());

      mutual5 = ((that5->readId() == fi) && (that5->read3p() == false));
    }

    if (this3->readId() != 0) {
      BestEdgeOverlap *that3 = getBestEdgeOverlap(this3->readId(), this3->read3p());

      mutual3 = ((that3->readId() == fi) && (that3->read3p() == true));
    }

    //  Compute spur, and mutual best

    if ((this5->readId() == 0) ||
        (this3->readId() == 0)) {
      nSpur++;
      nSpur1Mutual += (mutual5 || mutual3) ? 1 : 0;
      continue;
    }

    //  Otherwise, both edges exist

    nBoth++;
    nBoth1Mutual +=  (mutual5 != mutual3) ? 1 : 0;
    nBoth2Mutual += ((mutual5 == true) && (mutual3 == true)) ? 1 : 0;
  }

  writeLog("\n");
  writeLog("%s EDGES\n", label);
  writeLog("-------- ----------------------------------------\n");
  writeLog("%8u reads are contained\n", nContained);
  writeLog("%8u reads have no best edges (singleton)\n", nSingleton);
  writeLog("%8u reads have only one best edge (spur) \n", nSpur);
  writeLog("         %8u are mutual best\n", nSpur1Mutual);
  writeLog("%8u reads have two best edges \n", nBoth);
  writeLog("         %8u have one mutual best edge\n", nBoth1Mutual);
  writeLog("         %8u have two mutual best edges\n", nBoth2Mutual);
  writeLog("\n");
}



void
BestOverlapGraph::reportBestEdges(const char *prefix, const char *label) {
  char  N[FILENAME_MAX];

  //  Open output files.

  snprintf(N, FILENAME_MAX, "%s.%s.edges",     prefix, label);   FILE *BE = AS_UTL_openOutputFile(N);
  snprintf(N, FILENAME_MAX, "%s.%s.edges.gfa", prefix, label);   FILE *BG = AS_UTL_openOutputFile(N);

  //  Write best edges, singletons and suspicious edges.

  if (BE) {
    fprintf(BE, "readID   libID flags      5' edge M   length    eRate      5' edge M   length    eRate\n");
    fprintf(BE, "-------- ----- -----  ----------- - -------- --------  ----------- - -------- --------\n");

    for (uint32 id=1; id<RI->numReads() + 1; id++) {
      BestEdgeOverlap *e5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *e3 = getBestEdgeOverlap(id, true);

      BestEdgeOverlap  *e5back = getBestEdgeOverlap(e5->readId(), e5->read3p());   //  Get edges that are
      BestEdgeOverlap  *e3back = getBestEdgeOverlap(e3->readId(), e3->read3p());   //  potentially back to us.

      bool  e5e = (e5->readId() == 0) ? false : true;
      bool  e3e = (e3->readId() == 0) ? false : true;

      double e5err = AS_OVS_decodeEvalue(e5->evalue());
      double e3err = AS_OVS_decodeEvalue(e3->evalue());

      uint32 e5len = RI->overlapLength(id, e5->readId(), e5->ahang(), e5->bhang());
      uint32 e3len = RI->overlapLength(id, e3->readId(), e3->ahang(), e3->bhang());

      char  e5mutual = ((e5back->readId() == id) && (e5back->read3p() == false)) ? 'M' : ' ';
      char  e3mutual = ((e3back->readId() == id) && (e3back->read3p() ==  true)) ? 'M' : ' ';

      if (RI->readLength(id) == 0)
        continue;

      if      ((e5e == false) && (e3e == false)) {
        fprintf(BE, "%-8u %5u %c%c%c%c%c  -------- -- - -------- --------  -------- -- - -------- --------\n",
                id,
                RI->libraryIID(id),
                (isContained(id))   ? 'C' : '-',
                (isIgnored(id))     ? 'I' : '-',
                (isCoverageGap(id)) ? 'G' : '-',
                (isLopsided(id))    ? 'L' : '-',
                (isSpur(id))        ? 'S' : '-');
      }

      else if ((e5e == false) && (e3e ==  true)) {
        fprintf(BE, "%-8u %5u %c%c%c%c%c  -------- -- - -------- --------  %8u %c' %c %8u %8.6f\n",
                id,
                RI->libraryIID(id),
                (isContained(id))   ? 'C' : '-',
                (isIgnored(id))     ? 'I' : '-',
                (isCoverageGap(id)) ? 'G' : '-',
                (isLopsided(id))    ? 'L' : '-',
                (isSpur(id))        ? 'S' : '-',
                e3->readId(), e3->read3p() ? '3' : '5', e3mutual, e3len, e3err);
      }

      else if ((e5e ==  true) && (e3e == false)) {
        fprintf(BE, "%-8u %5u %c%c%c%c%c  %8u %c' %c %8u %8.6f  -------- -- - -------- --------\n",
                id,
                RI->libraryIID(id),
                (isContained(id))   ? 'C' : '-',
                (isIgnored(id))     ? 'I' : '-',
                (isCoverageGap(id)) ? 'G' : '-',
                (isLopsided(id))    ? 'L' : '-',
                (isSpur(id))        ? 'S' : '-',
                e5->readId(), e5->read3p() ? '3' : '5', e5mutual, e5len, e5err);
      }

      else if ((e5e ==  true) && (e3e ==  true)) {
        fprintf(BE, "%-8u %5u %c%c%c%c%c  %8u %c' %c %8u %8.6f  %8u %c' %c %8u %8.6f\n",
                id,
                RI->libraryIID(id),
                (isContained(id))   ? 'C' : '-',
                (isIgnored(id))     ? 'I' : '-',
                (isCoverageGap(id)) ? 'G' : '-',
                (isLopsided(id))    ? 'L' : '-',
                (isSpur(id))        ? 'S' : '-',
                e5->readId(), e5->read3p() ? '3' : '5', e5mutual, e5len, e5err,
                e3->readId(), e3->read3p() ? '3' : '5', e3mutual, e3len, e3err);
      }

      else {
        assert(0);
      }
    }
  }

  //  Write best edge graph.

  if (BG) {
    fprintf(BG, "H\tVN:Z:1.0\n");

    //  First, write the sequences used.  The sequence can be used as either
    //  a source node or a destination node (or both).

    set<uint32>  used;

    for (uint32 id=1; id<RI->numReads() + 1; id++) {
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      if ((bestedge5->readId() == 0) &&   //  Ignore singletons.
          (bestedge3->readId() == 0) &&
          (isContained(id) == false))
        continue;

      if (isContained(id) == true)        //  Ignore contained reads.
        continue;

      //if (isLopsided(id) == true)         //  Ignore lopsided.
      //  continue;

      //  Remember the source and destination of this edge.
      used.insert(id);
      used.insert(bestedge5->readId());
      used.insert(bestedge3->readId());
    }

    for (set<uint32>::iterator it=used.begin(); it != used.end(); it++)
      if (*it != 0)
        fprintf(BG, "S\tread%08u\t*\tLN:i:%u\n", *it, RI->readLength(*it));

    //  Now, report edges.  GFA wants edges in exactly this format:
    //
    //       -------------
    //             -------------
    //
    //  with read orientation given by +/-.  Conveniently, this is what we've saved (for the edges).

    for (uint32 id=1; id<RI->numReads() + 1; id++) {
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      if ((bestedge5->readId() == 0) &&   //  Ignore singletons.
          (bestedge3->readId() == 0) &&
          (isContained(id) == false))
        continue;

      if (isContained(id) == true)        //  Ignore contained reads.
        continue;

      //if (isLopsided(id) == true)         //  Ignore lopsided.
      //  continue;

      if (bestedge5->readId() != 0) {
        int32  ahang   = bestedge5->ahang();
        int32  bhang   = bestedge5->bhang();
        int32  olaplen = RI->overlapLength(id, bestedge5->readId(), bestedge5->ahang(), bestedge5->bhang());

        if ((ahang > 0) || (bhang > 0))
          fprintf(stderr, "BAD 5' overlap from read %u to read %u %c': hangs %d %d\n",
                  id,
                  bestedge3->readId(),
                  bestedge3->read3p() ? '3' : '5',
                  bestedge3->ahang(), bestedge3->bhang());
        assert((ahang <= 0) && (bhang <= 0));  //  ALL 5' edges should be this.

        fprintf(BG, "L\tread%08u\t-\tread%08u\t%c\t%uM\n",
                id,
                bestedge5->readId(), bestedge5->read3p() ? '-' : '+',
                olaplen);
      }

      if (bestedge3->readId() != 0) {
        int32  ahang   = bestedge3->ahang();
        int32  bhang   = bestedge3->bhang();
        int32  olaplen = RI->overlapLength(id, bestedge3->readId(), bestedge3->ahang(), bestedge3->bhang());

        if ((ahang < 0) || (bhang < 0))
          fprintf(stderr, "BAD 3' overlap from read %u to read %u %c': hangs %d %d\n",
                  id,
                  bestedge3->readId(),
                  bestedge3->read3p() ? '3' : '5',
                  bestedge3->ahang(), bestedge3->bhang());
        assert((ahang >= 0) && (bhang >= 0));  //  ALL 3' edges should be this.

        fprintf(BG, "L\tread%08u\t+\tread%08u\t%c\t%uM\n",
                id,
                bestedge3->readId(), bestedge3->read3p() ? '-' : '+',
                RI->overlapLength(id, bestedge3->readId(), bestedge3->ahang(), bestedge3->bhang()));
      }
    }
  }

  //  Close all the files.

  AS_UTL_closeFile(BE);
  AS_UTL_closeFile(BG);
}



inline
void
logEdgeScore(BAToverlap   &olap,
             const char   *message) {
  if ((logFileFlagSet(LOG_OVERLAP_SCORING)) &&    //  Report logging only if enabled, and only
      ((olap.a_iid != 0) ||                       //  for specific annoying reads.  (By default,
       (olap.a_iid == 97202) ||                   //  report for all reads; modify as needed).
       (olap.a_iid == 30701)))
    writeLog("scoreEdge()-- %6d %c' to %6d %c' -- hangs %6d %6d err %8.6f -- %s\n",
             olap.a_iid, olap.AEndIs3prime() ? '3' : '5',
             olap.b_iid, olap.BEndIs3prime() ? '3' : '5',
             olap.a_hang,
             olap.b_hang,
             olap.erate(),
             message);
}


void
BestOverlapGraph::scoreEdge(BAToverlap& olap) {

  if ((olap.isDovetail()         == false) ||   //  Ignore non-dovetail overlaps.
      (isContained(olap.b_iid)   == true) ||    //  Ignore edges into contained.
      (isCoverageGap(olap.b_iid) == true))      //  Ignore edges into coverage gap reads.
    return;

  if (isIgnored(olap.a_iid) == true) {          //  Ignore ignored reads.  This should
    assert(0);                                  //  never actually happen.
    logEdgeScore(olap, "ignored");
    return;
  }

  if (isIgnored(olap.b_iid) == true) {          //  Ignore ignored reads.  This could
    logEdgeScore(olap, "ignored");              //  happen; it's just easier to filter
    return;                                     //  them out here.
  }

  if (isOverlapBadQuality(olap) == true) {      //  Ignore the overlap if it is
    logEdgeScore(olap, "bad quality");          //  bad quality.
    return;
  }

  //  Compute the score for this overlap, and remember this overlap if the
  //  score is the best.

  uint64           newScr = scoreOverlap(olap);
  bool             a3p    = olap.AEndIs3prime();
  BestEdgeOverlap *best   = getBestEdgeOverlap(olap.a_iid, a3p);
  uint64          &score  = (a3p) ? (_best3score[olap.a_iid]) : (_best5score[olap.a_iid]);

  assert(newScr > 0);

  if (newScr <= score) {
    logEdgeScore(olap, "worse");
  }

  else {
    logEdgeScore(olap, "BEST");
    best->set(olap);
    score = newScr;
  }
}



bool
BestOverlapGraph::isOverlapBadQuality(BAToverlap& olap) {
  bool   isBad = true;
  bool   isIgn = false;

  if (olap.erate() <= _errorLimit)               //  Our only real test is on
    isBad = false;                               //  overlap error rate.

  if ((RI->isValid(olap.a_iid) == false) ||      //  But if either read is not valid,
      (RI->isValid(olap.b_iid) == false))        //  the overlap is bad.  This should
    isIgn = true;                                //  never occur.

  if ((isIgnored(olap.a_iid) == true) ||         //  But if either read is ignored,
      (isIgnored(olap.b_iid) == true))           //  the overlap is also bad.
    isIgn = true;

  olap.filtered = ((isBad == true) ||            //  The overlap is filtered out ("bad")
                   (isIgn == true));             //  if it's either Bad or Ignored.

  //  Now just a bunch of logging.

  if ((logFileFlagSet(LOG_OVERLAP_SCORING)) &&   //  Write the log only if enabled and the read
      ((olap.a_iid != 0) ||                      //  is specifically annoying (the default is to
       (olap.a_iid == 0) ||                      //  log for all reads (!= 0).
       (olap.a_iid == 0)))
    writeLog("isOverlapBadQuality()-- %6d %6d %c  hangs %6d %6d err %.3f -- %s\n",
             olap.a_iid, olap.b_iid,
             olap.flipped ? 'A' : 'N',
             olap.a_hang,
             olap.b_hang,
             olap.erate(),
             (isBad == true) ? ("REJECT!")
                             : ((isIgn == true) ? "good quality, but ignored"
                                                : "good quality"));

  return(olap.filtered);
}



uint64
BestOverlapGraph::scoreOverlap(BAToverlap& olap) {
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
  leng = RI->overlapLength(olap.a_iid, olap.b_iid, olap.a_hang, olap.b_hang);
#endif

  if (olap.a_hang > 0)
    leng = RI->readLength(olap.a_iid) - olap.a_hang;
  else
    leng = RI->readLength(olap.a_iid) + olap.b_hang;

  //  Convert the length into an expected number of matches.

#if 0
  assert(olap.erate() <= 1.0);
  leng  -= leng * olap.erate();
#endif

  //  And finally shift it to the correct place in the word.

  leng <<= AS_MAX_EVALUE_BITS;

  return(leng | rate);
}

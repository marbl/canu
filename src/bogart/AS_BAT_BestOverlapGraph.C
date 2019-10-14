
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
BestOverlapGraph::removeSuspicious(const char *UNUSED(prefix)) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32               no  = 0;
    BAToverlap          *ovl = OC->getOverlaps(fi, no);

    bool                 verified = false;
    intervalList<int32>  IL;

    uint32               fLen = RI->readLength(fi);

    if (fLen == 0)   //  Unused read, ignore.
      continue;

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

    if (verified == true)
      continue;

    //  Not verified.  Merge the intervals and declare it good if there
    //  is exactly one.

    IL.merge();

    //  If no intervals, it's a singleton.  Do nothing.
    if (IL.numberOfIntervals() == 0)
      continue;

    //  If one interval, it's good!  At worst, a spur, but we'll find that out later.
    if (IL.numberOfIntervals() == 1)
      continue;

    //  Uh oh.  Bad regions internal to the read!  Possible chimeric read.
#pragma omp critical (suspInsert)
    {
      writeLog("read %u has suspicious overlap pattern, ignored.\n", fi);
      _suspicious.insert(fi);
    }
  }
}



void
BestOverlapGraph::removeHighErrorBestEdges(void) {
  uint32  fiLimit    = RI->numReads();
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

    if (b5->readId() != 0)   edgeStats.insert(erates[eratesLen++] = b5->erate());
    if (b3->readId() != 0)   edgeStats.insert(erates[eratesLen++] = b3->erate());

    //  If there are NO best edges, find the overlap with the most matches and use that.

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
        edgeStats.insert(erates[eratesLen++] = bestE);
    }
  }

  _mean   = edgeStats.mean();
  _stddev = edgeStats.stddev();

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

  //  Compute an error limit based on the median or absolute deviation.

  double   Tmean = _mean   + _deviationGraph * _stddev;
  double   Tmad  = _median + _deviationGraph * 1.4826 * _mad;

  _errorLimit = min(_errorLimit, (_median > 1e-10) ? Tmad : Tmean);

  //  The real filtering is done on the next pass through findEdges().  Here, we're just collecting statistics.

  uint32  oneFiltered = 0;
  uint32  twoFiltered = 0;

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap *b5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap *b3 = getBestEdgeOverlap(fi, true);

    bool  b5filtered = (b5->erate() > _errorLimit);
    bool  b3filtered = (b3->erate() > _errorLimit);

    if      (b5filtered && b3filtered)
      _n2EdgeFiltered++;
    else if (b5filtered || b3filtered)
      _n1EdgeFiltered++;
  }

  writeLog("\n");
  writeLog("ERROR RATES (%u samples)\n", edgeStats.size());
  writeLog("-----------\n");
  writeLog("mean   %10.8f stddev %10.8f -> %10.8f fraction error = %10.6f%% error\n", _mean,   _stddev, Tmean, 100.0 * Tmean);
  writeLog("median %10.8f mad    %10.8f -> %10.8f fraction error = %10.6f%% error\n", _median, _mad,    Tmad,  100.0 * Tmad);
  writeLog("\n");
  writeLog("Using overlaps below %10.6f%% error\n", 100.0 * _errorLimit);
  writeLog("\n");
}



void
BestOverlapGraph::removeLopsidedEdges(const char *UNUSED(prefix)) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap *this5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap *this3 = getBestEdgeOverlap(fi, true);

    //  Ignore spurs and contains...and previously detected suspicious reads.  The suspicious reads
    //  do not have best edges back to them, and it's possible to find reads B where best edge A->B
    //  exists, yet no best edge from B exists.

    if ((isSuspicious(fi) == true) ||    //  Suspicious overlap pattern
        (isContained(fi)  == true) ||    //  Contained read (duh!)
        ((this5->readId() == 0) ||       //  Spur read
         (this3->readId() == 0)))
      continue;

    //  If there is a huge difference in error rates between the two best overlaps, that's a little
    //  suspicious.  This kind-of worked, but it is very sensitive to the 'limit', and was only
    //  tested on one bacteria.  It will also do very bad things in metagenomics.

#if 0
    double  this5erate  = this5->erate();
    double  this3erate  = this3->erate();

    double  limit       = 0.01;

    if (fabs(this5erate - this3erate) > limit) {
#pragma omp critical (suspInsert)
      {
        _suspicious.insert(fi);

        writeStatus("Incompatible error rates on best edges for read %u -- %.4f %.4f.\n", fi, this5erate, this3erate);

#warning NOT COUNTING ERATE DIFFS
        //_ERateIncompatible++;
      }
      continue;
    }
#endif

    //  Find the overlap for this5 and this3.

    int32   this5ovlLen = RI->overlapLength(fi, this5->readId(), this5->ahang(), this5->bhang());
    int32   this3ovlLen = RI->overlapLength(fi, this3->readId(), this3->ahang(), this3->bhang());

    //  Find the edges for our best overlaps.

    BestEdgeOverlap *that5 = getBestEdgeOverlap(this5->readId(), this5->read3p());
    BestEdgeOverlap *that3 = getBestEdgeOverlap(this3->readId(), this3->read3p());

    //  If both point back to us, we're done.  These must be symmetric, else overlapper is bonkers.

    if ((that5->readId() == fi) && (that5->read3p() == false) &&
        (that3->readId() == fi) && (that3->read3p() == true))
      continue;

    //  If there is an overlap to something with no overlaps out of it, that's a little suspicious.

    if ((that5->readId() == 0) ||
        (that3->readId() == 0)) {
      writeLog("WARNING: read %u has overlap to spur!  3' overlap to read %u back to read %u    5' overlap to read %u back to read %u\n",
               fi,
               this5->readId(), that5->readId(),
               this3->readId(), that3->readId());
#pragma omp critical (suspInsert)
      _suspicious.insert(fi);
      continue;
    }

    //  Something doesn't agree.  Find those overlaps...

    int32  that5ovlLen = RI->overlapLength(this5->readId(), that5->readId(), that5->ahang(), that5->bhang());
    int32  that3ovlLen = RI->overlapLength(this3->readId(), that3->readId(), that3->ahang(), that3->bhang());

    //  ...and compare.

    double  percDiff5 = 200.0 * abs(this5ovlLen - that5ovlLen) / (this5ovlLen + that5ovlLen);
    double  percDiff3 = 200.0 * abs(this3ovlLen - that3ovlLen) / (this3ovlLen + that3ovlLen);

    if ((percDiff5 <= 5.0) &&    //  Both good, keep 'em as is.
        (percDiff3 <= 5.0)) {
      //writeLog("fi %8u -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%% -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%% -- ACCEPTED\n",
      //         fi,
      //         this5->readId(), this5->read3p() ? '3' : '5', this5ovlLen, that5->readId(), that5->read3p() ? '3' : '5', that5ovlLen, percDiff5,
      //         this3->readId(), this3->read3p() ? '3' : '5', this3ovlLen, that3->readId(), that3->read3p() ? '3' : '5', that3ovlLen, percDiff3);
      continue;
    }

    //  Nope, one or both of the edges are too different.  Flag the read as suspicious.

    writeLog("lopsidedBest %8u -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%% -- %8u/%c' len %6u VS %8u/%c' len %6u %8.4f%%\n",
             fi,
             this5->readId(), this5->read3p() ? '3' : '5', this5ovlLen, that5->readId(), that5->read3p() ? '3' : '5', that5ovlLen, percDiff5,
             this3->readId(), this3->read3p() ? '3' : '5', this3ovlLen, that3->readId(), that3->read3p() ? '3' : '5', that3ovlLen, percDiff3);

#pragma omp critical (suspInsert)
    {
      _suspicious.insert(fi);

      if ((percDiff5 > 5.0) && (percDiff3 > 5.0))
        _n2EdgeIncompatible++;
      else
        _n1EdgeIncompatible++;
    }
  }
}



void
BestOverlapGraph::removeSimpleSpurs(const char *prefix) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  uint32  nSpur = 0;
  uint32  nSing = 0;

  FILE   *F = AS_UTL_openOutputFile(prefix, '.', "best.simpleSpurs");

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap  *edge5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap  *edge3 = getBestEdgeOverlap(fi, true);

    if (RI->readLength(fi) == 0)   //  Unused read, ignore.
      continue;

    //  Edge to a suspicious read makes us a spur!

    assert(isSuspicious(edge5->readId()) == false);
    assert(isSuspicious(edge3->readId()) == false);

    bool              spur5 = (edge5->readId() == 0);
    bool              spur3 = (edge3->readId() == 0);

    if (isContained(fi))      //  Contained, not a spur or singleton.
      continue;

    if ((spur5 == false) &&   //  Edges off of both ends.  Not a spur.
        (spur3 == false))
      continue;

    if ((spur5 == true) &&
        (spur3 == true)) {
      fprintf(F, F_U32" singleton\n", fi);
      nSing++;
    }

    else {
      fprintf(F, F_U32" %c' simle spur\n", fi, (spur5) ? '5' : '3');
      nSpur++;

      _spur.insert(fi);
    }
  }

  writeStatus("BestOverlapGraph()-- detected " F_SIZE_T " singleton reads.\n", nSing);
  writeStatus("BestOverlapGraph()-- detected " F_SIZE_T " spur reads.\n", nSpur);

  fclose(F);
}



uint32
BestOverlapGraph::spurDistance(BestEdgeOverlap *edge, uint32 limit, uint32 distance) {
  uint32  inid = edge->readId();
  bool    in3p = edge->read3p();
  bool    ot3p = (in3p == false) ? true : false;

  if (isSuspicious(inid) == true)   //  Edge to a suspicious read?
    return(distance);               //  Call it a spur.

  if (inid == 0)                    //  If no edge, we've terminated early,
    return(distance);               //  so we've followed a path to a spur.

  if (distance == limit)            //  If we've hit the limit, not a spur,
    return(UINT32_MAX);             //  return infinity.

  //  Otherwise, follow the path.

  return(spurDistance(getBestEdgeOverlap(inid, ot3p), limit, distance+1));
}



void
BestOverlapGraph::removeSpannedSpurs(const char *prefix, uint32 spurDepth) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  uint32  nSpur = 0;
  uint32  nSing = 0;

  FILE   *F = AS_UTL_openOutputFile(prefix, '.', "best.spannedSpurs");

  //  Find reads that terminate at a spur after some short traversal.

  set<uint32>  spurpath5;   //  spurpath is true if the edge out of this end
  set<uint32>  spurpath3;   //           leads to a dead-end spur.
  set<uint32>  spur;        //  spur     is true if this read is a spur.

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    if ((isContained(fi)  == true) ||   //  Contained read, ignore.
        (RI->readLength(fi) == 0))      //  Unused read, ignore.
      continue;

    uint32  dist5 = spurDistance(getBestEdgeOverlap(fi, false), spurDepth);
    uint32  dist3 = spurDistance(getBestEdgeOverlap(fi,  true), spurDepth);

    if ((dist5 == 0) ||
        (dist3 == 0)) {
      writeLog("SPUR read %u\n", fi);
      spur.insert(fi);
    }

    if (dist5 < spurDepth) {
      writeLog("SPUR path from read %u 5'\n", fi);
      spurpath5.insert(fi);
    }

    if (dist3 < spurDepth) {
      writeLog("SPUR path from read %u 3'\n", fi);
      spurpath3.insert(fi);
    }
  }

  writeStatus("Iter 0 - %d reads are spur.\n", spur.size());
  writeStatus("Iter 0 - %d 5' and %d 3' read ends leading to a spur.\n", spurpath5.size(), spurpath3.size());

  //  If any non-spur read has an edge to a spur flagged end,
  //  remove the flag on that end.

  for (uint32 iter=1; iter <= 2 * spurDepth; iter++) {
    for (uint32 fi=1; fi <= fiLimit; fi++) {
      if ((isContained(fi)  == true) ||   //  Skip contains (which should have no best edges anyway)
          (isSuspicious(fi) == true) ||   //  and suspicious (which DO have best edges).
          (isSpur(fi)       == true))
        continue;

      if (spur.count(fi) > 0)
        continue;

      BestEdgeOverlap  *edge5 = getBestEdgeOverlap(fi, false);
      BestEdgeOverlap  *edge3 = getBestEdgeOverlap(fi,  true);

#if 1

      //  If I'm not a spur, and there is a good path out of me, then the
      //  opposite-end edges out of this read can be used to remove the spur
      //  flag from the read ends they point to.

      if (spur.count(fi) == 0) {
        if (spurpath5.count(fi) == 0) {
          if (edge3->read3p() == false)   spurpath5.erase(edge3->readId());
          else                            spurpath3.erase(edge3->readId());
        }

        if (spurpath3.count(fi) == 0) {
          if (edge5->read3p() == false)   spurpath5.erase(edge5->readId());
          else                            spurpath3.erase(edge5->readId());
        }
      }

#else

      //  If I'm completely outside of a spurpath, then the edges out of this
      //  read can be used to remove the spur flag from the read ends they point
      //  to.

      if ((spur.count(fi)      == 0) &&
          (spurpath5.count(fi) == 0) &&
          (spurpath3.count(fi) == 0)) {
        if (edge5->read3p() == false)   spurpath5.erase(edge5->readId());
        else                            spurpath3.erase(edge5->readId());

        if (edge3->read3p() == false)   spurpath5.erase(edge3->readId());
        else                            spurpath3.erase(edge3->readId());
      }

#endif

    }

    //
    //  Compute new edges that don't point to a spur or spurpath, unless
    //  those are the only edges it has.
    //

    //  Reset just the scores.  Containment flags and best edges are not reset.
    memset(_scorA, 0, sizeof(BestScores)   * (fiLimit + 1));

    for (uint32 fi=1; fi <= fiLimit; fi++) {
      if ((isContained(fi)  == true) ||   //  Skip contains (which should have no best edges anyway)
          (isSuspicious(fi) == true) ||   //  and suspicious (which DO have best edges).
          (isSpur(fi)       == true))
        continue;

      uint32      no  = 0;
      BAToverlap *ovl = OC->getOverlaps(fi, no);

      BestEdgeOverlap  orig5 = *getBestEdgeOverlap(fi, false);
      BestEdgeOverlap  orig3 = *getBestEdgeOverlap(fi,  true);

      getBestEdgeOverlap(fi, false)->clear();
      getBestEdgeOverlap(fi,  true)->clear();

      //  Compute best edge, ignoring edges to spurs and to things that lead to spurs.

      for (uint32 ii=0; ii<no; ii++) {
        uint32  bid = ovl[ii].b_iid;

        if (ovl[ii].isDovetail() == false)
          continue;

        if (isOverlapBadQuality(ovl[ii]) == true)
          continue;

        if (isContained(ovl[ii].b_iid) == true)
          continue;

        if (isSuspicious(ovl[ii].b_iid) == true)
          continue;

        if (spur.count(ovl[ii].b_iid) != 0)
          continue;

        if (ovl[ii].AEndIs5prime()) {
          if (ovl[ii].BEndIs5prime() == true) {
            if (spurpath3.count(ovl[ii].b_iid) == 0)
              scoreEdge(ovl[ii]);
          } else {
            if (spurpath5.count(ovl[ii].b_iid) == 0)
              scoreEdge(ovl[ii]);
          }
        }

        if (ovl[ii].AEndIs3prime()) {
          if (ovl[ii].BEndIs5prime() == true) {
            if (spurpath3.count(ovl[ii].b_iid) == 0)
              scoreEdge(ovl[ii]);
          } else {
            if (spurpath5.count(ovl[ii].b_iid) == 0)
              scoreEdge(ovl[ii]);
          }
        }
      }

      //  If no edge found, reset to the previous edge.

      if (getBestEdgeOverlap(fi, false)->readId() == 0)   *getBestEdgeOverlap(fi, false) = orig5;
      if (getBestEdgeOverlap(fi,  true)->readId() == 0)   *getBestEdgeOverlap(fi,  true) = orig3;
    }

    writeStatus("Iter %u - %d 5' and %d 3' read ends leading to a spur.\n", iter, spurpath5.size(), spurpath3.size());
  }

  //  Now just flag any read with marks as a spur.

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    bool s5 = (spurpath5.count(fi) > 0);
    bool s3 = (spurpath3.count(fi) > 0);

    if ((isContained(fi)  == true) ||   //  Skip contains (which should have no best edges anyway)
        (isSuspicious(fi) == true) ||   //  and suspicious (which DO have best edges).
        (isSpur(fi)       == true))
      continue;

    if (spur.count(fi) > 0) {
      fprintf(F, "%u 5' is simple spur\n", fi);
      _spur.insert(fi);
      continue;
    }

    if ((s5 == true) && (s3 == true)) {
      fprintf(F, "%u    is non-spanned spur read\n", fi);
      _spur.insert(fi);
    }

    if ((s5 == true)  && (s3 == false)) {
      fprintf(F, "%u 5' is non-spanned spur end\n", fi);
      _spur.insert(fi);
    }

    if ((s5 == false) && (s3 == true))  {
      fprintf(F, "%u 3' is non-spanned spur end\n", fi);
      _spur.insert(fi);
    }
  }

  fclose(F);


  //  Remove edges that involve spurs.

  for (uint32 fi=1; fi <= fiLimit; fi++) {
    BestEdgeOverlap  *orig5 = getBestEdgeOverlap(fi, false);
    BestEdgeOverlap  *orig3 = getBestEdgeOverlap(fi,  true);

    //  Remove edges from spur to good reads.  Essentially, don't let spurs
    //  break contigs.

    if (_spur.count(fi) > 0) {
      if (_spur.count(orig5->readId()) == 0)   orig5->clear();
      if (_spur.count(orig3->readId()) == 0)   orig3->clear();

      orig5->clear();
      orig3->clear();
    }

    //  Remove edges from good to spurs.
#if 0
    if (_spur.count(orig5->readId()) > 0)
      orig5->clear();

    if (_spur.count(orig3->readId()) > 0)
      orig3->clear();
#endif
  }
}



//  Mark zombie masters.  Any read that has only contained overlaps (it is the container) and is the
//  smallest ID of those with no hangs, is a master.  These get promoted to unitigs.
//
void
BestOverlapGraph::findZombies(const char *prefix) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, no);
    uint32      nc = 0;

    if (no == 0)
      continue;

    for (uint32 ii=0; ii<no; ii++, nc++)       //  If any overlap makes A not
      if (ovl[ii].AisContainer() == false)     //  a container, it's not a zombie
        break;

    if (nc < no)
      continue;

    nc = UINT32_MAX;

    for (uint32 ii=0; ii<no; ii++)             //  Find the smallest ID
      if ((ovl[ii].a_hang == 0) &&             //  with no hangs.
          (ovl[ii].b_hang == 0) &&
          (ovl[ii].b_iid < nc))
        nc = ovl[ii].b_iid;

    if (fi < nc) {                             //  If we're smaller, we're a
#pragma omp critical (suspInsert)              //  Zombie Master!
      {
        writeLog("read %u is a zombie.\n", fi);
        _zombie.insert(fi);
      }
    }
  }

  writeStatus("BestOverlapGraph()-- detected " F_SIZE_T " zombie reads.\n", _zombie.size());
}



void
BestOverlapGraph::findEdges(void) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  memset(_bestA, 0, sizeof(BestOverlaps) * (fiLimit + 1));
  memset(_scorA, 0, sizeof(BestScores)   * (fiLimit + 1));

  //  One pass through all reads to flag any that are in a containment
  //  relationship.

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {
    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, no);

    for (uint32 ii=0; ii<no; ii++)
      markContainment(ovl[ii]);
  }

  //  A second pass to score and find the best edge for each end.
  //
  //  Contained reads can be skipped; we never use the edge here (in fact,
  //  it's just reset to zero at the end of graph construction).
  //
  //  Edges out of spurs are allowed, but edges into spurs are not.  This
  //  should prevent them from being incorporated into a promiscuous unitig,
  //  but still let them be popped as bubbles (but they shouldn't because
  //  they're spurs).

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi <= fiLimit; fi++) {

    if (isContained(fi) == true)
      continue;

    uint32      no  = 0;
    BAToverlap *ovl = OC->getOverlaps(fi, no);

    for (uint32 ii=0; ii<no; ii++) {
      uint32  bid = ovl[ii].b_iid;

      if ((isSpur(bid)       == false) &&    //  Edge into spur, ignore.
          (isSuspicious(bid) == false) &&    //  Edge into suspicious, ignore.
          (isContained(bid)  == false))
        scoreEdge(ovl[ii]);
    }
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



BestOverlapGraph::BestOverlapGraph(double        erateGraph,
                                   double        deviationGraph,
                                   const char   *prefix,
                                   bool          filterSuspicious,
                                   bool          filterHighError,
                                   bool          filterLopsided,
                                   bool          filterSpur,
                                   uint32        spurDepth) {

  writeStatus("\n");
  writeStatus("BestOverlapGraph()-- Computing Best Overlap Graph.\n");
  writeStatus("BestOverlapGraph()-- Allocating best edges (" F_SIZE_T "MB).\n",
           ((2 * sizeof(BestEdgeOverlap) * (RI->numReads() + 1)) >> 20));

  _bestA               = new BestOverlaps [RI->numReads() + 1];  //  Cleared in findEdges()
  _scorA               = new BestScores   [RI->numReads() + 1];

  _mean                = erateGraph;
  _stddev              = 0.0;

  _median              = erateGraph;
  _mad                 = 0.0;

  _errorLimit          = erateGraph;

  _n1EdgeFiltered      = 0;
  _n2EdgeFiltered      = 0;
  _n1EdgeIncompatible  = 0;
  _n2EdgeIncompatible  = 0;

  _suspicious.clear();

  _bestM.clear();
  _scorM.clear();

  _restrict            = NULL;
  _restrictEnabled     = false;

  _erateGraph          = erateGraph;
  _deviationGraph      = deviationGraph;

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

    removeHighErrorBestEdges();
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
    writeStatus("BestOverlapGraph()-- Filtering reads with suspicious overlap patterns.\n");

    removeSuspicious(prefix);
    findEdges();

    writeStatus("BestOverlapGraph()--   %u reads marked as suspicious.\n", _suspicious.size());
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

    writeStatus("BestOverlapGraph()--   %u reads marked as suspicious.\n", _suspicious.size());
  }

  else {
    writeStatus("BestOverlapGraph()-- NOT filtering reads with lopsided best edges.\n");
  }

  //
  //  Now that we know the max difference allowed in overlaps, mark reads as
  //  spurs if they're spurs.  Then don't find best edges to them.
  //

  if (filterSpur) {
    writeStatus("BestOverlapGraph()-- Filtering spur reads.\n");

    //  Simple spurs are essentially what was here before October 2019:
    //  flag any read with no overlaps on on end as a spur, and disallow
    //  best edges to them.
    //removeSimpleSpurs(prefix);
    //findEdges();

    //  Spanned spurs searches the graph for dead ends, then backs up
    //  to a previous branch.  It (supposedly) does this such that
    //  the spur path is preserved unless a branch is found.
    removeSpannedSpurs(prefix, spurDepth);

    //  This findEdges() should be redundant, and is probably harmful.
    //  It will reset our carefully adjusted branch edges back to paths
    //  into dead ends.
    //findEdges();
  }

  else {
    writeStatus("BestOverlapGraph()-- NOT filtering spur reads.\n");
  }

  //
  //  Do some final boring cleanup:
  //    flag zombie reads.
  //    remove best edges associated with contained reads.
  //

  findZombies(prefix);
  removeContainedDovetails();

  //
  //  Report filtering and final statistics.
  //

  writeLog("\n");
  writeLog("EDGE FILTERING\n");
  writeLog("-------- ------------------------------------------\n");
  writeLog("%8u reads have a suspicious overlap pattern\n", _suspicious.size());
  writeLog("%8u reads had edges filtered\n", _n1EdgeFiltered + _n2EdgeFiltered);
  writeLog("         %8u had one\n", _n1EdgeFiltered);
  writeLog("         %8u had two\n", _n2EdgeFiltered);
  writeLog("%8u reads have length incompatible edges\n", _n1EdgeIncompatible + _n2EdgeIncompatible);
  writeLog("         %8u have one\n", _n1EdgeIncompatible);
  writeLog("         %8u have two\n", _n2EdgeIncompatible);

  reportBestEdges(prefix, "best");

  reportEdgeStatistics(prefix, "FINAL");

  //  Done with scoring data.

  delete [] _scorA;
  _scorA = NULL;

  _spur.clear();

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
  FILE             *BCH = NULL;
  FILE *BE = NULL, *BEH = NULL, *BEG = NULL;
  FILE *BS = NULL;
  FILE *SS = NULL;

  //  Open output files.

  snprintf(N, FILENAME_MAX, "%s.%s.edges",               prefix, label);   BE = AS_UTL_openOutputFile(N);
  snprintf(N, FILENAME_MAX, "%s.%s.singletons",          prefix, label);   BS = AS_UTL_openOutputFile(N);
  snprintf(N, FILENAME_MAX, "%s.%s.edges.suspicious",    prefix, label);   SS = AS_UTL_openOutputFile(N);

  snprintf(N, FILENAME_MAX, "%s.%s.edges.gfa",           prefix, label);   BEG = AS_UTL_openOutputFile(N);

  snprintf(N, FILENAME_MAX, "%s.%s.contains.histogram",  prefix, label);   BCH = AS_UTL_openOutputFile(N);
  snprintf(N, FILENAME_MAX, "%s.%s.edges.histogram",     prefix, label);   BEH = AS_UTL_openOutputFile(N);

  //  Write best edges, singletons and suspicious edges.

  if ((BE) && (BS) && (SS)) {
    fprintf(BE, "#readId\tlibId\tbest5iid\tbest5end\tbest3iid\tbest3end\teRate5\teRate3\tbest5len\tbest3len\n");
    fprintf(BS, "#readId\tlibId\n");

    for (uint32 id=1; id<RI->numReads() + 1; id++) {
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      if ((bestedge5->readId() == 0) && (bestedge3->readId() == 0) && (isContained(id) == false)) {
        fprintf(BS, "%u\t%u\n", id, RI->libraryIID(id));
      }

      else if (isSuspicious(id) == true) {
        fprintf(SS, "%u\t%u\t%u\t%c'\t%u\t%c'\t%6.4f\t%6.4f\t%u\t%u%s\n", id, RI->libraryIID(id),
          bestedge5->readId(), bestedge5->read3p() ? '3' : '5',
                bestedge3->readId(), bestedge3->read3p() ? '3' : '5',
                AS_OVS_decodeEvalue(bestedge5->evalue()),
                AS_OVS_decodeEvalue(bestedge3->evalue()),
                (bestedge5->readId() == 0 ? 0 : RI->overlapLength(id, bestedge5->readId(), bestedge5->ahang(), bestedge5->bhang())),
                (bestedge3->readId() == 0 ? 0 : RI->overlapLength(id, bestedge3->readId(), bestedge3->ahang(), bestedge3->bhang())),
                isContained(id) ? "\tcontained" : "");
      }

      else {
        fprintf(BE, "%u\t%u\t%u\t%c'\t%u\t%c'\t%6.4f\t%6.4f\t%u\t%u%s\n", id, RI->libraryIID(id),
                bestedge5->readId(), bestedge5->read3p() ? '3' : '5',
                bestedge3->readId(), bestedge3->read3p() ? '3' : '5',
                AS_OVS_decodeEvalue(bestedge5->evalue()),
                AS_OVS_decodeEvalue(bestedge3->evalue()),
                (bestedge5->readId() == 0 ? 0 : RI->overlapLength(id, bestedge5->readId(), bestedge5->ahang(), bestedge5->bhang())),
                (bestedge3->readId() == 0 ? 0 : RI->overlapLength(id, bestedge3->readId(), bestedge3->ahang(), bestedge3->bhang())),
                isContained(id) ? "\tcontained" : "");
      }
    }
  }

  //  Write best edge graph.

  if (BEG) {
    fprintf(BEG, "H\tVN:Z:1.0\n");

    //  First, write the sequences used.  The sequence can be used as either
    //  a source node or a destination node (or both).

    set<uint32>  used;

    for (uint32 id=1; id<RI->numReads() + 1; id++) {
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      //  Ignore singletons, contained and suspicious reads.
      if (((bestedge5->readId() == 0) && (bestedge3->readId() == 0) && (isContained(id) == false)) ||
          (isContained(id) == true) ||
          (isSuspicious(id) == true))
        continue;

      //  Remember the source and destination of this edge.
      used.insert(id);
      used.insert(bestedge5->readId());
      used.insert(bestedge3->readId());
    }

    for (set<uint32>::iterator it=used.begin(); it != used.end(); it++)
      fprintf(BEG, "S\tread%08u\t*\tLN:i:%u\n", *it, RI->readLength(*it));

    //  Now, report edges.  GFA wants edges in exactly this format:
    //
    //       -------------
    //             -------------
    //
    //  with read orientation given by +/-.  Conveniently, this is what we've saved (for the edges).

    for (uint32 id=1; id<RI->numReads() + 1; id++) {
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      //  Ignore singletons, contained and suspicious reads.
      if (((bestedge5->readId() == 0) && (bestedge3->readId() == 0) && (isContained(id) == false)) ||
          (isContained(id) == true) ||
          (isSuspicious(id) == true))
        continue;

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

        fprintf(BEG, "L\tread%08u\t-\tread%08u\t%c\t%uM\n",
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

        fprintf(BEG, "L\tread%08u\t+\tread%08u\t%c\t%uM\n",
                id,
                bestedge3->readId(), bestedge3->read3p() ? '-' : '+',
                RI->overlapLength(id, bestedge3->readId(), bestedge3->ahang(), bestedge3->bhang()));
      }
    }
  }

  //  Write error rate histograms of best edges and contains.

  if ((BCH) && (BEH)) {
    double *bc = new double [RI->numReads() + 1 + RI->numReads() + 1];
    double *be = new double [RI->numReads() + 1 + RI->numReads() + 1];

    uint32  bcl = 0;
    uint32  bel = 0;

    for (uint32 id=1; id<RI->numReads() + 1; id++) {
      BestEdgeOverlap *bestedge5 = getBestEdgeOverlap(id, false);
      BestEdgeOverlap *bestedge3 = getBestEdgeOverlap(id, true);

      if (isContained(id)) {
        //bc[bcl++] = bestcont->erate();
#warning what is the error rate of the 'best contained' overlap?
        bc[bcl++] = bestedge5->erate();
        bc[bcl++] = bestedge3->erate();
      }
      else {
        if (bestedge5->readId() > 0)
          be[bel++] = bestedge5->erate();

        if (bestedge3->readId() > 0)
          be[bel++] = bestedge3->erate();
      }
    }

    sort(bc, bc+bcl);
    sort(be, be+bel);

    for (uint32 ii=0; ii<bcl; ii++)
      fprintf(BCH, "%f\n", bc[ii]);

    for (uint32 ii=0; ii<bel; ii++)
      fprintf(BEH, "%f\n", be[ii]);

    delete [] bc;
    delete [] be;
  }

  //  Close all the files.

  AS_UTL_closeFile(BE);
  AS_UTL_closeFile(BS);
  AS_UTL_closeFile(SS);

  AS_UTL_closeFile(BEG);

  AS_UTL_closeFile(BCH);
  AS_UTL_closeFile(BEH);
}



void
BestOverlapGraph::markContainment(BAToverlap& olap) {

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
    //  We only save if A is the contained read.
    return;

  setContained(olap.a_iid);
}



void
BestOverlapGraph::scoreEdge(BAToverlap& olap) {
  bool   enableLog = false;

  //if ((olap.a_iid == 97202) || (olap.a_iid == 30701))    //  Report logging only for
  //  enableLog = true;                                    //  specific annoying reads.

  //  We shold only be getting good clean dovetail overlaps here.  Except we don't.

  if (olap.isDovetail() == false) {
    //fprintf(stderr, "scoreEdge()-- OVERLAP CONT:     %6d %6d %c  hangs %6d %6d err %.3f -- not a dovetail\n",
    //        olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate());
    return;
  }

  assert(olap.isDovetail()        == true);

  assert(isSpur(olap.b_iid)       == false);
  assert(isSuspicious(olap.b_iid) == false);
  assert(isContained(olap.b_iid)  == false);

  //  Ignore the overlap if it's bad quality, or we aren't allowed to use it
  //  for this overlap graph.

  if (isOverlapBadQuality(olap)) {
    if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
      writeLog("scoreEdge()-- %6d %6d %c  hangs %6d %6d err %.3f -- bad quality\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate());
    return;
  }

  if (isOverlapRestricted(olap)) {
    if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
      writeLog("scoreEdge()-- %6d %6d %c  hangs %6d %6d err %.3f -- restricted\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate());
    return;
  }

  //  Compute the score for this overlap, and remember this overlap if the
  //  score is the best.

  uint64           newScr = scoreOverlap(olap);
  bool             a3p    = olap.AEndIs3prime();
  BestEdgeOverlap *best   = getBestEdgeOverlap(olap.a_iid, a3p);
  uint64          &score  = (a3p) ? (best3score(olap.a_iid)) : (best5score(olap.a_iid));

  assert(newScr > 0);

  if (newScr <= score) {
    if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
      writeLog("scoreEdge()-- %6d %6d %c  hangs %6d %6d err %.3f --      %c'\n",
               olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate(), (a3p) ? '3' : '5');
    return;
  }

  best->set(olap);

  score = newScr;

  if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
    writeLog("scoreEdge()-- %6d %6d %c  hangs %6d %6d err %.3f -- BEST %c'\n",
             olap.a_iid, olap.b_iid, olap.flipped ? 'A' : 'N', olap.a_hang, olap.b_hang, olap.erate(), (a3p) ? '3' : '5');
}



bool
BestOverlapGraph::isOverlapBadQuality(BAToverlap& olap) {
  bool   enableLog = false;

  //if ((olap.a_iid == 97202) || (olap.a_iid == 30701))    //  Report logging only for
  //  enableLog = true;                                    //  specific annoying reads.

  //  The overlap is bad if it involves deleted reads.  Shouldn't happen in a normal
  //  assembly, but sometimes us users want to delete reads after overlaps are generated.

  if ((RI->readLength(olap.a_iid) == 0) ||
      (RI->readLength(olap.b_iid) == 0)) {
    olap.filtered = true;
    return(true);
  }

  //  The overlap is GOOD (false == not bad) if the error rate is below the allowed erate.
  //  Initially, this is just the erate passed in.  After the first rount of finding edges,
  //  it is reset to the mean and stddev of selected best edges.

  if (olap.erate() <= _errorLimit) {
    if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
      writeLog("isOverlapBadQuality()-- %6d %6d %c  hangs %6d %6d err %.3f -- good quality\n",
               olap.a_iid, olap.b_iid,
               olap.flipped ? 'A' : 'N',
               olap.a_hang,
               olap.b_hang,
               olap.erate());

    return(false);
  }

  //  CA8.3 would further compute the (expected) number of errors in the alignment, and compare
  //  against another limit.  This was to allow very short overlaps where one error would push the
  //  error rate above a few percent.  canu doesn't do short overlaps.

  if ((enableLog == true) && (logFileFlagSet(LOG_OVERLAP_SCORING)))
    writeLog("isOverlapBadQuality()-- %6d %6d %c  hangs %6d %6d err %.3f -- REJECT!\n",
             olap.a_iid, olap.b_iid,
             olap.flipped ? 'A' : 'N',
             olap.a_hang,
             olap.b_hang,
             olap.erate());

  olap.filtered = true;
  return(true);
}



//  If no restrictions are known, this overlap is useful if both reads are not in a unitig
//  already.  Otherwise, we are restricted to just a specific set of reads (usually a single
//  unitig and all the mated reads).  The overlap is useful if both reads are in the set.
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

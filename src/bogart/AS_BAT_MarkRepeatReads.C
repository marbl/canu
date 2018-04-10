
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2016-MAR-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_AssemblyGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_MarkRepeatReads.H"

#include "intervalList.H"
#include "stddev.H"

#include <vector>

using namespace std;



//  Hack.
uint32 MIN_ANCHOR_HANG       = 500;   //  Require reads to be anchored by this many bases at boundaries of repeats.
int32  REPEAT_OVERLAP_MIN    = 50;

#define REPEAT_FRACTION   0.5

#undef  OLD_ANNOTATE
#undef  SHOW_ANNOTATE
#undef  SHOW_ANNOTATION_RAW             //  Show all overlaps used to annotate reads
#undef  SHOW_ANNOTATION_RAW_FILTERED    //  Show all overlaps filtered by high error rate

#undef  DUMP_READ_COVERAGE

//  Each evidence read picks its single best overlap to tig (based on overlaps to reads in the tig).
//  Filter out evidence that aligns at erate higher than expected.
//  Collapse to intervals on tig.
//  If still not significant and not spanned, break.



class olapDat {
public:
  olapDat(uint32 b, uint32 e, uint32 r, uint32 p) {
    tigbgn  = b;
    tigend  = e;
    eviRid  = r;
    eviPid  = p;
  };

  bool operator<(const olapDat &that)     const { return(tigbgn < that.tigbgn); };

  int32   tigbgn;   //  Location of the overlap on this tig
  int32   tigend;   //

  uint32  eviRid;   //  evidence read
  uint32  eviPid;   //  evidence read placeID
};



bool
olapDatByEviRid(const olapDat &A, const olapDat &B) {
  if (A.eviRid == B.eviRid)
    return(A.tigbgn < B.tigbgn);

  return(A.eviRid < B.eviRid);
}




class breakPointCoords {
public:
  breakPointCoords(int32 bgn, int32 end, bool rpt=false) {
    _bgn = bgn;
    _end = end;
    _rpt = rpt;
  };
  ~breakPointCoords() {
  };

  bool    operator<(breakPointCoords const &that) const {
    return(_bgn < that._bgn);
  };

  int32   _bgn;
  int32   _end;
  bool    _rpt;
};





//  Returns the coordinates the overlap intersects on the A read.
//
//        lo   hi
//        v    v
//  ------------
//        ----------
void
olapToReadCoords(ufNode *frg,
                 int32   ahang,  int32  bhang,
                 int32  &lo,     int32 &hi) {

  lo = 0;
  hi = RI->readLength(frg->ident);

  if (ahang > 0)
    lo += ahang;  //  Positive hang!

  if (bhang < 0)
    hi += bhang;  //  Negative hang!

  assert(0  <= lo);
  assert(0  <= hi);
  assert(lo <= hi);
  assert(lo <= RI->readLength(frg->ident));
  assert(hi <= RI->readLength(frg->ident));
}



void
findUnitigCoverage(Unitig               *tig,
                   intervalList<uint32> &coverage) {
  intervalList<uint32>  rawcoverage;

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode  frg = tig->ufpath[fi];

    if (frg.position.bgn < frg.position.end)
      rawcoverage.add(frg.position.bgn, frg.position.end - frg.position.bgn);
    else
      rawcoverage.add(frg.position.end, frg.position.bgn - frg.position.end);
  }

  coverage.clear();
  coverage.depth(rawcoverage);

#ifdef DUMP_READ_COVERAGE
  char  fn[FILENAME_MAX];
  snprintf(fn, FILENAME_MAX, "%08u.coverage", tig->id());
  FILE *F = AS_UTL_openOutputFile(fn);

  for (uint32 ii=0; ii<coverage.numberOfIntervals(); ii++)
    fprintf(F, "%u %u %u\n", coverage.lo(ii), coverage.hi(ii), coverage.depth(ii));

  AS_UTL_closeFile(F, fn);
#endif
}



uint32
splitTig(TigVector                &tigs,
         Unitig                   *tig,
         vector<breakPointCoords> &BP,
         Unitig                  **newTigs,
         int32                    *lowCoord,
         uint32                   *nRepeat,
         uint32                   *nUnique,
         bool                      doMove) {

  if (doMove == true) {
    memset(newTigs,  0, sizeof(Unitig *) * BP.size());
    memset(lowCoord, 0, sizeof(int32)    * BP.size());
  } else {
    memset(nRepeat,  0, sizeof(uint32)   * BP.size());
    memset(nUnique,  0, sizeof(uint32)   * BP.size());
  }

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode     &frg    = tig->ufpath[fi];
    int32       frgbgn = min(frg.position.bgn, frg.position.end);
    int32       frgend = max(frg.position.bgn, frg.position.end);

    //  Search for the region that matches the read.  BP's are sorted in increasing order.  It
    //  probably doesn't matter, but makes the logging a little easier to read.

    uint32      rid = UINT32_MAX;
    bool        rpt = false;

    //fprintf(stderr, "Searching for placement for read %u at %d-%d\n", frg.ident, frgbgn, frgend);

    for (uint32 ii=0; ii<BP.size(); ii++) {
      int32   rgnbgn = BP[ii]._bgn;
      int32   rgnend = BP[ii]._end;
      bool    repeat = BP[ii]._rpt;

      //  For repeats, the read must be contained fully.

      if ((repeat == true) && (rgnbgn <= frgbgn) && (frgend <= rgnend)) {
        rid = ii;
        rpt = true;
        break;
      }

      //  For non-repeat, the read just needs to intersect.

      if ((repeat == false) && (rgnbgn < frgend) && (frgbgn < rgnend)) {
        rid = ii;
        rpt = false;
        break;
      }
    }

    if (rid == UINT32_MAX) {
      fprintf(stderr, "Failed to place read %u at %d-%d\n", frg.ident, frgbgn, frgend);
      for (uint32 ii=0; ii<BP.size(); ii++)
        fprintf(stderr, "BP[%3u] at %8u-%8u repeat %u\n", ii, BP[ii]._bgn, BP[ii]._end, BP[ii]._rpt);
      flushLog();
    }
    assert(rid != UINT32_MAX);  //  We searched all the BP's, the read had better be placed!

    //  If moving reads, move the read!

    if (doMove) {
      if (newTigs[rid] == NULL) {
        lowCoord[rid] = frgbgn;

        newTigs[rid]  = tigs.newUnitig(true);

        if (nRepeat[rid] > nUnique[rid])
          newTigs[rid]->_isRepeat = true;
      }

      newTigs[rid]->addRead(frg, -lowCoord[rid], false);
    }

    //  Else, we're not moving, just count how many reads came from repeats or uniques.

    else {
      if (rpt)
        nRepeat[rid]++;
      else
        nUnique[rid]++;
    }
  }

  //  Return the number of tigs created.

  uint32  nTigsCreated = 0;

  for (uint32 ii=0; ii<BP.size(); ii++)
    if (nRepeat[ii] + nUnique[ii] > 0)
      nTigsCreated++;

  return(nTigsCreated);
}




//  Over all reads in tgA, return a vector of olapDat (tigBgn, tigEnd, eviRid)
//  for all reads that overlap into this tig.
//
//  The current AssemblyGraph is backwards to what we need.  It has, for each read, the
//  overlaps from that read that are compatible - but we need to the overlaps to each
//  read that are compatible, and the two are not symmetric.  A can be compatible in tig 1,
//  but the same overlapping read B can be incompatible with tig 2.
//
//  We can invert the graph at the start of repeat detection, making a list of
//    read B ---> overlaps to tig N position X-Y, with read A


void
annotateRepeatsOnRead(AssemblyGraph         *AG,
                      TigVector             &UNUSED(tigs),
                      Unitig                *tig,
                      double                 UNUSED(deviationRepeat),
                      vector<olapDat>       &repeats) {

  //  Over all reads in this tig,
  //  Grab pointers to all incoming edges.
  //  Push those locations onto our output list.

  for (uint32 ii=0; ii<tig->ufpath.size(); ii++) {
    ufNode               *read   = &tig->ufpath[ii];
    vector<BestReverse>  &rPlace = AG->getReverse(read->ident);

#if 0
    writeLog("annotateRepeatsOnRead()-- tig %u read #%u %u at %d-%d reverse %u items\n",
             tig->id(), ii, read->ident,
             read->position.bgn,
             read->position.end,
             rPlace.size());
#endif

    for (uint32 rr=0; rr<rPlace.size(); rr++) {
      uint32          rID    = rPlace[rr].readID;
      uint32          pID    = rPlace[rr].placeID;
      BestPlacement  &fPlace = AG->getForward(rID)[pID];

#ifdef SHOW_ANNOTATION_RAW
      writeLog("annotateRepeatsOnRead()-- tig %u read #%u %u place %u reverse read %u in tig %u placed %d-%d olap %d-%d%s\n",
               tig->id(), ii, read->ident, rr,
               rID,
               tig->inUnitig(rID),
               fPlace.placedBgn, fPlace.placedEnd,
               fPlace.olapBgn,   fPlace.olapEnd,
               (fPlace.isUnitig) ? " IN_UNITIG" : "");
#endif

      if ((fPlace.isUnitig == true) ||
          (fPlace.isContig == true))
        continue;

      repeats.push_back(olapDat(fPlace.olapBgn, fPlace.olapEnd, rID, pID));
    }
  }
}



void
mergeAnnotations(vector<olapDat> &repeatOlaps) {
  sort(repeatOlaps.begin(), repeatOlaps.end(), olapDatByEviRid);

#ifdef SHOW_ANNOTATE
  for (uint32 ii=0; ii<repeatOlaps.size(); ii++)
    if (repeatOlaps[ii].tigbgn < 1000000)
      writeLog("repeatOlaps[%u] %d-%d from read %u place %u RAW\n",
               ii,
               repeatOlaps[ii].tigbgn, repeatOlaps[ii].tigend,
               repeatOlaps[ii].eviRid, repeatOlaps[ii].eviPid);

  flushLog();
#endif

  for (uint32 dd=0, ss=1; ss<repeatOlaps.size(); ss++) {
    assert(repeatOlaps[dd].eviRid <= repeatOlaps[ss].eviRid);

    //  If different evidence reads, close the destination olap, set up
    //  for a new destination.

    if (repeatOlaps[dd].eviRid != repeatOlaps[ss].eviRid) {
      dd = ss;
      continue;
    }

    //  If the destination ends before the source begins, there is no overlap between the
    //  two regions.  Close dd, set up for a new dd.

    if (repeatOlaps[dd].tigend <= repeatOlaps[ss].tigbgn) {
      dd = ss;
      continue;
    }

    //  Otherwise, there must be an overlap.  Extend the destination region, erase the source
    //  region.

    repeatOlaps[dd].tigbgn = min(repeatOlaps[ss].tigbgn, repeatOlaps[dd].tigbgn);
    repeatOlaps[dd].tigend = max(repeatOlaps[ss].tigend, repeatOlaps[dd].tigend);

    repeatOlaps[ss].tigbgn = UINT32_MAX;
    repeatOlaps[ss].tigend = UINT32_MAX;
    repeatOlaps[ss].eviRid = UINT32_MAX;
    repeatOlaps[ss].eviPid = UINT32_MAX;
  }

  //  Sort overlaps again.  This pushes all those 'erased' regions to the end of the list, which
  //  we can then just pop off.

  sort(repeatOlaps.begin(), repeatOlaps.end(), olapDatByEviRid);

  for (uint32 ii=repeatOlaps.size(); ii--; )
    if (repeatOlaps[ii].eviRid == UINT32_MAX)
      repeatOlaps.pop_back();

  //  For logging, sort by coordinate

  sort(repeatOlaps.begin(), repeatOlaps.end());

#ifdef SHOW_ANNOTATE
  for (uint32 ii=0; ii<repeatOlaps.size(); ii++)
    if (repeatOlaps[ii].tigbgn < 1000000)
      writeLog("repeatOlaps[%d] %d-%d from tig %u read %u place %u MERGED\n",
               ii,
               repeatOlaps[ii].tigbgn, repeatOlaps[ii].tigend,
               repeatOlaps[ii].eviRid, repeatOlaps[ii].eviPid);
#endif
}



void
discardSpannedRepeats(Unitig              *tig,
                      intervalList<int32> &tigMarksR) {

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode     *frg       = &tig->ufpath[fi];
    bool        frgfwd    = (frg->position.bgn < frg->position.end);
    int32       frglo     = (frgfwd) ? frg->position.bgn : frg->position.end;
    int32       frghi     = (frgfwd) ? frg->position.end : frg->position.bgn;
    bool        discarded = false;

    for (uint32 ri=0; ri<tigMarksR.numberOfIntervals(); ri++) {
      bool   spanLo = false;
      bool   spanHi = false;

      //  The decision of 'spanned by a read' is broken into two pieces: does the read span the
      //  lower (higher) boundary of the region.  To be spanned, the boundary needs to be spanned
      //  by at least MIN_ANCHOR_HANG additional bases (to anchor the read to non-repeat
      //  sequence).
      //
      //  This is a problem at the start/end of the tig, beacuse no read will extend past the
      //  start/end of the tig.  Instead, if the repeat is contained within the first (last) read
      //  with no extension at the respective end, it is spanned.

      if ((frglo == 0) &&                                   //  Read at start of tig, spans off the high end
          (tigMarksR.hi(ri) + MIN_ANCHOR_HANG <= frghi))
        spanLo = spanHi = true;

      if ((frghi == tig->getLength()) &&                    //  Read at end of tig, spans off the low end
          (frglo + MIN_ANCHOR_HANG <= tigMarksR.lo(ri)))
        spanLo = spanHi = true;

      if (frglo + MIN_ANCHOR_HANG <= tigMarksR.lo(ri))      //  Read spanned off the low end
        spanLo = true;

      if (tigMarksR.hi(ri) + MIN_ANCHOR_HANG <= frghi)      //  Read spanned off the high end
        spanHi = true;

      if (spanLo && spanHi) {
        writeLog("discard region %8d:%-8d - contained in read %6u %8d-%8d\n",
                 tigMarksR.lo(ri), tigMarksR.hi(ri), frg->ident, frglo, frghi);

        tigMarksR.lo(ri) = 0;
        tigMarksR.hi(ri) = 0;

        discarded = true;
      }
    }

    if (discarded)
      tigMarksR.filterShort(1);
  }
}



void
reportThickestEdgesInRepeats(Unitig               *tig,
                             intervalList<int32>  &tigMarksR) {

  writeLog("thickest edges to the repeat regions:\n");

  for (uint32 ri=0; ri<tigMarksR.numberOfIntervals(); ri++) {
    uint32   t5 = UINT32_MAX, l5 = 0, t5bgn = 0, t5end = 0;
    uint32   t3 = UINT32_MAX, l3 = 0, t3bgn = 0, t3end = 0;

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode     *frg       = &tig->ufpath[fi];
      bool        frgfwd    = (frg->position.bgn < frg->position.end);
      int32       frglo     = (frgfwd) ? frg->position.bgn : frg->position.end;
      int32       frghi     = (frgfwd) ? frg->position.end : frg->position.bgn;
      bool        discarded = false;

      //  Overlap off the 5' end of the region.
      if (frglo <= tigMarksR.lo(ri) && (tigMarksR.lo(ri) <= frghi)) {
        uint32 olap = frghi - tigMarksR.lo(ri);
        if (l5 < olap) {
          l5    = olap;
          t5    = fi;
          t5bgn = frglo;  //  Easier than recomputing it later on...
          t5end = frghi;
        }
      }

      //  Overlap off the 3' end of the region.
      if (frglo <= tigMarksR.hi(ri) && (tigMarksR.hi(ri) <= frghi)) {
        uint32 olap = tigMarksR.hi(ri) - frglo;
        if (l3 < olap) {
          l3    = olap;
          t3    = fi;
          t3bgn = frglo;
          t3end = frghi;
        }
      }

      if (frglo <= tigMarksR.lo(ri) && (tigMarksR.hi(ri) <= frghi)) {
        writeLog("saved   region %8d:%-8d - closest    read %6u (%+6d) %8d:%-8d (%+6d) (contained)\n",
                 tigMarksR.lo(ri), tigMarksR.hi(ri),
                 frg->ident,
                 tigMarksR.lo(ri) - frglo, frglo,
                 frghi, frghi - tigMarksR.hi(ri));
      }
    }

    if (t5 != UINT32_MAX)
      writeLog("saved   region %8d:%-8d - closest 5' read %6u (%+6d) %8d:%-8d (%+6d)\n",
               tigMarksR.lo(ri), tigMarksR.hi(ri),
               tig->ufpath[t5].ident,
               tigMarksR.lo(ri) - t5bgn, t5bgn,
               t5end, t5end - tigMarksR.hi(ri));

    if (t3 != UINT32_MAX)
      writeLog("saved   region %8d:%-8d - closest 3' read %6u (%+6d) %8d:%-8d (%+6d)\n",
               tigMarksR.lo(ri), tigMarksR.hi(ri),
               tig->ufpath[t3].ident,
               tigMarksR.lo(ri) - t3bgn, t3bgn,
               t3end, t3end - tigMarksR.hi(ri));
  }
}



uint32 *
findConfusedEdges(TigVector            &tigs,
                  Unitig                *tig,
                  intervalList<int32>  &tigMarksR,
                  double                confusedAbsolute,
                  double                confusedPercent,
                  vector<confusedEdge> &confusedEdges) {

  uint32  *isConfused  = new uint32 [tigMarksR.numberOfIntervals()];

  memset(isConfused, 0, sizeof(uint32) * tigMarksR.numberOfIntervals());

  //  Examine every read in this tig.  If the read intersects a marked repeat, find the best edge
  //  that continues the tig in either direction.  If those reads are in the repeat region, scan all
  //  the overlaps of this read for any that are of comparable length.  If any are found, declare
  //  this repeat to be potentially confused.  If none are found - for the whole repeat region -
  //  then we can leave the repeat alone.

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode     *rdA       = &tig->ufpath[fi];
    uint32      rdAid     = rdA->ident;
    bool        rdAfwd    = (rdA->position.bgn < rdA->position.end);
    int32       rdAlo     = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
    int32       rdAhi     = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

    double      sc        = (rdAhi - rdAlo) / (double)RI->readLength(rdAid);

    if ((OG->isContained(rdAid)  == true) ||    //  Don't care about contained or suspicious
        (OG->isSuspicious(rdAid) == true))      //  reads; we'll use the container instead.
      continue;

    for (uint32 ri=0; ri<tigMarksR.numberOfIntervals(); ri++) {
      uint32  rMin = tigMarksR.lo(ri);
      uint32  rMax = tigMarksR.hi(ri);

      if ((rdAhi < rMin) ||   //  Read ends before the region
          (rMax  < rdAlo))    //  Read starts after the region
        continue;              //   -> don't care about this read!

      //  Compute the position (in the tig) of the best overlaps.

      int32  tig5bgn=0, tig5end=0;
      int32  tig3bgn=0, tig3end=0;

      //  Instead of using the best edge - which might not be the edge used in the unitig -
      //  we need to scan the layout to return the previous/next dovetail

      //  Put this in a function - what to return if no best overlap?

      BestEdgeOverlap   *b5 = OG->getBestEdgeOverlap(rdAid, false);
      BestEdgeOverlap   *b3 = OG->getBestEdgeOverlap(rdAid, true);

      //  If the best edge is to a read not in this tig, there is nothing to compare against.
      //  Is this confused by default?  Possibly.  The unitig was constructed somehow, and that
      //  must then be the edge coming into us.  We'll pick it up later.

      bool b5use = true;
      bool b3use = true;

      if (b5->readId() == 0)
        b5use = false;
      if (b3->readId() == 0)
        b3use = false;

      if ((b5use) && (tig->inUnitig(b5->readId()) != tig->id()))
        b5use = false;
      if ((b3use) && (tig->inUnitig(b3->readId()) != tig->id()))
        b3use = false;

      //  The best edge read is in this tig.  If they don't overlap, again, nothing to compare
      //  against.

      if (b5use) {
        ufNode     *rdB       = &tig->ufpath[tig->ufpathIdx(b5->readId())];
        uint32      rdBid     = rdB->ident;
        bool        rdBfwd    = (rdB->position.bgn < rdB->position.end);
        int32       rdBlo     = (rdBfwd) ? rdB->position.bgn : rdB->position.end;
        int32       rdBhi     = (rdBfwd) ? rdB->position.end : rdB->position.bgn;

        if ((rdAhi < rdBlo) ||
            (rdBhi < rdAlo))
          b5use = false;
      }

      if (b3use) {
        ufNode     *rdB       = &tig->ufpath[tig->ufpathIdx(b3->readId())];
        uint32      rdBid     = rdB->ident;
        bool        rdBfwd    = (rdB->position.bgn < rdB->position.end);
        int32       rdBlo     = (rdBfwd) ? rdB->position.bgn : rdB->position.end;
        int32       rdBhi     = (rdBfwd) ? rdB->position.end : rdB->position.bgn;

        if ((rdAhi < rdBlo) ||
            (rdBhi < rdAlo))
          b3use = false;
      }

      //  If we can use this edge, compute the placement of the overlap on the unitig.

      //  Call #1;

      if (b5use) {
        int32   bgn=0, end=0;

        olapToReadCoords(rdA,
                         b5->ahang(), b5->bhang(),
                         bgn, end);

        tig5bgn = (rdAfwd) ? (rdAlo + sc * bgn) : (rdAhi - sc * end);
        tig5end = (rdAfwd) ? (rdAlo + sc * end) : (rdAhi - sc * bgn);

        assert(tig5bgn < tig5end);

        if (tig5bgn < 0)                  tig5bgn = 0;
        if (tig5end > tig->getLength())   tig5end = tig->getLength();
      }

      //  Call #2

      if (b3use) {
        int32   bgn=0, end=0;

        olapToReadCoords(rdA,
                         b3->ahang(), b3->bhang(),
                         bgn, end);

        tig3bgn = (rdAfwd) ? (rdAlo + sc * bgn) : (rdAhi - sc * end);
        tig3end = (rdAfwd) ? (rdAlo + sc * end) : (rdAhi - sc * bgn);

        assert(tig3bgn < tig3end);

        if (tig3bgn < 0)                  tig3bgn = 0;
        if (tig3end > tig->getLength())   tig3end = tig->getLength();
      }

      //  If either of the 5' or 3' overlaps (or both!) are in the repeat region, we need to check for
      //  close overlaps on that end.

      uint32  len5 = 0;
      uint32  len3 = 0;

      if ((rMin    < tig5bgn) &&
          (tig5end < rMax) &&
          (b5use))
        len5 = RI->overlapLength(rdAid, b5->readId(), b5->ahang(), b5->bhang());
      else
        b5use = false;

      if ((rMin    < tig3bgn) &&
          (tig3end < rMax) &&
          (b3use))
        len3 = RI->overlapLength(rdAid, b3->readId(), b3->ahang(), b3->bhang());
      else
        b3use = false;

      double score5 = len5 * (1 - b5->erate());
      double score3 = len3 * (1 - b3->erate());

      //  Neither of the best edges are in the repeat region; move to the next region and/or read.
      if (len5 + len3 == 0)
        continue;

      //  At least one of the best edge overlaps is in the repeat region.  Scan for other edges
      //  that are of comparable length and quality.

      uint32        ovlLen   = 0;
      BAToverlap   *ovl      = OC->getOverlaps(rdAid, ovlLen);

      for (uint32 oo=0; oo<ovlLen; oo++) {
        uint32   rdBid    = ovl[oo].b_iid;
        uint32   tgBid    = tigs.inUnitig(rdBid);

        //  If the read is in a singleton, skip.  These are unassembled crud.
        if ((tgBid                      == 0) ||
            (tigs[tgBid]                == NULL) ||
            (tigs[tgBid]->ufpath.size() == 1))
          continue;

        //  Skip if this overlap is the best we're trying to match.
        //
        //  NOTE.  This doesn't care about potential duplicate overlaps between a pair of reads,
        //  as we're looking for thicker overlaps off only one end of the read.
        if ((rdBid == b5->readId()) ||
            (rdBid == b3->readId()))
          continue;

        //  Skip if this overlap is crappy quality
        if (OG->isOverlapBadQuality(ovl[oo]))
          continue;

        //  Skip if the read is contained or suspicious.
        if ((OG->isContained(rdBid)  == true) ||
            (OG->isSuspicious(rdBid) == true))
          continue;

        //  Skip if the overlap isn't dovetail.
        bool  ovl5 = ovl[oo].AEndIs5prime();
        bool  ovl3 = ovl[oo].AEndIs3prime();

        if ((ovl5 == false) &&
            (ovl3 == false))
          continue;

        //  Skip if we're not using this overlap
        if ((ovl5 == true) && (b5use == false))
          continue;

        if ((ovl3 == true) && (b3use == false))
          continue;


        uint32   rdBpos   =  tigs[tgBid]->ufpathIdx(rdBid);
        ufNode  *rdB      = &tigs[tgBid]->ufpath[rdBpos];

        bool     rdBfwd   = (rdB->position.bgn < rdB->position.end);
        int32    rdBlo    = (rdBfwd) ? rdB->position.bgn : rdB->position.end;
        int32    rdBhi    = (rdBfwd) ? rdB->position.end : rdB->position.bgn;

        //  If the overlap is to a read in a different tig, or
        //     the overlap is to a read in the same tig, but we don't overlap in the tig, check lengths.
        //  Otherwise, the overlap is present in the tig, and can't be confused.
        if ((tgBid == tig->id()) &&
            (rdBlo <= rdAhi) &&
            (rdAlo <= rdBhi))
          continue;

        uint32  len   = RI->overlapLength(rdAid, ovl[oo].b_iid, ovl[oo].a_hang, ovl[oo].b_hang);
        double  score = len * (1 - ovl[oo].erate());

        //  Compute percent difference.

        double  ad5 = fabs(score - score5);
        double  ad3 = fabs(score - score3);

        double  pd5 = 200 * ad5 / (score + score5);
        double  pd3 = 200 * ad3 / (score + score3);

        //  Skip if this overlap is vastly worse than the best.

        if ((ovl5 == true) && ((ad5 >= confusedAbsolute) || (pd5 > confusedPercent))) {
          writeLog("tig %7u read %8u pos %7u-%-7u NOT confused by 5' edge to read %8u - best edge read %8u len %6u erate %.4f score %8.2f - alt edge len %6u erate %.4f score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                   tig->id(), rdAid, rdAlo, rdAhi,
                   rdBid,
                   b5->readId(), len5, b5->erate(), score5,
                   len, ovl[oo].erate(), score,
                   ad5, pd5);
          continue;
        }

        if ((ovl3 == true) && ((ad3 >= confusedAbsolute) || (pd3 > confusedPercent))) {
          writeLog("tig %7u read %8u pos %7u-%-7u NOT confused by 3' edge to read %8u - best edge read %8u len %6u erate %.4f score %8.2f - alt edge len %6u erate %.4f score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                   tig->id(), rdAid, rdAlo, rdAhi,
                   rdBid,
                   b3->readId(), len3, b3->erate(), score3,
                   len, ovl[oo].erate(), score,
                   ad3, pd3);
          continue;
        }

        //  Potential confusion!

        if (ovl5 == true) {
          writeLog("tig %7u read %8u pos %7u-%-7u IS confused by 5' edge to read %8u - best edge read %8u len %6u erate %.4f score %8.2f - alt edge len %6u erate %.4f score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                   tig->id(), rdAid, rdAlo, rdAhi,
                   rdBid,
                   b5->readId(), len5, b5->erate(), score5,
                   len, ovl[oo].erate(), score,
                   ad5, pd5);

          confusedEdges.push_back(confusedEdge(rdAid, false, rdBid));
        }

        if (ovl3 == true) {
          writeLog("tig %7u read %8u pos %7u-%-7u IS confused by 3' edge to read %8u - best edge read %8u len %6u erate %.4f score %8.2f - alt edge len %6u erate %.4f score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                   tig->id(), rdAid, rdAlo, rdAhi,
                   rdBid,
                   b3->readId(), len3, b3->erate(), score3,
                   len, ovl[oo].erate(), score,
                   ad3, pd3);

          confusedEdges.push_back(confusedEdge(rdAid, true, rdBid));
        }

        isConfused[ri]++;
      }
    }  //  Over all marks (ri)
  }  //  Over all reads (fi)

  return(isConfused);
}



void
discardUnambiguousRepeats(TigVector            &tigs,
                          Unitig                *tig,
                          intervalList<int32>  &tigMarksR,
                          double                confusedAbsolute,
                          double                confusedPercent,
                          vector<confusedEdge> &confusedEdges) {

  uint32  *isConfused = findConfusedEdges(tigs, tig, tigMarksR, confusedAbsolute, confusedPercent, confusedEdges);

  //  Scan all the regions, and delete any that have no confusion.

  bool  discarded = false;

  for (uint32 ri=0; ri<tigMarksR.numberOfIntervals(); ri++) {
    if (isConfused[ri] == 0) {
      writeLog("discard region %8d:%-8d - no confusion in best edges\n",
               tigMarksR.lo(ri), tigMarksR.hi(ri));

      tigMarksR.lo(ri) = 0;
      tigMarksR.hi(ri) = 0;

      discarded = true;
    }

    else {
      writeLog("saved   region %8d:%-8d - %u best edges are potentially confused\n",
               tigMarksR.lo(ri), tigMarksR.hi(ri), isConfused[ri]);
    }
  }

  if (discarded)
    tigMarksR.filterShort(1);

  delete [] isConfused;
}



void
mergeAdjacentRegions(Unitig                *tig,
                     intervalList<int32>  &tigMarksR) {

  //  Extend, but don't extend past the end of the tig.

  for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++) {
    tigMarksR.lo(ii) = max<int32>(tigMarksR.lo(ii) - MIN_ANCHOR_HANG, 0);
    tigMarksR.hi(ii) = min<int32>(tigMarksR.hi(ii) + MIN_ANCHOR_HANG, tig->getLength());
  }

  //  Merge.

  bool  merged = false;

  for (uint32 ri=1; ri<tigMarksR.numberOfIntervals(); ri++) {
    uint32  rMin = min(tigMarksR.hi(ri-1), tigMarksR.lo(ri));
    uint32  rMax = max(tigMarksR.hi(ri-1), tigMarksR.lo(ri));

    if (tigMarksR.lo(ri) <= tigMarksR.hi(ri-1)) {
      writeLog("merge extended regions %8d:%-8d and %8d:%-8d\n",
               tigMarksR.lo(ri-1), tigMarksR.hi(ri-1),
               tigMarksR.lo(ri), tigMarksR.hi(ri));

      tigMarksR.lo(ri) = tigMarksR.lo(ri-1);

      tigMarksR.lo(ri-1) = 0;   //  CRITICAL to delete the ri-1 interval (and not ri) because the next
      tigMarksR.hi(ri-1) = 0;   //  iteration will be using ri (as its ri-1).  ri-1 here is never seen again.

      merged = true;
    }
  }

  if (merged)
    tigMarksR.filterShort(1);
}



void
reportTigsCreated(Unitig                    *tig,
                  vector<breakPointCoords>  &BP,
                  uint32                     nTigs,
                  Unitig                   **newTigs,
                  uint32                    *nRepeat,
                  uint32                    *nUnique) {

  for (uint32 ii=0; ii<BP.size(); ii++) {
    int32   rgnbgn = BP[ii]._bgn;
    int32   rgnend = BP[ii]._end;
    bool    repeat = BP[ii]._rpt;

    if      (nRepeat[ii] + nUnique[ii] == 0)
      writeLog("For tig %5u %s region %8d %8d - %6u/%6u repeat/unique reads - no new unitig created.\n",
               tig->id(), (repeat == true) ? "repeat" : "unique", rgnbgn, rgnend, nRepeat[ii], nUnique[ii]);

    else if (nTigs > 1)
      writeLog("For tig %5u %s region %8d %8d - %6u/%6u reads repeat/unique - unitig %5u created.\n",
               tig->id(), (repeat == true) ? "repeat" : "unique", rgnbgn, rgnend, nRepeat[ii], nUnique[ii], newTigs[ii]->id());

    else
      writeLog("For tig %5u %s region %8d %8d - %6u/%6u repeat/unique reads - unitig %5u remains unchanged.\n",
               tig->id(), (repeat == true) ? "repeat" : "unique", rgnbgn, rgnend, nRepeat[ii], nUnique[ii], tig->id());
  }
}





void
markRepeatReads(AssemblyGraph         *AG,
                TigVector             &tigs,
                double                 deviationRepeat,
                uint32                 confusedAbsolute,
                double                 confusedPercent,
                vector<confusedEdge>  &confusedEdges) {
  uint32  tiLimit = tigs.size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  writeLog("repeatDetect()-- working on " F_U32 " tigs, with " F_U32 " thread%s.\n", tiLimit, numThreads, (numThreads == 1) ? "" : "s");

  vector<olapDat>      repeatOlaps;   //  Overlaps to reads promoted to tig coords

  intervalList<int32>  tigMarksR;     //  Marked repeats based on reads, filtered by spanning reads
  intervalList<int32>  tigMarksU;     //  Non-repeat invervals, just the inversion of tigMarksR


  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = tigs[ti];

    if ((tig == NULL) ||                  //  Deleted, nothing to do.
        (tig->ufpath.size() == 1) ||      //  Singleton, nothing to do.
        (tig->_isUnassembled == true))    //  Unassembled, don't care.
      continue;

    writeLog("Annotating repeats in reads for tig %u/%u.\n", ti, tiLimit);

    //  Clear out all the existing marks.  They're not for this tig.


    //  Analyze overlaps for each read.  For each overlap to a read not in this tig, or not
    //  overlapping in this tig, and of acceptable error rate, add the overlap to repeatOlaps.

    repeatOlaps.clear();

    uint32  fiLimit    = tig->ufpath.size();
    uint32  numThreads = omp_get_max_threads();
    uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

    annotateRepeatsOnRead(AG, tigs, tig, deviationRepeat, repeatOlaps);

    writeLog("Annotated with %lu overlaps.\n", repeatOlaps.size());

    //  Merge marks for the same read into the largest possible.

    mergeAnnotations(repeatOlaps);

    //  Make a new set of intervals based on all the detected repeats.

    tigMarksR.clear();

    for (uint32 bb=0, ii=0; ii<repeatOlaps.size(); ii++)
      tigMarksR.add(repeatOlaps[ii].tigbgn, repeatOlaps[ii].tigend - repeatOlaps[ii].tigbgn);

    //  Collapse these markings Collapse all the read markings to intervals on the unitig, merging those that overlap
    //  significantly.

    tigMarksR.merge(REPEAT_OVERLAP_MIN);

    //  Scan reads, discard any mark that is contained in a read
    //
    //  We don't need to filterShort() after every one is removed, but it's simpler to do it Right Now than
    //  to track if it is needed.

    writeLog("Scan reads to discard spanned repeats.\n");

    discardSpannedRepeats(tig, tigMarksR);

    //  Run through again, looking for the thickest overlap(s) to the remaining regions.
    //  This isn't caring about the end effect noted above.

    reportThickestEdgesInRepeats(tig, tigMarksR);

    //  Scan reads.  If a read intersects a repeat interval, and the best edge for that read
    //  is entirely in the repeat region, decide if there is a near-best edge to something
    //  not in this tig.
    //
    //  A region with no such near-best edges is _probably_ correct.

    writeLog("search for confused edges:\n");

    discardUnambiguousRepeats(tigs, tig, tigMarksR, confusedAbsolute, confusedPercent, confusedEdges);


    //  Merge adjacent repeats.
    //
    //  When we split (later), we require a MIN_ANCHOR_HANG overlap to anchor a read in a unique
    //  region.  This is accomplished by extending the repeat regions on both ends.  For regions
    //  close together, this could leave a negative length unique region between them:
    //
    //   ---[-----]--[-----]---  before
    //   -[--------[]--------]-  after extending by MIN_ANCHOR_HANG (== two dashes)
    //
    //  To solve this, regions that were linked together by a single read (with sufficient overlaps
    //  to each) were merged.  However, there was no maximum imposed on the distance between the
    //  repeats, so (in theory) a 150kbp read could attach two repeats to a 149kbp unique unitig --
    //  and label that as a repeat.  After the merges were completed, the regions were extended.
    //
    //  This version will extend regions first, then merge repeats only if they intersect.  No need
    //  for a linking read.
    //
    //  The extension also serves to clean up the edges of tigs, where the repeat doesn't quite
    //  extend to the end of the tig, leaving a few hundred bases of non-repeat.

    mergeAdjacentRegions(tig, tigMarksR);


    //  Invert.  This finds the non-repeat intervals, which get turned into non-repeat tigs.

    tigMarksU = tigMarksR;
    tigMarksU.invert(0, tig->getLength());

    //  Create the list of intervals we'll use to make new tigs.

    vector<breakPointCoords>   BP;

    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++)
      BP.push_back(breakPointCoords(tigMarksR.lo(ii), tigMarksR.hi(ii), true));

    for (uint32 ii=0; ii<tigMarksU.numberOfIntervals(); ii++)
      BP.push_back(breakPointCoords(tigMarksU.lo(ii), tigMarksU.hi(ii), false));

    //  If there is only one BP, the tig is entirely resolved or entirely repeat.  Either case,
    //  there is nothing more for us to do.

    if (BP.size() == 1)
      continue;

    //  Report.

    sort(BP.begin(), BP.end());  //  Makes the report nice.  Doesn't impact splitting.

    writeLog("break tig %u into up to %u pieces:\n", ti, BP.size());
    for (uint32 ii=0; ii<BP.size(); ii++)
      writeLog("  %8d %8d %s (length %d)\n",
               BP[ii]._bgn, BP[ii]._end,
               BP[ii]._rpt ? "repeat" : "unique",
               BP[ii]._end - BP[ii]._bgn);

    //  Scan the reads, counting the number of reads that would be placed in each new tig.  This is done
    //  because there are a few 'splits' that don't move any reads around.

    Unitig **newTigs   = new Unitig * [BP.size()];
    int32   *lowCoord  = new int32    [BP.size()];
    uint32  *nRepeat   = new uint32   [BP.size()];
    uint32  *nUnique   = new uint32   [BP.size()];

    //  First call, count the number of tigs we would create if we let it create them.

    uint32  nTigs = splitTig(tigs, tig, BP, newTigs, lowCoord, nRepeat, nUnique, false);

    //  Second call, actually create the tigs, if anything would change.

    if (nTigs > 1)
      splitTig(tigs, tig, BP, newTigs, lowCoord, nRepeat, nUnique, true);

    //  Report the tigs created.

    reportTigsCreated(tig, BP, nTigs, newTigs, nRepeat, nUnique);

    //  Cleanup.

    delete [] newTigs;
    delete [] lowCoord;
    delete [] nRepeat;
    delete [] nUnique;

    //  Remove the old unitig....if we made new ones.

    if (nTigs > 1) {
      tigs[tig->id()] = NULL;
      delete tig;
    }
  }

#if 0
  FILE *F = AS_UTL_openOutputFile("junk.confusedEdges");
  for (uint32 ii=0; ii<confusedEdges.size(); ii++) {
    fprintf(F, "%7u %c' from read %7u\n",
            confusedEdges[ii].aid,
            confusedEdges[ii].a3p ? '3' : '5',
            confusedEdges[ii].bid);
  }
  AS_UTL_closeFile(F, "junk.confusedEdges");
#endif

  writeStatus("markRepeatReads()-- Found %u confused edges.\n", confusedEdges.size());
}

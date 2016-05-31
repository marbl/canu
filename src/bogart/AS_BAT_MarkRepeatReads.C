
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

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "intervalList.H"
#include "stddev.H"

#include <vector>

using namespace std;



//  Hack.
uint32 MIN_ANCHOR_HANG       = 500;   //  Require reads to be anchored by this many bases at boundaries of repeats.
int32  REPEAT_OVERLAP_MIN    = 50;

#define REPEAT_FRACTION   0.5

#undef  SHOW_ANNOTATION_RAW             //  Show all overlaps used to annotate reads
#undef  SHOW_ANNOTATION_RAW_FILTERED    //  Show all overlaps filtered by high error rate

#undef  DUMP_READ_COVERAGE

//  Each evidence read picks its single best overlap to tig (based on overlaps to reads in the tig).
//  Filter out evidence that aligns at erate higher than expected.
//  Collapse to intervals on tig.
//  If still not significant and not spanned, break.



class olapDat {
public:
  olapDat(uint32 b, uint32 e, uint32 t, uint32 r) {
    tigbgn  = b;
    tigend  = e;
    eviTid  = t;
    eviRid  = r;
  };

  bool operator<(const olapDat &that)     const { return(tigbgn < that.tigbgn); };

  uint32  tigbgn;   //  Location of the overlap on this tig
  uint32  tigend;   //

  uint32  eviTid;   //  tig that the evidence read came from
  uint32  eviRid;   //  evidence read
};


bool
olapDatByEviRid(const olapDat &A, const olapDat &B) {
  if (A.eviRid == B.eviRid)
    return(A.tigbgn < B.tigbgn);

  return(A.eviRid < B.eviRid);
}




class breakPointCoords {
public:
  breakPointCoords(uint32 tigID, int32 bgn, int32 end, bool rpt=false) {
    _tigID    = tigID;
    _bgn      = bgn;
    _end      = end;
    _isRepeat = rpt;
  };
  ~breakPointCoords() {
  };

  bool    operator<(breakPointCoords const &that) const {
    return(_bgn < that._bgn);
  };

  uint32  _tigID;
  int32   _bgn;
  int32   _end;
  bool    _isRepeat;
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
  hi = FI->fragmentLength(frg->ident);

  if (ahang > 0)
    lo += ahang;  //  Positive hang!

  if (bhang < 0)
    hi += bhang;  //  Negative hang!

  assert(0  <= lo);
  assert(0  <= hi);
  assert(lo <= hi);
  assert(lo <= FI->fragmentLength(frg->ident));
  assert(hi <= FI->fragmentLength(frg->ident));
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
  sprintf(fn, "%08u.coverage", tig->id());
  FILE *F = fopen(fn, "w");

  for (uint32 ii=0; ii<coverage.numberOfIntervals(); ii++)
    fprintf(F, "%u %u %u\n", coverage.lo(ii), coverage.hi(ii), coverage.depth(ii));

  fclose(F);
#endif
}








uint32
splitUnitigs(UnitigVector             &unitigs,
             Unitig                   *tig,
             vector<breakPointCoords> &BP,
             Unitig                  **newTigs,
             int32                    *lowCoord,
             uint32                   *nRepeat,
             uint32                   *nUnique,
             bool                      doMove) {
  uint32  nTigsCreated = 0;

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

    //fprintf(stderr, "Searching for placement for read %u at %u-%u\n", frg.ident, frgbgn, frgend);

    for (uint32 ii=0; ii<BP.size(); ii++) {
      int32   rgnbgn = BP[ii]._bgn;
      int32   rgnend = BP[ii]._end;
      bool    repeat = BP[ii]._isRepeat;

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
      fprintf(stderr, "Failed to place read %u at %u-%u\n", frg.ident, frgbgn, frgend);
      for (uint32 ii=0; ii<BP.size(); ii++)
        fprintf(stderr, "Breakpoints %2u %8u-%8u repeat %u\n", ii, BP[ii]._bgn, BP[ii]._end, BP[ii]._isRepeat);
    }
    assert(rid != UINT32_MAX);  //  We searched all the BP's, the read had better be placed!

    //  If moving reads, move the read!

    if (doMove) {
      if (newTigs[rid] == NULL) {
        lowCoord[rid] = frgbgn;

        newTigs[rid]  = unitigs.newUnitig(true);  // LOG_ADDUNITIG_BREAKING

        if (nRepeat[rid] > nUnique[rid])
          newTigs[rid]->_isRepeat = true;
      }

      newTigs[rid]->addFrag(frg, -lowCoord[rid], false);  //LOG_ADDFRAG_BREAKING);
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

  for (uint32 ii=0; ii<BP.size(); ii++)
    if (nRepeat[ii] + nUnique[ii] > 0)
      nTigsCreated++;

  return(nTigsCreated);
}






//  For each overlap, if the b-read is in this tig, ignore it.
//  Otherwise annotate the read with the overlap region.
//
//  Later, check if the two reads in this unitig overlap; if not, annotate also.
//
void
annotateRepeatsOnRead(UnitigVector          &unitigs,
                      Unitig                *tgA,
                      ufNode                *rdA,
                      double                 deviationRepeat,
                      vector<olapDat>       &repeats) {
  uint32               ovlLen   = 0;
  BAToverlap          *ovl      = OC->getOverlaps(rdA->ident, AS_MAX_ERATE, ovlLen);

  vector<olapDat>      readOlaps;        //  List of valid repeat overlaps to this read

  uint32               tgAid    = tgA->id();

  bool                 rdAfwd   = (rdA->position.bgn < rdA->position.end);
  int32                rdAlo    = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
  int32                rdAhi    = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

  assert(rdAlo < rdAhi);

  //  Beacuse the read is placed with a lot of fudging in the positions, we need
  //  to scale the coordinates we compute here.
  double               sc       = (rdAhi - rdAlo) / (double)FI->fragmentLength(rdA->ident);

  //  For all overlaps to this read, save the overlap if it is not represented in this tig.

  uint32  nOlaps = ovlLen;
  uint32  nDiff  = 0;       //  Overlap to different tig
  uint32  nSelf  = 0;       //  Overlap to same tig, different location
  uint32  nConf  = 0;       //  Overlap to same tig, confirmed good overlap

  for (uint32 oi=0; oi<ovlLen; oi++) {
    uint32   rdBid  = ovl[oi].b_iid;
    uint32   tgBid  = Unitig::fragIn(rdBid);

    int32    bgn    = 0;  //  Position in the read that
    int32    end    = 0;  //  the overlap covers

    //  If the read is in a singleton, skip.  These are unassembled crud.
    if ((tgBid                         == 0) ||
        (unitigs[tgBid]                == NULL) ||
        (unitigs[tgBid]->ufpath.size() == 1))
      continue;

    //  If the read is in an annotated bubble, skip.
    if (unitigs[tgBid]->_isBubble)
      continue;

    //  If the overlap is to a container read, skip it.
    if ((ovl[oi].a_hang < 0) && (ovl[oi].b_hang > 0))
      continue;

    //  If the overlap is to a contained read, skip it.
    if ((ovl[oi].a_hang > 0) && (ovl[oi].b_hang < 0))
      continue;

    uint32   rdBpos   =  unitigs[tgBid]->pathPosition(rdBid);
    ufNode  *rdB      = &unitigs[tgBid]->ufpath[rdBpos];

    bool     rdBfwd   = (rdB->position.bgn < rdB->position.end);
    int32    rdBlo    = (rdBfwd) ? rdB->position.bgn : rdB->position.end;
    int32    rdBhi    = (rdBfwd) ? rdB->position.end : rdB->position.bgn;

    assert(rdBlo < rdBhi);

    //  If the overlap is to a read in a different tig, save it.
    if (tgBid != tgAid) {
      nDiff++;
      olapToReadCoords(rdA, ovl[oi].a_hang, ovl[oi].b_hang, bgn, end);
    }

    // If the overlap is to a read in the same tig, but we don't overlap in the tig, save it.
    else if ((rdAhi < rdBlo) || (rdBhi < rdAlo)) {
      nSelf++;
      olapToReadCoords(rdA, ovl[oi].a_hang, ovl[oi].b_hang, bgn, end);
    }

    //  Otherwise, the overlap is present in the tig, and can't indicate a repeat.
    else {
      nConf++;
      continue;
    }

    //  Find the positions of the read that are covered by the overlap.

    int32  tigbgn = (rdAfwd) ? (rdAlo + sc * bgn) : (rdAhi - sc * end);
    int32  tigend = (rdAfwd) ? (rdAlo + sc * end) : (rdAhi - sc * bgn);

    assert(tigbgn < tigend);

    if (tigbgn < 0)                  tigbgn = 0;
    if (tigend > tgA->getLength())   tigend = tgA->getLength();

    //  Filter overlaps that are higher error than expected.

    double   consistent = tgA->overlapConsistentWithTig(deviationRepeat, tigbgn, tigend, ovl[oi].erate);

    if (consistent < REPEAT_FRACTION) {
#ifdef SHOW_ANNOTATION_RAW_FILTERED
      writeLog("tig %6u read %7u %8u-%8u OVERLAP from tig %6u read %7u %8u-%8u at tigpos %8u-%8u erate %.6f consistent %.4f FILTERED\n",
               tgAid, rdA->ident, rdAlo, rdAhi,
               tgBid, rdBid, rdBlo, rdBhi,
               tigbgn, tigend, ovl[oi].erate, consistent);
#endif
      continue;
    }

#ifdef SHOW_ANNOTATION_RAW
    writeLog("tig %6u read %7u %8u-%8u OVERLAP from tig %6u read %7u %8u-%8u at tigpos %8u-%8u erate %.6f consistent %.4f\n",
             tgAid, rdA->ident, rdAlo, rdAhi,
             tgBid, rdBid, rdBlo, rdBhi,
             tigbgn, tigend, ovl[oi].erate, consistent);
#endif

    readOlaps.push_back(olapDat(tigbgn, tigend, tgBid, rdBid));
  }

  //  All overlaps processed.  Save to the master list.

#pragma omp critical (repeatsPushBack)
  for (uint32 rr=0; rr<readOlaps.size(); rr++)
    repeats.push_back(readOlaps[rr]);
}




void
markRepeatReads(UnitigVector &unitigs,
                double        deviationRepeat,
                uint32        confusedAbsolute,
                double        confusedPercent) {
  uint32  tiLimit = unitigs.size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  writeLog("repeatDetect()-- working on "F_U32" unitigs, with "F_U32" threads.\n", tiLimit, numThreads);

  vector<olapDat>      repeatOlaps;   //  Overlaps to reads promoted to tig coords

  intervalList<int32>  tigMarksR;     //  Marked repeats based on reads, filtered by spanning reads
  intervalList<int32>  tigMarksU;     //  Non-repeat invervals, just the inversion of tigMarksR


  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = unitigs[ti];

    if (tig == NULL)
      continue;

    if (tig->ufpath.size() == 1)
      continue;

    vector<olapDat>   repeats;

    writeLog("Annotating repeats in reads for tig %u/%u.\n", ti, tiLimit);

    //  Clear out all the existing marks.  They're not for this tig.


    //  Analyze overlaps for each read.  For each overlap to a read not in this tig, or not
    //  overlapping in this tig, and of acceptable error rate, add the overlap to repeatOlaps.

    repeatOlaps.clear();

    uint32  fiLimit    = tig->ufpath.size();
    uint32  numThreads = omp_get_max_threads();
    uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

#pragma omp parallel for if(fiLimit > 100) schedule(dynamic, blockSize)
    for (uint32 fi=0; fi<fiLimit; fi++)
      annotateRepeatsOnRead(unitigs, tig, &tig->ufpath[fi], deviationRepeat, repeatOlaps);

    writeLog("Annotated with %lu overlaps.\n", repeatOlaps.size());

    //  Merge marks for the same read into the largest possible.

    sort(repeatOlaps.begin(), repeatOlaps.end(), olapDatByEviRid);

#ifdef SHOW_ANNOTATE
    for (uint32 ii=0; ii<repeatOlaps.size(); ii++)
      if (repeatOlaps[ii].tigbgn < 1000000)
        writeLog("repeatOlaps[%u] %u-%u from tig %u read %u RAW\n",
                 ii,
                 repeatOlaps[ii].tigbgn, repeatOlaps[ii].tigend,
                 repeatOlaps[ii].eviTid, repeatOlaps[ii].eviRid);

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
      repeatOlaps[ss].eviTid = UINT32_MAX;
      repeatOlaps[ss].eviRid = UINT32_MAX;
    }

    //  Sort overlaps again.  This pushes all those 'erased' regions to the end of the list, which
    //  we can then just pop off.

    sort(repeatOlaps.begin(), repeatOlaps.end(), olapDatByEviRid);

    for (uint32 ii=repeatOlaps.size(); ii--; )
      if (repeatOlaps[ii].eviTid == UINT32_MAX)
        repeatOlaps.pop_back();

    //  For logging, sort by coordinate

    sort(repeatOlaps.begin(), repeatOlaps.end());

#ifdef SHOW_ANNOTATE
    for (uint32 ii=0; ii<repeatOlaps.size(); ii++)
      writeLog("repeatOlaps[%d] %u-%u from tig %u read %u MERGED\n",
               ii,
               repeatOlaps[ii].tigbgn, repeatOlaps[ii].tigend,
               repeatOlaps[ii].eviTid, repeatOlaps[ii].eviRid);
#endif

    //  Make a new set of intervals based on all the detected repeats.

    tigMarksR.clear();

    for (uint32 bb=0, ii=0; ii<repeatOlaps.size(); ii++)
      tigMarksR.add(repeatOlaps[ii].tigbgn, repeatOlaps[ii].tigend - repeatOlaps[ii].tigbgn);

    //  Collapse these markings Collapse all the read markings to intervals on the unitig, merging those that overlap
    //  significantly.

    writeLog("Merge marks.\n");

    tigMarksR.merge(REPEAT_OVERLAP_MIN);

    //  Scan reads, discard any mark that is contained in a read
    //
    //  We don't need to filterShort() after every one is removed, but it's simpler to do it Right Now than
    //  to track if it is needed.

    writeLog("Scan reads to discard spanned repeats.\n");

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

    //  Run through again, looking for the thickest overlap(s) to the remaining regions.
    //  This isn't caring about the end effect noted above.

#if 1
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
#endif


    //  Scan reads.  If a read intersects a repeat interval, and the best edge for that read
    //  is entirely in the repeat region, decide if there is a near-best edge to something
    //  not in this tig.
    //
    //  A region with no such near-best edges is _probably_ correct.

    writeLog("search for confused edges:\n");

    uint32  *isConfused  = new uint32 [tigMarksR.numberOfIntervals()];

    memset(isConfused, 0, sizeof(uint32) * tigMarksR.numberOfIntervals());

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode     *rdA       = &tig->ufpath[fi];
      uint32      rdAid     = rdA->ident;
      bool        rdAfwd    = (rdA->position.bgn < rdA->position.end);
      int32       rdAlo     = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
      int32       rdAhi     = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

      double      sc        = (rdAhi - rdAlo) / (double)FI->fragmentLength(rdAid);

      if ((OG->isContained(rdAid)  == true) ||
          (OG->isSuspicious(rdAid) == true))
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

        if (b5->fragId() == 0)
          b5use = false;
        if (b3->fragId() == 0)
          b3use = false;

        if ((b5use) && (Unitig::fragIn(b5->fragId()) != tig->id()))
          b5use = false;
        if ((b3use) && (Unitig::fragIn(b3->fragId()) != tig->id()))
          b3use = false;

        //  The best edge read is in this tig.  If they don't overlap, again, nothing to compare
        //  against.

        if (b5use) {
          ufNode     *rdB       = &tig->ufpath[Unitig::pathPosition(b5->fragId())];
          uint32      rdBid     = rdB->ident;
          bool        rdBfwd    = (rdB->position.bgn < rdB->position.end);
          int32       rdBlo     = (rdBfwd) ? rdB->position.bgn : rdB->position.end;
          int32       rdBhi     = (rdBfwd) ? rdB->position.end : rdB->position.bgn;

          if ((rdAhi < rdBlo) ||
              (rdBhi < rdAlo))
            b5use = false;
        }

        if (b3use) {
          ufNode     *rdB       = &tig->ufpath[Unitig::pathPosition(b3->fragId())];
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
          len5 = FI->overlapLength(rdAid, b5->fragId(), b5->ahang(), b5->bhang());
        else
          b5use = false;

        if ((rMin    < tig3bgn) &&
            (tig3end < rMax) &&
            (b3use))
          len3 = FI->overlapLength(rdAid, b3->fragId(), b3->ahang(), b3->bhang());
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
        BAToverlap   *ovl      = OC->getOverlaps(rdAid, AS_MAX_ERATE, ovlLen);

        for (uint32 oo=0; oo<ovlLen; oo++) {
          uint32   rdBid    = ovl[oo].b_iid;
          uint32   tgBid    = Unitig::fragIn(rdBid);

          //  If the read is in a singleton, skip.  These are unassembled crud.
          if ((tgBid                         == 0) ||
              (unitigs[tgBid]                == NULL) ||
              (unitigs[tgBid]->ufpath.size() == 1))
            continue;

          //  If the read is in an annotated bubble, skip.
          if (unitigs[tgBid]->_isBubble)
            continue;

          //  Skip if this overlap is the best we're trying to match.
          if ((rdBid == b5->fragId()) ||
              (rdBid == b3->fragId()))
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


          uint32   rdBpos   =  unitigs[tgBid]->pathPosition(rdBid);
          ufNode  *rdB      = &unitigs[tgBid]->ufpath[rdBpos];

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

          uint32  len   = FI->overlapLength(rdAid, ovl[oo].b_iid, ovl[oo].a_hang, ovl[oo].b_hang);
          double  score = len * (1 - ovl[oo].erate);

          //  Compute percent difference.

          double  ad5 = fabs(score - score5);
          double  ad3 = fabs(score - score3);

          double  pd5 = 200 * ad5 / (score + score5);
          double  pd3 = 200 * ad3 / (score + score3);

          //  Skip if this overlap is vastly worse than the best.

          if ((ovl5 == true) && ((ad5 >= confusedAbsolute) || (pd3 > confusedPercent))) {
            writeLog("tig %7u read %8u pos %7u-%-7u NOT confused by 5' edge to read %8u - best edge read %8u len %6u erate %.4f score %8.2f - alt edge len %6u erate %.4f score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                     tig->id(), rdAid, rdAlo, rdAhi,
                     rdBid,
                     b5->fragId(), len5, b5->erate(), score5,
                     len, ovl[oo].erate, score,
                     ad5, pd5);
            continue;
          }

          if ((ovl3 == true) && ((ad3 >= confusedAbsolute) || (pd3 > confusedPercent))) {
            writeLog("tig %7u read %8u pos %7u-%-7u NOT confused by 3' edge to read %8u - best edge read %8u len %6u erate %.4f score %8.2f - alt edge len %6u erate %.4f score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                     tig->id(), rdAid, rdAlo, rdAhi,
                     rdBid,
                     b3->fragId(), len3, b3->erate(), score3,
                     len, ovl[oo].erate, score,
                     ad3, pd3);
            continue;
          }

          //  Potential confusion!

          if (ovl5 == true)
            writeLog("tig %7u read %8u pos %7u-%-7u IS confused by 5' edge to read %8u - best edge read %8u len %6u erate %.4f score %8.2f - alt edge len %6u erate %.4f score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                     tig->id(), rdAid, rdAlo, rdAhi,
                     rdBid,
                     b5->fragId(), len5, b5->erate(), score5,
                     len, ovl[oo].erate, score,
                     ad5, pd5);

          if (ovl3 == true)
            writeLog("tig %7u read %8u pos %7u-%-7u IS confused by 3' edge to read %8u - best edge read %8u len %6u erate %.4f score %8.2f - alt edge len %6u erate %.4f score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                     tig->id(), rdAid, rdAlo, rdAhi,
                     rdBid,
                     b3->fragId(), len3, b3->erate(), score3,
                     len, ovl[oo].erate, score,
                     ad3, pd3);

          isConfused[ri]++;
        }
      }  //  Over all marks (ri)
    }  //  Over all reads (fi)


    //  Scan all the regions, and delete any that have no confusion.

    {
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
    }

    delete [] isConfused;





    //  Scan reads, join any marks that have their junctions spanned by a sufficiently large amount.
    //
    //  If the read spans this junction be the usual amount, merge the intervals.
    //
    //  The intervals can be overlapping (by up to REPEAT_OVERLAP_MIN (x2?) bases.  For this junction
    //  to be spanned, the read must span from min-ROM to max+ROM, not just hi(ri-1) to lo(ri).
    //
    //  We DO need to filterShort() after every merge, otherwise, we'd have an empty bogus interval
    //  in the middle of our list, which could be preventing some other merge.  OK, we could
    //
    //  Anything that gets merged is now no longer a true repeat.  It's unique, just bordered by repeats.
    //  We can't track this through the indices (because we delete things).  We track it with a set of
    //  begin coordinates.

    set<int32>  nonRepeatIntervals;

    writeLog("Scan reads to merge repeat regions.\n");

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode     *frg       = &tig->ufpath[fi];
      bool        frgfwd    = (frg->position.bgn < frg->position.end);
      int32       frglo     = (frgfwd) ? frg->position.bgn : frg->position.end;
      int32       frghi     = (frgfwd) ? frg->position.end : frg->position.bgn;
      bool        merged    = false;

      for (uint32 ri=1; ri<tigMarksR.numberOfIntervals(); ri++) {
        uint32  rMin = min(tigMarksR.hi(ri-1), tigMarksR.lo(ri));
        uint32  rMax = max(tigMarksR.hi(ri-1), tigMarksR.lo(ri));

        if ((frglo + MIN_ANCHOR_HANG <= rMin) && (rMax + MIN_ANCHOR_HANG <= frghi)) {
          writeLog("merge regions %8d:%-8d and %8d:%-8d - junction contained in read %6u %5d-%5d\n",
                   tigMarksR.lo(ri-1), tigMarksR.hi(ri-1),
                   tigMarksR.lo(ri), tigMarksR.hi(ri),
                   frg->ident, frglo, frghi);

          tigMarksR.lo(ri) = tigMarksR.lo(ri-1);

          tigMarksR.lo(ri-1) = 0;   //  CRITICAL to delete this interval (and not ri) because the next
          tigMarksR.hi(ri-1) = 0;   //  iteration will be using ri-1 (== ri here) and ri (== ri+1).

          merged = true;

          nonRepeatIntervals.insert(tigMarksR.lo(ri));
        }
      }

      if (merged)
        tigMarksR.filterShort(1);
    }

    //  Extend the regions by MIN_ANCHOR_HANG.  This makes checking for reads that span and are
    //  anchored in the next region easier.  It also solved a quirk when the first/last repeat
    //  region doesn't extend to the end of the sequence:
    //    0-183     unique  (created from inversion below, but useless and incorrect)
    //    183-9942  repeat

    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++) {
      tigMarksR.lo(ii) = max<int32>(tigMarksR.lo(ii) - MIN_ANCHOR_HANG, 0);
      tigMarksR.hi(ii) = min<int32>(tigMarksR.hi(ii) + MIN_ANCHOR_HANG, tig->getLength());
    }

    //  Find the non-repeat intervals.

    tigMarksU = tigMarksR;
    tigMarksU.invert(0, tig->getLength());

    //  Create the list of intervals we'll use to make new unitigs.
    //
    //  The repeat intervals are extended by MIN_ANCHOR_HANG, and then any read fully contained in one of
    //  these is moved here.
    //
    //  The non-repeat intervals are shortened by the same amount, and any read that intersects one
    //  is moved there.
    //
    //  Does order matter?  Not sure.  The repeat intervals are first, then the formerly repeat
    //  merged intervals, then the unique intervals.  Splitting might depend on the repeats being
    //  first.

    writeLog("Make breakpoints.\n");

    vector<breakPointCoords>   BP;

    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++)
      if (nonRepeatIntervals.count(tigMarksR.lo(ii)) == 0)
        BP.push_back(breakPointCoords(ti, tigMarksR.lo(ii), tigMarksR.hi(ii), true));

    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++)
      if (nonRepeatIntervals.count(tigMarksR.lo(ii)) != 0)
        BP.push_back(breakPointCoords(ti, tigMarksR.lo(ii), tigMarksR.hi(ii), true));

    for (uint32 ii=0; ii<tigMarksU.numberOfIntervals(); ii++) {
      BP.push_back(breakPointCoords(ti, tigMarksU.lo(ii), tigMarksU.hi(ii), false));
    }

    //  If only one region, the whole unitig was declared repeat.  Nothing to do.

    if (BP.size() == 1)
      continue;

    sort(BP.begin(), BP.end());

    //  Report.

    writeLog("break tig %u into up to %u pieces:\n", ti, BP.size());
    for (uint32 ii=0; ii<BP.size(); ii++)
      writeLog("  %8d %8d %s (length %d)\n",
               BP[ii]._bgn, BP[ii]._end,
               BP[ii]._isRepeat ? "repeat" : "unique",
               BP[ii]._end - BP[ii]._bgn);

    //  Scan the reads, counting the number of reads that would be placed in each new tig.  This is done
    //  because there are a few 'splits' that don't move any reads around.

    Unitig **newTigs   = new Unitig * [BP.size()];
    int32   *lowCoord  = new int32    [BP.size()];
    uint32  *nRepeat   = new uint32   [BP.size()];
    uint32  *nUnique   = new uint32   [BP.size()];

    //  First call, count the number of tigs we would create if we let it create them.

    uint32  nTigs = splitUnitigs(unitigs, tig, BP, newTigs, lowCoord, nRepeat, nUnique, false);

    //  Second call, actually create the tigs, if anything would change.

    if (nTigs > 1)
      splitUnitigs(unitigs, tig, BP, newTigs, lowCoord, nRepeat, nUnique, true);

    //  Report the tigs created.

    for (uint32 ii=0; ii<BP.size(); ii++) {
      int32   rgnbgn = BP[ii]._bgn;
      int32   rgnend = BP[ii]._end;
      bool    repeat = BP[ii]._isRepeat;

      if      (nRepeat[ii] + nUnique[ii] == 0)
        writeLog("For tig %5u %s region %8d %8d - %6u/%6u repeat/unique reads - no new unitig created.\n",
                 ti, (repeat == true) ? "repeat" : "unique", rgnbgn, rgnend, nRepeat[ii], nUnique[ii]);

      else if (nTigs > 1)
        writeLog("For tig %5u %s region %8d %8d - %6u/%6u reads repeat/unique - unitig %5u created.\n",
                 ti, (repeat == true) ? "repeat" : "unique", rgnbgn, rgnend, nRepeat[ii], nUnique[ii], newTigs[ii]->id());

      else
        writeLog("For tig %5u %s region %8d %8d - %6u/%6u repeat/unique reads - unitig %5u remains unchanged.\n",
                 ti, (repeat == true) ? "repeat" : "unique", rgnbgn, rgnend, nRepeat[ii], nUnique[ii], tig->id());
    }

    //  Cleanup.

    delete [] newTigs;
    delete [] lowCoord;
    delete [] nRepeat;
    delete [] nUnique;

    //  Remove the old unitig....if we made new ones.

    if (nTigs > 1) {
      delete tig;
      unitigs[ti] = NULL;
    }
  }
}


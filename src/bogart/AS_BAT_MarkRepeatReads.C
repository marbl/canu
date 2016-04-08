
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
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_OverlapCache.H"

#include "intervalList.H"
#include "stddev.H"

#include <vector>

using namespace std;

//  Hack.
uint32 MIN_ANCHOR_HANG       = 500;   //  Require reads to be anchored by this many bases at boundaries of repeats.
int32  REPEAT_OVERLAP_MIN    = 50;

#define REPEAT_CONSISTENT 3.0

#undef  SHOW_ANNOTATION_RAW             //  Show all overlaps used to annotate reads
#undef  SHOW_ANNOTATION_RAW_FILTERED    //  Show all overlaps filtered by high error rate

#undef  DUMP_READ_COVERAGE
#undef  DUMP_ERROR_PROFILE

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
olapToReadCoords(ufNode *frg, BAToverlap *olap, int32 &lo, int32 &hi) {

  lo = 0;
  hi = FI->fragmentLength(frg->ident);

  if (olap->a_hang > 0)
    lo += olap->a_hang;  //  Postiive hang!

  if (olap->b_hang < 0)
    hi += olap->b_hang;  //  Negative hang!

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

  memset(newTigs,  0, sizeof(Unitig *) * BP.size());
  memset(lowCoord, 0, sizeof(int32)    * BP.size());
  memset(nRepeat,  0, sizeof(uint32)   * BP.size());
  memset(nUnique,  0, sizeof(uint32)   * BP.size());

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
      }

      newTigs[rid]->addFrag(frg, -lowCoord[rid], false);  //LOG_ADDFRAG_BREAKING);
    }

    //  Count how many reads cam from repeats or uniques.

    if (rpt)
      nRepeat[rid]++;
    else
      nUnique[rid]++;
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
                      double                 erateRepeat,
                      vector<olapDat>       &repeats) {
  uint32               ovlLen   = 0;
  BAToverlap          *ovl      = OC->getOverlaps(rdA->ident, erateRepeat, ovlLen);

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
      olapToReadCoords(rdA, ovl+oi, bgn, end);
    }

    // If the overlap is to a read in the same tig, but we don't overlap in the tig, save it.
    else if ((rdAhi < rdBlo) || (rdBhi < rdAlo)) {
      nSelf++;
      olapToReadCoords(rdA, ovl+oi, bgn, end);
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

    double   consistent = tgA->overlapConsistentWithTig(REPEAT_CONSISTENT, tigbgn, tigend, ovl[oi].erate);

    if (consistent < 0.5) {
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
                const char   *UNUSED(prefix),
                double        UNUSED(erateGraph),
                double        UNUSED(erateBubble),
                double        UNUSED(erateMerge),
                double        erateRepeat) {
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
      annotateRepeatsOnRead(unitigs, tig, &tig->ufpath[fi], erateRepeat, repeatOlaps);

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

#if 0
#ifdef SHOW_ANNOTATE
    writeLog("PRE-MERGE markings\n");
    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++) {
      writeLog("  %8d:%-8d size %7d (distance to next %7d)\n",
               tigMarksR.lo(ii), tigMarksR.hi(ii), tigMarksR.hi(ii) - tigMarksR.lo(ii),
               (ii < tigMarksR.numberOfIntervals()-1) ? (tigMarksR.lo(ii+1) - tigMarksR.hi(ii)) : (0));
    }
#endif
#endif

    tigMarksR.merge(REPEAT_OVERLAP_MIN);

#ifdef SHOW_ANNOTATE
    writeLog("POST-MERGE markings\n");
    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++) {
      writeLog("  %8d:%-8d size %7d (distance to next %7d)\n",
               tigMarksR.lo(ii), tigMarksR.hi(ii), tigMarksR.hi(ii) - tigMarksR.lo(ii),
               (ii < tigMarksR.numberOfIntervals()-1) ? (tigMarksR.lo(ii+1) - tigMarksR.hi(ii)) : (0));
    }
#endif

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

        if ((frglo == 0) &&
            (tigMarksR.hi(ri) + MIN_ANCHOR_HANG <= frghi))
          spanLo = spanHi = true;

        if ((frghi == tig->getLength()) &&
            (frglo + MIN_ANCHOR_HANG <= tigMarksR.lo(ri)))
          spanLo = spanHi = true;

        if (frglo + MIN_ANCHOR_HANG <= tigMarksR.lo(ri))
          spanLo = true;

        if (tigMarksR.hi(ri) + MIN_ANCHOR_HANG <= frghi)
          spanHi = true;

        if (spanLo && spanHi) {
          writeLog("discard region %8d:%-8d - contained in read %6u %5d-%5d\n",
                   tigMarksR.lo(ri), tigMarksR.hi(ri), frg->ident, frglo, frghi);

          tigMarksR.lo(ri) = 0;
          tigMarksR.hi(ri) = 0;

          discarded = true;
        }
      }

      if (discarded)
        tigMarksR.filterShort(1);
    }

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


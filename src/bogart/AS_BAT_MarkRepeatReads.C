
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
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
      fprintf(stderr, "Failed to place read %u at %d-%d with %lu BPs:\n", frg.ident, frgbgn, frgend, BP.size());
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



//  Build a vector of olapDat (tigBgn, tigEnd, eviRid) for all reads that
//  overlap into this tig.
void
annotateRepeatsOnRead(AssemblyGraph   *AG,
                      Unitig          *tig,
                      vector<olapDat> &repeats) {

  repeats.clear();

  for (uint32 ii=0; ii<tig->ufpath.size(); ii++) {
    ufNode               *read   = &tig->ufpath[ii];
    vector<BestReverse>  &rPlace = AG->getReverse(read->ident);

    for (uint32 rr=0; rr<rPlace.size(); rr++) {
      uint32          rID    = rPlace[rr].readID;
      uint32          pID    = rPlace[rr].placeID;
      BestPlacement  &fPlace = AG->getForward(rID)[pID];

      //  Also in AS_BAT_AssemblyGraph.C
      if (OG->isBubble(rID))              //  Ignore if the incoming overlap is from a bubble.
        continue;

      if (OG->isSpur(rID))               //  Ignore if the incoming overlap is from a spur.
        continue;

      //writeLog("annotateRepeatsOnRead()-- tig %u read #%u %u place %u reverse read %u in tig %u olap %d-%d%s\n",
      //         tig->id(), ii, read->ident, rr,
      //         rID,
      //         tig->inUnitig(rID),
      //         fPlace.olapMin,   fPlace.olapMax,
      //         (fPlace.isUnitig) ? " IN_UNITIG" : "");

      repeats.push_back(olapDat(fPlace.olapMin, fPlace.olapMax, rID, pID));
    }
  }

  writeLog("Annotated tig %u with %lu external overlaps.\n", tig->id(), repeats.size());
}



void
mergeAnnotations(vector<olapDat>      &repeatOlaps,
                 intervalList<int32>  &tigMarksR) {

  //  Sort the repeat markings by their evidence read id.

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

  //  Merge overlapping regions that come from the same evidence read.
  //
  //      -------------oo
  //          ---------oooooo
  //              ---------------
  //                 ---------------      This evidence read has two overlaps
  //                   ======   <-------  to the tig; we want to merge those
  //                         \            overlaps into one region.
  //                          |

  for (uint32 dd=0, ss=1; ss<repeatOlaps.size(); ss++) {
    assert(repeatOlaps[dd].eviRid <= repeatOlaps[ss].eviRid);

    //  If this olapDat is from a different evidence read, or if it has no
    //  intersection with the destination olapDat region, reset the
    //  destination olapDat pointer to the current olapDat.

    if ((repeatOlaps[dd].eviRid != repeatOlaps[ss].eviRid) ||
        (repeatOlaps[dd].tigend <= repeatOlaps[ss].tigbgn)) {
      dd = ss;
    }

    //  Otherwise, this olapDat overlaps with the region we're merging into.
    //  So merge it in too, then mark it as invalid.

    else {
      repeatOlaps[dd].tigbgn = min(repeatOlaps[ss].tigbgn, repeatOlaps[dd].tigbgn);
      repeatOlaps[dd].tigend = max(repeatOlaps[ss].tigend, repeatOlaps[dd].tigend);

      repeatOlaps[ss].tigbgn = UINT32_MAX;
      repeatOlaps[ss].tigend = UINT32_MAX;
      repeatOlaps[ss].eviRid = UINT32_MAX;
      repeatOlaps[ss].eviPid = UINT32_MAX;
    }
  }

  //  Sort overlaps again.  This pushes all those 'erased' regions to the end
  //  of the list, which we can then just pop off.

  sort(repeatOlaps.begin(), repeatOlaps.end(), olapDatByEviRid);

  for (uint32 ii=repeatOlaps.size(); ii--; )
    if (repeatOlaps[ii].eviRid == UINT32_MAX)
      repeatOlaps.pop_back();

  //  Sort by coordinate.

  sort(repeatOlaps.begin(), repeatOlaps.end());

#ifdef SHOW_ANNOTATE
  for (uint32 ii=0; ii<repeatOlaps.size(); ii++)
    if (repeatOlaps[ii].tigbgn < 1000000)
      writeLog("repeatOlaps[%d] %d-%d from tig %u read %u place %u MERGED\n",
               ii,
               repeatOlaps[ii].tigbgn, repeatOlaps[ii].tigend,
               repeatOlaps[ii].eviRid, repeatOlaps[ii].eviPid);
#endif

  //  Merge any adjcant regions if they overlap significantly.  If they
  //  don't, then this is just two independent repeats next to each other.

  tigMarksR.clear();
  for (uint32 bb=0, ii=0; ii<repeatOlaps.size(); ii++)
    tigMarksR.add(repeatOlaps[ii].tigbgn, repeatOlaps[ii].tigend - repeatOlaps[ii].tigbgn);

  tigMarksR.merge(REPEAT_OVERLAP_MIN);
}




void
discardSpannedRepeats(Unitig              *tig,
                      intervalList<int32> &tigMarksR) {

  //  Somewhat inefficiently, for each read in the tig, check if it covers
  //  each repeat region.

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode     *frg       = &tig->ufpath[fi];
    bool        frgfwd    = (frg->position.bgn < frg->position.end);
    int32       frglo     = (frgfwd) ? frg->position.bgn : frg->position.end;
    int32       frghi     = (frgfwd) ? frg->position.end : frg->position.bgn;
    bool        discarded = false;

    //  A read is spanning the repeat region if is at least MIN_ANCHOR_HANG
    //  bases larger than the region (on each end).  If the region is at
    //  the end of a tig, all we can ask for is that the read extends
    //  fully to the end of the tig.

    for (uint32 ri=0; ri<tigMarksR.numberOfIntervals(); ri++) {
      bool   spanLo = false;
      bool   spanHi = false;

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

    //  If any region was discarded, filter it out now, so that we don't
    //  falsely discard it again on later reads.

    if (discarded)
      tigMarksR.filterShort(1);
  }
}



//  This is purely logging.  For each repeat region, report
//  the thickest overlap between it and any read in the tig.
void
reportThickestEdgesInRepeats(Unitig               *tig,
                             intervalList<int32>  &tigMarksR) {

  writeLog("thickest edges to the repeat regions:\n");

  //  The variable names make no sense.
  //    t5    - read id with the thickest overlap covering the low end of the repeat region
  //    l5    - length of that overlap
  //    t5bgn - position of the read
  //    t5end - position of the read

  for (uint32 rr=0; rr<tigMarksR.numberOfIntervals(); rr++) {
    int32    rbgn = tigMarksR.lo(rr);
    int32    rend = tigMarksR.hi(rr);

    uint32   fi5 = UINT32_MAX, len5 = 0;
    uint32   fi3 = UINT32_MAX, len3 = 0;

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode     *frg       = &tig->ufpath[fi];
      int32       frglo     = frg->position.min();
      int32       frghi     = frg->position.max();

      //  Overlap off the 5' end of the region.
      if (frglo <= rbgn && (rbgn <= frghi) && (len5 < frghi - rbgn)) {
        fi5  = fi;
        len5 = frghi - rbgn;
      }

      //  Overlap off the 3' end of the region.
      if (frglo <= rend && (rend <= frghi) && (len3 < rend - frglo)) {
        fi3  = fi;
        len3 = rend - frglo;
      }

      //  Read cotains the overlap region.  This shouldn't happen frequently,
      //  but can if the read is just barely covering the repeat.

      if (frglo <= rbgn && (rend <= frghi))
        writeLog("saved   region %8d:%-8d - thickest    read %6u (anchor %7d) %8d:%-8d (anchor %7d) (contained)\n",
                 rbgn, rend, frg->ident, rbgn - frglo, frglo, frghi, frghi - rend);
    }

    //  Report the thickest overlap into the repeat.

    if (fi5 != UINT32_MAX) {
      ufNode     *frg       = &tig->ufpath[fi5];
      int32       frglo     = frg->position.min();
      int32       frghi     = frg->position.max();

      writeLog("saved   region %8d:%-8d - thickest 5' read %6u (anchor %7d) %8d:%-8d (repeat %7d)\n",
               rbgn, rend, frg->ident, rbgn - frglo, frglo, frghi, frghi - rbgn);
    }

    if (fi3 != UINT32_MAX) {
      ufNode     *frg       = &tig->ufpath[fi3];
      int32       frglo     = frg->position.min();
      int32       frghi     = frg->position.max();

      writeLog("saved   region %8d:%-8d - thickest 3' read %6u (repeat %7d) %8d:%-8d (anchor %7d)\n",
               rbgn, rend, frg->ident, rend - frglo, frglo, frghi, frghi - rend);
    }
  }
}



//  Returns true if coord is equal to or inside the region.
bool
isInside(int32 lo, int32 coord, int32 hi) {
  return((lo <= coord) && (coord <= hi));
}


uint32
findThickestPrevRead(Unitig *tig, uint32 fi, int32 rMin, int32 rMax) {
  ufNode     *rdA       = &tig->ufpath[fi];
  uint32      rdAid     = rdA->ident;
  bool        rdAfwd    = rdA->position.isForward();
  int32       rdAlo     = rdA->position.min();
  int32       rdAhi     = rdA->position.max();

  uint32   bestId  = 0;
  uint32   bestIdx = 0;
  uint32   bestLen = 0;
  uint32   olapMin = 0;
  uint32   olapMax = 0;;

  //  If rdA begins before rMin, we don't need to find the previous read.
  //  it's not confused.

  if (rdAlo < rMin)
    return(UINT32_MAX);

  //  Otherwise, search for the previous best read.

  for (int32 pi=fi-1; pi>0; pi--) {
    ufNode  *rdB   = &tig->ufpath[pi];
    int32    rdBlo = rdB->position.min();
    int32    rdBhi = rdB->position.max();

    if (OG->isContained(rdB->ident) == true)   //  Skip contained reads.
      continue;

    if (OG->isBackbone(rdB->ident) == false)   //  Skip non-backbone reads.
      continue;

    //if (rdAhi < rdBlo)                         //  We can stop if the read
    //  break;                                   //  starts after us.

    if ((rdBhi >= rdAlo) && (bestLen < rdBhi - rdAlo)) {
      bestId   = rdB->ident;
      bestIdx  = pi;
      bestLen  = rdBhi - rdAlo;

      olapMin  = rdAlo;
      olapMax  = rdBhi;

      assert(olapMin <= olapMax);
    }
  }

  //  Decide if the repeat region is to the left or right of this read.
  //  Then ignore this thickest read if it doesn't fall in the repeat correctly.
  //
  //  A prev that ends inside the repeat, and so could be confused:
  //
  //              ----[rrrrrrrrrrrrrrrrrrrrrrrrrrrrr]------------------
  //                                   ----------------------rdA
  //                    rdB------------------
  //
  //  A prev that ends outside the repeat, can't be confused.
  //
  //              ----[rrrrrrrrrrrrrrrrrrrrrrrrrrrrr]------------------
  //                                   ----------------------rdA
  //                    rdB----------------------------
  //

  assert(rdAlo <= rMax);
  assert(rMin <= rdAlo);

  int32    rdBlo = tig->ufpath[bestIdx].position.min();
  int32    rdBhi = tig->ufpath[bestIdx].position.max();

  //  If the high coordinate isn't inside the repeat, we don't care;
  //  rdB will be anchored to unique sequence.
  //
  //  For finding the previous read:
  //
  //        --------[repeat]]]--------
  //                   -------------         READ WE'RE SEARCHING FROM
  //                 ---------u              DON'T CARE, anchored in unique (u)
  //          -----                          DON'T CARE, shouldn't happen
  //
  //  For the next read:
  //
  //        ------[[[repeat]----------
  //       -------------                     READ WE'RE SEARCHING FROM
  //             u---------                  DON'T CARE, anchored in unique (u)
  //                          ------         DON'T CARE, shouldn't happen
  //
  //  The repeat is extended a bit so that the read will have some useful bit
  //  of unique anchoring.

  if (isInside(rMin, rdBhi, rMax + 50) == false) {
    bestId  = 0;
    bestIdx = UINT32_MAX;
    bestLen = 0;
  }

  //  Log.

  else {
    writeLog("find prev OLAP from read %u position %u %u to read %u position %u %u - olap %u %u (repeat at %u %u)\n",
             rdAid,  rdAlo, rdAhi,
             bestId, rdBlo, rdBhi,
             olapMin, olapMax, rMin, rMax);
  }

  return(bestIdx);
}


uint32
findPrevBestRead(Unitig *tig, uint32 fi, int32 rMin, int32 rMax) {
  ufNode     *rdA       = &tig->ufpath[fi];
  uint32      rdAid     = rdA->ident;
  bool        rdAfwd    = rdA->position.isForward();
  int32       rdAlo     = rdA->position.min();
  int32       rdAhi     = rdA->position.max();

  //  If rdA begins before rMin, we don't need to find the previous read.
  //  it's not confused.

  if (rdAlo < rMin)
    return(UINT32_MAX);

  //  If the read is reverse, we want to return the best edge from the 3'
  //  end, After that, just check that it's in the same tig, at an
  //  overlapping position and in the correct orientation.

  BestEdgeOverlap  *best = OG->getBestEdgeOverlap(rdAid, rdA->position.isReverse());

  if (best->readId() == 0)
    return(UINT32_MAX);

  if (tig->inUnitig(best->readId()) == tig->id()) {
    uint32   pi    =  tig->ufpathIdx(best->readId());
    ufNode  *rdB   = &tig->ufpath[pi];
    int32    rdBlo = rdB->position.min();
    int32    rdBhi = rdB->position.max();

    if (isInside(rMin, rdBhi, rMax + 50) == false)    //  Don't care if the read is outside the repeat.
      return(UINT32_MAX);                             //  See above for an example.

    if ((rdBlo <= rdAlo) &&
        (rdAlo <= rdBhi) &&
        (best->read3p() == rdB->isForward())) {
      writeLog("find prev BEST from read %u position %u %u to read %u position %u %u - olap %u %u (repeat at %u %u)\n",
               rdAid,  rdAlo, rdAhi,
               rdB->ident, rdBlo, rdBhi,
               rdAlo, rdBhi, rMin, rMax);
      return(pi);
    }
  }

  return(UINT32_MAX);
}



uint32
findThickestNextRead(Unitig *tig, uint32 fi, int32 rMin, int32 rMax) {
  ufNode     *rdA       = &tig->ufpath[fi];
  uint32      rdAid     = rdA->ident;
  bool        rdAfwd    = rdA->position.isForward();
  int32       rdAlo     = rdA->position.min();
  int32       rdAhi     = rdA->position.max();

  uint32   bestId  = 0;
  uint32   bestIdx = 0;
  uint32   bestLen = 0;
  uint32   olapMin = 0;
  uint32   olapMax = 0;;

  //  If rdA ends after rMax, we don't need to find the previous read.
  //  it's not confused.

  if (rMax < rdAhi)
    return(UINT32_MAX);

  //  Otherwise, search for the previous best read.

  for (int32 pi=fi+1; pi<tig->ufpath.size(); pi++) {
    ufNode  *rdB   = &tig->ufpath[pi];
    int32    rdBlo = rdB->position.min();
    int32    rdBhi = rdB->position.max();

    if (OG->isContained(rdB->ident) == true)   //  Skip contained reads.
      continue;

    if (OG->isBackbone(rdB->ident) == false)   //  Skip non-backbone reads.
      continue;

    if (rdAhi < rdBlo)                         //  We can stop if the read
      break;                                   //  starts after us.

    if ((rdAhi >= rdBlo) && (bestLen < rdAhi - rdBlo)) {
      bestId  = rdB->ident;
      bestIdx = pi;
      bestLen = rdAhi - rdBlo;

      olapMin = rdBlo;
      olapMax = rdAhi;

      assert(olapMin <= olapMax);
    }
  }

  //  If this thickest read ends after the repeat region, it's not confused.
  //  (see above)

  assert(rMin <= rdAhi);
  assert(rdAhi <= rMax);

  int32    rdBlo = tig->ufpath[bestIdx].position.min();
  int32    rdBhi = tig->ufpath[bestIdx].position.max();

  if (isInside(rMin - 50, rdBlo, rMax) == false) {    //  Don't care if the read is outside the repeat.
    bestId  = 0;                                      //  See above for an example.
    bestIdx = UINT32_MAX;
    bestLen = 0;
  }

  //  Log.

  else {
    writeLog("find next OLAP from read %u position %u %u to read %u position %u %u - olap %u %u (repeat at %u %u)\n",
             rdAid,  rdAlo, rdAhi,
             bestId, rdBlo, rdBhi,
             olapMin, olapMax, rMin, rMax);
  }

  return(bestIdx);
}

uint32
findNextBestRead(Unitig *tig, uint32 fi, int32 rMin, int32 rMax) {
  ufNode     *rdA       = &tig->ufpath[fi];
  uint32      rdAid     = rdA->ident;
  bool        rdAfwd    = rdA->position.isForward();
  int32       rdAlo     = rdA->position.min();
  int32       rdAhi     = rdA->position.max();

  //  If rdA ends after rMax, we don't need to find the previous read.
  //  it's not confused.

  if (rMax < rdAhi)
    return(UINT32_MAX);

  //  If the read is forward, we want to return the best edge from the 3'
  //  end, After that, just check that it's in the same tig, at an
  //  overlapping position and in the correct orientation.

  BestEdgeOverlap  *best = OG->getBestEdgeOverlap(rdAid, rdA->position.isForward());

  if (best->readId() == 0)
    return(UINT32_MAX);

  if (tig->inUnitig(best->readId()) == tig->id()) {
    uint32   pi    =  tig->ufpathIdx(best->readId());
    ufNode  *rdB   = &tig->ufpath[pi];
    int32    rdBlo = rdB->position.min();
    int32    rdBhi = rdB->position.max();

    if (isInside(rMin - 50, rdBlo, rMax) == false)    //  Don't care if the read is outside the repeat.
      return(UINT32_MAX);                             //  See above for an example.

    if ((rdBlo <= rdAhi) &&
        (rdAhi <= rdBhi) &&
        (best->read3p() == rdB->isReverse())) {
      writeLog("find next BEST from read %u position %u %u to read %u position %u %u - olap %u %u (repeat at %u %u)\n",
               rdAid,  rdAlo, rdAhi,
               rdB->ident, rdBlo, rdBhi,
               rdBlo, rdAhi, rMin, rMax);
      return(pi);
    }
  }

  return(UINT32_MAX);
}



//
//  Scan all of the rdA overlaps, searching for one to rdB on the correct end of rdA.
//

class bestSco {
public:
  double  score  = 0.0;
  uint32  readId = 0;
  uint32  tigId  = 0;
};


bestSco
scoreBestOverlap(TigVector &tigs, ufNode *rdA, ufNode *rdB, bool is3p, bool internal) {
  uint32        ovlLen = 0;
  BAToverlap   *ovl    = OC->getOverlaps(rdA->ident, ovlLen);

  bestSco       bestScore;

  // if we are searching for an internal read and have been given no id to search for, it means we're not confused so return score of 0 to signal that
  if (internal == true && rdB == NULL)
     return bestScore;

  //  is3p on input is telling us if we want overlaps on the low (false) or
  //  high (true) coordinate of rdA.  If the read is in the tig flipped, we
  //  need to flip this flag so it correctly refers to the 3' end of the
  //  read.

  if (rdA->position.isReverse() == true)
    is3p = !is3p;


  for (uint32 oo=0; oo<ovlLen; oo++) {
    BAToverlap *o      = ovl + oo;
    uint32      oBid   = o->b_iid;               //  Read id of the B read in the overlap.
    uint32      oTid   = tigs.inUnitig(oBid);    //  Tig id of the B read in the overlap.

    //  For simplicity, compute the score first.

    double score  = RI->overlapLength(o->a_iid, o->b_iid, o->a_hang, o->b_hang) * (1 - o->erate());

    //  Then do a bunch of tests to ignore overlaps we don't care about.

    if ((rdB != NULL) && (rdB->ident != oBid))   //  If we're looking for a specific read,
      continue;                                  //  ignore the others.

    if (o->AEndIs3prime() != is3p)               //  Always ignore overlaps off the wrong end.
      continue;

    //  If we're looking for a specific read, this must be it.
    //  Return the score.

    if (internal == true) {
      bestScore.score  = score;
      bestScore.readId = oBid;
      bestScore.tigId  = oTid;

      return(bestScore);
    }

    //  Otherwise, we're looking for the highest scoring overlap
    //  to a read in a real tig that isn't captured in this tig.

    if ((oTid == 0) ||                          //  Skip overlaps to reads not in a tig, or
        (tigs[oTid] == NULL) ||                 //  in a singleton tig.
        (tigs[oTid]->ufpath.size() == 1))
      continue;

    if (OG->isOverlapBadQuality(ovl[oo]))       //  Ignore overlaps that aren't
      continue;                                 //  of good quality.

    if (o->isDovetail() == false)               //  Skip containment overlaps.
      continue;

    //  Also in AS_BAT_AssemblyGraph.C
    if (OG->isBubble(oBid) == true)             //  Skip overlaps to bubble tigs.
      continue;

    if (OG->isSpur(oBid) == true)               //  Skip overlaps to spur tigs.
      continue;

    if (OG->isBackbone(oBid) == false)          //  Skip overlaps to non-backbone reads.
      continue;

    //  One last test.  We need to skip overlaps to reads at this location in the tig.
    //  For this, we need to get the reads.

    uint32      tgAid  = tigs.inUnitig(o->a_iid);
    uint32      tgBid  = tigs.inUnitig(o->b_iid);

    uint32      rdBidx =  tigs[tgBid]->ufpathIdx(o->b_iid);   //  The read is in a valid tig, so
    ufNode     *rdB    = &tigs[tgBid]->ufpath[rdBidx];        //  grab all the good bits about it.

    bool        rdAfwd = rdA->position.isForward();
    int32       rdAlo  = rdA->position.min();
    int32       rdAhi  = rdA->position.max();

    bool        rdBfwd = rdB->position.isForward();
    int32       rdBlo  = rdB->position.min();
    int32       rdBhi  = rdB->position.max();

    //  If the overlap is captured in the tig, skip it.

    if ((tgBid == tgAid) &&   //  The other read is in the same tig, and
        (rdBlo <= rdAhi) &&   //  overlaps the current read.  Ignore!
        (rdAlo <= rdBhi))
      continue;

    //  Remember the highest score.

    if (bestScore.score < score) {
      bestScore.score  = score;
      bestScore.readId = o->b_iid;
      bestScore.tigId  = tgBid;
    }
  }

  return(bestScore);
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

  //  Examine every read in this tig.  If the read intersects a marked
  //  repeat, find the best edge that continues the tig in either direction.
  //
  //  If those reads are in the repeat region, scan all the overlaps of this
  //  read for any that are of comparable length.  If any are found, declare
  //  this repeat to be potentially confused.  If none are found - for the
  //  whole repeat region - then we can leave the repeat alone.

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode     *rdA       = &tig->ufpath[fi];
    uint32      rdAid     = rdA->ident;
    bool        rdAfwd    = rdA->position.isForward();
    int32       rdAlo     = rdA->position.min();
    int32       rdAhi     = rdA->position.max();

    uint32      tgAid     = tig->id();

    double      sc        = (rdAhi - rdAlo) / (double)RI->readLength(rdAid);

    if ((OG->isContained(rdAid)   == true) ||   //  Don't care about contained or chimeric
        (OG->isCoverageGap(rdAid) == true) ||      //  or non-backbone reads; we'll use the container instead.
        (OG->isBackbone(rdAid) == false)) 
     continue;

    for (uint32 ri=0; ri<tigMarksR.numberOfIntervals(); ri++) {
      int32  rMin = tigMarksR.lo(ri);
      int32  rMax = tigMarksR.hi(ri);

      if ((rdAhi < rMin) ||     //  If the read doesn't intersect a repeat region,
          (rMax  < rdAlo)) {    //  there is no chance its overlapping read will be
        continue;               //  contained in the region.  Skip!
      }

/*
      if ((rMin  <= rdAlo) &&   //  If the read is contained in the repeat region,
          (rdAhi <= rMax)) {    //  it's useless for deciding if this is a
        continue;               //  confused repeat.
      }
*/

      //  This read intersects this repeat region.  Find the
      //  reads we used to construct the tig originally.

      uint32  best5idx = findThickestPrevRead(tig, fi, rMin, rMax);
      uint32  best3idx = findThickestNextRead(tig, fi, rMin, rMax);

      //uint32  best5idx = findPrevBestRead(tig, fi, rMin, rMax);
      //uint32  best3idx = findNextBestRead(tig, fi, rMin, rMax);

      ufNode *rdB5 = (best5idx < tig->ufpath.size()) ? &tig->ufpath[best5idx] : NULL;
      ufNode *rdB3 = (best3idx < tig->ufpath.size()) ? &tig->ufpath[best3idx] : NULL;

      //  If no overlaps, we're done with this read.

      if ((rdB5 == NULL) &&
          (rdB3 == NULL))
        continue;

      //  Find and score this best overlap.

      bestSco internal5sco = scoreBestOverlap(tigs, rdA, rdB5, false,  true);
      bestSco internal3sco = scoreBestOverlap(tigs, rdA, rdB3,  true,  true);

      bestSco external5sco = scoreBestOverlap(tigs, rdA, NULL, false, false);
      bestSco external3sco = scoreBestOverlap(tigs, rdA, NULL,  true, false);

      //
      //  Here, we have a read whose thickest overlap present in the tig
      //  extends into a repeat region, and the overlap region is entirely
      //  inside the repeat.
      //

      if (external5sco.score == 0.0) {
          writeLog("tig %7u read %8u pos %7u-%-7u 5' end NOT confused -- no external edge\n",
                   tgAid, rdAid, rdAlo, rdAhi);
      }

      if ((internal5sco.score > 0.0) &&
          (external5sco.score > 0.0)) {
        double  ad5 = internal5sco.score - external5sco.score;   //  Absolute difference.
        double  pd5 = 100 * ad5 / internal5sco.score;   //  Percent diffference.

        //  This read end is confused if the internal edge is worse than the
        //  external, or if the differences are small.

        if (internal5sco.score < external5sco.score ||
          (ad5 < confusedAbsolute && pd5 < confusedPercent)) {
          writeLog("tig %7u read %8u pos %7u-%-7u 5' end  IS confused by edge to tig %8u read %8u - internal edge score %8.2f external edge score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                   tgAid, rdAid, rdAlo, rdAhi,
                   external5sco.tigId, external5sco.readId,
                   internal5sco.score, external5sco.score, ad5, pd5);

          confusedEdges.push_back(confusedEdge(rdAid, false, external5sco.readId));
          isConfused[ri]++;
        } else {
          //  There is a second best edge, and we're better than it.
          writeLog("tig %7u read %8u pos %7u-%-7u 5' end NOT confused by edge to tig %8u read %8u - internal edge score %8.2f external edge score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                   tgAid, rdAid, rdAlo, rdAhi,
                   external5sco.tigId, external5sco.readId,
                   internal5sco.score, external5sco.score, ad5, pd5);
        }
      }

      if (external3sco.score == 0.0) {
          writeLog("tig %7u read %8u pos %7u-%-7u 3' end NOT confused -- no external edge\n",
                   tgAid, rdAid, rdAlo, rdAhi);
      }

      if ((internal3sco.score > 0.0) &&
          (external3sco.score > 0.0)) {
        double  ad3 = internal3sco.score - external3sco.score;   //  Absolute difference.
        double  pd3 = 100 * ad3 / internal3sco.score;   //  Percent diffference.

        if (internal3sco.score < external3sco.score ||
          (ad3 <= confusedAbsolute && pd3 < confusedPercent)) {
          writeLog("tig %7u read %8u pos %7u-%-7u 3' end  IS confused by edge to tig %8u read %8u - internal edge score %8.2f external edge score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                   tgAid, rdAid, rdAlo, rdAhi,
                   external3sco.tigId, external3sco.readId,
                   internal3sco.score, external3sco.score, ad3, pd3);

          confusedEdges.push_back(confusedEdge(rdAid, false, external3sco.readId));
          isConfused[ri]++;
        } else {
          //  There is a second best edge, and we're better than it.
          writeLog("tig %7u read %8u pos %7u-%-7u 3' end NOT confused by edge to tig %8u read %8u - internal edge score %8.2f external edge score %8.2f - absdiff %8.2f percdiff %8.4f\n",
                   tgAid, rdAid, rdAlo, rdAhi,
                   external3sco.tigId, external3sco.readId,
                   internal3sco.score, external3sco.score, ad3, pd3);
        }
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

  writeLog("search for confused edges:\n");

  //  For each repeat region, count the number of times we find a read
  //  external to the tig with an overlap more or less of the same strength
  //  as the overlap interal to the tig.

  uint32  *isConfused = findConfusedEdges(tigs, tig, tigMarksR, confusedAbsolute, confusedPercent, confusedEdges);

  //  Scan all the regions, delete any that have no confusion.

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

  //  Remove discarded regions.

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

    if ((tig == NULL) ||                  //  Ignore deleted and singleton tigs (nothing
        (tig->ufpath.size() == 1) ||      //  to do) and unassembled reads (don't care
        (tig->_isUnassembled == true))    //  about splitting them).
      continue;

    //  Copy overlaps from the AssemblyGraph to a list of OlapDat objects,
    //  then merge overlapping ones (from the same source read) into a single
    //  record.  This is thus a list of regions on each read that potentially
    //  contain repeats.
    //
    //  Finally, project that list of intervals into tig coordinates
    //  and merge any that overlap by a significant amount.
    //
    //  The end result is to have a list of repeat regions on this tig that
    //  have full support from reads not in this tig.  If two regions overlap
    //  but only a bit, then this indicates a location where two different
    //  repeats are next to each other, but this pair of repeats occurs only
    //  in this tig.

    annotateRepeatsOnRead(AG, tig, repeatOlaps);
    mergeAnnotations(repeatOlaps, tigMarksR);

    //  Scan reads, discard any region that is well-contained in a read.
    //  When done, report the thickest overlap between any remaining region
    //  and any read in the tig.

    discardSpannedRepeats(tig, tigMarksR);
    reportThickestEdgesInRepeats(tig, tigMarksR);

    //  Sacn reads.  If a read intersects a repeat interval, and the best
    //  edge for that read is entirely in the repeat region, decide if there
    //  is a near-best edge to something not in this tig.
    //
    //  A region with no such near-best edges is _probably_ correct.

    discardUnambiguousRepeats(tigs, tig, tigMarksR, confusedAbsolute, confusedPercent, confusedEdges);

    //  Merge adjacent repeats.
    //
    //  When we split (later), we require a MIN_ANCHOR_HANG overlap to anchor
    //  a read in a unique region.  This is accomplished by extending the
    //  repeat regions on both ends.  For regions close together, this could
    //  leave a negative length unique region between them:
    //
    //   ---[-----]--[-----]---  before
    //   -[--------[]--------]-  after extending by MIN_ANCHOR_HANG (== two dashes)
    //
    //  To solve this, regions that were linked together by a single read
    //  (with sufficient overlaps to each) were merged.  However, there was
    //  no maximum imposed on the distance between the repeats, so (in
    //  theory) a 150kbp read could attach two repeats to a 149kbp unique
    //  unitig -- and label that as a repeat.  After the merges were
    //  completed, the regions were extended.
    //
    //  This version will extend regions first, then merge repeats only if
    //  they intersect.  No need for a linking read.
    //
    //  The extension also serves to clean up the edges of tigs, where the
    //  repeat doesn't quite extend to the end of the tig, leaving a few
    //  hundred bases of non-repeat.

    mergeAdjacentRegions(tig, tigMarksR);

    //  Invert.  This finds the non-repeat intervals, which get turned into
    //  non-repeat tigs.

    tigMarksU = tigMarksR;
    tigMarksU.invert(0, tig->getLength());

#if 0
    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++)
      writeLog("tigMarksR[%2u] = %d %d\n", ii, tigMarksR.lo(ii), tigMarksR.hi(ii));
    for (uint32 ii=0; ii<tigMarksU.numberOfIntervals(); ii++)
      writeLog("tigMarksU[%2u] = %d %d\n", ii, tigMarksU.lo(ii), tigMarksU.hi(ii));
#endif

    //  Create the list of intervals we'll use to make new tigs.

    vector<breakPointCoords>   BP;

    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++)
      BP.push_back(breakPointCoords(tigMarksR.lo(ii), tigMarksR.hi(ii), true));

    for (uint32 ii=0; ii<tigMarksU.numberOfIntervals(); ii++)
      BP.push_back(breakPointCoords(tigMarksU.lo(ii), tigMarksU.hi(ii), false));

    //  If there is only one BP, the tig is entirely resolved or entirely
    //  repeat.  Either case, there is nothing more for us to do.

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

    //  Scan the reads, counting the number of reads that would be placed in
    //  each new tig.  This is done because there are a few 'splits' that
    //  don't move any reads around.

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

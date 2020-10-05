
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
#include "AS_BAT_SplitTig.H"

#include "AS_BAT_FindCircular.H"

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
  olapDat(int32 b, int32 e, uint32 r, uint32 p) : tigbgn(b), tigend(e),
                                                  eviRid(r), eviPid(p) {
  };

  int32   tigbgn;   //  Location of the overlap on this tig
  int32   tigend;   //

  uint32  eviRid;   //  evidence read
  uint32  eviPid;   //  evidence read placeID
};

auto byEviRid = [](olapDat const &A, olapDat const &B) { return(((A.eviRid == B.eviRid) && (A.tigbgn < B.tigbgn)) || (A.eviRid < B.eviRid)); };
auto byCoord  = [](olapDat const &A, olapDat const &B) { return(                            A.tigbgn < B.tigbgn);                            };



//  Build a vector of olapDat (tigBgn, tigEnd, eviRid) for all reads that
//  overlap into this tig.
void
annotateRepeatsOnRead(AssemblyGraph const *AG,
                      Unitig              *tig,
                      vector<olapDat>     &repeats) {

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

  writeLog("Annotated with %lu external overlaps.\n", repeats.size());
}



void
mergeAnnotations(vector<olapDat>      &repeatOlaps,
                 intervalList<int32>  &tigMarksR) {

  //  Sort the repeat markings by their evidence read id.

  sort(repeatOlaps.begin(), repeatOlaps.end(), byEviRid);

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

  sort(repeatOlaps.begin(), repeatOlaps.end(), byEviRid);

  for (uint32 ii=repeatOlaps.size(); ii--; )
    if (repeatOlaps[ii].eviRid == UINT32_MAX)
      repeatOlaps.pop_back();

  //  Sort by coordinate.

  sort(repeatOlaps.begin(), repeatOlaps.end(), byCoord);

#ifdef SHOW_ANNOTATE
  for (uint32 ii=0; ii<repeatOlaps.size(); ii++)
    if (repeatOlaps[ii].tigbgn < 1000000)
      writeLog("repeatOlaps[%d] %d-%d read %u place %u MERGED\n",
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
  bool   foundSpanners = false;

  //  Somewhat inefficiently, for each read in the tig, check if it covers
  //  each repeat region.

  writeLog("\n");
  writeLog("Dropping repeat regions contained in a read:\n");

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
        writeLog("  region %9d-%-9d inside read %7u %9d-%d\n",
                 tigMarksR.lo(ri), tigMarksR.hi(ri), frg->ident, frglo, frghi);

        tigMarksR.lo(ri) = 0;
        tigMarksR.hi(ri) = 0;

        discarded     = true;
        foundSpanners = true;
      }
    }

    //  If any region was discarded, filter it out now, so that we don't
    //  falsely discard it again on later reads.

    if (discarded)
      tigMarksR.filterShort(1);
  }

  if (foundSpanners == false)
    writeLog("  no regions contained in a read\n");
}



//  Returns true if coord is equal to or inside the region.
bool
isInside(int32 lo, int32 coord, int32 hi) {
  return((lo <= coord) && (coord <= hi));
}


//
//  General comments on findThickestPrevRead() (similar for findThickestNextRead())
//  that were distracting if inlined in the code below.
//
//
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
//


ufNode *
findThickestPrevRead(Unitig *tig, uint32 fi, int32 rMin, int32 rMax, char *logMsg) {
  ufNode  *rdA     = &tig->ufpath[fi];
  uint32   rdAid   = rdA->ident;
  bool     rdAfwd  = rdA->position.isForward();
  int32    rdAlo   = rdA->position.min();
  int32    rdAhi   = rdA->position.max();

  uint32   bestIdx = UINT32_MAX;
  uint32   bestLen = 0;
  uint32   olapMin = 0;
  uint32   olapMax = 0;

  //  If rdA begins before the repeat, it can't be confused; return a null
  //  previous read.

  if (rdAlo < rMin)
    return(nullptr);

  //  Otherwise, search for the previous best read.

  for (int32 pi=fi-1; pi >= 0; pi--) {
    ufNode  *rdB   = &tig->ufpath[pi];
    int32    rdBlo = rdB->position.min();
    int32    rdBhi = rdB->position.max();

    if (OG->isContained(rdB->ident) == true)   //  Skip contained reads.
      continue;

    if (OG->isBackbone(rdB->ident) == false)   //  Skip non-backbone reads.
      continue;

    if ((rdBhi >= rdAlo) && (bestLen < rdBhi - rdAlo)) {
      bestIdx  = pi;
      bestLen  = rdBhi - rdAlo;

      olapMin  = rdAlo;
      olapMax  = rdBhi;

      assert(olapMin <= olapMax);
    }
  }

  //  If this thickest overlap ends outside the repeat region, the overlap
  //  cannot be confused.

  ufNode *rdB   = nullptr;
  int32   rdBlo = 0;
  int32   rdBhi = 0;

  if (bestIdx != UINT32_MAX) {
    assert(rdAlo <= rMax);
    assert(rMin <= rdAlo);

    rdB   = &tig->ufpath[bestIdx];
    rdBlo =  tig->ufpath[bestIdx].position.min();
    rdBhi =  tig->ufpath[bestIdx].position.max();

    if (isInside(rMin, rdBhi, rMax) == false)
      rdB = nullptr;
  }

  if (rdB != nullptr) {
    writeLog(logMsg);
    writeLog("    prev-olap-to read %8u %9u-%-9u olap-at %6u-%u\n", rdB->ident, rdBlo, rdBhi, olapMin, olapMax);
    logMsg[0] = 0;   //  Clear it, so we don't report it again.
  }

  return(rdB);
}



ufNode *
findThickestNextRead(Unitig *tig, uint32 fi, int32 rMin, int32 rMax, char *logMsg) {
  ufNode  *rdA     = &tig->ufpath[fi];
  uint32   rdAid   = rdA->ident;
  bool     rdAfwd  = rdA->position.isForward();
  int32    rdAlo   = rdA->position.min();
  int32    rdAhi   = rdA->position.max();

  uint32   bestIdx = UINT32_MAX;
  uint32   bestLen = 0;
  uint32   olapMin = 0;
  uint32   olapMax = 0;

  //  If rdA ends after the repeat, it can't be confused; return a null
  //  next read.

  if (rMax < rdAhi)
    return(nullptr);

  //  Otherwise, search for the next best read.

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
      bestIdx = pi;
      bestLen = rdAhi - rdBlo;

      olapMin = rdBlo;
      olapMax = rdAhi;

      assert(olapMin <= olapMax);
    }
  }

  //  If this thickest overlap begins outside the repeat region, the overlap
  //  cannot be confused.

  ufNode *rdB   = nullptr;
  int32   rdBlo = 0;
  int32   rdBhi = 0;

  if (bestIdx != UINT32_MAX) {
    assert(rMin <= rdAhi);
    assert(rdAhi <= rMax);

    rdB   = &tig->ufpath[bestIdx];
    rdBlo =  tig->ufpath[bestIdx].position.min();
    rdBhi =  tig->ufpath[bestIdx].position.max();

    if (isInside(rMin, rdBlo, rMax) == false)        //  Don't care if the read is outside the repeat.
      rdB = nullptr;                                 //  See above for an example.
  }

  if (rdB != nullptr) {
    writeLog(logMsg);
    writeLog("    next-olap-to read %8u %9u-%-9u olap-at %6u-%u\n", rdB->ident, rdBlo, rdBhi, olapMin, olapMax);
    logMsg[0] = 0;       //  Clear it, so we don't report it again.
  }

  return(rdB);
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

  //  If we are searching for an internal read and have been given no id to
  //  search for, we cannot be confused, so return empty score.

  if ((internal == true) && (rdB == NULL))
    return(bestScore);

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



bool
isCircularizingEdge(Unitig   *tig,
                    ufNode   *rdA,
                    bestSco  &internalSco,
                    bestSco  &externalSco,
                    bool      rdAhi) {

  if ((externalSco.readId == 0) ||         //  Not circularizing if no read here,
      (externalSco.tigId  != tig->id()))   //  or if it's in a different tig.
    return(false);

  if ((internalSco.readId == 0) ||         //  Not circularizing if no read here,
      (internalSco.tigId  != tig->id()))   //  or if it's in a different tig.
    return(false);

  writeLog("      checkCircular from external read %6u [%4u] to confused %s end of read %6u [%4u]\n",
           externalSco.readId, tig->ufpathIdx(externalSco.readId),
           (rdAhi) ? "hi" : "lo",
           rdA->ident, tig->ufpathIdx(rdA->ident));

  return(isCircularizingEdge(tig, externalSco.readId, internalSco.readId));
}



void
checkConfusion(uint32                rdAid,
               bestSco const         internalSco,
               bestSco const         externalSco,
               bool                  isCircular,
               char const           *end,
               bool                  endFlag,
               double                confusedAbsolute,
               double                confusedPercent,
               vector<confusedEdge> &confusedEdges) {

  //  The read end is confused if the internal edge is worse than the
  //  external, or if the differences are small.
  //
  //  We'll call it confused if there is an overlap to some other read
  //  that is comparable in length to the overlap to the read in the tig.
  //    Confused if:
  //      extOlapLen > intOlapLen - 2500   (confusedAbsolute, -ca parameter)
  //      extOlapLen > intOlapLen * 0.85   (confusedPercent,  -cp parameter)
  //
  //  (if only it was implemented that way)

  if ((internalSco.score > 0.0) &&
      (externalSco.score > 0.0)) {
    double  ad = internalSco.score - externalSco.score;   //  Absolute difference.
    double  pd = 100 * ad / internalSco.score;            //  Percent diffference.

    if (isCircular == true) {
      writeLog("    %s end NOT confused by CIRCULAR edge to read     %8u - internal edge score %8.2f external edge score %8.2f - absdiff %8.2f percdiff %8.4f\n",
               end,
               externalSco.readId,
               internalSco.score, externalSco.score, ad, pd);
    }

    else if ((internalSco.score < externalSco.score) ||
             ((ad < confusedAbsolute) &&
              (pd < confusedPercent))) {
      writeLog("    %s end  IS confused by edge to tig %8u read %8u - internal edge score %8.2f external edge score %8.2f - absdiff %8.2f percdiff %8.4f\n",
               end,
               externalSco.tigId, externalSco.readId,
               internalSco.score, externalSco.score, ad, pd);

      confusedEdges.push_back(confusedEdge(rdAid, endFlag, externalSco.readId));
    }

    else {
      writeLog("    %s end NOT confused by edge to tig %8u read %8u - internal edge score %8.2f external edge score %8.2f - absdiff %8.2f percdiff %8.4f\n",
               end,
               externalSco.tigId, externalSco.readId,
               internalSco.score, externalSco.score, ad, pd);
    }
  }
  else if ((internalSco.score >  0.0) &&
           (externalSco.score == 0.0)) {
    writeLog("    %s end NOT confused -- no external edge\n",
             end);
  }
}



void
findConfusedEdges(TigVector            &tigs,
                  Unitig                *tig,
                  intervalList<int32>  &tigMarksR,
                  double                confusedAbsolute,
                  double                confusedPercent,
                  vector<confusedEdge> &confusedEdges) {

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
        (OG->isCoverageGap(rdAid) == true) ||   //  or non-backbone reads; we'll use the container instead.
        (OG->isBackbone(rdAid)    == false)) 
     continue;

    for (uint32 ri=0; ri<tigMarksR.numberOfIntervals(); ri++) {
      int32  rMin = tigMarksR.lo(ri);
      int32  rMax = tigMarksR.hi(ri);

      if ((rdAhi < rMin) ||     //  If the read doesn't intersect a repeat region,
          (rMax  < rdAlo)) {    //  there is no chance its overlapping read will be
        continue;               //  contained in the region.  Skip!
      }

      //  This read intersects this repeat region.  Find the reads we used to
      //  construct the tig originally.
      //
      //  Generate a log that we're checking this guy for confusion.  The log
      //  will be reported only if we find a thickest edge.

      char *logMsg = new char [1024];

      snprintf(logMsg, 1024, "\nCheck read %u position %u-%u for confusion; repeat is at %u-%u.\n",
               rdAid, rdAlo, rdAhi, rMin, rMax);

      ufNode *rdB5 = findThickestPrevRead(tig, fi, rMin, rMax, logMsg);
      ufNode *rdB3 = findThickestNextRead(tig, fi, rMin, rMax, logMsg);

      delete [] logMsg;

      //  If no overlaps, we're done with this read.

      if ((rdB5 == NULL) &&
          (rdB3 == NULL))
        continue;

      //  Find and score this best overlap.

      bestSco internal5sco = scoreBestOverlap(tigs, rdA, rdB5, false,  true);
      bestSco internal3sco = scoreBestOverlap(tigs, rdA, rdB3,  true,  true);

      bestSco external5sco = scoreBestOverlap(tigs, rdA, NULL, false, false);
      bestSco external3sco = scoreBestOverlap(tigs, rdA, NULL,  true, false);

      //  Decide if the external edge looks like it will circularize the tig.
      //  If it does, we don't call the read confused.

      bool isC5 = isCircularizingEdge(tig, rdA, internal5sco, external5sco, false);
      bool isC3 = isCircularizingEdge(tig, rdA, internal5sco, external3sco,  true);

      //  Now just check confusion, write a loely log, and add a confused
      //  edge to confusedEdges.

      checkConfusion(rdAid, internal5sco, external5sco, isC5, "lo", false, confusedAbsolute, confusedPercent, confusedEdges);
      checkConfusion(rdAid, internal3sco, external3sco, isC3, "hi",  true, confusedAbsolute, confusedPercent, confusedEdges);
    }  //  Over all marks (ri)
  }  //  Over all reads (fi)
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



vector<breakReadEnd>
buildBreakPoints(TigVector             &tigs,
                 Unitig                *tig,
                 intervalList<int32>   &tigMarksR,
                 intervalList<int32>   &tigMarksU,
                 vector<confusedEdge>  &confusedEdges) {
  vector<breakReadEnd>   BE;

  //  Iterate over the two lists of regions, in coordinate order, and:
  //   - report the region.
  //   - add any confused edges in that region to the output list of breakReadEnd
  //   - fail catastrophically if there is a break in a unique region
  //  

  writeLog("\n");
  writeLog("Region summary:\n");

  for (uint32 rr=0, uu=0; ((rr < tigMarksR.numberOfIntervals()) ||
                           (uu < tigMarksU.numberOfIntervals())); ) {
    bool    isRepeat  = false;
    int32   regionBgn = 0;
    int32   regionEnd = 0;

    //  Use the repeat interval if
    //    both R and U are valid, and R is first
    //    only R is valid

    if (((rr < tigMarksR.numberOfIntervals()) && (uu <  tigMarksU.numberOfIntervals()) && (tigMarksR.lo(rr) < tigMarksU.lo(uu))) ||
        ((rr < tigMarksR.numberOfIntervals()) && (uu >= tigMarksU.numberOfIntervals()))) {
      isRepeat  = true;
      regionBgn = tigMarksR.lo(rr);
      regionEnd = tigMarksR.hi(rr);

      //writeLog("tigMarksR[%2u] = %8d,%-8d (repeat)\n", rr, tigMarksR.lo(rr), tigMarksR.hi(rr));
      rr++;
    } else {
      isRepeat  = false;
      regionBgn = tigMarksU.lo(uu);
      regionEnd = tigMarksU.hi(uu);

      //writeLog("tigMarksU[%2u] = %8d,%-8d (unique)\n", uu, tigMarksU.lo(uu), tigMarksU.hi(uu));
      uu++;
    }

    writeLog("  %s interval %8d-%-8d\n", (isRepeat == true) ? "Repeat" : "Unique", regionBgn, regionEnd);

    //  Scan all the confused edges.  Remember the extents and count how many we have.

    int32   breakCount = 0;
    int32   breakBgn   = INT32_MAX;
    int32   breakEnd   = 0;

    for (uint32 ii=0; ii<confusedEdges.size(); ii++) {
      uint32  aid    = confusedEdges[ii].aid;
      uint32  a3p    = confusedEdges[ii].a3p;
      uint32  bid    = confusedEdges[ii].bid;

      uint32  atid   =  tigs.inUnitig(aid);
      Unitig *atig   =  tigs[atid];
      ufNode *aread  = &atig->ufpath[ tigs.ufpathIdx(aid) ];
      int32   abgn   =  aread->position.bgn;
      int32   aend   =  aread->position.end;
      int32   apoint =  (a3p == false) ? aread->position.min() : aread->position.max();

      uint32  btid   =  tigs.inUnitig(bid);
      Unitig *btig   =  tigs[btid];
      int32   bbgn   =  btig->ufpath[ tigs.ufpathIdx(bid) ].position.bgn;
      int32   bend   =  btig->ufpath[ tigs.ufpathIdx(bid) ].position.end;

      if (tig->id() != atid)   //  In a different tig.  (We're keeping a list of ALL confused edges, not just for this tig)
        continue;

      //  A bit of nastiness occurs at the junction between a repeat and a
      //  unique region: we don't know if the 'apoint' belongs to the repeat
      //  or the unique region.
      //
      //  Consider apoint=12637 with these regions:
      //
      //    Repeat interval        0-12637
      //       read  382236 -> at        0-11272    hi end <-- confused by read 1514978 in tig  6142 at     3585-12353
      //       read 2019594 <- at    11776-2420     hi end <-- confused by read 1514978 in tig  6142 at     3585-12353
      //       read 1390275 -> at     2580-12637    hi end <-- confused by read  807196 in tig  6142 at     5169-13423
      //    Unique interval    12637-13259
      //
      //  Is 'apoint' the last position in the repeat, or the first position
      //  in the unique?
      //
      //  For here, a point on the boundary is declared to be only in the
      //  repeat region.  This makes no actual difference, except to avoid
      //  the assert below for such points (and to not log the point in the
      //  unique region).

      if ((isRepeat ==  true) && ((apoint <  regionBgn) ||   //  Skip break points that are
                                  (regionEnd <  apoint)))    //  entirely outside of repeat
        continue;                                            //  regions.

      if ((isRepeat == false) && ((apoint <= regionBgn) ||   //  Skip break points that are
                                  (regionEnd <= apoint)))    //  outside or on the edge of
        continue;                                            //  unique regions.

      breakCount++;
      breakBgn = min(breakBgn, apoint);
      breakEnd = max(breakEnd, apoint);

      assert(isRepeat == true);   //  No breaks in unique regions!

      //  a3p isn't indicating the oriented end of the read, but rather if
      //  the coordinate we care about is the low or the high one.
      //
      //  When a3p is true:
      //    we want to split on the higher coordinate
      //    the 'repeat' region is to the left (before) the point
      //    reads that span the point need to go into the next region
      //
      //  When a3p is false:
      //    we want to split on the lower coordinate
      //    the 'repeat' region is to the right (after) the point
      //    reads that span the point need to go into the previous region

      BE.push_back(breakReadEnd(aid, a3p,
                                apoint,
                                regionBgn, regionEnd));

      if (a3p ==  true)   assert(apoint == aread->position.max());
      if (a3p == false)   assert(apoint == aread->position.min());

      writeLog("    read %7u %s at %8d-%-8d %s end <-- confused by read %7u in tig %5u at %8d-%d\n",
               aid, aread->position.isForward() ? "->" : "<-",
               abgn, aend, (apoint == aread->position.min()) ? "lo" : "hi",
               bid, btid, bbgn, bend);
    }

    //  If in a repeat, log if there are no confused edges, and
    //  extend the region to the end of the tig, if needed.

    if (isRepeat == true) {
      if (breakCount == 0) {
        writeLog("    no confused edges\n");
      }

      else {
        if (regionBgn == 0)                  breakBgn = 0;
        if (regionEnd == tig->getLength())   breakEnd = tig->getLength();
      }
    }
  }

  return(BE);
}



void
markRepeatReads(AssemblyGraph         *AG,
                TigVector             &tigs,
                double                 deviationRepeat,
                uint32                 confusedAbsolute,
                double                 confusedPercent,
                vector<confusedEdge>  &confusedEdgesGLOBAL) {
  uint32  tiLimit = tigs.size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  writeLog("repeatDetect()-- working on " F_U32 " tigs, with " F_U32 " thread%s.\n", tiLimit, numThreads, (numThreads == 1) ? "" : "s");

  vector<olapDat>      repeatOlaps;   //  Overlaps to reads promoted to tig coords

  intervalList<int32>  tigMarksR;     //  Marked repeats based on reads, filtered by spanning reads
  intervalList<int32>  tigMarksU;     //  Non-repeat invervals, just the inversion of tigMarksR

  vector<confusedEdge> confusedEdges;

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

    writeLog("\n");
    writeLog("----------------------------------------\n");
    writeLog("Working on tig %u.\n", ti);

    annotateRepeatsOnRead(AG, tig, repeatOlaps);
    mergeAnnotations(repeatOlaps, tigMarksR);

    //  Scan reads, discard any region that is well-contained in a read.
    //  When done, report the thickest overlap between any remaining region
    //  and any read in the tig.

    discardSpannedRepeats(tig, tigMarksR);

    //  Scan reads.  If a read intersects a repeat interval, and the best
    //  edge for that read is entirely in the repeat region, decide if there
    //  is a near-best edge to something not in this tig.
    //
    //  A region with no such near-best edges is _probably_ correct.

    confusedEdges.clear();

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

    //  For each repeat region, count the number of times we find a read
    //  external to the tig with an overlap more or less of the same strength
    //  as the overlap interal to the tig.
    //
    //  Prior to mid-June 2020 this was also removing any tigMarksR that had
    //  no confused edges in them.  With the new splitting introduced around
    //  then, this had the unintended consequence of mislabeling reads as
    //  unique when no confused edge was found in a region, which could lead
    //  to new 'repeat' tigs being flagged as unique when they were actually
    //  mostly repeat, for example: -------[rrrrr]--[rrrrrrrrrr]-[rrr]------
    //  If no confused edges were found in the middle repeat block, but were
    //  in the two outer blocks, the new tig created for the middle section
    //  would be called unique, even though it was mostly repeat.

    findConfusedEdges(tigs, tig, tigMarksR, confusedAbsolute, confusedPercent, confusedEdges);

    //  Invert.  This finds the non-repeat intervals, which get turned into
    //  non-repeat tigs.

    tigMarksU = tigMarksR;
    tigMarksU.invert(0, tig->getLength());

    //  Iterate over the marked intervals, in coordinate order.  Figure out
    //  which confused edges are in the interval.  If only one, all we can do
    //  is split the tig.  If multiple, we can split the tig AND flag the
    //  resulting pieces as either repeat or unique.

    vector<breakReadEnd> BE = buildBreakPoints(tigs, tig, tigMarksR, tigMarksU, confusedEdges);

    //  If there are breaks, split the tig.

    if (BE.size() > 0) {
      splitTigAtReadEnds(tigs, tig, BE, tigMarksR);

      tigs[ti] = nullptr;
      delete tig;
    }
  }


  writeStatus("markRepeatReads()-- Found %u confused edges.\n", confusedEdges.size());
}


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
uint32 MIN_COV_REPEAT        = 27;    //  Filter coverage regions with fewer than this many overlaps.
int32  MAX_COVERAGE_GAP      = 500;   //  Size of low coverage region we can span between two high coverage regions.  MUST be signed.
uint32 MIN_ANCHOR_HANG       = 500;   //  Require reads to be anchored by this many bases at boundaries of repeats.
int32  MAX_BLOCK_GAP         = 0;     //  Merge detected repeat blocks if they're closer than this; read must span merged blocks to discard repeat

int32  REPEAT_OVERLAP_MIN    = 50;

#undef  SHOW_OLAPS
#undef  SHOW_RAW_MARKS
#undef  DO_BREAKING

#define DUMP_READ_COVERAGE
#define DUMP_ERROR_PROFILE



class breakPointCoords {
public:
  breakPointCoords(int32 tigID, int32 bgn, int32 end, bool rpt=false) {
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

  int32  _tigID;
  int32  _bgn;
  int32  _end;
  bool   _isRepeat;
};




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




class olapDat {
public:
  olapDat(uint32 b, uint32 e, double er) {
    bgn    = b;
    end    = e;
    erate  = er;
  };

  olapDat(uint32 t, uint32 r, uint32 b, uint32 e, double er) {
    tigID  = t;
    readID = r;
    bgn    = b;
    end    = e;
    erate  = er;
  };

  bool operator<(const olapDat &that)     const { return(bgn < that.bgn); };

  uint32  tigID;
  uint32  readID;

  uint32  bgn;
  uint32  end;
  double  erate;
};




//  For each overlap, if the b-read is in this tig, ignore it.
//  Otherwise annotate the read with the overlap region.
//
//  Later, check if the two reads in this unitig overlap; if not, annotate also.
//
void
annotateRepeatsOnRead(Unitig                *tgA,
                      ufNode                *rdA,
                      double                 erateRepeat,
                      vector<olapDat>       &repeats) {
  uint32               ovlLen   = 0;
  BAToverlap          *ovl      = OC->getOverlaps(rdA->ident, erateRepeat, ovlLen);

  vector<olapDat>      readOlaps;
  intervalList<int32>  readRepeats;

  uint32               tgAid    = tgA->id();

  bool                 rdAfwd   = (rdA->position.bgn < rdA->position.end);
  int32                rdAlo    = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
  int32                rdAhi    = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

  //  Beacuse the read is placed with a lot of fudging in the positions, we need
  //  to scale the coordinates we compute here.
  double               sc       = (rdAhi - rdAlo) / (double)FI->fragmentLength(rdA->ident);

  //  For all overlaps to this read, save the overlap if it is not represented in this tig.

  uint32  nOlaps = ovlLen;
  uint32  nDiff  = 0;
  uint32  nSelf  = 0;
  uint32  nConf  = 0;

  for (uint32 oi=0; oi<ovlLen; oi++) {
    uint32   rdBid  = ovl[oi].b_iid;
    uint32   tgBid  = Unitig::fragIn(rdBid);

    int32    bgn    = 0;
    int32    end    = 0;

    //  If the read is in a different tig, save the overlap.

    if (tgBid != tgAid) {
      nDiff++;
      olapToReadCoords(rdA, ovl+oi, bgn, end);
    }

    //  Otherwise, the read is in the same tig.  If it doesn't intersect us, save the overlap.

    else {
      uint32   rdBpos =  Unitig::pathPosition(rdBid);
      ufNode  *rdB    = &tgA->ufpath[rdBpos];

      bool     rdBfwd   = (rdB->position.bgn < rdB->position.end);
      int32    rdBlo    = (rdBfwd) ? rdB->position.bgn : rdB->position.end;
      int32    rdBhi    = (rdBfwd) ? rdB->position.end : rdB->position.bgn;

      if ((rdAhi < rdBlo) || (rdBhi < rdAlo)) {
        nSelf++;
        olapToReadCoords(rdA, ovl+oi, bgn, end);
      } else {
        nConf++;
      }
    }

    //  If bgn != end, then we need to save this overlap.  bgn and end are the position in
    //  the read where the overlap hit.

    if (bgn != end) {
      int32  tigbgn = (rdAfwd) ? (rdAlo + sc * bgn) : (rdAhi - sc * end);
      int32  tigend = (rdAfwd) ? (rdAlo + sc * end) : (rdAhi - sc * bgn);

      assert(tigbgn < tigend);

      if (tigbgn < 0)                  tigbgn = 0;
      if (tigend > tgA->getLength())   tigend = tgA->getLength();

      fprintf(stderr, "tig %u read %u OVERLAP tig %u read %u at tigpos %u %u erate %f\n",
              tgAid, rdA->ident, tgBid, rdBid, tigbgn, tigend, ovl[oi].erate);

      readOlaps.push_back(olapDat(tgBid, rdBid, tigbgn, tigend, ovl[oi].erate));
      readRepeats.add(tigbgn, tigend - tigbgn);
    }
  }

  //  All overlaps processed.

  //  If there are no intervals saved, we're done.

  if (readRepeats.numberOfIntervals() == 0)
    return;

  //  If the number of overlaps from reads outside this tig is less than half the number of overlaps
  //  from reads inside this tig, we _probably_ got this correct.

  if (nDiff + nSelf < nConf / 2)
    return;

  //  Otherwise, squish the intervals together....

  readRepeats.merge();

  fprintf(stderr, "tig %u read %u nOlaps %u nDiff %u nSelf %u nConf %u -- %u regions %u-%u %u overlaps\n",
          tgAid, rdA->ident,
          nOlaps, nDiff, nSelf, nConf,
          readRepeats.numberOfIntervals(), readRepeats.lo(0), readRepeats.hi(0), readRepeats.count(0));

  //  ....and save overlaps for those intervals formed from a bunch of overlaps.
     
  for (uint32 rr=0; rr<readRepeats.numberOfIntervals(); rr++) {
    if (readRepeats.count(rr) < 5)
      continue;

    for (uint32 oo=0; oo<readOlaps.size(); oo++) {
      if ((readRepeats.hi(rr) < readOlaps[oo].bgn) ||
          (readOlaps[oo].end < readRepeats.lo(rr)))
        //  Repeat interval before this overlap, or 
        //  overlap before this repeat interval, skip.
        continue;

      //  Overlap intersects this interval, save.

      repeats.push_back(readOlaps[oo]);
    }
  }
}




void
annotateRepeatsInReads(UnitigVector          &unitigs, 
                       Unitig                *tig,
                       double                 erateRepeat,
                       intervalList<int32>   &tigMarksR) {     //  Output list of repeat regions in the unitig

  tigMarksR.clear();

  //  For each read, add overlaps from reads not in this unitig to the 'repeats' vector.  After all
  //  are added, it is sorted, and scanned against the errorProfile to weed out diverged repeats.
  vector<olapDat>  repeatOlaps;

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++)
    annotateRepeatsOnRead(tig, &tig->ufpath[fi], erateRepeat, repeatOlaps);

//
//  PROBLEM - we're filtering spurious based on unfiltered overlaps.  Diverged overlaps
//  should be removed THEN those with only a few overlaps can be ignored.
//

  fprintf(stderr, "Annotated with %lu overlaps.\n", repeatOlaps.size());

  //  Now that we have all the external reads overlapped, filter out diverged repeats.
  //  Any overlap that passes the filter is moved to a list of repeat intervals (tigMarksR).

  sort(repeatOlaps.begin(), repeatOlaps.end());

  for (uint32 bb=0, ii=0; ii<repeatOlaps.size(); ii++) {
    uint32  nBelow = 0;
    uint32  nAbove = 0;

    //  Find first region for this overlap.  All later overlaps will start at or after this region.
    for (; (bb < tig->errorProfile.size()) && (tig->errorProfile[bb].bgn < repeatOlaps[ii].bgn); bb++)
      ;

    //  Until the overlap stops overlapping, add its error rate to the regions.
    for (uint32 ee=bb; (ee < tig->errorProfile.size()) && (tig->errorProfile[ee].end <= repeatOlaps[ii].end); ee++)
#warning size compute is wrong
      if (repeatOlaps[ii].erate <= tig->errorProfile[ee].dev.mean() + 5 * tig->errorProfile[ee].dev.stddev())
        nBelow += tig->errorProfile[ee].end - tig->errorProfile[ee].bgn;
      else
        nAbove += tig->errorProfile[ee].end - tig->errorProfile[ee].bgn;

    if (nAbove > nBelow) {
      //writeLog("SKIP overlap iid %u iid %u at %.4f%% at unitig position %u %u - %u bases high error\n",
      //         0, 0, repeatOlaps[ii].erate, repeatOlaps[ii].bgn, repeatOlaps[ii].end, nAbove);
      continue;
    }

    writeLog("SAVE overlap to tig %u iid %u at %.4f%% at unitig position %u %u - %u bases high error\n",
             repeatOlaps[ii].tigID,
             repeatOlaps[ii].readID,
             repeatOlaps[ii].erate, repeatOlaps[ii].bgn, repeatOlaps[ii].end, nAbove);

    tigMarksR.add(repeatOlaps[ii].bgn, repeatOlaps[ii].end - repeatOlaps[ii].bgn);
  }
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





void
markRepeatReads(UnitigVector &unitigs,
                double        UNUSED(erateGraph),
                double        UNUSED(erateBubble),
                double        UNUSED(erateMerge),
                double        erateRepeat) {
  uint32  tiLimit = unitigs.size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  writeLog("repeatDetect()-- working on "F_U32" unitigs, with "F_U32" threads.\n", tiLimit, numThreads);

  intervalList<int32>  tigMarksR;  //  Marked repeats based on reads, filtered by spanning reads
  intervalList<int32>  tigMarksU;  //  Non-repeat invervals, just the inversion of tigMarksR

  //  Compute all error profiles.

  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = unitigs[ti];

    if (tig == NULL)
      continue;

    if (tig->ufpath.size() == 1)
      continue;

    fprintf(stderr, "Finding coverage and error profile for tig %u/%u.\n", ti, tiLimit);

    tig->computeErrorProfile();
  }



  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = unitigs[ti];

    if (tig == NULL)
      continue;

    if (tig->ufpath.size() == 1)
      continue;

    vector<olapDat>   repeats;

    fprintf(stderr, "Annotating repeats in reads for tig %u/%u.\n", ti, tiLimit);

    annotateRepeatsInReads(unitigs, tig, erateRepeat, tigMarksR);

    //  Collapse all the (scaled) read markings to intervals on the unitig, merging those that overlap
    //  significantly.

    fprintf(stderr, "Merge marks.\n");

#if 1
    writeLog("PRE-MERGE markings\n");
    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++) {
      writeLog("  %8d:%-8d size %7d (distance to next %7d)\n",
               tigMarksR.lo(ii), tigMarksR.hi(ii), tigMarksR.hi(ii) - tigMarksR.lo(ii),
               (ii < tigMarksR.numberOfIntervals()-1) ? (tigMarksR.lo(ii+1) - tigMarksR.hi(ii)) : (0));
    }
#endif

    tigMarksR.merge(REPEAT_OVERLAP_MIN);

#if 1
    writeLog("POST-MERGE markings\n");
    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++) {
      writeLog("  %8d:%-8d size %7d (distance to next %7d)\n",
               tigMarksR.lo(ii), tigMarksR.hi(ii), tigMarksR.hi(ii) - tigMarksR.lo(ii),
               (ii < tigMarksR.numberOfIntervals()-1) ? (tigMarksR.lo(ii+1) - tigMarksR.hi(ii)) : (0));
    }
#endif

    //  Scan reads, discard any mark that is contained in a read

    fprintf(stderr, "Scan reads.\n");

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode     *frg      = &tig->ufpath[fi];
      uint32      ovlLen   =  0;
      BAToverlap *ovl      =  OC->getOverlaps(frg->ident, erateRepeat, ovlLen);

      bool        frgfwd   = (frg->position.bgn < frg->position.end);
      int32       frglo    = (frgfwd) ? frg->position.bgn : frg->position.end;
      int32       frghi    = (frgfwd) ? frg->position.end : frg->position.bgn;

      for (uint32 ri=0; ri<tigMarksR.numberOfIntervals(); ri++)
        if ((frglo + MIN_ANCHOR_HANG <= tigMarksR.lo(ri)) && (tigMarksR.hi(ri) + MIN_ANCHOR_HANG <= frghi)) {
          writeLog("discard region %8d:%-8d - contained in read %6u %5d-%5d\n",
                   tigMarksR.lo(ri), tigMarksR.hi(ri), frg->ident, frglo, frghi);
          tigMarksR.lo(ri) = 0;
          tigMarksR.hi(ri) = 0;
        }

      tigMarksR.filterShort(1);
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

    fprintf(stderr, "Make breakpoints.\n");

    vector<breakPointCoords>   BP;

    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++)
      BP.push_back(breakPointCoords(ti,
                                    tigMarksR.lo(ii) - MIN_ANCHOR_HANG,
                                    tigMarksR.hi(ii) + MIN_ANCHOR_HANG,
                                    true));

    for (uint32 ii=0; ii<tigMarksU.numberOfIntervals(); ii++)
      BP.push_back(breakPointCoords(ti,
                                    tigMarksU.lo(ii) + MIN_ANCHOR_HANG,
                                    tigMarksU.hi(ii) - MIN_ANCHOR_HANG,
                                    false));

    //  If only one region, the whole unitig was declared repeat.  Nothing to do.

    if (BP.size() == 1)
      continue;

    sort(BP.begin(), BP.end());

    //  Report.

    fprintf(stderr, "Report.\n");

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


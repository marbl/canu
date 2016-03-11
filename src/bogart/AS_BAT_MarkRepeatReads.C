
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

//  Hack.
uint32 MIN_COV_REPEAT        = 27;    //  Filter coverage regions with fewer than this many overlaps.
int32  MAX_COVERAGE_GAP      = 500;   //  Size of low coverage region we can span between two high coverage regions.  MUST be signed.
uint32 MIN_ANCHOR_HANG       = 500;   //  Require reads to be anchored by this many bases at boundaries of repeats.
int32  MAX_BLOCK_GAP         = 0;     //  Merge detected repeat blocks if they're closer than this; read must span merged blocks to discard repeat

#undef  SHOW_OLAPS
#define SHOW_RAW_MARKS
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





class errorProfileVector {
  struct erp {
    double    sn;  //  Current variance
    double    mn;
    uint32    n;
  };

public:
  errorProfileVector(uint32 len) {
    erps = new erp [len];
    memset(erps, 0, sizeof(erp) * len);
  };

  ~errorProfileVector() {
    delete [] erps;
  };

  void   add(uint32 pos, double val) {
    double  m0 = erps[pos].mn;
    double  s0 = erps[pos].sn;
    uint32  n  = erps[pos].n + 1;

    erps[pos].mn = m0 + (val - m0) / n;
    erps[pos].sn = s0 + (val - m0) * (val - erps[pos].mn);
    erps[pos].n  = n;
  };

  double variance(uint32 pos) {
    if (erps[pos].n > 1)
      return(erps[pos].sn / (erps[pos].n - 1));
    return(0);
  };
  double stddev  (uint32 pos) {  return(sqrt(variance(pos)));   };

  double mean    (uint32 pos) {  return(erps[pos].mn);          };

private:
  erp    *erps;
};




void
findUnitigErrorProfile(UnitigVector                  &unitigs,
                       Unitig                        *tig,
                       errorProfileVector            &profile) {

  uint32      ti       = tig->id();

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode     *frg    = &tig->ufpath[fi];
    bool        frgfwd   = (frg->position.bgn < frg->position.end);
    int32       frglo    = (frgfwd) ? frg->position.bgn : frg->position.end;
    int32       frghi    = (frgfwd) ? frg->position.end : frg->position.bgn;

    uint32      ovlLen =  0;
    BAToverlap *ovl    =  OC->getOverlaps(frg->ident, AS_MAX_ERATE, ovlLen);

    fprintf(stderr, " profile - read %u/%u\r", fi, tig->ufpath.size());

    for (uint32 oi=0; oi<ovlLen; oi++) {
      uint32   tib  =  Unitig::fragIn(ovl[oi].b_iid);
      uint32   fib  =  Unitig::pathPosition(ovl[oi].b_iid);
      Unitig  *tigb =  unitigs[tib];
      ufNode  *frgb = &tigb->ufpath[fib];

      //  Skip if the two reads are in the different unitigs, or they don't overlap.

      if (ti != tib)
        continue;

      int32    lo  = frglo;
      int32    hi  = frghi;

      int32    lob = min(frgb->position.bgn, frgb->position.end);
      int32    hib = max(frgb->position.bgn, frgb->position.end);

      if ((hi < lob) || (hib < lo))
        continue;

      //  Now figure out what region is covered by the overlap.

      olapToReadCoords(frg, ovl+oi, lo, hi);

      //  Offset and adjust to tig coordinates

      //  Beacuse the read is placed with a lot of fudging in the positions, we need
      //  to scale the coordinates we compute here.
      double      sc = (frghi - frglo) / (double)FI->fragmentLength(frg->ident);

      lo = frglo + sc * lo;
      hi = frglo + sc * hi;

      //  Save the interval

      //rawcoverage.add(lo, hi-lo, ovl[oi].erate);

      for (uint32 ii=lo; ii<hi; ii++)
        profile.add(ii, ovl[oi].erate);
    }
  }

  fprintf(stderr, "\n");

  //  All reads and all overlaps added.  Squish down to coverage 

#ifdef DUMP_ERROR_PROFILE
  char  fn[FILENAME_MAX];
  sprintf(fn, "%08u.errors", tig->id());
  FILE *F = fopen(fn, "w");

  for (uint32 ii=0; ii<tig->getLength(); ii++)
    fprintf(F, "%u %f %f\n", ii, profile.mean(ii), profile.stddev(ii));

  fclose(F);
#endif
}







//  For each overlap, if the b-read is in this tig, ignore it.
//  Otherwise annotate the read with the overlap region.
//
//  Later, check if the two reads in this unitig overlap; if not, annotate also.
//
void
annotateRepeatsOnRead(UnitigVector          &unitigs,
                      Unitig                *tig,
                      ufNode                *frg,
                      double                 erateRepeat,
                      errorProfileVector    &errorProfile,
                      intervalList<int32>   &readMarks) {
  uint32      ovlLen   =  0;
  BAToverlap *ovl      =  OC->getOverlaps(frg->ident, erateRepeat, ovlLen);

  uint32      ti       = tig->id();

  bool        frgfwd   = (frg->position.bgn < frg->position.end);
  int32       frglo    = (frgfwd) ? frg->position.bgn : frg->position.end;
  int32       frghi    = (frgfwd) ? frg->position.end : frg->position.bgn;

  for (uint32 oi=0; oi<ovlLen; oi++) {
    uint32   tib  =  Unitig::fragIn(ovl[oi].b_iid);
    uint32   fib  =  Unitig::pathPosition(ovl[oi].b_iid);
    Unitig  *tigb =  unitigs[tib];
    ufNode  *frgb = &tigb->ufpath[fib];

    //  Skip if the overlap is garbage

    if (ovl[oi].erate > erateRepeat)
      continue;

    //  Skip if the two reads are in the same tig and they overlap

    int32    lo  = frglo;
    int32    hi  = frghi;

    int32    lob = min(frgb->position.bgn, frgb->position.end);
    int32    hib = max(frgb->position.bgn, frgb->position.end);

    if ((ti == tib) &&
        ((lo < hib) && (lob < hi)))
      continue;

    //  Now figure out what region is covered by the overlap.

    olapToReadCoords(frg, ovl+oi, lo, hi);

    //  Compare the overlap error against the error rates for this region.

    uint32  nBelow = 0;
    uint32  nAbove = 0;

    for (uint32 ii=lo; ii<hi; ii++)
      if (ovl[oi].erate < errorProfile.mean(ii) + 5 * errorProfile.stddev(ii))
        nBelow++;
      else
        nAbove++;

    if (nAbove > nBelow)
      writeLog("SKIP overlap iid %u iid %u at %.4f%% at unitig position %u %u - %u bases high error\n",
               ovl[oi].a_iid, ovl[oi].b_iid, ovl[oi].erate, lo, hi, nAbove);

    readMarks.add(lo, hi-lo);
  }
}





void
annotateRepeatsInReads(UnitigVector          &unitigs, 
                       Unitig                *tig,
                       double                 erateRepeat,
                       errorProfileVector    &errorProfile,
                       intervalList<int32>   &tigMarksR) {

  tigMarksR.clear();

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode     *frg      = &tig->ufpath[fi];

    bool        frgfwd   = (frg->position.bgn < frg->position.end);
    int32       frglo    = (frgfwd) ? frg->position.bgn : frg->position.end;
    int32       frghi    = (frgfwd) ? frg->position.end : frg->position.bgn;

    //  Beacuse the read is placed with a lot of fudging in the positions, we need
    //  to scale the coordinates we compute here.
    double      sc = (frghi - frglo) / (double)FI->fragmentLength(frg->ident);

    //  A list of the positions on this read with overlaps elsewhere.  One entry per overlap.
    intervalList<int32>  readOlaps;

    annotateRepeatsOnRead(unitigs, tig, frg, erateRepeat, errorProfile, readOlaps);

    uint32  nOlaps = readOlaps.numberOfIntervals();

    if (nOlaps == 0)
      continue;

    //  Compute the coverage intervals, and then squish the overlaps down to actual intervals of coverage.

    intervalList<int32>  readDepth(readOlaps);

    //  Filter out low coverage intervals.

    for (uint32 ii=0; ii<readDepth.numberOfIntervals(); ii++) {
      if (readDepth.depth(ii) < MIN_COV_REPEAT) {
        //writeLog("remove read %u at %d %d depth %d\n",
        //         frg->ident, readDepth.lo(ii), readDepth.hi(ii), readDepth.depth(ii));
        readDepth.lo(ii) = 0;
        readDepth.hi(ii) = 0;
      } else {
        //writeLog("retain read %u at %d %d depth %d\n",
        //         frg->ident, readDepth.lo(ii), readDepth.hi(ii), readDepth.depth(ii));
      }
    }
    readDepth.filterShort(1);

    if (readDepth.numberOfIntervals() == 0)
      continue;

    //  Merge the depths into markings, ignoring depth and merging those that are close.

    readDepth.merge(-MAX_COVERAGE_GAP);

    uint32  nRegions = readDepth.numberOfIntervals();

#ifdef SHOW_RAW_MARKS
    writeLog("tig %7u read %9u has %3u repeat markings from %3u overlaps:",
             tig->id(), frg->ident, nRegions, nOlaps);
#endif

    //  Map the (scaled) read markings to their (approximate) unitig position.

    for (uint32 ii=0; ii<nRegions; ii++) {
      int32  tigbgn = (frgfwd) ? (frglo + sc * readDepth.lo(ii)) : (frghi - sc * readDepth.hi(ii));
      int32  tigend = (frgfwd) ? (frglo + sc * readDepth.hi(ii)) : (frghi - sc * readDepth.lo(ii));

      assert(tigbgn < tigend);

      if (tigbgn < 0)                  tigbgn = 0;
      if (tigend > tig->getLength())   tigend = tig->getLength();

#ifdef SHOW_RAW_MARKS
      writeLog(" %6d-%-6d %8d:%-8d", readDepth.lo(ii), readDepth.hi(ii), tigbgn, tigend);
#endif

      tigMarksR.add(tigbgn, tigend - tigbgn);
    }

#ifdef SHOW_RAW_MARKS
    writeLog("\n");
#endif
  }  //  Annotate repeats in all reads.
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

  intervalList<int32>  tigMarksR;
  intervalList<int32>  tigMarksU;

  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = unitigs[ti];

    if (tig == NULL)
      continue;

    if (tig->ufpath.size() == 1)
      continue;

    fprintf(stderr, "Finding coverage and error profile for tig %u/%u.\n", ti, tiLimit);

    intervalList<uint32>         readCoverage;
    errorProfileVector           errorProfile(tig->getLength());

    fprintf(stderr, " coverage\n");
    findUnitigCoverage(tig, readCoverage);
    fprintf(stderr, " profile\n");
    findUnitigErrorProfile(unitigs, tig, errorProfile);
    fprintf(stderr, " continue\n");

    fprintf(stderr, "Annotating repeats in reads for tig %u/%u.\n", ti, tiLimit);

    //  Algorithm:
    //    map external olaps to read
    //    convert to coverage intervals
    //    drop any intervals less than X mappings - HARD, need heuristic to not drop low cov between high cov?
    //    convert to tig coords, add to tigMarksR
    //    collapse tigMarksR to intervals
    //
    annotateRepeatsInReads(unitigs, tig, erateRepeat, errorProfile, tigMarksR);

    //  Collapse all the (scaled) read markings to intervals on the unitig.

#if 1
    writeLog("PRE-MERGE markings\n");
    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++) {
      writeLog("  %8d:%-8d size %7d (distance to next %7d)\n",
               tigMarksR.lo(ii), tigMarksR.hi(ii), tigMarksR.hi(ii) - tigMarksR.lo(ii),
               (ii < tigMarksR.numberOfIntervals()-1) ? (tigMarksR.lo(ii+1) - tigMarksR.hi(ii)) : (0));
    }
#endif

    tigMarksR.merge();

#if 1
    writeLog("POST-MERGE markings\n");
    for (uint32 ii=0; ii<tigMarksR.numberOfIntervals(); ii++) {
      writeLog("  %8d:%-8d size %7d (distance to next %7d)\n",
               tigMarksR.lo(ii), tigMarksR.hi(ii), tigMarksR.hi(ii) - tigMarksR.lo(ii),
               (ii < tigMarksR.numberOfIntervals()-1) ? (tigMarksR.lo(ii+1) - tigMarksR.hi(ii)) : (0));
    }
#endif

    //  Scan reads, discard any mark that is contained in a read

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

    writeLog("break tig %u into up to %u pieces:\n", ti, BP.size());
    for (uint32 ii=0; ii<BP.size(); ii++)
      writeLog("  %8d %8d %s (length %d)\n",
               BP[ii]._bgn, BP[ii]._end,
               BP[ii]._isRepeat ? "repeat" : "unique",
               BP[ii]._end, BP[ii]._bgn);

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


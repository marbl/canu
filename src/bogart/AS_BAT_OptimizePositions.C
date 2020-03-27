
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

#include "runtime.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"
#include "stddev.H"
#include <vector>

//  History of this, as "remembered" by BPW in Feb 2020.
//
//  This was added 21 June 2017.  One week before this, I was working on GFA
//  and unitig creation, a pretty strong hint I was finding inconsistencies
//  in the layout when trying to split tigs.  After optimize was added,
//  mostly minor bug fixes, but overlap symmetry was added (meaning that A->B
//  was loaded but B->A was not).  So it's not terribly clear what propted
//  this to be added.
//
//  The original commit also only optimized after initial tigs were created
//  (step 003), not after contained reads are added.  Optimize after
//  contained reads are placed was added very quickly after this.  Most
//  of our failures are now after contained reads are placed.
//
//  My guess is that I was finding inconsistent locations to split a contig
//  into unitigs based on read-to-read overlaps, and that adding optimize was
//  enough to make contig-to-unitig work without additional changes.
//
//  Some commits
//    5ea89885c6e56db006ff3ae7154207c957ffc761 21 Jun 2017
//      Optimize read layout after creating initial contigs.
//    2077807173394cef13dbf8650ad9c58a1d915ea4  3 Jul 2017
//      Add a second pass for optimize() to place reads with missing overlaps.
//    47baa1d48e785d39ba0ef78ede0a8038314a2e97  7 Jul 2017
//      Drop optimize() iterations from 25 to 5, expand tigs so that short
//      reads aren't short any more.
//    f8bb08d9141f3a5c8b96c0ecac32507ad60919ef 17 Jul 2017
//      Move optimize() to own file.  (it was in AS_BAT_Unitig.C)
//    (I stopped here)

class optPos {
public:
  optPos() {
    ident = 0;
    min   = 0.0;
    max   = 0.0;
    fwd   = false;
  };

  void    set(ufNode &n) {
    ident = n.ident;
    min   = n.position.min();
    max   = n.position.max();
    fwd   = n.position.isForward();
  };

  uint32  ident;
  int32   min;
  int32   max;
  bool    fwd;
};



//  Decide if two reads in a tig are compatible with an overlap.
//
//  Returns true if the reads are overlapping in the tig, and are oriented
//  consistent with the overlap.
//
bool
Unitig::optimize_isCompatible(uint32       ii,
                              uint32       jj,
                              BAToverlap  &olap,
                              bool         inInit,
                              bool         firstPass,
                              bool         beVerbose) {

  SeqInterval  &ip  = ufpath[ii].position;   //  We're testing the ii'th read, to see if the
  SeqInterval  &jp  = ufpath[jj].position;   //  overlap to the jj'th is valid for this tig.

  uint32        iiid = ufpath[ii].ident;
  uint32        jjid = ufpath[jj].ident;

  assert(olap.a_iid == iiid);                //    The overlap has a_iid == ii.
  assert(olap.b_iid == jjid);

  assert(RI->readLength(iiid) > 0);
  assert(RI->readLength(jjid) > 0);

  //  We test four things:
  //
  //    Are the two reads overlapping in the current layout?  If not, then
  //    this overlap isn't correct for this positioning.  This can happen
  //    if the overlap orientation is flipped.
  //
  //    Does this overlap preserve the relative order of the two reads, e.g.,
  //    that read A starts before read B.  If not, the overlap is likely
  //    to reads at different loactions in the tig.
  //
  //    Does this overlap keep the read at _roughly_ the same location in
  //    the layout?  If not, it is definitely not the correct overlap.
  //
  //    Does the overlap keep the read in the same orientation?
  //
  bool  isOvl        = isOverlapping(ufpath[ii].position, ufpath[jj].position);
  bool  isOrdered    = true;
  bool  isPositioned = true;
  bool  isOriented   = true;

  //  Decide if the overlap preserves the ordering in the tig.
  //
  //  This is of questionable value, and might be incorrect.  We can have two
  //  overlaps between a pair of reads iff the orientation differs - the
  //  flipped flag.  That gets tested elsewhere.
  //
  //  isOvl tells the end of the read that should have the overlap (according
  //  to the tig layout), so check if the overlap actually is on that end.
  //  If not, flag this as 'mis-ordered'.
  //
  //  But if in a containment relationship, it cannot be mis-ordered but we rely on the next check instead.
  //
#if 1
  bool  isOvlLo = ((jp.min() <= ip.min()) && (ip.min() <= jp.max()) && (jp.max() <= ip.max()));
  bool  isOvlHi = ((ip.min() <= jp.min()) && (jp.min() <= ip.max()) && (ip.max() <= jp.max()));

  if (((isOvlLo == true) && (ip.isForward()) && (olap.AEndIs5prime() == false)) ||
      ((isOvlLo == true) && (ip.isReverse()) && (olap.AEndIs3prime() == false)) ||
      ((isOvlHi == true) && (ip.isForward()) && (olap.AEndIs3prime() == false)) ||
      ((isOvlHi == true) && (ip.isReverse()) && (olap.AEndIs5prime() == false)))
    isOrdered = false;

  if ((olap.AisContained() == true) ||
      (olap.AisContainer() == true))
    isOrdered = false;

#endif

  //  If the positions _roughly_ agree with the positions expected from the
  //  overlap, return true.
  //
  //  The overlap is from the ii read.  If that's forward, the hangs apply
  //  as-is.  If not, the hangs are swapped and inverted.
  //
  //       --------------> B          -B  <-------------
  //        A  --------------->     <------------  -A
  //
#if 1
  int32 expJbgn = ip.min() + (ip.isForward() ? olap.a_hang : -olap.b_hang);
  int32 expJend = ip.max() + (ip.isForward() ? olap.b_hang : -olap.a_hang);

  double  JbgnDiff = fabs(expJbgn - jp.min()) / (double)RI->readLength(jjid);
  double  JendDiff = fabs(expJend - jp.max()) / (double)RI->readLength(jjid);

  if ((JbgnDiff > 0.1) ||
      (JendDiff > 0.1))
    isPositioned = false;
#endif

  //  Are the reads in the layout in the same orientation as implied by the overlap?
  //
  if (((ip.isForward() == jp.isForward()) && (olap.flipped == true)) ||
      ((ip.isForward() != jp.isForward()) && (olap.flipped == false)))
    isOriented = false;


  //  Finally, dump some logging if in verbose mode or if in the second pass.
  //
  if (((beVerbose == true) || (firstPass == false)) &&
      (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
    writeLog("optimize_%s()-- tig %7u read %9u (at %9d %9d) olap to read %9u (at %9d %9d) - hangs %7d %7d - %s %s %s ovlLo %d ovlHi %d position %.4f %.4f flags %d%d%d%d%d\n",
             (inInit) ? "initPlace" : "recompute",
             id(),
             ufpath[ii].ident, ip.bgn, ip.end,
             ufpath[jj].ident, jp.bgn, jp.end,
             olap.a_hang, olap.b_hang,
             (isOvl        == true) ? "overlapping" : "not-overlapping",
             (isOrdered    == true) ? "ordered"     : "mis-ordered",
             (isPositioned == true) ? "positioned"  : "mis-positioned",
             isOvlLo, isOvlHi,
             JbgnDiff, JendDiff,
             ip.isForward(), olap.AisContained(), olap.AisContainer(), olap.AEndIs5prime(), olap.AEndIs3prime());

  return((isOvl         == true) &&    //  Good if the reads overlap in the current layout.
         ((isOrdered    == true) ||    //  Good if the reads are in the order implied by the overlap, OR
          (isPositioned == true)) &&   //       if the overlapping read is close enough.
         (isOriented    == true));     //  Good if the reads are in the orientation implied by the overlap.
}



void
Unitig::optimize_initPlace(uint32        ii,
                           optPos       *op,
                           optPos       *np,
                           bool          firstPass,
                           set<uint32>  &failed,
                           bool          beVerbose) {
  uint32   iid  = ufpath[ii].ident;
  vector<int32> hs;

  if ((firstPass == false) && (failed.count(iid) == 0))  //  If the second pass and not
    return;                                              //  failed, do nothing.

  if ((firstPass == false) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
    writeLog("optimize_initPlace()-- Second pass begins.\n");

  //  Then process all overlaps.

  if (ii > 0) {
    uint32       ovlLen  = 0;
    BAToverlap  *ovl     = OC->getOverlaps(iid, ovlLen);

    for (uint32 oo=0; oo<ovlLen; oo++) {
      uint32  jid = ovl[oo].b_iid;
      uint32  uu  = inUnitig (jid);
      uint32  jj  = ufpathIdx(jid);

      //  Skip if:
      //    the overlap is to a read in a different tig (note ufpath[jj] is invalid).
      //    the overlap is to a read after us in the layout, and we're in the first pass.
      //    the overlap isn't compatible with the current layout.

      if (uu != id())
        continue;

      if ((firstPass == true) && (ii < jj))
        continue;

      if (optimize_isCompatible(ii, jj, ovl[oo], true, !firstPass, beVerbose) == false)
        continue;

      //  Compute the position of the read using the overlap and the other read.
      int32 tmin = (op[iid].fwd) ? (op[jid].min - ovl[oo].a_hang) : (op[jid].min + ovl[oo].b_hang);

      if ((logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
         writeLog("optimize_initPlace()-- tig %7u read %9u placed at %d with overlap to %d\n", id(), iid, tmin, jid);
      hs.push_back(tmin);
    }  //  over all overlaps

    //  If no overlaps found, flag this read for a second pass.  If in the second pass,
    //  not much we can do.

    if ((firstPass == true) && (hs.size() == 0)) {
      if (logFileFlagSet(LOG_OPTIMIZE_POSITIONS))
        writeLog("optimize_initPlace()-- tig %7u read %9u FAILED TO FIND OVERLAPS (first pass)\n",
                 id(), iid);
      failed.insert(iid);
      return;
    }

    if ((firstPass == false) && (hs.size() == 0) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS))) {
      writeLog("optimize_initPlace()-- tig %7u read %9u FAILED TO FIND OVERLAPS (second pass)\n",
               id(), iid);
      flushLog();
    }

    assert(hs.size() > 0);
  }

  //  The initialization above does very little to enforce read lengths, and the optimization
  //  doesn't put enough weight in the read length to make it stable.  We simply force
  //  the correct read length here.

  op[iid].min = 0;
  if (hs.size() != 0)
     computeMedian(hs, op[iid].min);
  op[iid].max = op[iid].min + RI->readLength(ufpath[ii].ident);

  np[iid].min = 0;
  np[iid].max = 0;

  if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
    writeLog("optimize_initPlace()-- tig %7u read %9u initialized to position %7d %7d %s\n",
             id(), op[iid].ident, op[iid].min, op[iid].max, (firstPass == true) ? "" : " SECONDPASS");
}



void
Unitig::optimize_recompute(uint32        iid,
                           optPos       *op,
                           optPos       *np,
                           bool          beVerbose) {
  uint32       ii      = ufpathIdx(iid);

  int32        readLen = RI->readLength(iid);

  uint32       ovlLen  = 0;
  BAToverlap  *ovl     = OC->getOverlaps(iid, ovlLen);

  vector<int32> hsmin;
  vector<int32> hsmax;

  if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS))) {
    writeLog("optimize()--\n");
    writeLog("optimize()-- tig %8u read %9u previous  - %9d-%-9d length %9u\n", id(), iid, op[iid].min,           op[iid].max,           readLen);
    writeLog("optimize()-- tig %8u read %9u length    - %9d-%-9d length %9u\n", id(),  iid, op[iid].max - readLen, op[iid].min + readLen, readLen);
  }

  //  Process all overlaps.

  for (uint32 oo=0; oo<ovlLen; oo++) {
    uint32  jid = ovl[oo].b_iid;
    assert(RI->readLength(jid) > 0);
    double  jratio = (double)(op[jid].max-op[jid].min) / RI->readLength(jid);
    uint32  uu  = inUnitig (jid);
    uint32  jj  = ufpathIdx(jid);

    if (uu != id())   //  Skip if to a different tig.
      continue;       //  (otherwise, ufpath[jj] is invalid below)

    //  Compute the position of the read using the overlap and the other read.

    int32 tmin = (op[iid].fwd) ? (op[jid].min - (int32)(floor(ovl[oo].a_hang*jratio))) : (op[jid].min + (int32)(floor(ovl[oo].b_hang*jratio)));
    int32 tmax = (op[iid].fwd) ? (op[jid].max - (int32)(floor(ovl[oo].b_hang*jratio))) : (op[jid].max + (int32)(floor(ovl[oo].a_hang*jratio)));

    //  Skip if the overlap isn't compatible with the layout.
    if (optimize_isCompatible(ii, jj, ovl[oo], false, false, beVerbose) == false)
      continue;

   // if someone tried to position me to be negative vote for staying alone
   if (tmax - tmin <= 0) {
      writeLog("optimize()-- tig %8u read %9u to %d olap with hangs %d and %d currently at position %d %d set me to %d %d\n", id(), iid, jid, ovl[oo].a_hang, ovl[oo].b_hang, op[jid].min, op[jid].max, tmin, tmax);
      tmax=tmin+readLen;
    }
    assert(tmax - tmin >= 0);

    if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
      writeLog("optimize()-- tig %8u read %9u to %d olap %4u - %9d-%9d\n", id(), iid, jid, oo, tmin, tmax);

    //  Update estimate.

    hsmin.push_back(tmin);
    hsmax.push_back(tmax);
  }

  if (hsmin.size() == 0) {
    writeLog("Failed to optimize read %u in tig %u\n", iid, id());
    fprintf(stderr, "Failed to optimize read %u in tig %u\n", iid, id());
    flushLog();
  }
  assert(hsmin.size() > 0);

  //  Add in some evidence for the bases in the read.  We want higher weight than the overlaps,
  //  but not enough to swamp the hangs.

  hsmin.push_back(op[iid].max - readLen); hsmin.push_back(op[iid].max - readLen);
  hsmax.push_back(op[iid].min + readLen); hsmax.push_back(op[iid].min + readLen);

  //  Find the average and save.

  computeMedian(hsmin, np[iid].min);
  computeMedian(hsmax, np[iid].max);

  if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS))) {
    double dmin = 2.0 * (op[iid].min - np[iid].min) / (op[iid].min + np[iid].min);
    double dmax = 2.0 * (op[iid].max - np[iid].max) / (op[iid].max + np[iid].max);
    double npll = (double)np[iid].max - np[iid].min;

    writeLog("optimize()-- tig %8u read %9u           - %9d-%-9d length %9.2f/%-6d %7.2f%% posChange %+6.4f %+6.4f\n",
             id(), iid,
             np[iid].min, np[iid].max,
             npll, readLen,
             200.0 * (npll - readLen) / (npll + readLen),
             dmin, dmax);
  }
}



void
Unitig::optimize_expand(optPos  *op,
                       bool beVerbose) {

  for (uint32 ii=0; ii<ufpath.size(); ii++) {
    uint32       iid     = ufpath[ii].ident;

    int32        readLen = RI->readLength(iid);

    int32       opiimin = op[iid].min;             //  New start of this read, same as the old start
    int32       opiimax = op[iid].min + readLen;   //  New end of this read
    int32       opiilen = op[iid].max - op[iid].min;

    if (readLen <= opiilen)   //  This read is sufficiently long,
      continue;               //  do nothing.

    if (readLen <= 0 || opiilen <= 0)
       fprintf(stderr, "Error: tig %d negative read length %9u at positions %9d - %9d, should be %9u bases\n", id(), iid, op[iid].min, op[iid].max, readLen);
    assert(readLen > 0);
    assert(opiilen > 0);

    double       scale   = (double)readLen / opiilen;
    int32        expand  = opiimax - op[iid].max;   //  Amount we changed this read, bases
    if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
       writeLog("optimize_expand() -- read %9u at positions %9d - %9d length %9u and currently %9u scale %4.2f expanding by %9u bp\n", iid, op[iid].min, op[iid].max, readLen, opiilen, scale, expand);

    //  For each read, adjust positions based on how much they overlap with this read.

    for (uint32 jj=0; jj<ufpath.size(); jj++) {
      uint32 jid = ufpath[jj].ident;
      if (jid == iid) continue;		// don't update ourselves in the middle of the run, that happens at the end

      if      (op[jid].min < op[iid].min)
        ;
      else if (op[jid].min < op[iid].max) {
        uint32 old= op[jid].min;
        op[jid].min  = opiimin + (int32)(ceil((op[jid].min - op[iid].min) * scale));
        if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
          writeLog("optimize_expand1() -- read %9u being updated by %9u (%9d - %9d) min from %9d to %9d now it is at %9d - %9d and length %9u vs %9u\n", jid, iid, opiimin, opiimax, old, op[jid].min, op[jid].min, op[jid].max, op[jid].max-op[jid].min, RI->readLength(jid));
      } else {
         uint32 old= op[jid].min;
         if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
           writeLog("optimize_expand2() -- read %9u being updated by %9u (%9d - %9d) min from %9d to %9d now it is at %9d - %9d and length %9u vs %9u\n", jid, iid, opiimin, opiimax, old, op[jid].min, op[jid].min, op[jid].max, op[jid].max-op[jid].min, RI->readLength(jid));
        op[jid].min += expand;
      }

      if      (op[jid].max < op[iid].min)
        ;
      else if (op[jid].max < op[iid].max) {
        uint32 old= op[jid].max;
        op[jid].max  = opiimin + (int32)(ceil((op[jid].max - op[iid].min) * scale));
        if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
          writeLog("optimize_expand3() -- read %9u being updated by %9u (%9d - %9d) max from %9d to %9d now it is at %9d - %9d and length %9u vs %9u\n", jid, iid, opiimin, opiimax, old, op[jid].max, op[jid].min, op[jid].max, op[jid].max-op[jid].min, RI->readLength(jid));
      } else {
        uint32 old= op[jid].max;
        op[jid].max += expand;
        if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
          writeLog("optimize_expand4() -- read %9u being updated by %9u (%9d - %9d) max from %9d to %9d now it is at %9d - %9d and length %9u vs %9u\n", jid, iid, opiimin, opiimax, old, op[jid].max, op[jid].min, op[jid].max, op[jid].max-op[jid].min, RI->readLength(jid));
      }
    }

    //  Finally, actually shift us
    if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
      writeLog("optimize_expand() -- read %9u updated position from %9d - %9d to %9d - %9d\n", iid, op[iid].min, op[iid].max, opiimin, opiimax);

    op[iid].min = opiimin;
    op[iid].max = opiimax;
  }
}



void
Unitig::optimize_setPositions(optPos  *op,
                              bool     beVerbose) {

  for (uint32 ii=0; ii<ufpath.size(); ii++) {
    uint32  iid     = ufpath[ii].ident;

    int32   readLen = RI->readLength(iid);
    int32   opll    = (int32)op[iid].max - (int32)op[iid].min;
    double  opdd    = 200.0 * (opll - readLen) / (opll + readLen);

    if (op[iid].fwd) {
      if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
        writeLog("optimize()-- read %9u -> from %9d,%-9d %7d to %9d,%-9d %7d readLen %7d diff %7.4f%%\n",
                 iid,
                 ufpath[ii].position.bgn,
                 ufpath[ii].position.end,
                 ufpath[ii].position.end - ufpath[ii].position.bgn,
                 (int32)op[iid].min,
                 (int32)op[iid].max,
                 opll,
                 readLen,
                 opdd);

      ufpath[ii].position.bgn = (int32)op[iid].min;
      ufpath[ii].position.end = (int32)op[iid].max;
    } else {
      if ((beVerbose) && (logFileFlagSet(LOG_OPTIMIZE_POSITIONS)))
        writeLog("optimize()-- read %9u <- from %9d,%-9d %7d to %9d,%-9d %7d readLen %7d diff %7.4f%%\n",
                 iid,
                 ufpath[ii].position.bgn,
                 ufpath[ii].position.end,
                 ufpath[ii].position.bgn - ufpath[ii].position.end,
                 (int32)op[iid].max,
                 (int32)op[iid].min,
                 opll,
                 readLen,
                 opdd);

      ufpath[ii].position.bgn = (int32)op[iid].max;
      ufpath[ii].position.end = (int32)op[iid].min;
    }
  }
}



void
TigVector::optimizePositions(const char *prefix, const char *label) {
  uint32  numThreads  = omp_get_max_threads();

  uint32  tiLimit     = size();
  uint32  tiBlockSize = 10; //(tiLimit <   10 * numThreads) ? numThreads : tiLimit / 9;

  uint32  fiLimit     = RI->numReads() + 1;
  uint32  fiBlockSize = 100; //(fiLimit < 1000 * numThreads) ? numThreads : fiLimit / 999;

  bool    beVerbose   = false;

  writeStatus("optimizePositions()-- Optimizing read positions for %u reads in %u tigs, with %u thread%s.\n",
              fiLimit, tiLimit, numThreads, (numThreads == 1) ? "" : "s");

  //  Create work space and initialize to current read positions.

  writeStatus("optimizePositions()--   Allocating scratch space for %u reads (%u KB).\n", fiLimit, sizeof(optPos) * fiLimit * 2 / 1024);

  optPos *pp = NULL;
  optPos *op = new optPos [fiLimit];
  optPos *np = new optPos [fiLimit];

  for (uint32 fi=0; fi<fiLimit; fi++) {
    uint32    ti = inUnitig(fi);
    uint32    pp = ufpathIdx(fi);
    Unitig   *tig = operator[](ti);

    if (tig == NULL)
      continue;

    op[fi].set(tig->ufpath[pp]);
    np[fi].set(tig->ufpath[pp]);
  }

  //  Compute initial positions using previously placed reads and the read length.


  //
  //  Initialize positions using only reads before us.  If any reads fail to find overlaps, a second
  //  round will init positions using any read (before or after).
  //

  writeStatus("optimizePositions()--   Initializing positions with %u threads.\n", numThreads);

#pragma omp parallel for schedule(dynamic, tiBlockSize)
  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig       *tig = operator[](ti);
    set<uint32>   failed;

    if ((tig == NULL) || (tig->ufpath.size() == 1))
      continue;

    for (uint32 ii=0; ii<tig->ufpath.size(); ii++)
      tig->optimize_initPlace(ii, op, np, true,  failed, beVerbose);

    for (uint32 ii=0; ii<tig->ufpath.size(); ii++)
      tig->optimize_initPlace(ii, op, np, false, failed, beVerbose);
  }

  //
  //  Recompute positions using all overlaps and reads both before and after.  Do this for a handful of iterations
  //  so it somewhat stabilizes.
  //

  for (uint32 iter=0; iter<5; iter++) {

    //  Recompute positions

    writeStatus("optimizePositions()--   Recomputing positions, iteration %u, with %u threads.\n", iter+1, numThreads);

#pragma omp parallel for schedule(dynamic, fiBlockSize)
    for (uint32 fi=0; fi<fiLimit; fi++) {
      uint32        ti = inUnitig(fi);
      Unitig       *tig = operator[](ti);

      if ((tig == NULL) || (tig->ufpath.size() == 1))
        continue;

      tig->optimize_recompute(fi, op, np, beVerbose);
    }

    //  Reset zero

    writeStatus("optimizePositions()--     Reset zero.\n");

    for (uint32 ti=0; ti<tiLimit; ti++) {
      Unitig       *tig = operator[](ti);

      if ((tig == NULL) || (tig->ufpath.size() == 1))
        continue;

      int32  z = np[ tig->ufpath[0].ident ].min;

      for (uint32 ii=0; ii<tig->ufpath.size(); ii++) {
        uint32  iid = tig->ufpath[ii].ident;

        np[iid].min -= z;
        np[iid].max -= z;
      }
    }

    //  Decide if we've converged.  We used to compute percent difference in coordinates, but that is
    //  biased by the position of the read.  Just use percent difference from read length.

    writeStatus("optimizePositions()--     Checking convergence.\n");

    uint32  nConverged = 0;
    uint32  nChanged   = 0;

    for (uint32 fi=0; fi<fiLimit; fi++) {
      double  minp = 2.0 * (op[fi].min - np[fi].min) / (RI->readLength(fi));
      double  maxp = 2.0 * (op[fi].max - np[fi].max) / (RI->readLength(fi));

      if (minp < 0)  minp = -minp;
      if (maxp < 0)  maxp = -maxp;

      if ((minp < 0.005) && (maxp < 0.005))
        nConverged++;
      else
        nChanged++;
    }

    //  All reads processed, swap op and np for the next iteration.

    pp = op;
    op = np;
    np = pp;

    writeStatus("optimizePositions()--     converged: %6u reads\n", nConverged);
    writeStatus("optimizePositions()--     changed:   %6u reads\n", nChanged);

    if (nChanged == 0)
      break;
  }

  //
  //  Reset small reads.  If we've placed a read too small, expand it (and all reads that overlap)
  //  to make the length not smaller.
  //

  writeStatus("optimizePositions()--   Expanding short reads with %u threads.\n", numThreads);

#pragma omp parallel for schedule(dynamic, tiBlockSize)
  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig       *tig = operator[](ti);

    if ((tig == NULL) || (tig->ufpath.size() == 1))
      continue;

    tig->optimize_expand(op, beVerbose);
  }

  //
  //  Update the tig with new positions.  op[] is the result of the last iteration.
  //

  writeStatus("optimizePositions()--   Updating positions.\n");

  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig       *tig = operator[](ti);

    if ((tig == NULL) || (tig->ufpath.size() == 1))
      continue;

    tig->optimize_setPositions(op, beVerbose);
    tig->cleanUp();
  }

  //  Cleanup and finish.

  delete [] op;
  delete [] np;

  writeStatus("optimizePositions()--   Finished.\n");
}

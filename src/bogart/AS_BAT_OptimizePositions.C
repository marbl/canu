
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
 *    src/bogart/AS_BAT_Unitig.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2017-JUL-17
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"



class optPos {
public:
  optPos() {
  };
  ~optPos() {
  };

  void    set(ufNode &n) {
    ident = n.ident;
    min   = n.position.min();
    max   = n.position.max();
    fwd   = n.position.isForward();
  };

  uint32  ident;
  double  min;
  double  max;
  bool    fwd;
};



void
Unitig::optimize_initPlace(uint32        ii,
                           optPos       *op,
                           optPos       *np,
                           bool          firstPass,
                           set<uint32>  &failed,
                           bool          beVerbose) {
  uint32   fid  = ufpath[ii].ident;
  double   nmin = 0;
  int32    cnt  = 0;

  if ((firstPass == false) && (failed.count(fid) == 0))  //  If the second pass and not
    return;                                              //  failed, do nothing.

  //  Then process all overlaps.

  if (ii > 0) {
    uint32       ovlLen  = 0;
    BAToverlap  *ovl     = OC->getOverlaps(fid, ovlLen);

    for (uint32 oo=0; oo<ovlLen; oo++) {
      uint32  uu = inUnitig (ovl[oo].b_iid);
      uint32  jj = ufpathIdx(ovl[oo].b_iid);

      //  Probably overkill, but report ALL overlaps for the troubling reads.

      if ((beVerbose) || (firstPass == false))
        writeLog("olap %u a %u b %u hangs %d %d\n", oo, ovl[oo].a_iid, ovl[oo].b_iid, ovl[oo].a_hang, ovl[oo].b_hang);

      if (uu != id())   //  Skip if the overlap is to a different tig.
        continue;       //  (the ufpathIdx() call is valid, but using it isn't)

      //  Reads are in the same tig.  Decide if they overlap in position.

      bool  isOvl = isOverlapping(ufpath[ii].position, ufpath[jj].position);

      //  Log!  beVerbose should be true for the second pass, but just in case it isn't.

      if ((beVerbose) || (firstPass == false))
        writeLog("optimize_initPlace()-- olap %4u tig %7u read %8u (at %9d %9d) olap to read %8u (at %9d %9d) - hangs %7d %7d - %s %s\n",
                 oo, id(),
                 ovl[oo].a_iid, ufpath[ii].position.bgn, ufpath[ii].position.end,
                 ovl[oo].b_iid, ufpath[jj].position.bgn, ufpath[jj].position.end,
                 ovl[oo].a_hang, ovl[oo].b_hang,
                 (isOvl == true) ? "overlapping" : "not-overlapping",
                 (jj > ii)       ? "after"       : "before");

      if (isOvl == false)            //  Skip if the reads
        continue;                    //  don't overlap

      if ((firstPass) && (jj > ii))  //  We're setting initial positions, so overlaps to reads after
        continue;                    //  us aren't correct, unless we're in the 2nd pass

      //  Reads overlap.  Compute the position of the read using
      //  the overlap and the other read.

      nmin += (op[ii].fwd) ? (op[jj].min - ovl[oo].a_hang) : (op[jj].min + ovl[oo].b_hang);
      cnt  += 1;
    }  //  over all overlaps

    //  If no overlaps found, flag this read for a second pass.  If in the second pass,
    //  not much we can do.

    if ((firstPass == true) && (cnt == 0)) {
      writeLog("Failed to find overlaps for read %u in tig %u at %d-%d (first pass)\n",
               fid, id(), ufpath[ii].position.bgn, ufpath[ii].position.end);
      failed.insert(fid);
      return;
    }

    if ((firstPass == false) && (cnt == 0)) {
      writeLog("Failed to find overlaps for read %u in tig %u at %d-%d (second pass)\n",
               fid, id(), ufpath[ii].position.bgn, ufpath[ii].position.end);
      flushLog();
    }

    assert(cnt > 0);
  }

  //  The initialization above does very little to enforce read lengths, and the optimization
  //  doesn't put enough weight in the read length to make it stable.  We simply force
  //  the correct read length here.

  op[ii].min = (cnt == 0) ? 0 : (nmin / cnt);
  op[ii].max = op[ii].min + RI->readLength(ufpath[ii].ident);

  np[ii].min = 0;
  np[ii].max = 0;

  if (beVerbose)
    writeLog("optimize_initPlace()-- tig %7u read %9u initialized to position %9.2f %9.2f%s\n",
             id(), op[ii].ident,
             op[ii].min,
             op[ii].max,
             (firstPass == true) ? "" : " SECONDPASS");
}



void
Unitig::optimize(const char *prefix, const char *label) {
  bool   beVerbose = false;

  //if (id() == 43408)
  //  beVerbose = true;

  if (ufpath.size() == 1)
    return;

  if (beVerbose) {
    writeLog("\n");
    writeLog("optimize()-- tig %u\n", id());
    writeLog("\n");
  }

  optPos *pp = NULL;
  optPos *op = new optPos [ufpath.size()];
  optPos *np = new optPos [ufpath.size()];

  memset(op, 0, sizeof(optPos) * ufpath.size());
  memset(np, 0, sizeof(optPos) * ufpath.size());

  for (uint32 ii=0; ii<ufpath.size(); ii++) {
    op[ii].set(ufpath[ii]);
    np[ii].set(ufpath[ii]);
  }

  //
  //  Initialize positions using only reads before us.
  //  If any reads fail that, a second round will init positions using any read.
  //

  set<uint32>   failed;

  for (uint32 ii=0; ii<ufpath.size(); ii++)
    optimize_initPlace(ii, op, np, true,  failed, beVerbose);

  for (uint32 ii=0; ii<ufpath.size(); ii++)
    optimize_initPlace(ii, op, np, false, failed, true);

  //
  //  Optimize
  //

  for (uint32 iter=0; iter<5; iter++) {
    uint64  nOlapsTotal = 0;
    uint64  nOlapsUsed  = 0;

    //  Now run through all reads and compute positions based on overlaps.

    for (uint32 ii=0; ii<ufpath.size(); ii++) {
      uint32       fid     = op[ii].ident;

      int32        readLen = RI->readLength(fid);

      uint32       ovlLen  = 0;
      BAToverlap  *ovl     = OC->getOverlaps(fid, ovlLen);

      double       nmin = 0.0;
      double       nmax = 0.0;
      uint32       cnt  = 0;

      if (beVerbose) {
        writeLog("optimize()-- tig %8u read %8u previous  - %9.2f-%-9.2f\n", id(), fid, op[ii].min,           op[ii].max);
        writeLog("optimize()-- tig %8u read %8u length    - %9.2f-%-9.2f\n", id(), fid, op[ii].max - readLen, op[ii].min + readLen);
      }

      //  Process all overlaps.

      for (uint32 oo=0; oo<ovlLen; oo++) {
        uint32  uu = inUnitig (ovl[oo].b_iid);
        uint32  jj = ufpathIdx(ovl[oo].b_iid);

        if (uu != id())   //  Skip if the overlap is to a different tig.
          continue;       //  (the ufpathIdx() call is valid, but using it isn't)

        if (isOverlapping(ufpath[ii].position, ufpath[jj].position) == false)  //  Skip if the reads
          continue;                                                            //  don't overlap

        //  Reads overlap.  Compute the position of the read using
        //  the overlap and the other read.

        double tmin = (op[ii].fwd) ? (op[jj].min - ovl[oo].a_hang) : (op[jj].min + ovl[oo].b_hang);
        double tmax = (op[ii].fwd) ? (op[jj].max - ovl[oo].b_hang) : (op[jj].max + ovl[oo].a_hang);

        if (beVerbose)
          writeLog("optimize()-- tig %8u read %8u olap %4u - %9.2f-%-9.2f\n", id(), fid, oo, tmin, tmax);

        nmin += tmin;
        nmax += tmax;
        cnt  += 1;
      }  //  over all overlaps

      nOlapsTotal  += ovlLen;   //  Save some stats.
      nOlapsUsed   += cnt;

      //  Add in some evidence for the bases in the read.  We want higher weight than the overlaps,
      //  but not enough to swamp the hangs.

      nmin   += cnt/4 * (op[ii].max - readLen);
      nmax   += cnt/4 * (op[ii].min + readLen);
      cnt    += cnt/4;

      //  Find the average and save.

      np[ii].min = nmin / cnt;
      np[ii].max = nmax / cnt;

      double dmin = 2 * (op[ii].min - np[ii].min) / (op[ii].min + np[ii].min);
      double dmax = 2 * (op[ii].max - np[ii].max) / (op[ii].max + np[ii].max);
      double npll = np[ii].max - np[ii].min;

      if (beVerbose)
        writeLog("optimize()-- tig %8u read %8u           - %9.2f-%-9.2f length %9.2f/%-6d %7.2f%% posChange %+6.4f %+6.4f iter %2u\n",
                 id(), fid,
                 np[ii].min, np[ii].max,
                 npll, readLen,
                 200.0 * (npll - readLen) / (npll + readLen),
                 dmin, dmax,
                 iter);
    }  //  Over all reads in the tig

    //  Offset back to zero.

    int32  z = np[0].min;

    for (uint32 ii=0; ii<ufpath.size(); ii++) {
      np[ii].min -= z;
      np[ii].max -= z;
    }

    //  Decide if we've converged.  We used to compute percent difference in coordinates, but that is
    //  biased by the position of the read.  Just use percent difference from read length.

    uint32  nConverged = 0;
    uint32  nChanged   = 0;

    for (uint32 ii=0; ii<ufpath.size(); ii++) {
      double  minp = 2 * (op[ii].min - np[ii].min) / (RI->readLength(ufpath[ii].ident));
      double  maxp = 2 * (op[ii].max - np[ii].max) / (RI->readLength(ufpath[ii].ident));

      if (minp < 0)  minp = -minp;
      if (maxp < 0)  maxp = -maxp;

      if ((minp < 0.005) && (maxp < 0.005))
        nConverged++;
      else
        nChanged++;
    }

    //  All reads processed, swap op and np to compute the next iteration.

    pp = op;
    op = np;
    np = pp;

    if (beVerbose)
      writeLog("optimize()-- tig %8u iter %2u converged %6u changed %6u  olaps %9lu/%9lu\n",
               id(), iter, nConverged, nChanged, nOlapsUsed, nOlapsTotal);

    if (nChanged == 0)
      break;
  }

  //  If we've placed a read too small, expand it (and all reads that overlap) to make
  //  the length not smaller.

  for (uint32 ii=0; ii<ufpath.size(); ii++) {
    uint32       fid     = op[ii].ident;
    int32        readLen = RI->readLength(fid);

    double       opiimin = op[ii].min;             //  New start of this read, same as the old start
    double       opiimax = op[ii].min + readLen;   //  New end of this read
    double       opiilen = op[ii].max - op[ii].min;

    if (readLen <= opiilen)   //  This read is sufficiently long,
      continue;               //  do nothing.

    double       scale   = readLen / opiilen;
    double       expand  = opiimax - op[ii].max;   //  Amount we changed this read, bases

    //  For each read, adjust positions based on how much they overlap with this read.

    for (uint32 jj=0; jj<ufpath.size(); jj++) {

      if      (op[jj].min < op[ii].min)
        ;
      else if (op[jj].min < op[ii].max)
        op[jj].min  = opiimin + (op[jj].min - op[ii].min) * scale;
      else
        op[jj].min += expand;


      if      (op[jj].max < op[ii].min)
        ;
      else if (op[jj].max < op[ii].max)
        op[jj].max  = opiimin + (op[jj].max - op[ii].min) * scale;
      else
        op[jj].max += expand;
    }

    //  Finally, actually shift us

    op[ii].min = opiimin;
    op[ii].max = opiimax;
  }

  //  Update the tig with new positions.  op[] is the result of the last iteration.

  for (uint32 ii=0; ii<ufpath.size(); ii++) {
    int32   readLen = RI->readLength(ufpath[ii].ident);
    int32   opll    = (int32)op[ii].max - (int32)op[ii].min;
    double  opdd    = 200.0 * (opll - readLen) / (opll + readLen);

    if (op[ii].fwd) {
      if (beVerbose)
        writeLog("optimize()-- read %8u -> from %9d,%-9d %7d to %9d,%-9d %7d readLen %7d diff %7.4f%%\n",
                 ufpath[ii].ident,
                 ufpath[ii].position.bgn,
                 ufpath[ii].position.end,
                 ufpath[ii].position.end - ufpath[ii].position.bgn,
                 (int32)op[ii].min,
                 (int32)op[ii].max,
                 opll,
                 readLen,
                 opdd);

      ufpath[ii].position.bgn = (int32)op[ii].min;
      ufpath[ii].position.end = (int32)op[ii].max;
    } else {
      if (beVerbose)
        writeLog("optimize()-- read %8u <- from %9d,%-9d %7d to %9d,%-9d %7d readLen %7d diff %7.4f%%\n",
                 ufpath[ii].ident,
                 ufpath[ii].position.bgn,
                 ufpath[ii].position.end,
                 ufpath[ii].position.bgn - ufpath[ii].position.end,
                 (int32)op[ii].max,
                 (int32)op[ii].min,
                 opll,
                 readLen,
                 opdd);

      ufpath[ii].position.bgn = (int32)op[ii].max;
      ufpath[ii].position.end = (int32)op[ii].min;
    }
  }

  delete [] op;
  delete [] np;

  cleanUp();
}



void
TigVector::optimizePositions(const char *prefix, const char *label) {
  uint32  tiLimit = size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  writeStatus("optimizePositions()-- Optimizing read positions for %u tigs, with %u thread%s.\n", tiLimit, numThreads, (numThreads == 1) ? "" : "s");

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = operator[](ti);

    if (tig == NULL)
      continue;

    tig->optimize(prefix, label);
  }
}


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
 *    src/AS_CNS/tigStore.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2009-OCT-05 to 2014-MAR-31
 *      are Copyright 2009-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Michael Schatz from 2010-FEB-27 to 2010-MAY-17
 *      are Copyright 2010 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2014-APR-13
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2015-AUG-14
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-22
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"

#include "AS_UTL_decodeRange.H"
#include "intervalList.H"

#include "tgTigSizeAnalysis.H"

#undef  DEBUG_IGNORE

#define DUMP_UNSET               0
#define DUMP_STATUS              1
#define DUMP_TIGS                2
#define DUMP_CONSENSUS           3
#define DUMP_LAYOUT              4
#define DUMP_MULTIALIGN          5
#define DUMP_SIZES               6
#define DUMP_COVERAGE            7
#define DUMP_DEPTH_HISTOGRAM     8
#define DUMP_THIN_OVERLAP        9
#define DUMP_OVERLAP_HISTOGRAM  10


//  positions can be
//     layout (if no consensusExists()
//     gapped
//     ungapped
//
//
//  MISSING
//    masking and splitting consensus based on (low) coverage
//    coverage
//      single unitig coverage plot and dump of low coverage regions
//      all unitig histogram of coverage
//    reporting gc content (need to do this after consensus, while sequence is still in memory, then save in the tig itself)

class tgFilter {
public:
  tgFilter() {
    tigIDbgn        = 0;
    tigIDend        = UINT32_MAX;

    dumpAllClasses  = true;
    dumpUnassembled = false;
    dumpBubbles     = false;
    dumpContigs     = false;

    minNreads       = 0;
    maxNreads       = UINT32_MAX;

    minLength       = 0;
    maxLength       = UINT32_MAX;

    minCoverage     = 0.0;
    maxCoverage     = DBL_MAX;

    minGoodCov      = 0.0;
    maxGoodCov      = DBL_MAX;

    IL              = NULL;
    ID              = NULL;
  };

  ~tgFilter() {
    delete IL;
    delete ID;
  };

  bool          ignore(tgTig *tig, bool useGapped) {
#ifdef DEBUG_IGNORE
    bool   iI = ignoreID(tig);
    bool   iN = ignoreNreads(tig);
    bool   iL = ignoreLength(tig, useGapped);
    bool   iC = ignoreCoverage(tig, useGapped);
    bool   iS = ignoreClass(tig);

    fprintf(stderr, "ignore()--  tig %u - ignore id %s Nreads %s length %s coverage %s class %s\n",
            tig->tigID(),
            (iI) ? "true" : "false",
            (iN) ? "true" : "false",
            (iL) ? "true" : "false",
            (iC) ? "true" : "false",
            (iS) ? "true" : "false");
#endif

    return(ignoreID(tig) ||
           ignoreNreads(tig) ||
           ignoreLength(tig, useGapped) ||
           ignoreCoverage(tig, useGapped) ||
           ignoreClass(tig));
  };

  bool          ignoreID(tgTig *tig) {
    return((tig->tigID() < tigIDbgn) ||
           (tigIDend < tig->tigID()));
  };

  bool          ignoreClass(tgTig *tig) {
    if ((dumpAllClasses == true) ||
        ((tig->_class == tgTig_unassembled) && (dumpUnassembled == true)) ||
        ((tig->_class == tgTig_bubble)      && (dumpBubbles == true)) ||
        ((tig->_class == tgTig_contig)      && (dumpContigs == true)))
      return(false);

    return(true);
  };

  bool          ignoreNreads(tgTig *tig) {
    return((tig->numberOfChildren() < minNreads) ||
           (maxNreads < tig->numberOfChildren()));
  };

  bool          ignoreLength(tgTig *tig, bool useGapped) {
    uint32  length = tig->length(useGapped);

    return((length < minLength) ||
           (maxLength < length));
  };

  bool          ignoreCoverage(tgTig *tig, bool useGapped) {
    if ((minCoverage == 0) && (maxCoverage == UINT32_MAX))
      return(false);

    if (tig->consensusExists() == false)
      useGapped = true;

    delete IL;
    IL = new intervalList<int32>;

    for (uint32 i=0; i<tig->numberOfChildren(); i++) {
      tgPosition *pos = tig->getChild(i);

      int32  bgn = (useGapped) ? pos->min() : tig->mapGappedToUngapped(pos->min());
      int32  end = (useGapped) ? pos->max() : tig->mapGappedToUngapped(pos->max());

      IL->add(bgn, end - bgn);
    }

    delete ID;
    ID = new intervalList<int32>(*IL);

    uint32  goodCov = 0;
    uint32  badCov  = 0;

    for (uint32 ii=0; ii<ID->numberOfIntervals(); ii++)
      if ((minCoverage  <= ID->depth(ii)) &&
          (ID->depth(ii) <= maxCoverage))
        goodCov += ID->hi(ii) - ID->lo(ii);
      else
        badCov += ID->hi(ii) - ID->lo(ii);

    double fracGood = (goodCov) / (goodCov + badCov);

    return((fracGood < minGoodCov) ||
           (maxGoodCov < fracGood));
  };

  void          maskConsensus(char *cns) {
    for (uint32 ii=0; ii<ID->numberOfIntervals(); ii++) {
      if ((minCoverage  <= ID->depth(ii)) &&
          (ID->depth(ii) <= maxCoverage))
        for (uint32 pp=ID->lo(ii); pp<ID->hi(ii); pp++)
          cns[pp] = '_';
    }
  };

  uint32        tigIDbgn;
  uint32        tigIDend;

  bool          dumpAllClasses;
  bool          dumpUnassembled;
  bool          dumpBubbles;
  bool          dumpContigs;

  uint32        minNreads;
  uint32        maxNreads;

  uint32        minLength;
  uint32        maxLength;

  double        minCoverage;
  double        maxCoverage;

  double        minGoodCov;
  double        maxGoodCov;

  intervalList<int32>  *IL;
  intervalList<int32>  *ID;
};



void
dumpStatus(gkStore *UNUSED(gkpStore), tgStore *tigStore) {
  fprintf(stderr, "%u\n", tigStore->numTigs());
}



void
dumpTig(FILE *out, tgTig *tig, bool useGapped) {
  fprintf(out, F_U32"\t"F_U32"\t%s\t%.2f\t%s\t%s\t%s\t"F_U32"\n",
          tig->tigID(),
          tig->length(useGapped),
          tig->coordinateType(useGapped),
          tig->_coverageStat,
          toString(tig->_class),
          tig->_suggestRepeat ? "yes" : "no",
          tig->_suggestCircular ? "yes" : "no",
          tig->numberOfChildren());
}



void
dumpRead(FILE *out, tgTig *tig, tgPosition *read, bool useGapped) {
  fprintf(out, F_U32"\t"F_U32"\t%s\t"F_U32"\t"F_U32"\n",
          read->ident(),
          tig->tigID(),
          tig->coordinateType(useGapped),
          (useGapped) ? read->bgn() : tig->mapGappedToUngapped(read->bgn()),
          (useGapped) ? read->end() : tig->mapGappedToUngapped(read->end()));
}



void
dumpTigs(gkStore *UNUSED(gkpStore), tgStore *tigStore, tgFilter &filter, bool useGapped) {

  fprintf(stdout, "tigID\ttigLen\ttype\tcovStat\tsr\tsu\tsc\tsh\tnumChildren\n");

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (tig->consensusExists() == false)
      useGapped = true;

    if (filter.ignore(tig, useGapped) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    dumpTig(stdout, tig, useGapped);

    tigStore->unloadTig(ti);
  }
}



void
dumpConsensus(gkStore *UNUSED(gkpStore), tgStore *tigStore, tgFilter &filter, bool useGapped, char cnsFormat) {

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (tig->consensusExists() == false) {
      //fprintf(stderr, "dumpConsensus()-- tig %u has no consensus sequence.\n", ti);
      tigStore->unloadTig(ti);
      continue;
    }

    if (filter.ignore(tig, useGapped) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    switch (cnsFormat) {
      case 'A':
        tig->dumpFASTA(stdout, useGapped);
        break;

      case 'Q':
        tig->dumpFASTQ(stdout, useGapped);
        break;

      default:
        break;
    }

    tigStore->unloadTig(ti);
  }
}



void
dumpLayout(gkStore *UNUSED(gkpStore), tgStore *tigStore, tgFilter &filter, bool useGapped, char *outPrefix) {

  FILE *tigs   = NULL;    //  Length and flags of tigs, same as dumpTigs()
  FILE *reads  = NULL;    //  Length and flags of reads, mapping of read to tig
  FILE *layout = stdout;  //  Standard layout file

  if (outPrefix) {
    char T[FILENAME_MAX];  int32 Terr = 0;
    char R[FILENAME_MAX];  int32 Rerr = 0;
    char L[FILENAME_MAX];  int32 Lerr = 0;

    sprintf(T, "%s.layout.tigInfo",   outPrefix);
    sprintf(R, "%s.layout.readToTig", outPrefix);
    sprintf(L, "%s.layout",           outPrefix);

    errno = 0;

    tigs   = fopen(T, "w");  Terr = errno;
    reads  = fopen(R, "w");  Rerr = errno;
    layout = fopen(L, "w");  Lerr = errno;

    if (Terr)
      fprintf(stderr, "Failed to open '%s': %s\n", T, strerror(Terr));
    if (Rerr)
      fprintf(stderr, "Failed to open '%s': %s\n", R, strerror(Rerr));
    if (Lerr)
      fprintf(stderr, "Failed to open '%s': %s\n", L, strerror(Lerr));

    if (Terr + Rerr + Lerr > 0)
      exit(1);
  }

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (tig->consensusExists() == false)
      useGapped = true;

    if (filter.ignore(tig, useGapped) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    if (tigs)
      dumpTig(tigs, tig, useGapped);

    if (reads)
      for (uint32 ci=0; ci<tig->numberOfChildren(); ci++)
        dumpRead(reads, tig, tig->getChild(ci), useGapped);

    if (layout)
      tig->dumpLayout(layout);

    tigStore->unloadTig(ti);
  }

  if (outPrefix) {
    fclose(tigs);
    fclose(reads);
    fclose(layout);
  }
}



void
dumpMultialign(gkStore *gkpStore, tgStore *tigStore, tgFilter &filter, bool maWithQV, bool maWithDots, uint32 maDisplayWidth, uint32 maDisplaySpacing) {

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (filter.ignore(tig, true) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    tig->display(stdout, gkpStore, maDisplayWidth, maDisplaySpacing, maWithQV, maWithDots);

    tigStore->unloadTig(ti);
  }
}



void
dumpSizes(gkStore *UNUSED(gkpStore), tgStore *tigStore, tgFilter &filter, bool useGapped, uint64 genomeSize) {

  tgTigSizeAnalysis *siz = new tgTigSizeAnalysis(genomeSize);

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (tig->consensusExists() == false)
      useGapped = true;

    if (filter.ignore(tig, useGapped) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    siz->evaluateTig(tig, useGapped);

    tigStore->unloadTig(ti);
  }

  siz->finalize();
  siz->printSummary(stdout);

  delete siz;
}



void
plotDepthHistogram(char *N, uint64 *cov, uint32 covMax) {
  FILE  *F;

  //  Find the smallest and largest values with counts.

  uint32  minii = 0;
  uint32  maxii = covMax - 1;

  while ((minii < covMax) && (cov[minii] == 0))
    minii++;

  while ((minii < maxii) && (cov[maxii] == 0))
    maxii--;

  //  Extend a little bit, to give a little context (and to get zero).

  minii = (minii > 10) ? minii-10 : 0;
  maxii =                maxii+10;

  //  Dump the values from min to max

  errno = 0;
  F = fopen(N, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno));

  for (uint32 ii=minii; ii<=maxii; ii++)
    for (uint32 xx=0; xx<cov[ii]; xx++)
      fprintf(F, F_U32"\n", ii);

  fclose(F);

  //  Decide on a bucket size.  We want some even number, like 10, 100, 1000, 5000.

  uint32   boxsize  = 1;
  uint32   boxscale = 5;
  uint32   nboxes   = 50;

  while (boxsize * boxscale * nboxes < maxii - minii) {
  fprintf(stderr, "boxsize %u  boxscale %u  range %u-%u = %u  nboxes %u\n",
          boxsize, boxscale, maxii, minii, maxii-minii, (maxii-minii) / boxsize);
    boxsize  *= boxscale;
    boxscale  = (boxscale == 5) ? 2 : 5;
  }

  fprintf(stderr, "boxsize %u  boxscale %u  range %u-%u = %u  nboxes %u\n",
          boxsize, boxscale, maxii, minii, maxii-minii, (maxii-minii) / boxsize);


  //  Plot!

  F = popen("gnuplot > /dev/null 2>&1", "w");

  if (F) {
    fprintf(F, "set terminal 'png'\n");
    fprintf(F, "set output '%s.png'\n", N);
    fprintf(F, "set xlabel 'binwidth=%u'\n", boxsize);
    //fprintf(F, "set ylabel 'number of bases'\n");

    //  A histogram, assuming the data is sorted values
    fprintf(F, "binwidth=%u\n", boxsize);
    fprintf(F, "set boxwidth binwidth\n");
    fprintf(F, "bin(x,width) = width*floor(x/width) + binwidth/2.0\n");
    fprintf(F, "plot [%u:%u] '%s' using (bin($1,binwidth)):(1.0) smooth freq with boxes title '%s'\n", minii, maxii, N, N);

    //  A histogram, assuming the data is a histogram.  It's ugly though.
    //fprintf(F, "plot [%u:%u] [] '%s' using 1:2 with lines title '%s', \\\n",  //  Used to report 'tigs %u-%u' for the range
    //        minii, maxii, N, N);

    pclose(F);
  }

  //  Dump the values again, this time as a real histogram.

  errno = 0;
  F = fopen(N, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno));

  for (uint32 ii=minii; ii<=maxii; ii++)
    fprintf(F, F_U32"\t"F_U64"\n", ii, cov[ii]);

  fclose(F);
}



void
dumpDepthHistogram(gkStore *UNUSED(gkpStore), tgStore *tigStore, tgFilter &filter, bool useGapped, bool single, char *outPrefix) {
  char                  N[FILENAME_MAX];
  intervalList<uint32>  IL;

  int32     covMax = 1048576;
  uint64   *cov    = new uint64 [covMax];

  memset(cov, 0, sizeof(uint64) * covMax);

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (tig->consensusExists() == false)
      useGapped = true;

    if (filter.ignore(tig, useGapped) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    //  Save all the read intervals to the list.

    IL.clear();

    for (uint32 ci=0; ci<tig->numberOfChildren(); ci++) {
      tgPosition *read = tig->getChild(ci);
      uint32      bgn  = (useGapped) ? read->min() : tig->mapGappedToUngapped(read->min());
      uint32      end  = (useGapped) ? read->max() : tig->mapGappedToUngapped(read->max());

      IL.add(bgn, end - bgn);
    }

    //  Convert to depths.

    intervalList<uint32>  ID(IL);

    //  Add the depths to the histogram.

    for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++)
      cov[ID.depth(ii)] += ID.hi(ii) - ID.lo(ii);

    //  Maybe plot the histogram (and if so, clear it for the next tig).

    if (single == true) {
      sprintf(N, "%s.tig%06d.depthHistogram", outPrefix, tig->tigID());
      plotDepthHistogram(N, cov, covMax);

      memset(cov, 0, sizeof(uint64) * covMax);  //  Slight optimization if we do this in plotDepthHistogram of just the set values.
    }

    //  Repeat.

    tigStore->unloadTig(ti);
  }

  if (single == false) {
    sprintf(N, "%s.depthHistogram", outPrefix);
    plotDepthHistogram(N, cov, covMax);
  }

  delete cov;
}



void
dumpCoverage(gkStore *UNUSED(gkpStore), tgStore *tigStore, tgFilter &filter, bool useGapped, char *outPrefix) {
  uint32   covMax = 1024;
  uint64  *cov    = new uint64 [covMax];

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    tgTig    *tig    = tigStore->loadTig(ti);
    uint32    tigLen = tig->length(useGapped);

    if (tig->consensusExists() == false)
      useGapped = true;

    if (filter.ignore(tig, true) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    //  Do something.

    intervalList<int32>  allL;

    for (uint32 ci=0; ci<tig->numberOfChildren(); ci++) {
      tgPosition *read = tig->getChild(ci);
      uint32      bgn  = (useGapped) ? read->min() : tig->mapGappedToUngapped(read->min());
      uint32      end  = (useGapped) ? read->max() : tig->mapGappedToUngapped(read->max());

      allL.add(bgn, end - bgn);
    }

    intervalList<int32>   ID(allL);

    uint32  maxDepth    = 0;
    double  aveDepth    = 0;
    double  sdeDepth    = 0;

#if 0
    //  Report regions that have abnormally low or abnormally high coverage

    intervalList<int32>   minL;
    intervalList<int32>   maxL;

    for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
      if ((ID.depth(ii) < minCoverage) && (ID.lo(ii) != 0) && (ID.hi(ii) != tigLen)) {
        fprintf(stderr, "tig %d low coverage interval %ld %ld max %u coverage %u\n",
                tig->tigID(), ID.lo(ii), ID.hi(ii), tigLen, ID.depth(ii));
        minL.add(ID.lo(ii), ID.hi(ii) - ID.lo(ii) + 1);
      }

      if (maxCoverage <= ID.depth(ii)) {
        fprintf(stderr, "tig %d high coverage interval %ld %ld max %u coverage %u\n",
                tig->tigID(), ID.lo(ii), ID.hi(ii), tigLen, ID.depth(ii));
        maxL.add(ID.lo(ii), ID.hi(ii) - ID.lo(ii) + 1);
      }
    }
#endif

    //  Compute max and average depth, and save the depth in a histogram.
#warning replace this with genericStatistics

    for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
      if (ID.depth(ii) > maxDepth)
        maxDepth = ID.depth(ii);

      aveDepth += (ID.hi(ii) - ID.lo(ii) + 1) * ID.depth(ii);

      while (covMax <= ID.depth(ii))
        resizeArray(cov, covMax, covMax, covMax * 2);

      cov[ID.depth(ii)] += ID.hi(ii) - ID.lo(ii) + 1;
    }

    aveDepth /= tigLen;

    //  Now the std.dev

    for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++)
      sdeDepth += (ID.hi(ii) - ID.lo(ii) + 1) * (ID.depth(ii) - aveDepth) * (ID.depth(ii) - aveDepth);

    sdeDepth = sqrt(sdeDepth / tigLen);

    //  Merge the intervals to figure out what has coverage, or what is missing coverage.

#if 0
    allL.merge();
    minL.merge();
    maxL.merge();

    if      ((minL.numberOfIntervals() > 0) && (maxL.numberOfIntervals() > 0))
      fprintf(stderr, "tig %d has %u intervals, %u regions below %u coverage and %u regions at or above %u coverage\n",
              tig->tigID(),
              allL.numberOfIntervals(),
              minL.numberOfIntervals(), minCoverage,
              maxL.numberOfIntervals(), maxCoverage);
    else if (minL.numberOfIntervals() > 0)
      fprintf(stderr, "tig %d has %u intervals, %u regions below %u coverage\n",
              tig->tigID(),
              allL.numberOfIntervals(),
              minL.numberOfIntervals(), minCoverage);
    else if (maxL.numberOfIntervals() > 0)
      fprintf(stderr, "tig %d has %u intervals, %u regions at or above %u coverage\n",
              tig->tigID(),
              allL.numberOfIntervals(),
              maxL.numberOfIntervals(), maxCoverage);
    else
      fprintf(stderr, "tig %d has %u intervals\n",
              tig->tigID(),
              allL.numberOfIntervals());
#endif

    //  Plot the depth for each tig

    if (outPrefix) {
      char  outName[FILENAME_MAX];

      sprintf(outName, "%s.tig%08u.depth", outPrefix, tig->tigID());

      FILE *outFile = fopen(outName, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s': %s\n", outName, strerror(errno)), exit(1);

      for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
        fprintf(outFile, "%d\t%u\n", ID.lo(ii),     ID.depth(ii));
        fprintf(outFile, "%d\t%u\n", ID.hi(ii) - 1, ID.depth(ii));
      }

      fclose(outFile);

      FILE *gnuPlot = popen("gnuplot > /dev/null 2>&1", "w");

      if (gnuPlot) {
        fprintf(gnuPlot, "set terminal 'png'\n");
        fprintf(gnuPlot, "set output '%s.tig%08u.png'\n", outPrefix, tig->tigID());
        fprintf(gnuPlot, "set xlabel 'position'\n");
        fprintf(gnuPlot, "set ylabel 'coverage'\n");
        fprintf(gnuPlot, "set terminal 'png'\n");
        fprintf(gnuPlot, "plot '%s.tig%08u.depth' using 1:2 with lines title 'tig %u length %u', \\\n",
                outPrefix,
                tig->tigID(),
                tig->tigID(), tigLen);
        fprintf(gnuPlot, "     %f title 'mean %.2f +- %.2f', \\\n", aveDepth, aveDepth, sdeDepth);
        fprintf(gnuPlot, "     %f title '' lt 0 lc 2, \\\n", aveDepth - sdeDepth);
        fprintf(gnuPlot, "     %f title '' lt 0 lc 2\n",     aveDepth + sdeDepth);

        fclose(gnuPlot);
      }
    }

    //  Did something.

    tigStore->unloadTig(ti);
  }
}



void
dumpThinOverlap(gkStore *UNUSED(gkpStore), tgStore *tigStore, tgFilter &filter, bool useGapped, uint32 minOverlap) {

  fprintf(stderr, "reporting overlaps of at most %u bases\n", minOverlap);

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (tig->consensusExists() == false)
      useGapped = true;

    if (filter.ignore(tig, true) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    //  Do something.

    intervalList<int32>  allL;
    intervalList<int32>  ovlL;
    intervalList<int32>  badL;

    for (uint32 ri=0; ri<tig->numberOfChildren(); ri++) {
      tgPosition *read = tig->getChild(ri);
      uint32      bgn  = (useGapped) ? read->min() : tig->mapGappedToUngapped(read->min());
      uint32      end  = (useGapped) ? read->max() : tig->mapGappedToUngapped(read->max());

      allL.add(bgn, end - bgn);
      ovlL.add(bgn, end - bgn);
    }

    allL.merge();            //  Merge, requiring zero overlap (adjacent is OK) between pieces
    ovlL.merge(minOverlap);  //  Merge, requiring minOverlap overlap between pieces

    //  If there is more than one interval, make a list of the regions where we have thin overlaps.

    if (ovlL.numberOfIntervals() > 1)  //  Vertical space between tig reports
      fprintf(stderr, "\n");

    for (uint32 ii=1; ii<ovlL.numberOfIntervals(); ii++) {
      assert(ovlL.lo(ii) < ovlL.hi(ii-1));

      fprintf(stderr, "tig %d thin %u %u\n", tig->tigID(), ovlL.lo(ii), ovlL.hi(ii-1));

      badL.add(ovlL.lo(ii), ovlL.hi(ii-1) - ovlL.lo(ii));
    }

    //  Then report any reads that intersect that region.

    for (uint32 ri=0; ri<tig->numberOfChildren(); ri++) {
      tgPosition *read   = tig->getChild(ri);
      uint32      bgn    = (useGapped) ? read->min() : tig->mapGappedToUngapped(read->min());
      uint32      end    = (useGapped) ? read->max() : tig->mapGappedToUngapped(read->max());
      bool        report = false;

      for (uint32 oo=0; oo<badL.numberOfIntervals(); oo++)
        if ((badL.lo(oo) <= end) &&
            (bgn         <= badL.hi(oo))) {
          report = true;
          break;
        }

      if (report)
        fprintf(stderr, "tig %d read %u at %u %u\n",
                tig->tigID(),
                read->ident(),
                (useGapped) ? read->min() : tig->mapGappedToUngapped(read->min()),
                (useGapped) ? read->max() : tig->mapGappedToUngapped(read->max()));
    }

    if ((allL.numberOfIntervals() != 1) || (ovlL.numberOfIntervals() != 1))
      fprintf(stderr, "tig %d %s length %u has %u interval%s and %u interval%s after enforcing minimum overlap of %u\n",
              tig->tigID(), tig->coordinateType(useGapped), tig->length(),
              allL.numberOfIntervals(), (allL.numberOfIntervals() == 1) ? "" : "s",
              ovlL.numberOfIntervals(), (ovlL.numberOfIntervals() == 1) ? "" : "s",
              minOverlap);

    //  There, did something.

    tigStore->unloadTig(ti);
  }
}



void
dumpOverlapHistogram(gkStore *UNUSED(gkpStore), tgStore *tigStore, tgFilter &filter, bool useGapped, char *outPrefix) {
  uint32     histMax = AS_MAX_READLEN;
  uint64    *hist    = new uint64 [histMax];

  memset(hist, 0, sizeof(uint64) * histMax);

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    tgTig  *tig = tigStore->loadTig(ti);
    int32   tn  = tig->numberOfChildren();

    if (tig->consensusExists() == false)
      useGapped = true;

    if (filter.ignore(tig, true) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    //  Do something.  For each read, compute the thickest overlap off of each end.

    //  First, decide on positions for each read.  Store in an array for easier use later.

    uint32   *bgn = new uint32 [tn];
    uint32   *end = new uint32 [tn];

    for (uint32 ri=0; ri<tn; ri++) {
      tgPosition *read = tig->getChild(ri);

      bgn[ri] = (useGapped) ? read->min() : tig->mapGappedToUngapped(read->min());
      end[ri] = (useGapped) ? read->max() : tig->mapGappedToUngapped(read->max());
    }

    //  Scan these, marking contained reads.

    for (uint32 ri=0; ri<tn; ri++)
      for (uint32 ii=ri+1; ii<tn && bgn[ii] < end[ti]; ii++)
        if ((bgn[ri] <= bgn[ii]) && (end[ii] <= end[ti])) {
          bgn[ii] = UINT32_MAX;
          end[ii] = UINT32_MAX;
          break;
        }

    //  Now, scan the overlaps finding thickest.  There are no contained reads, and so we're guaranteed
    //  that as soon as we stop seeing overlaps, we'll see no more overlaps.

    for (uint32 ri=0; ri<tn; ri++) {
      uint32  thickest5 = 0;
      uint32  thickest3 = 0;

      if (bgn[ri] == UINT32_MAX)  //  Read is contained, no useful overlaps to report.
        continue;

      //  Off the 5' end, expect end[ii] < end[ri] and end[ii] > bgn[ri]
      for (int32 ii=ri-1; ii>0; ii--) {
        if (bgn[ii] == UINT32_MAX)
          continue;

        if (end[ii] < bgn[ri])  //  Read doesn't overlap, no more reads will.
          break;

        if (thickest5 < end[ii] - bgn[ri])
          thickest5 = end[ii] - bgn[ri];
      }

      //  Off the 3' end, expect bgn[ii] < end[ri] and bgn[ii] > bgn[ri]
      for (int32 ii=ri+1; ii<tn; ii++) {
        if (bgn[ii] == UINT32_MAX)
          continue;

        if (end[ri] < bgn[ii])  //  Read doesn't overlap, no more reads will.
          break;

        if (thickest5 < end[ri] - bgn[ii])
          thickest5 = end[ri] - bgn[ii];
      }

      //  Save those thickest (but not the boring zero cases).  Contained reads end up with no thickest overlaps.

      if (thickest5 > 0) {
        assert(thickest5 < histMax);
        hist[thickest5]++;
      }

      if (thickest3 > 0) {
        assert(thickest3 < histMax);
        hist[thickest3]++;
      }
    }

    delete [] bgn;
    delete [] end;

    //  There, did something.

    tigStore->unloadTig(ti);
  }

  //  All computed.  Dump the data and plot.

  char N[FILENAME_MAX];

  sprintf(N, "%s.thickestOverlapHistogram", outPrefix);

  plotDepthHistogram(N, hist, histMax);

  //  Cleanup and Bye!

  delete [] hist;
}





int
main (int argc, char **argv) {
  char         *gkpName           = NULL;
  char         *tigName           = NULL;
  int           tigVers           = -1;

  //  Tig Selection

  tgFilter      filter;

  //  Dump options

  uint32        dumpType          = DUMP_UNSET;

  bool          useGapped         = false;

  char          cnsFormat         = 'A';  //  Or 'Q' for FASTQ

  bool          maWithQV          = false;
  bool          maWithDots        = true;
  uint32        maDisplayWidth    = 100;
  uint32        maDisplaySpacing  = 3;

  uint64        genomeSize        = 0;

  char         *outPrefix         = NULL;

  bool          single            = false;

  uint32        minOverlap        = 0;


  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-tig") == 0) {
      AS_UTL_decodeRange(argv[++arg], filter.tigIDbgn, filter.tigIDend);
    }

    else if (strcmp(argv[arg], "-unassembled") == 0) {
      filter.dumpAllClasses  = false;
      filter.dumpUnassembled = true;
    }

    else if (strcmp(argv[arg], "-bubbles") == 0) {
      filter.dumpAllClasses  = false;
      filter.dumpBubbles     = true;
    }

    else if (strcmp(argv[arg], "-contigs") == 0) {
      filter.dumpAllClasses  = false;
      filter.dumpContigs     = true;
    }


    else if (strcmp(argv[arg], "-nreads") == 0) {
      filter.minNreads = atoi(argv[++arg]);
      filter.maxNreads = atoi(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-length") == 0) {
      filter.minLength = atoi(argv[++arg]);
      filter.maxLength = atoi(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-coverage") == 0) {
      if ((arg == argc-1) || (argv[arg+1][0] == '-')) {
        dumpType = DUMP_COVERAGE;
      }

      else if (arg + 4 < argc) {
        filter.minCoverage = atof(argv[++arg]);
        filter.maxCoverage = atof(argv[++arg]);
        filter.minGoodCov  = atof(argv[++arg]);
        filter.maxGoodCov  = atof(argv[++arg]);
      }

      else {
        fprintf(stderr, "ERROR: -coverage needs four values.\n");
        err++;
      }
    }

    //  Dump types.

    else if (strcmp(argv[arg], "-status") == 0)
      dumpType = DUMP_STATUS;
    else if (strcmp(argv[arg], "-tigs") == 0)
      dumpType = DUMP_TIGS;
    else if (strcmp(argv[arg], "-consensus") == 0)
      dumpType = DUMP_CONSENSUS;
    else if (strcmp(argv[arg], "-layout") == 0)
      dumpType = DUMP_LAYOUT;
    else if (strcmp(argv[arg], "-multialign") == 0)
      dumpType = DUMP_MULTIALIGN;
    else if (strcmp(argv[arg], "-sizes") == 0)
      dumpType = DUMP_SIZES;
    else if (strcmp(argv[arg], "-coverage") == 0)  //  NOTE!  Actually handled above.
      dumpType = DUMP_COVERAGE;
    else if (strcmp(argv[arg], "-depth") == 0)
      dumpType = DUMP_DEPTH_HISTOGRAM;
    else if (strcmp(argv[arg], "-overlap") == 0)
      dumpType = DUMP_THIN_OVERLAP;
    else if (strcmp(argv[arg], "-overlaphistogram") == 0)
      dumpType = DUMP_OVERLAP_HISTOGRAM;

    //  Options.

    else if (strcmp(argv[arg], "-gapped") == 0)
      useGapped = true;

    else if (strcmp(argv[arg], "-fasta") == 0)
      cnsFormat = 'A';
    else if (strcmp(argv[arg], "-fastq") == 0)
      cnsFormat = 'Q';

    else if (strcmp(argv[arg], "-w") == 0)
      maDisplayWidth = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-s") == 0)
      maDisplaySpacing = genomeSize = atol(argv[++arg]);

    else if (strcmp(argv[arg], "-o") == 0)
      outPrefix = argv[++arg];

    else if (strcmp(argv[arg], "-single") == 0)
      single = true;

    else if (strcmp(argv[arg], "-thin") == 0)
      minOverlap = atoi(argv[++arg]);

    //  Errors.

    else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }

  if (gkpName == NULL)
    err++;
  if (tigName == NULL)
    err++;
  if ((outPrefix == NULL) && (dumpType == DUMP_COVERAGE))
    err++;
  if (dumpType == DUMP_UNSET)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -G <gkpStore> -T <tigStore> <v> [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "STORE SELECTION (mandatory)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G <gkpStore>           path to the gatekeeper store\n");
    fprintf(stderr, "  -T <tigStore> <v>       path to the tigStore, version, to use\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "TIG SELECTION - if nothing specified, all tigs are reported\n");
    fprintf(stderr, "              - all ranges are inclusive.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -tig A[-B]              only dump tigs between ids A and B\n");
    fprintf(stderr, "  -unassembled            only dump tigs that are 'unassembled'\n");
    fprintf(stderr, "  -bubbles                only dump tigs that are 'bubbles'\n");
    fprintf(stderr, "  -contigs                only dump tigs that are 'contigs'\n");
    fprintf(stderr, "  -nreads min max         only dump tigs with between min and max reads\n");
    fprintf(stderr, "  -length min max         only dump tigs with length between 'min' and 'max' bases\n");
    fprintf(stderr, "  -coverage c C g G       only dump tigs with between fraction g and G at coverage between c and C\n");
    fprintf(stderr, "                            example:  -coverage 10 inf 0.5 1.0 would report tigs where half of the\n");
    fprintf(stderr, "                                      bases are at 10+ times coverage.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DUMP TYPE - all dumps, except status, report on tigs selected as above\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -status                 the number of tigs in the store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -tigs                   a list of tigs, and some information about them\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -consensus [opts]       the consensus sequence, with options:\n");
    fprintf(stderr, "                            -gapped           report the gapped (multialignment) consensus sequence\n");
    fprintf(stderr, "                            -fasta            report sequences in FASTA format (the default)\n");
    fprintf(stderr, "                            -fastq            report sequences in FASTQ format\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -layout [opts]          the layout of reads in each tig\n");
    fprintf(stderr, "                          if '-o' is supplied, three files are created, otherwise just the layout is printed to stdout\n");
    fprintf(stderr, "                            -gapped           report the gapped (multialignment) positions\n");
    fprintf(stderr, "                            -o outputPrefix   write plots to 'outputPrefix.*' in the current directory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -multialign [opts]      the full multialignment, output is to stdout\n");
    fprintf(stderr, "                            -w width          width of the page\n");
    fprintf(stderr, "                            -s spacing        spacing between reads on the same line\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -sizes [opts]           size statistics\n");
    fprintf(stderr, "                            -s genomesize     denominator to use for n50 computation\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -coverage [opts]        read coverage plots, one plot per tig\n");
    fprintf(stderr, "                            -o outputPrefix   write plots to 'outputPrefix.*' in the current directory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -depth [opts]           a histogram of depths\n");
    fprintf(stderr, "                            -single           one histogram per tig\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -overlap                read overlaps\n");
    fprintf(stderr, "                            -thin overlap     report regions where the (thickest) read overlap is less than 'overlap' bases\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -overlaphistogram       a histogram of the thickest overlaps used\n");
    fprintf(stderr, "                            -o outputPrefix   write plots to 'outputPrefix.*' in the current directory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

#if 0
    fprintf(stderr, "  -compress             Move tigs from earlier versions into the specified version.  This removes\n");
    fprintf(stderr, "                        historical versions of unitigs/contigs, and can save tremendous storage space,\n");
    fprintf(stderr, "                        but makes it impossible to back up the assembly past the specified versions\n");
#endif

    if (gkpName == NULL)
      err++;
    if (tigName == NULL)
      err++;
    if ((outPrefix == NULL) && (dumpType == DUMP_COVERAGE))
      err++;
    if (dumpType == DUMP_UNSET)
      err++;

    exit(1);
  }

  //  Open stores.

  gkStore *gkpStore = gkStore::gkStore_open(gkpName);
  tgStore *tigStore = new tgStore(tigName, tigVers);

  //  Check that the tig ID range is valid, and fix it if possible.

  uint32   nTigs = tigStore->numTigs();

  if (filter.tigIDend == UINT32_MAX)
    filter.tigIDend = nTigs-1;

  if (nTigs <= filter.tigIDend) {
    fprintf(stderr, "WARNING: adjusting tig ID range from "F_U32"-"F_U32" to "F_U32"-"F_U32" as there are only "F_U32" tigs in the store.\n",
            filter.tigIDbgn, filter.tigIDend, filter.tigIDbgn, nTigs-1, nTigs);
    filter.tigIDend = nTigs - 1;
  }

  if (filter.tigIDend < filter.tigIDbgn) {
    fprintf(stderr, "WARNING: adjusting inverted tig ID range -t "F_U32"-"F_U32"\n",
            filter.tigIDbgn, filter.tigIDend);
    uint32 x = filter.tigIDend;
    filter.tigIDend = filter.tigIDbgn;
    filter.tigIDbgn = x;
  }

  if (nTigs <= filter.tigIDbgn)
    fprintf(stderr, "ERROR: only "F_U32" tigs in the store (IDs 0-"F_U32" inclusive); can't dump requested range -t "F_U32"-"F_U32"\n",
            nTigs,
            nTigs-1,
            filter.tigIDbgn, filter.tigIDend), exit(1);

  //  Call the dump routine.

  switch (dumpType) {
    case DUMP_STATUS:
      dumpStatus(gkpStore, tigStore);
      break;
    case DUMP_TIGS:
      dumpTigs(gkpStore, tigStore, filter, useGapped);
      break;
    case DUMP_CONSENSUS:
      dumpConsensus(gkpStore, tigStore, filter, useGapped, cnsFormat);
      break;
    case DUMP_LAYOUT:
      dumpLayout(gkpStore, tigStore, filter, useGapped, outPrefix);
      break;
    case DUMP_MULTIALIGN:
      dumpMultialign(gkpStore, tigStore, filter, maWithQV, maWithDots, maDisplayWidth, maDisplaySpacing);
      break;
    case DUMP_SIZES:
      dumpSizes(gkpStore, tigStore, filter, useGapped, genomeSize);
      break;
    case DUMP_COVERAGE:
      dumpCoverage(gkpStore, tigStore, filter, useGapped, outPrefix);
      break;
    case DUMP_DEPTH_HISTOGRAM:
      dumpDepthHistogram(gkpStore, tigStore, filter, useGapped, single, outPrefix);
      break;
    case DUMP_THIN_OVERLAP:
      dumpThinOverlap(gkpStore, tigStore, filter, useGapped, minOverlap);
      break;
    case DUMP_OVERLAP_HISTOGRAM:
      dumpOverlapHistogram(gkpStore, tigStore, filter, useGapped, outPrefix);
      break;
    default:
      break;
  }

  //  Clean up.

  delete tigStore;

  gkpStore->gkStore_close();

  //  Bye.

  exit(0);
}

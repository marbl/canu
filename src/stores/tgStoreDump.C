
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

#include "strings.H"

#include "sqStore.H"
#include "tgStore.H"

#include "intervalList.H"

#include "tgTigSizeAnalysis.H"

#undef  DEBUG_IGNORE

#define DUMP_UNSET               0
#define DUMP_STATUS              1
#define DUMP_TIGS                2
#define DUMP_CONSENSUS           3
#define DUMP_LAYOUT              4
#define DUMP_INFO                5
#define DUMP_MULTIALIGN          6
#define DUMP_SIZES               7
#define DUMP_COVERAGE            8
#define DUMP_DEPTH_HISTOGRAM     9
#define DUMP_THIN_OVERLAP       10
#define DUMP_OVERLAP_HISTOGRAM  11


class tgFilter {
public:
  tgFilter() {
    tigIDbgn        = 0;
    tigIDend        = UINT32_MAX;

    dumpAllClasses  = true;
    dumpUnassembled = false;
    dumpContigs     = false;

    dumpRepeats     = false;
    dumpBubbles     = false;
    dumpCircular    = false;

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

  bool          ignore(uint32 id) {
    return(ignoreID(id));
  };

  bool          ignore(tgTig *tig) {
    return(ignoreID(tig->tigID()) ||
           ignoreNreads(tig)      ||
           ignoreLength(tig)      ||
           ignoreCoverage(tig)    ||
           ignoreClass(tig));
  };

  bool          ignoreID(uint32 id) {
    return((id < tigIDbgn) ||
           (tigIDend < id));
  };

  bool          ignoreClass(tgTig *tig) {
    if ((dumpAllClasses == true) ||
        ((tig->_class == tgTig_unassembled) && (dumpUnassembled == true)) ||
        ((tig->_class == tgTig_contig)      && (dumpContigs == true)))
      return(false);

    return(true);
  };

  bool          ignoreNreads(tgTig *tig) {
    return((tig->numberOfChildren() < minNreads) ||
           (maxNreads < tig->numberOfChildren()));
  };

  bool          ignoreLength(tgTig *tig) {
    uint32  length = tig->length();

    return((length < minLength) ||
           (maxLength < length));
  };

  bool          ignoreCoverage(tgTig *tig) {
    if ((minCoverage == 0) && (maxCoverage == UINT32_MAX))
      return(false);

    delete IL;
    IL = new intervalList<int32>;

    for (uint32 i=0; i<tig->numberOfChildren(); i++) {
      tgPosition *pos = tig->getChild(i);

      int32  bgn = pos->min();
      int32  end = pos->max();

      IL->add(bgn, end - bgn);
    }

    delete ID;
    ID = new intervalDepth<int32>(*IL);

    uint32  goodCov  = 0;
    uint32  badCov   = 0;
    double  fracGood = 0.0;

    for (uint32 ii=0; ii<ID->numberOfIntervals(); ii++)
      if ((minCoverage  <= ID->depth(ii)) &&
          (ID->depth(ii) <= maxCoverage))
        goodCov += ID->hi(ii) - ID->lo(ii);
      else
        badCov += ID->hi(ii) - ID->lo(ii);

    if (goodCov + badCov > 0)
      fracGood = (double)(goodCov) / (goodCov + badCov);

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
  bool          dumpContigs;

  bool          dumpRepeats;
  bool          dumpBubbles;
  bool          dumpCircular;

  uint32        minNreads;
  uint32        maxNreads;

  uint32        minLength;
  uint32        maxLength;

  double        minCoverage;
  double        maxCoverage;

  double        minGoodCov;
  double        maxGoodCov;

  intervalList<int32>  *IL;
  intervalDepth<int32> *ID;
};



void
dumpStatus(sqStore *UNUSED(seqStore), tgStore *tigStore) {
  fprintf(stderr, "%u\n", tigStore->numTigs());
}



void
dumpTigHeader(FILE *out) {
  fprintf(out, "#tigID\ttigLen\tcoverage\ttigClass\tsugRept\tsugBubb\tsugCirc\tnumChildren\n");
}

void
dumpTig(FILE *out, tgTig *tig) {
  fprintf(out, F_U32"\t" F_U32 "\t%.2f\t%s\t%s\t%s\t%s\t" F_U32 "\n",
          tig->tigID(),
          tig->length(),
          tig->computeCoverage(),
          toString(tig->_class),
          tig->_suggestRepeat ? "yes" : "no",
          tig->_suggestBubble ? "yes" : "no",
          tig->_suggestCircular ? "yes" : "no",
          tig->numberOfChildren());
}



void
dumpReadHeader(FILE *out) {
  fprintf(out, "#readID\ttigID\tbgn\tend\n");
}

void
dumpRead(FILE *out, tgTig *tig, tgPosition *read) {
  fprintf(out, F_U32"\t" F_U32 "\t" F_U32 "\t" F_U32 "\n",
          read->ident(),
          tig->tigID(),
          read->bgn(),
          read->end());
}



void
dumpTigs(sqStore *UNUSED(seqStore), tgStore *tigStore, tgFilter &filter) {

  dumpTigHeader(stdout);

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (filter.ignore(ti) == true)
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (filter.ignore(tig) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    dumpTig(stdout, tig);

    tigStore->unloadTig(ti);
  }
}



void
dumpConsensus(sqStore *UNUSED(seqStore), tgStore *tigStore, tgFilter &filter, bool useReverse, char cnsFormat) {

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (filter.ignore(ti) == true)
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (tig->consensusExists() == false) {
      //fprintf(stderr, "dumpConsensus()-- tig %u has no consensus sequence.\n", ti);
      tigStore->unloadTig(ti);
      continue;
    }

    if (filter.ignore(tig) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    if (useReverse)
      tig->reverseComplement();

    switch (cnsFormat) {
      case 'A':
        tig->dumpFASTA(stdout);
        break;

      case 'Q':
        tig->dumpFASTQ(stdout);
        break;

      default:
        break;
    }

    tigStore->unloadTig(ti);
  }
}



void
dumpInfo(sqStore *UNUSED(seqStore), tgStore *tigStore, tgFilter &filter, char *outPrefix) {
  FILE *tigs   = AS_UTL_openOutputFile(outPrefix, '.', "layout.tigInfo");
  FILE *reads  = AS_UTL_openOutputFile(outPrefix, '.', "layout.readToTig");

  dumpTigHeader(tigs);
  dumpReadHeader(reads);

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (filter.ignore(ti) == true)
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (filter.ignore(tig) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    if (tigs)
      dumpTig(tigs, tig);

    if (reads)
      for (uint32 ci=0; ci<tig->numberOfChildren(); ci++)
        dumpRead(reads, tig, tig->getChild(ci));

    tigStore->unloadTig(ti);
  }

  AS_UTL_closeFile(tigs,  outPrefix, '.', "layout.tigInfo");
  AS_UTL_closeFile(reads, outPrefix, '.', "layout.readToTig");
}



void
dumpLayout(sqStore *UNUSED(seqStore), tgStore *tigStore, tgFilter &filter, bool withSequence) {

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (filter.ignore(ti) == true)
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (filter.ignore(tig) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    tig->dumpLayout(stdout, withSequence);

    tigStore->unloadTig(ti);
  }
}



void
dumpMultialign(sqStore *seqStore, tgStore *tigStore, tgFilter &filter, bool maWithQV, bool maWithDots, uint32 maDisplayWidth, uint32 maDisplaySpacing) {

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (filter.ignore(ti) == true)
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (filter.ignore(tig) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    tig->display(stdout, seqStore, maDisplayWidth, maDisplaySpacing, maWithQV, maWithDots);

    tigStore->unloadTig(ti);
  }
}



void
dumpSizes(sqStore *UNUSED(seqStore), tgStore *tigStore, tgFilter &filter, uint64 genomeSize) {

  tgTigSizeAnalysis *siz = new tgTigSizeAnalysis(genomeSize);

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (filter.ignore(ti) == true)
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (filter.ignore(tig) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    siz->evaluateTig(tig);

    tigStore->unloadTig(ti);
  }

  siz->finalize();
  siz->printSummary(stdout);

  delete siz;
}



void
plotDepthHistogram(char *N, uint64 *cov, uint32 covMax) {

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

  FILE *F = AS_UTL_openOutputFile(N);

  for (uint32 ii=minii; ii<=maxii; ii++)
    for (uint32 xx=0; xx<cov[ii]; xx++)
      fprintf(F, F_U32"\n", ii);

  AS_UTL_closeFile(F, N);

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

  F = AS_UTL_openOutputFile(N);

  for (uint32 ii=minii; ii<=maxii; ii++)
    fprintf(F, F_U32"\t" F_U64 "\n", ii, cov[ii]);

  AS_UTL_closeFile(F, N);
}



void
dumpDepthHistogram(sqStore *UNUSED(seqStore), tgStore *tigStore, tgFilter &filter, bool single, char *outPrefix) {
  char                  N[FILENAME_MAX];
  intervalList<uint32>  IL;

  int32     covMax = 1048576;
  uint64   *cov    = new uint64 [covMax];

  memset(cov, 0, sizeof(uint64) * covMax);

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (filter.ignore(ti) == true)
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (filter.ignore(tig) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    //  Save all the read intervals to the list.

    IL.clear();

    for (uint32 ci=0; ci<tig->numberOfChildren(); ci++) {
      tgPosition *read = tig->getChild(ci);
      uint32      bgn  = read->min();
      uint32      end  = read->max();

      IL.add(bgn, end - bgn);
    }

    //  Convert to depths.

    intervalDepth<uint32> ID(IL);

    //  Add the depths to the histogram.

    for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++)
      cov[ID.depth(ii)] += ID.hi(ii) - ID.lo(ii);

    //  Maybe plot the histogram (and if so, clear it for the next tig).

    if (single == true) {
      snprintf(N, FILENAME_MAX, "%s.tig%06d.depthHistogram", outPrefix, tig->tigID());
      plotDepthHistogram(N, cov, covMax);

      memset(cov, 0, sizeof(uint64) * covMax);  //  Slight optimization if we do this in plotDepthHistogram of just the set values.
    }

    //  Repeat.

    tigStore->unloadTig(ti);
  }

  if (single == false) {
    snprintf(N, FILENAME_MAX, "%s.depthHistogram", outPrefix);
    plotDepthHistogram(N, cov, covMax);
  }

  delete [] cov;
}



void
dumpCoverage(sqStore *UNUSED(seqStore), tgStore *tigStore, tgFilter &filter, char *outPrefix) {
  uint32   covMax = 1024;
  uint64  *cov    = new uint64 [covMax];

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (filter.ignore(ti) == true)
      continue;

    tgTig    *tig    = tigStore->loadTig(ti);
    uint32    tigLen = tig->length();

    if (filter.ignore(tig) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    if (tigLen == 0) {
      tigStore->unloadTig(ti);
      continue;
    }

    //  Do something.

    intervalList<int32>  allL;

    for (uint32 ci=0; ci<tig->numberOfChildren(); ci++) {
      tgPosition *read = tig->getChild(ci);
      uint32      bgn  = read->min();
      uint32      end  = read->max();

      allL.add(bgn, end - bgn);
    }

    intervalDepth<int32>  ID(allL);

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

      snprintf(outName, FILENAME_MAX, "%s.tig%08u.depth", outPrefix, tig->tigID());

      FILE *outFile = AS_UTL_openOutputFile(outName);

      for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
        fprintf(outFile, "%d\t%u\n", ID.lo(ii),     ID.depth(ii));
        fprintf(outFile, "%d\t%u\n", ID.hi(ii) - 1, ID.depth(ii));
      }

      AS_UTL_closeFile(outFile, outName);

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

        pclose(gnuPlot);
      }
    }

    //  Did something.

    tigStore->unloadTig(ti);
  }

  delete [] cov;
}



void
dumpThinOverlap(sqStore *UNUSED(seqStore), tgStore *tigStore, tgFilter &filter, uint32 minOverlap) {

  fprintf(stderr, "reporting overlaps of at most %u bases\n", minOverlap);

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (filter.ignore(ti) == true)
      continue;

    tgTig  *tig = tigStore->loadTig(ti);

    if (filter.ignore(tig) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    //  Do something.

    intervalList<int32>  allL;
    intervalList<int32>  ovlL;
    intervalList<int32>  badL;

    for (uint32 ri=0; ri<tig->numberOfChildren(); ri++) {
      tgPosition *read = tig->getChild(ri);
      uint32      bgn  = read->min();
      uint32      end  = read->max();

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
      uint32      bgn    = read->min();
      uint32      end    = read->max();
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
                read->min(),
                read->max());
    }

    if ((allL.numberOfIntervals() != 1) || (ovlL.numberOfIntervals() != 1))
      fprintf(stderr, "tig %d length %u has %u interval%s and %u interval%s after enforcing minimum overlap of %u\n",
              tig->tigID(), tig->length(),
              allL.numberOfIntervals(), (allL.numberOfIntervals() == 1) ? "" : "s",
              ovlL.numberOfIntervals(), (ovlL.numberOfIntervals() == 1) ? "" : "s",
              minOverlap);

    //  There, did something.

    tigStore->unloadTig(ti);
  }
}



void
dumpOverlapHistogram(sqStore *UNUSED(seqStore), tgStore *tigStore, tgFilter &filter, char *outPrefix) {
  uint32     histMax = AS_MAX_READLEN;
  uint64    *hist    = new uint64 [histMax];

  memset(hist, 0, sizeof(uint64) * histMax);

  for (uint32 ti=0; ti<tigStore->numTigs(); ti++) {
    if (tigStore->isDeleted(ti))
      continue;

    if (filter.ignore(ti) == true)
      continue;

    tgTig  *tig = tigStore->loadTig(ti);
    int32   tn  = tig->numberOfChildren();

    if (filter.ignore(tig) == true) {
      tigStore->unloadTig(ti);
      continue;
    }

    //  Do something.  For each read, compute the thickest overlap off of each end.

    //  First, decide on positions for each read.  Store in an array for easier use later.

    uint32   *bgn = new uint32 [tn];
    uint32   *end = new uint32 [tn];

    for (uint32 ri=0; ri<tn; ri++) {
      tgPosition *read = tig->getChild(ri);

      bgn[ri] = read->min();
      end[ri] = read->max();
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

  snprintf(N, FILENAME_MAX, "%s.thickestOverlapHistogram", outPrefix);

  plotDepthHistogram(N, hist, histMax);

  //  Cleanup and Bye!

  delete [] hist;
}





int
main (int argc, char **argv) {
  char const   *seqName           = NULL;
  char const   *tigName           = NULL;
  int           tigVers           = -1;

  //  Tig Selection

  tgFilter      filter;

  //  Dump options

  uint32        dumpType          = DUMP_UNSET;

  bool          useReverse        = false;

  char          cnsFormat         = 'A';  //  Or 'Q' for FASTQ

  bool          maWithQV          = false;
  bool          maWithDots        = true;
  uint32        maDisplayWidth    = 100;
  uint32        maDisplaySpacing  = 3;

  bool          layWithSequence   = false;

  uint64        genomeSize        = 0;

  char         *outPrefix         = NULL;

  bool          single            = false;

  uint32        minOverlap        = 0;


  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;

  while (arg < argc) {
    if      (strcmp(argv[arg], "-S") == 0) {
      if (arg + 1 < argc)
        seqName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-T") == 0) {
      if (arg + 1 < argc)
        tigName = argv[++arg];
      if (arg + 1 < argc)
        tigVers = atoi(argv[++arg]);
    }

    else if ((strcmp(argv[arg], "-tig") == 0) ||
             (strcmp(argv[arg], "-t") == 0) ||    //  Deprecated!
             (strcmp(argv[arg], "-u") == 0)) {    //  Deprecated too!
      if (arg + 1 < argc)
        decodeRange(argv[++arg], filter.tigIDbgn, filter.tigIDend);
    }

    else if (strcmp(argv[arg], "-unassembled") == 0) {
      filter.dumpAllClasses  = false;
      filter.dumpUnassembled = true;
    }

    else if (strcmp(argv[arg], "-contigs") == 0) {
      filter.dumpAllClasses  = false;
      filter.dumpContigs     = true;
    }


    else if (strcmp(argv[arg], "-repeats") == 0) {
      filter.dumpRepeats     = true;
    }

    else if (strcmp(argv[arg], "-bubbles") == 0) {
      filter.dumpBubbles     = true;
    }

    else if (strcmp(argv[arg], "-circular") == 0) {
      filter.dumpCircular    = true;
    }


    else if (strcmp(argv[arg], "-nreads") == 0) {
      if (arg + 1 < argc)
        filter.minNreads = atoi(argv[++arg]);
      if (arg + 1 < argc)
        filter.maxNreads = atoi(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-length") == 0) {
      if (arg + 1 < argc)
        filter.minLength = atoi(argv[++arg]);
      if (arg + 1 < argc)
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
        char *s = new char [1024];
        snprintf(s, 1024, "ERROR: -coverage needs four values.\n");
        err.push_back(s);
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
      dumpType = DUMP_INFO;
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

    else if (strcmp(argv[arg], "-reverse") == 0)
      useReverse = true;

    else if (strcmp(argv[arg], "-fasta") == 0)
      cnsFormat = 'A';
    else if (strcmp(argv[arg], "-fastq") == 0)
      cnsFormat = 'Q';

    else if (strcmp(argv[arg], "-w") == 0) {
      if (arg + 1 < argc)
        maDisplayWidth = atoi(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-s") == 0) {
      if (arg + 1 < argc)
        maDisplaySpacing = genomeSize = atol(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-o") == 0) {
      if (arg + 1 < argc)
        outPrefix = argv[++arg];
    }

    else if (strcmp(argv[arg], "-single") == 0)
      single = true;

    else if (strcmp(argv[arg], "-thin") == 0) {
      if (arg + 1 < argc)
        minOverlap = atoi(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-sequence") == 0)
      layWithSequence = true;

    //  Errors.

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (seqName == NULL)
    err.push_back("No sequence store (-S option) supplied.\n");

  if (tigName == NULL)
    err.push_back("No tig store (-T option) supplied.\n");

  if (tigVers == -1)
    err.push_back("No tig store version (-T option) supplied.\n");

  if ((outPrefix == NULL) && (dumpType == DUMP_COVERAGE))
    err.push_back("-coverage needs and output prefix (-o option).\n");

  if (dumpType == DUMP_UNSET)
    err.push_back("No DUMP TYPE supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S <seqStore> -T <tigStore> <v> [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "STORE SELECTION (mandatory)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S <seqStore>           path to the sequence store\n");
    fprintf(stderr, "  -T <tigStore> <v>       path to the tigStore, version, to use\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "TIG SELECTION - if nothing specified, all tigs are reported\n");
    fprintf(stderr, "              - all ranges are inclusive.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -tig A[-B]              only dump tigs between ids A and B\n");
    fprintf(stderr, "  -unassembled            only dump tigs that are 'unassembled'\n");
    fprintf(stderr, "  -contigs                only dump tigs that are 'contigs'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -repeats                only dump tigs that are (probably) repeats\n");
    fprintf(stderr, "  -bubbles                only dump tigs that are (probably) bubbles\n");
    fprintf(stderr, "  -circular               only dump tigs that are (probably) circular\n");
    fprintf(stderr, "\n");
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
    fprintf(stderr, "                            -reverse          reverse complement the sequence\n");
    fprintf(stderr, "                            -fasta            report sequences in FASTA format (the default)\n");
    fprintf(stderr, "                            -fastq            report sequences in FASTQ format\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -layout -o name         metadata and read-to-tig position mapping.  -o is mandatory.\n");
    fprintf(stderr, "                          creates two files:\n");
    fprintf(stderr, "                            name.layout.readToTig - read to tig position\n");
    fprintf(stderr, "                            name.layout.tigInfo   - metadata for each tig\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -layout [opts]          the layout of reads in each tig and the consensus sequence if available\n");
    fprintf(stderr, "                          and enabled with option:\n");
    fprintf(stderr, "                            -sequence         also include consensus sequence and quality\n");
    fprintf(stderr, "\n");
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

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  //  Open stores.

  sqStore *seqStore = new sqStore(seqName);
  tgStore *tigStore = new tgStore(tigName, tigVers);

  //  Check that the tig ID range is valid, and fix it if possible.

  uint32   nTigs = tigStore->numTigs();

  if (filter.tigIDend == UINT32_MAX)
    filter.tigIDend = nTigs-1;

  if ((nTigs > 0) && (nTigs <= filter.tigIDend)) {
    fprintf(stderr, "WARNING: adjusting tig ID range from " F_U32 "-" F_U32 " to " F_U32 "-" F_U32 " as there are only " F_U32 " tigs in the store.\n",
            filter.tigIDbgn, filter.tigIDend, filter.tigIDbgn, nTigs-1, nTigs);
    filter.tigIDend = nTigs - 1;
  }

  if (filter.tigIDend < filter.tigIDbgn) {
    fprintf(stderr, "WARNING: adjusting inverted tig ID range -t " F_U32 "-" F_U32 "\n",
            filter.tigIDbgn, filter.tigIDend);
    uint32 x = filter.tigIDend;
    filter.tigIDend = filter.tigIDbgn;
    filter.tigIDbgn = x;
  }

  if ((nTigs > 0) && (nTigs <= filter.tigIDbgn))
    fprintf(stderr, "ERROR: only " F_U32 " tigs in the store (IDs 0-" F_U32 " inclusive); can't dump requested range -t " F_U32 "-" F_U32 "\n",
            nTigs,
            nTigs-1,
            filter.tigIDbgn, filter.tigIDend), exit(1);

  //  Call the dump routine.

  switch (dumpType) {
    case DUMP_STATUS:
      dumpStatus(seqStore, tigStore);
      break;
    case DUMP_TIGS:
      dumpTigs(seqStore, tigStore, filter);
      break;
    case DUMP_CONSENSUS:
      dumpConsensus(seqStore, tigStore, filter, useReverse, cnsFormat);
      break;
    case DUMP_INFO:
      if (outPrefix)
        dumpInfo(seqStore, tigStore, filter, outPrefix);
      else
        dumpLayout(seqStore, tigStore, filter, layWithSequence);
      break;
    case DUMP_MULTIALIGN:
      dumpMultialign(seqStore, tigStore, filter, maWithQV, maWithDots, maDisplayWidth, maDisplaySpacing);
      break;
    case DUMP_SIZES:
      dumpSizes(seqStore, tigStore, filter, genomeSize);
      break;
    case DUMP_COVERAGE:
      dumpCoverage(seqStore, tigStore, filter, outPrefix);
      break;
    case DUMP_DEPTH_HISTOGRAM:
      dumpDepthHistogram(seqStore, tigStore, filter, single, outPrefix);
      break;
    case DUMP_THIN_OVERLAP:
      dumpThinOverlap(seqStore, tigStore, filter, minOverlap);
      break;
    case DUMP_OVERLAP_HISTOGRAM:
      dumpOverlapHistogram(seqStore, tigStore, filter, outPrefix);
      break;
    default:
      break;
  }

  //  Clean up.

  delete tigStore;

  delete seqStore;

  //  Bye.

  exit(0);
}

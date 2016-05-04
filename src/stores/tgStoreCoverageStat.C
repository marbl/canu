
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
 *    src/AS_BAT/computeCoverageStat.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2011-DEC-20 to 2013-AUG-01
 *      are Copyright 2011-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2012-JUL-18
 *      are Copyright 2012 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2015-AUG-14
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"

#include <algorithm>

using namespace std;


//  This program will recompute the coverage statistic for all unitigs in the tigStore.
//  It replaces at least four implementations (AS_CGB, AS_BOG, AS_BAT, AS_CGW).
//
//  Notes from the AS_CGB version:
//
//  Rho is the number of bases in the chunk between the first fragment arrival and the last fragment
//  arrival.  It is the sum of the fragment overhangs in the chunk.  For intuitive purposes you can
//  think of it as the length of the chunk minus the length of the last fragment (if that last isn't
//  contained).  Thus a singleton chunk has a rho equal to zero.
//
//  A singleton chunk provides no information as to its local fragment arrival rate.  We need at
//  least two closely spaced fragments that are randomly sampled from the chunk to get a local
//  estimate of the fragment arrival rate.
//
//  The local arrival rate of fragments in the chunk is:
//    arrival_rate_local = (nfrag_randomly_sampled_in_chunk - 1) / rho
//
//  The arrival distance of fragments in the chunk is the reciprocal of the last formula:
//    arrival_distance_local = rho / (nfrag_randomly_sampled_in_chunk - 1)
//
//  Note a problem with this formula is that a singleton chunk has a coverage discriminator
//  statistic of 0/0.
//
//  The formula for the coverage discriminator statistic for the chunk is:
//    (arrival_rate_global / arrival_rate_local - ln(2)) * (nfrag_randomly_sampled_in_chunk - 1)
//
//  The division by zero singularity cancels out to give the formula:
//    (arrival_rate_global * rho - ln(2) * (nfrag_randomly_sampled_in_chunk - 1)
//
//  The coverage discriminator statistic should be positive for single coverage, negative for
//  multiple coverage, and near zero for indecisive.
//
//  ADJUST_FOR_PARTIAL_EXCESS: The standard statistic gives log likelihood ratio of expected depth
//  vs twice expected depth; but when enough fragments are present, we can actually test whether
//  depth exceeds expected even fractionally; in deeply sequenced datasets (e.g. bacterial genomes),
//  this has been observed for repetitive segments.
//



double    ln2 = 0.69314718055994530941723212145818;
double    globalArrivalRate = 0.0;

bool     *isNonRandom = NULL;
uint32   *readLength  = NULL;

bool      leniant = false;


//  No frags -> 1
//  One frag -> 1

double
computeRho(tgTig *tig) {
  int32  minBgn  = INT32_MAX;
  int32  maxEnd  = INT32_MIN;
  int32  fwdRho  = INT32_MIN;
  int32  revRho  = INT32_MAX;

  //  We compute the two rho's using the first definition above - distance between the first
  //  and last fragment arrival.  This changes based on the orientation of the unitig, so we
  //  return the average of those two.

  for (uint32 i=0; i<tig->numberOfChildren(); i++) {
    tgPosition  *pos = tig->getChild(i);

    minBgn = MIN(minBgn, pos->min());
    maxEnd = MAX(maxEnd, pos->max());

    fwdRho = MAX(fwdRho, pos->min());  //  largest begin coord
    revRho = MIN(revRho, pos->max());  //  smallest end coord
  }

  if ((leniant == false) && (minBgn != 0)) {
    fprintf(stderr, "tig %d doesn't begin at zero.  Layout:\n", tig->tigID());
    tig->dumpLayout(stderr);
  }
  if (leniant == false)
    assert(minBgn == 0);

  fwdRho = fwdRho - minBgn;
  revRho = maxEnd - revRho;

  assert(fwdRho >= 0);
  assert(revRho >= 0);

  //  AS_CGB is using the begin of the last fragment as rho

  return((fwdRho + revRho) / 2.0);
}


uint32
numRandomFragments(tgTig *tig) {
  uint32  numRand = 0;

  for (uint32 ii=0; ii<tig->numberOfChildren(); ii++)
    if (isNonRandom[tig->getChild(ii)->ident()] == false)
      numRand++;

  return(numRand);
}



double
getGlobalArrivalRate(tgStore         *tigStore,
                     FILE            *outSTA,
                     uint64           genomeSize,
                     bool             useN50) {
  double   globalRate  = 0;
  double   recalRate   = 0;

  double   sumRho      = 0;

  int32    arLen       = 0;
  double  *ar          = NULL;
  uint32   NF;
  uint64   totalRandom = 0;
  uint64   totalNF     = 0;
  int32    BIG_SPAN    = 10000;

  int32    big_spans_in_unitigs   = 0; // formerly arMax

  // Go through all the unitigs to sum rho and unitig arrival frags

  uint32 *allRho = new uint32 [tigStore->numTigs()];

  for (uint32 i=0; i<tigStore->numTigs(); i++) {
    tgTig  *tig = tigStore->loadTig(i);

    allRho[i] = 0;

    if (tig == NULL)
      continue;

    double rho       = computeRho(tig);
    int32  numRandom = numRandomFragments(tig);

    tigStore->unloadTig(i);

    sumRho                 += rho;
    big_spans_in_unitigs   += (int32) (rho / BIG_SPAN);  // Keep integral portion of fraction.
    totalRandom            += numRandom;
    totalNF                +=  (numRandom == 0) ? (0) : (numRandom - 1);

    allRho[i] = rho;
  }

  // Here is a rough estimate of arrival rate.
  // Use (number frags)/(unitig span) unless unitig span is zero; then use (reads)/(genome).
  // Here, (number frags) includes only random reads and omits the last read of each unitig.
  // Here, (sumRho) is total unitig span omitting last read of each unitig.

  if (genomeSize > 0) {
    globalRate = totalRandom / (double)genomeSize;

  } else {
    if (sumRho > 0)
      globalRate = totalNF / sumRho;
  }

  fprintf(outSTA, "BASED ON ALL UNITIGS:\n");
  fprintf(outSTA, "sumRho:                           %.0f\n", sumRho);
  fprintf(outSTA, "totalRandomFrags:                 "F_U64"\n", totalRandom);
  fprintf(outSTA, "Supplied genome size              "F_U64"\n", genomeSize);
  fprintf(outSTA, "Computed genome size:             %.2f\n", totalRandom / globalRate);
  fprintf(outSTA, "Calculated Global Arrival rate:   %f\n", globalRate);

  // Stop here and return the rough estimate under some circumstances.
  // *) If user suppled a genome size, we are done.
  // *) No unitigs.

  if (genomeSize > 0 || tigStore->numTigs()==0) {
    delete [] allRho;
    return(globalRate);
  }

  //  Calculate rho N50

  double rhoN50 = 0;
  if (useN50) {
    uint32 growUntil = sumRho / 2; // half is 50%, needed for N50
    uint64 growRho = 0;
    sort (allRho, allRho+tigStore->numTigs());
    for (uint32 i=tigStore->numTigs(); i>0; i--) { // from largest to smallest unitig...
      rhoN50 = allRho[i-1];
      growRho += rhoN50;
      if (growRho >= growUntil)
        break; // break when sum of rho > 50%
    }
  }
  delete [] allRho;

  //  Try for a better estimate based on just unitigs larger than N50.

  if (useN50) {
    double keepRho = 0;
    double keepNF = 0;
    for (uint32 i=0; i<tigStore->numTigs(); i++) {
      tgTig  *tig = tigStore->loadTig(i);

      if (tig == NULL)
        continue;

      double  rho = computeRho(tig);

      if (rho < rhoN50)
        continue; // keep only rho from unitigs > N50

      int32 numRandom =   numRandomFragments(tig);

      keepNF     +=  (numRandom == 0) ? (0) : (numRandom - 1);
      keepRho    +=  rho;

      tigStore->unloadTig(i);
    }

    fprintf(outSTA, "BASED ON UNITIGS > N50:\n");
    fprintf(outSTA, "rho N50:                          %.0f\n", rhoN50);
    if (keepRho > 1) {   // the cutoff 1 is arbitrary but larger than 0.0f
      globalRate = keepNF / keepRho;
      fprintf(outSTA, "sumRho:                           %.0f\n", keepRho);
      fprintf(outSTA, "totalRandomFrags:                 %.0f\n", keepNF);
      fprintf(outSTA, "Computed genome size:             %.2f\n", totalRandom / globalRate);
      fprintf(outSTA, "Calculated Global Arrival rate:   %f\n", globalRate);
      return (globalRate);
    } else {
      fprintf(outSTA, "It did not work to re-estimate using the N50 method.\n");
    }
  }

  //  Recompute based on just big unitigs. Big is 10Kbp.
  double BIG_THRESHOLD = 0.5;
  int32 big_spans_in_rho = (int32) (sumRho / BIG_SPAN);
  fprintf(outSTA, "Size of big spans is %d\n", BIG_SPAN);
  fprintf(outSTA, "Number of big spans in unitigs is %d\n", big_spans_in_unitigs);
  fprintf(outSTA, "Number of big spans in sum-of-rho is %d\n", big_spans_in_rho);
  fprintf(outSTA, "Ratio required for re-estimate is %f\n", BIG_THRESHOLD);
  if ( (big_spans_in_unitigs / big_spans_in_rho) <= BIG_THRESHOLD) {
    fprintf(outSTA, "Too few big spans to re-estimate using the big spans method.\n");
    return(globalRate);
  }
  //  The test above is a rewrite of the former version, where arMax=big_spans_in_unitigs...
  //  if (arMax <= sumRho / 20000)


  ar = new double [big_spans_in_unitigs];

  for (uint32 i=0; i<tigStore->numTigs(); i++) {
    tgTig  *tig = tigStore->loadTig(i);

    if (tig == NULL)
      continue;

    double  rho = computeRho(tig);

    if (rho <= BIG_SPAN)
      continue;

    int32   numRandom        = numRandomFragments(tig);
    double  localArrivalRate = numRandom / rho;
    uint32  rhoDiv10k        = rho / BIG_SPAN;

    assert(0 < rhoDiv10k);

    for (uint32 aa=0; aa<rhoDiv10k; aa++)
      ar[arLen++] = localArrivalRate;

    assert(arLen <= big_spans_in_unitigs);

    sort(ar, ar + arLen);

    double  maxDiff    = 0.0;
    uint32  maxDiffIdx = 0;

    uint32  idx = arLen / 10;

    double  arc = ar[idx];  //  ar[i]
    double  ard = 0;        //  Current difference in ar[i] - ar[i-1]
    double  arp = ar[idx];  //  ar[i-1]

    for (; idx < arLen / 2; idx++) {
      arc = ar[idx];
      ard = arc - arp;
      arp = arc;

      maxDiff = MAX(maxDiff, ard);
    }

    maxDiff    *= 2.0;
    maxDiffIdx  = arLen - 1;

    for (; idx < arLen; idx++) {
      arc = ar[idx];
      ard = arc - arp;
      arp = arc;

      if (ard > maxDiff) {
        maxDiff    = ard;
        maxDiffIdx = idx - 1;
        break;
      }
    }

    double recalRate = 0;

    recalRate  =                ar[arLen * 19 / 20];
    recalRate  = MIN(recalRate, ar[arLen *  1 / 10] * 2.0);
    recalRate  = MIN(recalRate, ar[arLen *  1 /  2] * 1.25);
    recalRate  = MIN(recalRate, ar[maxDiffIdx]);

    globalRate = MAX(globalRate, recalRate);

    tigStore->unloadTig(i);
  }

  delete [] ar;

  fprintf(outSTA, "BASED ON BIG SPANS IN UNITIGS:\n");
  fprintf(outSTA, "Computed genome size:             %.2f (reestimated)\n", totalRandom / globalRate);
  fprintf(outSTA, "Calculated Global Arrival rate:   %f (reestimated)\n", globalRate);

  return(globalRate);
}








int
main(int argc, char **argv) {
  char             *gkpName    = NULL;
  char             *tigName    = NULL;
  int32             tigVers    = -1;

  int64             genomeSize = 0;

  char             *outPrefix  = NULL;
  FILE             *outLOG     = NULL;
  FILE             *outSTA     = NULL;

  uint32            bgnID      = 0;
  uint32            endID      = 0;

  bool              doUpdate   = true;
  bool              use_N50    = true;

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-s") == 0) {
      genomeSize = atol(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-n") == 0) {
      doUpdate = false;

    } else if (strcmp(argv[arg], "-u") == 0) {
      use_N50 = false;

    } else if (strcmp(argv[arg], "-L") == 0) {
      leniant = true;

    } else {
      err++;
    }

    arg++;
  }

  if (gkpName == NULL)
    err++;
  if (tigName == NULL)
    err++;
  if (outPrefix == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -T tigStore version -o output-prefix [-s genomeSize] ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G <G>     Mandatory, path G to a gkpStore directory.\n");
    fprintf(stderr, "  -T <T> <v> Mandatory, path T to a tigStore, and version V.\n");
    fprintf(stderr, "  -o <name>  Mandatory, prefix for output files.\n");
    fprintf(stderr, "  -s <S>     Optional, assume genome size S.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -n         Do not update the tigStore (default = do update).\n");
    fprintf(stderr, "  -u         Do not estimate based on N50 (default = use N50).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L         Be leniant; don't require reads start at position zero.\n");
    fprintf(stderr, "\n");

    if (gkpName == NULL)
      fprintf(stderr, "No gatekeeper store (-G option) supplied.\n");

    if (tigName == NULL)
      fprintf(stderr, "No input tigStore (-T option) supplied.\n");

    if (outPrefix == NULL)
      fprintf(stderr, "No output prefix (-o option) supplied.\n");

    exit(1);
  }

  //  Open output files first, so we can fail before getting too far along.

  {
    char  outName[FILENAME_MAX];

    errno = 0;

    sprintf(outName, "%s.log", outPrefix);

    outLOG = fopen(outName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", outName, strerror(errno)), exit(1);

    sprintf(outName, "%s.stats", outPrefix);

    outSTA = fopen(outName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", outName, strerror(errno)), exit(1);
  }

  //
  //  Load fragment data
  //

  fprintf(stderr, "Opening gkpStore '%s'\n", gkpName);

  gkStore *gkpStore = gkStore::gkStore_open(gkpName, gkStore_readOnly);

  fprintf(stderr, "Reading read lengths and randomness for %u reads.\n",
          gkpStore->gkStore_getNumReads());

  isNonRandom = new bool   [gkpStore->gkStore_getNumReads() + 1];
  readLength  = new uint32 [gkpStore->gkStore_getNumReads() + 1];

  for (uint32 ii=0; ii<=gkpStore->gkStore_getNumReads(); ii++) {
    gkRead      *read = gkpStore->gkStore_getRead(ii);
    gkLibrary   *libr = gkpStore->gkStore_getLibrary(read->gkRead_libraryID());

    isNonRandom[ii] = libr->gkLibrary_isNonRandom();
    readLength[ii]  = read->gkRead_sequenceLength();
  }

  fprintf(stderr, "Closing gkpStore.\n");

  gkpStore->gkStore_close();
  gkpStore = NULL;

  //
  //  Open tigs.  Kind of important to do this.
  //

  fprintf(stderr, "Opening tigStore '%s'\n", tigName);

  tgStore *tigStore     = new tgStore(tigName, tigVers, tgStoreModify);

  if (endID == 0)
    endID = tigStore->numTigs();

  //
  //  Compute global arrival rate.  This ain't cheap.
  //

  fprintf(stderr, "Computing global arrival rate.\n");

  double  globalRate = getGlobalArrivalRate(tigStore, outSTA, genomeSize, use_N50);

  //
  //  Compute coverage stat for each unitig, populate histograms, write logging.
  //

  //  Three histograms were made, one for length, coverage stat and arrival distance.  The
  //  histograms included
  //
  //      columns: sum, cumulative_sum, cumulative_fraction, min, average, max
  //      rows:    reads, rs reads, nr reads, bases, rho, arrival, discriminator
  //
  //  Most of those are not defined or null (nr reads? currently zero), and the histograms
  //  in general aren't useful anymore.  They were used to help decide if the genome size
  //  was incorrect, causing too many non-unique unitigs.
  //
  //  They were removed 13 Aug 2015.

  fprintf(stderr, "Computing coverage stat for tigs %u-%u.\n", bgnID, endID-1);

  for (uint32 i=bgnID; i<endID; i++) {
    tgTig  *tig = tigStore->loadTig(i);

    if (tig == NULL)
      continue;

    int32   tigLength = tig->length(true);
    int32   numFrags  = tig->numberOfChildren();
    int32   numRandom = numRandomFragments(tig);

    double  rho       = computeRho(tig);

    double  covStat   = 0.0;
    double  arrDist   = 0.0;

    if (numRandom > 1)
      arrDist = rho / (numRandom - 1);

    if ((numRandom > 0) &&
        (globalRate > 0.0))
      covStat = (rho * globalRate) - (ln2 * (numRandom - 1));

    if (i == bgnID)
      fprintf(outLOG, "     tigID        rho    covStat    arrDist\n");
    fprintf(outLOG, "%10u %10.2f %10.2f %10.2f\n", tig->tigID(), rho, covStat, arrDist);

#undef ADJUST_FOR_PARTIAL_EXCESS
#ifdef ADJUST_FOR_PARTIAL_EXCESS
    //  Straight from AS_CGB_all.h, not refactored
    if(rho>0&&global_fragment_arrival_rate>0.f){
      double lambda = global_fragment_arrival_rate * rho;
      double zscore = ((number_of_randomly_sampled_fragments_in_chunk -1)-lambda) / sqrt(lambda);
      double p = .5 - erf(zscore/sqrt2)*.5;
      if(coverage_statistic>5 && p < .001){
        fprintf(outSTA, "Standard unitigger a-stat is %f, but only %e chance of this great an excess of fragments: obs = %d, expect = %g rho = " F_S64 " Will reset a-stat to 1.5\n",
                coverage_statistic,p,
                number_of_randomly_sampled_fragments_in_chunk-1,
                lambda,rho);
        covStat = 1.5;
      }
    }
#endif

    if (doUpdate)
      tigStore->setCoverageStat(tig->tigID(), covStat);

    tigStore->unloadTig(tig->tigID());
  }


  fclose(outLOG);

  delete [] isNonRandom;
  delete [] readLength;

  delete tigStore;

  exit(0);
}

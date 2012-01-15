
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id: computeCoverageStat.C,v 1.4 2012-01-15 23:49:34 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "MultiAlign.h"
#include "MultiAlignStore.h"

#include "AS_CGB_histo.h"

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
uint32   *fragLength  = NULL;


//  No frags -> 1
//  One frag -> 1

double
computeRho(MultiAlignT *ma) {
  int32  minBgn  = INT32_MAX;
  int32  maxEnd  = INT32_MIN;
  int32  fwdRho  = INT32_MIN;
  int32  revRho  = INT32_MAX;

  //  We compute the two rho's using the first definition above - distance between the first
  //  and last fragment arrival.  This changes based on the orientation of the unitig, so we
  //  return the average of those two.

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos  *frg = GetIntMultiPos(ma->f_list, i);

    if (frg->position.bgn < frg->position.end) {
      minBgn = MIN(minBgn, frg->position.bgn);
      maxEnd = MAX(maxEnd, frg->position.end);

      fwdRho = MAX(fwdRho, frg->position.bgn);  //  largest begin coord
      revRho = MIN(revRho, frg->position.end);  //  smallest end coord
    } else {
      minBgn = MIN(minBgn, frg->position.end);
      maxEnd = MAX(maxEnd, frg->position.bgn);
      
      fwdRho = MAX(fwdRho, frg->position.end);
      revRho = MIN(revRho, frg->position.bgn);
    }
  }

  if (minBgn != 0) {
    fprintf(stderr, "ERROR: unitig %d doesn't begin at zero.  Layout:\n", ma->maID);
    for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
      IntMultiPos  *frg = GetIntMultiPos(ma->f_list, i);
      fprintf(stderr, "  %10"F_IIDP" %5"F_U32P" %5"F_U32P"\n",
              frg->ident, frg->position.bgn, frg->position.end);
    }
  }
  assert(minBgn == 0);

  fwdRho = fwdRho - minBgn;
  revRho = maxEnd - revRho;

  assert(fwdRho >= 0);
  assert(revRho >= 0);

  //  AS_CGB is using the begin of the last fragment as rho

  return((fwdRho + revRho) / 2.0);
}


int32
numRandomFragments(MultiAlignT *ma) {
  int32  numRand = 0;

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos  *frg = GetIntMultiPos(ma->f_list, i);

    if (isNonRandom[frg->ident] == false)
      numRand++;
  }

  return(numRand);
}



double
getGlobalArrivalRate(MultiAlignStore *tigStore,
                     FILE            *outSTA,
                     uint64           genomeSize) {
  double   globalRate = 0;
  double   recalRate  = 0;

  double   sumRho     = 0;

  int32    arMax   = 0;
  int32    arLen   = 0;
  double  *ar      = NULL;

  uint64   totalRandom = 0;
  uint64   totalNF     = 0;

  // Go through all the unitigs to sum rho and unitig arrival frags

  for (uint32 i=0; i<tigStore->numUnitigs(); i++) {
    MultiAlignT  *ma = tigStore->loadMultiAlign(i, TRUE);

    if (ma == NULL)
      continue;

    double rho       = computeRho(ma);
    int32  numRandom = numRandomFragments(ma);

    sumRho  += rho;
    arMax   += rho / 10000;

    totalRandom += numRandom;
    totalNF     += (numRandom == 0) ? (0) : (numRandom - 1);
  }

  if (genomeSize > 0)
    globalRate = totalRandom / (double)genomeSize;

  if (sumRho > 0)
    globalRate = totalNF / sumRho;

  fprintf(outSTA, "\n");
  fprintf(outSTA, "Calculated Global Arrival rate:   %f\n", globalRate);
  fprintf(outSTA, "\n");
  fprintf(outSTA, "sumRho:                           %.0f\n", sumRho);
  fprintf(outSTA, "totalRandomFrags:                 "F_U64"\n", totalRandom);
  fprintf(outSTA, "\n");
  fprintf(outSTA, "Supplied genome size              "F_U64"\n", genomeSize);
  fprintf(outSTA, "Computed genome size:             %.2f\n", totalRandom / globalRate);
  fprintf(outSTA, "\n");

  if (genomeSize > 0)
    return(globalRate);

  if (arMax <= sumRho / 20000)
    return(globalRate);

  //  Recompute based on just big unitigs.

  ar = new double [arMax];

  for (uint32 i=0; i<tigStore->numUnitigs(); i++) {
    MultiAlignT  *ma = tigStore->loadMultiAlign(i, TRUE);

    if (ma == NULL)
      continue;

    double  rho = computeRho(ma);

    if (rho <= 10000)
      continue;

    int32   numRandom        = numRandomFragments(ma);
    double  localArrivalRate = numRandom / rho;
    uint32  rhoDiv10k        = rho / 10000;

    assert(0 < rhoDiv10k);

    for (uint32 i=0; i<rhoDiv10k; i++)
      ar[arLen++] = localArrivalRate;

    assert(arLen <= arMax);

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
  }

  delete [] ar;

  fprintf(outSTA, "Calculated Global Arrival rate:   %f (reestimated)\n", globalRate);

  return(globalRate);
}








int
main(int argc, char **argv) {
  char             *gkpName      = NULL;
  char             *tigName      = NULL;
  int32             tigVers      = -1;

  int64             genomeSize   = 0;

  char              outName[FILENAME_MAX];
  char             *outPrefix = NULL;
  FILE             *outCGA = NULL;
  FILE             *outLOG = NULL;
  FILE             *outSTA = NULL;

  AS_IID            bgnID = 0;
  AS_IID            endID = 0;

  MultiAlignStore  *tigStore = NULL;
  gkStore          *gkpStore = NULL;

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-s") == 0) {
      genomeSize = atol(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }

  if (gkpName == NULL)
    err++;
  if (tigName == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -g         Mandatory path to a gkpStore.\n");
    fprintf(stderr, "  -t         Mandatory path to a tigStore (can exist or not).\n");
    fprintf(stderr, "  -s         Optional genome size.\n");
    fprintf(stderr, "  -o         Output prefix; will create prefix.cga.0.\n");

    if (gkpName == NULL)
      fprintf(stderr, "No gatekeeper store (-g option) supplied.\n");

    if (tigName == NULL)
      fprintf(stderr, "No input tigStore (-t option) supplied.\n");

    if (outPrefix == NULL)
      fprintf(stderr, "No output prefix (-o option) supplied.\n");

    exit(1);
  }

  gkpStore     = new gkStore(gkpName, FALSE, FALSE);
  tigStore     = new MultiAlignStore(tigName, tigVers, 0, 0, TRUE, TRUE, FALSE);

  if (endID == 0)
    endID = tigStore->numUnitigs();

  errno = 0;

  sprintf(outName, "%s.cga.0", outPrefix);

  outCGA = fopen(outName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outName, strerror(errno)), exit(1);

  sprintf(outName, "%s.log", outPrefix);

  outLOG = fopen(outName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outName, strerror(errno)), exit(1);

  fprintf(outLOG, "tigID\trho\tcovStat\tarrDist\n");

  sprintf(outName, "%s.stats", outPrefix);

  outSTA = fopen(outName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outName, strerror(errno)), exit(1);



  //
  //  Load fragment data
  //

  double  globalRate     = 0;

  isNonRandom = new bool   [gkpStore->gkStore_getNumFragments() + 1];
  fragLength  = new uint32 [gkpStore->gkStore_getNumFragments() + 1];

  gkStream   *fs = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);
  gkFragment  fr;

  while(fs->next(&fr)) {
    uint32 iid = fr.gkFragment_getReadIID();

    isNonRandom[iid] = fr.gkFragment_getIsNonRandom();
    fragLength[iid]  = fr.gkFragment_getClearRegionLength(AS_READ_CLEAR_OBTCHIMERA);

    if ((iid % 10000000) == 0)
      fprintf(stderr, "Loading fragment information %9d out of %9d\n", iid, gkpStore->gkStore_getNumFragments());
  }
  
  delete fs;

  //
  //  Compute global arrival rate.  This ain't cheap.
  //

  globalRate = getGlobalArrivalRate(tigStore, outSTA, genomeSize);

  //
  //  Compute coverage stat for each unitig, populate histograms, write logging.
  //

  MyHistoDataType  z;

  memset(&z, 0, sizeof(MyHistoDataType));

  Histogram_t *len = create_histogram(500, 500, 0, TRUE);
  Histogram_t *cvg = create_histogram(500, 500, 0, TRUE);
  Histogram_t *arv = create_histogram(500, 500, 0, TRUE);

  extend_histogram(len, sizeof(MyHistoDataType), myindexdata, mysetdata, myaggregate, myprintdata);
  extend_histogram(cvg, sizeof(MyHistoDataType), myindexdata, mysetdata, myaggregate, myprintdata);
  extend_histogram(arv, sizeof(MyHistoDataType), myindexdata, mysetdata, myaggregate, myprintdata);

  for (uint32 i=bgnID; i<endID; i++) {
    MultiAlignT  *ma = tigStore->loadMultiAlign(i, TRUE);

    if (ma == NULL)
      continue;

    int32   tigLength = GetMultiAlignLength(ma);
    int32   numFrags  = GetNumIntMultiPoss(ma->f_list);
    int32   numRandom = numRandomFragments(ma);

    double  rho       = computeRho(ma);

    double  covStat   = 0.0;
    double  arrDist   = 0.0;

    if (numRandom > 1)
      arrDist = rho / (numRandom - 1);

    if ((numRandom > 0) &&
        (globalRate > 0.0))
      covStat = (rho * globalRate) - (ln2 * (numRandom - 1));

    fprintf(outLOG, F_U32"\t%.2f\t%.2f\t%.2f\n", ma->maID, rho, covStat, arrDist);

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

    z.nsamples     = 1;
    z.sum_frags    = z.min_frags    = z.max_frags    = numFrags;
    z.sum_bp       = z.min_bp       = z.max_bp       = tigLength;
    z.sum_rho      = z.min_rho      = z.max_rho      = rho;
    z.sum_discr    = z.min_discr    = z.max_discr    = covStat;
    z.sum_arrival  = z.min_arrival  = z.max_arrival  = arrDist;
    z.sum_rs_frags = z.min_rs_frags = z.max_rs_frags = numRandom;
    z.sum_nr_frags = z.min_nr_frags = z.max_nr_frags = 0;

    add_to_histogram(len, tigLength, &z);
    add_to_histogram(cvg, covStat,   &z);
    add_to_histogram(arv, arrDist,   &z);

    tigStore->setUnitigCoverageStat(ma->maID, covStat);
  }


  fprintf(outCGA, "\n");
  fprintf(outCGA, "\n");
  fprintf(outCGA, "Unitig Length\n");
  fprintf(outCGA, "label\tsum\tcummulative\tcummulative   min  average  max\n");
  fprintf(outCGA, "     \t   \t sum       \t fraction\n");
  print_histogram(outCGA, len, 0, 1);

  fprintf(outCGA, "\n");
  fprintf(outCGA, "\n");
  fprintf(outCGA, "Unitig Coverage Stat\n");
  fprintf(outCGA, "label\tsum\tcummulative\tcummulative   min  average  max\n");
  fprintf(outCGA, "     \t   \t sum       \t fraction\n");
  print_histogram(outCGA, cvg, 0, 1);

  fprintf(outCGA, "\n");
  fprintf(outCGA, "\n");
  fprintf(outCGA, "Unitig Arrival Distance\n");
  fprintf(outCGA, "label\tsum\tcummulative\tcummulative   min  average  max\n");
  fprintf(outCGA, "     \t   \t sum       \t fraction\n");
  print_histogram(outCGA, arv, 0, 1);

  free_histogram(len);
  free_histogram(cvg);
  free_histogram(arv);



  fclose(outCGA);
  fclose(outLOG);


  delete [] isNonRandom;
  delete [] fragLength;

  delete gkpStore;
  delete tigStore;
}

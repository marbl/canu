
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id$";

#include "AS_global.H"
#include "AS_OVS_overlap.H"
#include "AS_OVS_overlapStore.H"

#include "overlapStore.H"

#include <cmath>

//  Estimate the genome length based on overlap statistics.
//
//  Estimates genome length for a sequencing project from a low or intermediate coverage dataset,
//  discounting fragments whose coverage is high enough to suggest they may be repetitive and
//  avoiding near duplicate fragments in case of artifact.
//
//  Parameters
//
//    -K <int>              Exclude fragments with more than this many overlaps
//    -first <int>          First fragment to analyze; default=1
//    -last <int>           Last fragment to analyze; default=100000
//    -into <int>           Length of skipped prefix of analyzed fragments; default=50
//    -window <int>         Analysis window width; default=100
//    -minovl <int>         Minimum overlap length; default=40
//
//  Assumptions:
//
//    No deleted fragments in fragStore
//    K value at least twice the average depth of coverage in unique regions
//    Value of into+window is less than length of most fragments
//    Fragments in [first,last] are typical of entire dataset
//    Single genome with little contamination
//
//  Credit
//
//    Initial version described by Art Delcher (Kirkness et al 2004, Science).
//    Modifications and coding in perl by Aaron Halpern (1 May 2006) (AS_RUN/dataQC/glen_est_truncadjusted.pl)
//    Current implementation by Brian Walenz (24 Sep 2010).
//
//
//  Estimation of genome length (Kirkness et al 2004):
//
//  From samples of reads (200,000 or 1 million), the numbers and positions of overlaps that began 5
//  or more bases downstream from the 5' end of the read were computed. In order to eliminate reads
//  from repetitive regions, only reads with fewer than 5 overlaps beginning in this region were
//  considered. For a given window in this region, e.g., the first 100 bases of the region, the
//  number of overlaps beginning in that window is tabulated for each read (Table S1). From this,
//  the Poisson mean = 0 * 0.775082 + 1 * 0.196351 + ... + 4 * 0.000181 = 0.255094 is
//  calculated. Equating lambda to np, the mean of the binomial distribution with n = 6,225,420
//  reads and probability of a read beginning in a window of length 100 being p = 100/G , where G is
//  the genome length, yields G = 100 * n / lambda = 2.43 * 10^9.
//
//
//  Table S1. Overlap of sequence reads. Overlaps in a window of the first 100 bases, beginning 5 or
//  more bases downstream from the 5' end of the read were computed. In order to eliminate reads
//  from repetitive regions, only reads with fewer than 5 overlaps beginning in this region were
//  considered.
//
//  # Overlaps   # Reads   Fraction
//   0           687,875   0.775082
//   1           174,259   0.196351
//   2            23,228   0.026173
//   3             1,964   0.002213
//   4               161   0.000181
//  Total        887,487
//



uint32 *factorial = NULL;

//  Compute the mean overlap count after truncation given a poisson distrib with given mean.
//
//  Compute mean by considering full fragment from 1 to k overlaps.
//
//  (comments in the original were inlined, which made it tough to see)
//
//  For i fragments starting during the read, the probability if i is 'prob',
//  and the expected number occuring within the first 100bp is 'y = i / frac',
//  so the contribution to the expected value is 'L'.  But since we
//  exclude the cases affected trunction, need to compute the denominator.
//
double
estimate(double  mean,
         double  frac,
         uint32  overlapLimit) {

  double lambda = mean * frac;  //  this is the true mean in the full fragment
  double L = 0;                 //  this is the estimated mean in the first 100 bp
  double P = 0;                 //  normalizing probability ... see below

  for (uint32 i=0; i<overlapLimit; i++) {
    double  prob = exp(i * log(lambda)) * exp(-lambda) / factorial[i];

    L       += i / frac * prob;
    P       += prob;
  }

  return(L/P);
}


//  Need to adjust naive estimate of lambda to deal with the fact that if 'overlapLimit' is small
//  enough, we have filtered out a significant amount of upper tail of unique coverage.  The
//  following subroutine computes, by more or less Newton's method, the poisson mean that, when
//  truncated >= 'overlapLimit', would give (truncated) mean == naive lambda estimate
//
double
computeTruncatedMean(uint32  overlapLimit,
                     double  lambda,           //  goal value for fitted parameter
                     uint32  into,
                     uint32  avgFragLength,    //  unused
                     uint32  minOvl,
                     uint32  windowSize) {
  double precision = 0.00001;

  //  scale-factor for how much of read is evaluated for k and how much counts ovls
  double frac = (avgFragLength - minOvl - into) / windowSize;

  double low  = 0.00001;
  double high = lambda * 10;
  double mid  = (high + low) / 2;

  double highest = estimate(high, frac, overlapLimit);
  double lowest  = estimate(low,  frac, overlapLimit);
  double midest  = estimate(mid,  frac, overlapLimit);

  if (highest < lambda){
    fprintf(stderr, "Trouble getting initial bounds!\n");
    exit(1);
  }

  double err = fabs(midest - lambda);

  while (err > precision){
    if (midest > lambda){
      highest = midest;
      high    = mid;
    } else {
      lowest = midest;
      low    = mid;
    }

    mid    = (high + low) / 2;
    midest = estimate(mid, frac, overlapLimit);
    err    = abs(midest - lambda);
  }

  return(mid);
}






void
estimateGenomeLength(char    *ovlName,
                     char    *gkpName,
                     uint32   overlapLimit,  //  K
                     uint32   bgnIID,
                     uint32   endIID,
                     uint32   into,
                     uint32   windowSize,
                     uint32   minOvl) {

  if (overlapLimit < 1)
    fprintf(stderr, "error:  overlapLimit must be at least 1.\n"), exit(1);

  ////////////////////////////////////////

  factorial = new uint32 [overlapLimit];

  factorial[0] = 1;
  factorial[1] = 1;
  for (uint32 i=2; i<overlapLimit; i++)
    factorial[i] = factorial[i-1] * i;

  ////////////////////////////////////////

  uint64         avgFragLength = 0;
  uint32         numFrag       = 0;

  gkStore       *gkpStore = new gkStore(gkpName, FALSE, FALSE);

  if (bgnIID < 1)
    bgnIID = 1;
  if (gkpStore->gkStore_getNumLibraries() < endIID)
    endIID = gkpStore->gkStore_getNumFragments();

  gkStream      *gkpStream = new gkStream(gkpStore, bgnIID, endIID, GKFRAGMENT_INF);
  gkFragment     gkpFrag;

  while (gkpStream->next(&gkpFrag)) {
    if (gkpFrag.gkFragment_getIsDeleted())
      continue;

    avgFragLength += gkpFrag.gkFragment_getClearRegionLength();
    numFrag++;
  }

  avgFragLength /= numFrag;

  delete gkpStream;

  ////////////////////////////////////////

  OverlapStore  *ovlStore   = AS_OVS_openOverlapStore(ovlName);

  uint32         overlapLen = 0;
  uint32         overlapMax = 1024 * 1024;
  OVSoverlap    *overlap    = new OVSoverlap [overlapMax];

  uint32         readsWithOverlap  = 0;
  uint32         nonRepeatReads    = 0;
  uint32         nonRepeatOverlaps = 0;

  AS_OVS_setRangeOverlapStore(ovlStore, bgnIID, endIID);

  overlapLen = AS_OVS_readOverlapsFromStore(ovlStore, overlap, overlapMax, AS_OVS_TYPE_OVL);
  while (overlapLen > 0) {
    uint32  novlT = 0;  //  num overlaps total
    uint32  novlW = 0;  //  num overlaps in the window

    for (uint32 i=0; i<overlapLen; i++) {
      if (overlap[i].dat.ovl.a_hang <= into)
        //  Skip this overlap, too similar to another overlap, or negative hang
        continue;

      novlT++;

      if (overlap[i].dat.ovl.a_hang <= into + windowSize)
        novlW++;
    }

    if (novlT > 0)
      readsWithOverlap++;

    if (novlT < overlapLimit) {
      nonRepeatReads++;
      nonRepeatOverlaps += novlW;
    }

    overlapLen = AS_OVS_readOverlapsFromStore(ovlStore, overlap, overlapMax, AS_OVS_TYPE_OVL);
  }

  fprintf(stderr, "%0.2f%% of fragments marked repeat ("F_U32" out of "F_U32").\n\n",
          100.0 - 100.0 * (double)nonRepeatReads / (double)numFrag,
          numFrag - nonRepeatReads,
          numFrag);

  double lambda = nonRepeatOverlaps / (double)(numFrag - readsWithOverlap + nonRepeatReads);

  fprintf(stderr, "Final Estimate\n");
  fprintf(stderr, "  nFrags   "F_S32"\n", gkpStore->gkStore_getNumFragments());
  fprintf(stderr, "  lambda   %f\n",      lambda * avgFragLength / windowSize);
  fprintf(stderr, "  length   %f\n",      windowSize * numFrag / lambda);
  fprintf(stderr, "\n");

  lambda = computeTruncatedMean(overlapLimit,
                                lambda,
                                into,
                                avgFragLength,
                                minOvl,
                                windowSize);

  fprintf(stderr, "Final Estimate (adjusted)\n");
  fprintf(stderr, "  nFrags   "F_S32"\n", gkpStore->gkStore_getNumFragments());
  fprintf(stderr, "  lambda   %f\n",      lambda * avgFragLength / windowSize);
  fprintf(stderr, "  length   %f\n",      windowSize * numFrag / lambda);

  delete [] factorial;
  delete    gkpStore;
  delete    ovlStore;
}

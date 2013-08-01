
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

#ifndef ChiSquareTest_CGW_H
#define ChiSquareTest_CGW_H

static const char *rcsid_ChiSquareTest_CGW_H = "$Id$";

/* Data structures for clustering edges based on Chi Squared tests on
   combinations of edges */

typedef struct { // Stores the scores of Chi Squared tests on sets
  // of edges plus the resulting  statistics of the merged set
  double score; // Chi Squared score
  LengthT distance; // The statistics of the merged set
  unsigned int passed:1; // whether the score exceeded the defined threshold
}Chi2ResultT;

typedef struct { // Data structure to pass to Chi Squared Compute subroutine
  // to avoid recomputing the weighted mean and inverse variance
  LengthT distance; // The statistics of an input to the test
  double weightedMean;
  double inverseVariance;
}Chi2ComputeT;

typedef struct { // Stores the score of pairwise Chi Squared tests on pairs
  // of clusters plus the resulting cluster statistics if the pair is merged
  double score; // Chi Squared score
  LengthT distance; // The statistics of the merged cluster
  int rowIndex; // The smaller cluster index
  int colIndex; // The larger cluster index
  unsigned int passed:1; // whether the score exceeded the defined threshold
  unsigned int active:1; // whether both clusters are still active
}ClusterScoreChi2T;

typedef struct ClusterChi2{ // The statistics of a cluster
  ClusterScoreChi2T *pairwiseScores; // Pointer to pairwise Chi Squared scores
  // between this cluster and other clusters
  struct ClusterChi2 *replacedBy; // When clusters are merged this points to
  // the new cluster - if == NULL then this is an active cluster
  struct ClusterChi2 *replaced; // When clusters are merged this points to
  // the old cluster
  int numInCluster; // Number in cluster
  LengthT distance; // The statistics of a cluster
}ClusterChi2T;

static int PairwiseChiSquare(double mean1, double variance1,
                             double mean2, double variance2,
                             LengthT *distance,
                             double *chiSquaredValue,
                             double chiSquaredThreshold){

  //  Negative variance in GOS III was causing significant problems.  These used to be asserts, now
  //  we just fail the test.
  //
  if ((variance1 <= 0.0) ||
      (variance2 <= 0.0))
    return(FALSE);

  assert(variance1 > 0.0);
  assert(variance2 > 0.0);

  //  Kludge because some overlaps are conflicting with each other and we need to figure out why!!!
  if(variance1 < 25.0)
    variance1 = 25.0;
  if(variance2 < 25.0)
    variance2 = 25.0;

  variance1 = 1.0 / variance1;
  variance2 = 1.0 / variance2;

  double compVariance = 1.0 / (variance1 + variance2);
  double compMean     = ((mean1 * variance1) + (mean2 * variance2)) * compVariance;

  double chiSquared  = ((mean1 - compMean) * (mean1 - compMean) * variance1 +
                       (mean2 - compMean) * (mean2 - compMean) * variance2);

  *chiSquaredValue = chiSquared;

  if (distance != NULL) {
    distance->mean     = compMean;
    distance->variance = compVariance;
    assert(compVariance >= 0.0);
  }

  return(chiSquared < chiSquaredThreshold);
}

#if 1

//  These were the historical defaults (= 3 * 7.879).
#define PAIRWISECHI2THRESHOLD_500  0.455  //   50.0%
#define PAIRWISECHI2THRESHOLD_010 23.637  //    ?.?% of values higher than this; this was the default for all cases
#define PAIRWISECHI2THRESHOLD_001 23.637  //    ?.?% of values higher than this; this was the default for all cases

#else

#define PAIRWISECHI2THRESHOLD_995  0.000  //   99.5% chance these are the same
#define PAIRWISECHI2THRESHOLD_975  0.000  //   97.5%
#define PAIRWISECHI2THRESHOLD_900  0.016  //   90.0%
#define PAIRWISECHI2THRESHOLD_500  0.455  //   50.0%
#define PAIRWISECHI2THRESHOLD_100  2.706  //   10.0%
#define PAIRWISECHI2THRESHOLD_050  3.841  //    5.0%
#define PAIRWISECHI2THRESHOLD_025  5.024  //    2.5%
#define PAIRWISECHI2THRESHOLD_010  6.635  //    1.0%
#define PAIRWISECHI2THRESHOLD_005  7.879  //    0.5%
#define PAIRWISECHI2THRESHOLD_001 10.828  //    0.1% chance these are the same

#endif


#endif

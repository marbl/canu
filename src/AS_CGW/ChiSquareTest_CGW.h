
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
/* 	$Id: ChiSquareTest_CGW.h,v 1.1.1.1 2004-04-14 13:50:16 catmandew Exp $	 */
#ifndef ChiSquareTest_CGW_H
#define ChiSquareTest_CGW_H

/* Data structures for clustering edges based on Chi Squared tests on
   combinations of edges */

typedef struct { // Stores the scores of Chi Squared tests on sets
  // of edges plus the resulting  statistics of the merged set
  float score; // Chi Squared score
  LengthT distance; // The statistics of the merged set
  unsigned int passed:1; // whether the score exceeded the defined threshold
}Chi2ResultT;

typedef struct { // Data structure to pass to Chi Squared Compute subroutine
  // to avoid recomputing the weighted mean and inverse variance
  LengthT distance; // The statistics of an input to the test
  float weightedMean;
  float inverseVariance;
}Chi2ComputeT;

typedef struct { // Stores the score of pairwise Chi Squared tests on pairs
  // of clusters plus the resulting cluster statistics if the pair is merged
  float score; // Chi Squared score
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

static int PairwiseChiSquare(float mean1, float variance1, float mean2,
                             float variance2, LengthT *distance,
                             float *chiSquaredValue,
                             float chiSquaredThreshold){

  float compMean, compVariance;
  float chiSquared, chiTemp;

  assert(variance1 > 0.0);
  assert(variance2 > 0.0);
  /* Kludge because some overlaps are conflicting with each other
     and we need to figure out why!!! */
  if(variance1 < 25.0){
    variance1 = 25.0;
  }
  if(variance2 < 25.0){
    variance2 = 25.0;
  }
  compVariance = 1.0 / ((1.0 / variance1) + (1.0 / variance2));
  compMean = ((mean1 / variance1) + (mean2 / variance2)) * compVariance;
  chiTemp = mean1 - compMean;
  chiTemp *= chiTemp;
  chiSquared = chiTemp / variance1;
  chiTemp = mean2 - compMean;
  chiTemp *= chiTemp;
  chiSquared += chiTemp / variance2;
  *chiSquaredValue = chiSquared;
  if(distance != NULL){
    distance->mean = compMean;
    distance->variance = compVariance;
    assert(compVariance >= 0.0);
  }
  return(chiSquared < chiSquaredThreshold);
}

#define PAIRWISECHI2THRESHOLD_CGW 3.0 * 7.879 

#endif

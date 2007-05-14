
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
/*************************************************************************
 Module:  AS_UTL_SequenceBucket
 Description:
 These objects support collecting stats on the frequency of n-mers in 
 a collection of sequence.

 SequenceBucketT:
   An array of buckets to collect stats on n-mer frequency.

 SequenceBucketArrayT:
   An array of SequenceBucketT to collect stats on m-mer frequency
   for m = 1,n.

 Assumptions:
      None.
 Document:

 *************************************************************************/
#ifndef AS_UTL_SEQUENCE_BUCKET_H
#define AS_UTL_SEQUENCE_BUCKET_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "AS_global.h"

typedef struct{
  int bucketWidth; // Number of characters at the beginning of the string to use for bucketizing
  int32 numBuckets; // 4^bucketWidth
  int32 numSamples;
  int32 *buckets;
  float   *trate; // theoretical probability
  float   *arate; // actual probability
}SequenceBucketT;


/* Create a SequenceBucketT for collecting frequency data on n-mers where n = bucketWidth */
SequenceBucketT *CreateSequenceBucket(int bucketWidth);

/* Delete a SequenceBucketT */
void DeleteSequenceBucket(SequenceBucketT *sequenceBucket);

/* Given a bucket number, compute the n-mer the corresponds to the bucket */
void ComputeBucketString(SequenceBucketT *sequenceBucket, int bucketNum, char *buffer);

/* Given an n-mer, compute the corresponding bucket number */
int32 ComputeBucket(int bucketWidth, char *actgSequence);

/* Increment the appropriate bucket using the bucketWidth chars at the beginning
   of actgSequence */
void IncrementBucketPrefix(SequenceBucketT *sequenceBucket, char *actgSequence);

/* Increment the appropriate bucket for each of length - bucketWidth substrings
   over the entire length of actgSequence */
void IncrementBucketTotal(SequenceBucketT *sequenceBucket, char *actgSequence, int length);

/* Compute the actual n-mer frequencies based on the counts in the buckets */
void ComputeBucketActualRates(SequenceBucketT *sequenceBucket);

/* Compute the theoretical frequency of the n-mer in bucket bucketNum given  the
   assumptions:
     1) The probability of each char of the n-mer is independent
     2) The probability of an a,c,t,g is given by the values in singleCharProbs

   This computation assumes that the probability of seeing a given n-mer is:
   P(n-mer) = P(char 1 of n-mer)*.....*P(char n of n-mer)

   The mean and variance of the number of instances of this n-mer seen are according
   to the binomial distribution.
*/
void ComputeBucketProbability(SequenceBucketT *sequenceBucket, int bucketNum, float   *singleCharProbs);

/* Compute the theoretical probabilities for all buckets */
void ComputeBucketTheoreticalRates(SequenceBucketT *sequenceBucket, float   *singleCharProbs);


/**************************************************************/

typedef struct{
  int numSequenceBuckets;
  int numSamples;
  SequenceBucketT **sequenceBuckets;
}SequenceBucketArrayT;




/* Create an array of sequenceBucketTs, one for each width from 1-maxWidth */
SequenceBucketArrayT *CreateSequenceBucketArray(int32 maxWidth);

/* Delete an array of sequenceBucketTs */
void DeleteSequenceBucketArray(SequenceBucketArrayT *sba);

/* Increment the maxWidth buckets in the SequenceBucket based on the prefix to actgSequence */
void IncrementSequenceBucketArrayPrefix(SequenceBucketArrayT *sequenceBucketArray, char *actgSequence);

/* Look for statistical anomolies in all of the SequenceBuckets */
void CheckSequenceBucketArray(SequenceBucketArrayT *sequenceBucketArray, float   *actgProbabilities, float   num_sigma,
			      FILE *fout, char *label);

/* Look for statistical anomolies in all of the SequenceBuckets.  bucket array 2 is the control */
void CheckSequenceBucketArraySanity(SequenceBucketArrayT *sequenceBucketArray1, SequenceBucketArrayT *sequenceBucketArray2,
				    float   *actgProbabilities, float   num_sigma, FILE *fout, char *label);



/*** SequenceLengthHistogramT ***/

typedef struct{
  int bucketWidth;
  int32 numBuckets; // 4^bucketWidth
  FILE **bucketHistoFiles;
  FILE **bucketFragIDFiles;
  char name[256];
  int32 numActiveBuckets;
}SequenceLengthHistogramT;


void ActivateSequenceLengthHistogram(SequenceLengthHistogramT *sequenceLengthHistogram, char *sequence);

SequenceLengthHistogramT *CreateSequenceLengthHistogram(int bucketWidth, char *name);

void DeleteSequenceLengthHistogram(SequenceLengthHistogramT *sequenceLengthHistogram);

void IncrementSequenceLengthHistogram(SequenceLengthHistogramT *sequenceLengthHistogram, char *seq,  int32 clearRangeLength, uint64 fragID);


#endif

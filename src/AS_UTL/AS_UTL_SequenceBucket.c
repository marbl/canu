
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
static char CM_ID[] = "$Id: AS_UTL_SequenceBucket.c,v 1.2 2004-09-23 20:25:29 mcschatz Exp $";
/*************************************************************************
 Module:  AS_UTL_SequenceBucket
 Description:
 These objects support collecting stats on the frequency of n-mers in 
 a collection of sequence.
 See AS_GKP_checkFrag.c for usage examples.
 See AS_UTL_SequenceBucket.h for some descriptions.
 Assumptions:
      None.
 Document:

 *************************************************************************/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "AS_UTL_SequenceBucket.h"

#define BAD_SEQ_BUCKET_INDEX  -1

SequenceBucketT *CreateSequenceBucket(int bucketWidth){
  SequenceBucketT *sequenceBucket = (SequenceBucketT *)
    calloc(sizeof(SequenceBucketT),1);
  int i;

  sequenceBucket->bucketWidth = bucketWidth;
  sequenceBucket->numBuckets = 1<<(2 * bucketWidth);
#if 0
  fprintf(stderr,"* Creating SequenceBucket of width %d with %d buckets\n",
	  bucketWidth, sequenceBucket->numBuckets);
#endif
  sequenceBucket->buckets = (int32 *)
    calloc(sequenceBucket->numBuckets, sizeof(int32));
  sequenceBucket->trate = (float32 *)
    malloc(sequenceBucket->numBuckets * sizeof(float32));
  sequenceBucket->arate = (float32 *)
    malloc(sequenceBucket->numBuckets *sizeof(float32));
  for(i = 0; i < sequenceBucket->numBuckets; i++){
    sequenceBucket->trate[i] = 0.0;
    sequenceBucket->arate[i] = 0.0;
  }
  return sequenceBucket;
}
void DeleteSequenceBucket(SequenceBucketT *sequenceBucket){

  free(sequenceBucket->buckets);
  free(sequenceBucket->arate);
  free(sequenceBucket->trate);
  free(sequenceBucket);
}

void ComputeBucketString(SequenceBucketT *sequenceBucket, int bucketNum, char *buffer){
  char *cursor = buffer + sequenceBucket->bucketWidth ;
  int i;

  *cursor-- = '\0';
  for(i = 0; i <  sequenceBucket->bucketWidth; i++){
      char c = 'n';
    switch(bucketNum & 0x3){
    case 0:
      c = 'a';
      break;
    case 1:
      c = 'c';
      break;
    case 2:
      c = 'g';
      break;
    case 3:
      c = 't';
      break;
    }

    bucketNum>>= 2;
    *cursor-- = c;
  }

}



int32 ComputeBucket(int bucketWidth, char *actgSequence){
  int i;
  int32 bucketNum = 0;
  for(i = 0; i < bucketWidth;){
    int32 val ;
    switch(tolower(actgSequence[i])){
    case 'a':
      val = 0;
      break;
    case 'c':
      val = 1;
      break;
    case 'g':
      val = 2;
      break;
    case 't':
      val = 3;
      break;
    default:
      fprintf(stderr,"*At position %d encountered unexpected char %c\n", i, actgSequence[i]);
      return BAD_SEQ_BUCKET_INDEX;
    }
    bucketNum += val;

    // Don't shift on the last iteration
    if(++i < bucketWidth){
      bucketNum <<= 2;
    }
  }
  //  fprintf(stderr,"* seq:%s BucketWidth %d BucketNum = %d\n", actgSequence,bucketWidth, bucketNum);
  return bucketNum;
}

/* Increment the appropriate bucket using the bucketWidth chars at the beginning
   of actgSequence */
void IncrementBucketPrefix(SequenceBucketT *sequenceBucket, char *actgSequence){
  int32 bucketNum = ComputeBucket(sequenceBucket->bucketWidth, actgSequence);
  if(bucketNum != BAD_SEQ_BUCKET_INDEX)
  {
    sequenceBucket->numSamples++;
    sequenceBucket->buckets[bucketNum]++;
  }
}

/* Increment the appropriate bucket for each of length - bucketWidth substrings
   over the entire length of actgSequence */
void IncrementBucketTotal(SequenceBucketT *sequenceBucket, char *actgSequence, int length){
  int i;
  int max = length - sequenceBucket->bucketWidth + 1;
  for(i = 0; i < max; i++)
  {
    int32 bucketNum = ComputeBucket(sequenceBucket->bucketWidth, actgSequence + i);
    if(bucketNum != BAD_SEQ_BUCKET_INDEX)
    {
      sequenceBucket->numSamples++;
      sequenceBucket->buckets[bucketNum]++;
    }
  }
}
void ComputeBucketActualRates(SequenceBucketT *sequenceBucket){
  int i;
  for(i = 0; i < sequenceBucket->numBuckets; i++){
    if(sequenceBucket->numSamples > 0)
      sequenceBucket->arate[i] = (float32)(sequenceBucket->buckets[i])/(float32)(sequenceBucket->numSamples);
    else
      sequenceBucket->arate[i] = 0.0;
  }
}

void CheckBucketRates(SequenceBucketT *sequenceBucket, float32 nsigma, FILE *fout, char *label){
  int i;
  float32 delta, sigma;
  char buffer[256];
  for(i = 0; i < sequenceBucket->numBuckets; i++){
    // sigma = n * p (1-p)
    sigma = sequenceBucket->numSamples * (sequenceBucket->trate[i]) * (1.0 - sequenceBucket->trate[i]);
    delta = sequenceBucket->buckets[i] - sequenceBucket->numSamples * sequenceBucket->trate[i];
    ComputeBucketString(sequenceBucket, i, buffer);
    if(sequenceBucket->buckets[i] > 10 &&
       delta > nsigma * sigma){
      fprintf(fout,"%s* %1d-mer %-12s is overrepresented (%d/%d) theory = %8g  actual = %d sigma = %8g\n",
	      label,
	      sequenceBucket->bucketWidth, buffer, 
	      sequenceBucket->buckets[i], sequenceBucket->numSamples,
	      sequenceBucket->numSamples * sequenceBucket->trate[i], 
	      sequenceBucket->buckets[i],

	      sigma);
    }
#if 0
else{

      fprintf(fout,"%s* %d-mer %s is OK theory = %g  actual = %g sigma = %g\n",
	      label,
	      sequenceBucket->bucketWidth, buffer, 
	      sequenceBucket->trate[i],
	      sequenceBucket->arate[i],
	      sigma);
    }
#endif
  }
}

void CheckBucketRatesSanity(SequenceBucketT *sequenceBucket,
			    SequenceBucketT *sequenceBucket_sanity,
			    float32 nsigma, FILE *fout, char *label){
  int i;
  float32 delta, sigma;
  char buffer[256];
  for(i = 0; i < sequenceBucket->numBuckets; i++){
    // sigma = n * p (1-p)
    sigma = sequenceBucket->numSamples * (sequenceBucket->trate[i]) * (1.0 - sequenceBucket->trate[i]);
    delta = sequenceBucket->buckets[i] - sequenceBucket->numSamples * sequenceBucket->trate[i];
    ComputeBucketString(sequenceBucket, i, buffer);
    if(sequenceBucket->buckets[i] > 10 &&
       delta > nsigma * sigma){
      fprintf(fout,"%s* %1d-mer %-12s is overrepresented (%d/%d) theory = %8g  actual = %d sigma = %8g (sanity has %d)\n",
	      label,
	      sequenceBucket->bucketWidth, buffer, 
	      sequenceBucket->buckets[i], sequenceBucket->numSamples,
	      sequenceBucket->numSamples * sequenceBucket->trate[i], 
	      sequenceBucket->buckets[i],
	      sigma,
	      sequenceBucket_sanity->buckets[i]);
    }
#if 0
else{

      fprintf(fout,"%s* %d-mer %s is OK theory = %g  actual = %g sigma = %g\n",
	      label,
	      sequenceBucket->bucketWidth, buffer, 
	      sequenceBucket->trate[i],
	      sequenceBucket->arate[i],
	      sigma);
    }
#endif
  }
}


void ComputeBucketProbability(SequenceBucketT *sequenceBucket, int bucketNum, float32 *singleCharProbs){
  int i;
  float32 prob = 1.0;
  int32 scratch = bucketNum;
  for(i = 0; i <  sequenceBucket->bucketWidth; i++){
      int index = scratch & 0x3;
      
      assert(index >= 0 && index <= 3);
      prob *= singleCharProbs[index];

      scratch>>= 2;
  }

  sequenceBucket->trate[bucketNum] = prob;
}

void ComputeBucketTheoreticalRates(SequenceBucketT *sequenceBucket, float32 *singleCharProbs){
  int i;
  for(i = 0; i < sequenceBucket->numBuckets; i++){
    if(sequenceBucket->numSamples > 0)
      ComputeBucketProbability(sequenceBucket, i, singleCharProbs);
    else
      sequenceBucket->trate[i] = 0.0;
  }
}



/**************************************************************/

SequenceBucketArrayT *CreateSequenceBucketArray(int32 maxWidth){
  int i;
  SequenceBucketArrayT *sba = (SequenceBucketArrayT *)
    calloc(sizeof(SequenceBucketArrayT),1);

  sba->numSamples = 0;
  sba->numSequenceBuckets = maxWidth;
  sba->sequenceBuckets = (SequenceBucketT **)
    calloc(sizeof(SequenceBucketT *),maxWidth);
  for(i = 0; i < maxWidth;i++){
    sba->sequenceBuckets[i] = CreateSequenceBucket(i+1);
  }

  return sba;
}

void DeleteSequenceBucketArray(SequenceBucketArrayT *sba){
  int i;

  for(i = 0; i < sba->numSequenceBuckets;i++){
    DeleteSequenceBucket(sba->sequenceBuckets[i]);
  }

  free(sba->sequenceBuckets);
  free(sba);
}


void IncrementSequenceBucketArrayPrefix(SequenceBucketArrayT *sequenceBucketArray, char *actgSequence){
  int i;
  sequenceBucketArray->numSamples++;
  for(i = 0; i < sequenceBucketArray->numSequenceBuckets; i++){
    IncrementBucketPrefix(sequenceBucketArray->sequenceBuckets[i], actgSequence);
  }
}

void CheckSequenceBucketArray(SequenceBucketArrayT *sequenceBucketArray, float32 *actgProbabilities, float32 num_sigma,
			      FILE *fout, char *label){
  int i;
  for(i = 0; i < sequenceBucketArray->numSequenceBuckets; i++){
    ComputeBucketTheoreticalRates(sequenceBucketArray->sequenceBuckets[i], actgProbabilities);
    ComputeBucketActualRates(sequenceBucketArray->sequenceBuckets[i]);
    CheckBucketRates(sequenceBucketArray->sequenceBuckets[i], num_sigma, fout, label);
  }
}

void CheckSequenceBucketArraySanity(SequenceBucketArrayT *sequenceBucketArray1,
				    SequenceBucketArrayT *sequenceBucketArray2,
				    float32 *actgProbabilities, float32 num_sigma,
				    FILE *fout, char *label){
  int i;
  for(i = 0; i < sequenceBucketArray1->numSequenceBuckets; i++){
    ComputeBucketTheoreticalRates(sequenceBucketArray1->sequenceBuckets[i], actgProbabilities);
    ComputeBucketActualRates(sequenceBucketArray1->sequenceBuckets[i]);
    CheckBucketRatesSanity(sequenceBucketArray1->sequenceBuckets[i], 
			   sequenceBucketArray2->sequenceBuckets[i], 
			   num_sigma, fout, label);
  }
}






/*** SequenceLengthHistogramT ***/

void ActivateSequenceLengthHistogram(SequenceLengthHistogramT *sLH, char *sequence){
  char buffer[2048];
  int32 bucketNum = ComputeBucket(sLH->bucketWidth, sequence);
  FILE *fp = NULL;
  if(bucketNum == BAD_SEQ_BUCKET_INDEX){
    fprintf(stderr,"* ActiveSequenceLengthHistogram for %s ... bad sequence bucket index\n",
	    sLH->name);
    return;
  }
  if(sLH->bucketHistoFiles[bucketNum] != NULL){
    fprintf(stderr,"* ActiveSequenceLengthHistogram for %s ... bucket %d already active\n",
	    sLH->name, bucketNum);
    return;
  }
  sprintf(buffer,"%*s_%s.cgm", sLH->bucketWidth, sequence,sLH->name);
  fp = sLH->bucketHistoFiles[bucketNum] = fopen(buffer,"w");
  assert(NULL != fp);
  fprintf(fp,"Clear Range Lengths for %s %*s\n",sLH->name,sLH->bucketWidth, sequence);


  sprintf(buffer,"%*s_%s.frags", sLH->bucketWidth, sequence,sLH->name);
  fp = sLH->bucketFragIDFiles[bucketNum] = fopen(buffer,"w");
  assert(NULL != fp);

  sLH->numActiveBuckets++;

}

SequenceLengthHistogramT *CreateSequenceLengthHistogram(int bucketWidth, char *name){
  SequenceLengthHistogramT *sequenceLengthHistogram = (SequenceLengthHistogramT *)
    calloc(sizeof(SequenceLengthHistogramT),1);

  strncpy(sequenceLengthHistogram->name, name, 256);
  sequenceLengthHistogram->bucketWidth = bucketWidth;
  sequenceLengthHistogram->numBuckets = 1<<(2 * bucketWidth);
#if 0
  fprintf(stderr,"* Creating SequenceLengthHistogram of width %d with %d buckets\n",
	  bucketWidth, sequenceBucket->numBuckets);
#endif
  sequenceLengthHistogram->bucketHistoFiles = (FILE **)
    calloc(sequenceLengthHistogram->numBuckets, sizeof(FILE *));
  sequenceLengthHistogram->bucketFragIDFiles = (FILE **)
    calloc(sequenceLengthHistogram->numBuckets, sizeof(FILE *));

  return sequenceLengthHistogram;
}


void DeleteSequenceLengthHistogram(SequenceLengthHistogramT *sLH){
  int i;
  for(i = 0; i < sLH->numBuckets; i++){
    if(sLH->bucketHistoFiles[i])
    fclose(sLH->bucketHistoFiles[i]);
  }
  free(sLH->bucketHistoFiles);
}

void IncrementSequenceLengthHistogram(SequenceLengthHistogramT *sLH, char *seq, int32 clearRangeLength, uint64 fragID){
  int bucketNum = ComputeBucket(sLH->bucketWidth, seq);
  FILE *fp = NULL;

  if(bucketNum == BAD_SEQ_BUCKET_INDEX)
    return;
  
  fp = sLH->bucketHistoFiles[bucketNum];
  //  fprintf(stderr,"* Sequence %s has bucket %d and is %s active\n",
  //	  seq, bucketNum, (!fp?"NOT":""));
  if(fp){

    fprintf(fp,"%d\n", clearRangeLength);
    fflush(fp);


    fp = sLH->bucketFragIDFiles[bucketNum];
    assert(NULL != fp);
    fprintf(fp,F_UID "\n", fragID);
    fflush(fp);
  }
}

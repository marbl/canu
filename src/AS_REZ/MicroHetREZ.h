
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
/*********************************************************************
   CVS_ID:  $Id: MicroHetREZ.h,v 1.1.1.1 2004-04-14 13:53:23 catmandew Exp $
 *********************************************************************/
#ifndef MICROHETREZ_H
#define MICROHETREZ_H

//WHETHER OR NOT TO HAVE double *pval IN DISCRIMINATOR FUNCTIONS
#define RETURNPVALS

#include "AS_global.h"
#include "MicroHetPartitionsREZ.h"

//#include "MicroHetScoreREZ.h"
//#include "MicroHetPartitionsREZ.h"
//#include "MicroHetInterfaceREZ.h"

#include "UtilsREZ.h"
#include "AS_UTL_skiplist.h"
#include "AS_UTL_Var.h"
#include "PrimitiveVA_MSG.h"


#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_distStore.h"
#include "Array_CNS.h"

/* The order of the enum is important  */
/* UNITIG_IS_UNKOWN is greater than UNITIG_IS_SIMPLE */

typedef enum 
{
  UNITIG_IS_SHALLOW,
  UNITIG_IS_SIMPLE, 
  UNITIG_IS_UNKNOWN, 
  UNITIG_IS_REPETITIVE
} UnitigStatus_t;


/* This typdef describes a subalignment and its status */
typedef struct
{
  int start;
  int end;
  UnitigStatus_t simple;
  Partition_t *p;
  Marker_t    *m;
} TestSegment_t;


/* This is the maximal number for which the x-contributing attributes
   are stored in an Alignment_T */

#define MAX_X 4

/* The alignment structure. 
   The array in countA can be changed using col_count 
   with a specific marker */

typedef struct{
  char** ali;  // the actual alignment array
  int* countA; // an array containing the number of 'A' in each column
  int* countC; // an array containing the number of 'C' in each column
  int* countG; // an array containing the number of 'G' in each column
  int* countT; // an array containing the number of 'T' in each column
  int* countDash; // an array containing the number of '-' in each column
  int* countBlank; // an array containing the number of ' ' in each column
  int cols;
  int rows;
  int noOfSegs;
  double **seqErrArray; 
  int hasQuality;
} Alignment_t;


//AARON'S MODIFICATIONS
#define maxfive(v,w,x,y,z) max(v,max(max(w,x),max(y,z)))
//Slower, more accurate MP statistics
#define EXACT_EXPECTED_SAVEDSTEPS
typedef struct mpstat {
  int Obs;   /* Number of parsimony-style "steps" which could be saved by
                   the optimal partitioning of a column of data as compared to
		   the assumption that each mismatch is independent.
		   Allows each column to choose its optimal partition and so
		   is not guaranteed to be obtainable by a single tree for the
		   whole alignment */
  double Exp;   /* Expected number of columns with a pair of non-consensus 
		   matching characters; approximates the expected number that
		   corresponds to Obs */
  double pr;    /* Making a Poisson approximation, the probability of seeing 
		   Obs when expecting Exp */
} MPSTAT;

double ExpectedSavedSteps[200];

// END AARON

/* functions to allocate, free and print an alignment */
void        print_alignment(Alignment_t *p, int w);
Alignment_t *allocate_alignment(int c, int r);
void        free_alignment(Alignment_t *p);

/* functions to allocate, free and print an testSegment  */
TestSegment_t *allocate_testSegment(int l);
void        free_testSegment(TestSegment_t *t);


// The minimum length for which we test for microhet.
#define MIN_TEST_LENGTH_REZ 100

// The minimum length for which we test for Aarons MP score
#define MIN_MPTEST_LENGTH_REZ 100


// The minimum critical value
#define MIN_THRESHOLD_REZ 2
// We assume that an alignment with more rows than that is repetitive
#define TEST_UPPER_BOUND 40


UnitigStatus_t test_simple(Alignment_t *ali, double thresh, Marker_t* m,
			   int start, int end, int* critical
#ifdef RETURNPVALS
,double *pval
#endif
);

UnitigStatus_t test_MPsimple(Alignment_t *ali, double thresh, Marker_t* m, 
			     int start, int end
#ifdef RETURNPVALS
,double *pval
#endif
);

UnitigStatus_t test_PWsimple(Alignment_t *ali, double thresh, Marker_t* m, 
			     int start, int end
#ifdef RETURNPVALS
, double *pval
#endif
);


UnitigStatus_t is_IUM_simple(IntUnitigMesg* ium, FragStoreHandle handle,
			     Alignment_t **ali, double thresh, int variant
#ifdef RETURNPVALS
, double *pval
#endif
);

UnitigStatus_t is_IUM_MPsimple(IntUnitigMesg* ium, FragStoreHandle handle,
			     Alignment_t **ali, double thresh, int variant
#ifdef RETURNPVALS
, double *pval
#endif
);

UnitigStatus_t is_IUM_PWsimple(IntUnitigMesg* ium, FragStoreHandle handle,
			     Alignment_t **ali, double thresh, int variant
#ifdef RETURNPVALS
, double *pval
#endif
);

MPSTAT MP_score_alignment(Alignment_t *alignment,double erate, int s, int e);

void bipartition(Alignment_t* a, Marker_t* m, Partition_t* p, 
		 int start, int end, double alpha, int group);

void count_columns(Alignment_t* a, Marker_t* m);
TestSegment_t **get_segments(Alignment_t *ali, int *noOfSegs);

#endif









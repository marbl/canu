
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
   CVS_ID:  $Id: MicroHetREZ_test3.h,v 1.1.1.1 2004-04-14 13:53:22 catmandew Exp $
 *********************************************************************/
#ifndef AS_REZ_MICROHETREZ_H
#define AS_REZ_MICROHETREZ_H

#include "MicroHetPartitionsREZ_test3.h"

#include "UtilsREZ.h"
#include "AS_UTL_skiplist.h"
#include "AS_UTL_Var.h"
#include "PrimitiveVA_MSG.h"


#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_fragStorePartition.h"
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

// END AARON

/* functions to allocate, free and print an alignment */
void        AS_REZ_print_alignment(Alignment_t *p, int w);
Alignment_t *AS_REZ_allocate_alignment(int c, int r);
void        AS_REZ_free_alignment(Alignment_t *p);

// The minimum length for which we test for Aarons MP score
#define MIN_MPTEST_LENGTH_REZ 50


// We assume that an alignment with more rows than that is repetitive
#define TEST_UPPER_BOUND 40


UnitigStatus_t AS_REZ_test_MPsimple(Alignment_t *ali, double thresh, Marker_t* m, 
			     int start, int end,double *pval);

UnitigStatus_t AS_REZ_is_IUM_MPsimple(IntUnitigMesg* ium, FragStoreHandle handle, tFragStorePartition *phandle,
			     Alignment_t **ali, double thresh, int variant, double *pval);


MPSTAT AS_REZ_MP_score_alignment(Alignment_t *alignment,double erate, int s, int e);

void AS_REZ_count_columns(Alignment_t* a, Marker_t* m);

#endif

//Functions for Poisson upper tail probability
//  ... with calculation of lambda^k and k! by looping
double AS_REZ_Poisson(int Obs,double Exp);
//  ... with calculation of lambda^k and k! with exp(), log() and lgamma()
double AS_REZ_Poisson_prob(int Obs,double Exp);




/* AS_REZ_MP_MicroHet_prob() 

   RESULT: The function returns a (double) pvalue (probability) of an
   input unitig being SIMPLE -- meaning, having mismatches due to randomly 
   distributed sequencing errors.

   If the returned value is sufficiently small, the unitig should be treated
   as a likely repeat.

   A return value of 1.0 may indicate that the unitig was not deep enough for
   a meaningful test.

   Some false positives may be induced by polymorphisms; however, the 
   calculation should not be drastically misled by multibase indel 
   polymorphisms.

   INPUT:

   bqarray : an array of size [depth*2]*len of bases and quality values
             in alternative rows, giving a multialignment
   idarray : an array of size depth*len giving the fragment iid of each base
             in the multialignment
   handle  : the fragStore from which locale information for each fragment iid
             will be obtained  (-1 (NULLFRAGSTOREHANDLE) if paritioned store is used.)
   phandle  : the partitioned fragStore from which locale information for each fragment iid
             will be obtained (NULL if a traditional unpartitioned store is used);
   len     : number of columns in the multialignment
   depth   : number of rows in the multialignment
*/

double AS_REZ_MP_MicroHet_prob(char **bqarray,int **idarray,FragStoreHandle handle,tFragStorePartition *phandle,int len,int depth);


/* this is the main test function for a unitig. 
   It returns a pvalue (roughly, a probability) that the unitig is simple
   (mismatches are random errors).
   A return value of 1.0 may indicate that the unitig was not deep enough for meaningful test.
*/

double AS_REZ_prob_IUM_MPsimple(IntUnitigMesg* ium, FragStoreHandle handle, tFragStorePartition *phandle);




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
static char CM_ID[] = "$Id: MicroHetREZ.c,v 1.3 2005-03-22 19:07:44 jason_miller Exp $";

#include "MicroHetREZ.h"
#include "MicroHetScoreREZ.h"
#include "MicroHetPartitionsREZ.h"
#include "MicroHetInterfaceREZ.h"
#include "UtilsREZ.h"
#include "AS_UTL_skiplist.h"

#include "PrimitiveVA_MSG.h"

SL_DEF(Partition_t)

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
  
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_distStore.h"

//Stuff to control diagnostic output.
#define DEBUG 0
#define STATPRINT 0
static int doPrintMPobs=0;
int doPrintPWchisq=0;
 
SL_TYPE(Partition_t) *guess_partition(Alignment_t *a, int ct);

/* memory allocation and deallocation for TestSegment_t */
TestSegment_t *allocate_testSegment(int l)
{
  TestSegment_t* t = (TestSegment_t*) safe_malloc(sizeof(TestSegment_t));
  t->p             = allocate_partition(l);
  t->m             = allocate_marker(l);
  t->start         = -1;
  t->end           = -1;
  return t; 
}

void free_testSegment(TestSegment_t* t)
{
  free_partition(t->p);
  free_marker(t->m);
  free(t);
}

/* memory allocation and deallocation for Alignment_t */
Alignment_t *allocate_alignment(int c, int r)
{
  int i,j;
  Alignment_t *a = (Alignment_t*) safe_malloc(sizeof(Alignment_t));

  if( c > 0 && r > 0){
    a->cols = c;
    a->rows = r; 
    a->ali        = (char**) safe_calloc(sizeof(char*),c);
    a->countA     = (int*) safe_calloc(sizeof(int),c);
    a->countC     = (int*) safe_calloc(sizeof(int),c);
    a->countG     = (int*) safe_calloc(sizeof(int),c);
    a->countT     = (int*) safe_calloc(sizeof(int),c);
    a->countDash  = (int*) safe_calloc(sizeof(int),c);
    a->countBlank = (int*) safe_calloc(sizeof(int),c);
    for(i=0; i<c; i++)
      a->ali[i] = (char*) safe_calloc(sizeof(char),r);
    a->seqErrArray  = (double**) safe_calloc(sizeof(double*),c);
    for(i=0; i<c; i++)
      a->seqErrArray[i] = (double*) safe_calloc(sizeof(double),r);
    a->hasQuality = FALSE;
  }
  else
    {
      a->cols = 0;
      a->rows = 0;
      a->ali        = NULL;
      a->countA     = NULL;
      a->countC     = NULL;
      a->countG     = NULL;
      a->countT     = NULL;
      a->countDash  = NULL;
      a->countBlank = NULL;
      a->seqErrArray  = NULL;
    }
  return a;
}




void free_alignment(Alignment_t* a)
{
  int i;

  if(a->ali != NULL){
    for(i=0; i<a->cols; i++){
      assert(a->ali[i] != NULL);  
    free(a->ali[i]);
    }
    free(a->ali);

    for(i=0; i<a->cols; i++){
      assert(a->seqErrArray[i] != NULL);
      free(a->seqErrArray[i]);
    }
    free(a->seqErrArray);
    free(a->countA);
    free(a->countC);
    free(a->countG);
    free(a->countT);
    free(a->countDash);
    free(a->countBlank);
  }
  free(a);
}


/* a print function for alignments */

void print_alignment(Alignment_t *a,  int w)
{
  int i,j,l;
  int c = a->cols;
  int r = a->rows;
  int iter = 0;
  char consensus[a->cols];
  int count[6]; // A C G T Dash N
    
  for(i=0; i<a->cols; i++)
    {
      int max,k;
      for(k=0; k<6; k++)
	count[k] = 0;
      for(j=0; j<a->rows; j++)
	{
	  switch(a->ali[i][j])
	    {
	    case 'A' :
	      count[0]++;
	      break;	
	    case 'C' :
	      count[1]++;
	      break;
	    case 'G' :
	      count[2]++;
	      break;
	    case 'T' :
	      count[3]++;
	      break;
	    case '-' :
	      count[4]++;
	      break;
	    case 'N' :
	      count[5]++;
	      break;
	    }	
	}	
      max = 0;
      for(k=0; k<6; k++){
	if( count[max] < count[k] )
	  max = k;
      
	switch( max )
	  {
	  case 0 :
	    consensus[i] = 'A';
	    break;	
	  case 1 :
	    consensus[i] = 'C';
	    break;
	  case 2 :
	    consensus[i] = 'G';
	    break;
	  case 3 :
	    consensus[i] = 'T';
 	    break;
	  case 4 :
	    consensus[i] = '-';
	    break;
	  case 5 :
	    consensus[i] = 'N';
	    break;
	  }	
      }
    }
  do{
    int p = w*iter+w;
    l = ( p < c ? p : c );
    
    for(i=iter*w; i<l; i++)
      printf("%c",consensus[i]);
    printf(" %d\n",l);
    
    for(i=0; i<r; i++){
      for(j=iter*w; j<l; j++)
        if( a->ali[j][i] != consensus[j] )
          printf("%c",a->ali[j][i]);
        else
          printf("%c",'.');
      printf(" %d\n",l);
    }

    iter++;
    printf("\n");
  }while(iter*w < c);
  
}


void count_columns(Alignment_t* a, Marker_t* m)
{
  register int i,j;
  for(i=0; i<a->cols; i++)
    {
      a->countA[i] = 0;
      a->countC[i] = 0;
      a->countG[i] = 0;
      a->countT[i] = 0;
      a->countDash[i]  = 0;
      a->countBlank[i] = 0;

      for(j=0; j<a->rows; j++)
	if( m->set[j] == TRUE)
	  switch(a->ali[i][j])
	    {
	    case 'A' :
	      a->countA[i]++;
	      break;	
	    case 'C' :
	      a->countC[i]++;
	      break;
	    case 'G' :
	      a->countG[i]++;
	      break;
	    case 'T' :
	      a->countT[i]++;
	      break;
	    case '-' :
	      a->countDash[i]++;
	      break;
	    case ' ' :
	      a->countBlank[i]++;
	    case 'N' : // NOTE we count Ns as blanks
	      a->countBlank[i]++;
	      break;
	    }	
    }
}

/* returns TRUE is the two columns in the alignment differ
   in their arrangments of blanks, which means that a fragment
   enters or leaves the alignment*/

int is_different(Alignment_t *ali, int c1, int c2)
{
  int r;
  if( c2 >= ali->cols )
    return TRUE;
  for(r=0; r<ali->rows; r++)
    if( (ali->ali[c1][r] == ' ' && ali->ali[c2][r] != ' ') || 
	(ali->ali[c1][r] != ' ' && ali->ali[c2][r] == ' ') )
      return TRUE;

  return FALSE;
}
 

/* divides an alignment into testable segments, that is into consecutive
   runs of columns in such that the same fragments occur in the rows of the
   segment */
TestSegment_t **get_segments(Alignment_t *ali, int *noOfSegs)
{
  int i;
  int c1,c2;
  TestSegment_t **t;
  
  c1 = 0;
  c2 = 1;
  *noOfSegs = 0;

  while( c2 < ali->cols )
    {
      while( is_different(ali,c1,c2++) == FALSE )
	;
      c1 = c2-1;
      (*noOfSegs)++;
    }

  t = (TestSegment_t**) safe_calloc(sizeof(TestSegment_t*),*noOfSegs);

  for(i=0; i<*noOfSegs; i++)
    t[i] = allocate_testSegment(ali->rows);

  c1 = 0;
  c2 = 1;  
  *noOfSegs = 0;

  while( c2 < ali->cols )
    {
      while( is_different(ali,c1,c2++) == FALSE )
	;
      t[*noOfSegs]->start   = c1;
      t[(*noOfSegs)++]->end = c2-1;

#if DEBUG > 5
      printf("segment %d [%d,%d] length = %d\n",*noOfSegs,c1,c2-1,c2-c1);
#endif

      c1 = c2-1;
    }
  return t;
}




/* the main test function */

UnitigStatus_t test_simple(Alignment_t *ali, double thresh, Marker_t* m, int start, int end, int* critical
#ifdef RETURNPVALS
,double *pval
#endif
){
  int i,j,s;
  int ret = UNITIG_IS_SIMPLE;
  TestSegment_t* segmentList;
  int             noOfSegs;
  double p;  
  double exponent;
  double pThresh; 
  // pThresh is the critical value with which we have to conduct
  // the (rows choose 2) fixed rows tests such that the overall error
  // guarantee is thresh    

  int rows;

#ifdef RETURNPVALS
  *pval=1.0; //Sets up default return for cases when no test can be made.
#endif

  rows=0;
  // we count the number of marked rows
  for(i=0; i<m->len; i++)
    if( m->set[i] == TRUE )
      rows++;

  // if we only have three rows, we cannot separate them
  if( rows <= 3)
    return UNITIG_IS_SHALLOW;

  exponent = 1.0/(rows*(rows-1)/2);
  pThresh  = 1.0-pow(1.0-thresh,exponent);

  /* If the number of rows is way to big, we assume that 
     the alignment is repetitive. The test is here too costly */

  if( rows > TEST_UPPER_BOUND )
    {
      *critical = 0;
      return UNITIG_IS_REPETITIVE;
    }

  if( end-start >= MIN_TEST_LENGTH_REZ )
    {
      double heurSeqErr = guess_seqErr(ali,m,start,end);
      p = column_prob_two_fixed(rows,heurSeqErr);

#if DEBUG > 0
      printf("GuessedSeqErr=%lf\n",heurSeqErr);
#endif
      if( p > 0.0)
	{
	  *critical = critical_binomial_value(p,pThresh,end-start);
#if DEBUG > 1
	  printf("p=%e , pThresh=%e, end-start=%d, critical=%d\n",p,pThresh,end-start,*critical);
#endif      
	  
	  count_columns(ali,m);

	  for(i=0; i<ali->rows-1; i++)
	    {
	      if( m->set[i] != TRUE )
		continue;
	      for(j=i+1; j<ali->rows; j++)
		{
		  int cont;
		  if( m->set[j] != TRUE )
		    continue;
		  
		  cont = no_col_contributing_two_fixed(ali,i,j,start,end);

#ifdef RETURNPVALS
		  *pval=min(*pval,binomial_tail_prob(end-start,cont,p));
#endif
		  if( cont >= *critical )
		    {
		      if( cont <= MIN_THRESHOLD_REZ )
			{
			  if( ret == UNITIG_IS_SIMPLE)
			    ret = UNITIG_IS_UNKNOWN;
#if DEBUG > 0
			  printf("pair %d,%d test yields UNKNOWN in segment [%d,%d] with value %d. CV = %d\n",i,j,start,end,cont,*critical);
#endif
			}
		      else
			{
			  ret = UNITIG_IS_REPETITIVE;
#if DEBUG > 0
			  printf("pair %d,%d test yields REPETITIVE in segment [%d,%d] with value %d. CV = %d\n",i,j,start,end,cont,*critical);
#endif
			}
		    }
		} // j loop
	    } // i loop
	}
      else
	ret = UNITIG_IS_UNKNOWN;  // p is 0.0, hence we cannot test we cannot
    }
  else
    ret = UNITIG_IS_UNKNOWN;	 // The width of the segment is too small

  return ret;  
}




UnitigStatus_t test_MPsimple(Alignment_t *ali, double thresh, Marker_t* m, 
			     int start, int end
#ifdef RETURNPVALS
,double *pval
#endif
){
  int i,j,s;
  int ret = UNITIG_IS_SIMPLE;
  TestSegment_t* segmentList;
  int             noOfSegs;
  int rows;

#ifdef RETURNPVALS
  *pval=1.0; //Sets up default return for cases when no test can be made.
#endif

  // we count the number of marked rows
  rows=0;
  for(i=0; i<m->len; i++)
    if( m->set[i] == TRUE )
      rows++;

  // if we only have three rows, we cannot separate them
  if( rows <= 3)
    return UNITIG_IS_SHALLOW;

  /* If the number of rows is way to big, we assume that 
     the alignment is repetitive. The test is here too costly */

  if( rows > TEST_UPPER_BOUND ){
    return UNITIG_IS_REPETITIVE;
  }

#if DEBUG > 0
      printf("test_MPsimple start=%d, end=%d\n",start,end);
#endif

  if( end-start >= MIN_MPTEST_LENGTH_REZ ){
    double heurSeqErr = guess_seqErr(ali,m,start,end);
    MPSTAT mpresult;
#if DEBUG > 0
    printf("GuessedSeqErr=%lf for start=%d, end=%d\n",heurSeqErr,start,end);
#endif
    mpresult=MP_score_alignment(ali,heurSeqErr,start,end);

    if(doPrintMPobs)
      printf(">>>>> Aaron's MP evaluation: O = %d E = %f pr = %e\n",mpresult.Obs, mpresult.Exp,mpresult.pr);

    if(mpresult.pr<=thresh){
#if DEBUG > 0
      printf(" : NOT SIMPLE\n\n");
#endif
      ret = UNITIG_IS_REPETITIVE;
    }
    else {
#if DEBUG > 0
      printf(" : SIMPLE\n\n");
#endif
      ret = UNITIG_IS_SIMPLE;
    }
    *pval=mpresult.pr;
  }

  return ret;  
}


UnitigStatus_t test_PWsimple(Alignment_t *ali, double thresh, Marker_t* m, 
			     int start, int end
#ifdef RETURNPVALS
,double *pval
#endif
){
  int i,j,s;
  int ret = UNITIG_IS_SIMPLE;
  TestSegment_t* segmentList;
  int             noOfSegs;
  double PWpval;
  int rows;

#ifdef RETURNPVALS
  *pval=1.0; //Sets up default return for cases when no test can be made.
#endif

  // we count the number of marked rows
  rows=0;
  for(i=0; i<m->len; i++)
    if( m->set[i] == TRUE )
      rows++;

  // if we only have three rows, we cannot separate them
  if( rows <= 3)
    return UNITIG_IS_SHALLOW;

  /* If the number of rows is way to big, we assume that 
     the alignment is repetitive. The test is here too costly */

  if( rows > TEST_UPPER_BOUND ){
    return UNITIG_IS_REPETITIVE;
  }

#if DEBUG > 0
      printf("test_PWsimple start=%d, end=%d\n",start,end);
#endif

  if( end-start >= MIN_MPTEST_LENGTH_REZ ){
      PWpval=Pairwise_Only_Prob_Simple(ali,start,end, m);
#ifdef RETURNPVALS
      *pval=PWpval;
#endif
#if DEBUG > 0
      printf("\n\n>>>>> Aaron's PW evaluation gives p-value %e",PWpval);
#endif
      if(PWpval<=thresh){
#if DEBUG > 0
	printf(" : NOT SIMPLE\n\n");
#endif
	ret = UNITIG_IS_REPETITIVE;
      } else {
#if DEBUG > 0
	printf(" : SIMPLE\n\n");
#endif
	ret = UNITIG_IS_SIMPLE;
      }
  }
  return ret;  
}



/* this is the main test function for a unitig. 
   It returns 
   UNITIG_IS_SIMPLE     if the unitig is simple
   UNITIG_IS_UNKNOWN    if we cannot determine the status
                        since there is no big enough segment
		        or the critical value is below the
	   	        minimum critical value.
   UNITIG_IS_REPETITIVE if the test could determine that
                        the unitig is repetitive.
   UNITIG_IS_SHALLOW if the number of rows is smaller than 4
   The Alignment_t ali contains the alignment in which the list of segments
   can be inspected in order to find repetitive segments
*/


UnitigStatus_t is_IUM_simple(IntUnitigMesg* ium, FragStoreHandle handle,
			     Alignment_t **ali, double thresh, int variant
#ifdef RETURNPVALS
, double *pval
#endif
){
  int i;
  int rows;
  int noOfSegs;
  TestSegment_t **t;
  Marker_t *m;
  int critical;
  UnitigStatus_t simple;
  UnitigStatus_t ret = UNITIG_IS_SIMPLE;
  
  *ali = convert_IUM_to_alignment(ium,handle);

  /* for this test we allocate a marker that is by default TRUE
     for all columns */
  m = allocate_marker((*ali)->rows);

  count_columns(*ali,m);

  if( (*ali)->rows < 4)
    ret = UNITIG_IS_SHALLOW;
  else
    if( variant == 0 ){
#ifdef RETURNPVALS
      double segpval;
#endif
      /* divide the alignment into segments */
      t = get_segments(*ali,&noOfSegs);
      
      /* test each segment, whether it is repetitive. If one
	 segment is repetitive (or UNKNOWN) return this value */
      for(i=0; i<noOfSegs; i++)
	{
	  int start = t[i]->start;
          int end   = t[i]->end; 

#ifdef RETURNPVALS
	  simple = test_simple(*ali,thresh,m,start,end,&critical,&segpval);
	  *pval=min(*pval,segpval);
#else
	  simple = test_simple(*ali,thresh,m,start,end,&critical);
#endif

    
	  if( simple == UNITIG_IS_REPETITIVE )
	    {
	      ret = simple;
	      break;
	    }
	}
      for(i=0; i<noOfSegs; i++)
	free_testSegment(t[i]);
      free(t);
    }
    else
      if( variant == 1 ){
#ifdef RETURNPVALS
	double winpval;
#endif
	for(i=0; i<(*ali)->cols-601; i+=300){
#if DEBUG > 0
	  printf("start = %d, end = %d\n",i,i+600);
#endif
#ifdef RETURNPVALS
	  simple = test_simple(*ali,thresh,m,i,i+600,&critical,&winpval);
	  *pval=min(*pval,winpval);
#else
	  simple = test_simple(*ali,thresh,m,i,i+600,&critical);
#endif

	  
	  if( simple == UNITIG_IS_REPETITIVE ){
	    ret = simple;
	    break;
	  }
	}
	if( ret != UNITIG_IS_REPETITIVE ) {
#ifdef RETURNPVALS
	  ret = test_simple(*ali,thresh,m,i,(*ali)->cols,&critical,&winpval);
	  *pval=min(*pval,winpval);
#else
	  ret = test_simple(*ali,thresh,m,i,(*ali)->cols,&critical);
#endif
	}
      }
  /* clean up dynamically allocated stuff */
  free_marker(m);
 
  return ret;
}




/* this is the main test function for a unitig. 
   It returns 
   UNITIG_IS_SIMPLE     if the unitig is simple
   UNITIG_IS_UNKNOWN    if we cannot determine the status
                        since there is no big enough segment
		        or the critical value is below the
	   	        minimum critical value.
   UNITIG_IS_REPETITIVE if the test could determine that
                        the unitig is repetitive.
   UNITIG_IS_SHALLOW if the number of rows is smaller than 4
   The Alignment_t ali contains the alignment in which the list of segments
   can be inspected in order to find repetitive segments
*/


UnitigStatus_t is_IUM_MPsimple(IntUnitigMesg* ium, FragStoreHandle handle,
			       Alignment_t **ali, double thresh, int variant
#ifdef RETURNPVALS
, double *pval
#endif
){
  int i;
  int rows;
  int noOfSegs;
  TestSegment_t **t;
  Marker_t *m;
  UnitigStatus_t simple;
  UnitigStatus_t ret = UNITIG_IS_SIMPLE;
  
#ifdef RETURNPVALS
  *pval=1.0; //Sets up default return for cases when no test can be made.
#endif

  *ali = convert_IUM_to_alignment(ium,handle);

  /* for this test we allocate a marker that is by default TRUE
     for all columns */
  m = allocate_marker((*ali)->rows);
  count_columns(*ali,m);

  if( (*ali)->rows < 4)
    ret = UNITIG_IS_SHALLOW;
  else
    {
      if( variant == 2){
	/* first test on the whole alignment */

#if STATPRINT > -1
	doPrintMPobs=1;
#endif

	ret = test_MPsimple(*ali,thresh,m,0,(*ali)->cols-1
#ifdef RETURNPVALS
,pval
#endif
);
#if STATPRINT > -1
	doPrintMPobs=0;
#endif
      }
      else
	if( variant == 0 ){
#ifdef RETURNPVALS
	  double segpval;
#endif

	  /* divide the alignment into segments */
	  t = get_segments(*ali,&noOfSegs);
	
	  /* test each segment, whether it is repetitive. If one
	   segment is repetitive (or UNKNOWN) return this value */
#if DEBUG > 0
	  printf("No of Seqs = %d\n",noOfSegs);
#endif
	  for(i=0; i<noOfSegs; i++){
	    int start = t[i]->start;
	    int end   = t[i]->end; 
#if DEBUG > 0
	    printf("start = %d, end = %d\n",start,end);
#endif
#ifdef RETURNPVALS
	    simple = test_MPsimple(*ali,thresh,m,start,end,&segpval);
	    *pval=min(*pval,segpval);
#else
	    simple = test_MPsimple(*ali,thresh,m,start,end);
#endif
	    if( simple == UNITIG_IS_REPETITIVE ){
	      ret = simple;
	      break;
	    }
	  }
	  for(i=0; i<noOfSegs; i++)
	    free_testSegment(t[i]);
	  free(t);
	}
      else
	if( variant ==  1 ){
#ifdef RETURNPVALS
	  double winpval;
#endif
	  for(i=0; i<(*ali)->cols-601; i+=300){
#if DEBUG > 0
	    printf("start = %d, end = %d\n",i,i+600);
#endif

#ifdef RETURNPVALS
	    simple = test_MPsimple(*ali,thresh,m,i,i+600,&winpval);
	    *pval=min(*pval,winpval);
#else
	    simple = test_MPsimple(*ali,thresh,m,i,i+600);
#endif
	    
	    if( simple == UNITIG_IS_REPETITIVE ){
	      ret = simple;
	      break;
	    }
	  }
          if( ret != UNITIG_IS_REPETITIVE ){
#ifdef RETURNPVALS
	    ret = test_MPsimple(*ali,thresh,m,i,(*ali)->cols,&winpval);
  	    *pval=min(*pval,winpval);
#else
	    ret = test_MPsimple(*ali,thresh,m,i,(*ali)->cols);
#endif
          }
	}
    }
  /* clean up dynamically allocated stuff */
  free_marker(m);

  return ret;
}




UnitigStatus_t is_IUM_PWsimple(IntUnitigMesg* ium, FragStoreHandle handle,
			       Alignment_t **ali, double thresh, int variant
#ifdef RETURNPVALS
, double *pval
#endif
){
  int i;
  int rows;
  int noOfSegs;
  TestSegment_t **t;
  Marker_t *m;
  UnitigStatus_t simple;
  UnitigStatus_t ret = UNITIG_IS_SIMPLE;
  
#ifdef RETURNPVALS
  *pval=1.0; //Sets up default return for cases when no test can be made.
#endif

  *ali = convert_IUM_to_alignment(ium,handle);

  /* for this test we allocate a marker that is by default TRUE
     for all columns */
  m = allocate_marker((*ali)->rows);
  count_columns(*ali,m);

  if( (*ali)->rows < 4)
    ret = UNITIG_IS_SHALLOW;
  else
    {
      if( variant == 2){
	/* first test on the whole alignment */

#if STATPRINT > 0
	doPrintPWchisq=1;
#endif

#ifdef RETURNPVALS
	ret = test_PWsimple(*ali,thresh,m,0,(*ali)->cols-1,pval);
#else
	ret = test_PWsimple(*ali,thresh,m,0,(*ali)->cols-1);
#endif

#if STATPRINT > 0
	doPrintPWchisq=0;
#endif
      }
      else
	if( variant == 0 ){
#ifdef RETURNPVALS
	  double segpval;
#endif
	  /* divide the alignment into segments */
	  t = get_segments(*ali,&noOfSegs);
	
	  /* test each segment, whether it is repetitive. If one
	   segment is repetitive (or UNKNOWN) return this value */
#if DEBUG > 0
	  printf("No of Seqs = %d\n",noOfSegs);
#endif
	  for(i=0; i<noOfSegs; i++){
	    int start = t[i]->start;
	    int end   = t[i]->end; 
#if DEBUG > 0
	    printf("start = %d, end = %d\n",start,end);
#endif
#ifdef RETURNPVALS
	    simple = test_PWsimple(*ali,thresh,m,start,end,&segpval);
	    *pval=min(*pval,segpval);
#else
	    simple = test_PWsimple(*ali,thresh,m,start,end);
#endif
	    
	    if( simple == UNITIG_IS_REPETITIVE ){
	      ret = simple;
	      break;
	    }
	  }
	  for(i=0; i<noOfSegs; i++)
	    free_testSegment(t[i]);
	  free(t);
	}
      else
	if( variant ==  1){
#ifdef RETURNPVALS
	  double winpval;
#endif
	  for(i=0; i<(*ali)->cols-601; i+=300){
#if DEBUG > 0
	    printf("start = %d, end = %d\n",i,i+600);
#endif
#ifdef RETURNPVALS
	    simple = test_PWsimple(*ali,thresh,m,i,i+600,&winpval);
	    *pval=min(*pval,winpval);
#else
	    simple = test_PWsimple(*ali,thresh,m,i,i+600);
#endif
	    
	    if( simple == UNITIG_IS_REPETITIVE ){
	      ret = simple;
	      break;
	    }
	  }
          if( ret != UNITIG_IS_REPETITIVE ){

#ifdef RETURNPVALS
	    ret = test_PWsimple(*ali,thresh,m,i,(*ali)->cols,&winpval);
	    *pval=min(*pval,winpval);
#else
	    ret = test_PWsimple(*ali,thresh,m,i,(*ali)->cols);
#endif
         }
	}
    }
  /* clean up dynamically allocated stuff */
  free_marker(m);

  return ret;
}





// AARONS code

double choose(int n, int i){
  return(exp(lgamma(n+1.)-lgamma(i+1.)-lgamma(n-i+1.)));
}

double Poisson(int Obs,double Exp){
  int i;
  double complresult=0,complresult2=0,exact;
  #define MINP 1e-20
  for(i=Obs-1;i>=Exp;i--){
    exact=exp( i * log(Exp)- lgamma((double)(i+1)));
    complresult+=exact;
    //    printf("Subtracting Pr(%d|%f)= %e -> %e\n",i,Exp,exact/exp(Exp),1-complresult/exp(Exp));
  }

  for(i=0;i<Exp&&i<Obs;i++){
    exact=exp( i * log(Exp)- lgamma((double)(i+1)));
    complresult2+=exact;
    //    printf("Subtracting Pr(%d|%f)= %e -> %e\n",i,Exp,exact/exp(Exp),1-(complresult-complresult2)/exp(Exp));
  }

  complresult+=complresult2;

  complresult/=exp(Exp);
  // printf("Poisson final: %e\n",1-complresult);
  if(1-complresult<MINP){
    return(MINP);
  } else {
    return(1-complresult);
  }
}

MPSTAT MP_score_alignment(Alignment_t *alignment,double erate, int s, int e){
  MPSTAT result;
  double tmp;
  int col;
  int a,c,g,t,d,n,A,C,G,T,D,i;

  result.Obs=0;
  result.Exp=0;

  for(i=0;i<200;i++)
    ExpectedSavedSteps[i]=-1;

  for(col=s;col<e;col++){
    A=alignment->countA[col];
    C=alignment->countC[col];
    G=alignment->countG[col];
    T=alignment->countT[col];
    D=alignment->countDash[col];
    n=A+C+G+T+D;
    if(n<4)continue; /* No way to save a step with < 4 sequences */
    a=max(A-1,0);
    c=max(C-1,0);
    g=max(G-1,0);
    t=max(T-1,0);
    d=max(D-1,0);
    result.Obs+=a+c+g+t+d-maxfive(a,c,g,t,d);
#ifdef EXACT_EXPECTED_SAVEDSTEPS
    result.Exp+=expected_savedSteps(n,erate);
#else
    result.Exp+=choose(n,2)*pow(1.-erate/4.,n-2.);
#endif
  }
#ifndef EXACT_EXPECTED_SAVEDSTEPS
  result.Exp*=erate*erate/4.;
#endif
  //printf("n %d Obs %d Exp %e\n",n,result.Obs,result.Exp);
  result.pr=Poisson(result.Obs,result.Exp);
  //printf("--> pr %e\n",result.pr);
  return(result);
}
//END CHUNK OF AARON'S CODE

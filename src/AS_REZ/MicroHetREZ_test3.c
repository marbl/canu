
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
static char CM_ID[] = "$Id: MicroHetREZ_test3.c,v 1.1.1.1 2004-04-14 13:53:23 catmandew Exp $";

#include "MicroHetREZ_test3.h"
#include "MicroHetScoreREZ_test3.h"
#include "MicroHetPartitionsREZ_test3.h"
#include "MicroHetInterfaceREZ_test3.h"
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
#include <float.h>
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_distStore.h"

//Stuff to control diagnostic output.
#define DEBUG -1
#define STATPRINT -1
static int doPrintMPobs=0;


/* I don't think doPrintPWchisq is needed for "test3" and it is redundant with
   the definition in MicroHetREZ.c -- but I'm leaving it in, albeit commented
   out, in case I'm missing something!  ALH 9/20/2000 */
//int doPrintPWchisq=0;

//global: expected number of save steps (to reuse previously-calculated values)
double ExpectedSavedSteps[200];

//minimum p-value that Poisson probabilities can get right?
#define MINP 1e-16 

/* memory allocation and deallocation for Alignment_t */
Alignment_t *AS_REZ_allocate_alignment(int c, int r)
{
  int i;
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




void AS_REZ_free_alignment(Alignment_t* a)
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

void AS_REZ_print_alignment(Alignment_t *a,  int w)
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


void AS_REZ_count_columns(Alignment_t* a, Marker_t* m)
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




UnitigStatus_t AS_REZ_test_MPsimple(Alignment_t *ali, double thresh, Marker_t* m, 
			     int start, int end,double *pval)
{
  int i;
  int ret = UNITIG_IS_SIMPLE;
  int rows;

  *pval=1.0; //Sets up default return for cases when no test can be made.

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
#define PVAL_FOR_TOO_DEEP 1e-20
    *pval=PVAL_FOR_TOO_DEEP;
    return UNITIG_IS_REPETITIVE;
  }

#if DEBUG > 0
      printf("test_MPsimple start=%d, end=%d\n",start,end);
#endif

  if( end-start >= MIN_MPTEST_LENGTH_REZ ){
    double heurSeqErr = AS_REZ_guess_seqErr(ali,m,start,end);
    MPSTAT mpresult;
#if DEBUG > 0
    printf("GuessedSeqErr=%lf for start=%d, end=%d\n",heurSeqErr,start,end);
#endif
    mpresult=AS_REZ_MP_score_alignment(ali,heurSeqErr,start,end);

    if(doPrintMPobs)
      printf("O = %d E = %f pr = %e\n",mpresult.Obs, mpresult.Exp,mpresult.pr);

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


UnitigStatus_t AS_REZ_is_IUM_MPsimple(IntUnitigMesg* ium, FragStoreHandle handle,
                               tFragStorePartition *phandle,
			       Alignment_t **ali, double thresh, int variant, 
			       double *pval)
{
  Marker_t *m;
  UnitigStatus_t ret = UNITIG_IS_SIMPLE;
  
  *pval=1.0; //Sets up default return for cases when no test can be made.

  *ali = AS_REZ_convert_IUM_to_alignment(ium,handle,phandle);

  /* for this test we allocate a marker that is by default TRUE
     for all columns */

  if( (*ali)->rows < 4)
    ret = UNITIG_IS_SHALLOW;
  else
    {

      m = AS_REZ_allocate_marker((*ali)->rows);
      AS_REZ_count_columns(*ali,m);


#if STATPRINT > -1
	doPrintMPobs=1;
#endif

	ret = AS_REZ_test_MPsimple(*ali,thresh,m,0,(*ali)->cols-1,pval);

#if STATPRINT > -1
	doPrintMPobs=0;
#endif

	/* clean up dynamically allocated stuff */
	AS_REZ_free_marker(m);

    }

  return ret;
}




/* I don't think choose() is needed for "test3" and it is redundant with
   the definition in MicroHetREZ.c -- but I'm leaving it in, albeit commented
   out, in case I'm missing something!  ALH 9/20/2000 */
//double choose(int n, int i){
//  return(exp(lgamma(n+1.)-lgamma(i+1.)-lgamma(n-i+1.)));
//}

double AS_REZ_Poisson_prob(int Obs,double Exp){
  int i;
  double complresult=0,complresult2=0,exact;
  int overflow=0;
  double expfactor;
  double maxexpfactor;
  double logExp;

  maxexpfactor=log(DBL_MAX);

  logExp=log(Exp);

  // break the following for-loops into two stages so that we can sum 
  // increasing values--might minimize rounding errors?

  for(i=Obs-1;i>=Exp;i--){
    expfactor=i * logExp- lgamma((double)(i+1))-Exp;
    if (expfactor>maxexpfactor){
      overflow=1;
      break;
    }
    if(expfactor<-200)continue; // too small to have an appreciable effect ...
    exact=exp(expfactor);
    if(DBL_MAX-exact<complresult){
      overflow=1;
      break;
    }
    complresult+=exact;
  }

  for(i=0;i<Exp&&i<Obs;i++){
    expfactor=i * logExp- lgamma((double)(i+1))-Exp;
    if (expfactor>maxexpfactor){
      overflow=1;
      break;
    }
    if(expfactor<-200)continue; // too small to have an appreciable effect ...
    exact=exp(expfactor);
    if(DBL_MAX-exact<complresult2){
      overflow=1;
      break;
    }
    complresult2+=exact;
  }

  if(!overflow){
    complresult+=complresult2;
    if(1-complresult<MINP){
      return(MINP);
    } else {
      return(1-complresult);
    }
  } else {
    assert(0); // Sorry, I did not handle this overflow ...
               // The obvious way to do it would be to 
               // use the normal approximation to the Poisson.
  }
}

double AS_REZ_Poisson(int Obs,double Exp){
  int i;
  double complresult=0;
  double lambdaterm=1;
  double facterm=1;
  int inexact=0;

  for(i=0;i<Obs;){
    if(Exp>1){
      if(lambdaterm>DBL_MAX/Exp){
	inexact=1;
      }
    } else {
      if(lambdaterm<1e-200/Exp){
	inexact=1;
      }
    }
    if(facterm>DBL_MAX/(++i)){
      inexact=1;
    }
    if(inexact){
      #if DEBUG > -1
      fprintf(stderr,"Max'ed out Poisson probability calculation: obs: %d exp: %e\n",Obs,Exp);
      #endif
      return(AS_REZ_Poisson_prob(Obs,Exp));
    }
    complresult+=lambdaterm*exp(-Exp)/facterm;
    lambdaterm*=Exp;
    facterm*=i;
  }
  if(1-complresult<MINP){
    return(MINP);
  } else {
    return(1-complresult);
  }
}

MPSTAT AS_REZ_MP_score_alignment(Alignment_t *alignment,double erate, int s, int e){
  MPSTAT result;
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
#define USESNPMODEL
#ifdef USESNPMODEL
#  define SNPPROB .001
    result.Exp+=(1.-SNPPROB)*(double)AS_REZ_expected_savedSteps(n,erate)+
      SNPPROB*(double)(((int) (n/2))-1);
#else
    result.Exp+=AS_REZ_expected_savedSteps(n,erate);
#endif
  }
  result.pr=AS_REZ_Poisson(result.Obs,result.Exp);
  return(result);
}

/* AS_REZ_MP_MicroHet_prob 

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
             will be obtained (-1 (NULLFRAGSTOREHANDLE) if partitioned store is used)
   phandle  : the partitioned fragStore from which locale information for each fragment iid
             will be obtained (NULL if traditional non-partitioned store is used)
   len     : number of columns in the multialignment
   depth   : number of rows in the multialignment
*/

double AS_REZ_MP_MicroHet_prob(char **bqarray,int **idarray,FragStoreHandle handle,
                               tFragStorePartition *phandle,int len,int depth){
  double pvalue;
  UnitigStatus_t result;
  Marker_t *m;
  double thresh=1e-3; /* reasonable value; doesn't actually do anything
                        here except make AS_REZ_test_MPsimple() happy, 
			but cf. AS_REZ_is_IUM_MPsimple() */
  Alignment_t *ali = AS_REZ_convert_array_to_alignment(bqarray,len,depth);
  compress_shreds_and_null_indels(len,depth,handle,phandle,ali->ali,idarray,0);
  if(ali->rows<4){
    pvalue=1.0;
    result=UNITIG_IS_SHALLOW;
  } else {
    m= AS_REZ_allocate_marker(ali->rows);
    AS_REZ_count_columns(ali,m);
#if STATPRINT > -1
    doPrintMPobs=1;
#endif
    result=AS_REZ_test_MPsimple(ali,thresh,m,0,ali->cols-1,&pvalue);
#if STATPRINT > -1
    doPrintMPobs=0;
#endif
    AS_REZ_free_marker(m);
  }
  AS_REZ_free_alignment(ali);
  return(pvalue);
}

/* this is the main test function for a unitig. 
   It returns a pvalue (roughly, a probability) that the unitig is simple
   (mismatches are random errors).
   A return value of 1.0 may indicate that the unitig was not deep enough for meaningful test.
*/
double AS_REZ_prob_IUM_MPsimple(IntUnitigMesg* ium, FragStoreHandle handle,tFragStorePartition *phandle)
{
  int i;
  int rows;
  double ret;

  char **bqarray;
  int **idarray;
  int **oriarray;
  IMP2Array(ium->f_list,ium->num_frags,ium->length,
	    handle,phandle,NULLFRAGSTOREHANDLE,&rows,&bqarray,&idarray,&oriarray,0,READSTRUCT_LATEST);

  ret = AS_REZ_MP_MicroHet_prob(bqarray,idarray,handle,phandle, ium->length,rows);
  /* free the space that is allocated by IMP2Array */
  for(i=0; i<2*rows; i++)
    free(bqarray[i]);
  free(bqarray);

  for(i=0; i<rows; i++)
    free(idarray[i]);
  free(idarray);

  for(i=0; i<rows; i++)
    free(oriarray[i]);
  free(oriarray);

  return(ret);
}

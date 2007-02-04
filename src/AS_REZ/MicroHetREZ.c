
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
static char CM_ID[] = "$Id: MicroHetREZ.c,v 1.8 2007-02-04 09:30:46 brianwalenz Exp $";

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

#include "MicroHetREZ.h"

#include "UtilsREZ.h"
#include "PrimitiveVA_MSG.h"

// The minimum length for which we test for Aarons MP score
#define MIN_MPTEST_LENGTH_REZ 50

// We assume that an alignment with more rows than that is repetitive
#define TEST_UPPER_BOUND 40


#define maxfive(v,w,x,y,z) MAX(v,MAX(MAX(w,x),MAX(y,z)))


//minimum p-value that Poisson probabilities can get right?
#define MINP 1e-16 




void AS_REZ_free_marker(Marker_t *m) {
  free(m->set);
  free(m);
}

/* functions to allocate and free a marker of size l */
Marker_t *AS_REZ_allocate_marker(int l) {
  int i;
  Marker_t* m = (Marker_t*) safe_malloc(sizeof(Marker_t));
  m->set      = (int*) safe_calloc(sizeof(int),l);
  for(i=0; i<l; i++)
    m->set[i] = TRUE;
  m->len      = l;
  return m;
}

void AS_REZ_print_marker(Marker_t *m){
  int i;
  for(i=0; i<m->len; i++)
    printf("|%d|",m->set[i]);
  printf("\n");
}







/* memory allocation and deallocation for Alignment_t */
Alignment_t *AS_REZ_allocate_alignment(int c, int r) {
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


void AS_REZ_free_alignment(Alignment_t* a) {
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


void AS_REZ_print_alignment(Alignment_t *a,  int w) {
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



void AS_REZ_get_info(CDS_IID_t iid,
                     FragStoreHandle frag_store,
                     tFragStorePartition *pfrag_store,
                     CDS_UID_t *locale,
                     uint32 *beg, uint32 *end,
                     FragType *type,
                     VA_TYPE(CDS_UID_t) *locales,
                     VA_TYPE(uint32) *fragtype,
                     VA_TYPE(uint32) *locbeg,
                     VA_TYPE(uint32) *locend,
                     ReadStructp input){
  if( iid != 0 ){ // blank positions have IID 0

    uint32 *mytype = Getuint32(fragtype,iid);
    if ((mytype) && (*mytype == 0))
      mytype = NULL;

    if(!mytype ){ // we have not seen that iid before 
      if ( frag_store != NULLFRAGSTOREHANDLE ) {
        if(getFragStore(frag_store,iid,FRAG_S_FIXED,input) != 0)
          assert(0);
      } else {
        if(getFragStorePartition(pfrag_store,iid,FRAG_S_FIXED,input) != 0)
          assert(0);
      }
      getReadType_ReadStruct(input,type); 
      getLocID_ReadStruct(input,locale);
      getLocalePos_ReadStruct(input,beg,end);
      
      // we store the locales for later lookup
      SetCDS_UID_t(locales,iid,locale);
      Setuint32(locbeg,iid,beg);
      Setuint32(locend,iid,end);
      Setuint32(fragtype,iid,(uint32 *)type);
    }
    else{
      // We know these entries have been initialized, so we can use them directly
      // without going through the fetch routines that behave special for uninitialized entries.
      *locale = *GetCDS_UID_t(locales,iid);
      *beg    = *Getuint32(locbeg,iid);	
      *end    = *Getuint32(locend,iid);	
      *type   = *Getuint32(fragtype,iid);
    }
  }
  else{
    *beg = *end = *locale = 0;
    *type = 'X';
  }
}




// These statics accumulate the information for each
// fragment seen in this run of the Microhetrez tool
// they will grow to the length of the number of frags
// in the entire genome (not the number in this partition). SAK
//
static  VA_TYPE(CDS_UID_t)  *locales  = NULL;
static  VA_TYPE(uint32)  *locbeg   = NULL;
static  VA_TYPE(uint32)  *locend   = NULL;
static  VA_TYPE(uint32)  *fragtype = NULL;
static  ReadStructp input;

/* This function compresses shredded fragments from the same location 
 * into basically a 1x coverage, such that there are no aritfical microhets;
 * it also nulls out all but one position of a multibase gap, to reduce the
 * effect of multibase indel polymorphisms on microhet detection.
 */
void compress_shreds_and_null_indels(int c,
                                     int r,
                                     FragStoreHandle frag_store, 
                                     tFragStorePartition *pfrag_store,
                                     char **array,
                                     int **id_array,
                                     int verbose){
  int i,j,k;
  CDS_UID_t l1,l2;
  uint32 b1,b2;
  uint32 e1,e2;
  CDS_IID_t iid1,iid2;
  FragType t1,t2;

  if(locales == NULL){
    locales = CreateVA_CDS_UID_t(10000);
    locbeg = CreateVA_uint32(10000);
    locend = CreateVA_uint32(10000);
    fragtype = CreateVA_uint32(10000);
    input = new_ReadStruct();
  }

  if(verbose > 0)
    printf("Alignment has %d columns and %d rows\n",c,r);

  if( verbose > 0 )
    for(i=0; i<c; i++)
      for(j=0; j<r; j++)
	printf("first loop (%d,%d) iid1=%d\n",i,j,id_array[j][i]);


  for(i=0; i<c; i++){
    if(i==c-1){

      for(j=0; j<r-1; j++){
	iid1 = id_array[j][i];
	if( verbose > 0 )
	  printf("loop (%d,%d) iid1=" F_IID "\n",i,j,iid1);
	if( iid1 != 0 ){ // blank positions have IID 0
	  AS_REZ_get_info(iid1,frag_store,pfrag_store,&l1,&b1,&e1,&t1,locales,fragtype,locbeg,locend,input);
	  if( verbose > 0 )
	    printf("(id,loc,beg,end,type) = (" F_IID "," F_UID ",%d,%d,%c)\n",
                   iid1,l1,b1,e1,t1);
	  
	  if(AS_FA_SHREDDED(t1)){ 
	    for(k=j+1; k<r; k++){ // now eliminate the duplicates
	      if( verbose > 0 )
		printf("pos (%d,%d)\n",i,k);
	      iid2 = id_array[k][i];
	      AS_REZ_get_info(iid2,frag_store,pfrag_store,&l2,&b2,&e2,&t2,locales,fragtype,locbeg,locend,input);
	      if( verbose > 0)
		printf("(id2,loc2,beg2,end2,type2) = (" F_IID "," F_UID ",%d,%d,%c)\n",
                       iid2,l2,b2,e2,t2);
	      
	      if( t1 == t2  && l2 == l1 ){
		uint32 min2,max1;
		if( verbose > 0 )
		  printf("changed position (%d,%d) from %c to blank, types (%c,%c), locales (" F_UID "," F_UID ")\n",i,k,array[i][k],t1,t2,l1,l2);
		min2 = ( b2 < e2 ? b2 : e2);
		max1 = ( b1 > e1 ? b1 : e1);
		
		if(max1 < min2)
		  ;//		 fprintf(stdout,"WARNING begin of frag2 is less than end of second !!\n");
		else{
		  id_array[k][i] = 0;
		  array[i][k] = ' ';
		}
	      }
	    }
	  }
	}
      }
    } else {
      for(j=0; j<r; j++){

	if(array[i][j]=='-'&&array[i+1][j]=='-'){
	  //	  fprintf(stdout,"Multibase gap!\n");
	  array[i][j]=' ';
	}	

	if(j<r-1){

	  iid1 = id_array[j][i];
	  if( verbose > 0 )
	    printf("loop (%d,%d) iid1=" F_IID "\n",i,j,iid1);
	  if( iid1 != 0 ){ // blank positions have IID 0
	    AS_REZ_get_info(iid1,frag_store,pfrag_store,&l1,&b1,&e1,&t1,locales,fragtype,locbeg,locend,input);
	    if( verbose > 0 )
	      printf("(id,loc,beg,end,type) = (" F_IID "," F_UID ",%d,%d,%c)\n",
                     iid1,l1,b1,e1,t1);
	    
	    if(AS_FA_SHREDDED(t1)){ 
	      for(k=j+1; k<r; k++){ // now eliminate the duplicates
		if( verbose > 0 )
		  printf("pos (%d,%d)\n",i,k);
		iid2 = id_array[k][i];
		AS_REZ_get_info(iid2,frag_store,pfrag_store,&l2,&b2,&e2,&t2,locales,fragtype,locbeg,locend,input);
		if( verbose > 0)
		  printf("(id2,loc2,beg2,end2,type2) = (" F_IID "," F_UID ",%d,%d,%c)\n",
                         iid2,l2,b2,e2,t2);
		
		if( t1 == t2  && l2 == l1 ){
		  uint32 min2,max1;
		  if( verbose > 0 )
		    printf("changed position (%d,%d) from %c to blank, types (%c,%c), locales (" F_UID "," F_UID ")\n",
                           i,k,array[i][k],t1,t2,l1,l2);
		  min2 = ( b2 < e2 ? b2 : e2);
		  max1 = ( b1 > e1 ? b1 : e1);
		  
		  if(max1 < min2)
		    ;//		 fprintf(stdout,"WARNING begin of frag2 is less than end of second !!\n");
		  else{
		    id_array[k][i] = 0;
		    array[i][k] = ' ';
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}










/* this function converts an array computed by consensus into
 * an Alignment_t. The function that converts an IMP message 
 * into the array is called IMP2Array. The array contains in the odd
 * rows the sequence and in the even rows the quality values.
 */
Alignment_t *AS_REZ_convert_array_to_alignment(char **ar, int c, int r){
  int i,j;

  //  fprintf(stderr,"Working on an array of size rXc %d X %d\n",r,c);
  Alignment_t *a = AS_REZ_allocate_alignment(c,r);
  for(i=0; i<c; i++)
    for(j=0; j<r; j++)
      {
	assert(ar[2*j][i] == 'A' || ar[2*j][i] == 'C' || 
               ar[2*j][i] == 'G' || ar[2*j][i] == 'T' ||
               ar[2*j][i] == '-' || ar[2*j][i] == ' ' ||
               ar[2*j][i] == 'N');

	a->ali[i][j] = ar[2*j][i]; 
	
#define USE_QUALITY
#ifdef USE_QUALITY
	a->hasQuality = TRUE;
	if( ar[2*j][i] != ' '){
	  a->seqErrArray[i][j] = 1.0/pow(10,((double)ar[2*j+1][i]-'0')/10.0);
	  //	  printf("Seq = %c, Quality = %c, Seq err %f\n",ar[2*j][i],ar[2*j+1][i],a->seqErrArray[i][j]);
	}
	else
	  a->seqErrArray[i][j] = 0.0;
#else
	a->seqErrArray[i][j] = 0.0;
#endif

      } 

  return a;
}




#define FACLIMIT 160

double AS_REZ_fac(int n) {
  static double facREZ[FACLIMIT] = { 0.0 };

  assert(n < FACLIMIT);

  if( n < 0 )
    return 0.0;
  if( n == 0 || n == 1 )
    facREZ[n] = 1.0;
  if( facREZ[n] == 0.0 )
    facREZ[n] = n*AS_REZ_fac(n-1);
  return facREZ[n];
}
  




/* this computes the probability of a fournomial event 
   the assumption is, that c1,c2,c3 and c4 have the same probability
   of occuring */

double AS_REZ_fournomial(double seqErr, int c1, int c2, int c3, int c4, int n)
{
  register int i;
  double q = seqErr/4.0;
  double p = 1-seqErr;
  
  double prod = 1.0;
  int    sum  = 0;

  sum += c1;
  sum += c2;
  sum += c3;
  sum += c4;

  for(i=0; i<n-sum; i++)
    prod *= p;

  for(i=n-sum+1; i<=n; i++)
    prod *= i;

  for(i=0; i<sum; i++)
    prod *= q; 

  if( c1 >= FACLIMIT )
    {
      for(i=c1; i>=FACLIMIT; i--)
	prod /= i; 
      prod /= AS_REZ_fac(FACLIMIT-1);
    }
  else
    prod /= AS_REZ_fac(c1);

  if( c2 >= FACLIMIT )
    {
      for(i=c2; i>=FACLIMIT; i--)
	prod /= i; 
      prod /= AS_REZ_fac(FACLIMIT-1);
    }
  else
    prod /= AS_REZ_fac(c2);

  if( c3 >= FACLIMIT )
    {
      for(i=c3; i>=FACLIMIT; i--)
	prod /= i; 
      prod /= AS_REZ_fac(FACLIMIT-1);
    }
  else
    prod /= AS_REZ_fac(c3);

  if( c4 >= FACLIMIT )
    {
      for(i=c4; i>=FACLIMIT; i--)
	prod /= i; 
      prod /= AS_REZ_fac(FACLIMIT-1);
    }
  else
    prod /= AS_REZ_fac(c4);
  

  // The products get too big

  /*  prod *= AS_REZ_fac(n);
      prod /= AS_REZ_fac(c1);
      prod /= AS_REZ_fac(c2);
      prod /= AS_REZ_fac(c3);
      prod /= AS_REZ_fac(c4);
      prod /= AS_REZ_fac(n-sum); */

  return prod;
}
 



/* this function returns the expected number of "steps" which can be "saved"
   given the optimal partitioning of a column of data under the null hypothesis
   of a simple alignment */

//global: expected number of save steps (to reuse previously-calculated values)
double ExpectedSavedSteps[200];

double AS_REZ_expected_savedSteps(int r, double seqErr)
{
  double Exp = 0.0;
  int c1,c2,c3,c4,c5,a,c,g,t,d,saved;

  if(r<4)return(0);

  if(ExpectedSavedSteps[r]>=0)return(ExpectedSavedSteps[r]);

  for(c1=0; c1<=r; c1++)
    for(c2=0; c2<=r-c1; c2++)
      for(c3=0; c3<=r-c1-c2; c3++)
	for(c4=0; c4<=r-c1-c2-c3; c4++)
	  {
	    c5=r-c1-c2-c3-c4;

	    /* count the number of steps saved for these conditions */
	    a=MAX(c1-1,0);
	    c=MAX(c2-1,0);
	    g=MAX(c3-1,0);
	    t=MAX(c4-1,0);
	    d=MAX(c5-1,0);
	    saved=a+c+g+t+d-maxfive(a,c,g,t,d);

	    /* saved is the number of steps saved for the condition */
	    /* the probability of observing c1..c5 can be obtained by 
	       assuming that c5 is the true state without loss of generality */
	    /* So, expected is: */

	    Exp += saved * AS_REZ_fournomial(seqErr,c1,c2,c3,c4,r);

	  }

  ExpectedSavedSteps[r]=Exp;
  return Exp;
}


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
    a=MAX(A-1,0);
    c=MAX(C-1,0);
    g=MAX(G-1,0);
    t=MAX(T-1,0);
    d=MAX(D-1,0);
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



/* this function inspects a subalignment and -- assuming all changes are 
   due to sequencing error -- returns the mean of the most likely explanation 
   If we happen do have quality values, these are used to compute an exspected
   mean sequencing error
*/
double AS_REZ_guess_seqErr(Alignment_t *a, Marker_t* m, int s, int e)
{
  double seqErr = 0.0;
  int i,j,k;
  int rows = 0;
  int charcount = 0;
  int count[6]; // A C G T Dash (Blank or N)
  for(i=0; i<a->rows; i++)
    if( m->set[i] == TRUE )
      rows++;

  for(i=s; i<e; i++){
    int max;
    for(k=0; k<6; k++)
      count[k] = 0;
    for(j=0; j<a->rows; j++)
      {
        if( m->set[j] == TRUE )
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
	      case ' ' :
		count[5]++;
		break;
	      case 'N' :
		count[5]++;
		break;
            }	
      }	
      
    /* we counted the letters in the marked rows */
    /* now we determine the majority character   */
    max = 0;
    for(k=0; k<5; k++){
      if( count[max] < count[k] )
	max = k;
    }
      
    /* we add up the ratios between the counts of the 
       changed characters and the majority character */
    if( rows-count[5] > 0 )
      seqErr += ((double)(rows-count[5]-count[max]))/((double) (rows-count[5]));
  }
  seqErr /= (e-s+1);

  return seqErr;
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

  if( end-start >= MIN_MPTEST_LENGTH_REZ ){
    double heurSeqErr = AS_REZ_guess_seqErr(ali,m,start,end);
    MPSTAT mpresult;
    mpresult=AS_REZ_MP_score_alignment(ali,heurSeqErr,start,end);

    //printf("O = %d E = %f pr = %e\n",mpresult.Obs, mpresult.Exp,mpresult.pr);

    if(mpresult.pr<=thresh){
      ret = UNITIG_IS_REPETITIVE;
    }
    else {
      ret = UNITIG_IS_SIMPLE;
    }
    *pval=mpresult.pr;
  }

  return ret;  
}











/* AS_REZ_MP_MicroHet_prob 
 *
 * RESULT: The function returns a (double) pvalue (probability) of an
 * input unitig being SIMPLE -- meaning, having mismatches due to
 * randomly distributed sequencing errors.
 *
 * If the returned value is sufficiently small, the unitig should be
 * treated as a likely repeat.
 *
 * A return value of 1.0 may indicate that the unitig was not deep
 * enough for a meaningful test.
 *
 * Some false positives may be induced by polymorphisms; however, the
 * calculation should not be drastically misled by multibase indel
 * polymorphisms.
 *
 * INPUT:
 *
 * bqarray : an array of size [depth*2]*len of bases and quality
 *           values in alternative rows, giving a multialignment
 * idarray : an array of size depth*len giving the fragment iid of
 *           each base in the multialignment
 * handle : the fragStore from which locale information for each
 *          fragment iid will be obtained (-1 (NULLFRAGSTOREHANDLE)
 *          if partitioned store is used)
 * phandle : the partitioned fragStore from which locale information
 *           for each fragment iid will be obtained (NULL if
 *           traditional non-partitioned store is used)
 * len     : number of columns in the multialignment
 * depth   : number of rows in the multialignment
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
    m = AS_REZ_allocate_marker(ali->rows);
    AS_REZ_count_columns(ali,m);

    result=AS_REZ_test_MPsimple(ali,thresh,m,0,ali->cols-1,&pvalue);

    AS_REZ_free_marker(m);
  }
  AS_REZ_free_alignment(ali);
  return(pvalue);
}

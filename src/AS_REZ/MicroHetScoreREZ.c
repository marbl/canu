
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
static char CM_ID[] = "$Id: MicroHetScoreREZ.c,v 1.2 2004-09-23 20:25:28 mcschatz Exp $";

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "UtilsREZ.h"
#include "AS_UTL_rand.h"
#include "MicroHetScoreREZ.h"

#define DEBUG 0
extern double ExpectedSavedSteps[200];
extern int doPrintPWchisq;
/*
void set_xContributing(Alignment_t *a, int s, int e)
{
  int i,j;
  for(i=0; i<MAX_X; i++)
    for(j=s; j<e; j++)
      a->xContributing[i][j] = col_contributing(a,j,i+2); 

}

void set_tfContributing(Alignment_t *a, int r1, int r2, int s, int e)
{
  int i;
  for(i=s; i<e; i++)
    a->tfContributing[i] = col_contributing_two_fixed(a,i,r1,r2);
}
*/


static double fac(int n);

/* this function inspects a subalignment and -- assuming all changes are 
   due to sequencing error -- returns the mean of the most likely explanation 
   If we happen do have quality values, these are used to compute an exspected
   mean sequencing error
*/
double guess_seqErr(Alignment_t *a, Marker_t* m, int s, int e)
{
  double seqErr = 0.0;
  int i,j,k;
  int rows = 0;
  int charcount = 0;
  for(i=0; i<a->rows; i++)
    if( m->set[i] == TRUE )
      rows++;

  if( FALSE ){    
    for(i=s; i<e; i++)
      for(j=0; j<a->rows; j++)
	if( m->set[j] == TRUE )
	  if( a->ali[i][j] != ' '){
	    charcount++;
	    seqErr += a->seqErrArray[i][j];
	  }
    seqErr /= charcount;
  }
  else{ 
    int count[6]; // A C G T Dash (Blank or N)
        
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
  }
  return seqErr;
}



static char N[] = {'-','A','C','G','T' };

/* this function returns the probability of a sequencing
   error turning character at position (i,j) 
   given into character into.
   The sum of the conditional probabilities should be se.
   The probability of not changing the character is (1-se) */

double seq_err(Alignment_t* a, char given, int c, int r)
{	
  char into=a->ali[c][r];
  if( into == given )
    return (1.0-a->seqErrArray[c][r]);
  else
    return 0.25*a->seqErrArray[c][r];
}


/* This function returns the probability that the underlying character
   at position (i,j) is c */
static double prior(char c, int col){
  return 0.2;
}



/* This function returns the number of columns that have the xContributing
   predicate, meaning that there are at least two occurrences of non-majority
   characters that occur at least twice */
int no_col_contributing(Alignment_t *a, int t, int s, int e)
{
  int i;
  int count = 0;
  assert(t-2 < MAX_X);
  for(i=s; i<e; i++)
    if( col_contributing(a,i,t) == TRUE )
      count++;
  return count;
}


/* this function actually computes whether a ncolumn in an alignment is xContributing 
 precondition is that the counts are set according to the desired marker */
int col_contributing(Alignment_t *a, int c, int t)

{
  int observed=0;
 
  if( a->countA[c] >= t)
    observed++;
  if( a->countC[c] >= t)
    observed++;
  if( a->countG[c] >= t)
    observed++;
  if( a->countT[c] >= t)
    observed++;
  if( a->countDash[c] >= t)
    observed++;


  if( observed >= 2 )
    return TRUE;
  else
    return FALSE;    
}

/* this function returns the value k, for which the probability
   density becomes greater than 1.0-thresh, p is the success probability
   of the binomial distribution and n the number of trials */

int critical_binomial_value(double p, double thresh, int n)
{
  int j,k;
  int d;
  double sum = 0.0;

  for(k=0; k<n; k++)
    {
      double prod;

      if(sum > 1.0-thresh)
	break;
      
      prod  = pow(1-p,n-k);
      prod *= pow(p,k);
      
      d = n-k;
      if( d > k )
	{
	  for(j=n; j>d; j--)
	    prod *= j;
	  prod /= fac(k);
	}
      else
	{
	  for(j=n; j>k; j--)
	    prod *= j;
	    prod /= fac(d);
	}
      sum += prod;
      //      printf("SUM=%lf\n",sum);
    }

  return k-1;
}


int no_col_contributing_two_fixed(Alignment_t *a, int r1, int r2, int s, int e)
{
  int i;
  int count = 0;
  for(i=s; i<e; i++)
    if( col_contributing_two_fixed(a,i,r1,r2) == TRUE ){
      count++;
    }
  return count;
}


int col_contributing_two_fixed(Alignment_t *a, int c, 
			       int r1, int r2)
/* precondition is that the counts are set according to the desired marker */
{
  char c1 = a->ali[c][r1];
  char c2 = a->ali[c][r2];
  // we assume that rows r1 and r2 are in the same group
  if( c1 != c2 )
    return FALSE;

  if( c1 == ' ' )
    return FALSE;
    
  switch(c1)
    {
    case 'A' :
      if( a->countC[c] >= a->countA[c] )
	return TRUE;
      if( a->countG[c] >= a->countA[c] )
	return TRUE;
      if( a->countT[c] >= a->countA[c] )
	return TRUE;
      if( a->countDash[c] >= a->countA[c] )
	return TRUE;
      break;	
    case 'C' :
      if( a->countA[c] >= a->countC[c] )
	return TRUE;
      if( a->countG[c] >= a->countC[c] )
	return TRUE;
      if( a->countT[c] >= a->countC[c] )
	return TRUE;
      if( a->countDash[c] >= a->countC[c] )
	return TRUE;
      break;
    case 'G' :
      if( a->countC[c] >= a->countG[c] )
	return TRUE;
      if( a->countA[c] >= a->countG[c] )
	return TRUE;
      if( a->countT[c] >= a->countG[c] )
	return TRUE;
      if( a->countDash[c] >= a->countG[c] )
	return TRUE;
      break;
    case 'T' :
      if( a->countC[c] >= a->countT[c] )
	return TRUE;
      if( a->countG[c] >= a->countT[c] )
	return TRUE;
      if( a->countA[c] >= a->countT[c] )
	return TRUE;
      if( a->countDash[c] >= a->countT[c] )
	return TRUE;
      break;
    case '-' :
      if( a->countC[c] >= a->countDash[c] )
	return TRUE;
      if( a->countG[c] >= a->countDash[c] )
	return TRUE;
      if( a->countT[c] >= a->countDash[c] )
	return TRUE;
      if( a->countA[c] >= a->countDash[c] )
	return TRUE;
      break;
    }
  
  return FALSE;
}








int old_col_contributing_two_fixed(Alignment_t *a, int c, int r1, int r2,
			       Marker_t *m)
{
  int j,k,k1 = 4,k2 = 4,max;
  int rows = a->rows;

  int count[6]; // A C G T Dash Blank

  // we assume that rows r1 and r2 are in the same group
  switch(a->ali[c][r1])
    {
    case 'A' :
      k1 = 0;
      break;	
    case 'C' :
      k1 = 1;
      break;
    case 'G' :
      k1 = 2;
      break;
    case 'T' :
      k1 = 3;
      break;
    case '-' :
      k1 = 4;
      break;
    case ' ' :
      k1 = 5;
      break;
    case 'N' :
      k1 = 5;
      break;
    }	

  switch(a->ali[c][r2])
    {
    case 'A' :
      k2 = 0;
      break;	
    case 'C' :
      k2 = 1;
      break;
    case 'G' :
      k2 = 2;
      break;
    case 'T' :
      k2 = 3;
      break;
    case '-' :
      k2 = 4;
      break;
    case ' ' :
      k2 = 5;
      break;
    case 'N' :
      k2 = 5;
      break;
    }	
  
  if( k1 == 5 || k2 == 5 )
    return FALSE;

  if( k1 != k2 )
    return FALSE;

  for(k=0; k<6; k++)
    count[k] = 0;

  for(j=0; j<rows; j++)
    { 
      if( m->set[j] != TRUE )
	continue;

      switch(a->ali[c][j])
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

  /* assume that the character in the assigned rows is the majority character */
  max = k1;
 
  for(k=0; k<6; k++)
    {
      if( k == k1 )
	continue;
      if( count[max] <= count[k] ) /* if any of the other characters has */
	max = k;                   /* as many occurences as k1 this is fine */
    }
      
  if( max == 5 )
    return FALSE;  

  if( max == k1 )
    return FALSE;
  else
    return TRUE;
}





/* this function computes the log score of a colum */

double col_value(Partition_t* partition, Alignment_t *a, int col, Marker_t* m)
{ 
  int groups = partition->groups;
  int rows   = 0;
  double pValue[groups];
  double value = 1.0;
  int i,j,k;

  for(i=0; i<a->rows; i++)
    if( m->set[i] == TRUE )
      rows++;

  for(i=0; i<groups; i++){
    pValue[i] = 0.0;
  }

  for(i=0; i<groups; i++){  // Compute factor for each group number                                                // in the alignment
    for(j=0; j<5; j++){         // All possible underlying characters
      double prod=1.0;
      for(k=0; k<rows; k++){
	if( partition->part[k] == i )  // if the character in row k is in group i
	  prod *= seq_err(a,N[j],col,k);
      }
      pValue[i] += prior(N[j],col)*prod; // otherwise we add the contributions
    }
    value *= pValue[i];  // in this column
  }
  return -log(value);
}



/* this function computes the log score of all contributing columns of the
   alignment */
double value(Partition_t *partition, Alignment_t *a, int ct, Marker_t* m){
  int i;
  int contCols=0;
  double value = 0.0;

  for(i=0; i<a->cols; i++)
    if( col_contributing(a,i,ct) == TRUE ){
      contCols++;
      value += col_value(partition,a,i,m);
    }
  
  if( contCols > 0 )
    return value/contCols;
  else
    return 0.0;
}


/* this function computes the log score of all contributing columns of the
   alignment */
double value_two_fixed(Partition_t *partition, Alignment_t *a, int r1, int r2,
		       int start, int end, Marker_t* m){
  int i;
  int contCols=0;
  double value = 0.0;

  //  set_tfContributing(a,r1,r2,start,end);
  for(i=start; i<end; i++)
    if( col_contributing_two_fixed(a,i,r1,r2) == TRUE ){
      contCols++;
      value += col_value(partition,a,i,m);
    }
  
  if( contCols > 0 )
    return value/contCols;
  else
    return 0.0;
}


#define FACLIMIT 160

static double facREZ[FACLIMIT];

static double fac(int n)
{
  assert(n < FACLIMIT);
  if( n < 0 )
    return 0.0;
  if( n == 0 || n == 1 )
    facREZ[n] = 1.0;
  if( facREZ[n] == 0.0 )
    facREZ[n] = n*fac(n-1);
  return facREZ[n];
}
  

/* this computes the probability of a fournomial event 
   the assumption is, that c1,c2,c3 and c4 have the same probability
   of occuring */

double fournomial(double seqErr, int c1, int c2, int c3, int c4, int n)
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
      prod /= fac(FACLIMIT-1);
    }
  else
    prod /= fac(c1);

  if( c2 >= FACLIMIT )
    {
      for(i=c2; i>=FACLIMIT; i--)
	prod /= i; 
      prod /= fac(FACLIMIT-1);
    }
  else
    prod /= fac(c2);

  if( c3 >= FACLIMIT )
    {
      for(i=c3; i>=FACLIMIT; i--)
	prod /= i; 
      prod /= fac(FACLIMIT-1);
    }
  else
    prod /= fac(c3);

  if( c4 >= FACLIMIT )
    {
      for(i=c4; i>=FACLIMIT; i--)
	prod /= i; 
      prod /= fac(FACLIMIT-1);
    }
  else
    prod /= fac(c4);
  

  // The products get too big

  /*  prod *= fac(n);
  prod /= fac(c1);
  prod /= fac(c2);
  prod /= fac(c3);
  prod /= fac(c4);
  prod /= fac(n-sum); */

  return prod;
}
 


/* this function returns the 1-p, where p
   is the probability that there are at least two
   different characters in a column that occur 
   at least x times */
double column_prob(int x, int r, double seqErr)
{
  double prob = 0.0;
  int c1,c2,c3,c4;

  for(c1=0; c1<=r; c1++)
    for(c2=0; c2<=r-c1; c2++)
      for(c3=0; c3<=r-c1-c2; c3++)
	for(c4=0; c4<=r-c1-c2-c3; c4++)
	  {
	    /* now we check if at most one of the characters
	       occurs x times. This is the opposite event of what we
	       want. Hence we return 1-p 
	    */
	       
	    int count = 0;
	    if( c1 >= x )
	      count++;
	    if( c2 >= x )
	      count++;
	    if( c3 >= x )
	      count++;
	    if( c4 >= x )
	      count++;
	    if( r-c1-c2-c3-c4 >= x )
	      count++;
	    if( count < 2){
	      //	      printf("COL_PROB : %d,%d,%d,%d,%d\n",c1,c2,c3,c4,r-c1-c2-c3-c4);
	      prob += fournomial(seqErr,c1,c2,c3,c4,r);
	    }
	  }
  return prob;
 }




/* this function computes the probability that for two fixed rows
   in a column the two characters are the same, and that it is different
   from the majority character of that column */
double column_prob_two_fixed(int r, double seqErr)
{
  double prob1 = 0.0;
  double prob2 = 0.0;

  int c1,c2,c3,c4;

#if DEBUG > 0
  printf("calling column_prob_two_fixed with r=%d\n",r);
#endif

  /* case 1 : the character in the two fixed rows is mutated */

  /* w.l.o.g. we assume that c1 is the mutated character that
     occurs w.l.o.g in two rows. Since c1 cannot be the majority
     character it can occur at most r/2-2 times in the remaining rows */
  for(c1=0; c1<=r/2-2; c1++) 
    /* this is the highest numbers of c1s in the remaining rows without
       becoming the majority character */
    for(c2=0; c2<=r-2-c1; c2++) // c2 occurs at most r-2-c1 times and so on.
      for(c3=0; c3<=r-2-c1-c2; c3++) 
	for(c4=0; c4<=r-2-c1-c2-c3; c4++)
	  {
	    int max = TRUE;	
	    /* we check that c1 is not the majority character */
	    /* if this is the case we don not add the probability */
	    if( c2 >= c1+2 )
	      max = FALSE;
	    if( c3 >= c1+2 )
	      max = FALSE;
	    if( c4 >= c1+2 )
	      max = FALSE;
	    if( r-c1-c2-c3-c4 >= c1+2 )
	      max = FALSE;
	    if( max == FALSE )
	      prob1 += fournomial(seqErr,c1,c2,c3,c4,r-2);
	      
	  }

  /* for each of the four possible changed characters the probability
     of the two rows showing this character is (seqErr/4.0)^2 */
  prob1 = prob1*4.0*(seqErr/4.0)*(seqErr/4.0);


  /* case 2 : the character in the two fixed rows is  NOT mutated */

  for(c1=0; c1<=r-2; c1++) 
    for(c2=0; c2<=r-2-c1; c2++) // c2 occurs at most r-2-c1 times and so on.
      for(c3=0; c3<=r-2-c1-c2; c3++) 
	for(c4=0; c4<=r-2-c1-c2-c3; c4++)
	  {
	    int max = TRUE;	
	    int rest = r-2-c1-c2-c3-c4;
	    /* we check that the restt is not the majority character */
	    /* if this is the case we do not add the probability */
	    if( c1 >= rest+2 )
	      max = FALSE;
	    if( c2 >= rest+2 )
	      max = FALSE;
	    if( c3 >= rest+2 )
	      max = FALSE;
	    if( c4 >= rest+2 )
	      max = FALSE;
	    if( max == FALSE )
	      prob2 += fournomial(seqErr,c1,c2,c3,c4,r-2);
	      
	  }
 
  prob2 *= (1.0-seqErr)*(1.0-seqErr);
  // printf("%lf %lf\n",prob1,prob2);

  return prob1+prob2;

}


/* this function returns the expected number of "steps" which can be "saved"
   given the optimal partitioning of a column of data under the null hypothesis
   of a simple alignment */
double expected_savedSteps(int r, double seqErr)
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
	    a=max(c1-1,0);
	    c=max(c2-1,0);
	    g=max(c3-1,0);
	    t=max(c4-1,0);
	    d=max(c5-1,0);
	    saved=a+c+g+t+d-maxfive(a,c,g,t,d);

	    /* saved is the number of steps saved for the condition */
	    /* the probability of observing c1..c5 can be obtained by 
	       assuming that c5 is the true state without loss of generality */
	    /* So, expected is: */

	    Exp += saved*
	      fournomial(seqErr,c1,c2,c3,c4,r);

	    //AN ALTERNATIVE TO FOURNOMIAL THAT AARON THOUGHT SHOULD RUN
	    //FASTER BUT ACTUALLY SEEMS TO BE SLOWER:
	      //  pow(1-seqErr,c1)*pow(seqErr/4,r-c1)*
	      //  exp(lgamma(r+1)-lgamma(c1+1)-lgamma(c2+1)-lgamma(c3+1)-lgamma(c4+1)-lgamma(c5+1));

	  }

  ExpectedSavedSteps[r]=Exp;
  return Exp;
 }


//Give p-value associated with a chi-square statistic 
//with df degress of freedom
double chi_to_p(double chisq,int df)
{
  return(erfc(sqrt(chisq/(double)(df+1))));
}

double binomial_prob(int k,int n, double r){
  return(exp(lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))*pow(r,(double)k)*pow(1-r,(double)(n-k)));
}

double binomial_tail_prob(int n,int k,double p){
  double retval=1.0;
  int i;
  for(i=0;i<k;i++){
    retval-=binomial_prob(i,n,p);
  }
  return(retval);
}


//Function to give probability of alignment under null (simple) model
//based only on pairwise mismatch counts
double Pairwise_Only_Prob_Simple(Alignment_t *aln,int begincol,int endcol,Marker_t* m){
  int i,j,k,col;
  int Nmismatches_counts[20],currmm,currlen,ttlmm=0;
  double chisq=0,pval;
  double r;
  int pending,numsegs=0;
  int *lengths,availlengths=500,ttllen=0;
  double O,E,cum;
  int valid=0;


  /* Ensure that there is at least one pair with at least one column
     of overlap */
  for(i=0;i<aln->rows;i++){
    for(k=i+1;k<aln->rows;k++){
      for(j=begincol;j<endcol;j++){
	if(aln->ali[j][i]!=' '&&aln->ali[j][k]!=' '){
	  valid = 1;
	  break;
	}
      }
      if(valid)break;
    }
    if(valid)break;
  }
  
  //  printf("valid = %d\n",valid);
  if(valid==0)return(1.0);
  


  lengths=(int*)malloc(sizeof(int)*availlengths);
  if(lengths==NULL){
    fprintf(stderr,"Out of memory in Pairwise_Only_Prob_Simple()\n");
    exit(-1);
  }

  for(i=0;i<20;i++){
    Nmismatches_counts[i]=0;
  }


     
  for(i=0;i<aln->rows-1;i++){
    if((m->set)[i]==0)continue;
    for(j=i+1;j<aln->rows;j++){
      if((m->set)[j]==0)continue;
      currmm=pending=currlen=0;
      for(col=begincol;col<endcol;col++){
	if((aln->ali)[col][i]==' '||(aln->ali)[col][j]==' '){
	  if(pending==1){
	    if(currmm<=20-2){
	      Nmismatches_counts[currmm]++;
	    } else {
	      Nmismatches_counts[20-1]++;
	    }
	    while(numsegs>=availlengths){
	      availlengths*=2;
	      lengths=(int*)realloc(lengths,sizeof(int)*availlengths);
	      if(lengths==NULL){
		fprintf(stderr,"Out of memory in Pairwise_Only_Prob_Simple()\n");
		exit(-1);
	      }
	    }
	    ttllen+=currlen;
	    lengths[numsegs]=currlen;
	    pending=currlen=currmm=0;
	    numsegs++;
	  }
	} else {
	  pending=1;
	  currlen++;
	  if((aln->ali)[col][i]!=(aln->ali)[col][j]){
	    currmm++;
	    ttlmm++;
	  }
	}
      }
      if(pending==1){
	if(currmm<=20-2){
	  Nmismatches_counts[currmm]++;
	} else {
	  Nmismatches_counts[20-1]++;
	}

	while(numsegs>=availlengths){
	  availlengths*=2;
	  lengths=(int*)realloc(lengths,sizeof(int)*availlengths);
	  if(lengths==NULL){
	    fprintf(stderr,"Out of memory in Pairwise_Only_Prob_Simple()\n");
	    exit(-1);
	  }
	}
	lengths[numsegs]=currlen;
	ttllen+=currlen;
	pending=currlen=currmm=0;
	numsegs++;
      }
    }
  }

  if(ttlmm==0){
    free(lengths);
    return(1.0);
  }

  r=(double)ttlmm/(double)ttllen;
  cum=0;

  for(j=0;j<20-1;j++){
    O=(double)Nmismatches_counts[j];
    E=0;
    for(i=0;i<numsegs;i++){
      E+=binomial_prob(j,lengths[i],r);
    }
    if( E > 0.0 ){
      if(doPrintPWchisq)
	printf("%d: %d %e --> %e\n",j,(int)O,E,(O-E)*(O-E)/E);
      chisq+=(O-E)*(O-E)/E;
      cum+=E;
    }  else {
      if(doPrintPWchisq)
	printf("Non-positive E[%d]: %e ????; ttllen = %d, ttlmm = %d\n",j,E,ttllen,ttlmm);
    }
  }
  O=(double)Nmismatches_counts[20-1];
  E=numsegs-cum;
  if( E > 0.0 ){
    if(doPrintPWchisq)
      printf("%d: %d %e --> %e\n",j,(int)O,E,(O-E)*(O-E)/E);
    chisq+=(O-E)*(O-E)/E;
  } else {
    if(doPrintPWchisq)
      printf("Non-positive E[%d]: %e ????; ttllen = %d, ttlmm = %d\n",20-1,E,ttllen,ttlmm);
  }
  pval=chi_to_p(chisq,20-1);
  if(doPrintPWchisq)
    printf("Final chisq = %e (%d degrees of freedom; p-value = %e\n",chisq,20-1,pval);
  free(lengths);

  return(pval);

}


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
static char CM_ID[] = "$Id: ScoreREZ.c,v 1.2 2004-09-23 20:25:28 mcschatz Exp $";

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "UtilsREZ.h"
#include "AS_UTL_rand.h"
#include "ScoreREZ.h"

/* memory allocation for alignments */
/* note that we allocate first the columns, then the rows */
double **SeqErrArray;
double *MutErrArray;

Alignment_t *allocate_alignment(int c, int r)
{
  int i,j;
  Alignment_t *a = (Alignment_t*) safe_malloc(sizeof(Alignment_t));
  int **contributing = (int**) safe_calloc(sizeof(int*),MAX_X);
  for(i=0; i<MAX_X; i++)
    {
      contributing[i] = (int*) safe_calloc(sizeof(int),c);
      for(j=0; j<c; j++)
	contributing[i][j] = FALSE;
    }
  a->contributing = contributing;
  a->cols = c;
  a->rows = r;
  a->ali  = (char**) safe_calloc(sizeof(char*),c);
  for(i=0; i<c; i++)
    a->ali[i] = (char*) safe_calloc(sizeof(char),r);
  return a;
}



void free_alignment(Alignment_t* a)
{
  int i;
  for(i=0; i<MAX_X; i++)
    free(a->contributing[i]);
  free(a->contributing);

  for(i=0; i<a->cols; i++)
    free(a->ali[i]);
  free(a->ali);
  free(a);
}





static double fac(int n);


/* this function inspects an alignment and assuming all changes are sequencing
   error returns the median of the most likely explanation */
double guess_seqErr(Alignment_t *a)
{
  int i,j,k;
  int cols = a->cols;
  int rows = a->rows;

  int count[5]; // A C G T Dash
  double seqErr = 0.0;


  for(i=0; i<cols; i++){
    int max;
    for(k=0; k<5; k++)
      count[k] = 0;
    for(j=0; j<rows; j++)
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
	  }	
      }	
    max = 0;
    for(k=0; k<5; k++){
      if( count[max] < count[k] )
	max = k;
    }
    seqErr += ((double)(rows-count[max]))/((double) rows);
  }
  seqErr /= cols;
  return seqErr;
}


static char N[] = {'-','A','C','G','T' };

/* this function returns the probability of a sequencing
   error turning character at position (i,j) 
   given into character into.
   The sum of the conditional probabilities should be se.
   The probability of not changing the character is (1-se) */

double seq_err(char into, char given, int c, int r)
{	
  if( into == given )
    return (1.0-SeqErrArray[c][r]);
  else
    return 0.25*SeqErrArray[c][r];
}


/* This function returns the probability that the underlying character
   at position (i,j) is c */
static double prior(char c, int col){
  return 0.2;
}


int no_col_contributing(Alignment_t *a, int t)
{
  int i;
  int count = 0;
  for(i=0; i<a->cols; i++)
    if( a->contributing[t-2][i] == TRUE )
      count++;
  return count;
}


int col_contributing(Alignment_t *a, int c, int t)
{
  int j,k;
  int rows = a->rows;

  int count[5]; // A C G T Dash
  int observed=0;

  for(k=0; k<5; k++)
    count[k] = 0;
  for(j=0; j<rows; j++)
    {
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
	}	
    }	
  
  for(k=0; k<5; k++)
    if( count[k] >= t )
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


int no_col_contributing_two_fixed(Alignment_t *a, int r1, int r2)
{
  int i;
  int count = 0;
  for(i=0; i<a->cols; i++)
    if( col_contributing_two_fixed(a,i,r1,r2) ){
      //      printf("col %d is contributing\n",i);
      count++;
    }
  return count;
}


int col_contributing_two_fixed(Alignment_t *a, int c, int r1, int r2)
{
  int j,k,k1 = 4,k2 = 4,max;
  int rows = a->rows;

  int count[5]; // A C G T Dash


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
    }	
  
  if( k1 != k2 )
    return FALSE;


  for(k=0; k<5; k++)
    count[k] = 0;

  for(j=0; j<rows; j++)
    { 
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
	}	
    }	

  /* assume that the character in the assigned rows is the majority character */
  max = k1;
 
  for(k=0; k<5; k++)
    {
      if( k == k1 )
	continue;
      if( count[max] <= count[k] ) /* if any of the other characters has */
	max = k;                   /* as many occurences as k1 this is fine */
    }
      
  if( max == k1 )
    return FALSE;
  else
    return TRUE;
}







/* this function computes the log score of a colum */

double col_value(Partition_t* partition, Alignment_t *a, int col)
{ 
  int groups = partition->groups;
  int rows   = a->rows;
  double pValue[groups];
  double value = 1.0;
  int i,j,k;

  for(i=0; i<groups; i++){
    pValue[i] = 0.0;
  }

  for(i=0; i<groups; i++){  // Compute factor for each group number                                                // in the alignment
    for(j=0; j<5; j++){         // All possible underlying characters
      double prod=1.0;
      for(k=0; k<rows; k++){
	if( partition->part[k] == i )  // if the character in row k is in group i
	  prod *= seq_err(a->ali[col][k],N[j],col,k);
      }
      pValue[i] += prior(N[j],col)*prod; // otherwise we add the contributions
    }
    value *= pValue[i];  // in this column
  }
  return -log(value);
}





/* this function computes the log score of all contributing columns of the
   alignment */
double value(Partition_t *partition, Alignment_t *a, int ct){
  int i;
  int contCols=0;
  double value = 0.0;

  for(i=0; i<a->cols; i++)
    if( a->contributing[ct-2][i] == TRUE ){
      contCols++;
      value += col_value(partition,a,i);
    }
  
  if( contCols > 0 )
    return value/contCols;
  else
    return 0.0;
}







static double facREZ[300];

static double fac(int n)
{
  assert(n < 300);
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
  int i;
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

  for(i=0; i<sum; i++)
    prod *= q; 
  
  prod *= fac(n);
  prod /= fac(c1);
  prod /= fac(c2);
  prod /= fac(c3);
  prod /= fac(c4);
  prod /= fac(n-sum);

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
	    /* if this is the case we don not add the probability */
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

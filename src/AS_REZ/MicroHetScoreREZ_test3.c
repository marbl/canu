
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
static char CM_ID[] = "$Id: MicroHetScoreREZ_test3.c,v 1.4 2005-03-22 19:49:25 jason_miller Exp $";

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "UtilsREZ.h"
#include "AS_UTL_rand.h"
#include "MicroHetScoreREZ_test3.h"

#define DEBUG 0
extern double ExpectedSavedSteps[200];

static double AS_REZ_fac(int n);

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

#if 0
    for(i=s; i<e; i++)
      for(j=0; j<a->rows; j++)
	if( m->set[j] == TRUE )
	  if( a->ali[i][j] != ' '){
	    charcount++;
	    seqErr += a->seqErrArray[i][j];
	  }
    seqErr /= charcount;
#else
        
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
#endif
  return seqErr;
}

static char N[] = {'-','A','C','G','T' };



#define FACLIMIT 160

static double facREZ[FACLIMIT];

static double AS_REZ_fac(int n)
{
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
	      AS_REZ_fournomial(seqErr,c1,c2,c3,c4,r);

	  }

  ExpectedSavedSteps[r]=Exp;
  return Exp;
}


double AS_REZ_binomial_prob(int k,int n, double r){
  return(exp(lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))*pow(r,(double)k)*pow(1-r,(double)(n-k)));
}


/* this computes the probability of a binomial event */
double AS_REZ_binomial(int k,int n,double q)
{
  register int i;
  double p;
  
  double prod = 1;

  if(n-k > k){
    k=n-k;
    p=q;
    q=1-p;
  } else {
    p=1-q;
  }

  for(i=0; i<n-k; i++)
    prod *= p;

  for(i=0; i<k; i++)
    prod *= q; 

  for(i=k+1; i<=n; i++)
    prod *= i;

  if( k >= FACLIMIT )
    {
      for(i=k; i>=FACLIMIT; i--)
	prod /= i; 
      prod /= AS_REZ_fac(FACLIMIT-1);
    }
  else
    prod /= AS_REZ_fac(k);

  return prod;
}


double AS_REZ_binomial_tail_prob(int n,int k,double p){
  double retval=1.0;
  int i;
  for(i=0;i<k;i++){
    retval-=AS_REZ_binomial(i,n,p);
  }
  return(retval);
}

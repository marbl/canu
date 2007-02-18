
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
/**********************************************************************
 Module:      functions to to assign a score to overlaps
 Description:
 Assumptions:
**********************************************************************/

#include <assert.h>
#include <string.h>
#include <math.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_ALN_aligners.h"
#include "UtilsREZ.h"
#include "GreedyOverlapREZ.h"

#define GREEDYDEBUG 0

static FILE* LogFile=NULL;


int get_match_region(char* S1, char* Q1, char* S2, char* Q2, 
		     char* seq1, char* seq2, char* q1, char* q2,
		     int32 ia1, int32 ia2, OverlapMesg*  p);


/* this function returns the number of mismatches in the two null terminated
   strings. If the argument wdash is TRUE only substitutions are counted */
int  get_no_mismatches(char* S1, char* S2, char* Q1, char* Q2, 
		       int qt, int wdash)
{
  register int count=0,i=0;
  
  while(S1[i] != '\0')
    {
    if( S1[i] != S2[i] && (Q1[i]-'0' > qt) && (Q2[i]-'0' > qt) )
      {
#if GREEDYDEBUG > 0
	fprintf(LogFile,"(%c:%c:%d:%d) ",S1[i],S2[i],Q1[i]-'0',Q2[i]-'0'); 
#endif
	if( ! wdash )
	{
	  if( S1[i] != '-' && S2[i] != '-' )
	    count++;
	}
	else
	  count++;
      }
    i++;
    }
  return count;
}


/* this function returns the average mismatch probability
   given null terminated sequences and quality strings 
   If the argument wdash is TRUE, only substitutions are counted */
double get_average_mismatch_prob(char* S1, char* S2, char* Q1, char* Q2, 
				 int qt, int* cc, int wdash)
{
  double err1,err2,prob;
  double average_prob = 0.0;
  int i  = 0;
  int ct = 0;


  while( Q1[i] != '\0' )
    {
      if( Q1[i]-'0' > qt && Q2[i]-'0' > qt )
	{
	  err1 = pow(10,(Q1[i]-'0')/-10.0);
	  err2 = pow(10,(Q2[i]-'0')/-10.0);

	  prob = (1.0-err1)*(1.0-err2); 
	  // This is the probability that in both draws no error occured
      
	  prob += err1*err2/4.0;
	  // This is the probability that the underlying sequence was different
	  // the formula is 4*err1/4*err2/4;
      
	  if( wdash == FALSE )
	    {
	      if( S1[i] != '-' && S2[i] != '-' )
		{
		  average_prob += (1.0-prob);
		  // we add 1-p, where p is the probability 
		  // that there is a match 
		  ct++;
		}
	    }
	    else
	    {
	      average_prob += (1.0-prob);
	      // we add 1-p, where p is the probability 
	      // that there is a match 
	      ct++;
	    }  
	}
      i++;
    }

  if( ct > 0)
    average_prob /= ct;
  else
    average_prob = pow(10,-qt/10);

#if GREEDYDEBUG > 0
  fprintf(LogFile,"ct=%d\n",ct);
  fprintf(LogFile,"Average = %lf\n",average_prob);
#endif

  assert(average_prob <= 1.0);
  *cc = ct;

  return average_prob;
}


#define HIGH_QUALITY 10

/* this function returns a guess for the mutation error.
   It counts high quality mismatches and divides that number by
   the length of the alignment. IS NOT USED AT THE MOMENT. DID NOT SEEM
   TOO WORK VERY WELL */
double guess_mutation_prob(char* S1, char* S2, char* Q1, char* Q2, int* cc, int wdash)
{
  int i  = 0;
  int ct = 0;


  while( Q1[i] != '\0' )
    {
      if( Q1[i]-'0' > HIGH_QUALITY && Q2[i]-'0' > HIGH_QUALITY && S1[i] != S2[i] )
	{
	  ct++;
	}
      i++;
    }

  if( ct == 0)
    ct = 1;

  return (double)ct/(double)i;
}





/* This function computes a quality value for an overlap as follows:
 - First we infer the average probability that an alignment has a mismatch by using 
   the quality values of the matching region where the quality of a dash is arbitrarily
   set to the quality of a neighboring character.

 - Next we stipulate that the difference in between paralogous sequences is at least x 
   (x=MUTATIONERR) neglecting SNPs

-  We have now two modells : MS (the overlap is between fragments of the same location)
                           : MP (the overlap is between fragments of different locations)

-  If the sequences are from the same location in the genome we would expect the
   data to support modell MS better than if the data comes from different locations,
   meaning that P(MS|D) should be higher. Hence this probability can be taken as
   a quality measure. The variable we measure is the number of mismatches.
   Assuming modell MS we would expect DS=s*l mismatches on average 
   Assuming modell MP we would expect DP=(s+x)*l mismatches on average, where l
   is the length of the alignment. 

-  We modell the probability of seeing d mismatches under one of the two models
   by Poisson distributions with mean DS resp. DP.

-  The probability P(MS|D) can be computed using bayes' formula which results in
   P(MS|D) = P(D|MS)*P(MS) / (P(D|MS)*P(MP)+P(D|MP)*P(MP))
   Assuming that the priors are both 0.5 we have
   P(MS|D) = P(D|MS) / (P(D|MS)+P(D|MP))

-  filling in the formulas for the Poisson distributions
   we get P(MS|D) = 1 / 1 + ( (s+x)/s )^D* e^(-l*x) )

*/

#define MUTATION_ERR 0.04


int compute_bayesian_quality(InternalFragMesg* IFG1, InternalFragMesg* IFG2, 
			     OverlapMesg* olap, int qt, int *length, FILE* fp)
{
  int maxlength = 2*MIN(strlen(IFG1->sequence),strlen(IFG2->sequence))+2;
  {
  char  S1[maxlength];
  char  Q1[maxlength];
  char  S2[maxlength];
  char  Q2[maxlength];
  double s,x,p;
  int l,m;
  int cc; // number of contributing columns
  if( fp == NULL)
    LogFile = stderr;
  else
    LogFile = fp;


  // Get the matching region together with its quality values
  // this is in order to make the interface easier to change if
  // this function is not called with InternalFragMesgs

  /* if the first fragment is not the leftmost, change them */
  if( IFG1->iaccession == olap->bifrag)
    {
      InternalFragMesg* t;
      t = IFG1;
      IFG1 = IFG2;
      IFG2 = t;
    }
      

  /* we might have to complement one of the sequences */
  if(olap->orientation == AS_INNIE)
    Complement_Fragment_AS(IFG2);
  else 
    if(olap->orientation == AS_OUTTIE)
      Complement_Fragment_AS(IFG1);
  
  l = get_match_region(S1,Q1,S2,Q2,IFG1->sequence,IFG2->sequence,
		       IFG1->quality,IFG2->quality,
		       IFG1->iaccession,IFG2->iaccession,
		       olap);

  /* now we complement them back */

  if(olap->orientation == AS_INNIE)
    Complement_Fragment_AS(IFG2);
  else 
    if(olap->orientation == AS_OUTTIE)
      Complement_Fragment_AS(IFG1);

  // What is the average probability of a mismatch ?
  s = get_average_mismatch_prob(S1,S2,Q1,Q2,qt,&cc,TRUE);

  // we assume that the paralogous error rate (mutation rate)
  // is at least MUTATION_ERR
  x = MUTATION_ERR;
  //  x = 0.02+guess_mutation_prob(S1,S2,Q1,Q2,&cc,TRUE);;

  // count the number of mismatches
  m = get_no_mismatches(S1,S2,Q1,Q2,qt,TRUE);

  // compute the probability that this fits modell MS
  p  = 1.0/( 1.0+ (pow(((s+x)/s),m) / exp(cc*x) ));

#if GREEDYDEBUG > 0
  fprintf(LogFile,"\nAverage err        = %lf\n",s);
  fprintf(LogFile,"Guessed Mutation rate= %lf\n",x);
  fprintf(LogFile,"Number of mimatches  = %d\n",m);
  fprintf(LogFile,"Length               = %d\n",l);
  fprintf(LogFile,"Contributing columns = %d\n",cc);
  fprintf(LogFile,"quality              = %lf\n",1.0-p);
#endif

  olap->quality = 1.0-p;
  *length = l;

  return m;
  }
}



/* this function translates the special delta values into 'normal' ones using
   an array of ints instead of shorts */
int translate_delta(int* dn, signed char* d)
{
  int i,j;
  int dlen = strlen((char*)d);

  i = 0;
  j = 0;
  while( i < dlen )
    {
      if( d[i] == -127 )
	{
	  dn[j] = 126;
	  i++;
	  while( d[i] == -127 )
	    {
	      dn[j] += 126;
	      i++;
	    }
	  if( d[i] < 0 )
	    {
	      dn[j] *= -1;
	      dn[j] += d[i++];
	    }
	  else
	    dn[j] += d[i++];
	  j++;
	}
      else
	if( d[i] == -128 )
	  {
	    int k = abs(d[i]);
	    while(k >= 0)
	      if( d[i] > 0 )
		dn[j++] = 1;
	      else
		dn[j++] = -1;
	    i++;
	  }
      else
	dn[j++] = d[i++];
	
    }
  dn[j] = 0;
  return j;
}



/* This function fills into the preallocated arrays S1,Q1,S2,Q2
   the matching regions and the respective quality values
   of the match between seq1,seq2 that is given in p
 */

int get_match_region(char* S1, char* Q1, char* S2, char* Q2, 
		     char* seq1, char* seq2, char* q1, char* q2,
		     int32 ia1, int32 ia2, OverlapMesg*  p)
{
  int maxlength = 2*MIN(strlen(seq1),strlen(seq2))+2;
  {
  //  Display the alignment between strings  seq1  and  seq2  indicated
  //  in  (* p) .
  
  int   delta[2*maxlength];
  int  i, j; // indices in S and T
  int ks,kt,ns,nt,ct;
  char *S,*T,*SQ,*TQ;
  int slen,tlen;
  int dlen; 
  int k;
  char lastQ1 = 'I';
  char lastQ2 = 'I';

  S = seq1;
  T = seq2;
  SQ = q1;
  TQ = q2;

  k = translate_delta(delta,p->delta);

  slen = strlen(S);
  tlen = strlen(T);
  dlen = k;

  i = p->ahg;
  j = 0;

  ks = kt = 0;
  ns = nt = 1;
  
  if( p->bhg < 0)
    slen += p->bhg-1;
  if( p->bhg > 0)
    tlen -= p->bhg+1;


  while(i <= slen )
    {
      ct = 0;
      while(i <= slen)
	{
#if GREEDYDEBUG > 0
	  fprintf(stderr,"ct = %d\n", ct);
#endif
	  assert(ct < maxlength);
	  if(ks < dlen && ns == abs(delta[ks]))
	    {
	      if(delta[ks] < 0)
		{
		  S1[ct] = '-';
		  Q1[ct] = lastQ1;
		}
	      else
		{
		  S1[ct] = S[i];
		  Q1[ct] = SQ[i++];
		  lastQ1 = Q1[ct];
		}
	      ks++;
	      ns = 1;
	    }
	  else
	    {
              S1[ct] = S[i];
	      Q1[ct] = SQ[i++];
	      lastQ1 = Q1[ct];
	      ns ++; 
	    }
          ct ++;
	}
      S1[ct] = Q1[ct] = '\0';


      ct = 0;
      while(j <= tlen)
	{
	  if(kt < dlen && nt == abs(delta[kt]))
	    {
	      if(delta[kt] > 0)
		{
		  S2[ct] = '-';
		  Q2[ct] = lastQ2;
		  S2[ct] = '-';
		}
	      else
		{
		  S2[ct] = T[j];
		  Q2[ct] = TQ[j++];
		  lastQ2 = Q2[ct];
		}
	      kt++;
	      nt = 1;
	    }
	  else
	    {
              S2[ct] = T[j];
	      Q2[ct] = TQ[j++];
	      lastQ2 = Q2[ct];
              nt++;
	    }
	  ct++;
	}
      S2[ct] = Q2[ct] = '\0';


    }

  /* Now we assign the average of all mismatches that do not involve
     a dash */

  {
    double average1 = 0.0;
    double average2 = 0.0;
    int aq1;
    int aq2;

    int i  = 0;
    int ct1 = 0;
    int ct2 = 0;
    
    while(S1[i] != '\0')
      {
	if( S1[i] != '-' && S2[i] != '-' )
	  {
	    double err1 = pow(10,(Q1[i]-'0')/-10.0);
	    average1 += err1;
	    ct1++;
	  }
	if( S2[i] != '-' && S1[i] != '-' )
	  {
	    double err2 = pow(10,(Q2[i]-'0')/-10.0);
	    average2 += err2;
	    ct2++;
	  }	
	i++;
      }

    if( ct1 > 0)
      {
	average1 /= ct1;
#if GREEDYDEBUG > 0
	fprintf(LogFile,"average1=%lf\n",average1);
#endif	
    	aq1 = -10*log10(average1);

#if GREEDYDEBUG > 0
	fprintf(LogFile,"aq1=%d\n",aq1);
#endif	
	i=0;
	while( S1[i] != '\0')
	  {
	    if( S1[i] == '-')
	      Q1[i] = aq1+'0';
	    i++;
	  }
      }

    if( ct2 > 0)
      {
	average2 /= ct2;
#if GREEDYDEBUG > 0
	fprintf(LogFile,"average2=%lf\n",average2);
#endif	
	aq2 = -10*log10(average2);

#if GREEDYDEBUG > 0
	fprintf(LogFile,"aq2=%d\n",aq2);
#endif
	i=0;
	while( S2[i] != '\0')
	  {
	    if( S2[i] == '-')
	      Q2[i] = aq2+'0';
	    i++;
	  }
      } 
  }
#if GREEDYDEBUG > 0
      fprintf (LogFile,"\n%s\n%s\n",S1,Q1);
#endif
#if GREEDYDEBUG > 0
      fprintf (LogFile,"\n%s\n%s\n",S2,Q2);
#endif
      
  return(strlen(S1));
  }
}
  


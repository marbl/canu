
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
 Module:      AS_REZ
 Description: This file runs a number of trials to verify the theoretically
              computed distributions.
 Assumptions:
**********************************************************************/

static char CM_ID[] = "$Id: CheckDistributionREZ.c,v 1.3 2005-03-22 19:07:33 jason_miller Exp $";


#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "UtilsREZ.h"
#include "AS_UTL_rand.h"
#include "AS_UTL_skiplist.h"
#include "ScoreREZ.h"
#include "PartitionsREZ.h"



static int Trials = 100;
static int32 Cols    = 100;
static int32 Rows    = 4;
static int32 GenerateGroups = 1;
static double SeqErrMin = 0.0;
static double SeqErrMax = 0.0;
static double SeqErr=0.0;
static double MutErr  = 0.0;
static double MutErrMin = 0.0;
static double MutErrMax = 0.0;
static int ContributingColumns = 0;
static double MinVal;
static Alignment_t *Alignment;
static char* Seq;

static char N[] = {'-','A','C','G','T' };

static Alignment_t *generate_alignment(int *pmin, int* pmax,
				       int columns, int partitions);


static void print_alignment(Alignment_t *a, char* cor,  int w);

/* main routine */
int main (int argc, char *argv[]) {
  FILE *fp;
  int i,j;
  int rows;
  double seqErrGuessed = 0.0;
  char fileName[40];

  { /* Parse the argument list using "man 3 getopt". */
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "t:c:r:")) != EOF))
      switch(ch) 
	{
	case 't':
	  Trials = atoi(optarg);
	  break;
	case 'c':
	  Cols = atoi(optarg);
	  break;
	case 'r':
	  Rows = atoi(optarg);
	  break;
	case '?':
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	default :
	  errflg++;
	}
  }
 
  fprintf(stderr,"Running distribution checker for %d trials on alignments with %d columns\n",Trials,Cols);
  sprintf(fileName,"distribution-check.%d.%d",Trials,Cols);
  fp  = fopen(fileName,"w");
  
  fprintf(fp,"Number of trials %d, number of columns %d\n",Trials,Cols);
  fprintf(fp,"|  r |  SErr | SErrG |       2-FCC |       2-FCO |        2-CC |        2-CO |        3-CC |        3-CO |        4-CC |        4-CO |\n");
  fprintf(fp,"--------------------------------------------------------------------------------------------------------------------------------------\n");
  

  for(rows = 4; rows<=Rows; rows++)
    {
      int pmin = rows;
      int pmax = rows;
      //      fprintf(stderr,"Computing with : row %d\n ",rows);

      for(SeqErr = 0.005; SeqErr < 0.1; SeqErr += 0.005)
	{
	  double results2_CO[Trials];
	  double results3_CO[Trials];
	  double results4_CO[Trials];
	  double results2_FCO[Trials];
	  
	  double mean2_CO      = 0.0;
	  double variance2_CO  = 0.0;
	  double mean3_CO      = 0.0;
	  double variance3_CO  = 0.0;
	  double mean4_CO      = 0.0;
	  double variance4_CO  = 0.0;
	  double mean2_FCO     = 0.0;
	  double variance2_FCO = 0.0;

	  double p2_FCC;
	  double p2_CC;
	  double p3_CC;
	  double p4_CC;
	  
	  //  fprintf(stderr,"seqErr %lf \n",SeqErr);
     
	  GenerateGroups = FALSE; // We do not introduce mutation
	  SeqErrMin = SeqErrMax = SeqErr; 
	  MutErrMin = MutErrMax = MutErr = 0.0;
	  
	  for(j=0; j<Trials; j++)
	    {
	      /* generate an alignment with the specified 
		 sequencing error rate */
	      Alignment = generate_alignment(&pmin,&pmax,Cols,1);
	      //	      print_alignment(Alignment,Seq,100);

	      seqErrGuessed += guess_seqErr(Alignment);
	      
	      results2_FCO[j] = no_col_contributing_two_fixed(Alignment,0,1);
	      results2_CO[j]  = no_col_contributing(Alignment,2);
	      results3_CO[j]  = no_col_contributing(Alignment,3);
	      results4_CO[j]  = no_col_contributing(Alignment,4);


	      mean2_FCO += results2_FCO[j];
	      mean2_CO  += results2_CO[j];
	      mean3_CO  += results3_CO[j];
	      mean4_CO  += results4_CO[j];
	      
	      free_alignment(Alignment);
	      free(Seq);
	    }
	  mean2_FCO /= Trials;
	  mean2_CO  /= Trials;
	  mean3_CO  /= Trials;
	  mean4_CO  /= Trials;
	  seqErrGuessed /= Trials;
 
	  for(j=0;j<Trials;j++)
	    {
	      variance2_FCO += (mean2_FCO-results2_FCO[j])*(mean2_FCO-results2_FCO[j]);   
	      variance2_CO += (mean2_CO-results2_CO[j])*(mean2_CO-results2_CO[j]); 
	      variance3_CO += (mean3_CO-results3_CO[j])*(mean3_CO-results3_CO[j]); 
  	      variance4_CO += (mean4_CO-results4_CO[j])*(mean4_CO-results4_CO[j]); 
	    }
	  variance2_FCO /= Trials;
	  variance2_CO  /= Trials;
	  variance3_CO  /= Trials;
	  variance4_CO  /= Trials;    
	  
	  p2_FCC = column_prob_two_fixed(rows,seqErrGuessed);
	  p2_CC  = 1.0-column_prob(2,rows,seqErrGuessed);
	  p3_CC  = 1.0-column_prob(3,rows,seqErrGuessed);
	  p4_CC  = 1.0-column_prob(4,rows,seqErrGuessed);
	  //	  fprintf(stderr,"%lf %lf %lf %lf\n",p2_FCC,p2_CC,p3_CC,p4_CC);


	  fprintf(fp,"| %2d | %3.3lf | %3.3lf | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | \n",
		  rows,SeqErr,seqErrGuessed,
		  mean2_FCO,sqrt(variance2_FCO),p2_FCC*Cols,sqrt(p2_FCC*(1.0-p2_FCC)*Cols),
		  mean2_CO,sqrt(variance2_CO),p2_CC*Cols,sqrt(p2_CC*(1.0-p2_CC)*Cols),
		  mean3_CO,sqrt(variance3_CO),p3_CC*Cols,sqrt(p3_CC*(1.0-p3_CC)*Cols),
		  mean4_CO,sqrt(variance4_CO),p4_CC*Cols,sqrt(p4_CC*(1.0-p4_CC)*Cols)
		  );

	  fprintf(stderr,"| %2d | %3.3lf | %3.3lf | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | (%3.2lf,%3.2lf) | \n",
		  rows,SeqErr,seqErrGuessed,
		  mean2_FCO,sqrt(variance2_FCO),p2_FCC*Cols,sqrt(p2_FCC*(1.0-p2_FCC)*Cols),
		  mean2_CO,sqrt(variance2_CO),p2_CC*Cols,sqrt(p2_CC*(1.0-p2_CC)*Cols),
		  mean3_CO,sqrt(variance3_CO),p3_CC*Cols,sqrt(p3_CC*(1.0-p3_CC)*Cols),
		  mean4_CO,sqrt(variance4_CO),p4_CC*Cols,sqrt(p4_CC*(1.0-p4_CC)*Cols)
		  );

	  fflush(stderr);   
	  fflush(fp);
	}
    }
  fclose(fp);
}






static void print_alignment(Alignment_t *a, char* cor,  int w){
  int i,j,l;
  int c = a->cols;
  int r = a->rows;
  int iter = 0;

  do{
    int p = w*iter+w;
    l = ( p < c ? p : c );

    for(i=iter*w; i<l; i++)
      printf("%c",cor[i]);
    printf(" %d\n",l);
    
    for(i=0; i<r; i++){
      for(j=iter*w; j<l; j++)
	if( a->ali[j][i] != cor[j] )
	  printf("%c",a->ali[j][i]);
	else
	  printf("%c",'.');
      printf(" %d\n",l);
    }

    iter++;
    printf("\n");
  }while(iter*w < c);
  
}


 
Alignment_t *generate_alignment(int32 *pmin, int32 *pmax, 
				int32 cols, int32 partitions)
{
  int i,j;
  Alignment_t *a;
  Alignment_t *mut;
  int rows = 0;
  int r;
  int groups[partitions];

  int seed = getpid();
  /* set random seed */
  // srand48(seed);
  

  Seq = (char*) safe_calloc(sizeof(char),cols);
  /* Generate a random sequence */
  for(i=0; i<cols; i++)
    Seq[i] = N[GetRand_AS(0,4,TRUE)];


  /* throw a dice to determine how many rows are in each group */
  //  printf("\n");
  for(i=0; i<partitions; i++){
    groups[i] = GetRand_AS(pmin[i],pmax[i],TRUE);
    rows += groups[i];
  }

  /* Allocate and the sequencing error matrix */  


  if( SeqErrArray != NULL)
    {
      for(i=0; i<cols; i++)
	free(SeqErrArray[i]);
      free(SeqErrArray);
    }

  SeqErrArray = (double**) safe_calloc(sizeof(double*),cols);
  for(i=0; i<cols; i++)
    SeqErrArray[i] = (double*) safe_calloc(sizeof(double),Rows);

  for(i=0; i<cols; i++)
    for(j=0; j<Rows; j++)
      SeqErrArray[i][j] = GetDrand_AS(SeqErrMin,SeqErrMax); 
    
  /* Allocate and the mutation error array */  

  if( MutErrArray != NULL)
    free(MutErrArray);

  MutErrArray = (double*) safe_calloc(sizeof(double),partitions);

  for(i=0; i<partitions; i++)
    MutErrArray[i] = GetDrand_AS(MutErrMin,MutErrMax); 
    

  /* generate the alignment and the group sequences */
  a   = allocate_alignment(cols,rows);
  mut = allocate_alignment(cols,partitions);
  
  for(i=0; i<partitions; i++)
    for(j=0; j<cols; j++){
      double rv = drand48();
      
      /* Skip the mutation error */
      if( GenerateGroups == FALSE )
	rv = 1.0;
      
      if( MutErrArray[i] > rv ){
	char m;
	mut->ali[j][i] = m = N[GetRand_AS(0,4,TRUE)];
	while( m == Seq[j] )
	  mut->ali[j][i] = m = N[GetRand_AS(0,4,TRUE)];
      }
      else
	mut->ali[j][i] = Seq[j];
    }

  r = 0;
  for(i=0; i<partitions; i++)
    {
      int k;
      for(k=0; k<groups[i]; k++)
	{
	  for(j=0; j<cols; j++)
	    {
	      double rv = drand48();

	      if( SeqErrArray[j][r] > rv ){
		char s;
		a->ali[j][r] = s = N[GetRand_AS(0,4,TRUE)];
		while( s == mut->ali[j][i] )
		  a->ali[j][r] = s = N[GetRand_AS(0,4,TRUE)];
	      }
	      else
		a->ali[j][r] = mut->ali[j][i];
	    }
	  r++;
	}
    }
  
  a->hSeqErr = guess_seqErr(a);
  for(i=0; i<MAX_X; i++)
    for(j=0; j<cols; j++)
      a->contributing[i][j] = col_contributing(a,j,i+2);

  free_alignment(mut);

  return a;
}




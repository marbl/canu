
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
 Module:      The functions in this file are only for testing purposes
 Description:
 Assumptions:
**********************************************************************/

static char CM_ID[] = "$Id: MicroHetSimulatorREZ.c,v 1.3 2005-03-22 19:07:57 jason_miller Exp $";

//TO CONDUCT TESTS OF SIMPLE ONLY (NO ATTEMPT TO RESOLVE), DEFINE TESTSIMPLEONLY
#define TESTSIMPLEONLY

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
#include "MicroHetREZ.h"
#include "MicroHetScoreREZ.h"
#include "MicroHetPartitionsREZ.h"

/* The following global variables contain data to generate
   an alignment with a predefined mutation rate (may be an interval)
   and a predefined mutation rate */
static int32 Partitions     = 0; // the number of partitions in the 
static int32 Cols           = 0; // the number of columns
static int32 Rows           = 0; // the number of rows
static int32 GenerateGroups = TRUE; // if this is FALSE, no mutation is introduced
static double SeqErr        = 0.0;
static double SeqErrMin     = 0.0;
static double SeqErrMax     = 0.0;
static double MutErr        = 0.0;
static double MutErrMin     = 0.0;
static double TestThresh    = 1e-7;
static double MutErrMax     = 0.0;
static int *PMinRows;           // these fields hold the maximum respectively
static int *PMaxRows;           // resp. minimum number of elements in a group   
static Partition_t *CorrectPartition;
static Partition_t *MinPartition;
static Partition_t *SinglePartition;
static double *MutErrArray;
static double SingleVal;
static double CorrectVal;
static double MinVal;
static Alignment_t *Alignment;
static char* Seq;               // holds a randomly generated sequence
static uint64 NoOfPartitions = 0;
static int ContributingThreshold = 2; 

static char N[] = {'-','A','C','G','T' };

SL_DEF(Partition_t)
SL_TYPE(Partition_t) *guess_partition(Alignment_t *ali, int ct);

static Alignment_t *generate_alignment(int *pmin, int* pmax, int columns, 
				       int partitions);

static void read_conf(char* fileName);

static void stirling(Partition_t *set, int k, int g, Marker_t* m);
static void check_all_partitions(void);


/* main routine */
int main (int argc, char *argv[]) {

#ifdef RETURNPVALS
  double pval;
#endif

  /*Variable definitions for parsing the commandline */
  char  *fileName  = NULL;
  
  { /* Parse the argument list using "man 3 getopt". */
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "s:m:f:t:")) != EOF))
      switch(ch) 
	{
	case 's':
	  SeqErr = SeqErrMin = SeqErrMax = atof(optarg);
	  break;
	case 't':
	  TestThresh = atof(optarg);
	  break;
	case 'f':
	  fileName = strdup(optarg);
	  break;
	case 'm':
	  MutErr = MutErrMin = MutErrMax = atof(optarg);
	  break;
	case '?':
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	default :
	  errflg++;
	}
  }

  /* read the config file to get values for options 
     not set by the commandline */
  assert(fileName != NULL);
  read_conf(fileName);

  /* generate an alignment with the specified mutation and 
     sequencing error rate */
  Alignment = generate_alignment(PMinRows,PMaxRows,Cols,
				 Partitions);

  print_alignment(Alignment,100);

  //AARON'S STUFF
  { 
    double erate;
    MPSTAT mpresult;
    Marker_t    *all = allocate_marker(Alignment->rows);
    count_columns(Alignment,all);
    erate=guess_seqErr(Alignment,all,0,Alignment->cols-1);
    mpresult=MP_score_alignment(Alignment,erate,0,Alignment->cols-1);
    printf("\n\n>>>>> Aaron's MP evaluation: O = %d E = %f pr = %e",mpresult.Obs, mpresult.Exp,mpresult.pr);
    if(mpresult.pr<=TestThresh){
      printf(" : NOT SIMPLE\n\n");
    } else {
      printf(" : SIMPLE\n\n");
    }
    free_marker(all);
  }
  //END CHUNK OF AARON'S STUFF


  printf("Conducting 2 fixed test\n"); 
  {
    /* the 2 fixed contributing test */  
    UnitigStatus_t stat;
    int crit;
    Partition_t *p   = allocate_partition(Alignment->rows);
    Marker_t    *all = allocate_marker(Alignment->rows);

    printf("\nTesting with error probability %e \n",TestThresh);
    stat = test_simple(Alignment,TestThresh,all,0,Alignment->cols,&crit
#ifdef RETURNPVALS
,&pval
#endif
);

    switch(stat)
      {
      case UNITIG_IS_UNKNOWN:
	printf("Alignment is NOT SIMPLE \n");
#ifndef TESTSIMPLEONLY
	bipartition(Alignment,all,p,0,Alignment->cols,TestThresh,1);
	print_part(p);
#endif
	free_marker(all);
	free_partition(p);
	break;
      case UNITIG_IS_SIMPLE:
	printf("Alignment is SIMPLE \n");
	break;
      case UNITIG_IS_REPETITIVE:
	printf("Alignment is NOT SIMPLE \n");
#ifndef TESTSIMPLEONLY
	bipartition(Alignment,all,p,0,Alignment->cols,TestThresh,1);
	print_part(p);
#endif
	free_marker(all);
	free_partition(p);
	break;
      case UNITIG_IS_SHALLOW:
	printf("Alignment is SHALLOW \n");
	break;
	
      }
  }
  
#ifndef TESTSIMPLEONLY  
  //  MinPartition = allocate_partition(Rows);

  /* Now we guess a partition using all possible predicates */
  /*
  for(ContributingThreshold = 2; ContributingThreshold<=4; ContributingThreshold++)
    {
      printf("\nGuess partition with threshold %d\n",ContributingThreshold);
      sl = guess_partition(Alignment,ContributingThreshold);
      it = MinSL_Partition_t(sl);

      for(i=0; i<Rows; i++)
	MinPartition->part[i] = CorrectPartition->part[i];

      if( GenerateGroups == TRUE )
	{
	  CorrectVal = value(CorrectPartition,Alignment,ContributingThreshold);
	  printf("\nCorrect partition with value %3.3lf :\n",CorrectVal); 
	  print_part(CorrectPartition);
	  MinVal = CorrectVal;
	}
      else
	{
	  SinglePartition = allocate_partition(Rows);
	  for(i=0; i<Rows; i++)
	    SinglePartition->part[i] = 0;
	  SingleVal = value(SinglePartition,Alignment,ContributingThreshold);
	  printf("\nCorrect partition with value %3.3lf :\n",SingleVal); 
	  print_part(SinglePartition);
	  free_partition(SinglePartition);
	}

      printf("\nGuessed partition with score %3.3lf \n",it->key); 
      print_part(it->value);
      FreeSL_Partition_t(sl);

      if( Partitions < 4){
	check_all_partitions();
	printf("Best partition : \n"); 
	print_part(MinPartition);
	printf("Value of best partition = %3.3lf \n",MinVal);
      }
 
    }

  free_partition(CorrectPartition);
  free_partition(MinPartition);
  free(PMinRows);
  free(PMaxRows);
  */
#endif
  free_alignment(Alignment);
  return 0;
}


/* static functions */

static void read_conf(char* fileName)
{
  char buffer[100];

  FILE *confFile = file_open(fileName,"r");
  while( NULL != fgets(buffer,100,confFile))
    {
      switch(buffer[0]){
      case '#' :
	continue;
	break;
      case 'p':
	if( buffer[1] == ':')
	  {
	    sscanf(buffer,"p:%d\n",&Partitions);
	    printf("Number of partitions %d\n",Partitions);
	    PMinRows = (int*) safe_malloc(sizeof(int)*Partitions);
	    PMaxRows = (int*) safe_malloc(sizeof(int)*Partitions);
	  }
	else
	  {
	    int p;
	    int dummy;
	    sscanf(buffer,"p%d:",&p);
	    sscanf(buffer,"p%d:%d,%d\n",&dummy,&PMinRows[p],&PMaxRows[p]);
	  }
	break;
      case 'g':
	{
	  int g;
	  sscanf(buffer,"g:%d",&g);
	  if( g == 0 )
	    GenerateGroups = TRUE;
	  else
	    GenerateGroups = FALSE;
	}
	break;
      case 'c':
	sscanf(buffer,"c:%d\n",&Cols);
	printf("Number of cols       %d\n",Cols);
	break;
      case 's':
	if( SeqErr == 0.0 )
	  sscanf(buffer,"s:%lf:%lf\n",&SeqErrMin,&SeqErrMax);
	printf("Sequencing error between %0.3f and %0.3f\n",SeqErrMin,SeqErrMax);
	SeqErr = SeqErrMin;
	break;
      case 'm':
	if( MutErr == 0.0 )
	  sscanf(buffer,"m:%lf:%lf\n",&MutErrMin,&MutErrMax);
	printf("Mutation error between %0.3f and %0.3f\n",MutErrMin,MutErrMax);
	MutErr = MutErrMin;
	break;
      }
    }

}


 
Alignment_t *generate_alignment(int32 *pmin, int32 *pmax, int32 cols,
				int32 partitions)
{
  int i,j;
  Alignment_t *a;
  Alignment_t *mut;
  int rows = 0;
  int groups[partitions];

  /* set random seed */
  // int seed = getpid();
  // srand48(seed);
  

  Seq = (char*) safe_calloc(sizeof(char),cols);
  /* Generate a random sequence */
  for(i=0; i<cols; i++)
    Seq[i] = N[GetRand_AS(0,4,TRUE)];


  /* throw a dice to determine how many rows are in each group */
  //  printf("\n");
  for(i=0; i<partitions; i++){
    groups[i] = GetRand_AS(pmin[i],pmax[i],TRUE);
    //  printf("Partition %d has %d rows\n",i,groups[i]);
    rows += groups[i];
  }

  Rows = rows;

  /* generate the alignment and the group sequences */
  a   = allocate_alignment(cols,Rows);
  mut = allocate_alignment(cols,Partitions);

  /* Allocate the sequencing error matrix */  

  //  a->segments = (Test_Segment_t *) safe_malloc(sizeof(Test_Segment_t));  
  //  a->segments[0].start = 0;
  //  a->segments[0].end   = cols;
  //  a->noOfSegs = 1;
  for(i=0; i<cols; i++)
    for(j=0; j<Rows; j++)
      a->seqErrArray[i][j] = GetDrand_AS(SeqErrMin,SeqErrMax);
  

  a->hasQuality = TRUE;

  /* Allocate and the mutation error array */  
  MutErrArray = (double*) safe_calloc(sizeof(double),partitions);
  
  for(i=0; i<partitions; i++)
    MutErrArray[i] = GetDrand_AS(MutErrMin,MutErrMax); 
    
  CorrectPartition = allocate_partition(Rows);
  rows = 0;

  /* Now we know the correct partition */
  for(i=0; i<partitions; i++)
    for(j=0; j<groups[i]; j++)
      CorrectPartition->part[rows++] = i;

  CorrectPartition->groups = partitions;
  CorrectPartition->len    = Rows;

  
  for(i=0; i<partitions; i++)
    for(j=0; j<cols; j++){
      int rand;
      if( MutErrArray[i] > 0.0 )
	rand = GetRand_AS(0,(int)(1/MutErrArray[i]),TRUE);
      else
	rand = 1;
      /* Skip the mutation error */
      if( GenerateGroups == FALSE )
	rand = 1;
      
      if( rand == 0 ){
	char m;
	mut->ali[j][i] = m = N[GetRand_AS(0,4,TRUE)];
	while( m == Seq[j] )
	  mut->ali[j][i] = m = N[GetRand_AS(0,4,TRUE)];
      }
      else
	mut->ali[j][i] = Seq[j];
    }

  rows = 0;
  for(i=0; i<partitions; i++)
    {
      int k;
      for(k=0; k<groups[i]; k++)
	{
	  for(j=0; j<cols; j++)
	    {
	      int rand;
	      if( a->seqErrArray[j][rows] > 0.0 )
		rand = GetRand_AS(0,(int)(1/a->seqErrArray[j][rows]),TRUE);
	      else
		rand = 1;
	      if( rand == 0 ){
		char s;
		a->ali[j][rows] = s = N[GetRand_AS(0,4,TRUE)];
		while( s == mut->ali[j][i] )
		  a->ali[j][rows] = s = N[GetRand_AS(0,4,TRUE)];
	      }
	      else
		a->ali[j][rows] = mut->ali[j][i];
	    }
	  rows++;
	}
    }
  //  a->hSeqErr  = (double*) safe_calloc(sizeof(double),a->noOfSegs);
  //  for(i=0; i<a->noOfSegs; i++)
  //    a->hSeqErr[i] =  GetDrand_AS(SeqErrMin,SeqErrMax); 

  return a;
}




static void check_all_partitions(void){
  int i;
  Partition_t* start = allocate_partition(Rows);
  Marker_t*        m = allocate_marker(Rows);

  NoOfPartitions = 0;
  for(i=1; i<=Rows; i++){
    uint64 on = NoOfPartitions;
    start->len    = Rows;
    start->groups = i;
    stirling(start,Rows,i,m);
    printf("For (%d,%d) I checked " F_U64 " partitions\n",
           Rows,i,NoOfPartitions-on);

  }
  printf("Checked overall " F_U64 " partitions\n",NoOfPartitions);
  free_partition(start);
}


/* k is the number of elements that are left of inititally len
   which have to be placed in the remaining g of initially par groups */
static void stirling(Partition_t *set, int k, int g, Marker_t* m){
  int i;

  /* if we have only k==g elements they all come into
     a group of their own */
  if( k == g ){
    double val;
    
    for(i=0; i<k; i++)
      set->part[i] = i;
    
    val = value(set,Alignment,ContributingThreshold,m);

    if( val < CorrectVal ){
      //      print_part(set);
      //      printf("Partition has smaller value %lf than correct %lf \n",val,CorrectVal);
      if( val <= MinVal ){
	MinVal = val;
	for(i=0; i<set->len; i++)
	  MinPartition->part[i] = set->part[i];
	MinPartition->len    = set->len;
	MinPartition->groups = set->groups;
      }      
    }

    NoOfPartitions++;
    return;
  }
  
  /* if we have only one group, all elements come into that group */
  if( g == 1 ){
    double val;
    
    for(i=0; i<k; i++)
      set->part[i] = 0;

    val = value(set,Alignment,ContributingThreshold,m);

    if( val < CorrectVal ){
      //      print_part(set);
      //      printf("Partition has smaller value %lf than correct %lf \n",val,CorrectVal);
      if( val <= MinVal ){
	MinVal = val;
	for(i=0; i<set->len; i++)
	  MinPartition->part[i] = set->part[i];
	MinPartition->len    = set->len;
	MinPartition->groups = set->groups;
      }      
    }

    NoOfPartitions++;
    return;
  }
  
  /* part 1: the last element is in a group of its own */
  set->part[k-1] = g-1;
  stirling(set,k-1,g-1,m);
  
  /* part 2: the last element goes in one of the g groups 
     formed by the k-1 other elements */
  for(i=0; i<g; i++){
    set->part[k-1] = i;
    stirling(set,k-1,g,m);
  }
  return;
}


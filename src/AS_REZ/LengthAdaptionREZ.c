
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
static char CM_ID[] = "$Id: LengthAdaptionREZ.c,v 1.5 2005-09-15 15:20:16 eliv Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "GreedyOverlapREZ.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_ALN_aligners.h"
#include "UtilsREZ.h"


#define REZ_DP_ERATE  .1
#define REZ_DP_THRESH 1e-6
#define REZ_DP_MINLEN 30
#define DEBUG



int main(int argc, char *argv[]) 
{
  int i,j,m;
  char *fileName;
  /*Variable defintions for parsing the commandline */
  GenericMesg *pmesg    = NULL; 
  FILE        *fileInput= NULL;
  FILE **quotientFile;
  int **correct;
  int **false;

  int intQuality = 0;
  MesgReader   readerFn = NULL;
  ReadStructp  input;
  int noOfFrags = 0;
  int count;
  int length;
  int bucketLength=100;
  int probBucketLength=1;
  int olapLength;

  int qualityThreshold = 0;
  GenericMesg  **fragMesg;

  { /* Parse the argument list using "man 3 getopt". */
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "f:l:b:t:")) != EOF))
      switch(ch) 
	{
	case 'f':
	  fileName = strdup(optarg);
	  break;
	case 'l':
	  bucketLength = atoi(optarg);
	  break;
	case 'b':
	  probBucketLength = atoi(optarg);
	  break;
	case 't':
	  qualityThreshold = atoi(optarg);
	  break;
	case '?':
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	  break;
	}
  }
   
  if( (argc - optind != 0 )){
    fprintf (stderr,"USAGE:  length-adaption -f .inp file [-t quality threshold] [-l length bucket size] [-b prob bucket size * 100]\n");
    exit(1);
  }
  

  /* Test preconditions */
  assert( fileName != NULL);
  assert( strlen(fileName) != 0 );

  /* allocate the arrays */
  quotientFile = (FILE**) malloc(sizeof(FILE*)*600/bucketLength);
  correct      = (int**)  malloc(sizeof(int*)*600/bucketLength);
  false        = (int**)  malloc(sizeof(int*)*600/bucketLength);
  for(i=0; i<600; i+=bucketLength)
    {
      correct[i/bucketLength] = (int*) malloc(sizeof(int)*100/probBucketLength);
      false[i/bucketLength] =   (int*) malloc(sizeof(int)*100/probBucketLength);
    }
  for(i=0; i<600; i+=bucketLength)
     for(j=0; j<100; j+=probBucketLength)
       { 
	 correct[i/bucketLength][j/probBucketLength] = 0;
	 false[i/bucketLength][j/probBucketLength] = 0;
       }

  /* set up quotient files */ 
  for(olapLength = 0; olapLength < 600; olapLength += bucketLength)
      {
	char fileName[100];
	 
	sprintf(fileName,"quotient.%d.%d",olapLength,olapLength+bucketLength);
	quotientFile[olapLength/bucketLength] = file_open(fileName,"w");  
      }


  fileInput = file_open(fileName,"r");  
  readerFn  = (MesgReader)InputFileType_AS(fileInput);
  
  /* count the number of fragments */
  while( readerFn(fileInput,&pmesg) != EOF ) 
    {
      MessageType mesgtype = pmesg->t;
      if( mesgtype == MESG_IFG ) 
	noOfFrags++;
    }

  /* allocate the fragment array */
  fragMesg = (GenericMesg**) safe_calloc(sizeof(GenericMesg*),noOfFrags);

  fclose(fileInput);
  fileInput = file_open(fileName,"r");  

  /* read the fragment array */
  count = 0;
  while( readerFn(fileInput,&pmesg) != EOF ) 
    {
      MessageType mesgtype = pmesg->t;
      switch(mesgtype){
      case MESG_IFG : 
	fragMesg[count++] = DuplicateProtoMesg_AS(pmesg);
	break;
      }
    }

  /* we fill the lookup tables */
  //  MismatchTableREZ = fill_QV_Mismatch_table(prob_mismatch);
  //  MatchTableREZ    = fill_QV_Match_table(prob_match);

  
  /* we align all pairs of fragments and compute their quality */
  fprintf(stderr,"** number of fragments %d\n",noOfFrags);

  for(i=0; i<noOfFrags-1; i++)
    for(j=i+1; j<noOfFrags; j++)
      {
	OverlapMesg *olapp;
	OverlapMesg *olap_unit    = (OverlapMesg*) safe_malloc(sizeof(OverlapMesg));
	OverlapMesg *olap_quality = (OverlapMesg*) safe_malloc(sizeof(OverlapMesg));
	
	GenericMesg from,*to;
	int where;
	
	InternalFragMesg* IFG1 = (InternalFragMesg*) fragMesg[i]->m;
	InternalFragMesg* IFG2 = (InternalFragMesg*) fragMesg[j]->m;
	
	olapp = DP_Compare_AS(IFG1, IFG2,
			      -strlen(IFG2->sequence)+40,strlen(IFG1->sequence)-40,0,
			      REZ_DP_ERATE,REZ_DP_THRESH,REZ_DP_MINLEN, 
			      AS_FIND_ALIGN,&where);

	/* we copy the result for later use ? */
	if(olapp != NULL)
	  {
	    olap_unit->aifrag       = olapp->aifrag;
	    olap_unit->bifrag       = olapp->bifrag;
	    olap_unit->ahg          = olapp->ahg;
	    olap_unit->bhg          = olapp->bhg;
	    olap_unit->orientation  = olapp->orientation;
	    olap_unit->overlap_type = olapp->overlap_type;
	    olap_unit->quality      = olapp->quality;
	    olap_unit->min_offset   = olapp->min_offset;
	    olap_unit->max_offset   = olapp->max_offset;
	    olap_unit->polymorph_ct = olapp->polymorph_ct;
	    olap_unit->delta = (signed char*) safe_calloc(sizeof(signed char),strlen((char*)olapp->delta)+1);
	    strcpy((char*)olap_unit->delta,(char*)olapp->delta);
	    
	    /* Debugging : print the overlap and compute its quality */
	    
#ifdef DEBUG
	    fprintf(stdout,"DP orientation = %d\n",olap_unit->orientation);
	    fprintf(stdout,"afrag = %d\n",olap_unit->aifrag);
	    fprintf(stdout,"bfrag = %d\n",olap_unit->bifrag);
	    fprintf(stdout,"ahang = %d\n",olap_unit->ahg);
	    fprintf(stdout,"bhang = %d\n",olap_unit->bhg);
	    
	    Print_Overlap_AS(stdout, IFG1, IFG2, olap_unit);
#endif
	    
	    m  = compute_bayesian_quality(IFG1,IFG2,olap_unit,qualityThreshold,&length,NULL);
	    
	    if( length > 0)
	      {
		if( repetitive_overlap_sim(IFG1,IFG2,olap_unit) != TRUE )
		  false[length/bucketLength][(int)(olap_unit->quality*100)/probBucketLength]++;
		else
		  correct[length/bucketLength][(int)(olap_unit->quality*100)/probBucketLength]++;
	      }
	  }
	    


	/* no the other direction */
	olapp = DP_Compare_AS(IFG1, IFG2,
			      -strlen(IFG2->sequence)+40,strlen(IFG1->sequence)-40,1,
			      REZ_DP_ERATE,REZ_DP_THRESH,REZ_DP_MINLEN, 
			      AS_FIND_ALIGN,&where);

	
	if(olapp != NULL)
	  {
	    olap_unit->aifrag       = olapp->aifrag;
	    olap_unit->bifrag       = olapp->bifrag;
	    olap_unit->ahg          = olapp->ahg;
	    olap_unit->bhg          = olapp->bhg;
	    olap_unit->orientation  = olapp->orientation;
	    olap_unit->overlap_type = olapp->overlap_type;
	    olap_unit->quality      = olapp->quality;
	    olap_unit->min_offset   = olapp->min_offset;
	    olap_unit->max_offset   = olapp->max_offset;
	    olap_unit->polymorph_ct = olapp->polymorph_ct;
	    olap_unit->delta = (signed char*) safe_calloc(sizeof(signed char),strlen((char*)olapp->delta)+1);
	    strcpy((char*)olap_unit->delta,(char*)olapp->delta);
#ifdef DEBUG
	    fprintf(stdout,"DP orientation = %d\n",olap_unit->orientation);
	    fprintf(stdout,"afrag = %d\n",olap_unit->aifrag);
	    fprintf(stdout,"bfrag = %d\n",olap_unit->bifrag);
	    fprintf(stdout,"ahang = %d\n",olap_unit->ahg);
	    fprintf(stdout,"bhang = %d\n",olap_unit->bhg);
	    Print_Overlap_AS(stdout, IFG1, IFG2, olap_unit);
#endif
	    m  = compute_bayesian_quality(IFG1,IFG2,olap_unit,qualityThreshold,&length,NULL);
	    
	    if( length > 0)
	      {
		if( repetitive_overlap_sim(IFG1,IFG2,olap_unit) != TRUE )
		  false[length/bucketLength][(int)(olap_unit->quality*100)/probBucketLength]++;
		else
		  {
		    correct[length/bucketLength][(int)(olap_unit->quality*100)/probBucketLength]++;
		  }
	      }
	    
	    
	  }
      }


  for(i=0; i<600; i+=bucketLength)
    {
      for(j=1; j<100; j+=probBucketLength)
	{
	  int c1 = correct[i/bucketLength][j];
	  int f1 = false[i/bucketLength][j];
	  int c2 = correct[i/bucketLength][j-1];
	  int f2 = false[i/bucketLength][j-1];
	  float q1,q2;
	  if( c1+f1 > 0 )
	    q1 = (float)c1/(c1+f1);
	  else
	    q1 = 1.0;
 
	  if( c2+f2 > 0 )
	    q2 = (float)c2/(c2+f2);
	  else
	    q2 = 1.0;

	  fprintf(quotientFile[i/bucketLength],"%d %f \n",j,q1/q2);

	}
      //      fprintf(quotientFile[i/bucketLength],"\n");
      fclose(quotientFile[i/bucketLength]);
    }

  /* in the end we free all fragments */
  for(i=0; i<noOfFrags; i++)
    FreeProtoMesg_AS(fragMesg[i]);
 

  for(i=0; i<600; i+=bucketLength)
    {
      free( correct[i/bucketLength] );
      free( false[i/bucketLength] );
    }
  free(correct);
  free(false);

}
 

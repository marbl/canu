
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
static char CM_ID[] = "$Id: GreedyOverlapMainREZ.c,v 1.1.1.1 2004-04-14 13:53:18 catmandew Exp $";

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
  int i,j;
  char *fileName;
  /*Variable defintions for parsing the commandline */
  GenericMesg *pmesg    = NULL; 
  FILE        *fileInput= NULL;
  FILE *repFile[4][7],*normalFile[4][7];

  int intQuality = 0;
  MesgReader   readerFn = NULL;
  ReadStructp  input;
  int noOfFrags = 0;
  int count;
  int length;
  int olapMin,olapMax;
  int qualityThreshold = 0;
  GenericMesg  **fragMesg;

  { /* Parse the argument list using "man 3 getopt". */
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "f:s:e:t:")) != EOF))
      switch(ch) 
	{
	case 'f':
	  fileName = strdup(optarg);
	  break;
	case 's':
	  olapMin = atoi(optarg);
	  break;
	case 'e':
	  olapMax = atoi(optarg);
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
    fprintf (stderr,"USAGE:  greedy-overlap -f .inp file [-t quality threshold][-s min olap] [-e max olap]\n");
    exit(1);
  }
  

  /* Test preconditions */
  assert( fileName != NULL);
  assert( strlen(fileName) != 0 );

  /* set up histogram files */
  for(qualityThreshold=0; qualityThreshold <= 30; qualityThreshold += 5)
    for(olapMin = 0; olapMin <= 300; olapMin += 100)
      {
	char fileName1[100];
	char fileName2[100];
	int qt = qualityThreshold/5;

	if(olapMin == 300)
	  olapMax = 1000;
	else
	  olapMax = olapMin+100;
	
	sprintf(fileName1,"rep.hist.%d.%d.%d",olapMin,olapMax,qualityThreshold);		
	sprintf(fileName2,"normal.hist.%d.%d.%d",olapMin,olapMax,qualityThreshold);
	repFile[olapMin/100][qt]    = file_open(fileName1,"w");  
	normalFile[olapMin/100][qt] = file_open(fileName2,"w");  
	
	fprintf(repFile[olapMin/100][qt],"Quality of repeat overlaps of length %d-%d qt=%d\n",olapMin,olapMax,qualityThreshold);
	
	fprintf(normalFile[olapMin/100][qt],"Quality of normal overlaps of length  %d-%d qt=%d\n",olapMin,olapMax,qualityThreshold);

	fflush(repFile[olapMin/100][qt]);
	fflush(normalFile[olapMin/100][qt]);
      }


  fileInput = file_open(fileName,"r");  
  readerFn  = InputFileType_AS(fileInput);
  
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
  MismatchTableREZ = fill_QV_Mismatch_table(prob_mismatch);
  MatchTableREZ    = fill_QV_Match_table(prob_match);

  
  /* we align all pairs of fragments and compute their quality */
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
#endif
	    Print_Overlap_AS(stdout, IFG1, IFG2, olap_unit);

	    for(qualityThreshold=0; qualityThreshold<=30; qualityThreshold += 5)
	      {
		int qt = qualityThreshold/5;
		int m  = compute_bayesian_quality(IFG1,IFG2,olap_unit,qualityThreshold,&length,NULL);

		if( m >= 0 )
		  {
		    printf("Thresh = %d, Quality = %lf\n",qualityThreshold,olap_unit->quality);

		    if( repetitive_overlap_sim(IFG1,IFG2,olap_unit) != TRUE &&
			olap_unit->quality > 0.7)
		      printf("HREPEAT\n");
		    if( repetitive_overlap_sim(IFG1,IFG2,olap_unit) == TRUE &&
			olap_unit->quality < 0.2)
		      printf("LNORMAL\n");
		    
		    intQuality = (int) 100000*olap_unit->quality;

		    for(olapMin = 0; olapMin <= 300; olapMin += 100)
		      {
			int fi = olapMin/100;
			if( length > olapMin && length <= olapMax)
			  if( repetitive_overlap_sim(IFG1,IFG2,olap_unit) != TRUE )
			    fprintf(repFile[fi][qt],"%d ",intQuality);
			  else
			    fprintf(normalFile[fi][qt],"%d ",intQuality);
			fflush(repFile[fi][qt]);
		    fflush(normalFile[fi][qt]);
		      }
		  }
	      }

	    /* Realign the alignment. 
	       The macros SUB and DEL are defined in SUBDELREZ.h */
	    /*	    
	    olapp = QV_ReAligner_AS(IFG1,IFG2,olap_unit);

	    olap_quality->aifrag       = olapp->aifrag;
	    olap_quality->bifrag       = olapp->bifrag;
	    olap_quality->ahg          = olapp->ahg;
	    olap_quality->bhg          = olapp->bhg;
	    olap_quality->orientation  = olapp->orientation;
	    olap_quality->overlap_type = olapp->overlap_type;
	    olap_quality->quality      = olapp->quality;
	    olap_quality->min_offset   = olapp->min_offset;
	    olap_quality->max_offset   = olapp->max_offset;
	    olap_quality->polymorph_ct = olapp->polymorph_ct;
	    olap_quality->delta = (signed char*) safe_calloc(sizeof(signed char),strlen(olapp->delta)+1);
	    strcpy(olap_quality->delta,olapp->delta);
	    

#ifdef DEBUG
	    fprintf(stdout,"QV orientation = %d\n",olap_quality->orientation);
	    fprintf(stdout,"afrag = %d\n",olap_quality->aifrag);
	    fprintf(stdout,"bfrag = %d\n",olap_quality->bifrag);
	    fprintf(stdout,"ahang = %d\n",olap_quality->ahg);
	    fprintf(stdout,"bhang = %d\n",olap_quality->bhg);
#endif
	    if( repetitive_overlap_sim(IFG1,IFG2,olap_quality) != TRUE )
	      Print_Overlap_AS(stdout, IFG1, IFG2 ,olap_quality);

	    compute_bayesian_quality(IFG1,IFG2,olap_quality,&length,NULL);
	    printf("Quality = %lf\n",olap_quality->quality);

	    free(olap_quality->delta);
	    free(olap_quality);
	    */
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
#endif

	    Print_Overlap_AS(stdout, IFG1, IFG2, olap_unit);

	    for(qualityThreshold=0; qualityThreshold<=30; qualityThreshold += 5)
	      {
		int qt = qualityThreshold/5;
		int m = compute_bayesian_quality(IFG1,IFG2,olap_unit,qualityThreshold,&length,NULL);

		if( m >= 0)
		  {
		    printf("Thresh = %d, Quality = %lf\n",qualityThreshold,olap_unit->quality);

		    if( repetitive_overlap_sim(IFG1,IFG2,olap_unit) != TRUE &&
			olap_unit->quality > 0.7)
		      printf("HREPEAT\n");
		    if( repetitive_overlap_sim(IFG1,IFG2,olap_unit) == TRUE &&
			olap_unit->quality < 0.2)
		      printf("LNORMAL\n");
		    
		    intQuality = (int) 100000*olap_unit->quality;
		    
		    for(olapMin = 0; olapMin <= 300; olapMin += 100)
		      {
			int fi = olapMin/100;
			if( length > olapMin && length < olapMax)
			  if( repetitive_overlap_sim(IFG1,IFG2,olap_unit) != TRUE )
			fprintf(repFile[fi][qt],"%d ",intQuality);
			  else
			    fprintf(normalFile[fi][qt],"%d ",intQuality);
			fflush(repFile[fi][qt]);
			fflush(normalFile[fi][qt]);
		  }
		  }
	      }
	    /*	    show_alignment(IFG1->sequence,IFG2->sequence,
			   IFG1->iaccession,IFG2->iaccession,
			   &olap,40);*/
	    //	    free(olap_unit->delta);
	    
	    /*
	    olapp = QV_ReAligner_AS(IFG1,IFG2,olap_unit);

	    olap_quality->aifrag       = olapp->aifrag;
	    olap_quality->bifrag       = olapp->bifrag;
	    olap_quality->ahg          = olapp->ahg;
	    olap_quality->bhg          = olapp->bhg;
	    olap_quality->orientation  = olapp->orientation;
	    olap_quality->overlap_type = olapp->overlap_type;
	    olap_quality->quality      = olapp->quality;
	    olap_quality->min_offset   = olapp->min_offset;
	    olap_quality->max_offset   = olapp->max_offset;
	    olap_quality->polymorph_ct = olapp->polymorph_ct;
	    olap_quality->delta = (signed char*) safe_calloc(sizeof(signed char),strlen(olapp->delta)+1);
	    strcpy(olap_quality->delta,olapp->delta);

#ifdef DEBUG 
	    fprintf(stdout,"QV orientation = %d\n",olap_quality->orientation);
	    fprintf(stdout,"afrag = %d\n",olap_quality->aifrag);
	    fprintf(stdout,"bfrag = %d\n",olap_quality->bifrag);
	    fprintf(stdout,"ahang = %d\n",olap_quality->ahg);
	    fprintf(stdout,"bhang = %d\n",olap_quality->bhg);
#endif

	    if( repetitive_overlap_sim(IFG1,IFG2,olap_quality) != TRUE )
	      Print_Overlap_AS(stdout, IFG1, IFG2 ,olap_quality);
	    compute_bayesian_quality(IFG1,IFG2,olap_quality,&length,NULL);

	    printf("Quality = %lf\n",olap_quality->quality);

	    free(olap_quality->delta);
	    free(olap_quality);
	   */
	  }
      }


  /* in the end we free all fragments */
  for(i=0; i<noOfFrags; i++)
    FreeProtoMesg_AS(fragMesg[i]);
  for(qualityThreshold=0; qualityThreshold<=30; qualityThreshold += 5)
    for(olapMin = 0; olapMin <= 300; olapMin+=100)
      {
	fclose(repFile[olapMin/100][qualityThreshold/5]);
	fclose(normalFile[olapMin/100][qualityThreshold/5]);
      }
}
 

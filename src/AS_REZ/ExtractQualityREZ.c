
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
 Module: This function computes for a certain range of UOM's
         their quality value and filter certain classes of
	 them according to a threshold. In addition it writes
	 some statistics to files which are then concatenated
	 by the same perl script that concatenates the messages.
	 The resulting stat files are split into different overlap
	 types (O,M,I,Y and R for the rest) calling
	 celagram gives you the statistics.

 Description:
 Assumptions: 
**********************************************************************/

static char CM_ID[] = "$Id: ExtractQualityREZ.c,v 1.2 2004-09-23 20:25:27 mcschatz Exp $";

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
#include "UtilsREZ.h"
#include "GraphCGW_T.h"
#include "ScaffoldGraph_CGW.h"
#include "Input_CGW.h"
#include "AS_UTL_timer.h"
#include "StatisticsREZ.h"

#define DEBUG

void ProcessIUMsimple(IntUnitigMesg *ium_mesg, ScaffoldGraphT* sgraph){
  MultiAlignT *ma = CreateMultiAlignTFromIUM(ium_mesg, GetNumCIFragTs(sgraph->CIFrags), TRUE);
  int32 length = GetMultiAlignUngappedLength(ma);

  InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg->iaccession, TRUE, ma, TRUE);
  DuplicateEntryInSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg->iaccession, TRUE, ium_mesg->iaccession, FALSE, TRUE);
  //  SetMultiAlignInStore(sgraph->CIGraph->maStore, ium_mesg->iaccession, ma);
  //  SetMultiAlignInStore(sgraph->ContigGraph->maStore, ium_mesg->iaccession, ma);
  ProcessIUM_ScaffoldGraph(ium_mesg, length, TRUE);
}



int main(int argc, char *argv[]) 
{
  int i,j;
  char *fileName;
  char ofileName[200];
  char statFileName[200];
  char fileNameExt[200];
  ScaffoldGraphT *sgraph;
  TimerT timer;
  
  GenericMesg *pmesg      = NULL; 
  FILE        *fileInput  = NULL;
  FILE        *fileOutput = NULL;
  FILE        *statOutputO = NULL;
  FILE        *statOutputM = NULL;
  FILE        *statOutputI = NULL;
  FILE        *statOutputY = NULL;
  FILE        *statOutputX = NULL;
  FILE        *statOutputR = NULL;

  int intQuality = 0;
  MesgReader   readerFn = NULL;
  MesgWriter   writerFn = NULL;
  OutputType   outputType = AS_PROTO_OUTPUT; 
  // we have to use PROTO I/O. Concatenating binary files doe not work
                                               
  ReadStructp  input;
  float qualityThreshold = 0;
  int blockSize = 0;
  int blockNumber = 0;

  int notConfirmed = 0; // number of UOMs that did not have an overlap
  int tooBad       = 0; // number of UOMs whose overlap was too bad
  int goodEnough   = 0; // number of UOMs whose overlap was good enough
  int UOMcount     = 0; // number of UOMs
  int IUMcount     = 0; // number of IUMs
  int badKept      = 0; // number of UOMs of type O,M or I or Y 
                        // that have not sufficient quality but are
                        // kept nevertheless
  int processedUOMs = 0;
  int noOfNodes     = 2048;

  int processRest = FALSE;
  int timerStarted = FALSE;

  { /* Parse the argument list using "man 3 getopt". */
    int ch,errflg=0 ;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "rf:t:b:n:a:")) != EOF))
      switch(ch) 
	{
	case 'f':
	  fileName = strdup(optarg);
	  break;
	case 'r':
	  processRest = TRUE;
	  break;
	case 't':
	  qualityThreshold = atof(optarg);
	  break;
	case 'b':
	  blockSize = atoi(optarg);
	  break;
	case 'a':
	  noOfNodes = atoi(optarg);
	  break;
	case 'n':
	  blockNumber = atoi(optarg);
	  break;
	case '?':
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	  break;
	}
  }
   
  // make sure there is a uom directory
  init_uoms_dir();
  init_stat_dir();

#ifdef DEBUG
  fprintf(stderr,"\nProcessing UOMs from %d to %d\n",blockNumber*blockSize,(blockNumber+1)*blockSize);  
#endif


  if( (argc - optind != 0 )){
    fprintf (stderr,"USAGE:  extract-quality-uoms [-r (process rest of messages)]\n-f file (prefix of cgb file)\n-t quality_threshold (set threshold)\n-b blocksize (set blocksize for UOMs)\n-n block_num (set block number)\n[-a alloc (preallocate space for the graph)\n");
    exit(1);
  }
  
  sgraph = CreateScaffoldGraph(TRUE,fileName,noOfNodes,0);
  sgraph->RezGraph = sgraph->CIGraph;
  ScaffoldGraph = sgraph;

  /* Test preconditions */
  assert( fileName != NULL);
  assert( strlen(fileName) != 0 );

  // the input file has the extension cgb
  strcpy(fileNameExt,fileName);
  strcat(fileNameExt,".cgb");
  
  // the output file is in the uom directory and has the 
  // extension of the block processed
  sprintf(ofileName,"./uoms/%s.%d",fileName,blockNumber);


  fprintf(stderr,ofileName);


  fileInput  = file_open(fileNameExt,"r");  
  fileOutput = file_open(ofileName,"w");  


  /* we open the celagram files
     if the -r flag is set (== first batch) we write a header line 
     in the files  */  
  /* --- O --- */
  sprintf(statFileName,"./uoms/%s.O.stat.%d",fileName,blockNumber);
  statOutputO = file_open(statFileName,"w");  
  if( processRest == TRUE)
    fprintf(statOutputO,"Qualities of overlap type 'O' (AS_OVERLAP)\n");
  /* --- M --- */
  sprintf(statFileName,"./uoms/%s.M.stat.%d",fileName,blockNumber);
  statOutputM = file_open(statFileName,"w");  
  if( processRest == TRUE)
    fprintf(statOutputM,"Qualities of overlap type 'M' (AS_TOUCHES_CONTAINED_OVERLAP)\n");
  /* --- I --- */
  sprintf(statFileName,"./uoms/%s.I.stat.%d",fileName,blockNumber);
  statOutputI = file_open(statFileName,"w");  
  if( processRest == TRUE)
    fprintf(statOutputI,"Qualities of overlap type 'I' (AS_2_CONTAINS_1_OVERLAP)\n");
  /* --- Y --- */
  sprintf(statFileName,"./uoms/%s.Y.stat.%d",fileName,blockNumber);
  statOutputY = file_open(statFileName,"w");  
  if( processRest == TRUE)
    fprintf(statOutputY,"Qualities of overlap type 'Y' (AS_BETWEEN_CONTAINED_OVERLAP)\n");
  /* --- X --- */
  sprintf(statFileName,"./uoms/%s.X.stat.%d",fileName,blockNumber);
  statOutputX = file_open(statFileName,"w");  
  if( processRest == TRUE)
    fprintf(statOutputX,"Qualities of overlap type 'X' (AS_TRANSCHUNK_OVERLAP)\n");
  /* --- R (rest) --- */
  sprintf(statFileName,"./uoms/%s.R.stat.%d",fileName,blockNumber);
  statOutputR = file_open(statFileName,"w");  
  if( processRest == TRUE)
    fprintf(statOutputR,"Qualities of remaining overlap types\n");

  readerFn  = InputFileType_AS(fileInput);
  writerFn  = OutputFileType_AS(outputType);

  InitTimerT(&timer);

  /* read the fragment array */
  while( readerFn(fileInput,&pmesg) != EOF ) 
    {
      MessageType mesgtype = pmesg->t;

      if(mesgtype  == MESG_IUM)
	{
	  IntUnitigMesg  *iumMesg = (IntUnitigMesg*) pmesg->m;
	  ProcessIUMsimple(iumMesg,sgraph);
	  IUMcount++;
	  if( processRest == TRUE)
	    {
	      writerFn(fileOutput,pmesg);
	      fflush(fileOutput);
	    }
	}
      else
	if(mesgtype == MESG_UOM )
	  {
	    float quality;
	    ChunkOverlapCheckT olap;
	    int found;
	    UnitigOverlapMesg *uom_mesg = (UnitigOverlapMesg *) pmesg->m;

	    UOMcount++;
	    if( UOMcount < blockNumber*blockSize)
	      continue;
	    if( UOMcount >= blockNumber*blockSize+blockSize)
	      break;
	    {

	      if( ! timerStarted )
		{
		  StartTimerT(&timer);
		  timerStarted = TRUE;
		}

	      processedUOMs++;
	      found = ComputeUOMQualityOverlap(sgraph->RezGraph,
					       uom_mesg,
					       &olap,&quality);
	     

#ifdef DEBUG
	      if(UOMcount % 100000 == 0)
		fprintf(stderr,"*");
	      else
		if(UOMcount % 10000 == 0)
		  fprintf(stderr,".");
#endif
	      if( ! found )
		{
		  notConfirmed++;
		}
	      else
		{
		  UnitigOverlapType overlap_type = uom_mesg->overlap_type;
		  int intQual = (int) 10000*quality;
		  // adapt the best overlap field and the quality
		  uom_mesg->quality = quality;
		  uom_mesg->best_overlap_length = olap.overlap;

		  /* write the scores for each edge class into a celagram file */
		  switch( overlap_type ){
		  case AS_TOUCHES_CONTAINED_OVERLAP :
		    /* M */
		    fprintf(statOutputM,"%d ",intQual); 
		    break;
		  case AS_OVERLAP :
		    /* O */
		    fprintf(statOutputO,"%d ",intQual); 
		    break;
		  case AS_2_CONTAINS_1_OVERLAP :
		    /* I */
		    fprintf(statOutputI,"%d ",intQual); 
		    break;
		  case AS_BETWEEN_CONTAINED_OVERLAP :
		    /* Y */
		    fprintf(statOutputY,"%d ",intQual); 
		    break;
		  case AS_TRANSCHUNK_OVERLAP :
		    /* X */
		    fprintf(statOutputX,"%d ",intQual); 
		    break;
		  default:
		    /* R */
		    fprintf(statOutputR,"%d ",intQual); 
		    break;
		  }
		  if( quality > qualityThreshold)
		    {
		      // count the bad ones
		      tooBad++;
		      if( overlap_type == AS_TOUCHES_CONTAINED_OVERLAP || /* M */
			  overlap_type == AS_OVERLAP || /* O */
			  overlap_type == AS_2_CONTAINS_1_OVERLAP || /* I */
			  overlap_type == AS_BETWEEN_CONTAINED_OVERLAP) /* Y */
			{
			  /*if the type is one of the above we keep it anyway */ 
			  writerFn(fileOutput,pmesg);
			  badKept++;
			}
		    }
		  else
		    {
		      goodEnough++;
		      writerFn(fileOutput,pmesg);
		    }
		}
	    }
	  }
	else
	  {
	    if( processRest == TRUE)
	      {
		writerFn(fileOutput,pmesg);
		fflush(fileOutput);
	      }
	  }
    }
  
  {
    double time = StopTimerT(&timer);
    fprintf(stderr,"\nStatistics **************\n");
    fprintf(stderr,"*     UOMS %d\n",UOMcount);
    fprintf(stderr,"*    PUOMS %d\n",processedUOMs);
    fprintf(stderr,"*     IUMS %d\n",IUMcount);
    fprintf(stderr,"* NC  UOMS %d\n",notConfirmed);
    fprintf(stderr,"* ++  UOMS %d\n",goodEnough);
    fprintf(stderr,"* --  UOMS %d\n",tooBad);
    fprintf(stderr,"* Bad Kept %d\n",badKept);
    fprintf(stderr,"* =============================\n");
#ifdef DEBUG
    fprintf(stderr,"* Time elapsed %f\n",time);
    if( processedUOMs > 0 )
      fprintf(stderr,"* Processed %f UOMs/sec\n",processedUOMs/time);
#endif

  }

  /* close all files */
  fclose(fileInput);
  fclose(fileOutput);
  fclose(statOutputO);
  fclose(statOutputM);
  fclose(statOutputI);
  fclose(statOutputY);
  fclose(statOutputX);
  fclose(statOutputR);

  exit(0);
}
 

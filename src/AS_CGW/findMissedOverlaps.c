
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
static char CM_ID[] = "$Id: findMissedOverlaps.c,v 1.12 2007-08-03 20:45:03 brianwalenz Exp $";


/*********************************************************************/

//#define DEBUG 1
//#define DEBUG_BUCIS 1
//#define DEBUG_MERGE_SCAF 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "FbacREZ.h"
#include "PublicAPI_CNS.h"
#include "AS_ALN_aligners.h"
#include "AS_ALN_forcns.h"

#define NUM_STDDEV_CUTOFF 3.0

// externable variables for controlling use of Local_Overlap_AS_forCNS
extern int MaxGaps;		    // [initialized value is 12 -- no more than this many segments in the chain]
extern CDS_COORD_t MaxBegGap;       // [ init value is 200; this could be set to the amount you extend the clear 
// range of seq b, plus 10 for good measure]
extern CDS_COORD_t MaxEndGap;	    // [ init value is 200; this could be set to the amount you extend the 
// clear range of seq a, plus 10 for good measure]
extern CDS_COORD_t MaxInteriorGap;	// [ initial value is 1000 (should have almost no effect)
// and defines the largest gap between segments in the chain]
// Also: allowed size of gap within the alignment
// -- forcing relatively good alignments, compared to those
// allowed in bubble-smoothing where indel polymorphisms are expected
extern int asymmetricEnds;  // boolean to cause the size of an "end gap" to
                            // be evaluated with regard to the clear range extension
extern CDS_COORD_t max_indel_AS_ALN_LOCOLAP_GLOBAL;

static VA_TYPE(char) *lContigConsensus = NULL;
static VA_TYPE(char) *rContigConsensus = NULL;
static VA_TYPE(char) *lContigQuality = NULL;
static VA_TYPE(char) *rContigQuality = NULL;

LengthT FindGapLength( ChunkInstanceT * lchunk, ChunkInstanceT * rchunk, int verbose);

void SequenceComplement(char *sequence, char *quality);


void initVarArrays(void)
{
  if ( lContigConsensus == NULL )
    {
      lContigConsensus = CreateVA_char(1024);
      rContigConsensus = CreateVA_char(1024);
      lContigQuality = CreateVA_char(1024);
      rContigQuality = CreateVA_char(1024);
    }
}

int evalGap( CDS_CID_t sid,
	     ContigT *lcontig, 
	     ContigT *rcontig,
	     int gapNumber,int dumpSeqs)
{
  char *lSequence, *rSequence;
  NodeOrient lContigOrientation, rContigOrientation;
  int mean;
  int stddev;
  CDS_COORD_t beg, end, widebeg, wideend;
  int opposite = FALSE;
  double erate, thresh, minlen;
  CompareOptions what;
  Overlap *overlap;
  LengthT gapSize=FindGapLength( lcontig, rcontig, FALSE);
  CDS_COORD_t lLen, rLen;

  erate  = 0.06;
  thresh = 1e-6;
  minlen = 30;

  if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean){
    lContigOrientation = A_B;
  }else{
    lContigOrientation = B_A;
  }
  
  if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean){
    rContigOrientation = A_B;
  }else{
    rContigOrientation = B_A;
  }

  mean=gapSize.mean;
  stddev=sqrt(gapSize.variance);

  if(mean-NUM_STDDEV_CUTOFF*stddev>CNS_DP_MINLEN)return 0;

  // Get the consensus sequences for both chunks from the Store
  GetConsensus( ScaffoldGraph->ContigGraph, lcontig->id, lContigConsensus, lContigQuality);
  GetConsensus( ScaffoldGraph->ContigGraph, rcontig->id, rContigConsensus, rContigQuality);
  
  lSequence = Getchar( lContigConsensus, 0);
  lLen = (CDS_COORD_t) strlen(lSequence);

  rSequence = Getchar( rContigConsensus, 0);
  rLen = (CDS_COORD_t) strlen(rSequence);
  
  if (lContigOrientation == B_A)  //
    {
      SequenceComplement( lSequence, NULL);
    }
  
  if (rContigOrientation == B_A)  // the frag is oriented opposite to the contig in this case
    {
      SequenceComplement( rSequence, NULL);  // flip contig sequence to its orientation in scaffold
    }
  
  if(dumpSeqs){
    fprintf(stdout,"\n>ctg %d%s\n%s\n>ctg %d%s\n%s\n",
	    lcontig->id,(lContigOrientation == B_A ? " rc" : ""),lSequence,
	    rcontig->id,(rContigOrientation == B_A ? " rc" : ""),rSequence);
  }

  
  what       = AS_FIND_ALIGN;
  opposite   = 0;

  // mean is mean gap length: negative implies overlap
  beg        = MIN(lLen,lLen+mean-NUM_STDDEV_CUTOFF*stddev);
  end        = MIN(lLen,lLen+mean+NUM_STDDEV_CUTOFF*stddev);

  widebeg = MAX(-rLen,beg-10000);
  wideend = MIN(lLen,end+10000);


  fprintf(stdout,"testing .... Scaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d]  band " F_COORD " " F_COORD " wideband " F_COORD " " F_COORD,
	  sid,
	  lcontig->id,
	  rcontig->id,
	  gapNumber,
	  mean,stddev,
	  beg,
	  end,
	  widebeg,
	  wideend);


  overlap = DP_Compare( lSequence, rSequence,
			beg, end, opposite,
			erate, thresh, minlen,
			what);
	  

  fprintf(stdout,".");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (simple) , ahang " F_COORD " overlap " F_COORD " error %f\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    overlap->begpos, 
	    lLen-overlap->begpos,
	    overlap->diffs/(double)overlap->length);
    return 1;
  }


  overlap = DP_Compare( lSequence, rSequence,
			widebeg, wideend, opposite,
			erate, thresh, minlen,
			what);
	  
  fprintf(stdout,".");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (wide) , ahang " F_COORD " overlap " F_COORD " error %f\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    overlap->begpos, 
	    lLen-overlap->begpos,
	    overlap->diffs/(double)overlap->length);
    return 1;
  }


  overlap = DP_Compare( rSequence, lSequence,
			-wideend, -widebeg, opposite,
			erate, thresh, minlen,
			what);
	  
  fprintf(stdout,".");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (BAswap) , ahang " F_COORD " overlap " F_COORD " error %f\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    -overlap->begpos, 
	    lLen+overlap->begpos,
	    overlap->diffs/(double)overlap->length);
    return 1;
  }


  SequenceComplement(lSequence,0);
  SequenceComplement(rSequence,0);
  overlap = DP_Compare( lSequence, rSequence,
			-(rLen+mean+NUM_STDDEV_CUTOFF*stddev)-10000,
			-(rLen+mean-NUM_STDDEV_CUTOFF*stddev)+10000,
			opposite,
			erate, thresh, minlen,
			what);
	  
  fprintf(stdout,".");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (revOri) , ahang " F_COORD " overlap " F_COORD " error %f\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    -overlap->endpos, 
	    lLen-overlap->endpos,
	    overlap->diffs/(double)overlap->length);
    return 1;
  }

  SequenceComplement(lSequence,0);
  SequenceComplement(rSequence,0);


  opposite=1;

  overlap = DP_Compare( lSequence, rSequence,
			widebeg, wideend, opposite,
			erate, thresh, minlen,
			what);
	  
  fprintf(stdout,".");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (misorient1) , ahang " F_COORD " overlap " F_COORD " error %f\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    overlap->begpos, 
	    lLen-overlap->begpos,
	    overlap->diffs/(double)overlap->length);
    return 1;
  }

  overlap = DP_Compare( lSequence, rSequence,
			-rLen+lLen-wideend,-rLen+lLen-widebeg,opposite,
			erate, thresh, minlen,
			what);
	  
  fprintf(stdout,".");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (misorient2) , ahang " F_COORD " overlap " F_COORD " error %f\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    overlap->begpos, 
	    lLen-overlap->begpos,
	    overlap->diffs/(double)overlap->length);
    return 1;
  }


	
  //  overlapper = Local_Overlap_AS_forCNS;

  what       = AS_FIND_LOCAL_ALIGN;
  MaxGaps = 5;  // we don't want a lot of segments
  MaxEndGap = 200;
  MaxBegGap = 200;
  MaxInteriorGap = 400;
  asymmetricEnds = 0;

  opposite   = 0;


  overlap = Local_Overlap_AS_forCNS( lSequence, rSequence,
                                     beg, end, opposite,
                                     erate, thresh, minlen,
                                     what);
	  
  fprintf(stdout,",");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (simpleLOCO) , ahang " F_COORD " overlap " F_COORD " error %f maxGap " F_COORD "\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    overlap->begpos, 
	    lLen-overlap->begpos, 
	    overlap->diffs/(double)overlap->length,
	    max_indel_AS_ALN_LOCOLAP_GLOBAL);
    return 1;
  }


  overlap = Local_Overlap_AS_forCNS( lSequence, rSequence,
                                     widebeg, wideend, opposite,
                                     erate, thresh, minlen,
                                     what);
	  
  fprintf(stdout,",");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (wideLOCO) , ahang " F_COORD " overlap " F_COORD " error %f maxGap " F_COORD "\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    overlap->begpos, 
	    lLen-overlap->begpos, 
	    overlap->diffs/(double)overlap->length,
	    max_indel_AS_ALN_LOCOLAP_GLOBAL);
    return 1;
  }


  overlap = Local_Overlap_AS_forCNS( rSequence, lSequence,
                                     -wideend, -widebeg, opposite,
                                     erate, thresh, minlen,
                                     what);
	  
  fprintf(stdout,",");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (BAswapLOCO) , ahang " F_COORD " overlap " F_COORD " error %f maxGap " F_COORD "\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    -overlap->begpos, 
	    lLen+overlap->begpos, 
	    overlap->diffs/(double)overlap->length,
	    max_indel_AS_ALN_LOCOLAP_GLOBAL);
    return 1;
  }


  SequenceComplement(lSequence,0);
  SequenceComplement(rSequence,0);
  overlap = Local_Overlap_AS_forCNS( lSequence, rSequence,
                                     -(rLen+mean+NUM_STDDEV_CUTOFF*stddev)-10000,
                                     -(rLen+mean-NUM_STDDEV_CUTOFF*stddev)+10000,
                                     opposite,
                                     erate, thresh, minlen,
                                     what);
	  
  fprintf(stdout,",");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (revOriLOCO) , ahang " F_COORD " overlap " F_COORD " error %f maxGap " F_COORD "\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    -overlap->endpos, 
	    lLen+overlap->endpos, 
	    overlap->diffs/(double)overlap->length,
	    max_indel_AS_ALN_LOCOLAP_GLOBAL);
    return 1;
  }

  SequenceComplement(lSequence,0);
  SequenceComplement(rSequence,0);


  opposite=1;

  overlap = Local_Overlap_AS_forCNS( lSequence, rSequence,
                                     widebeg, wideend, opposite,
                                     erate, thresh, minlen,
                                     what);
	  
  fprintf(stdout,",");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (misorient1LOCO) , ahang " F_COORD " overlap " F_COORD " error %f maxGap " F_COORD "\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    overlap->begpos, 
	    lLen-overlap->begpos,
	    overlap->diffs/(double)overlap->length,
	    max_indel_AS_ALN_LOCOLAP_GLOBAL);
    return 1;
  }

  overlap = Local_Overlap_AS_forCNS( lSequence, rSequence,
				     -rLen+lLen-wideend,-rLen+lLen-widebeg,opposite,
				     erate, thresh, minlen,
				     what);
	  
  fprintf(stdout,",");
  if(overlap!=NULL){
    fprintf(stdout,"\nScaf " F_CID " Gap " F_CID " ^ " F_CID " %d [%d,%d] closed (misorient2LOCO) , ahang " F_COORD " overlap " F_COORD " error %f maxGap " F_COORD "\n",
            sid,lcontig->id,rcontig->id,
	    gapNumber,
	    mean,stddev,
	    overlap->begpos, 
	    overlap->length,
	    overlap->diffs/(double)overlap->length,
	    max_indel_AS_ALN_LOCOLAP_GLOBAL);
    return 1;
  }

  fprintf(stdout,"\n");
  return 0;
}





int main( int argc, char *argv[])
{
  Global_CGW *data;
  char *outputPath = NULL;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int setSingleSid = FALSE;
  CDS_CID_t singleSid = NULLINDEX;
  int ckptNum = NULLINDEX;
  CDS_CID_t sid;
  int startingGap = 0, setStartingGap = FALSE;
  int numGaps = 0, gapNumber = 0;
  int numGapsVarTooSmall = 0;
  int overlapFound = 0, noOverlapFound = 0;
  int dumpSeqs=0;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "c:C:df:g:m:n:s:")) != EOF)){
      switch(ch) {
        case 'c':
          {
            strcpy( data->File_Name_Prefix, argv[optind - 1]);
            setPrefixName = TRUE;		  
          }
          break;
        case 'C':
          startingGap = atoi(argv[optind - 1]);
          setStartingGap = TRUE;
          break;
        case 'd':
          {
            dumpSeqs=1;
          }
          break;
        case 'g':
          {
            strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
            setGatekeeperStore = TRUE;
          }
          break;	  
        case 'm':
          MaxInteriorGap = atoi(argv[optind - 1]);
          fprintf( stderr, "setting MaxInteriorGap to " F_COORD "\n",
                   MaxInteriorGap);
          break;
        case 'n':
          ckptNum = atoi(argv[optind - 1]);
          break;
        case 's':
          singleSid = atoi(argv[optind - 1]);
          setSingleSid = TRUE;
          fprintf( stderr, "setting singleSid to " F_CID "\n",
                   singleSid);
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if((setPrefixName == FALSE) || (setGatekeeperStore == 0))
      {
	fprintf(stderr,"* argc = %d optind = %d setGatekeeperStore = %d outputPath = %s\n",
		argc, optind,setGatekeeperStore, outputPath);
	fprintf (stderr, "USAGE:  %s [-d] -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum>\n\t-d option causes contig seqeunces to be dumped\n",argv[0]);
	exit (EXIT_FAILURE);
      }
  }

  if (setStartingGap == TRUE)
    fprintf( stderr, "set starting gap to %d\n", startingGap);
  
  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);
  GlobalData->aligner=Local_Overlap_AS_forCNS;
  initVarArrays();

  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++){

    CIScaffoldTIterator CIsTemp;
    CIScaffoldT * scaff;
    int icnt;
    CDS_CID_t rcontigID;
    ContigT *lcontig, *rcontig;
	
    scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
    // make sure the scaffold is there
    assert(scaff != NULL);
    
    // not interested in dead scaffold, not real scaffolds, or singleton scaffolds
    
    if ((isDeadCIScaffoldT(scaff)) ||
	(scaff->type != REAL_SCAFFOLD) ||
	(scaff->info.Scaffold.numElements < 2))
      {
        continue;
      }

    /*
      fprintf(stdout,"\n=====================================================================\n");
      fprintf(stdout,"=== examing scaffold " F_CID ", size %f\n",
      sid, scaff->bpLength.mean);
      fprintf(stderr,"Scaffold " F_CID "\n", sid);
    */
    
    // speed hack
    if ( setSingleSid == TRUE && sid != singleSid) continue;

	

    // make sure the scaffold is there
    assert(scaff != NULL);

    InitCIScaffoldTIterator( ScaffoldGraph, scaff, TRUE, FALSE, &CIsTemp);
    while ( NextCIScaffoldTIterator( &CIsTemp ) && 
	    CIsTemp.next != GetNumGraphNodes( ScaffoldGraph->ContigGraph ) - 1)
      {
        /*
          fprintf( stderr, "prev, curr, next: " F_CID ", " F_CID ", " F_CID "\n",
          CIsTemp.prev, CIsTemp.curr, CIsTemp.next);
        */
      }

  
    icnt = 0;

    lcontig = GetGraphNode( ScaffoldGraph->ContigGraph,
                            scaff->info.Scaffold.AEndCI);
    rcontigID = lcontig->BEndNext;
    while ( rcontigID != -1 ){
      NodeOrient lcontigOrientation, rcontigOrientation;
      LengthT gapSize;

      rcontig = GetGraphNode( ScaffoldGraph->ContigGraph, rcontigID);
	  
      assert(lcontig != NULL);
      assert(rcontig != NULL);		

      gapSize = FindGapLength( lcontig, rcontig, FALSE);
      numGaps++;

      if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
	lcontigOrientation = A_B;
      else
	lcontigOrientation = B_A;

      if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
	rcontigOrientation = A_B;
      else
	rcontigOrientation = B_A;

      /*
        fprintf( stdout, "\n\n\n---------------------------------------------------------------\n");
        fprintf( stdout, "examining gap %d from " F_CID " (orient: %c, pos: %f, %f) to " F_CID " (orient: %c, pos: %f, %f), size: %f \n", 
        gapNumber, 
        lcontig->id, lcontigOrientation,
        lcontig->offsetAEnd.mean, lcontig->offsetBEnd.mean,
        rcontig->id, rcontigOrientation,
        rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean,
        gapSize.mean);
      */
	  
      if ( setStartingGap == TRUE && gapNumber < startingGap ) continue;

      if ( gapSize.mean > NUM_STDDEV_CUTOFF * sqrt(gapSize.variance) && gapSize.mean > 100.0){
	numGapsVarTooSmall++;
      }else{ 
	if(evalGap( sid, lcontig, rcontig, gapNumber,dumpSeqs)) {
	  overlapFound++;
	}else{
	  noOverlapFound++;
	}

      }
      gapNumber++;
      lcontig = rcontig;
      rcontigID = lcontig->BEndNext;
    }
    fflush(stdout);
  }

  fprintf(stdout,"END OF findMissedOverlaps on %s ckp %d : %d gaps with %d unclosable %d unclosed and %d potentially closed\n",
	  data->File_Name_Prefix,
	  ckptNum,
	  gapNumber,
	  numGapsVarTooSmall,
	  noOverlapFound,
	  overlapFound);

  exit(0);
}

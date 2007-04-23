
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
#define FILTER_EDGES
static char CM_ID[] = "$Id: Input_CGW.c,v 1.33 2007-04-23 15:24:34 brianwalenz Exp $";

/*   THIS FILE CONTAINS ALL PROTO/IO INPUT ROUTINES */


//#define DEBUG 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_version.h"
#include "AS_CGW_dataTypes.h"
#include "AS_PER_gkpStore.h"
#include "ScaffoldGraph_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "Input_CGW.h"


static int32 DiscriminatorUniques = 0;
static int32 ShortDiscriminatorUniques = 0;

// Statistics Collection on Guides
static	int32 totalGuideFrags = 0;
static	int32 inUniqueGuideFrags = 0;
static	int32 inRepeatGuideFrags = 0;
static	int32 inTeenyUnitigGuideFrags = 0;
static	int32 inSingletonUnitigGuideFrags = 0;
static int32 onEndGuideFrags = 0;

// Statistics Collection on Reads
static	int32 totalReadFrags = 0;
static	int32 inUniqueReadFrags = 0;
static	int32 inRepeatReadFrags = 0;
static	int32 inTeenyUnitigReadFrags = 0;
static	int32 inSingletonUnitigReadFrags = 0;
static int32 onEndReadFrags = 0;

static int32 TouchesContained = 0;
static int32 TransChunk = 0;
static int32 Containment = 0;
static int32 DoveTail = 0;
static int32 Tandem = 0;
static int32 TotalUOMs = 0;
static int32 TotalUOMHandled = 0;
static int32 BetweenContained = 0;
static int32 ContainStack = 0;
static int32 BadQuality = 0;


static void ProcessUOM(UnitigOverlapMesg *uom_mesg, float transQualityCutoff);
static void ProcessIUM(IntUnitigMesg *ium_mesg);
static void ProcessFrags(void);
static void ProcessEdges(int32 maxDegree, int32 maxDegreeUnique);


static int32 msgCnt = 0;
static int32 numIDT = 0, numIUM = 0, numUOM = 0, numJNC = 0;

/****************************************************************************/
int ProcessInput(Global_CGW *data, int optind, int argc, char *argv[]){
  GenericMesg   *pmesg;
  FILE *infp;
  int i,j = 0;

  StartTimerT(&GlobalData->InputTimer);

  for(i = optind; i < argc; i++){
    fprintf(stderr,"* Opening file %d %s\n", j ++, argv[i]);
    infp = fopen(argv[i],"r");

    while(  (EOF != ReadProtoMesg_AS(infp, &pmesg))){

      switch(pmesg->t){

        case MESG_IDT:
          //  We load these from the gatekeeper store
          numIDT++;
          break;
        case MESG_UOM:  // UnitigOverlapMesg
          if (ScaffoldGraph->ignoreUOMs == 0)
            {
              if(numUOM == 0){
                fprintf(stderr,"* After reading all IUMs \n");
                ReportMemorySize(ScaffoldGraph, stderr);
              }
              numUOM++;
              ProcessUOM(pmesg->m, data->transQualityCutoff);
            }
          break;
        case MESG_IUM:  // IntUnitigMesg
          numIUM++;
          if(numIUM % 10000 == 0)
            fprintf(stderr,"Read %d unitigs\n",numIUM);
          ProcessIUM(pmesg->m);
          break;
        case MESG_ADT:
          // We already processes this
          break;
        case MESG_IBA:
        case MESG_FOM:
          break;

        default:
          fprintf(stderr,"* Oops: Read Message with type = %d\n", pmesg->t);
          if (data->cgwfp)
            WriteProtoMesg_AS(data->cgwfp,pmesg);
          break;
      }
    
    }
    fclose(infp);
  }
  
  ScaffoldGraph->numLiveCIs = ScaffoldGraph->numOriginalCIs = GetNumGraphNodes(ScaffoldGraph->CIGraph);
  fprintf(stderr,"* cgw read the following messages:\n");
  fprintf(stderr,"\tIUM:%d  (max IUM acc = %d) with %d fragments\n",
	  numIUM, (int) GetNumGraphNodes(ScaffoldGraph->CIGraph),
	  (int) GetNumCIFragTs(ScaffoldGraph->CIFrags));
  fprintf(stderr,"\tTotal UOM at Input:%d  BadQuality:%d Processed:%d\n",
	  numUOM,	BadQuality, TotalUOMHandled);
  fprintf(stderr,"\tUOMs by Type\n");
  fprintf(stderr,"\t\t containment: %d  overlap: %d tandem: %d\n",
	  Containment, DoveTail, Tandem);
  fprintf(stderr,"\t\t touchesContained: %d  betweenContained: %d ContainStack: %d\n",
	  TouchesContained, BetweenContained, ContainStack);
  fprintf(stderr,"\t\t transchunk: %d\n", TransChunk);
  if(numUOM != (Containment + DoveTail + Tandem + TouchesContained + BetweenContained + ContainStack + TransChunk)){
    fprintf(stderr,"Total edges read:%d != sum of types %d\n",
	    numUOM ,(Containment + DoveTail + Tandem + TouchesContained + BetweenContained + ContainStack + TransChunk));
  }
  if(numIDT)
    fprintf(stderr,"\tIDT:%d  (ignored...read from the gatekeeper store)\n",  numIDT);
  if(numJNC)
    fprintf(stderr,"\tJNC:%d  \n",  numJNC);
  
  /***vvvvv********  Process UOMs **************************vvvvvvvvvvv******/
  ProcessEdges(data->maxDegree, data->maxDegreeUnique);

  /***vvvvv********  Process Frags *************************vvvvvvvvvvv******/
  ProcessFrags();
  /***^^^^^********  End Process Frags ***************************^^^^^******/

  fprintf(stderr,"* Total Long Discriminator Uniques : %d   Short Uniques: %d\n",
	  DiscriminatorUniques, ShortDiscriminatorUniques);
  fprintf(stderr,"* Total Guides:%d in discriminator unique:%d in other:%d ; in teeny: %d in singles:%d on ends:%d\n",
	  totalGuideFrags, inUniqueGuideFrags, inRepeatGuideFrags, inTeenyUnitigGuideFrags, inSingletonUnitigGuideFrags,
	  onEndGuideFrags);
  fprintf(stderr,"* Total Reads:%d in discriminator unique:%d in other:%d ; in teeny: %d in singles:%d on ends:%d\n",
	  totalReadFrags, inUniqueReadFrags, inRepeatReadFrags, inTeenyUnitigReadFrags, inSingletonUnitigReadFrags,
	  onEndReadFrags);

  StopTimerT(&GlobalData->InputTimer);

  fprintf(stderr,"* CGW Input took %g seconds\n",
	  LapTimerT(&GlobalData->InputTimer));

  return(0);
}


/* ProcessInputADT:
   This is called when loading a checkpoint.  We still read the
   ADT stuff from the cgb file.
*/
int ProcessInputADT(Global_CGW *data, FILE *infp, int argc, char **argv){
  GenericMesg   *pmesg;
  int finished = FALSE;
  char buffer[16 * 1024];
  AuditMesg *adt_mesg = NULL;
  AuditMesg dummy_mesg;
  GenericMesg dummypmesg;

  int i;

  strcpy(buffer,"Invoked: ");
  for(i = 0; i < argc; i++){
    strcat(buffer,argv[i]);
    strcat(buffer," ");
    assert(strlen(buffer) < (16 * 1024));
  }
  ReadProtoMesg_AS(infp, &pmesg);  // Read 1 message

  switch(pmesg->t){

    case MESG_ADT:
      {
	adt_mesg = pmesg->m;
	pmesg->t = MESG_ADT;
      }
      finished = TRUE;
      break;


    default:
      fprintf(stderr,"* Oops: Read Message with type = %d ... expecting ADT message\n", pmesg->t);
      break;
  }
    
  if(!finished){
    // cook up an ADT message
    adt_mesg = &dummy_mesg;
    dummypmesg.m = adt_mesg;
    dummypmesg.t = MESG_ADT;
    pmesg = &dummypmesg;
    adt_mesg->list = NULL;
  }

  VersionStampADT(adt_mesg, argc, argv);

  if (data->cgwfp)
    WriteProtoMesg_AS(data->cgwfp,pmesg);

  return(!finished);
}



void ProcessUOM(UnitigOverlapMesg *uom_mesg, float transQualityCutoff){
  CDS_COORD_t range  = (uom_mesg->max_overlap_length - uom_mesg->min_overlap_length);
  int   hasTandemOverlap = (range > CGB_TANDEM_REPEAT_THRESHOLD);
  // For historical reasons, we prefer to insert edges in canonical form
  IntChunk_ID cidA = uom_mesg->chunk1;
  IntChunk_ID cidB = uom_mesg->chunk2;
  float32  quality = uom_mesg->quality;
  ChunkOrientationType orient = uom_mesg->orient;
  UnitigOverlapType overlap_type = uom_mesg->overlap_type;


  switch(overlap_type){
    case AS_TOUCHES_CONTAINED_OVERLAP: // M
      TouchesContained++;
      break;
    case AS_TANDEM_OVERLAP:   // T
      Tandem++;
      break;
    case AS_BETWEEN_CONTAINED_OVERLAP: // Y
      BetweenContained++;
      if(ScaffoldGraph->ignoreUOMBetweenContained)
        return;
      break;

    case AS_OVERLAP:          // O
      DoveTail++;
      break;
    case AS_TRANSCHUNK_OVERLAP: // X
    case AS_DOVETAIL_CHORD_OVERLAP: // d
      TransChunk++;
      if(ScaffoldGraph->ignoreUOMTranschunk)
        return;
      if(uom_mesg->quality > transQualityCutoff){
        fprintf(stderr,"* Bad Quality overlap (" F_CID "," F_CID ",%c,%c) overlap " F_COORD "  quality: %g  .... ignored\n",
                uom_mesg->chunk1, 
                uom_mesg->chunk2,
                uom_mesg->orient,
                uom_mesg->overlap_type,
                uom_mesg->best_overlap_length,
                uom_mesg->quality);
        BadQuality++;
        return;
      }
      break;

#if 0
    case AS_1_CONTAINS_2_STACK_OVERLAP: // Z
    case AS_CONTAINMENT_CHORD_OVERLAP: // c
      ContainStack++;
      if(ScaffoldGraph->ignoreUOMContains
         || ScaffoldGraph->ignoreUOMContainStack)
        return;
      break;
#endif
    
    case AS_1_CONTAINS_2_OVERLAP: // C
    case AS_2_CONTAINS_1_OVERLAP: // I
      Containment++;
      if(ScaffoldGraph->ignoreUOMContains)
        return;
      break;
    case AS_NO_OVERLAP:  // N
    default:
      assert(0);

  }



  TotalUOMHandled++;

  if((TotalUOMHandled % 100000) == 0){
    fprintf(stderr,"* Processed %d UOM messages\n",
	    TotalUOMHandled);
    ReportMemorySize(ScaffoldGraph, stderr);
  }

  if(uom_mesg->chunk1 > uom_mesg->chunk2){
    int exchange_chunks = FALSE;
    cidA = uom_mesg->chunk2;
    cidB = uom_mesg->chunk1;
    exchange_chunks = TRUE;
    if(orient == AB_AB){
      orient = BA_BA;
    }else if(orient == BA_BA){
      orient = AB_AB;
    }

    if( exchange_chunks ) {
      // The containment edges are anti-symmetric in the unitig
      // overlap type.
      switch(overlap_type) {
        case AS_1_CONTAINS_2_OVERLAP:
          overlap_type = AS_2_CONTAINS_1_OVERLAP; 
          break;
        case AS_2_CONTAINS_1_OVERLAP:
          overlap_type = AS_1_CONTAINS_2_OVERLAP; 
          break;
        default:
          break;
      }
    }
  }

  // NEW
  /* Create the raw overlap edges in the graph */
  /* We create the edges here but we DEFER their insertion in the graph
     until all of the ChunkInstances have been created
  */
  {
    LengthT distance;
    ChunkInstanceT *chunkA, *chunkB;
    int hasContributingOverlap = !hasTandemOverlap;
    int hasTransChunk = FALSE;
    int hasRepeatOverlap = FALSE;
    int aContainsB = FALSE;
    int bContainsA = FALSE;

    hasTransChunk = overlap_type == AS_TRANSCHUNK_OVERLAP;

    hasRepeatOverlap = FALSE;

    distance.mean = -uom_mesg->best_overlap_length;

    if(hasTandemOverlap){
      distance.variance = TANDEM_OVERLAP_VARIANCE;
    }else{
      double rangeVarianceEstimate = (double)(range * range) / (36.0);
      distance.variance = MAX(1.0, MAX(rangeVarianceEstimate, ComputeFudgeVariance(-distance.mean)));
    }
    assert(distance.variance > 0.0);

    chunkA = GetGraphNode(ScaffoldGraph->CIGraph, cidA);
    chunkB = GetGraphNode(ScaffoldGraph->CIGraph, cidB);

#if 0
    // This code, not yet deployed, tries to collect the high water mark of tandem contamination on
    // the end of each unitig.  We use the two inline instance fields for this purpose.

    if(hasTandemOverlap){
      int aEnd, bEnd;
      CDS_COORD_t maxOverlap = uom_mesg->max_overlap_length;

      switch(orient){
        case AB_AB:
          aEnd = B_END;
          bEnd = A_END;
          break;
        case AB_BA:
          aEnd = B_END;
          bEnd = B_END;
          break;
        case BA_AB:
          aEnd = A_END;
          bEnd = A_END;
          break;
        case BA_BA:
          aEnd = A_END;
          bEnd = B_END;
          break;
        default:
          assert(0);
      }

      if(aEnd == A_END){
	chunkA->info.CI.instances.in_line.instance1 = 	MAX(chunkA->info.CI.instances.in_line.instance1, maxOverlap);
      }else{
	chunkA->info.CI.instances.in_line.instance2 = 	MAX(chunkA->info.CI.instances.in_line.instance2, maxOverlap);
      }
      if(bEnd == A_END){
	chunkB->info.CI.instances.in_line.instance1 = 	MAX(chunkB->info.CI.instances.in_line.instance1, maxOverlap);
      }else{
	chunkB->info.CI.instances.in_line.instance2 = 	MAX(chunkB->info.CI.instances.in_line.instance2, maxOverlap);
      }

    }
#endif

    switch(overlap_type){
      /* These are all of the dovetails */
      case AS_OVERLAP:          // O
      case AS_TOUCHES_CONTAINED_OVERLAP: // M
      case AS_TANDEM_OVERLAP:   // T
      case AS_BETWEEN_CONTAINED_OVERLAP: // Y
      case AS_TRANSCHUNK_OVERLAP: // X
      case AS_DOVETAIL_CHORD_OVERLAP: // d
        if(uom_mesg->best_overlap_length >  chunkB->bpLength.mean){
          if(GlobalData->verbose > 0)
            fprintf(stderr,"* Warning: non-contain overlap is really A contains B:\n\tidA " F_CID " lenA %g idB " F_CID " lenB %g overlap " F_COORD "\n",
                    cidA, chunkA->bpLength.mean,
                    cidB, chunkB->bpLength.mean,
                    uom_mesg->best_overlap_length);
	      
          aContainsB = TRUE;
          hasContributingOverlap = FALSE;
        }
        else if(   uom_mesg->best_overlap_length >  chunkA->bpLength.mean){
          if(GlobalData->verbose > 0)
            fprintf(stderr,"* Warning: non-contain overlap is really B contains A:\n\tidA " F_CID " lenA %g idB " F_CID " lenB %g overlap " F_COORD "\n",
                    cidA, chunkA->bpLength.mean,
                    cidB, chunkB->bpLength.mean,
                    uom_mesg->best_overlap_length);
          bContainsA = TRUE;
          hasContributingOverlap = FALSE;
        }

        if(overlap_type == AS_TANDEM_OVERLAP)
          hasContributingOverlap = FALSE;

        break;
	    
#if 0
      case AS_1_CONTAINS_2_STACK_OVERLAP: // Z
      case AS_CONTAINMENT_CHORD_OVERLAP: // c
#endif
      case AS_1_CONTAINS_2_OVERLAP: // C
        // make sure that the containment can possibly hold
        //	    assert(chunkA->bpLength.mean >= chunkB->bpLength.mean);
        //	    assert(uom_mesg->best_overlap_length >= chunkB->bpLength.mean);
        // If we are ingoring these, just return, don't add any edges
        hasContributingOverlap = FALSE;
        aContainsB = TRUE;
        hasContributingOverlap = FALSE;
        if( chunkA->bpLength.mean < chunkB->bpLength.mean ){
          fprintf(stderr,"* Warning: containing chunk " F_CID " (%g) is shorter than contained chunk " F_CID " (%g)\n",
                  cidA, chunkA->bpLength.mean,
                  cidB, chunkB->bpLength.mean);
        }
        break;

      case AS_2_CONTAINS_1_OVERLAP: // I
        // make sure that the containment can possibly hold
        //	    assert(chunkA->bpLength.mean <= chunkB->bpLength.mean);
        //    assert(uom_mesg->best_overlap_length <= chunkB->bpLength.mean);
        // If we are ingoring these, just return, don't add any edges
        hasContributingOverlap = FALSE;
        bContainsA = TRUE;
        if( chunkA->bpLength.mean > chunkB->bpLength.mean ){
          fprintf(stderr,"* Warning: containing chunk " F_CID " (%g) is shorter than contained chunk " F_CID " (%g)\n",
                  chunkB->id, chunkB->bpLength.mean,
                  chunkA->id, chunkA->bpLength.mean);
        }
        break;

      case AS_NO_OVERLAP:  // N
      default:
        assert(0);

    }
	    
#ifdef DEBUG_DETAILED
    fprintf(stderr,"* Read UOM (" F_CID "," F_CID ",%c) and added CIedge (" F_CID "," F_CID ",%c)\n",
	    uom_mesg->chunk1,
	    uom_mesg->chunk2,
	    uom_mesg->orient,
	    cidA,
	    cidB,
	    orient);
#endif


    if(uom_mesg->chunk1 == uom_mesg->chunk2){
      fprintf(stderr,
	      "* ???? UOM mesg specifies ck1:" F_CID " ck2:" F_CID " orient:%c type:%c...ignoring\n",
	      uom_mesg->chunk1,
	      uom_mesg->chunk2,
	      uom_mesg->orient,
	      uom_mesg->overlap_type);
      return;

    }



    hasTransChunk = (overlap_type == AS_TRANSCHUNK_OVERLAP);


    if(!hasTandemOverlap && !aContainsB && !bContainsA){
      hasRepeatOverlap = IsRepeatOverlap(ScaffoldGraph->CIGraph, cidA, 
					 cidB, orient,
					 distance);

      /*      if(repeatOverlap)
	      fprintf(stderr,"* Edge (" F_CID "," F_CID ") is a repeat overlap\n",
	      edge->idA, edge->idB);
      */
      hasContributingOverlap = !hasRepeatOverlap;
    }


    AddGraphEdge(ScaffoldGraph->CIGraph, 
                 cidA, 
                 cidB, 
                 NULLINDEX, NULLINDEX, // frags
                 NULLINDEX,  // dist
                 distance,
                 quality,
                 range/2,  // fudge
                 orient,
                 FALSE, // inducedByUnknownOrient,
                 FALSE, // hasGuide,
                 FALSE, // hasMayJoin,
                 FALSE, // hasMustJoin,
                 hasContributingOverlap,
                 hasRepeatOverlap, //   isRepeat
                 hasTandemOverlap,
                 aContainsB,
                 bContainsA,
                 hasTransChunk,
                 FALSE, // isExtremalA
                 FALSE, // isExtremalB
                 UNKNOWN_EDGE_STATUS,
                 FALSE,  // collectOverlap
                 FALSE /* FALSE */); // do NOT insert
  }
}

/***************************************************************************/


void ProcessIUM_ScaffoldGraph(IntUnitigMesg *ium_mesg,
                              CDS_COORD_t length, int sequenceOnly){
  CDS_CID_t cfr;
  CDS_COORD_t simLength;
  // Data for the ScaffoldGraph construction
  ChunkInstanceT CI;

#ifdef DEBUG_DETAILED
  fprintf(stderr,"* Unitig " F_IID " has ungapped length " F_COORD "\n",
	  ium_mesg->iaccession, length);
#endif
  CI.id = ium_mesg->iaccession;
  CI.bpLength.mean = length;
  CI.bpLength.variance = MAX(1.0,ComputeFudgeVariance(CI.bpLength.mean));
  CI.edgeHead = NULLINDEX;
  CI.microhetScore = NULLINDEX;
  CI.setID = NULLINDEX;
  CI.scaffoldID = NULLINDEX;
  CI.indexInScaffold = NULLINDEX;
  CI.prevScaffoldID = NULLINDEX;
  CI.numEssentialA = 0;
  CI.numEssentialB = 0;
  CI.essentialEdgeA = NULLINDEX;
  CI.essentialEdgeB = NULLINDEX;
  CI.smoothExpectedCID = NULLINDEX;
  CI.BEndNext = CI.AEndNext = NULLINDEX;
  CI.info.CI.headOfFragments = GetNumCIFragTs(ScaffoldGraph->CIFrags);
  CI.info.CI.numFragments = ium_mesg->num_frags;
  CI.info.CI.coverageStat = (ium_mesg->coverage_stat < -1000.0? -1000:ium_mesg->coverage_stat);
  CI.info.CI.branchPointA = ium_mesg->a_branch_point;
  CI.info.CI.branchPointB = ium_mesg->b_branch_point;
  CI.info.CI.contigID = NULLINDEX;
  CI.info.CI.numInstances = 0;
  CI.info.CI.instances.in_line.instance1 = 0;
  CI.info.CI.instances.in_line.instance2 = 0;
  CI.info.CI.instances.va = NULL;
  CI.info.CI.source = NULLINDEX;
  CI.flags.all = 0;
  CI.offsetAEnd.mean = 0.0;
  CI.offsetAEnd.variance = 0.0;
  CI.offsetBEnd = CI.bpLength;

#ifdef AS_ENABLE_SOURCE
  if(ium_mesg->source){
    char *c = ium_mesg->source;
    CI.info.CI.source = GetNumchars(ScaffoldGraph->SourceFields);
    while(*c != '\0'){
      Appendchar(ScaffoldGraph->SourceFields,c++);
    }
    Appendchar(ScaffoldGraph->SourceFields,c);
  }else{
    CI.info.CI.source = NULLINDEX;
  }
#endif

  // Collect the microhetScore if available
  CI.microhetScore = 1.01;
#ifdef AS_ENABLE_SOURCE
  {
    char *mhp = strstr(ium_mesg->source,"mhp:");
    if(mhp)
      CI.microhetScore = atof(mhp+4);
      //fprintf(stderr,"* %s\n*  mhp:%g found *\n", ium_mesg->source, CI.microhetScore);
  }
#endif

  // See if this is a repeat, or we can pin it down to an interval
  {
    char *interval;
    char *type;
    int result;
    //	  fprintf(stderr,"* source = %s\n", ium_mesg->source);
	  
    CI.flags.bits.cgbType = XX_CGBTYPE;
    CI.aEndCoord = CI.bEndCoord = -1;
    simLength = CI.bpLength.mean;

    // See if this is a repeat, or we can pin it down to an interval
#ifdef AS_ENABLE_SOURCE
    type = strstr(ium_mesg->source,"gen> ");
    if(type){
      type += 5;
      if(!strncmp(type,"uu",2) || !strncmp(type,"@@",2)){
	CI.flags.bits.cgbType = (unsigned int)UU_CGBTYPE;
      }else if(!strncmp(type,"ru",2)){
	CI.flags.bits.cgbType = (unsigned int)RU_CGBTYPE;
      }else if(!strncmp(type,"rr",2)){
	CI.flags.bits.cgbType = (unsigned int)RR_CGBTYPE;
      }else if(!strncmp(type,"ur",2)){
	CI.flags.bits.cgbType = (unsigned int)UR_CGBTYPE;
      }

      if((interval = strstr(ium_mesg->source,"["))){
	//	    fprintf(stderr,"* interval = %s\n", interval);
	result = sscanf(interval + 1," " F_COORD "," F_COORD,
			&CI.aEndCoord, &CI.bEndCoord);
	simLength = abs(CI.aEndCoord - CI.bEndCoord);
      }else{
	CI.aEndCoord = CI.bEndCoord = -1;
	simLength = CI.bpLength.mean;
      }
    }
#endif
  }

  if(ium_mesg->coverage_stat >= GlobalData->cgbUniqueCutoff){
    if(length < CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH ||
       ium_mesg->num_frags < CGW_MIN_READS_IN_UNIQUE){
      ShortDiscriminatorUniques++;
    }else{
      DiscriminatorUniques++;
    }
  }

  
  {
    int isUnique = FALSE;
    if(ium_mesg->coverage_stat >= GlobalData->cgbUniqueCutoff &&
       length >= CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH &&
       ium_mesg->num_frags >= CGW_MIN_READS_IN_UNIQUE){
      // microhetScore is actually the probability of the sequence
      // being UNIQUE, based on microhet considerations.
      // Falling below threshhold makes something a repeat.
      if( CI.microhetScore < GlobalData->cgbMicrohetProb){
	if(ium_mesg->coverage_stat < GlobalData->cgbApplyMicrohetCutoff){
	  fprintf(stderr,"* CI " F_CID " with astat: %g classified as repeat based on microhet unique prob of %g < %g\n",
		  CI.id, ium_mesg->coverage_stat, CI.microhetScore, GlobalData->cgbMicrohetProb);
	  isUnique = FALSE;
	  if(CI.flags.bits.cgbType == XX_CGBTYPE)
	    CI.flags.bits.cgbType = RR_CGBTYPE;
	  CI.type = UNRESOLVEDCHUNK_CGW;
	}else{
	  isUnique = TRUE;
	  fprintf(stderr,"* WARNING: CI " F_CID " with coverage %g WOULD HAVE BEEN classified as repeat based on microhet unique prob of %g < %g\n",
		  CI.id, ium_mesg->coverage_stat, CI.microhetScore, GlobalData->cgbMicrohetProb);
	}
      }else{
	isUnique = TRUE;
      }
    }else{
      isUnique = FALSE;
    }

    if(isUnique){

      ScaffoldGraph->numDiscriminatorUniqueCIs++;
      //      fprintf(stderr,"* CI " F_CID " is unique\n",CI.id);
      CI.flags.bits.isUnique = 1;
      CI.type = DISCRIMINATORUNIQUECHUNK_CGW;
      if(CI.flags.bits.cgbType == XX_CGBTYPE)
	CI.flags.bits.cgbType = UU_CGBTYPE;
    }else{
      CI.flags.bits.isUnique = 0;
      CI.type = UNRESOLVEDCHUNK_CGW;
      if(CI.flags.bits.cgbType == XX_CGBTYPE)
	CI.flags.bits.cgbType = RR_CGBTYPE;
    }
  }

  CI.flags.bits.smoothSeenAlready = FALSE;
  CI.flags.bits.isCI = TRUE;
  CI.flags.bits.tandemOverlaps = NO_TANDEM_OVERLAP;
  CI.flags.bits.isChaff = FALSE;


  if( ! sequenceOnly )
    {
      CDS_CID_t extremalA = NULLINDEX;
      CDS_CID_t extremalB = NULLINDEX;
      CDS_COORD_t minOffset = CDS_COORD_MAX;
      CDS_COORD_t maxOffset = CDS_COORD_MIN;
	  
      /* Determine extremal fragments so we can label the fragments */
	  
	  
      for(cfr = 0; cfr < ium_mesg->num_frags; cfr++){
	IntMultiPos *cfr_mesg = ium_mesg->f_list + cfr;
	CDS_COORD_t end = MAX( cfr_mesg->position.end, cfr_mesg->position.bgn);
	CDS_COORD_t beg = MIN( cfr_mesg->position.end, cfr_mesg->position.bgn);
	    
	if(minOffset > beg){
	  minOffset = beg;
	  extremalA = cfr;
	}
	if(maxOffset < end){
	  maxOffset = end;
	  extremalB = cfr;
	}
      }
	  
	  
      for(cfr = 0; cfr < ium_mesg->num_frags; cfr++){
	CIFragT        cifrag;
	InfoByIID  info, *old_info;
	CDS_CID_t fragid = GetNumCIFragTs(ScaffoldGraph->CIFrags);
	IntMultiPos *cfr_mesg = ium_mesg->f_list + cfr;

	cifrag.iid      = cfr_mesg->ident;
        cifrag.mateOf   = NULLINDEX;
        cifrag.dist     = 0;
	cifrag.cid      = ium_mesg->iaccession;
	cifrag.CIid     = ium_mesg->iaccession;

	// These get set in UpdateNodeFragments, called below
        cifrag.offset5p.mean      = 0.0;
        cifrag.offset5p.variance  = 0.0;
        cifrag.offset3p.mean      = 0.0;
        cifrag.offset3p.variance  = 0.0;

	cifrag.contigID                 = NULLINDEX;
        cifrag.contigOffset5p.mean      = 0.0;
        cifrag.contigOffset5p.variance  = 0.0;
        cifrag.contigOffset3p.mean      = 0.0;
        cifrag.contigOffset3p.variance  = 0.0;

        cifrag.type      = cfr_mesg->type;
        cifrag.label     = AS_SINGLETON;

	cifrag.flags.bits.hasInternalOnlyCILinks     = FALSE; // set in CreateCIEdge
	cifrag.flags.bits.hasInternalOnlyContigLinks = FALSE; // set in CreateCIEdge
	cifrag.flags.bits.isPlaced                   = FALSE;
	cifrag.flags.bits.isSingleton                = FALSE;
	cifrag.flags.bits.isChaff                    = FALSE;
        cifrag.flags.bits.innieMate                  = FALSE;
        cifrag.flags.bits.hasMate                    = FALSE;
        cifrag.flags.bits.linkType                   = AS_UNKNOWN;
	cifrag.flags.bits.edgeStatus                 = INVALID_EDGE_STATUS;
        cifrag.flags.bits.mateDetail                 = UNASSIGNED_MATE;

        //  Singleton chunks are chaff; singleton frags are chaff
        //  unless proven otherwise
        //
        if (ium_mesg->num_frags < 2) {
          CI.flags.bits.isChaff         = TRUE;
          cifrag.flags.bits.isSingleton = TRUE;
          cifrag.flags.bits.isChaff     = TRUE;
	}

        cifrag.locale         = NULLINDEX;
        cifrag.localePos.bgn  = 0;
        cifrag.localePos.end  = 0;

	info.fragIndex   = fragid;
	info.set         = TRUE;

	// Check to see if we've already seen this fragment by IID!!!!
	old_info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, cifrag.iid);
	if(old_info && old_info->set){
	  CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, old_info->fragIndex);
	  fprintf(stderr,"*** FATAL ERROR:  Fragment with IID " F_CID " appears more than once with id " F_CID " and " F_CID "\n",
                  cifrag.iid, old_info->fragIndex, fragid);
	  fprintf(stderr,"***               First appearance was in unitig " F_CID ", currently found in unitig " F_CID "\n",
		  frag->cid, cifrag.cid);
	  exit(1);
	}

	SetInfoByIID(ScaffoldGraph->iidToFragIndex, cifrag.iid, &info);
	    
	// Collect guide stats
        if(AS_FA_GUIDE(cfr_mesg->type)){ 
	  totalGuideFrags++;
	  if(CI.flags.bits.isUnique)
	    inUniqueGuideFrags++;
          else
	    inRepeatGuideFrags++;
	  if(ium_mesg->num_frags <= 2)
	    inTeenyUnitigGuideFrags++;
	  if(ium_mesg->num_frags < 2)
	    inSingletonUnitigGuideFrags++;
	  if(cfr == extremalA || cfr == extremalB)
	    onEndGuideFrags++;
	}
	    
	// Collect read stats
        if (AS_FA_READ(cfr_mesg->type)) { 
	  totalReadFrags++;
	  if(CI.flags.bits.isUnique)
	    inUniqueReadFrags++;
	  else
	    inRepeatReadFrags++;
	  if(ium_mesg->num_frags <= 2)
            inTeenyUnitigReadFrags++;
          if(ium_mesg->num_frags < 2)
	    inSingletonUnitigReadFrags++;
	  if(cfr == extremalA || cfr == extremalB)
	    onEndReadFrags++;
	}
        //else if(AS_FA_SHREDDED(cfr_mesg->type)){ 
	//  CI.flags.bits.includesFinishedBacFragments = TRUE; 
        //}
	   
        AppendCIFragT(ScaffoldGraph->CIFrags, &cifrag);
      }
    }

  // Insert the Chunk Instance
  SetChunkInstanceT(ScaffoldGraph->CIGraph->nodes, CI.id, &CI);

  // Mark all frags as being members of this CI, and set their offsets within
  // the CI
  if( ! sequenceOnly )
    UpdateNodeFragments(ScaffoldGraph->CIGraph,CI.id, CI.type == DISCRIMINATORUNIQUECHUNK_CGW, TRUE ); // mark unitigs and contigs
}






/***************************************************************************/
void ProcessIUM(IntUnitigMesg *ium_mesg){
  MultiAlignT *ma = CreateMultiAlignTFromIUM(ium_mesg, GetNumCIFragTs(ScaffoldGraph->CIFrags),FALSE);
  CDS_COORD_t length = GetMultiAlignUngappedLength(ma);

  //  AddReferenceMultiAlignT(ma);

  //  SetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, ium_mesg->iaccession, ma);
  //  SetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, ium_mesg->iaccession, ma);


  // This adds a reference only if keepInCache is true...
  InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg->iaccession, TRUE, ma, TRUE);
  DuplicateEntryInSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg->iaccession, TRUE, ium_mesg->iaccession, FALSE, FALSE);


  ProcessIUM_ScaffoldGraph(ium_mesg, length, FALSE);

  UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg->iaccession, TRUE);
}


/****************************************************************************/



void LoadDistData(void){ // Load the distance record info from the gkpStore
  int32 numDists = getNumGateKeeperLibrarys(ScaffoldGraph->gkpStore->lib);
  CDS_CID_t i;
  
  for(i = 1; i <= numDists; i++){
    DistT dist;
    GateKeeperLibraryRecord gkpl;
    getGateKeeperLibraryStore(ScaffoldGraph->gkpStore->lib, i, &gkpl);
    
    if (gkpl.deleted)
      continue;

    dist.mu             = gkpl.mean;
    dist.sigma          = gkpl.stddev;
    dist.numSamples     = 0;
    dist.min            = CDS_COORD_MAX;
    dist.max            = CDS_COORD_MIN;
    dist.bnum           = 0;
    dist.bsize          = 0;
    dist.histogram      = NULL;
    dist.lower          = dist.mu - CGW_CUTOFF * dist.sigma;
    dist.upper          = dist.mu + CGW_CUTOFF * dist.sigma;
    dist.numReferences  = 0;
    dist.numBad         = 0;

    fprintf(GlobalData->stderrc,"* Loaded dist "F_UID","F_CID" (%g +/- %g)\n",
            gkpl.libraryUID, i, dist.mu, dist.sigma);

    SetDistT(ScaffoldGraph->Dists, i, &dist);
  }
}



void NullifyNodeEdges(NodeCGW_T * node)
{
  node->numEssentialA = node->numEssentialB = 0;
  node->essentialEdgeA = node->essentialEdgeB = NULLINDEX;
  node->edgeHead = NULLINDEX;
  if(node->flags.bits.isScaffold)
    {
      node->info.Scaffold.internalEdges =
        node->info.Scaffold.confirmedInternalEdges = 0;
    }
  else
    {
      node->scaffoldID = NULLINDEX;
      node->prevScaffoldID = NULLINDEX;
      node->indexInScaffold = NULLINDEX;
      if(node->flags.bits.isContig)
        {
          node->BEndNext = node->AEndNext = NULLINDEX;
          node->smoothExpectedCID = NULLINDEX;
        }
    }
}

void NullifyAllNodeEdges(GraphCGW_T * graph)
{
  GraphNodeIterator nodes;
  NodeCGW_T * node;
  
  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&nodes)))
    {
      NullifyNodeEdges(node);
    }
}


/****************************************************************************/
void ProcessFrags(void)
{
  CDS_CID_t i;
  int32 unmatedFrags = 0;
  GateKeeperFragmentRecord gkf;

  //  Do one pass through, reading from the gatekeeper store to fill
  //  out the cifrag info.

  for(i = 0; i < GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex); i++){
    InfoByIID *ciinfo = GetInfoByIID(ScaffoldGraph->iidToFragIndex, i);
    CIFragT   *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, ciinfo->fragIndex);

    if(!ciinfo->set)
      continue;

    assert(cifrag->iid == i);  //  If !set, this fails.

    getGateKeeperFragmentStore(ScaffoldGraph->gkpStore->frg, i, &gkf);

    if (gkf.mateIID != 0) {
      InfoByIID *miinfo = GetInfoByIID(ScaffoldGraph->iidToFragIndex, gkf.mateIID);

      if (miinfo->set) {
        cifrag->mateOf   = miinfo->fragIndex;
        cifrag->dist     = gkf.libraryIID;
        if (gkf.orientation == AS_READ_ORIENT_INNIE)
          cifrag->flags.bits.innieMate = TRUE;
        cifrag->flags.bits.linkType   = AS_MATE;
        cifrag->flags.bits.edgeStatus = UNKNOWN_EDGE_STATUS;
        cifrag->flags.bits.hasMate    = TRUE;
      } else {
        fprintf(stderr, "ProcessFrags()-- WARNING!  fragiid=%d,index=%d mateiid=%d,index=%d -- MATE DOESN'T EXIST!\n",
                i, ciinfo->fragIndex, gkf.mateIID, miinfo->fragIndex);
        //  This is not a critical failure, but does indicate
        //  something amiss with either the store or the unitigs.
        //
        assert(0);
      }

      //fprintf(stderr, "Frag: iid=%d,index=%d mateiid=%d,index=%d\n", i, ciinfo->fragIndex, gkf.mateIID, miinfo->fragIndex);
    }

    if (cifrag->flags.bits.hasMate == FALSE)
      unmatedFrags++;
  }

  fprintf(stderr,"* Unmated fragments %d\n", unmatedFrags);
  fprintf(stderr,"* Total IIDs:       %d\n", (int)GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex));

  //  Do a second pass to clean up mates.  We don't need to access the
  //  gatekeeper here, since everything is already populated.  If we
  //  put this stuff in the first loop, we'll be randomly accessing
  //  the gatekeeper.

  for(i = 0; i < GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex); i++){
    InfoByIID *ciinfo = GetInfoByIID(ScaffoldGraph->iidToFragIndex, i);
    CIFragT   *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, ciinfo->fragIndex);
    InfoByIID *miinfo = NULL;
    CIFragT   *mifrag = NULL;

    //  this frag not used, or no mate
    if ((!ciinfo->set) || (cifrag->mateOf == NULLINDEX))
      continue;

    mifrag = GetCIFragT(ScaffoldGraph->CIFrags, cifrag->mateOf);
    miinfo = GetInfoByIID(ScaffoldGraph->iidToFragIndex, mifrag->iid);

    if ((mifrag == NULL) || (mifrag->mateOf == NULLINDEX)) {
      //  We set up links to a dead frag, clean up...

      cifrag->mateOf   = NULLINDEX;
      cifrag->dist     = NULLINDEX;
      cifrag->flags.bits.linkType   = AS_UNKNOWN;
      cifrag->flags.bits.edgeStatus = INVALID_EDGE_STATUS;
      cifrag->flags.bits.hasMate    = FALSE;

      if (mifrag) {
        mifrag->mateOf   = NULLINDEX;
        mifrag->dist     = NULLINDEX;
        mifrag->flags.bits.linkType   = AS_UNKNOWN;
        mifrag->flags.bits.edgeStatus = INVALID_EDGE_STATUS;
        cifrag->flags.bits.hasMate    = FALSE;
      }
    } else {
      //  Both guys are alive, and we're mated.  Throw some asserts

      assert(ciinfo->set);
      assert(miinfo->set);
      assert(cifrag->dist   == mifrag->dist);
      assert(cifrag->mateOf == miinfo->fragIndex);
      assert(mifrag->mateOf == ciinfo->fragIndex);
    }

  }  //  for each frag

  {
    double numCIFrags = (double) GetNumCIFragTs(ScaffoldGraph->CIFrags);
    if (0 < numCIFrags) 
      fprintf(stderr,"* Found %d unmated frags (%g %%)\n",
  	      unmatedFrags, 100. * (double)unmatedFrags / numCIFrags);
    else
      fprintf(stderr, "* Found %d unmated frags\n", unmatedFrags);
  }
  
}


/****************************************************************************/

  
static int CompareEdgesByIdA_IdB_Orientation_Quality(const void *c1, const void *c2){
  EdgeCGW_T *e1 = (EdgeCGW_T *)c1;
  EdgeCGW_T *e2 = (EdgeCGW_T *)c2;
  int diff;
  float32 qualityDiff;

  diff = e1->idA - e2->idA;
  if(diff)
    return diff;

  diff = e1->idB - e2->idB;
  if(diff)
    return diff;
  
  diff = e1->orient - e2->orient;
  if(diff)
    return diff;
  
  qualityDiff = (e1->quality - e2->quality); // If e1 has larger quality, put it at the end
  if(qualityDiff >= 0)
    return 1;
  else
    return -1;
}

/********************************************/
/* This routine sorts the edges such that they are in reversed order
   by idA, idB, and orientation.  This should speed their insertion in the graph,
   and subsequent traversals of the graph to some degree */

static int CompareEdgesBeforeInsertion(const void *c1, const void *c2){
  EdgeCGW_T *e1 = (EdgeCGW_T *)c1;
  EdgeCGW_T *e2 = (EdgeCGW_T *)c2;
  int diff;
  
  diff = e2->idA - e1->idA;
  if(diff)
    return diff;

  diff = e2->idB - e1->idB;
  if(diff)
    return diff;

  diff = e2->orient - e1->orient;
  if(diff)
    return diff;

  return 0;
}

/********************************************/
static int CompareEdgesByIdA_Quality(const void *c1, const void *c2){
  EdgeCGW_T *e1 = (EdgeCGW_T *)c1;
  EdgeCGW_T *e2 = (EdgeCGW_T *)c2;
  int diff;
  float32 qualityDiff;

  diff = e1->idA - e2->idA;
  if(diff)
    return diff;
  
  qualityDiff = (e1->quality - e2->quality); // If e1 has larger quality, put it at the end
  if(qualityDiff >= 0)
    return 1;
  else
    return -1;
}

/********************************************/
static int CompareEdgesByIdB_Quality(const void *c1, const void *c2){
  EdgeCGW_T *e1 = (EdgeCGW_T *)c1;
  EdgeCGW_T *e2 = (EdgeCGW_T *)c2;
  int diff;
  float32 qualityDiff;


  diff = e1->idB - e2->idB;
  if(diff)
    return diff;
  
  qualityDiff = (e1->quality - e2->quality); // If e1 has larger quality, put it at the end
  if(qualityDiff >= 0)
    return 1;
  else
    return -1;
}

/********************************************/
int  EliminateDuplicates(void)
{
  int32 numEdges = 0;
  CDS_CID_t idA = NULLINDEX;
  CDS_CID_t idB = NULLINDEX;
  ChunkOrientationType orient = XX_XX;
  int32 numSeen = 0;
  CDS_CID_t i;
  EdgeCGW_T *edges = GetGraphEdge(ScaffoldGraph->CIGraph,0);
  EdgeCGW_T *headEdge = NULL;

  for(i = 0; i < GetNumGraphEdges(ScaffoldGraph->CIGraph);i++){
    EdgeCGW_T *edge = &edges[i];

    if(edge->idA != idA ||
       edge->idB != idB ||
       edge->orient != orient){
      headEdge = edge;
      numSeen = 1;
      idA = edge->idA;
      idB = edge->idB;
      orient = edge->orient;
    }else{
      numSeen++;
    }
    if(numSeen > 1){
      /* We are sorted in order of decreasing quality.  So the first edge is the best. */
      edge->flags.bits.isDeleted = TRUE;
#ifdef DEBUG_DETAILED
      fprintf(stderr,"* Deleted duplicate edge \n");
      PrintGraphEdge(stderr,ScaffoldGraph->CIGraph,"head ", headEdge, idA);
      PrintGraphEdge(stderr,ScaffoldGraph->CIGraph,"dupl ", edge, idA);
#endif	
      numEdges++;
    }
      
  }
  return numEdges;
}


int  MarkTopEdges(int32 maxDegree, int32 maxDegreeUnique, int useIDA)
{
  int32 numEdges = 0;
  int32 numAEnd = 0;
  int32 numBEnd = 0;
  int32 maxEndDegree = 0;
  NodeCGW_T *node = NULL;
  int nodeIsUnique = FALSE;
  int32 i;
  int32 maxEdges = (int32) GetNumGraphEdges(ScaffoldGraph->CIGraph);

  fprintf(stderr,"* MarkTopEdges maxDegree:%d useIDA:%d edges:%d\n",
          maxDegree, useIDA, maxEdges);

  for(i = 0; i < maxEdges;i++){
    EdgeCGW_T *edge = GetGraphEdge(ScaffoldGraph->CIGraph, i);
    CDS_CID_t nodeID = (useIDA?edge->idA:edge->idB);
    int markEdge = FALSE;

    if(edge->flags.bits.isDeleted)
      continue;

    //      PrintGraphEdge(stderr,ScaffoldGraph->CIGraph," ", edge, nodeID);

    if(!node || node->id != nodeID){
      assert(!node || nodeID > node->id);
      node = GetGraphNode(ScaffoldGraph->CIGraph, nodeID);
      numAEnd = 0;
      numBEnd = 0;
      nodeIsUnique = (node->info.CI.coverageStat > 0);
      maxEndDegree = maxDegree;
      if(nodeIsUnique){   /* Keep more edges for 'unique' nodes */
        maxEndDegree = maxDegreeUnique;
      }
#ifdef DEBUG_DETAILED
      fprintf(stderr,"* Node " F_CID " %s maxEndDegree = %d\n",
              node->id, (nodeIsUnique?" Unique ":" Not Unique "), maxEndDegree);
#endif
    }
    
    switch (GetEdgeOrientationWRT(edge,nodeID)){
      case AB_BA:
      case AB_AB:
        if(numBEnd++ < maxEndDegree)
          markEdge = TRUE;
        break;
      case BA_BA:
      case BA_AB:
        if(numAEnd++ < maxEndDegree)
          markEdge = TRUE;
        break;
      default:
        assert(0);
    }

    if(useIDA)
      edge->flags.bits.highQualityA = markEdge;
    else
      edge->flags.bits.highQualityB = markEdge;

    if(edge->flags.bits.highQualityA ||
       edge->flags.bits.highQualityB){
      numEdges++;
    }else{
#if 0
      if(!useIDA)
        PrintGraphEdge(stderr,ScaffoldGraph->CIGraph," Low Quality ", edge, nodeID);
#endif
    }
	  

  }
  return numEdges;
}

void  CopyActiveEdges(int32 activeEdges){
  VA_TYPE(EdgeCGW_T) *newEdges = CreateVA_EdgeCGW_T(activeEdges);
  EdgeCGW_T *edges;
  int32 i;
  int32 numInserted = 0;

  fprintf(stderr,"* CopyActiveEdges activeEdges = %d*\n", activeEdges);
  fflush(NULL);
  AssertPtr(newEdges);

  for(i = 0; i < GetNumGraphEdges(ScaffoldGraph->CIGraph); i++){
    EdgeCGW_T *edge = GetGraphEdge(ScaffoldGraph->CIGraph, i);
    if(!edge->flags.bits.isDeleted && 
       (edge->flags.bits.highQualityA ||
        edge->flags.bits.highQualityB)){
      AppendEdgeCGW_T(newEdges,edge);
      numInserted++;
    }
  }

  assert(numInserted == activeEdges);

  fprintf(stderr,"* Sorting *\n");
  fflush(NULL);
  // Sort edges by idA, idB, orient, quality, in reverse order to speed graph insertion
  edges = GetEdgeCGW_T(newEdges,0);

  if(edges){
    qsort((void *)edges, activeEdges, sizeof(EdgeCGW_T),CompareEdgesBeforeInsertion);
  }

  DeleteVA_EdgeCGW_T(ScaffoldGraph->CIGraph->edges);
  ScaffoldGraph->CIGraph->edges = newEdges;
  ScaffoldGraph->CIEdges = newEdges;

  fprintf(stderr,"* Inserting *\n");
  fflush(NULL);

  for(i = 0; i < GetNumGraphEdges(ScaffoldGraph->CIGraph); i++){
    EdgeCGW_T *edge = GetGraphEdge(ScaffoldGraph->CIGraph, i);
    assert(!edge->flags.bits.isDeleted);
    edge->topLevelEdge = i;
    InsertGraphEdge(ScaffoldGraph->CIGraph, i, FALSE);
    CreateChunkOverlapFromEdge(ScaffoldGraph->CIGraph, edge, TRUE);
  }

  fprintf(stderr,"* Inserted %d edges\n", numInserted);
  fflush(NULL);
}


void ProcessEdges(int32 maxDegree, int32 maxDegreeUnique){
  EdgeCGW_T *edges = GetGraphEdge(ScaffoldGraph->CIGraph,0);
  int32 numEdges;

  fprintf(stderr,"* Sort 1\n");
  fflush(NULL);

  // Sort edges by idA, idB, orient, quality
  qsort((void *)edges, GetNumGraphEdges(ScaffoldGraph->CIGraph), sizeof(EdgeCGW_T),CompareEdgesByIdA_IdB_Orientation_Quality);

  fprintf(stderr,"* Eliminate Duplicates\n");
  fflush(NULL);
  // Eliminate Duplicates..pick highest quality edge
  numEdges = EliminateDuplicates();
  fprintf(stderr,"* deleted %d duplicate edges\n", numEdges);

  fprintf(stderr,"* Sort 2\n");
  fflush(NULL);
  // Sort edges by idA, end, quality
  qsort((void *)edges, GetNumGraphEdges(ScaffoldGraph->CIGraph), sizeof(EdgeCGW_T), CompareEdgesByIdA_Quality);
  
  // Mark top N edges incident on non-unique ends
  numEdges = MarkTopEdges(maxDegree, maxDegreeUnique, TRUE /* idA */);
  fprintf(stderr," %d edges are marked to keep\n",
	  numEdges);

  // Sort edges by idB, end, quality
  qsort((void *)edges, GetNumGraphEdges(ScaffoldGraph->CIGraph), sizeof(EdgeCGW_T), CompareEdgesByIdB_Quality);

  // Mark top N edges incident on non-unique ends
  numEdges = MarkTopEdges(maxDegree, maxDegreeUnique, FALSE /* idB */);
  fprintf(stderr," %d edges are marked to keep\n",
	  numEdges);

  // Iterate through all edges, Copy edges that are not double marked to a new array
  // Delete original array
  // link new array to CIGraph->Edges
  // Insert all edges in graph
  // Add all edges to hashtable


  CopyActiveEdges(numEdges);



}

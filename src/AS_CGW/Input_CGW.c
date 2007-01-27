
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
static char CM_ID[] = "$Id: Input_CGW.c,v 1.16 2007-01-27 00:30:10 brianwalenz Exp $";

/*   THIS FILE CONTAINS ALL PROTO/IO INPUT ROUTINES */


//#define DEBUG 1
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
#include "AS_UTL_Var.h"
#include "AS_UTL_version.h"
#include "AS_CGW_dataTypes.h"
#include "AS_PER_gkpStore.h"
#include "ScaffoldGraph_CGW.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
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

/****  Temporary  ****/
//static VA_TYPE(InternalLinkMesg) *ILKs;


#if 0
static void ProcessIDT(InternalDistMesg *idt_mesg);
#endif
static void ProcessUOM(UnitigOverlapMesg *uom_mesg, float transQualityCutoff);
static void ProcessIUM(IntUnitigMesg *ium_mesg);
static void ProcessFrags(void);
static void ProcessEdges(int32 maxDegree, int32 maxDegreeUnique);


static int32 msgCnt = 0;
static int32 numIDT = 0, numILK = 0, numILKReread = 0, numIUM = 0, numUOM = 0, numJNC = 0;

/****************************************************************************/
int ProcessInput(Global_CGW *data, int optind, int argc, char *argv[]){
  GenericMesg   *pmesg;
  FILE *infp;
  int i,j = 0;

  StartTimerT(&GlobalData->InputTimer);

  for(i = optind; i < argc; i++){
    fprintf(stderr,"* Opening file %d %s\n", j ++, argv[i]);
    infp = fopen(argv[i],"r");
    data->reader = (MesgReader)InputFileType_AS(infp);

    while(  (EOF != (data->reader)(infp, &pmesg))){

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
        case MESG_ILK:
          numILK++;
          break;
        case MESG_IBA:
        case MESG_FOM:
          break;

        default:
          fprintf(stderr,"* Oops: Read Message with type = %d\n", pmesg->t);
          (data->writer)(data->cgwfp,pmesg);      // Echo to output
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
  if(numILK)
    fprintf(stderr,"\tILK:%d  (%d rereads ignored)(ignored)\n",  numILK, numILKReread);
  if(numIDT)
    fprintf(stderr,"\tIDT:%d  (ignored...read from the gatekeeper store)\n",  numIDT);
  if(numJNC)
    fprintf(stderr,"\tJNC:%d  \n",  numJNC);
  
  /***vvvvv********  Process UOMs **************************vvvvvvvvvvv******/
  ProcessEdges(data->maxDegree, data->maxDegreeUnique);

  /***vvvvv********  Process Frags *************************vvvvvvvvvvv******/
  ProcessFrags();
  //  DeleteVA_InternalLinkMesg(ILKs);
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
  (data->reader)(infp, &pmesg);  // Read 1 message

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
    
  (data->writer)(data->cgwfp,pmesg);        // Echo to output

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
                 FALSE, // hasSTSGuide,
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
#ifdef DEBUG_DATA
  CI.info.CI.source = NULLINDEX;
#endif
  CI.flags.all = 0;
  CI.offsetAEnd.mean = 0.0;
  CI.offsetAEnd.variance = 0.0;
  CI.offsetBEnd = CI.bpLength;

#ifdef DEBUG_DATA
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
  {
    char *mhp;
    char *score;

    mhp = strstr(ium_mesg->source,"mhp:");
    if(mhp){
      score = mhp+4;
      CI.microhetScore = atof(score);
#if 0
      fprintf(stderr,"* %s\n*  mhp:%g found *\n", ium_mesg->source, CI.microhetScore);
#endif

    }else{
      CI.microhetScore = 1.01;
    }

  }
  // See if this is a repeat, or we can pin it down to an interval
  {
    char *interval;
    char *type;
    int result;
    //	  fprintf(stderr,"* source = %s\n", ium_mesg->source);
	  
    // See if this is a repeat, or we can pin it down to an interval
    type = strstr(ium_mesg->source,"gen> ");
    CI.flags.bits.cgbType = (unsigned int)XX_CGBTYPE;
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
    }else{
      CI.aEndCoord = CI.bEndCoord = -1;
      simLength = CI.bpLength.mean;

    }
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
  CI.flags.bits.includesFinishedBacFragments = FALSE;
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
	cifrag.iid = cfr_mesg->ident;
	cifrag.cid = ium_mesg->iaccession;
	cifrag.CIid = ium_mesg->iaccession;
	cifrag.contigID = NULLINDEX;
	info.fragIndex = fragid;
	info.set = TRUE;
	    
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
	    
	//	cifrag.nextCIFrag = fragid + 1;   // obsolete
	cifrag.label = AS_SINGLETON; // interim


        if(ium_mesg->num_frags < 2){
          CI.flags.bits.isChaff = TRUE;
          if(GlobalData->debugLevel > 0)
            fprintf(stderr,"* Singleton chunk " F_CID " is CHAFF\n", CI.id);
        }



	// Collect guide stats
        if(AS_FA_GUIDE(cfr_mesg->type)){ 
	  totalGuideFrags++;
	  if(CI.flags.bits.isUnique){
	    inUniqueGuideFrags++;
	  }else{
	    inRepeatGuideFrags++;
	  }
	  if(ium_mesg->num_frags <= 2)
	    inTeenyUnitigGuideFrags++;
	      
	  if(ium_mesg->num_frags < 2)
	    inSingletonUnitigGuideFrags++;
	      
	      
	  if(cfr == extremalA || cfr == extremalB){
	    onEndGuideFrags++;
	  }
	}
	    
	// Collect read stats
        if (AS_FA_READ(cfr_mesg->type)) { 
	  totalReadFrags++;
	  if(CI.flags.bits.isUnique){
	    inUniqueReadFrags++;
	  }else{
	    inRepeatReadFrags++;
	  }
	  if(ium_mesg->num_frags <= 2)
	    inTeenyUnitigReadFrags++;
	      
	  if(ium_mesg->num_frags < 2){
	    inSingletonUnitigReadFrags++;
	  }
	  if(cfr == extremalA || cfr == extremalB){
	    onEndReadFrags++;
	  }
	      
	      
	}else if(AS_FA_SHREDDED(cfr_mesg->type)){ 
	  CI.flags.bits.includesFinishedBacFragments = TRUE; 
	}
	    
	cifrag.locale = NULLINDEX;
	cifrag.localePos.bgn = cifrag.localePos.end = 0;
	cifrag.mateOf = NULLINDEX;
	cifrag.flags.bits.mateStatus = MATE_NONE;
	cifrag.type  = cfr_mesg->type;
	cifrag.linkType = 0;
#ifdef DEBUG_DATA
	cifrag.source = NULLINDEX;
#endif
	cifrag.numLinks = -1; // not set yet
	cifrag.flags.bits.getLinksFromStore = TRUE;  // get links from store
	cifrag.flags.bits.hasInternalOnlyCILinks = FALSE; // set in CreateCIEdge
	cifrag.flags.bits.hasInternalOnlyContigLinks = FALSE; // set in CreateCIEdge
	cifrag.flags.bits.hasFalseMate = FALSE;
	cifrag.flags.bits.edgeStatus = INVALID_EDGE_STATUS;
	cifrag.flags.bits.isPlaced = FALSE;
	cifrag.flags.bits.isSingleton = FALSE;
	cifrag.flags.bits.isChaff = FALSE;
	// These get set in UpdateNodeFragments, called below
	cifrag.offset3p.mean  =  cifrag.offset5p.mean = 0.0;
	cifrag.offset3p.variance  =  cifrag.offset5p.variance = 0.0;
	    

	// Singleton frags are chaff unless proven otherwise
	if(ium_mesg->num_frags == 1){
	  cifrag.flags.bits.isSingleton = TRUE;
	  cifrag.flags.bits.isChaff = TRUE;
	}
	    
	// cifrag.fragOrient = getFragOrient(&cifrag);
	    
#ifdef DEBUG_DATA
	{
          cifrag.aEndCoord = cifrag.bEndCoord = -1;
	}
        cifrag.source = NULLINDEX;
#endif	    
	    
	AppendCIFragT(ScaffoldGraph->CIFrags, &cifrag);
	    
      }
#if 0
      // obsolete
      /* Terminate the linked list */
      {
	CIFragT *last = GetCIFragT(ScaffoldGraph->CIFrags,
                                   GetNumCIFragTs(ScaffoldGraph->CIFrags) - 1);
	last->nextCIFrag = NULLINDEX;
      }
#endif
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
#if 0
void ProcessIDT(InternalDistMesg *idt_mesg)
{


  switch(idt_mesg->action){
    case AS_DELETE:
      return;
      break; // Nothing to do...they won't be referenced
    case AS_ADD:
    case AS_REDEFINE:
      {
        DistT dist;
        dist.mean = idt_mesg->mean;
        dist.stddev = idt_mesg->stddev;
        dist.mu = dist.sigma = 0.0;
        dist.min = CDS_COORD_MAX;
        dist.max = CDS_COORD_MIN;
        dist.bsize = 0;
        dist.histogram = NULL;
        dist.lower = dist.upper = 0;
        dist.samples = CreateVA_CDS_COORD_t(1024);
        dist.numReferences = dist.numBad = 0;
        SetDistT(ScaffoldGraph->Dists, idt_mesg->iaccession, &dist);
      }
      break;
    default:
      assert(0);
  }
}
#endif



void LoadDistData(void){ // Load the distance record info from the gkpStore
  int32 numDists = getNumGateKeeperDistances(ScaffoldGraph->gkpStore.dstStore);
  CDS_CID_t i;
  
  for(i = 1; i <= numDists; i++){
    DistT dist;
    GateKeeperDistanceRecord gkpd;
    getGateKeeperDistanceStore(ScaffoldGraph->gkpStore.dstStore, i, &gkpd);
    
    if(gkpd.deleted)
      continue;

    dist.mean = gkpd.mean;
    dist.stddev = gkpd.stddev;
    dist.mu = dist.sigma = 0.0;
    dist.min = CDS_COORD_MAX;
    dist.max = CDS_COORD_MIN;
    dist.bsize = 0;
    dist.histogram = NULL;
    dist.lower = dist.upper = 0;
    dist.samples = CreateVA_CDS_COORD_t(1024);
    dist.numReferences = dist.numBad = 0;

    if(GlobalData->verbose)
      fprintf(GlobalData->stderrc,"* Loaded dist " F_CID " (%g +/- %g)\n", i, gkpd.mean, gkpd.stddev);

    SetDistT(ScaffoldGraph->Dists, i, &dist);
  }
}


void  LoadLocaleData(void){ // Load locale info from fragStore
  ReadStructp rsp = new_ReadStruct();
  CDS_CID_t i;
  
  fprintf(stderr,"* Loading Locale Data from FragStore *\n");

  for(i = 0; i < GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex); i++){
    InfoByIID *infobyiid = GetInfoByIID(ScaffoldGraph->iidToFragIndex, i);
    if(infobyiid->set){
      cds_uint32 bgn,end;

      CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, infobyiid->fragIndex);

      if(!AS_FA_SHREDDED(frag->type)) 
	continue;

      getFragStore(ScaffoldGraph->fragStore, i, FRAG_S_FIXED, rsp);
      getLocalePos_ReadStruct(rsp, &bgn, &end);
      frag->localePos.bgn = bgn;
      frag->localePos.end = end;

#ifdef DEBUG_DETAILED
      fprintf(stderr,"* Read locale " F_CID " [" F_COORD "," F_COORD "] for fragment id " F_CID " (iid " F_CID ")\n",
	      frag->locale, frag->localePos.bgn, frag->localePos.end,
	      infobyiid->fragIndex, i);
#endif
    }

  }

  delete_ReadStruct(rsp);

  fprintf(stderr,"* Done! *\n");
  fflush(NULL);
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

void ActivateLBACMatePairs(void)
{
  int32 i;
  int32 unMated = 0;
  int32 multiMated = 0;
  int32 singleMated = 0;
  int32 lBACFrags = 0;

  // loop to identify mated LBAC fragments
  for(i = 0; i < GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex); i++)
    {
      InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex,i);
      CIFragT *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
      GateKeeperFragmentRecord gkf;
    
      // assert(cifrag->iid == i); 
      getGateKeeperFragmentStore(ScaffoldGraph->gkpStore.frgStore, cifrag->iid, &gkf);
      if(gkf.type != AS_LBAC)
        continue;
      /*
        if(gkf.type != AS_EBAC)
        continue;
      */

      lBACFrags++;
      cifrag->linkHead = gkf.linkHead;
      cifrag->numLinks = gkf.numLinks;

      assert((cifrag->linkHead == NULL_LINK && gkf.numLinks == 0) ||
             (cifrag->linkHead != NULL_LINK && gkf.numLinks > 0));

      if(gkf.numLinks == 0)
        {
          unMated++;
        }
      else if(gkf.numLinks > 1)
        {
          cifrag->flags.bits.getLinksFromStore = TRUE;  // get links from store
          multiMated++;
        }
      else
        {
          cifrag->flags.bits.getLinksFromStore = FALSE;  // get links from store
          singleMated++;
        }
    }

  // Cycle through all of the fragments in iid order
  for(i = 0; i < GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex); i++)
    {
      InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex,i);
      CIFragT *cifrag;
      GateKeeperLinkRecordIterator GKPLinks;
      GateKeeperLinkRecord GKPLink;
      CIFragT *mcifrag;
      CDS_CID_t j;

      if(!info->set)
        continue;
    
      cifrag = GetCIFragT(ScaffoldGraph->CIFrags,info->fragIndex);
      if(cifrag->type != AS_LBAC || cifrag->linkHead == NULL_LINK)
        continue;
      /*
        if(cifrag->type != AS_EBAC || cifrag->linkHead == NULL_LINK)
        continue;
      */

      CreateGateKeeperLinkRecordIterator(ScaffoldGraph->gkpStore.lnkStore,
                                         cifrag->linkHead, cifrag->iid, &GKPLinks);
      while(NextGateKeeperLinkRecordIterator(&GKPLinks, &GKPLink))
        {
          if(GKPLink.type != AS_MATE &&
             GKPLink.type != AS_BAC_GUIDE)
            continue;

          assert(GKPLink.type == (unsigned int) AS_MATE &&
                 cifrag->type == AS_LBAC);
          /*
            assert(GKPLink.type == (unsigned int) AS_BAC_GUIDE &&
            cifrag->type == AS_EBAC);
          */
      
          assert(cifrag->iid == GKPLink.frag1 || cifrag->iid == GKPLink.frag2);

          if(cifrag->iid == GKPLink.frag2)
            {
              /* Avoid doing things i,j and j,i */
              InfoByIID *info2 =
                GetInfoByIID(ScaffoldGraph->iidToFragIndex,GKPLink.frag1);
              if(!info2->set)
                {
                  fprintf(stderr,"* Fragment with iid " F_CID ", mate of fragment with iid " F_CID " is NOT in assembly\n",
                          GKPLink.frag1, GKPLink.frag2);
                  cifrag->linkHead = NULL_LINK;
                  cifrag->numLinks = 0;
                }
              continue;
            }

          if(GKPLink.frag2 >= GetNumCIFragTs(ScaffoldGraph->CIFrags))
            continue;

          j = GKPLink.frag2;
      
          {
            InfoByIID *infoj = GetInfoByIID(ScaffoldGraph->iidToFragIndex,j);
	  
            if(!infoj->set)
              {
                // 2nd fragment is not part of assembly input
                fprintf(stderr,"* Fragment with iid " F_CID ", mate of fragment with iid " F_CID " is NOT in assembly\n",
                        GKPLink.frag2, GKPLink.frag1);
                continue;
              }

            mcifrag = GetCIFragT(ScaffoldGraph->CIFrags,infoj->fragIndex);
            cifrag->mateOf = infoj->fragIndex;
            mcifrag->mateOf = info->fragIndex;

            // If we saved the one and only link in the fragment
            // record, don't look in the store for it
            if(mcifrag->numLinks == 1 && mcifrag->mateOf != NULLINDEX)
              mcifrag->flags.bits.getLinksFromStore = FALSE;
            if(cifrag->numLinks == 1 && cifrag->mateOf != NULLINDEX)
              cifrag->flags.bits.getLinksFromStore = FALSE;

            mcifrag->dist = cifrag->dist = GKPLink.distance;
            cifrag->linkType = GKPLink.type;
            mcifrag->linkType = GKPLink.type;

            // The CIFragT data structure has been optimized for the case that
            // a fragment has one and only one mate/bac end link.  We need to know
            // whether the mate is a standard innie, or a funky outtie link.  We
            // use a flag bit to distinguish.  The main clients of this bit are
            // in ComputeMatePairStatistics and BuildGraphEdgesFromMultiAlign in GraphCGW_T.c
            //
            cifrag->flags.bits.innieMate = FALSE;
            mcifrag->flags.bits.innieMate = FALSE;
            // Mates can be either INNIE or OUTTIE
            if((GKPLink.type == AS_MATE &&
                GKPLink.orientation == AS_GKP_INNIE) ||
               (GKPLink.type == AS_BAC_GUIDE)){
              mcifrag->flags.bits.innieMate = TRUE;
              cifrag->flags.bits.innieMate = TRUE;
            }
            if(cifrag->flags.bits.hasFalseMate ||
               mcifrag->flags.bits.hasFalseMate ){
              // They must BOTH be false
              if (!(cifrag->flags.bits.hasFalseMate ==
                    mcifrag->flags.bits.hasFalseMate  )){
                fprintf(stderr,"* Frag " F_CID " and Frag " F_CID " have inconsistent mate status\n",
                        cifrag->iid, mcifrag->iid);
                assert(0);
              }
              cifrag->flags.bits.mateStatus = mcifrag->flags.bits.mateStatus = MATE_FALSE;
              cifrag->flags.bits.edgeStatus = UNKNOWN_EDGE_STATUS;
              mcifrag->flags.bits.edgeStatus = UNKNOWN_EDGE_STATUS;
            }else
              cifrag->flags.bits.mateStatus = mcifrag->flags.bits.mateStatus = MATE_OK;
            cifrag->flags.bits.edgeStatus = UNKNOWN_EDGE_STATUS;
            mcifrag->flags.bits.edgeStatus = UNKNOWN_EDGE_STATUS;
          }
        }
    }

  fprintf(GlobalData->stderrc, "Reactivated LBAC mate pairs.\n");
  fprintf(GlobalData->stderrc,
          "\t%d LBAC fragments: %d un-mated, %d singly-mated, %d multiply-mated\n",
          lBACFrags, unMated, singleMated, multiMated);
  /*
    fprintf(GlobalData->stderrc, "Reactivated EBAC mate pairs.\n");
    fprintf(GlobalData->stderrc,
    "\t%d EBAC fragments: %d un-mated, %d singly-mated, %d multiply-mated\n",
    lBACFrags, unMated, singleMated, multiMated);
  */
}

/****************************************************************************/
void ProcessFrags(void)
{
  CDS_CID_t i;
  int32 unmatedFrags = 0;
  GateKeeperFragmentRecord gkf;
  int32 numBACs = getNumGateKeeperLocales(ScaffoldGraph->gkpStore.locStore);
  VA_TYPE(CDS_CID_t) *bactigIDForFBAC = CreateVA_CDS_CID_t(numBACs + 1);

  for(i = 1; i <= numBACs; i++){
    GateKeeperLocaleRecord gkpLoc;
    CDS_CID_t id = NULLINDEX;
    getGateKeeperLocaleStore(ScaffoldGraph->gkpStore.locStore, i, &gkpLoc);
    if(gkpLoc.type == AS_FBAC){
      id = gkpLoc.firstBactig;
    }
    SetVA_CDS_CID_t(bactigIDForFBAC, i, &id);
  }

  // Cycle through all of the fragments in iid order
  for(i = 0; i < GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex); i++){
    InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex,i);
    CIFragT *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
    if(!info->set)
      continue;

    assert(cifrag->iid == i); 
    getGateKeeperFragmentStore(ScaffoldGraph->gkpStore.frgStore, i, &gkf);
    if(gkf.type == AS_FBAC){
      cifrag->locale = *GetCDS_CID_t(bactigIDForFBAC, gkf.localeID);
    }else if(gkf.type == AS_UBAC){
      cifrag->locale = gkf.bactigID;
    }else{
      cifrag->locale = NULLINDEX;
    }
       
    cifrag->linkHead = gkf.linkHead;
    cifrag->numLinks = gkf.numLinks;

#ifdef DEBUG_DETAILED
    fprintf(stderr,"* Frag " F_CID "  iid:" F_CID " numLinks:%d linkHead:" F_CID "\n",
	    info->fragIndex, i, cifrag->numLinks, cifrag->linkHead);
#endif
    
    assert((cifrag->linkHead == NULL_LINK && gkf.numLinks == 0) ||
           (cifrag->linkHead != NULL_LINK && gkf.numLinks > 0));

#ifdef RAT_RUN_1
    if(gkf.type == AS_LBAC)
      {
        fprintf(stderr, "Invalidating mate of AS_LBAC fragment " F_CID "\n", i);
        cifrag->linkHead = NULL_LINK;
        cifrag->numLinks = 0;
      }
    /*
      if(gkf.type == AS_EBAC)
      {
      fprintf(stderr, "Invalidating mate of AS_EBAC fragment " F_CID "\n", i);
      cifrag->linkHead = NULL_LINK;
      cifrag->numLinks = 0;
      }
    */
#endif
    
    if(gkf.numLinks == 0){
      unmatedFrags++;
    }else if(gkf.numLinks > 1){
      cifrag->flags.bits.getLinksFromStore = TRUE;  // get links from store
    }else{
      cifrag->flags.bits.getLinksFromStore = FALSE;  // get links from store
    }
  }
  fprintf(stderr,"* Unmated fragments %d\n", unmatedFrags);
      
  fprintf(stderr,"* Total IIDs: %d\n",
          (int) GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex));
  fflush(NULL);

  // Cycle through all of the fragments in iid order
  for(i = 0; i < GetNumInfoByIIDs(ScaffoldGraph->iidToFragIndex); i++){
    InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex,i);
    CIFragT *cifrag;
    GateKeeperLinkRecordIterator GKPLinks;
    GateKeeperLinkRecord GKPLink;
    CIFragT *mcifrag;
    CDS_CID_t j;


    if(!info->set){
      fprintf(stderr,"* Frag with iid " F_CID " is not set!\n",i);
      continue;
    }
    cifrag = GetCIFragT(ScaffoldGraph->CIFrags,info->fragIndex);
    if(cifrag->linkHead == NULL_LINK){
      //      fprintf(stderr,"* Frag with iid " F_CID " has NULL LINK!\n",i);
      continue;
    }

#ifdef DEBUG_DETAILED
    fprintf(stderr,"* Frag " F_CID "  iid:" F_CID " numLinks:%d linkHead:" F_CID "\n",
	    info->fragIndex, i, cifrag->numLinks, cifrag->linkHead);
#endif

    CreateGateKeeperLinkRecordIterator(ScaffoldGraph->gkpStore.lnkStore, cifrag->linkHead, i, &GKPLinks);
    while(NextGateKeeperLinkRecordIterator(&GKPLinks, &GKPLink)){
      if(GKPLink.type != AS_MATE &&
	 GKPLink.type != AS_BAC_GUIDE)
	continue;

#ifdef DEBUG_DETAILED
      fprintf(stderr,"* Link " F_CID "," F_CID " \n",
	      GKPLink.frag1, GKPLink.frag2);
#endif
      assert((GKPLink.type == (unsigned int)AS_MATE &&
              (cifrag->type == AS_READ ||
               cifrag->type == AS_EXTR ||
               cifrag->type == AS_TRNR ||
               cifrag->type == AS_LBAC)) ||
	     (GKPLink.type == (unsigned int)AS_BAC_GUIDE &&
              cifrag->type == AS_EBAC));
      
      assert(i == GKPLink.frag1 || i == GKPLink.frag2);

      if(i == GKPLink.frag2){ /* Avoid doing things i,j and j,i */
	//		  fprintf(stderr,"* Skipping (" F_CID "," F_CID ")\n",
	//			  GKPLink.frag1, GKPLink.frag2);
	InfoByIID *info2 = GetInfoByIID(ScaffoldGraph->iidToFragIndex,GKPLink.frag1);
	if(!info2->set){
	  // mate is not in assembly...this fragment effectively has no links
	  fprintf(stderr,"* Fragment with iid " F_CID ", mate of fragment with iid " F_CID " is NOT in assembly\n",
		  GKPLink.frag1, GKPLink.frag2);
	  cifrag->linkHead = NULL_LINK;
	  cifrag->numLinks = 0;
	}
	continue;
      }


#define WORRY_ABOUT_DELETED_FRAGS_NOT_INCREMENTAL_RUNS
#ifndef WORRY_ABOUT_DELETED_FRAGS_NOT_INCREMENTAL_RUNS
      if(GKPLink.frag2 > GetNumCIFragTs(ScaffoldGraph->CIFrags)){
	//		  fprintf(stderr,"* Skipping (" F_CID "," F_CID ")\n",
	//			  GKPLink.frag1, GKPLink.frag2);
	/* This can happen in an incremental run */
	continue;
      }
#endif

#if 0
      /* These asserts protect against the case that the input has
	 an ILK message, but the fragments that are referenced do not
	 appear in the input set */
      assert(GetInfoByIID(ScaffoldGraph->iidToFragIndex, GKPLink.frag2)->set);
      assert(GetInfoByIID(ScaffoldGraph->iidToFragIndex, GKPLink.frag1)->set);
      j = GetInfoByIID(ScaffoldGraph->iidToFragIndex, GKPLink.frag2)->fragIndex;
      assert(i == GetInfoByIID(ScaffoldGraph->iidToFragIndex, GKPLink.frag1)->fragIndex);
#endif

      j = GKPLink.frag2;
      {
	InfoByIID *infoj = GetInfoByIID(ScaffoldGraph->iidToFragIndex,j);

        //  BPW has seen one case where infoj was not defined -- it
        //  was the last frag in the store, and that frag was deleted.
        //  He suspected that the ScaffoldGraph was truncated, and
        //  returning NULL instead of an object with 'set' of 0.  He
        //  didn't check if that was true or not, but, hey, if it's
        //  not in the ScaffoldGraph, it's still not in the assembly!
        //
	if((!infoj) || (!infoj->set)){ // 2nd fragment is not part of assembly input
	  fprintf(stderr,"* Fragment with iid " F_CID ", mate of fragment with iid " F_CID " is NOT in assembly\n",
		  GKPLink.frag2, GKPLink.frag1);
	  continue;
	}

	mcifrag = GetCIFragT(ScaffoldGraph->CIFrags,infoj->fragIndex);
#ifdef DEBUG
	fprintf(stderr,"* cifrag->iid " F_CID " i " F_CID " mcifrag->iid " F_CID " j " F_CID " dist:" F_CID "\n",
		cifrag->iid,i,mcifrag->iid,j, GKPLink.distance);
        fprintf(stderr,"* Frag " F_CID " (" F_CID ")'s mate is " F_CID "(" F_CID ")\n",
                i, GKPLink.frag1, j, GKPLink.frag2);
#endif
	cifrag->mateOf = infoj->fragIndex;
	mcifrag->mateOf = info->fragIndex;

	// If we saved the one and only link in the fragment
	// record, don't look in the store for it
	if(mcifrag->numLinks == 1 && mcifrag->mateOf != NULLINDEX)
	  mcifrag->flags.bits.getLinksFromStore = FALSE;
	if(cifrag->numLinks == 1 && cifrag->mateOf != NULLINDEX)
	  cifrag->flags.bits.getLinksFromStore = FALSE;

	mcifrag->dist = cifrag->dist = GKPLink.distance;
	cifrag->linkType = GKPLink.type;
	mcifrag->linkType = GKPLink.type;

	// The CIFragT data structure has been optimized for the case that
	// a fragment has one and only one mate/bac end link.  We need to know
	// whether the mate is a standard innie, or a funky outtie link.  We
	// use a flag bit to distinguish.  The main clients of this bit are
	// in ComputeMatePairStatistics and BuildGraphEdgesFromMultiAlign in GraphCGW_T.c
	//
	cifrag->flags.bits.innieMate = FALSE;
	mcifrag->flags.bits.innieMate = FALSE;
	// Mates can be either INNIE or OUTTIE
	if((GKPLink.type == AS_MATE &&
            GKPLink.orientation == AS_GKP_INNIE) ||
	   (GKPLink.type == AS_BAC_GUIDE)){
	  mcifrag->flags.bits.innieMate = TRUE;
	  cifrag->flags.bits.innieMate = TRUE;
	}
	if(cifrag->flags.bits.hasFalseMate ||
	   mcifrag->flags.bits.hasFalseMate ){
	  // They must BOTH be false
	  if (!(cifrag->flags.bits.hasFalseMate ==
		mcifrag->flags.bits.hasFalseMate  )){
	    fprintf(stderr,"* Frag " F_CID " and Frag " F_CID " have inconsistent mate status\n",
		    cifrag->iid, mcifrag->iid);
	    assert(0);
	  }
	  cifrag->flags.bits.mateStatus = mcifrag->flags.bits.mateStatus = MATE_FALSE;
	  cifrag->flags.bits.edgeStatus = UNKNOWN_EDGE_STATUS;
	  mcifrag->flags.bits.edgeStatus = UNKNOWN_EDGE_STATUS;
	}else
          cifrag->flags.bits.mateStatus = mcifrag->flags.bits.mateStatus = MATE_OK;
	cifrag->flags.bits.edgeStatus = UNKNOWN_EDGE_STATUS;
	mcifrag->flags.bits.edgeStatus = UNKNOWN_EDGE_STATUS;
      }

    }
  }

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

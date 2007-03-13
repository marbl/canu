
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
/* All of the CGW celamy stuff is here */
static char CM_ID[] = "$Id: Celamy_CGW.c,v 1.14 2007-03-13 22:38:49 brianwalenz Exp $";

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
#include "AS_PER_gkpStore.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"

#define ANNOTATED_CELAMY_OUTPUT 1

int do_draw_frags_in_CelamyScaffold =0;
int do_compute_missing_overlaps=0;

// error rate (in parts per thousand) cutoff for counting overlaps
int CelamyOvlCutoff = 15; /* 1.5 percent, default for assembly with error correction */

/* The following is in support of defining a set of Celamy colors to draw with */

#define MAXCOLORS 18
#define COLOR_OFFSET 0
static char *Colors[MAXCOLORS] = {
  "CFF0000",
  "C00FF00",
  "C0000FF",
  "CFFFF00",
  "C00FFFF",
  "CFF00FF",

  "CFF0040",
  "C00FF40",
  "C4000FF",
  "CFF4000",
  "C40FF00",
  "C0040FF",

  "CFF4040",
  "C40FF40",
  "C4040FF",
  "CFF4040",
  "C40FF40",
  "C4040FF"
};


#define  NUM_COLOURS   18
#define  DUNIQUE_COLOUR         1
#define  CONSISTENT_COLOUR     2
#define  ONEFRAG_COLOUR        3
#define  REPEAT_COLOUR         4
#define  BADUNIQUE_COLOUR      5
#define  UNKNOWNUNIQUE_COLOUR      5
#define  CONT_BADUNIQUE_COLOUR 6
#define  ORPHAN_COLOUR         7
#define  LEFTBP_COLOUR         8
#define  RIGHTBP_COLOUR        9
#define  PUNIQUE_COLOUR        10
#define PROCK_COLOUR 11
#define PSTONE_COLOUR 12
#define PWALK_COLOUR 13
#define INVALIDCONTIG_COLOUR 14
#define INVALIDMISPLACEDCONTIG_COLOUR 15
#define CONTIG_COLOUR 16
#define MISPLACEDCONTIG_COLOUR 17

static  char  * Colour_String [NUM_COLOURS]
= {
  "C000000 T2 S  # Unused",
  "CFFFF00 T2 S  # DUnique",
  "CFF8040 T2 S  # Consistent",
  "C808000 T2 S  # OneFrag",
  "CFF0000 T2 S  # Repeat",
  "CFF00FF T2 S  # BadUnique",
  "CFF9A11 T2 S  # ContBadUnique",
  "C00FFFF T2 S  # OrphanFrag",  // Cyan
  "C00FF00 # LeftBP",
  "CFF0000 # RightBP",
  "CFF0077 T2 S  # PUnique",
  "CFF0000 T2 S  # RockCI",
  "C77EF77 T2 S  # StoneCI",
  "C8080FF T2 S  # WalkCI",
  "C00FFFF T2 S  # InvContig",
  "CFF0000 T2 S  # InvMispContig",
  "C0040FF T2 S  # Contig",
  "C8080FF T2 S  # MispContig"
};

#if 0
fprintf(file, "0ContigColor: %s T2 S # Contigs\n",
        Colors[11]);
fprintf(file, "0MisplacedContigColor: %s T2 S # MisplacedContigs\n",
        Colors[10]);
fprintf(file, "0InvalidContigColor: %s T2 S # InvalidContigs\n",
        Colors[4]);
fprintf(file, "0MisplacedInvalidContigColor: %s T2 S # MisplacedInvalidContigs\n",
        Colors[5]);
fprintf(file, "0ContigRealColor: %s T2 S # RealContigs\n",
        Colors[10]);
#endif

/* Forward Declarations */
static void CelamyOrderedScaffolds(FILE *fout,  FILE *fdregs,
                                   CIScaffoldT **scaffoldOrder,
                                   int64 *enda, int64 *endb);
static void OrderScaffoldsForOutput(CIScaffoldT **scaffoldOrder,
                                    int64 *enda, int64 *endb);

#define SCAFFOLD_ROW 1
#define CONTIG_ROW 1
#define DUNIQUE_ROW 3
#define PLACED_ROW 4
#define REALCONTIG_ROW 7
#define UNPLACED_ROW 6


static CIFragT * getFragByIID(ScaffoldGraphT * graph,
                              CDS_CID_t iid)
{
  InfoByIID * info = GetInfoByIID(graph->iidToFragIndex, iid);
  return(GetCIFragT(graph->CIFrags, info->fragIndex));
}


static char * ComputeCIUUCode(ChunkInstanceT *ci){
  switch(ci->flags.bits.cgbType){
    case UU_CGBTYPE:
      return "uu";
      break;
    case UR_CGBTYPE:
      return "ur";
      break;
    case RU_CGBTYPE:
      return "ru";
      break;
    case RR_CGBTYPE:
      return "rr";
      break;
    case XX_CGBTYPE:
      return "xx";
      break;
  }
  return "xx";
}
static int ComputeCIRow(ChunkInstanceT *ci, CIScaffoldT *scaffold){
  if(ci->flags.bits.isScaffold)
    return SCAFFOLD_ROW;
  if(ci->flags.bits.isContig)
    return CONTIG_ROW;
  if(ci->type == DISCRIMINATORUNIQUECHUNK_CGW)
    return DUNIQUE_ROW;
  if(scaffold && scaffold->type == REAL_SCAFFOLD)
    return PLACED_ROW;
  return UNPLACED_ROW;
}

static int ComputeContigColor(ContigT *CI, CIScaffoldT *scaffold){
  
  switch(CI->flags.bits.cgbType){
    case UU_CGBTYPE:
      return(CI->flags.bits.isMisplaced?MISPLACEDCONTIG_COLOUR:CONTIG_COLOUR);
    default:
      return(CI->flags.bits.isMisplaced?INVALIDMISPLACEDCONTIG_COLOUR:INVALIDCONTIG_COLOUR);
  }
}

static int ComputeCIColor(ChunkInstanceT *ci, CIScaffoldT *scaffold){

  int color = REPEAT_COLOUR;
  if(ci->type == DISCRIMINATORUNIQUECHUNK_CGW){
    switch(ci->flags.bits.cgbType){
      case UU_CGBTYPE:
        color = DUNIQUE_COLOUR;
        break;
      case UR_CGBTYPE:
        color = CONT_BADUNIQUE_COLOUR;
        break;
      case RU_CGBTYPE:
      case RR_CGBTYPE:
        color = BADUNIQUE_COLOUR;
        break;
      case XX_CGBTYPE:
        //assert(0);
        color = UNKNOWNUNIQUE_COLOUR;
        break;
    }
  }else{
    if(ci->scaffoldID != NULLINDEX){
      if(!ci->flags.bits.isSurrogate){
        if(ci->flags.bits.isStone){
          color = PSTONE_COLOUR;
        }else
          color = PROCK_COLOUR;
      }else if(ci->flags.bits.isStoneSurrogate){
        color = PSTONE_COLOUR;
      }else{
        color = PWALK_COLOUR;
      }
    }else{
      if(ci->info.CI.numFragments == 1)
        color = ONEFRAG_COLOUR;
      else if((ci->flags.bits.cgbType == UU_CGBTYPE) && (ci->info.CI.coverageStat > CGB_INVALID_CUTOFF))
        color = CONSISTENT_COLOUR;
      else
        color = REPEAT_COLOUR;
    }
  }
  return color;
}



/*******************************************************************************/
/* DumpCelamy Colors */
void DumpCelamyColors(FILE *file){
  { int icolour;
  for(icolour=0; icolour<NUM_COLOURS; icolour++) {
    fprintf(file,"%dCGBColor: %s\n",icolour,Colour_String[icolour]);
  }
  }

  {
    int i;
    for(i = 0; i < MAXCOLORS; i++){
      fprintf(file,"%dCGWColor: %s T2 S # C%d\n",
              i + COLOR_OFFSET, Colors[i],COLOR_OFFSET + i);
    }
    fprintf(file, "0ScaffoldColor: %s T2 S # Scaffolds\n",
            Colors[6]);
    fprintf(file, "0SingleScaffoldColor: %s T2 S # SingleScaffolds\n",
            Colors[9]);
    fprintf(file, "0ScaffoldEdgeColor: %s T2 S # ScaffoldEdges\n",
            Colors[7]);

  }

}


/*******************************************************************************/
/* print header/color info for mate pair indications (clone middle) */
void DumpCelamyMateColors(FILE *file){
  fprintf(file,"1CMColor: C00FF00 T2 S  # Satisfied\n");
  fprintf(file,"2CMColor: CFF0000 T2 S  # Anti\n");
  fprintf(file,"3CMColor: C0000FF T2 S  # TooClose\n");
  fprintf(file,"4CMColor: CFFFF00 T2 S  # TooFar\n");
  fprintf(file,"5CMColor: CFF00FF T2 S  # Normal\n");
  fprintf(file,"6CMColor: C00FFFF T2 S  # Transposed\n");
  fprintf(file,"7CMColor: C88FF88 T2 S  # ExtMateA_B\n");
  fprintf(file,"8CMColor: C880088 T2 S  # ExtMateB_A\n");
  fprintf(file,"9CMColor: CFFFFFF T2 S  # Unmated\n");
  fprintf(file,"0InterScfColor: CFFFFFF T3 S  # InterScaf\n");
  return;
}

/* print header/color info for fragments */
void DumpCelamyFragColors(FILE *file){
  fprintf(file,"0FragColor: C008080 T2 S  # ForwardFrg\n");
  fprintf(file,"1FragColor: C008000 T2 S  # ReverseFrg\n");
  fprintf(file,"2FragColor: C808000 T2 S  # ForwardSurro\n");
  fprintf(file,"3FragColor: C800080 T2 S  # ReverseSurro\n");
  return;
}

typedef enum {
  FWD_FRG_COLOR,
  REV_FRG_COLOR,
  FWD_SURRO_COLOR,
  REV_SURRO_COLOR
} fragColors;


/* CelamyAssembly
   Dumps a simulator-coordinate independent view of the assembly. Currently, scaffolds are drawn
   from left-right by their scaffold index.  Eventually, the ordering code can be refined.
*/
void CelamyAssembly(char *name){
  int numScaffolds = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
  FILE *camOut = NULL;
  FILE *dregsOut = NULL;
  char buffer[256];
  CIScaffoldT **scaffoldOrder = (CIScaffoldT **)safe_calloc((numScaffolds + 1) ,sizeof(CIScaffoldT *));
  int64 *scaffoldPositiona = (int64 *)safe_calloc((numScaffolds + 1) ,sizeof(int64));
  int64 *scaffoldPositionb = (int64 *)safe_calloc((numScaffolds + 1) ,sizeof(int64));
  sprintf(buffer,"%s.asm.cam", name);
  camOut = fopen(buffer,"w");
  sprintf(buffer,"%s.dregs.cam", name);
  dregsOut = fopen(buffer,"w");
  assert(camOut && dregsOut);
  DumpCelamyColors(camOut);
  DumpCelamyColors(dregsOut);
  OrderScaffoldsForOutput(scaffoldOrder, scaffoldPositiona, scaffoldPositionb);
  CelamyOrderedScaffolds(camOut, dregsOut, scaffoldOrder, scaffoldPositiona, scaffoldPositionb);
  safe_free(scaffoldOrder);
  safe_free(scaffoldPositiona);
  safe_free(scaffoldPositionb);
  fclose(camOut);
  fclose(dregsOut);
}

/* OrderScaffoldsForOutput
   Assign each Scaffold to a coordinate range, and intiialize the scaffoldOrder array.
   This is the routine to make more sophisticated to draw the scaffodls in a more
   intelligent order.
*/
void OrderScaffoldsForOutput(CIScaffoldT **scaffoldOrder,
                             int64 *scaffoldPositiona,
                             int64 *scaffoldPositionb){
  CDS_CID_t numScaffolds = 0;
  CIScaffoldT *scaffold;
  int64 aEndCoord = 0, bEndCoord = 0;
  int64 dregsAEndCoord = 0, dregsBEndCoord = 0;
  GraphNodeIterator scaffolds;
  
  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph,
			GRAPH_NODE_DEFAULT);
  while(NULL != (scaffold = NextGraphNodeIterator(&scaffolds))){
    assert(scaffold->bpLength.mean >= 0);
    if(scaffold->type == REAL_SCAFFOLD){
      bEndCoord = aEndCoord + (int64)scaffold->bpLength.mean;
      scaffoldPositiona[scaffold->id] = aEndCoord;
      scaffoldPositionb[scaffold->id] = bEndCoord;
      aEndCoord = bEndCoord + 10;
      assert(aEndCoord >= 0);
      assert(bEndCoord >= 0);
    }else{
      dregsBEndCoord = dregsAEndCoord + (int64)scaffold->bpLength.mean;
      scaffoldPositiona[scaffold->id] = dregsAEndCoord;
      scaffoldPositionb[scaffold->id] = dregsBEndCoord;
      dregsAEndCoord = dregsBEndCoord + 10;
      assert(dregsAEndCoord >= 0);
      assert(dregsBEndCoord >= 0);
    }
#if 0
    fprintf(stderr,"* Drawing scaffold " F_CID " at [" F_COORD "," F_COORD "] bpLength = (%g,%g)\n",
            numScaffolds, scaffold->aEndCoord, scaffold->bEndCoord, scaffold->bpLength.mean, scaffold->bpLength.variance);
#endif
    scaffoldOrder[numScaffolds++] = scaffold;
  }
}

/* CelamyOrderedScaffolds
   Iterate through the scaffolds and generate celamy output
*/
void CelamyOrderedScaffolds(FILE *fout,  FILE *fdregs,
                            CIScaffoldT **scaffoldOrder,
                            int64 *scaffoldPositiona,
                            int64 *scaffoldPositionb){
  int numScaffolds = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
  CDS_CID_t i;


  for(i = 0; i < numScaffolds ;i++){
    CIScaffoldT *scaffold = scaffoldOrder[i];
    if(!scaffold)
      break;
    if(isDeadCIScaffoldT(scaffold))
      continue;
#if 0
    fprintf(stderr,"* Drawing scaffold " F_CID " at [" F_COORD "," F_COORD "]\n",
            i,scaffold->aEndCoord, scaffold->bEndCoord);
#endif
    if(scaffold->type == REAL_SCAFFOLD)
      CelamyScaffold(fout, scaffold, scaffoldPositiona[scaffold->id], scaffoldPositionb[scaffold->id]);
    else
      CelamyScaffold(fdregs, scaffold, scaffoldPositiona[scaffold->id], scaffoldPositionb[scaffold->id]);

  }
}



void safelyAppendOvlInfo(char **ovlsString,OVSoverlap olap, int *lenString, int *lenUsed){
  char teststring[100];
  int testsize;
   
  testsize = snprintf(teststring,99," %d",
		      olap . b_iid);
  assert(testsize <= 100); /* test against truncation */
  assert(testsize >0); /* test against other error */
  if(*lenUsed+testsize>*lenString){
    *lenString+=1000;
    *ovlsString = (char *) safe_realloc((void*)*ovlsString, (*lenString) * sizeof(char));
  }
  strcat(*ovlsString,teststring);
  *lenUsed+=testsize-1; /* -1 because snprintf includes the '\0' in its return,
			   but strcat effectively adds one less since the '\0'
			   at the previous string end is overwritten */
}


void compute_overlaps_off_ends(int id, int *offAEnd, int *offBEnd,char **AEstr, char **BEstr){

  OVSoverlap            olap;
  int retval=0;
  static char *AEndString=NULL,*BEndString=NULL;
  static int lenAstring=0,lenBstring=0;
  int lenAused=1,lenBused=1;

  if(AEndString==NULL){
    lenAstring = 50000;
    AEndString = (char *) safe_malloc(lenAstring*sizeof(char));
  }
  AEndString[0]='\0';

  if(BEndString==NULL){
    lenBstring = 50000;
    BEndString = (char *) safe_malloc(lenBstring*sizeof(char));
  }
  BEndString[0]='\0';

  AS_OVS_setRangeOverlapStore(ScaffoldGraph->frgOvlStore, id, id);

  while  (AS_OVS_readOverlapFromStore(ScaffoldGraph->frgOvlStore, &olap)) {
    //    print_olap(olap);
    if(olap.dat.ovl.corr_erate>CelamyOvlCutoff)continue; /* skip overlaps missing the default conditions for unitigging */
    if  (olap . dat.ovl.a_hang < 0){
      (*offAEnd)++;
      safelyAppendOvlInfo(&AEndString,olap,&lenAstring,&lenAused);
    }
    if  (olap . dat.ovl.b_hang > 0){
      (*offBEnd)++;
      safelyAppendOvlInfo(&BEndString,olap,&lenBstring,&lenBused);
      //      fprintf(stderr,"BEndString x%x\n",BEndString);
    }
  }
  *AEstr = AEndString;
  *BEstr = BEndString;
}



void draw_surroFrags_in_contig_for_CelamyScaffold(FILE *fout, ContigT *ctg, int globallyReversed,int AEndCoord){
  static MultiAlignT *contig=NULL, *unitig=NULL;
  IntUnitigPos *u_list;
  int num_unitigs;
  int i,j;

  if(contig==NULL){
    contig = CreateEmptyMultiAlignT();
    unitig = CreateEmptyMultiAlignT();
  }
  ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,contig,ctg->id,FALSE);

  num_unitigs = GetNumIntUnitigPoss(contig->u_list);
  u_list = GetIntUnitigPos(contig->u_list,0);
  for (i=0;i<num_unitigs;i++) {
    cds_int32 utgID = u_list[i].ident;
    ChunkInstanceT *utg = GetChunkInstanceT(ScaffoldGraph->ChunkInstances,utgID);
    assert(utg!=NULL);
    if(utg->flags.bits.isStoneSurrogate ||
       utg->flags.bits.isWalkSurrogate){
      int utgAEnd,utgBEnd,num_frags,baseUtg;
      IntMultiPos *f_list;
      utg = GetGraphNode(ScaffoldGraph->CIGraph, utg->info.CI.baseID);
      utgAEnd = u_list[i].position.bgn;
      utgBEnd = u_list[i].position.end;
      ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,unitig,utg->id,TRUE);
      num_frags = GetNumIntMultiPoss(unitig->f_list);
      f_list = GetIntMultiPos(unitig->f_list,0);
      for(j=0;j<num_frags;j++){
	char *leftEndOvls,*rightEndOvls; // do not free these -- we don't own them
	int numLeftEndOvls=0,	numRightEndOvls=0;
	int frgAEnd = f_list[j].position.bgn;
	int frgBEnd = f_list[j].position.end;
	int surroColor;


	//	fprintf(stderr,"Surrogate idx %d frag %d [%d,%d]\n",
	//		j,f_list[j].ident,frgAEnd,frgBEnd);

	if(utgAEnd<utgBEnd){
	  frgAEnd += utgAEnd;
	  frgBEnd += utgAEnd;
	} else {
	  frgAEnd = utgAEnd - frgAEnd;
	  frgBEnd = utgAEnd - frgBEnd;
	}	  

	if(globallyReversed){
	  frgAEnd = AEndCoord - frgAEnd;
	  frgBEnd = AEndCoord - frgBEnd;
	} else {
	  frgAEnd += AEndCoord;
	  frgBEnd += AEndCoord;
	}

	if(do_compute_missing_overlaps){
	  compute_overlaps_off_ends(f_list[j].ident,&numLeftEndOvls,&numRightEndOvls,&leftEndOvls,&rightEndOvls);
	}

	if(frgAEnd < frgBEnd){
	  surroColor = FWD_SURRO_COLOR;
	} else {
	  int tmp;
	  char *tmpstr;
	  tmp=frgAEnd;
	  frgAEnd = frgBEnd;
	  frgBEnd=tmp;
	  surroColor = REV_SURRO_COLOR;
	  tmp=numLeftEndOvls;
	  numLeftEndOvls=numRightEndOvls;
	  numRightEndOvls=tmp;
	  tmpstr=leftEndOvls;
	  leftEndOvls = rightEndOvls;
	  rightEndOvls = tmpstr;
	}
	


	if(do_compute_missing_overlaps){

	  //CIFragT *f = getFragByIID(ScaffoldGraph,f_list[j].ident);
	  //	  if(f->dist>=0&&GetDistT(ScaffoldGraph->Dists,f->dist)->mean>15000)

	  fprintf(fout,"%dCtgSurro%d: %d A%dFragColor %d R10 # Contig %d Surrogate Frag %d Overlaps L/R %d/%d details: %s / %s\n",
                  ctg->id, f_list[j].ident,
                  frgAEnd,
                  surroColor,     
                  frgBEnd,
                  ctg->id,
                  f_list[j].ident,
                  numLeftEndOvls,
                  numRightEndOvls,
                  leftEndOvls,
                  rightEndOvls
                  );

	}else{
	  fprintf(fout,"%dCtgSurro%d: %d A%dFragColor %d R10 # Contig %d Surrogate Frag %d\n",
                  ctg->id, f_list[j].ident,
                  frgAEnd,
                  surroColor,     
                  frgBEnd,
                  ctg->id,
                  f_list[j].ident
                  );
	}
      }
    }
  }
}


void draw_frags_in_contig_for_CelamyScaffold(FILE *fout, ContigT *ctg, int globallyReversed,int AEndCoord){
  static MultiAlignT *contig=NULL;
  IntMultiPos *f_list;
  IntMultiPos *frag;
  int num_frags;
  int32 ci_leftcoord, ci_rightcoord;
  int32 t_leftcoord, t_rightcoord;
  int i;
  char buffer[64];
  SeqInterval sim;
  char *coord_start;
  fragColors frgcolor;

  // color is contig
  // output a line for the contig in the contig row with contig color
  // then, loop through the unitigs, outputting them too;

  if(contig==NULL){
    contig = CreateEmptyMultiAlignT();
  }

  ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,contig,ctg->id,FALSE);
  num_frags=GetNumIntMultiPoss(contig->f_list);

  f_list = GetIntMultiPos(contig->f_list,0);
  for (i=0;i<num_frags;i++) {
    uint64 fuid;
    int numLeftEndOvls=0,numRightEndOvls=0;
    char *leftEndOvls, *rightEndOvls; // do not free these -- we don't own them
    frag = &f_list[i];
    t_rightcoord = MAX(frag->position.bgn,frag->position.end);
    t_leftcoord =  MIN(frag->position.bgn,frag->position.end);
    if(globallyReversed){
      ci_leftcoord = AEndCoord - t_rightcoord;
      ci_rightcoord = AEndCoord - t_leftcoord;
    }else{
      ci_leftcoord = AEndCoord + t_leftcoord;
      ci_rightcoord = AEndCoord + t_rightcoord;
    }
#if 0
    if ( frag->source ) {
      coord_start = strchr(frag->source,'[');
    } else {
      coord_start = NULL;
    }
    if ( coord_start != NULL && (sscanf(coord_start,"[%d,%d]",&sim.bgn,&sim.end) == 2) ) {
      sprintf(buffer," [%d,%d]",sim.bgn,sim.end);
    } else {
      sprintf(buffer,"");
    }
    if (show_uids) {
      getFrag(frag_store,frag->ident,rsp,FRAG_S_FIXED);
      getAccID_ReadStruct(rsp, &fuid);
    } else {
      fuid = 0;
    }
#endif
    frgcolor = 
      globallyReversed ? 
      ( frag->position.bgn < frag->position.end ? REV_FRG_COLOR : FWD_FRG_COLOR ) :
      ( frag->position.bgn < frag->position.end ? FWD_FRG_COLOR : REV_FRG_COLOR );


    if(do_compute_missing_overlaps){
      compute_overlaps_off_ends(frag->ident,&numLeftEndOvls,&numRightEndOvls,&leftEndOvls,&rightEndOvls);
    }

    if(frgcolor == REV_FRG_COLOR){
      int tmp;
      char *tmpstr;
      tmp=numLeftEndOvls;
      numLeftEndOvls=numRightEndOvls;
      numRightEndOvls=tmp;
      tmpstr=leftEndOvls;
      leftEndOvls = rightEndOvls;
      rightEndOvls = tmpstr;
    }
       
    if(do_compute_missing_overlaps){
      //CIFragT *f = getFragByIID(ScaffoldGraph,frag->ident);
      //       if(f->dist>=0&&GetDistT(ScaffoldGraph->Dists,f->dist)->mean>15000)
      fprintf(fout,"%dCtgFrag%d: %d A%dFragColor %d R10 # Contig %d Frag %d Overlaps L/R %d/%d details: %s / %s\n",
              ctg->id, frag->ident,
              ci_leftcoord,
              frgcolor,     
              ci_rightcoord,
              ctg->id,
              frag->ident,
              numLeftEndOvls,
              numRightEndOvls,
              leftEndOvls,
              rightEndOvls
              );
    } else {
      fprintf(fout,"%dCtgFrag%d: %d A%dFragColor %d R10 # Contig %d Frag %d\n",
              ctg->id, frag->ident,
              ci_leftcoord,
              frgcolor,     
              ci_rightcoord,
              ctg->id,
              frag->ident
              );
    }
  }
  
  return;
}


/* Celamy Scaffold
   The workhorse routine for drawing a simulator-coordinate independent view of a scaffold.
*/
void CelamyScaffold(FILE *fout, CIScaffoldT *scaffold,
                    int64 scaffoldAEndCoord, int64 scaffoldBEndCoord){
  CIScaffoldTIterator CIs;
  ChunkInstanceT *CI;
  int scaffoldReversed = (scaffoldAEndCoord > scaffoldBEndCoord);
  int64 scaffoldMin = CDS_INT64_MAX;

  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);

  while(NULL != (CI = NextCIScaffoldTIterator(&CIs))){
    scaffoldMin = MIN(scaffoldMin, CI->offsetAEnd.mean);
    scaffoldMin = MIN(scaffoldMin, CI->offsetBEnd.mean);
  }

  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);

  while(NULL != (CI = NextCIScaffoldTIterator(&CIs))){
    ContigTIterator cis;
    ChunkInstanceT *ci;
    int64 CIaCoord, CIbCoord;
    int64 contigMin, contigMax;
    CDS_CID_t contigID = CI->id;
    int contigReversed = (CI->offsetAEnd.mean > CI->offsetBEnd.mean);

    if(CI->offsetAEnd.mean > scaffold->bpLength.mean ||
       CI->offsetBEnd.mean > scaffold->bpLength.mean ){
      fprintf(stderr,"* Contig " F_CID " in scaffold " F_CID " has offsets [%d,%d] outside of scaffold length [0,%d]\n",
              CI->id, scaffold->id,
              (int)CI->offsetAEnd.mean,
              (int)CI->offsetBEnd.mean,
              (int)scaffold->bpLength.mean);
    }

    UpdateContigSimCoordinates(CI);

    if(scaffoldReversed){
      CIaCoord = scaffoldBEndCoord - scaffoldMin + (int64)scaffold->bpLength.mean - (int64)MAX(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
      CIbCoord = scaffoldBEndCoord - scaffoldMin + (int64)scaffold->bpLength.mean - (int64)MIN(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
    }else{
      CIaCoord = scaffoldAEndCoord - scaffoldMin + (int64)MIN(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
      CIbCoord = scaffoldAEndCoord - scaffoldMin + (int64)MAX(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
    }

    contigMin = MIN(CIaCoord, CIbCoord); // seems like CIaCoord <= CIbCoord
    contigMax = MAX(CIaCoord, CIbCoord); // should be invariant?

    if(contigMin < MIN(scaffoldAEndCoord,scaffoldBEndCoord) ||
       contigMax > MAX(scaffoldAEndCoord,scaffoldBEndCoord)){
      fprintf(stderr,"* Contig " F_CID " in scaffold " F_CID " has drawing offsets [" F_S64 "," F_S64 "] outside of scaffold length [" F_S64 "," F_S64 "]\n",
              CI->id, scaffold->id,
              contigMax,contigMin,
              scaffoldAEndCoord, scaffoldBEndCoord);
    }

#if 0
    fprintf(stderr,"* CI " F_CID " CIaCoord = " F_S64 " CIbCoord = " F_S64 " contigMin = " F_S64 "\n",
            contigID,CIaCoord, CIbCoord, contigMin);
#endif
    CIaCoord = MAX(0,CIaCoord);
    CIbCoord = MAX(0,CIbCoord);
    if(scaffold->type == REAL_SCAFFOLD /* && (CIaCoord >= 0 && CIbCoord >= 0)*/ ){
      fprintf(fout,F_CID "ScaCtg" F_CID ": " F_S64 " A%dCGBColor " F_S64 " R%d # Scaffold " F_CID " Ctg " F_CID " %s\n",
              scaffold->id, contigID,
              CIaCoord,
              ComputeContigColor(CI, scaffold),
              CIbCoord,
              CONTIG_ROW,
              scaffold->id, contigID,(CI->flags.bits.isMisplaced?"MISPLACED":""));
    }
    if(contigMin == contigMax){
      fprintf(stderr,"* Warning: Zero Length Contig!!!! contigMin =  " F_S64 " contigMax = " F_S64 "\n", contigMin, contigMax);
      DumpContig(stderr,ScaffoldGraph, CI, FALSE);
    }
    InitContigTIterator(ScaffoldGraph, contigID, TRUE, FALSE, &cis);
    while(NULL != (ci = NextContigTIterator(&cis))){
      int64 ciACoord, ciBCoord;
      CDS_CID_t cid = ci->id;
      int color = ComputeCIColor(ci,scaffold);
#ifdef ANNOTATED_CELAMY_OUTPUT
      MultiAlignT *ci_ma = NULL;
#endif
      IntMultiPos *frag;
      CIFragT *ci_frag;
      int num_frags;
      int fi;
      

      num_frags = 0;
      if(scaffoldReversed ^ contigReversed){
        ciACoord = contigMax - MAX((int64)ci->offsetAEnd.mean, (int64)ci->offsetBEnd.mean);
        ciBCoord = contigMax - MIN((int64)ci->offsetAEnd.mean, (int64)ci->offsetBEnd.mean);
      }else{
        ciACoord = contigMin + MIN((int64)ci->offsetAEnd.mean, (int64)ci->offsetBEnd.mean);
        ciBCoord = contigMin + MAX((int64)ci->offsetAEnd.mean, (int64)ci->offsetBEnd.mean);

      }
      //	  assert(ciACoord >= 0 && ciBCoord >= 0);
      if(ciACoord < 0){
        fprintf(stderr,"* Warning ci " F_CID " has negative ciACoord " F_S64 " ==> 0",
                ci->id, ciACoord);
        ciACoord = 0;
      }
      if(ciBCoord < 0){
        fprintf(stderr,"* Warning ci " F_CID " has negative ciBCoord " F_S64 " ==> 0",
                ci->id, ciBCoord);
        ciBCoord = 0;
      }
         
      if(!ci->flags.bits.isSurrogate){
        fprintf(fout,F_CID "CtgCI" F_CID ": " F_S64 " A%dCGBColor " F_S64 " R%d # Contig " F_CID " CI " F_CID " %s cov:%d ",
                contigID, cid,
                ciACoord,
                color,
                ciBCoord,
                ComputeCIRow(ci, scaffold),
                contigID, cid,ComputeCIUUCode(ci), ci->info.CI.coverageStat);
      }else{
        NodeCGW_T *baseCI = GetGraphNode(ScaffoldGraph->CIGraph, ci->info.CI.baseID);
        fprintf(fout,F_CID "CtgCI" F_CID ": " F_S64 " A%dCGBColor " F_S64 " R%d # Contig " F_CID " CI " F_CID " (BaseCI " F_CID " copies %d) %s baseCov:%d ",
                contigID, cid,
                ciACoord,
                color,
                ciBCoord,
                ComputeCIRow(ci, scaffold),
                contigID, cid, 
                baseCI->id, baseCI->info.CI.numInstances, ComputeCIUUCode(ci), baseCI->info.CI.coverageStat);

      }
#ifdef ANNOTATED_CELAMY_OUTPUT
      if(GlobalData->annotateUnitigs && num_frags > 0 && ci_ma){
        frag = GetIntMultiPos(ci_ma->f_list,0); 
        if(ci_ma){ // if we loaded it, unload it...
          UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, cid, TRUE);
        }
      }
#endif
      fprintf(fout,"\n");
    }
    if(do_draw_frags_in_CelamyScaffold){
      draw_frags_in_contig_for_CelamyScaffold(
                                              fout,
                                              CI,
                                              scaffoldReversed ^ contigReversed,
                                              scaffoldReversed ^ contigReversed ? contigMax : contigMin);
      draw_surroFrags_in_contig_for_CelamyScaffold(
                                                   fout,
                                                   CI,
                                                   scaffoldReversed ^ contigReversed,
                                                   scaffoldReversed ^ contigReversed ? contigMax : contigMin);

	  
    }
  }
  
  if(scaffold->type == REAL_SCAFFOLD){
    if(scaffold->info.Scaffold.numElements > 1){
      fprintf(fout,"LNK: ");
      InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);
      while(NULL != (CI = NextCIScaffoldTIterator(&CIs))){
	fprintf(fout,F_CID "ScaCtg" F_CID " ",
		scaffold->id, CI->id);
      }
      fprintf(fout," A0ScaffoldColor \n");
    }else{
      assert(scaffold->info.Scaffold.numElements == 1);
      fprintf(fout,"LNK: " F_CID "ScaCtg" F_CID " A0SingleScaffoldColor\n",
              scaffold->id, scaffold->info.Scaffold.AEndCI);
    }
  }
}




/*******************************************************************************/
/* CelamyCIScaffolds
 *   Draw a simulator dependent view of the assembly
 */
void CelamyCIScaffolds(char *name, ScaffoldGraphT *graph){
  ChunkInstanceT *Ctga, *Ctgb, *CI;
  ChunkInstanceT *CIa;
  CIScaffoldTIterator CIs;
  CDS_CID_t i;
  int outputCalculatedOffsets = GlobalData->outputCalculatedOffsets;
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;
  FILE *fp;

  if(name){
    char temp[200];
    sprintf(temp,"%s.cam",name);
    fp = fopen(temp,"w");
  }else{
    assert(0);
  }

  DumpCelamyColors(fp);  

  InitGraphNodeIterator(&scaffolds, graph->ScaffoldGraph,
			GRAPH_NODE_DEFAULT);
  while(NULL != (scaffold = NextGraphNodeIterator(&scaffolds))){

    i = scaffold->id;
    
    Ctga = GetGraphNode(graph->RezGraph, scaffold->info.Scaffold.AEndCI);
    Ctgb = GetGraphNode(graph->RezGraph, scaffold->info.Scaffold.BEndCI);

    CIa = GetGraphNode(graph->RezGraph,Ctga->info.Contig.AEndCI);


    if(scaffold->type == REAL_SCAFFOLD ){
      int scaffoldReversed = (scaffold->aEndCoord > scaffold->bEndCoord);
      InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);

#if 0
      fprintf(stderr,"* Scaffold " F_CID " is %s [" F_COORD "," F_COORD "]\n",
	      i,
	      (scaffoldReversed?" REVERSED ":" NORMAL "),
	      scaffold->aEndCoord, scaffold->bEndCoord);
#endif
      while(NULL != (CI = NextCIScaffoldTIterator(&CIs))){
	ContigTIterator cis;
	ChunkInstanceT *ci;
	CDS_COORD_t CIaCoord, CIbCoord;
	CDS_COORD_t contigMin, contigMax;
	CDS_CID_t contigID = CI->id;
	int contigReversed = (CI->offsetAEnd.mean > CI->offsetBEnd.mean);

	UpdateContigSimCoordinates(CI);

	if(outputCalculatedOffsets){
	  if(scaffoldReversed){
	    CIaCoord = scaffold->bEndCoord + (CDS_COORD_t)scaffold->bpLength.mean - (CDS_COORD_t)MAX(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
	    CIbCoord = scaffold->bEndCoord + (CDS_COORD_t)scaffold->bpLength.mean - (CDS_COORD_t)MIN(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
	  }else{
	    CIaCoord = scaffold->aEndCoord + (CDS_COORD_t)MIN(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
	    CIbCoord = scaffold->aEndCoord + (CDS_COORD_t)MAX(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
	  }
	}else{
          CIaCoord = MIN(CI->aEndCoord,CI->bEndCoord);
          CIbCoord = MAX(CI->aEndCoord,CI->bEndCoord);
	}
	contigMin = MIN(CIaCoord, CIbCoord);
	contigMax = MAX(CIaCoord, CIbCoord);

	if(contigMin == contigMax){
	  fprintf(stderr,"* CelamyCIScaffolds... looks like there are no simulator coords...bye\n");
	  return;
	}

#if 0

	fprintf(stderr,"* Contig " F_CID " CIaCoord = " F_COORD " CIbCoord = " F_COORD " contigMin = " F_COORD "\n",
		contigID,CIaCoord, CIbCoord, contigMin);
#endif

        fprintf(fp,F_CID "ScaCtg" F_CID ": " F_COORD " A%dCGBColor " F_COORD " R%d # Scaffold " F_CID " Ctg " F_CID " %s\n",
                scaffold->id, contigID,
                CIaCoord,
                ComputeContigColor(CI, scaffold),
                //		  (CI->flags.bits.cgbType == UU_CGBTYPE?"A0ContigColor":"A0InvalidContigColor"),
                CIbCoord,
                CONTIG_ROW,
                scaffold->id, contigID,(CI->flags.bits.isMisplaced?"MISPLACED":""));

	InitContigTIterator(ScaffoldGraph, contigID, TRUE, FALSE, &cis);
	while(NULL != (ci = NextContigTIterator(&cis))){
	  CDS_COORD_t ciACoord, ciBCoord;
	  CDS_CID_t cid = ci->id;
	  int color = ComputeCIColor(ci,scaffold);

	  if(outputCalculatedOffsets){
            if(scaffoldReversed^contigReversed){
              ciACoord = contigMax - MAX((CDS_COORD_t)ci->offsetAEnd.mean, (CDS_COORD_t)ci->offsetBEnd.mean);
              ciBCoord = contigMax - MIN((CDS_COORD_t)ci->offsetAEnd.mean, (CDS_COORD_t)ci->offsetBEnd.mean);
            }else{
              ciACoord = contigMin + MIN((CDS_COORD_t)ci->offsetAEnd.mean, (CDS_COORD_t)ci->offsetBEnd.mean);
              ciBCoord = contigMin + MAX((CDS_COORD_t)ci->offsetAEnd.mean, (CDS_COORD_t)ci->offsetBEnd.mean);

            }
	  }else{
	    ciACoord = MIN(ci->aEndCoord, ci->bEndCoord);
	    ciBCoord = MAX(ci->aEndCoord, ci->bEndCoord);
	  }

	  if(!ci->flags.bits.isSurrogate){
            fprintf(fp,F_CID "CtgCI" F_CID ": " F_COORD " A%dCGBColor " F_COORD " R%d # Contig " F_CID " CI " F_CID " %s cov:%d\n",
                    contigID, cid,
                    ciACoord,
                    color,
                    ciBCoord,
                    ComputeCIRow(ci, scaffold),
                    contigID, cid, ComputeCIUUCode(ci), ci->info.CI.coverageStat);
	  }else{
	    NodeCGW_T *baseCI = GetGraphNode(ScaffoldGraph->CIGraph, ci->info.CI.baseID);
            fprintf(fp,F_CID "CtgCI" F_CID ": " F_COORD " A%dCGBColor " F_COORD " R%d # Contig " F_CID " CI " F_CID " (BaseCI " F_CID " copies %d) %s\n",
                    contigID, cid,
                    ciACoord,
                    color,
                    ciBCoord,
                    ComputeCIRow(ci, scaffold),
                    contigID, cid, 
                    baseCI->id, baseCI->info.CI.numInstances, ComputeCIUUCode(ci));

	  }
	}
      }



      if(scaffold->info.Scaffold.numElements > 1){
        fprintf(fp,"LNK: ");
        InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);

        while(NULL != (CI = NextCIScaffoldTIterator(&CIs))){
          //	if(CI->flags.bits.cgbType == UU_CGBTYPE)
          fprintf(fp,F_CID "ScaCtg" F_CID " ",
                  scaffold->id, CI->id);
        }
        fprintf(fp," A0ScaffoldColor \n");
      }else{
        assert(scaffold->info.Scaffold.numElements == 1);
	fprintf(fp,"LNK: " F_CID "ScaCtg" F_CID " A0SingleScaffoldColor\n",
		scaffold->id, scaffold->info.Scaffold.AEndCI);

      }
      if(outputCalculatedOffsets){
	InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);

	while(NULL != (CI = NextCIScaffoldTIterator(&CIs))){

          if(CI->aEndCoord > 0 && CI->bEndCoord > 0)
            fprintf(fp,F_CID "Ctg" F_CID ": " F_COORD " A0ContigRealColor " F_COORD " R%d # Scaffold " F_CID " Ctg " F_CID "\n",
                    scaffold->id, CI->id,
                    MIN(CI->aEndCoord, CI->bEndCoord),
                    MAX(CI->aEndCoord, CI->bEndCoord),
                    REALCONTIG_ROW,
                    scaffold->id, CI->id);
	}
      }
    }
  }

  {
    GraphNodeIterator contigs;
    ContigT *CI;
    InitGraphNodeIterator(&contigs, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
    while(NULL != (CI = NextGraphNodeIterator(&contigs))){
      if(CI->scaffoldID != NULLINDEX)
        continue;

      if(CI->type == CONTIG_CGW)
        continue;


      // Skip CIs that have instantiated surrogates
      if(CI->type == UNRESOLVEDCHUNK_CGW &&
         CI->info.CI.numInstances > 0)
        continue;

      if(CI->aEndCoord == CI->bEndCoord)
        continue;

      fprintf(fp,F_CID "CI: " F_COORD " A%dCGBColor " F_COORD " R%d # Unscaffolded CI " F_CID "\n",
              CI->id, MIN(CI->aEndCoord, CI->bEndCoord),
              ComputeCIColor(CI, NULL),
              MAX(CI->aEndCoord, CI->bEndCoord),
              ComputeCIRow(CI,NULL),
              CI->id);
    }
  }

  fflush(fp);

  if(name)
    fclose(fp);
}


void MarkMisplacedContigs(void){

  ChunkInstanceT *Ctga, *Ctgb, *CI;
  ChunkInstanceT *CIa;
  CIScaffoldTIterator CIs;
  CDS_CID_t i;
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;
  
  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph,
			GRAPH_NODE_DEFAULT);
  while(NULL != (scaffold = NextGraphNodeIterator(&scaffolds))){
    i = scaffold->id;

    Ctga = GetGraphNode(ScaffoldGraph->RezGraph, scaffold->info.Scaffold.AEndCI);
    Ctgb = GetGraphNode(ScaffoldGraph->RezGraph, scaffold->info.Scaffold.BEndCI);

    CIa = GetGraphNode(ScaffoldGraph->RezGraph,Ctga->info.Contig.AEndCI);


    if(scaffold->type == REAL_SCAFFOLD ){
      int scaffoldReversed = (scaffold->aEndCoord > scaffold->bEndCoord);
      CDS_COORD_t simAEndOffsetCurr = CDS_COORD_MIN;
      CDS_COORD_t simBEndOffsetCurr = CDS_COORD_MIN;
      CDS_COORD_t simAEndOffsetPrev = CDS_COORD_MIN;
      CDS_COORD_t calcAEndOffsetCurr = CDS_COORD_MIN;
      CDS_COORD_t calcAEndOffsetPrev = CDS_COORD_MIN;
      CDS_COORD_t calcBEndOffsetCurr = CDS_COORD_MIN;
      CDS_COORD_t calcBEndOffsetPrev = CDS_COORD_MIN;

#if 0
      fprintf(stderr,"* Scaffold " F_CID " is %s [" F_COORD "," F_COORD "]\n",
	      i,
	      (scaffoldReversed?" REVERSED ":" NORMAL "),
	      scaffold->aEndCoord, scaffold->bEndCoord);
#endif
      InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);
      while(NULL != (CI = NextCIScaffoldTIterator(&CIs))){
	CDS_COORD_t CIaCoord, CIbCoord;
	int contigReversed = (CI->offsetAEnd.mean > CI->offsetBEnd.mean);

	UpdateContigSimCoordinates(CI);
	CI->flags.bits.isMisplaced = FALSE;

	simAEndOffsetPrev = simAEndOffsetCurr;
	calcAEndOffsetPrev = calcAEndOffsetCurr;
	calcBEndOffsetPrev = calcBEndOffsetCurr;
	CIaCoord = CI->aEndCoord;
	CIbCoord = CI->bEndCoord;

	if(CIaCoord <= 0 &&
	   CIbCoord <= 0)
	  continue;
	   

	if(scaffoldReversed){
	  calcAEndOffsetCurr = (CDS_COORD_t)scaffold->bpLength.mean - (CDS_COORD_t)MAX(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
	  calcBEndOffsetCurr = (CDS_COORD_t)scaffold->bpLength.mean - (CDS_COORD_t)MIN(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
	  simAEndOffsetCurr = scaffold->aEndCoord - MAX(CIaCoord, CIbCoord) ;
	  simBEndOffsetCurr = scaffold->aEndCoord - MIN(CIaCoord, CIbCoord) ;
	}else{
	  calcAEndOffsetCurr =  (CDS_COORD_t)MIN(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
	  calcBEndOffsetCurr = (CDS_COORD_t)MAX(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
	  simAEndOffsetCurr = MIN(CIaCoord, CIbCoord) - scaffold->aEndCoord;
	  simBEndOffsetCurr = MAX(CIaCoord, CIbCoord) - scaffold->aEndCoord;
	}

	if(simAEndOffsetCurr == simBEndOffsetCurr)
	  continue;
	
	if(simAEndOffsetCurr <= simAEndOffsetPrev){
	  CI->flags.bits.isMisplaced = TRUE;
	  fprintf(stderr,"*@@@@@* Contig " F_CID " is misplaced should be [" F_COORD "," F_COORD "] is [" F_COORD "," F_COORD "]\n",
		  CI->id, simAEndOffsetCurr, simBEndOffsetCurr, calcAEndOffsetCurr,calcBEndOffsetCurr);
	  fprintf(stderr,"*#####* simAEndOffsetCurr = " F_COORD " simAEndOffsetPrev = " F_COORD "  contigReversed %d scaffold [" F_COORD "," F_COORD "]\n",
		  simAEndOffsetCurr, simAEndOffsetPrev, contigReversed, scaffold->aEndCoord, scaffold->bEndCoord);
	}
      }
    }
  }
}


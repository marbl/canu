
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

static char *rcsid = "$Id: Celamy_CGW.c,v 1.26 2009-10-05 22:49:42 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"

int do_draw_frags_in_CelamyScaffold  = 0;
int do_compute_missing_overlaps      = 0;

#define NUM_COLOURS                    18
#define DUNIQUE_COLOUR                  1
#define CONSISTENT_COLOUR               2
#define ONEFRAG_COLOUR                  3
#define REPEAT_COLOUR                   4
#define BADUNIQUE_COLOUR                5
#define UNKNOWNUNIQUE_COLOUR            5
#define CONT_BADUNIQUE_COLOUR           6
#define ORPHAN_COLOUR                   7
#define LEFTBP_COLOUR                   8
#define RIGHTBP_COLOUR                  9
#define PUNIQUE_COLOUR                 10
#define PROCK_COLOUR                   11
#define PSTONE_COLOUR                  12
#define PWALK_COLOUR                   13
#define INVALIDCONTIG_COLOUR           14
#define INVALIDMISPLACEDCONTIG_COLOUR  15
#define CONTIG_COLOUR                  16
#define MISPLACEDCONTIG_COLOUR         17

#define SCAFFOLD_ROW     1
#define CONTIG_ROW       1
#define DUNIQUE_ROW      3
#define PLACED_ROW       4
#define REALCONTIG_ROW   7
#define UNPLACED_ROW     6

#define FWD_FRG_COLOR    0
#define REV_FRG_COLOR    1
#define FWD_SURRO_COLOR  2
#define REV_SURRO_COLOR  3

static char *Colors[NUM_COLOURS] = {
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
  "C4040FF" };

static  char  * Colour_String [NUM_COLOURS] = {
  "C000000 T2 S  # Unused",
  "CFFFF00 T2 S  # DUnique",
  "CFF8040 T2 S  # Consistent",
  "C808000 T2 S  # OneFrag",
  "CFF0000 T2 S  # Repeat",
  "CFF00FF T2 S  # BadUnique",
  "CFF9A11 T2 S  # ContBadUnique",
  "C00FFFF T2 S  # OrphanFrag",  // Cyan
  "C00FF00       # LeftBP",
  "CFF0000       # RightBP",
  "CFF0077 T2 S  # PUnique",
  "CFF0000 T2 S  # RockCI",
  "C77EF77 T2 S  # StoneCI",
  "C8080FF T2 S  # WalkCI",
  "C00FFFF T2 S  # InvContig",
  "CFF0000 T2 S  # InvMispContig",
  "C0040FF T2 S  # Contig",
  "C8080FF T2 S  # MispContig"
};



static
int
ComputeCIRow(ChunkInstanceT *ci, CIScaffoldT *scaffold) {
  if (ci->flags.bits.isScaffold)
    return(SCAFFOLD_ROW);

  if (ci->flags.bits.isContig)
    return(CONTIG_ROW);

  if (ci->type == DISCRIMINATORUNIQUECHUNK_CGW)
    return(DUNIQUE_ROW);

  if (scaffold && scaffold->type == REAL_SCAFFOLD)
    return(PLACED_ROW);

  return(UNPLACED_ROW);
}


static
int
ComputeCIColor(ChunkInstanceT *ci, CIScaffoldT *scaffold) {
  if (ci->type == DISCRIMINATORUNIQUECHUNK_CGW)
    return(DUNIQUE_COLOUR);

  if ((ci->scaffoldID != NULLINDEX) && (ci->flags.bits.isSurrogate == 0))
    return((ci->flags.bits.isStone) ? PSTONE_COLOUR : PROCK_COLOUR);

  if (ci->scaffoldID != NULLINDEX)
    return((ci->flags.bits.isStoneSurrogate) ? PSTONE_COLOUR : PWALK_COLOUR);

  if (ScaffoldGraph->tigStore->getNumFrags(ci->id, TRUE) == 1)
    return(ONEFRAG_COLOUR);

  if (ScaffoldGraph->tigStore->getUnitigCoverageStat(ci->id) > CGB_INVALID_CUTOFF)
    return(CONSISTENT_COLOUR);

  return(REPEAT_COLOUR);
}



static
void
findOverlapsOffEnds(int    id,
                    int   *offAEnd, int   *offBEnd,
                    char **AEstr,   char **BEstr) {
  OVSoverlap   olap;

  AS_OVS_setRangeOverlapStore(ScaffoldGraph->frgOvlStore, id, id);

  uint64 nOvl = AS_OVS_numOverlapsInRange(ScaffoldGraph->frgOvlStore);

  char *AEpos = *AEstr = (char *)safe_malloc(sizeof(char) * 11 * nOvl);
  char *BEpos = *BEstr = (char *)safe_malloc(sizeof(char) * 11 * nOvl);

  while  (AS_OVS_readOverlapFromStore(ScaffoldGraph->frgOvlStore, &olap, AS_OVS_TYPE_OVL)) {
    if (AS_OVS_decodeQuality(olap.dat.ovl.corr_erate) > 0.015) //AS_UTG_ERROR_RATE
      //  skip overlaps missing the default conditions for unitigging
      continue;

    if (olap.dat.ovl.a_hang < 0) {
      (*offAEnd)++;
      sprintf(AEpos, " %d", olap.b_iid);
      while (*AEpos)  AEpos++;
    }

    if  (olap.dat.ovl.b_hang > 0) {
      (*offBEnd)++;
      sprintf(BEpos, " %d", olap.b_iid);
      while (*BEpos)  BEpos++;
    }
  }

  *AEpos = 0;
  *BEpos = 0;
}



static
void
drawSurrogateFrags(FILE *fout, ContigT *ctg, int globallyReversed,int AEndCoord) {
  int i,j;

  MultiAlignT  *contig = ScaffoldGraph->tigStore->loadMultiAlign(ctg->id, FALSE);
  IntUnitigPos *u_list = GetIntUnitigPos(contig->u_list,0);

  for (int i=0; i<GetNumIntUnitigPoss(contig->u_list); i++) {
    ChunkInstanceT *utgchk = GetChunkInstanceT(ScaffoldGraph->ChunkInstances, u_list[i].ident);

    if (utgchk->flags.bits.isStoneSurrogate ||
        utgchk->flags.bits.isWalkSurrogate) {

      ChunkInstanceT *utg = GetGraphNode(ScaffoldGraph->CIGraph, utg->info.CI.baseID);

      MultiAlignT *unitig = ScaffoldGraph->tigStore->loadMultiAlign(utg->id, TRUE);
      IntMultiPos *f_list = GetIntMultiPos(unitig->f_list,0);

      for (int j=0; j<GetNumIntMultiPoss(unitig->f_list); j++) {
	int frgAEnd    = f_list[j].position.bgn;
	int frgBEnd    = f_list[j].position.end;

	if (u_list[i].position.bgn < u_list[i].position.end) {
	  frgAEnd += u_list[i].position.bgn;
	  frgBEnd += u_list[i].position.bgn;
	} else {
	  frgAEnd = u_list[i].position.bgn - frgAEnd;
	  frgBEnd = u_list[i].position.bgn - frgBEnd;
	}

	if (globallyReversed) {
	  frgAEnd = AEndCoord - frgAEnd;
	  frgBEnd = AEndCoord - frgBEnd;
	} else {
	  frgAEnd += AEndCoord;
	  frgBEnd += AEndCoord;
	}

	int surroColor = (frgAEnd < frgBEnd) ? FWD_SURRO_COLOR : REV_SURRO_COLOR;

	if (do_compute_missing_overlaps) {
          int   numLeftEndOvls  = 0;
          int   numRightEndOvls = 0;
          char *leftEndOvls     = NULL;
          char *rightEndOvls    = NULL;

          if (frgAEnd < frgBEnd)
            findOverlapsOffEnds(f_list[j].ident, &numLeftEndOvls, &numRightEndOvls, &leftEndOvls, &rightEndOvls);
          else
            findOverlapsOffEnds(f_list[j].ident, &numRightEndOvls, &numLeftEndOvls, &rightEndOvls, &leftEndOvls);

          fprintf(fout,"%dCtgSurro%drand%d: %d A%dFragColor %d R10 # Contig %d Surrogate Frag %d Overlaps L/R %d/%d details: %s / %s\n",
                  ctg->id, f_list[j].ident, lrand48() % 9999,
                  frgAEnd,
                  surroColor,
                  frgBEnd,
                  ctg->id,
                  f_list[j].ident,
                  numLeftEndOvls,
                  numRightEndOvls,
                  leftEndOvls,
                  rightEndOvls);

          safe_free(leftEndOvls);
          safe_free(rightEndOvls);

	} else {
          fprintf(fout,"%dCtgSurro%drand%d: %d A%dFragColor %d R10 # Contig %d Surrogate Frag %d\n",
                  ctg->id, f_list[j].ident, lrand48() % 9999,
                  frgAEnd,
                  surroColor,
                  frgBEnd,
                  ctg->id,
                  f_list[j].ident);
	}
      }
    }
  }
}


static
void
drawFrags(FILE *fout, ContigT *ctg, int globallyReversed,int AEndCoord) {
  MultiAlignT *contig = ScaffoldGraph->tigStore->loadMultiAlign(ctg->id, FALSE);
  IntMultiPos *f_list = GetIntMultiPos(contig->f_list, 0);

  for (int i=0; i<GetNumIntMultiPoss(contig->f_list); i++) {
    int32 t_rightcoord = MAX(f_list[i].position.bgn,f_list[i].position.end);
    int32 t_leftcoord =  MIN(f_list[i].position.bgn,f_list[i].position.end);

    int32 ci_leftcoord  = 0;
    int32 ci_rightcoord = 0;
    int   frgcolor      = 0;

    if (globallyReversed) {
      ci_leftcoord  = AEndCoord - t_rightcoord;
      ci_rightcoord = AEndCoord - t_leftcoord;
      frgcolor      = (f_list[i].position.bgn < f_list[i].position.end) ? REV_FRG_COLOR : FWD_FRG_COLOR;
    } else {
      ci_leftcoord  = AEndCoord + t_leftcoord;
      ci_rightcoord = AEndCoord + t_rightcoord;
      frgcolor      = (f_list[i].position.bgn < f_list[i].position.end) ? FWD_FRG_COLOR : REV_FRG_COLOR;
    }

    if (do_compute_missing_overlaps) {
      int   numLeftEndOvls  = 0;
      int   numRightEndOvls = 0;
      char *leftEndOvls     = NULL;
      char *rightEndOvls    = NULL;

      if (frgcolor == FWD_FRG_COLOR)
        findOverlapsOffEnds(f_list[i].ident, &numLeftEndOvls, &numRightEndOvls, &leftEndOvls, &rightEndOvls);
      else
        findOverlapsOffEnds(f_list[i].ident, &numRightEndOvls, &numLeftEndOvls, &rightEndOvls, &leftEndOvls);

      fprintf(fout,"%dCtgFrag%d: %d A%dFragColor %d R10 # Contig %d Frag %d Overlaps L/R %d/%d details: %s / %s\n",
              ctg->id, f_list[i].ident,
              ci_leftcoord,
              frgcolor,
              ci_rightcoord,
              ctg->id,
              f_list[i].ident,
              numLeftEndOvls,
              numRightEndOvls,
              leftEndOvls,
              rightEndOvls);

      safe_free(leftEndOvls);
      safe_free(rightEndOvls);

    } else {
      fprintf(fout,"%dCtgFrag%d: %d A%dFragColor %d R10 # Contig %d Frag %d\n",
              ctg->id, f_list[i].ident,
              ci_leftcoord,
              frgcolor,
              ci_rightcoord,
              ctg->id,
              f_list[i].ident);
    }
  }
}





//static
void
DumpCelamyColors(FILE *file) {

  for (int i=0; i<NUM_COLOURS; i++)
    fprintf(file,"%dCGBColor: %s\n", i, Colour_String[i]);

  for(int i=0; i<NUM_COLOURS; i++)
    fprintf(file,"%dCGWColor: %s T2 S # C%d\n", i, Colors[i], i);

  fprintf(file, "0ScaffoldColor: %s T2 S # Scaffolds\n", Colors[6]);
  fprintf(file, "0SingleScaffoldColor: %s T2 S # SingleScaffolds\n", Colors[9]);
  fprintf(file, "0ScaffoldEdgeColor: %s T2 S # ScaffoldEdges\n", Colors[7]);
}


//static
void
DumpCelamyMateColors(FILE *file) {
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
}


//static
void
DumpCelamyFragColors(FILE *file) {
  fprintf(file,"0FragColor: C008080 T2 S  # ForwardFrg\n");
  fprintf(file,"1FragColor: C008000 T2 S  # ReverseFrg\n");
  fprintf(file,"2FragColor: C808000 T2 S  # ForwardSurro\n");
  fprintf(file,"3FragColor: C800080 T2 S  # ReverseSurro\n");
}


//  The workhorse routine for drawing a simulator-coordinate independent view of a scaffold.
//
//static
void
CelamyScaffold(FILE        *fout,
               CIScaffoldT *scaffold,
               int64        scaffoldAEndCoord,
               int64        scaffoldBEndCoord) {

  CIScaffoldTIterator CIs;
  ChunkInstanceT *CI;

  int scaffoldReversed = (scaffoldAEndCoord > scaffoldBEndCoord);
  int64 scaffoldMin = INT64_MAX;

  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    scaffoldMin = MIN(scaffoldMin, CI->offsetAEnd.mean);
    scaffoldMin = MIN(scaffoldMin, CI->offsetBEnd.mean);
  }

  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    ContigTIterator cis;
    ChunkInstanceT *ci;
    int64 CIaCoord;
    int64 CIbCoord;
    int64 contigMin;
    int64 contigMax;
    CDS_CID_t contigID = CI->id;
    int contigReversed = (CI->offsetAEnd.mean > CI->offsetBEnd.mean);

    if (CI->offsetAEnd.mean > scaffold->bpLength.mean ||
        CI->offsetBEnd.mean > scaffold->bpLength.mean ) {
      fprintf(stderr,"* Contig "F_CID" in scaffold "F_CID" has offsets [%d,%d] outside of scaffold length [0,%d]\n",
              CI->id, scaffold->id,
              (int)CI->offsetAEnd.mean,
              (int)CI->offsetBEnd.mean,
              (int)scaffold->bpLength.mean);
    }

    if (scaffoldReversed) {
      CIaCoord = scaffoldBEndCoord - scaffoldMin + (int64)scaffold->bpLength.mean - (int64)MAX(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
      CIbCoord = scaffoldBEndCoord - scaffoldMin + (int64)scaffold->bpLength.mean - (int64)MIN(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
    } else {
      CIaCoord = scaffoldAEndCoord - scaffoldMin + (int64)MIN(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
      CIbCoord = scaffoldAEndCoord - scaffoldMin + (int64)MAX(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
    }

    contigMin = MIN(CIaCoord, CIbCoord);
    contigMax = MAX(CIaCoord, CIbCoord);

    if (contigMin < MIN(scaffoldAEndCoord,scaffoldBEndCoord) ||
        contigMax > MAX(scaffoldAEndCoord,scaffoldBEndCoord)) {
      fprintf(stderr,"* Contig "F_CID" in scaffold "F_CID" has drawing offsets ["F_S64","F_S64"] outside of scaffold length ["F_S64","F_S64"]\n",
              CI->id, scaffold->id,
              contigMax,contigMin,
              scaffoldAEndCoord, scaffoldBEndCoord);
    }

    CIaCoord = MAX(0,CIaCoord);
    CIbCoord = MAX(0,CIbCoord);

    if (scaffold->type == REAL_SCAFFOLD) {
      fprintf(fout,F_CID "ScaCtg"F_CID": "F_S64" A%dCGBColor "F_S64" R%d # Scaffold "F_CID" Ctg "F_CID"\n",
              scaffold->id, contigID,
              CIaCoord,
              CONTIG_COLOUR,
              CIbCoord,
              CONTIG_ROW);
    }

    if (contigMin == contigMax) {
      fprintf(stderr,"* Warning: Zero Length Contig!!!! contigMin =  "F_S64" contigMax = "F_S64"\n", contigMin, contigMax);
      DumpContig(stderr,ScaffoldGraph, CI, FALSE);
    }

    InitContigTIterator(ScaffoldGraph, contigID, TRUE, FALSE, &cis);

    while(NULL != (ci = NextContigTIterator(&cis))){
      int64 ciACoord, ciBCoord;
      CDS_CID_t cid = ci->id;
      int color = ComputeCIColor(ci,scaffold);
      IntMultiPos *frag;
      CIFragT *ci_frag;
      int fi;

      if(scaffoldReversed ^ contigReversed){
        ciACoord = contigMax - MAX((int64)ci->offsetAEnd.mean, (int64)ci->offsetBEnd.mean);
        ciBCoord = contigMax - MIN((int64)ci->offsetAEnd.mean, (int64)ci->offsetBEnd.mean);
      }else{
        ciACoord = contigMin + MIN((int64)ci->offsetAEnd.mean, (int64)ci->offsetBEnd.mean);
        ciBCoord = contigMin + MAX((int64)ci->offsetAEnd.mean, (int64)ci->offsetBEnd.mean);

      }

      if (ciACoord < 0) {
        fprintf(stderr,"* Warning ci "F_CID" has negative ciACoord "F_S64" ==> 0", ci->id, ciACoord);
        ciACoord = 0;
      }

      if (ciBCoord < 0) {
        fprintf(stderr,"* Warning ci "F_CID" has negative ciBCoord "F_S64" ==> 0", ci->id, ciBCoord);
        ciBCoord = 0;
      }

      if (!ci->flags.bits.isSurrogate) {
        fprintf(fout, F_CID"CtgCI"F_CID": "F_S64" A%dCGBColor "F_S64" R%d # Contig "F_CID" CI "F_CID" cov:%d\n",
                contigID, cid,
                ciACoord,
                color,
                ciBCoord,
                ComputeCIRow(ci, scaffold),
                contigID, cid, ScaffoldGraph->tigStore->getUnitigCoverageStat(ci->id));
      } else {
        NodeCGW_T *baseCI = GetGraphNode(ScaffoldGraph->CIGraph, ci->info.CI.baseID);
        fprintf(fout, F_CID"CtgCI"F_CID": "F_S64" A%dCGBColor "F_S64" R%d # Contig "F_CID" CI "F_CID" (BaseCI "F_CID" copies %d) baseCov:%d\n",
                contigID, cid,
                ciACoord,
                color,
                ciBCoord,
                ComputeCIRow(ci, scaffold),
                contigID, cid,
                baseCI->id, baseCI->info.CI.numInstances, ScaffoldGraph->tigStore->getUnitigCoverageStat(baseCI->id));

      }
    }

    if (do_draw_frags_in_CelamyScaffold) {
      drawFrags(fout, CI,
                scaffoldReversed ^ contigReversed,
                scaffoldReversed ^ contigReversed ? contigMax : contigMin);

      drawSurrogateFrags(fout, CI,
                         scaffoldReversed ^ contigReversed,
                         scaffoldReversed ^ contigReversed ? contigMax : contigMin);
    }
  }

  if (scaffold->type == REAL_SCAFFOLD) {
    if (scaffold->info.Scaffold.numElements > 1) {
      fprintf(fout,"LNK: ");
      InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);
      while (NULL != (CI = NextCIScaffoldTIterator(&CIs)))
	fprintf(fout, F_CID"ScaCtg"F_CID" ", scaffold->id, CI->id);
      fprintf(fout," A0ScaffoldColor \n");
    } else {
      assert(scaffold->info.Scaffold.numElements == 1);
      fprintf(fout,"LNK: "F_CID"ScaCtg"F_CID" A0SingleScaffoldColor\n",
              scaffold->id, scaffold->info.Scaffold.AEndCI);
    }
  }
}



//  Dumps a simulator-coordinate independent view of the assembly. Currently, scaffolds are drawn
//  from left-right by their scaffold index.  Eventually, the ordering code can be refined.
//
//static
void
CelamyAssembly(char *name){
  char buffer[256];

  int numScaffolds = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);

  CIScaffoldT **scaffoldOrder     = (CIScaffoldT **)safe_calloc((numScaffolds + 1) ,sizeof(CIScaffoldT *));
  int64        *scaffoldPositiona = (int64 *)       safe_calloc((numScaffolds + 1) ,sizeof(int64));
  int64        *scaffoldPositionb = (int64 *)       safe_calloc((numScaffolds + 1) ,sizeof(int64));

  //  Assign each Scaffold to a coordinate range, and intiialize the scaffoldOrder array.  This is
  //  the routine to make more sophisticated to draw the scaffodls in a more intelligent order.

  numScaffolds = 0;

  GraphNodeIterator  scaffolds;
  CIScaffoldT       *scaffold;
  int64              aEndCoord = 0, dregsAEndCoord = 0;
  int64              bEndCoord = 0, dregsBEndCoord = 0;

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);

  while (NULL != (scaffold = NextGraphNodeIterator(&scaffolds))) {
    assert(scaffold->bpLength.mean >= 0);

    if (scaffold->type == REAL_SCAFFOLD) {
      bEndCoord = aEndCoord + (int64)scaffold->bpLength.mean;
      scaffoldPositiona[scaffold->id] = aEndCoord;
      scaffoldPositionb[scaffold->id] = bEndCoord;
      aEndCoord = bEndCoord + 10;
    } else {
      dregsBEndCoord = dregsAEndCoord + (int64)scaffold->bpLength.mean;
      scaffoldPositiona[scaffold->id] = dregsAEndCoord;
      scaffoldPositionb[scaffold->id] = dregsBEndCoord;
      dregsAEndCoord = dregsBEndCoord + 10;
    }

    scaffoldOrder[numScaffolds++] = scaffold;
  }

  //  Open files and write headers

  sprintf(buffer,"%s.asm.cam", name);
  FILE *camOut = fopen(buffer,"w");

  sprintf(buffer,"%s.dregs.cam", name);
  FILE *dregsOut = fopen(buffer,"w");

  DumpCelamyColors(camOut);
  DumpCelamyColors(dregsOut);

  //  Iterate through the scaffolds and generate celamy output

  for(int i=0; i<GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); i++) {
    CIScaffoldT *scaffold = scaffoldOrder[i];

    if (scaffold == NULL)
      break;

    if (isDeadCIScaffoldT(scaffold))
      continue;

    CelamyScaffold((scaffold->type == REAL_SCAFFOLD) ? camOut : dregsOut,
                   scaffold,
                   scaffoldPositiona[scaffold->id],
                   scaffoldPositionb[scaffold->id]);
  }

  safe_free(scaffoldOrder);
  safe_free(scaffoldPositiona);
  safe_free(scaffoldPositionb);

  fclose(camOut);
  fclose(dregsOut);
}

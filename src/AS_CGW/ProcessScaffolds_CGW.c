
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
/* $Id: ProcessScaffolds_CGW.c,v 1.2 2004-09-23 20:25:19 mcschatz Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <assert.h>
#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_ID_store.h"
#include "PrimitiveVA.h"
#include "PrimitiveVA_MSG.h"
#include "MultiAlignStore_CNS.h"

/* The following is in support of defining a set of Celamy colors to draw with */

#define MAXCOLORS 18
#define INTERSCAFFDIST 50
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


#define  DUNIQUE_COLOUR         1
#define  CONSISTENT_COLOUR      2
#define  ONEFRAG_COLOUR         3
#define  REPEAT_COLOUR          4
#define  BADUNIQUE_COLOUR       5
#define  CONT_BADUNIQUE_COLOUR  6
#define  ORPHAN_COLOUR          7
#define  LEFTBP_COLOUR          8
#define  RIGHTBP_COLOUR         9
#define  PUNIQUE_COLOUR        10
#define PROCK_COLOUR           11
#define PSTONE_COLOUR          12
#define PWALK_COLOUR           13
#define  NUM_COLOURS           14


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
	"C8080FF T2 S  # WalkCI"
     };


#define SCAFFOLD_ROW 1
#define CONTIG_ROW 1
#define DUNIQUE_ROW 3
#define PLACED_ROW 3
#define BACTIG_ROW 14
#define FRAG_ROW 16
FILE *fastaFile;
FILE *fastaDregsFile;
FILE *celamyFile;
char fastaFileName[FILENAME_MAX];
char fastaDregsFileName[FILENAME_MAX];
char celamyFileName[FILENAME_MAX];
int fasta;
int fastaDregs;
char *fastaIdent;
int show_stddev;
mode_t mode = S_IRWXU | S_IRWXG | S_IROTH;
VA_TYPE(CDS_IID_t) *scaff_index;
VA_TYPE(char) *scaffold_sequence=NULL;

VA_TYPE(CDS_IID_t) *cids=NULL;
VA_TYPE(int32) *reversedVA=NULL;
VA_TYPE(char) *ctmp=NULL;
VA_TYPE(char) *qtmp=NULL;



FragStoreHandle frag_store,bactig_store;
int show_uids;

void CleanExit(int rc) {
  char command[100+FILENAME_MAX];
  if( fastaFile != NULL ){
    fclose(fastaFile);
    sprintf(command,"rm -f %s",fastaFileName);
    fprintf(stderr,"%s\n",command);
    system(command);
  }
  if( celamyFile != NULL ) {
    fclose(celamyFile);
    sprintf(command,"rm -f %s",celamyFileName);
    fprintf(stderr,"%s\n",command);
    system(command);
  }
  exit(rc);
}

int HandleDir(char *filePathAndName, char *fileName) {
   mode_t mode = S_IRWXU | S_IRWXG | S_IROTH;
   char *suffix;
   char *DirName;
   char *FileName;
   DIR *Dir;
   suffix = strrchr(filePathAndName,(int)'/');
   if ( suffix != NULL ) {
      *suffix = '\0';
      DirName = filePathAndName; 
      if ( DirName != NULL ) {
        Dir = opendir(DirName);
        if ( Dir == NULL ) {
          if(mkdir(DirName,mode)){
            fprintf(stderr,"Failure to create directory %s\n", DirName);
            CleanExit(1);
          }
        }
      }
      *suffix = '/';
      FileName = filePathAndName;
    } else {
      FileName = filePathAndName;
    }
    strcpy(fileName,FileName);
    return 1;
}

/* DumpCelamy Colors */

static void DumpCelamyColors(FILE *file){
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
  fprintf(file, "0ContigColor: %s T2 S # Contigs\n",
	  Colors[11]);
  fprintf(file, "0InvalidContigColor: %s T2 S # InvalidContigs\n",
	  Colors[4]);
  fprintf(file, "0ContigRealColor: %s T2 S # RealContigs\n",
	  Colors[10]);
  fprintf(file, "0ScaffoldColor: %s T2 S # Scaffolds\n",
	  Colors[6]);
  fprintf(file, "0SingleScaffoldColor: %s T2 S # SingleScaffolds\n",
	  Colors[9]);
  fprintf(file, "0ContigLinkColor: %s T1 S # ContigLinks\n",
	  Colors[1]);
  fprintf(file, "0ScaffoldEdgeColor: %s T2 S # ScaffoldEdges\n",
	  Colors[7]);
  fprintf(file, "0FragColor: %s T2 S # Bactigs\n",
	  Colors[4]);
  fprintf(file, "1FragColor: %s T1 S # CeleraRead\n",
	  Colors[1]);
  fprintf(file, "2FragColor: %s T1 S # EBAC\n",
	  Colors[12]);
  fprintf(file, "3FragColor: %s T1 S # OtherGuides\n",
	  Colors[3]);
  fprintf(file, "4FragColor: %s T1 S # ShreddedBacFrag\n",
	  Colors[4]);

  }

}

static int ComputeCIRow(IntUnitigPos *ci){
  if(ci->type == AS_UNIQUE_UNITIG)
    return DUNIQUE_ROW;
  //else
  return PLACED_ROW;
}

static int ComputeFragRow(IntMultiPos *fi){
  if(fi->type == AS_BACTIG)
    return BACTIG_ROW;
  //else
  return FRAG_ROW;
}

static int ComputeCIColor(IntUnitigPos *ci) {
  int color;
  switch (ci->type) {
  case AS_UNIQUE_UNITIG:
     color = DUNIQUE_COLOUR;
     break;
  case AS_ROCK_UNITIG:
     color = PROCK_COLOUR;
     break;
  case AS_STONE_UNITIG:
     color = PSTONE_COLOUR;
     break;
  case AS_PEBBLE_UNITIG:
     color = PWALK_COLOUR;
     break;
  case AS_SINGLE_UNITIG:
     color = ONEFRAG_COLOUR;
     break;
  default:
     fprintf(stderr,"Invalid UnitigType %d\n",ci->type);
     assert(0);
  }
  return color;
}

static int ComputeFragColor(IntMultiPos *fi) {
  int color;
  switch (fi->type) {
  case AS_BACTIG:
     color = 0;
     break;
  case AS_READ:
  case AS_EXTR:
  case AS_TRNR:
     color = 1;
     break;
  case AS_EBAC:
     color = 2;
     break;
  case AS_LBAC:
  case AS_STS :
     color = 3;
     break;
  case AS_UBAC:
  case AS_FBAC:
     color = 4;
    break;
  default:
     fprintf(stderr,"Invalid FragType %d\n",fi->type);
     assert(0);
  }
  return color;
}

MultiAlignStoreT *cstore, *ustore;
CDS_COORD_t global_coord;

int IsForward(SeqInterval s) {
  return (s.bgn < s.end);
}

int getMultiAlignSimCoords(MultiAlignT *unitig,
                           CDS_COORD_t *bgn, CDS_COORD_t *end) {
  IntMultiPos *frag_0 = GetIntMultiPos(unitig->f_list,0); 
  IntMultiPos *frag_e = GetIntMultiPos(unitig->f_list,GetNumIntMultiPoss(unitig->f_list)-1); 
  SeqInterval sim_0, sim_e;
  CDS_COORD_t sim_left,sim_right;
  CDS_COORD_t unitig_len = (CDS_COORD_t) GetNumchars(unitig->consensus) - 1;
  char *coord_start;
  if ( frag_0 == NULL ) return 0;
  while ( frag_0 != frag_e ) {
     if ( min(frag_0->position.bgn,frag_0->position.end) == 0 ) break;
     frag_0++;
  }
  if ( max(frag_0->position.bgn,frag_0->position.end) == unitig_len ) { 
     frag_e = frag_0; // this is likely when first frag is a bactig
  } else {
     frag_e=GetIntMultiPos(unitig->f_list,GetNumIntMultiPoss(unitig->f_list)-1); 
     while ( frag_e != frag_0 ) {
       if ( max(frag_e->position.bgn,frag_e->position.end) == unitig_len ) break;
       frag_e--;
     }
  }
  coord_start = strchr(frag_0->source,'[');
  if ( coord_start == NULL || (sscanf(coord_start,"[" F_COORD "," F_COORD "]",&sim_0.bgn,&sim_0.end) != 2) ) return 0;
  coord_start = strchr(frag_e->source,'[');
  if ( coord_start == NULL || (sscanf(coord_start,"[" F_COORD "," F_COORD "]",&sim_e.bgn,&sim_e.end) != 2) ) return 0;
  sim_left = min(min(sim_0.bgn,sim_0.end),min(sim_e.bgn,sim_e.end));
  sim_right = max(max(sim_0.bgn,sim_0.end),max(sim_e.bgn,sim_e.end));
  if ( IsForward(sim_0) ^ IsForward(frag_0->position) ) {
     *bgn = sim_left;
     *end = sim_right;
  } else {
     *end = sim_left;
     *bgn = sim_right;
  }
  return 1;
}

int CelamyContig(FILE *out, CDS_IID_t scaffid, CDS_IID_t contigid, int reverse)  {
  MultiAlignT *contig=GetMultiAlignInStore(cstore,contigid);
  MultiAlignT *unitig;
  IntUnitigPos *u_list;
  IntMultiPos *f_list;
  IntMultiPos *frag;
  SeqInterval multialign_sim;
  CDS_IID_t num_unitigs = (CDS_IID_t) GetNumIntUnitigPoss(contig->u_list);
  CDS_IID_t num_frags = (CDS_IID_t) GetNumIntMultiPoss(contig->f_list);
  CDS_COORD_t leftcoord, rightcoord;
  CDS_COORD_t ci_leftcoord, ci_rightcoord;
  CDS_COORD_t t_leftcoord, t_rightcoord;
  int i;
  char buffer[64];
  SeqInterval sim;
  char *coord_start;
  ReadStructp rsp = new_ReadStruct();
  // color is contig
  // output a line for the contig in the contig row with contig color
  // then, loop through the unitigs, outputting them too;
  SetVA_CDS_IID_t(scaff_index,contigid, &scaffid);
  leftcoord = global_coord;
  rightcoord = global_coord+GetNumchars(contig->consensus);
  fprintf(out,F_IID "ScaCtg" F_IID ": " F_COORD " %s " F_COORD " R%d # Scaffold " F_IID " Ctg " F_IID "\n",
		  scaffid, contigid,
		  leftcoord,
		  "A0ContigColor",
		  rightcoord,
		  CONTIG_ROW,
		  scaffid, contigid);
  u_list = GetIntUnitigPos(contig->u_list,0);
  for (i=0;i<num_unitigs;i++) {
     unitig = GetMultiAlignInStore(ustore,u_list[i].ident);
     t_rightcoord = max(u_list[i].position.bgn,u_list[i].position.end);
     t_leftcoord =  min(u_list[i].position.bgn,u_list[i].position.end);
     if(reverse){
	ci_leftcoord = rightcoord - t_rightcoord;
	ci_rightcoord = rightcoord - t_leftcoord;
     }else{
        ci_leftcoord = leftcoord + t_leftcoord;
        ci_rightcoord = leftcoord + t_rightcoord;
     }
     if ( getMultiAlignSimCoords(unitig,&multialign_sim.bgn,&multialign_sim.end) ) {
        sprintf(buffer," [" F_COORD "," F_COORD "]",multialign_sim.bgn,multialign_sim.end);
     } else {
        buffer[0] = '\0';
     }
     fprintf(out,F_IID "CtgCI" F_IID ": " F_COORD " A%dCGBColor " F_COORD " R%d # Contig " F_IID " CI " F_IID "%s\n",
		  contigid, unitig->id,
		  ci_leftcoord,
		  ComputeCIColor(&u_list[i]),
		  ci_rightcoord,
		  ComputeCIRow(&u_list[i]),
		  contigid, unitig->id,buffer);
  }
  f_list = GetIntMultiPos(contig->f_list,0);
  for (i=0;i<num_frags;i++) {
     CDS_UID_t fuid;
     frag = &f_list[i];
     t_rightcoord = max(frag->position.bgn,frag->position.end);
     t_leftcoord =  min(frag->position.bgn,frag->position.end);
     if(reverse){
	ci_leftcoord = rightcoord - t_rightcoord;
	ci_rightcoord = rightcoord - t_leftcoord;
     }else{
        ci_leftcoord = leftcoord + t_leftcoord;
        ci_rightcoord = leftcoord + t_rightcoord;
     }
     if ( frag->source ) {
       coord_start = strchr(frag->source,'[');
     } else {
       coord_start = NULL;
     }
     if ( coord_start != NULL && (sscanf(coord_start,"[" F_COORD "," F_COORD "]",&sim.bgn,&sim.end) == 2) ) {
       sprintf(buffer," [" F_COORD "," F_COORD "]",sim.bgn,sim.end);
     } else {
       buffer[0] = '\0';
     }
     if (show_uids) {
      if (frag->type == AS_BACTIG ) {
       getFragStore(bactig_store,frag->ident,FRAG_S_FIXED,rsp);
      } else {
       getFragStore(frag_store,frag->ident,FRAG_S_FIXED,rsp);
      }
      getAccID_ReadStruct(rsp, &fuid);
     } else {
      fuid = 0;
     }
     if ( frag->type == AS_BACTIG ) {
       fprintf(out,F_IID "CtgBFrag" F_IID ": " F_COORD " A%dFragColor " F_COORD " R%d # Contig " F_IID " Frag " F_IID " (" F_UID ",%c)%s\n",
		  contigid, frag->ident,
		  ci_leftcoord,
		  ComputeFragColor(frag),
		  ci_rightcoord,
		  ComputeFragRow(frag),
		  contigid, frag->ident,fuid,frag->type,buffer);
     } else {
       fprintf(out,F_IID "CtgFrag" F_IID ": " F_COORD " A%dFragColor " F_COORD " R%d # Contig " F_IID " Frag " F_IID " (" F_UID ",%c)%s\n",
		  contigid, frag->ident,
		  ci_leftcoord,
		  ComputeFragColor(frag),
		  ci_rightcoord,
		  ComputeFragRow(frag),
		  contigid, frag->ident,fuid,frag->type,buffer);
     }
  }
  global_coord=rightcoord;
  delete_ReadStruct(rsp);
  return 1;
}

static void Complement(char *in, CDS_COORD_t len);

#define HENDERSON_GAP_SIZE (50)
#define ZLAI_NEG_GAP_SIZE (50)
#define ZLAI_SMALL_GAP_SIZE (51)
#define IR_GAP_SIZE (50)
#define ZLAI_GAP_SIZE (51)

CDS_COORD_t compute_gap(double gapsize, int32 alternate_gap_spec){
  switch(alternate_gap_spec){
  default:
  case 0:
    return (CDS_COORD_t) max(IR_GAP_SIZE, gapsize);

  case 1: // Henderson
    if(gapsize > 0){
      return (CDS_COORD_t)gapsize; // truncate and return
    }else{
      return HENDERSON_GAP_SIZE;
    }

  case 2:  // ZLai
    if(gapsize < 0){
       return (ZLAI_GAP_SIZE - 1);
    }else  if(gapsize < ZLAI_GAP_SIZE){
      return (ZLAI_GAP_SIZE );
    } else{
      return (CDS_COORD_t) gapsize;
    }
  }
}





int FastaScaffold(FILE *out, IntScaffoldMesg *scaff, int alternate_gap_spec)  {
  CDS_IID_t i,num_pairs=scaff->num_contig_pairs; 
  CDS_IID_t contigid;
  MultiAlignT *contig;
  IntContigPairs *cp=scaff->contig_pairs;
  int32 *reversed;
  char *fseq;
  char *cseq;
  CDS_COORD_t running_length=0;
  char nchar;
  CDS_COORD_t scaffold_length,contig_length,flen,ngaps;

  switch(alternate_gap_spec){
  case 0:
  default:
    nchar = 'n';
    break;
  case 1:
    nchar = 'N';
    break;
  case 2:
    nchar = 'S';
    break;
  }
  if(cids == NULL){
    cids=CreateVA_CDS_IID_t(1);
    reversedVA=CreateVA_int32(num_pairs + 1);
    ctmp=CreateVA_char(200000);
    qtmp=CreateVA_char(200000);
    scaffold_sequence = CreateVA_char(200000);
  }

    ResetVA_CDS_IID_t(cids);
    ResetVA_int32(reversedVA);
    ResetVA_char(ctmp);
    ResetVA_char(qtmp);
    ResetVA_char(scaffold_sequence);

    reversed = Getint32(reversedVA,0);
  contigid = cp[0].contig1;
  // output label line for this scaffold
  fprintf(stderr,"* alternate_gap_spec = %d\n", alternate_gap_spec);
  if(alternate_gap_spec == 2){
  /*The info line is as follows             
    >173000047796280 U       other       2 173000047796250 173000047796251
    scf UID             orient    class     #contigs    ctg1 UID      ctg2 UID
 
    so you should put U for orient, seed for class.
  */
    fprintf(out,">" F_IID " U seed " F_IID " ", scaff->iaccession, num_pairs + 1);
    fprintf(stderr,">" F_IID " U seed " F_IID " ", scaff->iaccession, num_pairs + 1);
    for(i = 0; i < (num_pairs>0? num_pairs:1); i++){
      fprintf(out,F_IID " ", cp[i].contig1);
      fprintf(stderr,F_IID " ", cp[i].contig1);
    }
    if(num_pairs > 0){
      fprintf(out,F_IID " ", cp[i].contig2);
      fprintf(stderr,F_IID " ", cp[i].contig2);
    }
    fprintf(out,"\n");
    fprintf(stderr,"\n");

  }else   if ( ! show_stddev ) {
    fprintf(out,">%s_Scaffold_" F_IID "\n",fastaIdent,scaff->iaccession);
  } else {
    fprintf(out,">%s_Scaffold_" F_IID " stddev: " F_IID,fastaIdent,scaff->iaccession,num_pairs);
    for ( i=0;i<num_pairs;i++ ) {
      fprintf(out," %.3f",cp[i].stddev);
    }
    fprintf(out,"\n");
  }

  scaffold_length = 0;
  // calculate length of sequence:
  scaffold_length += GetMultiAlignUngappedLength(GetMultiAlignInStore(cstore,contigid));
  for ( i=0;i<num_pairs;i++ ) {
    scaffold_length += compute_gap(cp[i].mean, alternate_gap_spec);
    //     scaffold_length += max(20,cp[i].mean);

     scaffold_length += GetMultiAlignUngappedLength(GetMultiAlignInStore(cstore,cp[i].contig2));
  }
  EnableRangeVA_char(scaffold_sequence,scaffold_length+1);
  EnableRangeVA_int32(reversedVA, num_pairs + 1);
  reversed = Getint32(reversedVA,0);
  // Compute the orientation of each contig in the scaffold
  if(num_pairs == 0){
    reversed[0] = FALSE;
  }else{
    switch(cp[0].orient){
    case AB_AB:
    case AB_BA:
      reversed[0] = FALSE;
      break;
    case BA_AB:
    case BA_BA:
      reversed[0] = TRUE;
      break;
    default:
      assert(0);
      break;
    }
  
    for ( i=0;i<num_pairs;i++ ) {
      switch(cp[i].orient){
      case AB_AB:
      case BA_BA:
	reversed[i+1] = reversed[i];
	break;
      case BA_AB:
      case AB_BA:
	reversed[i+1] = !reversed[i];
	break;
      default:
        assert(0);
        break;
      }
    }
  }


  // append first contig to scaffold_sequence
  contig = GetMultiAlignInStore(cstore,contigid);
  assert(contig!=NULL);
  contig_length = GetMultiAlignUngappedLength(contig);
  GetMultiAlignUngappedConsensus(contig,ctmp,qtmp);
  cseq = Getchar(ctmp,0);
  fseq = Getchar(scaffold_sequence,0);

  if(reversed[0]){
    Complement(cseq,contig_length);
  }
  running_length+=contig_length;
  if (running_length > GetNumchars(scaffold_sequence)) {
    fprintf(stderr,"FastaScaffold warning: unexpectedly long string in scaffold\n");
  }
  memcpy(fseq, Getchar(ctmp,0),contig_length);
  fseq+=contig_length;


  for ( i=0;i<num_pairs;i++ ) {
    ngaps = compute_gap(cp[i].mean, alternate_gap_spec);

     // append ngaps N's (or n's, depending on alternate_gap_spec)  to string as intercontig space
     ResetVA_char(ctmp);  
     EnableRangeVA_char(ctmp,ngaps+1);
     running_length+=ngaps;
     memset(Getchar(ctmp,0),nchar,ngaps);
     memset(fseq, nchar, ngaps);
     fseq+=ngaps;
     contigid = cp[i].contig2;
     contig = GetMultiAlignInStore(cstore,contigid);
     contig_length = GetMultiAlignUngappedLength(contig);
     GetMultiAlignUngappedConsensus(contig,ctmp,qtmp);
     if(reversed[i+1]){
       cseq = Getchar(ctmp,0);
       Complement(cseq,contig_length);
     }
     running_length+=contig_length;
     if (running_length > GetNumchars(scaffold_sequence)) {
       fprintf(stderr,"FastaScaffold warning: unexpectedly long string in scaffold\n");
     }
     memcpy(fseq, Getchar(ctmp,0),contig_length);
     fseq +=contig_length;
  }
  // now, output the scaffold
  fseq = Getchar(scaffold_sequence,0);
  flen = strlen(fseq);
  for (i = 0; i < flen; i += 70) { 
    fprintf(out,"%.*s\n",70,fseq);
    fseq += 70;
  }
  return 1;
}


int FastaDegenerateScaffold(FILE *out, IntDegenerateScaffoldMesg *scaff,VA_TYPE(CDS_IID_t) *is_placed)  {
  CDS_IID_t contigid = scaff->icontig;
  int i;
  MultiAlignT *contig;
  char *fseq;
  CDS_COORD_t running_length=0;
  CDS_COORD_t scaffold_length,contig_length,flen;

  if(cids == NULL){
    cids=CreateVA_CDS_IID_t(1);
    reversedVA=CreateVA_int32(1000);
    ctmp=CreateVA_char(200000);
    qtmp=CreateVA_char(200000);
    scaffold_sequence = CreateVA_char(200000);
  }
  ResetVA_CDS_IID_t(cids);
  ResetVA_int32(reversedVA);
  ResetVA_char(ctmp);
  ResetVA_char(qtmp);
  ResetVA_char(scaffold_sequence);

  // calculate length of sequence:
  scaffold_length = GetMultiAlignUngappedLength(GetMultiAlignInStore(cstore,contigid));

  EnableRangeVA_char(scaffold_sequence,scaffold_length+1);
  // append first contig to scaffold_sequence
  contig = GetMultiAlignInStore(cstore,contigid);
  contig_length = GetMultiAlignUngappedLength(contig);
  if ( GetNumIntUnitigPoss(contig->u_list) == 1 ) {
    CDS_IID_t uid= GetIntUnitigPos(contig->u_list,0)->ident;
    if ( *(GetCDS_IID_t(is_placed,uid)) ) return 0;
  }
  GetMultiAlignUngappedConsensus(contig,ctmp,qtmp);
  fseq = Getchar(scaffold_sequence,0);
  running_length+=contig_length;
  if (running_length > GetNumchars(scaffold_sequence)) {
    fprintf(stderr,"FastaScaffold warning: unexpectedly long string in scaffold\n");
  }
  // output label line for this scaffold
  fprintf(out,">%s_DegenerateScaffold_" F_IID "\n",fastaIdent,contigid);
  memcpy(fseq, Getchar(ctmp,0),contig_length);
  flen = strlen(fseq);
  for (i = 0; i < flen; i += 70) { 
    fprintf(out,"%.*s\n",70,fseq);
    fseq += 70;
  }
  return 1;
}

int CelamyProcessScaffold(FILE *out, IntScaffoldMesg *scaff)  {
  CDS_IID_t i,num_pairs=scaff->num_contig_pairs; 
  CDS_IID_t contigid;
  int reverse;
  IntContigPairs *cp=scaff->contig_pairs;
  contigid = cp[0].contig1;

  if(cids == NULL){
    cids=CreateVA_CDS_IID_t(1);
    reversedVA=CreateVA_int32(1000);
    ctmp=CreateVA_char(200000);
    qtmp=CreateVA_char(200000);
  }
  ResetVA_CDS_IID_t(cids);
  ResetVA_int32(reversedVA);
  ResetVA_char(ctmp);
  ResetVA_char(qtmp);

  if ( contigid == cp[0].contig2 && num_pairs > 0) {
    fprintf(stderr,"Warning: Singleton Scaffold " F_IID " has invalid num_contig_pairs "
                   "fields (should be 0)\n",scaff->iaccession);
    num_pairs = 0;
  }
  // draw first contig (forward oriented)
  CelamyContig(out,scaff->iaccession,contigid,0);
  AppendVA_CDS_IID_t(cids,&contigid);
  if ( num_pairs == 0 ) { // singleton scaffold
    global_coord += INTERSCAFFDIST;
    fprintf(out,"LNK: " F_IID "ScaCtg" F_IID " A0SingleScaffoldColor\n",
		scaff->iaccession, contigid);
  } else {
    for ( i=0;i<num_pairs;i++ ) {
	global_coord += cp[i].mean;
      contigid = cp[i].contig2;
      AppendVA_CDS_IID_t(cids,&contigid);
      reverse = (cp[i].orient == AB_AB || cp[i].orient == BA_AB)?0:1;
      CelamyContig(out,scaff->iaccession,contigid,reverse);
    }
    fprintf(out,"LNK: ");
    for (i=0;i<GetNumCDS_IID_ts(cids);i++) {
	fprintf(out,F_IID "ScaCtg" F_IID " ",
		scaff->iaccession, *GetCDS_IID_t(cids,i));
    }
    fprintf(out," A0ScaffoldColor \n");
  }
  global_coord += 10;
  fflush(out);
  return 1;
}

VA_DEF(IntContigLinkMesg)

int main(int argc, char *argv[])
{ GenericMesg *pmesg;
  MesgReader   reader;
  IntConConMesg *contig;
  IntUnitigMesg *unitig;
  IntContigLinkMesg *link;
  MultiAlignT *ma;
  char *suffix;
  char fastaNameBuffer[FILENAME_MAX];
  char fastaIdentifier[FILENAME_MAX];
  char celamyNameBuffer[FILENAME_MAX];
  VA_TYPE(char) *dummy_consensus;
  VA_TYPE(CDS_IID_t) *link_index;
  VA_TYPE(IntContigLinkMesg) *clinks;
  VA_TYPE(CDS_IID_t) *is_placed;
  uint32 placed=1;
  uint32 unplaced=0;
  CDS_IID_t linkid;
  CDS_IID_t num_links;
  int ch;
  int do_all = 0;
  int celamy = 0;
  int alternate_gap_spec = 0;
  show_stddev =0;
  show_uids = 0;
  do_all = 1;
  cstore = CreateMultiAlignStoreT(0);
  ustore = CreateMultiAlignStoreT(0);
  link_index = CreateVA_CDS_IID_t(0);
  scaff_index = CreateVA_CDS_IID_t(0);
  is_placed = CreateVA_CDS_IID_t(0);
  clinks = CreateVA_IntContigLinkMesg(0);
  global_coord = 0;
  optarg = NULL;
  fasta = 0;
  fastaDregs = 0;
  if ( argc < 2 ) {
     fprintf(stderr,"Try -h for usage\n");
     exit(1);
  }
  while ( ((ch = getopt(argc, argv, "h?df:c:sF:B:HZ")) != EOF)) {
        switch(ch) {
	case 'd':
	  fastaDregs = 1;
	  break;
        case 'f':
          fasta = 1;
          strcpy(fastaNameBuffer, optarg);
          HandleDir(fastaNameBuffer,fastaFileName);
          fastaFile = fopen(fastaFileName,"w");
          if (fastaFile == NULL ) {
            fprintf(stderr,"Failure to create fasta file %s\n", fastaFileName);
            CleanExit(1);
          }
	  if(fastaDregs){
	    fastaDregsFile = fopen(fastaDregsFileName,"w");
	    if (fastaDregsFile == NULL ) {
	      fprintf(stderr,"Failure to create fasta file %s\n", fastaDregsFileName);
	      CleanExit(1);
	    }
	  }
          strcpy(fastaIdentifier,fastaFileName);
          fastaIdent = strrchr(fastaIdentifier,'/');
          if ( fastaIdent == NULL ) { 
            fastaIdent = fastaIdentifier;
          } else {
            fastaIdent++;
          }
          suffix = strrchr(fastaIdentifier,(int)'.'); 
          if(suffix!=NULL) *suffix = '\0'; // this cuts off the ext, so filename root can be used
          while ( (suffix = strchr(fastaIdent,'.')) != NULL ) {
             *suffix = '_';
          }
          break;
        case 'c':
          celamy = 1;
          strcpy(celamyNameBuffer, optarg);
          HandleDir(celamyNameBuffer,celamyFileName);
          celamyFile = fopen(celamyFileName,"w");
          if (celamyFile == NULL ) {
            fprintf(stderr,"Failure to create celamy file %s\n", celamyFileName);
            CleanExit(1);
          }
          DumpCelamyColors(celamyFile);
          break;
        case 'F':  // specify FragStore for use to get fragment UIDs in Celamy file
          show_uids = 1;
          frag_store = openFragStore(optarg, "rb");
          break;
        case 'B':  // specify BactigStore for use to get fragment UIDs in Celamy file
          bactig_store = openFragStore(optarg, "rb");
          break;
       case 'H': // alternate gap spec as per Henderson
	  alternate_gap_spec = 1;
	  break;
       case 'Z': // alternate gap spec as per ZLai
	  alternate_gap_spec = 2;
	  fprintf(stderr,"* Z option invoked, alternate_gap_spec == 2\n");
	  break;
        case 's':
          show_stddev=1; 
          break;
        case 'h':
        case '?':
          {
            fprintf(stderr,"\n\nUsage: process_scaffolds [-H] [-d] [-h] [-c celamy.output.filename [-F Frag.Store -B Bactig.Store]] [-f fasta.output.filename [-s]]\n");
            fprintf(stderr,"\n The -f option produces a multi-fasta file, with one fasta record for each scaffold,\n");
            fprintf(stderr," and 'n's used as intercontig gap placeholders.  For positive gaps, the number of 'n's\n");
            fprintf(stderr," used represents the calculated mean gap size.  For negative gaps, 20 'n's are used.\n");
            fprintf(stderr," The -s suboption to -f (fasta output) directs the program to include gap standard\n");
            fprintf(stderr," deviations to the comment line of each scaffold.  The convention for this output is\n");
            fprintf(stderr,"\n   >UBAC_IID_UID_Scaffold_SID stddev: num_gaps fp[0] ... fp[num_gaps-1]\n");
            fprintf(stderr,"\n that is, num_gaps specifies how many floats to expect in remainder of line.\n");
            fprintf(stderr,"\n The -c option produces a celamy file for the assembly, down to the fragment level\n\n");
            fprintf(stderr,"\n        -F and -B can be used to open frag stores to capture fragment uids for comment field\n\n");
            fprintf(stderr,"\n The -d option produces a multi-fasta dregs file, with one fasta record for each degenerate scaffold,\n");
            fprintf(stderr,"\n The -d option is only meaningful when specified with -f\n");
            fprintf(stderr,"\n The -H option, handles the gaps as per Scott Henderson's spec \n");
            fprintf(stderr,"\n The -Z option, handles the headers and gaps as per Z Lai's  spec \n");
            CleanExit(1);
          }
        default:
          {
            fprintf(stderr,"Invalid option -%c, try -h for usage\n",ch);
            CleanExit(1);
          }
        }
  }

  if(fasta&&fastaDregs){
    sprintf(fastaDregsFileName,"%s.dregs",fastaFileName);
    fastaDregsFile = fopen(fastaDregsFileName,"w");
    if (fastaDregsFile == NULL ) {
      fprintf(stderr,"Failure to create fasta file %s\n", fastaDregsFileName);
      CleanExit(1);
    }
  }

  reader = InputFileType_AS( stdin );
#if 0
  sublist_file = argv[1];
  if ( sublist_file[0] == 'A' ) { do_all = 1;}
   
  if ( !do_all ) {
      char   string[1000];
      CDS_UID_t uid;
      CDS_UID_t  num_uids;
      sublist = fopen(sublist_file,"r");
      if( sublist == NULL )
      {
        fprintf( stderr, "Failed to open list file %s for reading.\n", argv[2] );
        CleanExit(1);
      }
      num_uids = 0;
      while( fgets( string, 1000, sublist ) )
      {
        num_uids++;
      }
      rewind( sublist );
      tig_iids = AllocateID_Array( num_uids );
      tig_iids_found = AllocateID_Array( num_uids );
      if( tig_iids == NULL || tig_iids_found == NULL ) return 1;
      for( this_id = 0; this_id < num_uids - 1; this_id++ )
      {
        fgets( string, 1000, sublist );
        sscanf(string, F_UID, &uid);
        AppendToID_Array( tig_iids, uid, 0 );
      }
      fgets( string, 1000, sublist );
      sscanf(string, F_UID, &uid);
      AppendToID_Array( tig_iids, uid, 1 );
  
      fclose( sublist );
  }
#endif
 dummy_consensus = CreateVA_char(200000);

 while (reader(stdin,&pmesg) != EOF){
    if (pmesg->t ==MESG_IUM)  {
      unitig = pmesg->m;
      if ( strlen(unitig->consensus) != unitig->length) {
         char *cptr;
#if 0
         if ( fasta ) {
           fprintf(stderr,"Input appears to be pre-consensus. Fasta information not available.\n");
           fprintf(stderr,"For fasta file, run consensus first.\n");
           CleanExit(1);
         }
#endif
         ResetVA_char(dummy_consensus);
         EnableRangeVA_char(dummy_consensus,unitig->length+1);
         cptr = Getchar(dummy_consensus,0); 
         memset(cptr,'N',unitig->length);
         unitig->consensus = cptr;
         unitig->quality = cptr;
      }
      SetVA_CDS_IID_t(is_placed, unitig->iaccession,  
                      ((unitig->status == AS_SEP) ? &placed : &unplaced));
      if(celamy){
	ma = CreateMultiAlignTFromIUM(unitig, -1,  0);
	SetMultiAlignInStore(ustore,ma->id,ma);
      }
    }
    if (pmesg->t ==MESG_ICM)  {
      contig = pmesg->m;
      if ( strlen(contig->consensus) != contig->length) {
         char *cptr;
#if 0
         if ( fasta ) {
           fprintf(stderr,"Input appears to be pre-consensus. Fasta information not available.\n");
           fprintf(stderr,"For fasta file, run consensus first.\n");
           CleanExit(1);
         }
#endif
         ResetVA_char(dummy_consensus);
         EnableRangeVA_char(dummy_consensus,contig->length+1);
         cptr = Getchar(dummy_consensus,0); 
         memset(cptr,'N',contig->length);
         contig->consensus = cptr;
         contig->quality = cptr;
      }
      ma = CreateMultiAlignTFromICM(contig, -1,  0);
      SetMultiAlignInStore(cstore,ma->id,ma);
    }
    if (pmesg->t ==MESG_ICL)  {
      link = pmesg->m;
      linkid = GetNumCDS_IID_ts(link_index);
      SetVA_CDS_IID_t(link_index, link->contig1, &linkid);
      AppendVA_IntContigLinkMesg(clinks,link);
    }
    if (pmesg->t ==MESG_ISF)  {
      if ( celamy ) {
        CelamyProcessScaffold(celamyFile, (IntScaffoldMesg *)pmesg->m);
      }
      if ( fasta ) {
        FastaScaffold(fastaFile, (IntScaffoldMesg *)pmesg->m, alternate_gap_spec);
        fflush(fastaFile);
      }
    }
    if (pmesg->t ==MESG_IDS)  {
      if ( fasta  && fastaDregs) {
        FastaDegenerateScaffold(fastaDregsFile, (IntDegenerateScaffoldMesg *)pmesg->m,is_placed);
        fflush(fastaDregsFile);
      }
    }
 }
 // Now, draw in the contig links:
 num_links = GetNumIntContigLinkMesgs(clinks);
#if 0
 for (i=0;i<num_links;i++ ) {
    link = GetIntContigLinkMesg(clinks,i);
    fprintf(stdout,"LNK: " F_IID "ScaCtg" F_IID " " F_IID "ScaCtg" F_IID " A0ContigLinkColor # %c %c (%.3f,%.3f)\n",
		*GetCDS_IID_t(scaff_index,link->contig1),link->contig1,
		*GetCDS_IID_t(scaff_index,link->contig2),link->contig2,
                link->orientation,link->overlap_type,link->mean_distance,link->std_deviation);
 } 
#endif
 if (fasta) {
   if(fastaDregs)
     fclose(fastaDregsFile);

   fclose(fastaFile);
 }
 exit (0);
}


/*** UTILITY ROUTINES ***/

/* Complement the sequence in fragment message a.  This include also
   revsersing the order of the quality values.  The operation does the
   complementation/reversal in place.  Calling it a second time on a
   given fragment restores it to its original state.                */

// Stolen from AS_ALN
static void Complement(char *in, CDS_COORD_t len)
{ static char WCinvert[256];
  static int Firstime = 1;

  if (Firstime)          /* Setup complementation array */
    { 
      int i;
      Firstime = 0;
      for(i = 0; i < 256;i++){
	WCinvert[i] = '?';
      }
      WCinvert['a'] = 't';
      WCinvert['c'] = 'g';
      WCinvert['g'] = 'c';
      WCinvert['t'] = 'a';
      WCinvert['n'] = 'n';
      WCinvert['A'] = 'T';
      WCinvert['C'] = 'G';
      WCinvert['G'] = 'C';
      WCinvert['T'] = 'A';
      WCinvert['N'] = 'N';
      WCinvert['-'] = '-'; // added this to enable alignment of gapped consensi
    }
      
  { /* Complement and reverse sequence */

    { register char *s, *t;
      int c;

      s = in;
      t = in + (len-1);
      while (s < t)
        { // Sanity Check!
	  assert(WCinvert[(int) *t] != '?' &&
		 WCinvert[(int) *s] != '?');

	  c = *s;
          *s++ = WCinvert[(int) *t];
          *t-- = WCinvert[c];
        }
      if (s == t)
        *s = WCinvert[(int) *s];
    }

  }
}

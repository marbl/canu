
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
/* $Id: CreateCelamy.c,v 1.2 2004-09-23 20:25:24 mcschatz Exp $ */

//  This program is intended to create a celamy .cam file for the 
//  object given by UID as the single command line argument.
//  The output is based on a .asm file, which has been pre-indexed.
//  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "AS_global.h"
#include "AS_UTL_HashCommon.h"
#include "AS_UTL_PHash.h"
#include "AS_UTL_Var.h"
#include <assert.h>

#define SEQ_NORMAL 0
#define SEQ_REVERSE 1
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
#define GAP_COLOUR             14
#define  NUM_COLOURS           15


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
	"C000000 T2 S  # IntraScaffoldGap"
     };


#define SCAFFOLD_ROW 1
#define CONTIG_ROW 1
#define GAP_ROW 2
#define DUNIQUE_ROW 3
#define PLACED_ROW 3
#define BACTIG_ROW 14
#define FRAG_ROW 16

/* DumpCelamy Colors */
static void DumpCelamyColors(FILE *file){
   { int icolour;
    for(icolour=0; icolour<NUM_COLOURS; icolour++) {
      fprintf(file,"%dCGBColor: %s\n",icolour,Colour_String[icolour]);
    }
    }

  {
  int32 cosIndex;
  int32 celamiIdNum = 0;
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

static int ComputeCIRow(UnitigPos *ci){
  if(ci->type == AS_UNIQUE_UNITIG)
    return DUNIQUE_ROW;
  //else
  return PLACED_ROW;
}

static int ComputeFragRow(SnapMultiPos *fi){
  if(fi->type == AS_BACTIG)
    return BACTIG_ROW;
  //else
  return FRAG_ROW;
}

static int ComputeCIColor(UnitigPos *ci) {
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

static int ComputeFragColor(SnapMultiPos *fi) {
  int color;
  switch (fi->type) {
  case AS_BACTIG:
     color = 0;
     break;
  case AS_READ:
  case AS_B_READ:
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

int CelamyGap(FILE *cam, CDS_UID_t scaffold_uid, int32 gap_no,SnapContigPairs *pair, int32 *global_coord) {
     int start_coord=*global_coord;
     int gap_len=pair->mean;
     int32 min_coord, max_coord;
     *global_coord = *global_coord+=pair->mean;
     if ( gap_len < 0 )  {
         min_coord = *global_coord;
         max_coord = start_coord;
     } else {
         max_coord = *global_coord;
         min_coord = start_coord;
     }
     fprintf(cam,F_UID "ScfGap" F_S32 ": " F_COORD " A%dCGBColor " F_COORD " R%d # Gap " F_S32 " mean: %f stddev: %f\n",
        scaffold_uid,gap_no, min_coord, GAP_COLOUR, max_coord, GAP_ROW, gap_no, pair->mean,pair->stddev);
     return 1;
}

int TestHash(PHashTable_AS *hash,CDS_UID_t uid){ 
     int lookup_rc,hash_rc; 
     PHashValue_AS value;
     lookup_rc = LookupInPHashTable_AS (hash, 1, uid, &value);
     if ( lookup_rc != HASH_SUCCESS) {
       return -1;
     } else {
       return value.IID;
     }
}

int MarkHash(PHashTable_AS *hash,CDS_UID_t uid){ 
     int lookup_rc,hash_rc; 
     PHashValue_AS value;
     lookup_rc = LookupInPHashTable_AS (hash, 1, 
	   uid, 
	   &value);
     value.IID+=1;
     DeleteFromPHashTable_AS(hash, 1, uid);
     hash_rc = InsertInPHashTable_AS(&hash,1, (CDS_UID_t) uid, &value, FALSE,FALSE);
     return value.IID;
}

int CelamyContig(FILE *cam, SnapConConMesg *contig, int32 *global_coord, int reverse, PHashTable_AS *occurences)  {
  UnitigPos *u_list;
  UnitigPos *unitig;
  SnapMultiPos *f_list;
  SnapMultiPos *frag;
  SeqInterval multialign_sim;
  int num_unitigs=contig->num_unitigs;
  int num_frags=contig->num_pieces;
  int32 leftcoord, rightcoord;
  int32 ci_leftcoord, ci_rightcoord;
  int32 t_leftcoord, t_rightcoord;
  int i;
  char buffer[64];
  SeqInterval sim;
  char *coord_start;
  PHashValue_AS value;
  int copy_number;
  // color is contig
  // output a line for the contig in the contig row with contig color
  // then, loop through the unitigs, outputting them too;
  leftcoord = *global_coord;
  rightcoord = leftcoord+contig->length;
  MarkHash(occurences,contig->eaccession); // mark this contig as having been seen and output to the cam file
  fprintf(cam,F_UID "Ctg: " F_COORD " %s " F_COORD " R%d # Ctg " F_UID " (" F_IID ")\n",
		  contig->eaccession,
		  leftcoord,
		  "A0ContigColor",
		  rightcoord,
		  CONTIG_ROW,
		  contig->eaccession,contig->iaccession);
  u_list = contig->unitigs;
  for (i=0;i<num_unitigs;i++) {
     unitig = &u_list[i];
     copy_number = MarkHash(occurences,unitig->eident);
     t_rightcoord = max(unitig->position.bgn,unitig->position.end);
     t_leftcoord =  min(unitig->position.bgn,unitig->position.end);
     if(reverse){
	ci_leftcoord = rightcoord - t_rightcoord;
	ci_rightcoord = rightcoord - t_leftcoord;
     }else{
        ci_leftcoord = leftcoord + t_leftcoord;
        ci_rightcoord = leftcoord + t_rightcoord;
     }
     fprintf(cam,F_UID "CtgCI%d: " F_COORD " A%dCGBColor " F_COORD " R%d # Contig " F_UID " CI " F_UID "\n",
		  unitig->eident,copy_number,
		  ci_leftcoord,
		  ComputeCIColor(unitig),
		  ci_rightcoord,
		  ComputeCIRow(unitig),
		  contig->eaccession, unitig->eident);
  }
  f_list = contig->pieces;
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
        sprintf(buffer,"");
     }
     fuid = frag->eident;
     if ( frag->type == AS_BACTIG ) {
       fprintf(cam,F_UID "CtgBFrag: " F_COORD " A%dFragColor " F_COORD " R%d # Contig " F_UID " Frag " F_UID " (%c)%s\n",
		  fuid,
		  ci_leftcoord,
		  ComputeFragColor(frag),
		  ci_rightcoord,
		  ComputeFragRow(frag),
		  contig->eaccession, fuid,frag->type,buffer);
     } else {
       fprintf(cam,F_UID "CtgFrag: " F_COORD " A%dFragColor " F_COORD " R%d # Contig " F_UID " Frag " F_UID " (%c)%s\n",
		  fuid,
		  ci_leftcoord,
		  ComputeFragColor(frag),
		  ci_rightcoord,
		  ComputeFragRow(frag),
		  contig->eaccession, fuid,frag->type,buffer);
     }
  }
  *global_coord=rightcoord;
  return 1;
}


GenericMesg *grab_message(FILE *mesg_file, long offset) {
  GenericMesg *pmesg;
  CDS_FSEEK(mesg_file,(off_t) offset,SEEK_SET);
  ReadProtoMesg_AS(mesg_file,&pmesg);
  return pmesg; 
}

typedef struct LinkAddress {
  CDS_UID_t uid;
  long offset;
} LinkAddress;

VA_DEF(LinkAddress);

int PrintLinks(FILE *cam,CDS_UID_t uid,PHashTable_AS *occurences,VA_TYPE(LinkAddress) *links, PHashTable_AS *linkptr, FILE *asm) {
      // lookup clks, and foreach, check whether paired contig is present
      //  if present, check whether its links have already been printed...
      //  if not, print the CLK
   GenericMesg *pmesg;
   LinkAddress *linkaddress;
   SnapContigLinkMesg *CLK;
   int hasclk=TestHash(linkptr, uid);
   int other_count;
   if ( hasclk < 0 ) return 0; // no clks recorded for this contig
   // hasclk now should point to spot in links with this contigs CLK offsets
   linkaddress = GetLinkAddress(links,hasclk++);
   while ( linkaddress->uid == uid ) {
      CDS_UID_t other;
      // Read protomesg from file
      pmesg = grab_message(asm,linkaddress->offset);
      CLK = (SnapContigLinkMesg *) pmesg->m;
      other = ( CLK->econtig1 == uid )? CLK->econtig2: CLK->econtig1;
      other_count = TestHash(occurences,other);
      switch(other_count) {
      case 0: // this contig isn't in the scaffold, don't refer to it in a link
             // later, it would be nice to add a non-scaffold contig for it showing where the link goes off to.
         break;
      case 1:
          // print the link
         {int i; 
         fprintf(cam,"LNK: " F_UID "Ctg " F_UID "Ctg A0ContigLinkColor # " F_UID " " F_UID " mea: %f stddev: %f num_contributing: " F_S32 " (",
          uid,other,uid,other,CLK->mean_distance,CLK->std_deviation,CLK->num_contributing);
         for (i=0;i<CLK->num_contributing;i++) {
           fprintf(cam," %c",CLK->jump_list[i].type);
         }
         fprintf(cam,")\n");
         break;
         }
      case 2:
         // the link was already reflected when "other" contig was examined
         break;
      default:
         break;
      }
      WriteProtoMesg_AS(stdout,pmesg);
      linkaddress = GetLinkAddress(links,hasclk++);
      if (linkaddress == NULL ) break; 
   }
   return 1;
}

int main(int c, char **argv)
{ GenericMesg *pmesg;
  int i;
  PHashTable_AS *asmindex;
  PHashTable_AS *asmcount;
  PHashTable_AS *linkptr;
  PHashValue_AS value;
  VA_TYPE(LinkAddress) *links= CreateVA_LinkAddress(200000);
  FILE *input;
  FILE *index;
  FILE *asm;
  FILE *celamy;
  FILE *clkfile;
  LinkAddress linkaddress;
  CDS_UID_t target_object;
  CDS_UID_t uid;
  int lookup_rc,hash_rc;
  long offset;
  input = stdin;
  asmindex = CreatePHashTable_AS(1000000,NULL);
  asmcount = CreatePHashTable_AS(1000000,NULL);
  linkptr = CreatePHashTable_AS(1000000,NULL);
  index = fopen("asm.index","r");
  asm = fopen("asm.asm","r");
  clkfile = fopen("clk.offsets","r");
  uid = 0;
  while ( fscanf(clkfile,F_UID " %ld",&linkaddress.uid,&linkaddress.offset) == 2 ) {
     if ( linkaddress.uid != uid ) {
        value.IID = GetNumLinkAddresss(links);
        hash_rc = InsertInPHashTable_AS(&linkptr,1, (CDS_UID_t) linkaddress.uid, &value, FALSE,FALSE);
        if ( hash_rc != HASH_SUCCESS) {
          fprintf(stderr,"Failure to insert ident " F_UID " in hashtable\n",linkaddress.uid); 
        }
     }
     AppendVA_LinkAddress(links,&linkaddress);
     uid = linkaddress.uid;
  }
  offset = CDS_FTELL(stdin); 
  while ( fscanf(index,F_UID " %ld",&uid,&offset) == 2 ) {
     value.IID = offset;
     hash_rc = InsertInPHashTable_AS(&asmindex,1, (CDS_UID_t) uid, &value, FALSE,FALSE);
     if ( hash_rc != HASH_SUCCESS) {
       fprintf(stderr,"Failure to insert ident " F_UID " in hashtable\n",uid); 
     }
     value.IID = 0;
     hash_rc = InsertInPHashTable_AS(&asmcount,1, (CDS_UID_t) uid, &value, FALSE,FALSE);
  }
  target_object = STR_TO_UID(argv[1],(char **) NULL,10);

  lookup_rc = LookupInPHashTable_AS (asmindex, 1, 
	   (CDS_UID_t) target_object, 
	   &value);
  offset = value.IID;
  pmesg = grab_message(asm,offset);
  //CDS_FSEEK(asm,(off_t) offset,SEEK_SET);
  
  //ReadProtoMesg_AS(asm,&pmesg);
  WriteProtoMesg_AS(stdout,pmesg);
  if (pmesg->t ==MESG_CCO) {
      SnapUnitigMesg *unitig = (SnapUnitigMesg *) pmesg->m;
  //    fprintf(stdout,F_UID " %ld\n",unitig->eaccession,offset);
  } 
  if (pmesg->t ==MESG_CCO) {
      SnapConConMesg  *contig = (SnapConConMesg *) pmesg->m;
      int scaffold_coord=0;
      celamy = fopen("scaffold.cam","w");
      DumpCelamyColors(celamy);
      CelamyContig(celamy,contig,&scaffold_coord,SEQ_REVERSE,asmcount);
  }     
  if (pmesg->t ==MESG_SCF) {
      SnapScaffoldMesg *scaff = (SnapScaffoldMesg *) pmesg->m;
      SnapContigPairs *pairs = (SnapContigPairs *) malloc(scaff->num_contig_pairs*sizeof(SnapContigPairs));
      SnapConConMesg  *contig;
      CDS_UID_t scaff_uid = scaff->eaccession;
      int reverse;
      int num_pairs = scaff->num_contig_pairs;
      int num_contigs = num_pairs+1;
      int scaffold_coord=0;
      for (i=0;i< num_pairs;i++){
          pairs[i] = scaff->contig_pairs[i];
//          fprintf(stdout,F_UID " " F_UID " scaffold_" F_UID " %f %f\n",scaff->contig_pairs[i].econtig1,scaff->contig_pairs[i].econtig2,
//                           scaff->eaccession,scaff->contig_pairs[i].mean, scaff->contig_pairs[i].stddev);
      }
      fprintf(stdout,F_UID " %ld\n",scaff->eaccession,offset);
      // copy the scaffold->contig_pairs so that subsequent proto read doesn't trash it
      //  then, lookup each contig and write it.
      celamy = fopen("scaffold.cam","w");
      DumpCelamyColors(celamy);
      target_object = pairs[0].econtig1;
      lookup_rc = LookupInPHashTable_AS (asmindex, 1, 
	   (CDS_UID_t) target_object, &value);
      offset = value.IID;
      pmesg = grab_message(asm,offset);
      if ( pmesg->t != MESG_CCO ) {
           fprintf(stderr,"Unanticipated message type at file location %ld (was expected CCO " F_UID ")\n",offset,target_object);
      }
      contig = pmesg->m;
      reverse = (pairs[0].orient == AS_NORMAL || pairs[0].orient == AS_INNIE )?0:1;
      CelamyContig(celamy,contig,&scaffold_coord,reverse,asmcount);
      for (i=0;i< scaff->num_contig_pairs;i++){
         SnapConConMesg  *contig;
         int reverse;
         CelamyGap(celamy,scaff_uid,i,&pairs[i],&scaffold_coord);
         target_object = pairs[i].econtig2;
         lookup_rc = LookupInPHashTable_AS (asmindex, 1, 
	   (CDS_UID_t) target_object, &value);
         offset = value.IID;
         pmesg = grab_message(asm,offset);
         if ( pmesg->t != MESG_CCO ) {
           fprintf(stderr,"Unanticipated message type at file location %ld (was expected CCO " F_UID ")\n",offset,target_object);
         }
         contig = pmesg->m;
         reverse = (pairs[i].orient == AS_NORMAL || pairs[i].orient == AS_OUTTIE )?0:1;
         CelamyContig(celamy,contig,&scaffold_coord,reverse,asmcount);
      } 
      // Print scaffold LNK message.
      fprintf(celamy,"LNK: ");
      fprintf(celamy,F_UID "Ctg ",pairs[0].econtig1);
      for (i=0;i<num_pairs;i++) {
	  fprintf(celamy,F_UID "Ctg ",pairs[i].econtig2);
      }
      fprintf(celamy," A0ScaffoldColor \n");

      // Print any pertinent contig link (CLK) messages
      PrintLinks(celamy,pairs[0].econtig1,asmcount,links,linkptr,asm);
      for (i=0;i<num_pairs;i++) {
          PrintLinks(celamy,pairs[i].econtig2,asmcount,links,linkptr,asm);
      }
  }
  fflush(celamy);
  exit (0);

}

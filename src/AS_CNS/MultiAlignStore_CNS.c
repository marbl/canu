
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
/*********************************************************************
   Module:       AS_CNS_MultiAlignStore_CNS.c
   Description:  MultiAlignT and MultiAlignStoreT
                 Data types for managing multi-alignments
   Assumptions:  libAS_UTL.a
 *********************************************************************/

static char CM_ID[] = "$Id: MultiAlignStore_CNS.c,v 1.26 2007-02-14 07:20:09 brianwalenz Exp $";


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_PER_SafeIO.h"
#include "UtilsREZ.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignStore_CNS.h"
#include "Array_CNS.h"
#include "AS_PER_encodeSequenceQuality.h"

VA_DEF(SnapMultiPos)
VA_DEF(UnitigPos)

MultiAlignT *
RevcomplMultiAlignT(MultiAlignT *ma)
{
  MultiAlignT *new_ma = CloneMultiAlignT(ma);
  int length = GetMultiAlignLength(ma);
  int num_frags=GetNumIntMultiPoss(ma->f_list);
  IntMultiPos *frags=GetIntMultiPos(ma->f_list,0);
  IntMultiPos *new_frags=GetIntMultiPos(new_ma->f_list,0);
  int num_unitigs=GetNumIntElementPoss(ma->u_list);
  IntUnitigPos *unitigs=GetIntUnitigPos(ma->u_list,0);
  IntUnitigPos *new_unitigs=GetIntUnitigPos(new_ma->u_list,0);
  int i,j;
  char *consensus = Getchar(new_ma->consensus,0);
  char *quality = Getchar(new_ma->quality,0);
  
  SequenceComplement(consensus, quality);
  new_ma->id = ma->id;
  for (i=0;i<num_frags;i++) {
    int jpos=frags[i].delta_length-1;
    int fraglen=frags[i].position.end - frags[i].position.bgn;
    if ( new_frags[i].ident != frags[i].ident ) {
        fprintf(stderr, "Expecting same fragment but have %d and %d\n",new_frags[i].ident, frags[i].ident);
    }
    if ( fraglen < 0 ) fraglen = -fraglen;
    fraglen-=(frags[i].delta_length-1);
    new_frags[i].position.bgn = length - frags[i].position.bgn;
    new_frags[i].position.end = length - frags[i].position.end;
    // handle special case of fragment endgaps...
    while ( new_frags[i].delta_length > 0 && frags[i].delta[jpos] >= fraglen ) {
       new_frags[i].delta_length--; jpos--;
    }
    if ( new_frags[i].delta_length > 0 ) {
      new_frags[i].delta = (int32 *)safe_malloc(frags[i].delta_length*sizeof(int32));
      for (j=0;j<new_frags[i].delta_length;j++) {
        new_frags[i].delta[j] = fraglen-frags[i].delta[jpos--]-1;
      }
    } else {
      new_frags[i].delta = NULL;
    }
  }
  for (i=0;i<num_unitigs;i++) {
    int jpos=unitigs[i].delta_length-1;
    int tiglen=unitigs[i].position.end - unitigs[i].position.bgn;
    if ( tiglen < 0 ) tiglen = -tiglen;
    tiglen-=(unitigs[i].delta_length-1);
    new_unitigs[i].position.bgn = length - unitigs[i].position.bgn;
    new_unitigs[i].position.end = length - unitigs[i].position.end;
    new_unitigs[i].delta = (int32 *)safe_malloc(unitigs[i].delta_length*sizeof(int32));
    for (j=0;j<unitigs[i].delta_length;j++) {
      new_unitigs[i].delta[j] = tiglen-unitigs[i].delta[jpos--]-1;
    }
  }
  return new_ma;
}

// Local checker function
static void 
CheckMAValidity(MultiAlignT *ma)
{
  char *consensus = Getchar(ma->consensus,0);
  char *quality = Getchar(ma->quality,0);
  char *c,*q;
  assert(strlen(consensus) == strlen(quality));
  for(c = consensus, q = quality; *c != '\0'; c++, q++){
    switch(tolower(*c)){
    case 'a':
    case 'c':
    case 't':
    case 'g':
    case '-':
    case 'n':
      break;
    default:
      assert(0);
    }
    assert(*q >= '0' && *q <= 'l');
  }
  assert(ma->source_alloc == 0 || ma->source_alloc == 1);
}

/**********************************************************************************************/
static int 
CompareUnitigPos (const void *c1, const void *c2)
{
  IntUnitigPos *u1 = (IntUnitigPos *)c1;
  IntUnitigPos *u2 = (IntUnitigPos *)c2;
  int diff;
  int32 bgn1 = MIN(u1->position.bgn, u1->position.end);
  int32 bgn2 = MIN(u2->position.bgn, u2->position.end);
  int32 end1, end2;

  diff = bgn1 - bgn2;
  if(diff)
    return diff;

  end1 = MAX(u1->position.bgn, u1->position.end);
  end2 = MAX(u2->position.bgn, u2->position.end);

  diff = end2 - end1;
  if(diff)
    return diff;

  return TRUE; // arbitrary

}

/**********************************************************************************************/
void 
MakeCanonicalMultiAlignT(MultiAlignT *ma)
{
  IntUnitigPos *unitigs = GetIntUnitigPos(ma->u_list,0);

#if 0
  int i;
  fprintf(stderr,"* Before sort *");
  for(i = 0; i < GetNumIntUnitigPoss(ma->u_list); i++){
    fprintf(stderr,"* unitig %d [%d,%d]\n",
	    unitigs[i].ident,
	    unitigs[i].position.bgn,
	    unitigs[i].position.end);
  }
#endif
  qsort((void *)unitigs, GetNumIntUnitigPoss(ma->u_list),
	sizeof(IntUnitigPos), CompareUnitigPos);

#if 0
  fprintf(stderr,"* After sort *");
  for(i = 0; i < GetNumIntUnitigPoss(ma->u_list); i++){
    fprintf(stderr,"* unitig %d [%d,%d]\n",
	    unitigs[i].ident,
	    unitigs[i].position.bgn,
	    unitigs[i].position.end);
  }
#endif

}


/*
 * Create a MultiAlign object from a protoio IntUnitigMesg, as received from CGB
 */

MultiAlignT *
CreateMultiAlignT(void)
{
  MultiAlignT *ma = (MultiAlignT *)safe_malloc(sizeof(MultiAlignT));
  ma->consensus = NULL;
  ma->quality = NULL;
  ma->delta = NULL;
  ma->f_list = NULL;
  ma->v_list = NULL;
  ma->udelta = NULL;
  ma->u_list = NULL;

  return ma;
}

MultiAlignT *
CreateEmptyMultiAlignT(void)
{
  MultiAlignT *ma = (MultiAlignT *)safe_malloc(sizeof(MultiAlignT));
  ma->consensus = CreateVA_char(0);
  ma->quality = CreateVA_char(0);
  ma->delta = CreateVA_int32(0);;
  ma->f_list = CreateVA_IntMultiPos(0);
  ma->v_list = CreateVA_IntMultiVar(0);
  ma->udelta = CreateVA_int32(0);
  ma->u_list = CreateVA_IntUnitigPos(0);

  return ma;
}

/* Create a clone of a MultiAlignT object */
MultiAlignT *
CloneMultiAlignT(MultiAlignT *ma)
{
  MultiAlignT *newma = CreateEmptyMultiAlignT();
  CopyMultiAlignT(newma, ma);
  return newma;
}

void 
CopyMultiAlignT(MultiAlignT *newma, MultiAlignT *ma)
{

  if(newma->consensus == NULL){
    newma->consensus = Clone_VA(ma->consensus);
    newma->quality = Clone_VA(ma->quality);
    // Save the delta pointers as offset from base of delta array
    newma->f_list = Clone_VA(ma->f_list);
    newma->v_list = Clone_VA(ma->v_list);
    newma->u_list = Clone_VA(ma->u_list);
    newma->delta = Clone_VA(ma->delta);
    newma->udelta = Clone_VA(ma->udelta);
  }else{
    ReuseClone_VA(newma->consensus,ma->consensus);
    ReuseClone_VA(newma->quality,ma->quality);
    // Save the delta pointers as offset from base of delta array
    ReuseClone_VA(newma->f_list,ma->f_list);
    ReuseClone_VA(newma->v_list,ma->v_list);
    ReuseClone_VA(newma->u_list,ma->u_list);
    ReuseClone_VA(newma->delta,ma->delta);
    ReuseClone_VA(newma->udelta, ma->udelta);
  }
  newma->forced = ma->forced;
  //  newma->id = ma->id;
  newma->refCnt = 0;
  newma->source_alloc = ma->source_alloc;
  {/* Adjust the delta pointers in the clone */
    int i;
    char *old_source, *old_var_seq, *old_nr_conf_alleles, *old_weights;
    int src_len;
    int32 *oldbase = Getint32(ma->delta, 0);
    int32 *newbase = Getint32(newma->delta, 0);
    int numf = GetNumIntMultiPoss(ma->f_list);
    int numv = GetNumIntMultiVars(ma->v_list);
    for(i = 0; i < numf; i++){
      IntMultiPos *npos = GetIntMultiPos(newma->f_list,i);
      int offset = (npos->delta - oldbase);
      npos->delta = newbase + offset;
    }
    for(i = 0; i < numv; i++)
    {
      IntMultiVar *nvar = GetIntMultiVar(newma->v_list,i);
      old_nr_conf_alleles = nvar->nr_conf_alleles;
      old_weights         = nvar->weights;
      old_var_seq         = nvar->var_seq;
      nvar->var_seq = (char *) safe_malloc((strlen(old_var_seq)+1)*sizeof(char));
      nvar->nr_conf_alleles = (char *) safe_malloc((strlen(old_nr_conf_alleles)+1)
           *sizeof(char));
      nvar->weights = (char *) safe_malloc((strlen(old_weights)+1)*sizeof(char));
      strcpy(nvar->nr_conf_alleles, old_nr_conf_alleles);
      strcpy(nvar->weights,         old_weights);
      strcpy(nvar->var_seq,         old_var_seq);
    }
  }
  {/* Adjust the delta pointers in the clone */
    int i;
    int32 *oldbase = Getint32(ma->udelta, 0);
    int32 *newbase = Getint32(newma->udelta, 0);
    int32 numu = GetNumIntUnitigPoss(ma->u_list);
    for(i = 0; i < numu; i++){
      IntUnitigPos *npos = GetIntUnitigPos(newma->u_list,i);
      int offset = (npos->delta - oldbase);
      npos->delta = newbase + offset;
    }
  }
}

// Create Surrogate
MultiAlignT *
CloneSurrogateOfMultiAlignT(MultiAlignT *oldMA, int32 newNodeID)
{
  MultiAlignT *newma = CreateMultiAlignT();
  IntUnitigPos *u;
  int32 oldLength = GetMultiAlignLength(oldMA);

  
 // We have a single unitig
   assert(GetNumIntUnitigPoss(oldMA->u_list) == 1);
#if 0
   // Surrogate has gapped consensus sequence
  newma->consensus = Clone_VA(oldMA->consensus);
  newma->quality = Clone_VA(oldMA->quality);
#else
  // Surrogate has UNGAPPED consensus sequence.  As fragments
  // get added, gaps will rematerialize

  newma->id = newNodeID;
  newma->consensus = CreateVA_char(oldLength);
  newma->quality = CreateVA_char(oldLength);
  GetMultiAlignUngappedConsensus(oldMA, newma->consensus, newma->quality);

  fprintf(stderr,"* oldMA has length:%d newma has length:%d\n",
	  oldLength, GetMultiAlignLength(newma));
#endif
  newma->delta = CreateVA_int32(0);
  newma->f_list = CreateVA_IntMultiPos(0);
  newma->v_list = CreateVA_IntMultiVar(0);
  newma->udelta = CreateVA_int32(0);
  newma->u_list = Clone_VA(oldMA->u_list);
  //  newma->id = newNodeID;
  newma->forced = 0;
  newma->refCnt = 0;
  newma->source_alloc = 0;
  u = GetIntUnitigPos(newma->u_list, 0);
  u->ident = newNodeID;
  u->position.end = GetMultiAlignLength(newma); // ungapped!
  u->delta_length = 0; // ungapped!
  u->delta = NULL; // ungapped!
  return newma;

}




IntUnitigPos *
GetAendUnitigPos(MultiAlignT *ma)
{
  return GetIntUnitigPos(ma->u_list,0);
}


IntUnitigPos *
GetBendUnitigPos(MultiAlignT *ma)
{
  int found = FALSE;
  int i;
  IntUnitigPos *result = NULL;
  long length = GetMultiAlignLength(ma); 
  int numu = GetNumIntUnitigPoss(ma->u_list);
  for(i = numu -1; i>= 0; i--)
    {
      IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);
      if( MAX(pos->position.bgn,pos->position.end) == length )
	{
	  result = pos;
	  found = TRUE;
	  break;
	}
    }
  assert(found);
  return result;
}


int32 
GetMultiAlignLength(MultiAlignT *ma)
{
  return (int32) GetNumchars(ma->consensus) - 1; // don't include the space for the null character
}


int32 
GetMultiAlignUngappedLength(MultiAlignT *ma)
{
  int32 ungappedLength = 0;
  char *consensus = Getchar(ma->consensus,0);
  char *c;

  for(c = consensus;
      *c != '\0';
      c++){

    if(*c != '-')
      ungappedLength++;
  }
  return ungappedLength;
}

MultiAlignT *
CreateMultiAlignTFromIUM(IntUnitigMesg *ium, int localID, int sequenceOnly)
{
/* if localID = -1, interpret the  frag source fields as strings , and copy them */
/* if localID = -2, preserve the special hijacked source fields in the frag messages */
/* if localID > 0, assign the frag source chars their special hijacked values, keyed from localID */
  int cfr, cvr, deltai;
  MultiAlignT *ma = (MultiAlignT *)safe_malloc(sizeof(MultiAlignT));
  char *ptr;
  IntUnitigPos unitigPos;
  long localFragID = localID;
  int delta_len=0;

  if (ium->length != strlen(ium->consensus))
  {
    fprintf(stderr, "Reported Length of IUM %d (%d) doesnt matches strlen (%d)\n",
            ium->iaccession, ium->length, strlen(ium->consensus));
  }

  assert(ium->length == strlen(ium->consensus));
  assert(ium->length == strlen(ium->quality));

  ma->id = ium->iaccession;
  if(ium->forced){
    SetMultiAlignForced(ma, TRUE);
    fprintf(stderr,"*** WARNING: IUM with accession %d has forced = TRUE\n", ium->iaccession);
  }

  /* We need real quality values...abort if we don't get them */
  {
    int ok = FALSE;
    for(ptr = ium->quality; *ptr != '\0'; ptr++){
      if(*ptr != '0'){
	ok = TRUE;
	break;
      }
    }
    if(!ok){ 
      fprintf(stderr,"* IUM with accession %d has bogus quality string...rerun consensus!\n",
	      ium->iaccession);
      if(!ium->forced){
	fprintf(stderr,"* -----> IUM DOES NOT have  forced flag set...exiting...");
	exit(1);
      }
      fprintf(stderr,"* -----> IUM has forced flag set...continuing...");
    }

  }

  ma->refCnt = 0; // set on insertion into a store
  if (localFragID == -2) {
    ma->source_alloc = 0;
  } 
  else if (localFragID < 0) {
    ma->source_alloc = 1;
  }
  else {
    ma->source_alloc = 0;
  }
  ma->consensus = CreateVA_char(ium->length + 1);
  EnableRangeVA_char(ma->consensus, ium->length + 1);

  ma->quality = CreateVA_char(ium->length + 1);
  EnableRangeVA_char(ma->quality, ium->length + 1);

  ma->forced = ium->forced;
  for(cfr = 0; cfr < ium->num_frags; cfr++){
    delta_len+=ium->f_list[cfr].delta_length;
  }

  if( ! sequenceOnly )
  {
      ma->delta = CreateVA_int32(delta_len);
      
      ma->f_list = CreateVA_IntMultiPos(ium->num_frags);
      ma->v_list = CreateVA_IntMultiVar(ium->num_vars); 
      ma->u_list = CreateVA_IntUnitigPos(1);
      

      for(cfr = 0,delta_len=0; cfr < ium->num_frags; cfr++)
      {
	IntMultiPos *cfr_mesg = ium->f_list + cfr;
	IntMultiPos tmp;
	

	
	tmp = *cfr_mesg; // Copy everything and then fix the pointers
	//	tmp.type = cfr_mesg->type;
	//	tmp.ident = cfr_mesg->ident;
	if (localFragID == -2) {
          tmp.sourceInt = cfr_mesg->sourceInt;
			 ma->source_alloc = 0;
	} else if (localFragID < 0) {
          tmp.sourceInt = INT_MAX;
          ma->source_alloc = 0;
	} else {
	  tmp.sourceInt = localFragID++;
	}
	//	tmp.position = cfr_mesg->position;
	//	tmp.delta_length = cfr_mesg->delta_length;
	for (deltai=0;deltai<cfr_mesg->delta_length;deltai++) {
	  //      fprintf(stderr,"* deltai = %d delta = %d\n",
	  //	      deltai, cfr_mesg->delta[deltai]);
	  AppendVA_int32(ma->delta,cfr_mesg->delta + deltai);
	} 
	tmp.delta = Getint32(ma->delta,delta_len);
	delta_len+=cfr_mesg->delta_length;
	SetIntMultiPos(ma->f_list, cfr, &tmp);
      }
      for(cvr = 0; cvr < ium->num_vars; cvr++)
      {
         IntMultiVar *cvr_mesg = ium->v_list + cvr; 
         IntMultiVar tmp;
         
         tmp = *cvr_mesg; 
         tmp.position = cvr_mesg->position;
         SetIntMultiVar(ma->v_list, cvr, &tmp);
      }
  }

  ptr = Getchar(ma->consensus,0);
  strcpy(ptr, ium->consensus);

  ptr = Getchar(ma->quality,0);
  strcpy(ptr, ium->quality);


  if( ! sequenceOnly )
    {
      ma->udelta = CreateVA_int32(0);
      { int32 ui,deltai;
      /* Add a multipos for this Unitig */

      unitigPos.type = AS_OTHER_UNITIG;  // Jason, 7/01.

      unitigPos.ident = ium->iaccession;
      unitigPos.position.bgn = 0;
      unitigPos.position.end = GetMultiAlignLength(ma);
      for (ui=0,deltai=0;ui<GetMultiAlignLength(ma);ui++) {
	if (ium->consensus[ui] == '-') {
	  Appendint32(ma->udelta,&deltai);
	} else {
	  deltai++;
	}
      }
      unitigPos.delta_length = GetNumint32s(ma->udelta);
      unitigPos.delta = Getint32(ma->udelta,0);
      delta_len+=GetNumint32s(ma->udelta);
      AppendIntUnitigPos(ma->u_list, &unitigPos);
      }
    }

#if 0
  fprintf(stderr,"* Added to u_list: ident:%d [%d,%d]\n",
	  unitigPos.ident, unitigPos.position.bgn, unitigPos.position.end);

  fprintf(stderr,"* Created MultiAlign for sequence of length %d with %d fragments\n",
	  ium->length, ium->num_frags);

  //  fprintf(stderr,"* MA %d ma->delta numElements = %ld size=%ld\n",
  //	  ium->iaccession, ma->delta->numElements, ma->delta->sizeofElement);
#endif
  if( ! sequenceOnly )
    assert(ium->num_frags == GetNumIntMultiPoss(ma->f_list));

  
  CheckMAValidity(ma);
  return ma;

}



MultiAlignT *
CreateMultiAlignTFromICM(IntConConMesg *icm, int localID, int sequenceOnly)
{
/* if localID is negative, use NULL source field, else, use source for localID */
  int cfr, cvr, deltai;
  MultiAlignT *ma = (MultiAlignT *)safe_malloc(sizeof(MultiAlignT));
  char *ptr;
  IntUnitigPos unitigPos;
  long localFragID = localID;
  int delta_len=0;

  assert(icm->length == strlen(icm->consensus));
  assert(icm->length == strlen(icm->quality));
  
  ma->id = icm->iaccession;

  ma->forced = icm->forced;
  ma->refCnt = 0; // set on insertion into a store
  if ( localID < 0 ) {
     ma->source_alloc = 1;
  } else {
     ma->source_alloc = 0;
  }
  ma->consensus = CreateVA_char(icm->length + 1);
  EnableRangeVA_char(ma->consensus, icm->length + 1);

  ma->quality = CreateVA_char(icm->length + 1);
  EnableRangeVA_char(ma->quality, icm->length + 1);

  for(cfr = 0; cfr < icm->num_pieces; cfr++){
    delta_len+=icm->pieces[cfr].delta_length;
  }

  if( ! sequenceOnly )
    {
      ma->delta = CreateVA_int32(delta_len);
      ma->f_list = CreateVA_IntMultiPos(icm->num_pieces);
      ma->udelta = CreateVA_int32(0);
      ma->u_list = CreateVA_IntUnitigPos(0);
      ma->v_list = CreateVA_IntMultiVar(icm->num_vars);
      

      for(cfr = 0,delta_len=0; cfr < icm->num_pieces; cfr++){
	IntMultiPos *cfr_mesg = icm->pieces + cfr;
	IntMultiPos tmp;
	

	
	tmp.type = cfr_mesg->type;
	tmp.ident = cfr_mesg->ident;
	/* if (localFragID == -2) {
             tmp.source = cfr_mesg->source;
			 } else */
	if (localFragID < 0) {
          tmp.sourceInt = cfr_mesg->sourceInt;
	} else {
	  tmp.sourceInt = localFragID++;
	}
	tmp.position = cfr_mesg->position;
	tmp.contained = cfr_mesg->contained;
	tmp.delta_length = cfr_mesg->delta_length;
	for (deltai=0;deltai<cfr_mesg->delta_length;deltai++) {
	  //      fprintf(stderr,"* deltai = %d delta = %d\n",
	  //	      deltai, cfr_mesg->delta[deltai]);
	  AppendVA_int32(ma->delta,cfr_mesg->delta + deltai);
	} 
	tmp.delta = Getint32(ma->delta,delta_len);
	delta_len+=cfr_mesg->delta_length;
	SetIntMultiPos(ma->f_list, cfr, &tmp);
      }
      for(cvr = 0; cvr < icm->num_vars; cvr++)
      {
         IntMultiVar *cvr_mesg = icm->v_list + cvr;
         IntMultiVar tmp;
         int         na  = (cvr_mesg->num_conf_alleles < 2) ? 2 : cvr_mesg->num_conf_alleles;

         tmp = *cvr_mesg;
         tmp.position         = cvr_mesg->position;
         tmp.num_reads        = cvr_mesg->num_reads;
         tmp.num_conf_alleles = cvr_mesg->num_conf_alleles;
         tmp.anchor_size      = cvr_mesg->anchor_size;
         tmp.var_length       = cvr_mesg->var_length;
         tmp.nr_conf_alleles  = (char *) safe_malloc(4*sizeof(char)*na + 1);
         tmp.weights          = (char *) safe_malloc(7*sizeof(char)*na + 1);
         tmp.var_seq          = (char *) safe_malloc((cvr_mesg->var_length+1)*sizeof(char)*na + 1);
         strcpy(tmp.nr_conf_alleles, cvr_mesg->nr_conf_alleles);
         strcpy(tmp.weights, cvr_mesg->weights);
         strcpy(tmp.var_seq, cvr_mesg->var_seq);
         SetIntMultiVar(ma->v_list, cvr, &tmp);
         safe_free(tmp.nr_conf_alleles);
         safe_free(tmp.weights);
         safe_free(tmp.var_seq);
      }
    }

  ptr = Getchar(ma->consensus,0);
  strcpy(ptr, icm->consensus);

  ptr = Getchar(ma->quality,0);
  strcpy(ptr, icm->quality);


  if( ! sequenceOnly )
    {
      { int32 ui;
      /* Add a unitigpos for each Unitig in icm->unitigs*/
      /* not authentic, since icm doesn't retain deltas */
      for (ui = 0;ui<icm->num_unitigs;ui++) {
        unitigPos.type = icm->unitigs[ui].type;
        unitigPos.ident = icm->unitigs[ui].ident;
        unitigPos.position = icm->unitigs[ui].position;
        unitigPos.delta_length = 0;
        unitigPos.delta = NULL;
        AppendIntUnitigPos(ma->u_list, &unitigPos);
      }
      }
    }

  if( ! sequenceOnly )
    assert(icm->num_pieces == GetNumIntMultiPoss(ma->f_list));

  
  CheckMAValidity(ma);
  return ma;
}

MultiAlignT *
CreateMultiAlignTFromCCO(SnapConConMesg *cco, int localID, int sequenceOnly)
{
/* if localID is negative, use NULL source field, else, use source for localID */
  int cfr, cvr, deltai;
  MultiAlignT *ma = (MultiAlignT *)safe_malloc(sizeof(MultiAlignT));
  char *ptr;
  UnitigPos unitigPos;
  long localFragID = localID;
  int delta_len=0;

  assert(cco->length == strlen(cco->consensus));
  assert(cco->length == strlen(cco->quality));
  
  ma->id = cco->iaccession;

  ma->forced = cco->forced;
  ma->refCnt = 0; // set on insertion into a store
  if ( localID < 0 ) {
     ma->source_alloc = 1;
  } else {
     ma->source_alloc = 0;
  }
  ma->consensus = CreateVA_char(cco->length + 1);
  EnableRangeVA_char(ma->consensus, cco->length + 1);

  ma->quality = CreateVA_char(cco->length + 1);
  EnableRangeVA_char(ma->quality, cco->length + 1);

  for(cfr = 0; cfr < cco->num_pieces; cfr++){
    delta_len+=cco->pieces[cfr].delta_length;
  }

  if( ! sequenceOnly )
    {
      ma->delta = CreateVA_int32(delta_len);
      ma->f_list = CreateVA_SnapMultiPos(cco->num_pieces);
      ma->v_list = CreateVA_IntMultiVar(cco->num_vars);  
      ma->udelta = CreateVA_int32(0);
      ma->u_list = CreateVA_UnitigPos(cco->num_unitigs);


      for(cfr = 0,delta_len=0; cfr < cco->num_pieces; cfr++){
	SnapMultiPos *cfr_mesg = cco->pieces + cfr;
	SnapMultiPos tmp;
	
	int length = cfr_mesg->delta_length * sizeof(int32);
	
	tmp.type = cfr_mesg->type;
	tmp.eident = cfr_mesg->eident;
	/* if (localFragID == -2) {
             tmp.source = cfr_mesg->source;
			 } else */
	if (localFragID < 0) {
          int32 src_len;
          if (cfr_mesg->source) {
             src_len =  strlen(cfr_mesg->source);
             tmp.source = (char *) safe_malloc((src_len+1)*sizeof(char));
             strcpy(tmp.source,cfr_mesg->source);
	  } else {
             tmp.source = cfr_mesg->source;
          }
	} else {
	  tmp.source = (char *)localFragID++;
	}
	tmp.position = cfr_mesg->position;
	//	tmp.contained = cfr_mesg->contained;
	tmp.delta_length = cfr_mesg->delta_length;
	for (deltai=0;deltai<cfr_mesg->delta_length;deltai++) {
	  //      fprintf(stderr,"* deltai = %d delta = %d\n",
	  //	      deltai, cfr_mesg->delta[deltai]);
	  AppendVA_int32(ma->delta,cfr_mesg->delta + deltai);
	} 
	tmp.delta = Getint32(ma->delta,delta_len);
	delta_len+=cfr_mesg->delta_length;
	SetSnapMultiPos(ma->f_list, cfr, &tmp);
      }
      for(cvr = 0; cvr < cco->num_vars; cvr++)
      {
         IntMultiVar *cvr_mesg = cco->vars + cvr;
         IntMultiVar tmp;
         int         na  = (cvr_mesg->num_conf_alleles < 2) ? 2 : cvr_mesg->num_conf_alleles;

         tmp = *cvr_mesg;
         tmp.position         = cvr_mesg->position;
         tmp.num_reads        = cvr_mesg->num_reads;
         tmp.num_conf_alleles = cvr_mesg->num_conf_alleles;
         tmp.anchor_size      = cvr_mesg->anchor_size;
         tmp.var_length       = cvr_mesg->var_length;
         tmp.nr_conf_alleles  = (char *) safe_malloc(4*sizeof(char)*na + 1);
         tmp.weights          = (char *) safe_malloc(7*sizeof(char)*na + 1);
         tmp.var_seq          = (char *) safe_malloc((cvr_mesg->var_length+1)*sizeof(char)*na + 1);
         strcpy(tmp.nr_conf_alleles, cvr_mesg->nr_conf_alleles);
         strcpy(tmp.weights, cvr_mesg->weights);
         strcpy(tmp.var_seq, cvr_mesg->var_seq);
         SetIntMultiVar(ma->v_list, cvr, &tmp);
      }
    }

  ptr = Getchar(ma->consensus,0);
  strcpy(ptr, cco->consensus);

  ptr = Getchar(ma->quality,0);
  strcpy(ptr, cco->quality);


  if( ! sequenceOnly )
    {
      { int32 ui,deltai;
      /* Add a unitigpos for each Unitig in cco->unitigs*/
      /* not authentic, since cco doesn't retain deltas */
      for (ui = 0;ui<cco->num_unitigs;ui++) {
        unitigPos.type = cco->unitigs[ui].type;
        unitigPos.eident = cco->unitigs[ui].eident;
        unitigPos.position = cco->unitigs[ui].position;
        unitigPos.delta_length = 0;
        unitigPos.delta = NULL;
        AppendUnitigPos(ma->u_list, &unitigPos);
      }
      }
    }

  if( ! sequenceOnly )
    assert(cco->num_pieces == GetNumIntMultiPoss(ma->f_list));

  
  CheckMAValidity(ma);
  return ma;
}



/********************************************************************************/
int32 
AddReferenceMultiAlignT(MultiAlignT *ma)
{
  ma->refCnt++;
  return ma->refCnt;
}

int32 
RemoveReferenceMultiAlignT(MultiAlignT *ma)
{
  ma->refCnt--;
  return ma->refCnt;
}



void 
DeleteMultiAlignT(MultiAlignT *ma)
{
  int i;

  RemoveReferenceMultiAlignT(ma);
  if(ma->refCnt > 0)
    return;
  //  fprintf(stderr,"* Freeing ma at 0x%x source_alloc = %d\n",
  //	  ma, ma->source_alloc);

  if (ma->source_alloc) {
    // make sure space alloced to hold source is freed
    IntMultiPos *t = NULL;
    IntMultiVar *v = NULL;
    int n_frags=GetNumIntMultiPoss(ma->f_list);
    int n_vars=GetNumIntMultiVars(ma->v_list);
    if (n_frags > 0) t=GetIntMultiPos(ma->f_list,0);
    if (n_vars > 0) v=GetIntMultiVar(ma->v_list, 0);
    for (i=0;i<n_vars;i++){
       if (v->nr_conf_alleles) safe_free(v->nr_conf_alleles);
       if (v->weights) safe_free(v->weights);
       if (v->var_seq) safe_free(v->var_seq);
       v++;
    }
  }   
  DeleteVA_char(ma->consensus);
  DeleteVA_char(ma->quality);
  DeleteVA_int32(ma->udelta);
  DeleteVA_IntUnitigPos(ma->u_list);
  DeleteVA_int32(ma->delta);
  DeleteVA_IntMultiPos(ma->f_list);
  DeleteVA_IntMultiVar(ma->v_list);
  safe_free(ma);
}

/********************************************************************************/
// Persistence
/********************************************************************************/
static void 
SaveReferenceMultiAlignTToStream(MultiAlignT *ma, FILE *stream)
{
  //char reference = TRUE;
  int32 reference = ma->id;
  int status;
  char isPresent = (ma != NULL);

  // Sentinel to say this is non-null
    status = safeWrite(stream, &isPresent, sizeof(char));
    assert(status == FALSE);

    if(!isPresent)
      return;

    // Sentinel to say this is a reference
    //    fprintf(stderr,"* Saving reference to ma %d\n", reference);

    status = safeWrite(stream, &reference, sizeof(int32));
    assert(status == FALSE);


}


size_t 
SaveMultiAlignTToStream(MultiAlignT *ma, FILE *stream)
{
  int i;
  int status;
  size_t totalSize = 0;
  int32 reference = NULLINDEX;
  char isPresent = (ma != NULL);
  //  if(!isPresent)
    //  fprintf(stderr,"* SaveMultiAlignTToStream   NULL MultiAlignT!!!\n");

  // Sentinel to say this is non-null
  totalSize++;
    status = safeWrite(stream, &isPresent, sizeof(char));
    assert(status == FALSE);

    if(!isPresent)
      return(sizeof(char));;

  // Sentinel to say this is a real one
    totalSize += sizeof(int32);
    status = safeWrite(stream, &reference, sizeof(int32));
    assert(status == FALSE);

    //  CheckMAValidity(ma);

  // Save the delta pointers as offset from base of delta array
  {
    int32 *base = Getint32(ma->delta, 0);
    //    fprintf(stderr,"* base = 0x%x\n",base);
    int numf = GetNumIntMultiPoss(ma->f_list);
    for(i = 0; i < numf; i++){
      IntMultiPos *pos = GetIntMultiPos(ma->f_list,i);
      long offset = (pos->delta - base);
      //      fprintf(stderr,"* %d delta:%x 0x%x \n", i, offset, pos->delta);
      pos->delta = (int32 *)offset;
    }
  }
  // Save the delta pointers as offset from base of delta array
  {
    int32 *base = Getint32(ma->udelta, 0);
    int numu = GetNumIntUnitigPoss(ma->u_list);
    for(i = 0; i < numu; i++){
      IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);
      long offset = (pos->delta - base);
      pos->delta = (int32 *)offset;
    }
  }
  totalSize += CopyToFileVA_char(ma->consensus, stream);
  totalSize += CopyToFileVA_char(ma->quality, stream);
  totalSize += CopyToFileVA_int32(ma->delta, stream);
  totalSize += CopyToFileVA_IntMultiPos(ma->f_list, stream);
  totalSize += CopyToFileVA_IntMultiVar(ma->v_list, stream);
  totalSize += CopyToFileVA_int32(ma->udelta, stream);
  totalSize += CopyToFileVA_IntMultiPos(ma->u_list, stream);
  totalSize += (3 * sizeof(int32));
  //  fprintf(stderr,"*totalSize is %d\n", totalSize);
  status = safeWrite(stream, &ma->forced, sizeof(int32));
  status = safeWrite(stream, &ma->id, sizeof(int32));
  status = safeWrite(stream, &ma->source_alloc, sizeof(int32));
  assert(status == FALSE);
 //  fprintf(stderr,"* ma %d start:%ld total:%ld\n",
  //	  ma->id, size, CDS_FTELL(stream) - size);
  // Restore the delta pointers since they were saved as offset from base of delta array
  //  fprintf(stderr,"* ma->delta = 0x%x\n", ma->delta);
  {
    int32 *base = Getint32(ma->delta, 0);
    //    fprintf(stderr,"* base = 0x%x\n",base);
    int32 numf = GetNumIntMultiPoss(ma->f_list);
    for(i = 0; i < numf; i++){
      IntMultiPos *pos = GetIntMultiPos(ma->f_list,i);
      //      fprintf(stderr,"* %d delta:%d", i,pos->delta);
      pos->delta = base + (long)(pos->delta);
      //      fprintf(stderr," after 0x%x\n", pos->delta);
    }
  }
  // Restore the delta pointers since they were saved as offset from base of delta array
  //  fprintf(stderr,"* ma->udelta = 0x%x\n", ma->udelta);
  {
    int32 *base = Getint32(ma->udelta, 0);
    int32 numu = GetNumIntUnitigPoss(ma->u_list);
    for(i = 0; i < numu; i++){
      IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);
      pos->delta = base + (long)(pos->delta);
    }
  }
  return totalSize;

}

/********************************************************************************/
MultiAlignT *
LoadMultiAlignTFromStream(FILE *stream, int32 *reference)
{
  int status;
  int i;
  MultiAlignT *ma;
  char isPresent; 

  // Sentinel to say this is non-null
    status = safeRead(stream, &isPresent, sizeof(char));
    assert(status == FALSE);

    if(!isPresent){
      //      fprintf(stderr,"* Read NULL MultiAlignT...returning\n");
    *reference = NULLINDEX;
      return NULL;
    }
  status = safeRead(stream, reference, sizeof(int32));
  assert(status == FALSE);

  if(*reference != NULLINDEX){
    //    fprintf(stderr,"* ma ?? start:%ld total:%ld\n",
    //	    size, CDS_FTELL(stream) - size);
    return NULL;
  }

  
  ma = CreateMultiAlignT();

  ma->consensus = CreateFromFileVA_char(stream,0);
  ma->quality = CreateFromFileVA_char(stream,0);
  ma->delta = CreateFromFileVA_int32(stream,0);
  ma->f_list = CreateFromFileVA_IntMultiPos(stream,0);
  ma->v_list = CreateFromFileVA_IntMultiVar(stream,0);
  ma->udelta = CreateFromFileVA_int32(stream,0);
  ma->u_list = CreateFromFileVA_IntUnitigPos(stream,0);
  status = safeRead(stream, &ma->forced, sizeof(int32));
  status = safeRead(stream, &ma->id, sizeof(int32));
  status = safeRead(stream, &ma->source_alloc, sizeof(int32));
  assert(status == FALSE);

  // Restore the delta pointers since they were saved as offset from base of delta array
  //  fprintf(stderr,"* ma->delta = 0x%x\n", ma->delta);
  {
    int32 *base = Getint32(ma->delta, 0);
    int numf = GetNumIntMultiPoss(ma->f_list);
    for(i = 0; i < numf; i++){
      IntMultiPos *pos = GetIntMultiPos(ma->f_list,i);
      //      fprintf(stderr,"* %d delta:%d\n", i,pos->delta);
      pos->delta = base + (long)(pos->delta);
    }
  }
  // Restore the udelta pointers since they were saved as offset from base of delta array
  //  fprintf(stderr,"* ma->udelta = 0x%x\n", ma->udelta);
  {
    int32 *base = Getint32(ma->udelta, 0);
    int32 numu = GetNumIntUnitigPoss(ma->u_list);
    for(i = 0; i < numu; i++){
      IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);
      pos->delta = base + (long)(pos->delta);
    }
  }

  ma->refCnt = 0;
  CheckMAValidity(ma);
  return (ma);
}
/******************************************************************/
void 
ReLoadMultiAlignTFromStream(FILE *stream, MultiAlignT *ma, int32 *reference)
{
  int status;
  int i;
  char isPresent; 

  AssertPtr(ma);
  ResetVA_char(ma->consensus);
  ResetVA_char(ma->quality);
  ResetVA_int32(ma->delta);
  ResetVA_int32(ma->udelta);
  ResetVA_IntMultiPos(ma->f_list);
  ResetVA_IntMultiVar(ma->v_list);
  ResetVA_IntUnitigPos(ma->u_list);

  // Sentinel to say this is non-null
    status = safeRead(stream, &isPresent, sizeof(char));
    assert(status == FALSE);

    if(!isPresent){
      //      fprintf(stderr,"* Read NULL MultiAlignT...returning\n");
    *reference = NULLINDEX;
      return;
    }
  status = safeRead(stream, reference, sizeof(int32));
  assert(status == FALSE);

  if(*reference != NULLINDEX){
    //    fprintf(stderr,"* ma ?? start:%ld total:%ld\n",
    //	    size, CDS_FTELL(stream) - size);
    return;
  }

  
  LoadFromFileVA_char(stream,ma->consensus,0);
  LoadFromFileVA_char(stream,ma->quality,0);
  LoadFromFileVA_int32(stream,ma->delta,0);
  LoadFromFileVA_IntMultiPos(stream,ma->f_list,0);
  LoadFromFileVA_IntMultiVar(stream,ma->v_list,0);
  LoadFromFileVA_int32(stream,ma->udelta,0);
  LoadFromFileVA_IntUnitigPos(stream,ma->u_list,0);
  status = safeRead(stream, &ma->forced, sizeof(int32));
  status = safeRead(stream, &ma->id, sizeof(int32));
  status = safeRead(stream, &ma->source_alloc, sizeof(int32));
  assert(status == FALSE);

  // Restore the delta pointers since they were saved as offset from base of delta array
  //  fprintf(stderr,"* ma->delta = 0x%x\n", ma->delta);
  {
    int32 *base = Getint32(ma->delta, 0);
    int32 numf = GetNumIntMultiPoss(ma->f_list);
    for(i = 0; i < numf; i++){
      IntMultiPos *pos = GetIntMultiPos(ma->f_list,i);
      //      fprintf(stderr,"* %d delta:%d\n", i,pos->delta);
      pos->delta = base + (long)(pos->delta);
    }
  }
  // Restore the udelta pointers since they were saved as offset from base of delta array
  //  fprintf(stderr,"* ma->udelta = 0x%x\n", ma->udelta);
  {
    int32 *base = Getint32(ma->udelta, 0);
    int32 numu = GetNumIntUnitigPoss(ma->u_list);
    for(i = 0; i < numu; i++){
      IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);
      pos->delta = base + (long)(pos->delta);
    }
  }

  ma->refCnt = 0;
  CheckMAValidity(ma);
  return;
}

size_t 
GetMemorySize(MultiAlignT *ma)
{
  size_t size;


  if(!ma)
    return 0;

  size = GetMemorySize_VA(ma->consensus) * 2 +
         GetMemorySize_VA(ma->delta) +
         GetMemorySize_VA(ma->f_list) +
         GetMemorySize_VA(ma->v_list) +
         GetMemorySize_VA(ma->udelta) +
         GetMemorySize_VA(ma->u_list);

//{  int i;
//  for (i = 0; i < GetNumIntMultiPoss(ma->f_list); i++){
//    IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
//    size += mp->delta_length;
//  }}

  return size;
}

/**************************************************/
int 
CompareMultiAlignT(MultiAlignT *thisMAT, MultiAlignT *otherMAT)
{
  int diff;

  diff = strcmp(
		Getchar(thisMAT->consensus,0), 
		Getchar(otherMAT->consensus,0));
  if(diff)
    return diff;

  diff = strcmp(
		Getchar(thisMAT->quality,0), 
		Getchar(otherMAT->quality,0));
  if(diff)
    return diff;

  return 0;
}



/* GetMultiAlignUngappedConsensus */
void 
GetMultiAlignUngappedConsensus(MultiAlignT *ma, VA_TYPE(char) *ungappedConsensus, 
  VA_TYPE(char) *ungappedQuality)
{
  char *consensus = Getchar(ma->consensus,0);
  char *quality = Getchar(ma->quality,0);
  char *c, *q;
  char nullChar = '\0';

  CheckMAValidity(ma);

  ResetVA_char(ungappedConsensus);
  ResetVA_char(ungappedQuality);

  for(c = consensus, q = quality;
      *c != '\0';
      c++, q++){

    if(*c == '-')
      continue;

    Appendchar(ungappedConsensus, c);
    Appendchar(ungappedQuality, q);
  }

  // Make sure we finish with a null char
    Appendchar(ungappedConsensus, &nullChar);
    Appendchar(ungappedQuality, &nullChar);
}



void 
GetMultiAlignUngappedConsensusFromInterval(MultiAlignT *ma, 
  SeqInterval gappedInterval,  VA_TYPE(char) *ungappedConsensus, 
  VA_TYPE(char) *ungappedQuality)
{
  int offset = gappedInterval.bgn;
  char *consensus = Getchar(ma->consensus,offset);
  char *quality = Getchar(ma->quality,offset);
  char *c, *q;
  char nullChar = '\0';

  CheckMAValidity(ma);

  ResetVA_char(ungappedConsensus);
  ResetVA_char(ungappedQuality);

  for(c = consensus, q = quality;
      offset < gappedInterval.end && *c != '\0';
      c++, q++, offset++){

    if(*c == '-')
      continue;

    Appendchar(ungappedConsensus, c);
    Appendchar(ungappedQuality, q);
  }

  // Make sure we finish with a null char
    Appendchar(ungappedConsensus, &nullChar);
    Appendchar(ungappedQuality, &nullChar);
}


/* GetMultiAlignUngappedOffsets */
void 
GetMultiAlignUngappedOffsets(MultiAlignT *ma, VA_TYPE(int32) *ungappedOffsets)
{
  char *consensus = Getchar(ma->consensus,0);
  // char *quality = Getchar(ma->quality,0);
  char *c;
  int ungapped = 0;

  Resetint32(ungappedOffsets);

  for(c = consensus;
      *c != '\0';
      c++){

    Appendint32(ungappedOffsets, &ungapped);

    if(*c != '-')
      ungapped++;
  }
  Appendint32(ungappedOffsets, &ungapped);

  return;
}

/****************************************************************************************************/
/*        MultiAlignStore                                                                           */
/****************************************************************************************************/


MultiAlignStoreT *
CreateMultiAlignStoreT(int32 size)
{
  MultiAlignStoreT *mas = (MultiAlignStoreT *)safe_malloc(sizeof(MultiAlignStoreT));
  mas->multiAligns = CreateVA_PtrT(size);
  return mas;
}

// Empty the MultiAlignStore
size_t 
ClearMultiAlignStoreT(MultiAlignStoreT *multiAlignStore)
{
  int i;
  size_t redeemed = 0;
  int32 numMultiAligns = GetNumMultiAlignTs(multiAlignStore->multiAligns);
  MultiAlignT **map = (MultiAlignT **) GetPtrT(multiAlignStore->multiAligns, 0); 
  for(i = 0; i < numMultiAligns; i++, map++){
    void *dummy = NULL;
    MultiAlignT *ma = *map;
    if(ma){
      if( ma->refCnt <= 1)
	redeemed += GetMemorySize(ma);
      DeleteMultiAlignT(ma);
      SetPtrT(multiAlignStore->multiAligns,i,&dummy);
    }
  }
  
  return redeemed;
}

// Delete the multiAlignStore and all of its referenced data
void 
DeleteMultiAlignStoreT(MultiAlignStoreT *multiAlignStore)
{
  int i;
  int32 numMultiAligns = GetNumMultiAlignTs(multiAlignStore->multiAligns);
  for(i = 0; i < numMultiAligns; i++){
    MultiAlignT *ma = (MultiAlignT *) *GetPtrT(multiAlignStore->multiAligns, i);
    if(ma){
      DeleteMultiAlignT(ma);
    }
  }
  
  DeleteVA_PtrT(multiAlignStore->multiAligns);
  safe_free(multiAlignStore);
}

// Persistence
void 
SaveMultiAlignStoreTToStream(MultiAlignStoreT *mas, FILE *stream, 
  int withReferences)
{
  int i;
  int status;
  int32 size = GetNumPtrTs(mas->multiAligns);
  status = safeWrite(stream, &size, sizeof(int32));
  assert(status == FALSE);
  for(i = 0; i < size; i++){
    MultiAlignT *ma = (MultiAlignT *) *GetPtrT(mas->multiAligns, i);
    //    fprintf(stderr,"* i = %d ma = 0x%x\n", i, ma);
    if(!ma || withReferences || GetReferenceCountMultiAlignT(ma) == 1){
      //      fprintf(stderr,"Saving ma %d \n",i);
      SaveMultiAlignTToStream(ma,stream);
    }else{
      //      fprintf(stderr,"Saving ma %d as reference (cnt = %d) to ma %d\n", 
      //	      i,GetReferenceCountMultiAlignT(ma),ma->id );
      SaveReferenceMultiAlignTToStream(ma,stream);
    }
  }
}


MultiAlignStoreT *
LoadMultiAlignStoreTFromStream(FILE *stream)
{
  MultiAlignStoreT *mas = NULL;
  int i;
  int status;
  int32 size;
  int32 reference;
  status = safeRead(stream, &size, sizeof(int32));
  mas = CreateMultiAlignStoreT(size);
  assert(status == FALSE);
  for(i = 0; i < size; i++){
    MultiAlignT *ma =(MultiAlignT *)
         LoadMultiAlignTFromStream(stream, &reference);

    //    fprintf(stderr,"* i = %d ma = 0x%x\n", i, ma);
    //    if(ma)assert(i == ma->id);
    assert(reference == NULLINDEX);
    SetPtrT(mas->multiAligns, i, (const void *) &ma);
    if(ma)
      AddReferenceMultiAlignT(ma);
  }
  return mas;
}

MultiAlignStoreT *
LoadMultiAlignStoreTFromStreamWithReferences(FILE *stream, 
  MultiAlignStoreT *original)
{
  MultiAlignStoreT *mas;
  int i;
  int status;
  int32 size;
  int32 reference;
  status = safeRead(stream, &size, sizeof(int32));
  mas = CreateMultiAlignStoreT(size);
  assert(status == FALSE);
  for(i = 0; i < size; i++){
    MultiAlignT *ma = (MultiAlignT *) 
        LoadMultiAlignTFromStream(stream, &reference);

    if(reference != NULLINDEX){
      //      fprintf(stderr,"* MultiAlign %d is reference to multiAlign %d in original store\n",
      //	      i, reference);
      ma = (MultiAlignT *) 
           GetMultiAlignInStore(original, reference);
    }
    //    fprintf(stderr,"* i = %d ma = 0x%x ref:%d\n", i, ma, reference);
    if(ma){
      AddReferenceMultiAlignT(ma);
    }
    SetPtrT(mas->multiAligns, i, (const void *)&ma);
  }
  return mas;
}


// Accessors
void 
SetMultiAlignInStore(MultiAlignStoreT *mas, int index, MultiAlignT *ma)
{
#ifdef DEBUG_MULTIALIGN
  fprintf(stderr,"* Inserting multiAlign with length %ld and %ld frags at index %d\n",
	  GetNumchars(ma->consensus), GetNumIntMultiPoss(ma->f_list), index);
#endif
    SetPtrT(mas->multiAligns, index, (const void *)&ma);
    if(ma)
      AddReferenceMultiAlignT(ma);
}

MultiAlignT *
GetMultiAlignInStore(MultiAlignStoreT *mas,  int index)
{
  MultiAlignT **ptrRetValue = (MultiAlignT **)GetPtrT(mas->multiAligns, index);
  if(ptrRetValue)
    return *ptrRetValue;
  return (MultiAlignT *)NULL;
}

size_t 
RemoveMultiAlignFromStore(MultiAlignStoreT *mas, int index)
{
  MultiAlignT *ma = GetMultiAlignInStore(mas, index);
  size_t redeemed = 0;
  if(!ma)
    return redeemed;

  if(ma->refCnt <= 1){
      redeemed += GetMemorySize(ma);
  }
  DeleteMultiAlignT(ma);
  SetMultiAlignInStore(mas, index, NULL);
  return redeemed;
}


int64 
StatsMultiAlignStore(MultiAlignStoreT *maStore, FILE *fout, int owner)
{
  size_t totalMemorySize = 0;
  size_t maSize = 0;
  int32 numMultiAligns = GetNumMultiAlignsInStore(maStore);
  int i;
  totalMemorySize += ReportMemorySize_VA(maStore->multiAligns,"MultiAligns",fout);
  for(i = 0; i < numMultiAligns; i++){
    MultiAlignT *ma = *(MultiAlignT **)GetPtrT(maStore->multiAligns, i);
    if(ma && (owner || ma->refCnt == 1))
      maSize += GetMemorySize(ma);
  }
  totalMemorySize += maSize;
  fprintf(fout,"* MultiAlignStore has %d multiAligns occupying " F_SIZE_T " total size\n",
	  numMultiAligns, totalMemorySize);
  return totalMemorySize;
}


// Clone
MultiAlignStoreT *
CloneMultiAlignStoreT(MultiAlignStoreT *original)
{
  MultiAlignStoreT *mas;
  int i;
  int32 size = GetNumMultiAlignsInStore(original);
  mas = CreateMultiAlignStoreT(size);
  for(i = 0; i < size; i++){
    MultiAlignT *ma = GetMultiAlignInStore(original, i);
    AddReferenceMultiAlignT(ma);
    SetPtrT(mas->multiAligns, i, (const void *)&ma);
  }
  return mas;
}

int 
GetCoverageInMultiAlignT(MultiAlignT *ma, SeqInterval range,
                VA_TYPE(int) *covinput, int includeExternal) 
{

/* Returns 1 if completely covered by fragment sequence, else returns 0 */
/*    if ( ! includeExternal ), then only CeleraRead data is counted */
/*    otherwise, all fragments are counted */
    int32 left;
    int32 right;
    int *cov;
    int rc;
    IntMultiPos *reads=GetIntMultiPos(ma->f_list,0);
    int i; // tracks reads
    int j; // tracks columns
    int num_reads = GetNumIntMultiPoss(ma->f_list);
    int range_width = range.end-range.bgn;
    
    cov = (int *) safe_malloc(range_width*sizeof(int));
    for (i=0;i<range_width;i++) { cov[i] = 0;}
   
    for (i=0;i<num_reads;i++) {
      left = (reads[i].position.bgn < reads[i].position.end)? reads[i].position.bgn : reads[i].position.end;
      right = (reads[i].position.bgn > reads[i].position.end)?reads[i].position.bgn:reads[i].position.end;
      if ( left > range.end ) break; // beyond end of range
      if ( left <= range.end && right >= range.bgn ) {
       j = ( left < range.bgn )?range.bgn:left;
       if ( j >= left && j<right) {
        while ( j<right && j <range.end) {
          if ( includeExternal ||
               ( reads[i].type == AS_READ ) ||
               ( reads[i].type == AS_EXTR ) ||
               ( reads[i].type == AS_TRNR ) ){
             cov[j-range.bgn]++; 
          }
          j++;
        }
       }
      }
    }
    ResetVA_int(covinput);
    rc = 1;
    for (j=0;j<range_width;j++) {
       if (cov[j] == 0) rc = 0;
       SetVA_int(covinput, j, &cov[j]);
    }
    safe_free(cov);
    return rc;
}

typedef enum {
    CNS_INSERT = (int) 'I',
    CNS_DELETE = (int) 'D',    
    CNS_SUBSTITUTE = (int) 'S'
}ErrorType;

typedef struct{
    int position;
    ErrorType type;
}ErrorStruct;
 
VA_DEF(ErrorStruct)

void 
CollectStats(MultiAlignT *ma,
             GateKeeperStore *frag_store, 
             FILE *column_stats, 
             FILE *frag_stats,
             uint32 clrrng_flag)
{
/*  
    Need to append to column_stats and frag_stats the following:
    To column_stats:
       Foreach column in multialignment, print contigID, column index, coverage, quality value

    To frag_stats:
       Foreach fragment in multialignment, print fragIID, fragUID,  clr_bgn, clr_end, errors (in apos,type pairs)
*/
    int32 readptr;
    int32 delptr;
    int32 left;
    int32 right;
    int32 flen;
    int32 ungapped=0;
    uint clrbgn;
    uint clrend;
    CDS_UID_t accession;
    IntMultiPos *reads=GetIntMultiPos(ma->f_list,0);
    int i; // tracks reads
    int j; // tracks columns
    int32 ma_len = GetMultiAlignLength(ma);
    int num_reads = GetNumIntMultiPoss(ma->f_list);
    int *column_cov;
    int *column_mm;
    char column_call;
    int num_errors;
    VA_TYPE(ErrorStruct) *errors; 
    ErrorStruct frag_error;
    char tmpseq[AS_READ_MAX_LEN+2];
    char tmpqv[AS_READ_MAX_LEN+2];
    char seqdata[AS_READ_MAX_LEN+2];
    char qvdata[AS_READ_MAX_LEN+2];
    ReadStructp rsp = new_ReadStruct();
    
    column_cov = (int *) safe_malloc(ma_len*sizeof(int));
    column_mm = (int *) safe_malloc(ma_len*sizeof(int));
    errors = CreateVA_ErrorStruct(250);

    assert(column_cov && column_mm);

    // special case for singletons
    if (num_reads == 1) {
       getFrag(frag_store,reads[0].ident,rsp,FRAG_S_ALL);
       getClearRegion_ReadStruct(rsp, &clrbgn,&clrend, clrrng_flag);
       fprintf(frag_stats,F_IID "  " F_UID " %c %d %d\n",
               reads[0].ident,accession,
               reads[0].type,(int) clrbgn,(int) clrend);
      flen = clrend - clrbgn;
      if(getSequence_ReadStruct(rsp, tmpseq, tmpqv, AS_READ_MAX_LEN+1) != 0)
        assert(0);
      // capture only the clear range for analysis
      // reverse complement if necessary:
      memcpy(seqdata, &tmpseq[clrbgn], (flen+1)*sizeof(char));
      memcpy(qvdata, &tmpqv[clrbgn], (flen+1)*sizeof(char));
      seqdata[flen] = '\0';
      qvdata[flen] = '\0';
       for (j=0;j<ma_len;j++) {
         fprintf(column_stats,"%d %d %d %d %d %c %d %d %d\n",ma->id,j,1,
              0,
              (reads[0].type != AS_READ &&
               reads[0].type != AS_EXTR &&
               reads[0].type != AS_TRNR)?1:0, 
             seqdata[j],qvdata[j] - '0',0,j);
       }
    } else {
        
     // initialize column coverage to zero
     for ( j=0;j<ma_len;j++) {
       column_cov[j]=0;
       column_mm[j]=0;
     }

    for(i=0;i<num_reads;i++) {
      left = (reads[i].position.bgn < reads[i].position.end)? reads[i].position.bgn : reads[i].position.end;
      right= (reads[i].position.bgn > reads[i].position.end)?reads[i].position.bgn:reads[i].position.end;
      getFrag(frag_store,reads[i].ident,rsp,FRAG_S_ALL);
      getClearRegion_ReadStruct(rsp, &clrbgn,&clrend, clrrng_flag);
      flen = clrend - clrbgn;
      assert(flen < AS_READ_MAX_LEN);
      assert(flen > 0);
      if(getSequence_ReadStruct(rsp, tmpseq, tmpqv, AS_READ_MAX_LEN+1) != 0)
        assert(0);
      // capture only the clear range for analysis
      // reverse complement if necessary:
      memcpy(seqdata, &tmpseq[clrbgn], (flen+1)*sizeof(char));
      memcpy(qvdata, &tmpqv[clrbgn], (flen+1)*sizeof(char));
      seqdata[flen] = '\0';
      qvdata[flen] = '\0';
      if (reads[i].position.bgn > reads[i].position.end) {
        SequenceComplement(seqdata,qvdata);
      }
      getAccID_ReadStruct(rsp, &accession);
      ResetErrorStruct(errors);
      
      readptr= 0;
      delptr = 0;
       for ( j=left;j<right;j++) {
         assert (j < ma_len );
         column_call=*Getchar(ma->consensus,j);
         if ( delptr < reads[i].delta_length ) {
             if ( readptr != *(reads[i].delta + delptr)) {
               // non gap coverage for this fragment in this column
               // compare base at readptr[i] in frag sequence to column_call
               if ( seqdata[readptr] != column_call ) {
                  // record the error
                  frag_error.position = readptr;
                  if (column_call == '-') { // insertion
                     frag_error.type = CNS_INSERT; 
                  } else {
                     frag_error.type = CNS_SUBSTITUTE; 
                  } 
                  AppendErrorStruct(errors,&frag_error);
                  column_mm[j]+=1;
               }
               readptr++;  column_cov[j]+=1;
             } else {
               // gap for this fragment in this column
               if ( '-' != column_call ) {
                  // record the error
                  frag_error.position = readptr;
                  frag_error.type = CNS_DELETE;
                  AppendErrorStruct(errors,&frag_error);
                  column_mm[j]+=1;
                  column_cov[j]+=1; //adding this so that intra-fragment gaps count as coverage
               }
               delptr++;
             }
           } else {
             // non gap coverage for this fragment in this column
             // compare base at readptr[i] in frag sequence to column_call
             if ( seqdata[readptr] != column_call ) {
                  // record the error
                  frag_error.position = readptr;
                  if (column_call == '-') { // insertion
                     frag_error.type = CNS_INSERT; 
                  } else {
                     frag_error.type = CNS_SUBSTITUTE; 
                  } 
                  AppendErrorStruct(errors,&frag_error);
                  column_mm[j]+=1;
             }
             readptr++;  column_cov[j]+=1;
           }
       }
       fprintf(frag_stats,F_IID " " F_UID " %c %d %d",
               reads[i].ident,accession,
               reads[i].type,(int) clrbgn,(int) clrend);
       num_errors = 0;
       if ( GetNumErrorStructs(errors) > 75 ) {
          fprintf(frag_stats," misaligned fragment with %d mismatches\n",
                  (int) GetNumErrorStructs(errors));
       } else {
       while ( GetErrorStruct(errors,num_errors) ) {
          frag_error = *GetErrorStruct(errors,num_errors);
          fprintf(frag_stats," %d %c",frag_error.position, frag_error.type);
          num_errors++;
       }
       fprintf(frag_stats,"\n"); 
       }
    }
    for (j=0;j<ma_len;j++) {
      fprintf(column_stats,"%d %d %d %c %d %d %d\n",ma->id,j,
              column_cov[j],
              *Getchar(ma->consensus,j),
              (int) *Getchar(ma->quality,j) - '0',
              column_mm[j],ungapped);
      if (*Getchar(ma->consensus,j) != '-') ungapped++;
    }
    }
    fflush(column_stats);
    fflush(frag_stats);
    safe_free(column_cov);
    safe_free(column_mm);
    DeleteVA_ErrorStruct(errors);
    delete_ReadStruct(rsp);
}


// Format the fragment type for display.

static char 
getFragTypeDisplay (FragType fragType) 
{
    char dispType;
    switch (fragType) {
      // AS_READ is normally 'R'
      // But for reports, we were asked to show ' '
      // for these perponderant normal reads.
    case AS_READ: dispType  = ' ';  // Celera Read
        break;
    case AS_EXTR: dispType  = 'X';  //External WGS read
        break;
    case AS_TRNR: dispType  = 'T';  //Transposon library read
        break;
    default: dispType = '?';
        break;
    }
    return dispType;
}



// Print a character representation of alignment.

int 
PrintMultiAlignT(FILE *out,
	         MultiAlignT *ma,
	         GateKeeperStore *frag_store, 
	         tFragStorePartition *pfrag_store,
	         int show_qv, 
	         int dots,
                 uint32 clrrng_flag) 
{
  char frgTypeDisplay;
  FragType frgTypeData;
  int depth;
  int rc,i;
  int window,length;
  char **multia=NULL; 
  int **idarray;
  int **oriarray;
  char *consensus = Getchar(ma->consensus,0);
  char *quality = Getchar(ma->quality,0);
  char *nonblank;
  static char *sep0="___________________________________________________________________________________________________________________________________";
  static char *sep1="         |         |         |         |         |         |         |         |         |         |";
  static  ReadStructp rsp=NULL;
  int partitioned=0;
  length = strlen(consensus);
  if (rsp==NULL) {
     rsp  = new_ReadStruct();
  }
  if ( frag_store == NULL ) {
   partitioned = 1;
  }
   
  rc = MultiAlignT2Array(ma, frag_store, pfrag_store,
                         &depth, &multia, &idarray, &oriarray, clrrng_flag);
  if (rc) {
    fprintf(out,"<<< begin Contig %d >>>",ma->id);;
       int ungapped=1;
       int tick=1;
       for (window=0;window<length;) {
          CDS_UID_t uid=0;
          int row_id=0;
          int rowind=0;
          int orient=0;
          int rowlen=0;
          int labelchars=0;
          char *rowchar=consensus+window;
          labelchars = fprintf(out,"\n\n                                        <<<  Contig %d, gapped length: %d  >>>\n",ma->id,length);
          fprintf(out,"%d gapped\n",window+1);
          fprintf(out,"%-100.100s\n",sep1);
          fprintf(out,"%d ungapped\n",ungapped+tick-1);
          rowlen = (window+100 < length)?100:length-window;
          for (rowind=0;rowind<rowlen;rowind++,rowchar++){
             if ( tick==10 ) {
               ungapped+=10;
               tick=0;
             }
             if ( tick==0 && *rowchar!='-') {
                 fprintf(out,"u");
             } else {
                 fprintf(out," ");
             }
             if (*rowchar!='-') {
               tick++;
             } 
          }     
          fprintf(out,"\n");
          fprintf(out,"%-100.100s  cns  (uid,iid) type\n",consensus+window);
          if (show_qv) fprintf(out,"%-100.100s  qlt\n",quality+window);
          fprintf(out,"%-130.130s\n",sep0);
          for (i=0;i<depth;i++) {
           if (multia[2*i] == NULL) continue;
           nonblank = strpbrk(multia[2*i]+window,"ACGT");
           if ( nonblank == NULL || nonblank-(multia[2*i]+window) > 100 ) continue;
           {
             int j;
             for (j=0;j<100;j++) {
                if ( window+j> length) break;
                if ( dots && *(multia[2*i]+window+j) == *(consensus+window+j) ) {
                   *(multia[2*i]+window+j) = '.';
                   *(multia[2*i+1]+window+j) = ' ';
                } else {
                   *(multia[2*i]+window+j) = tolower(*(multia[2*i]+window+j));
                }
             }
           }
           {int last = (window+99< length-1)?window+99:length-1;
           if ( *(idarray[i]+last) == 0 ) {
               row_id = *(idarray[i]+window);
               orient = *(oriarray[i]+window);
           } else {
               row_id = *(idarray[i]+last);
               orient = *(oriarray[i]+last);
           }
           }
           // Look up UID for row_id
           if ( row_id > 0 ) {
             if ( partitioned ) {
                  getFragStorePartition(pfrag_store,
					row_id,
					FRAG_S_INF,
                                        rsp);
             } else {
                  getFrag(frag_store,
			       row_id,
                               rsp,
                               FRAG_S_INF);
             }
             frgTypeDisplay = ' ';
             getAccID_ReadStruct(rsp, &uid); 
             //getReadType_ReadStruct(rsp, &frgTypeData);
             frgTypeData = AS_READ;

	     frgTypeDisplay = getFragTypeDisplay(frgTypeData);
             ///if ( type == AS_READ) type = ' ';

             fprintf(out,
		     "%-100.100s   %c   (" F_UID ",%d) %c\n",
		     multia[2*i]+window,
                     (orient>0)?'>':'<',
		     uid,
		     row_id,
		     frgTypeDisplay);
	             ///type);
             if (show_qv) {
               fprintf(out,
		       "%-100.100s   %c   (" F_UID ",%d) %c\n",
		       multia[2*i+1]+window,
                       (orient>0)?'>':'<',
		       uid,
		       row_id,
		       frgTypeDisplay);
	               ///type);
             }
           }
          }
          window+=100;
       }
       fprintf(out,"\n<<< end Contig %d >>>\n",ma->id);
  } else {
       fprintf(stderr,"Error returned from MultiAlignT2Array.\n");
  }
  if (multia) {
       for (i=0;i<2*depth;i++) {
         safe_free((char *)multia[i]);
       }
       safe_free(multia);
       for (i=0;i<depth;i++) {
         safe_free((int *)idarray[i]);
         safe_free((int *)oriarray[i]);
       }
       safe_free(idarray);
       safe_free(oriarray);
  }
  return 1;
}

int 
PrintMultiAlignTSNPs(
		     FILE *out,
		     MultiAlignT *ma,
		     GateKeeperStore *frag_store, 
		     tFragStorePartition *pfrag_store,
		     int show_qv, 
		     int dots,uint32 clrrng_flag) 
{
  int depth;
  int rc,i,j;
  int length;
  char **multia; 
  int **idarray;
  int **oriarray;
  char *consensus = Getchar(ma->consensus,0);
  static char *sep0="_______________________________________________________________________________________________________________";
  static char *sep1="         |         |         |         |         |         |         |         |         |         |";
  int partitioned=0;

  
  length = strlen(consensus);
  if ( frag_store == NULL ) {
   partitioned = 1;
  }
   
  rc = MultiAlignT2Array(ma, frag_store, pfrag_store,
                         &depth, &multia, &idarray,&oriarray,clrrng_flag);
  if (rc) {
       BaseCount profile;
       // First, run through all colunms and use the "oriarray" to store whether a column has
       // a confirmed mismatch
       for (j=0;j<length;j++) {
         char cns=consensus[j];
         char mm;
         ResetBaseCount(&profile);
         for (i=0;i<depth;i++) {
           if (multia[2*i] == NULL) continue;
           {
             char b=*(multia[2*i]+j);
             if ( b != ' ') {
               IncBaseCount(&profile,b);
               if ( dots && b == *(consensus+j) ) {
                   *(multia[2*i]+j) = '.';
                   *(multia[2*i+1]+j) = ' ';
               } else {
                   *(multia[2*i]+j) = tolower(*(multia[2*i]+j));
               }
             }
           }
         }
         mm = GetConfMM(&profile,BaseToInt(cns));
         if (mm == cns) { 
           *(oriarray[0]+j)=0; 
         } else {
           *(oriarray[0]+j)=1; 
         } 
      }
     // at this point, all of the columns with confirmed mismatches are marked (via oriarray[0])
     // now, go through these and output them in compressed form, say 100/line
     // outputting also the coordinate of the first in the row
     {
       int prev_ids[depth];
       int win_start=0;
       int lcase[depth];
       while ( *(oriarray[0]+win_start) == 0 ) win_start++;
       for (i=0;i<depth;i++) { prev_ids[i] = -1; lcase[i]=1;}
       while ( win_start < length ) {
          int snp_cnt=0;
          int j=win_start;
          int tick_cnt=1;
          fprintf(out,"MultiAlign offset: %d\n",win_start);
          while (j<length && snp_cnt < 100) { // output the consensus line
                if ( *(oriarray[0]+j) == 1 ) { 
                  // this is a SNP column 
                  char oc= *(consensus+j);
                  fprintf(out,"%c %d",oc,j);
                  if ( tick_cnt==10 ) {
                     fprintf(out," -- %d",snp_cnt+1);
                     tick_cnt=1;
                  }  else {
                     tick_cnt++;
                  }
                  fprintf(out,"\n");
                  snp_cnt++;
                }
                j++;
          }
          snp_cnt=0;
          j=win_start;
          fprintf(out,"\n%-110.110s\n",sep1);
          while (j<length && snp_cnt < 100) { // output the consensus line
                if ( *(oriarray[0]+j) == 1 ) { 
                  // this is a SNP column 
                  char oc= *(consensus+j);
                  fprintf(out,"%c",oc);
                  snp_cnt++;
                }
                j++;
          }
          while ( snp_cnt < 100 ) {
                  fprintf(out," ");
                  snp_cnt++;
          }
          fprintf(out,"  cns\n");
          fprintf(out,"%-110.110s\n",sep0);
          for (i=0;i<depth;i++)  {
             snp_cnt=0;
             j=win_start;
             while (j<length && snp_cnt < 100) {
                if ( *(oriarray[0]+j) == 1 ) { 
                  // this is a SNP column 
                  char oc= *(multia[2*i]+j);
                  if ( oc != ' ') {
                   if  ( *(idarray[i]+j) != prev_ids[i] ) {
                     // new fragment, change case
                     if (lcase[i]==0) lcase[i]++;
                       else lcase[i]--;
                     prev_ids[i]= *(idarray[i]+j);
                   }
                   if ( ! lcase[i] && oc == '.' ) oc=',';
                   if ( ! lcase[i] && oc == '-' ) oc='~';
                   else if ( ! lcase[i] ) oc = toupper(oc);
                  }
                  fprintf(out,"%c",oc);
                  snp_cnt++;
                }
                j++;
             }
             fprintf(out,"\n");
          }
          win_start=j;
          while ( *(oriarray[0]+win_start) == 0 ) win_start++;
       }
    }
  }
  if (multia) {
       for (i=0;i<2*depth;i++) {
         safe_free((char *)multia[i]);
       }
       safe_free(multia);
       for (i=0;i<depth;i++) {
         safe_free((int *)idarray[i]);
         safe_free((int *)oriarray[i]);
       }
       safe_free(idarray);
       safe_free(oriarray);
  }
  return 1;
}

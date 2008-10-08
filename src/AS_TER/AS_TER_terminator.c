
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

const char *mainid = "$Id: AS_TER_terminator.c,v 1.29 2008-10-08 22:03:00 brianwalenz Exp $";

//  Assembly terminator module. It is the backend of the assembly
//  pipeline and replaces internal accession numbers by external
//  accession numbers.

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <unistd.h>
#include  <assert.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "SYS_UIDclient.h"


//  Info loaded from gatekeeper.
typedef struct {
  int32       loaded:1;
  int32       deleted:1;
  int32       clearBeg:15;
  int32       clearEnd:15;
  AS_UID      uid;
} fragInfo_t;

fragInfo_t    *fragInfo;

VA_DEF(AS_UID);
VA_DEF(AS_IID);

VA_TYPE(AS_IID) *ICMinISF1;
VA_TYPE(AS_IID) *ICMinISF2;

VA_TYPE(AS_UID) *IUMmap;
VA_TYPE(AS_UID) *ICMmap;
VA_TYPE(AS_UID) *ISFmap;
VA_TYPE(AS_UID) *FRGmap;
VA_TYPE(AS_UID) *DSTmap;

UIDserver       *uids;

void
DumpIID2UIDmap(VA_TYPE(AS_UID) *map,
               char *label,
               FILE *F) {
  int      i;

  for(i=0; i<GetNumAS_UIDs(map); i++){
    AS_UID  di = *GetAS_UID(map,i);
    if (AS_UID_isDefined(di))
      fprintf(F,"%s\t%d\t%s\n", label, i, AS_UID_toString(di));
  }
}

AS_UID
lookupUID(VA_TYPE(AS_UID) *map, AS_IID    iid){
  AS_UID    *ret = GetAS_UID(map, iid);
  if ((ret == NULL) || (AS_UID_isDefined(*ret) == FALSE))
    return(AS_UID_undefined());
  return(*ret);
}

int
existsUID(VA_TYPE(AS_UID) *map, AS_IID    iid){
  AS_UID    *ret = GetAS_UID(map, iid);
  if ((ret == NULL) || (AS_UID_isDefined(*ret) == FALSE))
    return(0);
  return(1);
}

AS_IID
lookupIID(VA_TYPE(AS_IID) *map, AS_IID    iid){
  AS_IID    *ret = GetAS_IID(map, iid);
  if ((ret == NULL) || (AS_IID_isDefined(*ret) == FALSE))
    return(0);
  return(*ret);
}


void
convertIAF(GenericMesg *pmesg,
           FILE        *fileOutput) {
  IntAugFragMesg *iafMesg = (IntAugFragMesg *)pmesg->m;
  AugFragMesg     afgMesg;

  // This assertion was added by Jason, Oct 2001, while adding
  // modified clear range fields to the frag store.  If we reach
  // this point in the code, then some program upstream of
  // Terminator has written a modified clear range in the IAF
  // message.  That violates my assumption that all clear range
  // modifications get written to the frag store.
  //
  assert(iafMesg->clear_rng.bgn == -1);
  assert(iafMesg->clear_rng.end == -1);

  if (existsUID(FRGmap, iafMesg->iaccession)) {
    fprintf(stderr, "IAF: Spotted FRG internal ID "F_IID" second time\n", iafMesg->iaccession);
    exit(1);
  }

  afgMesg.eaccession    = fragInfo[iafMesg->iaccession].uid;
  afgMesg.iaccession    = iafMesg->iaccession;
  afgMesg.mate_status   = iafMesg->mate_status;
  afgMesg.chaff         = iafMesg->chaff;
  afgMesg.clear_rng.bgn = fragInfo[iafMesg->iaccession].clearBeg;
  afgMesg.clear_rng.end = fragInfo[iafMesg->iaccession].clearEnd;

  SetAS_UID(FRGmap, iafMesg->iaccession, &afgMesg.eaccession);

  pmesg->m = &afgMesg;
  pmesg->t = MESG_AFG;

  WriteProtoMesg_AS(fileOutput, pmesg);
}


void
convertIAM(GenericMesg *pmesg,
           FILE        *fileOutput) {
  IntAugMatePairMesg *iamMesg = (IntAugMatePairMesg *)pmesg->m;
  AugMatePairMesg     ampMesg;

  if (!existsUID(FRGmap, iamMesg->fragment1)) {
    fprintf(stderr, "IAM:  reference before definition error for fragment ID "F_IID"\n", iamMesg->fragment1);
    exit(1);
  }
  if (!existsUID(FRGmap, iamMesg->fragment2)) {
    fprintf(stderr, "IAM:  reference before definition error for fragment ID "F_IID"\n", iamMesg->fragment2);
    exit(1);
  }

  ampMesg.fragment1   = fragInfo[iamMesg->fragment1].uid;
  ampMesg.fragment2   = fragInfo[iamMesg->fragment2].uid;
  ampMesg.mate_status = iamMesg->mate_status;

  pmesg->m = &ampMesg;
  pmesg->t = MESG_AMP;

  WriteProtoMesg_AS(fileOutput, pmesg);
}


void
convertIUM(GenericMesg *pmesg,
           FILE        *fileOutput) {

  IntUnitigMesg  *iumMesg = (IntUnitigMesg *)pmesg->m;
  SnapUnitigMesg  utgMesg;
  int             i;

  if (existsUID(IUMmap, iumMesg->iaccession)) {
    fprintf(stderr, "IUM: Spotted IUM internal ID "F_IID" second time\n", iumMesg->iaccession);
    exit(1);
  }

  utgMesg.eaccession     = AS_UID_fromInteger(getUID(uids));
  utgMesg.iaccession     = iumMesg->iaccession;
#ifdef AS_ENABLE_SOURCE
  utgMesg.source         = strdup(iumMesg->source);
#endif
  utgMesg.coverage_stat  = iumMesg->coverage_stat;
  utgMesg.status         = iumMesg->status;
  utgMesg.length         = iumMesg->length;
  utgMesg.consensus      = strdup(iumMesg->consensus);
  utgMesg.quality        = strdup(iumMesg->quality);
  utgMesg.forced         = iumMesg->forced;
  utgMesg.num_frags      = iumMesg->num_frags;
  utgMesg.num_vars       = 0;
  utgMesg.f_list         = NULL;
  utgMesg.v_list         = NULL;

  SetAS_UID(IUMmap, iumMesg->iaccession, &utgMesg.eaccession);

  if (iumMesg->num_frags > 0) {
    utgMesg.f_list = (SnapMultiPos*)safe_malloc(iumMesg->num_frags * sizeof(SnapMultiPos));

    for(i=0; i<iumMesg->num_frags; i++){
      if (existsUID(FRGmap, iumMesg->f_list[i].ident) == 0) {
        fprintf(stderr,"IUM: Reference before definition for fragment ID "F_IID"\n", iumMesg->f_list[i].ident);
        exit(1);
      }

      utgMesg.f_list[i].type         = iumMesg->f_list[i].type;
#ifdef AS_ENABLE_SOURCE
      utgMesg.f_list[i].source       = NULL;
#endif
      utgMesg.f_list[i].eident       = lookupUID(FRGmap, iumMesg->f_list[i].ident);
      utgMesg.f_list[i].delta_length = iumMesg->f_list[i].delta_length;
      utgMesg.f_list[i].position     = iumMesg->f_list[i].position;
      utgMesg.f_list[i].delta        = NULL;

      if (utgMesg.f_list[i].delta_length > 0)
	utgMesg.f_list[i].delta = iumMesg->f_list[i].delta;
    }
  }

  pmesg->m = &utgMesg;
  pmesg->t = MESG_UTG;

  WriteProtoMesg_AS(fileOutput,pmesg);

  safe_free(utgMesg.f_list);
}


void
convertIUL(GenericMesg *pmesg,
           FILE        *fileOutput) {

  IntUnitigLinkMesg  *iulMesg = (IntUnitigLinkMesg*) pmesg->m;
  SnapUnitigLinkMesg  ulkMesg;
  int32               jumplistLength;
  int                 i;

  if (existsUID(IUMmap, iulMesg->unitig1) == 0) {
    fprintf(stderr,"IUL: reference before definition error for unitig ID "F_IID"\n",iulMesg->unitig1);
    exit(1);
  }
  if (existsUID(IUMmap, iulMesg->unitig2) == 0) {
    fprintf(stderr,"IUL: reference before definition error for unitig ID "F_IID"\n",iulMesg->unitig2);
    exit(1);
  }

  ulkMesg.eunitig1            = lookupUID(IUMmap, iulMesg->unitig1);
  ulkMesg.eunitig2            = lookupUID(IUMmap, iulMesg->unitig2);

  ulkMesg.orientation         = iulMesg->orientation;
  ulkMesg.overlap_type        = iulMesg->overlap_type;
  ulkMesg.is_possible_chimera = iulMesg->is_possible_chimera;
  ulkMesg.includes_guide      = iulMesg->includes_guide;
  ulkMesg.mean_distance       = iulMesg->mean_distance;
  ulkMesg.std_deviation       = iulMesg->std_deviation;
  ulkMesg.num_contributing    = iulMesg->num_contributing;
  ulkMesg.status              = iulMesg->status;

  jumplistLength = ulkMesg.num_contributing;

  //  a case distinction to find out the number of elements in jump_list
  if (iulMesg->overlap_type != AS_NO_OVERLAP)
    jumplistLength--;

  ulkMesg.jump_list = NULL;
  if (jumplistLength > 0)
    ulkMesg.jump_list = (SnapMate_Pairs *)safe_malloc(jumplistLength * sizeof(SnapMate_Pairs));

  //  traverse the jump list and replace in the MatePairs the internal
  //  IDs by the external IDs

  for(i=0; i<jumplistLength; i++) {

    if (existsUID(FRGmap, iulMesg->jump_list[i].in1) == 0) {
      fprintf(stderr,"IUL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", iulMesg->jump_list[i].in1);
      exit(1);
    }
    if (existsUID(FRGmap, iulMesg->jump_list[i].in2) == 0) {
      fprintf(stderr,"IUL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", iulMesg->jump_list[i].in2);
      exit(1);
    }

    ulkMesg.jump_list[i].in1 = lookupUID(FRGmap, iulMesg->jump_list[i].in1);
    ulkMesg.jump_list[i].in2 = lookupUID(FRGmap, iulMesg->jump_list[i].in2);

    ulkMesg.jump_list[i].type = iulMesg->jump_list[i].type;
  }

  pmesg->m = &ulkMesg;
  pmesg->t = MESG_ULK;

  WriteProtoMesg_AS(fileOutput,pmesg);

  safe_free(ulkMesg.jump_list);
}


void
convertICM(GenericMesg *pmesg,
           FILE        *fileOutput) {

  IntConConMesg   *icmMesg = (IntConConMesg*) pmesg->m;
  SnapConConMesg   ccoMesg;
  int              i;

  if (existsUID(ICMmap, icmMesg->iaccession)) {
    fprintf(stderr, "ICM: Spotted ICM internal ID "F_IID" second time\n", icmMesg->iaccession);
    exit(1);
  }

  ccoMesg.eaccession = AS_UID_fromInteger(getUID(uids));
  ccoMesg.iaccession = icmMesg->iaccession;
  ccoMesg.placed     = icmMesg->placed;
  ccoMesg.length     = icmMesg->length;
  ccoMesg.consensus  = icmMesg->consensus;
  ccoMesg.quality    = icmMesg->quality;
  ccoMesg.forced     = icmMesg->forced;
  ccoMesg.num_pieces = icmMesg->num_pieces;
  ccoMesg.num_unitigs= icmMesg->num_unitigs;
  ccoMesg.num_vars   = icmMesg->num_vars;
  ccoMesg.pieces     = NULL;
  ccoMesg.vars       = NULL;
  ccoMesg.unitigs    = NULL;

  SetAS_UID(ICMmap, icmMesg->iaccession, &ccoMesg.eaccession);

  if (icmMesg->num_vars > 0) {
    ccoMesg.vars = (IntMultiVar*) safe_malloc(icmMesg->num_vars * sizeof(IntMultiVar));

    for(i=0; i<icmMesg->num_vars; i++) {
      ccoMesg.vars[i].position         = icmMesg->v_list[i].position;
      ccoMesg.vars[i].num_reads        = icmMesg->v_list[i].num_reads;
      ccoMesg.vars[i].num_conf_alleles = icmMesg->v_list[i].num_conf_alleles;
      ccoMesg.vars[i].min_anchor_size  = icmMesg->v_list[i].min_anchor_size;
      ccoMesg.vars[i].var_length       = icmMesg->v_list[i].var_length ;
      ccoMesg.vars[i].curr_var_id      = icmMesg->v_list[i].curr_var_id;
      ccoMesg.vars[i].phased_var_id    = icmMesg->v_list[i].phased_var_id;
      ccoMesg.vars[i].nr_conf_alleles  = icmMesg->v_list[i].nr_conf_alleles;
      ccoMesg.vars[i].weights          = icmMesg->v_list[i].weights;
      ccoMesg.vars[i].var_seq          = icmMesg->v_list[i].var_seq;
      ccoMesg.vars[i].conf_read_iids   = icmMesg->v_list[i].conf_read_iids;
    }
  }

  if( ccoMesg.num_pieces > 0 ){
    ccoMesg.pieces = (SnapMultiPos*) safe_malloc(icmMesg->num_pieces*sizeof(SnapMultiPos));

    for(i=0; i<icmMesg->num_pieces; i++) {
      ccoMesg.pieces[i].type = icmMesg->pieces[i].type;
#ifdef AS_ENABLE_SOURCE
      ccoMesg.pieces[i].source = NULL;
#endif

      if (existsUID(FRGmap, icmMesg->pieces[i].ident) == 0) {
	fprintf(stderr,"ICM: Reference before definition for fragment ID "F_IID"\n", icmMesg->pieces[i].ident);
	exit(1);
      }
      ccoMesg.pieces[i].eident       = lookupUID(FRGmap, icmMesg->pieces[i].ident);
      ccoMesg.pieces[i].delta_length = icmMesg->pieces[i].delta_length;
      ccoMesg.pieces[i].position     = icmMesg->pieces[i].position;

      ccoMesg.pieces[i].delta = NULL;

      if (ccoMesg.pieces[i].delta_length > 0)
	ccoMesg.pieces[i].delta = icmMesg->pieces[i].delta;
    }
  }

  if( ccoMesg.num_unitigs > 0 ){
    ccoMesg.unitigs = (UnitigPos*) safe_malloc(icmMesg->num_unitigs*sizeof(UnitigPos));

    for(i=0; i<icmMesg->num_unitigs; i++){
      ccoMesg.unitigs[i].type  = icmMesg->unitigs[i].type;

      if (existsUID(IUMmap, icmMesg->unitigs[i].ident) == 0) {
	fprintf(stderr,"ICM: Reference before definition for unitig ID "F_IID"\n", icmMesg->unitigs[i].ident);
	exit(1);
      }
      ccoMesg.unitigs[i].eident       = lookupUID(IUMmap, icmMesg->unitigs[i].ident);
      ccoMesg.unitigs[i].position     = icmMesg->unitigs[i].position;
      ccoMesg.unitigs[i].delta        = icmMesg->unitigs[i].delta;
      ccoMesg.unitigs[i].delta_length = icmMesg->unitigs[i].delta_length;
    }
  }

  pmesg->m = &ccoMesg;
  pmesg->t = MESG_CCO;

  WriteProtoMesg_AS(fileOutput,pmesg);

  safe_free(ccoMesg.vars);
  safe_free(ccoMesg.pieces);
  safe_free(ccoMesg.unitigs);
}


void
convertICL(GenericMesg *pmesg,
           FILE        *fileOutput) {

  IntContigLinkMesg  *iclMesg = (IntContigLinkMesg*) pmesg->m;
  SnapContigLinkMesg  clkMesg;
  int32               jumplistLength;
  int                 i;

  if (existsUID(ICMmap, iclMesg->contig1) == 0) {
    fprintf(stderr,"ICL: reference before definition error for contig ID "F_IID"\n",iclMesg->contig1);
    exit(1);
  }
  if (existsUID(ICMmap, iclMesg->contig2) == 0) {
    fprintf(stderr,"ICL: reference before definition error for contig ID "F_IID"\n",iclMesg->contig2);
    exit(1);
  }

  clkMesg.econtig1            = lookupUID(ICMmap, iclMesg->contig1);
  clkMesg.econtig2            = lookupUID(ICMmap, iclMesg->contig2);

  clkMesg.orientation         = iclMesg->orientation;
  clkMesg.overlap_type        = iclMesg->overlap_type;
  clkMesg.is_possible_chimera = iclMesg->is_possible_chimera;
  clkMesg.includes_guide      = iclMesg->includes_guide;
  clkMesg.mean_distance       = iclMesg->mean_distance;
  clkMesg.std_deviation       = iclMesg->std_deviation;
  clkMesg.num_contributing    = iclMesg->num_contributing;
  clkMesg.status              = iclMesg->status;

  jumplistLength = clkMesg.num_contributing;

  //  a case distinction to determine the length of the jump_list
  if( iclMesg->overlap_type != AS_NO_OVERLAP )
    jumplistLength--;

  clkMesg.jump_list = NULL;
  if (jumplistLength > 0)
    clkMesg.jump_list = (SnapMate_Pairs*) safe_malloc(jumplistLength * sizeof(SnapMate_Pairs));

  //  traverse the jump_list and and replace the internal fragment IDs
  //  by external fragment IDs

  for (i=0; i<jumplistLength; i++) {
    if (existsUID(FRGmap, iclMesg->jump_list[i].in1) == 0) {
      fprintf(stderr,"ICL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", iclMesg->jump_list[i].in1);
      exit(1);
    }
    if (existsUID(FRGmap, iclMesg->jump_list[i].in2) == 0) {
      fprintf(stderr,"ICL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", iclMesg->jump_list[i].in2);
      exit(1);
    }

    clkMesg.jump_list[i].in1 = lookupUID(FRGmap, iclMesg->jump_list[i].in1);
    clkMesg.jump_list[i].in2 = lookupUID(FRGmap, iclMesg->jump_list[i].in2);

    clkMesg.jump_list[i].type = iclMesg->jump_list[i].type;
  }

  pmesg->m = &clkMesg;
  pmesg->t = MESG_CLK;

  WriteProtoMesg_AS(fileOutput,pmesg);

  safe_free(clkMesg.jump_list);
}


void
convertISL(GenericMesg *pmesg,
           FILE        *fileOutput) {

  InternalScaffoldLinkMesg *islMesg = (InternalScaffoldLinkMesg*) pmesg->m;
  SnapScaffoldLinkMesg      slkMesg;
  int32                     jumplistLength;
  int                       i;

  if (existsUID(ISFmap, islMesg->iscaffold1) == 0) {
    fprintf(stderr,"ISL: reference before definition error for scaffold ID "F_IID"\n", islMesg->iscaffold1);
    exit(1);
  }
  if (existsUID(ISFmap, islMesg->iscaffold2) == 0) {
    fprintf(stderr,"ISL: reference before definition error for scaffold ID "F_IID"\n", islMesg->iscaffold2);
    exit(1);
  }

  slkMesg.escaffold1 = lookupUID(ISFmap, islMesg->iscaffold1);
  slkMesg.escaffold2 = lookupUID(ISFmap, islMesg->iscaffold2);

  slkMesg.orientation    = islMesg->orientation;
  slkMesg.includes_guide = islMesg->includes_guide;
  slkMesg.mean_distance  = islMesg->mean_distance;
  slkMesg.std_deviation  = islMesg->std_deviation;
  slkMesg.num_contributing = islMesg->num_contributing;

  jumplistLength = slkMesg.num_contributing;

  slkMesg.jump_list = NULL;
  if( jumplistLength > 0 )
    slkMesg.jump_list = (SnapMate_Pairs*) safe_malloc(jumplistLength*sizeof(SnapMate_Pairs));

  //  traverse the jump_list and and replace the internal fragment IDs
  //  by external fragment IDs

  for (i=0; i<jumplistLength; i++) {
    if (existsUID(FRGmap, islMesg->jump_list[i].in1) == 0) {
      fprintf(stderr,"ISL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", islMesg->jump_list[i].in1);
      exit(1);
    }
      if (existsUID(FRGmap, islMesg->jump_list[i].in2) == 0) {
      fprintf(stderr,"ISL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", islMesg->jump_list[i].in2);
      exit(1);
    }

    slkMesg.jump_list[i].in1 = lookupUID(FRGmap, islMesg->jump_list[i].in1);
    slkMesg.jump_list[i].in2 = lookupUID(FRGmap, islMesg->jump_list[i].in2);

    slkMesg.jump_list[i].type = islMesg->jump_list[i].type;
  }

  pmesg->m = &slkMesg;
  pmesg->t = MESG_SLK;

  WriteProtoMesg_AS(fileOutput,pmesg);
}


void
convertISF(GenericMesg *pmesg,
           FILE        *fileOutput) {
  IntScaffoldMesg   *isfMesg = (IntScaffoldMesg *)pmesg->m;
  SnapScaffoldMesg   scfMesg;
  int                i;

  if (existsUID(ISFmap, isfMesg->iaccession)) {
    fprintf(stderr,"ISF: duplicate definition for scaffold ID "F_IID"\n", isfMesg->iaccession);
    exit(1);
  }

  scfMesg.iaccession       = isfMesg->iaccession;
  scfMesg.eaccession       = AS_UID_fromInteger(getUID(uids));

  scfMesg.num_contig_pairs = isfMesg->num_contig_pairs;
  scfMesg.contig_pairs     = NULL;

  if (scfMesg.num_contig_pairs > 0) {
    scfMesg.contig_pairs = (SnapContigPairs *)safe_malloc(scfMesg.num_contig_pairs * sizeof(SnapContigPairs));

    for(i=0; i<scfMesg.num_contig_pairs; i++) {
      if (existsUID(ICMmap, isfMesg->contig_pairs[i].contig1) == 0) {
	fprintf(stderr,"ISF: reference before definition for contig ID "F_IID"\n", isfMesg->contig_pairs[i].contig1);
	exit(1);
      }
      if (existsUID(ICMmap, isfMesg->contig_pairs[i].contig2) == 0) {
	fprintf(stderr,"ISF: reference before definition for contig ID "F_IID"\n", isfMesg->contig_pairs[i].contig2);
	exit(1);
      }

      //  decide if we've already put this contig in a scaffold

      if (lookupIID(ICMinISF1, isfMesg->contig_pairs[i].contig1)) {
	fprintf(stderr,"ISF: Contig1 ID "F_IID" already used in scaffold uid:%s\n",
                isfMesg->contig_pairs[i].contig1,
                AS_UID_toString(lookupUID(ICMinISF1, isfMesg->contig_pairs[i].contig1)));
        fprintf(stderr, "tried in scaffold iid="F_IID" uid=%s\n",
                scfMesg.iaccession,
                AS_UID_toString(scfMesg.eaccession));
	exit(1);
      }
      if (lookupIID(ICMinISF2, isfMesg->contig_pairs[i].contig2)) {
	fprintf(stderr,"ISF: Contig2 ID "F_IID" already used in scaffold uid:%s\n",
                isfMesg->contig_pairs[i].contig2,
                AS_UID_toString(lookupUID(ICMinISF2, isfMesg->contig_pairs[i].contig2)));
        fprintf(stderr, "tried in scaffold iid="F_IID" uid=%s\n",
                scfMesg.iaccession,
                AS_UID_toString(scfMesg.eaccession));
	exit(1);
      }

      SetAS_IID(ICMinISF1, isfMesg->contig_pairs[i].contig1, &isfMesg->iaccession);
      SetAS_IID(ICMinISF2, isfMesg->contig_pairs[i].contig2, &isfMesg->iaccession);

      scfMesg.contig_pairs[i].econtig1 = lookupUID(ICMmap, isfMesg->contig_pairs[i].contig1);
      scfMesg.contig_pairs[i].econtig2 = lookupUID(ICMmap, isfMesg->contig_pairs[i].contig2);
      scfMesg.contig_pairs[i].mean     = isfMesg->contig_pairs[i].mean;
      scfMesg.contig_pairs[i].stddev   = isfMesg->contig_pairs[i].stddev;
      scfMesg.contig_pairs[i].orient   = isfMesg->contig_pairs[i].orient;
    }

  } else {

    //  special case if there are zero contig pairs, then there is ONE
    //  with the second contig id -1

    scfMesg.contig_pairs = (SnapContigPairs*) safe_malloc(1 * sizeof(SnapContigPairs));

    if (existsUID(ICMmap, isfMesg->contig_pairs[0].contig1) == 0) {
      fprintf(stderr,"ISF: reference before definition for contig ID "F_IID"\n", isfMesg->contig_pairs[0].contig1);
      exit(1);
    }
    if (lookupIID(ICMinISF1, isfMesg->contig_pairs[0].contig1)) {
      fprintf(stderr,"ISF: Contig1 ID "F_IID" already used in scaffold uid:%s\n",
              isfMesg->contig_pairs[0].contig1,
              AS_UID_toString(lookupUID(ICMinISF1, isfMesg->contig_pairs[0].contig1)));
      fprintf(stderr, "tried in scaffold iid="F_IID" uid="F_IID"\n", scfMesg.iaccession, scfMesg.eaccession);
      exit(1);
    }
    if (lookupIID(ICMinISF2, isfMesg->contig_pairs[0].contig2)) {
      fprintf(stderr,"ISF: Contig2 ID "F_IID" already used in scaffold uid:%s\n",
              isfMesg->contig_pairs[0].contig2,
              AS_UID_toString(lookupUID(ICMinISF2, isfMesg->contig_pairs[0].contig2)));
      fprintf(stderr, "tried in scaffold iid="F_IID" uid="F_IID"\n", scfMesg.iaccession, scfMesg.eaccession);
      exit(1);
    }

    SetAS_IID(ICMinISF1, isfMesg->contig_pairs[0].contig1, &isfMesg->iaccession);
    SetAS_IID(ICMinISF2, isfMesg->contig_pairs[0].contig2, &isfMesg->iaccession);

    scfMesg.contig_pairs[0].econtig1 = lookupUID(ICMmap, isfMesg->contig_pairs[0].contig1);
    scfMesg.contig_pairs[0].econtig2 = lookupUID(ICMmap, isfMesg->contig_pairs[0].contig1);
    scfMesg.contig_pairs[0].mean     = isfMesg->contig_pairs[0].mean;
    scfMesg.contig_pairs[0].stddev   = isfMesg->contig_pairs[0].stddev;
    scfMesg.contig_pairs[0].orient   = isfMesg->contig_pairs[0].orient;
  }

  SetAS_UID(ISFmap, scfMesg.iaccession, &scfMesg.eaccession);

  pmesg->m = &scfMesg;
  pmesg->t = MESG_SCF;

  WriteProtoMesg_AS(fileOutput,pmesg);

  safe_free(scfMesg.contig_pairs);
}


void
convertIMD(GenericMesg     *pmesg,
           GateKeeperStore *gkpStore,
           FILE            *fileOutput) {
  IntMateDistMesg  *imdMesg = (IntMateDistMesg*) pmesg->m;
  SnapMateDistMesg  mdiMesg;

  mdiMesg.erefines    = getGateKeeperLibrary(gkpStore, imdMesg->refines)->libraryUID;
  mdiMesg.irefines    = imdMesg->refines;
  mdiMesg.min         = imdMesg->min;
  mdiMesg.max         = imdMesg->max;
  mdiMesg.mean        = imdMesg->mean;
  mdiMesg.stddev      = imdMesg->stddev;
  mdiMesg.num_buckets = imdMesg->num_buckets;
  mdiMesg.histogram   = imdMesg->histogram;

  pmesg->m = &mdiMesg;
  pmesg->t = MESG_MDI;

  WriteProtoMesg_AS(fileOutput,pmesg);
}



int main (int argc, char *argv[]) {
  char *outputPrefix       = NULL;
  char *gkpStoreName       = NULL;

  uint64      uidStart     = 0;

  GenericMesg *pmesg       = NULL;
  FILE        *fileOutput  = NULL;

  int numMSG = 0;
  int numIAF = 0;
  int numIAM = 0;
  int numIUM = 0;
  int numIUL = 0;
  int numICM = 0;
  int numICL = 0;
  int numISL = 0;
  int numISF = 0;
  int numIMD = 0;

  GateKeeperStore *gkpStore;
  FragStream      *fs;
  fragRecord       fr;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-s") == 0) {
        uidStart = strtoul(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-n") == 0) {
      SYS_UIDset_euid_namespace(argv[++arg]);
    } else if (strcmp(argv[arg], "-E") == 0) {
      SYS_UIDset_euid_server(argv[++arg]);

    } else if (strcmp(argv[arg], "-h") == 0) {
      err++;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((gkpStoreName == NULL) || (err)) {
    fprintf(stderr, "usage: %s -g gkpStore [-o prefix] [-s firstUID] [-n namespace] [-E server] [-h]\n", argv[0]);
    fprintf(stderr, "  -g gkpStore      mandatory path to the gatekeeper store\n");
    fprintf(stderr, "  -o prefix        write the output here, otherwise to stdout (.asm appended)\n");
    fprintf(stderr, "  -s firstUID      don't use real UIDs, but start counting from here\n");
    fprintf(stderr, "  -n namespace     use this UID namespace\n");
    fprintf(stderr, "  -E server        use this UID server\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Reads internal assembly messages (from scaffolder and consensus) on stdin,\n");
    fprintf(stderr, "  converts them to external assembly messages (mostly by just attaching a\n");
    fprintf(stderr, "  UID to each object) and writes them to a single output file.\n");
    exit(1);
  }

  ICMinISF1     = CreateVA_AS_IID(8192);
  ICMinISF2     = CreateVA_AS_IID(8192);

  IUMmap        = CreateVA_AS_UID(8192);
  ICMmap        = CreateVA_AS_UID(8192);
  ISFmap        = CreateVA_AS_UID(8192);
  FRGmap        = CreateVA_AS_UID(65536);
  DSTmap        = CreateVA_AS_UID(64);

  fprintf(stderr, "Reading gatekeeper store\n");

  gkpStore = openGateKeeperStore(gkpStoreName, FALSE);
  fs       = openFragStream(gkpStore, FRAG_S_INF);
  fragInfo = (fragInfo_t *)safe_calloc(getLastElemFragStore(gkpStore) + 1, sizeof(fragInfo_t));

  while (nextFragStream(fs, &fr)) {
    AS_IID    iid = getFragRecordIID(&fr);

    fragInfo[iid].loaded   = 1;
    fragInfo[iid].deleted  = getFragRecordIsDeleted(&fr);
    fragInfo[iid].clearBeg = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_LATEST);
    fragInfo[iid].clearEnd = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_LATEST);
    fragInfo[iid].uid      = getFragRecordUID(&fr);

    if ((uidStart > 0) &&
        (AS_UID_isString(fragInfo[iid].uid) == FALSE) &&
        (uidStart < AS_UID_toInteger(fragInfo[iid].uid)))
      uidStart = AS_UID_toInteger(fragInfo[iid].uid) + 1;
  }

  closeFragStream(fs);

  //  We still use the gkpStore for getting library info, so leave it open for now.

  if ((outputPrefix == NULL) || (strcmp(outputPrefix, "-") == 0)) {
    outputPrefix = NULL;
    fileOutput = stdout;
  } else {
    char N[FILENAME_MAX] = {0};
    sprintf(N, "%s.asm", outputPrefix);
    errno = 0;
    fileOutput = fopen(N, "w");
    if (errno) {
      fprintf(stderr, "%s: Couldn't open '%s' for write: %s\n", N, strerror(errno));
      exit(1);
    }
  }

  uids     = UIDserverInitialize(256, uidStart);

  fprintf(stderr, "Writing assembly file\n");

  while(ReadProtoMesg_AS(stdin,&pmesg) != EOF){
    switch(pmesg->t){
      case MESG_IAF :
        convertIAF(pmesg, fileOutput);
        numMSG += 1;
        numIAF++;
        break;
      case MESG_IAM :
        convertIAM(pmesg, fileOutput);
        numMSG += 1;
        numIAM++;
        break;
      case MESG_IUM :
        convertIUM(pmesg, fileOutput);
        numMSG += 19;
        numIUM++;
        break;
      case MESG_IUL :
        convertIUL(pmesg, fileOutput);
        numMSG += 7;
        numIUL++;
        break;
      case MESG_ICM :
        convertICM(pmesg, fileOutput);
        numMSG += 79;
        numICM++;
        break;
      case MESG_ICL :
        convertICL(pmesg, fileOutput);
        numMSG += 19;
        numICL++;
        break;
      case MESG_ISL :
        convertISL(pmesg, fileOutput);
        numMSG += 19;
        numISL++;
        break;
      case MESG_ISF :
        convertISF(pmesg, fileOutput);
        numMSG += 59;
        numISF++;
        break;
      case MESG_IMD :
        convertIMD(pmesg, gkpStore, fileOutput);
        numMSG += 102;
        numIMD++;
        break;
      default:
        break;
    }

    if (numMSG > 462583) {
      numMSG = 0;
      fprintf(stderr, "numIAF:%d numIAM:%d numIUM:%d numIUL:%d numICM:%d numICL:%d numISL:%d numISF:%d numIMD:%d\n",
              numIAF, numIAM, numIUM, numIUL, numICM, numICL, numISL, numISF, numIMD);
    }
  }

  fprintf(stderr, "numIAF:%d numIAM:%d numIUM:%d numIUL:%d numICM:%d numICL:%d numISL:%d numISF:%d numIMD:%d\n",
          numIAF, numIAM, numIUM, numIUL, numICM, numICL, numISL, numISF, numIMD);

  if (outputPrefix)
    fclose(fileOutput);

  fprintf(stderr, "Assembly file complete.\n");

  fprintf(stderr, "Writing IID to UID mapping files.\n");

  {
    FILE    *F = NULL;
    char     N[1024];
    sprintf(N, "%s.iidtouid", outputPrefix);
    errno = 0;
    F = fopen(N, "w");
    if (errno) {
      fprintf(stderr, "%s: Couldn't open '%s' for write: %s\n", argv[0], N, strerror(errno));
      exit(1);
    }

    DumpIID2UIDmap(FRGmap, "FRG", F);
    DumpIID2UIDmap(IUMmap, "UTG", F);
    DumpIID2UIDmap(ICMmap, "CTG", F);
    DumpIID2UIDmap(ISFmap, "SCF", F);
    DumpIID2UIDmap(DSTmap, "LIB", F);

    fclose(F);

    fprintf(stderr, "IID to UID mapping files complete.\n");
  }

  //  CANNOT close this before dumping IID to UID maps -- needed for
  //  UID lookups!
  closeGateKeeperStore(gkpStore);

  DeleteVA_AS_IID(ICMinISF1);
  DeleteVA_AS_IID(ICMinISF2);

  DeleteVA_AS_UID(IUMmap);
  DeleteVA_AS_UID(ICMmap);
  DeleteVA_AS_UID(ISFmap);
  DeleteVA_AS_UID(FRGmap);
  DeleteVA_AS_UID(DSTmap);

  return(0);
}

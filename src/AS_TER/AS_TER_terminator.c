
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

static const char CM_ID[] = "$Id: AS_TER_terminator.c,v 1.19 2007-05-23 15:32:56 skoren Exp $";

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
#include "AS_UTL_version.h"
#include "SYS_UIDclient.h"


typedef struct {
  int32       loaded:1;
  int32       deleted:1;
  int32       clearBeg:15;
  int32       clearEnd:15;
  CDS_UID_t   uid;
} fragInfo_t;

fragInfo_t    *fragInfo;

VA_DEF(CDS_UID_t);
VA_DEF(CDS_IID_t);

VA_TYPE(CDS_IID_t) *ICMinISF1;
VA_TYPE(CDS_IID_t) *ICMinISF2;

VA_TYPE(CDS_UID_t) *IUMmap;
VA_TYPE(CDS_UID_t) *ICMmap;
VA_TYPE(CDS_UID_t) *ISFmap;
VA_TYPE(CDS_UID_t) *FRGmap;
VA_TYPE(CDS_UID_t) *DSCmap;
VA_TYPE(CDS_UID_t) *DSTmap;

UIDserver       *uids;

void
DumpIID2UIDmap(VA_TYPE(CDS_UID_t) *map,
               char *execname,
               char *prefix,
               char *name,
               char *label) {
  FILE    *F = NULL;
  char     N[1024];
  int      i;

  sprintf(N, name, prefix);
  errno = 0;
  F = fopen(N, "w");
  if (errno) {
    fprintf(stderr, "%s: Couldn't open '%s' for write: %s\n", execname, N, strerror(errno));
    exit(1);
  }

  fprintf(F, "%s\n", label);

  for(i=0; i<GetNumCDS_UID_ts(map); i++){
    CDS_UID_t *di = GetCDS_UID_t(map,i);
    if (*di != 0)
      fprintf(F,"%d\t" F_U64 "\n",i,*di);
  }

  fclose(F);
}

CDS_UID_t
lookupUID(VA_TYPE(CDS_UID_t) *map, CDS_IID_t iid){
  CDS_UID_t *ret = GetCDS_UID_t(map, iid);
  if ((ret == NULL) || (*ret == 0))
    return(0);
  return(*ret);
}

CDS_IID_t
lookupIID(VA_TYPE(CDS_IID_t) *map, CDS_IID_t iid){
  CDS_IID_t *ret = GetCDS_IID_t(map, iid);
  if ((ret == NULL) || (*ret == 0))
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

  if (lookupUID(IUMmap, iafMesg->iaccession)) {
    fprintf(stderr, "IAF: Spotted FRG internal ID "F_IID" second time\n", iafMesg->iaccession);
    exit(1);
  }

  afgMesg.eaccession    = fragInfo[iafMesg->iaccession].uid;
  afgMesg.iaccession    = iafMesg->iaccession;
  afgMesg.mate_status   = iafMesg->mate_status;
  afgMesg.chimeric      = iafMesg->chimeric;
  afgMesg.chaff         = iafMesg->chaff;
  afgMesg.clear_rng.bgn = fragInfo[iafMesg->iaccession].clearBeg;
  afgMesg.clear_rng.end = fragInfo[iafMesg->iaccession].clearEnd;

  SetCDS_UID_t(FRGmap, iafMesg->iaccession, &afgMesg.eaccession);

  pmesg->m = &afgMesg;
  pmesg->t = MESG_AFG;

  WriteProtoMesg_AS(fileOutput, pmesg);
}


void
convertIUM(GenericMesg *pmesg,
           FILE        *fileOutput) {

  IntUnitigMesg  *iumMesg = (IntUnitigMesg *)pmesg->m;
  SnapUnitigMesg  utgMesg;
  int             i;

  if (lookupUID(IUMmap, iumMesg->iaccession)) {
    fprintf(stderr, "IUM: Spotted IUM internal ID "F_IID" second time\n", iumMesg->iaccession);
    exit(1);
  }

  utgMesg.eaccession     = getUID(uids);
  utgMesg.iaccession     = iumMesg->iaccession;
#ifdef AS_ENABLE_SOURCE
  utgMesg.source         = strdup(iumMesg->source);
#endif
  utgMesg.coverage_stat  = iumMesg->coverage_stat;
  utgMesg.status         = iumMesg->status;
  utgMesg.a_branch_point = iumMesg->a_branch_point;
  utgMesg.b_branch_point = iumMesg->b_branch_point;
  utgMesg.length         = iumMesg->length;
  utgMesg.consensus      = strdup(iumMesg->consensus);
  utgMesg.quality        = strdup(iumMesg->quality);
  utgMesg.forced         = iumMesg->forced;
  utgMesg.num_frags      = iumMesg->num_frags;
  utgMesg.num_vars       = 0;
  utgMesg.f_list         = NULL;
  utgMesg.v_list         = NULL;

  SetCDS_UID_t(IUMmap, iumMesg->iaccession, &utgMesg.eaccession);

  if (iumMesg->num_frags > 0) {
    utgMesg.f_list = (SnapMultiPos*)safe_malloc(iumMesg->num_frags * sizeof(SnapMultiPos));

    for(i=0; i<iumMesg->num_frags; i++){
      if (lookupUID(FRGmap, iumMesg->f_list[i].ident) == 0) {
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

  if (lookupUID(IUMmap, iulMesg->unitig1) == 0) {
    fprintf(stderr,"IUL: reference before definition error for unitig ID "F_IID"\n",iulMesg->unitig1);
    exit(1);
  }
  if (lookupUID(IUMmap, iulMesg->unitig2) == 0) {
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
    ulkMesg.jump_list[i].in1 = lookupUID(FRGmap, iulMesg->jump_list[i].in1);
    ulkMesg.jump_list[i].in2 = lookupUID(FRGmap, iulMesg->jump_list[i].in2);

    if (ulkMesg.jump_list[i].in1 == 0) {
      fprintf(stderr,"IUL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", iulMesg->jump_list[i].in1);
      exit(1);
    }
    if (ulkMesg.jump_list[i].in2 == 0) {
      fprintf(stderr,"IUL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", iulMesg->jump_list[i].in2);
      exit(1);
    }

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

  if (lookupUID(ICMmap, icmMesg->iaccession)) {
    fprintf(stderr, "ICM: Spotted ICM internal ID "F_IID" second time\n", icmMesg->iaccession);
    exit(1);
  }

  ccoMesg.eaccession = getUID(uids);
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

  SetCDS_UID_t(ICMmap, icmMesg->iaccession, &ccoMesg.eaccession);

  if (icmMesg->num_vars > 0) {
    ccoMesg.vars = (IntMultiVar*) safe_malloc(icmMesg->num_vars * sizeof(IntMultiVar));

    for(i=0; i<icmMesg->num_vars; i++) {
      ccoMesg.vars[i].position         = icmMesg->v_list[i].position;
      ccoMesg.vars[i].num_reads        = icmMesg->v_list[i].num_reads;
      ccoMesg.vars[i].num_conf_alleles = icmMesg->v_list[i].num_conf_alleles;
      ccoMesg.vars[i].anchor_size      = icmMesg->v_list[i].anchor_size;
      ccoMesg.vars[i].var_length       = icmMesg->v_list[i].var_length ;
      ccoMesg.vars[i].curr_var_id      = icmMesg->v_list[i].curr_var_id;
      ccoMesg.vars[i].phased_var_id    = icmMesg->v_list[i].phased_var_id;
      ccoMesg.vars[i].nr_conf_alleles  = icmMesg->v_list[i].nr_conf_alleles;
      ccoMesg.vars[i].weights          = icmMesg->v_list[i].weights;
      ccoMesg.vars[i].var_seq          = icmMesg->v_list[i].var_seq;
    }
  }

  if( ccoMesg.num_pieces > 0 ){
    ccoMesg.pieces = (SnapMultiPos*) safe_malloc(icmMesg->num_pieces*sizeof(SnapMultiPos));

    for(i=0; i<icmMesg->num_pieces; i++) {
      ccoMesg.pieces[i].type = icmMesg->pieces[i].type;
#ifdef AS_ENABLE_SOURCE
      ccoMesg.pieces[i].source = NULL;
#endif

      if (lookupUID(FRGmap, icmMesg->pieces[i].ident) == 0) {
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

      if (lookupUID(IUMmap, icmMesg->unitigs[i].ident) == 0) {
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

  if (lookupUID(ICMmap, iclMesg->contig1) == 0) {
    fprintf(stderr,"ICL: reference before definition error for contig ID "F_IID"\n",iclMesg->contig1);
    exit(1);
  }
  if (lookupUID(ICMmap, iclMesg->contig2) == 0) {
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
    clkMesg.jump_list[i].in1 = lookupUID(FRGmap, iclMesg->jump_list[i].in1);
    clkMesg.jump_list[i].in2 = lookupUID(FRGmap, iclMesg->jump_list[i].in2);

    if (clkMesg.jump_list[i].in1 == 0) {
      fprintf(stderr,"ICL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", iclMesg->jump_list[i].in1);
      exit(1);
    }
    if (clkMesg.jump_list[i].in2 == 0) {
      fprintf(stderr,"ICL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", iclMesg->jump_list[i].in2);
      exit(1);
    }

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

  if (lookupUID(ISFmap, islMesg->iscaffold1) == 0) {
    fprintf(stderr,"ISL: reference before definition error for scaffold ID "F_IID"\n", islMesg->iscaffold1);
    exit(1);
  }
  if (lookupUID(ISFmap, islMesg->iscaffold2) == 0) {
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
    slkMesg.jump_list[i].in1 = lookupUID(FRGmap, islMesg->jump_list[i].in1);
    slkMesg.jump_list[i].in2 = lookupUID(FRGmap, islMesg->jump_list[i].in2);

    if (slkMesg.jump_list[i].in1 == 0) {
      fprintf(stderr,"ISL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", islMesg->jump_list[i].in1);
      exit(1);
    }
    if (slkMesg.jump_list[i].in2 == 0) {
      fprintf(stderr,"ISL: Internal Fragment ID "F_IID" does not exist in Fragstore\n", islMesg->jump_list[i].in2);
      exit(1);
    }

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

  if (lookupUID(ISFmap, isfMesg->iaccession)) {
    fprintf(stderr,"ISF: duplicate definition for scaffold ID "F_IID"\n", isfMesg->iaccession);
    exit(1);
  }

  scfMesg.iaccession       = isfMesg->iaccession;
  scfMesg.eaccession       = getUID(uids);

  scfMesg.num_contig_pairs = isfMesg->num_contig_pairs;
  scfMesg.contig_pairs     = NULL;

  if (scfMesg.num_contig_pairs > 0) {
    scfMesg.contig_pairs = (SnapContigPairs *)safe_malloc(scfMesg.num_contig_pairs * sizeof(SnapContigPairs));

    for(i=0; i<scfMesg.num_contig_pairs; i++) {
      if (lookupUID(ICMmap, isfMesg->contig_pairs[i].contig1) == 0) {
	fprintf(stderr,"ISF: reference before definition for contig ID "F_IID"\n", isfMesg->contig_pairs[i].contig1);
	exit(1);
      }
      if (lookupUID(ICMmap, isfMesg->contig_pairs[i].contig2) == 0) {
	fprintf(stderr,"ISF: reference before definition for contig ID "F_IID"\n", isfMesg->contig_pairs[i].contig2);
	exit(1);
      }

      //  decide if we've already put this contig in a scaffold

      if (lookupIID(ICMinISF1, isfMesg->contig_pairs[i].contig1)) {
	fprintf(stderr,"ISF: Contig1 ID "F_IID" already used in scaffold uid:"F_UID"\n",
                isfMesg->contig_pairs[i].contig1, lookupUID(ICMinISF1, isfMesg->contig_pairs[i].contig1));
        fprintf(stderr, "tried in scaffold iid="F_IID" uid="F_IID"\n", scfMesg.iaccession, scfMesg.eaccession);
	exit(1);
      }
      if (lookupIID(ICMinISF2, isfMesg->contig_pairs[i].contig2)) {
	fprintf(stderr,"ISF: Contig2 ID "F_IID" already used in scaffold uid:"F_UID"\n",
                isfMesg->contig_pairs[i].contig2, lookupUID(ICMinISF2, isfMesg->contig_pairs[i].contig2));
        fprintf(stderr, "tried in scaffold iid="F_IID" uid="F_IID"\n", scfMesg.iaccession, scfMesg.eaccession);
	exit(1);
      }

      SetCDS_IID_t(ICMinISF1, isfMesg->contig_pairs[i].contig1, &isfMesg->iaccession);
      SetCDS_IID_t(ICMinISF2, isfMesg->contig_pairs[i].contig2, &isfMesg->iaccession);

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

    if (lookupUID(ICMmap, isfMesg->contig_pairs[0].contig1) == 0) {
      fprintf(stderr,"ISF: reference before definition for contig ID "F_IID"\n", isfMesg->contig_pairs[0].contig1);
      exit(1);
    }
    if (lookupIID(ICMinISF1, isfMesg->contig_pairs[0].contig1)) {
      fprintf(stderr,"ISF: Contig1 ID "F_IID" already used in scaffold uid:"F_UID"\n",
              isfMesg->contig_pairs[0].contig1, lookupUID(ICMinISF1, isfMesg->contig_pairs[0].contig1));
      fprintf(stderr, "tried in scaffold iid="F_IID" uid="F_IID"\n", scfMesg.iaccession, scfMesg.eaccession);
      exit(1);
    }
    if (lookupIID(ICMinISF2, isfMesg->contig_pairs[0].contig2)) {
      fprintf(stderr,"ISF: Contig2 ID "F_IID" already used in scaffold uid:"F_UID"\n",
              isfMesg->contig_pairs[0].contig2, lookupUID(ICMinISF2, isfMesg->contig_pairs[0].contig2));
      fprintf(stderr, "tried in scaffold iid="F_IID" uid="F_IID"\n", scfMesg.iaccession, scfMesg.eaccession);
      exit(1);
    }

    SetCDS_IID_t(ICMinISF1, isfMesg->contig_pairs[0].contig1, &isfMesg->iaccession);
    SetCDS_IID_t(ICMinISF2, isfMesg->contig_pairs[0].contig2, &isfMesg->iaccession);

    scfMesg.contig_pairs[0].econtig1 = lookupUID(ICMmap, isfMesg->contig_pairs[0].contig1);
    scfMesg.contig_pairs[0].econtig2 = lookupUID(ICMmap, isfMesg->contig_pairs[0].contig1);
    scfMesg.contig_pairs[0].mean     = isfMesg->contig_pairs[0].mean;
    scfMesg.contig_pairs[0].stddev   = isfMesg->contig_pairs[0].stddev;
    scfMesg.contig_pairs[0].orient   = isfMesg->contig_pairs[0].orient;
  }

  SetCDS_UID_t(ISFmap, scfMesg.iaccession, &scfMesg.eaccession);

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


void
convertIDS(GenericMesg *pmesg,
           FILE        *fileOutput) {

  IntDegenerateScaffoldMesg   *idsMesg = (IntDegenerateScaffoldMesg*) pmesg->m;
  SnapDegenerateScaffoldMesg   dscMesg;

  if (lookupUID(ICMmap, idsMesg->icontig) == 0) {
    fprintf(stderr,"IDS: reference before definition error for contig ID "F_IID"\n", idsMesg->icontig);
    exit(1);
  }
  if (lookupUID(DSCmap, idsMesg->icontig)) {
    fprintf(stderr,"IDS: duplicate definition for contig ID "F_IID"\n", idsMesg->icontig);
    exit(1);
  }

  dscMesg.eaccession = getUID(uids);
  dscMesg.econtig    = lookupUID(ICMmap, idsMesg->icontig);

  SetCDS_UID_t(DSCmap, idsMesg->icontig, &dscMesg.eaccession);

  pmesg->m = &dscMesg;
  pmesg->t = MESG_DSC;

  WriteProtoMesg_AS(fileOutput,pmesg);
}



int main (int argc, char *argv[]) {
  char *outputFileName     = NULL;
  char *mapFileName        = NULL;
  char *gkpStoreName       = NULL;

  uint64      uidStart     = 0;

  GenericMesg *pmesg       = NULL;
  FILE        *fileOutput  = NULL;

  int numMSG = 0;
  int numIAF = 0;
  int numIUM = 0;
  int numIUL = 0;
  int numICM = 0;
  int numICL = 0;
  int numISL = 0;
  int numISF = 0;
  int numIMD = 0;
  int numIDS = 0;

  GateKeeperStore *gkpStore;
  FragStream      *fs;
  fragRecord       fr;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];
    } else if (strcmp(argv[arg], "-m") == 0) {
      mapFileName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputFileName = argv[++arg];

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
    fprintf(stderr, "usage: %s -g gkpStore [-o output.asm] [-m mapprefix] [-s firstUID] [-n namespace] [-E server] [-h]\n", argv[0]);
    fprintf(stderr, "  -g gkpStore      mandatory path to the gatekeeper store\n");
    fprintf(stderr, "  -o output.asm    write the output here, otherwise to stdout\n");
    fprintf(stderr, "  -m mapprefix     write IID to UID mappings to files with this prefix\n");
    fprintf(stderr, "  -s firstUID      don't use real UIDs, but start counting from here\n");
    fprintf(stderr, "  -n namespace     use this UID namespace\n");
    fprintf(stderr, "  -E server        use this UID server\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Reads internal assembly messages (from scaffolder and consensus) on stdin,\n");
    fprintf(stderr, "  converts them to external assembly messages (mostly by just attaching a\n");
    fprintf(stderr, "  UID to each object) and writes them to a single output file.\n");
    exit(1);
  }

  ICMinISF1     = CreateVA_CDS_IID_t(8192);
  ICMinISF2     = CreateVA_CDS_IID_t(8192);

  IUMmap        = CreateVA_CDS_UID_t(8192);
  ICMmap        = CreateVA_CDS_UID_t(8192);
  ISFmap        = CreateVA_CDS_UID_t(8192);
  FRGmap        = CreateVA_CDS_UID_t(65536);
  DSCmap        = CreateVA_CDS_UID_t(8192);
  DSTmap        = CreateVA_CDS_UID_t(64);

  fprintf(stderr, "Reading gatekeeper store\n");

  gkpStore = openGateKeeperStore(gkpStoreName, FALSE);
  fs       = openFragStream(gkpStore, FRAG_S_INF);
  fragInfo = (fragInfo_t *)safe_calloc(getLastElemFragStore(gkpStore) + 1, sizeof(fragInfo_t));

  while (nextFragStream(fs, &fr)) {
    CDS_IID_t iid = getFragRecordIID(&fr);

    fragInfo[iid].loaded   = 1;
    fragInfo[iid].deleted  = getFragRecordIsDeleted(&fr);
    fragInfo[iid].clearBeg = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_LATEST);
    fragInfo[iid].clearEnd = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_LATEST);
    fragInfo[iid].uid      = getFragRecordUID(&fr);

    if ((uidStart > 0) && (uidStart < fragInfo[iid].uid))
      uidStart = fragInfo[iid].uid + 1;
  }

  closeFragStream(fs);

  //  We still use the gkpStore for getting library info, so leave it open for now.

  if ((outputFileName == NULL) || (strcmp(outputFileName, "-") == 0)) {
    outputFileName = NULL;
    fileOutput = stdout;
  } else {
    errno = 0;
    fileOutput = fopen(outputFileName, "w");
    if (errno) {
      fprintf(stderr, "%s: Couldn't open '%s' for write: %s\n", outputFileName, strerror(errno));
      exit(1);
    }
  }

  uids     = UIDserverInitialize(256, uidStart);

  fprintf(stderr, "Writing assembly file\n");

  while(ReadProtoMesg_AS(stdin,&pmesg) != EOF){
    switch(pmesg->t){
      case MESG_ADT :
        VersionStampADT((AuditMesg *)pmesg->m, argc, argv);
        WriteProtoMesg_AS(fileOutput,pmesg);
        break;
      case MESG_IAF :
        convertIAF(pmesg, fileOutput);
        numMSG += 1;
        numIAF++;
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
      case MESG_IDS :
        convertIDS(pmesg, fileOutput);
        numMSG += 89;
        numIDS++;
        break;
      default:
        break;
    }

    if (numMSG > 462583) {
      numMSG = 0;
      fprintf(stderr, "numIAF:%d numIUM:%d numIUL:%d numICM:%d numICL:%d numISL:%d numISF:%d numIMD:%d numIDS:%d\n",
              numIAF, numIUM, numIUL, numICM, numICL, numISL, numISF, numIMD, numIDS);
    }
  }

  fprintf(stderr, "numIAF:%d numIUM:%d numIUL:%d numICM:%d numICL:%d numISL:%d numISF:%d numIMD:%d numIDS:%d\n",
          numIAF, numIUM, numIUL, numICM, numICL, numISL, numISF, numIMD, numIDS);

  if (outputFileName)
    fclose(fileOutput);

  fprintf(stderr, "Assembly file complete.\n");

  closeGateKeeperStore(gkpStore);

  fprintf(stderr, "Writing IID to UID mapping files.\n");

  DumpIID2UIDmap(FRGmap, argv[0], mapFileName, "%s.fragment.iidtouid", "Fragment IID2UID map");
  DumpIID2UIDmap(IUMmap, argv[0], mapFileName, "%s.unitig.iidtouid",   "Unitig IID2UID map");
  DumpIID2UIDmap(ICMmap, argv[0], mapFileName, "%s.contig.iidtouid",   "Contig IID2UID map");
  DumpIID2UIDmap(ISFmap, argv[0], mapFileName, "%s.scaffold.iidtouid", "Scaffold IID2UID map");
  DumpIID2UIDmap(DSTmap, argv[0], mapFileName, "%s.distrib.iidtouid",  "Distrib IID2UID map");
  DumpIID2UIDmap(DSCmap, argv[0], mapFileName, "%s.degeneratecontig.iidtouid", "Degenerate Contig IID to Scaffold UID map");

  fprintf(stderr, "IID to UID mapping files complete.\n");

  DeleteVA_CDS_IID_t(ICMinISF1);
  DeleteVA_CDS_IID_t(ICMinISF2);

  DeleteVA_CDS_UID_t(IUMmap);
  DeleteVA_CDS_UID_t(ICMmap);
  DeleteVA_CDS_UID_t(ISFmap);
  DeleteVA_CDS_UID_t(FRGmap);
  DeleteVA_CDS_UID_t(DSCmap);
  DeleteVA_CDS_UID_t(DSTmap);

  return(0);
}

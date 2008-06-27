
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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignStore_CNS.h"

MultiAlignT *
CreateMultiAlignT(void) {
  MultiAlignT *ma = (MultiAlignT *)safe_calloc(1, sizeof(MultiAlignT));
  ma->maID      = -1;
  ma->consensus = NULL;
  ma->quality   = NULL;
  ma->fdelta    = NULL;
  ma->f_list    = NULL;
  ma->v_list    = NULL;
  ma->udelta    = NULL;
  ma->u_list    = NULL;
  //fprintf(stderr, "CreateMultiAlignT()--  ma 0x%016p created.\n", ma);
  return(ma);
}

MultiAlignT *
CreateEmptyMultiAlignT(void) {
  MultiAlignT *ma = (MultiAlignT *)safe_calloc(1, sizeof(MultiAlignT));
  ma->maID      = -1;
  ma->consensus = CreateVA_char(0);
  ma->quality   = CreateVA_char(0);
  ma->fdelta    = CreateVA_int32(0);;
  ma->f_list    = CreateVA_IntMultiPos(0);
  ma->v_list    = CreateVA_IntMultiVar(0);
  ma->udelta    = CreateVA_int32(0);
  ma->u_list    = CreateVA_IntUnitigPos(0);
  //fprintf(stderr, "CreateMultiAlignT()--  ma 0x%016p created.\n", ma);
  return(ma);
}

void
ClearMultiAlignT(MultiAlignT *ma) {
  ma->maID      = -1;
  ResetVA_char(ma->consensus);
  ResetVA_char(ma->quality);
  ResetVA_int32(ma->fdelta);
  ResetVA_IntMultiPos(ma->f_list);
  ResetVA_int32(ma->udelta);
  ResetVA_IntUnitigPos(ma->u_list);
  ResetVA_IntMultiVar(ma->v_list);
}

void
DeleteMultiAlignTWorker(MultiAlignT *ma) {

  if (ma == NULL)
    return;

  //fprintf(stderr, "DeleteMultiAlignTWorker()--  deleting ma 0x%016p\n", ma);

  if (ma->v_list) {
    int i;
    for (i=0; i<GetNumIntMultiVars(ma->v_list); i++){
      IntMultiVar *v = GetIntMultiVar(ma->v_list, i);
      safe_free(v->nr_conf_alleles);
      safe_free(v->weights);
      safe_free(v->var_seq);
      safe_free(v->conf_read_iids);
    }
  }

  DeleteVA_char(ma->consensus);
  DeleteVA_char(ma->quality);
  DeleteVA_int32(ma->fdelta);
  DeleteVA_IntMultiPos(ma->f_list);
  DeleteVA_int32(ma->udelta);
  DeleteVA_IntUnitigPos(ma->u_list);
  DeleteVA_IntMultiVar(ma->v_list);

  //  But first, trash it.
  ma->maID      = -1;
  ma->consensus = (void *)0xdeadbeef;
  ma->quality   = (void *)0xdeadbeef;
  ma->fdelta    = (void *)0xdeadbeef;
  ma->f_list    = (void *)0xdeadbeef;
  ma->udelta    = (void *)0xdeadbeef;
  ma->u_list    = (void *)0xdeadbeef;
  ma->v_list    = (void *)0xdeadbeef;

  safe_free(ma);
}

MultiAlignT *
CopyMultiAlignT(MultiAlignT *newma, MultiAlignT *oldma) {

  if (newma == NULL)
    newma = CreateMultiAlignT();

  if (newma->consensus == NULL) {
    newma->consensus = Clone_VA(oldma->consensus);
    newma->quality   = Clone_VA(oldma->quality);
    newma->fdelta    = Clone_VA(oldma->fdelta);
    newma->f_list    = Clone_VA(oldma->f_list);
    newma->udelta    = Clone_VA(oldma->udelta);
    newma->u_list    = Clone_VA(oldma->u_list);
    newma->v_list    = Clone_VA(oldma->v_list);
  } else {
    ReuseClone_VA(newma->consensus,oldma->consensus);
    ReuseClone_VA(newma->quality,  oldma->quality);
    ReuseClone_VA(newma->fdelta,   oldma->fdelta);
    ReuseClone_VA(newma->f_list,   oldma->f_list);
    ReuseClone_VA(newma->udelta,   oldma->udelta);
    ReuseClone_VA(newma->u_list,   oldma->u_list);
    ReuseClone_VA(newma->v_list,   oldma->v_list);
  }

  assert(oldma->maID != -1);
  newma->maID = oldma->maID;

  //  Adjust the delta pointers in the clone

  int32 *oldbase, *newbase, i;

  oldbase = Getint32(oldma->fdelta, 0);
  newbase = Getint32(newma->fdelta, 0);
  for (i=0; i<GetNumIntMultiPoss(oldma->f_list); i++){
    IntMultiPos *npos = GetIntMultiPos(newma->f_list,i);
    IntMultiPos *opos = GetIntMultiPos(oldma->f_list,i);
    if (opos->delta)
      npos->delta = newbase + (opos->delta - oldbase);
  }

  oldbase = Getint32(oldma->udelta, 0);
  newbase = Getint32(newma->udelta, 0);
  for (i=0; i<GetNumIntUnitigPoss(oldma->u_list); i++){
    IntUnitigPos *npos = GetIntUnitigPos(newma->u_list,i);
    IntUnitigPos *opos = GetIntUnitigPos(oldma->u_list,i);
    if (opos->delta)
      npos->delta = newbase + (opos->delta - oldbase);
  }

  for (i=0; i<GetNumIntMultiVars(oldma->v_list); i++) {
    IntMultiVar *ovar = GetIntMultiVar(oldma->v_list,i);
    IntMultiVar *nvar = GetIntMultiVar(newma->v_list,i);

    nvar->nr_conf_alleles = strdup(ovar->nr_conf_alleles);
    nvar->weights         = strdup(ovar->weights);
    nvar->var_seq         = strdup(ovar->var_seq);
    nvar->conf_read_iids  = strdup(ovar->conf_read_iids);
  }

  return(newma);
}



MultiAlignT *
CloneSurrogateOfMultiAlignT(MultiAlignT *oldma, int32 newNodeID) {
  MultiAlignT *newma     = CreateEmptyMultiAlignT();

  assert(GetNumIntUnitigPoss(oldma->u_list) == 1);

  // Surrogate has UNGAPPED consensus sequence.  As fragments
  // get added, gaps will rematerialize

  assert(newNodeID != -1);
  newma->maID = newNodeID;

  //  Copy the consensus and quality from the old multialign, also
  //  clone the unitig list.
  //
  GetMultiAlignUngappedConsensus(oldma, newma->consensus, newma->quality);
  ReuseClone_VA(newma->u_list, oldma->u_list);

  assert(GetNumIntUnitigPoss(newma->u_list) == 1);

  IntUnitigPos *u = GetIntUnitigPos(newma->u_list, 0);

  u->ident        = newNodeID;
  u->position.end = GetMultiAlignLength(newma);
  u->delta_length = 0;
  u->delta        = NULL;

  return(newma);
}


static
void
saveDeltaPointers(MultiAlignT *ma) {
  int32 *base, i;

  base = Getint32(ma->fdelta, 0);
  for (i=0; i<GetNumIntMultiPoss(ma->f_list); i++){
    IntMultiPos *pos = GetIntMultiPos(ma->f_list,i);

    if (pos->delta == NULL) {
      assert(pos->delta_length == 0);
    } else {
      assert(pos->delta >= base);
      assert(pos->delta_length > 0);
      pos->delta = (int32 *)(pos->delta - base);
    }
  }

  base = Getint32(ma->udelta, 0);
  for (i=0; i<GetNumIntUnitigPoss(ma->u_list); i++){
    IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);

    if (pos->delta == NULL) {
      assert(pos->delta_length == 0);
    } else {
      assert(pos->delta >= base);
      assert(pos->delta_length > 0);
      pos->delta = (int32 *)(pos->delta - base);
    }
  }
}

static
void
restoreDeltaPointers(MultiAlignT *ma) {
  int32 *base, i;

  base = Getint32(ma->fdelta, 0);
  for (i=0; i<GetNumIntMultiPoss(ma->f_list); i++){
    IntMultiPos *pos = GetIntMultiPos(ma->f_list,i);
    if (pos->delta_length > 0) {
      assert((long)pos->delta >= 0);
      pos->delta = base + (long)pos->delta;
    }
  }

  base = Getint32(ma->udelta, 0);
  for (i=0; i<GetNumIntUnitigPoss(ma->u_list); i++){
    IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);
    if (pos->delta_length > 0) {
      assert((long)pos->delta >= 0);
      pos->delta = base + (long)pos->delta;
    }
  }
}


void
SaveMultiAlignTToStream(MultiAlignT *ma, FILE *stream) {
  int32  maID = -1;

  if (ma) {
    assert(ma->maID != -1);
    maID = ma->maID;
  }

  //  If maID == -1, then the multialign isn't here.
  AS_UTL_safeWrite(stream, &maID, "SaveMultiAlignTToStream", sizeof(int32), 1);

  if (ma) {
    saveDeltaPointers(ma);

    CopyToFileVA_char(ma->consensus, stream);
    CopyToFileVA_char(ma->quality, stream);
    CopyToFileVA_int32(ma->fdelta, stream);
    CopyToFileVA_IntMultiPos(ma->f_list, stream);
    CopyToFileVA_int32(ma->udelta, stream);
    CopyToFileVA_IntMultiPos(ma->u_list, stream);
    CopyToFileVA_IntMultiVar(ma->v_list, stream);

    restoreDeltaPointers(ma);
  }
}

void
ReLoadMultiAlignTFromStream(FILE *stream, MultiAlignT *ma) {
  int32  maID   = 0;
  int    status = 0;

  assert(ma != NULL);
  ClearMultiAlignT(ma);

  status = AS_UTL_safeRead(stream, &maID, "ReLoadMultiAlignTFromStream", sizeof(int32), 1);
  assert(status == 1);

  //  If maID == -1, then the multialign isn't here.
  if (maID != -1) {
    ma->maID = maID;

    LoadFromFileVA_char(stream,ma->consensus);
    LoadFromFileVA_char(stream,ma->quality);
    LoadFromFileVA_int32(stream,ma->fdelta);
    LoadFromFileVA_IntMultiPos(stream,ma->f_list);
    LoadFromFileVA_int32(stream,ma->udelta);
    LoadFromFileVA_IntUnitigPos(stream,ma->u_list);
    LoadFromFileVA_IntMultiVar(stream,ma->v_list);

    restoreDeltaPointers(ma);

    CheckMAValidity(ma);
  }
}



//  used to return NULL if (!isPresent) or (*reference != NULLINDEX)

MultiAlignT *
LoadMultiAlignTFromStream(FILE *stream) {
  MultiAlignT *ma = CreateEmptyMultiAlignT();
  ReLoadMultiAlignTFromStream(stream, ma);
  return(ma);
}






void
CheckMAValidity(MultiAlignT *ma) {
  char *c      = Getchar(ma->consensus,0);
  char *q      = Getchar(ma->quality,0);
  char  v[256] = {0};

  assert(strlen(c) == strlen(q));

  v['a'] = v['c'] = v['g'] = v['t'] = v['-'] = v['n'] = 1;
  v['A'] = v['C'] = v['G'] = v['T'] = v['-'] = v['N'] = 1;

  for (; *c; c++, q++){
    assert(v[*c]);
    assert(*q >= '0' && *q <= 'l');
  }
}





void
GetMultiAlignUngappedConsensus(MultiAlignT *ma,
                               VA_TYPE(char) *ungappedConsensus,
                               VA_TYPE(char) *ungappedQuality) {

  ResetVA_char(ungappedConsensus);
  ResetVA_char(ungappedQuality);

  if (ma) {
    char *c = Getchar(ma->consensus,0);
    char *q = Getchar(ma->quality,0);

    CheckMAValidity(ma);

    while (*c != 0) {
      if (*c != '-') {
        Appendchar(ungappedConsensus, c);
        Appendchar(ungappedQuality, q);
      }
      c++;
      q++;
    }
  }

  char n  = '\0';
  Appendchar(ungappedConsensus, &n);
  Appendchar(ungappedQuality,   &n);
}


void
GetMultiAlignUngappedOffsets(MultiAlignT *ma,
                             VA_TYPE(int32) *ungappedOffsets) {

  Resetint32(ungappedOffsets);

  if (ma) {
    char *c = Getchar(ma->consensus, 0);
    int   u = 0;

    while (*c != 0) {
      Appendint32(ungappedOffsets, &u);
      if (*c != '-')
        u++;
      c++;
    }
    Appendint32(ungappedOffsets, &u);
  }
}



static int
CompareUnitigPos (const void *c1, const void *c2) {
  IntUnitigPos *u1 = (IntUnitigPos *)c1;
  IntUnitigPos *u2 = (IntUnitigPos *)c2;

  int diff = (MIN(u1->position.bgn, u1->position.end) -
              MIN(u2->position.bgn, u2->position.end));
  if (diff)
    return diff;

  diff = (MAX(u2->position.bgn, u2->position.end) -
          MAX(u1->position.bgn, u1->position.end));
  return(diff);
}

void
MakeCanonicalMultiAlignT(MultiAlignT *ma) {
  qsort(GetIntUnitigPos(ma->u_list,0),
        GetNumIntUnitigPoss(ma->u_list),
        sizeof(IntUnitigPos),
        CompareUnitigPos);
}



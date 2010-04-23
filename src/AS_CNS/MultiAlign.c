
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

static const char *rcsid = "$Id: MultiAlign.c,v 1.15 2010-04-23 15:11:06 brianwalenz Exp $";

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "MultiAlignment_CNS.h"

#include "splitToWords.H"

#undef  DEBUG_CREATE
#undef  DEBUG_FILES

MultiAlignT *
CreateMultiAlignT(void) {
  MultiAlignT *ma = (MultiAlignT *)safe_calloc(1, sizeof(MultiAlignT));

  ma->maID                       = -1;

  ma->data.unitig_coverage_stat  = 0.0;
  ma->data.unitig_microhet_prob  = 1.0;

  ma->data.unitig_status         = AS_UNASSIGNED;
  ma->data.unitig_unique_rept    = AS_FORCED_NONE;

  ma->data.contig_status         = AS_UNPLACED;

  ma->data.num_frags             = 0;
  ma->data.num_unitigs           = 0;

  ma->consensus = NULL;
  ma->quality   = NULL;

  ma->f_list    = NULL;
  ma->u_list    = NULL;
  ma->v_list    = NULL;

  ma->fdelta    = NULL;
  ma->udelta    = NULL;

#ifdef DEBUG_CREATE
  fprintf(stderr, "CreateMultiAlignT()--  ma 0x%016p created.\n", ma);
#endif
  return(ma);
}

MultiAlignT *
CreateEmptyMultiAlignT(void) {
  MultiAlignT *ma = (MultiAlignT *)safe_calloc(1, sizeof(MultiAlignT));

  ma->maID                       = -1;

  ma->data.unitig_coverage_stat  = 0.0;
  ma->data.unitig_microhet_prob  = 1.0;

  ma->data.unitig_status         = AS_UNASSIGNED;
  ma->data.unitig_unique_rept    = AS_FORCED_NONE;

  ma->data.contig_status         = AS_UNPLACED;

  ma->data.num_frags             = 0;
  ma->data.num_unitigs           = 0;

  ma->consensus = CreateVA_char(0);
  ma->quality   = CreateVA_char(0);

  ma->f_list    = CreateVA_IntMultiPos(0);
  ma->u_list    = CreateVA_IntUnitigPos(0);
  ma->v_list    = CreateVA_IntMultiVar(0);

  ma->fdelta    = CreateVA_int32(0);;
  ma->udelta    = CreateVA_int32(0);

#ifdef DEBUG_CREATE
  fprintf(stderr, "CreateEmptyMultiAlignT()--  ma 0x%016p created.\n", ma);
#endif
  return(ma);
}

void
ClearMultiAlignT(MultiAlignT *ma) {

  ma->maID                       = -1;

  ma->data.unitig_coverage_stat  = 0.0;
  ma->data.unitig_microhet_prob  = 1.0;

  ma->data.unitig_status         = AS_UNASSIGNED;
  ma->data.unitig_unique_rept    = AS_FORCED_NONE;

  ma->data.contig_status         = AS_UNPLACED;

  ma->data.num_frags             = 0;
  ma->data.num_unitigs           = 0;

  if (ma->v_list) {
    for (uint32 i=0; i<GetNumIntMultiVars(ma->v_list); i++) {
      IntMultiVar *v = GetIntMultiVar(ma->v_list, i);
      safe_free(v->alleles);
      safe_free(v->var_seq_memory);
      safe_free(v->read_id_memory);
    }
  }

  ResetVA_char(ma->consensus);
  ResetVA_char(ma->quality);

  ResetVA_IntMultiPos(ma->f_list);
  ResetVA_IntUnitigPos(ma->u_list);
  ResetVA_IntMultiVar(ma->v_list);

  ResetVA_int32(ma->fdelta);
  ResetVA_int32(ma->udelta);
}

void
DeleteMultiAlignTWorker(MultiAlignT *ma) {

  if (ma == NULL)
    return;

#ifdef DEBUG_CREATE
  fprintf(stderr, "DeleteMultiAlignTWorker()--  deleting ma 0x%016p\n", ma);
#endif

  if (ma->v_list) {
    for (uint32 i=0; i<GetNumIntMultiVars(ma->v_list); i++) {
      IntMultiVar *v = GetIntMultiVar(ma->v_list, i);
      safe_free(v->alleles);
      safe_free(v->var_seq_memory);
      safe_free(v->read_id_memory);
    }
  }

  DeleteVA_char(ma->consensus);
  DeleteVA_char(ma->quality);

  DeleteVA_IntMultiPos(ma->f_list);
  DeleteVA_IntUnitigPos(ma->u_list);
  DeleteVA_IntMultiVar(ma->v_list);

  DeleteVA_int32(ma->fdelta);
  DeleteVA_int32(ma->udelta);

  //  But first, trash it.
  memset(ma, 0xee, sizeof(MultiAlignT));

  safe_free(ma);
}

MultiAlignT *
CopyMultiAlignT(MultiAlignT *newma, MultiAlignT *oldma) {

  assert(oldma->maID != -1);

  if (newma == NULL)
    newma = CreateEmptyMultiAlignT();
  else
    ClearMultiAlignT(newma);

#ifdef DEBUG_CREATE
  fprintf(stderr, "CopyMultiAlignT()--  copy from ma 0x%016p to ma 0x%016p\n", ma, newma);
#endif

  //  Shallow copy first.

  newma->maID = oldma->maID;
  newma->data = oldma->data;

  //  Then a little deeper.

  ReuseClone_VA(newma->consensus,oldma->consensus);
  ReuseClone_VA(newma->quality,  oldma->quality);

  ReuseClone_VA(newma->f_list,   oldma->f_list);
  ReuseClone_VA(newma->u_list,   oldma->u_list);
  ReuseClone_VA(newma->v_list,   oldma->v_list);

  ReuseClone_VA(newma->fdelta,   oldma->fdelta);
  ReuseClone_VA(newma->udelta,   oldma->udelta);

  //  And rob the graves.

  int32 *oldbase, *newbase;

  oldbase = Getint32(oldma->fdelta, 0);
  newbase = Getint32(newma->fdelta, 0);
  for (uint32 i=0; i<GetNumIntMultiPoss(oldma->f_list); i++){
    IntMultiPos *npos = GetIntMultiPos(newma->f_list,i);
    IntMultiPos *opos = GetIntMultiPos(oldma->f_list,i);
    if (opos->delta)
      npos->delta = newbase + (opos->delta - oldbase);
  }

  oldbase = Getint32(oldma->udelta, 0);
  newbase = Getint32(newma->udelta, 0);
  for (uint32 i=0; i<GetNumIntUnitigPoss(oldma->u_list); i++){
    IntUnitigPos *npos = GetIntUnitigPos(newma->u_list,i);
    IntUnitigPos *opos = GetIntUnitigPos(oldma->u_list,i);
    if (opos->delta)
      npos->delta = newbase + (opos->delta - oldbase);
  }

  for (uint32 i=0; i<GetNumIntMultiVars(oldma->v_list); i++) {
    IntMultiVar *ovar = GetIntMultiVar(oldma->v_list,i);
    IntMultiVar *nvar = GetIntMultiVar(newma->v_list,i);

    int32 aSize = sizeof(IntVarAllele) * ovar->num_alleles;
    int32 vSize = sizeof(char)         * ovar->num_alleles * ovar->var_length + ovar->num_alleles;
    int32 rSize = sizeof(int32)        * ovar->num_reads;

    nvar->alleles        = (IntVarAllele *)safe_malloc(aSize);
    nvar->var_seq_memory = (char         *)safe_malloc(vSize);
    nvar->read_id_memory = (int32        *)safe_malloc(rSize);

    memcpy(nvar->alleles,        nvar->alleles,        aSize);
    memcpy(nvar->var_seq_memory, nvar->var_seq_memory, vSize);
    memcpy(nvar->read_id_memory, nvar->read_id_memory, rSize);
  }

  return(newma);
}



MultiAlignT *
CloneSurrogateOfMultiAlignT(MultiAlignT *oldma, int32 newNodeID) {
  MultiAlignT *newma     = CreateEmptyMultiAlignT();

#ifdef DEBUG_CREATE
  fprintf(stderr, "CloneSurrogateOfMultiAlignT()--  clone from ma 0x%016p to ma 0x%016p\n", oldma, newma);
#endif

  assert(GetNumIntUnitigPoss(oldma->u_list) == 1);

  // Surrogate has UNGAPPED consensus sequence.  As fragments
  // get added, gaps will rematerialize

  assert(newNodeID != -1);
  newma->maID = newNodeID;

  newma->data = oldma->data;

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
  int32 *base;

  base = Getint32(ma->fdelta, 0);
  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++){
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
  for (uint32 i=0; i<GetNumIntUnitigPoss(ma->u_list); i++){
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
  int32 *base;

  base = Getint32(ma->fdelta, 0);
  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++){
    IntMultiPos *pos = GetIntMultiPos(ma->f_list,i);
    if (pos->delta_length > 0) {
      assert((long)pos->delta >= 0);
      pos->delta = base + (long)pos->delta;
    }
  }

  base = Getint32(ma->udelta, 0);
  for (uint32 i=0; i<GetNumIntUnitigPoss(ma->u_list); i++){
    IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);
    if (pos->delta_length > 0) {
      assert((long)pos->delta >= 0);
      pos->delta = base + (long)pos->delta;
    }
  }
}


static
size_t
saveVARData(MultiAlignT *ma, char *&memory) {
  size_t  s = 0;

  //  If no memory pointer, we're just requesting the size needed
  //  Otherwise, copy into pre-allocated memory

  for (int32 i=0; i<GetNumIntMultiVars(ma->v_list); i++) {
    IntMultiVar  *imv = GetIntMultiVar(ma->v_list, i);

    s += sizeof(IntVarAllele) * imv->num_alleles;
    s += sizeof(char)         * imv->num_alleles * (imv->var_length + 1);
    s += sizeof(int32)        * imv->num_reads;

    if (memory) {
      memcpy(memory, imv->alleles, sizeof(IntVarAllele) * imv->num_alleles);
      memory += sizeof(IntVarAllele) * imv->num_alleles;

      memcpy(memory, imv->var_seq_memory, sizeof(char) * imv->num_alleles * (imv->var_length + 1));
      memory += sizeof(char) * imv->num_alleles * (imv->var_length + 1);

      memcpy(memory, imv->read_id_memory, sizeof(int32) * imv->num_reads);
      memory += sizeof(int32) * imv->num_reads;
    }
  }

  return(s);
}


static
void
restoreVARData(char *&memory, MultiAlignT *ma) {

  //  Allocate our memory.

  for (int32 i=0; i<GetNumIntMultiVars(ma->v_list); i++) {
    IntMultiVar  *imv = GetIntMultiVar(ma->v_list, i);

    //  The IMVs should already be empty, but our pointers are invalid and non-NULL.
    //safe_free(imv->alleles);
    //safe_free(imv->var_seq_memory);
    //safe_free(imv->read_id_memory);

    imv->alleles        = (IntVarAllele *)safe_malloc(sizeof(IntVarAllele) * imv->num_alleles);
    imv->var_seq_memory = (char         *)safe_malloc(sizeof(char)         * imv->num_alleles * (imv->var_length + 1));
    imv->read_id_memory = (int32        *)safe_malloc(sizeof(int32)        * imv->num_reads);
  }

  //  And copy data back into the correct spots in each VAR entry.

  for (int32 i=0; i<GetNumIntMultiVars(ma->v_list); i++) {
    IntMultiVar  *imv = GetIntMultiVar(ma->v_list, i);

    memcpy(imv->alleles, memory, sizeof(IntVarAllele) * imv->num_alleles);
    memory += sizeof(IntVarAllele) * imv->num_alleles;

    memcpy(imv->var_seq_memory, memory, sizeof(char) * imv->num_alleles * (imv->var_length + 1));
    memory += sizeof(char) * imv->num_alleles * (imv->var_length + 1);

    memcpy(imv->read_id_memory, memory, sizeof(int32) * imv->num_reads);
    memory += sizeof(int32) * imv->num_reads;
  }
}


void
SaveMultiAlignTToStream(MultiAlignT *ma, FILE *stream) {
  size_t   memorySize = 0;
  char    *memory     = NULL;
  char    *memoryBase = NULL;

  if (ma == NULL) {
    AS_UTL_safeWrite(stream, &memorySize, "SaveMultiAlignTToStream0", sizeof(size_t), 1);
    return;
  }

  assert(ma->maID != -1);

  saveDeltaPointers(ma);

  memorySize += sizeof(int32);
  memorySize += sizeof(MultiAlignD);
  memorySize += CopyToMemory_VA(ma->consensus, memory);
  memorySize += CopyToMemory_VA(ma->quality,   memory);

  memorySize += CopyToMemory_VA(ma->fdelta, memory);
  memorySize += CopyToMemory_VA(ma->udelta, memory);

  memorySize += CopyToMemory_VA(ma->f_list, memory);
  memorySize += CopyToMemory_VA(ma->u_list, memory);
  memorySize += CopyToMemory_VA(ma->v_list, memory);

  memorySize += saveVARData(ma, memory);


  memoryBase = memory = (char *)safe_malloc(sizeof(char) * memorySize);


  memcpy(memory, &ma->maID, sizeof(int32));
  memory += sizeof(int32);

  memcpy(memory, &ma->data, sizeof(MultiAlignD));
  memory += sizeof(MultiAlignD);

  CopyToMemoryVA_char(ma->consensus, memory);
  CopyToMemoryVA_char(ma->quality, memory);

  CopyToMemoryVA_int32(ma->fdelta, memory);
  CopyToMemoryVA_int32(ma->udelta, memory);

  CopyToMemoryVA_IntMultiPos(ma->f_list, memory);
  CopyToMemoryVA_IntUnitigPos(ma->u_list, memory);
  CopyToMemoryVA_IntMultiVar(ma->v_list, memory);

  restoreDeltaPointers(ma);

  saveVARData(ma, memory);

  AS_UTL_safeWrite(stream, &memorySize, "SaveMultiAlignTToStream0", sizeof(size_t), 1);
  AS_UTL_safeWrite(stream,  memoryBase, "SaveMultiAlignTToStream1", sizeof(char),   memorySize);

  safe_free(memoryBase);
}

void
ReLoadMultiAlignTFromStream(FILE *stream, MultiAlignT *ma) {
  size_t   memorySize = 0;
  char    *memory     = NULL;
  char    *memoryBase = NULL;
  size_t   status     = 0;

  assert(ma != NULL);

  ClearMultiAlignT(ma);

  status = AS_UTL_safeRead(stream, &memorySize, "ReLoadMultiAlignTFromStream0", sizeof(size_t), 1);
  assert(status == 1);

  if (memorySize == 0)
    return;

  memoryBase = memory = (char *)safe_malloc(sizeof(char) * memorySize);

  status = AS_UTL_safeRead(stream,  memoryBase, "ReLoadMultiAlignTFromStream1", sizeof(char), memorySize);
  assert(status == memorySize);

  memcpy(&ma->maID, memory, sizeof(int32));
  memory += sizeof(int32);

  memcpy(&ma->data, memory, sizeof(MultiAlignD));
  memory += sizeof(MultiAlignD);

  LoadFromMemoryVA_char(memory, ma->consensus);
  LoadFromMemoryVA_char(memory, ma->quality);

  LoadFromMemoryVA_int32(memory, ma->fdelta);
  LoadFromMemoryVA_int32(memory, ma->udelta);

  LoadFromMemoryVA_IntMultiPos(memory, ma->f_list);
  LoadFromMemoryVA_IntUnitigPos(memory, ma->u_list);
  LoadFromMemoryVA_IntMultiVar(memory, ma->v_list);

  restoreDeltaPointers(ma);

  restoreVARData(memory, ma);

  safe_free(memoryBase);
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

  if ((c == NULL) || (q == NULL))
    return;

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






void
DumpMultiAlignForHuman(FILE *out, MultiAlignT *ma, bool isUnitig) {

  char *cns = ma->consensus ? Getchar(ma->consensus, 0) : NULL;
  char *qlt = ma->quality   ? Getchar(ma->quality, 0)   : NULL;

  fprintf(out, "%s %d\n", (isUnitig) ? "unitig" : "contig", ma->maID);
  fprintf(out, "len %d\n", (cns) ? strlen(cns) : 0);
  fprintf(out, "cns %s\n", (cns) ? cns : "");
  fprintf(out, "qlt %s\n", (qlt) ? qlt : "");
  fprintf(out, "data.unitig_coverage_stat %f\n", ma->data.unitig_coverage_stat);
  fprintf(out, "data.unitig_microhet_prob %f\n", ma->data.unitig_microhet_prob);
  fprintf(out, "data.unitig_status        %c\n", ma->data.unitig_status);
  fprintf(out, "data.unitig_unique_rept   %c\n", ma->data.unitig_unique_rept);
  fprintf(out, "data.contig_status        %c\n", ma->data.contig_status);
  fprintf(out, "data.num_frags            %u\n", ma->data.num_frags);
  fprintf(out, "data.num_unitigs          %u\n", ma->data.num_unitigs);

  for (int32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    fprintf(stdout, "FRG type %c ident %9d container %9d parent %9d hang %6d %6d position %6d %6d\n",
            imp->type,
            imp->ident,
            imp->contained,
            imp->parent,
            imp->ahang, imp->bhang,
            imp->position.bgn, imp->position.end);
  }

  for (int32 i=0; i<GetNumIntUnitigPoss(ma->u_list); i++) {
    IntUnitigPos *iup = GetIntUnitigPos(ma->u_list, i);

    fprintf(stdout, "UTG type %c ident %9d position %6d %6d num_instances %d\n",
            iup->type,
            iup->ident,
            iup->position.bgn, iup->position.end,
            iup->num_instances);
  }
}


bool
LoadMultiAlignFromHumanGetLine(FILE *in, int32 LINElen, int32 LINEmax, char *LINE, splitToWords &W) {

  LINE[0] = 0;

 LoadMultiAlignFromHumanGetLineAgain:
  while ((fgets(LINE, LINEmax, in) != NULL) &&
         (LINE[0] == '#'))
    ;

  chomp(LINE);
  W.split(LINE);

  if (feof(in))
    //  Nothing more to read!
    return(false);

  if (LINE[0] == 0)
    //  Blank or whitespace only line.  Recurse or just goto, hmmmm.
    goto LoadMultiAlignFromHumanGetLineAgain;

  return(W.numWords() > 0);
}


bool
LoadMultiAlignFromHuman(MultiAlignT *ma, bool &isUnitig, FILE *in) {

  int32  LINElen = 0;
  int32  LINEmax = 1048576;
  char  *LINE    = new char [LINEmax];

  splitToWords  W;

  ClearMultiAlignT(ma);

  if (LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W) == 0)
    //  Nothing loaded.
    return(false);

  if        (strcmp(W[0], "unitig") == 0) {
    isUnitig = true;
    ma->maID = atoi(W[1]);
  } else if (strcmp(W[0], "contig") == 0) {
    isUnitig = false;
    ma->maID = atoi(W[1]);
  } else {
    fprintf(stderr, "MultiAlign not loaded:  Unknown MultiAlign type (expected 'unitig' or 'contig') in '%s'\n", LINE);
    exit(1);
  }

  if (LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W) == 0)
    fprintf(stderr, "LoadMultiAlignFromHuman()-- Expecting 'len' line, got empty line.\n"), exit(1);

  if (strcmp(W[0], "len") == 0) {
    int32 l = atoi(W[1]) + 1024;
    if (LINEmax < l) {
      delete [] LINE;
      LINEmax = l;
      LINE = new char [LINEmax];
    }
  } else {
    fprintf(stderr, "MultiAlign %d not loaded:  Unknown length in '%s'\n", ma->maID, LINE);
    exit(1);
  }

  if (LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W) == 0)
    fprintf(stderr, "MultiAlign %d not loaded:  Expecting 'cns' line, got empty line.\n", ma->maID), exit(1);

  if (strcmp(W[0], "cns") == 0) {
    if (W[1])
      SetRangeVA_char(ma->consensus, 0, strlen(W[1]) + 1, W[1]);
  } else {
    fprintf(stderr, "MultiAlign %d not loaded:  Unknown cns in '%s'\n", ma->maID, LINE);
    exit(1);
  }

  if (LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W) == 0)
    fprintf(stderr, "MultiAlign %d not loaded:  Expecting 'qlt' line, got empty line.\n", ma->maID), exit(1);

  if (strcmp(W[0], "qlt") == 0) {
    if (W[1])
      SetRangeVA_char(ma->quality, 0, strlen(W[1]) + 1, W[1]);
  } else {
    fprintf(stderr, "MultiAlign %d not loaded:  Unknown qlt in '%s'\n", ma->maID, LINE);
    exit(1);
  }

  ////////////////////////////////////////

  if (LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W) == 0)
    fprintf(stderr, "MultiAlign %d not loaded:  Expecting 'data.unitig_coverage_stat' line, got empty line.\n", ma->maID), exit(1);

  if (strcmp(W[0], "data.unitig_coverage_stat") == 0) {
    ma->data.unitig_coverage_stat = atof(W[1]);
  } else {
    fprintf(stderr, "MultiAlign %d not loaded:  Unknown data.unitig_converage_stat in '%s'\n", ma->maID, LINE);
    exit(1);
  }

  if (LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W) == 0)
    fprintf(stderr, "MultiAlign %d not loaded:  Expecting 'data.unitig_microhet_prob' line, got empty line.\n", ma->maID), exit(1);

  if (strcmp(W[0], "data.unitig_microhet_prob") == 0) {
    ma->data.unitig_microhet_prob = atof(W[1]);
  } else {
    fprintf(stderr, "MultiAlign %d not loaded:  Unknown data.unitig_microhet_prob in '%s'\n", ma->maID, LINE);
    exit(1);
  }

  if (LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W) == 0)
    fprintf(stderr, "MultiAlign %d not loaded:  Expecting 'data.unitig_status' line, got empty line.\n", ma->maID), exit(1);

  if (strcmp(W[0], "data.unitig_status") == 0) {
    switch (W[1][0]) {
      case AS_UNIQUE:
        ma->data.unitig_status = AS_UNIQUE;
        break;
      case AS_NOTREZ:
        ma->data.unitig_status = AS_NOTREZ;
        break;
      case AS_SEP:
        ma->data.unitig_status = AS_SEP;
        break;
      case AS_UNASSIGNED:
        ma->data.unitig_status = AS_UNASSIGNED;
        break;
      default:
        fprintf(stderr, "MultiAlign %d not loaded:  Unknown data.unitig_status in '%s'\n", ma->maID, LINE);
        exit(1);
        break;
    }
  } else {
    fprintf(stderr, "MultiAlign %d not loaded:  Unknown data.unitig_status in '%s'\n", ma->maID, ma->maID, LINE);
    exit(1);
  }

  if (LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W) == 0)
    fprintf(stderr, "MultiAlign %d not loaded:  Expecting 'data.unitig_unique_rept' line, got empty line.\n"), exit(1);

  if (strcmp(W[0], "data.unitig_unique_rept") == 0) {
    switch (W[1][0]) {
      case AS_FORCED_NONE:
        ma->data.unitig_unique_rept = AS_FORCED_NONE;
        break;
      case AS_FORCED_UNIQUE:
        ma->data.unitig_unique_rept = AS_FORCED_UNIQUE;
        break;
      case AS_FORCED_REPEAT:
        ma->data.unitig_unique_rept = AS_FORCED_REPEAT;
        break;
      default:
        fprintf(stderr, "MultiAlign %d not loaded:  Unknown data.unitig_unique_rept in '%s'\n", ma->maID, ma->maID, LINE);
        exit(1);
        break;
    }
  } else {
    fprintf(stderr, "MultiAlign %d not loaded:  Unknown data.unitig_unique_rept in '%s'\n", ma->maID, LINE);
    exit(1);
  }

  if (LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W) == 0)
    fprintf(stderr, "MultiAlign %d not loaded:  Expecting 'data.contig_status' line, got empty line.\n", ma->maID), exit(1);

  if (strcmp(W[0], "data.contig_status") == 0) {
    switch (W[1][0]) {
      case AS_PLACED:
        ma->data.contig_status = AS_PLACED;
        break;
      case AS_UNPLACED:
        ma->data.contig_status = AS_UNPLACED;
        break;
      default:
        fprintf(stderr, "MultiAlign %d not loaded:  Unknown data.contig_status in '%s'\n", ma->maID, LINE);
        exit(1);
        break;
    }
  } else {
    fprintf(stderr, "Unknown data.contig_status in '%s'\n", ma->maID, LINE);
    exit(1);
  }

  if (LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W) == 0)
    fprintf(stderr, "MultiAlign %d not loaded:  Expecting 'data.num_frags' line, got empty line.\n"), exit(1);

  if (strcmp(W[0], "data.num_frags") == 0) {
    ma->data.num_frags = atoi(W[1]);
  } else {
    fprintf(stderr, "MultiAlign %d not loaded:  Unknown data.num_frags in '%s'\n", ma->maID, ma->maID, LINE);
    exit(1);
  }

  if (LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W) == 0)
    fprintf(stderr, "MultiAlign %d not loaded:  Expecting 'data.num_unitigs' line, got empty line.\n"), exit(1);

  if (strcmp(W[0], "data.num_unitigs") == 0) {
    ma->data.num_unitigs = atoi(W[1]);
  } else {
    fprintf(stderr, "MultiAlign %d not loaded:  Unknown data.num_unitigs in '%s'\n", ma->maID, ma->maID, LINE);
    exit(1);
  }

  ////////////////////////////////////////

  ResetToRange_VA(ma->f_list, ma->data.num_frags);
  ResetToRange_VA(ma->u_list, ma->data.num_unitigs);

  ////////////////////////////////////////

  for (int32 i=0; i<ma->data.num_frags; i++) {
    IntMultiPos  *imp = GetIntMultiPos(ma->f_list, i);

    while ((LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W)) &&
           (W[0][0] == '#'))
      ;

    if (W.numWords() == 0) {
      fprintf(stderr, "MultiAlign %d not loaded:  Too few FRG lines; read %d, expected %d\n", ma->maID, i, ma->data.num_frags);
      exit(1);
    }

    if ((W.numWords() != 15) ||
        (strcmp(W[ 0], "FRG")       != 0) ||
        (strcmp(W[ 1], "type")      != 0) ||
        (strcmp(W[ 3], "ident")     != 0) ||
        (strcmp(W[ 5], "container") != 0) ||
        (strcmp(W[ 7], "parent")    != 0) ||
        (strcmp(W[ 9], "hang")      != 0) ||
        (strcmp(W[12], "position")  != 0)) {
      fprintf(stderr, "MultiAlign %d not loaded:  Unknown FRG line in '%s'\n", ma->maID, LINE);
      exit(1);
    }

    switch (W[2][0]) {
      case AS_READ:
        imp->type = AS_READ;
        break;
      default:
        fprintf(stderr, "MultiAlign %d not loaded:  Unknown FRG type %c/%d in '%s'\n", ma->maID, W[2][0], W[2][0], ma->maID, LINE);
        exit(1);
        break;
    }

    imp->ident        = atoi(W[4]);
    imp->contained    = atoi(W[6]);
    imp->parent       = atoi(W[8]);
    imp->ahang        = atoi(W[10]);
    imp->bhang        = atoi(W[11]);
    imp->position.bgn = atoi(W[13]);
    imp->position.end = atoi(W[14]);
  }

  ////////////////////////////////////////

  for (int32 i=0; i<ma->data.num_unitigs; i++) {
    IntUnitigPos *iup = GetIntUnitigPos(ma->u_list, i);

    while ((LoadMultiAlignFromHumanGetLine(in, LINElen, LINEmax, LINE, W)) &&
           (W[0][0] == '#'))
      ;

    if (W.numWords() == 0) {
      fprintf(stderr, "MultiAlign %d not loaded:  Too few UTG lines; read %d, expected %d\n", ma->maID, i, ma->data.num_unitigs);
      exit(1);
    }

    if ((W.numWords() != 10) ||
        (strcmp(W[0], "UTG")           != 0) ||
        (strcmp(W[1], "type")          != 0) ||
        (strcmp(W[3], "ident")         != 0) ||
        (strcmp(W[5], "position")      != 0) ||
        (strcmp(W[8], "num_instances") != 0)) {
      fprintf(stderr, "MultiAlign not loaded:  Unknown UTG line in '%s'\n", ma->maID, LINE);
      exit(1);
    }

    switch (W[2][0]) {
      case AS_UNIQUE_UNITIG:
        iup->type = AS_UNIQUE_UNITIG;
        break;
      case AS_ROCK_UNITIG:
        iup->type = AS_ROCK_UNITIG;
        break;
      case AS_STONE_UNITIG:
        iup->type = AS_STONE_UNITIG;
        break;
      case AS_PEBBLE_UNITIG:
        iup->type = AS_PEBBLE_UNITIG;
        break;
      case AS_SINGLE_UNITIG:
        iup->type = AS_SINGLE_UNITIG;
        break;
      case AS_OTHER_UNITIG:
        iup->type = AS_OTHER_UNITIG;
        break;
      default:
        fprintf(stderr, "MultiAlign %d not loaded:  Unknown UTG type %c/%d in '%s'\n", ma->maID, W[2][0], W[2][0], ma->maID, LINE);
        exit(1);
        break;
    }

    iup->ident         = atoi(W[4]);
    iup->position.bgn  = atoi(W[6]);
    iup->position.end  = atoi(W[7]);
    iup->num_instances = atoi(W[9]);
  }

  return(true);
}


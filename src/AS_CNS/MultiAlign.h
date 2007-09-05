
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

#ifndef MULTIALIGN_H
#define MULTIALIGN_H

#include "AS_MSG_pmesg.h"
#include "AS_UTL_Var.h"
#include "AS_PER_gkpStore.h"

VA_DEF(char);
VA_DEF(int32);

typedef struct {
  int32                   maID;
  VA_TYPE(char)          *consensus;  // gapped consensus
  VA_TYPE(char)          *quality;    // quality
  VA_TYPE(int32)         *fdelta;     // deltas for all fragments in f_list
  VA_TYPE(IntMultiPos)   *f_list;     // positions of fragments
  VA_TYPE(int32)         *udelta;     // deltas for all unitigs in u_list
  VA_TYPE(IntUnitigPos)  *u_list;     // positions of unitigs
  VA_TYPE(IntMultiVar)   *v_list;     // variations                  
} MultiAlignT;
VA_DEF(MultiAlignT)


MultiAlignT *CreateMultiAlignT(void);
MultiAlignT *CreateEmptyMultiAlignT(void);
void         ClearMultiAlignT(MultiAlignT *multiAlign);

#define      DeleteMultiAlignT(M) { DeleteMultiAlignTWorker(M); (M) = NULL; }
void         DeleteMultiAlignTWorker(MultiAlignT *multiAlign);

MultiAlignT *CreateMultiAlignTFromIUM(IntUnitigMesg *ium, int32 localFragID, int sequenceOnly);
MultiAlignT *CreateMultiAlignTFromICM(IntConConMesg *ium, int32 localFragID, int sequenceOnly);
MultiAlignT *CreateMultiAlignTFromCCO(SnapConConMesg *ium, int32 localFragID, int sequenceOnly);

//  Copies oldma into newma.  If newma is NULL, a new one is allocated.
//  Both cases return the copy.
MultiAlignT *CopyMultiAlignT(MultiAlignT *newma, MultiAlignT *oldma);

MultiAlignT *CloneSurrogateOfMultiAlignT(MultiAlignT *oldMA, int32 newNodeID);

void         SaveMultiAlignTToStream(MultiAlignT *ma, FILE *stream);
MultiAlignT *LoadMultiAlignTFromStream(FILE *stream);
void         ReLoadMultiAlignTFromStream(FILE *stream, MultiAlignT *ma);

void         CheckMAValidity(MultiAlignT *ma);

void         GetMultiAlignUngappedConsensus(MultiAlignT *ma, VA_TYPE(char) *ungappedSequence, VA_TYPE(char) *ungappedQuality);
void         GetMultiAlignUngappedOffsets(MultiAlignT *ma, VA_TYPE(int32) *ungappedOffsets);

void         MakeCanonicalMultiAlignT(MultiAlignT *ma);

void         PrintMultiAlignT(FILE *out,
                              MultiAlignT *ma,
                              GateKeeperStore *gkp_store,
                              int show_qv,
                              int dots,
                              uint32 clrrng_flag);

static
size_t
GetMultiAlignTMemorySize(MultiAlignT *ma) {
  size_t size = 0;
  if (ma)
    size = (GetMemorySize_VA(ma->consensus) * 2 +
            GetMemorySize_VA(ma->fdelta) +
            GetMemorySize_VA(ma->f_list) +
            GetMemorySize_VA(ma->udelta) +
            GetMemorySize_VA(ma->u_list) +
            GetMemorySize_VA(ma->v_list));
  return(size);
}



static
int32 
GetMultiAlignLength(MultiAlignT *ma) {
  // don't include the space for the null character
  return((int32)GetNumchars(ma->consensus) - 1);
}

static
int32 
GetMultiAlignUngappedLength(MultiAlignT *ma) {
  int32   u = 0;
  char   *c = Getchar(ma->consensus,0);

  while (*c) {
    if (*c != '-')
      u++;
    c++;
  }
  return(u);
}



static
IntUnitigPos *
GetAendUnitigPos(MultiAlignT *ma) {
  return(GetIntUnitigPos(ma->u_list,0));
}

static
IntUnitigPos *
GetBendUnitigPos(MultiAlignT *ma) {
  long length = GetMultiAlignLength(ma); 
  int  i;

  for(i = GetNumIntUnitigPoss(ma->u_list) - 1; i>=0; i--) {
    IntUnitigPos *pos = GetIntUnitigPos(ma->u_list,i);
    if (MAX(pos->position.bgn, pos->position.end) == length)
      return(pos);
  }
  return(NULL);
}

#endif //  MULTIALIGN_H


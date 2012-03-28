
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

static const char *rcsid_MULTIALIGN_H = "$Id: MultiAlign.h,v 1.12 2012-03-28 06:10:55 brianwalenz Exp $";

#include "AS_MSG_pmesg.h"
#include "AS_UTL_Var.h"
#include "AS_PER_gkpStore.h"

VA_DEF(char)
VA_DEF(int32)

typedef struct {
  double                     unitig_coverage_stat;
  double                     unitig_microhet_prob;

  UnitigStatus               unitig_status;
  UnitigFUR                  unitig_unique_rept;

  ContigStatus               contig_status;

  uint32                     num_frags;
  uint32                     num_unitigs;
} MultiAlignD;

typedef struct {
  int32                      maID;
  MultiAlignD                data;

  VA_TYPE(char)             *consensus;  // gapped consensus
  VA_TYPE(char)             *quality;    // gapped quality

  VA_TYPE(IntMultiPos)      *f_list;     // positions of fragments
  VA_TYPE(IntUnitigPos)     *u_list;     // positions of unitigs
  VA_TYPE(IntMultiVar)      *v_list;     // variations

  VA_TYPE(int32)            *fdelta;     // deltas for all fragments in f_list
  VA_TYPE(int32)            *udelta;     // deltas for all unitigs in u_list
} MultiAlignT;


MultiAlignT *CreateMultiAlignT(void);
MultiAlignT *CreateEmptyMultiAlignT(void);
void         ClearMultiAlignT(MultiAlignT *multiAlign);

#define      DeleteMultiAlignT(M) do { DeleteMultiAlignTWorker(M); (M) = NULL; } while (0)
void         DeleteMultiAlignTWorker(MultiAlignT *multiAlign);

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

void         DumpMultiAlignForHuman(FILE *out, MultiAlignT *ma, bool isUnitig);
bool         LoadMultiAlignFromHuman(MultiAlignT *ma, bool &isUnitig, FILE *in);

void         PrintMultiAlignT(FILE *out,
                              MultiAlignT *ma,
                              gkStore *gkp_store,
                              int32 show_qv,
                              int32 dots,
                              uint32 clrrng_flag);

int32        GetMultiAlignLength(MultiAlignT *ma);
int32        GetMultiAlignUngappedLength(MultiAlignT *ma);

#endif //  MULTIALIGN_H


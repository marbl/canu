
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
#ifndef MULTIALIGNMENT_CNS_INCLUDE
#define MULTIALIGNMENT_CNS_INCLUDE

static const char *rcsid_MULTIALIGNMENT_CNS_INCLUDE = "$Id: MultiAlignment_CNS.h,v 1.59 2011-12-04 23:46:58 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_MSG_pmesg.h"
#include "AS_ALN_aligners.h"

#include "MultiAlign.h"
#include "MultiAlignStore.h"

extern int32 DUMP_UNITIGS_IN_MULTIALIGNCONTIG;
extern int32 VERBOSE_MULTIALIGN_OUTPUT;

#define CNS_OPTIONS_SPLIT_ALLELES_DEFAULT  1
#define CNS_OPTIONS_MIN_ANCHOR_DEFAULT    11
#define CNS_OPTIONS_DO_PHASING_DEFAULT     1

typedef struct {
  int32 split_alleles;
  int32 smooth_win;
  int32 do_phasing;
} CNS_Options;

typedef enum {
  CNS_SMOOTH = 1, // only eliminate pairwise construction artifacts
  CNS_POLYX  = 2, // align poly-X regions
  CNS_INDEL  = 4  // push apart mushed block indels
}  CNS_RefineLevel;


MultiAlignT *MergeMultiAlignsFast_new(VA_TYPE(IntElementPos) *,
                                      CNS_Options *opp);

MultiAlignT *ReplaceEndUnitigInContig(uint32,
                                      uint32,
                                      int,
                                      CNS_Options *opp);

bool
MultiAlignUnitig(MultiAlignT   *ma,
                 gkStore       *fragStore,
                 CNS_Options   *opp,
                 int32         *failed);


bool
MultiAlignContig(MultiAlignT   *ma,
                 gkStore       *fragStore,
                 CNS_Options   *opp);


//  Options to things in MultiAligment_CNS.c

extern int32 allow_neg_hang;

#endif

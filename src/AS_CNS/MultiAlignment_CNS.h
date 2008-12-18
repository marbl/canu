
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

static const char *rcsid_MULTIALIGNMENT_CNS_INCLUDE = "$Id: MultiAlignment_CNS.h,v 1.43 2008-12-18 07:13:22 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_MSG_pmesg.h"
#include "AS_SDB_SequenceDB.h"
#include "MultiAlignStore_CNS.h"
#include "AS_ALN_aligners.h"

#define CNS_MIN_QV 0
#define CNS_MAX_QV 60
// 1/28/2008 - include the N character in the allowed alphabet
#define CNS_NALPHABET 7
#define CNS_NALPHABET_NULL_ARRAY {0, 0, 0, 0, 0, 0, 0}
#define CNS_NP 32

#define MIN_SIZE_OF_MANODE 10000
#define BC_MAX(a,b)  (((a)>(b))?(a):(b))
#define BC_MIN(a,b)  (((a)<(b))?(a):(b))
#define MIN_ALLOCATED_DEPTH 100

//  This is probably broken, or extremely inefficient, as of Nov 4 2007.
#undef PRINTUIDS

extern int DUMP_UNITIGS_IN_MULTIALIGNCONTIG;
extern int VERBOSE_MULTIALIGN_OUTPUT;

#define CNS_OPTIONS_SPLIT_ALLELES_DEFAULT  1
#define CNS_OPTIONS_MIN_ANCHOR_DEFAULT    11

typedef struct {
  int split_alleles;
  int smooth_win;
} CNS_Options;

typedef enum {
  CNS_QUIET       = (int)'Q', // quiet,  print nothing
  CNS_STATS_ONLY  = (int)'S', // print only 1-line statistic summary
  CNS_ALIGNMENT   = (int)'A', // print the multialignment, sans CNS
  CNS_CONSENSUS   = (int)'C', // print the multialignment, with CNS
  CNS_DOTS        = (int)'D', // print the multialignment, dot format
  CNS_NODOTS      = (int)'N', // print the multialignment, "nodot" format
  CNS_EDIT_SCORE  = (int)'E', // print the edit score column by column
  CNS_VIEW_UNITIG = (int)'U',  // show the unitigs in the contig alignment
  CNS_VERBOSE     = (int)'V'  // verbose pre-post refinment output
} CNS_PrintKey;   // determine the format for PrintAlignment

typedef enum {
  CNS_SMOOTH = 1, // only eliminate pairwise construction artifacts
  CNS_POLYX  = 2, // align poly-X regions
  CNS_INDEL  = 4  // push apart mushed block indels
}  CNS_RefineLevel;

typedef enum  {
  AS_CONSENSUS   = (int) 'C',
  AS_OVERLAY     = (int) 'O',
  AS_MERGE       = (int) 'M',
  AS_NOALIGNTYPE = (int) 'U'
} ALIGN_TYPE;

extern ALIGN_TYPE ALIGNMENT_CONTEXT;



MultiAlignT *MergeMultiAlignsFast_new(tSequenceDB *,
                                      GateKeeperStore *,
                                      VA_TYPE(IntElementPos) *,
                                      int,
                                      int,
                                      AS_ALN_Aligner *COMPARE_FUNC,
                                      CNS_Options *opp);

MultiAlignT *ReplaceEndUnitigInContig(tSequenceDB *,
                                      GateKeeperStore * ,
                                      uint32,
                                      uint32,
                                      int,
                                      AS_ALN_Aligner *COMPARE_FUNC,
                                      CNS_Options *opp);


int MultiAlignUnitig(IntUnitigMesg *,
                     GateKeeperStore *,
                     VA_TYPE(char) *,
                     VA_TYPE(char) *,
                     VA_TYPE(int32) *,
                     CNS_PrintKey,
                     AS_ALN_Aligner *COMPARE_FUNC,
                     CNS_Options *opp);

int MultiAlignContig(IntConConMesg *,
                     VA_TYPE(char) *,
                     VA_TYPE(char) *,
                     VA_TYPE(int32) *,
                     CNS_PrintKey ,
                     AS_ALN_Aligner *COMPARE_FUNC,
                     CNS_Options *opp);

int MultiAlignContig_ReBasecall(MultiAlignT *,
                                VA_TYPE(char) *,
                                VA_TYPE(char) *,
                                CNS_Options *);

void SequenceComplement(char *sequence, char *quality);

//  Options to things in MultiAligment_CNS.c

extern int allow_forced_frags;
extern int allow_neg_hang;

#endif

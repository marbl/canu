
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

static const char *rcsid_MULTIALIGN_H = "$Id$";

#include "AS_UTL_Var.H"
#include "AS_PER_gkpStore.H"

#include <vector>

using namespace std;

VA_DEF(char)
VA_DEF(int32)

typedef enum {
  AS_READ    = (int)'R',  //  Celera Read
  AS_EXTR    = (int)'X',  //  External WGS read
  AS_TRNR    = (int)'T',  //  Transposon library read
  AS_UNITIG  = (int)'U',  //  Assembled unitig
  AS_CONTIG  = (int)'C'   //  Assembled contig
} FragType;

typedef enum {
  AS_UNIQUE =     (int)'U',
  AS_NOTREZ =     (int)'N',
  AS_SEP =        (int)'S',
  AS_UNASSIGNED = (int)'X'
} UnitigStatus;

typedef enum {
  AS_PLACED	= (int)'P',
  AS_UNPLACED   = (int)'U'
} ContigStatus;

typedef enum {
  AS_UNIQUE_UNITIG   = (int)'U',  // U-Unique
  AS_ROCK_UNITIG     = (int)'R',  // Rock
  AS_STONE_UNITIG    = (int)'S',  // Stone
  AS_PEBBLE_UNITIG   = (int)'P',  // Pebble
  AS_SINGLE_UNITIG   = (int)'s',  // Singleton Unitig in Unplaced Contig
  AS_OTHER_UNITIG    = (int)'X'   // Unspecified surrogate unitig
} UnitigType;

typedef struct IntMultiPos {
  FragType        type;
  AS_IID          ident;
  AS_IID          contained;
  AS_IID          parent;     //  IID of the fragment we align to

  int32           ahang;      //  If parent defined, these are relative
  int32           bhang;      //  that fragment

  SeqInterval     position;
  int32           delta_length;
  int32          *delta;
} IntMultiPos;  //  IMP

typedef struct IntVarAllele {
  int32           num_reads;
  int32           weight;
  int32           var_seq_offset;  //  offset into var_seq_memory for this allele
  int32           read_id_offset;  //  offset into read_id_memory for this allele
} IntVarAllele;

typedef struct IntMultiVar {
  int32           var_id;                 // id of current VAR record
  int32           phased_id;              // id of the VAR record phased with the current one
  SeqInterval     position;               // position of the var region
  int32           num_reads;              // num reads total
  int32           num_alleles;            // num alleles total
  int32           num_alleles_confirmed;  // num confirmed alleles
  int32           min_anchor_size;        //
  int32           var_length;             // bases in the var region
  IntVarAllele   *alleles;                // list of num_alleles alleles

  char           *var_seq_memory;         // single allocation for all memory
  int32          *read_id_memory;         // single allocation for all memory

  char           *enc_num_reads;          //  the nra field
  char           *enc_weights;            //  the wgt field
  char           *enc_var_seq;            //  the seq field
  char           *enc_read_ids;           //  the rid field
} IntMultiVar;  //  IMV

typedef struct {
  UnitigType    type;
  AS_IID        ident;
  SeqInterval   position;
  int32         num_instances;
  int32         delta_length;
  int32        *delta;
} IntUnitigPos;  //  IUP

VA_DEF(IntMultiPos)
VA_DEF(IntMultiVar)
VA_DEF(IntUnitigPos)



typedef struct {
  double                     unitig_coverage_stat;
  double                     unitig_microhet_prob;

  UnitigStatus               unitig_status;

  uint32                     unitig_suggest_repeat:1;
  uint32                     unitig_suggest_unique:1;
  uint32                     unitig_force_repeat:1;
  uint32                     unitig_force_unique:1;
  uint32                     spare:28;

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

void         GetMultiAlignUngappedConsensus(MultiAlignT *ma, char *ungappedSequence, char *ungappedQuality);
void         GetMultiAlignUngappedConsensus(MultiAlignT *ma, VA_TYPE(char) *ungappedSequence, VA_TYPE(char) *ungappedQuality);

void         GetMultiAlignGapToUngap(MultiAlignT *ma, vector<int32> &map);
void         GetMultiAlignUngapToGap(MultiAlignT *ma, vector<int32> &map);

void         MakeCanonicalMultiAlignT(MultiAlignT *ma);

void         DumpMultiAlignForHuman(FILE *out, MultiAlignT *ma, bool isUnitig);
bool         LoadMultiAlignFromHuman(MultiAlignT *ma, bool &isUnitig, FILE *in);

void         PrintMultiAlignT(FILE *out,
                              MultiAlignT *ma,
                              gkStore *gkp_store,
                              int32 show_qv,
                              int32 dots,
                              uint32 clrrng_flag);

int32        GetMultiAlignLength(MultiAlignT *ma, bool force=false);
int32        GetMultiAlignUngappedLength(MultiAlignT *ma);

#endif //  MULTIALIGN_H


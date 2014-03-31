
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
/* $Id$   */

#ifndef AS_MSG_PMESG_INCLUDE_H
#define AS_MSG_PMESG_INCLUDE_H

static const char *rcsid_AS_MSG_PMESG_INCLUDE_H = "$Id$";

#include <cstdio>
#include <ctime>

using namespace std;

#include "AS_global.H"
#include "AS_MSG_types.H"
#include "AS_UTL_UID.H"
#include "AS_UTL_Var.H"

//
//  Generic message object handle
//

typedef enum {
  AS_ADD      = (int)'A',
  AS_DELETE   = (int)'D',
  AS_IGNORE   = (int)'I',
  AS_UPDATE   = (int)'U'
} ActionType;

typedef enum {
  MESG_NUL = 0,
  MESG_BAT, MESG_VER, MESG_DST, MESG_LIB, MESG_FRG, MESG_LKG, MESG_PLC,
  MESG_OVL,
  MESG_UOM,
  MESG_IMD, MESG_IAF, MESG_IAM, MESG_IUM, MESG_IUL, MESG_ICM, MESG_ICL, MESG_ISF, MESG_ISL,
  MESG_MDI, MESG_AFG, MESG_AMP, MESG_UTG, MESG_ULK, MESG_CCO, MESG_CLK, MESG_SCF, MESG_SLK,
  MESG_EOF
} MessageType;

#define NUM_OF_REC_TYPES MESG_EOF

static const char  *MessageTypeName[NUM_OF_REC_TYPES + 1] = {
  "NUL",
  "BAT", "VER", "DST", "LIB", "FRG", "LKG", "PLC",
  "OVL",
  "UOM",
  "IMD", "IAF", "IAM", "IUM", "IUL", "ICM", "ICL", "ISF", "ISL",
  "MDI", "AFG", "AMP", "UTG", "ULK", "CCO", "CLK", "SCF", "SLK",
  "EOF"
};

typedef struct {
  void         *m;         // A pointer to the message object
  MessageType  t;          // The message type
} GenericMesg;


//
//  Common Data Structures
//

typedef struct {
  int32 bgn;
  int32 end;
} SeqInterval;

//
//  Internal Messages
//

typedef struct InternalBatchMesgTag {
  char         *name;
  AS_UID        eaccession;
  char         *comment;
} BatchMesg;  //  BAT

typedef struct {
  uint32     version;
} VersionMesg;  //  VER

typedef struct {
  ActionType      action;
  LinkType        type;
  PairOrient      link_orient;
  AS_UID          frag1;
  AS_UID          frag2;
  AS_UID          distance;
} LinkMesg;  //  LKG

typedef struct {
  ActionType      action;
  AS_UID          frag;
  AS_UID          bound1;
  AS_UID          bound2;
} PlacementMesg;  //  PLC, only version 2

typedef struct {
  ActionType   action;
  AS_UID       eaccession;
  double       mean;
  double       stddev;
  char        *source;
  PairOrient   link_orient;
  uint32       num_features;
  char       **features;
  char       **values;
} LibraryMesg;  //  LIB, only version 2

typedef struct {
  ActionType   action;
  AS_UID       eaccession;
  double       mean;
  double       stddev;
} DistanceMesg;  //  DST, only version 1

typedef struct {
  ActionType       action;
  uint32                version;
  AS_UID         eaccession;
  AS_UID                library_uid;     //  only version 2
  IntLibrary_ID         library_iid;     //  only version 2
  AS_UID                plate_uid;       //  only version 2
  uint32                plate_location;  //  only version 2
  FragType         type;            //  only version 1
  uint32                is_random;       //  only version 2
  char                  status_code;     //  only version 2
  SeqInterval      clear_rng;
  SeqInterval      clear_vec;       //  only version 2
  SeqInterval      clear_max;       //  only version 2
  SeqInterval           contamination;   //  only version 2
  char            *source;
  char            *sequence;
  char            *quality;
  char                  *hps;            //  only version 2
  IntFragment_ID     iaccession;
} FragMesg;  //  FRG

typedef FragMesg InternalFragMesg;

typedef struct {
  IntFragment_ID   aifrag, bifrag;
  int32            ahg, bhg;
  PairOrient       orientation;
  OverlapType      overlap_type;
  double           quality;
  int32            min_offset, max_offset;
  int32            polymorph_ct;
  int32            *alignment_trace;
#ifdef AS_MSG_USE_OVL_DELTA
  //  This is left in Just In Case someone really wants to output
  //  an overlap delta from overlapper.
  signed char      *alignment_delta;
#endif
} OverlapMesg;  //  OVL

typedef struct {
  IntChunk_ID    chunk1;
  IntChunk_ID       chunk2;
  PairOrient            orient;
  UnitigOverlapType  overlap_type;
  int32                 best_overlap_length;
  int32                 min_overlap_length;
  int32                 max_overlap_length;
  double                quality;
} UnitigOverlapMesg;  //  UOM

typedef struct {
  IntChunk_ID     iaccession;
  int32           bp_length;
  double          coverage_stat;
  int32           num_frags;
  int32           a_degree;
  int32           b_degree;
} ChunkMesg;

typedef struct IntMultiPos {
  FragType        type;
  IntFragment_ID  ident;
  IntFragment_ID  contained;
  IntFragment_ID  parent;     //  IID of the fragment we align to

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
  IntUnitig_ID  ident;
  SeqInterval   position;
  int32         num_instances;
  int32         delta_length;
  int32        *delta;
} IntUnitigPos;  //  IUP

typedef struct {
  UnitigType   type;
  SeqInterval  position;
  int32        delta_length;
  int32        *delta;
  AS_UID        eident;
} UnitigPos;

typedef struct {
  FragType        type;
  IntFragment_ID  ident;
  SeqInterval     position;
} IntElementPos;  //  IEP

typedef struct {
  IntChunk_ID     iaccession;
  double          coverage_stat;
  double          microhet_prob;
  UnitigStatus    status;
  char            unique_rept;   //  NOT USED?
  int32           length;
  char            *consensus;
  char            *quality;
  int32           forced;
  int32           num_frags;
  IntMultiPos    *f_list;
} IntUnitigMesg;  //  IUM

typedef struct {
  AS_UID          eaccession;
  UnitigStatus    status;
  int32           num_occurences;
  AS_UID          *occurences;
  int32           length;
  char            *consensus;
  char            *quality;
  int32           num_reads;
  int32           num_guides;
  int32           sum_delta_lengths;
} UnitigMesg;

typedef struct {
  IntContig_ID               iaccession;
  ContigStatus               placed;
  int32                      length;
  char                       *consensus;
  char                       *quality;
  int32                 forced;
  int32                      num_pieces;
  int32                      num_unitigs;
  int32                      num_vars;
  IntMultiPos               *pieces;
  IntUnitigPos              *unitigs;
  IntMultiVar               *v_list;
} IntConConMesg;  //  ICM

typedef struct {
  IntFragment_ID  in1, in2;
  LinkType        type;
} IntMate_Pairs;

typedef struct {
  IntChunk_ID    unitig1;
  IntChunk_ID    unitig2;
  PairOrient            orientation;
  UnitigOverlapType  overlap_type;
  int32      is_possible_chimera;
  double     mean_distance;
  double     std_deviation;
  int32      num_contributing;
  PlacementStatusType  status;
  IntMate_Pairs    *jump_list;
} IntUnitigLinkMesg;  //  IUL

typedef struct {
  IntChunk_ID    contig1;
  IntChunk_ID    contig2;
  PairOrient            orientation;
  UnitigOverlapType  overlap_type;
  int32      is_possible_chimera;
  double     mean_distance;
  double     std_deviation;
  int32      num_contributing;
  PlacementStatusType  status;
  IntMate_Pairs    *jump_list;
} IntContigLinkMesg;  //  ICL

typedef struct {
  IntScaffold_ID  iscaffold1;
  IntScaffold_ID  iscaffold2;
  PairOrient            orientation;
  double     mean_distance;
  double     std_deviation;
  int32      num_contributing;
  IntMate_Pairs    *jump_list;
} InternalScaffoldLinkMesg;  //  ISL

typedef struct {
  IntFragment_ID  iaccession;
  FragType        type;
  int32           chimeric_NOTUSED;
  int32           chaff;
  SeqInterval     clear_rng;
  MateStatType    mate_status;
} IntAugFragMesg;  //  IAF

typedef struct {
  AS_UID          eaccession;
  IntFragment_ID  iaccession;
  MateStatType    mate_status;
  int32           chimeric_NOTUSED;
  int32           chaff;
  SeqInterval     clear_rng;
} AugFragMesg;  //  AFG

typedef struct {
  AS_IID          fragment1;
  AS_IID          fragment2;
  MateStatType    mate_status;
} IntAugMatePairMesg;  //  IAM

typedef struct {
  AS_UID          fragment1;
  AS_UID          fragment2;
  MateStatType    mate_status;
} AugMatePairMesg;  //  AMP

typedef struct {
  IntDist_ID  refines;
  double      mean;
  double      stddev;
  int32        min;
  int32        max;
  int32        num_buckets;
  int32       *histogram;
} IntMateDistMesg;  //  IMD

typedef struct {
  IntContig_ID    contig1;
  IntContig_ID    contig2;
  double     mean;
  double     stddev;
  PairOrient            orient;
} IntContigPairs;  //  ICP

typedef struct {
  IntScaffold_ID  iaccession;
  int32      num_contig_pairs;
  IntContigPairs  *contig_pairs;
} IntScaffoldMesg;  //  ISF

//
//  Genome Snapshot Messages
//

typedef struct {
  FragType       type;
  AS_UID         eident;
  SeqInterval    position;
  int32          delta_length;
  int32         *delta;
} SnapMultiPos;  //  MPS

typedef struct {
  AS_UID          eaccession;  // changed in comparison to internal message
  IntChunk_ID     iaccession;
  double          coverage_stat;
  double          microhet_prob;
  UnitigStatus    status;
  int32           length;
  char            *consensus;
  char            *quality;
  int32      forced;
  int32           num_frags;
  int32           num_vars;
  SnapMultiPos    *f_list;// changed in comparison to internal message
  IntMultiVar     *v_list;
} SnapUnitigMesg;  //  UTG

typedef struct {
  AS_UID       in1, in2;
  LinkType     type;
} SnapMate_Pairs;

typedef struct {
  AS_UID       eunitig1;
  AS_UID       eunitig2;
  PairOrient            orientation;
  UnitigOverlapType  overlap_type;
  int32      is_possible_chimera;
  double     mean_distance;
  double     std_deviation;
  int32      num_contributing;
  PlacementStatusType  status;
  SnapMate_Pairs       *jump_list; // changed in comparison to internal message
} SnapUnitigLinkMesg;  //  ULK

typedef struct {
  AS_UID                      eaccession;
  IntContig_ID                iaccession;
  ContigStatus                placed;
  int32                       length;
  char                       *consensus;
  char                       *quality;
  int32                       forced;
  int32                       num_pieces;
  int32                       num_unitigs;
  int32                       num_vars;
  SnapMultiPos               *pieces; // changed in comparison to internal message
  IntMultiVar                *vars;
  UnitigPos                  *unitigs;// changed in comparison to internal message
} SnapConConMesg;  //  CCO

typedef struct {
  AS_UID       econtig1; // changed in comparison to internal message
  AS_UID       econtig2; // changed in comparison to internal message
  PairOrient            orientation;
  UnitigOverlapType  overlap_type;
  int32      is_possible_chimera;
  double     mean_distance;
  double     std_deviation;
  int32      num_contributing;
  PlacementStatusType  status;
  SnapMate_Pairs  *jump_list; // changed in comparison to internal message
} SnapContigLinkMesg;  //  CLK

typedef struct {
  AS_UID                escaffold1;
  AS_UID                escaffold2;
  PairOrient            orientation;
  double     mean_distance;
  double     std_deviation;
  int32      num_contributing;
  SnapMate_Pairs  *jump_list;
} SnapScaffoldLinkMesg;  //  SLK

typedef struct {
  AS_UID       econtig1; // changed in comparison to internal message
  AS_UID       econtig2; // changed in comparison to internal message
  double     mean;
  double     stddev;
  PairOrient            orient;
} SnapContigPairs;  //  CTP

typedef struct {
  AS_UID                eaccession;
  IntScaffold_ID        iaccession;
  int32      num_contig_pairs;
  SnapContigPairs     *contig_pairs; // changed in comparison to internal message
} SnapScaffoldMesg;  //  SCF

typedef struct {
  AS_UID       erefines; // changed in comparison to internal message
  IntDist_ID    irefines; // changed in comparison to internal message
  double    mean;
  double    stddev;
  int32      min;
  int32      max;
  int32      num_buckets;
  int32      *histogram;
} SnapMateDistMesg;  //  MDI

typedef struct EndOfFileMesgTag {
  int32   status;
  char    *comment;
} EndOfFileMesg;  //  EOF



VA_DEF(IntMultiPos)
VA_DEF(IntMultiVar)
VA_DEF(IntUnitigPos)
VA_DEF(IntElementPos)
VA_DEF(IntUnitigMesg)  //  Used by unitigger.


//  Returns a number in the range [1, NUM_OF_REC_TYPES -1]
//  as a function of the first 3 characters of the passed string.
//
int GetMessageType(char *string);

//   Returns a string as a function of message type
//
const char  *GetMessageName(int type);

//  External Routines

void       AS_MSG_resetProtoLineNum(void);
void       AS_MSG_setFormatVersion(int format);

int        ReadProtoMesg_AS(FILE *fin, GenericMesg **pmesg);
void       WriteProtoMesg_AS(FILE *fout, GenericMesg *mesg);

#endif  /* AS_MSG_PMESG_INCLUDE */

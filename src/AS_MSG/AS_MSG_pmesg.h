
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
/* $Id: AS_MSG_pmesg.h,v 1.8 2005-08-02 02:37:39 gdenisov Exp $   */

#ifndef AS_MSG_PMESG_INCLUDE
#define AS_MSG_PMESG_INCLUDE

#include <stdio.h>
#include <time.h>

#include "AS_global.h"

/* Exported constant definitions; Macro definitions; type definitions */

// Defining the following enables internal source fields for testing
#define AS_ENABLE_SOURCE

/* All externally created accession numbers are 64 bits. */

/* This supports a backward compatibility mode for DROS data */
typedef enum{
  AS_DROS_MODE =  (int)'D',
  AS_HUMAN_MODE = (int)'H'
}ProtoIOMode;

void SetProtoMode_AS(ProtoIOMode mode);
ProtoIOMode GetProtoMode_AS(void);


#define AS_FRAG_MAX_LEN (2048)
#define AS_FRAG_MIN_LEN (64)

#define AS_BACTIG_MAX_LEN (500000)
#define AS_BACTIG_MIN_LEN (150)

#define AS_BAC_MAX_LEN (500000)
#define AS_BAC_MIN_LEN (150)

#ifdef GENERIC_STORE_USE_LONG_STRINGS
#define AS_READ_MAX_LEN AS_BACTIG_MAX_LEN
#else
#define AS_READ_MAX_LEN AS_FRAG_MAX_LEN
#endif
#define AS_READ_MIN_LEN AS_FRAG_MIN_LEN
#define DEFINE_IDs(type)\
typedef CDS_UID_t type##_ID;\
typedef CDS_IID_t Int##type##_ID;

// Traditional UID types
DEFINE_IDs(Fragment)
DEFINE_IDs(Locale)
DEFINE_IDs(Distance)
DEFINE_IDs(ScreenItem)
DEFINE_IDs(Chunk)
DEFINE_IDs(Unitig)
DEFINE_IDs(Contig)
DEFINE_IDs(Dist)
DEFINE_IDs(Scaffold)

// New UID types for human
DEFINE_IDs(Repeat)
DEFINE_IDs(Batch)
DEFINE_IDs(Bactig)
DEFINE_IDs(Bac)
DEFINE_IDs(Bin)
DEFINE_IDs(Sequence)
DEFINE_IDs(Plate)
DEFINE_IDs(Library)
DEFINE_IDs(Donor)
// Special case for Well IDs, since external numbers refer to plate position
// Internal numbers are internally unique
typedef uint16 Well_ID;
typedef uint32 IntWell_ID;

typedef enum {
 AS_ADD    = (int)'A',
 AS_DELETE = (int)'D',
 AS_UPDATE = (int)'U',
 AS_REDEFINE = (int)'R'
} ActionType;

typedef struct { CDS_COORD_t bgn, end; } SeqInterval;

  /* Keep this in sync with recordtypes array */

typedef enum {
  MESG_NONE = 0,
  MESG_ADT,
  MESG_FRG,
  MESG_IFG,
  MESG_SFG,
  MESG_OFG, // 5
  MESG_LKG,
  MESG_ILK,
  MESG_DST,
  MESG_IDT,
  MESG_SCN,// 10
  MESG_ISN,
  MESG_RPT,
  MESG_OVL,
  MESG_BRC,
  MESG_UOM,// 15
  MESG_IUM,
  MESG_IUL,
  MESG_ICL,
  MESG_AFG,
  MESG_ISF,// 20
  MESG_IMD,
  MESG_IAF,
  MESG_UTG,
  MESG_ULK,
  MESG_ICM,// 25
  MESG_CCO,
  MESG_CLK,
  MESG_SCF,
  MESG_MDI,
  // NEW FOR HUMAN
  MESG_BAT,// 30  
  MESG_IBA,
  MESG_BAC,
  MESG_IBC,
  MESG_BIN,
  MESG_IBI,// 35
  MESG_PLA,
  MESG_LKP,
  MESG_SP0, // spare 0
  MESG_SP1, // spare 1
  MESG_IRP,// 40
  MESG_IDS,
  MESG_DSC,
  MESG_SLK,
  MESG_ISL,
  MESG_FOM,// 45
  MESG_OFR,
  MESG_BUG,
  MESG_LIB,
  MESG_SP2, // spare 2
  MESG_EOF // 50
} MessageType;


/**** NOTE: If you add a message type, update the following!!!!! ****/

#define NUM_OF_REC_TYPES MESG_EOF

static char  *MessageTypeName[NUM_OF_REC_TYPES + 1] = {
  "NONE",
  "ADT",
  "FRG",
  "IFG",
  "SFG",
  "OFG", // 5
  "LKG",
  "ILK",
  "DST",
  "IDT",
  "SCN", // 10
  "ISN",
  "RPT", 
  "OVL", 
  "BRC",
  /* "CHK",               */
  /*  "RTG",     Obsolete */
  /* "CTG",      Obsolete */
  "UOM",// 15
  "IUM", 
  "IUL", 
  "ICL",
  "AFG",
  "ISF",// 20  
  "IMD", 
  "IAF", 
  "UTG", 
  "ULK",
  "ICM",// 25 
  "CCO", 
  "CLK",
  "SCF", 
  "MDI",
  // NEW FOR HUMAN
  "BAT",// 30  
  "IBA",
  "BAC",
  "IBC", 
  "BIN",
  "IBI",// 35
  "PLA", 
  "LKP",
  "SP0", // spare 0
  "SP1", // spare 1
  "IRP",// 40
  "IDS",
  "DSC",
  "SLK",
  "ISL",
  "FOM", // 45
  "OFR",
  "BUG",
  "LIB",
  "SP2", // spare 2
  "EOF"  // 50
};

  /*Generic message object handle */

typedef struct { 
  void         *m;          /* A pointer to a message. */
  MessageType  t;          /* The message type discriminator. */
  int32        s;          /* The message size in bytes. */
} GenericMesg;



/* BAT && IBA record */
typedef struct InternalBatchMesgTag {
  char         *name;
  time_t       created;
  Batch_ID     eaccession;
  char         *comment;
  IntBatch_ID  iaccession;
}InternalBatchMesg;

typedef InternalBatchMesg BatchMesg;

  /* ADL record */

typedef struct AuditLineTag {
  struct AuditLineTag  *next;
  char                 *name;
  time_t               complete;
  char                 *version;
  char                 *comment;
} AuditLine;

  /* ADT message */

typedef struct {
  AuditLine *list;
} AuditMesg;


/* BIN message */
typedef struct BinMesgTag{
  ActionType      action;
  Bin_ID          eaccession;
  time_t          entry_time;
  char            *source;  
  IntBin_ID       iaccession;
}BinMesg;

typedef BinMesg InternalBinMesg;



  /* LKG message */

typedef enum {
  AS_MATE       = (int)'M', // Mate
  AS_STS_GUIDE  = (int)'S', // STS
  AS_BAC_GUIDE  = (int)'B', // BAC
  AS_REREAD     = (int)'R', // Reread
  AS_MAY_JOIN   = (int)'Y', // maY
  AS_MUST_JOIN  = (int)'T',  // musT
  AS_B_MATE =     (int)'G'   // BGLii mate
} LinkType;

typedef enum {
  AS_NORMAL     = (int)'N',
  AS_INNIE      = (int)'I', 
  AS_OUTTIE     = (int)'O',
  AS_ANTI	= (int)'A',
  AS_UNKNOWN	= (int)'U'
} OrientType;

typedef struct {
  ActionType      action;
  LinkType        type;
  time_t          entry_time;
  OrientType      link_orient;
  Fragment_ID     frag1;
  Fragment_ID     frag2;
  Distance_ID     distance;
} LinkMesg;

  /* ILK message */
 
typedef struct {
  ActionType      action;
  LinkType        type;
  time_t          entry_time;
  OrientType      link_orient;
  IntFragment_ID  ifrag1;
  IntFragment_ID  ifrag2;
  IntDist_ID      idistance;
} InternalLinkMesg;

  /* DST message */

typedef struct {
  ActionType   action;
  Distance_ID  eaccession;
  float32      mean;
  float32      stddev;
} DistanceMesg;

  /* IDT message */

typedef struct {
  ActionType   action;
  Distance_ID  eaccession;
  float32      mean;
  float32      stddev;
  IntDist_ID   iaccession;
} InternalDistMesg;

  /* SCN message */

typedef enum {
  AS_UBIQREP       = (int)'U', 
  AS_CONTAMINANT   = (int)'C',
  AS_NO_SCREENTYPE = (int) 'N'
} ScreenType;

/* bit-vector values for relevance field */
#define AS_OVL_HEED_RPT     (1)
#define AS_URT_IS_VECTOR    (2)
#define AS_URT_IS_SATELLITE (4)
#define AS_URT_IS_SIMPLE    (8)

typedef struct {
  ActionType        action;
  ScreenType        type;
  ScreenItem_ID     eaccession;
  Repeat_ID         erepeat_id;
  int32             relevance;
  char              *source;
  char              *sequence;
  float             variation;
  CDS_COORD_t    min_length;
  IntScreenItem_ID  iaccession;
  IntRepeat_ID      irepeat_id;
} InternalScreenItemMesg;

#if 0
  /* ISN message */

typedef struct {
  ActionType 	  action;
  ScreenType      type;
  ScreenItem_ID   eaccession;
  Repeat_ID       repeat_id;
  int32           relevance;
  char            *source;
  char            *sequence;
  float           variation;
  CDS_COORD_t  min_length;
  IntScreen_ID    iaccession;
} InternalScreenItemMesg;
#endif

typedef InternalScreenItemMesg ScreenItemMesg;

  /* RPT message */

typedef struct {
  Repeat_ID       erepeat_id;
  char            *which;
  CDS_COORD_t  length;
  IntRepeat_ID    irepeat_id;
} InternalRepeatItemMesg;


typedef InternalRepeatItemMesg RepeatItemMesg;

typedef enum {
  AS_READ    = (int)'R',  //Celera Read
  AS_EXTR    = (int)'X',  //External WGS read
  AS_TRNR    = (int)'T',  //Transposon library read
  AS_EBAC    = (int)'E',  //End of BAC
  AS_LBAC    = (int)'L',  //Lightly shotgunned
  AS_UBAC    = (int)'U',  //Unfinished
  AS_FBAC    = (int)'F',  //Finished
  AS_STS     = (int)'S',  //Sts
  AS_UNITIG  = (int)'u',  //Unitig
  AS_CONTIG  = (int)'c',   //Contig
  AS_BACTIG  = (int) 'B',   // BacTig
  AS_FULLBAC = (int)'C',   // Full Bac C = Complete)
  AS_B_READ  = (int)'G' // BGLII read
} FragType;

typedef enum {
  AS_UNIQUE_UNITIG   = (int)'U',  // U-Unique
  AS_ROCK_UNITIG     = (int)'R',  // Rock
  AS_STONE_UNITIG    = (int)'S',  // Stone
  AS_PEBBLE_UNITIG   = (int)'P',  // Pebble
  AS_SINGLE_UNITIG   = (int)'s', // Singleton Unitig in Unplaced Contig
  AS_OTHER_UNITIG    = (int)'X'  // Unspecified surrogate unitig
} UnitigType;

  /* ISM & SMA record */

typedef enum {
 AS_FORWARD   = (int)'F',
 AS_REVERSE = (int)'R'
} DirectionType;

typedef struct iScreenMatchTag {
  struct iScreenMatchTag  *next;
  SeqInterval             where;
  CDS_IID_t            iwhat;
  IntRepeat_ID            repeat_id;
  int32                   relevance;
  SeqInterval             portion_of;
  DirectionType           direction;
} IntScreenMatch;

typedef struct ScreenMatchTag {
  struct ScreenMatchTag  *next;
  SeqInterval            where;
  CDS_UID_t           what;
  Repeat_ID              repeat_id;
  int32                  relevance;
  SeqInterval            portion_of;
  DirectionType          direction;
} ScreenMatch;

  /* SFG message */


typedef struct {
  ActionType   		action;
  Fragment_ID  		eaccession;
  FragType     		type;
  Locale_ID    		elocale;
  Sequence_ID  		eseq_id;
  Bactig_ID  		ebactig_id;
  SeqInterval  		locale_pos;
  time_t       		entry_time;
  SeqInterval  		clear_rng;
  char        		*source;
  char        		*sequence;
  char        		*quality;
  IntFragment_ID   	iaccession;
  IntScreenMatch	*screened;
  IntLocale_ID          ilocale;
  IntSequence_ID          iseq_id;
  IntBactig_ID          ibactig_id;
} ScreenedFragMesg;

/* OFG, FRG, IFG, OFR messages */

typedef ScreenedFragMesg OFGMesg;
typedef ScreenedFragMesg FragMesg;
typedef ScreenedFragMesg InternalFragMesg;
typedef ScreenedFragMesg OFRMesg;

  /*OVL message*/

#define AS_LONGEST_DELTA     (126)
#define AS_LONG_DELTA_CODE  (-127)
#define AS_POLY_DELTA_CODE  (-128)
#define AS_ENDOF_DELTA_CODE    (0)

typedef enum {
  AS_DOVETAIL    = (int)'D',
  AS_CONTAINMENT = (int)'C', 
  AS_SUPERREPEAT = (int)'S',
  
  AS_DOVETAIL_TRANS = (int) 'X',
  // A dovetail overlap that is transitively inferrable by a path in
  // the fragment overlap graph, for example using a path of dovetail
  // overlaps followed by a path of containment overlaps.

  /* The following are created during de-chording the fragment overlap
     graph. They are not by default present in a path compressed
     fragment overlap graph but instead removed since they are easily
     inferred. */
  AS_DOVETAIL_CHORD  = (int)'d',
  // This dovetail overlap is known to be a chord in the dovetail
  // overlap sub-graph.
  AS_CONTAINMENT_CHORD = (int)'c'
  // This containment overlap is known to be a chord in the
  // containment overlap sub-graph.

} OverlapType;


typedef struct {
  IntFragment_ID   aifrag, bifrag;
  CDS_COORD_t   ahg,bhg;
  OrientType       orientation;
  OverlapType      overlap_type;
  float32          quality;
  CDS_COORD_t   min_offset, max_offset;
  int32            polymorph_ct;
  signed char      *delta;
} OverlapMesg;

  /* BRC message */

typedef struct {
  ActionType        action;
  IntFragment_ID    ifrag;
  CDS_COORD_t    pre_br, suf_br, pre_end, suf_end;
} BranchMesg;


typedef enum {
  AS_A_END = (int)'A',
  AS_B_END = (int)'B'
} ChunkOrientType;

typedef enum {
  AS_SINGLETON    = (int)'S',
  AS_INTERCHUNK_A = (int)'A',
  AS_INTERCHUNK_B = (int)'B',
  AS_INTRACHUNK   = (int)'I',
  AS_CONTCHUNK    = (int)'C'
} LabelType;

typedef enum {
  AS_INTO_REPEAT = (int) 'R',
  AS_INTO_UNIQUE = (int) 'U',
  AS_NO_BPOINT   = (int) 'N'
} BranchType;

typedef struct  {
  IntFragment_ID  ifrag;
  CDS_COORD_t  offset5p;
  CDS_COORD_t  offset3p;
  FragType        type;
  LabelType       label;
  char            *source;
} ChunkFrag;

typedef struct {
  IntChunk_ID     chunk;
  ChunkOrientType orient;
  CDS_COORD_t  best_overlap_length;
  CDS_COORD_t  min_overlap_length;
  CDS_COORD_t  max_overlap_length;
} ChunkOverlap;

/* UOM message */

/*
 * AS_NORMAL		AB_AB
 * AS_ANTINORMAL	BA_BA
 * AS_OUTIE		BA_AB
 * AS_INNIE		AB_BA
 *
 */

// ChunkOrientationType discontinued by Jason 7/01 
// because code intermingled it with OrientType.
// Both enums used the same integer values.
//typedef enum {
//  AB_AB		= (int) 'N',
//  BA_BA		= (int) 'A',
//  BA_AB		= (int) 'O',
//  AB_BA		= (int) 'I',
//  XX_XX         = (int) 'U'    // unknown relative orientation
//} ChunkOrientationType;
#define ChunkOrientationType OrientType
#define AB_AB AS_NORMAL  
#define BA_BA AS_ANTI  
#define BA_AB AS_OUTTIE  
#define AB_BA AS_INNIE 
#define XX_XX AS_UNKNOWN 

typedef enum {
  /* The following are CGW overlap output classifications: */
  AS_NO_OVERLAP	= (int) 'N', // Not used by the unitigger.
  AS_OVERLAP    = (int) 'O', // A dovetail overlap between unitigs.
  AS_TANDEM_OVERLAP = (int) 'T', // Not used by the unitigger.
    
  
  /* The following are unitigger overlap output classifications: */
  AS_1_CONTAINS_2_OVERLAP = (int) 'C', 
  AS_2_CONTAINS_1_OVERLAP = (int) 'I',
  // Two types of containment overlap types.

  AS_TOUCHES_CONTAINED_OVERLAP   = (int) 'M',
  // A dovetail overlap to a singleton unitig composed of a contained
  // fragment from a non-contained unitig. This overlaps touch
  // multiply contained fragments and orphan contained fragments. This
  // is another kind of overlap necessary when processing containment
  // overlaps.

  AS_TRANSCHUNK_OVERLAP   = (int) 'X',
  // A dovetail overlap between unitigs that was marked for removal by
  // because it was transitively inferrable in the fragment overlap
  // graph but was not a chord in the dovetail overlap sub-graph.
  
  AS_DOVETAIL_CHORD_OVERLAP  = (int)'d',
  // This dovetail overlap is known to be a chord in the dovetail
  // overlap sub-graph.  This clasification is created during
  // de-chording the fragment overlap graph. They are by default not
  // present in a path compressed fragment overlap graph but instead
  // removed since they are easily inferred.

  AS_1_CONTAINS_2_CHORD_OVERLAP = (int)'c',
  // This containment overlap is known to be a chord in the
  // containment overlap sub-graph.  This clasification is created
  // during de-chording the fragment overlap graph. They are by
  // default not present in a path compressed fragment overlap graph
  // but instead removed since they are easily inferred.

  AS_2_CONTAINS_1_CHORD_OVERLAP = (int)'i',
  // This containment overlap is known to be a chord in the
  // containment overlap sub-graph.  This clasification is created
  // during de-chording the fragment overlap graph. They are by
  // default not present in a path compressed fragment overlap graph
  // but instead removed since they are easily inferred.

  AS_BETWEEN_CONTAINED_OVERLAP   = (int) 'Y',
  // A dovetail overlap between unitigs each spanned by a contained
  // fragment. Thus, this overlap is between multiply contained
  // fragments and/or orphan contained fragments.

  AS_1_CONTAINS_2_STACK_OVERLAP   = (int) 'Z'
  // A containment overlap between two globally contained fragments
  // that is not a chord in the containment overlap graph.

} UnitigOverlapType;

typedef struct {
  IntChunk_ID		chunk1;
  IntChunk_ID     	chunk2;
  ChunkOrientationType	orient;
  UnitigOverlapType	overlap_type;
# ifdef AS_ENABLE_SOURCE
  char			*source;
# endif
  CDS_COORD_t        best_overlap_length;
  CDS_COORD_t        min_overlap_length;
  CDS_COORD_t        max_overlap_length;
  float32               quality;
} UnitigOverlapMesg;

/* FOM */
// This is the alternative overlap message.  The FOM message includes
// the same information as the OVL message in the absence of alignment
// deltas.
typedef struct {
  IntFragment_ID	afrag;
  IntFragment_ID     	bfrag;
  ChunkOrientationType	orient;
  UnitigOverlapType	overlap_type;
# ifdef AS_ENABLE_SOURCE
  char			*source;
# endif
  CDS_COORD_t        best_overlap_length;
  CDS_COORD_t        min_overlap_length;
  CDS_COORD_t        max_overlap_length;
  float32               quality;
} FragOverlapMesg;

typedef struct {
  IntChunk_ID     iaccession;
  CDS_COORD_t  bp_length;
  float           coverage_stat;
  BranchType      a_branch_type;
  BranchType      b_branch_type;
  CDS_COORD_t  a_branch_point;
  CDS_COORD_t  b_branch_point;
  int32           num_frags;
  int32           a_degree;
  int32           b_degree;
  ChunkFrag      *f_list;
  ChunkOverlap   *a_list;
  ChunkOverlap   *b_list;
  char           *source;
} ChunkMesg;

  /* ICT messages */

typedef enum {
  AS_KEEP   = (int) 'K',
  AS_NOKEEP = (int) 'N'
} ResolveType;

typedef struct {
  FragType             type;
  IntFragment_ID       ident;
  ResolveType          label;
  DirectionType        orientation;
  CDS_COORD_t       position;
} LayoutPos;



typedef enum {
  AS_UNIQUE =     (int)'U',
  AS_CHIMER =     (int)'C',
  AS_NOTREZ =     (int)'N',
  AS_SEP =        (int)'S',
  AS_UNASSIGNED = (int) 'X'
} UnitigStatus;

/* This will eventually be Replaced by IntMultiPos */
typedef struct MultiPos {
  FragType        type;
  IntFragment_ID  ident;
  SeqInterval     position;
  int32           delta_length;
  int32           *delta;
} MultiPos;

/* IMP message */

typedef struct IntMultiPos {
  FragType        type;
  IntFragment_ID  ident;

# ifdef AS_ENABLE_SOURCE
  char		  *source;
#ifdef i386
  int32 ptrPad1;
#endif
# endif

  SeqInterval     position;

  IntFragment_ID       contained;
  int32           delta_length;

  int32           *delta;
#ifdef i386
  int32 ptrPad2;
#endif

} IntMultiPos;

/* IMV message */

typedef struct IntMultiVar {
  IntFragment_ID  ident;
#ifdef i386
  int32           ptrPad1;
#endif
  SeqInterval     position;
  int32           nreads;
  int32           nreads_best;
  float           ratio;     
  int32           var_length;
  char           *var_sequence;
  int32           window;  
  IntFragment_ID *aindent;
#ifdef i386
  int32           ptrPad2;
#endif
} IntMultiVar;

/* This is a variant of IntMultiPos to handle deltas in a longer (unitig) sequence */
typedef struct {
  UnitigType    type;
  IntUnitig_ID  ident;
  SeqInterval   position;
  int32         delta_length;
  int32                 wordPad;
  int32         *delta;
#ifdef i386
  int32 ptrPad;
#endif
} IntUnitigPos;


typedef struct {
  UnitigType   type;
  SeqInterval  position;
  int32        delta_length;
  int32        *delta;
  Unitig_ID    eident;
} UnitigPos;

/* IEP messages */

typedef struct {
  FragType        type;
  IntFragment_ID  ident;
  SeqInterval     position;
} IntElementPos;

/* IUM */

typedef struct {
  IntChunk_ID     iaccession;
# ifdef AS_ENABLE_SOURCE
  char		  *source;
# endif
  float32         coverage_stat;
  UnitigStatus    status;
  CDS_COORD_t  a_branch_point;
  CDS_COORD_t  b_branch_point;
  CDS_COORD_t  length;
  char            *consensus;
  char            *quality;
  int32		  forced;
  int32           num_frags;
  IntMultiPos    *f_list;
  int32           num_vars; 
  IntMultiVar    *v_list;
} IntUnitigMesg;


/* The following message type will eventually be Removed */
typedef struct {
  Chunk_ID        eaccession;
  UnitigStatus    status;
  int32           num_occurences;
  Contig_ID       *occurences;
  CDS_COORD_t  length;
  char            *consensus;
  char            *quality;
  int32           num_reads;
  int32           num_guides;
  int32           sum_delta_lengths;
  MultiPos        *reads;
  MultiPos        *guides;
} UnitigMesg;


typedef enum {
  AS_PLACED	= (int) 'P',
  AS_UNPLACED   = (int) 'U'
} ContigPlacementStatusType;

/* ICM */

typedef struct {
  IntContig_ID               iaccession;
  ContigPlacementStatusType  placed;
  CDS_COORD_t             length;
  char                       *consensus;
  char                       *quality;
  int32		             forced;
  int32                      num_pieces;
  int32                      num_unitigs;
  int32                      num_vars;
  IntMultiPos                *pieces;
  IntUnitigPos               *unitigs;
  IntMultiVar                *v_list;
} IntConConMesg;


/* IUL Mesage */

typedef struct {
  IntFragment_ID  in1, in2; 
  LinkType        type;
} IntMate_Pairs;


typedef enum {
  AS_IN_ASSEMBLY	= (int) 'A',
  AS_POLYMORPHIRSM	= (int) 'P',
  AS_BAD		= (int) 'B',
  AS_CHIMERA		= (int) 'C',
  AS_UNKNOWN_IN_ASSEMBLY= (int) 'U'
} PlacementStatusType;

typedef struct {
  IntChunk_ID		unitig1;
  IntChunk_ID		unitig2;
  ChunkOrientationType	orientation;
  UnitigOverlapType	overlap_type;
  int32			is_possible_chimera;
  int32			includes_guide;
  float32		mean_distance;
  float32		std_deviation;
  int32			num_contributing;
  PlacementStatusType	status;
  IntMate_Pairs		*jump_list;
} IntUnitigLinkMesg;


/* ICL message */

typedef struct {
  IntChunk_ID		contig1;
  IntChunk_ID		contig2;
  ChunkOrientationType	orientation;
  UnitigOverlapType	overlap_type;
  int32			is_possible_chimera;
  int32			includes_guide;
  float32		mean_distance;
  float32		std_deviation;
  int32			num_contributing;
  PlacementStatusType	status;
  IntMate_Pairs		*jump_list;
} IntContigLinkMesg;

/* ISL message */

typedef struct {
  IntScaffold_ID	iscaffold1;
  IntScaffold_ID	iscaffold2;
  ChunkOrientationType	orientation;
  int32			includes_guide;
  float32		mean_distance;
  float32		std_deviation;
  int32			num_contributing;
  IntMate_Pairs		*jump_list;
} InternalScaffoldLinkMesg;


/* IAF message */

typedef enum {
  GOOD_MATE		= (int) 'G',
  BAD_MATE		= (int) 'B',
  NO_MATE		= (int) 'N',
  UNRESOLVED_MATE	= (int) 'U'
} MateStatType;

typedef struct {
  IntFragment_ID  iaccession;
  FragType        type;
  int32           chimeric;
  int32           chaff;
  SeqInterval     clear_rng;
  MateStatType    mate_status;
} IntAugFragMesg;

/* AFG message */

typedef struct {
  Fragment_ID     eaccession;
  IntFragment_ID  iaccession;              
  ScreenMatch	  *screened;
  MateStatType    mate_status;
  int32           chimeric;
  int32           chaff;
  SeqInterval     clear_rng;
} AugFragMesg;

/* IMD message */

typedef struct {
  IntDist_ID  refines;
  float32     mean;
  float32     stddev;
  int32	      min;
  int32	      max;
  int32	      num_buckets;
  int32       *histogram;
} IntMateDistMesg;

/* ISF message */

typedef struct {
  IntContig_ID		contig1;
  IntContig_ID		contig2;
  float32		mean;
  float32		stddev;
  ChunkOrientationType	orient;
} IntContigPairs;

typedef struct {
  IntScaffold_ID  iaccession;
  int32		  num_contig_pairs;
  IntContigPairs  *contig_pairs;
} IntScaffoldMesg;

/* IDS message */

typedef struct {
  IntContig_ID icontig;
} IntDegenerateScaffoldMesg;


/* Genome Snapshot typedefs */
/****************************/

/* MPS message */
typedef struct {
  FragType      type;
  Fragment_ID   eident;
# ifdef AS_ENABLE_SOURCE
  char		*source;
# endif
  SeqInterval   position;
  int32         delta_length;
  int32         *delta;
} SnapMultiPos;

/* EPS messages */
typedef struct {
  FragType     type;
  Fragment_ID  eident;
  SeqInterval  position;
} SnapElementPos;


/* UTG Message */
typedef struct {
  Chunk_ID        eaccession;  // changed in comparison to internal message
  IntChunk_ID     iaccession;
# ifdef AS_ENABLE_SOURCE
  char		  *source;
# endif
  float32         coverage_stat;
  UnitigStatus    status;
  CDS_COORD_t  a_branch_point;
  CDS_COORD_t  b_branch_point;
  CDS_COORD_t  length;
  char            *consensus;
  char            *quality;
  int32		  forced;
  int32           num_frags;
  SnapMultiPos    *f_list;// changed in comparison to internal message
} SnapUnitigMesg;


typedef struct {
  Fragment_ID  in1, in2; 
  LinkType     type;
} SnapMate_Pairs;

/* ULK message */
typedef struct {
  Chunk_ID		eunitig1;
  Chunk_ID		eunitig2;
  ChunkOrientationType	orientation;
  UnitigOverlapType	overlap_type;
  int32			is_possible_chimera;
  int32			includes_guide;
  float32		mean_distance;
  float32		std_deviation;
  int32			num_contributing;
  PlacementStatusType	status;
  SnapMate_Pairs	*jump_list; // changed in comparison to internal message
} SnapUnitigLinkMesg;


/* CCO message */
typedef struct {
  Contig_ID                   eaccession;
  IntContig_ID                iaccession;
  ContigPlacementStatusType   placed;
  CDS_COORD_t              length;
  char                        *consensus;
  char                        *quality;
  int32                       forced;
  int32                       num_pieces;
  int32                       num_unitigs;
  SnapMultiPos                *pieces; // changed in comparison to internal message
  UnitigPos                   *unitigs;// changed in comparison to internal message
} SnapConConMesg;


/* CLK message */
typedef struct {
  Chunk_ID		econtig1; // changed in comparison to internal message
  Chunk_ID		econtig2; // changed in comparison to internal message
  ChunkOrientationType	orientation;
  UnitigOverlapType	overlap_type;
  int32			is_possible_chimera;
  int32			includes_guide;
  float32		mean_distance;
  float32		std_deviation;
  int32			num_contributing;
  PlacementStatusType	status;
  SnapMate_Pairs	*jump_list; // changed in comparison to internal message
} SnapContigLinkMesg;

/* SLK message */
typedef struct {
  Scaffold_ID           escaffold1;
  Scaffold_ID           escaffold2;
  ChunkOrientationType	orientation;
  int32			includes_guide;
  float32		mean_distance;
  float32		std_deviation;
  int32			num_contributing;
  SnapMate_Pairs	*jump_list;
} SnapScaffoldLinkMesg;


/* SCF and CTP message */
typedef struct {
  Contig_ID		econtig1; // changed in comparison to internal message
  Contig_ID		econtig2; // changed in comparison to internal message
  float32		mean;
  float32		stddev;
  ChunkOrientationType	orient;
} SnapContigPairs;

typedef struct {
  Scaffold_ID           eaccession;
  IntScaffold_ID        iaccession;
  int32			num_contig_pairs;
  SnapContigPairs   	*contig_pairs; // changed in comparison to internal message
} SnapScaffoldMesg;


/* DSC message */
typedef struct {
  Scaffold_ID      eaccession;
  Contig_ID        econtig;
} SnapDegenerateScaffoldMesg;



/* MDI message */
typedef struct {
  Dist_ID		erefines; // changed in comparison to internal message
  IntDist_ID		irefines; // changed in comparison to internal message
  float			mean;
  float			stddev;
  int32			min;
  int32			max;
  int32			num_buckets;
  int32			*histogram;
} SnapMateDistMesg;



/* BTG & IBT subrecord */
typedef struct InternalBactigMesgTag{
  Bactig_ID       eaccession;
  CDS_COORD_t  length;
  IntBactig_ID    iaccession;
}InternalBactigMesg;

typedef InternalBactigMesg BactigMesg;

typedef enum {
  AS_ENDS            = (int)'E',  //End of BAC
  AS_LIGHT_SHOTGUN   = (int)'L',  //Lightly shotgunned
  AS_UNFINISHED      = (int)'U',  //Unfinished
  AS_FINISHED        = (int)'F'   //Finished
} BACType;


/* BAC && IBC record */
typedef struct InternalBacMesgTag{
  ActionType          action;
  Bac_ID              ebac_id;
  Sequence_ID         eseq_id;
  BACType             type;
  time_t              entry_time;
  int16               num_bactigs;
  InternalBactigMesg  *bactig_list;  // array of length num_bactigs
  Distance_ID         elength;
  char                *source;
  IntBac_ID           ibac_id;
  IntSequence_ID      iseq_id;
  IntDistance_ID      ilength;
}InternalBacMesg;

typedef InternalBacMesg BacMesg;

/* LIB */
typedef struct LibDonorMesgTag {
  ActionType  action;
  Library_ID  eaccession;
  Donor_ID    donor;
  char        *source;
} LibDonorMesg;

/* WEL */
typedef struct WellMesgTag{
  Fragment_ID  efrag;
  Well_ID      ewell;  // NOTE this external ID is uint16
  Library_ID   elibrary;
}WellMesg;

/* PLA */
typedef struct PlateMesgTag{
  ActionType  action;
  Plate_ID    eaccession;
  int16       num_wells;
  WellMesg    *well_list; // array of length num_wells
}PlateMesg;

/* LKP */
typedef struct LinkPlateMesgTag{
  ActionType  action;
  Plate_ID    eplate_for;
  Plate_ID    eplate_rev;
}LinkPlateMesg;

/* BSP */
typedef struct BugSplitPosTag {
  CDS_COORD_t  position;
  Fragment_ID     bactig_eaccession;
} BugSplitPos;


/* BUG */
typedef enum {
  AS_CREATE_BAC =       (int) 'C',
  AS_REMOVE_BACTIG =    (int) 'R',
  AS_ADD_BACTIG =       (int) 'A',
  AS_SPLIT_BACTIG =     (int) 'S',
  AS_CHANGE_CLR_RANGE = (int) 'G'
} BugMesgType;    

typedef struct BugMesgTag {
  BugMesgType  type;
  char         *source;
  Locale_ID    bac_eaccession;
  Locale_ID    seq_eaccession;
  Locale_ID    src_bac_eaccession;
  Locale_ID    src_seq_eaccession;
  Locale_ID    dst_bac_eaccession;
  Locale_ID    dst_seq_eaccession;
  Locale_ID    bactig_eaccession;
  Locale_ID    new_bactig_eaccession;
  Locale_ID    old_frag_eaccession;
  Locale_ID    new_frag_eaccession;
  int32        num_pos;
  BugSplitPos  *split_array;
  SeqInterval  clear_rng;
} BugMesg;

/* EOF */
typedef struct EndOfFileMesgTag {
  int32   status;
  time_t  created;
  char    *comment;
} EndOfFileMesg;

typedef enum {
  AS_BINARY_OUTPUT  = (int) 'B',
  AS_PROTO_OUTPUT   = (int) 'P'
} OutputType;

/* External Routines */

#undef TODELETE
#ifdef TODELETE
#endif

typedef int (*MesgReader)(FILE *, GenericMesg **);
typedef int (*MesgWriter)(FILE *, GenericMesg *);

#if 0

/* Function: DuplicateProtoMesg_AS 

   Description: Copies the generic message structure including all of
   its variable length data. The variable length data must be released
   by free_mesg().
   
   Return Value: A pointer to a copy of the generic message that has a
   lifetime determined by the client.

   Inputs: omesg - A pointer to a generic message.

*/

extern GenericMesg *DuplicateProtoMesg_AS(GenericMesg *omesg);

/* Function: FreeProtoMesg_AS 

   Description: Frees the memory for a generic message structure.

   Input/Outputs: omesg - A pointer to a generic message.

*/
   
extern void FreeProtoMesg_AS(GenericMesg *omesg);

/* Function: InputFileType_AS

   Description: The file cursor for fin must be at the start of the file
   when this routine is called.  The routine examines the first byte of
   the file and decides whether it is a binary or proto file.  It then
   returns a pointer to either ReadProtoMesg_AS or ReadBinaryMesg_AS
   accordingly.
*/

extern MesgReader InputFileType_AS(FILE *fin);


/* Function: OutputFileType_AS
   Description: This routine
   returns a pointer to either WriteProtoMesg_AS or WriteBinaryMesg_AS
   as a function of its input.
*/

extern MesgWriter OutputFileType_AS(OutputType type);




/* Functions: ReadProtoMesg_AS 
              ReadBinaryMesg_AS

   Description: Reads the next message from the file "fin" and a
   returns the memory location of a generic message.  This memory is
   managed by the routine only. This memory location and all of its
   variable data will go out of scope or be trashed by the next call
   to the routine.  The first routine is for reading proto ASCI input
   and the second for reading binary encoded input.
   
   Return Value: The return value is EOF if an end of file is
   encountered occurs, otherwise the return value is zero indicating
   success.
   
   Outputs: mesg - A handle to a generic message.
   
   Input/Outputs: fin - A file openned for text reading. 
*/

#endif
extern int ReadProtoMesg_AS(FILE *fin, GenericMesg **pmesg);
extern int ReadBinaryMesg_AS(FILE *fin, GenericMesg **pmesg);
#if 0

/* Functions:  WriteProtoMesg_AS 
               WriteBinaryMesg_AS

   Description: Writes a generic message to the file "fout" in either
   ASCII or binary mode depending on the routine. 

   Return Value: The return value is negative if an error occured.

   Inputs: mesg - A pointer to the generic message to be output. 

   Input/Outputs: fout - A file openned for text writing.
*/

#endif
extern int WriteProtoMesg_AS(FILE *fout, GenericMesg *mesg);
extern int WriteBinaryMesg_AS(FILE *fout, GenericMesg *mesg);
#if 0
/* Function: Transfer_XXX_to_YYY

   Description: Transfers the fields of an XXX message to that of a
   YYY message.  Note carefully that second level memory is not duplicated,
   so that such items (e.g. sequence) are *shared* between the two structures.
   A kludge to allow items to be throughput without assignment of internal
   IDs, should no longer be in use!
*/

extern void Transfer_FRG_to_IFG_AS(FragMesg         *frg_mesg,
                                   InternalFragMesg *ifg_mesg);

extern void Transfer_IFG_to_SFG_AS(InternalFragMesg *ifg_mesg,
                                   ScreenedFragMesg *sfg_mesg);

extern void Transfer_SFG_to_OFG_AS(ScreenedFragMesg *sfg_mesg,
                                   OFGMesg *ofg_mesg);

extern void Transfer_SFG_to_OFR_AS(ScreenedFragMesg *sfg_mesg,
                                   OFRMesg *ofr_mesg);

extern void Transfer_DST_to_IDT_AS(DistanceMesg     *dst_mesg,
                                   InternalDistMesg *idt_mesg);

extern void Transfer_LKG_to_ILK_AS(LinkMesg         *lkg_mesg,
                                   InternalLinkMesg *ilk_mesg);

extern void Transfer_SCN_to_ISN_AS(ScreenItemMesg         *scn_mesg,
                                   InternalScreenItemMesg *isn_mesg);

extern void AppendAuditLine_AS(AuditMesg *adt_mesg,
                               AuditLine *auditLine,
                               time_t t, char *name,
                               char *version, char *comment);

/* Function: GetProtoLineNum_AS

   Description: When reading in proto mode, this function will return
   the line number the input is currently on.  This function operates
   correctly only when a single input is being read.  The return value
   in all other cases is the sum of the number of lines read in all
   proto files thus far.
*/

extern int GetProtoLineNum_AS(void);

/* GetMessageType:
   Returns a number in the range [1, NUM_OF_REC_TYPES -1]
   as a function of the first 3 characters of the passed string.
*/
int GetMessageType(char *string);

/* GetMessageName:
   Returns a string as a function of message type
*/
const char  *GetMessageName(int type);

/* Free all memory allocated by the Proto IO Package */

extern void ResetProto_AS(void);
extern void ResetBinary_AS(void);

#endif
#endif /* AS_MSG_PMESG_INCLUDE */

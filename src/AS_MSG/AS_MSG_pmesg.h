
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
/* $Id: AS_MSG_pmesg.h,v 1.30 2007-01-29 20:41:13 brianwalenz Exp $   */

#ifndef AS_MSG_PMESG_INCLUDE
#define AS_MSG_PMESG_INCLUDE

#include <stdio.h>
#include <time.h>

#include "AS_global.h"

/* Exported constant definitions; Macro definitions; type definitions */

// Defining the following enables internal source fields for testing
#define AS_ENABLE_SOURCE

/* All externally created accession numbers are 64 bits. */

#define AS_BACTIG_MAX_LEN (500000)
#define AS_BACTIG_MIN_LEN (150)

#define AS_BAC_MAX_LEN (500000)
#define AS_BAC_MIN_LEN (150)

extern int novar;

#ifdef GENERIC_STORE_USE_LONG_STRINGS
#define AS_READ_MAX_LEN AS_BACTIG_MAX_LEN
#else
#define AS_READ_MAX_LEN AS_FRAG_MAX_LEN
#endif
#define AS_READ_MIN_LEN AS_FRAG_MIN_LEN

//   #define NEW_UNITIGGER_INTERFACE

#define DEFINE_IDs(type)\
typedef CDS_UID_t type##_ID;\
typedef CDS_IID_t Int##type##_ID;

DEFINE_IDs(Fragment);
DEFINE_IDs(Locale);
DEFINE_IDs(Distance);
DEFINE_IDs(Chunk);
DEFINE_IDs(Unitig);
DEFINE_IDs(Contig);
DEFINE_IDs(Dist);
DEFINE_IDs(Scaffold);
DEFINE_IDs(Batch);
DEFINE_IDs(Bactig);
DEFINE_IDs(Bac);
DEFINE_IDs(Sequence);

typedef enum {
  AS_ADD      = (int)'A',
  AS_DELETE   = (int)'D',
  AS_UPDATE   = (int)'U',
  AS_REDEFINE = (int)'R'
} ActionType;

typedef struct {
  CDS_COORD_t bgn;
  CDS_COORD_t end;
} SeqInterval;


typedef enum {
  MESG_NONE = 0,
  MESG_ADT,
  MESG_FRG,
  MESG_IFG,
  MESG_SPe,
  MESG_OFG, // 5
  MESG_LKG,
  MESG_ILK,
  MESG_DST,
  MESG_IDT,
  MESG_SPb, // 10
  MESG_SPc,
  MESG_SP1,
  MESG_OVL,
  MESG_SPf,
  MESG_UOM, // 15
  MESG_IUM,
  MESG_IUL,
  MESG_ICL,
  MESG_AFG,
  MESG_ISF, // 20
  MESG_IMD,
  MESG_IAF,
  MESG_UTG,
  MESG_ULK,
  MESG_ICM, // 25
  MESG_CCO,
  MESG_CLK,
  MESG_SCF,
  MESG_MDI,
  MESG_BAT, // 30  
  MESG_IBA,
  MESG_BAC,
  MESG_IBC,
  MESG_SP2,
  MESG_IBI, // 35
  MESG_SP3,
  MESG_SP4,
  MESG_SP5,
  MESG_SP6,
  MESG_SP7, // 40
  MESG_IDS,
  MESG_DSC,
  MESG_SLK,
  MESG_ISL,
  MESG_FOM, // 45
  MESG_SPd,
  MESG_SP8,
  MESG_SP9,
  MESG_SPa,
  MESG_EOF  // 50
} MessageType;


#define NUM_OF_REC_TYPES MESG_EOF

static char  *MessageTypeName[NUM_OF_REC_TYPES + 1] = {
  "NONE",
  "ADT",
  "FRG",
  "IFG",
  "SPe",
  "OFG", // 5
  "LKG",
  "ILK",
  "DST",
  "IDT",
  "SPb", // 10
  "SPc",
  "SP1", 
  "OVL", 
  "SPf",
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
  "BAT",// 30  
  "IBA",
  "BAC",
  "IBC", 
  "SP2",
  "IBI",// 35
  "SP3", 
  "SP4",
  "SP5",
  "SP6",
  "SP7",// 40
  "IDS",
  "DSC",
  "SLK",
  "ISL",
  "FOM", // 45
  "SPd",
  "SP8",
  "SP9",
  "SPa",
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


/* Extended Granger's comment to fragment types

AS_READ    Celera Read
This should be a .trusted. randomly generated shotgun read from the entire
DNA target sequence . usually a whole genome. Features: randomly sampled, trusted for consensus,
can have a mate pair.

AS_EXTR    External WGS read
This should be an .untrusted. randomly generated shotgun read . same as above. Features:
randomly sampled, untrusted for consensus, can have a untrusted mate pair

AS_TRNR    Transposon library read
Read was generated from a subclone of the target sequence using transposon .bombing..
Features: nonrandom, might include information about subclone, trusted for consensus, can have
mate pair but oriented in .outtie. rather .innie..

AS_EBAC    End of BAC
Read is a BAC end . that is a sequence from one end or the other of a BAC subclone.
Assumed to be randomly sampled. Variance is assumed to be high for size of insert so not
used in initial scaffolding. Features: random sampling, high variance insert size, trusted for consensus,
can have mate pair.

AS_LBAC    Lightly shotgunned
Read is shotgun sampled from a subclone (BAC or otherwise). Features: nonrandom, trusted,
can have a mate pair.

AS_UBAC    Unfinished
Read is an .artificial. shredding of an unfinished subclone (unfinished means at least one gap
in the sequence . so multiple contigs). Features: nonrandom, untrusted, no mate pair.

AS_FBAC    Finished
Same as above but subclone is a single contig.

AS_STS     Sts
Never used or implemented. The idea is to allow constraints to be input to the assembly
based on some sequence or clone tag. For instance, the sequence .ACTGTGTGT.. should be
on the same chromosome as .GCGACGAGAT.. no orientation known. This is only a read or
fragment in the sense that there is an associated sequence. Features: nonrandom, not to be
used for consensus, mate pair relationship is present only by analogy.

AS_UNITIG  Unitig
This was used for an old flavor of the assembler . no longer used.

AS_CONTIG  Contig
This was used for an old flavor of the assembler . no longer used.

AS_BACTIG  BacTig
This was used for an old flavor of the assembler . no longer used.

AS_FULLBAC Full Bac C = Complete)
This was used for an old flavor of the assembler . no longer used.

AS_B_READ  BGLII read (should not be used)
This read is sequenced out from an internal spacer sequence added to a clone to
replace a section cut out by a restriction enzyme to make the size of the insert stable. Only
used once with crappy results . not sure code was fully implemented. Features: nonrandom,
trusted, multiple mate pair relationships . outtie, innie, normal, antinormal.
*/

/*

type/attribute table for used types only 

AS_READ	read			random	trusted			can have mate 
AS_EXTR	read			random		external read	can have mate 
AS_TRNR	read				trusted			can have mate 
AS_EBAC   	guide		random	trusted			can have mate 
AS_LBAC   				trusted			can have mate 
AS_UBAC   		shredded 
AS_FBAC   		shredded 
AS_STS     

*/

#define AS_FA_READ(type) 		((type == AS_READ) || (type == AS_EXTR) || (type == AS_EXTR) || (type == AS_B_READ))
#define AS_FA_RANDOM(type) 		((type == AS_READ) || (type == AS_EXTR) || (type == AS_EBAC))
#define AS_FA_SHREDDED(type) 		((type == AS_FBAC) || (type == AS_UBAC))
#define AS_FA_CAN_HAVE_MATE(type) 	((type == AS_READ) || (type == AS_EXTR) || (type == AS_TRNR) || (type == AS_LBAC) || (type == AS_EBAC)) 
#define AS_FA_GUIDE(type)         	(type == AS_EBAC)
#define AS_FA_TRUSTED_FOR_CNS(type) 	((type == AS_READ) || (type == AS_TRNR) || (type == AS_EBAC) || (type == AS_LBAC) || (type == AS_B_READ))
#define AS_FA_EXTERNAL_READ(type)	(type == AS_EXTR)

typedef enum {
  AS_UNIQUE_UNITIG   = (int)'U',  // U-Unique
  AS_ROCK_UNITIG     = (int)'R',  // Rock
  AS_STONE_UNITIG    = (int)'S',  // Stone
  AS_PEBBLE_UNITIG   = (int)'P',  // Pebble
  AS_SINGLE_UNITIG   = (int)'s',  // Singleton Unitig in Unplaced Contig
  AS_OTHER_UNITIG    = (int)'X'   // Unspecified surrogate unitig
} UnitigType;


/* Fragment messages, FRG, IFG, OFG */

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
  IntLocale_ID          ilocale;
  IntSequence_ID        iseq_id;
  IntBactig_ID          ibactig_id;
} FragMesg;

typedef FragMesg OFGMesg;
typedef FragMesg InternalFragMesg;

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


/* UOM message */

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
  AS_NO_OVERLAP	    = (int) 'N', // Not used by the unitigger.
  AS_OVERLAP        = (int) 'O', // A dovetail overlap between unitigs.
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
#ifdef AS_ENABLE_SOURCE
  char			*source;
#endif
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
#ifdef AS_ENABLE_SOURCE
  char			*source;
#endif
  CDS_COORD_t           best_overlap_length;
  CDS_COORD_t           min_overlap_length;
  CDS_COORD_t           max_overlap_length;
  float32               quality;
} FragOverlapMesg;

typedef struct {
  IntChunk_ID     iaccession;
  CDS_COORD_t     bp_length;
  float           coverage_stat;
  BranchType      a_branch_type;
  BranchType      b_branch_type;
  CDS_COORD_t     a_branch_point;
  CDS_COORD_t     b_branch_point;
  int32           num_frags;
  int32           a_degree;
  int32           b_degree;
  char           *source;
} ChunkMesg;


//  Though DirectionType is not explicitly referenced, AS_FORWARD and
//  AS_REVERSE are used in the code.
//
typedef enum {
  AS_FORWARD = (int)'F',
  AS_REVERSE = (int)'R'
} DirectionType;

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
#ifdef NEW_UNITIGGER_INTERFACE
  IntFragment_ID  ident2; // iid of the fragment that will align with current one
#endif
  
#ifdef AS_ENABLE_SOURCE
  int32        sourceInt;
#endif
#ifdef NEW_UNITIGGER_INTERFACE
  int32           ahang;
  int32           bhang;
#endif

  SeqInterval     position;
  IntFragment_ID  contained;
  int32           delta_length;
  int32          *delta;
#ifdef i386
  int32           ptrPad1; // 4 for pointer
#ifndef NEW_UNITIGGER_INTERFACE
  int32           ptrPad2; // 4 for 8-byte alignment
#endif
#endif
} IntMultiPos;

/* IMV message */

typedef struct IntMultiVar {
  SeqInterval     position;
  int32           num_reads;
  int32           num_conf_alleles;
  int32           anchor_size;
  int32           var_length;
  char           *nr_conf_alleles;
  char           *weights;
  char           *var_seq;
} IntMultiVar;

/* This is a variant of IntMultiPos to handle deltas in a longer (unitig) sequence */
typedef struct {
  UnitigType    type;
  IntUnitig_ID  ident;
  SeqInterval   position;
  int32         delta_length;
  int32         wordPad;
  int32        *delta;
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
#ifdef AS_ENABLE_SOURCE
  char		  *source;
#endif
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
  IntMultiPos               *pieces;
  IntUnitigPos              *unitigs;
  IntMultiVar               *v_list;
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

/* Grangers new list, roughly in order of precedence */
typedef enum {
  UNASSIGNED         = 'Z',
  GOOD_MATE          = 'G',
  BAD_SHORT_MATE     = 'C',
  BAD_LONG_MATE      = 'L',
  SAME_ORIENT_MATE   = 'S',
  OUTTIE_ORIENT_MATE = 'O',
  NO_MATE            = 'N',
  BOTH_CHAFF_MATE    = 'H',
  CHAFF_MATE         = 'A',
  BOTH_DEGEN_MATE    = 'D',
  DEGEN_MATE         = 'E',
  BOTH_SURR_MATE     = 'U',
  SURR_MATE          = 'R',
  DIFF_SCAFF_MATE    = 'F'
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
#ifdef AS_ENABLE_SOURCE
  char		*source;
#endif
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
#ifdef AS_ENABLE_SOURCE
  char		  *source;
#endif
  float32         coverage_stat;
  UnitigStatus    status;
  CDS_COORD_t     a_branch_point;
  CDS_COORD_t     b_branch_point;
  CDS_COORD_t     length;
  char            *consensus;
  char            *quality;
  int32		  forced;
  int32           num_frags;
  int32           num_vars;
  SnapMultiPos    *f_list;// changed in comparison to internal message
  IntMultiVar     *v_list;
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
  SnapMate_Pairs       *jump_list; // changed in comparison to internal message
} SnapUnitigLinkMesg;

/* CCO message */
typedef struct {
  Contig_ID                   eaccession;
  IntContig_ID                iaccession;
  ContigPlacementStatusType   placed;
  CDS_COORD_t                 length;
  char                       *consensus;
  char                       *quality;
  int32                       forced;
  int32                       num_pieces;
  int32                       num_unitigs;
  int32                       num_vars;
  SnapMultiPos               *pieces; // changed in comparison to internal message
  IntMultiVar                *vars;   
  UnitigPos                  *unitigs;// changed in comparison to internal message
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

/* CTP message */
typedef struct {
  Contig_ID		econtig1; // changed in comparison to internal message
  Contig_ID		econtig2; // changed in comparison to internal message
  float32		mean;
  float32		stddev;
  ChunkOrientationType	orient;
} SnapContigPairs;

/* SCF message */
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


/* EOF */
typedef struct EndOfFileMesgTag {
  int32   status;
  time_t  created;
  char    *comment;
} EndOfFileMesg;


//  External Routines


//  Function: DuplicateProtoMesg_AS 
//
//  Description: Copies the generic message structure including all of
//  its variable length data. The variable length data must be released
//  by free_mesg().
//
//  Return Value: A pointer to a copy of the generic message that has a
//  lifetime determined by the client.
//
//  Inputs: omesg - A pointer to a generic message.
//
//  This function is no longer supported.
//
//GenericMesg *DuplicateProtoMesg_AS(GenericMesg *omesg);


//  Function: FreeProtoMesg_AS 
//
//  Description: Frees the memory for a generic message structure.
//
//  Input/Outputs: omesg - A pointer to a generic message.
//
//  This function is no longer supported.
//
//void FreeProtoMesg_AS(GenericMesg *omesg);


//  Functions: ReadProtoMesg_AS 
//
//  Description: Reads the next message from the file "fin" and a
//  returns the memory location of a generic message.  This memory is
//  managed by the routine only. This memory location and all of its
//  variable data will go out of scope or be trashed by the next call
//  to the routine.
//  
//  Return Value: The return value is EOF if an end of file is
//  encountered occurs, otherwise the return value is zero indicating
//  success.
//  
//  Outputs: mesg - A handle to a generic message.
//  
//  Input/Outputs: fin - A file openned for text reading. 
//
int ReadProtoMesg_AS(FILE *fin, GenericMesg **pmesg);


//  Functions:  WriteProtoMesg_AS 
//
//  Description: Writes a generic message to the file "fout"
//
//  Return Value: The return value is negative if an error occured.
//
//  Inputs: mesg - A pointer to the generic message to be output. 
//
//  Input/Outputs: fout - A file openned for text writing.
//
int WriteProtoMesg_AS(FILE *fout, GenericMesg *mesg);


//  Function: Transfer_XXX_to_YYY
//
//  Description: Transfers the fields of an XXX message to that of a
//  YYY message.  Note carefully that second level memory is not duplicated,
//  so that such items (e.g. sequence) are *shared* between the two structures.
//  A kludge to allow items to be throughput without assignment of internal
//  IDs, should no longer be in use!
//
void Transfer_FRG_to_IFG_AS(FragMesg         *frg_mesg,
                            InternalFragMesg *ifg_mesg);

void Transfer_IFG_to_OFG_AS(InternalFragMesg *ifg_mesg,
                            OFGMesg          *ofg_mesg);

void Transfer_DST_to_IDT_AS(DistanceMesg     *dst_mesg,
                            InternalDistMesg *idt_mesg);

void Transfer_LKG_to_ILK_AS(LinkMesg         *lkg_mesg,
                            InternalLinkMesg *ilk_mesg);



void AppendAuditLine_AS(AuditMesg *adt_mesg,
                        AuditLine *auditLine,
                        time_t t, char *name,
                        char *version, char *comment);


//  Description: When reading in proto mode, this function will return
//  the line number the input is currently on.  This function operates
//  correctly only when a single input is being read.  The return value
//  in all other cases is the sum of the number of lines read in all
//  proto files thus far.
//
int GetProtoLineNum_AS(void);


//  Returns a number in the range [1, NUM_OF_REC_TYPES -1]
//  as a function of the first 3 characters of the passed string.
//
int GetMessageType(char *string);


//   Returns a string as a function of message type
//
const char  *GetMessageName(int type);


#endif  /* AS_MSG_PMESG_INCLUDE */

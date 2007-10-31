
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
/* $Id: AS_MSG_pmesg.h,v 1.57 2007-10-31 17:26:43 eliv Exp $   */

#ifndef AS_MSG_PMESG_INCLUDE_H
#define AS_MSG_PMESG_INCLUDE_H

#include <stdio.h>
#include <time.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"

// Defining the following enables internal source fields for testing
#define AS_ENABLE_SOURCE

//   #define NEW_UNITIGGER_INTERFACE

typedef CDS_IID_t   IntLibrary_ID;
typedef CDS_IID_t   IntDist_ID;
typedef CDS_IID_t   IntFragment_ID;
typedef CDS_IID_t   IntChunk_ID;
typedef CDS_IID_t   IntUnitig_ID;
typedef CDS_IID_t   IntContig_ID;
typedef CDS_IID_t   IntScaffold_ID;

typedef enum {
  AS_ADD      = (int)'A',
  AS_DELETE   = (int)'D',
  AS_IGNORE   = (int)'I',
  AS_UPDATE   = (int)'U'
} ActionType;

typedef struct {
  CDS_COORD_t bgn;
  CDS_COORD_t end;
} SeqInterval;

typedef enum {
  MESG_NUL = 0,
  MESG_ADT, MESG_VER, MESG_FRG, MESG_IFG, MESG_SPo, // 5
  MESG_LKG, MESG_SPg, MESG_DST, MESG_IDT, MESG_LIB, // 10
  MESG_SPc, MESG_SP1, MESG_OVL, MESG_SPf, MESG_UOM, // 15
  MESG_IUM, MESG_IUL, MESG_ICL, MESG_AFG, MESG_ISF, // 20
  MESG_IMD, MESG_IAF, MESG_UTG, MESG_ULK, MESG_ICM, // 25
  MESG_CCO, MESG_CLK, MESG_SCF, MESG_MDI, MESG_BAT, // 30  
  MESG_SPl, MESG_SPn, MESG_SPm, MESG_SP2, MESG_IBI, // 35
  MESG_SP3, MESG_SP4, MESG_SP5, MESG_SP6, MESG_SP7, // 40
  MESG_IDS, MESG_DSC, MESG_SLK, MESG_ISL, MESG_SPk, // 45
  MESG_SPd, MESG_SP8, MESG_SP9, MESG_SPa, MESG_EOF  // 50
} MessageType;

#define NUM_OF_REC_TYPES MESG_EOF

static const char  *MessageTypeName[NUM_OF_REC_TYPES + 1] = {
  "NUL",
  "ADT", "VER", "FRG", "IFG", "SPo", // 5
  "LKG", "SPg", "DST", "IDT", "LIB", // 10
  "SPc", "SP1", "OVL", "SPf", "UOM", // 15
  "IUM", "IUL", "ICL", "AFG", "ISF", // 20  
  "IMD", "IAF", "UTG", "ULK", "ICM", // 25 
  "CCO", "CLK", "SCF", "MDI", "BAT", // 30  
  "SPl", "SPn", "SPm", "SP2", "IBI", // 35
  "SP3", "SP4", "SP5", "SP6", "SP7", // 40
  "IDS", "DSC", "SLK", "ISL", "SPk", // 45
  "SPd", "SP8", "SP9", "SPa", "EOF"  // 50
};

/*Generic message object handle */

typedef struct { 
  void         *m;         // A pointer to the message object
  MessageType  t;          // The message type
  int32        s;          // The message size in bytes
} GenericMesg;

/* BAT record */

typedef struct InternalBatchMesgTag {
  char         *name;
  CDS_UID_t     eaccession;
  char         *comment;
}BatchMesg;

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

/* VER message */

typedef struct {
  uint32     version;
} VersionMesg;

/* LKG message */

typedef enum {
  AS_MATE       = (int)'M', // Mate
  AS_REREAD     = (int)'R', // Reread
  AS_MAY_JOIN   = (int)'Y', // maY
  AS_MUST_JOIN  = (int)'T'  // musT
  //AS_UNKNOWN	= (int)'U'  // the initial value, can't be defined here as it is in OrientType
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
  OrientType      link_orient;
  CDS_UID_t       frag1;
  CDS_UID_t       frag2;
  CDS_UID_t       distance;
} LinkMesg;

/* LIB message -- only for version 2 */

typedef struct {
  ActionType   action;
  CDS_UID_t    eaccession;
  float        mean;
  float        stddev;
#ifdef AS_ENABLE_SOURCE
  char        *source;
#endif
  OrientType   link_orient;
  uint32       num_features;
  char       **features;
  char       **values;
} LibraryMesg;

/* DST message  --  only for version 1 */

typedef struct {
  ActionType   action;
  CDS_UID_t    eaccession;
  float        mean;
  float        stddev;
} DistanceMesg;



//
//  XXX When more types are added, or when the library starts telling
//  what type the reads are, it will be VERY useful to remove ALL of
//  these types to catch the numerous places where the type of the
//  read is assumed to be AS_READ.
//

typedef enum {
  AS_READ    = (int)'R',  //  Celera Read
  AS_EXTR    = (int)'X',  //  External WGS read
  AS_TRNR    = (int)'T',  //  Transposon library read

  //  These aren't FragTypes exactly, but are used in CGW/eCR/CNS like
  //  they are.  Someone should figure out what they exactly are, and
  //  fix this.
  //
  AS_UNITIG  = (int)'U',  //  Assembled unitig
  AS_CONTIG  = (int)'C'   //  Assembled contig
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

type/attribute table for used types only 

AS_READ	read			random	trusted			can have mate 
AS_EXTR	read			random		external read	can have mate 
AS_TRNR	read				trusted			can have mate 

*/

#define AS_FA_READ(type) 		((type == AS_READ) || (type == AS_EXTR))
#define AS_FA_RANDOM(type) 		((type == AS_READ) || (type == AS_EXTR))
#define AS_FA_SHREDDED(type) 		(0)
#define AS_FA_CAN_HAVE_MATE(type) 	((type == AS_READ) || (type == AS_EXTR) || (type == AS_TRNR))
#define AS_FA_GUIDE(type)         	(0)
#define AS_FA_TRUSTED_FOR_CNS(type) 	((type == AS_READ) || (type == AS_TRNR))
#define AS_FA_EXTERNAL_READ(type)	((type == AS_EXTR))

typedef enum {
  AS_UNIQUE_UNITIG   = (int)'U',  // U-Unique
  AS_ROCK_UNITIG     = (int)'R',  // Rock
  AS_STONE_UNITIG    = (int)'S',  // Stone
  AS_PEBBLE_UNITIG   = (int)'P',  // Pebble
  AS_SINGLE_UNITIG   = (int)'s',  // Singleton Unitig in Unplaced Contig
  AS_OTHER_UNITIG    = (int)'X'   // Unspecified surrogate unitig
} UnitigType;

/* Fragment messages, FRG, IFG */

typedef struct {
  ActionType   		action;
  uint32                version;
  CDS_UID_t  		eaccession;
  CDS_UID_t             library_uid;     //  only version 2
  IntLibrary_ID         library_iid;     //  only version 2
  CDS_UID_t             plate_uid;       //  only version 2
  uint32                plate_location;  //  only version 2
  FragType     		type;            //  only version 1
  uint32                is_random;       //  only version 2
  char                  status_code;     //  only version 2
  SeqInterval  		clear_rng;
  SeqInterval  		clear_vec;       //  only version 2
  SeqInterval  		clear_qlt;       //  only version 2
  char        		*source;
  char        		*sequence;
  char        		*quality;
  char                  *hps;            //  only version 2
  IntFragment_ID   	iaccession;
} FragMesg;

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
  CDS_COORD_t      ahg,bhg;
  OrientType       orientation;
  OverlapType      overlap_type;
  float            quality;
  CDS_COORD_t      min_offset, max_offset;
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
  CDS_COORD_t           best_overlap_length;
  CDS_COORD_t           min_overlap_length;
  CDS_COORD_t           max_overlap_length;
  float                 quality;
} UnitigOverlapMesg;



typedef struct {
  IntChunk_ID     iaccession;
  CDS_COORD_t     bp_length;
  float           coverage_stat;
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
  AS_UNASSIGNED = (int)'X'
} UnitigStatus;

typedef enum {
  AS_FORCED_NONE    = (int)'X',
  AS_FORCED_UNIQUE  = (int)'U',
  AS_FORCED_REPEAT  = (int)'R'
} UnitigFUR;

/* IMP message */

typedef struct IntMultiPos {
  FragType        type;
  IntFragment_ID  ident;
#ifdef NEW_UNITIGGER_INTERFACE
  IntFragment_ID  ident2; // iid of the fragment that will align with current one
#endif
  
  //  This should probably be called frgSource
  IntFragment_ID  sourceInt;

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

VA_DEF(IntMultiPos);

/* IMV message */

typedef struct IntMultiVar {
  SeqInterval     position;
  int32           num_reads;
  int32           num_conf_alleles;
  int32           min_anchor_size;
  int32           var_length;
  int32           curr_var_id; // id of current VAR record
  int32           phased_var_id;  // id of the VAR record phased with the current one
  char           *nr_conf_alleles;
  char           *weights;
  char           *var_seq;
  char           *conf_read_iids; // iids of phased reads
#ifdef i386
  int32           ptrPad[4];
#endif
} IntMultiVar;

VA_DEF(IntMultiVar);

/* This is a variant of IntMultiPos to handle deltas in a longer (unitig) sequence */
typedef struct {
  UnitigType    type;
  IntUnitig_ID  ident;
  SeqInterval   position;
  int32         delta_length;
  int32         wordPad;
  int32        *delta;
#ifdef i386
  int32         ptrPad;
#endif
} IntUnitigPos;

VA_DEF(IntUnitigPos);

typedef struct {
  UnitigType   type;
  SeqInterval  position;
  int32        delta_length;
  int32        *delta;
  CDS_UID_t     eident;
} UnitigPos;

/* IEP messages */

typedef struct {
  FragType        type;
  IntFragment_ID  ident;
  SeqInterval     position;
} IntElementPos;

VA_DEF(IntElementPos);

/* IUM */

typedef struct {
  IntChunk_ID     iaccession;
#ifdef AS_ENABLE_SOURCE
  char		  *source;
#endif
  float           coverage_stat;
  UnitigStatus    status;
  UnitigFUR       unique_rept;
  CDS_COORD_t     length;
  char      *consensus;
  const char      *quality;
  int32		  forced;
  int32           num_frags;
  IntMultiPos    *f_list;
} IntUnitigMesg;

VA_DEF(IntUnitigMesg);  //  Used by unitigger.

typedef struct {
  CDS_UID_t       eaccession;
  UnitigStatus    status;
  int32           num_occurences;
  CDS_UID_t       *occurences;
  CDS_COORD_t     length;
  char            *consensus;
  char            *quality;
  int32           num_reads;
  int32           num_guides;
  int32           sum_delta_lengths;
} UnitigMesg;


typedef enum {
  AS_PLACED	= (int)'P',
  AS_UNPLACED   = (int)'U'
} ContigPlacementStatusType;

/* ICM */

typedef struct {
  IntContig_ID               iaccession;
  ContigPlacementStatusType  placed;
  CDS_COORD_t                length;
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
  float  		mean_distance;
  float  		std_deviation;
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
  float  		mean_distance;
  float  		std_deviation;
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
  float  		mean_distance;
  float  		std_deviation;
  int32			num_contributing;
  IntMate_Pairs		*jump_list;
} InternalScaffoldLinkMesg;


/* IAF message */

/* Grangers new list, roughly in order of precedence */
typedef enum {
  UNASSIGNED_MATE    = 'Z',
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
  CDS_UID_t       eaccession;
  IntFragment_ID  iaccession;              
  MateStatType    mate_status;
  int32           chimeric;
  int32           chaff;
  SeqInterval     clear_rng;
} AugFragMesg;

/* IMD message */

typedef struct {
  IntDist_ID  refines;
  float       mean;
  float       stddev;
  int32	      min;
  int32	      max;
  int32	      num_buckets;
  int32       *histogram;
} IntMateDistMesg;

/* ISF message */

typedef struct {
  IntContig_ID		contig1;
  IntContig_ID		contig2;
  float  		mean;
  float  		stddev;
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
  CDS_UID_t     eident;
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
  CDS_UID_t    eident;
  SeqInterval  position;
} SnapElementPos;

/* UTG Message */
typedef struct {
  CDS_UID_t       eaccession;  // changed in comparison to internal message
  IntChunk_ID     iaccession;
#ifdef AS_ENABLE_SOURCE
  char		  *source;
#endif
  float           coverage_stat;
  UnitigStatus    status;
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
  CDS_UID_t    in1, in2; 
  LinkType     type;
} SnapMate_Pairs;

/* ULK message */
typedef struct {
  CDS_UID_t		eunitig1;
  CDS_UID_t		eunitig2;
  ChunkOrientationType	orientation;
  UnitigOverlapType	overlap_type;
  int32			is_possible_chimera;
  int32			includes_guide;
  float  		mean_distance;
  float  		std_deviation;
  int32			num_contributing;
  PlacementStatusType	status;
  SnapMate_Pairs       *jump_list; // changed in comparison to internal message
} SnapUnitigLinkMesg;

/* CCO message */
typedef struct {
  CDS_UID_t                   eaccession;
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
  CDS_UID_t		econtig1; // changed in comparison to internal message
  CDS_UID_t		econtig2; // changed in comparison to internal message
  ChunkOrientationType	orientation;
  UnitigOverlapType	overlap_type;
  int32			is_possible_chimera;
  int32			includes_guide;
  float  		mean_distance;
  float  		std_deviation;
  int32			num_contributing;
  PlacementStatusType	status;
  SnapMate_Pairs	*jump_list; // changed in comparison to internal message
} SnapContigLinkMesg;

/* SLK message */
typedef struct {
  CDS_UID_t             escaffold1;
  CDS_UID_t             escaffold2;
  ChunkOrientationType	orientation;
  int32			includes_guide;
  float  		mean_distance;
  float  		std_deviation;
  int32			num_contributing;
  SnapMate_Pairs	*jump_list;
} SnapScaffoldLinkMesg;

/* CTP message */
typedef struct {
  CDS_UID_t		econtig1; // changed in comparison to internal message
  CDS_UID_t		econtig2; // changed in comparison to internal message
  float  		mean;
  float  		stddev;
  ChunkOrientationType	orient;
} SnapContigPairs;

/* SCF message */
typedef struct {
  CDS_UID_t             eaccession;
  IntScaffold_ID        iaccession;
  int32			num_contig_pairs;
  SnapContigPairs   	*contig_pairs; // changed in comparison to internal message
} SnapScaffoldMesg;

/* DSC message */
typedef struct {
  CDS_UID_t             eaccession;
  CDS_UID_t             econtig;
} SnapDegenerateScaffoldMesg;

/* MDI message */
typedef struct {
  CDS_UID_t		erefines; // changed in comparison to internal message
  IntDist_ID		irefines; // changed in comparison to internal message
  float			mean;
  float			stddev;
  int32			min;
  int32			max;
  int32			num_buckets;
  int32			*histogram;
} SnapMateDistMesg;



/* EOF */
typedef struct EndOfFileMesgTag {
  int32   status;
  char    *comment;
} EndOfFileMesg;




//  Semi-External Routines -- need to be exported for use in AS_MSG,
//  might be used and/or useful outside AS_MSG.

//  This function will return the line number the input is currently
//  on.  This function operates correctly only when a single input is
//  being read.  The return value in all other cases is the sum of the
//  number of lines read in all proto files thus far.
//
int GetProtoLineNum_AS(void);


//  Returns a number in the range [1, NUM_OF_REC_TYPES -1]
//  as a function of the first 3 characters of the passed string.
//
int GetMessageType(char *string);

//   Returns a string as a function of message type
//
const char  *GetMessageName(int type);

//  External Routines

void       AS_MSG_setFormatVersion(int format);

int        ReadProtoMesg_AS(FILE *fin, GenericMesg **pmesg);
int        WriteProtoMesg_AS(FILE *fout, GenericMesg *mesg);

void       AppendAuditLine_AS(AuditMesg *adt_mesg,
                              AuditLine *auditLine,
                              time_t t, char *name,
                              char *version, char *comment);


#endif  /* AS_MSG_PMESG_INCLUDE */

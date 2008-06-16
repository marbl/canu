
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
/* $Id: AS_MSG_pmesg.h,v 1.66 2008-06-16 22:53:26 brianwalenz Exp $   */

#ifndef AS_MSG_PMESG_INCLUDE_H
#define AS_MSG_PMESG_INCLUDE_H

#include <stdio.h>
#include <time.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"

// Defining the following enables internal source fields for testing
#define AS_ENABLE_SOURCE

//   #define NEW_UNITIGGER_INTERFACE

typedef AS_IID      IntLibrary_ID;
typedef AS_IID      IntDist_ID;
typedef AS_IID      IntFragment_ID;
typedef AS_IID      IntChunk_ID;
typedef AS_IID      IntUnitig_ID;
typedef AS_IID      IntContig_ID;
typedef AS_IID      IntScaffold_ID;

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
  MESG_BAT, MESG_VER, MESG_DST, MESG_LIB, MESG_FRG, MESG_LKG,
  MESG_OVL,
  MESG_UOM,
  MESG_IMD, MESG_IAF, MESG_IAM, MESG_IUM, MESG_IUL, MESG_ICM, MESG_ICL, MESG_ISF, MESG_ISL,
  MESG_MDI, MESG_AFG, MESG_AMP, MESG_UTG, MESG_ULK, MESG_CCO, MESG_CLK, MESG_SCF, MESG_SLK, 
  MESG_EOF
} MessageType;

#define NUM_OF_REC_TYPES MESG_EOF

static const char  *MessageTypeName[NUM_OF_REC_TYPES + 1] = {
  "NUL",
  "BAT", "VER", "DST", "LIB", "FRG", "LKG",
  "OVL",
  "UOM",
  "IMD", "IAF", "IAM", "IUM", "IUL", "ICM", "ICL", "ISF", "ISL",
  "MDI", "AFG", "AMP", "UTG", "ULK", "CCO", "CLK", "SCF", "SLK", 
  "EOF"
};

/*Generic message object handle */

typedef struct { 
  void         *m;         // A pointer to the message object
  MessageType  t;          // The message type
} GenericMesg;

/* BAT record */

typedef struct InternalBatchMesgTag {
  char         *name;
  AS_UID        eaccession;
  char         *comment;
}BatchMesg;

/* VER message */

typedef struct {
  uint32     version;
} VersionMesg;

/* LKG message */

typedef enum {
  AS_MATE       = (int)'M', // Mate
  AS_REREAD     = (int)'R', // Reread
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
  AS_UID          frag1;
  AS_UID          frag2;
  AS_UID          distance;
} LinkMesg;

/* LIB message -- only for version 2 */

typedef struct {
  ActionType   action;
  AS_UID       eaccession;
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
  AS_UID       eaccession;
  float        mean;
  float        stddev;
} DistanceMesg;



//
//  XXX When more types are added, or when the library starts telling
//  what type the reads are, it will be VERY useful to remove ALL of
//  these types to catch the numerous places where the type of the
//  read is assumed to be AS_READ.
//
//  After wasting a couple hours, BPW figured out why AS_UNITIG and
//  AS_CONTIG are FragTypes.  Consensus is using the IMP to store both
//  fragments and unitigs.
//

typedef enum {
  AS_READ    = (int)'R',  //  Celera Read
  AS_EXTR    = (int)'X',  //  External WGS read
  AS_TRNR    = (int)'T',  //  Transposon library read
  AS_UNITIG  = (int)'U',  //  Assembled unitig
  AS_CONTIG  = (int)'C'   //  Assembled contig
} FragType;

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

/* Fragment messages, FRG */

typedef struct {
  ActionType   		action;
  uint32                version;
  AS_UID     		eaccession;
  AS_UID                library_uid;     //  only version 2
  IntLibrary_ID         library_iid;     //  only version 2
  AS_UID                plate_uid;       //  only version 2
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
  AS_UID        eident;
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
  char            *consensus;
  char            *quality;
  int32		  forced;
  int32           num_frags;
  IntMultiPos    *f_list;
} IntUnitigMesg;

VA_DEF(IntUnitigMesg);  //  Used by unitigger.

typedef struct {
  AS_UID          eaccession;
  UnitigStatus    status;
  int32           num_occurences;
  AS_UID          *occurences;
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
  int32           chimeric_NOTUSED;
  int32           chaff;
  SeqInterval     clear_rng;
  MateStatType    mate_status;
} IntAugFragMesg;

/* AFG message */

typedef struct {
  AS_UID          eaccession;
  IntFragment_ID  iaccession;              
  MateStatType    mate_status;
  int32           chimeric_NOTUSED;
  int32           chaff;
  SeqInterval     clear_rng;
} AugFragMesg;

/* IAM message */

typedef struct {
  AS_IID          fragment1;
  AS_IID          fragment2;
  MateStatType    mate_status;
} IntAugMatePairMesg;

/* AMP message */

typedef struct {
  AS_UID          fragment1;
  AS_UID          fragment2;
  MateStatType    mate_status;
} AugMatePairMesg;

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


/* Genome Snapshot typedefs */
/****************************/

/* MPS message */
typedef struct {
  FragType      type;
  AS_UID        eident;
#ifdef AS_ENABLE_SOURCE
  char		*source;
#endif
  SeqInterval   position;
  int32         delta_length;
  int32         *delta;
} SnapMultiPos;

/* UTG Message */
typedef struct {
  AS_UID          eaccession;  // changed in comparison to internal message
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
  AS_UID       in1, in2; 
  LinkType     type;
} SnapMate_Pairs;

/* ULK message */
typedef struct {
  AS_UID   		eunitig1;
  AS_UID   		eunitig2;
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
  AS_UID                      eaccession;
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
  AS_UID   		econtig1; // changed in comparison to internal message
  AS_UID   		econtig2; // changed in comparison to internal message
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
  AS_UID                escaffold1;
  AS_UID                escaffold2;
  ChunkOrientationType	orientation;
  int32			includes_guide;
  float  		mean_distance;
  float  		std_deviation;
  int32			num_contributing;
  SnapMate_Pairs	*jump_list;
} SnapScaffoldLinkMesg;

/* CTP message */
typedef struct {
  AS_UID   		econtig1; // changed in comparison to internal message
  AS_UID   		econtig2; // changed in comparison to internal message
  float  		mean;
  float  		stddev;
  ChunkOrientationType	orient;
} SnapContigPairs;

/* SCF message */
typedef struct {
  AS_UID                eaccession;
  IntScaffold_ID        iaccession;
  int32			num_contig_pairs;
  SnapContigPairs   	*contig_pairs; // changed in comparison to internal message
} SnapScaffoldMesg;

/* MDI message */
typedef struct {
  AS_UID   		erefines; // changed in comparison to internal message
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

#endif  /* AS_MSG_PMESG_INCLUDE */

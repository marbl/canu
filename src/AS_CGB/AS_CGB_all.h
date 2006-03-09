
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
/*********************************************************************
 $Id: AS_CGB_all.h,v 1.6 2006-03-09 17:42:34 brianwalenz Exp $
 Module: Chunk Graph Builder
 Description: A catch-all include file for the Chunk Graph Builder
 Assumptions:
 Author: Clark M. Mobarry
 *********************************************************************/

#ifndef AS_CGB_ALL_INCLUDE
#define AS_CGB_ALL_INCLUDE

/* Condition Compilation */
#undef SAVE_FRGSRC_IN_VA

#undef STORE_BRANCH_POINTS_AT_FRAGMENT

#undef CONTAINMENT_STACKING
#undef USE_CROSS_POINTERS

#undef DONT_RUN_IN_SYMMETRIC_MODE    

#define USE_CORRECTED_ERROR_RATE

// Allow un-mated dovetail edges in CGB.

/* System include files */
//#define _LINT_
//#define _ANSI_C_SOURCE

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h> /* man 3 getopt */
#include <sys/types.h>

#include <fcntl.h>
#include <dirent.h>
#include <sys/stat.h>
#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Celera Assembler include files */
#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_UTL_Var.h"
#include "PrimitiveVA.h"
#include "AS_UTL_param_proc.h"

VA_DEF(BranchMesg)
VA_DEF(IntMultiPos)
VA_DEF(IntMultiVar);
VA_DEF(IntUnitigMesg)
VA_DEF(OFGMesg)
VA_DEF(OverlapMesg)

/**************************************************************/
/* Utility functions */

//  See also a copy in AS_FGB_FragmentHash.c

#define  SAFE_MALLOC(the_name, the_type, length) \
  assert(NULL == the_name); \
  the_name = (the_type *) malloc(length*sizeof(the_type)); \
  assert(NULL != the_name)
#define  SAFE_MALLOC_NOISY(the_name, the_type, length) \
  fprintf(stderr, "SAFE_MALLOC " F_SIZE_T " " #the_name " " #the_type  " " F_SIZE_T " " F_SIZE_T "\n", length * sizeof(the_type), length, sizeof(the_type) ); \
  SAFE_MALLOC(the_name, the_type, length)

#define  SAFE_CALLOC(the_name, the_type, length) \
  assert(NULL == the_name); \
  the_name = (the_type *) calloc(length,sizeof(the_type)); \
  assert(NULL != the_name)
#define  SAFE_CALLOC_NOISY(the_name, the_type, length) \
  fprintf(stderr, "SAFE_CALLOC " F_SIZE_T " " #the_name " " #the_type  " " F_SIZE_T " " F_SIZE_T "\n", length * sizeof(the_type), length, sizeof(the_type) ); \
  SAFE_CALLOC(the_name, the_type, length)

#define  SAFE_REALLOC(the_name, the_type, length) \
  assert(NULL != the_name); \
  the_name = (the_type *) realloc(the_name,length * sizeof(the_type)); \
  assert(NULL != the_name)
#define  SAFE_REALLOC_NOISY(the_name, the_type, length) \
  fprintf(stderr, "SAFE_REALLOC " F_SIZE_T " " #the_name " " #the_type  " " F_SIZE_T " " F_SIZE_T "\n", length * sizeof(the_type), length, sizeof(the_type) ); \
  SAFE_REALLOC(the_name, the_type, length)

#define SAFE_FREE(the_name) \
  assert(NULL != the_name); \
  free(the_name); \
  the_name = NULL
#define SAFE_FREE_NOISY(the_name) \
  fprintf(stderr, "SAFE_FREE " #the_name  "\n" ); \
  SAFE_FREE(the_name)


#define SAFE_FOPEN( file_handle, file_name, the_mode) \
assert(NULL == file_handle); \
assert(NULL != file_name); \
assert(NULL != the_mode); \
if(NULL == (file_handle = fopen(file_name, the_mode))) \
   { fprintf(stderr,"WARNING: Can not open file %s\n", file_name); exit(1); }

#define  SAFE_FCLOSE(file_handle) \
assert(NULL != file_handle); fclose(file_handle); file_handle = NULL;
  
#define warning_AS(expr) ((void) ((expr) ? 0 : (fprintf(stderr, "WARNING file=" __FILE__ " line=%d : %s\n",  __LINE__ , #expr ))))

#define ABS(a) ((a)>= 0 ?(a):-(a))

static void system_date(void) {
#if 1
  system("`date 1>&2`");
#endif  
}

static void system_top(void) {
  system_date();
}

/**************************************************************/
/* Global Variables */
#if 0
extern int TIMINGS;
extern int TIMINGS1;
#endif


/**************************************************************/
/* Global Preprocessor Variables */

#if 1
typedef cds_int64 BPTYPE;
#endif

# define BPFORMAT    F_S64
# define BPFORMATS10 "% 10" F_S64P
# define BPFORMAT6   "% 6" F_S64P
# define BPFORMAT15  "% 15" F_S64P

/* The Fancy formats */
#define FRAG_FORMAT "%15" F_IIDP

/* Exported constant definitions; Macro definitions; type definitions */
#define AS_CGB_NO_MATE 0
#define FRAGMENT_NOT_VISITED CDS_INT32_MAX
#define CHUNK_NOT_VISITED    CDS_INT32_MAX
typedef size_t IntRank;
typedef uint32 IntEdge_ID;
#define CMD_BUFFER_SIZE 1024


/* The transitive overlap inference slop is 
   
   abs(
   (length of the middle fragment) 
   + (length of the opposite overlap)
   - (length of the righthand overlap)
   - (length of the lefthand overlap)
   )
   <=  AS_CGB_TRANSITIVE_SLOP_ALPHA 
   + AS_CGB_TRANSITIVE_SLOP_EPSILON*(overlap in bp)
*/ 
#define AS_CGB_TRANSITIVE_SLOP_ALPHA   20     // bp
#define AS_CGB_TRANSITIVE_SLOP_EPSILON 0.07  // fraction 
// These values allow for 40 bp slop in an average 512 bp
// containment overlap.

#undef UNDIRECTED_DEBUG

typedef enum {

  AS_CGB_UNLABELED_FRAG=0,
  
  /* Initial non-contained fragment types */
  AS_CGB_SOLO_FRAG=1, 
  // A non-contained fragment with no raw dovetail overlaps at
  // all. This could be because of crappy sequencing, hard repetitive
  // sequence screening on both fragment ends, or a physical
  // sequencing gap.
  AS_CGB_HANGING_FRAG=2, 
  // A non-contained fragment with no raw dovetail overlaps on one
  // side. This could be because of crappy sequencing, hard repetitive
  // sequence screening on one fragment end, or a physical sequencing
  // gap.
  AS_CGB_THRU_FRAG=3, 
  // A non-contained fragment with raw dovetail overlaps on both
  // fragment ends. 

  /* Modified non-contained fragment types */
  AS_CGB_HANGING_CRAPPY_FRAG=4,    
  // A formerly hanging fragment that needs to be treated like a
  // contained unitig so that it does not inhibit chunking.
  AS_CGB_HANGING_CHUNK_FRAG=5, 
  // A formerly hanging fragment that terminates a non-singleton light
  // chunk. This could be due to crappy sequencing or hard repetitive
  // sequence screening.
  AS_CGB_INTERCHUNK_FRAG=6, 
  // A formerly thru fragment that terminates a light chunk and
  // has essentail non-containment overlaps to more than one light
  // chunk.
  AS_CGB_INTRACHUNK_FRAG=7, 
  // A formerly thru fragment that is interior to a light chunk. 

  /* Initial contained fragment type */
  AS_CGB_ORPHANEDCONT_FRAG=8,
  // A contained fragment that needs to be singleton unitig because it
  // was not placed into a unitig by a containment edge.

  /* Modified contained fragment types */
  AS_CGB_MULTICONT_FRAG=9, 
  // A contained fragment that needs to be singleton unitig because it
  // could be placed in multiple unitigs by containment edges.
  AS_CGB_SINGLECONT_FRAG=10, 
  // A contained fragment that has essential containment overlaps to
  // only one light chunk.

  /* spur & chimera & other bad fragments */
  AS_CGB_MARKED_BREAKER_FRAG=11,
  // soft - marked
  AS_CGB_REMOVED_BREAKER_FRAG=12,
  // hard - removed
  
  AS_CGB_UNPLACEDCONT_FRAG=13,
  // A contained fragment that has not yet been placed into a chunk.
  AS_CGB_BRANCHMULTICONT_FRAG=14,
  // A contained fragment that even after transitive overlap removal
  // still has a dovetail overlap to a non-contained fragment. Many of
  // these contained fragments are in repetitive regions near a
  // branch-point.
  AS_CGB_ESSENTIAL_CONT_FRAG=15,
  
  /* Other fragment types */
  AS_CGB_DELETED_FRAG=127

} Tlab;


//           dovetail 
// \superset nonchordal
// \superset nontransitive
// \superset nonblizzard
// \superset thickest (a directed edge concept)
// \superset interchunk (essential edges for chunking.)
// \superset buddychunk
// \superset intrachunk

// nonchordal-nonblizzard = nonchordal overlaps between two
// fragment-ends that are not buddies of each other but do have buddy
// overlaps.

// Most of the repeat separation and terrible unitigging errors occurs
// when the blizzard overlaps are discarded in reaping.


typedef enum {
  AS_CGB_UNUSED_EDGE=0,

  AS_CGB_DOVETAIL_EDGE=1,
  
  AS_CGB_THICKEST_EDGE=2,
  // A dovetail overlap edge that is the thickest from the proximal
  // fragment-end

  AS_CGB_BUDDY_EDGE=3,
  // A dovetail overlap that is (mutually) thickest from both
  // fragment-ends.

  AS_CGB_BETWEEN_CONTAINED_EDGE=4,
  // A dovetail overlap between globally contained fragments.

  AS_CGB_TOUCHES_CONTAINED_EDGE=6,
  // A dovetail overlap touching a globally contained fragment and a
  // non-contained fragment.

  /* Containment Overlap types (asymmetric edges) */
  AS_CGB_CONTAINED_EDGE=12,

  /* Dovetail Overlap types (symmetric edges) */
  AS_CGB_INTERCHUNK_EDGE=21,
  // A dovetail overlap exterior to a chunk.

  AS_CGB_BUDDYCHUNK_EDGE=22,
  // A dovetail overlap that qualifies as interchunk and is buddy.
  
  AS_CGB_INTRACHUNK_EDGE=23,
  // A dovetail overlap interior to a chunk.


  AS_CGB_TOUCHES_CRAPPY_DVT=32,
  // A dovetail overlap touching a crappy fragment and a non-crappy
  // fragment.
  AS_CGB_BETWEEN_CRAPPY_DVT=33,
  // A dovetail overlap between crappy fragments.
  AS_CGB_TOUCHES_CRAPPY_CON=38,
  AS_CGB_BETWEEN_CRAPPY_CON=39,


  AS_CGB_MARKED_BY_BRANCH_DVT=57,
  // An dovetail overlap removed by being on the repeat side of a
  // branch point.

  AS_CGB_MARKED_BY_BREAKER=58,
  // An edge to spur, chimera & other bad fragment
  
  AS_CGB_MARKED_BY_DELETED_DVT=61,
  AS_CGB_MARKED_BY_DELETED_CON=64,

  AS_CGB_REMOVED_BY_TRANSITIVITY_DVT=101,
  // An dovetail overlap that is inferrable by the dovetail chain overlap transivity
  // rules.
  AS_CGB_REMOVED_BY_TRANSITIVITY_CON=104,
  // An containment overlap that is inferrable by the overlap
  // transivity rules.

  AS_CGB_REMOVED_BY_THRESHOLD_DVT=111,
  // An dovetail overlap that is beyond the adjacency degree
  // threshold.
  AS_CGB_REMOVED_BY_THRESHOLD_CON=114,
  // An containment overlap that is beyond the adjacency degree
  // threshold.

  AS_CGB_REMOVED_BY_BREAKER=114,
  // An edge to spur, chimera & other bad fragment
  
  AS_CGB_REMOVED_BY_DUPLICATE_DVT=116,
  AS_CGB_REMOVED_BY_DUPLICATE_CON=119
  // An edge removed because it was a duplicate.

#if 0  
  // A dovetail overlap touching a deleted fragment.
  AS_CGB_REMOVED_BY_DELETED_DVT=121,
  AS_CGB_REMOVED_BY_DELETED_CON=124
#endif  

} Tnes;


/* Consult AS_UTL_Var.h. */

// The methods to access the vertex and edge data store.
#include "AS_CGB_methods.h"
#include "AS_CGB_store.h"
#include "AS_FGB_FragmentHash.h"

/* Prototypes for exported function definitions */
#include "AS_CGB_walk.h"
//#include "AS_CGB_branch_walk.h"
//#include "AS_CGB_triplets.h"
//#include "AS_CGB_contained.h"
#include "AS_CGB_edgemate.h"
#include "AS_CGB_fgb.h"
#include "AS_CGB_io.h"
#include "AS_CGB_cgb.h"
//#include "AS_CGB_branchpts.h"
#include "AS_CGB_fga.h"
#include "AS_CGB_cga.h"
#include "AS_CGB_histo.h"
#include "AS_CGB_traversal.h"
//#include "AS_CGB_edgetrimmer.h"
#include "AS_CGB_classify.h"
#include "AS_CGB_chimeras.h"
#include "AS_CGB_count_fragment_and_edge_labels.h"
#include "AS_FGB_io.h"

#endif /*AS_CGB_ALL_INCLUDE*/

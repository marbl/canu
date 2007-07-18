
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
 $Id: AS_CGB_all.h,v 1.17 2007-07-18 15:19:55 brianwalenz Exp $
 Module: Chunk Graph Builder
 Description: A catch-all include file for the Chunk Graph Builder
 Assumptions:
 Author: Clark M. Mobarry
 *********************************************************************/

#ifndef AS_CGB_ALL_INCLUDE
#define AS_CGB_ALL_INCLUDE

#undef CONTAINMENT_STACKING
#undef DONT_RUN_IN_SYMMETRIC_MODE    

// Allow un-mated dovetail edges in CGB.

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h> /* man 3 getopt */

#include <fcntl.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_param_proc.h"


#define ABS(a) ((a)>= 0 ?(a):-(a))

typedef int64 BPTYPE;

# define BPFORMAT    F_S64
# define BPFORMATS10 "% 10" F_S64P
# define BPFORMAT6   "% 6" F_S64P
# define BPFORMAT15  "% 15" F_S64P

/* The Fancy formats */
#define FRAG_FORMAT "%15" F_IIDP

/* Exported constant definitions; Macro definitions; type definitions */
#define AS_CGB_NO_MATE        0
#define FRAGMENT_NOT_VISITED INT32_MAX
#define CHUNK_NOT_VISITED    INT32_MAX

typedef size_t IntRank;
typedef uint32 IntEdge_ID;


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
// \superset thickest (a directed edge concept)
// \superset interchunk (essential edges for chunking.)
// \superset intrachunk


typedef enum {
  AS_CGB_UNUSED_EDGE=0,

  AS_CGB_DOVETAIL_EDGE=1,
  
  AS_CGB_THICKEST_EDGE=2,
  // A dovetail overlap edge that is the thickest from the proximal
  // fragment-end

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

} Tnes;



// The methods to access the vertex and edge data store.
#include "AS_CGB_methods.h"
#include "AS_CGB_store.h"
#include "AS_FGB_FragmentHash.h"

/* Prototypes for exported function definitions */
#include "AS_CGB_walk.h"
#include "AS_CGB_edgemate.h"
#include "AS_CGB_fgb.h"
#include "AS_CGB_cgb.h"
#include "AS_CGB_fga.h"
#include "AS_CGB_cga.h"
#include "AS_CGB_histo.h"
#include "AS_CGB_traversal.h"
#include "AS_CGB_classify.h"
#include "AS_CGB_chimeras.h"
#include "AS_CGB_count_fragment_and_edge_labels.h"

#endif /*AS_CGB_ALL_INCLUDE*/


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
/*************************************************
* Module:  GapFillREZ.c
* Description:
*   Assig  non-unique chunks to positions between unique-chunks
*   in (or between) scaffolds.
* 
*    Programmer:  A. Delcher
*       Written:  11 May 99
* 
* Assumptions:
*  Called from CGW main program.  Uses its global data structures.
* 
* Notes:
*
*************************************************/

static char fileID[] = "$Id: GapFillREZ.c,v 1.18 2007-02-04 09:30:45 brianwalenz Exp $";


#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_timer.h"
#include "dpc_CNS.h"

//
// AS_CGW
//
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ChiSquareTest_CGW.h"

//
// AS_REZ
//
#include "DataTypesREZ.h"
#include "UtilsREZ.h"
#include "CommonREZ.h"
#include "ConsistencyChecksREZ.h"
#include "GapWalkerREZ.h"
#include "SubgraphREZ.h"
#include "GWDriversREZ.h"
#include "RepeatRez.h"

// Temporary to show memory allocation

#if  0
#define  PALLOC(x)  fprintf (stderr, "  ALLOC %8ld bytes at line %5d\n", (long) (x), __LINE__);
#define  PRALLOC(x)  fprintf (stderr, "REALLOC %8ld bytes at line %5d\n", (long) (x), __LINE__);
#else
#define  PALLOC(x)  // Nothing
#define  PRALLOC(x)  // Nothing
#endif
                            
//#define  REF(i)   (Ref [i])   // Macro to reference chunk_ref data
#define  REF(i)   (Ref_Data [Ref_Index [i]])   // Macro to reference chunk_ref data

// #define  GAPS_TO_ADJUST  1;  j < fill_chunks [scaff_id] . num_gaps - 1;  j ++
#define  GAPS_TO_ADJUST  0;  j < fill_chunks [scaff_id] . num_gaps;  j ++

#define  DEBUG1                 0

#define  ALLOW_LOOSE_END_ROCKS     0
    // If  1  do usual rocks off ends of scaffolds; otherwise,
    // don't throw a rock off the end unless it's:  path-confirmed
    // or indicates a scaffold join.  Only allow ONE rock for
    // each indicated scaffold join.
#define  MAX_DIFF_SCAFFS           3
    // Most different scaffolds that a unitig can have links to without being
    // rejected.
#define  REQUIRE_SCAFFOLD_EDGE_FOR_JOIN  0
    // If  1  don't add a rock that appears to join two scaffolds
    // unless there is a scaffold edge that indicates the scaffolds
    // should be joined.
#define  SHOW_CALC_COORDS          0
    // If  1  then use calculated coordinates relative to celamy start
    // of scaffold for unique and placed chunks
#define  SHOW_FRAG_DETAILS         0
    // If  1  show individual fragment mate link info in
    // Print_Potential_Fill_Chunks
#define  SHOW_OLAP_DETAILS         0
    // If  1  then print info about overlaps that are computed
#define  SHOW_STONE_CONFIRM        0
    // If  1  print info about paths to confirm stones
#define  SIMPLE_ADJUST             1
    // If  1  adjust chunk positions to values from confirming overlap
    // path.  Otherwise, combine that value with the existing
    // position estimates.
#define  SKIP_CONTAINED_STONES     1
    // If  1  stones that are contained will not be used
#define  TEST_HOPELESS_SCAFFS      1
    // Ignore scaffolds not modified in previous round of rocks
#define  UNIQUES_CAN_BE_STONES     0
    // If  1  allow singleton discriminator-unique contigs to be
    // thrown as stones.
#define  USE_JOIN_SCAFFS           1
    // Whether or not to insert chunks that appear to separate 2 (or
    // more) scaffolds
#define  USE_MY_INSERT             0
    // If  1  use simple insert function for rocks instead of
    // old  Update_Scaffold_Graph -- this might not do anything,
    // look for the comments below (search for USE_MY_INSERT)
#define  VERBOSE                   0
    // If  1  print out lots of extra stuff


const double  CALC_SCAFF_SEPARATION = 100.0;
    // The number of bases to separate scaffolds in calc-cam files
const double  EPSILON = 0.5;
    // Added to variances to make sure they are in ascending order
    // for chunks in scaffolds.  Should be big enough to cover
    // any floating-point roundoff errors.
const double  EXPANSION_FACTOR = 1.3;
    // Amount to grow dynamically allocated arrays by.
const int  FILL_GAPS_THRESHOLD = 2;
    // If more than this many chunks inserted by level-1
    // then don't do stone throwing
const char  GAP_END_CHAR = '$';
const int  GOOD_LINKS_IF_BAD = 3;
    // Number of good links needed if there is a bad link present
const int  INITIAL_GAP_ENTRIES = 100;
    // Number of entries initially allocated for doublechecking
    // positions.
const int  INITIAL_REF_DATA_SIZE = 1000;
    // Number of Ref_Data entries to allocate initially
const int  INITIAL_SCAFF_JOIN_SIZE = 100;
    // Number of entries initially allocated for  Scaff_Join  array
const int  INITIAL_TREE_SIZE = 100;
    // Number of nodes in overlap-path dfs tree initially
const char  JOINER_ROCK_CHAR = '~';
const int32  MAX_MATE_DISTANCE = 50000;
    // Used to set end of gaps at ends of scaffolds
const int  MAX_SCAFFOLD_SIZE = 1000;
    // Most chunks in a scaffold that celamy can handle
const int  MAX_STRING_LEN = 1000;
    // Maximum length of annotation string in celamy
const int  MIN_GOOD_LINKS = 2;
    // The minimum number of "good" mate links to scaffold chunks
    // in order to add an unresolved chunk to the scaffold
const int  MIN_OLAP_LEN = 25;
    // The minimum number of bases of overlap needed to expect
    // overlap to be found in the hash lookup.
const int  MIN_ROCK_COVER_STAT = -2;      // Was  0
    // Only place rocks whose discriminator score is at least this large
const int  MIN_STONE_COVER_STAT = - 1000;
    // Only place stones whose discriminator score is at least this large
const int  MIN_TARGET_SLOP = 30.0;
    // The minimum number of bases either way that a target position
    // can be constrained to.
const double  MIN_VARIANCE = 3.0;
    // Smallest variance we're willing to set for a calculated position
const int  STONES_PER_CHECKPOINT = 25000;
    // The number of stones to throw before writing a checkpoint
const int  STACK_SIZE = 1000;
    // *** Make this dynamic eventually ***
    // Needs to be large enough to hold the most edges that any
    // chunk can can have to another
//const double  VARIANCE_FUDGE_FACTOR = CGW_FUDGE_FACTOR;
const double  VARIANCE_FUDGE_FACTOR = 3.0 * CGW_FUDGE_FACTOR;
    // Multiply this by total overlap length to estimate variance

#define  TorF(x)  ((x) ? 'T' : 'F')


char  * Colour_String [NUM_COLOURS]
    = {
       "C000000 T2 S  # Unused",
       "CFFFF00 T2 S  # Unique",
       "CFF8080 T2 S  # Inserted",
       "C7F7F7F T2 S  # Unconnected", 
       "C0000FF T2 S  # Unresolved",
       "C49F8F8 T2 S  # Consistent",
       "CFF0000 T2 S  # Placed",
       "CFF00FF T3 S  # Misplaced",
       "C7FFF00 T3 S  # Rejected",
       "CDFBF4F T2 S  # Scaffold",
       "CA0A0FF T2 S  # LoCoverStat",
       "CFF8000 T2 S  # Stone",
       "C80C080 T3 S  # rr_ru",
       "C00FFFF T3 S  # ur"
      };
FILE  * Cam_File;
#if  MAKE_CAM_FILE
#if  SHOW_CALC_COORDS
FILE  * Calc_Cam_File;
#endif
#endif


typedef  enum
  {
   ALL_FALSE,
   FALSE_IFF_SINGLETON
  }  Set_Split_t;

typedef  struct
  {
#if  CHECK_CELSIM_COORDS
   int  celsim_left, celsim_right;
#endif
   char  * annotation;
   int  calc_left, calc_right;
   int  scaff_id;
   signed int  colour : 6;
   unsigned int  flipped : 1;
   unsigned int  keep : 1;
   unsigned int  cgb_type : 4;
   signed int  use_ct : 16;
  }  Chunk_Info_t;


typedef  struct
  {
   int  scaff_id;
   int  rel_pos;
   LengthT  a_end, b_end;
   unsigned int  is_singleton : 1;       // true if in a scaffold by itself
   unsigned int  is_unthrowable : 1;     // true if chunk is in a scaffold
                                         //   that's been modified
  }  Chunk_Ref_t;


//  DFS_Info_t  contains nodes in the depth-first-search tree created
//  by  Find_Olap_Path

typedef  struct

  {
   int  id;
   double  distance;         // from the root in DNA bases
   int  max_hits;            // most target nodes reachable from this node
                             //   not counting this node itself
   int  max_first;           // subscript of first of those targets in max_hits
                             //   not counting this node itself
   int  max_first_dist;      // number of DNA bases from this node to the max_first
                             //   under it.  Done to accommodate different path
                             //   lengths of cross-edges in DFS tree.
   double  max_first_total;  // total length of fragments in overlaps that got
                             //   to  max_first .  Used to estimate variance.
   int  first;               // subscript of closest target node reachable from this node
                             //   not counting this node itself
   int  target_sub;          // subscript in target array of this node if it's
                             //   one of the targets
   unsigned  is_hit : 1;     // is this node a target of the search that
                             //   was hit within the specified distance range
   unsigned  succeeded : 1;  // did this node get to the destination (if there was one)
  }  DFS_Info_t;


//  Node_Index_t  is used to map nodes in the repeat-rez graph to
//  a depth-first-search tree used to find overlap paths

typedef  struct
  {
   unsigned  sub : 30;      // subscript in array of dfs-tree info
   unsigned  finished : 1;  // whether this node has finished exploring its subtree
   unsigned  visited : 1;   // whether this node has been visited or not
  }  Node_Index_t;


//  Placement_t  stores the tentative position of a chunk in a scaffold.
//  Used to check overlaps and adjust positions.

typedef  struct
  {
   LengthT  A_end, B_end;
   Gap_Chunk_t  * chunk;
   int  id;
   unsigned  keep : 1;
   unsigned  flipped : 1;
  }  Placement_t;


//  Interval_t  stores an interval of doubles (e.g., representing the
//  position of a chunk

typedef  struct
  {
   double  lo, hi;
  }  Interval_t;


typedef  struct
  {
   int  chunk_id;
       // Of a unique chunk.
   int  scaff_id;
       // Scaffold the unique chunk is in
   int  num_good_mates;
       // The number of mate links included in this edge, discounting
       // those flagged as suspicious in some way.
   int  num_tlaps;
   int  num_ulaps;
   int  num_rlaps;
   float  source_variance;
       // Variance of the end of the unique chunk that this edge
       // passes over, and hence the variance that induces the variance
       // on  left_end  and  right_end  below.
   unsigned int  flipped : 1;
       // Chunk relative to scaffold orientation
   unsigned int  left_link : 1;
       // If true this edge goes left from this chunk to an
       // earlier position in the scaffold, i.e., the scaffold chunk is on the
       // left of this chunk in scaffold order.
   unsigned int  is_bad : 1;
       // If true this edge should be ignored in computing chunk end
       // offset positions
   CIEdgeT  * edge;
   LengthT  left_end, right_end;
       // Of chunk in scaffold coordinates induced by this edge
   int  celsim_offset;
       // Difference between celsim coordinate and scaffold coordinate of
       // end of scaffold chunk passed over by this edge
   int  partition;
       // Indicates which group of edges this one belongs to when
       // creating multiple copies of a stone.
  }  Stack_Entry_t;


typedef  struct
  {
   int  scaff1, scaff2;
   int  chunk_id;
   int  m;
   double  b;
       // Coefficients of linear relation  y = mx + b
       // where  x  is scaff1 coordinate and  y  is scaff2 coordinate
   double  variance;
   int  insert_scaff;
   LengthT  left_end, right_end;
   int  gap;
   float  edge_quality;
   int  cover_stat;
   int  link_ct;
   Stack_Entry_t  stack_val;
   unsigned int  is_bad : 1;
   unsigned int  flipped : 1;
   unsigned int  violated : 1;
  }  Scaff_Join_t;


VA_DEF(Scaff_Join_t)


static Chunk_Info_t  * Chunk_Info = NULL;
    // Global array of info about each chunk
static int  Contained_Only_Switch = FALSE;
    // If  TRUE  then only insert rocks/stones that are
    // contained within a previously scaffolded contig.
char  * Filename_Prefix;
    // the prefix of the file name in the cgw command line
long int  Global_Fill_Info_Size;
    // To sum the number of bytes in the fill data structure
#if  TEST_HOPELESS_SCAFFS
static char  * Is_Hopeless_Scaff = NULL;
static int  Is_Hopeless_Size = 0;
static char  Hopeless_False_Mask;
static char  Hopeless_True_Mask;
#endif
static int  Num_Chunks;
    // Number of nodes in the rez graph at the start of repeat-rez
    // Number can grow as a result of insertions made.
static int  Num_Scaffolds;
    // Number of scaffolds in the rez graph at the start of repeat-rez
    // Number can grow/shrink as a result of insertions made.
#if  0
static Chunk_Ref_t  * Ref = NULL;
    // Global array of info about unique chunks
#else
static Chunk_Ref_t  * Ref_Data = NULL;
    // Global array of info about unique chunks
static int  * Ref_Index = NULL;
    // Global array of subscripts into above array to reduce wasted space.
static int  Ref_Data_Size = 0;
#endif
static VA_TYPE(Scaff_Join_t)  * Scaff_Join;
    // Global array of indications that scaffolds should merge
static int64  * Scaffold_End = NULL;
    // Global array of celsim start coordinates of each scaffold.
    // Used to make cam file showing calculated coordinates
static char  * Scaffold_Flipped = NULL;
    // Global array of booleans for each scaffold indicating its
    // orientation relative to celsim coordinates
static int64  * Scaffold_Start = NULL;
    // Global array of celsim start coordinates of each scaffold.
    // Used to make cam file showing calculated coordinates
static int  Single_Fragment_Only = FALSE;
    // When true, stones that contain more than one fragment will be
    // ignored
static int  Use_Partial_Stone_Paths = FALSE;
    // If true, then when throwing stones, allow paths that do
    // not span the entire gap


static void  Add_Gap_Ends
    (Scaffold_Fill_t * fill_chunks);
static void  Add_Join_Entry
    (int cid, int scaff1, int scaff2,
     LengthT a_end1, LengthT b_end1, LengthT a_end2, LengthT b_end2,
     LengthT left_end, LengthT right_end, int flipped,
     int gap, float edge_quality, int cover_stat, int link_ct,
     Stack_Entry_t stack_val);
static void  Adjust_By_Ref_Variance
    (Scaffold_Fill_t * fill_chunks);
static void  Adjust_By_Ref_Variance_One_Scaffold
    (Scaffold_Fill_t * fill_chunks, int scaff_id);
static void  Adjust_Positions
    (Gap_Fill_t * this_gap, int num_targets, Target_Info_t target [],
     int gap_sub [], int max_hits, int max_first, double start_coord,
     DirectionType direction, double * high_variance, FILE * fp);
static void  Adjust_Stone_Positions
     (int list [], int num_list, Gap_Chunk_t * node [], double ref_position,
      double factor, Path_Info_t path_info [], int target_sub,
      LengthT * target_position);
static void  Analyze_Placement
    (FILE * fp, Scaffold_Fill_t * fill);
void  Analyze_Rock_Fill
    (FILE * fp, Scaffold_Fill_t * fill_chunks);
void  Analyze_Stone_Fill
    (FILE * fp, Scaffold_Fill_t * fill_chunks);
static int  Ascending_Positions
    (const void * a, const void * b);
static int  Assign_To_Gap
    (int cid, LengthT left_end, LengthT right_end, int gap, int scaff_id,
     int flipped, Scaffold_Fill_t * fill_chunks, float edge_quality,
     int cover_stat, int link_ct, char id);
static int  Between
    (double a, double b, double lo, double hi);
static void  Build_Path_Subgraph
    (int start_sub, int target_sub, Gap_Chunk_t * node [], int num_nodes,
     int forward_edge [], int reverse_edge [], Stone_Edge_t pool [],
     int sorted [], int * num_sorted);
static int  By_Index
    (const void * a, const void * b);
static int  By_High_Placement
    (const void * a, const void * b);
static int  By_High_Position
    (const void * a, const void * b);
static int  By_Interval_Lo
    (const void * a, const void * b);
static int  By_Keep_And_Low_Position
    (const void * a, const void * b);
static int  By_Low_Position
    (const void * a, const void * b);
static int  By_Scaff_And_Flipped
    (const void * a, const void * b);
static int  By_Scaff_Flipped_And_Left_End
    (const void * a, const void * b);
static void  Check_Olaps
    (Gap_Fill_t * gap);
static void  Check_Other_Links
    (Scaffold_Fill_t * fill_chunks);
static void  Check_Other_Links_One_Scaffold
    (Scaffold_Fill_t * fill_chunks, int scaff_id,
     int * have_bad_ct, int * no_bad_ct, int * keep_ct, int * reject_ct);
static void  Check_Rock_Olaps
    (Gap_Fill_t * gap);
static void  Check_Rocks
    (FILE * fp, Scaffold_Fill_t * fill_chunks);
static int  Check_Scaffold_and_Orientation
    (int cid, Stack_Entry_t * stack, int stack_top, int * good_total,
     Scaffold_Fill_t * fill_chunks, int * bad_links, int min_good_links);
static void  Check_Scaffold_Join
    (int cid, Stack_Entry_t * stack, int stack_top, int scaff_id [],
     int scaff_links [], Scaffold_Fill_t * fill_chunks, int bad_links);
static int  Choose_Best_Stones
    (int start_id, int target_id, Gap_Chunk_t * node [], int num_nodes,
     int edge [], Stone_Edge_t pool [], double ref_position, double factor,
     LengthT * target_position);
static void  Choose_Safe_Chunks
    (Scaffold_Fill_t * fill_chunks, int min_good_links, int min_cover_stat);
static void  Choose_Stones
    (Scaffold_Fill_t * fill_chunks, int min_good_links, int min_cover_stat,
     int allow_bogus_edges);
static int  Chunk_Contained_In_Chunk
    (Gap_Chunk_t * A, Gap_Chunk_t * B);
static int  Chunk_Contained_In_Scaff
    (Gap_Chunk_t * A, int cid);
static void  Clear_Keep_Flags
    (Scaffold_Fill_t * fill_chunks, int except_num);
static void  Confirm_Contained
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int use_all);
static void  Confirm_Stones
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int use_all);
static int  Contained_In
    (Placement_t * p, Placement_t * q);
static int  Depth_First_Visit
    (ChunkInstanceT * from, int from_end, ChunkInstanceT * to,
     int num_targets, Target_Info_t target [], int bound, double so_far,
     int * first, int * max_hits, int * max_first, int * max_first_dist,
     double * max_first_total, int * tree_sub, int * tree_size,
     DFS_Info_t * * dfs_tree, Node_Index_t * node, int level,
     unsigned int edge_mask);
static int  Determine_Components
    (int list [], int num_list, Gap_Chunk_t * node [], int num_nodes,
     int start_sub, int target_sub, int edge [], Stone_Edge_t pool [],
     double ref_position, double factor, LengthT * target_position,
     int * complete_path, Gap_Fill_t * gap);
static void  DFS_Stone_Visit
    (int i, Gap_Chunk_t * node [], int edge [], Stone_Edge_t pool [],
     int * bad_edge);
static void  Disqualify_Scaff_Chunks
    (Scaffold_Fill_t * fill_chunks);
static void  Disqualify_Scaff_Chunks_One_Scaffold
    (Scaffold_Fill_t * fill_chunks, int scaff_id);
static void  Doublecheck_Positions
    (Scaffold_Fill_t * fill_chunks, int make_adjustments);
static void  Doublecheck_Positions_One_Scaffold
    (Scaffold_Fill_t * fill_chunks, int make_adjustments, int scaff_id);
static void  Eliminate_Encumbered_Uniques
    (Scaffold_Fill_t * fill);
static void  Eliminate_Encumbered_Uniques_One_Scaffold
    (Scaffold_Fill_t * fill, int scaff_id);
static int  Estimate_Chunk_Ends
    (Stack_Entry_t * stack, int stack_top,
     LengthT * left_end, LengthT * right_end, ChunkInstanceT * chunk,
     float * edge_quality, Scaffold_Fill_t * fill_chunks,
     int * gap, int * scaff_id, int * allowed_bad_links);
static void  Fasta_Print
    (FILE * fp, char * s, char * label);
static void  Fixup_Chunk_End_Variances
    (LengthT * left_end, LengthT * right_end, double diff);
void  Force_Increasing_Variances
    (void);
static void  Free_Global_Arrays
    (void);
static Overlap *  Get_Chunk_Overlap
    (Gap_Chunk_t * a, Gap_Chunk_t * b, char * * a_seq,
     char * * b_seq, FILE * fp);
static char *  Get_Contig_Sequence
    (int id);
static void  Identify_Best_Rocks
    (Scaffold_Fill_t * fill_chunks, int check_keep);
static void  Include_Good_Joins
    (Scaffold_Fill_t * fill_chunks);
static int  Insert_Chunks_In_Graph
    (ScaffoldGraphT * sgraph, Scaffold_Fill_t * fill, Kind_Of_Fill_t kind);
static int  Insert_Chunks_In_Graph_One_Scaffold
    (ScaffoldGraphT * sgraph, Scaffold_Fill_t * fill, int scaff_id,
     Kind_Of_Fill_t kind);
static int  Is_Good_Scaff_Edge
    (SEdgeT * edge);
static void UnJigglePositions(void);
static void  Jiggle_Positions
    (Scaffold_Fill_t * fill_chunks);
static void  Jiggle_Positions_One_Scaffold
    (Scaffold_Fill_t * fill_chunks, int scaff_id);
static int  Just_True
    (ContigT * chunk);
static void  Kill_Duplicate_Stones
    (Scaffold_Fill_t * fill);
static void  Kill_Duplicate_Stones_One_Scaffold
    (Scaffold_Fill_t * fill, int scaff_id);
static int  Maybe_Rock
    (ContigT * chunk);
static int  Maybe_Stone
    (ContigT * chunk);
static int  Might_Overlap
    (double a_frag_start, double a_frag_end,
     double b_frag_start, double b_frag_end,
     double slop, ChunkOrientationType * orient,
     int * min_ahang, int * max_ahang);
static void  New_Confirm_Stones
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int use_all);
static void  New_Confirm_Stones_One_Scaffold
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int use_all, int scaff_id);
static int  Num_Keep_Entries
    (Scaffold_Fill_t * fill, int scaff_id);
static void  Output_Cam_Files
    (Scaffold_Fill_t * fill);
static void  Partition_Edges
    (int cid, Stack_Entry_t * stack, int stack_top, int min_good_links);
#if  CHECK_CELSIM_COORDS
static void  Print_Coverage
    (FILE * fp);
#endif
void  Print_Fill_Info
    (FILE * fp, Scaffold_Fill_t * fill_chunks);
void  Print_Fill_Info_One_Scaffold
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int scaff_id,
     int * total_chunks, int * total_keep);
static void  Print_Frag_Info
    (FILE * fp, int cid);
static void  Print_Potential_Fill_Chunks
    (FILE * fp, int (* has_potential) (ContigT *), int allow_bogus_edges);
static void  Print_Scaffolds
    (FILE * fp);
static void  Print_Unique_Chunks
    (FILE * fp);
static int  Prior_Olaps_OK
    (int  from, int to, int to_low, Path_Info_t path_info [], int edge [],
     Stone_Edge_t pool []);
static void  Re_Check_Inserted_Rocks
    (Scaffold_Fill_t * fill, int min_good_links);
static void  Reject_Non_Reachable
    (int id, Gap_Chunk_t * node [], int n, int edge [], Stone_Edge_t pool []);
static void  Remove_From_Scaffold
    (Gap_Chunk_t * cp);
static int  Repeat_Colour
    (unsigned int cgb_type);
static void  Requalify_Scaff_Chunks
    (Scaffold_Fill_t * fill_chunks);
static void  Restore_Best_Rocks
    (Scaffold_Fill_t * fill_chunks);
static ChunkOrientationType  Reverse_Orientation
    (ChunkOrientationType orient);
static void  Reverse_Positions
    (Gap_Fill_t * this_gap);
static int  Scaff_Join_Cmp
    (const void *, const void *);
Scaffold_Fill_t *  Scan_Gaps
    (void);
static int  Select_Good_Edges
    (Stack_Entry_t * stack, int stack_top, ChunkInstanceT * chunk);
static void  Set_Is_Hopeless
    (Scaffold_Fill_t * fill);
static int  Set_Longest_Path
    (int list [], int num_list, Gap_Chunk_t * node [], int num_nodes,
     int target_sub, int edge [], Stone_Edge_t pool [], double ref_position,
     double factor, LengthT * target_position);
static void  Set_Position_From_Left_Olap
    (int left_cid, Gap_Chunk_t * this_chunk, Overlap * olap);
static void  Set_Position_From_Right_Olap
    (Gap_Chunk_t * this_chunk, int right_cid, Overlap * olap);
static void  Set_Split_Flags
    (Scaffold_Fill_t * fill, Set_Split_t set);
static void  Set_Split_Flags_One_Scaffold
    (Scaffold_Fill_t * fill, Set_Split_t set, int scaff_id);
static int  Should_Overlap
    (Placement_t * left, Placement_t * right, ChunkOrientationType * orient,
     double * how_much);
void  Show_Gap_Reads_One_Scaff
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int scaff_id);
static void  Show_Read_Info
    (FILE * fp, int cid);
static void  Sort_Insertions
    (Scaffold_Fill_t * fill_chunks,
     int (* cmp) (const void *, const void *));
static void  Sort_Insertions_One_Scaffold
    (Scaffold_Fill_t * fill_chunks,
     int (* cmp) (const void *, const void *), int scaff_id);
int  Throw_Stones
    (char *, int, int);
static void  Topologically_Sort
    (Gap_Chunk_t * node [], int num_nodes, int start_sub, int edge [],
     Stone_Edge_t pool [], int sorted [], int * num_sorted, int sort_all);
static void  Top_Sort_Visit
    (int i, Gap_Chunk_t * node [], int edge [], Stone_Edge_t pool [],
     int sorted [], int * num_sorted);
static void  Update_Colours
    (Scaffold_Fill_t * fill_chunks);
static void   Verify_Single_Placement
    (Scaffold_Fill_t * fill);
static int  Violates_Scaff_Edges
    (Scaff_Join_t  * p);
static void  Adjust_By_Ref_Variance_One_Scaffold(Scaffold_Fill_t * fill_chunks, int scaff_id);


int  Global_Debug_Flag = FALSE;


static void  Check_Loads
    (void)

//  Load multialignments of all chunks not yet in scaffolds.
//  Just for debugging--was seg faulting on this.

  {
   GraphNodeIterator  contig_iterator;
   ContigT  * chunk;

   InitGraphNodeIterator (& contig_iterator, ScaffoldGraph -> RezGraph,
                          GRAPH_NODE_DEFAULT);
   while  ((chunk = NextGraphNodeIterator (& contig_iterator)) != NULL)
     {
      if  (! Is_Unique (chunk))
          {
           MultiAlignT  * ma;
   
           fprintf (stderr, "### Check LoadMultiAlignTFromSequenceDB chunk = %d\n",
                    chunk -> id);
           ma = LoadMultiAlignTFromSequenceDB
                    (ScaffoldGraph -> sequenceDB, chunk -> id,
                     ScaffoldGraph -> RezGraph -> type == CI_GRAPH);
           fprintf (stderr, "### LoadMultiAlignTFromSequenceDB succeeded\n");
          }
     }

   return;
  }


static void  Add_Gap_Ends
    (Scaffold_Fill_t * fill_chunks)

//  Add scaffold chunks at ends of each gap to  fill_chunks
//  in any gap that already has a chunk

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;

         if  (this_gap -> num_chunks > 0)
             {
              int  cid, flipped, cover_stat, assign_succeeded;
              LengthT  left_end, right_end;

              if  (j > 0)
                  {
                   cid =  this_gap -> left_cid;
                   if  (REF (cid) . a_end . mean
                          <= REF (cid) . b_end . mean)
                       {
                        left_end = REF (cid) . a_end;
                        right_end = REF (cid) . b_end;
                        flipped = FALSE;

                        left_end . variance -= this_gap -> ref_variance;
                        if  (left_end . variance < MIN_VARIANCE)
                            left_end . variance = MIN_VARIANCE;

                        right_end . variance -= this_gap -> ref_variance;
                        if  (right_end . variance < MIN_VARIANCE)
                            right_end . variance = MIN_VARIANCE;
                       }
                     else
                       {
                        left_end = REF (cid) . b_end;
                        right_end = REF (cid) . a_end;
                        flipped = TRUE;

                        left_end . variance -= this_gap -> ref_variance;
                        if  (left_end . variance < MIN_VARIANCE)
                            left_end . variance = MIN_VARIANCE;

                        right_end . variance -= this_gap -> ref_variance;
                        if  (right_end . variance < MIN_VARIANCE)
                            right_end . variance = MIN_VARIANCE;
                       }
                   cover_stat = GetCoverageStat
                                    (GetGraphNode (ScaffoldGraph->RezGraph,
                                                   cid));
                   assign_succeeded
                       = Assign_To_Gap (cid, left_end, right_end,
                                        j, scaff_id, flipped, fill_chunks, 0.0,
                                        cover_stat, 0, GAP_END_CHAR);
                  }

              if  (j < fill_chunks [scaff_id] . num_gaps - 1)
                  {
                   cid =  this_gap -> right_cid;
                   if  (REF (cid) . a_end . mean
                          <= REF (cid) . b_end . mean)
                       {
                        left_end = REF (cid) . a_end;
                        right_end = REF (cid) . b_end;
                        flipped = FALSE;

                        left_end . variance -= this_gap -> ref_variance;
                        if  (left_end . variance < MIN_VARIANCE)
                            left_end . variance = MIN_VARIANCE;

                        right_end . variance -= this_gap -> ref_variance;
                        if  (right_end . variance < MIN_VARIANCE)
                            right_end . variance = MIN_VARIANCE;
                       }
                     else
                       {
                        left_end = REF (cid) . b_end;
                        right_end = REF (cid) . a_end;
                        flipped = TRUE;

                        left_end . variance -= this_gap -> ref_variance;
                        if  (left_end . variance < MIN_VARIANCE)
                            left_end . variance = MIN_VARIANCE;

                        right_end . variance -= this_gap -> ref_variance;
                        if  (right_end . variance < MIN_VARIANCE)
                            right_end . variance = MIN_VARIANCE;
                       }
                   cover_stat = GetCoverageStat
                                    (GetGraphNode (ScaffoldGraph->RezGraph,
                                                   cid));
                   assign_succeeded
                       = Assign_To_Gap (cid, left_end, right_end,
                                        j, scaff_id, flipped, fill_chunks, 0.0,
                                        cover_stat, 0, GAP_END_CHAR);
                  }
             }
        }
     }

   return;
  }



static void  Add_Join_Entry
    (int cid, int scaff1, int scaff2,
     LengthT a_end1, LengthT b_end1, LengthT a_end2, LengthT b_end2,
     LengthT left_end, LengthT right_end, int flipped, int gap,
     float edge_quality, int cover_stat, int link_ct,
     Stack_Entry_t stack_val)

//  Append to global  Scaff_Join  an entry that chunk  cid  had evidence
//  to join  scaff1  and  scaff2 .   a_end1 ,  b_end1 ,  a_end2  and
//  b_end2  are computed coordinates relative to these scaffolds.
//  left_end  and  right_end  are the extreme ends of the chunk in
//  scaff1 's coordinates.
//  flipped  indicates the orientation of the chunk's connection to
//  scaff1 .   gap  is the gap in that scaffold that the chunk would
//  be inserted in.   edge_quality  is the average quality of the connections
//  to that scaffold.   stack_val  is info about an edge that connects
//  this chunk to  scaff1 , where it will be inserted.
//  cover_stat  is the coverage statistic for this chunk.
//  link_ct  is the number of good mate links to *both* of the scaffolds
//  this chunk joins.

  {
   Scaff_Join_t  join;
   int  insert_scaff;
   double  m;

   assert (scaff1 != scaff2);

   insert_scaff = scaff1;
   if  (scaff2 < scaff1)
       {
        // Swap them so entries are consistent
        int  save;
        LengthT  save_len;

        save = scaff1;
        scaff1 = scaff2;
        scaff2 = save;

        save_len = a_end1;
        a_end1 = a_end2;
        a_end2 = save_len;

        save_len = b_end1;
        b_end1 = b_end2;
        b_end2 = save_len;
       }

   join . chunk_id = cid;
   join . scaff1 = scaff1;
   join . scaff2 = scaff2;
   join . variance = Max_double (a_end1 . variance, b_end1 . variance)
                       + Max_double (a_end2 . variance, b_end2 . variance);
   
   assert (b_end1 . mean != a_end1 . mean);
   assert (b_end2 . mean != a_end2 . mean);

   m = (b_end2 . mean - a_end2 . mean)
             / (b_end1 . mean - a_end1 . mean);
   if  (fabs (m - 1.0) <= 0.1)
       join . m = 1;
   else if  (fabs (m + 1.0) <= 0.1)
       join . m = -1;
     else
       {
        fprintf (stderr, "ERROR:  Bad m in Scaffold Join\n");
        fprintf (stderr, "  cid = %d  scaff1 = %d  scaff2 = %d  m = %.1f\n",
                 cid, scaff1, scaff2, m);
        return;
       }
   join . b = a_end2 . mean - m * a_end1 . mean;
   join . insert_scaff = insert_scaff;
   join . left_end = left_end;
   join . right_end = right_end;
   join . gap = gap;
   join . edge_quality = edge_quality;
   join . cover_stat = cover_stat;
   join . link_ct = link_ct;
   join . stack_val = stack_val;
   join . is_bad = FALSE;
   join . flipped = flipped;
   join . violated = FALSE;

   AppendScaff_Join_t (Scaff_Join, & join);

   return;
  }



static void  Adjust_Stone_Positions
     (int list [], int num_list, Gap_Chunk_t * node [], double ref_position,
      double factor, Path_Info_t path_info [], int target_sub,
      LengthT * target_position)

//  Adjust positions of entries  list [0 .. (num_list - 1)]  which are
//  subscripts of entries in  Node [] .
//  Use the corresponding information in  path_info [] .   ref_position
//  indicates the start position of the path and  factor  is either
//  +1.0  or  -1.0  indicating the direction of the path.
//  target_sub  is the subscript of the path destination if there is one;
//  otherwise, it's -1.  Set  (* target_position)  to the position
//  the stones indicate the end of gap should be.

  {
   LengthT  * far_end, * near_end, * max_position = NULL;
   double  prev_near_end = 0;
   int  found_target = FALSE;
   int  i;

   for  (i = 0;  i < num_list;  i ++)
     {
      double  old_len, new_len;
      int  sub = list [i];

      if  ((node [sub] -> flipped && factor > 0.0)
           || (! node [sub] -> flipped && factor < 0.0))
          {
           far_end = & (node [sub] -> start);
           near_end = & (node [sub] -> end);
          }
        else
          {
           far_end = & (node [sub] -> end);
           near_end = & (node [sub] -> start);
          }
      old_len = fabs (far_end -> mean - near_end -> mean);
      far_end -> mean
          = factor * path_info [sub] . hi_position + ref_position;
      far_end -> variance
          = Max_double
                (ComputeFudgeVariance (path_info [sub] . total_olap),
                 MIN_VARIANCE);
      if  (i == 0)
          near_end -> mean
              = far_end -> mean - factor * node [sub] -> len;
        else
          near_end -> mean = prev_near_end + factor * path_info [sub] . a_hang;
      prev_near_end = near_end -> mean;

      near_end -> variance
          = Max_double
              (ComputeFudgeVariance
                   (path_info [sub] . total_olap - node [sub] -> len),
               MIN_VARIANCE);

      if  (sub == target_sub)
          {
           (* target_position) = (* near_end);
           found_target = TRUE;
          }
      else if  (max_position == NULL
                  || far_end -> variance > max_position -> variance)
          {
           max_position = far_end;
          }
      new_len = fabs (far_end -> mean - near_end -> mean);
      if  (fabs (old_len - new_len) >= 20.0)
          {
           fprintf (stderr, "YIKES!!  old_len = %.0f  new_len = %.0f for chunk %d\n",
                    old_len, new_len, node [sub] -> chunk_id);
          }
     }

   if  (! found_target)
       (* target_position) = (* max_position);

   return;
  }



static void  Analyze_Placement
    (FILE * fp, Scaffold_Fill_t * fill)

//  List information about positions of entries in  fill  with
//  respect to their celsim coordinates.  Output goes to  fp .

  {
#if  CHECK_CELSIM_COORDS
   int  scaff_id;
   
   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
         ContigT  * contig;
         double  ref_celsim, ref_calc, slope, intercept;
         double  delta, denom;
         double  calc_gap_start, calc_gap_end;
         double  celsim_gap_start, celsim_gap_end;
         double  hi, lo, next;
         int  num_kept = 0;
         int  k, m;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                num_kept ++;
           }

         if  (num_kept == 0)
             continue;

         if  (j < fill [scaff_id] . num_gaps - 1)
             {
              contig
                  = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> right_cid);

              assert (contig != NULL);
              delta = contig -> offsetBEnd . mean
                        - contig -> offsetAEnd . mean;
              if  (delta >= 0.0)
                  {
                   ref_calc = contig -> offsetAEnd . mean;
                   ref_celsim = contig -> aEndCoord;
                  }
                else
                  {
                   ref_calc = contig -> offsetBEnd . mean;
                   ref_celsim = contig -> bEndCoord;
                  }
              calc_gap_end = ref_calc;
              celsim_gap_end = ref_celsim;
             }
         if  (j > 0)
             {
              contig
                  = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> left_cid);

              assert (contig != NULL);
              delta = contig -> offsetBEnd . mean
                        - contig -> offsetAEnd . mean;
              if  (delta >= 0.0)
                  {
                   ref_calc = contig -> offsetBEnd . mean;
                   ref_celsim = contig -> bEndCoord;
                  }
                else
                  {
                   ref_calc = contig -> offsetAEnd . mean;
                   ref_celsim = contig -> aEndCoord;
                  }
              calc_gap_start = ref_calc;
              celsim_gap_start = ref_celsim;
             }
         denom = contig -> bEndCoord - contig -> aEndCoord;
         assert (denom != 0.0);
         slope = rint (delta / denom);
         if  (fabs (slope) != 1.0)
             {
              fprintf (stderr, "ERROR:  Unexpected slope = %.2f  scaff %d  gap %d\n",
                       slope, scaff_id, j);
              fprintf (stderr, "        delta = %.1f  denom = %.1f\n",
                       delta, denom);
              if  (slope >= 0.0)
                  slope = 1.0;
                else
                  slope = -1.0;
             }
         intercept = ref_calc - slope * ref_celsim;

         fprintf (fp, "Scaff %4d  Gap %3d", scaff_id, j);
         if  (j > 0 && j < fill [scaff_id] . num_gaps - 1)
             fprintf (fp,
                      "  celsim: = [%8.0f,%8.0f] = %6.0f  calc: [%8.0f,%8.0f] = %6.0f",
                      celsim_gap_start, celsim_gap_end, 
                      slope * (celsim_gap_end - celsim_gap_start),
                      calc_gap_start, calc_gap_end,
                      calc_gap_end - calc_gap_start);
         fprintf (fp, "\n");

         hi = -1.0;
         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                {
                 contig
                     = GetGraphNode(ScaffoldGraph->RezGraph, this_chunk -> chunk_id);

                 assert (contig != NULL);

                 lo = Min_double (this_chunk -> start . mean,
                                  this_chunk -> end . mean);
                 if  (hi >= 0
                        && lo >= hi - 5.0)
                     fprintf (fp, "*** break ***\n");

                 next = Max_double (this_chunk -> start . mean,
                                    this_chunk -> end . mean);
                 if  (next > hi)
                     hi = next;
                 
                 fprintf (fp,
                 "  %6d  celsim: = [%8d,%8d]  calc: [%8.0f,%8.0f]"
                 "  delta: [%8.0f,%8.0f]\n",
                          this_chunk -> chunk_id,
                          contig -> aEndCoord, contig -> bEndCoord,
                          this_chunk -> start . mean,
                          this_chunk -> end . mean,
                          slope * contig -> aEndCoord + intercept
                            - this_chunk -> start . mean,
                          slope * contig -> bEndCoord + intercept
                            - this_chunk -> end . mean);
                }
           }
        }
     }
#endif

   return;
  }



void  Analyze_Rock_Fill
    (FILE * fp, Scaffold_Fill_t * fill_chunks)

//  List information about rocks in  fill_chunks .
//  Output goes to  fp .

  {
   int  scaff_id;
   int  inserted_chunks = 0;
   int  internal_rock_ct = 0;
   int  end_rock_ct = 0;
   int  links3_path_ct = 0;
   int  links3_nopath_ct = 0;
   int  links2_path_ct = 0;
   int  links2_nopath_ct = 0;
   int  scaffold_ct = 0;
   int  internal_gap_ct = 0;
   int  end_gap_ct = 0;
   int  internal_with_fill = 0;
   int  end_with_fill = 0;
   double  scaffold_bases = 0.0;
   double  internal_bases_filled = 0.0;
   double  end_bases_filled = 0.0;

   fprintf (fp, "%7s  %7s  %5s  %4s\n",
            "ID", "Length", "Links", "Conf");
   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      if  (fill_chunks [scaff_id] . num_gaps > 0)
          {
           scaffold_ct ++;
           if  (fill_chunks [scaff_id] . num_gaps < 2)
               fprintf (fp, "### Scaff id = %d has only 1 gap\n",
                        scaff_id);
             else
               {
                internal_gap_ct += fill_chunks [scaff_id] . num_gaps - 2;
                end_gap_ct += 2;
               }
          }

      for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
         ChunkInstanceT  * scaff_chunk;
         double  prev_left, target_right;
         int  has_internal_fill = FALSE;
         int  has_end_fill = FALSE;
         int  k;

         if  (j > 0)
             {
              scaff_chunk
                  = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> left_cid);
              scaffold_bases += fabs (scaff_chunk -> offsetAEnd . mean
                                        - scaff_chunk -> offsetBEnd . mean);
              prev_left = Max_double (scaff_chunk -> offsetAEnd . mean,
                                      scaff_chunk -> offsetBEnd . mean);
             }
           else
             {
              scaff_chunk
                  = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> right_cid);
              prev_left = - DBL_MAX;
             }
         if  (j < fill_chunks [scaff_id] . num_gaps - 1)
             {
              ChunkInstanceT  * right_chunk;

              right_chunk
                  = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> right_cid);
              target_right = Min_double (right_chunk -> offsetAEnd . mean,
                                         right_chunk -> offsetBEnd . mean);
             }
           else
             target_right = DBL_MAX;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                {
                 double  lo, hi;

                 lo = Min_double (this_chunk -> start .mean,
                                  this_chunk -> end .mean);
                 hi = Max_double (this_chunk -> start .mean,
                                  this_chunk -> end .mean);
                                  
                 inserted_chunks ++;
                 if  (j == 0 || j == fill_chunks [scaff_id] . num_gaps - 1)
                     end_rock_ct ++;
                   else
                     internal_rock_ct ++;
                 fprintf (fp, "%7d  %7.0f  %5d  %-4s\n",
                          this_chunk -> chunk_id,
                          hi - lo,
                          this_chunk -> link_ct,
                          this_chunk -> path_confirmed ? "Path" : "None");

                 if  (lo < prev_left)
                     lo = prev_left;
                 if  (hi > target_right)
                     hi = target_right;
                 if  (hi < lo)
                     hi = lo;
                 
                 if  (j == 0 || j == fill_chunks [scaff_id] . num_gaps - 1)
                     {
                      end_bases_filled += hi - lo;
                      has_end_fill = TRUE;
                     }
                   else
                     {
                      internal_bases_filled += hi - lo;
                      has_internal_fill = TRUE;
                     }
                 prev_left = hi;

                 if  (this_chunk -> link_ct > 2)
                     {
                      if  (this_chunk -> path_confirmed)
                          links3_path_ct ++;
                        else
                          links3_nopath_ct ++;
                     }
                   else
                     {
                      if  (this_chunk -> path_confirmed)
                          links2_path_ct ++;
                        else
                          links2_nopath_ct ++;
                     }
                }
           }
         if  (has_internal_fill)
             internal_with_fill ++;
         if  (has_end_fill)
             end_with_fill ++;
        }
     }

   fprintf (fp, "\n");
   fprintf (fp, "       Inserted rocks = %6d\n", inserted_chunks);
   fprintf (fp, "  (>=3)-links w/ path = %6d\n", links3_path_ct);
   fprintf (fp, "  (>=3)-links no path = %6d\n", links3_nopath_ct);
   fprintf (fp, "      2-links w/ path = %6d\n", links2_path_ct);
   fprintf (fp, "      2-links no path = %6d\n", links2_nopath_ct);
   fprintf (fp, "       scaffold count = %6d\n", scaffold_ct);
   fprintf (fp, "        internal gaps = %6d\n", internal_gap_ct);
   fprintf (fp, "             end gaps = %6d\n", end_gap_ct);
   fprintf (fp, "  internal rock count = %6d\n", internal_rock_ct);
   fprintf (fp, "       end rock count = %6d\n", end_rock_ct);
   fprintf (fp, " internal gaps w/fill = %6d\n", internal_with_fill);
   fprintf (fp, "      end gaps w/fill = %6d\n", end_with_fill);
   fprintf (fp, "internal bases filled = %8.0f\n", internal_bases_filled);
   fprintf (fp, "     end bases filled = %8.0f\n", end_bases_filled);
   fprintf (fp, "   bases in scaffolds = %8.0f\n", scaffold_bases);

#if  0
   fprintf (fp, "\n\n");
   fprintf (fp, "%7s  %7s  %7s  %7s  %8s\n",
            "Scaff", "Gap", "Len", "Inserts", "CoverLen");
   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 1;  j < fill_chunks [scaff_id] . num_gaps - 1;  j ++)
        {
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
         ChunkInstanceT  * left_scaff_chunk, * right_scaff_chunk;
         double  start_position, gap_len;
         static Interval_t  * place = NULL;
         static int  place_size = 0;
         int  k, num_inserts = 0;

         left_scaff_chunk
             = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> left_cid);
         right_scaff_chunk
             = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> right_cid);

         start_position = Max_double (left_scaff_chunk -> offsetAEnd . mean,
                                        left_scaff_chunk -> offsetBEnd . mean);
         gap_len = Min_double (right_scaff_chunk -> offsetAEnd . mean,
                                        right_scaff_chunk -> offsetBEnd . mean)
                     - start_position;

         if  (this_gap -> num_chunks > place_size)
             {
              place_size = 2 * this_gap -> num_chunks;
PRALLOC (place_size * sizeof (Interval_t));
              place = (Interval_t *) safe_realloc
                          (place, place_size * sizeof (Interval_t));
             }

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                {
                 place [num_inserts] . lo
                     = Min_double (this_chunk -> start . mean,
                                   this_chunk -> end . mean)
                         - start_position;
                 if  (place [num_inserts] . lo < 0)
                     place [num_inserts] . lo = 0;
                 place [num_inserts] . hi
                     = Max_double (this_chunk -> start . mean,
                                   this_chunk -> end . mean)
                         - start_position;
                 if  (place [num_inserts] . hi > gap_len)
                     place [num_inserts] . hi = gap_len;
                 if  (place [num_inserts] . hi < place [num_inserts] . lo)
                     place [num_inserts] . hi = place [num_inserts] . lo;
                 num_inserts ++;
                }
           }

         if  (num_inserts > 0)
             {
              double  left, right;
              double  combined_len = 0.0;
              int  i;

              qsort (place, num_inserts, sizeof (Interval_t), By_Interval_Lo);

              left = place [0] . lo;
              right = place [0] . hi;

              for  (i = 1;  i < num_inserts;  i ++)
                if  (place [i] . lo > right)
                    {
                     combined_len += right - left;
                     left = place [i] . lo;
                     right = place [i] . hi;
                    }
                else if  (place [i] . hi > right)
                    {
                     right = place [i] . hi;
                    }
              combined_len += right - left;

              fprintf (fp, "%7d  %7d  %7.0f  %7d  %8.0f\n",
                       scaff_id, j,
                       gap_len,
                       num_inserts, combined_len);
             }
        }
     }
#else
   fprintf (fp, "\n\n");
   fprintf (fp, "%7s  %7s  %7s  %7s  %8s  %7s\n",
            "Scaff", "Gap", "Thrown", "Placed", "CoverLen", "GapLen");
   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
         ChunkInstanceT  * left_scaff_chunk, * right_scaff_chunk;
         double  start_position, gap_len;
         double  combined_len;
         static Interval_t  * place = NULL;
         static int  place_size = 0;
         int  num_inserts = 0;
         int  num_thrown = 0;
         int  k;

         if  (j > 0)
             left_scaff_chunk
                 = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> left_cid);
           else
             left_scaff_chunk = NULL;
         if  (j < fill_chunks [scaff_id] . num_gaps - 1)
             right_scaff_chunk
                 = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> right_cid);
           else
             right_scaff_chunk = NULL;

         if  (left_scaff_chunk == NULL)
             start_position = 0.0;
           else
             start_position = Max_double (left_scaff_chunk -> offsetAEnd . mean,
                                          left_scaff_chunk -> offsetBEnd . mean);
         if  (right_scaff_chunk == NULL)
             gap_len = DBL_MAX;
           else
             gap_len = Min_double (right_scaff_chunk -> offsetAEnd . mean,
                                   right_scaff_chunk -> offsetBEnd . mean)
                          - start_position;

         if  (this_gap -> num_chunks > place_size)
             {
              place_size = 2 * this_gap -> num_chunks;
PRALLOC (place_size * sizeof (Interval_t));
              place = (Interval_t *) safe_realloc
                          (place, place_size * sizeof (Interval_t));
             }

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> copy_letter != GAP_END_CHAR)
                num_thrown ++;
            if  (this_chunk -> keep)
                {
                 place [num_inserts] . lo
                     = Min_double (this_chunk -> start . mean,
                                   this_chunk -> end . mean)
                         - start_position;
                 if  (place [num_inserts] . lo < 0)
                     place [num_inserts] . lo = 0;
                 place [num_inserts] . hi
                     = Max_double (this_chunk -> start . mean,
                                   this_chunk -> end . mean)
                         - start_position;
                 if  (place [num_inserts] . hi > gap_len)
                     place [num_inserts] . hi = gap_len;
                 if  (place [num_inserts] . hi < place [num_inserts] . lo)
                     place [num_inserts] . hi = place [num_inserts] . lo;
                 num_inserts ++;
                }
           }

         combined_len = 0.0;
         if  (num_inserts > 0)
             {
              double  left, right;
              int  i;

              qsort (place, num_inserts, sizeof (Interval_t), By_Interval_Lo);

              left = place [0] . lo;
              right = place [0] . hi;

              for  (i = 1;  i < num_inserts;  i ++)
                if  (place [i] . lo > right)
                    {
                     combined_len += right - left;
                     left = place [i] . lo;
                     right = place [i] . hi;
                    }
                else if  (place [i] . hi > right)
                    {
                     right = place [i] . hi;
                    }
              combined_len += right - left;

             }

         fprintf (fp, "%7d  %7d  %7d  %7d  %8.0f",
                  scaff_id, j, num_thrown, num_inserts, combined_len);
         if  (left_scaff_chunk != NULL && right_scaff_chunk != NULL)
             fprintf (fp, "  %7.0f", gap_len);
         fprintf (fp, "\n");
        }
     }
#endif

   return;
  }



void  Analyze_Stone_Fill
    (FILE * fp, Scaffold_Fill_t * fill_chunks)

//  List information about stones in  fill_chunks .
//  Output goes to  fp .

  {
   int  scaff_id;
   int  inserted_chunks = 0;
   int  internal_stone_ct = 0;
   int  end_stone_ct = 0;
   int  scaffold_ct = 0;
   int  internal_gap_ct = 0;
   int  end_gap_ct = 0;
   int  internal_with_fill = 0;
   int  end_with_fill = 0;
   double  scaffold_bases = 0.0;
   double  internal_bases_filled = 0.0;
   double  end_bases_filled = 0.0;

   fprintf (fp, "%7s  %7s  %5s\n",
            "ID", "Length", "Links");
   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      if  (fill_chunks [scaff_id] . num_gaps > 0)
          {
           scaffold_ct ++;
           if  (fill_chunks [scaff_id] . num_gaps < 2)
               fprintf (fp, "### Scaff id = %d has only 1 gap\n",
                        scaff_id);
             else
               {
                internal_gap_ct += fill_chunks [scaff_id] . num_gaps - 2;
                end_gap_ct += 2;
               }
          }

      for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
         ChunkInstanceT  * scaff_chunk;
         double  prev_left, target_right;
         int  has_internal_fill = FALSE;
         int  has_end_fill = FALSE;
         int  k;

         if  (j > 0)
             {
              scaff_chunk
                  = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> left_cid);
              scaffold_bases += fabs (scaff_chunk -> offsetAEnd . mean
                                        - scaff_chunk -> offsetBEnd . mean);
              prev_left = Max_double (scaff_chunk -> offsetAEnd . mean,
                                      scaff_chunk -> offsetBEnd . mean);
             }
           else
             {
              scaff_chunk
                  = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> right_cid);
              prev_left = - DBL_MAX;
             }
         if  (j < fill_chunks [scaff_id] . num_gaps - 1)
             {
              ChunkInstanceT  * right_chunk;

              right_chunk
                  = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> right_cid);
              target_right = Min_double (right_chunk -> offsetAEnd . mean,
                                         right_chunk -> offsetBEnd . mean);
             }
           else
             target_right = DBL_MAX;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                {
                 double  lo, hi;

                 lo = Min_double (this_chunk -> start .mean,
                                  this_chunk -> end .mean);
                 hi = Max_double (this_chunk -> start .mean,
                                  this_chunk -> end .mean);
                                  
                 inserted_chunks ++;
                 if  (j == 0 || j == fill_chunks [scaff_id] . num_gaps - 1)
                     end_stone_ct ++;
                   else
                     internal_stone_ct ++;
                 fprintf (fp, "%7d  %7.0f  %5d\n",
                          this_chunk -> chunk_id,
                          hi - lo,
                          this_chunk -> link_ct);
                 if  (lo < prev_left)
                     lo = prev_left;
                 if  (hi > target_right)
                     hi = target_right;
                 if  (hi < lo)
                     hi = lo;
                 
                 if  (j == 0 || j == fill_chunks [scaff_id] . num_gaps - 1)
                     {
                      end_bases_filled += hi - lo;
                      has_end_fill = TRUE;
                     }
                   else
                     {
                      internal_bases_filled += hi - lo;
                      has_internal_fill = TRUE;
                     }
                 prev_left = hi;
                }
           }
         if  (has_internal_fill)
             internal_with_fill ++;
         if  (has_end_fill)
             end_with_fill ++;
        }
     }

   fprintf (fp, "\n");
   fprintf (fp, "      Inserted stones = %6d\n", inserted_chunks);
   fprintf (fp, "       scaffold count = %6d\n", scaffold_ct);
   fprintf (fp, "        internal gaps = %6d\n", internal_gap_ct);
   fprintf (fp, "             end gaps = %6d\n", end_gap_ct);
   fprintf (fp, " internal stone count = %6d\n", internal_stone_ct);
   fprintf (fp, "      end stone count = %6d\n", end_stone_ct);
   fprintf (fp, " internal gaps w/fill = %6d\n", internal_with_fill);
   fprintf (fp, "      end gaps w/fill = %6d\n", end_with_fill);
   fprintf (fp, "internal bases filled = %8.0f\n", internal_bases_filled);
   fprintf (fp, "     end bases filled = %8.0f\n", end_bases_filled);
   fprintf (fp, "   bases in scaffolds = %8.0f\n", scaffold_bases);


   fprintf (fp, "\n\n");
   fprintf (fp, "%7s  %7s  %7s  %7s  %8s  %7s\n",
            "Scaff", "Gap", "Thrown", "Placed", "CoverLen", "GapLen");
   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
         ChunkInstanceT  * left_scaff_chunk, * right_scaff_chunk;
         double  start_position, gap_len;
         double  combined_len;
         static Interval_t  * place = NULL;
         static int  place_size = 0;
         int  num_inserts = 0;
         int  num_thrown = 0;
         int  k;

         if  (j > 0)
             left_scaff_chunk
                 = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> left_cid);
           else
             left_scaff_chunk = NULL;
         if  (j < fill_chunks [scaff_id] . num_gaps - 1)
             right_scaff_chunk
                 = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> right_cid);
           else
             right_scaff_chunk = NULL;

         if  (left_scaff_chunk == NULL)
             start_position = 0.0;
           else
             start_position = Max_double (left_scaff_chunk -> offsetAEnd . mean,
                                          left_scaff_chunk -> offsetBEnd . mean);
         if  (right_scaff_chunk == NULL)
             gap_len = DBL_MAX;
           else
             gap_len = Min_double (right_scaff_chunk -> offsetAEnd . mean,
                                   right_scaff_chunk -> offsetBEnd . mean)
                          - start_position;

         if  (this_gap -> num_chunks > place_size)
             {
              place_size = 2 * this_gap -> num_chunks;
PRALLOC (place_size * sizeof (Interval_t));
              place = (Interval_t *) safe_realloc
                          (place, place_size * sizeof (Interval_t));
             }

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> copy_letter != GAP_END_CHAR)
                num_thrown ++;
            if  (this_chunk -> keep)
                {
                 place [num_inserts] . lo
                     = Min_double (this_chunk -> start . mean,
                                   this_chunk -> end . mean)
                         - start_position;
                 if  (place [num_inserts] . lo < 0)
                     place [num_inserts] . lo = 0;
                 place [num_inserts] . hi
                     = Max_double (this_chunk -> start . mean,
                                   this_chunk -> end . mean)
                         - start_position;
                 if  (place [num_inserts] . hi > gap_len)
                     place [num_inserts] . hi = gap_len;
                 if  (place [num_inserts] . hi < place [num_inserts] . lo)
                     place [num_inserts] . hi = place [num_inserts] . lo;
                 num_inserts ++;
                }
           }

         combined_len = 0.0;
         if  (num_inserts > 0)
             {
              double  left, right;
              int  i;

              qsort (place, num_inserts, sizeof (Interval_t), By_Interval_Lo);

              left = place [0] . lo;
              right = place [0] . hi;

              for  (i = 1;  i < num_inserts;  i ++)
                if  (place [i] . lo > right)
                    {
                     combined_len += right - left;
                     left = place [i] . lo;
                     right = place [i] . hi;
                    }
                else if  (place [i] . hi > right)
                    {
                     right = place [i] . hi;
                    }
              combined_len += right - left;

             }

         fprintf (fp, "%7d  %7d  %7d  %7d  %8.0f",
                  scaff_id, j, num_thrown, num_inserts, combined_len);
         if  (left_scaff_chunk != NULL && right_scaff_chunk != NULL)
             fprintf (fp, "  %7.0f", gap_len);
         fprintf (fp, "\n");
        }
     }

   return;
  }



static int  Ascending_Positions
    (const void * a, const void * b)

//  First compare  keep  flags of  a  and  b .  If not equal
//  return the value that indicates the  keep = TRUE  value
//  is first.
//  Otherwise, compare scaffold positions of  Gap_Chunk_t 's  a
//  and  b  and return  -1  if  a  is before  b ,  0  if  a  and
//  b  are in the same place , and  1  if  b  is before  a .
//  Used for  qsort .

  {
   double  x_left, x_right, y_left, y_right;
   Gap_Chunk_t  * x, * y;

   x = (Gap_Chunk_t *) a;
   y = (Gap_Chunk_t *) b;

   if  (x -> keep && ! y -> keep)
       return  -1;
   else if  (y -> keep && ! x -> keep)
       return  1;

   if  (x -> start . mean < x -> end . mean)
       {
        x_left = x -> start . mean;
        x_right = x -> end . mean;
       }
     else
       {
        x_left = x -> end . mean;
        x_right = x -> start . mean;
       }

   if  (y -> start . mean < y -> end . mean)
       {
        y_left = y -> start . mean;
        y_right = y -> end . mean;
       }
     else
       {
        y_left = y -> end . mean;
        y_right = y -> start . mean;
       }

   if  (x_left < y_left)
       return  -1;
   else if  (x_left > y_left)
       return  1;
   else if  (x_right < y_right)
       return  -1;
   else if  (x_right > y_right)
       return  1;

   return  0;
  }



static int  Assign_To_Gap
    (int cid, LengthT left_end, LengthT right_end, int gap, int scaff_id,
     int flipped, Scaffold_Fill_t * fill_chunks, float edge_quality,
     int cover_stat, int link_ct, char id)

//  Assign this chunk with id =  cid  to gap number  gap
//  in scaffold  scaff_id .  The end coordinates of this chunk
//  are  left_end  and  right_end .  
//  Add it to  fill_chunks  structure and set the  copy_letter
//  field there to  id .
//  flipped  indicates if this chunk is reversed
//  in the scaffold.
//  cover_stat  is the discriminator-unique coverage stat for this chunk
//  link_ct  is the number of edge mate links determining the postion
//  of this chunk in the scaffold.  Return  TRUE  if the assignment
//  succeeded;  FALSE, otherwise.

  {
   Gap_Chunk_t  * this_chunk;
   Gap_Fill_t  * g;
   double  left_extreme, right_extreme;

#if  0
if  (scaff_id > 0)
    return  FALSE;
#endif
#if  0
fprintf (stderr, ">>> Adding cid = %d to gap %d in scaff %d\n",
         cid, gap, scaff_id);
#endif

#if  TEST_HOPELESS_SCAFFS
   if  (Is_Hopeless_Scaff [scaff_id] & Hopeless_True_Mask)
       return  FALSE;
#endif

   Chunk_Info [cid] . flipped = flipped;

   assert (gap < fill_chunks [scaff_id] . num_gaps);
   g = fill_chunks [scaff_id] . gap + gap;

   assert (g -> start . mean <= g -> end . mean);

   if  (Contained_Only_Switch)
       {
        left_extreme = left_end . mean + 5.0 * sqrt (left_end . variance);
        right_extreme = right_end . mean - 5.0 * sqrt (right_end . variance);

        if  (g -> len > 0 && right_extreme > g -> start . mean
                  && left_extreme < g -> end . mean)
           {
#if  VERBOSE
            fprintf (stderr, ">>> Chunk %d not a contained at gap %d of scaff %d\n",
                     cid, gap, scaff_id);
#endif
            return  FALSE;
           }
       }
     else
       {
        left_extreme = left_end . mean - 5.0 * sqrt (left_end . variance);
        right_extreme = right_end . mean + 5.0 * sqrt (right_end . variance);
        if  (! Interval_Intersection (left_extreme, right_extreme,
                                      g -> start . mean, g -> end . mean))
           {
#if  VERBOSE
            fprintf (stderr, ">>> Chunk %d was contained at gap %d of scaff %d\n",
                     cid, gap, scaff_id);
#endif
            return  FALSE;
           }
       }

   g -> num_chunks ++;
PRALLOC (g -> num_chunks * sizeof (Gap_Chunk_t));
   g -> chunk = (Gap_Chunk_t *) safe_realloc
                    (g -> chunk,
                     g -> num_chunks
                       * sizeof (Gap_Chunk_t));

   this_chunk = g -> chunk + (g -> num_chunks - 1);
   this_chunk -> chunk_id = cid;
   this_chunk -> scaff_id = scaff_id;
   this_chunk -> gap = gap;
   this_chunk -> cover_stat = cover_stat;
   this_chunk -> link_ct = link_ct;
   this_chunk -> avg_edge_quality = edge_quality;
   this_chunk -> index = -1;
   this_chunk -> copy_letter = id;
   this_chunk -> keep = FALSE;
   this_chunk -> flipped = flipped;
   this_chunk -> split = FALSE;
   this_chunk -> best = FALSE;
   this_chunk -> candidate = FALSE;
   this_chunk -> path_confirmed = FALSE;
   this_chunk -> len = (int32) (right_end . mean - left_end . mean);
   if  (flipped)
       {
        this_chunk -> start . mean = right_end . mean;
        this_chunk -> start . variance = right_end . variance;
        this_chunk -> end . mean = left_end . mean;
        this_chunk -> end . variance = left_end . variance;
       }
     else
       {
        this_chunk -> start . mean = left_end . mean;
        this_chunk -> start . variance = left_end . variance;
        this_chunk -> end . mean = right_end . mean;
        this_chunk -> end . variance = right_end . variance;
       }

#if  CHECK_CELSIM_COORDS
   if  ((Scaffold_Start [scaff_id] <= Scaffold_End [scaff_id] && ! flipped)
            || (Scaffold_Start [scaff_id] > Scaffold_End [scaff_id] && flipped))
       {
        this_chunk -> sim_start = Chunk_Info [cid] . celsim_left;
        this_chunk -> sim_end = Chunk_Info [cid] . celsim_right;
       }
     else
       {
        this_chunk -> sim_start = Chunk_Info [cid] . celsim_right;
        this_chunk -> sim_end = Chunk_Info [cid] . celsim_left;
       }
   if  (Scaffold_Start [scaff_id] <= Scaffold_End [scaff_id])
       {
        this_chunk -> sim_start -= Scaffold_Start [scaff_id];
        this_chunk -> sim_end -= Scaffold_Start [scaff_id];
       }
     else
       {
        this_chunk -> sim_start
            = Scaffold_Start [scaff_id] - this_chunk -> sim_start;
        this_chunk -> sim_end
            = Scaffold_Start [scaff_id] - this_chunk -> sim_end;
       }
#endif

   return  TRUE;
  }



static int  Between
    (double a, double b, double lo, double hi)

//  Return  TRUE  iff both  a  and  b  are in the closed interval
//   [lo, hi] .

  {
   return  (lo <= a && a <= hi && lo <= b && b <= hi);
  }



static void  Build_Path_Subgraph
    (int start_sub, int target_sub, Gap_Chunk_t * node [], int num_nodes,
     int forward_edge [], int reverse_edge [], Stone_Edge_t pool [],
     int sorted [], int * num_sorted)

//  Extract a topologically sorted set of nodes from  node [start_sub]  to
//  node [target_sub]  based on edges in  forward_edge []  and  reverse_edge []
//  (which refer to  pool [] ).  Put subscripts of resulting nodes,
//  in order, into  sorted []  and set  (* num_sorted)  to the number
//  of entries placed.  All references are subscripts into  node []
//  which has  num_nodes  entries.

  {
   Topologically_Sort (node, num_nodes, start_sub, forward_edge, pool,
                       sorted, num_sorted, TRUE);

#if  SHOW_STONE_CONFIRM
{
 int  i;

 fprintf (stderr,
          "Build_Path_Subgraph Topological Sort:  num_nodes = %d  num_sorted = %d\n",
          num_nodes, (* num_sorted));
 for  (i = 0;  i < (* num_sorted);  i ++)
   {
    int  cid = node [sorted [i]] -> chunk_id;

    fprintf (stderr, "%5d  %s\n", cid,
             sorted [i] != start_sub && sorted [i] != target_sub
               && REF (cid) . is_singleton ? "*UNIQUE*" : "");
   }
}
#endif

   return;
  }



static int  By_Index
    (const void * a, const void * b)

//  First compare  keep  flags of  a  and  b .  If not equal
//  return the value that indicates the  keep = TRUE  value
//  is first.
//  Otherwise, compare  index  field  in  a  and  b  (as  Gap_Chunk_t 's)
//  and return  -1  if  a  is before  b ,  0  if  a  and
//  b  are in the same place , and  1  if  b  is before  a .
//  Used for  qsort .

  {
   Gap_Chunk_t  * x, * y;

   x = (Gap_Chunk_t *) a;
   y = (Gap_Chunk_t *) b;

   if  (x -> keep && ! y -> keep)
       return  -1;
   else if  (y -> keep && ! x -> keep)
       return  1;
   else if  (x -> index < y -> index)
       return  -1;
   else if  (x -> index > y -> index)
       return  1;

   return  0;
  }



static int  By_High_Placement
    (const void * a, const void * b)

//  Compare the higher of the two positions ( start  and  end )
//  in  a  and  b  as  (Placement_t *) 's and return  -1  if  a  is
//  before  b ,  0  if  a  and  b  are in the same place , and  1
//  if  b  is before  a .
//  Used for  qsort .

  {
   Placement_t  * x, * y;
   double  x_hi, y_hi;

   x = (Placement_t *) a;
   y = (Placement_t *) b;

   x_hi = Max_double (x -> A_end . mean, x -> B_end . mean);
   y_hi = Max_double (y -> A_end . mean, y -> B_end . mean);

   if  (x_hi < y_hi)
       return  -1;
   else if  (x_hi > y_hi)
       return  1;

   return  0;
  }



static int  By_High_Position
    (const void * a, const void * b)

//  Compare the higher of the two positions ( start  and  end )
//  in  a  and  b  as  (Gap_Chunk_t * *) 's and return  -1  if  a  is
//  before  b ,  0  if  a  and  b  are in the same place , and  1
//  if  b  is before  a .
//  Used for  qsort .

  {
   Gap_Chunk_t  * * x, * * y;
   double  x_hi, y_hi;

   x = (Gap_Chunk_t * *) a;
   y = (Gap_Chunk_t * *) b;

   x_hi = Max_double ((* x) -> start . mean, (* x) -> end . mean);
   y_hi = Max_double ((* y) -> start . mean, (* y) -> end . mean);

   if  (x_hi < y_hi)
       return  -1;
   else if  (x_hi > y_hi)
       return  1;

   return  0;
  }



static int  By_Interval_Lo
    (const void * a, const void * b)

//  Compare the  lo  fields  in  a  and  b  as  (Interval_t *) 's and
//  return  -1  if  a  is before  b ,  0  if  a  and  b  are in the same
//  place , and  1  if  b  is before  a .
//  Used for  qsort .

  {
   Interval_t  * x, * y;

   x = (Interval_t *) a;
   y = (Interval_t *) b;

   if  (x -> lo < y -> lo)
       return  -1;
   else if  (x -> lo > y -> lo)
       return  1;

   return  0;
  }



static int  By_Keep_And_Low_Position
    (const void * a, const void * b)

//  First check the  keep  flags in  a  and  b  as  (Gap_Chunk_t *) 's
//  and return  -1  if  a  is  keep  and  b  is not, or  +1
//  if  b  is  keep  and  a  is not.  Otherwise,
//  compare the lower of the two positions ( start  and  end )
//  in  a  and  b  and return  -1  if  a  is
//  before  b ,  0  if  a  and  b  are in the same place , and  1
//  if  b  is before  a .
//  Used for  qsort .

  {
   Gap_Chunk_t  * x, * y;
   double  x_lo, y_lo;

   x = (Gap_Chunk_t *) a;
   y = (Gap_Chunk_t *) b;

   if  (x -> keep && ! y -> keep)
       return  -1;
   if  (! x -> keep && y -> keep)
       return  1;

   if  (x -> start . mean < x -> end . mean)
       x_lo = x -> start . mean;
     else
       x_lo = x -> end . mean;

   if  (y -> start . mean < y -> end . mean)
       y_lo = y -> start . mean;
     else
       y_lo = y -> end . mean;

   if  (x_lo < y_lo)
       return  -1;
   else if  (x_lo > y_lo)
       return  1;

   return  0;
  }



static int  By_Low_Position
    (const void * a, const void * b)

//  Compare the lower of the two positions ( start  and  end )
//  in  a  and  b  as  (Gap_Chunk_t * *) 's and return  -1  if  a  is
//  before  b ,  0  if  a  and  b  are in the same place , and  1
//  if  b  is before  a .  Use  variance  to break ties.
//  Used for  qsort .

  {
   Gap_Chunk_t  * * x, * * y;
   LengthT  * x_lo, * y_lo;

   x = (Gap_Chunk_t * *) a;
   y = (Gap_Chunk_t * *) b;

   if  ((* x) -> start . mean < (* x) -> end . mean)
       x_lo = & ((* x) -> start);
     else
       x_lo = & ((* x) -> end);

   if  ((* y) -> start . mean < (* y) -> end . mean)
       y_lo = & ((* y) -> start);
     else
       y_lo = & ((* y) -> end);

   if  (x_lo -> mean < y_lo -> mean)
       return  -1;
   else if  (x_lo -> mean > y_lo -> mean)
       return  1;
   else if  (x_lo -> variance < y_lo -> variance)
       return  -1;
   else if  (x_lo -> variance > y_lo -> variance)
       return  1;

   return  0;
  }



static int  By_Scaff_And_Flipped
    (const void * a, const void * b)

//  Regard  (* a)  and  (* b)  as  Stack_Entry_t 's  and
//  return their order with  scaff_id  as the major key
//  and  flipped  as the minor key.
//  Used for  qsort .

  {
   Stack_Entry_t  * x, * y;

   x = (Stack_Entry_t *) a;
   y = (Stack_Entry_t *) b;

   if  (x -> scaff_id < y -> scaff_id)
       return  -1;
   else if  (x -> scaff_id > y -> scaff_id)
       return  1;
   else if  (x -> flipped < y -> flipped)
       return  -1;
   else if  (x -> flipped > y -> flipped)
       return  1;

   return  0;
  }



static int  By_Scaff_Flipped_And_Left_End
    (const void * a, const void * b)

//  Regard  (* a)  and  (* b)  as  Stack_Entry_t 's  and
//  return their order with  scaff_id  as the major key,
//  and  flipped  as the second key, and  left_end . mean
//  (less 3.0 time std dev) as the third key.
//  Used for  qsort .

  {
   Stack_Entry_t  * x, * y;
   double  x_left, y_left;

   x = (Stack_Entry_t *) a;
   y = (Stack_Entry_t *) b;
   x_left = x -> left_end . mean - 3.0 * sqrt (x -> left_end . variance);
   y_left = y -> left_end . mean - 3.0 * sqrt (y -> left_end . variance);

   if  (x -> scaff_id < y -> scaff_id)
       return  -1;
   else if  (x -> scaff_id > y -> scaff_id)
       return  1;
   else if  (x -> flipped < y -> flipped)
       return  -1;
   else if  (x -> flipped > y -> flipped)
       return  1;
   else if  (x_left < y_left)
       return  -1;
   else if  (x_left > y_left)
       return  1;

   return  0;
  }



static void  Calc_End_Coords
    (Stack_Entry_t * stack, int stack_top,
     LengthT * left_end, LengthT * right_end, ChunkInstanceT * chunk,
     double ref_variance)

//  Set  (* left_end)  and  (* right_end)  to best estimates of
//  the ends  (* chunk)  based on the edges in
//  stack [0 .. (stack_top - 1)]  with variances scaled down to
//  ref_variance .

  {
   double  left_numerator_sum, left_denom_sum;
   double  right_numerator_sum, right_denom_sum;
   double  var;
   int  i;

   left_numerator_sum = 0.0;
   left_denom_sum = 0.0;
   right_numerator_sum = 0.0;
   right_denom_sum = 0.0;

   assert (stack_top < STACK_SIZE);
   for  (i = 0;  i < stack_top;  i ++)
     {
      if  (stack [i] . is_bad)
          continue;

      var = fabs (stack [i] . source_variance - ref_variance)
                + stack [i] . edge -> distance . variance;
      assert (var > 0.0);
      if  (stack [i] . left_link)
          {
           left_denom_sum += 1.0 / var;
           left_numerator_sum += stack [i] . left_end . mean / var;
           stack [i] . left_end . variance = var;
           var += chunk -> bpLength . variance;
           right_denom_sum += 1.0 / var;
           right_numerator_sum += stack [i] . right_end . mean / var;
           stack [i] . right_end . variance = var;
          }
        else
          {
           right_denom_sum += 1.0 / var;
           right_numerator_sum += stack [i] . right_end . mean / var;
           stack [i] . right_end . variance = var;
           var += chunk -> bpLength . variance;
           left_denom_sum += 1.0 / var;
           left_numerator_sum += stack [i] . left_end . mean / var;
           stack [i] . left_end . variance = var;
          }
if  (Global_Debug_Flag)
    fprintf (stderr, "source_var = %.0f  ref_var = %.0f  edge_var = %.0f  var = %.0f\n",
             stack [i] . source_variance, ref_variance,
             stack [i] . edge -> distance . variance, var);
     }
if  (Global_Debug_Flag)
    fprintf (stderr, "left_end_var = %.0f\n", 1.0 / left_denom_sum);

   left_end -> mean = left_numerator_sum / left_denom_sum;
   left_end -> variance = 1.0 / left_denom_sum;
   right_end -> mean = right_numerator_sum / right_denom_sum;
   right_end -> variance = 1.0 / right_denom_sum;

   return;
  }



char *  CGB_Type_As_String
    (unsigned int t)

//  Return string equivalent of CGB_Type  t

  {
   switch  (t)
     {
      case  UU_CGBTYPE :
        return  "uu";
      case  UR_CGBTYPE :
        return  "ur";
      case  RU_CGBTYPE :
        return  "ru";
      case  RR_CGBTYPE :
        return  "rr";
      case  XX_CGBTYPE :
        return  "xx";
     }

   return  "*???*";
  }



static void  Adjust_By_Ref_Variance
    (Scaffold_Fill_t * fill_chunks)

//  Add the  ref_variance  of each gap to the variances
//  of all chunks with the  keep  flag true for each gap
//  in  fill_chunks .

  {
   int  scaff_id;

   fprintf (stderr, "### Adjust_By_Ref_Variance ###\n");

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     Adjust_By_Ref_Variance_One_Scaffold (fill_chunks, scaff_id);

   return;
  }



static void  Adjust_By_Ref_Variance_One_Scaffold
    (Scaffold_Fill_t * fill_chunks, int scaff_id)

//  Add the  ref_variance  of each gap to the variances
//  of all chunks with the  keep  flag true for each gap
//  in  fill_chunks [scaff_id] .

  {
   int  j;

   for  (j = GAPS_TO_ADJUST)
     {
      int  k;
      Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;


      for  (k = 0;  k < this_gap -> num_chunks;  k ++)
        {
         Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

//  Adjust them all so all entries printed out will be same scale
//            if  (this_chunk -> keep)
             {
              this_chunk -> start . variance
                  += this_gap -> ref_variance;
              this_chunk -> end . variance
                  += this_gap -> ref_variance;
             }
        }
     }

   return;
  }



static void  Adjust_Positions
    (Gap_Fill_t * this_gap, int num_targets, Target_Info_t target [],
     int gap_sub [], int max_hits, int max_first, double start_coord,
     DirectionType direction, double * high_variance, FILE * fp)

//  Adjust positions of chunks in  this_gap  that were confirmed
//  by an overlap path.  Use the length of the overlap path
//  to make the adjustment.
//  The path information is in  target .
//  num_targets  is the total number of entries in  target .
//  gap_sub  connects targets to gap entries, i.e.,  target [i]  is
//  this_gap -> chunk [gap_sub [i]] .
//  max_hits  is the number of targets to which an overlap path was found
//  max_first  is the target subscript of the first of those max hits.
//  start_coord  is the position from which the path originated.
//  direction  is the direction the walk was taken w.r.t. scaffold coordinates.
//  Thus, for  AS_FORWARD  direction the scaffold coordinates should be going
//  the same direction as the  where  values; for  AS_REVERSE  direction they
//  should be in opposite directions.
//  Set  (* high_variance)  to the highest adjust variance of any target
//  on the path.
//  fp  is the log file.
//
//  In general, the  where  coordinates  should be ascending, but in
//  the case of a containment overlap, they may decrease.

  {
   (* high_variance) = 0.0;

   if  (max_hits > 0)
       {
        int  i;
        int  sub = max_first;

        for  (i = 0;  i < max_hits;  i ++)
          {
           Gap_Chunk_t  * gap_chunk;
           ChunkInstanceT  * chunk;
           double  chunk_len, factor;
           LengthT  new_start, new_end;

           gap_chunk = this_gap -> chunk + gap_sub [sub];
           gap_chunk -> index = i;

           chunk = GetGraphNode (ScaffoldGraph -> RezGraph, target [sub] . id);
           chunk_len = chunk -> bpLength . mean;

fprintf (fp, "Adjust:  %d  was (%.0f [%.0f], %.0f [%.0f])",
         chunk -> id,
         gap_chunk -> start . mean,
         gap_chunk -> start . variance,
         gap_chunk -> end . mean,
         gap_chunk -> end . variance);

           switch  (direction)
             {
              case  AS_FORWARD :
                factor = 1.0;
                break;
              case  AS_REVERSE :
                factor = -1.0;
                break;
              default :
                fprintf (stderr, "ERROR:  Bad direction\n");
                assert (FALSE);
             }
           switch  (target [sub] . orient)
             {
              case  AB_AB :
              case  BA_AB :
                new_end . mean = start_coord + factor * target [sub] . where;
                new_end . variance = ComputeFudgeVariance(target [sub] . total); //  * VARIANCE_FUDGE_FACTOR;
                new_start . mean = new_end . mean - factor * chunk_len;
                new_start . variance
		  = ComputeFudgeVariance(target [sub] . total - chunk_len);//  * VARIANCE_FUDGE_FACTOR;
                break;

              case  AB_BA :
              case  BA_BA :
                new_start . mean = start_coord + factor * target [sub] . where;
                new_start . variance = ComputeFudgeVariance(target [sub] . total); // * VARIANCE_FUDGE_FACTOR;
                new_end . mean = new_start . mean - factor * chunk_len;
                new_end . variance
		  = ComputeFudgeVariance(target [sub] . total - chunk_len); //  * VARIANCE_FUDGE_FACTOR;
                break;

              default :
               fprintf (stderr, "ERROR:  Bad target orientation\n");
               assert (FALSE);
             }

#if  SIMPLE_ADJUST
           gap_chunk -> start = new_start;
           gap_chunk -> end = new_end;
#else
           {
            double  denom;

            if  (new_start . variance < MIN_VARIANCE)
                new_start . variance = MIN_VARIANCE;
            if  (gap_chunk -> start . variance < MIN_VARIANCE)
                gap_chunk -> start . variance = MIN_VARIANCE;

            denom = 1.0 / new_start . variance
                      + 1.0 / gap_chunk -> start . variance;
            gap_chunk -> start . mean
                = (gap_chunk -> start . mean
                        / gap_chunk -> start . variance
                    + new_start . mean / new_start . variance) / denom;
            gap_chunk -> start . variance = 1.0 / denom;

            if  (new_end . variance < MIN_VARIANCE)
                new_end . variance = MIN_VARIANCE;
            if  (gap_chunk -> end . variance < MIN_VARIANCE)
                gap_chunk -> end . variance = MIN_VARIANCE;

            denom = 1.0 / new_end . variance
                      + 1.0 / gap_chunk -> end . variance;
            gap_chunk -> end . mean
                = (gap_chunk -> end . mean
                        / gap_chunk -> end . variance
                    + new_end . mean / new_end . variance) / denom;
            gap_chunk -> end . variance = 1.0 / denom;
           }
#endif
           if  (gap_chunk -> start . variance < MIN_VARIANCE)
               gap_chunk -> start . variance = MIN_VARIANCE;
           if  (gap_chunk -> end . variance < MIN_VARIANCE)
               gap_chunk -> end . variance = MIN_VARIANCE;
           if  (gap_chunk -> start . variance > (* high_variance))
               (* high_variance) = gap_chunk -> start . variance;
           if  (gap_chunk -> end . variance > (* high_variance))
               (* high_variance) = gap_chunk -> end . variance;

fprintf (fp, "  now (%.0f [%.0f], %.0f [%.0f])\n",
         gap_chunk -> start . mean,
         gap_chunk -> start . variance,
         gap_chunk -> end . mean,
         gap_chunk -> end . variance);

           sub = target [sub] . next;
          }
       }

   return;
  }



static void  Check_Olaps
    (Gap_Fill_t * gap)

//  Check if there are overlaps in the hash table to support all
//  overlaps indicated by positions of chunks in  gap .  Only
//  do chunks whose  keep  flag  is true.  Adjust positions if
//  the overlap indicates chunks are in the wrong relative order.

  {
   Placement_t  * place;
   ChunkInstanceT  * chunk;
   int  i, j, ct;

PALLOC ((2 + gap -> num_chunks) * sizeof (Placement_t));
   place = (Placement_t *) safe_malloc
             ((2 + gap -> num_chunks) * sizeof (Placement_t));
   ct = 0;
   if  (gap -> left_cid > 0)
       {
        chunk = GetGraphNode (ScaffoldGraph -> RezGraph,
                              gap -> left_cid);
        place [0] . A_end = chunk -> offsetAEnd;
        place [0] . B_end = chunk -> offsetBEnd;
        place [0] . A_end . variance
            = fabs (chunk -> offsetAEnd . variance
                      - gap -> ref_variance);
        place [0] . B_end . variance
            = fabs (chunk -> offsetBEnd . variance
                      - gap -> ref_variance);
        place [0] . id = gap -> left_cid;
        place [0] . keep = TRUE;
        if  (place [0] . B_end . mean < place [0] . A_end . mean)
            place [0] . flipped = TRUE;
          else
            place [0] . flipped = FALSE;
        place [0] . chunk = NULL;
        ct ++;
       }

   for  (i = 0;  i < gap -> num_chunks;  i ++)
     if  (gap -> chunk [i] . keep)
         {
          place [ct] . A_end = gap -> chunk [i] . start;
          place [ct] . B_end = gap -> chunk [i] . end;
          place [ct] . id = gap -> chunk [i] . chunk_id;
          place [ct] . keep = TRUE;
          if  (place [ct] . B_end . mean < place [ct] . A_end . mean)
              place [ct] . flipped = TRUE;
            else
              place [ct] . flipped = FALSE;
          place [ct] . chunk = gap -> chunk + i;
          ct ++;
         }

   if  (gap -> right_cid > 0)
       {
        chunk = GetGraphNode (ScaffoldGraph -> RezGraph,
                              gap -> right_cid);
        place [ct] . A_end = chunk -> offsetAEnd;
        place [ct] . B_end = chunk -> offsetBEnd;
        place [ct] . A_end . variance
            = fabs (chunk -> offsetAEnd . variance
                      - gap -> ref_variance);
        place [ct] . B_end . variance
            = fabs (chunk -> offsetBEnd . variance
                      - gap -> ref_variance);
        place [ct] . id = gap -> right_cid;
        place [ct] . keep = TRUE;
        if  (place [ct] . B_end . mean < place [ct] . A_end . mean)
            place [ct] . flipped = TRUE;
          else
            place [ct] . flipped = FALSE;
        place [ct] . chunk = NULL;
        ct ++;
       }

   qsort (place, ct, sizeof (Placement_t), By_High_Placement);

   for  (i = 1;  i < ct;  i ++)
     for  (j = i - 1;  j >= 0 && place [i] . keep;  j --)
       {
        ChunkOrientationType  orient;
        double  how_much, allowed_error;
        int  min_ahang, max_ahang;

        how_much = Interval_Intersection
                       ((int) Min_double (place [j] . A_end . mean,
                                          place [j] . B_end . mean),
                        (int) Max_double (place [j] . A_end . mean,
                                          place [j] . B_end . mean),
                        (int) Min_double (place [i] . A_end . mean,
                                          place [i] . B_end . mean),
                        (int) Max_double (place [i] . A_end . mean,
                                          place [i] . B_end . mean));
        if  (gap -> has_path)
            allowed_error = 30.0 + CGW_FUDGE_FACTOR * how_much;
          else
            allowed_error = 3.0 * sqrt (Max_double (place [j] . A_end . variance,
                           place [i] . A_end . variance));

        if  (place [j] . keep
#if  1
               && Might_Overlap
                      (place [j] . A_end . mean, place [j] . B_end . mean,
                       place [i] . A_end . mean, place [i] . B_end . mean,
                       - allowed_error,     
                       & orient, & min_ahang, & max_ahang))
                            // - allowed_error really makes this "must overlap"
#else
               && Should_Overlap (place + j, place + i, & orient, & how_much))
#endif
            {
#if 0
             ChunkOrientationType  new_orient;
             int  olap_found;
#endif
//             ChunkOverlapCheckT  olap;
             Overlap  * olap;
             char  * i_seq, * j_seq;

#if  1
             i_seq = Get_Contig_Sequence (place [i] . id);
             j_seq = Get_Contig_Sequence (place [j] . id);

             olap = OverlapSequences
                        (j_seq, i_seq, orient, min_ahang - (int) (3.0 * allowed_error),
                         max_ahang + (int) (3.0 * allowed_error),
                         CGW_DP_ERATE, CGW_DP_THRESH, CGW_DP_MINLEN,
                         AS_FIND_ALIGN);			 
#else
             olap = OverlapChunks                 // debug code does NOT handle suspicious
                        (ScaffoldGraph -> RezGraph,
                         place [j] . id, place [i] . id,
                         orient,
                         (int) (how_much - allowed_error),
                         (int) (how_much + allowed_error),
                         CGW_DP_ERATE, FALSE);

             olap_found = (olap . overlap > 0);
#endif

#ifdef DEBUG_DETAILED
         fprintf (stderr, ">>> Check olap %d/%d %s  exp: [%.0f,%.0f]  %s",
         place [j] . id, place [i] . id, Orientation_As_String (orient),
         how_much - allowed_error, how_much + allowed_error,
         gap -> has_path ? "Has Path" : "No Path");
#endif
         
#if  VERBOSE
fprintf (stderr,
         "olapping j = %d  [%.0f,%.0f]   i = %d  [%.0f,%.0f]\n"
         "  min_ahang = %d  max_ahang = %d  allowed_error = %d\n"
         "  orient = %s  olap = %p\n",
         place [j] . id, 
         place [j] . A_end . mean, place [j] . B_end . mean,
         place [i] . id,
         place [i] . A_end . mean, place [i] . B_end . mean,
         min_ahang, max_ahang,
         (int) allowed_error, Orientation_As_String (orient),
         olap);
if  (olap == NULL)
    fprintf (stderr, "  Not found\n");
  else
    fprintf (stderr, "  begpos = %d  endpos = %d  length = %d\n",
             olap -> begpos, olap -> endpos, olap -> length);
#endif

#if  0
if  (olap_found)
    {
#ifdef DEBUG_DETAILED
     fprintf (stderr, "  Found  overlap = %d\n", olap . overlap);
#endif
    }
  else
    {
     ChunkOverlapCheckT  rev_olap;
     int  rev_olap_found;
     int  i_len, j_len;

     switch  (orient)
       {
        case  AB_AB :
        case  BA_BA :
          new_orient = orient;
          break;
        case  AB_BA :
          new_orient = BA_AB;
          break;
        case  BA_AB :
          new_orient = AB_BA;
          break;
        default :
          fprintf (stderr, "YIKES:  Bad orientation = %d\n", (int) orient);
          assert (FALSE);
       }

     i_len = (int) fabs (place [i] . A_end . mean - place [i] . B_end . mean);
     j_len = (int) fabs (place [j] . A_end . mean - place [j] . B_end . mean);
     rev_olap = OverlapChunks      // debug code does NOT handle suspicious
                (ScaffoldGraph -> RezGraph,
                 place [i] . id, place [j] . id,
                 new_orient,
                 0, Max_int (i_len, j_len),
                 CGW_DP_ERATE, FALSE);

     rev_olap_found = (rev_olap . overlap > 0);
     fprintf (stderr, "  Didn't find expected overlap between %d and %d",
              place [j] . id, place [i] . id);
     if  (rev_olap_found)
         fprintf (stderr, "  Reverse found  overlap = %d\n", rev_olap . overlap);
       else
         fprintf (stderr, "  Reverse not found\n");
    }
#ifdef DEBUG_DETAILED
fprintf (stderr, "     AContainsB = %c  BContainsA = %c\n",
         olap . AContainsB ? 'T' : 'F', olap . BContainsA ? 'T' : 'F');
#endif
#endif

             if  (olap == NULL)
                 {
                  if  (place [i] . chunk == NULL)
                      place [j] . keep = FALSE;
                    else
                      place [i] . keep = FALSE;
                 }
             free (i_seq);
             free (j_seq);
            }
       }

   for  (i = 0;  i < ct;  i ++)
     if  (! place [i] . keep && place [i] . chunk != NULL)
         {
fprintf (stderr, ">>> Rejected chunk %d\n", place [i] . id);
          place [i] . chunk -> keep = FALSE;
         }

   free (place);

   return;
  }



static void  Check_Other_Links
    (Scaffold_Fill_t * fill_chunks)

//  Check if entries in  fill_chunks  have any clone mate links to
//  chunks that are neither unique nor in  fill_chunks.

  {
   int  have_bad_ct = 0, no_bad_ct = 0;
   int  keep_ct = 0, reject_ct = 0;
   double  denom;
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     Check_Other_Links_One_Scaffold
         (fill_chunks, scaff_id, & have_bad_ct, & no_bad_ct, & keep_ct, & reject_ct);

   denom = have_bad_ct + no_bad_ct;
   if  (denom == 0.0)
       denom = 1.0;
     
   fprintf (stderr, "        Removed by check links: %7d\n", reject_ct);
   fprintf (stderr, "     Sent to consistency check: %7d\n", keep_ct);
   fprintf (stderr,
            "                Have bad links: %7d (%.1f%%)\n"
            "                  No bad links: %7d (%.1f%%)\n",
            have_bad_ct, (100.0 * have_bad_ct) / denom,
            no_bad_ct, (100.0 * no_bad_ct) / denom);

   return;
  }



static void  Check_Other_Links_One_Scaffold
    (Scaffold_Fill_t * fill_chunks, int scaff_id,
     int * have_bad_ct, int * no_bad_ct, int * keep_ct, int * reject_ct)

//  Check if entries in  fill_chunks [scaff_id]  have any clone mate links to
//  chunks that are neither unique nor in  fill_chunks.  Increment values
//  of  (* have_bad_ct) ,  (* no_bad_ct) ,  (* keep_ct)  and  (* reject_ct)
//  for entries with (resp.) at least one bad link, no bad links,
//  that will be kept, and that will be rejected.

  {
   int  j;

   for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
     {
      int  k, sub;

      for  (k = sub = 0;  k < fill_chunks [scaff_id] . gap [j] . num_chunks;  k ++)
        {
         ChunkInstanceT  * chunk;
         GraphEdgeIterator  ci_edges;
         CIEdgeT  * edge;
         int  cid, good_links, bad_links;

         cid = fill_chunks [scaff_id] . gap [j] . chunk [k] . chunk_id;
         chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);

         InitGraphEdgeIterator (ScaffoldGraph->RezGraph, cid, ALL_END, ALL_EDGES,
                              GRAPH_EDGE_DEFAULT, & ci_edges);

         good_links = bad_links = 0;
         while  ((edge = NextGraphEdgeIterator (& ci_edges)) != NULL)
           {
            int  clone_mates;

            if  (isProbablyBogusEdge (edge)
                   || isSloppyEdge (edge))
                continue;

            clone_mates = edge -> edgesContributing;
            if  (edge -> flags . bits . isPossibleChimera)
                clone_mates = 1;
              else
                {
                  if  (isOverlapEdge(edge)){
                     clone_mates --;
                  }
                }

            if  (clone_mates > 0)
                {
                 int  other_cid;

                 if  (edge -> idA == cid)
                     other_cid = edge -> idB;
                   else
                     other_cid = edge -> idA;

                 if  (Chunk_Info [other_cid] . colour != UNIQUE_COLOUR
                        && Chunk_Info [other_cid] . colour != PLACED_COLOUR
                        && Chunk_Info [other_cid] . colour != CONSISTENT_COLOUR)
                     bad_links += clone_mates;
                   else
                     good_links += clone_mates;
                }
           }

#if  0
         printf ("Chunk %6d:  good links = %3d  bad links = %3d\n",
                 cid, good_links, bad_links);
#endif
         if  (bad_links == 0)
             (* no_bad_ct) ++;
           else
             (* have_bad_ct) ++;

         if  (bad_links == 0
                || (bad_links == 1 && good_links >= GOOD_LINKS_IF_BAD)
             )
             {
              fill_chunks [scaff_id] . gap [j] . chunk [sub ++]
                = fill_chunks [scaff_id] . gap [j] . chunk [k];
              (* keep_ct) ++;
             }
           else
             {
              Chunk_Info [cid] . colour = CONSISTENT_COLOUR;
              (* reject_ct) ++;
             }
        }

      fill_chunks [scaff_id] . gap [j] . num_chunks = sub;
     }

   return;
  }



static void  Check_Rock_Olaps
    (Gap_Fill_t * gap)

//  Check for dubious overlaps of chunks in  gap .

  {
   Placement_t  place1, place2;
   ChunkInstanceT  * chunk;
   int  i;

   if  (gap -> left_cid > 0)
       {
        chunk = GetGraphNode (ScaffoldGraph -> RezGraph,
                              gap -> left_cid);
        place1 . A_end = chunk -> offsetAEnd;
        place1 . B_end = chunk -> offsetBEnd;
        place1 . id = gap -> left_cid;
        place1 . keep = TRUE;
        if  (place1 . B_end . mean < place1 . A_end . mean)
            place1 . flipped = TRUE;
          else
            place1 . flipped = FALSE;
        place1 . chunk = NULL;

        for  (i = 0;  i < gap -> num_chunks;  i ++)
          if  (gap -> chunk [i] . keep)
              {
               ChunkOrientationType  orient;
               double  how_much;

               place2 . A_end = gap -> chunk [i] . start;
               place2 . B_end = gap -> chunk [i] . end;
               place2 . id = gap -> chunk [i] . chunk_id;
               place2 . keep = TRUE;
               if  (place2 . B_end . mean < place2 . A_end . mean)
                   place2 . flipped = TRUE;
                 else
                   place2 . flipped = FALSE;
               place2 . chunk = gap -> chunk + i;

               if  (Max_double (place2 . A_end . mean, place2 . B_end . mean)
                      <= Max_double (place1 . A_end . mean, place1 . B_end . mean))
                   gap -> chunk [i] . keep = FALSE;
               else if  (Should_Overlap (& place1, & place2, & orient, & how_much))
                   {
                    ChunkOverlapCheckT  olap;
                    int  olap_found;
                    double  allowed_error;

                    allowed_error = 30.0 + CGW_FUDGE_FACTOR * how_much;

                    olap_found = LookupOverlap (ScaffoldGraph -> RezGraph,
                                                place1 . id,
                                                place2 . id,
                                                orient, & olap);
                    if  (olap_found)
                        {
                         gap -> chunk [i] . start . mean -= olap . overlap - how_much;
                         gap -> chunk [i] . end . mean -= olap . overlap - how_much;
                        }
                      else
                        gap -> chunk [i] . keep = FALSE;
                   }
              }
       }

   if  (gap -> right_cid > 0)
       {
        chunk = GetGraphNode (ScaffoldGraph -> RezGraph,
                              gap -> right_cid);
        place1 . A_end = chunk -> offsetAEnd;
        place1 . B_end = chunk -> offsetBEnd;
        place1 . id = gap -> right_cid;
        place1 . keep = TRUE;
        if  (place1 . B_end . mean < place1 . A_end . mean)
            place1 . flipped = TRUE;
          else
            place1 . flipped = FALSE;
        place1 . chunk = NULL;

        for  (i = 0;  i < gap -> num_chunks;  i ++)
          if  (gap -> chunk [i] . keep)
              {
               ChunkOrientationType  orient;
               double  how_much;

               place2 . A_end = gap -> chunk [i] . start;
               place2 . B_end = gap -> chunk [i] . end;
               place2 . id = gap -> chunk [i] . chunk_id;
               place2 . keep = TRUE;
               if  (place2 . B_end . mean < place2 . A_end . mean)
                   place2 . flipped = TRUE;
                 else
                   place2 . flipped = FALSE;
               place2 . chunk = gap -> chunk + i;

               if  (Min_double (place2 . A_end . mean, place2 . B_end . mean)
                      >= Min_double (place1 . A_end . mean, place1 . B_end . mean))
                   gap -> chunk [i] . keep = FALSE;
               else if  (Should_Overlap (& place2, & place1, & orient, & how_much))
                   {
                    ChunkOverlapCheckT  olap;
                    int  olap_found;
                    double  allowed_error;

                    allowed_error = 30.0 + CGW_FUDGE_FACTOR * how_much;

                    olap_found = LookupOverlap (ScaffoldGraph -> RezGraph,
                                                place2 . id,
                                                place1 . id,
                                                orient, & olap);
                    if  (olap_found)
                        {
                         gap -> chunk [i] . start . mean += olap . overlap - how_much;
                         gap -> chunk [i] . end . mean += olap . overlap - how_much;
                        }
                      else
                        gap -> chunk [i] . keep = FALSE;
                   }
              }
       }

   return;
  }


static void  Check_Rocks
    (FILE * fp, Scaffold_Fill_t * fill_chunks)

//  Check all rocks in  fill_chunks  for:
//    - overlaps with scaffold chunks on ends of gaps.  If have and
//      compatible with position range, then change position to that
//      indicated by the overlap
//    - overlaps with each other that indicate out of order.  If any
//      change positions or discard both for now.

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = GAPS_TO_ADJUST)
        {
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;

         Check_Rock_Olaps (this_gap);
        }
     }

   return;
  }




static int  Check_Scaffold_and_Orientation
    (int cid, Stack_Entry_t * stack, int stack_top, int * good_total,
     Scaffold_Fill_t * fill_chunks, int * bad_links,
     int  min_good_links)

//  Check if edges in  stack [0 .. (stack_top - 1)] all go to
//  same scaffold and have same orientation.  Return  TRUE  if they
//  do,  FALSE  otherwise.  Set  (* good_total)  to the number
//  of mate_links included in all the edges.
//  cid  is the id of the current chunk.
//  Set  (* bad_links) to the number of bad links found.
//  fill_chunks  contains scaffold information.
//  min_good_links  is the minimum number of link mates to scaffold
//  chunks needed to approve this chunk.
//  If there are links to multiple scaffolds, save them in global
//  Scaff_Join  for later consistency checking.

  {
   int  scaff_id [MAX_DIFF_SCAFFS], scaff_links [MAX_DIFF_SCAFFS] = {0};
   int  flipped_ct, non_flipped_ct;
   int  i, j, num_scaffs;

   (* good_total) = 0;
   (* bad_links) = 0;


   // Count number of links to each scaffold

   num_scaffs = 0;
   for  (i = 0;  i < stack_top;  i ++)
     {
      for  (j = 0;  j < num_scaffs;  j ++)
        if  (REF (stack [i] . chunk_id) . scaff_id == scaff_id [j])
            {
             scaff_links [j] += stack [i] . num_good_mates;
             break;
            }
      if  (j == num_scaffs)
          {
           if  (num_scaffs == MAX_DIFF_SCAFFS)
               return  FALSE;
           scaff_id [j] = REF (stack [i] . chunk_id) . scaff_id;
           scaff_links [j] = stack [i] . num_good_mates;
           num_scaffs ++;
          }
     }

   // Sort to find scaffold with most links

   for  (i = 0;  i < num_scaffs - 1;  i ++)
     for  (j = i + 1;  j < num_scaffs;  j ++)
       if  (scaff_links [i] < scaff_links [j])
           {
            int  save;

            save = scaff_id [i];
            scaff_id [i] = scaff_id [j];
            scaff_id [j] = save;
            save = scaff_links [i];
            scaff_links [i] = scaff_links [j];
            scaff_links [j] = save;
           }

   switch  (num_scaffs)
     {
      case  0 :
        fprintf (stderr, "ERROR:  Check_Scaffold_and_Orientation has empty stack\n");
        assert (FALSE);

      case  1 :
        break;

      case  2 :
        if  (scaff_links [1] == 1 && scaff_links [0] >= GOOD_LINKS_IF_BAD)
            {
             (* bad_links) = 1;
             break;
            }

#if  USE_JOIN_SCAFFS
        //  Check if this chunk "joins" two scaffolds
        if  (! Contained_Only_Switch
               && scaff_links [0] >= min_good_links
               && scaff_links [1] >= min_good_links)
            {
             Check_Scaffold_Join
                 (cid, stack, stack_top, scaff_id, scaff_links,
                  fill_chunks, (* bad_links));
            }
#endif
        return  FALSE;

      case  3 :
#if  USE_JOIN_SCAFFS
        //  Check if this chunk "joins" two scaffolds
        if  (! Contained_Only_Switch
               && scaff_links [0] >= min_good_links
               && scaff_links [1] >= min_good_links
               && scaff_links [2] == 1)
            {
             (* bad_links) = 0;
             Check_Scaffold_Join
                 (cid, stack, stack_top, scaff_id, scaff_links,
                  fill_chunks, (* bad_links));
            }
#endif
        return  FALSE;

      default :
        assert (num_scaffs > 0);
        //  Might be separating two chunks--will check for that later.
        return  FALSE;
     }

   (* good_total) = scaff_links [0];

   // Mark good entries in stack
   for  (i = 0;  i < stack_top;  i ++)
     if  (REF (stack [i] . chunk_id) . scaff_id != scaff_id [0])
         stack [i] . is_bad = TRUE;

#if  0
   //**ALD
   // If multiple links to same contig, mark links
   // with fewer than max good mates as bad
   // Must be same contig if is good and rel_pos is same
   // This should help prevent collapsed tandems.
   for  (i = 0;  i < stack_top - 1;  i ++)
     {
      if  (stack [i] . is_bad)
          continue;
      for  (j = i + 1;  j < stack_top;  j ++)
        if  (! stack [j] . is_bad
                && REF (stack [i] . chunk_id) . rel_pos
                      == REF (stack [i] . chunk_id) . rel_pos)
            {
             if  (stack [i] . num_good_mates <= stack [j] . num_good_mates)
                 {
                  stack [i] . is_bad = TRUE;
                  (* bad_links) += stack [i] . num_good_mates;
                 }
               else
                 {
                  stack [j] . is_bad = TRUE;
                  (* bad_links) += stack [j] . num_good_mates;
                 }
            }
     }
#endif
   
   // Check for orientation discrepancy
   flipped_ct = non_flipped_ct = 0;
   for  (i = 0;  i < stack_top;  i ++)
     {
      if  (stack [i] . is_bad)
          continue;
      if  (stack [i] . flipped)
          flipped_ct += stack [i] . num_good_mates;
        else
          non_flipped_ct += stack [i] . num_good_mates;
     }

   if  (flipped_ct > 1 && non_flipped_ct > 1)
       return  FALSE;

   if  (flipped_ct == 1 && non_flipped_ct >= GOOD_LINKS_IF_BAD)
       {
        if  ((* bad_links) > 0)
            return  FALSE;
        for  (i = 0;  i < stack_top;  i ++)
          if  (stack [i] . flipped)
              stack [i] . is_bad = TRUE;
        (* bad_links) = 1;
        (* good_total) --;
       }
   else if  (non_flipped_ct == 1 && flipped_ct >= GOOD_LINKS_IF_BAD)
       {
        if  ((* bad_links) > 0)
            return  FALSE;
        for  (i = 0;  i < stack_top;  i ++)
          if  (! stack [i] . flipped)
              stack [i] . is_bad = TRUE;
        (* bad_links) = 1;
        (* good_total) --;
       }
   else if  (flipped_ct != 0 && non_flipped_ct != 0)
       return  FALSE;
   
   return  TRUE;
  }



static void  Check_Scaffold_Join
    (int cid, Stack_Entry_t * stack, int stack_top, int scaff_id [],
     int scaff_links [], Scaffold_Fill_t * fill_chunks, int bad_links)

//  See if chunk  cid  has evidence of joining scaffolds in
//  scaff_id [0 .. 1] .  Evidence is edges on  stack [0 .. (stack_top - 1)] .
//  scaff_links [0 .. 1]  are the number of mate links to respective
//  scaffolds.  If there is enough evidence, store the connection to
//  allow checking for later inconsistencies, i.e., if some other
//  chunk suggests a conflicting join.
//  fill_chunks  contains scaffold information.
//  bad_links is the number of faulty links already encountered--in
//  general only 1 is tolerated

  {
   Stack_Entry_t  save;
   ContigT  * contig;
   float  edge_quality1, edge_quality2;
   int  bad_allowed, consistent;
   LengthT  left_end1, right_end1, left_end2, right_end2;
   LengthT  a_end1, b_end1, a_end2, b_end2;
   int  gap1, gap2, new_scaff1, new_scaff2;
   int  flipped, not_flipped, min, max;
   int  group1, group2, cover_stat, link_ct;
   int  i, j;

   // Move  scaff_id [0]  edges to front of stack

   for  (i = j = 0;  j < stack_top;  j ++)
     {
      if  (REF (stack [j] . chunk_id) . scaff_id == scaff_id [0])
          {
           if  (i != j)
               {
                save  = stack [i];
                stack [i] = stack [j];
                stack [j] = save;
               }
           i ++;
          }
     }
   group1 = i;

   // Move  scaff_id [1]  edges in back of them

   for  (j = i;  j < stack_top;  j ++)
     {
      if  (REF (stack [j] . chunk_id) . scaff_id == scaff_id [1])
          {
           if  (i != j)
               {
                save  = stack [i];
                stack [i] = stack [j];
                stack [j] = save;
               }
           i ++;
          }
     }
   group2 = i;


   // Check both groups of edges for consistent orientation

   flipped = not_flipped = 0;
   for  (i = 0;  i < group1;  i ++)
     if  (stack [i] . flipped)
         flipped += stack [i] . num_good_mates;
       else
         not_flipped += stack [i] . num_good_mates;
   min = Min_int (flipped, not_flipped);
   max = Max_int (flipped, not_flipped);
   if  (min + bad_links > 1 || (min > 0 && max < GOOD_LINKS_IF_BAD))
       return;
   if  (min > 0)
       {
        bad_links += min;
        for  (i = 0;  i < group1;  i ++)
          if  ((min == flipped && stack [i] . flipped)
                 || (min == not_flipped && ! stack [i] . flipped))
              stack [i] . is_bad = TRUE;
       }
   link_ct = max;

   flipped = not_flipped = 0;
   for  (i = group1;  i < group2;  i ++)
     if  (stack [i] . flipped)
         flipped += stack [i] . num_good_mates;
       else
         not_flipped += stack [i] . num_good_mates;
   min = Min_int (flipped, not_flipped);
   max = Max_int (flipped, not_flipped);
   if  (min + bad_links > 1 || (min > 0 && max < GOOD_LINKS_IF_BAD))
       return;
   if  (min > 0)
       {
        bad_links += min;
        for  (i = group1;  i < group2;  i ++)
          if  ((min == flipped && stack [i] . flipped)
                 || (min == not_flipped && ! stack [i] . flipped))
              stack [i] . is_bad = TRUE;
       }
   link_ct += max;

   contig = GetGraphNode(ScaffoldGraph->RezGraph, cid);
   cover_stat = GetCoverageStat (contig);

   bad_allowed = Max_int (0, 1 - bad_links);
   consistent = Estimate_Chunk_Ends
                    (stack, group1, & left_end1,
                     & right_end1, contig, & edge_quality1,
                     fill_chunks, & gap1, & new_scaff1, & bad_allowed);
   if  (! consistent)
       return;

   if  (stack [0] . flipped)
       {
        a_end1 = right_end1;
        b_end1 = left_end1;
       }
     else
       {
        a_end1 = left_end1;
        b_end1 = right_end1;
       }

   bad_allowed = Max_int (0, bad_allowed);
   consistent = Estimate_Chunk_Ends
                    (stack + group1, group2 - group1, & left_end2,
                     & right_end2, contig, & edge_quality2,
                     fill_chunks, & gap2, & new_scaff2, & bad_allowed);
   if  (! consistent)
       return;

   if  (stack [group1] . flipped)
       {
        a_end2 = right_end2;
        b_end2 = left_end2;
       }
     else
       {
        a_end2 = left_end2;
        b_end2 = right_end2;
       }

#if  0
fprintf (stderr, "### chunk %d joins\n", cid);
fprintf (stderr, "    scaff %d  gap %d  a_end = (%.0f, %.0f)  b_end = (%.0f, %.0f)\n",
         new_scaff1, gap1, a_end1 . mean, sqrt (a_end1 . variance),
         b_end1 . mean, sqrt (b_end1 . variance));
fprintf (stderr, "    scaff %d  gap %d  a_end = (%.0f, %.0f)  b_end = (%.0f, %.0f)\n",
         new_scaff2, gap2, a_end2 . mean, sqrt (a_end2 . variance),
         b_end2 . mean, sqrt (b_end2 . variance));
#endif

   Add_Join_Entry (cid, new_scaff1, new_scaff2, a_end1, b_end1, a_end2, b_end2,
                   left_end1, right_end1, stack [0] . flipped, gap1, edge_quality1,
                   cover_stat, link_ct, stack [0]);

   return;
  }



static int  Choose_Best_Stones
    (int start_id, int target_id, Gap_Chunk_t * node [], int num_nodes,
     int edge [], Stone_Edge_t pool [], double ref_position, double factor,
     LengthT * target_position)

//  Select best stones in gap from chunk  start_id  to chunk
//  target_id  using chunks in  node [0 .. (num_nodes - 1)]
//  and edges in  edge and  pool .  Only use edges not marked
//  as  bad .   ref_position  is the scaffold position of the
//  start of the search.   factor  is  +1.0  for forward searches
//  and  -1.0  for reverse.  Return the number of stones marked
//  kept that are left in the gap.  Set  (* target_position)
//  to the position of the end of the gap indicated by the
//  stones.

  {
   int  sorted [num_nodes];    // non-standard run-time allocation
   int  start_sub = 0, target_sub = 0, num_sorted;
   int  num_kept;
   int  i;

   for  (i = 0;  i < num_nodes;  i ++)
     {
      if  (node [i] -> chunk_id == start_id)
          start_sub = i;
      if  (node [i] -> chunk_id == target_id)
          target_sub = i;
     }
   if  (target_id == -1)
       target_sub = -1;

   Topologically_Sort (node, num_nodes, start_sub, edge, pool,
                       sorted, & num_sorted, FALSE);

#if  SHOW_STONE_CONFIRM
{
 int  i;

 fprintf (stderr, "Topological Sort:  num_nodes = %d  num_sorted = %d\n",
          num_nodes, num_sorted);
 for  (i = 0;  i < num_sorted;  i ++)
   {
    int  cid = node [sorted [i]] -> chunk_id;

    fprintf (stderr, "%5d  %s\n", cid,
             cid != start_id && cid != target_id
               && REF (cid) . is_singleton ? "*UNIQUE*" : "");
    
   }
}
#endif

   for  (i = 0;  i < num_nodes;  i ++)
     node [i] -> keep = FALSE;
   for  (i = 0;  i < num_sorted;  i ++)
     node [sorted [i]] -> keep = TRUE;

   if  (num_sorted <= 0)
       num_kept = 0;
     else
       {
        assert (node [sorted [0]] -> chunk_id == start_id);

        num_kept = Set_Longest_Path (sorted, num_sorted, node, num_nodes,
                                     target_sub, edge, pool, ref_position,
                                     factor, target_position);
       }

// Special temporary hack for Lactobacillus
#if  0
{
 for  (i = 0;  i < num_nodes;  i ++)
   if  (! node [i] -> keep
          && node [i] -> link_ct >= 10)
       {
        node [i] -> keep = TRUE;
        num_kept ++;
       }
}
#endif

   return  num_kept;
  }



static void  Choose_Safe_Chunks
    (Scaffold_Fill_t * fill_chunks, int min_good_links, int min_cover_stat)

//  Choose unresolved chunks that have clear evidence of belonging
//  to a gap and assign them to that gap in the  fill_chunks  structure.
//  Chunks must have at least  min_good_links  to a scaffold in order
//  to be selected.
//  Chunks whose coverage statistic is below  min_cover_stat  are not
//  selected.
//  Global  Ref  has information about positions of unique chunks in
//  scaffolds.
//
//
// The strategy is the following:
// 
//  1.  Identify chunks with edge mates that incorporate at least
//      two mate links to unique chunks in the same scaffold, and no
//      edge mates with links to any other scaffold.  (Later we
//      handle chunks that appear to be between two scaffolds.)  Edge
//      mates that are flagged as anomalous or that are based
//      solely on repeat overlaps are ignored.
// 
//  2.  For each selected chunk, combine the information from the
//      selected edge mates to estimate the start and end positions
//      of the chunk in the scaffold.  Calculations are done with
//      variances adjusted so that the nearest left
//      position) unique chunk in the scaffold (we'll call it the
//      "cornerstone") has left position with variance zero.
// 
//  3.  Check the calculated start and end positions for
//      consistency with each edge mate involved in their
//      calculation.  We'll say they're consistent if the
//      3-stddev interval around the calculated position intersects
//      the 3-stddev interval specified by the edge mate.  If
//      any edge mate is violated, the chunk is excluded from
//      further processing in this phase.
//
// 4.  For each consistent chunk, assign it to the nearest gap
//     in the scaffold structure.
//     We also will assign chunks to "end gaps", i.e., before the
//     first chunk in a scaffold or after the last one.

  {
   int  cid, scaff_id, prev_cid = -1;
   int  cover_stat;
   int  non_unique_ct = 0;
   int  unique_connect_ct = 0;
   int  num_consistent = 0;
   int  num_placed = 0;
   ContigT  * contig, * prev_contig = NULL;
   GraphNodeIterator  contig_iterator;
#if  MAKE_CAM_FILE
   int   cam_colour;
   char  annotation_string [MAX_STRING_LEN];
#endif

   InitGraphNodeIterator (& contig_iterator, ScaffoldGraph -> RezGraph,
                          GRAPH_NODE_DEFAULT);
   while  ((contig = NextGraphNodeIterator (& contig_iterator)) != NULL)
     {
      if  (contig == prev_contig)
          {
           fprintf (stderr, "YIKES:  got contig %p (cid %d) twice\n",
                    contig, contig -> id);
           exit (-1);
          }
      cid = contig -> id;
      if  (cid == prev_cid)
          {
           fprintf (stderr, "YIKES:  got cid %d twice from different contigs\n",
                    cid);
           exit (-1);
          }
      prev_contig = contig;
      prev_cid = cid;

#if  CHECK_CELSIM_COORDS
      if  (contig -> aEndCoord >=0 && contig -> bEndCoord >= 0)
          {
           if  (contig -> aEndCoord < contig -> bEndCoord)
               {
                Chunk_Info [cid] . celsim_left = contig -> aEndCoord;
                Chunk_Info [cid] . celsim_right = contig -> bEndCoord;
               }
             else
               {
                Chunk_Info [cid] . celsim_left = contig -> bEndCoord;
                Chunk_Info [cid] . celsim_right = contig -> aEndCoord;
               }
          }
        else
          Chunk_Info [cid] . celsim_left = Chunk_Info [cid] . celsim_right = 0;
#endif
      Chunk_Info [cid] . calc_left = -1;
      Chunk_Info [cid] . calc_right = -1;
      cover_stat = GetCoverageStat (contig);

      if  (Is_Unique (contig) || IsSurrogate(contig))
          {
#if  MAKE_CAM_FILE
           int  rel_pos, left, right;
           
           cam_colour = UNIQUE_COLOUR;
           scaff_id = REF (cid) . scaff_id;
           rel_pos = REF (cid) . rel_pos;

           sprintf (annotation_string,
"  Scaff #%d  rel pos #%d  start = <%.0f,%.0f>  end = <%.0f,%.0f>  cov = %d  typ = %s",
                    scaff_id, rel_pos,
                    contig -> offsetAEnd . mean, sqrt (contig -> offsetAEnd . variance),
                    contig -> offsetBEnd . mean, sqrt (contig -> offsetBEnd . variance),
                    cover_stat,
                    CGB_Type_As_String (contig -> flags . bits . cgbType) 
                   );
#if  SHOW_CALC_COORDS
	   assert(contig->scaffoldID > NULLINDEX);
           if  (contig -> offsetAEnd . mean < contig -> offsetBEnd . mean)
               {
                left = contig -> offsetAEnd . mean;
                right = contig -> offsetBEnd . mean;
               }
             else
               {
                left = contig -> offsetBEnd . mean;
                right = contig -> offsetAEnd . mean;
               }
                
           Chunk_Info [cid] . scaff_id = scaff_id;
           if  (Scaffold_Start [scaff_id] <= Scaffold_End [scaff_id])
               {
                Chunk_Info [cid] . calc_left = left + Scaffold_Start [scaff_id];
                Chunk_Info [cid] . calc_right = right + Scaffold_Start [scaff_id];
               }
             else
               {
                Chunk_Info [cid] . calc_left = Scaffold_Start [scaff_id] - right;
                Chunk_Info [cid] . calc_right = Scaffold_Start [scaff_id] - left;
               }
#endif
#endif
          }
#if  0
      else if  (! Maybe_Rock (contig))
          {
           non_unique_ct ++;
#if  MAKE_CAM_FILE
           cam_colour = NO_CONNECT_COLOUR;
           sprintf (annotation_string,
                    "  No connections to uniques  cov = %d  typ = %s",
                    cover_stat,
                    CGB_Type_As_String (contig -> flags . bits . cgbType));
#endif
          }
#endif
      else if  (contig -> info . CI . numInstances > 0)
          {
           fprintf (stderr, "SURPRISE:  contig %d has surrogates...skipping it\n",
                    cid);
          }
        else
          {
           int  gap;
           Stack_Entry_t  stack [STACK_SIZE];
           int  stack_top = 0;
           int  is_lacto_rock = FALSE;
           GraphEdgeIterator  ci_edges;
           ChunkInstanceT  * unitig;
           CIEdgeT  * edge;
           
           non_unique_ct ++;
#if  MAKE_CAM_FILE
           cam_colour = NO_CONNECT_COLOUR;
           sprintf (annotation_string,
                    "  No connections to uniques  cov = %d  typ = %s",
                    cover_stat,
                    CGB_Type_As_String (contig -> flags . bits . cgbType));
#endif

           // Put edges from chunk to a unique chunk onto a stack

           InitGraphEdgeIterator (ScaffoldGraph->RezGraph, cid, ALL_END,
                                  ALL_EDGES, GRAPH_EDGE_DEFAULT,
                                  & ci_edges);
           while  ((edge = NextGraphEdgeIterator (& ci_edges)) != NULL)
             {
              ChunkInstanceT  * other_chunk;

              if  (isProbablyBogusEdge (edge)
                     || isSloppyEdge (edge))
                  continue;

              if  (edge -> idA == cid)
                  other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                   edge -> idB);
                else
                  other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                   edge -> idA);
              
              if  (Is_Unique (other_chunk) && !IsSurrogate(other_chunk))
                  {
                   assert (stack_top < STACK_SIZE);
		   assert (other_chunk -> scaffoldID > NULLINDEX);
                   stack [stack_top] . chunk_id = other_chunk -> id;
                   stack [stack_top] . edge = edge;
                   stack_top ++;
                  }
             }

           // Move "good" edges to front of stack

           if  (stack_top > 0)
               {
                unique_connect_ct ++;
#if  MAKE_CAM_FILE
                cam_colour = CONNECT_COLOUR;
                sprintf (annotation_string,
                         "  %d mate edges to uniques  cov = %d  typ = %s",
                         stack_top, cover_stat,
                         CGB_Type_As_String (contig -> flags . bits . cgbType));
#endif
                stack_top = Select_Good_Edges (stack, stack_top, contig);
               }

// Special temporary hack for Lactobacillus
#if  0
{
 int  total_links = 0;
 int  hi_ct = -1, best_lo, best_hi;
 int  i, j;

 Partition_Edges  (cid, stack, stack_top, 0);

 for  (i = 0;  i < stack_top;  i = j)
   {
    int  link_ct;

    link_ct = stack [i] . num_good_mates;
    for  (j = i + 1;
            j < stack_top
              && stack [j] . partition == stack [i] . partition;
            j ++)
      link_ct += stack [j] . num_good_mates;

    if  (stack [i] . is_bad)
        continue;

    total_links += link_ct;
    if  (link_ct > hi_ct)
        {
         hi_ct = link_ct;
         best_lo = i;
         best_hi = j;
        }
   }

 if  (cover_stat < min_cover_stat
        && hi_ct >= 10
        && (double) hi_ct / total_links >= 0.80)
     {
      int  k;

      fprintf (stderr, "Lacto Rock %4d  cov = %d  hi = %d  tot = %d\n",
               cid, cover_stat, hi_ct, total_links);
      is_lacto_rock = TRUE;

      if  (best_lo == 0)
          stack_top = best_hi;
        else
          {
           for  (k = best_lo;  k < best_hi;  k ++)
             stack [k - best_lo] = stack [k];
           stack_top = k - best_lo;
          }
     }

}
#endif

           unitig = GetGraphNode (ScaffoldGraph -> CIGraph,
                                  contig -> info . Contig . AEndCI);
           if  (contig -> info . Contig . AEndCI == contig -> info . Contig . BEndCI
                  && contig -> info . Contig . AEndCI == cid
                  && unitig -> info . CI . numInstances > 0)
               {
                fprintf (stderr, "SURPRISE:  contig %d has surrogates...skipping it\n",
                         contig -> info . Contig . AEndCI);
               }
           else if  (cover_stat < min_cover_stat
                  && ! is_lacto_rock)
               {
#if  MAKE_CAM_FILE
                cam_colour = LO_COVERSTAT_COLOUR;
#endif
               }
           else if  (stack_top > 0)
               {
                int  bad_allowed, bad_links, consistent, good_total;
                LengthT  left_end, right_end;
                float  edge_quality;

                consistent = Check_Scaffold_and_Orientation
                                 (cid, stack, stack_top, & good_total,
                                  fill_chunks, & bad_links, min_good_links);

                if  (consistent)
                    {
                     num_consistent ++;
#if  MAKE_CAM_FILE
                     cam_colour = CONSISTENT_COLOUR;
                     sprintf (annotation_string,
"  %d mate links (in %d edges) to uniques  low = <%.0f, %.0f>  high = <%.0f, %.0f>  cov = %d  typ = %s",
                              good_total, stack_top,
                              left_end . mean, sqrt (left_end . variance),
                              right_end . mean, sqrt (right_end . variance),
                              cover_stat,
                              CGB_Type_As_String (contig -> flags . bits . cgbType));
#endif
                    }

                bad_allowed = Max_int (0, 1 - bad_links);
                if  (consistent && good_total >= min_good_links)
                    consistent = Estimate_Chunk_Ends
                                     (stack, stack_top, & left_end,
                                      & right_end, contig, & edge_quality,
                                      fill_chunks, & gap, & scaff_id,
                                      & bad_allowed);

                if  (consistent && good_total >= min_good_links)
                    {
#if  CHECK_CELSIM_COORDS
                     int  diff, this_offset;
                     double  cutoff;
#endif
#if  MAKE_CAM_FILE && SHOW_CALC_COORDS
                     int64  left_coord, right_coord;
#endif
                     int assign_succeeded;
                     assign_succeeded
                         = Assign_To_Gap (cid, left_end, right_end,
                                          gap, scaff_id,
                                          stack [0] . flipped, fill_chunks,
                                          edge_quality, cover_stat, good_total,
                                          ' ');

                     num_placed ++;

#if  MAKE_CAM_FILE
                     if  (assign_succeeded)
                         cam_colour = PLACED_COLOUR;
                     sprintf (annotation_string,
"  %d mate links (in %d edges) to uniques  low = <%.0f, %.0f>  high = <%.0f, %.0f>  cov = %d  typ = %s",
                              good_total, stack_top,
                              left_end . mean, sqrt (left_end . variance),
                              right_end . mean, sqrt (right_end . variance),
                              cover_stat,
                              CGB_Type_As_String (contig -> flags . bits . cgbType));
#endif
#if  CHECK_CELSIM_COORDS
                     if  (Scaffold_Flipped [REF (stack [0] . chunk_id) . scaff_id])
                         this_offset
                             = Max_int (contig -> aEndCoord, contig -> bEndCoord)
                                   + left_end . mean;
                       else
                         this_offset
                             = Min_int (contig -> aEndCoord, contig -> bEndCoord)
                                   - left_end . mean;
                     diff = abs (this_offset - stack [0] . celsim_offset);
                     cutoff = Max_int (6 * sqrt (stack [0] . edge
                                                   -> distance . variance
                                                   + stack [0] . source_variance),
                                       100);

                     if  (contig -> flags . bits . cgbType != UU_CGBTYPE)
                         fprintf (stderr,
                                  "### Placed confused chunk #%d  type = %s\n",
                                  cid,
                                  CGB_Type_As_String
                                      (contig -> flags . bits . cgbType));
                     else if  (diff > cutoff )
                         {
                          fprintf (stderr, "### Misplaced chunk #%d  diff = %d\n",
                                   cid, diff);
                          cam_colour = MISPLACED_COLOUR;
                         }
#endif
#if  MAKE_CAM_FILE && SHOW_CALC_COORDS
                     scaff_id = REF (stack [0] . chunk_id) . scaff_id;
                     if  (Scaffold_Start [scaff_id] < Scaffold_End [scaff_id])
                         {
                          left_coord = left_end . mean + Scaffold_Start [scaff_id];
                          right_coord = right_end . mean + Scaffold_Start [scaff_id];
                         }
                       else
                         {
                          left_coord = Scaffold_Start [scaff_id] - right_end . mean;
                          right_coord = Scaffold_Start [scaff_id] - left_end . mean;
                         }
                     if  (left_coord < 0)
                         left_coord = 0;
                     if  (right_coord < 0)
                         right_coord = 0;
                     Chunk_Info [cid] . scaff_id = scaff_id;
                     Chunk_Info [cid] . calc_left = left_coord;
                     Chunk_Info [cid] . calc_right = right_coord;
#endif
                    }
               }
          }

#if  MAKE_CAM_FILE
      Chunk_Info [cid] . colour = cam_colour;
      Chunk_Info [cid] . annotation = strdup (annotation_string);
      if  (contig -> aEndCoord >= 0 && contig -> bEndCoord >= 0)
          {
           if  (contig -> aEndCoord <= contig -> bEndCoord)
               {
#if  0
                fprintf (Cam_File,
                         "%dCHUNKREZ: %d A%dREZ %d # chunk %d forward %s\n",
                         cid, contig -> aEndCoord,
                         cam_colour, contig -> bEndCoord,
                         cid, annotation_string);
#endif
#if  SHOW_CALC_COORDS
                if  (cam_colour != UNIQUE_COLOUR
                       && cam_colour != PLACED_COLOUR
                       && cam_colour != MISPLACED_COLOUR)
                    {
                     Chunk_Info [cid] . scaff_id = scaff_id;
                     Chunk_Info [cid] . calc_left = contig -> aEndCoord;
                     Chunk_Info [cid] . calc_right = contig -> bEndCoord;
                    }
#endif
                Chunk_Info [cid] . colour = cam_colour;
               }
             else
               {
#if  0
                fprintf (Cam_File,
                         "%dCHUNKREZ: %d A%dREZ %d # chunk %d reverse %s\n",
                         cid, contig -> bEndCoord, cam_colour,
                         contig -> aEndCoord,
                         cid, annotation_string);
#endif
#if  SHOW_CALC_COORDS
                if  (cam_colour != UNIQUE_COLOUR
                       && cam_colour != PLACED_COLOUR
                       && cam_colour != MISPLACED_COLOUR)
                    {
                     Chunk_Info [cid] . scaff_id = scaff_id;
                     Chunk_Info [cid] . calc_left = contig -> bEndCoord;
                     Chunk_Info [cid] . calc_right = contig -> aEndCoord;
                    }
#endif
                Chunk_Info [cid] . colour = cam_colour;
               }
          }
#endif
     }

   fprintf (stderr, "             Non-unique chunks: %7d\n", non_unique_ct);
   fprintf (stderr, "         With edges to uniques: %7d\n", unique_connect_ct);
   fprintf (stderr, "                    Consistent: %7d\n", num_consistent);
   fprintf (stderr, "       Placed in scaffold gaps: %7d\n", num_placed);
   fprintf (stderr, "                  Total chunks: %7d\n", Num_Chunks);
   fprintf (stderr, "   Discriminator unique chunks: %7d\n",
            ScaffoldGraph -> numDiscriminatorUniqueCIs);
   fprintf (stderr, "             Scaffolded chunks: %7d\n",
            Num_Chunks - non_unique_ct);

   return;
  }



static void  Choose_Stones
    (Scaffold_Fill_t * fill_chunks, int min_good_links, int min_cover_stat,
     int allow_bogus_edges)

//  Choose unresolved chunks that have at least  min_good_links
//  to a scaffold and assign them to gaps in the  fill_chunks  structure.
//  Chunks may be assigned to multiple gaps if the links are
//  inconsistent.
//  Chunks whose coverage statistic is below  min_cover_stat  are not
//  selected.
//  If  allow_bogus_edges  is true, can use mate links that are marked
//   isProbablyBogusEdge () ; otherwise, ignore those edges.
//  Global  Ref  has information about positions of unique chunks in
//  scaffolds.
//
//
// The strategy is the following:
// 
//  1.  Find chunks that have at least  min_good_links  to the
//      same scaffold.  Edge
//      mates that are flagged as anomalous or that are based
//      solely on repeat overlaps are ignored.
// 
//  2.  For each selected chunk, combine the information from the
//      selected edge mates to estimate the start and end positions
//      of the chunk in the scaffold.  Calculations are done with
//      variances adjusted so that the nearest left
//      position) unique chunk in the scaffold (we'll call it the
//      "cornerstone") has left position with variance zero.
// 
//  3.  Check the calculated start and end positions for
//      consistency with each edge mate involved in their
//      calculation.  We'll say they're consistent if the
//      3-stddev interval around the calculated position intersects
//      the 3-stddev interval specified by the edge mate.  If
//      any edge mate is violated, it is put in a separate group
//      to imply a different position relative to the scaffold.
//
// 4.  For each chunk position assign it to all the gaps in the
//     scaffold that it intersects.
//     We also will assign chunks to "end gaps", i.e., before the
//     first chunk in a scaffold or after the last one.
//     [Should we prevent the same chunk being confirmed off 2
//     different scaffold ends?  It might really be the same position.]

  {
   int  cid, scaff_id;
   int  cover_stat;
   int  non_unique_ct = 0;
   int  unique_connect_ct = 0;
   int  num_stones = 0;
   int  cam_colour;
   ContigT  * chunk;
   GraphNodeIterator  contig_iterator;
#if  MAKE_CAM_FILE
   char  annotation_string [MAX_STRING_LEN];
#endif

   InitGraphNodeIterator (& contig_iterator, ScaffoldGraph -> RezGraph,
                          GRAPH_NODE_DEFAULT);
   while  ((chunk = NextGraphNodeIterator (& contig_iterator)) != NULL)
     {
      cid = chunk->id;

#if  CHECK_CELSIM_COORDS
      if  (chunk -> aEndCoord >=0 && chunk -> bEndCoord >= 0)
          {
           if  (chunk -> aEndCoord < chunk -> bEndCoord)
               {
                Chunk_Info [cid] . celsim_left = chunk -> aEndCoord;
                Chunk_Info [cid] . celsim_right = chunk -> bEndCoord;
               }
             else
               {
                Chunk_Info [cid] . celsim_left = chunk -> bEndCoord;
                Chunk_Info [cid] . celsim_right = chunk -> aEndCoord;
               }
          }
        else
          Chunk_Info [cid] . celsim_left = Chunk_Info [cid] . celsim_right = 0;
#endif
      Chunk_Info [cid] . calc_left = -1;
      Chunk_Info [cid] . calc_right = -1;
      cover_stat = GetCoverageStat (chunk);

      if  (Is_Unique (chunk)
             && ! (REF (cid) . is_singleton && UNIQUES_CAN_BE_STONES))
          {
#if  MAKE_CAM_FILE
           int  rel_pos, left, right;
           
           cam_colour = UNIQUE_COLOUR;
           scaff_id = REF (cid) . scaff_id;
           rel_pos = REF (cid) . rel_pos;

           sprintf (annotation_string,
"  Scaff #%d  rel pos #%d  start = <%.0f,%.0f>  end = <%.0f,%.0f>  cov = %d  typ = %s",
                    scaff_id, rel_pos,
                    chunk -> offsetAEnd . mean, sqrt (chunk -> offsetAEnd . variance),
                    chunk -> offsetBEnd . mean, sqrt (chunk -> offsetBEnd . variance),
                    cover_stat,
                    CGB_Type_As_String (chunk -> flags . bits . cgbType) 
                   );
#if  SHOW_CALC_COORDS
	   assert(chunk->scaffoldID > NULLINDEX);
           if  (chunk -> offsetAEnd . mean < chunk -> offsetBEnd . mean)
               {
                left = chunk -> offsetAEnd . mean;
                right = chunk -> offsetBEnd . mean;
               }
             else
               {
                left = chunk -> offsetBEnd . mean;
                right = chunk -> offsetAEnd . mean;
               }
                
           Chunk_Info [cid] . scaff_id = scaff_id;
           if  (Scaffold_Start [scaff_id] <= Scaffold_End [scaff_id])
               {
                Chunk_Info [cid] . calc_left = left + Scaffold_Start [scaff_id];
                Chunk_Info [cid] . calc_right = right + Scaffold_Start [scaff_id];
               }
             else
               {
                Chunk_Info [cid] . calc_left = Scaffold_Start [scaff_id] - right;
                Chunk_Info [cid] . calc_right = Scaffold_Start [scaff_id] - left;
               }
#endif
#endif
          }
        else
          {
           int  gap;
           Stack_Entry_t  stack [STACK_SIZE];
           int  stack_top = 0;
           GraphEdgeIterator  ci_edges;
           ChunkInstanceT  * unitig;
           int  problem;
           CIEdgeT  * edge;
           
           non_unique_ct ++;
#if  MAKE_CAM_FILE
           cam_colour = NO_CONNECT_COLOUR;
           sprintf (annotation_string,
                    "  No connections to uniques  cov = %d  typ = %s",
                    cover_stat,
                    CGB_Type_As_String (chunk -> flags . bits . cgbType));
#endif

           if  (Single_Fragment_Only)
               {
                ChunkInstanceT  * first_chunk;
                ContigT  * contig;

                contig = GetGraphNode
                           (ScaffoldGraph -> RezGraph, cid);
                assert (contig != NULL);
                first_chunk
                    = GetGraphNode
                          (ScaffoldGraph -> CIGraph, contig -> info . Contig . AEndCI);
                problem = (contig -> info . Contig . AEndCI
                             != contig -> info . Contig . BEndCI);
                if  (problem)
                    fprintf (stderr,
                             "YOWZA!! Contig %d not unique but has two unitigs\n",
                             cid);
                if  (problem || first_chunk -> info . CI . numFragments != 1)
                    continue;
               }

           if  (chunk -> info . Contig . AEndCI != chunk -> info . Contig . BEndCI)
               {
                fprintf (stderr, "SURPRISE:  contig %d has > 1 unitigs...skipping it\n",
                         cid);
                continue;
               }
           if  (chunk -> info . Contig . AEndCI != cid)
               {
                fprintf (stderr, "SURPRISE:  contig %d has unitig id %d...skipping it\n",
                         cid, chunk -> info . Contig . AEndCI);
                continue;
               }
           unitig = GetGraphNode (ScaffoldGraph -> CIGraph,
                                  chunk -> info . Contig . AEndCI);
           if  (unitig -> type != UNRESOLVEDCHUNK_CGW)
               {
                fprintf (stderr, "SURPRISE:  unitig %d type = %d...skipping it\n",
                         chunk -> info . Contig . AEndCI, (int) (unitig -> type));
                continue;
               }
           if  (Contained_Only_Switch && unitig -> info . CI . numInstances > 0)
               {
                fprintf (stderr, "SURPRISE:  contig %d has surrogates...skipping it\n",
                         chunk -> info . Contig . AEndCI);
                continue;
               }

           // Put edges from chunk to a unique chunk onto a stack

           InitGraphEdgeIterator (ScaffoldGraph->RezGraph, cid, ALL_END,
                                  ALL_EDGES, GRAPH_EDGE_DEFAULT,
                                  & ci_edges);
           while  ((edge = NextGraphEdgeIterator (& ci_edges)) != NULL)
             {
              ChunkInstanceT  * other_chunk;

              if  ((isProbablyBogusEdge (edge) && ! allow_bogus_edges)
                     || isSloppyEdge (edge))
                  continue;

              if  (edge -> idA == cid)
                  other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                   edge -> idB);
                else
                  other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                   edge -> idA);

              if  (Is_Unique (other_chunk))
                  {
                   assert (stack_top < STACK_SIZE);
                   assert(other_chunk->scaffoldID > NULLINDEX);
                   stack [stack_top] . chunk_id = other_chunk -> id;
                   stack [stack_top] . edge = edge;
                   stack_top ++;
                  }
             }
#if  VERBOSE
           if  (! Maybe_Stone (unitig) && stack_top > 0)
               fprintf (stderr, "cid = %d  Not potential stone but stack = %d\n",
                        cid, stack_top);
#endif

           // Move "good" edges to front of stack

           if  (stack_top > 0)
               {
                unique_connect_ct ++;
#if  MAKE_CAM_FILE
                cam_colour = CONNECT_COLOUR;
                sprintf (annotation_string,
                         "  %d mate edges to uniques  cov = %d  typ = %s",
                         stack_top, 
			 cover_stat,
                         CGB_Type_As_String (chunk -> flags . bits . cgbType));
#endif
                stack_top = Select_Good_Edges (stack, stack_top, chunk);
               }

           if  (cover_stat < min_cover_stat)
               {
#if  MAKE_CAM_FILE
                cam_colour = LO_COVERSTAT_COLOUR;
#endif
               }
           else if  (stack_top > 0)
               {
                LengthT  left_end, right_end;
                float  edge_quality;
                char  copy_letter = 'a';
                int  bad_allowed, consistent;
                int  i, j;

                Partition_Edges  (cid, stack, stack_top, min_good_links);

                for  (i = 0;  i < stack_top;  i = j)
                  {
                   int  link_ct, assign_succeeded;
#if  MAKE_CAM_FILE && SHOW_CALC_COORDS
                   int64  left_coord, right_coord;
#endif

                   link_ct = stack [i] . num_good_mates;
                   for  (j = i + 1;
                           j < stack_top
                             && stack [j] . partition == stack [i] . partition;
                           j ++)
                     link_ct += stack [j] . num_good_mates;

                   if  (stack [i] . is_bad)
                       continue;

                   bad_allowed = 0;
                   consistent = Estimate_Chunk_Ends
                                    (stack + i, j - i, & left_end,
                                     & right_end, chunk, & edge_quality,
                                     fill_chunks, & gap, & scaff_id,
                                     & bad_allowed);

                   if  (! consistent)   // Come back here and add partition refinement
                       continue;        //   for case where there are alternative
                                        //   positions in the same scaffold with the
                                        //   same orientation.

                   num_stones ++;

#if  MAKE_CAM_FILE
                   cam_colour = STONE_COLOUR;
                   sprintf (annotation_string,
"  stone copy %c  low = <%.0f, %.0f>  high = <%.0f, %.0f>  cov = %d  typ = %s",
                            copy_letter,
                            left_end . mean, sqrt (left_end . variance),
                            right_end . mean, sqrt (right_end . variance),
                            cover_stat,
                            CGB_Type_As_String (chunk -> flags . bits . cgbType));
#endif

                   assign_succeeded
                       = Assign_To_Gap (cid, left_end, right_end,
                                        gap, scaff_id,
                                        stack [i] . flipped, fill_chunks,
                                        edge_quality, cover_stat, link_ct,
                                        copy_letter);
                   if  (! assign_succeeded)
                       cam_colour = REJECT_COLOUR;

#if  MAKE_CAM_FILE && SHOW_CALC_COORDS
                   scaff_id = REF (stack [0] . chunk_id) . scaff_id;
                   if  (Scaffold_Start [scaff_id] < Scaffold_End [scaff_id])
                       {
                        left_coord = left_end . mean + Scaffold_Start [scaff_id];
                        right_coord = right_end . mean + Scaffold_Start [scaff_id];
                       }
                     else
                       {
                        left_coord = Scaffold_Start [scaff_id] - right_end . mean;
                        right_coord = Scaffold_Start [scaff_id] - left_end . mean;
                       }
                   if  (left_coord < 0)
                       left_coord = 0;
                   if  (right_coord < 0)
                       right_coord = 0;
                   Chunk_Info [cid] . scaff_id = scaff_id;
                   Chunk_Info [cid] . calc_left = left_coord;
                   Chunk_Info [cid] . calc_right = right_coord;
#endif
                   copy_letter ++;
                  }
               }
          }

#if  MAKE_CAM_FILE
      Chunk_Info [cid] . colour = cam_colour;
      Chunk_Info [cid] . annotation = strdup (annotation_string);
      if  (chunk -> aEndCoord >= 0 && chunk -> bEndCoord >= 0)
          {
           if  (chunk -> aEndCoord <= chunk -> bEndCoord)
               {
#if  SHOW_CALC_COORDS
                if  (cam_colour != UNIQUE_COLOUR
                       && cam_colour != PLACED_COLOUR
                       && cam_colour != MISPLACED_COLOUR)
                    {
                     Chunk_Info [cid] . scaff_id = scaff_id;
                     Chunk_Info [cid] . calc_left = chunk -> aEndCoord;
                     Chunk_Info [cid] . calc_right = chunk -> bEndCoord;
                    }
#endif
                Chunk_Info [cid] . colour = cam_colour;
               }
             else
               {
#if  SHOW_CALC_COORDS
                if  (cam_colour != UNIQUE_COLOUR
                       && cam_colour != PLACED_COLOUR
                       && cam_colour != MISPLACED_COLOUR)
                    {
                     Chunk_Info [cid] . scaff_id = scaff_id;
                     Chunk_Info [cid] . calc_left = chunk -> bEndCoord;
                     Chunk_Info [cid] . calc_right = chunk -> aEndCoord;
                    }
#endif
                Chunk_Info [cid] . colour = cam_colour;
               }
          }    
#endif
     }

   fprintf (stderr, "             Non-unique chunks: %7d\n", non_unique_ct);
   fprintf (stderr, "         With edges to uniques: %7d\n", unique_connect_ct);
   fprintf (stderr, "                        Stones: %7d\n", num_stones);
   fprintf (stderr, "                  Total chunks: %7d\n", Num_Chunks);
   fprintf (stderr, "   Discriminator unique chunks: %7d\n",
            ScaffoldGraph -> numDiscriminatorUniqueCIs);
   fprintf (stderr, "             Scaffolded chunks: %7d\n",
            Num_Chunks - non_unique_ct);

   return;
  }



float  CIEdge_Quality
    (CIEdgeT  * edge)

//  Return a value that indicates the reliability of the edge mate
//  (* edge) .
  {
   float  val;

   if  (isProbablyBogusEdge (edge))
       return  0.001;

   val = 3.0 * edge -> edgesContributing;
   if  (edge -> flags . bits . isPossibleChimera)
       val -= 0.2;
   if  (edge -> flags . bits . hasContributingOverlap)
       val -= 1.0;
   if  (edge -> flags . bits . hasRepeatOverlap)
       val -= 2.0;
   if  (edge -> flags . bits . hasTandemOverlap)
       val -= 3.0;
   if  (edge -> flags . bits . hasGuide)
       val -= 0.5;

   return  val;
  }



static int  Chunk_Contained_In_Chunk
    (Gap_Chunk_t * A, Gap_Chunk_t * B)

//  Return  TRUE  iff chunk  A  has an overlap in the hash table
//  that indicates it is contained within chunk  B .

  {
   ChunkOrientationType  orient;
   ChunkOverlapCheckT  olap;
   int  olap_found;

   // Try with  A  on the left of  B
   if  (A -> start. mean <= A -> end . mean)
       {
        if  (B -> start . mean <= B -> end . mean)
            orient = AB_AB;
          else
            orient = AB_BA;
       }
     else
       {
        if  (B -> start . mean <= B -> end . mean)
            orient = BA_AB;
          else
            orient = BA_BA;
       }
   olap_found = LookupOverlap (ScaffoldGraph -> RezGraph,
                               A -> chunk_id,
                               B -> chunk_id,
                               orient, & olap);
   if  (olap_found && olap . BContainsA)
       return  TRUE;

   switch  (orient)
     {
      case  AB_AB :
      case  BA_BA :
        break;              // No change
      case  AB_BA :
        orient = BA_AB;
        break;
      case  BA_AB :
        orient = AB_BA;
        break;
      default :
        fprintf (stderr, "YIKES:  Bad orientation = %d\n", (int) orient);
        assert (FALSE);
     }
   olap_found = LookupOverlap (ScaffoldGraph -> RezGraph,
                               A -> chunk_id,
                               B -> chunk_id,
                               orient, & olap);
   if  (olap_found && olap . BContainsA)
       return  TRUE;
   
   return  FALSE;
  }



static int  Chunk_Contained_In_Scaff
    (Gap_Chunk_t * A, int cid)

//  Return  TRUE  iff chunk  A  has an overlap in the hash table
//  that indicates it is contained within scaffold chunk with ID  cid .

  {
   ChunkInstanceT  * B;
   ChunkOrientationType  orient;
   ChunkOverlapCheckT  olap;
   int  olap_found;

   B = GetGraphNode (ScaffoldGraph -> RezGraph, cid);

   // Try with  A  on the left of  B
   if  (A -> start. mean <= A -> end . mean)
       {
        if  (B -> offsetAEnd . mean <= B -> offsetBEnd . mean)
            orient = AB_AB;
          else
            orient = AB_BA;
       }
     else
       {
        if  (B -> offsetAEnd . mean <= B -> offsetBEnd . mean)
            orient = BA_AB;
          else
            orient = BA_BA;
       }
   olap_found = LookupOverlap (ScaffoldGraph -> RezGraph,
                               A -> chunk_id,
                               cid,
                               orient, & olap);
   if  (olap_found && olap . BContainsA)
       return  TRUE;

   switch  (orient)
     {
      case  AB_AB :
      case  BA_BA :
        break;              // No change
      case  AB_BA :
        orient = BA_AB;
        break;
      case  BA_AB :
        orient = AB_BA;
        break;
      default :
        fprintf (stderr, "YIKES:  Bad orientation = %d\n", (int) orient);
        assert (FALSE);
     }
   olap_found = LookupOverlap (ScaffoldGraph -> RezGraph,
                               A -> chunk_id,
                               cid,
                               orient, & olap);
   if  (olap_found && olap . BContainsA)
       return  TRUE;
   
   return  FALSE;
  }


static void  Clear_Keep_Flags
    (Scaffold_Fill_t * fill_chunks, int except_num)

//  Set all the keep flags of the entries in  fill_chunks  to  FALSE ,
//  except leave the first  except_num  of them that are  TRUE
//  unchanged.

  {
   int  scaff_id;
   int  total_chunks = 0;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
         int  k;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep
                   && (total_chunks ++) >= except_num)
                this_chunk -> keep = FALSE;
           }
        }
     }

   return;
  }



static void  Confirm_Contained
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int use_all)

//  Check whether each entry in  fill_chunks  overlaps one
//  of the scaffold contigs bounding the gap to be considered
//  contained in that contig.  If so mark its  keep  flag
//  and set its position; otherwise, set its  keep  flag
//  false.  Send output log to  fp .  If  fp  is  NULL , do
//  no actual writing.
//  If  use_all  is true, use all the chunks in  fill_chunks ;
//  otherwise, use only the ones whose  keep  flag is already
//  true.

  {
   int  scaff_id;

   if  (fp != NULL)
       fprintf (fp, "\n Confirm_Contained:\n");

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = GAPS_TO_ADJUST)
        {
         int  k, ct;
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
         Gap_Chunk_t  * left_scaff_contig = NULL;
         Gap_Chunk_t  * right_scaff_contig = NULL;
         char  * left_scaff_sequence = NULL, * right_scaff_sequence = NULL;
         int  avail_extend, used_left_extend, used_right_extend;

         if  (this_gap -> num_chunks == 0)
             continue;

#if  VERBOSE
if  (fp != NULL)
    fprintf (fp, "\nConfirm_Contained Scaff %d  Gap %d\n", scaff_id, j);
#endif

         if  (this_gap -> left_cid >= 0)
             {
              for  (k = 0;  k < this_gap -> num_chunks;  k ++)
                if  (this_gap -> chunk [k] . chunk_id == this_gap -> left_cid)
                    {
                     left_scaff_contig = this_gap -> chunk + k;
                     break;
                    }
              if  (k >= this_gap -> num_chunks)
                  {
                   fprintf (stderr, "didn't find left cid = %d  scaff = %d  gap = %d\n",
                            this_gap -> left_cid, scaff_id, j);
                  }
             }
         if  (this_gap -> right_cid >= 0)
             {
              for  (k = 0;  k < this_gap -> num_chunks;  k ++)
                if  (this_gap -> chunk [k] . chunk_id == this_gap -> right_cid)
                    {
                     right_scaff_contig = this_gap -> chunk + k;
                     break;
                    }
              if  (k >= this_gap -> num_chunks)
                  {
                   fprintf (stderr, "didn't find right cid = %d\n",
                            this_gap -> right_cid);
                  }
             }

         // Don't let rock/stones stick out past the end of the the containing contig
         // enough to cause the gap to become less than - CGW_DP_MINLEN
         //  avail_extend  will have the number of bases
         // still available for sticking out into the gap on either side.
         //  used_(left|right)_extend  will have
         // the number of bases that have been used by earlier stones
         // sticking out on the respective side.
         avail_extend = (int) (this_gap -> len + CGW_DP_MINLEN);
         used_left_extend = used_right_extend = 0;
         ct = 0;
         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Overlap  * left_olap = NULL, * right_olap = NULL;
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;
            char  * this_sequence = NULL;
            int  trim;

            if  ((use_all || this_chunk -> keep)
                   && this_chunk != left_scaff_contig
                   && this_chunk != right_scaff_contig)
                {
                 this_sequence = NULL;
                 if  (left_scaff_contig != NULL)
                     left_olap
                         = Get_Chunk_Overlap
                               (left_scaff_contig, this_chunk,
                                & left_scaff_sequence, & this_sequence, fp);
                 if  (right_scaff_contig != NULL)
                     right_olap
                         = Get_Chunk_Overlap
                               (this_chunk, right_scaff_contig,
                                & this_sequence, & right_scaff_sequence, fp);

                 // allow up to 3 bases of overhang and still regard
                 // as contained
                 if  (right_olap == NULL
                        && left_olap != NULL
                        && left_olap -> begpos >= 0 && left_olap -> endpos <= 3)
                     {
                      if  (fp != NULL)
                          fprintf (fp, "cid %5d Left olap  %6d %6d %6d\n",
                               this_chunk -> chunk_id,
                               left_olap -> begpos, left_olap -> endpos,
                               left_olap -> length);
                      trim = left_olap -> endpos - used_left_extend;
                      if  (trim > 0)
                          {
                           if  (trim > avail_extend)
                               {
                                trim -= avail_extend;
                                used_left_extend += avail_extend;
                                avail_extend = 0;
                                left_olap -> endpos -= trim;
                               }
                             else
                               {
                                used_left_extend += trim;
                                avail_extend -= trim;
                               }
                          }
                      Set_Position_From_Left_Olap
                          (this_gap -> left_cid, this_chunk, left_olap);
                      this_chunk -> keep = TRUE;
                      ct ++;
                     }
                 else if  (left_olap == NULL
                        && right_olap != NULL
                        && right_olap -> begpos <= 3 && right_olap -> endpos >= 0)
                     {
                      if  (fp != NULL)
                          fprintf (fp, "cid %5d Right olap  %6d %6d %6d\n",
                               this_chunk -> chunk_id,
                               right_olap -> begpos, right_olap -> endpos,
                               right_olap -> length);
                      trim = right_olap -> begpos - used_right_extend;
                      if  (trim > 0)
                          {
                           if  (trim > avail_extend)
                               {
                                trim -= avail_extend;
                                used_right_extend += avail_extend;
                                avail_extend = 0;
                                right_olap -> begpos -= trim;
                               }
                             else
                               {
                                used_right_extend += trim;
                                avail_extend -= trim;
                               }
                          }
                      Set_Position_From_Right_Olap
                          (this_chunk, this_gap -> right_cid, right_olap);
                      this_chunk -> keep = TRUE;
                      ct ++;
                     }
                   else
                     this_chunk -> keep = FALSE;

                 if  (this_sequence != NULL)
                     free (this_sequence);
                }
              else
                this_chunk -> keep = FALSE;
           }

         if  (left_scaff_sequence != NULL)
             free (left_scaff_sequence);
         if  (right_scaff_sequence != NULL)
             free (right_scaff_sequence);
        }
     }

   return;
  }




static void  Confirm_Stones
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int use_all)

//  Try to confirm stones in  fill_chunks  by finding overlap paths
//  through them.  Send output log to  fp .
//  If  use_all  is true, use all the chunks in  fill_chunks ;
//  otherwise, use only the ones whose  keep  flag is already
//  true.

  {
   Target_Info_t  * target = NULL;
   int  * fill_sub = NULL;
   int  target_size = 0;
   int  scaff_id;

#if  SHOW_STONE_CONFIRM
   fprintf (fp, "\n*** Confirm Stones ***\n");
#endif

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

#if  SHOW_STONE_CONFIRM
      fprintf (fp, "\n Scaffold #%d  num_gaps = %d:\n",
               scaff_id, fill_chunks [scaff_id] . num_gaps);
#endif
      for  (j = GAPS_TO_ADJUST)
        {
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
         ChunkInstanceT  * from, * to;
         double  bound;
         LengthT  start_position;
         int  first, max_hits, max_first;
         int  found, from_end, is_forward;
         int  num_targets;
         int  k;

#if  SHOW_STONE_CONFIRM
         fprintf (fp, " Gap %3d:  (%8.0f,%7.0f)  (%8.0f,%7.0f)  <%6d,%6d> %3d\n",
                  j,
                  this_gap -> start . mean,
                  this_gap -> start . variance,
                  this_gap -> end . mean,
                  this_gap -> end . variance,
                  this_gap -> left_cid,
                  this_gap -> right_cid,
                  this_gap -> num_chunks);
         fprintf (fp, "   ref_variance = %.0f\n",
                  this_gap -> ref_variance);
#endif

         if  (this_gap -> num_chunks > target_size)
             {
              target_size = this_gap -> num_chunks;
PRALLOC (target_size * sizeof (Target_Info_t));
              target = (Target_Info_t *) safe_realloc
                           (target, target_size * sizeof (Target_Info_t));
PRALLOC (target_size * sizeof (int));
              fill_sub = (int *) safe_realloc
                             (fill_sub, target_size * sizeof (int));
             }

         if  (j == 0)
             from = GetGraphNode (ScaffoldGraph -> RezGraph,
                                  this_gap -> right_cid);
           else
             from = GetGraphNode (ScaffoldGraph -> RezGraph,
                                  this_gap -> left_cid);
         if  (j == 0 || j == fill_chunks [scaff_id] . num_gaps - 1)
             to = NULL;
           else
             to = GetGraphNode (ScaffoldGraph -> RezGraph,
                                this_gap -> right_cid);

         is_forward = (from -> offsetAEnd . mean <= from -> offsetBEnd . mean);
         if  (j == 0)
             is_forward = ! is_forward;
         if  (is_forward)
             {
              from_end = B_END;
              start_position = from -> offsetBEnd;
             }
           else
             {
              from_end = A_END;
              start_position = from -> offsetAEnd;
             }

         if  (to == NULL)
             bound = MAX_MATE_DISTANCE;
           else
             bound = this_gap -> end . mean
                        - this_gap -> start . mean
                        + 3.0 * sqrt (this_gap -> end . variance
                                       - this_gap -> start . variance);

         this_gap -> adjustment . mean = this_gap -> adjustment . variance
             = 0.0;

         num_targets = 0;
         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            double  delta, place;
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (use_all || this_chunk -> keep)
                {
                 delta = 3.0 * sqrt (this_chunk -> end . variance);
                 if  (delta < MIN_TARGET_SLOP)
                     delta = MIN_TARGET_SLOP;
                 target [num_targets] . id = this_chunk ->  chunk_id;

                 is_forward = (this_chunk -> start . mean <= this_chunk -> end . mean);
                 if  (j == 0)
                     is_forward = ! is_forward;
                 if  (is_forward)
                     {
                      place = this_chunk -> end . mean;
                      target [num_targets] . orient = AB_AB;
                                                 // only the 2nd part matters
                     }
                   else
                     {
                      place = this_chunk -> start . mean;
                      target [num_targets] . orient = AB_BA;
                                                 // only the 2nd part matters
                     }

                 if  (j == 0)
                     {
                      target [num_targets] . lo
                          = start_position . mean - place - delta;
                      target [num_targets] . hi
                          = start_position . mean - place + delta;
                     }
                   else
                     {
                      target [num_targets] . lo
                          = place - delta - start_position . mean;
                      target [num_targets] . hi
                          = place + delta - start_position . mean;
                     }
                 fill_sub [num_targets] = k;
                 num_targets ++;

                 this_chunk -> keep = FALSE;
                }
           }

         if  (num_targets > 0)
             {
              LengthT  to_position;

              found = Find_Olap_Path
                        (from, from_end, to, num_targets, target, bound,
                         & first, & max_hits, & max_first, & to_position,
//                         SKIP_TANDEM_OLAPS);
                         SKIP_TANDEM_OLAPS | SKIP_CONTAINMENT_OLAPS);
              if  (found && max_hits > 0)
                  {
                   LengthT  far_position, near_position;
                   double  min_var_adjustment, high_target_variance;
                   DirectionType  direction;
                   int  sub = max_first;

                   this_gap -> has_path = TRUE;

                   for  (k = 0;  k < max_hits;  k ++)
                     {
                      this_gap -> chunk [fill_sub [sub]] . keep = TRUE;
                      this_gap -> chunk [fill_sub [sub]] . path_confirmed = TRUE;
                      sub = target [sub] . next;
                     }

                   if  (j == 0)
                       direction = AS_REVERSE;
                     else
                       direction = AS_FORWARD;
		   Adjust_Positions (this_gap, num_targets, target, fill_sub,
                                     max_hits, max_first, start_position . mean,
                                     direction, & high_target_variance, fp);
		   
                   if  (j == 0)
                       Reverse_Positions (this_gap);

                   if  (to != NULL)
                       {
                        if  (to -> offsetAEnd . mean <= to -> offsetBEnd . mean)
                            {
                             far_position = to -> offsetBEnd;
                             near_position = to -> offsetAEnd;
                            }
                          else
                            {
                             far_position = to -> offsetAEnd;
                             near_position = to -> offsetBEnd;
                            }
                        this_gap -> adjustment . mean
                            = to_position . mean
                                - (far_position . mean - start_position . mean);
                        this_gap -> adjustment . variance
                            = to_position . variance
                                - (far_position . variance - start_position . variance);
#ifdef DEBUG_DETAILED
fprintf (stderr, "### Scaff %d  gap %d  to_pv = %.1f  far_pv = %.1f  st_pv = %.1f\n",
         scaff_id, j,
         to_position . variance, far_position . variance, start_position . variance);
#endif

                        //  Ensure that, after the insertion, the variances
                        //  of chunks in the scaffolds will be strictly
                        //  ascending.  EPSILON should be big enough to
                        //  cover any floating-point roundoff errors.

                        min_var_adjustment = start_position . variance
                                               - near_position . variance
                                               + EPSILON;
                        if  (high_target_variance > 0.0)
                            min_var_adjustment += high_target_variance;
                        if  (this_gap -> adjustment . variance
                              < min_var_adjustment)
                            this_gap -> adjustment . variance = min_var_adjustment;

#ifdef DEBUG_DETAILED
fprintf (stderr, "     adjv = %.1f  minvadj = %.1f\n",
         this_gap -> adjustment . variance, min_var_adjustment);
#endif

                       }                     
                  }
#if  SHOW_STONE_CONFIRM
              fprintf (fp, "found = %c  first = %3d  max_hits = %3d  max_first = %3d\n",
                       found ? 'T' : 'F', first, max_hits, max_first);
              if  (found && max_hits > 0)
                  {
                   int  sub = max_first;

                   for  (k = 0;  k < max_hits;  k ++)
                     {
                      if  (k > 0)
                          fprintf (fp, " --> ");
                      fprintf (fp, "%d", target [sub] . id);
                      sub = target [sub] . next;
                     }
                   fprintf (fp, "\n");
                  }

              fprintf (fp, "Targets:\n");
              for  (k = 0;  k < num_targets;  k ++)
                fprintf (fp, "%3d: %6d %8.0f %8.0f %8.0f %8.0f %3d\n",
                         k, target [k] . id, target [k] . lo, target [k] . hi,
                         target [k] . where,
                         target [k] . total,
                         target [k] . next);
              if  (found && max_hits > 0)
                  fprintf (fp, "  New to_position = [%.0f, %.0f]\n",
                           to_position . mean, to_position . variance);
#endif
             }
        }
     }

   free (target);
   free (fill_sub);

   return;
  }



static int  Contained_In
    (Placement_t * p, Placement_t * q)

//  Return  TRUE  iff the positions of  p  and  q  indicates that
//  p  is contained within  q .

  {
   double  a, b, c, d;

   a = Min_double (p -> A_end . mean, p -> B_end . mean);
   b = Max_double (p -> A_end . mean, p -> B_end . mean);
   c = Min_double (q -> A_end . mean, q -> B_end . mean);
   d = Max_double (q -> A_end . mean, q -> B_end . mean);

   return  (c <= a && a <= d && c <= b && b <= d);
  }



static int  Depth_First_Visit
    (ChunkInstanceT * from, int from_end, ChunkInstanceT * to,
     int num_targets, Target_Info_t target [], int bound, double so_far,
     int * first, int * max_hits, int * max_first, int * max_first_dist,
     double * max_first_total, int * tree_sub, int * tree_size,
     DFS_Info_t * * dfs_tree, Node_Index_t * node, int level,
     unsigned int edge_mask)

//  Continue a depth-first search from node  from , using overlaps
//  off its  from_end  (A_END or B_END) looking for an overlap
//  path to node  to  hitting as many of the nodes in
//  target [0 .. (num_targets - 1)]  as possible.
//  Do not pursue any path longer than  bound  DNA bases.
//  so_far  is the number of DNA bases in the path to  from .
//  Set  (* first)  to the subscript of the first target hit.
//  Set  (* max_hits)  to the number of targets on the path with
//  the most targets (break ties arbitrarily).
//  Set  (* max_first)  to the subscript of the first target hit on the
//  max_hits path.
//  Set  (* max_first_dist)  to the number of bases to  (* max_first) .
//  Set  (* max_first_total)  to the sum of fragment lengths on
//    the path to  (* max_first) .
//  Save  from  in  (* dfs_tree) [* tree_sub]  and increment
//  (* tree_sub) .  Increase  (* tree_size)  if necessary and
//  realloc  dfs_tree  to be larger.  Use  node  to track the
//  entries in  dfs_tree .
//  level  is used for debugging to keep track of the depth of the recursion.
//  edge_mask  indicates which edges can be used in the path.
//  Return  TRUE  if a successful path was found.  Successful
//  means hit any target nodes and reached  to .  If  to  is
//  NULL  the path needn't end anywhere special.  If  (num_targets == 0)
//  a successful path just reaches  to .

  {
   GraphEdgeIterator  edge_iterator;
   CIEdgeT  * edge;
   double  curr_max_first_total;
   int  curr_first, curr_max_hits, curr_max_first, curr_max_first_dist;
   int  i, sub, succeeded = FALSE;

   sub = (* tree_sub) ++;
   if  (sub >= (* tree_size) - 1)
       {
        (* tree_size) *= EXPANSION_FACTOR;
PRALLOC ((* tree_size) * sizeof (DFS_Info_t));
        (* dfs_tree) = (DFS_Info_t *) safe_realloc ((* dfs_tree),
                           (* tree_size) * sizeof (DFS_Info_t));
       }

if  (Global_Debug_Flag)
    fprintf (stderr, "%*s From %d  so_far = %.0f:\n",
             4 * level + 1, "", from -> id, so_far);

   node [from -> id] . visited = TRUE;
   node [from -> id] . sub = sub;
   (* dfs_tree) [sub] . id = from -> id;
   (* dfs_tree) [sub] . distance = so_far;

   (* dfs_tree) [sub] . is_hit = FALSE;
   for  (i = 0;  i < num_targets;  i ++)
     if  (from -> id == target [i] . id)
         {
          int  hit_correct_end;

          switch  (target [i] . orient)
            {
             case  AB_AB :
             case  BA_AB :
               hit_correct_end = (from_end == B_END);
               break;
             case  AB_BA :
             case  BA_BA :
               hit_correct_end = (from_end == A_END);
               break;
             default :
               fprintf (stderr, "ERROR:  Bad target orientation\n");
               assert (FALSE);
            }

if  (Global_Debug_Flag)
    fprintf (stderr, "%*s   hit  lo = %.0f  hi = %.0f  so_far = %.0f  end = %c\n",
             4 * level + 1, "", target [i] . lo, target [i] . hi,
             so_far, hit_correct_end ? 'T' : 'F');

          if  (target [i] . lo <= so_far && so_far <= target [i] . hi
                 && hit_correct_end)
              {
               (* dfs_tree) [sub] . is_hit = TRUE;
               if  (to == NULL)
                   succeeded = TRUE;
              }
          target [i] . where = so_far;
          (* dfs_tree) [sub] . target_sub = i;
          break;
         }

   // Temporarily only allow target nodes in search
   if  (! (* dfs_tree) [sub] . is_hit && from != to && level > 0)
       return  FALSE;

   curr_max_hits = 0;
   curr_first = curr_max_first = curr_max_first_dist = -1;
   curr_max_first_total = 0.0;

   if  (from == to)
       {
        if  ((* dfs_tree) [sub] . is_hit)
            {
             (* first) = (* max_first) = sub;
             (* max_hits) = 1;
            }
          else
            {
             (* first) = (* max_first) = -1;
             (* max_hits) = 0;
            }

        (* max_first_dist) = 0;
        (* max_first_total) = curr_max_first_total;

        (* dfs_tree) [sub] . max_hits = curr_max_hits;
        (* dfs_tree) [sub] . first = curr_first;
        (* dfs_tree) [sub] . max_first = curr_max_first;
        (* dfs_tree) [sub] . max_first_dist = curr_max_first_dist;
        (* dfs_tree) [sub] . max_first_total = curr_max_first_total;
        (* dfs_tree) [sub] . succeeded = TRUE;
        node [from -> id] . finished = TRUE;

        return  TRUE;
       }
         
//  This node automatically goes in the tree.  If it is past the
//  distance bound, then don't pursue any of its edges.
//  bound  measures the distance to the *low* end of the  to  chunk
//  so that we don't give up too soon in the case of a containment.

   if  (so_far - from -> bpLength . mean > bound)
       {
        if  (to == NULL && (* dfs_tree) [sub] . is_hit)
            {
             (* first) = (* max_first) = sub;
             (* max_hits) = 1;
            }
          else
            {
             (* first) = (* max_first) = -1;
             (* max_hits) = 0;
            }

        (* max_first_dist) = curr_max_first_dist;
        (* max_first_total) = curr_max_first_total;

        (* dfs_tree) [sub] . max_hits = curr_max_hits;
        (* dfs_tree) [sub] . first = curr_first;
        (* dfs_tree) [sub] . max_first = curr_max_first;
        (* dfs_tree) [sub] . max_first_dist = curr_max_first_dist;
        (* dfs_tree) [sub] . max_first_total = curr_max_first_total;
        (* dfs_tree) [sub] . succeeded = succeeded;
        node [from -> id] . finished = TRUE;

        return  succeeded;
       }

   InitGraphEdgeIterator (ScaffoldGraph->RezGraph, from -> id, from_end,
                          ALL_EDGES, GRAPH_EDGE_DEFAULT, & edge_iterator);

   while  ((edge = NextGraphEdgeIterator (& edge_iterator)) != NULL)
     {
      ChunkInstanceT  * next_node;
      double  frag_len, new_max_first_total;
      int  new_first, new_max_hits, new_max_first, new_max_first_dist;
      int  next_id, progress;

if  (Global_Debug_Flag)
    fprintf (stderr, ">>> Trying edge idA = %d  idB = %d  mean = %.0f  %s %s %s %s\n",
             edge -> idA, edge -> idB,
             edge -> distance . mean,
             isOverlapEdge (edge) ? "" : "NotOlap",
             isProbablyBogusEdge (edge) ? "Bogus" : "",
             edge -> flags . bits . hasTandemOverlap ? "Tandem" : "",
             edge -> flags . bits . hasContainmentOverlap ? "Contain" : "");

      if  (! isOverlapEdge (edge)
             || isProbablyBogusEdge (edge))
          continue;
      if  ((edge_mask & SKIP_TANDEM_OLAPS)
             && (edge -> flags . bits . hasTandemOverlap
                 || edge -> distance . mean > 0))  // Probably a tandem repeat overlap
          continue;

      if  (edge -> idA == from -> id)
          next_id = edge -> idB;
        else
          next_id = edge -> idA;

if  (Global_Debug_Flag)
    fprintf (stderr, "%*s   next %d\n", 4 * level + 1, "", next_id);

      if  (node [next_id] . visited && ! node [next_id] . finished)
          continue;          // A back-edge, ignore it.

      next_node = GetGraphNode (ScaffoldGraph -> RezGraph, next_id);
      progress = next_node -> bpLength . mean + edge -> distance . mean;
      frag_len = next_node -> bpLength . mean;


      // Skip edges marked as contains if SKIP_CONTAINMENT_OLAPS indicates
      // so and if it's not the first or last step

      if  ((edge_mask & SKIP_CONTAINMENT_OLAPS)
             && edge -> flags . bits . hasContainmentOverlap
             && (to == NULL
                   || next_node -> id != to -> id
                   || (next_node -> id == to -> id
                         && next_node -> bpLength . mean
                              >= from -> bpLength . mean)
                   || level != 0))
          continue;


      // Require forward progress on each step, except possibly
      // for final step to destination

      if  ((to == NULL || next_node -> id != to -> id)
              && progress <= 0)
          continue;


      if  (node [next_id] . finished)
          {
           int  j;

           j = node [next_id] . sub;

           if  (! (* dfs_tree) [j] . succeeded)
               continue;

           succeeded = TRUE;
           if  ((* dfs_tree) [j] . is_hit)
               {
                new_first = j;
                new_max_hits = 1 + (* dfs_tree) [j] . max_hits;
                new_max_first = j;
                new_max_first_dist = progress;
                new_max_first_total = frag_len - edge -> distance . mean;
               }
             else
               {
                new_first = (* dfs_tree) [j] . first;
                new_max_hits = (* dfs_tree) [j] . max_hits;
                new_max_first = (* dfs_tree) [j] . max_first;
                new_max_first_dist
                    = progress + (* dfs_tree) [j] . max_first_dist;
                new_max_first_total
                    = frag_len - edge -> distance . mean
                        + (* dfs_tree) [j] . max_first_total;
               }
           if  (new_max_hits > curr_max_hits)
               {
                curr_max_hits = new_max_hits;
                curr_max_first = new_max_first;
                curr_max_first_dist = new_max_first_dist;
                curr_max_first_total = new_max_first_total;
               }
           else if  (curr_max_first_total == 0.0)
               {
                curr_max_first_dist = new_max_first_dist;
                curr_max_first_total = new_max_first_total;
               }
           if  (curr_first == -1
                  || (new_first != -1
                       && (* dfs_tree) [new_first] . distance
                             < (* dfs_tree) [curr_first] . distance))
               curr_first = new_first;
          }
        else
          {
           int  found, next_end;

           switch  (GetEdgeOrientationWRT (edge, from -> id))
             {
              case  AB_AB :
                next_end = B_END;
                break;
              case  AB_BA :
                next_end = A_END;
                break;
              case  BA_AB :
                next_end = B_END;
                break;
              case  BA_BA :
                next_end = A_END;
                break;
              default :
                assert (FALSE);
             }

           found = Depth_First_Visit
                      (next_node, next_end, to, num_targets, target, bound,
                       so_far + progress, & new_first, & new_max_hits,
                       & new_max_first, & new_max_first_dist,
                       & new_max_first_total, tree_sub,
                       tree_size, dfs_tree, node, level + 1,
                       edge_mask);

           if  (found)
               {
                if  (new_max_hits > curr_max_hits)
                    {
                     curr_max_hits = new_max_hits;
                     curr_max_first = new_max_first;
                     curr_max_first_dist = new_max_first_dist + progress;
                     curr_max_first_total = new_max_first_total + frag_len
                                              - edge -> distance . mean;
                    }
                else if  (curr_max_first_total == 0.0)
                    {
                     curr_max_first_dist = new_max_first_dist + progress;
                     curr_max_first_total = new_max_first_total + frag_len
                                              - edge -> distance . mean;
                    }
                if  (curr_first == -1
                       || (new_first != -1
                            && (* dfs_tree) [new_first] . distance
                                  < (* dfs_tree) [curr_first] . distance))
                    curr_first = new_first;
               }
           succeeded = (succeeded || found);
          }
     }

   if  ((* dfs_tree) [sub] . is_hit)
       {
        (* max_hits) = 1 + curr_max_hits;
        (* max_first) = sub;
        (* max_first_dist) = 0;
        (* max_first_total) = 0.0;
        (* first) = sub;
        if  (curr_max_hits > 0)
            target [(* dfs_tree) [sub] . target_sub] . next
                = (* dfs_tree) [curr_max_first] . target_sub;
       }
     else
       {
        (* max_hits) = curr_max_hits;
        (* max_first) = curr_max_first;
        (* max_first_dist) = curr_max_first_dist;
        (* max_first_total) = curr_max_first_total;
        (* first) = curr_first;
       }

   (* dfs_tree) [sub] . max_hits = curr_max_hits;
   (* dfs_tree) [sub] . first = curr_first;
   (* dfs_tree) [sub] . max_first = curr_max_first;
   (* dfs_tree) [sub] . max_first_dist = curr_max_first_dist;
   (* dfs_tree) [sub] . max_first_total = curr_max_first_total;
   (* dfs_tree) [sub] . succeeded = succeeded;
   node [from -> id] . finished = TRUE;

   return  succeeded;
  }



static int  Determine_Components
    (int list [], int num_list, Gap_Chunk_t * node [], int num_nodes,
     int start_sub, int target_sub, int edge [], Stone_Edge_t pool [],
     double ref_position, double factor, LengthT * target_position,
     int * complete_path, Gap_Fill_t * gap)

//  Set  keep  flag true for the longest sequence of nodes in
//  list [0 .. (num_list - 1)]  starting with  list [start_sub]  and going
//  to  list [target_sub] , unless  target_sub  is negative
//  in which case the path can go anywhere.  Use edges in
//  edge  and  pool .  Values in  list  are subscripts of entries
//  in  node [0 .. (num_nodes - 1)] .   ref_position  is the scaffold
//  position corresponding to the zero path position.   factor
//  is  +1.0  for forward paths and  -1.0  for reverse.
//  Set  (* target_position)  to the place where the end of the gap
//  should be, based on the path.
//  Set  (* complete_path)  TRUE  iff there is a single component from
//  start_sub  to  target_sub  or if one of those sub's is  -1  and
//  there is only a single component.
//  gap  is a pointer to the information for the gap containing these
//  nodes.

  {
   const int  LARGE_FALSE_VAL = 10000;
   const double  NUM_STD_DEVS = 5.0;
   static Component_Info_t  * comp_info = NULL;
   static Path_Info_t  * path_info = NULL;
   static int  path_info_size = 0;
   double  ref_pos, ref_var, hi_delta, olap_diff;
   double  needed_expansion, avail_expansion, gap_avail;
   int  a_hang_sum;
   int  num_good_components = 0;
   int  num_kept = 0, index, comp_num, total_kept;
   int  hit_start = FALSE, hit_target = FALSE;
   int  i, j;

   assert (num_list > 0);
   target_position -> mean = 0.0;
   (* complete_path) = FALSE;

   if  (num_nodes > path_info_size)
       {
PRALLOC (num_nodes * sizeof (Path_Info_t));
        path_info = (Path_Info_t *) safe_realloc
                        (path_info, num_nodes * sizeof (Path_Info_t));
PRALLOC (num_nodes * sizeof (Component_Info_t));
        comp_info = (Component_Info_t *) safe_realloc
                        (comp_info, num_nodes * sizeof (Component_Info_t));
        path_info_size = num_nodes;
       }

   for  (i = 0;  i < num_nodes;  i ++)
     {
      path_info [i] . path_len = node [i] -> link_ct;
      path_info [i] . from = -1;
      path_info [i] . hi_position = 0;
      path_info [i] . lo_position = 0;
      path_info [i] . a_hang = 0;
      path_info [i] . total_olap = 0;
      path_info [i] . hit = FALSE;
      path_info [i] . done = FALSE;
     }
   if  (start_sub >= 0)
       path_info [start_sub] . path_len = LARGE_FALSE_VAL;
           // Force paths to use start_sub if possible
           // Remember to subtract off this value if it matters

#if  VERBOSE
fprintf (stderr, "Determine_Components:\n");
#endif
   for  (i = 0;  i < num_list;  i ++)
     {
      for  (j = edge [list [i]];  j >= 0;  j = pool [j] . next)
        {
         assert (pool [j] . from == list [i]);

         if  (node [pool [j] . to] -> keep && ! pool [j] . bad)
             {
#if  SHOW_STONE_CONFIRM
              fprintf (stderr,
                       "From %4d (%5d) to %4d (%5d)  ahg = %6d  prog = %6d  qual = %6.4f",
                       pool [j] . from, node [pool [j] . from] -> chunk_id,
                       pool [j] . to, node [pool [j] . to] -> chunk_id,
                       pool [j] . a_hang, pool [j] . progress,
                       pool [j] . quality);
#endif
              if  (path_info [pool [j] . to] . done)
                  {
#if  SHOW_STONE_CONFIRM
                   fprintf (stderr, "  forms a cycle, ignored\n");
#endif
                   continue;
                  }

              if  (path_info [pool [j] . from] . path_len
                     + node [pool [j] . to] -> link_ct
                         >= path_info [pool [j] . to] . path_len)
                  {
                   LengthT  from_near_end, to_near_end;
                   double  desired_a_hang, error_band, std_dev;

                   if  ((node [pool [j] . from] -> flipped && factor >= 0.0)
                         || (! node [pool [j] . from] -> flipped && factor < 0.0))
                       from_near_end = node [pool [j] . from] -> end;
                     else
                       from_near_end = node [pool [j] . from] -> start;
                   if  ((node [pool [j] . to] -> flipped && factor >= 0.0)
                         || (! node [pool [j] . to] -> flipped && factor < 0.0))
                       to_near_end = node [pool [j] . to] -> end;
                     else
                       to_near_end = node [pool [j] . to] -> start;
                   desired_a_hang
                       = factor * (to_near_end . mean - from_near_end . mean);

                   std_dev = sqrt (from_near_end . variance)
                               + sqrt (to_near_end . variance);
                   error_band = NUM_STD_DEVS * std_dev;
                   if  (pool [j] . a_hang < desired_a_hang - error_band
                          || pool [j] . a_hang > desired_a_hang + error_band)
                       {
#if  SHOW_STONE_CONFIRM
                        fprintf (stderr,
                                 "  missed  %6d not in [%6.0f,%6.0f]  %4.2f stdvs\n",
                                 pool [j] . a_hang,
                                 desired_a_hang - error_band,
                                 desired_a_hang + error_band,
                                 (pool [j] . a_hang - desired_a_hang) / std_dev);
#endif
                        continue;
                       }

#if  SHOW_STONE_CONFIRM
                   fprintf (stderr,
                            "  %6d in [%6.0f,%6.0f]  %4.2f stdvs\n",
                            pool [j] . a_hang,
                            desired_a_hang - error_band,
                            desired_a_hang + error_band,
                            (pool [j] . a_hang - desired_a_hang) / std_dev);
#endif

                   assert (pool [j] . progress >= 0);
                   if  (! Prior_Olaps_OK (pool [j] . from, pool [j] . to,
                                          path_info [pool [j] . from] . hi_position
                                            - pool [j] . length,
                                          path_info, edge, pool))
                       {
#if  SHOW_STONE_CONFIRM
                        fprintf (stderr, "  failed prior olaps\n");
#endif
                        continue;
                       }

                   path_info [pool [j] . to] . path_len
                       = path_info [pool [j] . from] . path_len
                           + node [pool [j] . to] -> link_ct;
                   path_info [pool [j] . to] . from = pool [j] . from;
                   path_info [pool [j] . to] . a_hang = pool [j] . a_hang;
                   path_info [pool [j] . to] . hi_position
                       = path_info [pool [j] . from] . hi_position
                           + pool [j] . progress;
                   path_info [pool [j] . to] . lo_position
                       = path_info [pool [j] . from] . hi_position
                           - pool [j] . length;
                   path_info [pool [j] . to] . total_olap
                       = path_info [pool [j] . from] . total_olap
                           + 2 * node [pool [j] . to] -> len
                           - pool [j] . progress;
                  }
#if  SHOW_STONE_CONFIRM
              fprintf (stderr, "  path_len = %2d  from = %2d\n",
                       path_info [pool [j] . to] . path_len,
                       path_info [pool [j] . to] . from);
#endif
             }
        }
      path_info [list [i]] . done = TRUE;
     }   

// list all components
{
 int  hi_unhit;
 int  ct;

 if  (target_sub >= 0)
     {
      ct = 0;
      for  (i = target_sub;  i >= 0;  i = path_info [i] . from)
        {
         ct ++;
         path_info [i] . hit = TRUE;
#if  VERBOSE
         fprintf (stderr,
            "%3d (%5d)  hi_pos = %6d  AEnd = (%6.0f,%6.0f)  BEnd = (%6.0f,%6.0f)\n",
                  i, node [i] -> chunk_id, path_info [i] . hi_position,
                  node [i] -> start . mean, sqrt (node [i] -> start . variance),
                  node [i] -> end . mean, sqrt (node [i] -> end . variance));
#endif
         if  (i == start_sub)
             hit_start = TRUE;
        }

      if  (hit_start)
          {
           (* complete_path) = TRUE;
           return  ct - 2;    // not counting start and target
          }

      if  (ct > 1)
          {
           num_kept += ct - 1;     // not counting target
           comp_info [num_good_components] . sub = target_sub;
           comp_info [num_good_components] . ct = ct;   // including target
           comp_info [num_good_components] . has_target = TRUE;
           comp_info [num_good_components] . has_start = FALSE;
           num_good_components ++;
           hit_target = TRUE;
          }
     }

 hi_unhit = -1;
 for  (i = 0;  i < num_list;  i ++)
   if  (! path_info [list [i]] . hit
          && (hi_unhit < 0 || path_info [hi_unhit] . path_len
                                < path_info [list [i]] . path_len))
       hi_unhit = list [i];

 while  (hi_unhit >= 0)
   {
    int  is_start_component = FALSE;

    for  (i = hi_unhit;  i >= 0;  i = path_info [i] . from)
      if  (path_info [i] . hit)
          break;
    if  (i < 0)
        {
#if  VERBOSE
         fprintf (stderr, "Component:  links = %d\n",
                  path_info [hi_unhit] . path_len);
#endif

         ct = 0;
         for  (i = hi_unhit;  i >= 0;  i = path_info [i] . from)
           {
#if  VERBOSE
            fprintf (stderr,
               "%3d (%5d)  hi_pos = %6d  AEnd = (%6.0f,%6.0f)  BEnd = (%6.0f,%6.0f)\n",
                     i, node [i] -> chunk_id, path_info [i] . hi_position,
                     node [i] -> start . mean, sqrt (node [i] -> start . variance),
                     node [i] -> end . mean, sqrt (node [i] -> end . variance));
#endif
            if  (i == start_sub && i != hi_unhit)
                is_start_component = TRUE;
              else
                ct ++;
           }
         if  (ct > 1 || is_start_component
                || (hi_unhit != start_sub && path_info [hi_unhit] . path_len > 1))
             {
              num_kept += ct;
              comp_info [num_good_components] . sub = hi_unhit;
              comp_info [num_good_components] . ct = ct;
              if  (is_start_component)
                  comp_info [num_good_components] . ct ++;   // count start, too
              comp_info [num_good_components] . has_target = FALSE;
              comp_info [num_good_components] . has_start = is_start_component;
              num_good_components ++;
              if  (is_start_component)
                  hit_start = TRUE;
             }
        }
    assert (! path_info [hi_unhit] . hit);
    for  (i = hi_unhit;  i >= 0 && ! path_info [i] . hit;  i = path_info [i] . from)
      path_info [i] . hit = TRUE;

    hi_unhit = -1;
    for  (i = 0;  i < num_list;  i ++)
      if  (! path_info [list [i]] . hit
             && (hi_unhit < 0 || path_info [hi_unhit] . path_len
                                   < path_info [list [i]] . path_len))
          hi_unhit = list [i];
   }

 if  (num_good_components == 1
        && (hit_start || hit_target))
     (* complete_path) = TRUE;

#if  VERBOSE
 fprintf (stderr, "num_good_components = %d  num_kept = %d  complete_path = %c\n",
          num_good_components, num_kept, (* complete_path) ? 'T' : 'F');
#endif

 if  ((* complete_path) || start_sub == -1 || target_sub == -1
        || num_kept <= 0)
     return  num_kept;    // handle usual way for now
}

   for  (i = 0;  i < num_nodes;  i ++)
     node [i] -> keep = FALSE;


   total_kept = 0;
   avail_expansion = 5.0 * sqrt (fabs (gap -> end . variance
                                         - gap -> start . variance));

   gap_avail = gap -> len;
   gap -> adjustment . mean = 0.0;
   gap -> adjustment . variance = 0.0;
   for  (comp_num = 0;  comp_num < num_good_components;  comp_num ++)
     {
      num_kept = 0;
      for  (i = comp_info [comp_num] . sub;  i >= 0;  i = path_info [i] . from)
        {
         num_kept ++;
#if  SHOW_STONE_CONFIRM
         fprintf (stderr,
         "%3d (%5d)  hi_position = %6d  AEnd = (%6.0f,%6.0f)  BEnd = (%6.0f,%6.0f)\n",
                  i, node [i] -> chunk_id, path_info [i] . hi_position,
                  node [i] -> start . mean, sqrt (node [i] -> start . variance),
                  node [i] -> end . mean, sqrt (node [i] -> end . variance));
#endif
        }

      assert  (num_kept >= 1);

//  Create list of nodes in order by low position and set
//  that value in  path_info  based on  a_hang  fields in edges.

      index = num_kept;
      for  (i = comp_info [comp_num] . sub;  i >= 0;  i = path_info [i] . from)
        list [-- index] = i;
      assert (index == 0);

      for  (i = 1;  i < num_kept;  i ++)
        {
         int  save = list [i];
         int  last_neg_ahang = 0;
         int  missing_olap = FALSE;
         int  k = 0;

         for  (j = i - 1;  j >= 0;  j --)
           {
            int  m;

            // Find edge from  list [j]  to  list [i]
            for  (k = edge [list [j]];  k >= 0;  k = pool [k] . next)
              if  (pool [k] . to == save)
                  break;
            if  (k < 0)
                {
                 int  sum = 0;

                 fprintf (stderr, "YIKES!  No edge from %d to %d\n",
                          node [list [j]] -> chunk_id, node [save] -> chunk_id);
                 if  (j != i - 1)
                     {
                      fprintf (stderr,
                         "YOUCH! Inconsistent overlaps.  Leave stones alone.\n");
                      for  (m = 0;  m < num_kept;  m ++)
                        node [list [m]] -> keep = FALSE;
                      return  0;
                     }
                 for  (m = j;  m >= 0 && list [m] != path_info [save] . from;
                         m --)
                   sum += path_info [list [m]] . a_hang;
                 assert (m >= 0);
                 path_info [save] . a_hang -= sum;
                 assert (path_info [save] . a_hang >= 0);
                 missing_olap = TRUE;
                 break;
                }

            // if  a_hang is positive or back to first chunk,  stop
            // and put  i  after  j  in  list; otherwise, keep going

            if  (j == 0 || pool [k] . a_hang >= 0)
                break;

            last_neg_ahang = pool [k] . a_hang;
            list [j + 1] = list [j];
           }

         // when stop, save best edge from  j  to  i  (now at  j + 1)
         // and  i  to previous  j + 1
         // this is like an insertion sort

         if  (! missing_olap)
             {
              list [j + 1] = save;
              path_info [list [j + 1]] . a_hang = pool [k] . a_hang;
              if  (j < i - 1)
                  path_info [list [j + 2]] . a_hang = - last_neg_ahang;
             }
        }

#if  SHOW_STONE_CONFIRM
{
 int  i;

 fprintf (stderr, "Component revised path:\n");
 for  (i = 0;  i < num_kept;  i ++)
   {
    fprintf (stderr,
             "%3d %3d (%5d)  a_hang = %5d  hi_position = %6d\n",
               i, list [i], node [list [i]] -> chunk_id,
               path_info [list [i]] . a_hang,
               path_info [list [i]] . hi_position);
   }
}
#endif

      if  (comp_info [comp_num] . has_target)
          {
           // adjust from high position back to low based on overlaps

           a_hang_sum = 0;
           for  (i = 0;  i < num_kept;  i ++)
             a_hang_sum += path_info [list [i]] . a_hang;
           needed_expansion = Max_double (0.0, a_hang_sum - gap_avail);
#if  SHOW_STONE_CONFIRM
           fprintf (stderr,
               "target comp = %d  needed_expansion = %.0f  avail_expansion = %.0f"
               "  gap_avail = %.0f  used = %d\n",
               comp_num, needed_expansion, avail_expansion, gap_avail,
               a_hang_sum);
#endif
           if  (needed_expansion > avail_expansion)
               continue;      // doesn't fit

           avail_expansion -= needed_expansion;
           gap_avail -= a_hang_sum - needed_expansion;

           gap -> adjustment . mean += needed_expansion;
           if  (gap -> len >= 0)
               {
                ref_pos = gap -> end . mean + needed_expansion
                            - path_info [list [num_kept - 1]] . a_hang;
                ref_var = gap -> end . variance - gap -> ref_variance;
               }
             else
               {
                ref_pos = gap -> start . mean + needed_expansion
                            - path_info [list [num_kept - 1]] . a_hang;
                ref_var = gap -> start . variance - gap -> ref_variance;
               }
           hi_delta = Max_double (node [list [num_kept - 1]] -> start . mean,
                                  node [list [num_kept - 1]] -> end . mean)
                        + needed_expansion
                        - path_info [list [num_kept - 1]] . hi_position;
           for  (i = num_kept - 2;  i >= 0;  i --)
             {
              olap_diff = path_info [list [i]] . total_olap;
              if  (i > 0)
                  olap_diff -= path_info [list [i - 1]] . total_olap;
              ref_var -= ComputeFudgeVariance (olap_diff);

              if  (node [list [i]] -> flipped)
                  {
                   node [list [i]] -> end . mean = ref_pos;
                   node [list [i]] -> end . variance = ref_var;
                   node [list [i]] -> start . mean
                       = path_info [list [i]] . hi_position
                           + hi_delta;
                   node [list [i]] -> start . variance
                       = node [list [i]] -> end . variance
                           + ComputeFudgeVariance (node [list [i]] -> len);
                  }
                else
                  {
                   node [list [i]] -> start . mean = ref_pos;
                   node [list [i]] -> start . variance = ref_var;
                   node [list [i]] -> end . mean
                       = path_info [list [i]] . hi_position
                           + hi_delta;
                   node [list [i]] -> end . variance
                       = node [list [i]] -> start . variance
                           + ComputeFudgeVariance (node [list [i]] -> len);
                  }
              node [list [i]] -> keep = TRUE;
              ref_pos -= path_info [list [i]] . a_hang;
             }
          }
      else if  (comp_info [comp_num] . has_start)
          {
           needed_expansion
               = Max_double (0.0,
                             path_info [comp_info [comp_num] . sub] . hi_position
                               - gap_avail);
#if  SHOW_STONE_CONFIRM
           fprintf (stderr,
                    "start comp = %d  needed_expansion = %.0f  avail_expansion = %.0f"
                    "  gap_avail = %.0f  used = %d\n",
                    comp_num, needed_expansion, avail_expansion, gap_avail,
                    path_info [comp_info [comp_num] . sub] . hi_position);
#endif
           if  (needed_expansion > avail_expansion)
               continue;      // doesn't fit

           avail_expansion -= needed_expansion;
           gap_avail -= path_info [comp_info [comp_num] . sub] . hi_position
                          - needed_expansion;

           gap -> adjustment . mean += needed_expansion;

           // adjust positions of target component nodes if necessary
           if  (needed_expansion > 0.0 && comp_num > 0)
               for  (i = comp_info [0] . sub;  i >= 0;
                         i = path_info [i] . from)
                 {
                  node [i] -> start . mean += needed_expansion;
                  node [i] -> end . mean += needed_expansion;
                 }

           Adjust_Stone_Positions
               (list, num_kept, node, ref_position, factor, path_info,
                -1, target_position);

           for  (i = 1;  i < num_kept;  i ++)
             node [list [i]] -> keep = TRUE;
          }
        else
          {
           // set left position
#if  VERBOSE
           for  (i = 0;  i < num_kept;  i ++)
             fprintf (stderr, "left = [%6.0f,%6.0f]  lo_position = %d\n",
                      node [list [i]] -> flipped ?
                        node [list [i]] -> end . mean :
                        node [list [i]] -> start . mean,
                      sqrt (Min_double (node [list [i]] -> start . variance,
                                        node [list [i]] -> start . variance)),
                      path_info [list [i]] . lo_position);
#endif
           num_kept = 0;  // temporary until do other cases
          }

#if  SHOW_STONE_CONFIRM
{
 int  i;

 fprintf (stderr, "New positions:\n");
 for  (i = 0;  i < num_kept;  i ++)
   {
    fprintf (stderr,
             "%3d %3d (%5d)  start = (%6.0f,%6.0f)  end = (%6.0f,%6.0f)\n",
               i, list [i], node [list [i]] -> chunk_id,
               node [list [i]] -> start . mean,
               sqrt (node [list [i]] -> start . variance),
               node [list [i]] -> end . mean,
               sqrt (node [list [i]] -> end . variance));
   }
}
#endif

      total_kept += num_kept;

//      break;  // temporary so just does back from target ones
     }

   return  total_kept;
  }



static void  DFS_Stone_Visit
    (int i, Gap_Chunk_t * node [], int edge [], Stone_Edge_t pool [],
     int * bad_edge)

//  Visit all nodes marked  keep  and not yet visited starting at
//  node [i]  using edges in  edge  and  pool .
//  Set  (* bad_edge)  true if the edge to this node is a backedge
//  in the DFS tree.

  {
   int  j;

   (* bad_edge) = node [i] -> visited && ! node [i] -> finished;

   if  (node [i] -> visited || ! (node [i] -> keep))
       return;

   node [i] -> visited = TRUE;
   for  (j = edge [i];  j >= 0;  j = pool [j] . next)
     {
      int  mark_edge;

      assert (i == pool [j] . from);
      DFS_Stone_Visit (pool [j] . to, node, edge, pool, & mark_edge);
      if  (mark_edge)
          pool [j] . bad = TRUE;
     }
   node [i] -> finished = TRUE;

   return;
  }



static void  Disqualify_Scaff_Chunks
    (Scaffold_Fill_t * fill_chunks)

//  Set the  keep  flag to  FALSE  for any scaffold chunks
//  in  fill_chunks .

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     Disqualify_Scaff_Chunks_One_Scaffold (fill_chunks, scaff_id);

   return;
  }



static void  Disqualify_Scaff_Chunks_One_Scaffold
    (Scaffold_Fill_t * fill_chunks, int scaff_id)

//  Set the  keep  flag to  FALSE  for any scaffold chunks
//  in  fill_chunks [scaff_id] .

  {
   int  j;

   for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
     {
      int  k;
      Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;

      for  (k = 0;  k < this_gap -> num_chunks;  k ++)
        {
         Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

//            if  (REF (this_chunk -> chunk_id) . scaff_id != NULLINDEX)
         if  (this_chunk -> copy_letter == GAP_END_CHAR)
             this_chunk -> keep = FALSE;
        }
     }

   return;
  }



static void  Doublecheck_Positions
    (Scaffold_Fill_t * fill_chunks, int make_adjustments)

//  Make sure the variances of all entries in  fill_chunks
//  are non-decreasing, and wherever an overlap is indicated,
//  confirm it in hash table.  If  make_adjustments  is true, actually
//  make the changes; otherwise, just report anomalies.

  {
   int  scaff_id;

   fprintf (stderr, "### Doublecheck_Positions ###\n");

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     Doublecheck_Positions_One_Scaffold
         (fill_chunks, make_adjustments, scaff_id);

   return;
  }



static void  Doublecheck_Positions_One_Scaffold
    (Scaffold_Fill_t * fill_chunks, int make_adjustments, int scaff_id)

//  Make sure the variances of all entries in  fill_chunks [scaff_id]
//  are non-decreasing, and wherever an overlap is indicated,
//  confirm it in hash table.  If  make_adjustments  is true, actually
//  make the changes; otherwise, just report anomalies.

  {
   static Gap_Chunk_t  * * check = NULL;
   static double  * adjustment = NULL;
   static int  check_size = 0;
   double  var_adjust, prev_adjust;
   int  j;

   if  (check_size == 0)
       {
        check_size = INITIAL_GAP_ENTRIES;
PALLOC (check_size * sizeof (Gap_Chunk_t *));
        check = (Gap_Chunk_t * *)
                     safe_malloc (check_size * sizeof (Gap_Chunk_t *));
PALLOC (check_size * sizeof (double));
        adjustment = (double *)
                     safe_malloc (check_size * sizeof (double));
       }


   // First make sure any overlaps implied by computed positions
   // are really there. 

   for  (j = GAPS_TO_ADJUST)
     {
      Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;

      Check_Olaps (this_gap);
     }


   // Now check for dips in variance and negative variances

   prev_adjust = var_adjust = 0.0;

   for  (j = GAPS_TO_ADJUST)
     {
      ChunkInstanceT  * left_chunk, * right_chunk;
      Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
      double  next_lo, prev_hi, lo_var, hi_var;
      int  k, ct;

      ct = 0;

      if  (this_gap -> num_chunks > 0)
          {
           if  (this_gap -> num_chunks > check_size)
               {
                check_size *= 2;
PRALLOC (check_size * sizeof (Gap_Chunk_t *));
                check = (Gap_Chunk_t * *)
                            safe_realloc
                                (check,
                                 check_size * sizeof (Gap_Chunk_t *));
PRALLOC (check_size * sizeof (double));
                adjustment = (double *)
                            safe_realloc
                                (adjustment,
                                 check_size * sizeof (double));
               }

           for  (k = 0;  k < this_gap -> num_chunks;  k ++)
             {
              Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

              if  (this_chunk -> keep)
                  {
                   adjustment [ct] = 0.0;
                   check [ct ++] = this_chunk;
                  }
             }

           qsort (check, ct, sizeof (Gap_Chunk_t *), By_Low_Position);
//              qsort (check, ct, sizeof (Gap_Chunk_t *), By_High_Position);
          }


      if  (j <= 0)
          {
           prev_hi = 0.0;
           left_chunk = NULL;
          }
        else
          {
           left_chunk = GetGraphNode (ScaffoldGraph -> RezGraph,
                                       this_gap -> left_cid);
           prev_hi = Max_double (left_chunk -> offsetAEnd . variance,
                                 left_chunk -> offsetBEnd . variance);
          }

      hi_var = 0.0;
      for  (k = 0;  k < ct;  k ++)
        {
         hi_var = Max_double (check [k] -> start . variance,
                              check [k] -> end . variance);
         lo_var = Min_double (check [k] -> start . variance,
                              check [k] -> end . variance);
         if  (lo_var < 0.0)
             {
              fprintf (stderr,
              "Doublecheck ERROR:  variance = %f  chunk %d  scaff %d\n",
                       lo_var, check [k] -> chunk_id, check [k] -> scaff_id);
              var_adjust = Max_double (var_adjust, - lo_var);
             }
         if  (lo_var < prev_hi)
             {
              fprintf (stderr,
              "Doublecheck ERROR:  variance dip = %f  chunk %d  scaff %d\n",
                       lo_var - prev_hi, check [k] -> chunk_id,
                       check [k] -> scaff_id);
              var_adjust +=  prev_hi - lo_var;
             }
         adjustment [k] = var_adjust;
         prev_hi = hi_var;
        }

      if  (j < fill_chunks [scaff_id] . num_gaps - 1)
          {
           right_chunk = GetGraphNode (ScaffoldGraph->RezGraph,
                                       this_gap -> right_cid);
           next_lo = Min_double (right_chunk -> offsetAEnd . variance,
                                 right_chunk -> offsetBEnd . variance);
           if  (next_lo < hi_var)
               {
                fprintf (stderr,
                "Doublecheck ERROR:  variance dip = %f  chunk %d  scaff %d\n",
                         next_lo - hi_var,  right_chunk -> id,
                         scaff_id);
                var_adjust +=  hi_var - next_lo;
               }
          }

      if  (var_adjust > 0.0 && make_adjustments)
          {
           //fprintf (stderr, ">>> Adjusting insertion variances <<<\n");
           for  (k = 0;  k < ct;  k ++)
             {
              check [k] -> start . variance += adjustment [k];
              check [k] -> end . variance += adjustment [k];
             }
          }

      if  (prev_adjust > 0.0 && make_adjustments && left_chunk != NULL)
          {
           //fprintf (stderr, ">>> Adjusting scaffold variances <<<\n");
           left_chunk -> offsetAEnd . variance += prev_adjust;
           left_chunk -> offsetBEnd . variance += prev_adjust;
          }

      prev_adjust = var_adjust;
     }

   return;
  }



static void  Eliminate_Encumbered_Uniques
    (Scaffold_Fill_t * fill)

//  Set  keep  flag to false of any unique (i.e. previously scaffolded)
//  entry in  fill  that has entries scheduled to be added to it.

  {
   int  scaff_id;
   int  change_made;


   // Count the number of entries to be added to each scaffold
   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     fill [scaff_id] . keep_ct = fill [scaff_id] . added_to_ct = 0;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
         int  k;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                fill [scaff_id] . keep_ct ++;
           }
        }
     }


   // Eliminate from insertion any scaffold that has at least
   // two entries to be inserted into it
   do
     {
      change_made = FALSE;

      for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
        {
         int  j;

         for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
           {
            Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
            int  k;

            for  (k = 0;  k < this_gap -> num_chunks;  k ++)
              {
               Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

               if  (this_chunk -> keep)
                   {
                    int  sid = REF (this_chunk -> chunk_id) . scaff_id;

                    if  (sid >= 0
                           && fill [sid] . keep_ct >= 2)
                        {
                         this_chunk -> keep = FALSE;
                         fill [scaff_id] . keep_ct --;
fprintf (stderr,
         "Removed chunk %d from insertion in scaff %d gap %d\n"
         "  was scaff %d with %d items\n",
         this_chunk -> chunk_id, scaff_id, j, sid, fill [sid] . keep_ct);
                         change_made = TRUE;
                        }
                   }
              }
           }
        }
     }  while  (change_made);


   // Eliminate one of any pair of scaffolds that each is trying
   // to insert itself in the other.
   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
         int  k;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                {
                 int  sid = REF (this_chunk -> chunk_id) . scaff_id;

                 if  (sid >= 0)
                     {
                      fill [sid] . added_to_ct ++;
                      fill [sid] . added_to_id = scaff_id;
                     }
                }
           }
        }
     }

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
         int  k;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                {
                 int  sid = REF (this_chunk -> chunk_id) . scaff_id;

                 if  (sid >= 0)
                     {
                      if  (fill [sid] . added_to_ct > 1
                             || fill [sid] . keep_ct > 1
                             || (fill [sid] . keep_ct == 1
                                  && (fill [scaff_id] . added_to_id != sid
                                        || fill [scaff_id] . added_to_ct != 1
                                        || scaff_id < sid)))
                          this_chunk -> keep = FALSE;
                     }
                }
           }
        }
     }

   return;
  }



static void  Eliminate_Encumbered_Uniques_One_Scaffold
    (Scaffold_Fill_t * fill, int scaff_id)

//  Set  keep  flag to false of any unique (i.e. previously scaffolded)
//  entry in  fill [scaff_id]  that has entries scheduled to be added to it.
//  If don't set  keep  flag false, mark that unique as unavailable
//  for future additions.
//  Note:  this is weaker than the versions that does all of  fill
//    at the same time.

  {
   int  j;

   for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
     {
      Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
      int  k;

      for  (k = 0;  k < this_gap -> num_chunks;  k ++)
        {
         Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

         if  (this_chunk -> keep)
             {
              int  sid = REF (this_chunk -> chunk_id) . scaff_id;

              if  (fill [scaff_id] . added_to_ct > 0)
                  this_chunk -> keep = FALSE;

              if  (sid >= 0)
                  {
                   if  (Num_Keep_Entries (fill, sid) > 0)
                       this_chunk -> keep = FALSE;
                     else
                       {
                        fill [sid] . added_to_ct ++;

                        //  uniques can only be inserted in one place (and hence shouldn't
                        //  use a surrogate)
                        if  (fill [sid] . added_to_ct > 1)
                            this_chunk -> keep = FALSE;
                       }
                  }
             }
        }
     }

   return;
  }



static int  Estimate_Chunk_Ends
    (Stack_Entry_t * stack, int stack_top,
     LengthT * left_end, LengthT * right_end, ChunkInstanceT * chunk,
     float * edge_quality, Scaffold_Fill_t * fill_chunks,
     int * gap, int * scaff_id, int * allowed_bad_links)

//  Set  (* left_end)  and  (* right_end)  to best estimates of
//  the ends  (* chunk)  based on the edges in  stack [0 .. (stack_top - 1)] .
//  Variances and standard deviations are relative to the highest relative
//  position unique chunk in  stack  that is to the left of this chunk.
//  If none, then variances are relative to the lowest relative position
//  unique chunk in  stack .  Return whether the estimated positions
//  are consistent with the positions implied by all the edges in  stack .
//  They are regarded as consistent iff the 3-std-deviation intervals
//  around the mean intersect.  Set  (* edge_quality)  to the average
//  quality of the edges that determine the end positions.
//  Use  fill_chunks  to extract scaffold positions and set  (* gap)
//  and  (* scaff_id)  to the relative gap position and scaffold number
//  into which this chunk should be placed if the computed coordinates
//  are consistent with the edges on the stack.
//  allowed_bad_links  is the number of bad links that can be tolerated
//  (generally 0 or 1) given a sufficient number of good links.

  {
   Gap_Fill_t  * g;
   Gap_Chunk_t  new_chunk;
   double  max_left_variance, min_right_variance, ref_variance;
   double  center, min_dist, min_variance, max_variance, new_ref_variance;
   int  got_new_bad_links;
   int  i, is_OK, bad_links, closest, edge_ct, good_links;

   max_left_variance = - 1.0;
   min_right_variance = DBL_MAX;
   min_variance = DBL_MAX;
   max_variance = - 1.0;

   assert (0 < stack_top && stack_top < STACK_SIZE);

   bad_links = 0;
   do
     {
      int  good_edge = -1;

      edge_ct = 0;
      for  (i = 0;  i < stack_top;  i ++)
        {
         if  (stack [i] . is_bad)
             continue;
         if  (stack [i] . source_variance < min_variance)
             min_variance = stack [i] . source_variance;
         if  (stack [i] . source_variance > max_variance)
             max_variance = stack [i] . source_variance;
         if  (stack [i] . left_link
                  && stack [i] . source_variance > max_left_variance)
             max_left_variance = stack [i] . source_variance;
         else if  (! stack [i] . left_link
                       && stack [i] . source_variance < min_right_variance
                       && max_left_variance == - 1.0)
             min_right_variance = stack [i] . source_variance;
         edge_ct ++;
         good_edge = i;
        }

if  (good_edge < 0)
    {
     fprintf (stderr, "will die on chunk %d\n", chunk -> id);
    }
      assert (good_edge >= 0);

      // Move the good edge to slot 0 in the stack for future reference
      if  (good_edge != 0)
          {
           Stack_Entry_t  save;

           save = stack [0];
           stack [0] = stack [good_edge];
           stack [good_edge] = save;
          }

      assert (min_right_variance >= 0.0);
      if  (edge_ct == 0)
          {
           (* edge_quality) = 0.0;
	   fprintf(stderr,"* Chunk %d: no edges\n", chunk->id);
           return  FALSE;
          }

      if  (max_left_variance > - 1.0)
          ref_variance = max_left_variance;
        else
          ref_variance = min_right_variance;

      assert (0.0 <= ref_variance && ref_variance <= FLT_MAX);

      Calc_End_Coords (stack, stack_top, left_end, right_end, chunk,
                       ref_variance);

      (* edge_quality) = 0.0;
      is_OK = TRUE;

      new_chunk . chunk_id = chunk -> id;
      assert(new_chunk.chunk_id != NULLINDEX);

      new_chunk . scaff_id = REF (stack [0] . chunk_id) . scaff_id;

      if  (stack [0] . flipped)
          {
           new_chunk . start = * right_end;
           new_chunk . end = * left_end;
          }
        else
          {
           new_chunk . start = * left_end;
           new_chunk . end = * right_end;
          }

      got_new_bad_links = FALSE;
      good_links = 0;
      for  (i = 0;  i < stack_top;  i ++)
        {
         Gap_Chunk_t  scaff_chunk;

         if  (stack [i] . is_bad)
             continue;
         scaff_chunk . chunk_id = stack [i] . chunk_id;
         scaff_chunk . scaff_id = REF (stack [i] . chunk_id) . scaff_id;
         scaff_chunk . start = REF (stack [i] . chunk_id) . a_end;
         scaff_chunk . end = REF (stack [i] . chunk_id) . b_end;

         if  (Single_Fragment_Only
                  || Is_Edge_Consistent (stack [i] . edge, & new_chunk, & scaff_chunk))
             {
              (* edge_quality) += CIEdge_Quality (stack [i] . edge);
              good_links += stack [i] . num_good_mates;
             }
           else
             {
              bad_links += stack [i] . num_good_mates;
              stack [i] . is_bad = TRUE;
              got_new_bad_links = TRUE;
             }
        }

      if  (bad_links > (* allowed_bad_links)
             || (bad_links > 0 && good_links < GOOD_LINKS_IF_BAD))
          {
           (* edge_quality) = 0.0;
           (* allowed_bad_links) = 0;
	   /*	   fprintf(stderr,"* Chunk %d: Too many bad links %d or not enough good links %d\n",
		   chunk->id, bad_links, good_links);
	   */
           return  FALSE;
          }
     }  while  (got_new_bad_links);


#if  TEST_HOPELESS_SCAFFS
   if  (Is_Hopeless_Scaff [new_chunk . scaff_id] & Hopeless_True_Mask)
       return  FALSE;
#endif

   (* edge_quality) /= edge_ct;
   g = fill_chunks [new_chunk . scaff_id] . gap;

   closest = 0;
   center = (left_end -> mean + right_end -> mean) / 2.0;
   min_dist = fabs (center);  // gap zero must start at coord zero

   for  (i = 1;  i < fill_chunks [new_chunk . scaff_id] . num_gaps;  i ++)
     {
      double  dist;

      dist = fabs (center - g [i] . start . mean);
      if  (dist < min_dist)
          {
           min_dist = dist;
           closest = i;
          }

      dist = fabs (center - g [i] . end . mean);
      if  (dist < min_dist)
          {
           min_dist = dist;
           closest = i;
          }
     }
   new_ref_variance = g [closest] . ref_variance;

   if  (new_ref_variance != ref_variance)
       Calc_End_Coords (stack, stack_top, left_end, right_end, chunk,
                        new_ref_variance);

   Fixup_Chunk_End_Variances (left_end, right_end, chunk -> bpLength . variance);

   (* gap) = closest;
   (* scaff_id) = new_chunk . scaff_id;

   return  TRUE;
  }



static void  Fasta_Print
    (FILE * fp, char * s, char * label)

//  Print sequence  s  in FASTA format to file  fp  with
//  string  label  on the '>' line at the beginning.

  {
   int  i = 0;
   fprintf (fp, ">%s\n", label);

   while  (* s != '\0')
     {
      fputc (* s, fp);
      s ++;
      i ++;
      if  (i == 60)
          {
           fputc ('\n', fp);
           i = 0;
          }
     }

   if  (i != 0)
       fputc ('\n', fp);
   
   return;
  }



int  Fill_Gaps
    (Global_CGW * data, char * prefix, int level, int redo_index)

//  Assign unresolved fragments to scaffolds.  This is the main entry
//  point from the chunk-graph-walker module.
//  (* data)  is no longer
//  used--the information needed comes from global variables.
//  prefix  is the prefix of the name to use for output files
//  level  indicates which steps to perform.
//  redo_index  is the number of previous calls to this routine without
//  renumbering or merging scaffolds

  {
   FILE  * log_file;
   char  filename [1000], iter_string [20];
   Scaffold_Fill_t  * fill_chunks;
   static int  iteration = 0;
   clock_t  start_time, stop_time;
   time_t  now;
   int  i;
   int  inserted = 0;

   StartTimerT(&GlobalData->GapFillTimer);

   Num_Scaffolds = GetNumGraphNodes (ScaffoldGraph -> ScaffoldGraph);
   if  (Num_Scaffolds == 0)
     return 0;
     
#if  TEST_HOPELESS_SCAFFS
   Hopeless_False_Mask = '\376';
   Hopeless_True_Mask = '\001';
   if  (Is_Hopeless_Scaff == NULL)
       Is_Hopeless_Scaff
           = (char *) safe_calloc (Num_Scaffolds, sizeof (char));
     else
       {
        Is_Hopeless_Scaff
            = (char *) safe_realloc (Is_Hopeless_Scaff,
                                     Num_Scaffolds * sizeof (char));
        if  (redo_index <= 0)
            for  (i = 0;  i < Num_Scaffolds;  i ++)
              Is_Hopeless_Scaff [i] &= Hopeless_False_Mask;
        for  (i = Is_Hopeless_Size;  i < Num_Scaffolds;  i ++)
          Is_Hopeless_Scaff [i] = '\0';
       }
   Is_Hopeless_Size = Num_Scaffolds;
#endif

   Filename_Prefix = prefix;

   fprintf (stderr, "### Fill_Gaps iteration #%d\n", iteration);
   
   sprintf (iter_string, "%d", iteration ++);

   now = time (NULL);
   fprintf (stderr, "### Start Rocks iteration #%s   %s\n",
            iter_string, ctime (& now));
   start_time = clock ();

   strcpy (filename, prefix);
   strcat (filename, ".rez.i");
   strcat (filename, iter_string);
   strcat (filename, ".log");
   log_file = file_open (filename, "w");

#if  MAKE_CAM_FILE
   strcpy (filename, prefix);
   strcat (filename, ".rez.i");
   strcat (filename, iter_string);
   strcat (filename, ".cam");
   Cam_File = file_open (filename, "w");

   for  (i = 0;  i < NUM_COLOURS;  i ++)
     fprintf (Cam_File, "%dREZ: %s\n", i, Colour_String [i]);
#if  SHOW_CALC_COORDS
   strcpy (filename, prefix);
   strcat (filename, ".calcrez.i");
   strcat (filename, iter_string);
   strcat (filename, ".cam");
   Calc_Cam_File = file_open (filename, "w");

   for  (i = 0;  i < NUM_COLOURS;  i ++)
     fprintf (Calc_Cam_File, "%dREZ: %s\n", i, Colour_String [i]);
#endif
#endif

PALLOC (Num_Scaffolds * sizeof (int64));
   Scaffold_Start = (int64 *) safe_calloc
                      (Num_Scaffolds, sizeof (int64));
PALLOC (Num_Scaffolds * sizeof (int64));
   Scaffold_End = (int64 *) safe_calloc
                    (Num_Scaffolds, sizeof (int64));
PALLOC (Num_Scaffolds * sizeof (char));
   Scaffold_Flipped = (char *) safe_calloc
                        (Num_Scaffolds, sizeof (char));

   for  (i = 0;  i < Num_Scaffolds;  i ++)
     Scaffold_Start [i] = Scaffold_End [i] = -1;

   Scaff_Join = CreateVA_Scaff_Join_t (INITIAL_SCAFF_JOIN_SIZE);

fprintf (stderr, ">>> Before  Print_Unique_Chunks\n");
   Print_Unique_Chunks (log_file);

fprintf (stderr, ">>> Before  Print_Scaffolds\n");
   Print_Scaffolds (log_file);

fprintf (stderr, ">>> Before  Print_Potential_Fill_Chunks\n");
   Print_Potential_Fill_Chunks (log_file, Maybe_Rock, FALSE);

   StartTimerT(&GlobalData->ChooseChunksTimer);

fprintf (stderr, ">>> Before  Scan_Gaps\n");
   fill_chunks = Scan_Gaps ();

fprintf (stderr, ">>> Before  Choose_Safe_Chunks\n");
   Choose_Safe_Chunks (fill_chunks, MIN_GOOD_LINKS, MIN_ROCK_COVER_STAT);

#if  VERBOSE
{
 fprintf (log_file, "\n>>> Fill after Choose_Safe_Chunks <<<\n");
 Print_Fill_Info (log_file, fill_chunks);
}
#endif

   Include_Good_Joins (fill_chunks);
fprintf (stderr, "Back from Include_Good_Joins\n");

#if  0
{
 fprintf (log_file, "\n>>> Fill after Include_Good_Joins <<<\n");
 Print_Fill_Info (log_file, fill_chunks);
}
#endif

   if (level <= 2)
     Check_Other_Links (fill_chunks);
   
   StopTimerT(&GlobalData->ChooseChunksTimer);

   if  (level == 5)
       {
fprintf (stderr, "Before  Identify_Best_Rocks\n");
        Identify_Best_Rocks (fill_chunks, FALSE);
       }

#if  0
{
 fprintf (log_file, "\n>>> Fill before  check_consistency <<<\n");
 Print_Fill_Info (log_file, fill_chunks);
}
#endif

   StartTimerT(&GlobalData->ConsistencyCheckTimer);
   {
    int  passed_consistency_check;

fprintf (stderr, "Before  Add_Gap_Ends\n");
    Add_Gap_Ends (fill_chunks);
fprintf (stderr, "Before  check_consistency\n");
    passed_consistency_check
        = check_consistency (fill_chunks, Num_Scaffolds, iteration - 1);
fprintf (stderr, "Before  Disqualify_Scaff_Chunks\n");
    Disqualify_Scaff_Chunks (fill_chunks);

    fprintf (stderr, "      Passed consistency check: %7d\n",
             passed_consistency_check);
   }
   StopTimerT(&GlobalData->ConsistencyCheckTimer);

#if  0
{
 fprintf (log_file, "\n>>> Fill after  check_consistency <<<\n");
 Print_Fill_Info (log_file, fill_chunks);
}
#endif

   if  (level == 4)
       {
fprintf (stderr, "Before  Identify_Best_Rocks\n");
        Identify_Best_Rocks (fill_chunks, TRUE);   // Original
//        Identify_Best_Rocks (fill_chunks, FALSE);
                    // Better until fix consistency check
       }

   if  (level <= 3)
       {
        // Do overlap path walking to confirm inserted chunks, same as for
        // stones.  Reject those not confirmed by walk.

fprintf (stderr, "Before Confirm_Stones\n");
#if  0
        Confirm_Stones (log_file, fill_chunks, FALSE);
#else
        Requalify_Scaff_Chunks (fill_chunks);
        New_Confirm_Stones (log_file, fill_chunks, FALSE);
        Disqualify_Scaff_Chunks (fill_chunks);
#endif
fprintf (stderr, "Before Jiggle_Positions\n");
        Jiggle_Positions (fill_chunks);
fprintf (stderr, "Before Adjust_By_Ref_Variance\n");
        Adjust_By_Ref_Variance (fill_chunks);
fprintf (stderr, "Before Sort_Insertions\n");
        Sort_Insertions (fill_chunks, By_Index);

fprintf (stderr, "Before Doublecheck_positions\n");
        Doublecheck_Positions (fill_chunks, TRUE);
       }
   else if  (level == 4 || level == 5)
       {
        // Do overlap path walking to confirm inserted chunks, same as for
        // stones.  Reject those not confirmed by walk.

fprintf (stderr, "Before Confirm_Stones\n");
#if  0
        Confirm_Stones (log_file, fill_chunks, FALSE);
#else
        Requalify_Scaff_Chunks (fill_chunks);
#if  VERBOSE
   fprintf (log_file, "\n\n>>>> Fill BEFORE New_Confirm_Stones <<<<\n");
   Print_Fill_Info (log_file, fill_chunks);
   fflush (log_file);
#endif

        Use_Partial_Stone_Paths = TRUE;
        New_Confirm_Stones (log_file, fill_chunks, TRUE);
        Use_Partial_Stone_Paths = FALSE;

#if  VERBOSE
   fprintf (log_file, "\n\n>>>> Fill AFTER New_Confirm_Stones <<<<\n");
   Print_Fill_Info (log_file, fill_chunks);
   fflush (log_file);
#endif
        Disqualify_Scaff_Chunks (fill_chunks);
#endif

fprintf (stderr, "Before Restore_Best_Rocks\n");
        Restore_Best_Rocks (fill_chunks);

#if  0
fprintf (stderr, "\n### Check_Loads after Restore_Best_Rocks\n");
Check_Loads ();
fprintf (stderr, "Flushing CacheSequenceDB\n");
ClearCacheSequenceDB(ScaffoldGraph->sequenceDB, TRUE);
ClearCacheSequenceDB(ScaffoldGraph->sequenceDB, FALSE);
Check_Loads ();
#endif

#if  VERBOSE
fprintf (log_file, "\n>>> Fill after  Restore_Best_Rocks <<<\n");
Print_Fill_Info (log_file, fill_chunks);
#endif
fprintf (stderr, "Before Jiggle_Positions\n");
        Jiggle_Positions (fill_chunks);
fprintf (stderr, "After Jiggle\n");
        Adjust_By_Ref_Variance (fill_chunks);
fprintf (stderr, "After Adjust_By_Ref_Variance\n");
        Sort_Insertions (fill_chunks, By_Index);
fprintf (stderr, "After Sort_Insertions\n");

#if  0
fprintf (stderr, "Before Doublecheck_Positions\n");
fprintf (log_file, "\n>>> Fill before  Doublecheck_Positions <<<\n");
Print_Fill_Info (log_file, fill_chunks);
#endif
        Doublecheck_Positions (fill_chunks, TRUE);
#if  0
fprintf (stderr, "After Doublecheck_Positions\n");
fprintf (log_file, "\n>>> Fill after  Doublecheck_Positions <<<\n");
Print_Fill_Info (log_file, fill_chunks);
#endif
       }
   else if  (level == 6)
       {
        Adjust_By_Ref_Variance (fill_chunks);
        Check_Rocks (log_file, fill_chunks);
       }
       

#if  MAKE_CAM_FILE
   Update_Colours (fill_chunks);
   Output_Cam_Files (fill_chunks);

   fclose (Cam_File);
#if  SHOW_CALC_COORDS
   fclose (Calc_Cam_File);
#endif
#endif

fprintf (stderr, "Set_Split_Flags\n");
   Set_Split_Flags (fill_chunks, ALL_FALSE);
fprintf (stderr, "After Set_Split_Flags\n");
fprintf (log_file, "\n>>> Fill before  Update_Scaffold_Graph <<<\n");
   Print_Fill_Info (log_file, fill_chunks);

   fclose (log_file);

   strcpy (filename, prefix);
   strcat (filename, ".rez.i");
   strcat (filename, iter_string);
   strcat (filename, ".analysis");
   log_file = file_open (filename, "w");
   Analyze_Rock_Fill (log_file, fill_chunks);
   fclose (log_file);


#if  CHECK_CELSIM_COORDS
   strcpy (filename, prefix);
   strcat (filename, ".rez.i");
   strcat (filename, iter_string);
   strcat (filename, ".place");
   log_file = file_open (filename, "w");
   Analyze_Placement (log_file, fill_chunks);
   fclose (log_file);
#endif

   StartTimerT(&GlobalData->UpdateTimer);
   if  (level > 1)
       {
#if  0
Clear_Keep_Flags (fill_chunks, 1);
#endif
#if  USE_MY_INSERT
        inserted = Insert_Chunks_In_Graph (ScaffoldGraph, fill_chunks, ROCKS);
#else
        inserted
            = Update_Scaffold_Graph
                  (ScaffoldGraph, fill_chunks, FALSE, FALSE, TRUE,
                   TRUE /* copyAllOverlaps...not used */, -1, ROCKS);
#endif
        fprintf (stderr, "             Actually inserted: %7d\n", inserted);
       }

#if  0     // Not ready yet
   Re_Check_Inserted_Rocks (fill_chunks, MIN_GOOD_LINKS);
#endif

   StopTimerT(&GlobalData->UpdateTimer);


#if  TEST_HOPELESS_SCAFFS
   Set_Is_Hopeless (fill_chunks);
#endif

   Free_Fill_Array (fill_chunks);
   Free_Global_Arrays ();

   UnJigglePositions();

   now = time (NULL);
   fprintf (stderr, "### Finish Rocks iteration #%s   %s\n",
            iter_string, ctime (& now));
   stop_time = clock ();
   fprintf (stderr, "### cpu time = %.1f sec\n",
               (double) (stop_time - start_time) / CLOCKS_PER_SEC);

   StopTimerT(&GlobalData->GapFillTimer);

   return inserted;     // <*** Change this back ***>
  }



static void  Fine_Tune_Positions
    (Scaffold_Fill_t * fill)

//  Check chunks in each gap of  fill  for overlaps and
//  adjust their positions where necessary.

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 1;  j < fill [scaff_id] . num_gaps - 1;  j ++)
        {
         int  k;
         Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                {
                 ;  // Fill in later
                }
           }
        }
     }

   return;
  }



void  Free_Fill_Array
    (Scaffold_Fill_t * fill_chunks)

//  Free all the memory in  fill_chunks .

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
        free (fill_chunks [scaff_id] . gap [j] . chunk);

      free (fill_chunks [scaff_id] . gap);
     }

   free (fill_chunks);

   return;
  }



void  Free_Fill_Array_Scaffold
    (Scaffold_Fill_t * fill_chunks)

//  Free all the memory in  fill_chunks as allocated by Scan_Gaps_In_Scaffold.

{
  int  j;
  
  for  (j = 0;  j < fill_chunks -> num_gaps;  j ++)
	free (fill_chunks -> gap [j] . chunk);
  
  free (fill_chunks -> gap);
  
  return;
}



int  Find_Olap_Path
    (ChunkInstanceT * from, int from_end, ChunkInstanceT * to,
     int num_targets, Target_Info_t target [], double bound,
     int * first, int * max_hits, int * max_first, LengthT * to_position,
     unsigned int edge_mask)

//  Look for a path of overlaps starting at vertex  from , stopping at
//  vertex  to , and trying to hit as many of the vertices in
//  target [0 .. (num_targets)]  as possible.  If  to  is  NULL ,
//  then the path need not terminate any particular place.
//  The path must issue from the  from_end  (A_END or B_END) of  from .
//  Do not consider any path longer (in DNA bases) than  bound .
//  Set  (* first)  to the subscript of the closest (measured in base-pair distance)
//  target to which a path can be found.
//  Set  (* max_hits)  to the most targets that can be hit on a path.
//  Set  (* max_first)  to the first (closest) target subscript on the
//  max-hit path.
//  Set  (* to_position)  to an estimate of the mean and variance
//  of the leading position of  to  based on the overlap path found.
//  edge_mask  indicates which edges can be used in the path.
//  Return  TRUE  if a successful path was found.  Successful
//  means hit any target nodes and reached  to .
//  If  (num_targets == 0)  a successful path just reaches  to .

  {
   static Node_Index_t  * node = NULL;
   static DFS_Info_t  * dfs_tree = NULL;
   static int  tree_size;
   double  max_first_total;
   int  max_first_dist;
   int  found, num_nodes, tree_sub;
   int  i;

   if  (to == NULL && num_targets == 0)
       return  TRUE;                        // Nothing to do

   num_nodes = GetNumGraphNodes (ScaffoldGraph -> RezGraph);
PRALLOC (num_nodes * sizeof (Node_Index_t));
   node = (Node_Index_t *) safe_realloc (node, num_nodes * sizeof (Node_Index_t));
   for  (i = 0;  i < num_nodes;  i ++)
     node [i] . visited = node [i] . finished = FALSE;

   if  (dfs_tree == NULL)
       {
        tree_size = INITIAL_TREE_SIZE;
PALLOC (tree_size * sizeof (DFS_Info_t));
        dfs_tree = (DFS_Info_t *) safe_malloc (tree_size * sizeof (DFS_Info_t));
       }
   tree_sub = 0;

   for  (i = 0;  i < num_targets;  i ++)
     {
      target [i] . found = 0;
      target [i] . next = -1;
      target [i] . where = -1.0;
      target [i] . total = 0.0;
     }

#if  0
if  (from -> id == 842)
    Global_Debug_Flag = TRUE;
#endif

if  (Global_Debug_Flag)
    {
     int  k;

     fprintf (stderr, "Targets\n");
     fprintf (stderr, "%3s %5s %5s %5s\n",
              "sub", "id", "lo", "hi");
     for  (k = 0;  k < num_targets;  k ++)
       fprintf (stderr, "%2d: %5d %5.0f %5.0f\n",
                k, target [k] . id, target [k] . lo,
                target [k] . hi);
    }

   found = Depth_First_Visit
       (from, from_end, to, num_targets, target, bound, 0,
        first, max_hits, max_first, & max_first_dist, & max_first_total,
        & tree_sub, & tree_size, & dfs_tree, node, 0, edge_mask);

if  (Global_Debug_Flag)
    {
     int  k;

     fprintf (stderr, "found = %s\n", found ? "True" : "False");
     fprintf (stderr, "%3s %5s %5s %5s %6s %6s %3s %3s\n",
              "sub", "id", "hits", "first", "dist", "total", "hit", "suc");
     for  (k = 0;  k < tree_sub;  k ++)
       fprintf (stderr, "%2d: %5d %5d %5d %6d %6.0f %3c %3c\n",
                k, dfs_tree [k] . id, dfs_tree [k] . max_hits,
                dfs_tree [k] . max_first,
                dfs_tree [k] . max_first_dist,
                dfs_tree [k] . max_first_total,
                dfs_tree [k] . is_hit ? 'T' : 'F',
                dfs_tree [k] . succeeded ? 'T' : 'F');
    }

Global_Debug_Flag = FALSE;

   if  ((* max_hits) > 0)
       {
        double  total = 0.0;
        int  place = 0;
        int  sub = 0;

        for  (i = 0;  i < (* max_hits);  i ++)
          {
           place += dfs_tree [sub] . max_first_dist;
           total += dfs_tree [sub] . max_first_total;
           sub = dfs_tree [sub] . max_first;
           target [dfs_tree [sub] . target_sub] . where = place;
           target [dfs_tree [sub] . target_sub] . total = total;
          }

        place += dfs_tree [sub] . max_first_dist;
        total += dfs_tree [sub] . max_first_total;
        to_position -> mean = place;
        to_position -> variance = ComputeFudgeVariance(total);

        (* first) = dfs_tree [(* first)] . target_sub;
        (* max_first) = dfs_tree [(* max_first)] . target_sub;
       }
   else if  (found)
      {
       to_position -> mean = dfs_tree [0] . max_first_dist;
       to_position -> variance = ComputeFudgeVariance(dfs_tree [0] . max_first_total);

      }

   return  (found && (num_targets == 0 || (* max_hits) > 0));
  }



static void  Fixup_Chunk_End_Variances
    (LengthT * left_end, LengthT * right_end, double diff)

//  Make difference in variances in  (* left_end)  and  (* right_end)
//  at least  diff .

  {
   double  delta;

   delta = diff - fabs (left_end -> variance - right_end -> variance);
   if  (delta > 0)
       {
        if  (left_end -> variance < right_end -> variance)
            right_end -> variance += delta;
          else
            left_end -> variance += delta;
       }

   return;
  }



void  Force_Increasing_Variances
    (void)

//  Make variances be non-decreasing in all scaffolds

  {
   GraphNodeIterator scaffolds;
   CIScaffoldT  * scaffold;

   InitGraphNodeIterator (& scaffolds, ScaffoldGraph -> ScaffoldGraph,
                          GRAPH_NODE_DEFAULT);
   while  ((scaffold = NextGraphNodeIterator (& scaffolds)) != NULL)
     {

       Force_Increasing_Variances_One_Scaffold (scaffold->id);
     }
   return;
  }



void  Force_Increasing_Variances_One_Scaffold
    (int scaff_id)

//  Make variances be non-decreasing in scaffold  scaff_id .

  {
   CIScaffoldT  * scaff;
   CIScaffoldTIterator  scaff_iterator;
   ChunkInstanceT  * chunk, * prev_chunk;
   double prev_variance, incr;

   scaff = GetGraphNode (ScaffoldGraph -> ScaffoldGraph, scaff_id);

   InitCIScaffoldTIterator (ScaffoldGraph, scaff, TRUE, FALSE,
                            & scaff_iterator);
   incr = prev_variance = 0.0;

   prev_chunk = NULL;

   while  ((chunk = NextCIScaffoldTIterator (& scaff_iterator)) != NULL)
     {
      double  dip, min_var;

      min_var = Min_double (chunk -> offsetAEnd . variance,
                            chunk -> offsetBEnd . variance);
      if  (min_var >= prev_variance)
          dip = 0.0;
        else
          {
           dip = prev_variance - min_var + EPSILON;

	   /*
           fprintf (stderr,
                    ">>> Forced variance up by %.1f at chunk %d in scaff %d\n",
                    dip, chunk -> id, scaff -> id);
	   */
          }

      prev_variance = Max_double (chunk -> offsetAEnd . variance,
                                  chunk -> offsetBEnd . variance);
      if  (prev_chunk != NULL)
          {
           prev_chunk -> offsetAEnd . variance += incr;
           prev_chunk -> offsetBEnd . variance += incr;
          }

      incr += dip;
      prev_chunk = chunk;
     }

   if  (prev_chunk != NULL)
       {
        prev_chunk -> offsetAEnd . variance += incr;
        prev_chunk -> offsetBEnd . variance += incr;
		scaff->bpLength.variance = MAX( prev_chunk->offsetAEnd.variance, 
										prev_chunk->offsetBEnd.variance);
       }

#if  1
   // Check to make sure it worked

   InitCIScaffoldTIterator
       (ScaffoldGraph, scaff, TRUE, FALSE, & scaff_iterator);
   prev_variance = 0.0;

   while  ((chunk = NextCIScaffoldTIterator (& scaff_iterator)) != NULL)
     {
      double  min_var;

      min_var = Min_double (chunk -> offsetAEnd . variance,
                            chunk -> offsetBEnd . variance);
      if  (min_var < prev_variance)
          {
           fprintf (stderr,
                    ">>> Force_Increasing failed by %.1f at chunk %d in scaff %d\n",
                    prev_variance - min_var, chunk -> id, scaff -> id);
          }

      prev_variance = Max_double (chunk -> offsetAEnd . variance,
                                  chunk -> offsetBEnd . variance);
     }
#endif

   return;
  }



static void   Free_Global_Arrays
    (void)

//  Free dynamically allocated global arrays.

  {
   int  cid;

   for  (cid = 0;  cid < Num_Chunks;  cid ++)
     if  (Chunk_Info [cid] . annotation != NULL)
         free (Chunk_Info [cid] . annotation);

   free (Chunk_Info);
   free (Ref_Data);
   free (Ref_Index);

   free (Scaffold_Start);
   free (Scaffold_End);
   free (Scaffold_Flipped);

   DeleteVA_Scaff_Join_t (Scaff_Join);

   return;
  }



static Overlap *  Get_Chunk_Overlap
    (Gap_Chunk_t * a, Gap_Chunk_t * b, char * * a_seq, char * * b_seq, FILE * fp)

//  Determine if contigs  a  and  b  might overlap, and if so
//  look for the overlap and return it if found.
//  Print log info to  fp .   (* a_seq)  and  (* b_seq)  contain the
//  sequence of contigs  a  and  b , or if they are NULL, then
//  the sequence will be loaded from the multialignment store..

  {
   double  max_var, min_var, var_diff, slop;
   ChunkOrientationType  orient;
   int  min_ahang, max_ahang;

   max_var = Max_double (a -> start . variance, a -> end . variance);
   max_var = Max_double (max_var, b -> start . variance);
   max_var = Max_double (max_var, b -> end . variance);

   min_var = Min_double (a -> start . variance, a -> end . variance);
   min_var = Min_double (min_var, b -> start . variance);
   min_var = Min_double (min_var, b -> end . variance);

   var_diff = max_var - min_var;

   if  (var_diff < 0.0)
       {
        fprintf (stderr, "ERROR:  Negative variance diff = %.1f\n", var_diff);
        slop = 0.0;
       }
     else
       slop = 5.0 * sqrt (var_diff);

#if  SHOW_OLAP_DETAILS
        fprintf (fp, "a_id = %d  b_id = %d  slop = %.0f",
                 a -> chunk_id, b -> chunk_id, slop);
#endif

   if  (Might_Overlap (a -> start . mean, a -> end . mean,
                       b -> start . mean, b -> end . mean,
                       slop, & orient, & min_ahang, & max_ahang))
       {
        Overlap  * result;

        if  ((* a_seq) == NULL)
            (* a_seq) = Get_Contig_Sequence (a -> chunk_id);
        if  ((* b_seq) == NULL)
            (* b_seq) = Get_Contig_Sequence (b -> chunk_id);
        result = OverlapSequences
                    ((* a_seq), (* b_seq), orient, min_ahang, max_ahang,
                     CGW_DP_ERATE, CGW_DP_THRESH, CGW_DP_MINLEN,
                     AS_FIND_ALIGN);
#if  SHOW_OLAP_DETAILS
        fprintf (fp, "  min_ahang = %4d  max_ahang = %4d  orient = %s\n",
                 min_ahang, max_ahang, Orientation_As_String (orient));
        if  (result == NULL)
            fprintf (fp, "   No overlap\n");
          else
            fprintf (fp, "   begpos = %d  endpos = %d  length = %d\n",
                     result -> begpos, result -> endpos, result -> length);
#endif
        return  result;
       }
     else
       {
#if  SHOW_OLAP_DETAILS
        fprintf (fp, "  Can't overlap\n");
#endif
       }

   return  NULL;
  }



static char *  Get_Contig_Sequence
    (int id)

//  Extract the ungapped sequence for the contig with  id
//  and return a pointer to it.  Allocates memory that must
//  be freed by the customer.

  {
   static MultiAlignT  * ma = NULL;
   char  * p, * gapped_seq, * ungapped_seq;
   int  ct, len;

   if  (ma == NULL)
       ma = CreateEmptyMultiAlignT ();

   ReLoadMultiAlignTFromSequenceDB
            (ScaffoldGraph -> sequenceDB, ma, id,
             ScaffoldGraph -> RezGraph -> type == CI_GRAPH);

   gapped_seq = Getchar (ma -> consensus, 0);
   len = strlen (gapped_seq);

   if  (len <= 0)
       {
        fprintf (stderr, "ERROR:  No multialign for contig %d\n", id);
        exit (EXIT_FAILURE);
       }

PALLOC (len + 1);
   ungapped_seq = (char *) safe_malloc (len + 1);

   ct = 0;
   for  (p = gapped_seq;  * p != '\0';  p ++)
     if  (isalpha (* p))
         ungapped_seq [ct ++] = * p;
   ungapped_seq [ct] = '\0';

//   free (gapped_seq);

   return  ungapped_seq;
  }



int  Hurl_Contained_Rocks
    (char * prefix, int level, int redo_index)

//  Find rocks that are contained within a contig already in a scaffold.
//  This is a main entry point from the chunk-graph-walker module.
//  prefix  is the prefix of the name to use for output files
//  level  indicates which steps to perform.
//  redo_index  is the number of previous calls to this routine without
//  renumbering or merging scaffolds

  {
   FILE  * log_file;
   char  filename [1000], iter_string [20];
   Scaffold_Fill_t  * fill_chunks;
   static int  iteration = 0;
   clock_t  start_time, stop_time;
   time_t  now;
   int  i;
   int  inserted = 0;
     
   StartTimerT(&GlobalData->GapFillTimer);

   Num_Scaffolds = GetNumGraphNodes (ScaffoldGraph -> ScaffoldGraph);
   if (Num_Scaffolds == 0)
     return 0;
     
#if  TEST_HOPELESS_SCAFFS
   Hopeless_False_Mask = '\375';
   Hopeless_True_Mask = '\002';
   if  (Is_Hopeless_Scaff == NULL)
       Is_Hopeless_Scaff
           = (char *) safe_calloc (Num_Scaffolds, sizeof (char));
     else
       {
        Is_Hopeless_Scaff
            = (char *) safe_realloc (Is_Hopeless_Scaff,
                                     Num_Scaffolds * sizeof (char));
        if  (redo_index <= 0)
            for  (i = 0;  i < Num_Scaffolds;  i ++)
              Is_Hopeless_Scaff [i] &= Hopeless_False_Mask;
        for  (i = Is_Hopeless_Size;  i < Num_Scaffolds;  i ++)
          Is_Hopeless_Scaff [i] = '\0';
       }
   Is_Hopeless_Size = Num_Scaffolds;
#endif

   Filename_Prefix = prefix;

   fprintf (stderr, "### Hurl_Contained_Rocks iteration #%d\n", iteration);
   Contained_Only_Switch = TRUE;
   
   sprintf (iter_string, "%d", iteration ++);

   now = time (NULL);
   fprintf (stderr, "### Start Contained Rocks iteration #%s   %s\n",
            iter_string, ctime (& now));
   start_time = clock ();

   strcpy (filename, prefix);
   strcat (filename, ".crocks.i");
   strcat (filename, iter_string);
   strcat (filename, ".log");
   log_file = file_open (filename, "w");

#if  MAKE_CAM_FILE
   strcpy (filename, prefix);
   strcat (filename, ".crocks.i");
   strcat (filename, iter_string);
   strcat (filename, ".cam");
   Cam_File = file_open (filename, "w");

   for  (i = 0;  i < NUM_COLOURS;  i ++)
     fprintf (Cam_File, "%dREZ: %s\n", i, Colour_String [i]);
#if  SHOW_CALC_COORDS
   strcpy (filename, prefix);
   strcat (filename, ".calccr.i");
   strcat (filename, iter_string);
   strcat (filename, ".cam");
   Calc_Cam_File = file_open (filename, "w");

   for  (i = 0;  i < NUM_COLOURS;  i ++)
     fprintf (Calc_Cam_File, "%dREZ: %s\n", i, Colour_String [i]);
#endif
#endif

PALLOC (Num_Scaffolds * sizeof (int64));
   Scaffold_Start = (int64 *) safe_calloc
                      (Num_Scaffolds, sizeof (int64));
PALLOC (Num_Scaffolds * sizeof (int64));
   Scaffold_End = (int64 *) safe_calloc
                    (Num_Scaffolds, sizeof (int64));
PALLOC (Num_Scaffolds * sizeof (char));
   Scaffold_Flipped = (char *) safe_calloc
                        (Num_Scaffolds, sizeof (char));

   for  (i = 0;  i < Num_Scaffolds;  i ++)
     Scaffold_Start [i] = Scaffold_End [i] = -1;

   Scaff_Join = CreateVA_Scaff_Join_t (INITIAL_SCAFF_JOIN_SIZE);
       // Not used, but put in to be compatible with regular rocks code

fprintf (stderr, ">>> Before  Print_Unique_Chunks\n");
   Print_Unique_Chunks (log_file);

fprintf (stderr, ">>> Before  Print_Scaffolds\n");
   Print_Scaffolds (log_file);

fprintf (stderr, ">>> Before  Print_Potential_Fill_Chunks\n");
   Print_Potential_Fill_Chunks (log_file, Maybe_Rock, TRUE);

   StartTimerT(&GlobalData->ChooseChunksTimer);

fprintf (stderr, ">>> Before  Scan_Gaps\n");
   fill_chunks = Scan_Gaps ();

fprintf (stderr, ">>> Before  Choose_Safe_Chunks\n");
   Choose_Safe_Chunks (fill_chunks, MIN_GOOD_LINKS, MIN_ROCK_COVER_STAT);

#if  VERBOSE
{
 fprintf (log_file, "\n>>> Fill after Choose_Safe_Chunks <<<\n");
 Print_Fill_Info (log_file, fill_chunks);
}
#endif

   if (level <= 2)
     Check_Other_Links (fill_chunks);
   
   StopTimerT(&GlobalData->ChooseChunksTimer);

   Add_Gap_Ends (fill_chunks);

#if  VERBOSE
{
 fprintf (log_file, "\n>>> Fill after Add_Gap_Ends <<<\n");
 Print_Fill_Info (log_file, fill_chunks);
}
#endif

   Confirm_Contained (log_file, fill_chunks, TRUE);

   Disqualify_Scaff_Chunks (fill_chunks);
   Sort_Insertions (fill_chunks, By_Keep_And_Low_Position);

#if  MAKE_CAM_FILE
   Update_Colours (fill_chunks);
   Output_Cam_Files (fill_chunks);

   fclose (Cam_File);
#if  SHOW_CALC_COORDS
   fclose (Calc_Cam_File);
#endif
#endif

fprintf (stderr, "Set_Split_Flags\n");
   Set_Split_Flags (fill_chunks, ALL_FALSE);
fprintf (stderr, "After Set_Split_Flags\n");
fprintf (log_file, "\n>>> Fill before  Update_Scaffold_Graph <<<\n");
   Print_Fill_Info (log_file, fill_chunks);

   fclose (log_file);

   strcpy (filename, prefix);
   strcat (filename, ".crocks.i");
   strcat (filename, iter_string);
   strcat (filename, ".analysis");
   log_file = file_open (filename, "w");
   Analyze_Rock_Fill (log_file, fill_chunks);
   fclose (log_file);


   StartTimerT(&GlobalData->UpdateTimer);
   if (level > 1) {
#if  0
Clear_Keep_Flags (fill_chunks, 1);
#endif
#if  USE_MY_INSERT
     inserted = Insert_Chunks_In_Graph (ScaffoldGraph, fill_chunks, ROCKS);
#else
     inserted = Update_Scaffold_Graph
       (ScaffoldGraph, fill_chunks, FALSE, FALSE, TRUE,
        TRUE /* copyAllOverlaps...not used */, -1, ROCKS);
#endif
     fprintf (stderr, "             Actually inserted: %7d\n", inserted);
   }
   StopTimerT(&GlobalData->UpdateTimer);

   Contained_Only_Switch = FALSE;

#if  TEST_HOPELESS_SCAFFS
   Set_Is_Hopeless (fill_chunks);
#endif

   Free_Fill_Array (fill_chunks);
   Free_Global_Arrays ();

   now = time (NULL);
   fprintf (stderr, "### Finish Contained Rocks iteration #%s   %s\n",
            iter_string, ctime (& now));
   stop_time = clock ();
   fprintf (stderr, "### cpu time = %.1f sec\n",
               (double) (stop_time - start_time) / CLOCKS_PER_SEC);

   StopTimerT(&GlobalData->GapFillTimer);

   return inserted;
  }



static void  Identify_Best_Rocks
    (Scaffold_Fill_t * fill_chunks, int check_keep)

//  Mark the rock in each gap of  fill_chunks  that has the smallest
//  variance if its position could place it in only one gap.
//  If  check_keep  is true, only consider chunks with  keep  flag
//  set;  otherwise, ignore the  keep  flag.

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = GAPS_TO_ADJUST)
        {
         int  k;
         double  min_variance, prev_gap_end, next_gap_begin;
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
         Gap_Chunk_t  * best_chunk = NULL;

         if  (j <= 0)
             prev_gap_end = -DBL_MAX;
           else
             prev_gap_end = fill_chunks [scaff_id] . gap [j - 1] . end . mean;
         if  (j >= fill_chunks [scaff_id] . num_gaps - 1)
             next_gap_begin = DBL_MAX;
           else
             next_gap_begin = fill_chunks [scaff_id] . gap [j + 1] . start . mean;


         // First identify chunks that could only go in this gap and
         // whose mean position indicates it extends into the gap sufficiently

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;
            double  left_extent, right_extent;
            LengthT  left, right;

            this_chunk -> best = FALSE;
            this_chunk -> candidate = FALSE;

            if  (this_chunk -> start . mean <= this_chunk -> end . mean)
                {
                 left = this_chunk -> start;
                 right = this_chunk -> end;
                }
              else
                {
                 left = this_chunk -> end;
                 right = this_chunk -> start;
                }
            left_extent = left . mean - 3.0 * sqrt (left . variance);
            right_extent = right . mean + 3.0 * sqrt (right . variance);

            if  ((! check_keep || this_chunk -> keep)
                   && Between (left_extent, right_extent,
                               prev_gap_end, next_gap_begin)
                   && right . mean > this_gap -> start . mean
                   && left . mean < this_gap -> end . mean
                   && this_chunk -> chunk_id != this_gap -> left_cid
                   && this_chunk -> chunk_id != this_gap -> right_cid)
                {
                 this_chunk -> candidate = TRUE;
                }
           }

         if  (j % 51 == 50 || this_gap -> num_chunks >= 50)
             CheckScaffoldGraphCache (ScaffoldGraph);

         // Now find the best chunk among the candidates.  Eliminate
         // any chunk that is contained in another candidate or
         // contained in the scaffold chunks on the ends of this gap.

         min_variance = DBL_MAX;
         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;
            int  m, is_contained;
            if  (! this_chunk -> candidate
                    || (this_gap -> left_cid >= 0
                        && Chunk_Contained_In_Scaff (this_chunk,
                                                     this_gap -> left_cid))
                    || (this_gap -> right_cid >= 0
                             && Chunk_Contained_In_Scaff (this_chunk,
                                                         this_gap -> right_cid)))
                continue;

            is_contained = FALSE;
            for  (m = 0;  m < this_gap -> num_chunks && ! is_contained;  m ++)
              {
               Gap_Chunk_t  * other_chunk = this_gap -> chunk + m;

               if  (m == k || ! other_chunk -> candidate)
                   continue;

               if  (Chunk_Contained_In_Chunk (this_chunk, other_chunk))
                   is_contained = TRUE;
              }

            if  (! is_contained
                     && this_chunk -> start . variance < min_variance)
                {
                 best_chunk = this_chunk;
                 min_variance = this_chunk -> start . variance;
                }
           }

         if  (best_chunk != NULL)
             best_chunk -> best = TRUE;
        }

      CheckScaffoldGraphCache (ScaffoldGraph);
     }

   return;
  }



static void  Include_Good_Joins
    (Scaffold_Fill_t * fill_chunks)

//  Check consistency of joins in global  Scaff_Join  and add
//  consistent ones to  fill_chunks .

  {
   Scaff_Join_t  * * p;
   int  i, j, n;

   n = GetNumVA_Scaff_Join_t (Scaff_Join);
   if  (n == 0)
       return;

PALLOC (n * sizeof (Scaff_Join_t *));
   p = (Scaff_Join_t * *) safe_calloc (n, sizeof (Scaff_Join_t *));

   for  (i = 0;  i < n;  i ++)
     p [i] = GetVA_Scaff_Join_t (Scaff_Join, i);

   qsort (p, n, sizeof (Scaff_Join_t *), Scaff_Join_Cmp);

   for  (i = 0;  i < n;  i = j + 1)
     {
      CIScaffoldT  * scaffold;
      ChunkInstanceT  * chunk;
      double  tail, lo1, hi1, lo2, hi2;
      int  r, s;
      int  is_bad;

      for  (j = i + 1;  j < n && Scaff_Join_Cmp (p + i, p + j) == 0;  j ++)
        ;
      j --;

      is_bad = FALSE;
      for  (r = 0;  r < i;  r ++)
        if  (p [r] -> scaff1 == p [i] -> scaff1)
            {
             scaffold = GetCIScaffoldT (ScaffoldGraph -> CIScaffolds, p [r] -> scaff2);
             chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                        scaffold -> info . Scaffold . BEndCI);
             tail = Max_double (chunk -> offsetAEnd . mean, chunk -> offsetBEnd . mean);
             if  (p [r] -> m > 0)
                 {
                  lo1 = (0.0 - p [r] -> b) / p [r] -> m;
                  hi1 = (tail - p [r] -> b) / p [r] -> m;
                 }
               else
                 {
                  hi1 = (0.0 - p [r] -> b) / p [r] -> m;
                  lo1 = (tail - p [r] -> b) / p [r] -> m;
                 }
             scaffold = GetCIScaffoldT (ScaffoldGraph -> CIScaffolds, p [i] -> scaff2);
             chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                        scaffold -> info . Scaffold . BEndCI);
             tail = Max_double (chunk -> offsetAEnd . mean, chunk -> offsetBEnd . mean);
             if  (p [i] -> m > 0)
                 {
                  lo2 = (0.0 - p [i] -> b) / p [i] -> m;
                  hi2 = (tail - p [i] -> b) / p [i] -> m;
                 }
               else
                 {
                  hi2 = (0.0 - p [i] -> b) / p [i] -> m;
                  lo2 = (tail - p [i] -> b) / p [i] -> m;
                 }
             if  (Interval_Intersection (lo1, hi1, lo2, hi2))
                 {
                  p [r] -> is_bad = TRUE;
                  is_bad = TRUE;
                 }
            }
        else if  (p [r] -> scaff2 == p [i] -> scaff2)
            {
             scaffold = GetCIScaffoldT (ScaffoldGraph -> CIScaffolds, p [r] -> scaff1);
             chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                        scaffold -> info . Scaffold . BEndCI);
             tail = Max_double (chunk -> offsetAEnd . mean, chunk -> offsetBEnd . mean);
             if  (p [r] -> m > 0)
                 {
                  lo1 = 0.0 + p [r] -> b;
                  hi1 = p [r] -> m * tail + p [r] -> b;
                 }
               else
                 {
                  hi1 = 0.0 + p [r] -> b;
                  lo1 = p [r] -> m * tail + p [r] -> b;
                 }
             scaffold = GetCIScaffoldT (ScaffoldGraph -> CIScaffolds, p [i] -> scaff1);
             chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                        scaffold -> info . Scaffold . BEndCI);
             tail = Max_double (chunk -> offsetAEnd . mean, chunk -> offsetBEnd . mean);
             if  (p [i] -> m > 0)
                 {
                  lo2 = 0.0 + p [i] -> b;
                  hi2 = p [i] -> m * tail + p [i] -> b;
                 }
               else
                 {
                  hi2 = 0.0 + p [i] -> b;
                  lo2 = p [i] -> m * tail + p [i] -> b;
                 }
             if  (Interval_Intersection (lo1, hi1, lo2, hi2))
                 {
                  p [r] -> is_bad = TRUE;
                  is_bad = TRUE;
                 }
            }
        else if  (p [r] -> scaff1 == p [i] -> scaff2)
            {
             scaffold = GetCIScaffoldT (ScaffoldGraph -> CIScaffolds, p [r] -> scaff2);
             chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                        scaffold -> info . Scaffold . BEndCI);
             tail = Max_double (chunk -> offsetAEnd . mean, chunk -> offsetBEnd . mean);
             if  (p [r] -> m > 0)
                 {
                  lo1 = (0.0 - p [r] -> b) / p [r] -> m;
                  hi1 = (tail - p [r] -> b) / p [r] -> m;
                 }
               else
                 {
                  hi1 = (0.0 - p [r] -> b) / p [r] -> m;
                  lo1 = (tail - p [r] -> b) / p [r] -> m;
                 }
             scaffold = GetCIScaffoldT (ScaffoldGraph -> CIScaffolds, p [i] -> scaff1);
             chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                        scaffold -> info . Scaffold . BEndCI);
             tail = Max_double (chunk -> offsetAEnd . mean, chunk -> offsetBEnd . mean);
             if  (p [i] -> m > 0)
                 {
                  lo2 = 0.0 + p [i] -> b;
                  hi2 = p [i] -> m * tail + p [i] -> b;
                 }
               else
                 {
                  hi2 = 0.0 + p [i] -> b;
                  lo2 = p [i] -> m * tail + p [i] -> b;
                 }
             if  (Interval_Intersection (lo1, hi1, lo2, hi2))
                 {
                  p [r] -> is_bad = TRUE;
                  is_bad = TRUE;
                 }
            }
        else if  (p [r] -> scaff2 == p [i] -> scaff1)
            {
             scaffold = GetCIScaffoldT (ScaffoldGraph -> CIScaffolds, p [r] -> scaff1);
             chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                        scaffold -> info . Scaffold . BEndCI);
             tail = Max_double (chunk -> offsetAEnd . mean, chunk -> offsetBEnd . mean);
             if  (p [r] -> m > 0)
                 {
                  lo1 = 0.0 + p [r] -> b;
                  hi1 = p [r] -> m * tail + p [r] -> b;
                 }
               else
                 {
                  hi1 = 0.0 + p [r] -> b;
                  lo1 = p [r] -> m * tail + p [r] -> b;
                 }
             scaffold = GetCIScaffoldT (ScaffoldGraph -> CIScaffolds, p [i] -> scaff2);
             chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                        scaffold -> info . Scaffold . BEndCI);
             tail = Max_double (chunk -> offsetAEnd . mean, chunk -> offsetBEnd . mean);
             if  (p [i] -> m > 0)
                 {
                  lo2 = (0.0 - p [i] -> b) / p [i] -> m;
                  hi2 = (tail - p [i] -> b) / p [i] -> m;
                 }
               else
                 {
                  hi2 = (0.0 - p [i] -> b) / p [i] -> m;
                  lo2 = (tail - p [i] -> b) / p [i] -> m;
                 }
             if  (Interval_Intersection (lo1, hi1, lo2, hi2))
                 {
                  p [r] -> is_bad = TRUE;
                  is_bad = TRUE;
                 }
            }

      for  (r = i;  r < j && ! is_bad;  r ++)
        {
         double  r_del = 3 * sqrt (p [r] -> variance);

         for  (s = i + 1;  s <= j && ! is_bad;  s ++)
           {
            double  s_del = 3 * sqrt (p [s] -> variance);

            is_bad = (! Interval_Intersection
                        ((int) (p [r] -> b - r_del),
                         (int) (p [r] -> b + r_del),
                         (int) (p [s] -> b - s_del),
                         (int) (p [s] -> b + s_del)));
           }
        }

      if  (! is_bad)
          {
           p [i] -> violated = Violates_Scaff_Edges (p [i]);
           if  (p [i] -> violated)
               is_bad = TRUE;
          }

      if  (is_bad)
          for  (r = i;  r <= j;  r ++)
            p [r] -> is_bad = TRUE;
      else if  (! ALLOW_LOOSE_END_ROCKS)
          {
           int  best = -1;

           for  (r = i;  r <= j;  r ++)
             if  (! p [r] -> is_bad
                    && (best == -1
                          || p [best] -> link_ct < p [r] -> link_ct))
                 best = r;
           for  (r = i;  r <= j;  r ++)
             if  (r != best)
                 p [r] -> is_bad = TRUE;
          }
     }

   fprintf (stderr, "%6s %6s %5s %10s %7s %5s %5s %7s\n",
            "scaff1", "scaff2", "m", "b", "stdev", "bad", "vio", "chunk");
   for  (i = 0;  i < n;  i ++)
     {
      fprintf (stderr, "%6d %6d %5d %10.1f %7.1f %5s %5s %7d\n",
               p [i] -> scaff1, p [i] -> scaff2, p [i] -> m, p [i] -> b,
               sqrt (p [i] -> variance), p [i] -> is_bad ? "y" : "n",
               p [i] -> violated ? "y" : "n",
               p [i] -> chunk_id);
     }

   for  (i = 0;  i < n;  i ++)
     if  (! p [i] -> is_bad)
         {
          ChunkInstanceT  * chunk;
          int  assign_succeeded;
#if  CHECK_CELSIM_COORDS
          double  this_offset, diff, cutoff;
          int  cam_colour;
#endif
#if  MAKE_CAM_FILE && SHOW_CALC_COORDS
          char  annotation_string [MAX_STRING_LEN];
#if SHOW_CALC_COORDS
          int left_coord, right_coord;
#endif
#endif

          chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                     p [i] -> chunk_id);
          assign_succeeded
              = Assign_To_Gap (p [i] -> chunk_id, p [i] -> left_end,
                               p [i] -> right_end, p [i] -> gap,
                               p [i] -> insert_scaff, p [i] -> flipped,
                               fill_chunks, p [i] -> edge_quality,
                               p [i] -> cover_stat, p [i] -> link_ct,
                               JOINER_ROCK_CHAR);
#if  MAKE_CAM_FILE
          if  (assign_succeeded)
              Chunk_Info [p [i] -> chunk_id] . colour = PLACED_COLOUR;
          sprintf (annotation_string,
"  joins scaffs  low = <%.0f, %.0f>  high = <%.0f, %.0f>  cov = %d  typ = %s",
                   p [i] -> left_end . mean, sqrt (p [i] -> left_end . variance),
                   p [i] -> right_end . mean, sqrt (p [i] -> right_end . variance),
                   p [i] -> cover_stat,
                   CGB_Type_As_String (chunk -> flags . bits . cgbType));
          if  (Chunk_Info [p [i] -> chunk_id] . annotation != NULL)
              free (Chunk_Info [p [i] -> chunk_id] . annotation);
          Chunk_Info [p [i] -> chunk_id] . annotation
              = strdup (annotation_string);
#endif
#if  CHECK_CELSIM_COORDS
          if  (Scaffold_Flipped [p [i] -> insert_scaff])
              this_offset
                  = Max_int (chunk -> aEndCoord, chunk -> bEndCoord)
                        + p [i] -> left_end . mean;
            else
              this_offset
                  = Min_int (chunk -> aEndCoord, chunk -> bEndCoord)
                        - p [i] -> left_end . mean;
          diff = fabs (this_offset - p [i] -> stack_val . celsim_offset);
          cutoff = Max_double (6 * sqrt (p [i] -> stack_val . edge
                                        -> distance . variance
                                        + p [i] -> stack_val . source_variance),
                            100);

          if  (chunk -> flags . bits . cgbType != UU_CGBTYPE)
              fprintf (stderr,
                       "### Placed confused join chunk #%d  type = %s\n",
                       p [i] -> chunk_id,
                       CGB_Type_As_String
                           (chunk -> flags . bits . cgbType));
          else if  (diff > cutoff )
              {
               fprintf (stderr, "### Misplaced join chunk #%d  diff = %.0f\n",
                        p [i] -> chunk_id, diff);
               cam_colour = MISPLACED_COLOUR;
              }
#endif
#if  MAKE_CAM_FILE && SHOW_CALC_COORDS
          if  (Scaffold_Start [p [i] -> insert_scaff]
                 < Scaffold_End [p [i] -> insert_scaff])
              {
               left_coord = p [i] -> left_end . mean
                              + Scaffold_Start [p [i] -> insert_scaff];
               right_coord = p [i] -> right_end . mean
                               + Scaffold_Start [p [i] -> insert_scaff];
              }
            else
              {
               left_coord = Scaffold_Start [p [i] -> insert_scaff]
                              - p [i] -> right_end . mean;
               right_coord = Scaffold_Start [p [i] -> insert_scaff]
                               - p [i] -> left_end . mean;
              }
          if  (left_coord < 0)
              left_coord = 0;
          if  (right_coord < 0)
              right_coord = 0;
          Chunk_Info [p [i] -> chunk_id] . scaff_id = p [i] -> insert_scaff;
          Chunk_Info [p [i] -> chunk_id] . calc_left = left_coord;
          Chunk_Info [p [i] -> chunk_id] . calc_right = right_coord;
#endif
         }

   free (p);

   return;
  }


// return TRUE iff (* chunk) has 0 instances of unitig

int IsSurrogate(ChunkInstanceT * chunk)
{
  if(chunk->flags.bits.isContig)
    return(chunk->info.Contig.numCI == 1 &&
           chunk->info.Contig.AEndCI != NULLINDEX &&
           (GetGraphNode(ScaffoldGraph->CIGraph,
               chunk->info.Contig.AEndCI))->info.CI.numInstances != 0);
  else
    return(chunk->info.CI.numInstances != 0);
  
}

int  Is_Unique
    (ChunkInstanceT * chunk)

//  Return  TRUE  iff  (* chunk)  has been resolved as not being
//  repetitive.

 {
   if  (chunk -> flags . bits. isContig == 1)
       return chunk->scaffoldID > NULLINDEX;

  return  (chunk -> type == DISCRIMINATORUNIQUECHUNK_CGW
	      || chunk -> type == UNIQUECHUNK_CGW);
 }



static int  Insert_Chunks_In_Graph
    (ScaffoldGraphT * sgraph, Scaffold_Fill_t * fill, Kind_Of_Fill_t kind)

//  Insert all entries marked  keep  in  fill  into  sgraph .
//  kind  indicates type of entries.
//  Return the number of entries inserted.

  {
   int  inserted = 0;
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     inserted += Insert_Chunks_In_Graph_One_Scaffold (sgraph, fill, scaff_id, kind);

   return  inserted;
  }



static int  Insert_Chunks_In_Graph_One_Scaffold
    (ScaffoldGraphT * sgraph, Scaffold_Fill_t * fill, int scaff_id,
     Kind_Of_Fill_t kind)

//  Insert all entries marked  keep  in  fill [scaff_id]  into  sgraph .
//  Use  split  flag  to determine whether to  split off a surrogate or
//  not.   kind  indicates type of entries.
//  Return the number of entries inserted.

  {
   int  inserted = 0;
   int  j;

#if  TEST_HOPELESS_SCAFFS
   if  (Is_Hopeless_Scaff [scaff_id] & Hopeless_True_Mask)
       {
        fprintf (stderr, "HOPELESS  scaff = %d\n", scaff_id);
        return  0;
       }
#endif

   for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
     {
      Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
      int  k;

      if  (j > 0)
          REF (this_gap -> left_cid) . is_unthrowable = TRUE;

      for  (k = 0;  k < this_gap -> num_chunks;  k ++)
        {
         Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

         if  (this_chunk -> keep)
             {
              ContigT  * contig;
              int  insert_id;

              Remove_From_Scaffold (this_chunk);

              contig = GetGraphNode (ScaffoldGraph -> RezGraph,
                                     this_chunk -> chunk_id);
              if  (contig -> info . Contig . numCI == 1)
                  {
                   ChunkInstanceT  * unitig;

                   assert (contig -> info . Contig . AEndCI
                             == contig -> info . Contig . BEndCI);
                   unitig
                       = GetGraphNode (ScaffoldGraph -> CIGraph,
                                       contig -> info . Contig . AEndCI);
                   switch  (kind)
                     {
                      case  ROCKS :
                        unitig -> flags . bits . isRock = TRUE;
                        break;
                      case  STONES :
                        unitig -> flags . bits . isStone = TRUE;
                        break;
                      case  WALKS :
                        unitig -> flags . bits . isWalk = TRUE;
                        break;
                      default :
                        fprintf (stderr, "ERROR:  Unexpected insert type = %d\n",
                                 (int) kind);
                     }
                  }
                  

              if  (this_chunk -> split)
                  insert_id = SplitUnresolvedContig
                                  (ScaffoldGraph -> RezGraph,
                                   this_chunk -> chunk_id, NULL, TRUE);
                else
                  insert_id = this_chunk -> chunk_id;

              contig = GetGraphNode (ScaffoldGraph -> RezGraph, insert_id);

              switch  (kind)
                {
                 case  ROCKS :
                   contig -> flags . bits . isRock = TRUE;
                   break;
                 case  STONES :
                   contig -> flags . bits . isStone = TRUE;
                   break;
                 case  WALKS :
                   contig -> flags . bits . isWalk = TRUE;
                   break;
                 default :
                   fprintf (stderr, "ERROR:  Unexpected insert type = %d\n",
                            (int) kind);
                }

              InsertCIInScaffold (sgraph, insert_id, scaff_id,
                                  this_chunk -> start, this_chunk -> end,
                                  TRUE, NO_CONTIGGING);
              inserted ++;
             }
        }
     }

   return  inserted;
  }



static int  Is_Good_Scaff_Edge
    (SEdgeT * edge)

//  Return whether scaffold edge  edge  should be used.

  {
   if  (edge -> flags . bits . isProbablyBogus
          || edge -> flags . bits . isMarkedForDeletion
          || (edge -> edgesContributing == 1
                && edge -> distance . mean < -10.000))
       return  FALSE;

   return  TRUE;
  }



/*
  this function is needed to undo Jiggle_Positions() changes
  to singleton scaffolds where the contig offset may not be 0
  imd 07/11/02
*/
static void UnJigglePositions(void)
{
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;

  // iterate over scaffolds
  InitGraphNodeIterator(&scaffolds,
                        ScaffoldGraph->ScaffoldGraph,
                        GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL)
  {
    // work with singleton, real, live scaffolds
    if(scaffold->type == REAL_SCAFFOLD &&
       !scaffold->flags.bits.isDead)
    {
      NodeCGW_T * ci;
      ci = GetGraphNode(ScaffoldGraph->RezGraph,
                        scaffold->info.Scaffold.AEndCI);

      // work on scaffolds where cis are off
      if(MIN(ci->offsetAEnd.mean, ci->offsetBEnd.mean) != 0.0)
      {
        LengthT delta;

        fprintf(stderr, "\n\nSQUAWK!\n"
                "scaffold %d's first contig (%d) doesn't start at offset 0.0!\n",
                scaffold->id, ci->id);
        fprintf(stderr, "\toffsets: (%.3f,%.3f) - (%.3f,%.3f)\n",
                ci->offsetAEnd.mean, ci->offsetAEnd.variance,
                ci->offsetBEnd.mean, ci->offsetBEnd.variance);
        fprintf(stderr,
                "\tlengths: scaffold (%.3f,%.3f), contig (%.3f, %.3f)\n",
                scaffold->bpLength.mean, scaffold->bpLength.variance,
                ci->bpLength.mean, ci->bpLength.variance);
        
        fprintf(stderr, "*** BEFORE UNJIGGLING POSITIONS:\n");
        DumpCIScaffold(stderr, ScaffoldGraph, scaffold, 0);
        
        if(ci->offsetAEnd.mean < ci->offsetBEnd.mean)
        {
          delta.mean = - ci->offsetAEnd.mean;
          delta.variance = - fabs(ci->offsetAEnd.variance);
        }
        else
        {
          delta.mean = - ci->offsetBEnd.mean;
          delta.variance = - fabs(ci->offsetBEnd.variance);
        }
        AddDeltaToScaffoldOffsets(ScaffoldGraph, scaffold->id,
                                  ci->id, TRUE, FALSE, delta);
        fprintf(stderr, "*** AFTER UNJIGGLING POSITIONS:\n");
        DumpCIScaffold(stderr, ScaffoldGraph, scaffold, 0);

      }
    }
  }
}


static void  Jiggle_Positions
    (Scaffold_Fill_t * fill_chunks)

//  Apply the gap adjustments in  fill_chunks  to scaffolds
//  and to the chunks to be inserted in them.
//  This moves positions of scaffold elements to make them
//  more compatible with subsequent insertions.

  {
   int  scaff_id;

   fprintf (stderr, "### Jiggle_Positions ###\n");

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     Jiggle_Positions_One_Scaffold
         (fill_chunks, scaff_id);

   return;
  }



static void  Jiggle_Positions_One_Scaffold
    (Scaffold_Fill_t * fill_chunks, int scaff_id)

//  Apply the gap adjustments in  fill_chunks [scaff_id]  to
//  that scaffold and to the chunks to be inserted in it.
//  This moves positions of scaffold elements to make them
//  more compatible with subsequent insertions.

  {
   LengthT  cum_adjust;
   int  j;

   cum_adjust . mean = cum_adjust . variance = 0.0;

#ifdef DEBUG_DETAILED
fprintf (stderr, "### Scaffold %d:\n", scaff_id);
#endif

   for  (j = GAPS_TO_ADJUST)
     {
      Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;

      if  (j < fill_chunks [scaff_id] . num_gaps - 1
             && (this_gap -> adjustment . mean != 0.0
                   || this_gap -> adjustment . variance != 0.0))
          AddDeltaToScaffoldOffsets
              (ScaffoldGraph, scaff_id,
               this_gap -> right_cid, TRUE, FALSE,
               this_gap -> adjustment);

#ifdef DEBUG_DETAILED
fprintf (stderr, "### Gap %d  gapadjv = %.1f  cumadjv = %.1f  oldrefv = %.1f",
      j, this_gap -> adjustment . variance, cum_adjust . variance,
      this_gap -> ref_variance);
#endif

      this_gap -> ref_variance += cum_adjust . variance;

#ifdef DEBUG_DETAILED
fprintf (stderr, "  newrefv = %.1f\n",
      this_gap -> ref_variance);
#endif

      if  (j > 0)
          {
           int  k;

           for  (k = 0;  k < this_gap -> num_chunks;  k ++)
             {
              Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

//                 if  (this_chunk -> keep)
                  {
                   this_chunk -> start . mean += cum_adjust . mean;
                   this_chunk -> end . mean += cum_adjust . mean;
                  }
             }
           this_gap -> start . mean += cum_adjust . mean;
           this_gap -> start . variance += cum_adjust . variance;
          }

      cum_adjust . mean += this_gap -> adjustment . mean;
      cum_adjust . variance += this_gap -> adjustment . variance;
      if  (j < fill_chunks [scaff_id] . num_gaps - 1)
          {
           this_gap -> end . mean += cum_adjust . mean;
           this_gap -> end . variance += cum_adjust . variance;
          }
     }

   return;
  }



static int  Just_True
    (ContigT * chunk)

//  Always return  TRUE;

  {
   return  TRUE;
  }



static void  Kill_Duplicate_Stones
    (Scaffold_Fill_t * fill)

//  Set  keep  flag to false of duplicate stones in the same gap
//  anywhere in  fill  if at (essentially) the same position.

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     Kill_Duplicate_Stones_One_Scaffold (fill, scaff_id);

   return;
  }



#define  TOLERANCE    30

static void  Kill_Duplicate_Stones_One_Scaffold
    (Scaffold_Fill_t * fill, int scaff_id)

//  Set  keep  flag to false of duplicate stones in the same gap of
//  fill [scaff_id]  at (essentially) the same position.

  {
   int  j;

   for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
     {
      Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
      int  k;

      for  (k = 0;  k < this_gap -> num_chunks;  k ++)
        {
         Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

         if  (this_chunk -> keep)
             {
              int  m;

              for  (m = 0;  m < k;  m ++)
                {
                 Gap_Chunk_t  * prev_chunk = this_gap -> chunk + m;

                 if  (prev_chunk -> chunk_id == this_chunk -> chunk_id
                        && fabs (prev_chunk -> start . mean
                                   - this_chunk -> start . mean) <= TOLERANCE
                        && fabs (prev_chunk -> end . mean
                                   - this_chunk -> end . mean) <= TOLERANCE)
                     {
                      this_chunk -> keep = FALSE;
                      break;
                     }
                }
             }
        }
     }

   return;
  }



static int  Maybe_Rock
    (ContigT * chunk)

//  Return  TRUE  iff  chunk  has a chance of being a
//  rock.

  {
   return  chunk -> flags . bits . isPotentialRock
             && chunk -> info . CI . numInstances == 0;
  }



static int  Maybe_Stone
    (ContigT * chunk)

//  Return  TRUE  iff  chunk  has a chance of being a
//  rock.

  {
   return  chunk -> flags . bits . isPotentialStone;
  }



static int  Might_Overlap
    (double a_frag_start, double a_frag_end,
     double b_frag_start, double b_frag_end,
     double slop, ChunkOrientationType * orient,
     int * min_ahang, int * max_ahang)

//  Return  TRUE  iff the sequence from  a_frag_start  to
//  a_frag_end  could overlap the sequence from  b_frag_start
//  to  b_frag_end  if we allow one of the sequences to shift
//  left or right by  slop  positions.  Set (* orient)  to
//  the orientation of the sequences and  (* min_ahang)  and
//  (* max_ahang)  to the range of distances the A fragment
//  could extend to the left of the B fragment in the possible
//  overlaps (including negative values).

  {
   double  a_lo, a_hi, b_lo, b_hi;
   double  a_len, b_len;

   a_lo = Min_double (a_frag_start, a_frag_end);
   a_hi = Max_double (a_frag_start, a_frag_end);
   b_lo = Min_double (b_frag_start, b_frag_end);
   b_hi = Max_double (b_frag_start, b_frag_end);
   a_len = a_hi - a_lo;
   b_len = b_hi - b_lo;

   if  (a_frag_start < a_frag_end)
       {
        if  (b_frag_start < b_frag_end)
            (* orient) = AB_AB;
          else
            (* orient) = AB_BA;
       }
     else
       {
        if  (b_frag_start < b_frag_end)
            (* orient) = BA_AB;
          else
            (* orient) = BA_BA;
       }

//  If  slop  is negative, the intersection represents a "must overlap"
//  condition, but  min_ahang  and  max_ahang  are still set for
//  "might overlap" to be used in a call to  OverlapSequences, for
//  instance.
   if  (Interval_Intersection (a_lo - slop, a_hi + slop,
                               b_lo, b_hi))
       {
        (* max_ahang) = (int) rint (Min_double (a_len, b_lo - a_lo + abs (slop)));
        (* min_ahang) = (int) rint (Max_double (- b_len, b_lo - a_lo - abs (slop)));
        return  TRUE;
       }

   return  FALSE;
  }



static void  New_Confirm_Stones
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int use_all)

//  Try to confirm stones in  fill_chunks  by checking overlaps
//  among them.  Send output log to  fp .
//  If  use_all  is true, use all the chunks in  fill_chunks ;
//  otherwise, use only the ones whose  keep  flag is already
//  true.

  {
   int  scaff_id;

   fprintf (fp, "\n New_Confirm_Stones:\n");

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     New_Confirm_Stones_One_Scaffold
         (fp, fill_chunks, use_all, scaff_id);

   return;
  }



static void  New_Confirm_Stones_One_Scaffold
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int use_all, int scaff_id)

//  Try to confirm stones in scaffold  scaff_id  of  fill_chunks  by
//  checking overlaps among them.  Send output log to  fp .
//  If  use_all  is true, use all the chunks in  fill_chunks ;
//  otherwise, use only the ones whose  keep  flag is already
//  true.

  {
   static Gap_Chunk_t  * * check = NULL;
   static Stone_Edge_t  * edge_pool = NULL;
   static int  * forward_edge = NULL;
   static int  * reverse_edge = NULL;
   static int  * sorted = NULL;
   static char  * * sequence = NULL;
   static int  check_size = 0;
   static int  edge_pool_size = 0;
   int  num_sorted;
   LengthT  target_position;
   int  j;

   if  (check_size == 0 || edge_pool_size == 0)
       {
        check_size = INITIAL_GAP_ENTRIES;
        edge_pool_size = INITIAL_GAP_ENTRIES;
PALLOC (check_size * sizeof (Gap_Chunk_t *));
        check = (Gap_Chunk_t * *)
                     safe_malloc (check_size * sizeof (Gap_Chunk_t *));
PALLOC (edge_pool_size * sizeof (Stone_Edge_t));
        edge_pool = (Stone_Edge_t *)
                     safe_malloc (edge_pool_size * sizeof (Stone_Edge_t));
PALLOC (check_size * sizeof (int));
        forward_edge = (int *)
                     safe_malloc (check_size * sizeof (int));
PALLOC (check_size * sizeof (int));
        reverse_edge = (int *)
                     safe_malloc (check_size * sizeof (int));
PALLOC (check_size * sizeof (int));
        sorted = (int *)
                     safe_malloc (check_size * sizeof (int));
PALLOC (check_size * sizeof (char *));
        sequence = (char * *)
                     safe_malloc (check_size * sizeof (char *));
       }

   for  (j = GAPS_TO_ADJUST)
     {
      //  Find all overlaps
      //  Find all edges
      //  Get score for each stone
      //  Sort by scores
      //  Put in stones if sufficient score and consistent

      int  k, p, q, ct;
      int  next_edge, num_kept;
      Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;

#if  VERBOSE
fprintf (fp, "\nNew_Confirm Scaff %d  Gap %d\n", scaff_id, j);
fprintf (stderr, "\nScaff %d  Gap %d\n", scaff_id, j);
#endif

      ct = 0;
      if  (this_gap -> num_chunks > check_size)
          {
           check_size *= 2;
           if  (check_size < this_gap -> num_chunks)
               check_size = this_gap -> num_chunks;
PRALLOC (check_size * sizeof (Gap_Chunk_t *));
           check = (Gap_Chunk_t * *)
                       safe_realloc
                           (check,
                            check_size * sizeof (Gap_Chunk_t *));
PRALLOC (check_size * sizeof (int));
           forward_edge = (int *)
                              safe_realloc
                                  (forward_edge,
                                   check_size * sizeof (int));
PRALLOC (check_size * sizeof (int));
           reverse_edge = (int *)
                              safe_realloc
                                  (reverse_edge,
                                   check_size * sizeof (int));
PRALLOC (check_size * sizeof (int));
           sorted = (int *)
                        safe_realloc
                            (sorted, check_size * sizeof (int));
PRALLOC (check_size * sizeof (char *));
           sequence = (char * *)
                          safe_realloc (sequence,
                                        check_size * sizeof (char *));
          }

      for  (k = 0;  k < this_gap -> num_chunks;  k ++)
        {
         Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

         if  (REF (this_chunk -> chunk_id) . scaff_id != NULLINDEX
                && REF (this_chunk -> chunk_id) . is_unthrowable)
             this_chunk -> keep = FALSE;
         else if  (use_all || this_chunk -> keep)
             {
              forward_edge [ct] = reverse_edge [ct] = -1;
              this_chunk -> keep = TRUE;
              check [ct] = this_chunk;
              ct ++;
             }
        }

//      if  (ct == 0)
//          continue;

      if  (ct <= 1
             || (ct == 2 && this_gap -> left_cid >= 0
                         && this_gap -> right_cid >= 0))
          continue;    // only gap ends (already in scaffold)

      qsort (check, ct, sizeof (Gap_Chunk_t *), By_Low_Position);

      for  (p = 0;  p < ct;  p ++)
        sequence [p] = NULL;

      next_edge = 0;
      for  (p = 0;  p < ct - 1;  p ++)
        {
         for  (q = p + 1;  q < ct;  q ++)
           {
            Overlap  * olap;

            olap = Get_Chunk_Overlap
                       (check [p], check [q],
                        sequence + p, sequence + q, fp);
            if  (olap != NULL)
                {
                 if  (next_edge >= edge_pool_size - 1)
                     {
                      edge_pool_size *= 2;
PRALLOC (edge_pool_size * sizeof (Stone_Edge_t));
                      edge_pool
                          = (Stone_Edge_t *)
                              safe_realloc
                                (edge_pool,
                                 edge_pool_size * sizeof (Stone_Edge_t));
                     }
                 if  (olap -> endpos > 0)
                     {
                      edge_pool [next_edge] . from = p;
                      edge_pool [next_edge] . to = q;
                      edge_pool [next_edge] . progress = olap -> endpos;
                      edge_pool [next_edge] . a_hang = olap -> begpos;
                      edge_pool [next_edge] . length = olap -> length;
                      edge_pool [next_edge] . quality
                          = (double) olap -> diffs / olap -> length;
                      edge_pool [next_edge] . kind = OVERLAP;
                      edge_pool [next_edge] . next = forward_edge [p];
                      edge_pool [next_edge] . bad = FALSE;
if  (SKIP_CONTAINED_STONES && edge_pool [next_edge] . a_hang < 0)
 {
  continue;
 }
                      forward_edge [p] = next_edge;
                      next_edge ++;
                     }
                   else
                     {
                      edge_pool [next_edge] . from = q;
                      edge_pool [next_edge] . to = p;
                      edge_pool [next_edge] . progress = - (olap -> endpos);
                      edge_pool [next_edge] . a_hang = - (olap -> begpos);
                      edge_pool [next_edge] . length = olap -> length;
                      edge_pool [next_edge] . quality
                          = (double) olap -> diffs / olap -> length;
                      edge_pool [next_edge] . kind = OVERLAP;
                      edge_pool [next_edge] . next = forward_edge [q];
                      edge_pool [next_edge] . bad = FALSE;
if  (SKIP_CONTAINED_STONES && edge_pool [next_edge] . a_hang < 0)
 {
  continue;
 }
                      forward_edge [q] = next_edge;
                      next_edge ++;
                     }
                 if  (olap -> begpos > 0)
                     {
                      edge_pool [next_edge] . from = q;
                      edge_pool [next_edge] . to = p;
                      edge_pool [next_edge] . progress = olap -> begpos;
                      edge_pool [next_edge] . a_hang = olap -> endpos;
                      edge_pool [next_edge] . length = olap -> length;
                      edge_pool [next_edge] . quality
                          = (double) olap -> diffs / olap -> length;
                      edge_pool [next_edge] . kind = OVERLAP;
                      edge_pool [next_edge] . next = reverse_edge [q];
                      edge_pool [next_edge] . bad = FALSE;
if  (SKIP_CONTAINED_STONES && edge_pool [next_edge] . a_hang < 0)
 {
  continue;
 }
                      reverse_edge [q] = next_edge;
                      next_edge ++;
                     }
                   else
                     {
                      edge_pool [next_edge] . from = p;
                      edge_pool [next_edge] . to = q;
                      edge_pool [next_edge] . progress = - (olap -> begpos);
                      edge_pool [next_edge] . a_hang = - (olap -> endpos);
                      edge_pool [next_edge] . length = olap -> length;
                      edge_pool [next_edge] . quality
                          = (double) olap -> diffs / olap -> length;
                      edge_pool [next_edge] . kind = OVERLAP;
                      edge_pool [next_edge] . next = reverse_edge [p];
                      edge_pool [next_edge] . bad = FALSE;
if  (SKIP_CONTAINED_STONES && edge_pool [next_edge] . a_hang < 0)
 {
  continue;
 }
                      reverse_edge [p] = next_edge;
                      next_edge ++;
                     }
                }
              else
                {
                 // search for link from p to q
                 GraphEdgeIterator  ci_edges;
                 CIEdgeT  * edge;

                 InitGraphEdgeIterator
                     (ScaffoldGraph->RezGraph, check [p] -> chunk_id,
//                         ALL_END,
                      check [p] -> flipped ? A_END : B_END,
                      ALL_EDGES, GRAPH_EDGE_DEFAULT,
                      & ci_edges);

                 while  ((edge = NextGraphEdgeIterator (& ci_edges)) != NULL)
                   {
                    ChunkInstanceT  * other_chunk;

                    if  (isProbablyBogusEdge (edge)
                           || isSloppyEdge (edge))
                        continue;

                    if  (edge -> idA == check [p] -> chunk_id)
                        other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                         edge -> idB);
                      else
                        other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                         edge -> idA);

                    if  (other_chunk -> id == check [q] -> chunk_id)
                        {
                         // add to edges--not done yet
                        }
                   }
                }
           }
        }

      for  (p = 0;  p < ct;  p ++)
        if  (sequence [p] != NULL)
            free (sequence [p]);

if  (use_all && Use_Partial_Stone_Paths)   // for stones only
{
 int  start_sub = -1, target_sub = -1;
 int  num_kept, complete_path;
 int  i;

 for  (i = 0;  i < ct;  i ++)
   {
    if  (check [i] -> chunk_id == this_gap -> left_cid)
        start_sub = i;
    if  (check [i] -> chunk_id == this_gap -> right_cid)
        target_sub = i;
   }

#if  VERBOSE
 fprintf (stderr, "\nScaff %d  Gap %d\n", scaff_id, j);
 fprintf (stderr, "Before Build_Path_Subgraph\n");
 fprintf (stderr, "  Gap start = (%.0f,%.0f)  end = (%.0f,%.0f)\n",
          this_gap -> start . mean, sqrt (this_gap -> start . variance),
          this_gap -> end . mean, sqrt (this_gap -> end . variance));
#endif

 Build_Path_Subgraph
     (start_sub, target_sub,
      check, ct, forward_edge, reverse_edge, edge_pool,
      sorted, & num_sorted);

 num_kept = Determine_Components
                (sorted, num_sorted, check, ct,
                 start_sub, target_sub, forward_edge, edge_pool,
                 this_gap -> start . mean, 1.0, & target_position,
                 & complete_path, this_gap);

 if  (num_kept == 0)
     {
      for  (i = 0;  i < ct;  i ++)
        check [i] -> keep = FALSE;
      continue;
     }
 if  (start_sub != -1 && target_sub != -1 && ! complete_path)
     continue;
}


      if  (this_gap -> left_cid >= 0)
          Reject_Non_Reachable (this_gap -> left_cid, check, ct,
                                forward_edge, edge_pool);
      if  (this_gap -> right_cid >= 0)
          Reject_Non_Reachable (this_gap -> right_cid, check, ct,
                                reverse_edge, edge_pool);

#if  VERBOSE
fprintf (fp, "Forward edges:  ct = %d\n", ct);
for  (p = 0;  p < ct;  p ++)
{
fprintf (fp, "  From %d\n", check [p] -> chunk_id);
for  (q = forward_edge [p];  q != -1;  q = edge_pool [q] . next)
  fprintf (fp, "    to %4d  ahg = %5d  prog = %5d  qual = %6.4f  bad = %c\n",
           check [edge_pool [q] . to] -> chunk_id,
           edge_pool [q] . a_hang,
           edge_pool [q] . progress,
           edge_pool [q] . quality,
           edge_pool [q] . bad ? 'T' : 'F');
}
fflush (fp);
#endif

#if  VERBOSE
{
 int  i;

 fprintf (fp, "Keep:\n");
 for  (i = 0;  i < ct;  i ++)
   if  (check [i] -> keep)
       fprintf (fp, "  %5d\n", check [i] -> chunk_id);
}
#endif

      if  (this_gap -> left_cid >= 0)
          {
           double  start_coord;

           if  (this_gap -> len >= 0)
               start_coord = this_gap -> start . mean;
             else
               start_coord = this_gap -> end . mean;
           num_kept
               = Choose_Best_Stones
                    (this_gap -> left_cid, this_gap -> right_cid,
                     check, ct, forward_edge, edge_pool,
                     start_coord, 1.0,
                     & target_position);
          }
        else
          num_kept
              = Choose_Best_Stones
                    (this_gap -> right_cid, -1,
                     check, ct, reverse_edge, edge_pool,
                     this_gap -> end . mean, -1.0,
                     & target_position);

      if  (num_kept < 2)
          continue;

      if  (this_gap -> left_cid < 0)           // Left end gap
          Reverse_Positions (this_gap);
      else if  (this_gap -> right_cid >= 0     // Not right end gap
                  && num_kept > 2)
          {
           if  (this_gap -> len >= 0)
               {
                this_gap -> adjustment . mean
                    = target_position . mean - this_gap -> end . mean;
                this_gap -> adjustment . variance
                    = target_position . variance
                        - this_gap -> end . variance
                            + this_gap -> start . variance
                                + EPSILON;
               }
             else
               {
                this_gap -> adjustment . mean
                    = target_position . mean - this_gap -> start . mean;
                this_gap -> adjustment . variance
                    = target_position . variance
                        - this_gap -> start . variance
                            + this_gap -> end . variance
                                + EPSILON;
               }
          }
     }

   CheckScaffoldGraphCache (ScaffoldGraph);

   return;
  }



static int  Num_Keep_Entries
    (Scaffold_Fill_t * fill, int scaff_id)

//  Return the number of entries in all the gaps of  fill [scaff_id]
//  that have their  keep  flag  set.

  {
   int  j, ct = 0;

   for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
     {
      Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
      int  k;

      for  (k = 0;  k < this_gap -> num_chunks;  k ++)
        {
         Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

         if  (this_chunk -> keep)
             ct ++;
        }
     }

   return  ct;
  }



static void  Output_Cam_Files
    (Scaffold_Fill_t * fill)

//  Output information in global  Chunk_Info  to  Cam_File  and
//  info in  fill  to  Calc_Cam_File .

  {
   char  annotation_string [MAX_STRING_LEN];
   int  scaff_id;
   int  j;
#if  MAKE_CAM_FILE
   double  scaff_start = 0.0, scaff_end = 0.0;
#endif
#if  CHECK_CELSIM_COORDS
   ChunkInstanceT  * chunk;
   GraphNodeIterator  nodes;
   int cid;
#endif

   if  (! MAKE_CAM_FILE)
       return;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
       {
        Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
        int  k;

        for  (k = 0;  k < this_gap -> num_chunks;  k ++)
          {
           Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

           // Use  keep  so that stop overwriting annotation
           // after first  keep  entry in fill.

           if  (Chunk_Info [this_chunk -> chunk_id] . keep)
               continue;

           sprintf (annotation_string,
                    "copy %c  start = <%.0f,%.0f>  end = <%.0f,%.0f>",
                    this_chunk -> copy_letter,
                    this_chunk -> start . mean,
                    sqrt (this_chunk -> start . variance),
                    this_chunk -> end . mean,
                    sqrt (this_chunk -> end . variance));
           Chunk_Info [this_chunk -> chunk_id] . annotation
               = strdup (annotation_string);
           Chunk_Info [this_chunk -> chunk_id] . keep
               = this_chunk -> keep;
          }
       }
   

#if  CHECK_CELSIM_COORDS
   InitGraphNodeIterator(&nodes,ScaffoldGraph->RezGraph, GRAPH_NODE_DEFAULT);

   fprintf(stderr,"* Output_Cam_Files: num_chunks = %d\n", Num_Chunks);

   while  ((chunk = NextGraphNodeIterator (& nodes)) != NULL)
     {
      ContigT  * contig;
      int  row, mid;

      cid = chunk->id;

      if  (chunk -> flags . bits . isStoneSurrogate)
          continue;

      if  (cid >= Num_Chunks)
          {
           fprintf (stderr, "cid = %d  too big, skipped\n", cid);
           continue;
          }

      if  (Chunk_Info [cid] . celsim_left < 0
             || Chunk_Info [cid] . celsim_right < 0)
          continue;


      contig = GetGraphNode (ScaffoldGraph -> RezGraph, cid);
      assert (contig != NULL);

      mid = (Chunk_Info [cid] . celsim_left
              + Chunk_Info [cid] . celsim_right) / 2;
      
      switch  (Chunk_Info [cid] . colour)
        {
         case  UNIQUE_COLOUR :
           row = 1;
#if  0
           chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);
           if  (chunk -> type != DISCRIMINATORUNIQUECHUNK_CGW)
               Chunk_Info [cid] . colour = INSERTED_COLOUR;
#endif
           fprintf (Cam_File, "%dCHUNKREZ: %d",
                    cid, Chunk_Info [cid] . celsim_left);
           if  (Chunk_Info [cid] . cgb_type != UU_CGBTYPE)
               fprintf (Cam_File, " A%dREZ M%d",
                        Repeat_Colour (Chunk_Info [cid] . cgb_type),
                        mid);
           fprintf (Cam_File,
                    " A%dREZ %d R%d # chunk %d"
                    "  Scaff #%d  rel pos #%d  start = <%.0f,%.0f>"
                    "  end = <%.0f,%.0f>  cov = %d  typ = %s\n",
                    Chunk_Info [cid] . colour,
                    Chunk_Info [cid] . celsim_right, row, cid,
                    REF (cid) . scaff_id, REF (cid) . rel_pos,
                    contig -> offsetAEnd . mean,
                    sqrt (contig -> offsetAEnd . variance),
                    contig -> offsetBEnd . mean,
                    sqrt (contig -> offsetBEnd . variance),
                    GetCoverageStat (contig),
                    CGB_Type_As_String (contig -> flags . bits . cgbType));
           break;
         case  PLACED_COLOUR :
         case  MISPLACED_COLOUR :
         case  REJECT_COLOUR :
           row = 3;
           fprintf (Cam_File, "%dCHUNKREZ: %d",
                    cid, Chunk_Info [cid] . celsim_left);
           if  (Chunk_Info [cid] . cgb_type != UU_CGBTYPE)
               fprintf (Cam_File, " A%dREZ M%d",
                        Repeat_Colour (Chunk_Info [cid] . cgb_type),
                        mid);
           fprintf (Cam_File,
                    "A%dREZ %d R%d # chunk %d"
                    "  %s  cov = %d  typ = %s\n",
                    Chunk_Info [cid] . colour,
                    Chunk_Info [cid] . celsim_right, row, cid,
                    Chunk_Info [cid] . annotation,
                    GetCoverageStat (contig),
                    CGB_Type_As_String (contig -> flags . bits . cgbType));
           break;
         default :
           row = 4;
           fprintf (Cam_File, "%dCHUNKREZ: %d",
                    cid, Chunk_Info [cid] . celsim_left);
           if  (Chunk_Info [cid] . cgb_type != UU_CGBTYPE)
               fprintf (Cam_File, " A%dREZ M%d",
                        Repeat_Colour (Chunk_Info [cid] . cgb_type),
                        mid);
           fprintf (Cam_File,
                    "A%dREZ %d R%d # chunk %d"
                    "  %s  cov = %d  typ = %s\n",
                    Chunk_Info [cid] . colour,
                    Chunk_Info [cid] . celsim_right, row, cid,
                    Chunk_Info [cid] . annotation,
                    GetCoverageStat (contig),
                    CGB_Type_As_String (contig -> flags . bits . cgbType));
        }
     }
#endif


#if  SHOW_CALC_COORDS
   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      double  scaff_contig_lo, scaff_contig_hi;
      int  j;

      for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
         int  k;

         if  (j > 0)
             {
              if  (this_gap -> len >= 0)
                  scaff_contig_hi = this_gap -> start . mean;
                else
                  scaff_contig_hi = this_gap -> end . mean;
#if  MAKE_CAM_FILE
              fprintf (Calc_Cam_File,
                       "%dCHUNKREZ: %.0f A%dREZ %.0f R%d # contig %d  %s\n",
                       this_gap -> left_cid,
                       scaff_start + scaff_contig_lo,
                       UNIQUE_COLOUR,
                       scaff_start + scaff_contig_hi,
                       1, this_gap -> left_cid,
                       Chunk_Info [this_gap -> left_cid] . annotation);
#endif
             }
         if  (this_gap -> len >= 0)
             scaff_contig_lo = this_gap -> end . mean;
           else
             scaff_contig_lo = this_gap -> start . mean;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;
            int  colour, row;
            double  lo, hi;

            if  (this_chunk -> copy_letter == GAP_END_CHAR)
                continue;
            else if  (this_chunk -> keep)
                {
                 colour = INSERTED_COLOUR;
                 row = 3;
                }
              else
                {
                 colour = REJECT_COLOUR;
                 row = 4;
                }

            lo = Min_double (this_chunk -> start . mean,
                             this_chunk -> end . mean);
            hi = Max_double (this_chunk -> start . mean,
                             this_chunk -> end . mean);
            if  (hi > scaff_end)
                scaff_end = hi;

#if  MAKE_CAM_FILE
            fprintf (Calc_Cam_File,
                     "%dCHUNK%dg%d: %.0f A%dREZ %.0f R%d # contig %d  %s\n",
                     this_chunk -> chunk_id, scaff_id, j,
                     scaff_start + lo,
                     colour, scaff_start + hi, row,
                     this_chunk -> chunk_id,
                     Chunk_Info [this_chunk -> chunk_id] . annotation);
#endif
           }
        }
      if  (scaff_end > scaff_start)
          scaff_start = scaff_end + CALC_SCAFF_SEPARATION;
     }
#endif

   return;
  }



static void  Print_Scaffolds
    (FILE * fp)

//  Print list of all scaffolds and the chunks in them, sending output
//  to  fp .  Set global  Ref  to array of references of scaffold and relative
//  position in it for each scaffold.  If  fp  is  NULL  do no
//  actual printing.

  {
   GraphNodeIterator scaffolds;
   CIScaffoldT  * scaffold;
   int  scaff_id;
   int  last_coord_used = 0;

   if  (fp != NULL)
       fprintf (fp, "\nScaffolds\n");

   InitGraphNodeIterator (& scaffolds, ScaffoldGraph -> ScaffoldGraph,
                          GRAPH_NODE_DEFAULT);
   while  ((scaffold = NextGraphNodeIterator (& scaffolds)) != NULL)
     {
      CIScaffoldTIterator  scaff_iterator;
      ChunkInstanceT  * chunk;
      double  prev_variance;
      int  scaff_index;
#if  MAKE_CAM_FILE && SHOW_CALC_COORDS
      int  seg_num = 0;
#endif
#if CHECK_CELSIM_COORDS
      double  slope = 0.0, y_intercept;
#endif

      scaff_id = scaffold -> id;
      if  (fp != NULL)
          {
           fprintf (fp, "Scaffold %3d:\n", scaff_id);
           fprintf (fp, "  %3s  %7s    (%8s, %8s)  (%8s, %8s)"
#if  CHECK_CELSIM_COORDS
           "  [%8s  %-8s]"
#endif
           "  %7s\n",
                    "Idx",  "Chunk", "Start", "Variance",
                    "End", "Variance",
#if  CHECK_CELSIM_COORDS
                    "Celsim", "Coords",
#endif
                    "Len");
          }

#if  MAKE_CAM_FILE
      fprintf (Cam_File, "LNK: ");
#if  SHOW_CALC_COORDS
      fprintf (Calc_Cam_File, "LNK: ");
#endif
#endif

#if  CHECK_CELSIM_COORDS
      chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                 scaffold -> info.Scaffold.AEndCI);

      assert(chunk->scaffoldID == scaff_id);

      if  (chunk -> offsetAEnd . mean < chunk -> offsetBEnd . mean)
          Scaffold_Start [scaff_id] = chunk -> aEndCoord;
        else
          Scaffold_Start [scaff_id] = chunk -> bEndCoord;
      chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                 scaffold -> info.Scaffold.BEndCI);
      if  (chunk -> offsetAEnd . mean < chunk -> offsetBEnd . mean)
          Scaffold_End [scaff_id] = chunk -> bEndCoord;
        else
          Scaffold_End [scaff_id] = chunk -> aEndCoord;
#else
      last_coord_used += MAX_MATE_DISTANCE;
      Scaffold_Start [scaff_id] = last_coord_used;
      chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                 scaffold -> info.Scaffold.BEndCI);
      if  (chunk -> offsetAEnd . mean < chunk -> offsetBEnd . mean)
          last_coord_used += (int) chunk -> offsetBEnd . mean;
        else
          last_coord_used += (int) chunk -> offsetAEnd . mean;
      Scaffold_End [scaff_id] = last_coord_used;
#endif

      InitCIScaffoldTIterator (ScaffoldGraph, scaffold, TRUE, FALSE,
                              & scaff_iterator);
      prev_variance = -1.0;

      for  (scaff_index = 0;
            (chunk = NextCIScaffoldTIterator (& scaff_iterator)) != NULL;
            scaff_index ++)
        {
         double  dip;
#if CHECK_CELSIM_COORDS
         double  denom;
#endif

         if  (fp != NULL)
             fprintf (fp,
                  "  %3d: %7d    (%8.0f, %8.0f)  (%8.0f, %8.0f)"
#if  CHECK_CELSIM_COORDS
                  "  [%8ld, %8ld]"
#endif
                  "  %7.0f\n",
                  scaff_index,
                  chunk -> id,
                  chunk -> offsetAEnd . mean,
                  chunk -> offsetAEnd . variance,
                  chunk -> offsetBEnd . mean,
                  chunk -> offsetBEnd . variance,
#if  CHECK_CELSIM_COORDS
                  chunk -> aEndCoord,
                  chunk -> bEndCoord,
#endif
                  chunk -> bpLength . mean
                 );

         if  (chunk -> offsetAEnd . mean <= chunk -> offsetBEnd . mean)
             Chunk_Info [chunk -> id] . flipped = TRUE;
           else
             Chunk_Info [chunk -> id] . flipped = FALSE;

         if  (chunk -> offsetAEnd . variance < 0
                || chunk -> offsetBEnd . variance < 0)
             {
              fprintf (stderr,
                       "ERROR:  Negative variance for chunk %d in scaffold %d\n",
                       chunk -> id, scaff_id);
              if  (fp != NULL)
                  fclose (fp);
              assert (FALSE);
             }
         if  ((dip = chunk -> offsetAEnd . variance - prev_variance) < 0.0
                || (dip = chunk -> offsetBEnd . variance - prev_variance) < 0.0)
             {
              fprintf (stderr,
                       "ERROR:  Variance dip by %f for chunk %d in scaffold %d\n",
                       dip, chunk -> id, scaff_id);
              if  (fp != NULL)
                  fclose (fp);
              assert (FALSE);
             }
#if CHECK_CELSIM_COORDS
         if  (scaff_index > 0)
             {
              double  delta, target, diff, tolerance;

              delta = 3.0 * sqrt (chunk -> offsetAEnd . variance
                                    - prev_variance);
              target = slope * chunk -> offsetAEnd . mean + y_intercept;
              diff = fabs (chunk -> aEndCoord - target);
              tolerance = fabs (slope * delta);
              if  (diff > tolerance && fp != NULL)
                  fprintf (fp,
                  "OOPS  AEnd:  dif = %.0f tol = %.0f slope = %.1f y_int = %.1f\n",
                           diff, tolerance, slope, y_intercept);

              delta = 5.0 * sqrt (chunk -> offsetBEnd . variance
                                    - prev_variance);
              target = slope * chunk -> offsetBEnd . mean + y_intercept;
              diff = fabs (chunk -> bEndCoord - target);
              tolerance = fabs (slope * delta);
              if  (tolerance < 50.0)
                  tolerance = 50.0;
              if  (diff > tolerance && fp != NULL)
                  fprintf (fp,
                  "OOPS  BEnd:  dif = %.0f tol = %.0f slope = %.1f y_int = %.1f\n",
                           diff, tolerance, slope, y_intercept);
             }

         prev_variance = Max_double (chunk -> offsetAEnd . variance,
                                     chunk -> offsetBEnd . variance);

         denom = chunk -> offsetBEnd . mean - chunk -> offsetAEnd . mean;
         if  (fabs (denom) < 1e-6)
             {
              fprintf (stderr, "\nERROR:  Zero denominator for chunk %d\n",
                       chunk -> id);
              fprintf (stderr, "        offsetA = %.0f  offsetB = %.0f\n",
                       chunk -> offsetAEnd . mean,
                       chunk -> offsetBEnd . mean);
             }
           else
             {
              slope = (chunk -> bEndCoord - chunk -> aEndCoord)
                        / denom;
              if  (fabs (fabs (slope) - 1.0) > 1.0)
                  {
                   if  (fp != NULL)
                       fprintf (fp, "OOPS:  Bad slope = %.3f\n", slope);
                  }
              else if  (slope > 0.0)
                  slope = 1.0;
                else
                  slope = -1.0;
              y_intercept = chunk -> aEndCoord - slope * chunk -> offsetAEnd . mean;
             }
#endif

#if  MAKE_CAM_FILE
#if  SHOW_CALC_COORDS
         fprintf (Calc_Cam_File, " %dCHUNKREZ", chunk -> id);
         if  ((scaff_index + 1) % MAX_SCAFFOLD_SIZE == 0
                && scaff_index != scaffold -> info.Scaffold.numElements - 1)
             {
              fprintf (Calc_Cam_File, " A%dREZ  # Scaffold #%d segment %d\n",
                       SCAFFOLD_COLOUR, scaff_id, seg_num);
              fprintf (Calc_Cam_File, "LNK: ");
              fprintf (Calc_Cam_File, " %dCHUNKREZ", chunk -> id);
             }
#endif
         fprintf (Cam_File, " %dCHUNKREZ", chunk -> id);
         if  ((scaff_index + 1) % MAX_SCAFFOLD_SIZE == 0
                && scaff_index != scaffold -> info.Scaffold.numElements - 1)
             {
              fprintf (Cam_File, " A%dREZ  # Scaffold #%d segment %d\n",
                       SCAFFOLD_COLOUR, scaff_id, seg_num ++);
              fprintf (Cam_File, "LNK: ");
              fprintf (Cam_File, " %dCHUNKREZ", chunk -> id);
             }
#endif


         if  (REF (chunk -> id) . scaff_id != scaff_id)
             {
              fprintf (stderr, "chunk -> id = %d\n", chunk -> id);
              fprintf (stderr, "REF . scaff_id = %d\n", REF (chunk -> id) . scaff_id);
              fprintf (stderr, "scaff_id = %d\n", scaff_id);
              if  (fp != NULL)
                  fflush (fp);
             }
         assert (REF (chunk -> id) . scaff_id == scaff_id);
         REF (chunk -> id) . rel_pos = scaff_index;
        }

      Scaffold_Flipped [scaff_id] = Scaffold_Start [scaff_id] > Scaffold_End [scaff_id];


#if  MAKE_CAM_FILE
      fprintf (Cam_File, " A%dREZ  # Scaffold #%d", SCAFFOLD_COLOUR, scaff_id);
      if  (seg_num == 0)
          fprintf (Cam_File, "\n");
        else
          fprintf (Cam_File, " segment %d\n", seg_num);
#if  SHOW_CALC_COORDS
      fprintf (Calc_Cam_File, " A%dREZ  # Scaffold #%d", SCAFFOLD_COLOUR, scaff_id);
      if  (seg_num == 0)
          fprintf (Calc_Cam_File, "\n");
        else
          fprintf (Calc_Cam_File, " segment %d\n", seg_num);
#endif
#endif
     }

   // List all scaffold edges

   if  (fp != NULL)
       fprintf (fp, "\nScaffold Edges\n");

   InitGraphNodeIterator (& scaffolds, ScaffoldGraph -> ScaffoldGraph,
                          GRAPH_NODE_DEFAULT);
   while  ((scaffold = NextGraphNodeIterator (& scaffolds)) != NULL)
     {
      GraphEdgeIterator SEdges;
      SEdgeT  * edge;

      scaff_id = scaffold -> id;
      if  (fp != NULL)
          {
           fprintf (fp, "Scaffold %3d:\n", scaff_id);

           fprintf (fp, "A Edges\n");
           InitGraphEdgeIterator(ScaffoldGraph->ScaffoldGraph, scaff_id, A_END, ALL_EDGES, 
                               GRAPH_EDGE_DEFAULT,   &SEdges);

           //      InitSEdgeTIterator (ScaffoldGraph, scaff_id, FALSE, FALSE,
           //                          A_END, FALSE,  & SEdges);
           while  ((edge = NextGraphEdgeIterator (& SEdges)) != NULL)
             {
              fprintf (fp, " %3d -> %3d  [%8.0f,%8.0f]  %5s  %3d  %4s  %5s\n",
                       edge -> idA, edge -> idB,
                       edge -> distance . mean,
                       sqrt (edge -> distance . variance),
                       Orientation_As_String (edge -> orient),
                       edge -> edgesContributing,
                       Is_Good_Scaff_Edge (edge) ? "good" : "bad",
                       edge -> flags . bits . isBogus ? "bogus" : "valid");
             }

           fprintf (fp, "B Edges\n");
           InitGraphEdgeIterator(ScaffoldGraph->ScaffoldGraph, scaff_id, B_END, ALL_EDGES, 
                               GRAPH_EDGE_DEFAULT,   &SEdges);

           //InitSEdgeTIterator (ScaffoldGraph, scaff_id, FALSE, FALSE,
           //                          B_END, FALSE,  & SEdges);
           while  ((edge = NextGraphEdgeIterator (& SEdges)) != NULL)
             {
              fprintf (fp, " %3d -> %3d  [%8.0f,%8.0f]  %5s  %3d  %4s  %5s\n",
                       edge -> idA, edge -> idB,
                       edge -> distance . mean,
                       sqrt (edge -> distance . variance),
                       Orientation_As_String (edge -> orient),
                       edge -> edgesContributing,
                       Is_Good_Scaff_Edge (edge) ? "good" : "bad",
                       edge -> flags . bits . isBogus ? "bogus" : "valid");
             }
          }
     }

#if  CHECK_CELSIM_COORDS
   if  (fp != NULL)
       {
        fprintf (fp, "\nCelsim coordinates of scaffolds\n");

        for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
          fprintf (fp, "Scaffold %2d:  start = %7ld  end = %7ld\n",
                   scaff_id, Scaffold_Start [scaff_id], Scaffold_End [scaff_id]);
       }
#endif

   return;
  }



static void  Partition_Edges
    (int cid, Stack_Entry_t * stack, int stack_top, int min_good_links)

//  Partition entries in  stack [0 .. (stack_top - 1)]  based on
//  which scaffold and orientation they connect with.  Mark as
//  bad all edges whose partition set collectively does not contain
//  at least min_good_links  mate links .

  {
   int  i, j;
   int  partition;

   // Sort by scaffold id and  whether flipped w.r.t. that
   // scaffold

//   qsort (stack, stack_top, sizeof (Stack_Entry_t), By_Scaff_And_Flipped);
   qsort (stack, stack_top, sizeof (Stack_Entry_t), By_Scaff_Flipped_And_Left_End);

   // Mark good entries in stack

#if  VERBOSE
fprintf (stderr, ">>> Partition #%d:\n", cid);
#endif

   partition = 0;
   for  (i = 0;  i < stack_top;  i = j)
     {
      int  k, total_mates = 0;
      double  right_extent;
      double  new_left, new_right;

      new_left = stack [i] . left_end . mean
                    - 3.0 * sqrt (stack [i]  . left_end . variance);
      right_extent = stack [i] . right_end . mean
                       + 3.0 * sqrt (stack [i] . right_end . variance);
#if  VERBOSE
fprintf (stderr, "   1st_left = %8.0f  1st_right = %8.0f\n",
         new_left, right_extent);
#endif

      for  (j = i + 1;
              j < stack_top
                && stack [j] . scaff_id == stack [i] . scaff_id
                && stack [j] . flipped == stack [i] . flipped;  j ++)
           {
            new_left = stack [j] . left_end . mean
                          - 3.0 * sqrt (stack [j]  . left_end . variance);
            if  (new_left > right_extent)
                break;

            new_right = stack [j] . right_end . mean
                          + 3.0 * sqrt (stack [j]  . right_end . variance);
            if  (new_right > right_extent)
                right_extent = new_right;
#if  VERBOSE
fprintf (stderr, "   new_left = %8.0f  new_right = %8.0f\n",
         new_left, new_right);
#endif
           }

      for  (k = i;  k < j;  k ++)
        {
         total_mates += stack [k] . num_good_mates;
         stack [k] . partition = partition;
        }
#if  VERBOSE
fprintf (stderr, "  %3d  %3d .. %3d  total_mates = %3d\n",
         partition, i, j - 1, total_mates);
#endif
      partition ++;

      if  (total_mates < min_good_links)
          {
           for  (k = i;  k < j;  k ++)
             stack [k] . is_bad = TRUE;
          }
     }

   return;
  }



#if  CHECK_CELSIM_COORDS
static void  Print_Coverage
    (FILE * fp)

//  Print to  fp  total coverage of genome by uniques, and by uniques
//  together with placed's.

  {
   int64  * lo, * hi, sum;
   int  i, j, n;

PALLOC (Num_Chunks * sizeof (int64));
   lo = (int64 *) safe_malloc (Num_Chunks * sizeof (int64));
PALLOC (Num_Chunks * sizeof (int64));
   hi = (int64 *) safe_malloc (Num_Chunks * sizeof (int64));

   for  (i = n = 0;  i < Num_Chunks;  i ++)
     if  (Chunk_Info [i] . colour == UNIQUE_COLOUR)
         {
          lo [n] = Chunk_Info [i] . celsim_left;
          hi [n] = Chunk_Info [i] . celsim_right;
          n ++;
         }

   for  (i = 0;  i < n - 1;  i ++)
     for  (j = i + 1;  j < n;  j ++)
       if  (lo [i] > lo [j])
           {
            int64  save;

            save = lo [i];
            lo [i] = lo [j];
            lo [j] = save;

            save = hi [i];
            hi [i] = hi [j];
            hi [j] = save;
           }

   j = 0;
   for  (i = 1;  i < n;  i ++)
     if  (lo [i] > hi [j])
         {
          j ++;
          lo [j] = lo [i];
          hi [j] = hi [i];
         }
     else if  (hi [j] < hi [i])
          hi [j] = hi [i];

   sum = 0;
   for  (i = 0;  i <= j;  i ++)
     sum += hi [i] - lo [i];

   fprintf (fp, "Coverage by uniques = %ld bases\n", sum);

   for  (i = 0, n = 0;  i < Num_Chunks;  i ++)
     if  (Chunk_Info [i] . colour == UNIQUE_COLOUR
             || Chunk_Info [i] . colour == PLACED_COLOUR)
         {
          lo [n] = Chunk_Info [i] . celsim_left;
          hi [n] = Chunk_Info [i] . celsim_right;
          n ++;
         }

   for  (i = 0;  i < n - 1;  i ++)
     for  (j = i + 1;  j < n;  j ++)
       if  (lo [i] > lo [j])
           {
            int64  save;

            save = lo [i];
            lo [i] = lo [j];
            lo [j] = save;

            save = hi [i];
            hi [i] = hi [j];
            hi [j] = save;
           }

   j = 0;
   for  (i = 1;  i < n;  i ++)
     if  (lo [i] > hi [j])
         {
          j ++;
          lo [j] = lo [i];
          hi [j] = hi [i];
         }
     else if  (hi [j] < hi [i])
          hi [j] = hi [i];

   sum = 0;
   for  (i = 0;  i <= j;  i ++)
     sum += hi [i] - lo [i];

   fprintf (fp, "Coverage by uniques & placed = %ld bases\n", sum);

   free (lo);
   free (hi);

   return;
  }
#endif



void  Print_Fill_Info
    (FILE * fp, Scaffold_Fill_t * fill_chunks)

//  List the locations of all gaps in scaffolds in  fill_chunks .  Output
//  goes to  fp .

  {
   int  scaff_id;
   int  total_chunks = 0;
   int  total_keep = 0;

   fprintf (fp, "\n*** Chunks to fill in gaps ***\n");

   Global_Fill_Info_Size = Num_Scaffolds * sizeof (Scaffold_Fill_t);

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     Print_Fill_Info_One_Scaffold
         (fp, fill_chunks, scaff_id, & total_chunks, & total_keep);

   fprintf (fp, "\n*** Total chunks in fill_info = %d ***\n", total_chunks);
   fprintf (fp, "***               Keep chunks = %d ***\n", total_keep);
   fprintf (fp, "***        Bytes in fill_info = %ld ***\n",
            Global_Fill_Info_Size);

   return;
  }


void  Print_Fill_Info_One_Scaffold
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int scaff_id,
     int * total_chunks, int * total_keep)

//  List the locations of all gaps in scaffold  fill_chunks [scaff_id] .
//  Output goes to  fp .  Add to  (* total_chunks)  and  (* total_keep)
//  the number of chunks (marked  keep , resp.) in the gaps.

  {
   int  num_entries;
   int  j;

   fprintf (fp, "\n Scaffold #%d  num_gaps = %d:\n",
            scaff_id, fill_chunks [scaff_id] . num_gaps);

#if  TEST_HOPELESS_SCAFFS
   if  (Is_Hopeless_Scaff [scaff_id] & Hopeless_True_Mask)
       {
        fprintf (fp, "   HOPELESS\n");
        return;
       }
#endif

   Global_Fill_Info_Size
       += fill_chunks [scaff_id] . num_gaps * sizeof (Gap_Fill_t);

   num_entries = 0;
   for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
     num_entries += fill_chunks [scaff_id] . gap [j] . num_chunks;
   if  (num_entries == 0)
       {
        fprintf (fp, "   EMPTY\n");
        return;
       }


   for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
     {
      Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
      ChunkInstanceT  * scaff_chunk;
      int  k;
      double ref_var;
#if  CHECK_CELSIM_COORDS
      double  slope = 0.0, intercept = 0.0;
#endif
      
      if  (j > 0)
          {
           scaff_chunk
               = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> left_cid);
#if  0
           fprintf (fp,
                    " Scaffold Chunk %8d:  (%8.0f,%7.0f)  (%8.0f,%7.0f)\n",
                    scaff_chunk -> id,
                    scaff_chunk -> offsetAEnd . mean,
                    scaff_chunk -> offsetAEnd . variance,
                    scaff_chunk -> offsetBEnd . mean,
                    scaff_chunk -> offsetBEnd . variance);
#endif
           ref_var = Max_double (scaff_chunk -> offsetAEnd . variance,
                                 scaff_chunk -> offsetBEnd . variance);
          }
        else
          {
           scaff_chunk
               = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> right_cid);
           ref_var = Min_double (scaff_chunk -> offsetAEnd . variance,
                                 scaff_chunk -> offsetBEnd . variance);
          }

#if  CHECK_CELSIM_COORDS
      if  (scaff_chunk -> aEndCoord >= 0
             && scaff_chunk -> bEndCoord >= 0)
          {
           assert (scaff_chunk -> offsetAEnd . mean
                     != scaff_chunk -> offsetBEnd . mean);
           slope = (scaff_chunk -> bEndCoord - scaff_chunk -> aEndCoord)
                     / (scaff_chunk -> offsetBEnd . mean
                          - scaff_chunk -> offsetAEnd . mean);
           if  (fabs (fabs (slope) - 1.0) > 0.1)
               fprintf (fp, "YIKES:  Bad slope = %.3f  [%7d,%7d]\n",
                        slope, scaff_chunk -> aEndCoord,
                        scaff_chunk -> bEndCoord);
           else if  (slope > 0.0)
               slope = 1.0;
             else
               slope = -1.0;

           intercept = scaff_chunk -> bEndCoord
                         - slope * scaff_chunk -> offsetBEnd . mean;
          }
#endif

      fprintf (fp,
               " Gap %3d:  (%8.0f,%7.0f)  (%8.0f,%7.0f)  <%6d,%6d>  len = %.0f"
               "  ref_var = %.0f  adj = (%.0f,%.0f)\n",
               j,
               this_gap -> start . mean,
               this_gap -> start . variance,
               this_gap -> end . mean,
               this_gap -> end . variance,
               this_gap -> left_cid,
               this_gap -> right_cid,
               this_gap -> len,
               this_gap -> ref_variance,
               this_gap -> adjustment . mean,
               this_gap -> adjustment . variance);
{
ContigT  * contig;

if  (this_gap -> left_cid >= 0)
  {
   contig = GetGraphNode (ScaffoldGraph -> RezGraph, this_gap -> left_cid);
   fprintf (fp, "   Left end at [%7.0f,%7.0f]",
            contig -> offsetAEnd . mean, contig -> offsetBEnd . mean);
  }
if  (this_gap -> right_cid >= 0)
  {
   contig = GetGraphNode (ScaffoldGraph -> RezGraph, this_gap -> right_cid);
   fprintf (fp, "   Right end at [%7.0f,%7.0f]",
            contig -> offsetAEnd . mean, contig -> offsetBEnd . mean);
  }
fprintf (fp, "\n");
}

      Global_Fill_Info_Size
          += this_gap -> num_chunks * sizeof (Gap_Chunk_t);

      for  (k = 0;  k < this_gap -> num_chunks;  k ++)
        {
         Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;
         ChunkInstanceT  * contig;
         ChunkInstanceT  * chunk;
//            MultiAlignT  * ma;
//            int  ma_count;

         if  (REF (this_chunk -> chunk_id) . scaff_id != NULLINDEX
                && this_chunk -> copy_letter != GAP_END_CHAR
                && REF (this_chunk -> chunk_id) . is_unthrowable)
             {
              fprintf (fp, "   Chunk %6d%c is unthrowable\n",
                       this_chunk -> chunk_id, this_chunk -> copy_letter);
              break;
             }
             

         contig = GetGraphNode
                    (ScaffoldGraph -> RezGraph, this_chunk -> chunk_id);
         if  (contig -> flags . bits . isDead)
             {
              fprintf (fp,
                       "   Chunk %6d%c  <%7.0f,%7.0f>  <%7.0f,%7.0f>  %3d,%2d"
                       " %s %s %s %s **Previously Inserted Unique**\n",
                       this_chunk -> chunk_id,
                       this_chunk -> copy_letter,
                       this_chunk -> start . mean,
                       this_chunk -> start . variance,
                       this_chunk -> end . mean,
                       this_chunk -> end . variance,
                       this_chunk -> cover_stat,
                       this_chunk -> link_ct,
                       this_chunk -> keep ? "Keep" : "Rej",
                       this_chunk -> best ? "Best" : " ",
                       this_chunk -> path_confirmed ? "Path" : " ",
                       REF (this_chunk -> chunk_id) . is_singleton
                         && this_chunk -> copy_letter != GAP_END_CHAR ?
                              "*UNIQUE*" : ""
                      );
              this_chunk -> keep = FALSE;
              fprintf (stderr,
              "*** Chunk %d  Scaff %d  Gap %d is previously inserted unique ***\n",
                       this_chunk -> chunk_id, scaff_id, j);
              continue;
             }

         chunk = GetGraphNode
                   (ScaffoldGraph -> CIGraph, contig -> info . Contig . AEndCI);
#if  0
         ma = LoadMultiAlignTFromSequenceDB
                  (ScaffoldGraph -> sequenceDB, this_chunk -> chunk_id,
                   ScaffoldGraph -> RezGraph -> type == CI_GRAPH);
//            ma = GetMultiAlignInStore
//                   (ScaffoldGraph -> RezGraph -> maStore, this_chunk -> chunk_id);
         ma_count = GetNumIntMultiPoss (ma -> f_list);
#endif

         fprintf (fp,
#if  0
"   Chunk %6d%c  <%7.0f,%7.0f>  <%7.0f,%7.0f>  %3d,%2d  [%+7.0f] %2d %-4s %-4s %-4s\n",
"   Chunk %6d%c  <%7.0f,%7.0f>  <%7.0f,%7.0f>  %3d,%2d,%7.0f,%3d[%3d,%3d]  %2d %s %s %s %s\n",
#endif
"   Chunk %6d%c  <%7.0f,%7.0f>  <%7.0f,%7.0f>  %3d,%2d,%7.0f,%3d,%3d  %2d %s %s %s %s\n",
                  this_chunk -> chunk_id,
                  this_chunk -> copy_letter,
                  this_chunk -> start . mean,
                  this_chunk -> start . variance,
                  this_chunk -> end . mean,
                  this_chunk -> end . variance,
//                     this_chunk -> avg_edge_quality,
                  this_chunk -> cover_stat,
                  this_chunk -> link_ct,
#if  0
                  this_chunk -> start . mean
                    - this_chunk -> sim_start,
#else
                  contig -> bpLength . mean,
                  chunk -> info . CI . numFragments,
                  contig -> info . Contig . numCI,
//                     ma_count,
#endif
                  this_chunk -> index,
                  this_chunk -> keep ? "Keep" : "Rej",
                  this_chunk -> best ? "Best" : " ",
                  this_chunk -> path_confirmed ? "Path" : " ",
                  REF (this_chunk -> chunk_id) . is_singleton
                    && this_chunk -> copy_letter != GAP_END_CHAR ?
                         "*UNIQUE*" : ""
                 );
         if  (this_chunk -> keep)
             (* total_keep) ++;

#if  CHECK_CELSIM_COORDS
         if  (this_chunk -> keep)
             {
              ChunkInstanceT  * chunk
                  = GetGraphNode (ScaffoldGraph->RezGraph,
                                  this_chunk -> chunk_id);

              if  (chunk -> aEndCoord >= 0 && chunk -> bEndCoord >= 0)
                  {
                    double  mid, delta;
                   delta = 3.0 * sqrt (fabs (this_chunk -> start . variance
                                         - ref_var))
                             * fabs (slope);
                   mid = slope * (this_chunk -> start . mean) + intercept;
                   if  (fabs (chunk -> aEndCoord - mid) > delta)
                       fprintf (fp,
                           "YIKES:  A_End at %d  Should be in [%.0f, %.0f]\n",
                           chunk -> aEndCoord, mid - delta, mid + delta);

                   delta = 3.0 * sqrt (fabs (this_chunk -> end . variance
                                         - ref_var))
                             * fabs (slope);
                   mid = slope * (this_chunk -> end . mean) + intercept;
                   if  (fabs (chunk -> bEndCoord - mid) > delta)
                       fprintf (fp,
                           "YIKES:  B_End at %d  Should be in [%.0f, %.0f]\n",
                           chunk -> bEndCoord, mid - delta, mid + delta);
                  }
             }
#endif
         (* total_chunks) ++;
        }
     }

   return;
  }



static void  Print_Frag_Info
    (FILE * fp, int cid)

//  Print to  fp  info about the fragments in chunk  cid .

  {
   ChunkInstanceT  * contig
       = GetGraphNode (ScaffoldGraph -> RezGraph, cid);
   ChunkInstanceT  * chunk;
   MultiAlignT  * ma;
   int  cover_stat;
   int  i, chunk_id, num_frags;

   chunk_id = contig -> info . Contig . AEndCI;
   chunk = GetGraphNode (ScaffoldGraph -> CIGraph, chunk_id);
   cover_stat = GetCoverageStat (chunk);

//   assert (contig -> info . Contig . AEndCI
//             == contig -> info . Contig . BEndCI);

   ma = LoadMultiAlignTFromSequenceDB
            (ScaffoldGraph -> sequenceDB, cid,
             ScaffoldGraph -> RezGraph -> type == CI_GRAPH);
   assert (ma != NULL);

   // cycle through fragments 
   num_frags = GetNumIntMultiPoss (ma -> f_list);

   if  (contig -> info . Contig . numCI > 1)
       fprintf (fp,
                "  Contig = %6d  UNIQUE  frags = %d\n",
                cid, num_frags);
     else
       fprintf (fp,
                "  Contig = %6d  cov = %d  frags = %d\n",
                cid, cover_stat, num_frags);
   fprintf (fp, "    %17s %13s %28s %22s %6s %6s\n",
            "Frag ID/iid/type", "Chunk Coords",
            "Mate Frag/Chunk & Coords", "Contig ID & Coords ",
            "Scaff", "Dist CI/Ctg");

   for  (i = 0;  i < num_frags;  i ++)
     {
      IntMultiPos  * mp = GetIntMultiPos (ma -> f_list, i);
      CDS_CID_t  fragID = (CDS_CID_t) mp -> sourceInt;
        // This is an internal-data-structure ID
      CDS_CID_t  ident = (CDS_CID_t) mp -> ident;
        // This is the read's IID
      CIFragT  * frag = GetCIFragT (ScaffoldGraph -> CIFrags,
                            fragID);

      assert (ident == frag -> iid);
      fprintf (fp,
               "    %7" F_CIDP "/%7" F_CIDP "/%c [%5.0f,%5.0f]",
               fragID, ident, (char) frag -> type,
               frag -> offset5p . mean, frag -> offset3p . mean);

      if  (frag -> numLinks == 1 && frag -> mateOf != NULLINDEX)
          {
           CIFragT  * mateFrag
               = GetCIFragT (ScaffoldGraph -> CIFrags, frag -> mateOf);
           ChunkInstanceT  * mateChunk
               = GetGraphNode (ScaffoldGraph -> CIGraph,
                     mateFrag -> CIid);
           DistT  * fragDist
               = GetDistT (ScaffoldGraph -> Dists, frag -> dist);
           assert (mateChunk != NULL);
           fprintf (fp,
                    " %7d/%6d [%5.0f,%5.0f] %6d [%6.0f,%6.0f] %6d %6.0f  %c %c",
                    frag -> mateOf,
#if  1
                    mateChunk -> id,
                    mateFrag -> offset5p . mean,
                    mateFrag -> offset3p . mean,
#endif
                    mateChunk -> info . CI . contigID,
                    mateFrag -> contigOffset5p . mean,
                    mateFrag -> contigOffset3p . mean,
                    mateChunk -> scaffoldID,
                    fragDist -> mean,
                    mateFrag -> flags . bits . hasInternalOnlyCILinks
                      ? 'T' : 'F',
                    mateFrag -> flags . bits . hasInternalOnlyContigLinks
                      ? 'T' : 'F');
          }
      fprintf (fp, "\n");
     }

   return;
  }




static void  Print_Potential_Fill_Chunks
    (FILE * fp, int (* has_potential) (ContigT *), int allow_bogus_edges)

//  Print list of all chunks that have edges to chunks already in
//  scaffolds.  Output goes to  fp .   has_potential ()  is a predicate
//  to screen entries to process.  Only chunks for which it is
//  true are printed.
//  If  allow_bogus_edges  is true, can use mate links that are marked
//   isProbablyBogusEdge () ; otherwise, ignore those edges.
//  If  fp  is NULL, do no actual printing.


  {
   GraphNodeIterator  contig_iterator;
   ContigT  * chunk;

   if  (fp != NULL)
       fprintf (fp, "\nConnections from unresolved chunks to uniques\n");

   assert (Num_Chunks == GetNumGraphNodes (ScaffoldGraph -> RezGraph));

   InitGraphNodeIterator (& contig_iterator, ScaffoldGraph -> RezGraph,
                          GRAPH_NODE_DEFAULT);
   while  ((chunk = NextGraphNodeIterator (& contig_iterator)) != NULL)
     {
      int  cid;

      cid = chunk -> id;

      if  (has_potential (chunk)
             && (! Is_Unique (chunk)
                  || (REF (cid) . is_singleton && UNIQUES_CAN_BE_STONES)))
          {  // show mate edges
           Stack_Entry_t  stack [STACK_SIZE];
           int  stack_top = 0;
           int  other_links = 0;
           GraphEdgeIterator  ci_edges;
           CIEdgeT  * edge;
           
           InitGraphEdgeIterator (ScaffoldGraph->RezGraph, cid,
                                ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, & ci_edges);
           while  ((edge = NextGraphEdgeIterator (& ci_edges)) != NULL)
             {
              ChunkInstanceT  * other_chunk;

              if  ((isProbablyBogusEdge (edge) && ! allow_bogus_edges)
                     || isSloppyEdge (edge))
                  continue;

              if  (edge -> idA == cid)
                  other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                   edge -> idB);
                else
                  other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                   edge -> idA);

              if  (Is_Unique (other_chunk))
                  {
                   assert (stack_top < STACK_SIZE);
		   assert(other_chunk->scaffoldID > NULLINDEX);
                   stack [stack_top] . chunk_id = other_chunk -> id;
                   stack [stack_top] . edge = edge;
                   stack_top ++;
                  }
                else
                  {
                   other_links += edge -> edgesContributing;
                   if  (isOverlapEdge (edge))
                       other_links --;
                  }
                  
             }

           if  (fp != NULL && stack_top > 0)
               {
                ChunkInstanceT  * ci = GetGraphNode (ScaffoldGraph -> CIGraph,
                                                     chunk -> info . Contig . AEndCI);
                int  i;

                fprintf (fp,
                    "Non-unique #%d (len = %.0f  frags = %d)"
#if  CHECK_CELSIM_COORDS
                    " [%7d,%7d] CGB_Type = %s"
#endif
                    " Links to non-unique = %d  %s\n",
                    cid, chunk -> bpLength . mean,
                    ci -> info . CI . numFragments,
#if  CHECK_CELSIM_COORDS
                    chunk -> aEndCoord, chunk -> bEndCoord,
                    CGB_Type_As_String (chunk -> flags . bits . cgbType)
#endif
                    other_links,
                    Is_Unique (chunk) ? "*UNIQUE*" : "");
                fprintf (fp,
                    "  %8s %6s %6s %8s %6s %6s %7s  %5s\n",
                    "To Chunk", "Scaff", "RelPos", "Gap Len", "Orient",
                    "Items", "Quality", "Flags");

		assert (stack_top < STACK_SIZE);
                for  (i = 0;  i < stack_top;  i ++)
                  {
                   fprintf (fp, "  %8d %6d %6d %8.0f %6s %6d %7.1f  ",
                        stack [i] . chunk_id,
                        REF (stack [i] . chunk_id) . scaff_id,
                        REF (stack [i] . chunk_id) . rel_pos,
                        stack [i] . edge -> distance . mean,
                        Orientation_As_String
                            (GetEdgeOrientationWRT (stack [i] . edge, cid)),
                        stack [i] . edge -> edgesContributing,
                        CIEdge_Quality (stack [i] . edge));
                   if  (isOverlapEdge (stack [i] . edge))
                       fprintf (fp, " Olap");
                   if  (stack [i] . edge -> flags . bits . isPossibleChimera)
                       fprintf (fp, " Chim");
                   if  (stack [i] . edge -> flags . bits . hasContributingOverlap)
                       fprintf (fp, " U-lap");
                   if  (stack [i] . edge -> flags . bits . hasTandemOverlap)
                       fprintf (fp, " T-lap");
                   if  (stack [i] . edge -> flags . bits . hasRepeatOverlap)
                       fprintf (fp, " r-lap");
                   if  (stack [i] . edge -> flags . bits . rangeTruncated)
                       fprintf (fp, " Trunc");
                   if  (isProbablyBogusEdge (stack [i] . edge))
                       fprintf (fp, " Bogus");
                   if  (isSloppyEdge (stack [i] . edge))
                       fprintf (fp, " Sloppy");
#if  VERBOSE
                   fprintf (fp, " numInstances=%d", ci -> info . CI . numInstances);
#endif
                   fprintf (fp, "\n");
                  }
#if  SHOW_FRAG_DETAILS
                Print_Frag_Info (fp, cid);
#endif
               }
          }
     }

   return;
  }



static void  Print_Unique_Chunks
    (FILE * fp)

//  Print listing of all unique chunks, sending output to  fp .
//  Allocate memory for global reference array  Ref  and fill in
//  chunk id's for unique chunks.
//  Also sets global  Num_Chunks  to the number of nodes in
//  the  RezGraph .  Do no actual printing if  fp  is  NULL .

  {
   int  i, cid, unique_ct;
   GraphNodeIterator  chunk_iterator;
   ChunkInstanceT  * chunk;
   int  ref_entry_ct;
   
   Num_Chunks = GetNumGraphNodes (ScaffoldGraph -> RezGraph);

   if  (fp != NULL)
       fprintf (fp, "### Contigs in graph = %d\n", Num_Chunks);

PALLOC (Num_Chunks * sizeof (Chunk_Ref_t));
#if  0
   Ref = (Chunk_Ref_t *) safe_calloc
           (Num_Chunks, sizeof (Chunk_Ref_t));
   for  (i = 0;  i < Num_Chunks;  i ++)
     REF (i) . scaff_id = REF (i) . rel_pos = NULLINDEX;
#else
   Ref_Index = (int *) safe_calloc
                   (Num_Chunks, sizeof (int));
   for  (i = 0;  i < Num_Chunks;  i ++)
     Ref_Index [i] = 0;
#endif

   assert (INITIAL_REF_DATA_SIZE >= 1);
   Ref_Data_Size = INITIAL_REF_DATA_SIZE;
   Ref_Data = (Chunk_Ref_t *) safe_malloc (Ref_Data_Size * sizeof (Chunk_Ref_t));
   Ref_Data [0] . scaff_id = Ref_Data [0] . rel_pos = NULLINDEX;
   Ref_Data [0] . is_singleton = Ref_Data [0] . is_unthrowable = FALSE;
        // Use 0 entry as dummy for non-scaffolded contigs
   ref_entry_ct = 1;


PALLOC (Num_Chunks * sizeof (Chunk_Info_t));
   Chunk_Info = (Chunk_Info_t *) safe_calloc
                  (Num_Chunks, sizeof (Chunk_Info_t));

fprintf (stderr, "Size of Ref = " F_SIZE_T "\n",
         Num_Chunks * sizeof (Chunk_Ref_t));
fprintf (stderr, "Size of Chunk_Info = " F_SIZE_T "\n",
         Num_Chunks * sizeof (Chunk_Info_t));

   if  (fp != NULL)
       {
        fprintf (fp, "\nUnique Chunks\n");
        fprintf (fp, "%6s %8s %8s %8s (Celsim coords)\n",
                 "Chunk#", "Scaff#",
                 "A End", "B End");
       }

   InitGraphNodeIterator(& chunk_iterator, ScaffoldGraph -> RezGraph,
                         GRAPH_NODE_DEFAULT);

   unique_ct = 0;
   while  ((chunk = NextGraphNodeIterator (& chunk_iterator)) != NULL)
     {
      cid = chunk -> id;
      Chunk_Info [cid] . cgb_type = chunk -> flags . bits . cgbType;
      if  (Is_Unique (chunk))
          {
           CIScaffoldT  * scaff;

           assert (chunk -> scaffoldID > NULLINDEX);

           if  (ref_entry_ct >= Ref_Data_Size)
               {
                Ref_Data_Size += INITIAL_REF_DATA_SIZE;
                Ref_Data = (Chunk_Ref_t *) safe_realloc
                              (Ref_Data, Ref_Data_Size * sizeof (Chunk_Ref_t));
               }
           Ref_Index [chunk -> id] = ref_entry_ct;
           REF (chunk -> id) . scaff_id = chunk -> scaffoldID;
           REF (chunk -> id) . a_end = chunk -> offsetAEnd;
           REF (chunk -> id) . b_end = chunk -> offsetBEnd;
           REF (chunk -> id) . is_singleton = FALSE;
           REF (chunk -> id) . is_unthrowable = FALSE;
           ref_entry_ct ++;

           scaff = GetGraphNode (ScaffoldGraph -> ScaffoldGraph,
                                  chunk -> scaffoldID);
           if  (scaff -> info. Scaffold . numElements == 1
//                  && chunk -> info . Contig . AEndCI
//                       == chunk -> info . Contig . BEndCI
               )
               REF (cid) . is_singleton = TRUE;

           if  (chunk -> flags . bits . isContig)
               UpdateContigSimCoordinates (chunk);

           if  (fp != NULL)
               fprintf (fp, "%6d %8" F_CIDP " %8" F_COORDP " %8" F_COORDP "\n",
                        cid, chunk -> scaffoldID,
                        chunk -> aEndCoord, chunk -> bEndCoord);
           unique_ct ++;
          }
     }

   if  (fp != NULL)
       fprintf (fp, "\n### total unique chunks = %d\n", unique_ct);
   fprintf (stderr, "### actual Ref entries = %d  size= " F_SIZE_T "\n",
            ref_entry_ct, ref_entry_ct * sizeof (Chunk_Ref_t));

   return;
  }



static int  Prior_Olaps_OK
    (int  from, int to, int to_low, Path_Info_t path_info [], int edge [],
     Stone_Edge_t pool [])

//  Check edges in  pool []  to confirm that transitive overlaps that
//  should be in  path_info  really are.  Overlaps are on the path
//  up to  from  and overlap  to .  Return  TRUE  if overlaps are all OK;
//  FALSE , otherwise.

  {
   int  i;

   for  (i = path_info [from] . from;  i >= 0;  i = path_info [i] . from)
     {
      if  (path_info [i] . hi_position >= to_low + CGW_DP_MINLEN)
          {
           int  j;

           for  (j = edge [i];  j >= 0;  j = pool [j] . next)
             if  (pool [j] . to == to)
                 break;

           if  (j < 0)
               return  FALSE;
          }
        else
          break;
     }

   return  TRUE;
  }



static void  Re_Check_Inserted_Rocks
    (Scaffold_Fill_t * fill, int min_good_links)

//  Retest links of entries marked  keep  in  fill  to see
//  if they still meet the link-consistency criteria to be rocks,
//  including the links to newly inserted rocks.  If an entry
//  fails, it is removed from its scaffold, and its  keep
//  flag is set  FALSE .
//  min_good_links  is the minimum number of good links a rock
//  needs to be kept.

  {
   GraphNodeIterator  chunk_iterator;
   ChunkInstanceT  * chunk;
   int  scaff_id, j;

   // First need to update  REF  info

   InitGraphNodeIterator(& chunk_iterator, ScaffoldGraph -> RezGraph,
                         GRAPH_NODE_DEFAULT);

   while  ((chunk = NextGraphNodeIterator (& chunk_iterator)) != NULL)
     {
      if  (Is_Unique (chunk))
          {
           assert (chunk -> scaffoldID > NULLINDEX);

           REF (chunk -> id) . scaff_id = chunk -> scaffoldID;
           REF (chunk -> id) . a_end = chunk -> offsetAEnd;
           REF (chunk -> id) . b_end = chunk -> offsetBEnd;
          }
     }


   // Now check all entries in  fill  that are marked  Keep

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
       {
        Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
        int  k;

        for  (k = 0;  k < this_gap -> num_chunks;  k ++)
          {
           Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

           if  (this_chunk -> keep
//                  && this_chunk -> copy_letter != JOINER_ROCK_CHAR
               )
               {
                Stack_Entry_t  stack [STACK_SIZE];
                int  cid;
                int  stack_top = 0;
                GraphEdgeIterator  ci_edges;
                CIEdgeT  * edge;

                cid = this_chunk -> chunk_id;

                InitGraphEdgeIterator
                    (ScaffoldGraph -> RezGraph, cid,
                     ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, & ci_edges);
                while  ((edge = NextGraphEdgeIterator (& ci_edges)) != NULL)
                  {
                   ChunkInstanceT  * other_chunk;

                   if  (isProbablyBogusEdge (edge)
                          || isSloppyEdge (edge))
                       continue;

                   if  (edge -> idA == cid)
                       other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                        edge -> idB);
                     else
                       other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                        edge -> idA);

                   if  (Is_Unique (other_chunk))
                       {
                        assert (stack_top < STACK_SIZE);
                        assert (other_chunk -> scaffoldID > NULLINDEX);
                        stack [stack_top] . chunk_id = other_chunk -> id;
                        stack [stack_top] . edge = edge;
                        stack_top ++;
                       }
                  }

                if  (stack_top == 0)
                    fprintf (stderr, "WHOOAH! no links for rock %d\n",
                             cid);
                  else
                    {
                     ContigT  * contig;
                     LengthT  left_end, right_end;
                     float  edge_quality;
                     int  bad_allowed, gap;
                     int  consistent, good_total, bad_links;

                     contig
                         = GetGraphNode
                               (ScaffoldGraph -> RezGraph, cid);
                     stack_top = Select_Good_Edges (stack, stack_top, contig);
                     consistent = Check_Scaffold_and_Orientation
                                      (cid, stack, stack_top, & good_total,
                                       fill, & bad_links, min_good_links);
                     bad_allowed = Max_int (0, 1 - bad_links);
                     if  (consistent && good_total >= min_good_links)
                         consistent = Estimate_Chunk_Ends
                                          (stack, stack_top, & left_end,
                                           & right_end, contig, & edge_quality,
                                           fill, & gap, & scaff_id,
                                           & bad_allowed);

                     if  (! consistent || good_total < min_good_links)
                         {
                          fprintf (stderr,
"### Upon further review, rock %d (copy = \'%c\') is removed from scaffold %d gap %d\n",
                                   cid, this_chunk -> copy_letter, scaff_id, j);

                          Remove_From_Scaffold (this_chunk);

                          this_chunk -> keep = FALSE;
                         }
                    }
               }
          }
       }

   return;
  }



static void  Reject_Non_Reachable
    (int id, Gap_Chunk_t * node [], int n, int edge [], Stone_Edge_t pool [])

//  Remove  keep  flag from every node in  node [0 .. (n - 1)]  that is
//  not reachable from the node with id  id  using the edges in  edge
//  that are linked through  pool [] .  Use only nodes that are marked  keep .

  {
   int  bad_edge;
   int  i;

   for  (i = 0;  i < n;  i ++)
     node [i] -> visited = node [i] -> finished = FALSE;

   for  (i = 0;  node [i] -> chunk_id != id && i < n;  i ++)
     ;

   if  (i >= n)
       {
        fprintf (stderr, "ERROR:  Node %d not in list\n", id);
        exit (EXIT_FAILURE);
       }

   node [i] -> keep = TRUE;          // Automatically keep the start node

   DFS_Stone_Visit (i, node, edge, pool, & bad_edge);

   for  (i = 0;  i < n;  i ++)
     if  (! node [i] -> visited)
         node [i] -> keep = FALSE;

   return;
  }



static void  Remove_From_Scaffold
    (Gap_Chunk_t * cp)

// See if the contig pointed to by  cp  is currently in a scaffold.
// If so, remove it.  If it was the last chunk in the scaffold,
// delete the scaffold, too.

  {
   ContigT  * chunk = GetGraphNode (ScaffoldGraph -> RezGraph, cp -> chunk_id);

   if  (chunk -> scaffoldID != NULLINDEX)
       {
        CIScaffoldT  * scaff;

        scaff = GetGraphNode (ScaffoldGraph -> ScaffoldGraph,
                              chunk -> scaffoldID);
        fprintf (stderr,
                 "Removing chunk %d from scaff %d\n",
                 cp -> chunk_id, chunk -> scaffoldID);
        RemoveCIFromScaffold (ScaffoldGraph, scaff, chunk, TRUE);
        if  (scaff -> info . Scaffold . numElements == 0)
            {
             fprintf (stderr, "Deleting scaffold %d\n", scaff -> id);
             DeleteGraphNode (ScaffoldGraph -> ScaffoldGraph, scaff);
            }

        if  (cp -> split)
            {
             fprintf (stderr,
                 "ERROR:  chunk %d was in scaffold and had split = TRUE\n",
                 cp -> chunk_id);
             cp -> split = FALSE;
            }
       }

   return;
  }



static int  Repeat_Colour
    (unsigned int cgb_type)

//  Return the celamy colour number for mark for  cgb_type .

  {
   switch  (cgb_type)
     {
      case  RU_CGBTYPE :
      case  RR_CGBTYPE :
        return  RU_RR_COLOUR;
      case  UR_CGBTYPE :
        return  UR_COLOUR;
     }

   return  0;
  }



static void  Requalify_Scaff_Chunks
    (Scaffold_Fill_t * fill_chunks)

//  Set the  keep  flag to  TRUE  for any scaffold chunks
//  in  fill_chunks .

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
        {
         int  k;
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (REF (this_chunk -> chunk_id) . scaff_id != NULLINDEX)
                this_chunk -> keep = TRUE;
           }
        }
     }

   return;
  }



static void  Restore_Best_Rocks
    (Scaffold_Fill_t * fill_chunks)

//  Set the  keep  flag true for the rock with the  best  flag set
//  in each gap in  fill_chunks  that has no chunks with the
//   keep  flag already set.

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = GAPS_TO_ADJUST)
        {
         int  k, keep_any;
         Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
         Gap_Chunk_t  * best_chunk = NULL;

         keep_any = FALSE;
         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                {
                 keep_any = TRUE;
                 this_chunk -> path_confirmed = TRUE;
                }
            else if  (this_chunk -> best)
                best_chunk = this_chunk;
           }

         if  (! keep_any && best_chunk != NULL
                && (ALLOW_LOOSE_END_ROCKS
                      || (this_gap -> left_cid >= 0
                           && this_gap -> right_cid >= 0)
                      || best_chunk -> copy_letter == JOINER_ROCK_CHAR))
             {
              ChunkOrientationType  orient;
              int  min_ahang, max_ahang;
#if 0
              int  how_much;
#endif
              int  had_left_overlap = FALSE;
//              ChunkOverlapCheckT  olap;
              char  * scaff_seq, * rock_seq;
              double  delta, implied_olap;
              double  left_hi_mean = 0., right_lo_mean = 0.;
              LengthT  * best_lo_pos, * best_hi_pos;
              double  max_gap_expansion, max_left_shift = DBL_MAX;
              double  gap_adjustment = 0.0;
              Overlap  * olap;

              best_chunk -> keep = TRUE;
              best_chunk -> path_confirmed = FALSE;
              if  (best_chunk -> flipped)
                  {
                   best_lo_pos = & (best_chunk -> end);
                   best_hi_pos = & (best_chunk -> start);
                  }
                else
                  {
                   best_lo_pos = & (best_chunk -> start);
                   best_hi_pos = & (best_chunk -> end);
                  }
              max_gap_expansion
                  = 3.0 * sqrt (fabs (this_gap -> end . variance
                                        - this_gap -> start . variance));

              rock_seq = Get_Contig_Sequence (best_chunk -> chunk_id);

              if  (this_gap -> left_cid >= 0)
                  {
                   left_hi_mean
                       = Max_double
                             (REF (this_gap -> left_cid) . a_end . mean,
                              REF (this_gap -> left_cid) . b_end . mean);
                   max_left_shift = best_lo_pos -> mean - left_hi_mean
                                      + CGW_MISSED_OVERLAP;
                  }
                else
                  max_left_shift = MAX_MATE_DISTANCE;
              if  (this_gap -> left_cid >= 0
                     && Might_Overlap
                          (REF (this_gap -> left_cid) . a_end . mean,
                           REF (this_gap -> left_cid) . b_end . mean,
                           best_chunk -> start . mean,
                           best_chunk -> end . mean,
                           5.0 * sqrt (best_chunk -> start . variance),
                           & orient, & min_ahang, & max_ahang))
                  {
#if  0
                   how_much
                     = (int) rint (Min_double
                                    (fabs (REF (this_gap -> left_cid) . a_end . mean
                                            - REF (this_gap -> left_cid) . b_end . mean),
                                     fabs (best_chunk -> start . mean
                                            - best_chunk -> end . mean)));
                                           
                   olap = OverlapChunks      // debug code does NOT handle suspicious
                              (ScaffoldGraph -> RezGraph,
                               this_gap -> left_cid, best_chunk -> chunk_id,
                               orient, 0, how_much,
                               CGW_DP_ERATE, FALSE);
#else
//  Replace with call to Mike's function.

                   scaff_seq = Get_Contig_Sequence (this_gap -> left_cid);
                   olap = OverlapSequences
                              (scaff_seq, rock_seq, orient, min_ahang, max_ahang,
                               CGW_DP_ERATE, CGW_DP_THRESH, CGW_DP_MINLEN,
                               AS_FIND_ALIGN);
#if  VERBOSE
fprintf (stderr, "Restore_Best olapping rock %d with left scaff_contig %d\n"
                 "  min_ahang = %d  max_ahang = %d  orient = %s\n",
         best_chunk -> chunk_id, this_gap -> left_cid,
         min_ahang, max_ahang, Orientation_As_String (orient));
if  (olap == NULL)
    fprintf (stderr, "  Not found\n");
  else
    fprintf (stderr, "  Found  begpos = %d  endpos = %d  length = %d\n",
             olap -> begpos, olap -> endpos, olap -> length);
#endif
#endif

                   if  (olap != NULL)
                       {
                        if  (olap -> endpos <= 0)
                            best_chunk -> keep = FALSE;
                          else
                            {
                             best_lo_pos -> mean
                                 = left_hi_mean - olap -> length;
                             if  (olap -> begpos < 0)
                                 {
                                  fprintf (stderr,
                                           "WOW!! Rock %d contains scaffold chunk %d\n",
                                           best_chunk -> chunk_id,
                                           this_gap -> left_cid);
                                  best_lo_pos -> mean += olap -> begpos;
                                 }
                             best_hi_pos -> mean
                                 = left_hi_mean + olap -> endpos;
                             best_lo_pos -> variance
                                 = ComputeFudgeVariance (olap -> length);
                             best_hi_pos -> variance
                                 = ComputeFudgeVariance
                                     (olap -> length
                                        + best_chunk -> len);
                             this_gap -> ref_variance
                                 = this_gap -> start . variance;
                             max_left_shift = 0;
                             had_left_overlap = TRUE;
                            }
                       }
                     else
                       {
                        implied_olap
                            = left_hi_mean
                                - best_lo_pos -> mean;
                        if  (implied_olap > CGW_MISSED_OVERLAP)
                            {
                             delta = implied_olap - CGW_MISSED_OVERLAP;
                             best_lo_pos -> mean += delta;
                             best_hi_pos -> mean += delta;
                             max_left_shift = 0;
                            }
                          else
                            max_left_shift = CGW_MISSED_OVERLAP - implied_olap;
                       }

                   free (scaff_seq);
                  }

              if  (best_chunk -> keep
                     && this_gap -> right_cid >= 0
                     && Might_Overlap
                          (best_chunk -> start . mean,
                           best_chunk -> end . mean,
                           REF (this_gap -> right_cid) . a_end . mean,
                           REF (this_gap -> right_cid) . b_end . mean,
                           5.0 * sqrt (best_chunk -> start . variance),
                           & orient, & min_ahang, & max_ahang))
                  {
//  Using Mike's function.
                   scaff_seq = Get_Contig_Sequence (this_gap -> right_cid);
                   olap = OverlapSequences
                              (rock_seq, scaff_seq, orient, min_ahang, max_ahang,
                               CGW_DP_ERATE, CGW_DP_THRESH, CGW_DP_MINLEN,
                               AS_FIND_ALIGN);
#if  VERBOSE
fprintf (stderr, "Restore_Best olapping rock %d with right scaff_contig %d\n"
                 "  min_ahang = %d  max_ahang = %d  orient = %s\n",
         best_chunk -> chunk_id, this_gap -> right_cid,
         min_ahang, max_ahang, Orientation_As_String (orient));
if  (olap == NULL)
    fprintf (stderr, "  Not found\n");
  else
    fprintf (stderr, "  Found  begpos = %d  endpos = %d  length = %d\n",
             olap -> begpos, olap -> endpos, olap -> length);
#endif

                   right_lo_mean
                       = Min_double
                             (REF (this_gap -> right_cid) . a_end . mean,
                              REF (this_gap -> right_cid) . b_end . mean);

                   if  (olap != NULL)
                       {
                        if  (olap -> begpos <= 0)
                            best_chunk -> keep = FALSE;
                        else if  (had_left_overlap)
                            best_chunk -> keep = FALSE;  // Should have been
                                                         // path confirmed
                          else
                            {
                             delta = right_lo_mean - olap->begpos - best_lo_pos->mean;
                             if  (delta < - max_left_shift)
                                 {
                                  gap_adjustment = - max_left_shift - delta;
                                  delta = - max_left_shift;
                                 }
                               else
                                 gap_adjustment = 0.0;

                             if  (gap_adjustment > max_gap_expansion)
                                 best_chunk -> keep = FALSE;
                               else
                                 {
                                  best_lo_pos -> mean
                                      = right_lo_mean - olap -> begpos
                                          + gap_adjustment;
                                  best_hi_pos -> mean += delta;
                                  if  (olap -> endpos < 0)
                                      {
                                       fprintf (stderr,
                                       "WOW!! Rock %d contains scaffold chunk %d\n",
                                                best_chunk -> chunk_id,
                                                this_gap -> right_cid);
                                       best_hi_pos -> mean -= olap -> endpos;
                                      }

                                  this_gap -> ref_variance
                                      = this_gap -> start . variance;
                                  best_hi_pos -> variance
                                      = this_gap -> end . variance
                                          - this_gap -> start . variance
                                          + ComputeFudgeVariance (olap -> length);
                                  if  (best_hi_pos -> variance < 0.0)
                                      best_hi_pos -> variance = 0.0;
                                  best_lo_pos -> variance
                                      = best_hi_pos -> variance
                                          - ComputeFudgeVariance
                                              (best_chunk -> len);
                                  if  (best_lo_pos -> variance < 0.0)
                                      best_lo_pos -> variance = 0.0;
                                 }
                           }
                       }
                     else
                       {
                        implied_olap
                            = best_hi_pos -> mean
                                - right_lo_mean;

                        if  (implied_olap > CGW_MISSED_OVERLAP)
                            {
                             delta = implied_olap - CGW_MISSED_OVERLAP;
                             if  (delta > max_left_shift)
                                 {
                                  gap_adjustment = delta - max_left_shift;
                                  delta = max_left_shift;
                                 }
                               else
                                 gap_adjustment = 0.0;

                             if  (gap_adjustment > max_gap_expansion)
                                 best_chunk -> keep = FALSE;
                               else
                                 {
                                  best_lo_pos -> mean -= delta;
                                  best_hi_pos -> mean -= delta;
                                 }
                            }
#if  VERBOSE
fprintf (stderr,
         "Restore_Best  chunk = %d  no rightolap   implied_olap = %.0f\n"
         "  delta = %.0f  gap_adjustment = %.0f  max_gap_exp = %.0f\n",
         best_chunk -> chunk_id, implied_olap, delta, gap_adjustment,
         max_gap_expansion);
#endif
                       }
                   free (scaff_seq);
                  }

              if  (best_chunk -> keep)
                  {
                   if  (j == 0)
                       Reverse_Positions (this_gap);
                     else
                       {
                        double  hi_var;

                        if  ((best_chunk -> start . mean <= best_chunk -> end . mean
                                && best_chunk -> start . variance
                                     > best_chunk -> end . variance)
                             ||
                               (best_chunk -> start . mean > best_chunk -> end . mean
                                  && best_chunk -> start . variance
                                     < best_chunk -> end . variance))
                            {
                             double  save;

                             save = best_chunk -> start . variance;
                             best_chunk -> start . variance
                                 = best_chunk -> end . variance;
                             best_chunk -> end . variance = save;
                            }
                        hi_var = Max_double (best_chunk -> start . variance,
                                             best_chunk -> end . variance)
                                   + this_gap -> ref_variance;
                        if  (hi_var > this_gap -> end . variance)
                            this_gap -> adjustment . variance
                                += hi_var - this_gap -> end . variance + EPSILON;
                        this_gap -> adjustment . mean = gap_adjustment;
                       }
                  }
              free (rock_seq);
             }
        }
      CheckScaffoldGraphCache (ScaffoldGraph);
     }

   return;
  }



static ChunkOrientationType  Reverse_Orientation
    (ChunkOrientationType orient)

//  Return the flipped edge orientation of  orient .

  {
   switch  (orient)
     {
      case  AB_AB :
        return  BA_BA;
      case  AB_BA :
        return  AB_BA;
      case  BA_AB :
        return  BA_AB;
      case  BA_BA :
        return  AB_AB;
      default :
        fprintf (stderr, "ERROR:  Bad edge orientation\n");
        assert (FALSE);
     }

   return  AB_AB;
  }



static void  Reverse_Positions
    (Gap_Fill_t * this_gap)

//  Reverse the order of positions of chunks in this gap and
//  make the lowest position and variance both 0.  Set
//  this_gap -> adjustment  to the value that must be added
//  to the position of all chunks after this one in the scaffold.

  {
   double  lo_position = DBL_MAX;
   double  hi_variance = -1.0;
   int  i, found = 0;

   for  (i = 0;  i < this_gap -> num_chunks;  i ++)
     {
      Gap_Chunk_t  * this_chunk = this_gap -> chunk + i;

      if  (this_chunk -> keep)
          {
           if  (this_chunk -> start . mean < lo_position)
               lo_position = this_chunk -> start . mean;
           if  (this_chunk -> start . variance > hi_variance)
               hi_variance = this_chunk -> start . variance;
           if  (this_chunk -> end . mean < lo_position)
               lo_position = this_chunk -> end . mean;
           if  (this_chunk -> end . variance > hi_variance)
               hi_variance = this_chunk -> end . variance;
           found ++;
          }
     }

   if  (found == 0)
       this_gap -> adjustment . mean = this_gap -> adjustment . mean = 0.0;
     else
       {
        this_gap -> adjustment . mean = - lo_position;
        this_gap -> adjustment . variance = hi_variance;

        for  (i = 0;  i < this_gap -> num_chunks;  i ++)
          {
           Gap_Chunk_t  * this_chunk = this_gap -> chunk + i;

           if  (this_chunk -> keep)
               {
                this_chunk -> start . mean -= lo_position;
                this_chunk -> end . mean -= lo_position;
                this_chunk -> start . variance
                    = hi_variance - this_chunk -> start . variance;
                this_chunk -> end . variance
                    = hi_variance - this_chunk -> end . variance;
               }
          }
       }

   return;
  }



static int  Scaff_Join_Cmp
    (const void * a, const void * b)

//  Return the order of  a  and  b  as  (Scaff_Join_t *) 's
//  based on  scaff1  (primary key) and  scaff2  (secondary key)
//  fields of what they point to.  Used for  qsort .

  {
   Scaff_Join_t  * * x, * * y;

   x = (Scaff_Join_t * *) a;
   y = (Scaff_Join_t * *) b;

   if  ((* x) -> scaff1 == (* y) -> scaff1
          && (* x) -> scaff2 == (* y) -> scaff2)
       return  0;
   else if  ((* x) -> scaff1 < (* y) -> scaff1
               || ((* x) -> scaff1 == (* y) -> scaff1
                     && (* x) -> scaff2 < (* y) -> scaff2))
       return  -1;

   return  1;
  }



Scaffold_Fill_t *  Scan_Gaps
    (void)

//  Find all gaps in  ScaffoldGraph  and allocate and return a data
//  structure to store them.

  {
   GraphNodeIterator  scaffolds;
   CIScaffoldT  * scaff;
   Scaffold_Fill_t  * fill;
   int  i;

   Num_Scaffolds = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);

PALLOC (Num_Scaffolds * sizeof (Scaffold_Fill_t));
   fill = (Scaffold_Fill_t *) safe_calloc (Num_Scaffolds, sizeof (Scaffold_Fill_t));

   for  (i = 0;  i < Num_Scaffolds;  i ++)
     {
      fill [i] . num_gaps = 0;
      fill [i] . added_to_ct = 0;
     }

   InitGraphNodeIterator (& scaffolds, ScaffoldGraph -> ScaffoldGraph,
                          GRAPH_NODE_DEFAULT);
   while  ((scaff = NextGraphNodeIterator (& scaffolds)) != NULL)
     {
      CIScaffoldTIterator  scaff_iterator;
      ChunkInstanceT  * chunk;
      int  scaff_id;
      LengthT  prev_end;
      int  scaff_index, num_chunks;

      scaff_id = scaff -> id;
      fill [scaff_id] . scaff_id = scaff_id;

#if  TEST_HOPELESS_SCAFFS
      if  (Is_Hopeless_Scaff [scaff_id] & Hopeless_True_Mask)
          continue;
#endif

      num_chunks = scaff -> info.Scaffold.numElements;
      fill [scaff_id] . num_gaps = 1 + num_chunks;   // including ends
PALLOC ((1 + num_chunks) * sizeof (Gap_Fill_t));
      fill [scaff_id] . gap
          = (Gap_Fill_t *) safe_calloc (1 + num_chunks, sizeof (Gap_Fill_t));
      fill [scaff_id] . gap [0] . gap = 0;
      fill [scaff_id] . gap [0] . num_chunks = 0;
      fill [scaff_id] . gap [0] . chunk = NULL;

      fill [scaff_id] . gap [0] . start . mean = - MAX_MATE_DISTANCE;
      fill [scaff_id] . gap [0] . start . variance = 0.0;
      fill [scaff_id] . gap [0] . ref_variance = 0.0;
      fill [scaff_id] . gap [0] . has_path = FALSE;

      InitCIScaffoldTIterator (ScaffoldGraph, scaff, TRUE, FALSE,
                              & scaff_iterator);

      chunk = NextCIScaffoldTIterator (& scaff_iterator);
      assert(chunk->scaffoldID == scaff_id);
      if  (chunk -> offsetAEnd . mean <= chunk -> offsetBEnd . mean)
          {
           fill [scaff_id] . gap [0] . end
               = chunk -> offsetAEnd;
           prev_end = chunk -> offsetBEnd;
          }
        else
          {
           fill [scaff_id] . gap [0] . end
               = chunk -> offsetBEnd;
           prev_end = chunk -> offsetAEnd;
          }
      fill [scaff_id] . gap [0] . left_cid = -1;
      fill [scaff_id] . gap [0] . right_cid = chunk -> id;
      fill [scaff_id] . gap [0] . len = 0;

      scaff_index = 1;
      for  (;  (chunk = NextCIScaffoldTIterator (& scaff_iterator)) != NULL;
              scaff_index ++)
        {
         Gap_Fill_t  * g = fill [scaff_id] . gap + scaff_index;

	 assert(chunk->scaffoldID == scaff_id);
         g -> gap = scaff_index;
         g -> num_chunks = 0;
         g -> chunk = NULL;
         g -> ref_variance = prev_end . variance;
         g -> has_path = FALSE;

         g -> left_cid = g [-1] . right_cid;
         g -> right_cid = chunk -> id;

         if  (chunk -> offsetAEnd . mean <= chunk -> offsetBEnd . mean)
             {
              if  (prev_end . mean <= chunk -> offsetAEnd . mean)
                  {
                   g -> start = prev_end;
                   g -> end = chunk -> offsetAEnd;
                   g -> len = g -> end . mean - g -> start . mean;
                  }
                else
                  {
                   g -> start = chunk -> offsetAEnd;
                   g -> end = prev_end;
                   g -> len = g -> start . mean - g -> end . mean;
                  }
              prev_end = chunk -> offsetBEnd;
             }
           else
             {
              if  (prev_end . mean <= chunk -> offsetBEnd . mean)
                  {
                   g -> start = prev_end;
                   g -> end = chunk -> offsetBEnd;
                   g -> len = g -> end . mean - g -> start . mean;
                  }
                else
                  {
                   g -> start = chunk -> offsetBEnd;
                   g -> end = prev_end;
                   g -> len = g -> start . mean - g -> end . mean;
                  }
              prev_end = chunk -> offsetAEnd;
             }
        }

      fill [scaff_id] . gap [scaff_index] . start
          = prev_end;
      fill [scaff_id] . gap [scaff_index] . end . mean
          = prev_end . mean + MAX_MATE_DISTANCE;
      fill [scaff_id] . gap [scaff_index] . end . variance
          = prev_end . variance;
      fill [scaff_id] . gap [scaff_index] . ref_variance
          = prev_end . variance;
      fill [scaff_id] . gap [scaff_index] . has_path = FALSE;
      fill [scaff_id] . gap [scaff_index] . left_cid
          = fill [scaff_id] . gap [scaff_index - 1] . right_cid;
      fill [scaff_id] . gap [scaff_index] . right_cid
          = -1;
      fill [scaff_id] . gap [scaff_index] . len = 0;

#if  0
// No longer needed since we now use the  ref_variance  field
// in the gap structure to adjust variances to a local value.

      // Now shrink variances to what's induced by previous 3 chunks
      // so that large values at end of scaffolds aren't used.  Should be
      // OK since calculations are done based on local chunks.

      for  (scaff_index = num_chunks;  scaff_index >= 3;  scaff_index --)
        {
         double  smaller;
         Gap_Fill_t  * g = fill [scaff_id] . gap + scaff_index;

         if  ((g - 3) -> end . variance < (g - 3) -> start . variance)
             smaller = (g - 3) -> end . variance;
           else
             smaller = (g - 3) -> start . variance;
         g -> end . variance -= smaller;
         g -> start . variance -= smaller;
        }
#endif

     }


   return  fill;
  }



Scaffold_Fill_t *  Scan_Gaps_In_Scaffold (int target_scaffold)

//  Find all gaps in  ScaffoldGraph  and allocate and return a data
//  structure to store them.

  {
   GraphNodeIterator  scaffolds;
   CIScaffoldT  * scaff;
   Scaffold_Fill_t  * fill;

PALLOC (1 * sizeof (Scaffold_Fill_t));
   fill = (Scaffold_Fill_t *) safe_calloc (1, sizeof (Scaffold_Fill_t));

   fill -> num_gaps = 0;

   InitGraphNodeIterator (& scaffolds, ScaffoldGraph -> ScaffoldGraph,
                          GRAPH_NODE_DEFAULT);
   while  ((scaff = NextGraphNodeIterator (& scaffolds)) != NULL)
     {
      CIScaffoldTIterator  scaff_iterator;
      ChunkInstanceT  * chunk;
      int  scaff_id;
      LengthT  prev_end;
      int  scaff_index, num_chunks;

	  if (scaff->id != target_scaffold)
		break;
	  
      scaff_id = scaff -> id;
      fill -> scaff_id = scaff_id;

      num_chunks = scaff -> info.Scaffold.numElements;
      fill -> num_gaps = 1 + num_chunks;   // including ends
PALLOC ((1 + num_chunks) * sizeof (Gap_Fill_t));
      fill -> gap
          = (Gap_Fill_t *) safe_calloc (1 + num_chunks, sizeof (Gap_Fill_t));
      fill -> gap [0] . gap = 0;
      fill -> gap [0] . num_chunks = 0;
      fill -> gap [0] . chunk = NULL;

      fill -> gap [0] . start . mean = - MAX_MATE_DISTANCE;
      fill -> gap [0] . start . variance = 0.0;
      fill -> gap [0] . ref_variance = 0.0;
      fill -> gap [0] . has_path = FALSE;

      InitCIScaffoldTIterator (ScaffoldGraph, scaff, TRUE, FALSE,
                              & scaff_iterator);

      chunk = NextCIScaffoldTIterator (& scaff_iterator);
      assert(chunk->scaffoldID == scaff_id);
      if  (chunk -> offsetAEnd . mean <= chunk -> offsetBEnd . mean)
          {
           fill -> gap [0] . end
               = chunk -> offsetAEnd;
           prev_end = chunk -> offsetBEnd;
          }
        else
          {
           fill -> gap [0] . end
               = chunk -> offsetBEnd;
           prev_end = chunk -> offsetAEnd;
          }
      fill -> gap [0] . left_cid = -1;
      fill -> gap [0] . right_cid = chunk -> id;

      scaff_index = 1;
      for  (;  (chunk = NextCIScaffoldTIterator (& scaff_iterator)) != NULL;
              scaff_index ++)
        {
         Gap_Fill_t  * g = fill -> gap + scaff_index;

	 assert(chunk->scaffoldID == scaff_id);
         g -> gap = scaff_index;
         g -> num_chunks = 0;
         g -> chunk = NULL;
         g -> ref_variance = prev_end . variance;
         g -> has_path = FALSE;

         g -> left_cid = g [-1] . right_cid;
         g -> right_cid = chunk -> id;

         if  (chunk -> offsetAEnd . mean <= chunk -> offsetBEnd . mean)
             {
              if  (prev_end . mean <= chunk -> offsetAEnd . mean)
                  {
                   g -> start = prev_end;
                   g -> end = chunk -> offsetAEnd;
                  }
                else
                  {
                   g -> start = chunk -> offsetAEnd;
                   g -> end = prev_end;
                  }
              prev_end = chunk -> offsetBEnd;
             }
           else
             {
              if  (prev_end . mean <= chunk -> offsetBEnd . mean)
                  {
                   g -> start = prev_end;
                   g -> end = chunk -> offsetBEnd;
                  }
                else
                  {
                   g -> start = chunk -> offsetBEnd;
                   g -> end = prev_end;
                  }
              prev_end = chunk -> offsetAEnd;
             }
        }

      fill -> gap [scaff_index] . start
          = prev_end;
      fill -> gap [scaff_index] . end . mean
          = prev_end . mean + MAX_MATE_DISTANCE;
      fill -> gap [scaff_index] . end . variance
          = prev_end . variance;
      fill -> gap [scaff_index] . ref_variance
          = prev_end . variance;
      fill -> gap [scaff_index] . has_path = FALSE;
      fill -> gap [scaff_index] . left_cid
          = fill -> gap [scaff_index - 1] . right_cid;
      fill -> gap [scaff_index] . right_cid
          = -1;

#if  0
// No longer needed since we now use the  ref_variance  field
// in the gap structure to adjust variances to a local value.

      // Now shrink variances to what's induced by previous 3 chunks
      // so that large values at end of scaffolds aren't used.  Should be
      // OK since calculations are done based on local chunks.

      for  (scaff_index = num_chunks;  scaff_index >= 3;  scaff_index --)
        {
         double  smaller;
         Gap_Fill_t  * g = fill -> gap + scaff_index;

         if  ((g - 3) -> end . variance < (g - 3) -> start . variance)
             smaller = (g - 3) -> end . variance;
           else
             smaller = (g - 3) -> start . variance;
         g -> end . variance -= smaller;
         g -> start . variance -= smaller;
        }
#endif

     }


   return  fill;
  }



static int  Select_Good_Edges
    (Stack_Entry_t * stack, int stack_top, ChunkInstanceT * chunk)

//  Move the "good" edges in  stack [0 .. (stack_top - 1)]  to the
//  front of stack and return how many there are.  Edges go from the
//  chunk pointed to by  chunk .  Also calculate and store on the stack
//  the orientation and scaffold position information of the chunk
//  pointed to.
//  
//  "Good" means including at least one mate pair of fragments
//  Other bad flags like "probablyBogus" or "hasGuide" should not
//  have been put on the stack in the first place.

  {
   int  good_ct, i;

   good_ct = 0;
   assert (stack_top < STACK_SIZE);	   
   for  (i = 0;  i < stack_top;  i ++)
     {
      stack [i] . num_good_mates
          = stack [i] . edge -> edgesContributing;
      stack [i] . num_ulaps = stack [i] . num_tlaps
          = stack [i] . num_rlaps = 0 ;
      if  (stack [i] . edge -> flags . bits . isPossibleChimera)
          stack [i] . num_good_mates = 1;   // Discount the overlap but
                                            //   count the mate
        else
          {
           if  (stack [i] . edge -> flags . bits . hasContributingOverlap)
               stack [i] . num_ulaps = 1;
           if  (stack [i] . edge -> flags . bits . hasTandemOverlap)
               stack [i] . num_tlaps = 1;
           if  (stack [i] . edge -> flags . bits . hasRepeatOverlap)
               stack [i] . num_rlaps = 1;
           if  (isOverlapEdge (stack [i] . edge))
               stack [i] . num_good_mates --;
          }

      if  (stack [i] . num_good_mates > 0)
          {
           int  scaff_id;
           ChunkInstanceT  * scaffold_chunk;
               // Chunk in scaffold at other end of edge
           LengthT  * a_side, * b_side;
               // Of chunk in scaffold
           ChunkOrientationType  orientation;

           scaffold_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                  stack [i] . chunk_id);
           a_side = & (REF (stack [i] . chunk_id) . a_end);
           b_side = & (REF (stack [i] . chunk_id) . b_end);
           scaff_id = REF (stack [i] . chunk_id) . scaff_id;

           stack [i] . scaff_id = scaff_id;

           if  (stack [i] . edge -> idB == stack [i] . chunk_id)
               orientation = stack [i] . edge -> orient;
             else
               orientation = Reverse_Orientation (stack [i] . edge -> orient);

           // Find position in scaffold
           switch  (orientation)
             {
              case  AB_AB :
              case  BA_AB :
                if  (Scaffold_Flipped [scaff_id])
                    stack [i] . celsim_offset
                        = scaffold_chunk -> aEndCoord + a_side -> mean;
                  else
                    stack [i] . celsim_offset
                        = scaffold_chunk -> aEndCoord - a_side -> mean;
                if  (a_side -> mean
                       <= b_side -> mean)
                    {
                     stack [i] . right_end . mean
                         = a_side -> mean
                               - stack [i] . edge -> distance . mean;
                     stack [i] . right_end . variance
                         = a_side -> variance
                               + stack [i] . edge -> distance . variance;
                     stack [i] . left_end . mean
                         = stack [i] . right_end . mean
                               - chunk -> bpLength . mean;
                     stack [i] . left_end . variance
                         = stack [i] . right_end . variance
                               + chunk -> bpLength . variance;
                     stack [i] . flipped
                         = (orientation == BA_AB);
                     stack [i] . source_variance = a_side -> variance;
                     stack [i] . left_link = FALSE;
                    }
                  else
                    {
                     stack [i] . left_end . mean
                         = a_side -> mean
                               + stack [i] . edge -> distance . mean;
                     stack [i] . left_end . variance
                         = a_side -> variance
                               + stack [i] . edge -> distance . variance;
                     stack [i] . right_end . mean
                         = stack [i] . left_end . mean
                               + chunk -> bpLength . mean;
                     stack [i] . right_end . variance
                         = stack [i] . left_end . variance
                               + chunk -> bpLength . variance;
                     stack [i] . flipped
                         = (orientation == AB_AB);
                     stack [i] . source_variance = a_side -> variance;
                     stack [i] . left_link = TRUE;
                    }
                break;
              case  AB_BA :
              case  BA_BA :
                if  (Scaffold_Flipped [scaff_id])
                    stack [i] . celsim_offset
                        = scaffold_chunk -> bEndCoord + b_side -> mean;
                  else
                    stack [i] . celsim_offset
                        = scaffold_chunk -> bEndCoord - b_side -> mean;
                if  (a_side -> mean
                       <= b_side -> mean)
                    {
                     stack [i] . left_end . mean
                         = b_side -> mean
                               + stack [i] . edge -> distance . mean;
                     stack [i] . left_end . variance
                         = b_side -> variance
                               + stack [i] . edge -> distance . variance;
                     stack [i] . right_end . mean
                         = stack [i] . left_end . mean
                               + chunk -> bpLength . mean;
                     stack [i] . right_end . variance
                         = stack [i] . left_end . variance
                               + chunk -> bpLength . variance;
                     stack [i] . flipped
                         = (orientation == AB_BA);
                     stack [i] . source_variance = b_side -> variance;
                     stack [i] . left_link = TRUE;
                    }
                  else
                    {
                     stack [i] . right_end . mean
                         = b_side -> mean
                               - stack [i] . edge -> distance . mean;
                     stack [i] . right_end . variance
                         = b_side -> variance
                               + stack [i] . edge -> distance . variance;
                     stack [i] . left_end . mean
                         = stack [i] . right_end . mean
                               - chunk -> bpLength . mean;
                     stack [i] . left_end . variance
                         = stack [i] . right_end . variance
                               + chunk -> bpLength . variance;
                     stack [i] . flipped
                         = (orientation == BA_BA);
                     stack [i] . source_variance = b_side -> variance;
                     stack [i] . left_link = FALSE;
                    }
                break;
              case  XX_XX :
              default :
                fprintf (stderr, "ERROR:  Unexpected edge orientation, line %d\n",
                         __LINE__);
                assert (FALSE);
             }
           stack [i] . is_bad = FALSE;
	   assert (good_ct < STACK_SIZE);
           stack [good_ct ++] = stack [i];
          }
     }

   return  good_ct;
  }



#if  TEST_HOPELESS_SCAFFS
static void  Set_Is_Hopeless
    (Scaffold_Fill_t * fill)

//  Set global  Is_Hopeless_Scaff  entries based on contents of  fill .
//  Any scaffold that has no  keep  entries has its hopeless entry set true.
//  Scaffolds already hopeless are ignored.

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  keep_ct;
      int  j;

      if  (Is_Hopeless_Scaff [scaff_id] & Hopeless_True_Mask)
          continue;

      keep_ct = 0;

      for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
         int  k;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                keep_ct ++;
           }
        }

      if  (keep_ct == 0)
          Is_Hopeless_Scaff [scaff_id] |= Hopeless_True_Mask;
     }

   return;
  }
#endif



static int  Set_Longest_Path
    (int list [], int num_list, Gap_Chunk_t * node [], int num_nodes,
     int target_sub, int edge [], Stone_Edge_t pool [],
     double ref_position, double factor, LengthT * target_position)

//  Set  keep  flag true for the longest sequence of nodes in
//  list [0 .. (num_list - 1)]  starting with  list [0]  and going
//  to  list [target_sub] , unless  target_sub  is negative
//  in which case the path can go anywhere.  Use edges in
//  edge  and  pool .  Values in  list  are subscripts of entries
//  in  node [0 .. (num_nodes - 1)] .   ref_position  is the scaffold
//  position corresponding to the zero path position.   factor
//  is  +1.0  for forward paths and  -1.0  for reverse.
//  Set  (* target_position)  to the place where the end of the gap
//  should be, based on the path.

  {
   static const double  NUM_STD_DEVS = 5.0;
   static Path_Info_t  * path_info = NULL;
   static int  path_info_size = 0;
   int  longest_sub = 0, longest_path = 0;
   int  num_kept = 0, index;
   int  i, j;

   assert (num_list > 0);
   target_position -> mean = 0.0;

   if  (num_nodes > path_info_size)
       {
PRALLOC (num_nodes * sizeof (Path_Info_t));
        path_info = (Path_Info_t *) safe_realloc
                        (path_info, num_nodes * sizeof (Path_Info_t));
        path_info_size = num_nodes;
       }

   for  (i = 0;  i < num_nodes;  i ++)
     {
      path_info [i] . path_len = 0;
      path_info [i] . from = -1;
      path_info [i] . hi_position = 0;
      path_info [i] . lo_position = 0;
      path_info [i] . a_hang = 0;
      path_info [i] . total_olap = 0;
      path_info [i] . hit = FALSE;
     }
   path_info [list [0]] . hit = TRUE;

#if  VERBOSE
fprintf (stderr, "Set_Longest_Path:\n");
#endif
   for  (i = 0;  i < num_list;  i ++)
     for  (j = edge [list [i]];  j >= 0;  j = pool [j] . next)
       {
        assert (pool [j] . from == list [i]);

        if  (node [pool [j] . to] -> keep && ! pool [j] . bad)
            {
#if  SHOW_STONE_CONFIRM
             fprintf (stderr,
                      "From %4d (%5d) to %4d (%5d)  ahg = %6d  prog = %6d  qual = %6.4f",
                      pool [j] . from, node [pool [j] . from] -> chunk_id,
                      pool [j] . to, node [pool [j] . to] -> chunk_id,
                      pool [j] . a_hang, pool [j] . progress,
                      pool [j] . quality);
#endif

             if  (path_info [pool [j] . from] . hit
                    && path_info [pool [j] . from] . path_len
                         >= path_info [pool [j] . to] . path_len)
                 {
                  double  test_position, error_band, std_dev;
                  LengthT  desired_position;

                  test_position
                      = (path_info [pool [j] . from] . hi_position
                          + pool [j] . progress) * factor
                              + ref_position;
                  if  ((node [pool [j] . to] -> flipped && factor >= 0.0)
                        || (! node [pool [j] . to] -> flipped && factor < 0.0))
                      desired_position = node [pool [j] . to] -> start;
                    else
                      desired_position = node [pool [j] . to] -> end;
                  std_dev = sqrt (desired_position . variance)
                                      + sqrt (ComputeFudgeVariance
                                         (path_info [pool [j] . from] . total_olap
                                           + 2.0 * node [pool [j] . to] -> len
                                           - pool [j] . progress));
                  error_band = NUM_STD_DEVS * std_dev;
                  if  (test_position < desired_position . mean - error_band
                         || test_position > desired_position . mean + error_band)
                      {
#if  SHOW_STONE_CONFIRM
                       fprintf (stderr,
                                "  missed  %6.0f not in [%6.0f,%6.0f]  %4.2f stdvs\n",
                                test_position,
                                desired_position . mean - error_band,
                                desired_position . mean + error_band,
                                (test_position - desired_position . mean) / std_dev);
#endif
                       continue;
                      }

#if  SHOW_STONE_CONFIRM
                  fprintf (stderr,
                           "  %6.0f in [%6.0f,%6.0f]  %4.2f stdvs\n",
                           test_position,
                           desired_position . mean - error_band,
                           desired_position . mean + error_band,
                           (test_position - desired_position . mean) / std_dev);
#endif

                  assert (pool [j] . progress >= 0);
                  if  (! Prior_Olaps_OK (pool [j] . from, pool [j] . to,
                                         path_info [pool [j] . from] . hi_position
                                           - pool [j] . length,
                                         path_info, edge, pool))
                      {
#if  SHOW_STONE_CONFIRM
                       fprintf (stderr, "  failed prior olaps\n");
#endif
                       continue;
                      }

                  path_info [pool [j] . to] . hit = TRUE;
                  path_info [pool [j] . to] . path_len
                      = path_info [pool [j] . from] . path_len + 1;
                  path_info [pool [j] . to] . from = pool [j] . from;
                  path_info [pool [j] . to] . a_hang = pool [j] . a_hang;
                  path_info [pool [j] . to] . hi_position
                      = path_info [pool [j] . from] . hi_position
                          + pool [j] . progress;
                  path_info [pool [j] . to] . lo_position
                      = path_info [pool [j] . from] . hi_position
                          - pool [j] . length;
                  path_info [pool [j] . to] . total_olap
                      = path_info [pool [j] . from] . total_olap
                          + 2 * node [pool [j] . to] -> len
                          - pool [j] . progress;
                  if  (path_info [pool [j] . to] . path_len > longest_path ||
                       (path_info [pool [j] . to] . path_len == longest_path
                        && path_info [pool [j] . to] . hi_position
                        > path_info [longest_sub] . hi_position))
                      {
                       longest_path = path_info [pool [j] . to] . path_len;
                       longest_sub = pool [j] . to;
                      }
                 }
#if  SHOW_STONE_CONFIRM
             fprintf (stderr, "  path_len = %2d  from = %2d\n",
                      path_info [pool [j] . to] . path_len,
                      path_info [pool [j] . to] . from);
#endif
            }
       }

   if  (target_sub < 0)
       target_sub = longest_sub;

   for  (i = 0;  i < num_list;  i ++)
     node [list [i]] -> keep = FALSE;

#if  SHOW_STONE_CONFIRM
fprintf (stderr, "target_sub = %d  hit = %c\n",
         target_sub, path_info [target_sub] . hit ? 'T' : 'F');
#endif

   if  (! path_info [target_sub] . hit)
       return  0;

#if  SHOW_STONE_CONFIRM
fprintf (stderr, "Reverse path:\n");
#endif

   num_kept = 0;
   for  (i = target_sub;  i >= 0;  i = path_info [i] . from)
     {
      node [i] -> keep = TRUE;
      num_kept ++;
#if  SHOW_STONE_CONFIRM
      fprintf (stderr,
      "%3d (%5d)  hi_position = %6d  AEnd = (%6.0f,%6.0f)  BEnd = (%6.0f,%6.0f)\n",
               i, node [i] -> chunk_id, path_info [i] . hi_position,
               node [i] -> start . mean, sqrt (node [i] -> start . variance),
               node [i] -> end . mean, sqrt (node [i] -> end . variance));
#endif
     }

   if  (num_kept <= 1)
       return  num_kept;

//  Create list of nodes in order by low position and set
//  that value in  path_info  based on  a_hang  fields in edges.

   index = num_kept;
   for  (i = target_sub;  i >= 0;  i = path_info [i] . from)
     list [-- index] = i;
   assert (index == 0);

   for  (i = 1;  i < num_kept;  i ++)
     {
      int  save = list [i];
      int  last_neg_ahang = 0;
      int  missing_olap = FALSE;
      int  k = 0;

      for  (j = i - 1;  j >= 0;  j --)
        {
         int  m;

         // Find edge from  list [j]  to  list [i]
         for  (k = edge [list [j]];  k >= 0;  k = pool [k] . next)
           if  (pool [k] . to == save)
               break;
         if  (k < 0)
             {
              int  sum = 0;

              fprintf (stderr, "YIKES!  No edge from %d to %d\n",
                       node [list [j]] -> chunk_id, node [save] -> chunk_id);
              if  (j != i - 1)
                  {
                   fprintf (stderr,
                      "YOUCH! Inconsistent overlaps.  Leave stones alone.\n");
                   for  (m = 0;  m < num_kept;  m ++)
                     node [list [m]] -> keep = FALSE;
                   return  0;
                  }
              for  (m = j;  m >= 0 && list [m] != path_info [save] . from;
                      m --)
                sum += path_info [list [m]] . a_hang;
              assert (m >= 0);
              path_info [save] . a_hang -= sum;
              assert (path_info [save] . a_hang >= 0);
              missing_olap = TRUE;
              break;
             }

         // if  a_hang is positive or back to first chunk,  stop
         // and put  i  after  j  in  list; otherwise, keep going

         if  (j == 0 || pool [k] . a_hang >= 0)
             break;

         last_neg_ahang = pool [k] . a_hang;
         list [j + 1] = list [j];
        }

      // when stop, save best edge from  j  to  i  (now at  j + 1)
      // and  i  to previous  j + 1
      // this is like an insertion sort

      if  (! missing_olap)
          {
           list [j + 1] = save;
           path_info [list [j + 1]] . a_hang = pool [k] . a_hang;
           if  (j < i - 1)
               path_info [list [j + 2]] . a_hang = - last_neg_ahang;
          }
     }

#if  SHOW_STONE_CONFIRM
{
 int  i;

 fprintf (stderr, "Revised path:\n");
 for  (i = 0;  i < num_kept;  i ++)
   {
    fprintf (stderr,
             "%3d %3d (%5d)  a_hang = %5d  hi_position = %6d\n",
               i, list [i], node [list [i]] -> chunk_id,
               path_info [list [i]] . a_hang,
               path_info [list [i]] . hi_position);
   }
}
#endif

   //  Use edges and order to set  a_hang  in  path_info
   //  Use  a_hang  in  Adjust_Stone_Positions

   Adjust_Stone_Positions (list, num_kept, node, ref_position, factor, path_info,
                           target_sub, target_position);

   return  num_kept;
  }



static void  Set_Position_From_Left_Olap
    (int left_cid, Gap_Chunk_t * this_chunk, Overlap * olap)

//  Set  start  and  end  positions in  (* this_chunk)  based on
//  the overlap between it and the contig  left_cid  recorded in
//  olap .

  {
   LengthT  * scaff_hi, * scaff_lo, * lo, * hi;

#if  VERBOSE
   fprintf (stderr,
            "left_scaff_contig:  start = <%8.0f,%8.0f>  end = <%8.0f,%8.0f>  id = %d\n",
            REF (left_cid) . a_end . mean, REF (left_cid) . a_end . variance,
            REF (left_cid) . b_end . mean, REF (left_cid) . b_end . variance,
            left_cid);
  
   fprintf (stderr,
            "       this_chunk:  start = <%8.0f,%8.0f>  end = <%8.0f,%8.0f>  id = %d\n",
            this_chunk -> start . mean, this_chunk -> start . variance,
            this_chunk -> end . mean, this_chunk -> end . variance,
            this_chunk -> chunk_id);
   fprintf (stderr,
            "             olap:  begpos = %7d  endpos = %7d  length = %7d\n",
            olap -> begpos, olap -> endpos, olap -> length);
#endif
   
   if  (REF (left_cid) . a_end . mean <= REF (left_cid) . b_end . mean)
       {
        scaff_lo = & (REF (left_cid) . a_end);
        scaff_hi = & (REF (left_cid) . b_end);
       }
     else
       {
        scaff_lo = & (REF (left_cid) . b_end);
        scaff_hi = & (REF (left_cid) . a_end);
       }
   if  (this_chunk -> flipped)
       {
        lo = & (this_chunk -> end);
        hi = & (this_chunk -> start);
       }
     else
       {
        lo = & (this_chunk -> start);
        hi = & (this_chunk -> end);
       }

   lo -> mean = scaff_lo -> mean + olap -> begpos;
   hi -> mean = scaff_hi -> mean + olap -> endpos;

   lo -> variance = scaff_hi -> variance + ComputeFudgeVariance (olap -> length);
   hi -> variance = lo -> variance
                      + ComputeFudgeVariance (hi -> mean - lo -> mean);

#if  VERBOSE
   fprintf (stderr,
            "   new this_chunk:  start = <%8.0f,%8.0f>  end = <%8.0f,%8.0f>\n",
            this_chunk -> start . mean, this_chunk -> start . variance,
            this_chunk -> end . mean, this_chunk -> end . variance);
#endif
   
   return;
  }



static void  Set_Position_From_Right_Olap
    (Gap_Chunk_t * this_chunk, int right_cid, Overlap * olap)

//  Set  start  and  end  positions in  (* this_chunk)  based on
//  the overlap between it and the contig  right_cid  recorded in
//  olap .

  {
   LengthT  * scaff_hi, * scaff_lo, * lo, * hi;

#if  VERBOSE
   fprintf (stderr,
            "right_scaff_contig:  start = <%8.0f,%8.0f>  end = <%8.0f,%8.0f>  id = %d\n",
            REF (right_cid) . a_end . mean, REF (right_cid) . a_end . variance,
            REF (right_cid) . b_end . mean, REF (right_cid) . b_end . variance,
            right_cid);
  
   fprintf (stderr,
            "       this_chunk:  start = <%8.0f,%8.0f>  end = <%8.0f,%8.0f>  id = %d\n",
            this_chunk -> start . mean, this_chunk -> start . variance,
            this_chunk -> end . mean, this_chunk -> end . variance,
            this_chunk -> chunk_id);
   fprintf (stderr,
            "             olap:  begpos = %7d  endpos = %7d  length = %7d\n",
            olap -> begpos, olap -> endpos, olap -> length);
#endif
   
   if  (REF (right_cid) . a_end . mean <= REF (right_cid) . b_end . mean)
       {
        scaff_lo = & (REF (right_cid) . a_end);
        scaff_hi = & (REF (right_cid) . b_end);
       }
     else
       {
        scaff_lo = & (REF (right_cid) . b_end);
        scaff_hi = & (REF (right_cid) . a_end);
       }
   if  (this_chunk -> flipped)
       {
        lo = & (this_chunk -> end);
        hi = & (this_chunk -> start);
       }
     else
       {
        lo = & (this_chunk -> start);
        hi = & (this_chunk -> end);
       }

   lo -> mean = scaff_lo -> mean - olap -> begpos;
   hi -> mean = scaff_hi -> mean - olap -> endpos;

   lo -> variance = scaff_lo -> variance + ComputeFudgeVariance (olap -> length);
   hi -> variance = lo -> variance
                      + ComputeFudgeVariance (hi -> mean - lo -> mean);

#if  VERBOSE
   fprintf (stderr,
            "   new this_chunk:  start = <%8.0f,%8.0f>  end = <%8.0f,%8.0f>\n",
            this_chunk -> start . mean, this_chunk -> start . variance,
            this_chunk -> end . mean, this_chunk -> end . variance);
#endif
   
   return;
  }



static void  Set_Split_Flags
    (Scaffold_Fill_t * fill, Set_Split_t set)

//  Set the  split  bit  of all entries in  fill  whose  keep
//  bit is true.  Set as specified by  set .

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     Set_Split_Flags_One_Scaffold (fill, set, scaff_id);

   return;
  }



static void  Set_Split_Flags_One_Scaffold
    (Scaffold_Fill_t * fill, Set_Split_t set, int scaff_id)

//  Set the  split  bit  of all entries in  fill [scaff_id]  whose  keep
//  bit is true.  Set as specified by  set .

  {
   int  j;

   for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
     {
      Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
      int  k;

      for  (k = 0;  k < this_gap -> num_chunks;  k ++)
        {
         Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

         if  (this_chunk -> keep)
             {
              switch  (set)
                {
                 case  ALL_FALSE :
                   this_chunk -> split = FALSE;
                   break;
                 case  FALSE_IFF_SINGLETON :
                   {
                    MultiAlignT  * ma
                      = LoadMultiAlignTFromSequenceDB
                            (ScaffoldGraph -> sequenceDB,
                             this_chunk -> chunk_id,
                             ScaffoldGraph -> RezGraph -> type == CI_GRAPH);
                      //                            = GetMultiAlignInStore
                      //                                  (ScaffoldGraph -> RezGraph -> maStore,
                      //        this_chunk -> chunk_id);


                    this_chunk -> split = (GetNumIntMultiPoss (ma -> f_list) != 1);
                    break;
                   }
                 default :
                   fprintf (stderr, "ERROR:  Unexpected split type\n");
                   assert (FALSE);
                }

              // if this chunk is currently in a scaffold (i.e., a unique)
              // its split flag is automatically false

              if  (REF (this_chunk -> chunk_id) . is_singleton)
                  this_chunk -> split = FALSE;
             }
        }
     }

   return;
  }



static int  Should_Overlap
    (Placement_t * left, Placement_t * right, ChunkOrientationType * orient,
     double * how_much)

//  Check positions of  left  and  right  and return whether or not
//  they indicate an overlap.  Set  (* orient)  to the orientation of
//  the overlap (whether or not it exists) and set  (* how_much)  to
//  the number of bases the overlap should be.

  {
   double  a, b, c, d;

   a = Min_double (left -> A_end . mean, left -> B_end . mean);
   b = Max_double (left -> A_end . mean, left -> B_end . mean);
   c = Min_double (right -> A_end . mean, right -> B_end . mean);
   d = Max_double (right -> A_end . mean, right -> B_end . mean);

   assert (b <= d);

   if  (d < a || b < c)
     (* how_much) = 0.0;
   else
     (* how_much) = Min_double (b, d) - Max_double (a, c);

   if  (left -> flipped)
       {
        if  (right -> flipped)
            (* orient) = BA_BA;
          else
            (* orient) = BA_AB;
       }
     else
       {
        if  (right -> flipped)
            (* orient) = AB_BA;
          else
            (* orient) = AB_AB;
       }

   return  ((* how_much) >= MIN_OLAP_LEN);
  }



void  Show_Gap_Reads_One_Scaff
    (FILE * fp, Scaffold_Fill_t * fill_chunks, int scaff_id)

//  Show the unitigs and reads in the gaps in scaffold  fill_chunks [scaff_id] .
//  Output goes to  fp .

  {
   int  num_entries;
   int  j;

   num_entries = 0;
   for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
     num_entries += fill_chunks [scaff_id] . gap [j] . num_chunks;
   if  (num_entries == 0)
       return;

   fprintf (fp, "\nScaff %d  num_gaps= %d\n",
            scaff_id, fill_chunks [scaff_id] . num_gaps);

   for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
     {
      Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;
      ChunkInstanceT  * scaff_chunk;
      int  k;
      double ref_var;
      
      if  (j > 0)
          {
           scaff_chunk
               = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> left_cid);
           ref_var = Max_double (scaff_chunk -> offsetAEnd . variance,
                                 scaff_chunk -> offsetBEnd . variance);
          }
        else
          {
           scaff_chunk
               = GetGraphNode(ScaffoldGraph->RezGraph, this_gap -> right_cid);
           ref_var = Min_double (scaff_chunk -> offsetAEnd . variance,
                                 scaff_chunk -> offsetBEnd . variance);
          }

      fprintf (fp,
               "Gap %3d  %8.0f .. %-8.0f  len= %.0f  L/Rnbrs= %d %d\n",
               j,
               this_gap -> start . mean,
               this_gap -> end . mean,
               this_gap -> len,
               this_gap -> left_cid,
               this_gap -> right_cid);

      for  (k = 0;  k < this_gap -> num_chunks;  k ++)
        {
         Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;
         ChunkInstanceT  * contig;
         ChunkInstanceT  * chunk;

         if  (REF (this_chunk -> chunk_id) . scaff_id != NULLINDEX
                && this_chunk -> copy_letter != GAP_END_CHAR
                && REF (this_chunk -> chunk_id) . is_unthrowable)
             {
              fprintf (fp, "   Chunk %6d%c is unthrowable\n",
                       this_chunk -> chunk_id, this_chunk -> copy_letter);
              break;
             }
             

         contig = GetGraphNode
                    (ScaffoldGraph -> RezGraph, this_chunk -> chunk_id);
         if  (contig -> flags . bits . isDead)
             continue;

         chunk = GetGraphNode
                   (ScaffoldGraph -> CIGraph, contig -> info . Contig . AEndCI);
         fprintf (fp,
"Uni %6d %c  %7.0f .. %-7.0f  cov= %3d links= %2d len= %7.0f nfr= %3d\n",
                  this_chunk -> chunk_id,
                  this_chunk -> copy_letter,
                  this_chunk -> start . mean,
                  this_chunk -> end . mean,
                  this_chunk -> cover_stat,
                  this_chunk -> link_ct,
                  contig -> bpLength . mean,
                  chunk -> info . CI . numFragments);
         Show_Read_Info (fp, this_chunk -> chunk_id);
        }
     }

   return;
  }



static void  Show_Read_Info
    (FILE * fp, int cid)

//  Print to  fp  info about the reads in unitig  cid .

  {
   ChunkInstanceT  * contig
       = GetGraphNode (ScaffoldGraph -> RezGraph, cid);
   ChunkInstanceT  * chunk;
   MultiAlignT  * ma;
   int  cover_stat;
   int  i, chunk_id, num_frags;

   chunk_id = contig -> info . Contig . AEndCI;
   chunk = GetGraphNode (ScaffoldGraph -> CIGraph, chunk_id);
   cover_stat = GetCoverageStat (chunk);

   ma = LoadMultiAlignTFromSequenceDB
            (ScaffoldGraph -> sequenceDB, cid,
             ScaffoldGraph -> RezGraph -> type == CI_GRAPH);
   assert (ma != NULL);

   // cycle through fragments 
   num_frags = GetNumIntMultiPoss (ma -> f_list);

   fprintf (fp, " %7s %13s %13s %6s %4s %6s %6s\n",
            "ReadIID", "Place", "Mate: Place", "Contig", "nfr",
            "Scaff", "Dist");

   for  (i = 0;  i < num_frags;  i ++)
     {
      IntMultiPos  * mp = GetIntMultiPos (ma -> f_list, i);
      CDS_CID_t  fragID = (CDS_CID_t) mp -> sourceInt;
        // This is an internal-data-structure ID
      CDS_CID_t  ident = (CDS_CID_t) mp -> ident;
        // This is the read's IID
      CIFragT  * frag = GetCIFragT (ScaffoldGraph -> CIFrags,
                            fragID);

      fprintf (fp,
               " %7" F_CIDP " %6.0f %6.0f",
               ident, frag -> offset5p . mean, frag -> offset3p . mean);

      if  (frag -> numLinks == 1 && frag -> mateOf != NULLINDEX)
          {
           CIFragT  * mateFrag
               = GetCIFragT (ScaffoldGraph -> CIFrags, frag -> mateOf);
           ChunkInstanceT  * mateChunk
               = GetGraphNode (ScaffoldGraph -> CIGraph,
                     mateFrag -> CIid);
           DistT  * fragDist
               = GetDistT (ScaffoldGraph -> Dists, frag -> dist);
           assert (mateChunk != NULL);
           fprintf (fp, " %6.0f %6.0f %6d %4d %6d %6.0f",
                    mateFrag -> offset5p . mean,
                    mateFrag -> offset3p . mean,
                    mateChunk -> info . CI . contigID,
                    mateChunk -> info . CI . numFragments,
                    mateChunk -> scaffoldID,
                    fragDist -> mean);
          }
      fprintf (fp, "\n");
     }

   return;
  }




int  Show_Reads_In_Gaps
    (char * prefix)

//  Print unitigs and reads in them that could go in gaps of
//  current scaffolds.  Output goes to file <prefix>.gapreads

  {
   FILE  * fp;
   char  filename [1000];
   Scaffold_Fill_t  * fill_stones;
   clock_t  start_time, stop_time;
   time_t  now;
   int  i, scaff_id;
     
   Num_Scaffolds = GetNumGraphNodes (ScaffoldGraph -> ScaffoldGraph);
   if (Num_Scaffolds == 0)
     return 0;

   now = time (NULL);
   fprintf (stderr, "### Start Show_Reads_In_Gaps at %s\n",
            ctime (& now));
   start_time = clock ();

#if  TEST_HOPELESS_SCAFFS
   Hopeless_False_Mask = '\373';
   Hopeless_True_Mask = '\004';
   if  (Is_Hopeless_Scaff == NULL)
       Is_Hopeless_Scaff
           = (char *) safe_calloc (Num_Scaffolds, sizeof (char));
     else
       {
        Is_Hopeless_Scaff
            = (char *) safe_realloc (Is_Hopeless_Scaff,
                                     Num_Scaffolds * sizeof (char));
        for  (i = 0;  i < Num_Scaffolds;  i ++)
          Is_Hopeless_Scaff [i] &= Hopeless_False_Mask;
        for  (i = Is_Hopeless_Size;  i < Num_Scaffolds;  i ++)
          Is_Hopeless_Scaff [i] = '\0';
       }
   Is_Hopeless_Size = Num_Scaffolds;
#endif

   strcpy (filename, prefix);
   strcat (filename, ".gapreads");
   fp = file_open (filename, "w");

PALLOC (Num_Scaffolds * sizeof (int64));
   Scaffold_Start = (int64 *) safe_calloc
                      (Num_Scaffolds, sizeof (int64));
PALLOC (Num_Scaffolds * sizeof (int64));
   Scaffold_End = (int64 *) safe_calloc
                    (Num_Scaffolds, sizeof (int64));
PALLOC (Num_Scaffolds * sizeof (char));
   Scaffold_Flipped = (char *) safe_calloc
                        (Num_Scaffolds, sizeof (char));

   Scaff_Join = CreateVA_Scaff_Join_t (INITIAL_SCAFF_JOIN_SIZE);

   // Need these to allocate global arrays
   Print_Unique_Chunks (NULL);
   Print_Scaffolds (NULL);
   Print_Potential_Fill_Chunks (NULL, Just_True, TRUE);

   fill_stones = Scan_Gaps ();

   Choose_Stones (fill_stones, 1, -1e6, FALSE);

   Clear_Keep_Flags (fill_stones, 0);

   Kill_Duplicate_Stones (fill_stones);

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     Show_Gap_Reads_One_Scaff
         (fp, fill_stones, scaff_id);

   fclose (fp);

   Free_Fill_Array (fill_stones);
   Free_Global_Arrays ();

   now = time (NULL);
   fprintf (stderr, "### Finish Show_Reads_In_Gaps at %s\n",
            ctime (& now));
   stop_time = clock ();
   fprintf (stderr, "### cpu time = %.1f sec\n",
               (double) (stop_time - start_time) / CLOCKS_PER_SEC);

   return  0;
  }



static void  Sort_Insertions
    (Scaffold_Fill_t * fill_chunks,
     int (* cmp) (const void *, const void *))

//  Sort the entries in each gap of  fill_chunks  using the function
//  cmp .

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     Sort_Insertions_One_Scaffold (fill_chunks, cmp, scaff_id);

   return;
  }


static void  Sort_Insertions_One_Scaffold
    (Scaffold_Fill_t * fill_chunks,
     int (* cmp) (const void *, const void *), int scaff_id)

//  Sort the entries in each gap of  fill_chunks [scaff_id]  using
//  the function  cmp .

  {
   int  j;

   for  (j = 0;  j < fill_chunks [scaff_id] . num_gaps;  j ++)
     {
      Gap_Fill_t  * this_gap = fill_chunks [scaff_id] . gap + j;

      if  (this_gap -> num_chunks > 0)
          qsort (this_gap -> chunk, this_gap -> num_chunks,
                 sizeof (Gap_Chunk_t), cmp);
     }

   return;
  }



#if VERBOSE
void
Throw_Stones_Debug(char *msg,
                   FILE *log_file,
                   Scaffold_Fill_t *fill_stones,
                   int scaff_id,
                   int *total_entries,
                   int *keep_entries) {
  fprintf (log_file, "\n\n>>>> fill_stones %s <<<<\n", msg);
  *total_entries = *keep_entries = 0;
  Print_Fill_Info_One_Scaffold (log_file, fill_stones, scaff_id, total_entries, keep_entries);
  fflush (log_file);
}
#else
#define Throw_Stones_Debug(A,B,C,D,E,F) ;
#endif


int Throw_Stones
    (char * prefix, int level, int use_partial)

//  Assign unresolved contigs to scaffolds, considering any contig with
//  any mate link to a scaffold--we call these "stones".
//  Stones are inserted if there is an overlap path that confirms
//  the position estimate from the link(s).
//  prefix  is the prefix of the name to use for output files
//  level  indicates which steps to perform.
//  If  use_partial  is true, then stones do not need to span
//  the entire gap.

  {
   FILE  * log_file;
   char  filename [1000], iter_string [20];
   Scaffold_Fill_t  * fill_stones;
   static int  iteration = 0;
   clock_t  start_time, stop_time;
   time_t  now;
   int  inserted = 0, total_stones = 0;
   int  stones_last_chkpt = 0;
   int  i, scaff_id;

   int  splitscaffolds = 0;

   Num_Scaffolds = GetNumGraphNodes (ScaffoldGraph -> ScaffoldGraph);
   if (Num_Scaffolds == 0)
     return 0;

   Use_Partial_Stone_Paths = use_partial;
     
#if  TEST_HOPELESS_SCAFFS
   Hopeless_False_Mask = '\367';
   Hopeless_True_Mask = '\010';
   if  (Is_Hopeless_Scaff == NULL)
       Is_Hopeless_Scaff = (char *) safe_calloc (Num_Scaffolds, sizeof (char));
     else
       {
        Is_Hopeless_Scaff
            = (char *) safe_realloc (Is_Hopeless_Scaff,
                                     Num_Scaffolds * sizeof (char));
        for  (i = 0;  i < Num_Scaffolds;  i ++)
          Is_Hopeless_Scaff [i] &= Hopeless_False_Mask;
       }
   Is_Hopeless_Size = Num_Scaffolds;
#endif

   Filename_Prefix = prefix;

   fprintf (stderr, "### Throw_Stones iteration #%d\n", iteration);

   sprintf (iter_string, "%d", iteration ++);

   now = time (NULL);
   fprintf (stderr, "### Start Stones iteration #%s   %s\n",
            iter_string, ctime (& now));
   start_time = clock ();

   strcpy (filename, prefix);
   strcat (filename, ".stone.i");
   strcat (filename, iter_string);
   strcat (filename, ".log");
   log_file = file_open (filename, "w");

#if  MAKE_CAM_FILE
   strcpy (filename, prefix);
   strcat (filename, ".stone.i");
   strcat (filename, iter_string);
   strcat (filename, ".cam");
   Cam_File = file_open (filename, "w");

   for  (i = 0;  i < NUM_COLOURS;  i ++)
     fprintf (Cam_File, "%dREZ: %s\n", i, Colour_String [i]);
#if  SHOW_CALC_COORDS
   strcpy (filename, prefix);
   strcat (filename, ".calcstone.i");
   strcat (filename, iter_string);
   strcat (filename, ".cam");
   Calc_Cam_File = file_open (filename, "w");

   for  (i = 0;  i < NUM_COLOURS;  i ++)
     fprintf (Calc_Cam_File, "%dREZ: %s\n", i, Colour_String [i]);

#endif
#endif

PALLOC (Num_Scaffolds * sizeof (int64));
   Scaffold_Start = (int64 *) safe_calloc
                      (Num_Scaffolds, sizeof (int64));
PALLOC (Num_Scaffolds * sizeof (int64));
   Scaffold_End = (int64 *) safe_calloc
                    (Num_Scaffolds, sizeof (int64));
PALLOC (Num_Scaffolds * sizeof (char));
   Scaffold_Flipped = (char *) safe_calloc
                        (Num_Scaffolds, sizeof (char));

   Scaff_Join = CreateVA_Scaff_Join_t (INITIAL_SCAFF_JOIN_SIZE);

   Print_Unique_Chunks (log_file);

   Print_Scaffolds (log_file);

//   Print_Potential_Fill_Chunks (log_file, Maybe_Stone);
   Print_Potential_Fill_Chunks (log_file, Just_True, FALSE);

   fill_stones = Scan_Gaps ();

   Choose_Stones (fill_stones, 1, MIN_STONE_COVER_STAT, FALSE);

#if  0
   {
    int  passed_consistency_check;

    Add_Gap_Ends (fill_stones);
    passed_consistency_check
        = check_consistency (fill_stones, Num_Scaffolds, iteration - 1);
    Disqualify_Scaff_Chunks (fill_stones);

    fprintf (stderr, "      Passed consistency check: %7d\n",
             passed_consistency_check);
   }
#else
   Clear_Keep_Flags (fill_stones, 0);
#endif

   Add_Gap_Ends (fill_stones);


   //  The original, very old, version did each step "sequentially" --
   //  we would do New_Confirm_Stones() on all scaffolds, then
   //  Disqualify_Scaff_Chunks(), etc.  Sometime before version 1.1 in
   //  sourceforge that original version was replaced with one that
   //  does each scaffold separately, and lets us checkpoint.
   //
   //  Removing the original version allows us to cleanup (ha, ha)
   //  AS_CGW_main by removing a duplicate CleanupScaffolds() call,
   //  without danger that somebody would reenable the historical
   //  code.

   fprintf (stderr, "### starting_stone_scaffold = %d\n",
            GlobalData -> starting_stone_scaffold);

   for  (scaff_id = GlobalData -> starting_stone_scaffold;
           scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  total_entries, keep_entries;

      // skip any scaffold that has been thrown as a stone into
      // another scaffold
      if  (fill_stones [scaff_id] . added_to_ct > 0)
          continue;

      Throw_Stones_Debug("Before New_Confirm_Stones", log_file, fill_stones, scaff_id, &total_entries, &keep_entries);
      New_Confirm_Stones_One_Scaffold (log_file, fill_stones, TRUE, scaff_id);
      Throw_Stones_Debug("After New_Confirm_Stones", log_file, fill_stones, scaff_id, &total_entries, &keep_entries);

      Disqualify_Scaff_Chunks_One_Scaffold (fill_stones, scaff_id);

      Set_Split_Flags_One_Scaffold (fill_stones, FALSE_IFF_SINGLETON, scaff_id);

      Throw_Stones_Debug("Before Jiggle_Positions", log_file, fill_stones, scaff_id, &total_entries, &keep_entries);
      Jiggle_Positions_One_Scaffold (fill_stones, scaff_id);
      Throw_Stones_Debug("After Jiggle_Positions", log_file, fill_stones, scaff_id, &total_entries, &keep_entries);

      Adjust_By_Ref_Variance_One_Scaffold (fill_stones, scaff_id);
      Throw_Stones_Debug("After Adjust_By_Ref_Variance", log_file, fill_stones, scaff_id, &total_entries, &keep_entries);

      Sort_Insertions_One_Scaffold (fill_stones, By_Keep_And_Low_Position, scaff_id);

      Throw_Stones_Debug("Before Doublecheck_Positions", log_file, fill_stones, scaff_id, &total_entries, &keep_entries);
      Doublecheck_Positions_One_Scaffold (fill_stones, TRUE, scaff_id);

      if  (level == 1)
          {
           int  have_bad_ct = 0, no_bad_ct = 0;
           int  keep_ct = 0, reject_ct = 0;

           Check_Other_Links_One_Scaffold
               (fill_stones, scaff_id, & have_bad_ct, & no_bad_ct,
                & keep_ct, & reject_ct);
          }

#if  UNIQUES_CAN_BE_STONES
          Eliminate_Encumbered_Uniques_One_Scaffold (fill_stones, scaff_id);
#endif

      Kill_Duplicate_Stones_One_Scaffold (fill_stones, scaff_id);

      fprintf (log_file, "\n\n>>>> Stones <<<<\n");

      total_entries = keep_entries = 0;
      Print_Fill_Info_One_Scaffold(log_file, fill_stones, scaff_id, & total_entries, & keep_entries);

      if  (keep_entries > 0) {
        int  components0 = 0;
        int  components1 = 0;

#if 1
        components0 = IsScaffoldInternallyConnected(ScaffoldGraph,
                                                    GetGraphNode(ScaffoldGraph->ScaffoldGraph, scaff_id),
                                                    ALL_EDGES);
#endif

        //  XXXXX: Even though USE_MY_INSERT is not defined, we still
        //  use the "my insert" function here.  On human, this made no
        //  difference in scaffold lengths....

        inserted = Insert_Chunks_In_Graph_One_Scaffold(ScaffoldGraph,
                                                       fill_stones,
                                                       scaff_id,
                                                       STONES);

        fprintf (stderr, ">>> Threw %d stones into scaffold %d\n", inserted, scaff_id);

        //  Occasionally we insert stones and get bogus edges
        //  included.  So, we remark all the edges.

        if (inserted > 0) {
          NodeCGW_T  *scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, scaff_id);
          int         nc    = 0;

          //  Someone includes edges to other scaffolds in this scaffold, so we
          //  should re-mark all the edges.
          //
          MarkInternalEdgeStatus(ScaffoldGraph, scaff, PAIRWISECHI2THRESHOLD_CGW,
                                 SLOPPY_EDGE_VARIANCE_THRESHHOLD,
                                 TRUE, TRUE, 0, TRUE);

          CleanupAScaffold(ScaffoldGraph,
                           scaff,
                           FALSE, NULLINDEX, FALSE);

          Force_Increasing_Variances_One_Scaffold(scaff_id);

#if 1
          //  A bad idea to use ALL_TRUSTED_EDGES.  We destroy our nice big scaffolds.
          //
          components1 = CheckScaffoldConnectivityAndSplit(ScaffoldGraph, scaff_id, ALL_EDGES, FALSE);  //  last arg is verbose

          if (components1 > 1) {
            splitscaffolds++;

            fprintf(stderr, "Throw_Stones()-- Scaffold %d components: ALL_EDGES=%d (%d before stones); ALL_TRUSTED_EDGES=%d components.\n",
                    components1,
                    components0,
                    IsScaffoldInternallyConnected(ScaffoldGraph,
                                                  GetGraphNode(ScaffoldGraph->ScaffoldGraph, scaff_id),
                                                  ALL_TRUSTED_EDGES));
          }
#endif


          total_stones += inserted;
        }
      }  //  stones to throw

     if  (total_stones - stones_last_chkpt >= STONES_PER_CHECKPOINT)
         {
          fprintf (stderr,
                   "* Stones CleanupScaffolds through scaffold %d\n",
                   scaff_id);

          fprintf (stderr,
                   "Writing Stone Checkpoint after %d stones at scaffold %d\n",
                   total_stones, scaff_id);
          fprintf (GlobalData -> timefp,
                   "\n\nCheckpoint %d written during Stones CleanupScaffolds"
                   " after scaffold %d\n",
                   ScaffoldGraph -> checkPointIteration, scaff_id);
          CheckpointScaffoldGraph (ScaffoldGraph, 0);

          // Clear cache if too large
          CheckScaffoldGraphCache (ScaffoldGraph);

          stones_last_chkpt = total_stones;
         }
     }  //  Main loop

   fprintf (stderr, "             Actually inserted: %7d\n", total_stones);
   fprintf (stderr, "             Scaffolds split because they were disconnected: %d\n", splitscaffolds);

   fclose (log_file);

    //  But since Throw_Stones can now split scaffolds, we should be
    //  rebuilding scaffold edges.
   //
   if (splitscaffolds > 0)
     BuildNewScaffoldEdges(ScaffoldGraph, 0);


#if  MAKE_CAM_FILE
   Update_Colours (fill_stones);
   Output_Cam_Files (fill_stones);

   fclose (Cam_File);
#if  SHOW_CALC_COORDS
   fclose (Calc_Cam_File);
#endif
#endif

   strcpy (filename, prefix);
   strcat (filename, ".stone.i");
   strcat (filename, iter_string);
   strcat (filename, ".analysis");
   log_file = file_open (filename, "w");
   Analyze_Stone_Fill (log_file, fill_stones);
   fclose (log_file);



   Free_Fill_Array (fill_stones);
   Free_Global_Arrays ();

   now = time (NULL);
   fprintf (stderr, "### Finish Stones iteration #%s   %s\n",
            iter_string, ctime (& now));
   stop_time = clock ();
   fprintf (stderr, "### cpu time = %.1f sec\n",
               (double) (stop_time - start_time) / CLOCKS_PER_SEC);

   Use_Partial_Stone_Paths = FALSE;

   return inserted;
  }



static void  Topologically_Sort
    (Gap_Chunk_t * node [], int num_nodes, int start_sub, int edge [],
     Stone_Edge_t pool [], int sorted [], int * num_sorted, int sort_all)

//  Topologically sort the nodes in  node [0 .. (num_nodes - 1)]
//  that are marked  keep  starting with  start_sub  and using
//  non-  bad  edges in  edge and  pool .  Put result in  sorted
//  in order, and set  (* num_sorted)  to the number put there.
//  If  sort_all  is true, sort all nodes, otherwise, only sort those
//  reachable from  start_sub .

  {
   int  i, j;

   for  (i = 0;  i < num_nodes;  i ++)
     node [i] -> visited = node [i] -> finished = FALSE;

#if  SHOW_STONE_CONFIRM
{
 int  i;

 fprintf (stderr, "Start Topologically_Sort:  start_sub = %d\n",
          start_sub);
 for  (i = 0;  i < num_nodes;  i ++)
   fprintf (stderr, "%3d: %5d %s\n",
            i, node [i] -> chunk_id, node [i] -> keep ? "keep" : "");
}
#endif

   (* num_sorted) = 0;
   if  (sort_all)
       {
        for  (i = 0;  i < num_nodes;  i ++)
          if  (! node [i] -> visited)
              Top_Sort_Visit (i, node, edge, pool, sorted, num_sorted);
       }
     else
       Top_Sort_Visit (start_sub, node, edge, pool, sorted, num_sorted);

   // Reverse finish time order

   for  (i = 0, j = (* num_sorted) - 1;  i < j;  i ++, j --)
     {
      int  save = sorted [i];

      sorted [i] = sorted [j];
      sorted [j] = save;
     }

   return;
  }



static void  Top_Sort_Visit
    (int i, Gap_Chunk_t * node [], int edge [], Stone_Edge_t pool [],
     int sorted [], int * num_sorted)

//  Recursively visit nodes in  node [i] 's subtree using edges in
//  edge  and  pool .  Only visit nodes marked  keep .  As nodes
//  are finished add them to  sorted  and increment  (* num_sorted) .

  {
   int  j;

   if  (node [i] -> visited || ! (node [i] -> keep))
       return;

   node [i] -> visited = TRUE;
   for  (j = edge [i];  j >= 0;  j = pool [j] . next)
     {
      assert (i == pool [j] . from);

      if  (! pool [j] . bad)
          Top_Sort_Visit (pool [j] . to, node, edge, pool, sorted,
                          num_sorted);
     }
   node [i] -> finished = TRUE;
   sorted [(* num_sorted)] = i;
   (* num_sorted) ++;

   return;
  }



int  Toss_Contained_Stones
    (char * prefix, int level, int redo_index)

//  Find stones that are contained within a contig already in a scaffold.
//  This is a main entry point from the chunk-graph-walker module.
//  prefix  is the prefix of the name to use for output files
//  level  indicates which steps to perform.
//  Note:  These stones are inserted containing their fragments--not
//  as surrogates.  This allows any other links they might have to
//  pull in additional chunks and/or join scaffolds.  Because of
//  this, they are inserted only if there is just one location for them.
//  redo_index  is the number of previous calls to this routine without
//  renumbering or merging scaffolds

  {
   FILE  * log_file = NULL;
   char  filename [1000], iter_string [20];
   Scaffold_Fill_t  * fill_stones;
   static int  iteration = 0;
   clock_t  start_time, stop_time;
   time_t  now;
   int  i;
   int  inserted = 0;
     
   StartTimerT(&GlobalData->GapFillTimer);

   Num_Scaffolds = GetNumGraphNodes (ScaffoldGraph -> ScaffoldGraph);
   if (Num_Scaffolds == 0)
     return 0;

   Single_Fragment_Only = TRUE;
     
#if  TEST_HOPELESS_SCAFFS
   Hopeless_False_Mask = '\373';
   Hopeless_True_Mask = '\004';
   if  (Is_Hopeless_Scaff == NULL)
       Is_Hopeless_Scaff
           = (char *) safe_calloc (Num_Scaffolds, sizeof (char));
     else
       {
        Is_Hopeless_Scaff
            = (char *) safe_realloc (Is_Hopeless_Scaff,
                                     Num_Scaffolds * sizeof (char));
        if  (redo_index <= 0)
            for  (i = 0;  i < Num_Scaffolds;  i ++)
              Is_Hopeless_Scaff [i] &= Hopeless_False_Mask;
        for  (i = Is_Hopeless_Size;  i < Num_Scaffolds;  i ++)
          Is_Hopeless_Scaff [i] = '\0';
       }
   Is_Hopeless_Size = Num_Scaffolds;
#endif

   Filename_Prefix = prefix;

   fprintf (stderr, "### Toss_Contained_Stones iteration #%d\n", iteration);
   Contained_Only_Switch = TRUE;
   
   sprintf (iter_string, "%d", iteration ++);

   now = time (NULL);
   fprintf (stderr, "### Start Contained Stones iteration #%s   %s\n",
            iter_string, ctime (& now));
   start_time = clock ();

#if VERBOSE
   strcpy (filename, prefix);
   strcat (filename, ".cstones.i");
   strcat (filename, iter_string);
   strcat (filename, ".log");
   log_file = file_open (filename, "w");
#endif

#if  MAKE_CAM_FILE
   strcpy (filename, prefix);
   strcat (filename, ".cstones.i");
   strcat (filename, iter_string);
   strcat (filename, ".cam");
   Cam_File = file_open (filename, "w");

   for  (i = 0;  i < NUM_COLOURS;  i ++)
     fprintf (Cam_File, "%dREZ: %s\n", i, Colour_String [i]);
#if  SHOW_CALC_COORDS
   strcpy (filename, prefix);
   strcat (filename, ".calccs.i");
   strcat (filename, iter_string);
   strcat (filename, ".cam");
   Calc_Cam_File = file_open (filename, "w");

   for  (i = 0;  i < NUM_COLOURS;  i ++)
     fprintf (Calc_Cam_File, "%dREZ: %s\n", i, Colour_String [i]);
#endif
#endif

PALLOC (Num_Scaffolds * sizeof (int64));
   Scaffold_Start = (int64 *) safe_calloc
                      (Num_Scaffolds, sizeof (int64));
PALLOC (Num_Scaffolds * sizeof (int64));
   Scaffold_End = (int64 *) safe_calloc
                    (Num_Scaffolds, sizeof (int64));
PALLOC (Num_Scaffolds * sizeof (char));
   Scaffold_Flipped = (char *) safe_calloc
                        (Num_Scaffolds, sizeof (char));

   Scaff_Join = CreateVA_Scaff_Join_t (INITIAL_SCAFF_JOIN_SIZE);

   for  (i = 0;  i < Num_Scaffolds;  i ++)
     Scaffold_Start [i] = Scaffold_End [i] = -1;

   fprintf (stderr, ">>> Before  Print_Unique_Chunks\n");
   Print_Unique_Chunks (log_file);

   fprintf (stderr, ">>> Before  Print_Scaffolds\n");
   Print_Scaffolds (log_file);

   fprintf (stderr, ">>> Before  Print_Potential_Fill_Chunks\n");
   //   Print_Potential_Fill_Chunks (log_file, Maybe_Stone);
   Print_Potential_Fill_Chunks (log_file, Just_True, TRUE);

   StartTimerT(&GlobalData->ChooseChunksTimer);

   fprintf (stderr, ">>> Before  Scan_Gaps\n");
   fill_stones = Scan_Gaps ();

   fprintf (stderr, ">>> Before  Choose_Stones\n");
   Choose_Stones (fill_stones, 1, MIN_STONE_COVER_STAT, TRUE);

#if  VERBOSE
   fprintf (log_file, "\n>>> Fill after Choose_Stones <<<\n");
   Print_Fill_Info (log_file, fill_stones);
#endif

   StopTimerT(&GlobalData->ChooseChunksTimer);

   Add_Gap_Ends (fill_stones);

#if  VERBOSE
   fprintf (log_file, "\n>>> Fill after Add_Gap_Ends <<<\n");
   Print_Fill_Info (log_file, fill_stones);
#endif

   Confirm_Contained (log_file, fill_stones, TRUE);

   Disqualify_Scaff_Chunks (fill_stones);

#if  UNIQUES_CAN_BE_STONES
       Eliminate_Encumbered_Uniques (fill_stones);
#endif

   Verify_Single_Placement (fill_stones);

   Sort_Insertions (fill_stones, By_Keep_And_Low_Position);

#if  MAKE_CAM_FILE
   Update_Colours (fill_stones);
   Output_Cam_Files (fill_stones);

   fclose (Cam_File);
#if  SHOW_CALC_COORDS
   fclose (Calc_Cam_File);
#endif
#endif

   fprintf (stderr, "Set_Split_Flags\n");
   Set_Split_Flags (fill_stones, ALL_FALSE);
   fprintf (stderr, "After Set_Split_Flags\n");

#if  VERBOSE
   fprintf (log_file, "\n>>> Fill before  Update_Scaffold_Graph <<<\n");
   Print_Fill_Info (log_file, fill_stones);

   fclose (log_file);

   strcpy (filename, prefix);
   strcat (filename, ".cstones.i");
   strcat (filename, iter_string);
   strcat (filename, ".analysis");
   log_file = file_open (filename, "w");
   Analyze_Rock_Fill (log_file, fill_stones);

   fclose (log_file);
#endif


   StartTimerT(&GlobalData->UpdateTimer);
   if  (level > 1)
       {
#if  USE_MY_INSERT
        inserted = Insert_Chunks_In_Graph (ScaffoldGraph, fill_stones, STONES);
#else
        // Not using surrogates on contains.
        inserted
            = Update_Scaffold_Graph
                  (ScaffoldGraph, fill_stones, FALSE, TRUE,
                   FALSE, TRUE /* copyAllOverlaps...not used */, -1, STONES);
#endif
        fprintf (stderr, "             Actually inserted: %7d\n", inserted);
       }
   StopTimerT(&GlobalData->UpdateTimer);

   Contained_Only_Switch = FALSE;

#if  TEST_HOPELESS_SCAFFS
   Set_Is_Hopeless (fill_stones);
#endif

   Free_Fill_Array (fill_stones);
   Free_Global_Arrays ();

   now = time (NULL);
   fprintf (stderr, "### Finish Contained Stones iteration #%s   %s\n",
            iter_string, ctime (& now));
   stop_time = clock ();
   fprintf (stderr, "### cpu time = %.1f sec\n",
               (double) (stop_time - start_time) / CLOCKS_PER_SEC);

   StopTimerT(&GlobalData->GapFillTimer);

   Single_Fragment_Only = FALSE;

   return inserted;
  }



static void  Update_Colours
    (Scaffold_Fill_t * fill_chunks)

//  Update the  colour  field in global  Chunk_Info  based on information
//  returned in  fill_chunks .

  {
   int  scaff_id;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      Gap_Fill_t  * g;
      int  i;

      g = fill_chunks [scaff_id] . gap;

      for  (i = 0;  i < fill_chunks [scaff_id] . num_gaps;  i ++)
        {
         int  j, num_chunks;
         Gap_Chunk_t  * chunk;

         num_chunks = g [i] . num_chunks;
         chunk = g [i] . chunk;

         for  (j = 0;  j < num_chunks;  j ++)
           if  (! chunk [j] . keep
                  && Chunk_Info [chunk [j] . chunk_id] . colour != UNIQUE_COLOUR)
                Chunk_Info [chunk [j] . chunk_id] . colour = REJECT_COLOUR;
        }
     }

   return;
  }



static void   Verify_Single_Placement
    (Scaffold_Fill_t * fill)

//  Set the  keep  flag to false of any chunk in  fill  that has more than
//  one occurrence (with  keep  set true).

  {
   int  scaff_id;
   int  i;

   for  (i = 0;  i < Num_Chunks;  i ++)
     Chunk_Info [i] . use_ct = 0;

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
         int  k;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (this_chunk -> keep)
                Chunk_Info [this_chunk -> chunk_id] . use_ct ++;
           }
        }
     }

   for  (scaff_id = 0;  scaff_id < Num_Scaffolds;  scaff_id ++)
     {
      int  j;

      for  (j = 0;  j < fill [scaff_id] . num_gaps;  j ++)
        {
         Gap_Fill_t  * this_gap = fill [scaff_id] . gap + j;
         int  k;

         for  (k = 0;  k < this_gap -> num_chunks;  k ++)
           {
            Gap_Chunk_t  * this_chunk = this_gap -> chunk + k;

            if  (Chunk_Info [this_chunk -> chunk_id] . use_ct > 1)
                this_chunk -> keep = FALSE;
           }
        }
     }

   return;
  }



static int  Violates_Scaff_Edges
    (Scaff_Join_t  * p)

//  Return whether the scaffold join pointed to by  p  is inconsistent
//  with any scaffold edges in the scaffold graph.

  {
   CIScaffoldT  * scaffold1, * scaffold2;
   GraphEdgeIterator  SEdges;
   SEdgeT  * edge;
   LengthT  gap;
   int  found_supporting_edge, trust_edge;

   // Check if there is an edge consistent with p's info
   // If there is any trusted edge that is inconsistent with the join
   // or there is no edge (trusted or not) that is consistent with it
   // (but only if  REQUIRE_SCAFFOLD_EDGE_FOR_JOIN  is true),
   // then TRUE (i.e., violated) is returned.

   scaffold1 = GetCIScaffoldT (ScaffoldGraph -> CIScaffolds, p -> scaff1);
   scaffold2 = GetCIScaffoldT (ScaffoldGraph -> CIScaffolds, p -> scaff2);
   InitGraphEdgeIterator(ScaffoldGraph->ScaffoldGraph, p->scaff1, ALL_END, ALL_EDGES, 
			  GRAPH_EDGE_DEFAULT,   &SEdges);
      

   //   InitSEdgeTIterator (ScaffoldGraph, p -> scaff1, FALSE, FALSE,
   //                       ALL_END, FALSE,  & SEdges);

   found_supporting_edge = FALSE;
   while  ((edge = NextGraphEdgeIterator (& SEdges)) != NULL)
     {
      int  edge_lo, edge_hi, gap_lo, gap_hi;
      double  edge_delta, gap_delta;

      assert (edge -> idA < edge -> idB);
      trust_edge = Is_Good_Scaff_Edge (edge);

      if  (edge -> idA == p -> scaff1 && edge -> idB == p -> scaff2)
          {
           switch  (edge -> orient)
             {
              case  AB_AB :
                if  (p -> m != 1 && trust_edge)
                    return  TRUE;
                gap . mean = - (scaffold1 -> bpLength . mean + p -> b);
                gap . variance = scaffold1 -> bpLength . variance + p -> variance;
                break;
              case  AB_BA :
                if  (p -> m != -1 && trust_edge)
                    return  TRUE;
                gap . mean = p -> b - scaffold1 -> bpLength . mean
                               - scaffold2 -> bpLength . mean;
                gap . variance = scaffold1 -> bpLength . variance
                              + scaffold2 -> bpLength . variance
                              + p -> variance;
                break;
              case  BA_AB :
                if  (p -> m != -1 && trust_edge)
                    return  TRUE;
                gap . mean = - (p -> b);
                gap . variance = p -> variance;
                break;
              case  BA_BA :
                if  (p -> m != 1 && trust_edge)
                    return  TRUE;
                gap . mean = p -> b - scaffold2 -> bpLength . mean;
                gap . variance = scaffold2 -> bpLength . variance + p -> variance;
                break;
              default :
                return  TRUE;   // edge without orientation kills join
             }

           edge_delta = 3.0 * sqrt (edge -> distance . variance);
           edge_lo = (int) (edge -> distance . mean - edge_delta);
           edge_hi = (int) (edge -> distance . mean + edge_delta);
           gap_delta = 3.0 * sqrt (gap . variance);
           gap_lo = (int) (gap . mean - gap_delta);
           gap_hi = (int) (gap . mean + gap_delta);

           if  (Interval_Intersection (edge_lo, edge_hi, gap_lo, gap_hi))
               found_supporting_edge = TRUE;
           else if  (trust_edge)
               {
                return  TRUE;
               }
          }
      else if  (edge -> idA == p -> scaff2 && edge -> idB == p -> scaff1)
          {
           switch  (edge -> orient)
             {
              case  AB_AB :
                if  (p -> m != 1 && trust_edge)
                    return  TRUE;
                gap . mean = p -> b - scaffold2 -> bpLength . mean;
                gap . variance = scaffold2 -> bpLength . variance + p -> variance;
                break;
              case  AB_BA :
                if  (p -> m != -1 && trust_edge)
                    return  TRUE;
                gap . mean = p -> b - scaffold1 -> bpLength . mean
                               - scaffold2 -> bpLength . mean;
                gap . variance = scaffold1 -> bpLength . variance
                              + scaffold2 -> bpLength . variance
                              + p -> variance;
                break;
              case  BA_AB :
                if  (p -> m != -1 && trust_edge)
                    return  TRUE;
                gap . mean = - (p -> b);
                gap . variance = p -> variance;
                break;
              case  BA_BA :
                if  (p -> m != 1 && trust_edge)
                    return  TRUE;
                gap . mean = - (p -> b + scaffold1 -> bpLength . mean);
                gap . variance = scaffold1 -> bpLength . variance + p -> variance;
                break;
              default :
                return  TRUE;   // edge without orientation kills join
             }

           edge_delta = 3.0 * sqrt (edge -> distance . variance);
           edge_lo = (int) (edge -> distance . mean - edge_delta);
           edge_hi = (int) (edge -> distance . mean + edge_delta);
           gap_delta = 3.0 * sqrt (gap . variance);
           gap_lo = (int) (gap . mean - gap_delta);
           gap_hi = (int) (gap . mean + gap_delta);

           if  (Interval_Intersection (edge_lo, edge_hi, gap_lo, gap_hi))
               found_supporting_edge = TRUE;
           else if  (trust_edge)
               {
                return  TRUE;
               }
          }
     }

   if  (found_supporting_edge || ! REQUIRE_SCAFFOLD_EDGE_FOR_JOIN)
       return  FALSE;

   return  TRUE;
  }





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
/**********************************************************************
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/AS_ORA_main.c,v $
$Revision: 1.12 $
**********************************************************************/

/**********************************************************************
Module: AS_ORA (Assembler Overlap Detector Diagnostic Module)

Description: AS_ORA reads in simulated fragment messages and associated
             overlap messages (generaged by AS_OVL), performs a series
             correctness checks on the detected overlaps & generates
             several statistics

Assumptions: Fragment messages MUST be "truthed" - contain start & stop
             indices into the source sequence.

**********************************************************************/


/*********************************************************************/
// headers
// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include  <assert.h>
#include  <ctype.h>

// project headers
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_distStore.h"
#include "AS_ORA_fragments.h"
#include "AS_ORA_overlaps.h"
#include "AS_ORA_statistics.h"
#include "AS_ORA_inlines.h"
/*********************************************************************/


/*********************************************************************/
// defines
#ifndef MAX_SEQUENCE_LENGTH
#define MAX_SEQUENCE_LENGTH AS_READ_MAX_LEN
#endif

#ifndef MAX_SOURCE_LENGTH
#define MAX_SOURCE_LENGTH 512
#endif
/*********************************************************************/

static int  Global_Debug = FALSE;

/*********************************************************************/
// structures
/*********************************************************************/


/*********************************************************************/
// forward function declarations
int ComputeTrueOverlaps( FragmentArrayp fragments,
                         OverlapStatisticsp stats,
                         uint32 heap_overlaps );
int CompareOverlaps( FragmentArrayp fragments,
                     char * overlap_filename,
                     OverlapStatisticsp stats );
int AddOverlap( FragmentArrayp fragments, uint64 i1, uint64 i2 );
void InsertionSortFragmentOverlap( FragmentObjectp frag,
                                   OverlapObjectp overlap );
int SetOverlapBestMatch( FragmentArrayp fragments,
                         uint32 with_id,
                         OverlapObjectp foop,
                         OverlapMesg * osp );
void PopulateOverlapObject( OverlapObjectp oop,
                            OverlapMesg * osp,
                            uint32 a_id,
                            int32 hds );
int GenerateOverlapFile( char * output_ovl_filename,
                         char * fragstore_name,
                         FragmentArrayp fragments );
int32 MinHangDeltaSum( FragmentArrayp fragments,
                       OverlapObjectp oop,
                       uint32 a_id,
                       OverlapMesg * osp );
int ComputeGaps( FragmentArrayp fragments,
                 OverlapStatisticsp stats );

// New stuff added by Art.
static void  Print_Canonical_Overlap
    (uint32 a_id, int a_len, char a_comp,
     uint32 b_id, int b_len, char b_comp, 
     int a_hang, int b_hang, int errors, int min_overlap);
static void  Output_Actual_Overlaps
    (FragmentArrayp fragments, int min_overlap,
     char * fragstore_name);
static void  Calculate_Overlap
    (char a [], char b [], int * a_hang, int * b_hang,
     int * errors, double * score);
static void  Rev_Complement
    (char * S);
static char  Complement
    (char Ch);
static void  Extract_Clear_Sequence
    (ReadStructp rs, char seq []);

/*********************************************************************/


/*********************************************************************/
/* Function:
     main()
   Description:
     top-level function for AS_ORA
   Return Value:
     0 if ok
   Parameters:
     <filename>.???    Binary frag store directory name
     <filename>.ovl    ASCII overlap message filename
                       (NOTE: This will change to a binary file when
                              the I/O routines (and files) are ready)
     min_interval      Minimum overlap interval to pay attention to
*/
int main( int argc, char ** argv )
{
  char              * fragstore_name = NULL;
  char              * input_ovl_filename = NULL;
  char              * output_ovl_filename = NULL;
  uint32            heap_overlaps;
  FragmentArrayp    fragments;
  OverlapStatistics statistics;
  int               output_actual = FALSE;

  // initialize structures
  InitializeStatistics( &statistics );
  
  // parse the command line parameters
  // use getopt(): see "man 3 getopt"
  {
    int ch, errflg = 0;
    optarg = NULL;
    while( !errflg && ((ch = getopt( argc, argv, "as:i:o:l:" )) != EOF) )
    {
      switch( ch )
      {
        case 'a' :
          output_actual = TRUE;
          break;
        case 's':
          fragstore_name = optarg;
          break;
        case 'i':
          input_ovl_filename = optarg;
          break;
        case 'o':
          output_ovl_filename = optarg;
          break;
        case 'l':
          statistics.min_overlap = atoi( optarg );
          break;
        case '?':
          fprintf( stderr, "Unrecognized option -%c\n", optopt );
        default:
          errflg++;
          break;
      }
    }

    // need fragstore_name & min_overlap and one or both of
    // input and output ovl filenames
    if( errflg != 0 ||
        fragstore_name == NULL ||
        statistics.min_overlap == 0 ||
        (input_ovl_filename == NULL && output_ovl_filename == NULL) )
    {
      fprintf( stderr, "Usage: %s\n"
               "       -s fragstorename\n"
               "       -l min-valid-overlap-size\n"
               "       [-i input-overlap-filename |"
               " -o output-overlap-filename]\n",
               argv[0] );
      fprintf (stderr,
               "   -a  to compute and print overlaps using full dyn pgmming\n");
    return 1;
    }
  }

  // read in the fragment store
  fprintf( stderr, "reading fragments...\n" );
  fragments = ReadFragments( fragstore_name );
  if( fragments == NULL )
  {
    fprintf( stderr,
             "Failed to read fragment data from store %s into array.\n",
             fragstore_name );
    return 1;
  }

  // sort the fragments array by min & max sequence interface
  fprintf( stderr, "sorting fragments by start & end indices...\n" );
  SortFragmentsByInterval( fragments );

  // get a number of overlaps to allocate an initial heap
  if( input_ovl_filename != NULL )
  {
    heap_overlaps = GetNumFoundOverlaps( input_ovl_filename );
    if( heap_overlaps == 0 )
    {
      fprintf( stderr, "No overlaps in overlaps file %s.\n",
               input_ovl_filename );
      return 1;
    }
  }
  else
  {
    heap_overlaps = 10 * fragments->num_items;
  }
  
  // generate true overlaps
  fprintf( stderr, "computing 'true' overlaps...\n" );

  if  (output_actual)
      {
       // Do full dynamic-programming alignment on the actual fragments
       // including sampling errors.

       Output_Actual_Overlaps (fragments, statistics . min_overlap, fragstore_name);

       return  0;
      }

  if( ComputeTrueOverlaps( fragments,
                           &statistics,
                           heap_overlaps ) )
  {
    fprintf( stderr, "Failed to compute true overlaps.\n" );
    return 1;
  }

  // Compare computed-true overlaps with found overlaps in overlaps file
  if( input_ovl_filename != NULL )
  {
    fprintf( stderr,
             "reading 'found' overlaps and comparing with 'true' ones...\n" );
    if( CompareOverlaps( fragments,
                         input_ovl_filename,
                         &statistics ) )
    {
      fprintf( stderr, "Failed to compare true and found overlaps.\n" );
      return 1;
    }

    if( ComputeGaps( fragments, &statistics ) )
    {
      fprintf( stderr, "Failed to compute gap statistics.\n" );
      return 1;
    }
    
    fprintf( stderr, "printing results...\n" );
    PrintStatistics( fragments,
                     fragstore_name,
                     input_ovl_filename,
                     &statistics );
  }

  if( output_ovl_filename != NULL )
  {
    fprintf( stderr, "generating overlap file with 'true' overlaps...\n" );
    if( GenerateOverlapFile( output_ovl_filename,
                             fragstore_name,
                             fragments ) )
    {
      fprintf( stderr,
               "Failed to write overlap file %s.\n",
               output_ovl_filename );
      return 1;
    }
  }
  
  // Clean up
  FreeFragmentArray( fragments );
  FreeStatistics( &statistics );

  fprintf( stderr, "Done.\n" );
  return 0;
}



static void  Output_Actual_Overlaps
    (FragmentArrayp fragments, int min_overlap,
     char * fragstore_name)

//  Output to file  actual.olaps  all real overlaps
//  between fragments in  fragments  by doing full dynamic
//  programming on the actual fragment sequence.
//  The sequence is obtained from the fragment store named
//  fragstore_name .   min_overlap  is the minimum length of the
//  overlap region.

  {
   FragStoreHandle  store_handle;
   ReadStructp  read_struct;
   int  left;

   // open the store
   store_handle = openFragStore (fragstore_name, "r");
   if  (store_handle < 0)
       {
        fprintf (stderr,
                 "Failed to open fragment store %s.\n",
                 fragstore_name );
        exit (-1);
       }

   // get a reusable read structure
   read_struct = new_ReadStruct ();
   if  (read_struct == NULL)
       {
        fprintf (stderr, "Failed to allocate read structure.\n");
        closeFragStore( store_handle );
        exit (-1);
       }

   for  (left = 0;  left < fragments -> num_items - 1;  left ++)
     {
      char  left_seq [2 * AS_READ_MAX_LEN];
      int  right;

      getFragStore (store_handle, (fragments -> ptrs [left]) -> internal_id,
                    FRAG_S_SEQUENCE, read_struct);

      Extract_Clear_Sequence (read_struct, left_seq);

      if  ((fragments -> ptrs [left]) -> wc_complement)
          Rev_Complement (left_seq);

      for  (right = left + 1;
            right < fragments -> num_items
                && (fragments -> ptrs [right]) -> begin + min_overlap
                     <= (fragments -> ptrs [left]) -> end;
            right ++)
        {
         char  right_seq [2 * AS_READ_MAX_LEN];
         int  a_hang, b_hang, errors;
         double  score;

         getFragStore (store_handle, (fragments -> ptrs [right]) -> internal_id,
                       FRAG_S_SEQUENCE, read_struct);

         Extract_Clear_Sequence (read_struct, right_seq);

         if  ((fragments -> ptrs [right]) -> wc_complement)
             Rev_Complement (right_seq);

         Calculate_Overlap (left_seq, right_seq, & a_hang, & b_hang,
                            & errors, & score);
         Print_Canonical_Overlap
             ((fragments -> ptrs [left]) -> internal_id,
              strlen (left_seq),
              (fragments -> ptrs [left]) -> wc_complement,
              (fragments -> ptrs [right]) -> internal_id,
              strlen (right_seq),
              (fragments -> ptrs [right]) -> wc_complement,
              a_hang, b_hang, errors, min_overlap);
        }
     }


   // clean up
   delete_ReadStruct (read_struct);
   closeFragStore (store_handle);

   return;
  }



static void  Print_Canonical_Overlap
    (uint32 a_id, int a_len, char a_comp,
     uint32 b_id, int b_len, char b_comp, 
     int a_hang, int b_hang, int errors, int min_overlap)

//  Print the overlap indicated by the parameters in a canonical
//  form for comparison with overlapper overlaps.   a_id  and  b_id
//  are the fragment id's;  a_len  and  b_len  are their respective
//  lengths;  a_comp  and  b_comp  if true indicate the respective fragments
//  came from the reverse strand of the celsim genome.   a_hang  and
//  b_hang  are the amounts the fragments extend past the other fragment
//  (these may be negative).   errors  is the number of errors on the
//  alignment.
//
//  Canonical form is:
//     low id first, then high id
//     low id is automatically in the forward direction
//     high id is followed by  'f'  if it's forward,  'r'  if it's reversed
//     then the "low hang":  the number of bases the low fragment "hangs
//     off the left end" of the high fragment--negative if high is to
//     the left of low
//     then the error rate:  the number of errors in the alignment
//     divided by the average of the a-frag and b-frag lengths in the
//     alignment region.
//  Print nothing is average alignment-region length < min_overlap .

  {
   uint32  lo_id, hi_id;
   int  output_hang;
   int  a_olap_len, b_olap_len;
   char  output_comp;
   double  olap_len;
   double  error_rate;

   assert (a_id != b_id);

   if  (a_id < b_id)
       {
        lo_id = a_id;
        hi_id = b_id;
        if  (a_comp)
            {
             output_hang = - b_hang;
             output_comp = ! b_comp;
            }
          else
            {
             output_hang = a_hang;
             output_comp = b_comp;
            }
       }
     else
       {
        lo_id = b_id;
        hi_id = a_id;
        if  (b_comp)
            {
             output_hang = b_hang;
             output_comp = ! a_comp;
            }
          else
            {
             output_hang = - a_hang;
             output_comp = a_comp;
            }
       }

   if  (a_hang >= 0)
       {
        if  (b_hang >= 0)
            {
             a_olap_len = a_len - a_hang;
             b_olap_len = b_len - b_hang;
            }
          else
            {
             a_olap_len = a_len - a_hang + b_hang;
             b_olap_len = b_len;
            }
       }
     else
       {
        if  (b_hang >= 0)
            {
             a_olap_len = a_len;
             b_olap_len = b_len - b_hang + a_hang;
            }
          else
            {
             a_olap_len = a_len + b_hang;
             b_olap_len = b_len + a_hang;
            }
       }
   olap_len = (a_olap_len + b_olap_len) / 2.0;

   if  (olap_len < min_overlap)
       return;

   assert (olap_len > 0.0);
   error_rate = errors / olap_len;
   printf ("%8u %8u %c %4d %5.2f %3d %3d\n",
           lo_id, hi_id, output_comp ? 'r' : 'f', output_hang,
           100.0 * error_rate, a_olap_len, b_olap_len);

   return;
  }


typedef  struct
  {
   float  score;
   int  a_hang;
   int  errors;
  }  Align_t;

#define  MATCH_SCORE  1.0
#define  MISMATCH_SCORE  -2.0
#define  INSERT_SCORE  -2.0
#define  DELETE_SCORE  -2.0

static void  Calculate_Overlap
    (char a [], char b [], int * a_hang, int * b_hang,
     int * errors, double * score)

//  Find the best overlap between sequences  a  and  b .  Set
//  a_hang  and  b_hang  to the number of bases that  a  extends
//  to the left of  b , and that  b  extends to the right of  a
//  respectively.  Both numbers can be negative.  Set  errors
//  to the number of insert, delete and mismatch errors on the
//  alignment.

  {
   static Align_t  matrix [2 * AS_READ_MAX_LEN];
   double  best_score;
   int  a_len, b_len;
   int  i, j;

   a_len = strlen (a);
   b_len = strlen (b);

   assert (b_len < 2 * AS_READ_MAX_LEN - 1);

if  (Global_Debug)
    printf ("a_len = %d  b_len = %d\n", a_len, b_len);

   // initialize row 0 of the alignment scoring matrix
   for  (j = 0;  j <= b_len;  j ++)
     {
      matrix [j] . score = 0.0;
      matrix [j] . a_hang = - j;
      matrix [j] . errors = 0;
     }

   best_score = -1.0;
   * a_hang = a_len;
   * b_hang = b_len;
   * errors = 0;
   * score = best_score;

   for  (i = 1;  i <= a_len;  i ++)
     {
      Align_t  newa, prev;

      prev . score = 0.0;
      prev . a_hang = i;
      prev . errors = 0;

      for  (j = 1;  j <= b_len;  j ++)
        {
         double  tmp;

         newa = matrix [j - 1];
         if  (a [i - 1] == b [j - 1])
             newa . score += MATCH_SCORE;
           else
             {
              newa . score += MISMATCH_SCORE;
              newa . errors ++;
             }

         tmp = matrix [j] . score + DELETE_SCORE;
         if  (tmp > newa . score)
             {
              newa = matrix [j];
              newa . score = tmp;
              newa . errors ++;
             }

         tmp = prev . score + INSERT_SCORE;
         if  (tmp > newa . score)
             {
              newa = prev;
              newa . score = tmp;
              newa . errors ++;
             }

         matrix [j - 1] = prev;
         prev = newa;
        }

      matrix [b_len] = prev;
      if  (prev . score > best_score)
          {
           best_score = prev . score;
           * a_hang = prev . a_hang;
           * b_hang = i - a_len;
           * errors = prev . errors;
           * score = best_score;
          }
     }

   for  (j = 1;  j < b_len;  j ++)
     if  (matrix [j] . score > best_score)
         {
          best_score = matrix [j] . score;
          * a_hang = matrix [j] . a_hang;
          * b_hang = b_len - j;
          * errors = matrix [j] . errors;
          * score = best_score;
         }

   return;
  }



static void  Rev_Complement
    (char * S)

/* Set  S [0 .. ]  to its DNA reverse complement. */

  {
   int  i, j, len;

   len = strlen (S);

   for  (i = 0, j = len - 1;  i < j;  i ++, j --)
     {
      char  ch;

      ch = Complement (S [i]);
      S [i] = Complement (S [j]);
      S [j] = ch;
     }

   if  (i == j)
       S [i] = Complement (S [i]);

   return;
  }



static char  Complement
    (char Ch)

/*  Return the DNA complement of  Ch . */

  {
   switch  (tolower ((int) Ch))
     {
      case  'a' :
        return  't';
      case  'c' :
        return  'g';
      case  'g' :
        return  'c';
      case  't' :
        return  'a';
      default :
        fprintf (stderr, "ERROR(complement):  Unexpected character `%c\'\n", Ch);
        exit (-1);
     }
  }



static void  Extract_Clear_Sequence
    (ReadStructp rs, char seq [])

//  Extract the clear portion of the sequence in  rs , make it
//  all lower-case, and put it in  seq .

  {
   char  quality [2 * AS_READ_MAX_LEN];
   int  result;
   unsigned int  clear_start, clear_end;
   int  i, len;

   result = getSequence_ReadStruct
                (rs, seq, quality, 2 * AS_READ_MAX_LEN);
   if  (result != 0)
       {
        fprintf (stderr, "Error reading frag store\n");
        exit (EXIT_FAILURE);
       }

   len = strlen (seq);
   for  (i = 0;  i < len;  i ++)
     seq [i] = tolower (seq [i]);

   result = getClearRegion_ReadStruct
                (rs, & clear_start, & clear_end, READSTRUCT_LATEST);
   if  (result != 0)
       {
        fprintf (stderr, "Error reading frag store\n");
        exit (EXIT_FAILURE);
       }

   len = clear_end - clear_start;
   if  (clear_start > 0)
       memmove (seq, seq + clear_start, len);
   seq [len] = '\0';

   return;
  }




/* Function:
     ComputeTrueOverlaps
   Description:
     Walks through the fragments array & determines which ones overlap in
     terms of start & stop indices
   Return Value:
     0 if ok
   Parameters:
     FragmentArrayp fragments: array structure of fragments
     OverlapStatisticsp stats: pointer to statistics structure
     uint32 heap_overlaps: number of overlaps to allocate initially
*/
int ComputeTrueOverlaps( FragmentArrayp fragments,
                         OverlapStatisticsp stats,
                         uint32 heap_overlaps )
{
  uint64 focus;
  uint64 relative;
  int gap_flag;

  // allocate a heap of overlaps to draw from
  fragments->overlap_heap = AllocateOverlapHeap( heap_overlaps );
  if( fragments->overlap_heap == NULL )
  {
    fprintf( stderr,
             "Failed to allocate overlap heap of %d overlaps\n",
             heap_overlaps );
    return 1;
  }

  stats->max_overlap = 0;

  // loop over sorted fragments to find overlaps
  // main loop steps through all fragments but last
  for( focus = 0, relative = 1;
       focus < fragments->num_items - 1;
       focus++, relative = focus + 1 )
  {
    gap_flag = 1;
    // inner loop steps through set of "relative" fragments
    // to compare with "focus" fragment
    while( relative < fragments->num_items &&
           (fragments->ptrs[relative])->begin + stats->min_overlap <=
           (fragments->ptrs[focus])->end )
    {
      // if overlaps array was passed in, populate an overlap structure
      if( AddOverlap( fragments, focus, relative ) )
      {
        fprintf( stderr, "Failed to add overlap of fragments.\n" );
        return 1;
      }

      // increment the number of true overlaps
      stats->num_true++;
      
      // compute maximum overlap interval for histogram allocation later
      stats->max_overlap = GetMaxOverlap( stats->max_overlap,
                                          fragments->ptrs[focus],
                                          fragments->ptrs[relative] );
      
      // increment overlap counter & relative fragment index
      relative++;
      gap_flag = 0;
    }
  }

  return 0;
}

/* Function:
     AddOverlap
   Description:
     Fills an OverlapObject given two overlapping fragments
   Return Value:
     0 if ok
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     uint64 i1: index of first fragment (in fragments->ptrs array)
     uint64 i2: index of second fragment (in fragments->ptrs array)
   NOTE:
     This function assumes that i1's begin <= i2's begin
*/
int AddOverlap( FragmentArrayp fragments, uint64 i1, uint64 i2 )
{
  OverlapObjectp overlap = GetOverlapObject( fragments->overlap_heap );
  if( overlap == NULL )
  {
    fprintf( stderr, "Failed to get overlap object from heap.\n" );
    return 1;
  }

  // set the min & max ids for easy sorting & walking
  overlap->frag_id = fragments->ptrs[i2]->internal_id;

  // determine containment or dovetail & compute a_hang & b_hang
  if( fragments->ptrs[i1]->begin <= fragments->ptrs[i2]->begin &&
      fragments->ptrs[i1]->end >= fragments->ptrs[i2]->end )
  {
    overlap->type = AS_CONTAINMENT;
  }
  else
  {
    overlap->type = AS_DOVETAIL;
  }
  overlap->a_hang = fragments->ptrs[i2]->begin - fragments->ptrs[i1]->begin;
  overlap->b_hang = (int32) ((int64) fragments->ptrs[i2]->end -
                             (int64) fragments->ptrs[i1]->end);

  // determine orientation
  if( fragments->ptrs[i1]->wc_complement ==
      fragments->ptrs[i2]->wc_complement )
  {
    overlap->orientation = AS_NORMAL;
  }
  else if( fragments->ptrs[i1]->wc_complement == 1 )
  {
    overlap->orientation = AS_OUTTIE;
  }
  else
  {
    overlap->orientation = AS_INNIE;
  }

  // Not used for true overlaps
  overlap->min_offset = overlap->a_hang;
  overlap->max_offset = overlap->a_hang;

  InsertionSortOverlap( &(fragments->ptrs[i1]->overlaps), overlap );

  return 0;
}

/* Function:
     ComputeGaps
   Description:
     computes sequence & overlap gaps
   Return Value:
     0 if ok
   Parameters:
     FragmentArrayp fragments: structure holding fragment array & computed
                               overlaps
     OverlapStatisticsp stats: pointer to statistics structure
*/
int ComputeGaps( FragmentArrayp fragments, OverlapStatisticsp stats )
{
  uint64 focus;
  uint64 end_base = fragments->ptrs[0]->end;
  uint64 end_frag = 0;
  OverlapObjectp ovl = NULL;

  // loop over sorted fragments to find overlaps
  // main loop steps through all fragments but last
  for( focus = 0; focus < fragments->num_items - 1; focus++ )
  {
    if( end_base - stats->min_overlap < fragments->ptrs[focus]->begin )
    {
      stats->num_ovl_gaps++;
      stats->ovl_gap_bases += fragments->ptrs[focus]->begin -
        (end_base - stats->min_overlap);
      if( end_base < fragments->ptrs[focus]->begin )
      {
        stats->num_seq_gaps++;
        stats->seq_gap_bases += fragments->ptrs[focus]->begin - end_base;
      }
    }
    end_base = MAX( end_base, fragments->ptrs[focus]->end );

    if( end_frag < focus )
      stats->num_chunks++;
      
    ovl = fragments->ptrs[focus]->overlaps;
    while( ovl != NULL )
    {
      end_frag = MAX( end_frag,
                      fragments->data[ovl->frag_id - fragments->min_id].ptr_index );
      ovl = ovl->next;
    }
  }

  return 0;
}

/* Function:
     CompareOverlaps
   Description:
     compares computed-true overlaps with found overlaps from file
   Return Value:
     0 if ok
   Parameters:
     FragmentArrayp fragments: structure holding fragment array & computed
                               overlaps
     char * input_ovl_filename: name of overlap message file
     OverlapStatisticsp stats: pointer to statistics structure
*/
int CompareOverlaps( FragmentArrayp fragments,
                     char * input_ovl_filename,
                     OverlapStatisticsp stats )
{
  FILE * fp = fopen( input_ovl_filename, "r" );
  OverlapObjectp foop = NULL;
  GenericMesg * gmesg = NULL;
  OverlapMesg * osp = NULL;
  uint32        with_id;
  uint32        overlap_size;

  if( fp == NULL )
  {
    fprintf( stderr,
             "Failed to open overlap store %s for reading.\n",
             input_ovl_filename );
    return 1;
  }

  // check that it makes sense to collect statistics
  if( stats->max_overlap <= stats->min_overlap )
  {
    fprintf( stderr,
             "Max overlap interval (%d) is less than min specified (%d).\n",
             stats->max_overlap, stats->min_overlap );
    return 1;
  }
  // allocate the false_negatives interval histogram
  stats->histo_false_negatives =
    AllocateHistogram( (stats->max_overlap - stats->min_overlap) + 1 );
  stats->histo_fn_critical =
    AllocateHistogram( (stats->max_overlap - stats->min_overlap) + 1 );
  stats->histo_false_positives = AllocateHistogram( 1500 );
  if( stats->histo_false_negatives == NULL ||
      stats->histo_fn_critical == NULL ||
      stats->histo_false_positives == NULL )
  {
    fprintf( stderr, "Failed to allocate one or more histograms.\n" );
    return 1;
  }

  // read the found overlaps in one-at-a-time
  while( ReadProtoMesg_AS( fp, &gmesg ) != EOF )
  {
    if( gmesg && gmesg->t == MESG_OVL )
    {
      // deal with the generic message as an overlap message
      osp = (OverlapMesg *) gmesg->m;
      
      // increment the number of found overlaps
      stats->num_found++;
      
      // search for the overlap among the fragments
      foop = FindTrueOverlap( fragments, osp->aifrag, osp->bifrag, &with_id );
      if( foop == NULL )
      {
        // false positive
        stats->num_false_positives++;
        UpdateRepeatStats( stats, fragments, osp );
        overlap_size = fragments->data[osp->aifrag - fragments->min_id].end -
          fragments->data[osp->aifrag - fragments->min_id].begin - osp->ahg;
        overlap_size = MIN( overlap_size, overlap_size - osp->bhg );
        overlap_size =
          MAX( overlap_size, stats->min_overlap ) - stats->min_overlap;
        stats->histo_false_positives->bins[overlap_size]++;
      }
      else
      {
        if( MinHangDeltaSum( fragments, foop, with_id, osp ) <
            MAX_OVERLAP_MISMATCH )
        {
          foop->found_count++;
          SetOverlapBestMatch( fragments, with_id, foop, osp );
        }
        else
        {
          stats->num_near_misses++;
        }
      }
    }
  }

  // clean up
  fclose( fp );

  // gather statistics in a final loop through fragments
  FinalizeStatistics( fragments, stats );
  
  return 0;
}

/* Function:
     SetOverlapBestMatch
   Description:
     given a 'true' overlap and a corresponding 'found' overlap, set the
     true overlap's best_match structure member to hold the best of the
     matching found overlaps
   Return Value:
     0 if ok
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     uint32 a_id: internal id of the fragment referencing the true overlap
     OverlapObjectp foop: pointer to true overlap object
     OverlapMesg * osp: pointer to matching found overlap message
*/
int SetOverlapBestMatch( FragmentArrayp fragments,
                         uint32 a_id,
                         OverlapObjectp foop,
                         OverlapMesg * osp )
{
  int32 hds = MinHangDeltaSum( fragments, foop, a_id, osp );
  
  if( foop->best_match == NULL )
  {
    foop->best_match = GetOverlapObject( fragments->overlap_heap );
    if( foop->best_match == NULL )
    {
      fprintf( stderr, "Failed to get overlap object from heap.\n" );
      return 1;
    }
    PopulateOverlapObject( foop->best_match, osp, a_id, hds );
  }
  else
  {
    // pick the better match of the two matching found overlaps to store
    if( hds < foop->best_match->hds )
    {
      PopulateOverlapObject( foop->best_match, osp, a_id, hds );
    }
  }
  
  return 0;
}

/* Function:
     MinHangDeltaSum
   Description:
     Computes a minimum composite a_hang + b_hang comparison value given that
     two overlaps (one true & one found) are already similar
     (i.e., have compatible types & orientations)
   Return Value:
     returns the best possible legitiamte match up composite a_hang & b_hang
     comparison value between two overlaps.
     if overlaps are not a legitimate match, returns MAX_HANG_DELTA_SUM
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     OverlapObjectp oop: pointer to overlap object to populate
     uint32 a_id: id of A fragment in true overlap
     OverlapMesg * osp: pointer to overlap message containing data
*/
int32 MinHangDeltaSum( FragmentArrayp fragments,
                       OverlapObjectp oop,
                       uint32 a_id,
                       OverlapMesg * osp )
{
  // NOTE: comparing overlaps is more intricate than it might seem
  // The regressor uses indices that do not take sequencing errors
  // into account. The regressor also knows which fragments have been
  // reported as reverse complements.
  //
  // So, overlaps may be a good match even when:
  //   1. overlap types don't match but a_hangs are very small
  //   2. orientations don't match for containment overlaps
  //   3. fragment A is listed as fragment B for dovetails
  //   4. fragment A is listed as fragment B for containments
  //      but both a_hangs & b_hangs are very small

  // check to see if fragments are nearly the same
  if( abs( (int64) fragments->data[a_id - fragments->min_id].begin -
           (int64) fragments->data[oop->frag_id - fragments->min_id].begin ) <
      MAX_OVERLAP_MISMATCH / 2 &&
      abs( (int64) fragments->data[a_id - fragments->min_id].end -
           (int64) fragments->data[oop->frag_id - fragments->min_id].end ) <
      MAX_OVERLAP_MISMATCH / 2 )
  {
    // fragments are nearly the same
    // compute favorable comparison for overlapper
    return( MIN( abs( oop->a_hang - osp->ahg ) +
                 abs( oop->b_hang - osp->bhg ),
                 MIN( abs( oop->b_hang - osp->ahg ) +
                      abs( oop->a_hang - osp->bhg ),
                      abs( oop->b_hang + osp->ahg ) +
                      abs( oop->a_hang + osp->bhg ) ) ) );
  }
  
  if( oop->type == AS_DOVETAIL && osp->overlap_type == AS_DOVETAIL )
  {
    // dovetail type match
    if( oop->orientation == osp->orientation )
    {
      // orientation match
      if( a_id == osp->aifrag )
      {
        // a == a
        // hangs must match up
        return( abs( oop->a_hang - osp->ahg ) +
                abs( oop->b_hang - osp->bhg ) );
      }
      else
      {
        // a == b
        // hangs are switched
        return( abs( oop->b_hang - osp->ahg ) +
                abs( oop->a_hang - osp->bhg ) );
      }
    }
    else
    {
      // dovetail match with orientation mismatch
      // this is only ok if fragments are nearly the same - handled above
      return MAX_HANG_DELTA_SUM;
    }
  }
  else if( oop->type == AS_CONTAINMENT && osp->overlap_type == AS_CONTAINMENT )
  {
    // containment type match
    if( oop->orientation == osp->orientation )
    {
      // orientation match
      if( a_id == osp->aifrag )
      {
        // a == a
        // check to see if a or b were complemented
        switch( fragments->data[a_id - fragments->min_id].wc_complement +
              fragments->data[oop->frag_id - fragments->min_id].wc_complement )
        {
          case 0:
          case 2:
            // neither were reversed
            return( abs( oop->a_hang - osp->ahg ) +
                    abs( oop->b_hang - osp->bhg ) );
            break;
          case 1:
            // exactly one was reversed - give overlapper benefit of doubt
            return( MIN( abs( oop->a_hang - osp->ahg ) +
                         abs( oop->b_hang - osp->bhg ),
                         abs( oop->a_hang + osp->bhg ) +
                         abs( oop->b_hang + osp->ahg ) ) );
            break;
            /*
          case 2:
            // both were reversed
            return( abs( oop->a_hang + osp->bhg ) +
                    abs( oop->b_hang + osp->ahg ) );
            break;
            */
        }
      }
      else
      {
        // a == b
        // this is only ok if fragments are nearly the same - handled above
        return MAX_HANG_DELTA_SUM;
      }
    }
    else
    {
      // containment match with orientation mismatch
      // acceptable provided neither orientation is AS_NORMAL
      if( oop->orientation != AS_NORMAL && osp->orientation != AS_NORMAL )
      {
        // hangs are switched
        return( abs( oop->a_hang + osp->bhg ) +
                abs( oop->b_hang + osp->ahg ) );
      }
      else
      {
        // bad match
        return MAX_HANG_DELTA_SUM;
      }
    }
  }
  else 
  {
    // overlap type mismatch
    // NOTE: handle the case where my a_hang == 0 & i call it a containment
    // but overlapper calls it a dovetail with b_hang == 0
    // see GenerateOverlapFile() section where OVL messages are generated
    if( oop->a_hang == 0 &&
        oop->type == AS_DOVETAIL )
    {
      // punt for now - give overlapper benefit of the doubt
      return( MIN( abs( oop->a_hang - osp->ahg ) +
                   abs( oop->b_hang + osp->bhg ),
                   abs( oop->a_hang - osp->bhg ) +
                   abs( oop->b_hang - osp->ahg ) ) );
    }
    else
    {
      // this is only ok if fragments are nearly the same - handled above
      return MAX_HANG_DELTA_SUM;
    }
  }

  // if here, something went wrong
  return MAX_HANG_DELTA_SUM;
}

/* Function:
     PopulateOverlapObject
   Description:
     populates an overlap object with overlap message data
   Return Value:
     none
   Parameters:
     OverlapObjectp oop: pointer to overlap object to populate
     OverlapMesg * osp: pointer to overlap message containing data
     uint32 a_id: id of fragment holding oop - used to identify the other id
     int32 hds: hang delta sum of overlap comparison
*/
void PopulateOverlapObject( OverlapObjectp oop,
                            OverlapMesg * osp,
                            uint32 a_id,
                            int32 hds )
{
  if( osp->aifrag == a_id )
    oop->frag_id = osp->bifrag;
  else
    oop->frag_id = osp->aifrag;
  oop->hds = hds;
  oop->a_hang = osp->ahg;
  oop->b_hang = osp->bhg;
  oop->min_offset = osp->min_offset;
  oop->max_offset = osp->max_offset;
}

/* Function:
     GenerateOverlapFile
   Description:
     creates & populates a .ovl file containing 'true' overlaps
   Return Value:
     0 if ok
   Parameters:
     char * output_ovl_filename: name of .ovl file to create
     char * fragstore_name: name of fragment store
     FragmentArrayp fragments: pointer to fragment array structure
*/
int GenerateOverlapFile( char * output_ovl_filename,
                         char * fragstore_name,
                         FragmentArrayp fragments )
{
  GenericMesg outMesg;
  FILE *fp;
  
  // open the .ovl file to write
  fp = fopen( output_ovl_filename, "w" );
  if( fp == NULL )
  {
    fprintf( stderr,
             "Failed to open file %s for writing.\n",
             output_ovl_filename );
    return 1;
  }
  
  // create & write ADT messages
  {
    AuditMesg   auditMesg;
    AuditLine   auditLine1;
    AuditLine   auditLine2;
    char comment[256];
    
    // bogus Celsim ADL message
    auditLine1.complete = time(0);
    auditLine1.name = "Celsim";
    auditLine1.version = "1.29 1999/01/07";
    sprintf( comment, "Genome Length is " F_U64, fragments->sequence_size );
    auditLine1.comment = comment;
    auditLine1.next = &auditLine2;
    
    // regressor ADL message
    auditLine2.complete = time(0);
    auditLine2.name = "overlap_regressor";
    auditLine2.version = "$Revision: 1.12 $ $Date: 2007-04-26 14:07:04 $";
    auditLine2.comment = "(empty)";
    auditLine2.next = NULL;

    auditMesg.list = &auditLine1;
    outMesg.m = &auditMesg;
    outMesg.t = MESG_ADT;

    if( WriteProtoMesg_AS( fp, &outMesg ) < 0 )
    {
      fprintf( stderr, "Failed to write ADT message.\n" );
      fclose( fp );
      return 1;
    }
  } // write ADT messages

  // create & write IDT messages
  {
    char store_name[256];
    DistStore dstore_handle;
    DistRecord dist_record;
    InternalDistMesg idt_mesg;
    int index;

    sprintf( store_name, "%s/db.dst", fragstore_name );
    // open the distance store
    dstore_handle = openDistStore( store_name, "r" );
    // NOTE: is this a valid check?
    if( dstore_handle < 0 )
    {
      fprintf( stderr,
               "Failed to open distance store %s.\n",
               fragstore_name );
      fclose( fp );
      return 1;
    }

    // read, convert, & write all IDT messages
    outMesg.m = &idt_mesg;
    outMesg.t = MESG_IDT;
    idt_mesg.action = AS_ADD;
    for( index = getFirstElemStore( dstore_handle );
         index < getLastElemStore( dstore_handle ) + 1;
         index++ )
    {
      getDistStore( dstore_handle, index, &dist_record );
      idt_mesg.eaccession = dist_record.UID;
      idt_mesg.iaccession = dist_record.IID;
      idt_mesg.mean = dist_record.mean;
      idt_mesg.stddev = dist_record.stddev;
      
      if( WriteProtoMesg_AS( fp, &outMesg ) < 0 )
      {
        fprintf( stderr, "Failed to write IDT message.\n" );
        closeDistStore( dstore_handle );
        fclose( fp );
        return 1;
      }
    } // loop to read/write ADT/IDT messages
    
    closeDistStore( dstore_handle );
  } // write IDT messages

  // create & write OFG messages
  {
    FragStoreHandle fstore_handle;
    ReadStructp read_struct;
    OFGMesg ofg_mesg;
    // char sequence[MAX_SEQUENCE_LENGTH];
    // char quality[MAX_SEQUENCE_LENGTH];
    char source[MAX_SOURCE_LENGTH];
    uint32 index;
    
    // open the fragment store
    fstore_handle = openFragStore( fragstore_name, "r" );
    // NOTE: is this a valid check?
    if( fstore_handle < 0 )
    {
      fprintf( stderr,
               "Failed to open fragment store %s.\n",
               fragstore_name );
      fclose( fp );
      return 1;
    }
    
    // get a reusable read structure
    read_struct = new_ReadStruct();
    if( read_struct == NULL )
    {
      fprintf( stderr, "Failed to allocate read structure.\n" );
      closeFragStore( fstore_handle );
      fclose( fp );
      return 1;
    }
    
    // read, create, & write OFG messages
    outMesg.m = &ofg_mesg;
    outMesg.t = MESG_OFG;
    ofg_mesg.action = AS_ADD;
    // write fragments out in begin/end lexicographic order
    for( index = 0; index < fragments->num_items; index++ )
    {
      if( getFragStore( fstore_handle,
                        fragments->ptrs[index]->internal_id,
                        FRAG_S_ALL,
                        read_struct ) == 0 )
      {
        getAccID_ReadStruct( read_struct, &ofg_mesg.eaccession );
        getReadIndex_ReadStruct( read_struct, &ofg_mesg.iaccession );
        getReadType_ReadStruct( read_struct, &ofg_mesg.type );
        //ofg_mesg.type = AS_READ;
        getSource_ReadStruct( read_struct, source, MAX_SOURCE_LENGTH );
        ofg_mesg.source = source;
        getClearRegion_ReadStruct( read_struct,
                                   (unsigned int *) &ofg_mesg.clear_rng.bgn,
                                   (unsigned int *) &ofg_mesg.clear_rng.end,
				   READSTRUCT_LATEST);
/* OFG messages don't have sequence & quality fields
        getSequence_ReadStruct( read_struct, sequence, quality,
                                MAX_SEQUENCE_LENGTH );
        ofg_mesg.sequence = sequence;
        ofg_mesg.quality = quality;
*/
        if( WriteProtoMesg_AS( fp, &outMesg ) < 0 )
        {
          fprintf( stderr, "Failed to write OFG message.\n" );
          delete_ReadStruct( read_struct );
          closeFragStore( fstore_handle );
          fclose( fp );
          return 1;
        }
      }
      else
      {
        fprintf( stderr, "Failed to read fragment index %d\n", index );
      }
    } // read/write in all SFG/OFG messages
    delete_ReadStruct( read_struct );
    closeFragStore( fstore_handle );
  } // write OFG messages

  // create & write OVL messages
  {
    OverlapMesg ovl_mesg;
    uint32 f_index;
    OverlapObjectp ovl = NULL;

    outMesg.m = &ovl_mesg;
    outMesg.t = MESG_OVL;
    // NOTE: the following is not available from fragstore or regressor
    ovl_mesg.quality = 1.0f;
    ovl_mesg.polymorph_ct = 0;
    ovl_mesg.delta = (signed char *) "\0";
    for( f_index = 0; f_index < fragments->num_items; f_index++ )
    {
      ovl = fragments->data[f_index].overlaps;
      while( ovl != NULL )
      {
        // set meaningful items
        // avoid an assert failure on reading in overlaps in AS_MSG_pmesg.c
        // that doesn't allow dovetails when a_hang == 0
        if( ovl->a_hang == 0 )
        {
          // NOTE: b_hang should be positive, since fragments were sorted
          // on begin & end indices
          ovl_mesg.ahg = abs( ovl->b_hang );
          ovl_mesg.aifrag = ovl->frag_id;
          ovl_mesg.bhg = 0;
          ovl_mesg.bifrag = fragments->data[f_index].internal_id;
          ovl_mesg.overlap_type = AS_CONTAINMENT;
          if( ovl->orientation == AS_INNIE )
            ovl_mesg.orientation = AS_OUTTIE;
          else if( ovl->orientation == AS_OUTTIE )
            ovl_mesg.orientation = AS_INNIE;
          else
            ovl_mesg.orientation = ovl->orientation;
        }
        else
        {
          ovl_mesg.ahg = ovl->a_hang;
          ovl_mesg.aifrag = fragments->data[f_index].internal_id;
          ovl_mesg.bhg = ovl->b_hang;
          ovl_mesg.bifrag = ovl->frag_id;
          ovl_mesg.overlap_type = ovl->type;
          ovl_mesg.orientation = ovl->orientation;
        }
          
        // set meaningless items
        ovl_mesg.min_offset = ovl_mesg.ahg;
        ovl_mesg.max_offset = ovl_mesg.ahg;
        
        if( WriteProtoMesg_AS( fp, &outMesg ) < 0 )
        {
          fprintf( stderr, "Failed to write OVL message.\n" );
          fclose( fp );
          return 1;
        }
        ovl = ovl->next;
      } // loop over overlaps of fragment
    } // loop over fragments
  } // write overlaps to .ovl file
  
  // clean up
  fclose( fp );
  
  return 0;
}

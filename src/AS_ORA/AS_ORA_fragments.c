
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/AS_ORA_fragments.c,v $
$Revision: 1.5 $
**********************************************************************/
#define DEBUG_ORA
/*********************************************************************/
// headers
// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// project headers
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_ORA_repeats.h"
#include "AS_ORA_fragments.h"
#include "AS_ORA_inlines.h"
/*********************************************************************/


/*********************************************************************/
// constants
#define INDEX_STRING "\n["
#define ANNOT_STRING "] ["
/*********************************************************************/


/*********************************************************************/
// local function declarations
static int PopulateFragmentObject( FragmentArrayp fragments,
                                   ReadStructp read_struct );
static int ParseStringForAnnotations( FragmentArrayp fragments,
                                      uint32 f_index,
                                      char * string );
static int ParseStringForIndices( FragmentArrayp fragments,
                                  uint32 f_index,
                                  char * src_string );
static int ParseSourceString( FragmentArrayp fragments,
                              uint32 f_index,
                              char * src_string );
static int GetCommaDelimitedIntegers( char * string,
                                      uint64 * first,
                                      uint64 * second );
/*********************************************************************/


/*********************************************************************/
/* Function:
     AllocateFragmentArray
   Description:
     Allocates an array of fragments
   Return Value:
     pointer to valid FragmentArray structure if ok
   Parameters:
     uint32 num_fragments: number of fragment objects to allocate
*/
FragmentArrayp AllocateFragmentArray( uint32 num_fragments )
{
  uint32 i;
  FragmentArrayp fragments;

  // allocate array object
  fragments = (FragmentArrayp) calloc( 1, sizeof( FragmentArray ) );
  if( fragments == NULL )
  {
    fprintf( stderr, "Failed to allocate fragment array object.\n" );
    return NULL;
  }

  // allocate data array
  fragments->data =
    (FragmentObjectp) calloc( num_fragments, sizeof( FragmentObject ) );
  if( fragments->data == NULL )
  {
    fprintf( stderr,
             "Failed to allocate array of %d fragments.\n",
             num_fragments );
    FreeFragmentArray( fragments );
    return NULL;
  }

  // allocate pointer array
  fragments->ptrs =
    (FragmentObjectp *) malloc( num_fragments * sizeof( FragmentObjectp ) );
  if( fragments->ptrs == NULL )
  {
    fprintf( stderr,
             "Failed to allocate array of %d fragment pointers.\n",
             num_fragments );
    FreeFragmentArray( fragments );
    return NULL;
  }

  // get the repeat heap initialized & allocated
  fragments->repeat_heap = AllocateRepeatHeap( num_fragments / 2 );
  if( fragments->repeat_heap == NULL )
  {
    fprintf( stderr, "Failed to allocate repeat heap.\n" );
    FreeFragmentArray( fragments );
    return NULL;
  }

  // make pointers point to data
  for( i = 0; i < num_fragments; i++ )
  {
    fragments->ptrs[i] = &(fragments->data[i]);
  }

  fragments->num_items = num_fragments;
  return fragments;
}

/* Function:
     FreeFragmentArray
   Description:
     Frees a FragmentArray array structure
   Return Value:
     none
   Parameters:
     FragmentArrayp fragments: pointer to FragmentArray structure
*/
void FreeFragmentArray( FragmentArrayp fragments )
{
  if( fragments != NULL )
  {
    if( fragments->data != (FragmentObjectp) NULL )
      free( fragments->data );
    if( fragments->ptrs != (FragmentObjectp *) NULL )
      free( fragments->ptrs );
    FreeRepeatHeap( fragments->repeat_heap );
    FreeOverlapHeap( fragments->overlap_heap );
    free( fragments );
  }
    
}

/* Function:
     ReadFragments
   Description:
     allocates array and reads elements of a fragment store from a file
   Return Value:
     pointer to valid allocated & populated FragmentArray structure if ok
   Parameters:
     char * fragstore_name:  name of fragment store
*/
FragmentArrayp ReadFragments( char * fragstore_name )
{
  FragmentArrayp fragments;
  FragStoreHandle store_handle;
  FragStreamHandle stream_handle;
  ReadStructp read_struct;
  uint32 first_fragment;
  uint32 num_fragments;

  // open the store
  store_handle = openFragStore( fragstore_name, "r" );
  // NOTE: is this a valid check?
  if( store_handle < 0 )
  {
    fprintf( stderr,
             "Failed to open fragment store %s.\n",
             fragstore_name );
    return NULL;
  }

  // determine the number of fragments
  first_fragment = (uint32) getFirstElemFragStore( store_handle );
  num_fragments = ((uint32) getLastElemFragStore( store_handle ) -
                   first_fragment) + 1;
  if( num_fragments == 0 )
  {
    fprintf( stderr,
             "No fragments in fragment store file %s.\n",
             fragstore_name );
    closeFragStore( store_handle );
    return NULL;
  }

  // allocate the memory to hold fragment data
  if( (fragments = AllocateFragmentArray( num_fragments )) == NULL )
  {
    fprintf( stderr, "Failed to allocate fragment array.\n" );
    closeFragStore( store_handle );
    return NULL;
  }
  fragments->min_id = first_fragment;

  // get a stream handle to iterate through fragments
  stream_handle = openFragStream( store_handle, NULL, 0 );
  // NOTE: is this a valid check?
  if( stream_handle < 0 )
  {
    fprintf( stderr, "Failed to open fragment stream on store %s.\n",
             fragstore_name );
    FreeFragmentArray( fragments );
    closeFragStore( store_handle );
    return NULL;
  }

  // get a reusable read structure
  read_struct = new_ReadStruct();
  if( read_struct == NULL )
  {
    fprintf( stderr, "Failed to allocate read structure.\n" );
    FreeFragmentArray( fragments );
    closeFragStream( stream_handle );
    closeFragStore( store_handle );
    return NULL;
  }

  // read in the fragments & copy pertinent data to the fragments array
  while( nextFragStream( stream_handle, read_struct, FRAG_S_SOURCE ) )
  {
    // populate a fragment item
    if( PopulateFragmentObject( fragments, read_struct ) )
    {
      fprintf( stderr, "Failed to set fragment object.\n" );
      delete_ReadStruct( read_struct );
      FreeFragmentArray( fragments );
      closeFragStream( stream_handle );
      closeFragStore( store_handle );
      return NULL;
    }
  }

  // clean up
  delete_ReadStruct( read_struct );
  closeFragStream( stream_handle );
  closeFragStore( store_handle );

  return fragments;
}

/* Function:
     FindTrueOverlap
   Description:
     returns a pointer to a 'true' overlap given a 'found' overlap
   Return Value:
     pointer to true overlap object, if one is found
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     uint32 a: internal ID of first fragment
     uint32 b: internal ID of second fragment
     uint32 * with_id: internal ID of fragment with which overlap was found
*/
OverlapObjectp FindTrueOverlap( FragmentArrayp fragments,
                                uint32 a,
                                uint32 b,
                                uint32 * with_id )
{
  OverlapObjectp oop = fragments->data[a - fragments->min_id].overlaps;
  *with_id = a;
  
  // look at the a fragment for b
  while( oop && oop->frag_id < b )
  {
    oop = oop->next;
  }

  // check the outcome
  // perhaps look at the b fragment for a
  if( oop == (OverlapObjectp) NULL ||
      oop->frag_id != b )
  {
    *with_id = b;
    oop = fragments->data[b - fragments->min_id].overlaps;

    // iterate through ordered list of overlaps of b
    while( oop && oop->frag_id < a )
    {
      oop = oop->next;
    }
    if( oop == (OverlapObjectp) NULL ||
        oop->frag_id != a )
    {
      return( (OverlapObjectp) NULL );
    }
  }

  return( oop );
}

/* Function:
     SortFragmentsByInterval
   Description:
     sorts elements in a fragment array by min & max beginning index
     uses double bucket sort
   Return Value:
     none
   Parameters:
     FragmentArrayp fragments: pointer to structure of fragment array
*/
void SortFragmentsByInterval( FragmentArrayp fragments )
{
  FragmentObjectp * arrayA;
  FragmentObjectp * arrayB;
  FragmentObjectp   tempP;
  uint32 bucket;
  int64 i;
  uint32 sss = (uint32) (sqrt( (double) fragments->sequence_size ) + 1.5);
  
  if( fragments->num_items <= 1 )
    return;

  // allocate arrays for double bucket sort
  arrayA = (FragmentObjectp *) calloc( sss, sizeof( FragmentObjectp ) );
  if( arrayA == NULL )
  {
    fprintf( stderr, "Failed to allocate bucket sort array A.\n" );
    return;
  }
  arrayB = (FragmentObjectp *) calloc( sss, sizeof( FragmentObjectp ) );
  if( arrayB == NULL )
  {
    fprintf( stderr, "Failed to allocate bucket sort array B.\n" );
    free( arrayA );
    return;
  }

  /*******************************************/
  // first sort, on end
  // first half of double bucket sort
  // populate arrayA with end % sss
  // fill forwards
  for( i = 0; i < fragments->num_items; i++ )
  {
    bucket = fragments->data[i].end % sss;
    fragments->data[i].next = arrayA[bucket];
    arrayA[bucket] = &(fragments->data[i]);
  }

  // second half of double bucket sort
  // move arrayA items to arrayB with end / sss
  // remove from back to front
  for( i = sss - 1; i >= 0; i-- )
  {
    while( arrayA[i] != NULL )
    {
      tempP = arrayA[i];
      arrayA[i] = arrayA[i]->next;
      bucket = tempP->end / sss;
      tempP->next = arrayB[bucket];
      arrayB[bucket] = tempP;
    }
  }

  /*******************************************/
  // Second sort, on begin
  // first half of double bucket sort
  // populate arrayA with begin % sss
  // fill forwards
  for( i = 0; i < sss; i++ )
  {
    while( arrayB[i] != NULL )
    {
      tempP = arrayB[i];
      arrayB[i] = arrayB[i]->next;
      bucket = tempP->begin % sss;
      tempP->next = arrayA[bucket];
      arrayA[bucket] = tempP;
    }
  }
  
  /*
  // first half of double bucket sort
  // populate arrayA with begin % sss
  // fill forwards
  for( i = 0; i < fragments->num_items; i++ )
  {
    bucket = fragments->data[i].begin % sss;
    fragments->data[i].next = arrayA[bucket];
    arrayA[bucket] = &(fragments->data[i]);
  }
  */
  
  // second half of double bucket sort
  // move arrayA items to arrayB with begin / sss
  // remove from back to front
  for( i = sss - 1; i >= 0; i-- )
  {
    while( arrayA[i] != NULL )
    {
      tempP = arrayA[i];
      arrayA[i] = arrayA[i]->next;
      bucket = tempP->begin / sss;
      tempP->next = arrayB[bucket];
      arrayB[bucket] = tempP;
    }
  }

  // move arrayB items into fragments->ptrs
  // remove from front to back
  for( bucket = i = 0; bucket < sss; bucket++ )
  {
    while( arrayB[bucket] != NULL )
    {
      fragments->ptrs[i] = arrayB[bucket];
      fragments->ptrs[i]->ptr_index = i;
      arrayB[bucket] = arrayB[bucket]->next;
      i++;
    }
  }

#ifdef DEBUG_ORA
  // check the sort
  for( i = 0; i < fragments->num_items - 1; i++ )
  {
    if( fragments->ptrs[i]->begin > fragments->ptrs[i+1]->begin )
      fprintf( stderr, "Out of order fragments after sort!.\n" );
    if( fragments->ptrs[i]->begin == fragments->ptrs[i+1]->begin &&
        fragments->ptrs[i]->end > fragments->ptrs[i+1]->end )
    {
      fprintf( stderr, "%8d: %12" F_U64P ", %12" F_U64P "\n",
               fragments->ptrs[i]->internal_id,
               fragments->ptrs[i]->begin,
               fragments->ptrs[i]->end );
      fprintf( stderr, "%8d: %12" F_U64P ", %12" F_U64P "\n\n",
               fragments->ptrs[i+1]->internal_id,
               fragments->ptrs[i+1]->begin,
               fragments->ptrs[i+1]->end );
    }
  }
#endif // DEBUG_ORA
  
  // clean up
  free( arrayA );
  free( arrayB );
}

/* Function:
     GetOverlapInterval
   Description:
     returns the size of overlap between two fragments
   Return Value:
     overlap size
   Parameters:
     FragmentObjectp a: pointer to first fragment
     FragmentObjectp b: pointer to second fragment
*/
uint32 GetOverlapInterval( FragmentObjectp a, FragmentObjectp b )
{
  // NOTE: this assumes a->begin < b->begin
  return( (uint32) min( (int64) a->end - (int64) b->begin,
                        (int64) b->end - (int64) b->begin ) );
}

/* Function:
     GetMaxOverlap
   Description:
     returns the max of two fragments' overlap size and a current max size
   Return Value:
     new max overlap size
   Parameters:
     uint32 curr_max: the current maximum overlap size
     FragmentObjectp a: pointer to first fragment
     FragmentObjectp b: pointer to second fragment
*/
uint32 GetMaxOverlap( uint32 curr_max,
                      FragmentObjectp a,
                      FragmentObjectp b )
{
  return( MAX( curr_max, GetOverlapInterval( a, b ) ) );
}

/* Function:
     UpdateFragmentSetRepeats
   Description:
     Updates the fragment set's repeats based on a single fragment's repeats
   Return Value:
     0 if ok
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     uint32 f_index: index into data array of the fragment
*/
int UpdateFragmentSetRepeats( FragmentArrayp fragments, uint32 f_index )
{
  RepeatObjectp f_repeat = fragments->data[f_index].repeats;
  RepeatObjectp s_repeat = fragments->repeats;
  RepeatObjectp new_repeat = NULL;
  RepeatObjectp o_repeat = NULL;

  // loop over fragment's annotations
  while( f_repeat != NULL )
  {
    // find it in the fragment set's repeats
    while( s_repeat != NULL && s_repeat->repeat_id < f_repeat->repeat_id )
      s_repeat = s_repeat->next;

    // if not found
    if( s_repeat == NULL || s_repeat->repeat_id != f_repeat->repeat_id )
    {
      // establish a new repeat
      new_repeat = GetRepeatObject( fragments->repeat_heap );
      if( new_repeat == NULL )
      {
        fprintf( stderr, "Failed to get repeat object from heap.\n" );
        return 1;
      }

      // set its values
      new_repeat->repeat_id = f_repeat->repeat_id;
      new_repeat->rpt_end = f_repeat->rpt_end - f_repeat->rpt_begin;
      new_repeat->count = 1;

      // insert it into the list
      InsertionSortRepeat( &(fragments->repeats), new_repeat );
    }
    else
    {
      // update the size of the repeat
      s_repeat->rpt_end = MAX( s_repeat->rpt_end,
                               f_repeat->rpt_end - f_repeat->rpt_begin );
      // count repeats only once
      if( o_repeat == NULL || o_repeat->repeat_id != f_repeat->repeat_id )
        s_repeat->count++;
    }
    o_repeat = f_repeat;
    f_repeat = f_repeat->next;
  }

  return 0;
}

/* Function:
     IsCriticalOverlap
   Description:
     determines whether or not a false negative overlap is critical in a
     first-order sense - whether or not there is a single other fragment
     that is true & found that overlaps both frags in the false neg overlap
   Return Value:
     0 if not critical, 1 if critical
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     uint32 adi: fragment index in data array
     uint32 bdi: fragment index in data array
*/
int IsCriticalOverlap( FragmentArrayp fragments, uint32 adi, uint32 bdi )
{
  uint32 min_pi = min( fragments->data[adi].ptr_index,
                       fragments->data[bdi].ptr_index );
  uint32 max_pi = MAX( fragments->data[adi].ptr_index,
                       fragments->data[bdi].ptr_index );
  uint32 i;
  uint32 with_id;
  OverlapObjectp foop1;
  
  // first search from min + 1 to max - 1
  for( i = min_pi + 1; i < max_pi; i++ )
  {
    foop1 = FindTrueOverlap( fragments,
                             fragments->ptrs[min_pi]->internal_id,
                             fragments->ptrs[i]->internal_id,
                             &with_id );
    if( foop1 != NULL && foop1->found_count > 0 )
    {
      foop1 =  FindTrueOverlap( fragments,
                             fragments->ptrs[i]->internal_id,
                             fragments->ptrs[max_pi]->internal_id,
                             &with_id );
      if( foop1 != NULL && foop1->found_count > 0 )
        return 0;
    }
  }

  // search forward from max + 1
  i = max_pi + 1;
  while( i < fragments->num_items &&
         (foop1 = FindTrueOverlap( fragments,
                                   fragments->ptrs[min_pi]->internal_id,
                                   fragments->ptrs[i]->internal_id,
                                   &with_id )) != NULL )
  {
    if( foop1->found_count > 0 )
    {
      foop1 =  FindTrueOverlap( fragments,
                                fragments->ptrs[max_pi]->internal_id,
                                fragments->ptrs[i]->internal_id,
                                &with_id );
      if( foop1 != NULL && foop1->found_count > 0 )
        return 0;
    }
    i++;
  }

  // search backward from min - 1
  i = min_pi - 1;
  while( i != 0 &&
         (foop1 = FindTrueOverlap( fragments,
                                   fragments->ptrs[i]->internal_id,
                                   fragments->ptrs[min_pi]->internal_id,
                                   &with_id )) != NULL )
  {
    if( foop1->found_count > 0 )
    {
      foop1 =  FindTrueOverlap( fragments,
                                fragments->ptrs[i]->internal_id,
                                fragments->ptrs[max_pi]->internal_id,
                                &with_id );
      if( foop1 != NULL && foop1->found_count > 0 )
        return 0;
    }
    i--;
  }

  
  return 1;
}

/* Function:
     PopulateFragmentObject
   Description:
     Fills elements of a FragmentObject based on ReadStruct elements
     and the src string
   Return Value:
     0 if ok
   Parameters:
     FragmentObjectp fragment: pointer to FragmentObject to fill
     ReadStruct * read_struct: fragment/read structure to read from
*/
static int PopulateFragmentObject( FragmentArrayp fragments,
                                   ReadStructp read_struct )
{
  char src_string[SOURCE_STRING_LENGTH];
  uint64 temp_swap;
  uint32 f_index;

  // get the fragment id
  getReadIndex_ReadStruct( read_struct, &f_index );
  fragments->data[f_index - fragments->min_id].internal_id = f_index;
  f_index -= fragments->min_id;

  // get src field
  getSource_ReadStruct( read_struct, src_string, SOURCE_STRING_LENGTH );

  // parse src string for repeats & sequence indices
  if( ParseSourceString( fragments, f_index, src_string ) )
  {
    fprintf( stderr,
             "Failed to parse src string of fragment %d.\n",
             fragments->data[f_index].internal_id );
    return 1;
  }

  // set min, max, & complement flag correctly
  if( fragments->data[f_index].begin > fragments->data[f_index].end )
  {
    temp_swap = fragments->data[f_index].end;
    fragments->data[f_index].end = fragments->data[f_index].begin;
    fragments->data[f_index].begin = temp_swap;
    fragments->data[f_index].wc_complement = (char) 1;
  }
  else
  {
    fragments->data[f_index].wc_complement = (char) 0;
  }

  // update the sequence size (for sorting purposes later)
  fragments->sequence_size = MAX( fragments->sequence_size,
                                  fragments->data[f_index].end );

  return 0;
}

/* Function:
     ParseStringForIndices
   Description:
     parses string to extract begin & end indices
   Return Value:
     0 if ok
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     uint32 f_index: index into fragments->data of fragment
     char * string: string to parse
*/
static int ParseStringForIndices( FragmentArrayp fragments,
                                  uint32 f_index,
                                  char * string )
{
  char * temp_char = strstr( string, INDEX_STRING );
  if( !temp_char )
  {
    fprintf( stderr,
             "Failed to access start index for fragment ID %d.\n",
             fragments->data[f_index].internal_id );
    return 1;
  }
  temp_char += strlen( INDEX_STRING );
  
  if( GetCommaDelimitedIntegers( temp_char,
                                 &(fragments->data[f_index].begin),
                                 &(fragments->data[f_index].end) ) )
  {
    fprintf( stderr, "Failed to parse bracketed numbers in string.\n" );
    return 1;
  }
  return 0;
}

/* Function:
     ParseStringForAnnotations
   Description:
     parses string to extract annotation information
   Return Value:
     0 if ok
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     uint32 f_index: index into fragments->data of fragment
     char * string: string to parse
   NOTE: this assumes repeats identifiers are in [A,Z]
*/
static int ParseStringForAnnotations( FragmentArrayp fragments,
                                      uint32 f_index,
                                      char * string )
{
  char * temp_char1 = string;
  char * temp_char2;
  uint64 temp_uint64;
  RepeatObjectp new_repeat = NULL;

  fragments->data[f_index].repeat_id = 0;
  
  // search for all annotations
  while( (temp_char1 = strstr( temp_char1, ANNOT_STRING )) )
  {
    // new annotation
    new_repeat = GetRepeatObject( fragments->repeat_heap );
    if( new_repeat == NULL )
    {
      fprintf( stderr, "Failed to get repeat object from heap.\n" );
      return 1;
    }
    
    // back-up to \n character
    temp_char2 = temp_char1;
    while( temp_char2 > string && temp_char2[0] != '\n' )
    {
      temp_char2--;
    }

    // check that temp_char2 is positioned in a parsable location
    if( temp_char2[0] != '\n' || temp_char2[1] < 'A' || temp_char2[1] > 'Z' )
    {
      fprintf( stderr, "Failed to find repeat identifier.\n" );
      fprintf( stderr,
               "Failed to parse annotations for fragment %d. Continuing.\n",
               fragments->data[f_index].internal_id );
      return 0;
    }
    
    // character after \n is the annotation identifier
    new_repeat->repeat_id = 1 << (int) (temp_char2[1] - 'A');

    // update the fragment's repeat id composite bit vector
    fragments->data[f_index].repeat_id |= new_repeat->repeat_id;

    // get the size of the repeat
    temp_char2 = index( temp_char2, '[' );
    if( temp_char2 == NULL )
    {
      fprintf( stderr, "Failed to find repeat indices.\n" );
      fprintf( stderr,
               "Failed to parse annotations for fragment %d. Continuing.\n",
               fragments->data[f_index].internal_id );
      return 0;
    }
    temp_char2++;
    if( GetCommaDelimitedIntegers( temp_char2,
                                   &(new_repeat->rpt_begin),
                                   &(new_repeat->rpt_end) ) )
    {
      fprintf( stderr,
               "Failed to parse repeat indices in annotation string.\n" );
      fprintf( stderr,
               "Failed to parse annotations for fragment %d. Continuing.\n",
               fragments->data[f_index].internal_id );
      return 0;
    }
    
    // set the repeat indices & wc_complement flag
    if( new_repeat->rpt_begin > new_repeat->rpt_end  )
    {
      temp_uint64 = new_repeat->rpt_end;
      new_repeat->rpt_end = new_repeat->rpt_begin;
      new_repeat->rpt_begin = temp_uint64;
      new_repeat->wc_complement = 1;
    }
    
    // get the indices of the repeat in the fragment
    temp_char2 = index( temp_char2, '[' );
    if( temp_char2 == NULL )
    {
      fprintf( stderr, "Failed to find repeat location in fragment.\n" );
      fprintf( stderr,
               "Failed to parse annotations for fragment %d. Continuing.\n",
               fragments->data[f_index].internal_id );
      return 0;
    }
    temp_char2++;
    if( GetCommaDelimitedIntegers( temp_char2,
                                   &(new_repeat->frag_begin),
                                   &(new_repeat->frag_end) ) )
    {
      fprintf( stderr, "Failed to parse repeat location"
               " in annotation string.\n" );
      fprintf( stderr,
               "Failed to parse annotations for fragment %d. Continuing.\n",
               fragments->data[f_index].internal_id );
      return 0;
    }
    
    // insertion sort the repeat into the list of fragment repeats
    InsertionSortRepeat( &(fragments->data[f_index].repeats), new_repeat );

    // move temp_char1 forward
    temp_char1 = temp_char2;
  }

  return 0;
}

/* Function:
     ParseSourceString
   Description:
     parses fragment src string to extract begin & end indices & annotations
   Return Value:
     0 if ok
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     uint32 f_index: index into fragments->data of fragment
     char * src_string: string to parse
*/
static int ParseSourceString( FragmentArrayp fragments,
                              uint32 f_index,
                              char * src_string )
{
  // parse the src string for begin & end indices
  if( ParseStringForIndices( fragments, f_index, src_string ) )
  {
    fprintf( stderr, "Failed to parse src string for indices.\n" );
    return 1;
  }

  // parse the src string for annotations
  if( ParseStringForAnnotations( fragments, f_index, src_string ) )
  {
    fprintf( stderr, "Failed to parse src string for annotations.\n" );
    return 1;
  }

  // update fragment set's annotations
  if( UpdateFragmentSetRepeats( fragments, f_index ) )
  {
    fprintf( stderr, "Failed to update fragment set annotations.\n" );
    return 1;
  }

  return 0;
}

/* Function:
     GetCommaDelimitedIntegers
   Description:
     parses a string to find front-end comma-delimited integers
   Return Value:
     0 if ok
   Parameters:
     char * string: the string to parse
     uint64 * first: pointer to the first number to return
     uint64 * second: pointer to the second number to return
*/
static int GetCommaDelimitedIntegers( char * string,
                                      uint64 * first,
                                      uint64 * second )
{
  char * temp_char = string;
  
  *first = (uint64) atoi( temp_char );

  temp_char = index( temp_char, ',' );
  if( !temp_char )
  {
    fprintf( stderr,
             "Failed to find comma in delimited number pair [%s\n",
             string );
    return 1;
  }
  temp_char++;
  *second = (uint64) atoi( temp_char );

  return 0;
}

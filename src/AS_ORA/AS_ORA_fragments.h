
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
#ifndef AS_ORA_FRAGMENTS_H
#define AS_ORA_FRAGMENTS_H
/**********************************************************************
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/AS_ORA_fragments.h,v $
$Revision: 1.1.1.1 $
**********************************************************************/


/*********************************************************************/
// headers
// project headers
#include "AS_global.h"
#include "AS_ORA_repeats.h"
#include "AS_ORA_overlaps.h"
/*********************************************************************/


/*********************************************************************/
// defines
#define URC_EXTENSION        "urc"

#ifdef MAX_SOURCE_LENGTH
#define SOURCE_STRING_LENGTH  MAX_SOURCE_LENGTH
#else
#define SOURCE_STRING_LENGTH  512
#endif
/*********************************************************************/


/*********************************************************************/
// types & structures
// FragmentObject
typedef struct FOStruct
{
  uint32            internal_id;     // internal id of fragment
  uint32            ptr_index;       // index of frag in ptrs array
  uint64            begin;           // lower index of fragment in sequence
  uint64            end;             // higher index of fragment in sequence
  char              wc_complement;   // Watson-Crick complement sentinel
  RepeatID          repeat_id;       // bit-vector composite of repeat IDs
  RepeatObjectp     repeats;         // linked list of repeats in this fragment
  OverlapObjectp    overlaps;        // linked list of overlaps for this frag
  struct FOStruct * next;            // linked list used in sorting
} FragmentObject;
typedef FragmentObject * FragmentObjectp;

// structure for array of fragments
typedef struct
{
  FragmentObjectp   data;           // array to hold data
  FragmentObjectp * ptrs;           // array of pointers to point to data
  uint32            num_items;      // number of items in arrays
  uint32            min_id;         // minimum internal id in fragstore
  uint64            sequence_size;  // size of sequence
  RepeatObjectp     repeats;        // list of repeats in fragment set
  RepeatHeapp       repeat_heap;    // heap for getting repeat objects
  OverlapHeapp      overlap_heap;   // heap for getting overlap objects
} FragmentArray;
typedef FragmentArray * FragmentArrayp;
/*********************************************************************/


/*********************************************************************/
// function declarations
/* Function:
     AllocateFragmentArray
   Description:
     Allocates an array of fragments
   Return Value:
     pointer to valid FragmentArray structure if ok
   Parameters:
     uint32 num_fragments: number of fragment objects to allocate
*/
FragmentArrayp AllocateFragmentArray( uint32 num_fragments );

/* Function:
     FreeFragmentArray
   Description:
     Frees a FragmentArray array structure
   Return Value:
     none
   Parameters:
     FragmentArrayp fragments: pointer to FragmentArray structure
*/
void FreeFragmentArray( FragmentArrayp fragments );

/* Function:
     ReadFragments
   Description:
     allocates array and reads elements of a fragment store from a file
   Return Value:
     pointer to valid allocated & populated FragmentArray structure if ok
   Parameters:
     char * fragstore_name:  name of fragment store
*/
FragmentArrayp ReadFragments( char * frag_filename );

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
                                uint32 *with_id );

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
void SortFragmentsByInterval( FragmentArrayp fragments );

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
uint32 GetOverlapInterval( FragmentObjectp a, FragmentObjectp b );

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
                      FragmentObjectp b );

/* Function:
     UpdateFragmentSetRepeats
   Description:
     Updates the fragment set's repeats based on a single fragment's repeats
   Return Value:
     0 if ok
   Parameters:
     FragmentArrayp fragments: pointer to fragment array structure
     uint32 index: index into data array of the fragment
*/
int UpdateFragmentSetRepeats( FragmentArrayp fragments, uint32 index );

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
int IsCriticalOverlap( FragmentArrayp fragments, uint32 adi, uint32 bdi );
/*********************************************************************/
  
#endif // #ifndef AS_ORA_FRAGMENTS_H


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
#ifndef AS_ORA_OVERLAPS_H
#define AS_ORA_OVERLAPS_H
/**********************************************************************
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/AS_ORA_overlaps.h,v $
$Revision: 1.3 $
**********************************************************************/


/*********************************************************************/
// headers
// project headers
#include "AS_global.h"
#include "AS_UTL_heap.h"
#include "AS_MSG_pmesg.h"
/*********************************************************************/


/*********************************************************************/
// constants
#define OVERLAP_EXTENSION     "ovl"
#define OVERLAP_3_CODE        "OVL"
#define AS_ORA_TEMPFILE       "/tmp/AS_ORA_tempfile"
#define MAX_OVERLAP_MISMATCH  10
#define MAX_HANG_DELTA_SUM    0xffff
/*********************************************************************/


/*********************************************************************/
// types & structures
// OverlapObject
typedef struct OOStruct
{
  uint32             frag_id;      // ID of other fragment
  OrientType         orientation;  // orientation of overlap
  OverlapType        type;         // overlap type
  int32              hds;          // hang delta sum
  int32              a_hang;       // a_hang of overlap
  int32              b_hang;       // b_hang of overlap
  int32              min_offset;   // min_offset of overlap
  int32              max_offset;   // max_offset of overlap
  int32              found_count;  // counter for matching found overlaps
  struct OOStruct *  best_match;   // best match among "found" overlaps
  struct OOStruct *  prev;         // pointer to prev overlap in linked list
  struct OOStruct *  next;         // pointer to next overlap in linked list
} OverlapObject;
typedef OverlapObject * OverlapObjectp;

// Generics
typedef GenericArray  OverlapArray;
typedef GenericArrayp OverlapArrayp;
typedef GenericHeap   OverlapHeap;
typedef GenericHeapp  OverlapHeapp;
/*********************************************************************/


/*********************************************************************/
// function declarations
/* Function:
     AllocateOverlapHeap
   Description:
     wrapper around AllocateGenericHeap for overlap objects
   Return Value:
     valid OverlapHeapp pointer if ok
   Parameters:
     uint32 num_overlaps: number of overlaps to allocate initially
*/
OverlapHeapp AllocateOverlapHeap( uint32 num_overlaps );

/* Function:
     GrowOverlapHeap
   Description:
     wrapper around GrowGenericHeap to add memory to overlap heap
   Return Value:
     0 if ok
   Parameters:
     OverlapHeapp overlap_heap: pointer to existing overlap heap
*/
int GrowOverlapHeap( OverlapHeapp overlap_heap );

/* Function:
     FreeOverlapHeap
   Description:
     wrapper around FreeGenericHeap - frees an overlap heap
   Return Value:
     none
   Parameters:
     OverlapHeapp overlap_heap: pointer to existing overlap heap
*/
void FreeOverlapHeap( OverlapHeapp overlap_heap );

/* Function:
     GetOverlapObject
   Description:
     wrapper around GetGenericHeapItem
     gets an OverlapObject from an Overlap heap
   Return Value:
     valid pointer to OverlapObject if ok
   Parameters:
     OverlapHeapp overlap_heap: pointer to existing overlap heap
*/
OverlapObjectp GetOverlapObject( OverlapHeapp overlap_heap );

/* Function:
     AllocateOverlapArray
   Description:
     wrapper around AllocateGenericArray
     initial allocation of an overlap array
   Return Value:
     valid OverlapArrayp if ok
   Parameters:
     uint32 num_overlaps: number of Overlap objects to allocate
*/
OverlapArrayp AllocateOverlapArray( uint32 num_overlaps );

/* Function:
     FreeOverlapArray
   Description:
     wrapper around FreeGenericArray
     Frees a OverlapArray array & resets structure values
   Return Value:
     none
   Parameters:
     OverlapArrayp overlaps:  pointer to OverlapArray structure
*/
void FreeOverlapArray( OverlapArrayp overlaps );

/* Function:
     GetNumFoundOverlaps
   Description:
     performs system commands to count number of OVL records in overlap file
   Return Value:
     non-zero number of overlaps in overlap file if ok
   Parameters:
     char * filename: name of overlap message file
*/
uint32 GetNumFoundOverlaps( char * filename );

/* Function:
     InsertionSortOverlap
   Description:
     inserts a new overlap object into a linked list of overlap objects in
     increasing order of frag_id
   Return Value:
     none
   Parameters:
     OverlapObjectp * start: pointer to first pointer in linked list
     OverlapObjectp overlap: pointer to overlap object to insert
*/
void InsertionSortOverlap( OverlapObjectp * start, OverlapObjectp overlap );

/*********************************************************************/
  
#endif // #ifndef AS_ORA_OVERLAPS_H


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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/AS_ORA_overlaps.c,v $
$Revision: 1.3 $
**********************************************************************/


/*********************************************************************/
// headers
// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// project headers
#include "AS_UTL_heap.h"
#include "AS_ORA_overlaps.h"
/*********************************************************************/


/*********************************************************************/
/* Function:
     AllocateOverlapHeap
   Description:
     wrapper around AllocateGenericHeap for overlap objects
   Return Value:
     valid OverlapHeapp pointer if ok
   Parameters:
     uint32 num_overlaps: number of overlaps to allocate initially
*/
OverlapHeapp AllocateOverlapHeap( uint32 num_overlaps )
{
  return( AllocateGenericHeap( num_overlaps, sizeof( OverlapObject ) ) );
}

/* Function:
     GrowOverlapHeap
   Description:
     wrapper around GrowGenericHeap to add memory to overlap heap
   Return Value:
     0 if ok
   Parameters:
     OverlapHeapp overlap_heap: pointer to existing overlap heap
*/
int GrowOverlapHeap( OverlapHeapp overlap_heap )
{
  return( GrowGenericHeap( overlap_heap ) );
}

/* Function:
     FreeOverlapHeap
   Description:
     wrapper around FreeGenericHeap - frees an overlap heap
   Return Value:
     none
   Parameters:
     OverlapHeapp overlap_heap: pointer to existing overlap heap
*/
void FreeOverlapHeap( OverlapHeapp overlap_heap )
{
  FreeGenericHeap( overlap_heap );
}

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
OverlapObjectp GetOverlapObject( OverlapHeapp overlap_heap )
{
  return( (OverlapObjectp) GetGenericHeapItem( overlap_heap ) );
}

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
OverlapArrayp AllocateOverlapArray( uint32 num_overlaps )
{
  return( AllocateGenericArray( num_overlaps, sizeof( OverlapObject ) ) );
}

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
void FreeOverlapArray( OverlapArrayp overlaps )
{
  FreeGenericArray( overlaps );
}

/* Function:
     GetNumFoundOverlaps
   Description:
     performs system commands to count number of OVL records in overlap file
   Return Value:
     non-zero number of overlaps in overlap file if ok
   Parameters:
     char * filename: name of overlap message file
*/
uint32 GetNumFoundOverlaps( char * filename )
{
  FILE * fp = fopen( filename, "r" );
  uint32 num_messages = 0;
  GenericMesg * gmesg = NULL;
  MesgReader reader;

  if( fp == NULL )
  {
    fprintf( stderr,
             "Failed to open overlap store %s for reading.\n",
             filename );
    return 1;
  }
  reader = InputFileType_AS( fp );

  while( reader( fp, &gmesg ) != EOF )
  {
    if( gmesg && gmesg->t == MESG_OVL )
    {
      num_messages++;
    }
  }

  // close up & return
  fclose( fp );
  return num_messages;
}

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
void InsertionSortOverlap( OverlapObjectp * start, OverlapObjectp overlap )
{
  OverlapObjectp curr = *start;
  OverlapObjectp prev = NULL;

  // if the list is empty
  if( curr == NULL )
  {
    *start = overlap;
    return;
  }
  
  while( curr && curr->frag_id < overlap->frag_id )
  {
    prev = curr;
    curr = curr->next;
  }

  // check for prepending
  if( prev == NULL )
  {
    overlap->next = curr;
    *start = curr->prev = overlap;
  }
  else
  {
    // check for appending
    if( curr == NULL )
    {
      overlap->prev = prev;
      prev->next = overlap;
    }
    else  // otherwise inserting between prev & curr
    {
      overlap->next = curr;
      overlap->prev = prev;
      prev->next = curr->prev = overlap;
    }
  }
}

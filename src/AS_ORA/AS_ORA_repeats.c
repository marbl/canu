
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/AS_ORA_repeats.c,v $
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
#include "AS_ORA_repeats.h"
/*********************************************************************/


/*********************************************************************/
/* Function:
     AllocateRepeatHeap
   Description:
     wrapper around AllocateGenericHeap to allocate a heap of RepeatObjects
   Return Value:
     valid RepeatHeap pointer if ok
   Parameters:
     uint32 num_repeats: number of repeats to allocate initially
*/
RepeatHeapp AllocateRepeatHeap( uint32 num_repeats )
{
  return( AllocateGenericHeap( num_repeats, sizeof( RepeatObject ) ) );
}

/* Function:
     GrowRepeatHeap
   Description:
     wrapper around GrowGenericHeap
     increases size of existing repeat heap
   Return Value:
     0 if ok
   Parameters:
     RepeatHeapp repeat_heap: pointer to repeat heap
*/
int GrowRepeatHeap( RepeatHeapp repeat_heap )
{
  return( GrowGenericHeap( repeat_heap ) );
}

/* Function:
     FreeRepeatHeap
   Description:
     wrapper around FreeGenericHeap
   Return Value:
     none
   Parameters:
     RepeatHeapp repeat_heap: pointer to repeat heap
*/
void FreeRepeatHeap( RepeatHeapp repeat_heap )
{
  FreeGenericHeap( repeat_heap );
}

/* Function:
     GetRepeatObject
   Description:
     wrapper around GetGenericHeapItem
     gets a RepeatObject from a heap
   Return Value:
     valid RepeatObject pointer if ok
   Parameters:
     RepeatHeapp repeat_heap: pointer to repeat heap
*/
RepeatObjectp GetRepeatObject( RepeatHeapp repeat_heap )
{
  return( (RepeatObjectp) GetGenericHeapItem( repeat_heap ) );
}

/* Function:
     AllocateRepeatArray
   Description:
     wrapper around AllocateGenericArray
     initial allocation of an repeat array
   Return Value:
     valid RepeatArrayp if ok
   Parameters:
     uint32 num_repeats: number of Repeat objects to allocate
*/
RepeatArrayp AllocateRepeatArray( uint32 num_repeats )
{
  return( AllocateGenericArray( num_repeats, sizeof( RepeatObject ) ) );
}

/* Function:
     FreeRepeatArray
   Description:
     wrapper around FreeGenericArray
     frees a RepeatArray array & resets structure values
   Return Value:
     none
   Parameters:
     RepeatArrayp repeats:  pointer to RepeatArray structure
*/
void FreeRepeatArray( RepeatArrayp repeats )
{
  FreeGenericArray( repeats );
}

/* Function:
     InsertionSortRepeat
   Description:
     inserts a new repeat object into a linked list of repeat objects in
     increasing order of repeat_id
   Return Value:
     none
   Parameters:
     RepeatObjectp * start: pointer to first pointer in linked list
     RepeatObjectp repeat: pointer to repeat object to insert
*/
void InsertionSortRepeat( RepeatObjectp * start, RepeatObjectp repeat )
{
  RepeatObjectp curr = *start;
  RepeatObjectp prev = NULL;

  // if the list is empty
  if( curr == NULL )
  {
    *start = repeat;
    return;
  }
  
  while( curr && curr->repeat_id < repeat->repeat_id )
  {
    prev = curr;
    curr = curr->next;
  }

  // check for prepending
  if( prev == NULL )
  {
    repeat->next = curr;
    *start = repeat;
  }
  else
  {
    // check for appending
    if( curr == NULL )
    {
      prev->next = repeat;
    }
    else  // otherwise inserting between prev & curr
    {
      repeat->next = curr;
      prev->next = repeat;
    }
  }
}

/* Function:
     GetRepeatID
   Description:
     Returns the first or next annotation character from a repeat bit vector
   Return Value:
     char: the character representing the repeat
   Parameters:
     RepeatID id: the repeat id bit vector
*/
char GetRepeatID( RepeatID id )
{
  char ret_val = 'A';
  while( id && !(id & 1) )
  {
    id >>= 1;
    ret_val ++;
  }
  return ret_val;
}

/* Function:
     GetRepeatString
   Description:
     Produces a string listing all repeats in a repeat id bit vector
   Return Value:
     none
   Parameters:
     RepeatID id: the repeat id bit vector
     char * string: string to hold result
*/
void GetRepeatString( RepeatID id, char * string )
{
  char repeat_char = 'A';
  int str_index = 0;

  string[0] = '\0';
  // get the first
  while( id && !(id & 1) )
  {
    id >>= 1;
    repeat_char ++;
  }

  if( id != 0 )
  {
    string[str_index++] = repeat_char;
    id >>= 1;
    repeat_char++;
    while( id != 0 )
    {
      if( id & 1 )
      {
        string[str_index++] = ',';
        string[str_index++] = ' ';
        string[str_index++] = repeat_char;
      }
      id >>= 1;
      repeat_char++;
    }
  }
  string[str_index] = '\0';
}

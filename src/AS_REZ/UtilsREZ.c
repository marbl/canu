
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

 Module: UtilsREZ.c 

 Description: Contains small utility functions

 Programmer: K. Reinert
             S. Lonardi (stelo@cs.purdue.edu)

 Assumptions: none

**********************************************************************/

static char CM_ID[] = "$Id: UtilsREZ.c,v 1.2 2004-09-23 20:25:28 mcschatz Exp $";

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "UtilsREZ.h"
#include "DataTypesREZ.h"
#include "dpc_CNS.h"

// -------------------------------------
// bit manipulation (for chunk subgraph)
// -------------------------------------

void Clear_All_Path_Bit(chunk_subgraph * s) {
  //
  // clear all path bits
  //
  int
    i;

  assert(s != NULL);
  for (i = 0; i < s->size; i++)
    s->node[i].path_bit = FALSE;
}



void Clear_Path_Bit(chunk_subgraph * s, int32 cid) {
  //
  // clear the path bit of <cid>
  //
  assert(s != NULL);
  assert(cid < s->max);
  assert(s->table[cid] != NULL);
  s->table[cid]->path_bit = FALSE;  
}



void Set_Path_Bit(chunk_subgraph * s, int32 cid) {
  //
  // set the path bit of <cid>
  //
  assert(s != NULL);
  assert(cid < s->max);
  assert(s->table[cid] != NULL);
  s->table[cid]->path_bit = TRUE;  
}



void Clear_All_Visited_Bit(chunk_subgraph * s) {
  //
  // clear all visited bits
  //
  int
    i;

  assert(s != NULL);
  for (i = 0; i < s->size; i++)
    s->node[i].visited = FALSE;
}



void Clear_Visited_Bit(chunk_subgraph * s, int32 cid) {
  //
  // clear the visited bit of <cid>
  //
  assert(s != NULL);
  assert(cid < s->max);
  assert(s->table[cid] != NULL);
  s->table[cid]->visited = FALSE;  
}



void Set_Visited_Bit(chunk_subgraph * s, int32 cid) {
  //
  // set the visited bit of <cid>
  //
  assert(s != NULL);
  assert(cid < s->max);
  assert(s->table[cid] != NULL);
  s->table[cid]->visited = TRUE;
}



void Clear_All_Done_Bit(chunk_subgraph * s) {
  //
  // clear all done bits
  //
  int
    i;

  assert(s != NULL);
  for (i = 0; i < s->size; i++)
    s->node[i].done = FALSE;
}



void Set_Done_Bit(chunk_subgraph * s, int32 cid) {
  //
  // set the done bit of <cid>
  //
  assert(s != NULL);
  assert(cid < s->max);
  assert(s->table[cid] != NULL);
  s->table[cid]->done = TRUE;
}


// ---------------------
// <nodes_stack> methods
// ---------------------


void Push_Node(nodes_stack * s, int v) {
  //
  // push a node
  //
  assert(s != NULL);
  assert(s->top < s->max_size);
  s->nodes[s->top] = v;
  (s->top)++;
}



int Pop_Node(nodes_stack * s) {
  //
  // pop a node
  //
  assert(s != NULL);
  assert(s->top);
  return s->nodes[--(s->top)];
}



int Top(nodes_stack * s) {
  //
  // get the value of Top without popping
  //
  assert(s != NULL);
  assert(s->top);
  return s->nodes[(s->top - 1)];
}



nodes_stack * Create_Stack(int no_elements) {
  //
  // crate a stack of <no_element> size
  //
  nodes_stack
    * s = (nodes_stack *)safe_calloc(1, sizeof(nodes_stack));

  assert(no_elements);
  s->top = 0;
  s->max_size = no_elements;
  s->nodes = (int *)safe_calloc(no_elements, sizeof(int));

  return s;
}



void Free_Stack(nodes_stack * s) {
  //
  // do you want the memory back?
  //
  assert(s != NULL);
  free(s->nodes);
  free(s);
}

// ------------
// CIEdge stuff
// ------------

char * Orientation_As_String (ChunkOrientationType orient) {
  //
  //  Return string equivalent of orient 
  //
  switch  (orient) {
  case  AB_AB :
    return  "AB_AB";
  case  AB_BA :
    return  "AB_BA";
  case  BA_BA :
    return  "BA_BA";
  case  BA_AB :
    return  "BA_AB";
  default :
    return  "*???*";
  }
}



int or2num(ChunkOrientationType o) {
  //
  // convert a <ChunkOrientationType> to a number
  // between 0 and 3
  //
  switch (o) {
  case AB_AB :
    return OR2NUM_AB_AB;
  case AB_BA :
    return OR2NUM_AB_BA;
  case BA_AB :
    return OR2NUM_BA_AB;
  case BA_BA :
    return OR2NUM_BA_BA;
  default :
    assert(0);
  }
  return -1;
}

/*--------------------------------------------------------------------*/
/* Interval Math */
/*--------------------------------------------------------------------*/


int Intersection (LengthT * a, LengthT * b) {
  //  Return the number of bases by which the closed interval  [mean - 3stdDev, mean + 3stdDev]
  //  intersects the closed interval  [mean - 3stdDev, mean + 3stdDev]
  return Interval_Intersection (a->mean - 3.0 * sqrt(a->variance), a->mean + 3.0 * sqrt(a->variance),
				b->mean - 3.0 * sqrt(b->variance), b->mean + 3.0 * sqrt(b->variance)); 
}



int Interval_Intersection (int a, int b, int c, int d) {
  //  Return the number of bases by which the closed interval [a, b]
  //  intersects the closed interval [c, d] 
  if  (d < a || b < c)
    return  0;
  else
    return  1 + Min_int (b, d) - Max_int (a, c);
}



double Max_double (double a, double b) {
  //  Return the larger of a and  b 
  if  (a < b)
    return  b;
  else
    return  a;
}



double  Min_double (double a, double b) {
  //  Return the smaller of a and  b 
  if  (a < b)
    return  a;
  else
    return  b;
}


float Max_float (float a, float b) {
  //  Return the larger of a and  b 
  if  (a < b)
    return  b;
  else
    return  a;
}



float  Min_float(float a, float b) {
  //  Return the smaller of a and  b 
  if  (a < b)
    return  a;
  else
    return  b;
}



int Max_int (int a, int b) {
  //  Return the larger of a and  b 
  if  (a < b)
    return  b;
  else
    return  a;
}



int  Min_int (int a, int b) {
  //  Return the smaller of a and  b 
  if  (a < b)
    return  a;
  else
    return  b;
}


/*--------------------------------------------------------------------*/
/* Error Handling */
/*--------------------------------------------------------------------*/

static char ErrorString[NumOfErrorsREZ][30] = 
{
  "FILE ERROR",
  "MEMORY ERROR",
  "PRECONDITION ERROR"
};


void error(ErrorCodeREZ err, const char* errorMesg,  ExitStatusREZ ex, 
	   const char* file, int line)
{
  /************************************************************************/ 
  /* The function checks whether the error code err is valid. If not it   */
  /* issues an error message. Otherwise it issues the error message given */
  /* by errorMesg and exits with status ex.                               */
  /************************************************************************/

  int intErr = (int)err;
  if( intErr < 0 || intErr >= NumOfErrorsREZ )
    {
    fprintf(stderr,"ERROR in error function : Invalid ErrorCode %d ",intErr);
    exit(EXIT_FAILURE_REZ);
  }
  fprintf(stderr,"%s in file %s, line %d : %s\n",
	  ErrorString[intErr],file,line,errorMesg);
  
  exit(ex);
}



/*--------------------------------------------------------------------*/
/*  File Handling routines */
/*--------------------------------------------------------------------*/


FILE*  file_open(const char* fileName, const char* mode)
{
  /*****************************************************************/
  /* Open  Filename  in  Mode  and return a pointer to its control */
  /* block.  If fail, print a message and exit.                    */
  /*****************************************************************/
   FILE* fp;

   fp = fopen (fileName,mode);
   if(fp == NULL)
     {
        char dummy[40];
	sprintf(dummy,"Could not open file  %s in mode %s\n",fileName, mode);
        error(FILE_ERROR_REZ,dummy,EXIT_FAILURE_REZ,__FILE__,__LINE__);
   }
   return  fp;
}


FileStatusREZ file_exists (const char * fileName)
{
  /*****************************************************************/
  /* Test for filename existence.                                  */
  /*****************************************************************/
  FILE* fp;

  fp = fopen (fileName, "r+");
  if(fp)
    {
      fclose(fp);
      return FILE_EXISTS_REZ;
    }
  
  return FILE_EXISTS_NOT_REZ;
}


/*--------------------------------------------------------------------*/
/*  Memory Managment routines */
/*--------------------------------------------------------------------*/


void *safe_calloc(size_t num, size_t len)
{
  /*****************************************************************/
  /* Allocate and return a pointer to an array of  num  elements of*/
  /* len  bytes each.  All are set to 0.  Exit if fai.             */
  /*****************************************************************/
  void  *p;

   p = calloc (num, len);
   if  (p == NULL)
     {
       char dummy[100];
       sprintf(dummy,"Could not calloc memory (%d number * %d bytes) \n",
               (int) num,(int) len);
       assert(0);
       error(MEMORY_ERROR_REZ,dummy,EXIT_FAILURE_REZ,__FILE__,__LINE__);
     }
   return  p;
}



void *safe_malloc(size_t len)
{
  /*****************************************************************/
  /* Allocate and return a pointer to len bytes of memory.         */
  /* Len  bytes each.  Exit if fail.                               */
  /*****************************************************************/
  void  *p;

  p = malloc (len);
  if(p == NULL)
    {
      char dummy[40];
      sprintf(dummy,"Could not malloc memory (%d bytes) \n",(int) len);
      error(MEMORY_ERROR_REZ,dummy,EXIT_FAILURE_REZ,__FILE__,__LINE__);
    }
  
  return  p;
}



void *safe_realloc(void *q, size_t len)
{
  /*****************************************************************/
  /* Reallocate memory for q to len  bytes and return a pointer    */
  /* to the new memory.  Exit if fail.                             */     
  /*****************************************************************/
  void  *p;

  if  (q == NULL)
      p = malloc (len);        // to prevent some compilers' warnings
    else
      p = realloc (q, len);

  if  (p == NULL)
    {
      char dummy[100];
      sprintf(dummy,"Could not realloc memory (%d bytes) \n",(int) len);
      error(MEMORY_ERROR_REZ,dummy,EXIT_FAILURE_REZ,__FILE__,__LINE__);
    }
  
  return  p;
}




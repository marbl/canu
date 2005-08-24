
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
   CVS_ID:  $Id: UtilsREZ.h,v 1.5 2005-08-24 07:47:15 brianwalenz Exp $
 *********************************************************************/
#ifndef UTILSREZ_H
#define UTILSREZ_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "AS_CGW_dataTypes.h"
#include "DataTypesREZ.h"
#include "dpc_CNS.h"

#define OR2NUM_AB_AB         0
#define OR2NUM_AB_BA         1
#define OR2NUM_BA_AB         2
#define OR2NUM_BA_BA         3

// ----------------
// bit manipulation
// ----------------

void Clear_All_Path_Bit(chunk_subgraph *);

void Clear_All_Visited_Bit(chunk_subgraph *);

void Clear_All_Done_Bit(chunk_subgraph *);

void Set_Path_Bit(chunk_subgraph *, int32);

void Set_Visited_Bit(chunk_subgraph *, int32);

void Set_Done_Bit(chunk_subgraph *, int32);

void Clear_Path_Bit(chunk_subgraph *, int32);

void Clear_Visited_Bit(chunk_subgraph *, int32);

// -----
// stack
// -----

nodes_stack * Create_Stack(int);

void Push_Node(nodes_stack *, int);

int Top(nodes_stack *);

int Pop_Node(nodes_stack *);

void Free_Stack(nodes_stack *);

// ------------
// CIEdge stuff
// ------------

char * Orientation_As_String (ChunkOrientationType);

int or2num(ChunkOrientationType);

/*--------------------------------------------------------------------*/
/* Interval Math */
/*--------------------------------------------------------------------*/

int Intersection (LengthT *, LengthT *);
int Interval_Intersection (int, int, int, int);
double Max_double (double, double);
double Min_double (double, double);
float Max_float (float, float);
float Min_float (float, float);
int Max_int (int, int);
int Min_int (int, int);

/*---------------------------------------------------------*/
/* #defines */
/*---------------------------------------------------------*/

#define NumOfErrorsREZ 3

/*---------------------------------------------------------*/
/* typedefs */
/*---------------------------------------------------------*/

/* this typedef gives a list of possible error codes
   that index into a table used by the error Handler 
   Edit this table and accordingly update the number
   of error messages as well as the actual error messages
   in the array ErrorString in the file UtilsREZ.c   */

typedef enum {
  FILE_ERROR_REZ,
  MEMORY_ERROR_REZ,
  PRECONDITION_ERROR_REZ
} ErrorCodeREZ; 

typedef enum {
  EXIT_SUCCESS_REZ = 0,
  EXIT_FAILURE_REZ = 1
} ExitStatusREZ;


typedef enum {
  FILE_EXISTS_REZ     = 0,
  FILE_EXISTS_NOT_REZ = 1
} FileStatusREZ;



/*---------------------------------------------------------*/
/* functions */
/*---------------------------------------------------------*/


void error(ErrorCodeREZ err, const char* errorMesg,  ExitStatusREZ ex, 
	   const char* file, int line);
/************************************************************************/ 
/* The function checks whether the error code err is valid. If not it   */
/* issues an error message. Otherwise it issues the error message given */
/* by errorMesg and exits with status ex.                               */
/************************************************************************/


FileStatusREZ file_exists (const char* filename);
/*****************************************************************/
/* Test for filename existence. Returns AS_TER_FILE_EXISTS if    */
/* Filename exists and AS_TER_FILE_EXISTS_NOT o.w.               */
/*****************************************************************/


FILE* file_open(const char* filename, const char* mode);
/*****************************************************************/
/* Open  Filename  in  Mode  and return a pointer to its control */
/* block.  If it fails it issues an error message and exits.     */
/*****************************************************************/


#endif









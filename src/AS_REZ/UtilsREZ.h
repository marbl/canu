
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

#ifndef UTILSREZ_H
#define UTILSREZ_H

static const char *rcsid_UTILSREZ_H = "$Id: UtilsREZ.h,v 1.8 2008-10-08 22:03:00 brianwalenz Exp $";

#include "DataTypesREZ.h"

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

//--------------
// Interval Math
//--------------

int Intersection (LengthT *, LengthT *);
int Interval_Intersection (int, int, int, int);

//-----
// misc
//-----

FILE* file_open(const char* filename, const char* mode);

#endif









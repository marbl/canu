
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
   Module:  AS_OVL
   Description:  Assembly Overlap Module--computes overlaps between
      pairs of DNA strings.
 *********************************************************************/

/* RCS info
 * $Id: MatchListOVL.h,v 1.1.1.1 2004-04-14 13:52:42 catmandew Exp $
 * $Revision: 1.1.1.1 $
*/


#ifndef MATCH_LIST_OVL_H
#define MATCH_LIST_OVL_H

// Component:
//
//   MatchListOVL.h
//
//   Last revised:  10 July 2001
//
// Description:
// 
//   These routines implement the data structure that holds matches
//   between fragments that could potentially overlap.  A match is
//   an exact-match substring between two fragments.
// 
// Design:
// 
// 
// Limitations:
// 
// 
// Status:
// 
// 
// Architecture and Dependencies:
// 
//   MatchListOVL.h   This file
//


typedef  struct Match_Node
  {
   int32  Offset;
     // To start of exact match in hash-table frag
   int32  Len;
     // Of exact match
   int32  Start;
     // Of exact match in current (streaming) frag
   int32  Next;
     // Subscript of next match in list
  }  Match_Node_t;

typedef  struct Match_Set
  {
   int32  first;
     // Subscript of first match in set
   Match_Node_t  * match_node_space;
     // Array to which subscripts refer
  }  Match_Set_t;


int  Is_Empty
    (const Match_Set_t * matches);


#endif

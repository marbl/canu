
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
 * $Id: AS_CGB_edgemate.h,v 1.3 2005-03-22 19:02:07 jason_miller Exp $
 *
 * Module: AS_CGB_edgemate.h
 *
 * Description: These routines find and access the mate directed edge
 * for a given edge of an overlap.
 *
 * Assumptions: 
 *
 * Author: Clark Mobarry
 *
 *******************************************************************/

#ifndef AS_CGB_EDGEMATE_INCLUDE
#define AS_CGB_EDGEMATE_INCLUDE

#define AS_CGB_EDGE_NOT_FOUND (~((IntEdge_ID)(0)))

void reflect_Aedge( Aedge *new_edge, Aedge *old_edge);
/*  fragment overlaps:
    The overlapper connects two fragment-ends in an overlap
    relationship:
    
    A    ---------------->
    B          -------------->
    
    The direction mate edge preserves the which fragment-ends are
    in the overlap:
    
    B^c  <----------------
    A^c       <---------------- 
*/


void granger_Aedge( Aedge *new_edge, Aedge *old_edge);
/*
  To-contained edges have ahg>0
  A    ---------------->
  B          ------->...
  
  So create the granger mate edge:
  A^c  <----------------
  B^c     <-------......

  
  From-contained edges have bhg <0
  A    ...------->      
  B    ---------------->
  
  So create the granger mate edge:
  A^c  .......<-------
  B^c  <-----------------
*/


typedef int (*QsortCompare)(const void *, const void *);

void set_compare_edge_function( QsortCompare compare_edge);
//  int (*compare_edge)(const void *,const void *)

QsortCompare get_compare_edge_function(void);

int verify_that_the_edges_are_in_order(Tedge edges[]);

IntEdge_ID find_overlap_edge_mate
(/* Input Only */
 const Tfragment frags[], 
 const Tedge edges[],
 const IntEdge_ID ie0
);

void fix_overlap_edge_mate
(/* Input Only */
 const Tfragment frags[], 
 Tedge edges[],
 const IntEdge_ID ie0);

void append_the_edge_mates
(
 Tfragment frags[],
 Tedge edges[],
 TIntEdge_ID * next_edge_obj 
);

IntEdge_ID check_symmetry_of_the_edge_mates
(
 Tfragment frags[],
 Tedge edges[],
 TIntEdge_ID * next_edge_obj
);

#endif // AS_CGB_EDGEMATE_INCLUDE


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

        Module:  BccREZ.h

   Description:  Declaration of common data types used in BccREZ.c

    Programmer:  S. Lonardi (stelo@cs.purdue.edu)

       Written:  30 June 99

 **********************************************************************/

/*********************************************************************
   CVS_ID: $Id: BccREZ.h,v 1.2 2004-09-23 20:25:26 mcschatz Exp $
 *********************************************************************/

#ifndef BCC_REZ_H
#define BCC_REZ_H

#define MAX_LNK_CHUNKS      30
#define BCC_NUM_COLOURS     27
#define BCC_MIN_SIZE         3
#define DELTA_THRE           5.0 // admissible error in the placement
#define MIN_UNIQUES_PROP     0.0 // minimum proportion in a bcc

#define CHECK_COMPONENTS     1   // set to 1 to run a test on the edges of each component
#define CREATE_DOT_FILE      1   // set to 1 to produce a .dot file for really bad components

// ------
// protos
// ------

int Scan_Components(chunk_subgraph *,
		    bcc_array *,
		    float (*)(CIEdgeT *, int32),
		    int (* filter)(CIEdgeT *));

//
// do a DFS of the graph s
// avoid edges such that filter(e) == FALSE
// 
void DFS(chunk_subgraph *,
	 int (*)(CIEdgeT *));

//
// do a DFS of the subgraph from the <this> node
//
// assign start_time dfs, end_time, bcc_id (see CLR)
//
// avoid edges such that filter(e) == FALSE
//
void DFS_Visit(chunk_subgraph *,
               int32,
               int32,
               nodes_stack *, int,
               int (*)(CIEdgeT *));

//
// return a bcc_array that can be indexed by the component
// number and get all the chunks that re in that component
//
// note that one chunk could be in more than one component
// (articulation point)
//
// (the size of the array is <bcc_time_stamp>)
//
bcc_array * Compute_Bcc_Array(chunk_subgraph *);

void Print_Cam_BCC(bcc_array *,
		   chunk_subgraph *);

void Print_Cam_BCC_Links(bcc_array *);

int Edge_Feasible_Bcc(CIEdgeT *);

#endif


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
/* 	$Id: UnionFind_AS.h,v 1.4 2005-03-22 19:49:30 jason_miller Exp $	 */
#ifndef UNION_FIND_H
#define UNION_FIND_H

/* Union Find code based on Aho Hopcroft Ullman pg 131

  Conceptually this code is easy to use.

  First, create a UFData object for the problem you are solving, where
  problemSize is the number of nodes in your graph:

      UFDataT *UFData = UFCreateSets(<problemSize>); 

  Second, associate each set with one of your nodes:

      for(i = 0; i < <problemSize>; i++){
          UFSetT *set = UFGetSet(UFData,i);
          node = <Get your graph node i>;
	  set->data = (void *)node;  // associate your node with the set
      }

  Third, for each edge in your graph, perform a Union operation:

     for(edge = <iterate over all your graph edges>){
       int setA = <setID of the 1st node on which this edge is incident>;
       int setB = <setID of the 2nd node on which this edge is incident>;

       UFUnion(UFData, setA, setB);
     }


  Fourth, renumber the sets to achieve a dense encoding:

     int numSets = UFRenumberSets(UFData);


  Lastly, mark your nodes with their setid:

      for(i = 0; i < <problemSize>; i++){
          UFSetT *set = UFGetSet(UFData,i);
          node = <Get your graph node i>;
	  node->setID = set->component;
      }


  Congratulations, you're done.
     UFFreeSets(UFData);

 */

/* UFSetT:
   Data structure for a set in the Union Find Algorithm.
*/
typedef struct {
  void *data;        // Client uses this 
  int count;         // Number of elements in a set. Zero if not a definition
  int parent;        // Used for set tree structure
  int component;     // Final numbering of this set
  int isDefinition;  // used by renumber components only
}UFSetT;


/* UFDataT
   Data allocated to do a Union Find on numSets sets
*/
#define MAX_UNIONFIND_SIZE_AS 100000
typedef struct {
  int    numSets;
  UFSetT sets[MAX_UNIONFIND_SIZE_AS];  /* use a large number so bounds checking doesn't catch us
					  in practice, allocate less */
}UFDataT;

/* UFGetSet
   Range checking accessor to get data for set setID
*/
static UFSetT *UFGetSet(UFDataT *data, int setID){
  if(setID < data->numSets)
    return data->sets + setID;
  //else
    return NULL;
}

// Create data structures to do UF on numSets sets
//
UFDataT *UFCreateSets(int numSets);

// Free resources
//
void UFFreeSets(UFDataT *data);


/* UFFind
   Find the root setID for set index setID
   This employs the path compression algorithm.
*/
int UFFind(UFDataT *data, int setID);

/* UFUnion
   Performs the Union operator between the two sets 
*/
void UFUnion(UFDataT *data, int setA, int setB);

/* UFRenumberSets
   Renumber the sets so that they are densely encoded.
   Returns the number of sets.
*/
int UFRenumberSets(UFDataT *data);



#endif

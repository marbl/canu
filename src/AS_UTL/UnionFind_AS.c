
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
static char CM_ID[] = "$Id: UnionFind_AS.c,v 1.1.1.1 2004-04-14 13:53:51 catmandew Exp $";

#include <stdlib.h>
#include "UnionFind_AS.h"
#include <stdio.h>

// Union Find code based on Aho Hopcroft Ullman pg 131

static int isRoot(UFSetT *set){
  return set->count != 0;
}

static UFSetT *getParent(UFDataT *data, UFSetT *set){
  return data->sets + set->parent;
}

static void UFMark(UFDataT *data, int setID, int component);
static int UFFindRoot(UFDataT *data, int setID);


// AHU pg 4.19 
// A tree based Union-Find with path compression.
//
void UFMark(UFDataT *data, int setID, int component){
  UFSetT *parent, *child;
  UFSetT *set = data->sets + setID;
  int end;

  parent = getParent(data, set);
  child = set;
  end = 0;
  while(!end){
    end = isRoot(child);
    child->parent = component;
    child = parent;
    parent = getParent(data, parent);
  }
}


static  int UFFindRoot(UFDataT *data, int setID){
  UFSetT *set = data->sets + setID;
  UFSetT *parent;

    for(parent = set; 
	!isRoot(parent);
	parent = getParent(data, parent));
      
      return parent->parent;
}

int UFFind(UFDataT *data, int setId){
  int root = UFFindRoot(data, setId);
  UFMark(data, setId, root); // Path compression

    return root;
}


void UFUnion(UFDataT *data, int setA, int setB){
  UFSetT *rootA = data->sets + UFFind(data,setA);
  UFSetT *rootB = data->sets + UFFind(data,setB);
  UFSetT *smaller, *larger;

  // Don't do anything, since the two sets are already merged!
  if(rootA == rootB)
    return;

  // Always merge the smaller set into the larger one.
  if(rootA->count <= rootB->count){
    smaller = rootA;
    larger = rootB;
  }else{
    smaller = rootB;
    larger = rootA;
  }

  smaller->parent = larger->parent;
  larger->count += smaller->count;
  smaller->count = 0;

  return;
}



// Renumber components in range 0-numComponents
//
int UFRenumberSets(UFDataT *data){
  int i;
  int components = 0;
  for(i = 0; i < data->numSets; i++){
    UFSetT *compI = data->sets + i;
    if(isRoot(compI)){
      compI->component = components++;
    }
  }
  for(i = 0; i < data->numSets; i++){
    UFSetT *compI = data->sets + i;
      int newComponent = (data->sets + UFFindRoot(data, i))->component;
      if(!isRoot(compI)){
	compI->component = newComponent;
      }
  }
  return components;
}


//
UFDataT *UFCreateSets(int numSets){
  UFDataT *data = (UFDataT *)malloc(sizeof(UFDataT) + numSets * sizeof(UFSetT)); // This is a bit too large...OK
  int i;
  data->numSets = numSets;
  for(i = 0; i < numSets; i++){
    UFSetT *set = data->sets + i;
    set->data = NULL;
    set->component = i;
    set->parent = i;
    set->count = 1;
  }
  return data;
}

 void UFFreeSets(UFDataT *data){
   free(data);
 }

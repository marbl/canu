
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
static char CM_ID[] = "$Id: PartitionsREZ.c,v 1.4 2005-03-22 19:49:25 jason_miller Exp $";

#include "AS_UTL_skiplist.h"
#include "UtilsREZ.h"
#include "ScoreREZ.h"
#include "PartitionsREZ.h"


SL_DEF(Partition_t)


/* functions to allocate and free a partition of size l */
Partition_t *allocate_partition(int l)
{
  Partition_t* p = (Partition_t*) safe_malloc(sizeof(Partition_t));
  p->part        = (int*) safe_calloc(sizeof(int),l);
  p->groups      = 1;
  p->len         = l;
  return p;
}


void free_partition(Partition_t *p)
{
  free(p->part);
  free(p);
}


/* get_induced_partition looks at a column of an alignment and returns
   the partition induced by tat colmn */

Partition_t* get_induced_partition(int col, Alignment_t *a)
{
  int i;
  int rows = a->rows;
  int g = 0;
  Partition_t *par = allocate_partition(rows);

  int group_dash = -1;
  int group_A    = -1;
  int group_C    = -1;
  int group_G    = -1;
  int group_T    = -1;
  
  for(i=0; i<rows; i++)
    {
      switch(a->ali[col][i]){
      case '-' :
	if( group_dash < 0 )
	  group_dash = g++;
	par->part[i] = group_dash;
	break;	
      case 'A' :
	if( group_A < 0 )
	  group_A = g++;
	par->part[i] = group_A;
	break;
      case 'C' :
	if( group_C < 0 )
	  group_C = g++;
	par->part[i] = group_C;
	break;
      case 'G' :
	if( group_G < 0 )
	  group_G = g++;
	par->part[i] = group_G;
	break;
      case 'T' :
	if( group_T < 0 )
	  group_T = g++;
	par->part[i] = group_T;
	break;
      }
    }
  par->groups = g;
  par->len    = rows;

  return par;
}



/* The function combine_partitions takes as argument two partitions
   of the same length and combines them to a 'finer' partition which 
   is returned */

Partition_t *combine_partitions(Partition_t* part1, Partition_t* part2)
{
  int i,j,k;
  int rows = part1->len;
  Partition_t *ret;
  int g;

  int newGroup[rows][rows]; // mapping of the group combination in part1
                            //      and part2 to the new groups

  /* we only combine partition of the same length */
  assert(part1->len == part2->len);
  ret =  allocate_partition(rows);

  for(j=0; j<rows; j++)
    for(k=0; k<rows; k++)
      newGroup[j][k] = -1;


  g = 0; // the current number of new groups
 
  for(i=0; i<rows; i++)
    {         
      if( newGroup[part1->part[i]][part2->part[i]] < 0)
	newGroup[part1->part[i]][part2->part[i]] = g++;

      ret->part[i] = newGroup[part1->part[i]][part2->part[i]];
    }
  ret->groups = g;
  ret->len    = rows;
 
  return ret;  
}




SL_TYPE(Partition_t) *guess_partition(Alignment_t *a, int ct)
{
  SL_TYPE(Partition_t) *sl = CreateSL_Partition_t(TRUE,(SLF) free_partition);
  int i;
  int cols = a->cols;
  Partition_t *si = allocate_partition(a->rows);

  for(i=0; i<cols; i++)
    if(col_contributing(a,i,ct))
      {
	Partition_t *par = get_induced_partition(i,a);
	InsertSL_Partition_t(value(par,a,ct),par,sl);
       }

  /* In addition we insert the singleton partition */	
  InsertSL_Partition_t(value(si,a,ct),si,sl);


  /* Now we iterate over the columns in the skiplist */
  if( GetNumSL_Partition_t(sl) > 1 )
    {
      int improved;
      do{
	int i,j;
	int noOfCurrent = GetNumSL_Partition_t(sl);
	sl_item it = MinSL_Partition_t(sl);
	keyType currentMin = it->key;
	sl_item current[noOfCurrent];
	
	for(i=0; i<noOfCurrent; i++){
	  current[i] = it;
	  it = it->succ;
	}

	improved = FALSE;

	for(i=0; i<noOfCurrent-1; i++){
	  for(j=i+1; j<noOfCurrent; j++){
	    Partition_t *comb = combine_partitions(current[i]->value,current[j]->value);
	    double combValue = value(comb,a,ct);
	    if( combValue < current[i]->key && combValue < current[j]->key)
	      InsertSL_Partition_t(combValue,comb,sl);
	    else
	      free_partition(comb);
	  }
	}

	if( MinSL_Partition_t(sl)->key < currentMin ){
	  improved = TRUE;
	  //	  printf("Improved to %lf \n",MinSL_int(sl)->key);
	}

      }while( improved == TRUE );
    }
  return sl;	 
}






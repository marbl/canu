
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
#include "AS_global.h"
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "AS_UTL_HashCommon.h"
#include "AS_UTL_PHash.h"

int main(int argc, char **argv){
  int r;
  PHashTable_AS *hashtable = CreatePHashTable_AS(10, NULL); /* In Memory */
  PHashTable_AS *hashtable1 = CreatePHashTable_AS(10, "harry"); /* Memory Mapped File */
  PHashValue_AS value;
  PHashTable_Iterator_AS iterator;
  int rvalue;
  uint64 input;
  // uint64 output;
  int i;
  char nameSpace = 'A';
  char n;

  for(i = 1; i <= 1000; i++){
    input = i*1000 + i * 100 + i * 10 + i;

    r = InsertInPHashTable_AS(&hashtable, nameSpace, input, &value, FALSE);
    assert(r == HASH_SUCCESS);
  }

  fprintf(stderr,"* 1st round of lookups \n");
  for(i = 1; i <= 1000; i++){
    input = i*1000 + i * 100 + i * 10 + i;
    rvalue = LookupInPHashTable_AS(hashtable, nameSpace, input, &value);
    assert(rvalue == HASH_SUCCESS &&
           value.IID == i);
  }

  fprintf(stderr,"* TestHashTable:  Before Concat \n");
  InitializePHashTable_Iterator_AS(hashtable, &iterator);
    while(HASH_SUCCESS == NextPHashTable_Iterator_AS(&iterator, &n, &input, &value)){
      fprintf(stderr,"\t* ns = %c key = " F_S64
	      " value = %d\n",n, input, value.IID);
    }
	      

    fprintf(stderr,"\n\n* BEFORE: TestHashTable hashtable numNodes = %d htbl1->numNodes = %d\n",
	    hashtable->numNodes, hashtable1->numNodes);

  ConcatPHashTable_AS(&hashtable1, hashtable);

    fprintf(stderr,"\n\n* AFTER: TestHashTable hashtable numNodes = %d htbl1->numNodes = %d\n",
	    hashtable->numNodes, hashtable1->numNodes);
  fprintf(stderr,"\n\n* TestHashTableAfter Concat HashTable1 \n");
  InitializePHashTable_Iterator_AS(hashtable1, &iterator);
    while(HASH_SUCCESS == NextPHashTable_Iterator_AS(&iterator, &n, &input, &value)){
      fprintf(stderr,"* ns = %c key = " 
	      F_S64 " value = %d\n",n, input, value.IID);
    }

  fprintf(stderr,"\n\n* TestHashTable After Concat HashTable \n");
  InitializePHashTable_Iterator_AS(hashtable, &iterator);
    while(HASH_SUCCESS == NextPHashTable_Iterator_AS(&iterator, &n, &input, &value)){
      fprintf(stderr,"\t* ns = %c key = "
	      F_S64 " value = %d\n",n, input, value.IID);
    }
	      
  ClosePHashTable_AS(hashtable);

  fprintf(stderr,"* 2nd round of lookups \n");
  for(i = 1; i <= 1000; i++){
    input = i*1000 + i * 100 + i * 10 + i;
    rvalue = LookupInPHashTable_AS(hashtable1, nameSpace, input, &value);
    if(rvalue != HASH_SUCCESS ||
       value.IID != i){
      fprintf(stderr,"* Lookup failure i = %d ns = %c input = "
	      F_S64 " rvalue = %d value.IID = %d\n",
	      i, nameSpace, input, rvalue, value.IID);
      exit(1);
    }
  }

  ClosePHashTable_AS(hashtable);
  ClosePHashTable_AS(hashtable1);

  return 0;
}


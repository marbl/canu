
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "AS_global.h"
#include "AS_UTL_Hash.h"

#define NUM_TESTS     1000
#define VAL_START   900000

int main(int argc, char **argv){
  int r;
  HashTable_AS *hashtable = CreateHashTable_uint64_AS(10);
  int i;
  uint64 inputs[NUM_TESTS];
  uint64 * output;

  for(i = 0; i < 1000; i++){
    inputs[i] = VAL_START + i;
    r = InsertInHashTable_AS(hashtable, (void *) &(inputs[i]),
                             sizeof(uint64), (void *) &(inputs[i]));
    fprintf(stdout,"Inserted " F_U64 " into hashtable with result %d\n",
	    inputs[i], r);
  }

  for(i = 0; i < 1000; i++){
    output = (uint64 *)LookupInHashTable_AS(hashtable, 
                                            (void *) &(inputs[i]),
                                            sizeof(uint64));
    fprintf(stdout,"Looked up " F_U64 " in Hashtable with result " F_U64 "\n",
            inputs[i], *output);
  }
  DeleteHashTable_AS(hashtable);

return 0;
}


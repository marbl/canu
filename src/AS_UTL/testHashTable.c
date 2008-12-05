
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

static const char *rcsid = "$Id: testHashTable.c,v 1.5 2008-12-05 19:06:12 brianwalenz Exp $";

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "AS_global.h"
#include "AS_UTL_Hash.h"

//  cc -o test -I.. -I. testHashTable.c AS_UTL_Hash.c AS_UTL_heap.c AS_UTL_alloc.c AS_UTL_fileIO.c -lm

#define NUM_ENTRIES  130000000

int
main(int argc, char **argv) {

  HashTable_AS *hashtable;
  uint64       *inputs;
  uint64        output;

  FILE *fp;

  int i;

  hashtable = CreateScalarHashTable_AS();
  inputs    = (uint64 *)safe_malloc(NUM_ENTRIES * sizeof(uint64));

  srand48(time(NULL));

  for (i=0; i<NUM_ENTRIES; i++){
    inputs[i]   = lrand48();
    inputs[i] <<= 32;
    inputs[i]  |= lrand48();

    InsertInHashTable_AS(hashtable, inputs[i], 0, inputs[i], 0);

    if ((i % 1000000) == 0)
      fprintf(stderr, "inserting %d.\n", i);
  }

  fprintf(stderr, "testing.\n");
  for (i=0; i<NUM_ENTRIES; i++) {
    if (LookupValueInHashTable_AS(hashtable, inputs[i], 0) != inputs[i]) {
      fprintf(stderr, "hash error for "F_U64"\n", inputs[i]);
    }
  }

  fprintf(stderr, "writing.\n");
  SaveHashTable_AS("test.hashtable", hashtable);

  fprintf(stderr, "deleting.\n");
  DeleteHashTable_AS(hashtable);

  fprintf(stderr, "reading.\n");
  hashtable = LoadUIDtoIIDHashTable_AS("test.hashtable");

  fprintf(stderr, "testing.\n");
  for (i=0; i<NUM_ENTRIES; i++) {
    if (LookupValueInHashTable_AS(hashtable, inputs[i], 0) != inputs[i]) {
      fprintf(stderr, "hash error for "F_U64"\n", inputs[i]);
    }
  }

  fprintf(stderr, "deleting.\n");
  DeleteHashTable_AS(hashtable);

  fprintf(stderr, "all done.\n");

  unlink("test.hashtable");
  safe_free(inputs);

  return(0);
}


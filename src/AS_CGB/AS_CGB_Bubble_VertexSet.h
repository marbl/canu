
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
#ifndef _AS_CGB_BUBBLE_LABELSET_H_
#define _AS_CGB_BUBBLE_LABELSET_H_

/*
 * CONSTANTS AND GLOBAL VARIABLES
 */

#define AS_CGB_BUBBLE_DEFAULT_SET_SIZE      10
#define AS_CGB_BUBBLE_DEFAULT_MAX_AGE       -1
#define AS_CGB_BUBBLE_DEFAULT_MAX_OUTDEGREE  4

extern int AS_CGB_BUBBLE_set_size_G;
extern int AS_CGB_BUBBLE_max_age_G;
extern int AS_CGB_BUBBLE_max_outdegree_G;

/*
 * VertexSet Type
 */

typedef struct BubVertexEntry {
  IntFragment_ID id;
  int age;
} BubVertexEntry;

typedef struct BubVertexSet {
  int numEntries;
  BubVertexEntry *entries;
} BubVertexSet;

typedef BubVertexSet * BubVertexSet_t;

/*
 * Operations on vertex sets.
 */

/*
 * IMPORTANT: Before any BVS objects are created or destroyed, the BVS system (specifically, memory allocation
 * subsystem) must be initialized.  It should later be destroyed to free memory.
 */
void BVS_sysInit(void);

void BVS_sysDone(void);

/* Initialize and deallocate BVS objects. */

void BVS_initialize(BubVertexSet_t bvs);
void BVS_destroy(BubVertexSet_t bvs);

/* Get a hash value from a BVS. */
uint64 BVS_hash(BubVertexSet_t bvs);

/* Tests to see if the first BVS is "less than" the second, using
   lexicographical ordering.  Returns -1 if less than, 0 if equal,
   or 1 if greater than. */
int BVS_compare(BubVertexSet_t bvs1, BubVertexSet_t bvs2);

/* Returns the number of entries in the BVS. */
int BVS_numEntries(BubVertexSet_t bvs);
int BVS_empty(BubVertexSet_t bvs);
int BVS_full(BubVertexSet_t bvs);

/* Computes the intersection of the first argument with the second,
   leaving the result in the first argument. */
void BVS_intersect(BubVertexSet_t bvs1, BubVertexSet_t bvs2);

/* Copies the first argument into the second, which is assumed to be
   an initialized, empty BVS. */
void BVS_copy(BubVertexSet_t src_bvs, BubVertexSet_t dst_bvs);

/* Adds amount inc to the age in each entry, and deletes those greater
   than the max age. */
void BVS_age(BubVertexSet_t bvs, int inc);

/* Inserts a new entry into the BVS with age 0.  If the BVS is full,
   the oldest entry is removed to make space. */
void BVS_insert(BubVertexSet_t bvs, IntFragment_ID new_entry);

/* Returns true if the second argument is a subset of the first. */
int BVS_subset(BubVertexSet_t bvs1, BubVertexSet_t bvs2);

/* Returns a BVS to its empty state. */
void BVS_makeEmpty(BubVertexSet_t bvs);

/* Prints a BVS vertex set to the specific file. */
void BVS_print(BubVertexSet_t bvs, FILE *of);

#endif

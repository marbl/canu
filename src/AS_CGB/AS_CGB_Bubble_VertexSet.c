
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

static char CM_ID[] = "$Id: AS_CGB_Bubble_VertexSet.c,v 1.9 2008-06-27 06:29:13 brianwalenz Exp $";

#include <string.h>
#include "AS_CGB_all.h"
#include "AS_CGB_methods.h"
#include "AS_CGB_Bubble_VertexSet.h"

// CMM: Is this an arbitrary size that needs to be revisited?
#define BVS_MEMSTACK_SIZE (2000 * CGB_MULTIPLIER)

int AS_CGB_BUBBLE_set_size_G = AS_CGB_BUBBLE_DEFAULT_SET_SIZE;
int AS_CGB_BUBBLE_max_age_G = AS_CGB_BUBBLE_DEFAULT_MAX_AGE;
int AS_CGB_BUBBLE_max_outdegree_G = AS_CGB_BUBBLE_DEFAULT_MAX_OUTDEGREE;

static BubVertexEntry **BVS_memstack_G;
static int BVS_memstack_top_G;

static BubVertexEntry *
BVS__alloc(void)
{
  if (BVS_memstack_top_G == 0) {
    for (BVS_memstack_top_G = 0; BVS_memstack_top_G < BVS_MEMSTACK_SIZE / 2; ++BVS_memstack_top_G) {
      BVS_memstack_G[BVS_memstack_top_G] = safe_malloc(sizeof(BubVertexEntry) * AS_CGB_BUBBLE_set_size_G );
    }
  }
  return (BVS_memstack_G[--BVS_memstack_top_G]);
}


static void
BVS__free(BubVertexEntry *block)
{
  if (!block)
    return;

  if (BVS_memstack_top_G < BVS_MEMSTACK_SIZE) {
    BVS_memstack_G[BVS_memstack_top_G++] = block;
  } else {
    safe_free(block);
  }
}



static void
BVS__removeIndex(BubVertexSet_t bvs, int idx)
{
  while ((idx + 1) < bvs->numEntries) {
    bvs->entries[idx] = bvs->entries[idx + 1];
    idx++;
  }

  bvs->numEntries--;
}


static void
BVS__removeOldest(BubVertexSet_t bvs)
{
  int max_idx;
  int max_age;
  int i;

  if (bvs->entries == 0)
    return;

  max_idx = 0;
  max_age = bvs->entries[0].age;

  for (i = 1; i < bvs->numEntries; ++i)
    if (bvs->entries[i].age > max_age) {
      max_age = bvs->entries[i].age;
      max_idx = i;
    }

  BVS__removeIndex(bvs, max_idx);
}


static void
BVS__insertWithAge(BubVertexSet_t bvs, IntFragment_ID id, int age)
{
  int ip = 0;
  int sp;

  if (!bvs->entries) {
    bvs->entries = BVS__alloc();
  }
  else if (BVS_full(bvs))
    BVS__removeOldest(bvs);

  while ((ip < bvs->numEntries) && (id > bvs->entries[ip].id))
    ip++;

  sp = bvs->numEntries;
  while (sp > ip) {
    bvs->entries[sp] = bvs->entries[sp - 1];
    sp--;
  }

  bvs->entries[ip].id = id;
  bvs->entries[ip].age = age;
  bvs->numEntries++;
}


/*
 * BVS PUBLIC METHOD IMPLEMENTATIONS
 */

void
BVS_sysInit(void)
{
  BVS_memstack_G = (BubVertexEntry **)safe_calloc(sizeof(BubVertexEntry *), BVS_MEMSTACK_SIZE);
  for (BVS_memstack_top_G = 0; BVS_memstack_top_G < BVS_MEMSTACK_SIZE / 2; ++BVS_memstack_top_G) {
    BVS_memstack_G[BVS_memstack_top_G] = (BubVertexEntry *)safe_malloc(sizeof(BubVertexEntry) * AS_CGB_BUBBLE_set_size_G );
  }
}


void
BVS_sysDone(void)
{
  int t;

  for (t = 0; t < BVS_memstack_top_G; ++t) {
    safe_free(BVS_memstack_G[t]);
  }
  safe_free(BVS_memstack_G);
}


void
BVS_initialize(BubVertexSet_t bvs)
{
  bvs->numEntries = 0;
  bvs->entries = NULL;
}


void
BVS_destroy(BubVertexSet_t bvs)
{
  if (bvs->entries) {
    safe_free(bvs->entries);
  }
}


uint64
BVS_hash(BubVertexSet_t bvs)
{
  uint64 result = 0;
  int i;

  for (i = bvs->numEntries -1 ; i >= 0; --i)
    result = (result << 2) + bvs->entries[i].id;

  return result;
}


int
BVS_compare(BubVertexSet_t bvs1, BubVertexSet_t bvs2)
{
  int p1 = 0, p2 = 0;

  while ((p1 < bvs1->numEntries) && (p2 < bvs2->numEntries)) {
    if (bvs1->entries[p1].id < bvs2->entries[p2].id) {
      return -1;
    }
    else if (bvs1->entries[p1].id == bvs2->entries[p2].id) {
      p1++; p2++;
    }
    else /* bvs1->entries[p1].id > bvs2->entries[p2].id */
      return 1;
  }


  if (bvs1->numEntries < bvs2->numEntries)
    return -1;
  else if (bvs1->numEntries == bvs2->numEntries)
    return 0;
  else /* (bvs1->numEntries > bvs2->numEntries) */
    return 1;
}


int
BVS_numEntries(BubVertexSet_t bvs)
{
  return bvs->numEntries;
}


int
BVS_empty(BubVertexSet_t bvs)
{
  return (bvs->numEntries == 0);
}


int
BVS_full(BubVertexSet_t bvs)
{
  return (bvs->numEntries == AS_CGB_BUBBLE_set_size_G);
}


void
BVS_intersect(BubVertexSet_t bvs1, BubVertexSet_t bvs2)
{
  int p1_rd = 0, p2 = 0, p1_wr = 0;

  if (BVS_empty(bvs1))
    return;

  while ((p1_rd < bvs1->numEntries) && (p2 < bvs2->numEntries)) {
    if (bvs1->entries[p1_rd].id == bvs2->entries[p2].id) {
      bvs1->entries[p1_wr].id = bvs1->entries[p1_rd].id;
      bvs1->entries[p1_wr].age =
	MAX(bvs1->entries[p1_rd].age, bvs2->entries[p2].age);
      p1_rd++; p1_wr++; p2++;
    }
    else if (bvs1->entries[p1_rd].id < bvs2->entries[p2].id)
      p1_rd++;
    else /* bvs1->entries[p1].id > bvs2->entries[p2].id */
      p2++;
  }

  bvs1->numEntries = p1_wr;
  if (bvs1->numEntries == 0) {
    BVS__free(bvs1->entries);
    bvs1->entries = NULL;
  }
}


void
BVS_copy(BubVertexSet_t src_bvs, BubVertexSet_t dst_bvs)
{
  dst_bvs->numEntries = src_bvs->numEntries;
  if (src_bvs->numEntries > 0) {
    if (!dst_bvs->entries)
      dst_bvs->entries = BVS__alloc();
    memcpy(dst_bvs->entries, src_bvs->entries,
	   src_bvs->numEntries * sizeof(BubVertexEntry));
  }
  else
    dst_bvs->entries = NULL;
}


void BVS_age(BubVertexSet_t bvs, int inc)
{
  int i = 0;

  while (i < bvs->numEntries) {
    bvs->entries[i].age += inc;
    if ((AS_CGB_BUBBLE_max_age_G > 0) &&
	(bvs->entries[i].age > AS_CGB_BUBBLE_max_age_G))
      BVS__removeIndex(bvs, i);
    else
      i++;
  }

  if (bvs->numEntries == 0) {
    BVS__free(bvs->entries);
    bvs->entries = NULL;
  }
}


void
BVS_insert(BubVertexSet_t bvs, IntFragment_ID id)
{
  BVS__insertWithAge(bvs, id, 0);
}


int
BVS_subset(BubVertexSet_t bvs1, BubVertexSet_t bvs2)
{
  int p1 = 0, p2 = 0;

  while ((p1 < bvs1->numEntries) && (p2 < bvs2->numEntries)) {
    if (bvs1->entries[p1].id == bvs2->entries[p2].id) {
      p1++; p2++;
    }
    else if (bvs1->entries[p1].id < bvs2->entries[p2].id)
      p1++;
    else /* bvs1->entries[p1].id > bvs2->entries[p2].id */
      return FALSE;
  }

  return (p2 == bvs2->numEntries);
}


void
BVS_makeEmpty(BubVertexSet_t bvs)
{
  bvs->numEntries = 0;
  if (bvs->entries) {
    BVS__free(bvs->entries);
    bvs->entries = NULL;
  }
}


void
BVS_print(BubVertexSet_t bvs, FILE *of)
{
  int i;

  for (i = 0; i < bvs->numEntries; ++i)
    fprintf(of, F_IID " (%d)  ", bvs->entries[i].id, bvs->entries[i].age);
}

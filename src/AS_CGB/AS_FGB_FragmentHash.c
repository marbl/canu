
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
 * $Id: AS_FGB_FragmentHash.c,v 1.3 2005-03-22 19:02:37 jason_miller Exp $
 *
 * Module:
 * Description:
 * Assumptions: Too many to count.
 * Author: Clark Mobarry
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "AS_FGB_FragmentHash.h"

#define  SAFE_MALLOC_QUIET(the_name, the_type, length) \
  assert(NULL == the_name); \
  the_name = (the_type *) malloc(length*sizeof(the_type)); \
  assert(NULL != the_name);
#define  SAFE_MALLOC(the_name, the_type, length) \
  fprintf(stderr, "SAFE_MALLOC " F_SIZE_T " " #the_name " " #the_type  " " F_SIZE_T " " F_SIZE_T "\n", length * sizeof(the_type), length, sizeof(the_type) ); \
  SAFE_MALLOC_QUIET(the_name, the_type, length)

#define  SAFE_CALLOC(the_name, the_type, length) \
  fprintf(stderr, "SAFE_CALLOC " F_SIZE_T " " #the_name " " #the_type  " " F_SIZE_T " " F_SIZE_T "\n", length * sizeof(the_type), length, sizeof(the_type) ); \
  assert(NULL == the_name); \
  the_name = (the_type *) calloc(length,sizeof(the_type)); \
  assert(NULL != the_name);

#define  SAFE_REALLOC(the_name, the_type, length) \
  fprintf(stderr, "SAFE_REALLOC " F_SIZE_T " " #the_name " " #the_type  " " F_SIZE_T " " F_SIZE_T "\n", length * sizeof(the_type), length, sizeof(the_type) ); \
  assert(NULL != the_name); \
  the_name = (the_type *) realloc(the_name,length * sizeof(the_type)); \
  assert(NULL != the_name);

#define SAFE_FREE_QUIET(the_name) \
  assert(NULL != the_name); free(the_name); the_name = NULL;
#define SAFE_FREE(the_name) \
  fprintf(stderr, "SAFE_FREE " #the_name  "\n" ); \
  SAFE_FREE_QUIET(the_name)

#if 0
#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)
#define FALSE 0
#undef DEBUGGING1
#undef DEBUGGING2
#endif

typedef struct {
  IntFragment_ID * iid_to_vid;
  IntFragment_ID max_iid;
  IntFragment_ID min_iid;
  IntFragment_ID allocated;
} FragmentHashStruct;
// Hide the implementation of the object.  I encapsulate the object
// (1) data with a void pointer and (2) methods by extern functions.

#if 0
static void set_max_iid_FragmentHash(FragmentHashObject * _self, IntFragment_ID max_iid)
{
  FragmentHashStruct * self = (FragmentHashStruct *) _self ;
  assert(NULL != self);
  self->max_iid = max_iid;
}

static void set_min_iid_FragmentHash(FragmentHashObject * _self, IntFragment_ID min_iid)
{
  FragmentHashStruct * self = (FragmentHashStruct *) _self ;
  assert(NULL != self);
  self->min_iid = min_iid;
}
#endif

void set_vid_FragmentHash
( FragmentHashObject * _self, 
  IntFragment_ID iid,
  IntFragment_ID vid)
{ 
  FragmentHashStruct * self = (FragmentHashStruct *) _self ;
  //if(AS_CGB_NOT_SEEN_YET != vid)
  assert(NULL != self);
  assert(NULL != self->iid_to_vid);
  self->iid_to_vid[iid] = vid;
}

IntFragment_ID get_vid_FragmentHash
( FragmentHashObject * _self, IntFragment_ID iid)
{
  FragmentHashStruct * self = (FragmentHashStruct *) _self ;
  assert(NULL != self);
  assert(NULL != self->iid_to_vid);
  return self->iid_to_vid[iid];
}

static void clear_FragmentHash(FragmentHashObject * _self)
{
  FragmentHashStruct * self = (FragmentHashStruct *) _self ;
  IntFragment_ID iid;
  assert(NULL != self);
  for(iid=0; iid < self->allocated; iid++) {
    set_vid_FragmentHash(self,iid,AS_CGB_NOT_SEEN_YET);
  }
}


FragmentHashObject * create_FragmentHash(IntFragment_ID max_iid) {
  // Create an object. The return value is non-NULL when successful.
  // The maxiid is a pre-allocation amount for the maximum number of
  // fragments.  You can specify zero for the pre-allocation amounts.
  FragmentHashStruct * self = 
    (FragmentHashStruct *)malloc(sizeof(FragmentHashStruct));
  assert(NULL != self);
  memset(self,0,sizeof(self));
  assert(max_iid < AS_CGB_NOT_SEEN_YET);
  self->allocated = max_iid+1;
  SAFE_MALLOC(self->iid_to_vid, IntFragment_ID, self->allocated);
  clear_FragmentHash(self);
  return (FragmentHashObject *)self;
}

int destroy_FragmentHash(FragmentHashObject * _self)
{
  // Destroy an object. The return value is zero when successful.
  FragmentHashStruct * self = (FragmentHashStruct *) _self ;
  assert(NULL != self->iid_to_vid); 
  free( self->iid_to_vid ); 
  self->iid_to_vid = NULL;
  SAFE_FREE(self);
  return 0;
}

//////////////////////////////////////////////////////////////////////

// Re-hash the fragment IID to fragment VID mapping using the
// fragments already in the store.


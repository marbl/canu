
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
 * $Id: AS_FGB_FragmentHash.c,v 1.7 2007-02-14 07:20:05 brianwalenz Exp $
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

//  Copied from AS_CGB_all.h

typedef struct {
  IntFragment_ID * iid_to_vid;
  IntFragment_ID max_iid;
  IntFragment_ID min_iid;
  IntFragment_ID allocated;
} FragmentHashStruct;
// Hide the implementation of the object.  I encapsulate the object
// (1) data with a void pointer and (2) methods by extern functions.

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
    (FragmentHashStruct *)safe_malloc(sizeof(FragmentHashStruct));
  memset(self,0,sizeof(self));
  assert(max_iid < AS_CGB_NOT_SEEN_YET);
  self->allocated = max_iid+1;
  self->iid_to_vid = safe_malloc(sizeof(IntFragment_ID) * self->allocated);
  clear_FragmentHash(self);
  return (FragmentHashObject *)self;
}

int destroy_FragmentHash(FragmentHashObject * _self)
{
  // Destroy an object. The return value is zero when successful.
  FragmentHashStruct * self = (FragmentHashStruct *) _self ;
  assert(NULL != self->iid_to_vid); 
  safe_free( self->iid_to_vid ); 
  self->iid_to_vid = NULL;
  safe_free(self);
  return 0;
}

//////////////////////////////////////////////////////////////////////

// Re-hash the fragment IID to fragment VID mapping using the
// fragments already in the store.


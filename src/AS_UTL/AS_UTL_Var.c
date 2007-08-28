
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

static char CM_ID[] = "$Id: AS_UTL_Var.c,v 1.19 2007-08-28 22:49:14 brianwalenz Exp $";

/********************************************************************/
/* Variable Length C Array Package 
 * 
 *     Saul A. Kravitz
 *     January 1999
 *
 * This package is meant to simplify the coding and manipulation of
 * variable length, auto-resizing arrays.
 * It defines a basic set of operations, and provides a set of
 * macros that expand to support typesafe manipulation of the
 * arrays.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "AS_global.h"
#include "math_AS.h"
#include "AS_UTL_Var.h"

#include "AS_UTL_fileIO.h"


//  As a debugging aid (and to torture developers) we can force the VA
//  to always invalidate pointers.
//
#undef ALWAYS_MOVE_VA_ON_MAKEROOM


//  We write a VarArrayType as a known-size structure -- The size of a
//  VarArrayType changes between 32-bit and 64-bit platforms.
//
typedef struct {
  uint64 Elements;
  uint64 sizeofElement;
  uint64 numElements;
  uint64 allocatedElements;
  char   typeofElement[VA_TYPENAMELEN];
} FileVarArrayType;



int
MakeRoom_VA(VarArrayType * const va,
            const size_t         maxElements,
            const int            pad_to_a_power_of_two) {

  size_t newElements, newSize, tentativeNewSize, oldSize;
  char *mem = NULL;

#ifdef DEBUG
  fprintf(stderr,"* MakeRoom_VA for handle %p\n",va);
  fprintf(stderr,"* Elements              = %p\n", va->Elements);
  fprintf(stderr,"* numElements           = " F_SIZE_T "\n", va->numElements);
  fprintf(stderr,"* allocatedElements     = " F_SIZE_T "\n", va->allocatedElements);
  fprintf(stderr,"* sizeofElement         = " F_SIZE_T "\n", va->sizeofElement);
  fprintf(stderr,"* typeofElement         = %s\n", va->typeofElement);
  fprintf(stderr,"* requested maxElements = " F_SIZE_T "\n", maxElements);
#endif


#ifndef ALWAYS_MOVE_VA_ON_MAKEROOM
  /* If we have enough space, then return now. */
  if(maxElements <= (va->allocatedElements)){
#ifdef DEBUG
    fprintf(stderr,"* MakeRoom_VA called with "
	    "maxElements " F_SIZE_T " < allocElements " F_SIZE_T "...returning\n",
	    maxElements, va->allocatedElements);
#endif
    return FALSE;
  }
#endif

  oldSize = (va->allocatedElements)*(va->sizeofElement);
  
#ifdef DEBUG
  fprintf(stderr,"* MakeRoom_VA oldSize=" F_SIZE_T "\n",oldSize);
#endif

  // Minimimum allocation is one element;
  newSize = MAX(maxElements, 1)*(va->sizeofElement);
    
  if(pad_to_a_power_of_two) {
    /* Compute a power-of-two allocation size for the va */
    
    // Only allocate a power of 2 number of bytes.
    tentativeNewSize = (((size_t)1) << ceil_log2(newSize));

    // Cap alloc'd size at 512MB to decrease failure rate
    if ( tentativeNewSize - newSize > (2 << 28) )
        tentativeNewSize = oldSize + (2 << 28);
    
    // If we need to use the end of the block, do it
    newSize = MAX(newSize, tentativeNewSize);
  }

#ifdef DEBUG
  fprintf(stderr,"* MakeRoom_VA newSize=" F_SIZE_T "\n",newSize);
#endif

#ifdef ALWAYS_MOVE_VA_ON_MAKEROOM
  if (newSize < oldSize)
    newSize = oldSize;
#endif
  
  assert(oldSize <= newSize);
  newElements = (newSize)/va->sizeofElement;
  assert( va->allocatedElements <= newElements);
  assert( maxElements <= newElements);
  
#ifndef ALWAYS_MOVE_VA_ON_MAKEROOM
  //  Do not need to do anything.
  if (newSize <= oldSize)
    return FALSE;
#endif

  if( NULL == va->Elements ) {
    mem = (char *)safe_calloc(newSize, sizeof(char));
  } else {
#ifdef ALWAYS_MOVE_VA_ON_MAKEROOM
    mem = (char *)safe_calloc(newSize, sizeof(char));
    memcpy(mem, va->Elements, oldSize);
    memset(va->Elements, 0xff, oldSize);
    safe_free(va->Elements);
#else
    mem = (char *)safe_realloc(va->Elements, newSize);
    memset(mem + oldSize, 0, newSize - oldSize);

    if (0) {
      int i=0, j=0;
      for (i=0; i<newSize; i++)
        j += mem[i];
      fprintf(stderr, "realloc'd VA (i=%d) with sum %d\n", i, j);
    }
#endif
  }

  va->Elements = mem;
  va->allocatedElements = newElements;

#ifdef ALWAYS_MOVE_VA_ON_MAKEROOM
  if (oldSize > 0)
    fprintf(stderr, "* MakeRoom_VA reallocated '%s' from "F_SIZE_T" bytes to "F_SIZE_T" bytes.\n", va->typeofElement, oldSize, newSize);
#endif

  return TRUE;
}


VarArrayType *
Create_VA(const size_t numElements,
          const size_t sizeofElement,
          const char * const thetype) {
  VarArrayType * const va = (VarArrayType *)safe_calloc(1, sizeof(VarArrayType));

  va->Elements          = NULL;
  va->sizeofElement     = sizeofElement;
  va->numElements       = 0;
  va->allocatedElements = 0;

  strncpy(va->typeofElement,thetype,VA_TYPENAMELEN);
  va->typeofElement[VA_TYPENAMELEN-1] = (char)0;

  MakeRoom_VA(va, numElements, FALSE);

  return va;
}


void
Clear_VA(VarArrayType * const va){
  if(NULL == va)
    return;
  safe_free(va->Elements);
  memset(va, 0, sizeof(VarArrayType));
}


//  Append vb's data onto va
void
Concat_VA(VarArrayType * const va,
          const VarArrayType * const vb){

  assert(NULL != va);
  assert(NULL != vb);
  assert(va != vb);

  size_t asize = va->numElements * va->sizeofElement;
  size_t bsize = vb->numElements * vb->sizeofElement;

  if ((asize + bsize == 0) || (bsize == 0))
    return;

  MakeRoom_VA(va, va->numElements + vb->numElements, FALSE);
  va->numElements += vb->numElements;
  memcpy(va->Elements + asize, vb->Elements, bsize);
}


void
ResetToRange_VA(VarArrayType * const va, const size_t indx){

  //  Resetting to a larger array is equivalent to EnableRange
  //
  if (indx > va->numElements) {
    EnableRange_VA(va,indx);
    return;
  }

  if (indx == va->numElements)
    return;

  // Resetting to a smaller array, zeros out the unused elements and
  // resets numElements.

  memset(va->Elements + (va->sizeofElement * indx), 0, va->sizeofElement * (va->numElements - indx));

  va->numElements = indx;
}


void
EnableRange_VA(VarArrayType * const va, const size_t maxElements){

  if (maxElements > va->allocatedElements)
    MakeRoom_VA(va,maxElements,TRUE);

  assert(maxElements == 0 || va->Elements != NULL);

  if(maxElements >= va->numElements)
    va->numElements = maxElements;
}


void
SetElement_VA(VarArrayType * const va,
              const size_t indx, 
              const void * const data){
  EnableRange_VA(va, (indx+1));
  memcpy(va->Elements + indx * va->sizeofElement, (char *)data, va->sizeofElement);
}


void
SetRange_VA(VarArrayType * const va,
            const size_t indx, 
            const size_t numElements, 
            const void * const data){
  EnableRange_VA(va, indx + numElements);
  memcpy(va->Elements + indx * va->sizeofElement, (char *)data, va->sizeofElement * numElements);
}


VarArrayType *
Clone_VA(const VarArrayType *fr){
  VarArrayType *to = (VarArrayType *)safe_calloc(1, sizeof(VarArrayType));

  to->Elements          = NULL;
  to->sizeofElement     = fr->sizeofElement;
  to->numElements       = 0;
  to->allocatedElements = 0;

  strncpy(to->typeofElement, fr->typeofElement, VA_TYPENAMELEN);

  MakeRoom_VA(to, fr->numElements, FALSE);
  EnableRange_VA(to, fr->numElements);

  if (fr->numElements > 0)
    memcpy(to->Elements, fr->Elements, fr->sizeofElement * fr->numElements);

  return(to);
}

void
ReuseClone_VA(VarArrayType *to, const VarArrayType *fr){

  if ((fr->sizeofElement != to->sizeofElement) ||
      (strcmp(fr->typeofElement, to->typeofElement) != 0)) {

    safe_free(to->Elements);

    to->sizeofElement      = fr->sizeofElement;
    to->numElements        = 0;
    to->allocatedElements  = 0;

    strncpy(to->typeofElement, fr->typeofElement, VA_TYPENAMELEN);
  } else {
    ResetToRange_VA(to, 0);
  }

  MakeRoom_VA(to, fr->numElements, FALSE);
  EnableRange_VA(to, fr->numElements);

  if (fr->numElements > 0)
    memcpy(to->Elements, fr->Elements, fr->sizeofElement * fr->numElements);
}

void
LoadFromFile_VA(FILE *fp,
                VarArrayType *va) {

  FileVarArrayType    vat = {0};

  AS_UTL_safeRead(fp, &vat, "CreateFromFile_VA (vat)", sizeof(FileVarArrayType), 1);

  assert(vat.numElements <= vat.allocatedElements);

  if(strncmp(va->typeofElement, vat.typeofElement, VA_TYPENAMELEN)){
    fprintf(stderr,"* Expecting array of type <%s> but read array of type <%s>\n",
            va->typeofElement, vat.typeofElement);
    assert(0);
  }

  MakeRoom_VA(va, vat.numElements, FALSE);
  EnableRange_VA(va, vat.numElements);

  assert((vat.numElements == 0) || (va->Elements != NULL));
  assert(vat.numElements == va->numElements);
  assert(vat.sizeofElement == va->sizeofElement);
  assert(vat.sizeofElement > 0);

  AS_UTL_safeRead(fp, va->Elements, "CreateFromFile_VA", va->sizeofElement, va->numElements);
}


VarArrayType *
CreateFromFile_VA(FILE *fp,
                  const char * const thetype,
                  size_t additional_elements) {

  FileVarArrayType    vat = {0};
  VarArrayType       *va  = (VarArrayType *)safe_calloc(1, sizeof(VarArrayType));

  AS_UTL_safeRead(fp, &vat, "CreateFromFile_VA (vat)", sizeof(FileVarArrayType), 1);

  assert(vat.numElements <= vat.allocatedElements);

  if(strncmp(vat.typeofElement,thetype,VA_TYPENAMELEN)){
    fprintf(stderr,"* Expecting array of type <%s> but read array of type <%s>\n",
	    thetype, vat.typeofElement);
    assert(0);
  }

  // We construct a VA big enough to hold numElements.  This will
  // 'compress' the VA, if it was initially allocated much larger than
  // necessary.

  va->Elements          = NULL;
  va->sizeofElement     = vat.sizeofElement;
  va->numElements       = 0;
  va->allocatedElements = 0;

  strncpy(va->typeofElement, vat.typeofElement, VA_TYPENAMELEN);
  va->typeofElement[VA_TYPENAMELEN-1] = 0;

  MakeRoom_VA(va, vat.numElements + additional_elements, FALSE);

  EnableRange_VA(va, vat.numElements);

  assert(vat.numElements == 0 || (va->Elements != NULL));
  assert(vat.numElements == va->numElements);
  assert(vat.sizeofElement == va->sizeofElement);
  assert(vat.sizeofElement > 0);

  AS_UTL_safeRead(fp, va->Elements, "CreateFromFile_VA", va->sizeofElement, va->numElements);

  return(va);
}



size_t CopyToFile_VA(const VarArrayType * const va,FILE *fp){
  FileVarArrayType vat = {0};

  assert(fp != NULL);
  assert(va != NULL);
  assert(va->numElements == 0 || va->Elements != NULL);
  assert(va->sizeofElement > 0);

  vat.Elements           = 0;
  vat.sizeofElement      = va->sizeofElement;
  vat.numElements        = va->numElements;
  vat.allocatedElements  = va->allocatedElements;

  strncpy(vat.typeofElement, va->typeofElement, VA_TYPENAMELEN);

  AS_UTL_safeWrite(fp, &vat,         "CopyToFile_VA (vat)", sizeof(FileVarArrayType), 1);
  AS_UTL_safeWrite(fp, va->Elements, "CopyToFile_VA (dat)", va->sizeofElement,        va->numElements);

  return(sizeof(FileVarArrayType) + va->sizeofElement * va->numElements);
}

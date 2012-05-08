
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

static char *rcsid = "$Id: AS_UTL_Var.c,v 1.35 2012-05-08 23:17:55 brianwalenz Exp $";

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

#include "AS_UTL_Var.h"
#include "AS_UTL_fileIO.h"

//  As a debugging aid (and to torture developers) we can force the VA
//  to always invalidate pointers.
//
#undef ALWAYS_MOVE_VA_ON_MAKEROOM

//  Trash the memory used by a deleted VA.  Very, very slow.
//
#undef TRASH_DELETED_VA


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
MakeRoom_VA(VarArrayType *va,
            size_t         maxElements) {

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
Create_VA(size_t numElements,
          size_t sizeofElement,
          const char *thetype) {
  VarArrayType *va = (VarArrayType *)safe_calloc(1, sizeof(VarArrayType));

  va->Elements          = NULL;
  va->sizeofElement     = sizeofElement;
  va->numElements       = 0;
  va->allocatedElements = 0;

  strncpy(va->typeofElement,thetype,VA_TYPENAMELEN);
  va->typeofElement[VA_TYPENAMELEN-1] = (char)0;

  MakeRoom_VA(va, numElements);

  return va;
}


void
Clear_VA(VarArrayType *va){
  if (NULL == va)
    return;
  safe_free(va->Elements);
  memset(va, 0, sizeof(VarArrayType));
}


void
Trash_VA(VarArrayType *va){
  if (NULL == va)
    return;
#ifdef TRASH_DELETED_VA
  if (va->allocatedElements * va->sizeofElement < 1024) {
    int i;
    for (i=0; i<va->allocatedElements * va->sizeofElement; i++)
      va->Elements[i] = 0xff;
  } else {
    memset(va->Elements, 0xff, va->allocatedElements * va->sizeofElement);
  }
#endif
  safe_free(va->Elements);
  safe_free(va);
}


//  Append vb's data onto va
void
Concat_VA(VarArrayType *va,
          VarArrayType *vb){

  assert(NULL != va);
  assert(NULL != vb);
  assert(va != vb);

  size_t asize = va->numElements * va->sizeofElement;
  size_t bsize = vb->numElements * vb->sizeofElement;

  if ((asize + bsize == 0) || (bsize == 0))
    return;

  MakeRoom_VA(va, va->numElements + vb->numElements);
  va->numElements += vb->numElements;
  assert(va->Elements + asize != vb->Elements);
  memcpy(va->Elements + asize, vb->Elements, bsize);
}


void
ResetToRange_VA(VarArrayType *va, size_t indx){

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

  assert(va->numElements > 0);
  assert(va->Elements != NULL);

  memset(va->Elements + (va->sizeofElement * indx), 0, va->sizeofElement * (va->numElements - indx));

  va->numElements = indx;
}


void
EnableRange_VA(VarArrayType *va, size_t maxElements){

  if (maxElements > va->allocatedElements)
    MakeRoom_VA(va, maxElements);  //  Was allocating a power of two

  assert(maxElements == 0 || va->Elements != NULL);

  if(maxElements >= va->numElements)
    va->numElements = maxElements;
}


void
SetElements_VA(VarArrayType *va,
               size_t        indx,
               void         *data,
               size_t       nume){
  EnableRange_VA(va, (indx+nume));
  if (va->Elements + indx * va->sizeofElement != data)
    memcpy(va->Elements + indx * va->sizeofElement, data, va->sizeofElement * nume);
}


VarArrayType *
Clone_VA(VarArrayType *fr){
  VarArrayType *to = (VarArrayType *)safe_calloc(1, sizeof(VarArrayType));

  to->Elements          = NULL;
  to->sizeofElement     = fr->sizeofElement;
  to->numElements       = 0;
  to->allocatedElements = 0;

  strncpy(to->typeofElement, fr->typeofElement, VA_TYPENAMELEN);

  MakeRoom_VA(to, fr->numElements);
  EnableRange_VA(to, fr->numElements);

  if (fr->numElements > 0)
    memcpy(to->Elements, fr->Elements, fr->sizeofElement * fr->numElements);

  return(to);
}

void
ReuseClone_VA(VarArrayType *to, VarArrayType *fr){

  assert(to != fr);

  if ((fr->sizeofElement != to->sizeofElement) ||
      (strcmp(fr->typeofElement, to->typeofElement) != 0)) {

    safe_free(to->Elements);

    to->Elements           = NULL;
    to->sizeofElement      = fr->sizeofElement;
    to->numElements        = 0;
    to->allocatedElements  = 0;

    strncpy(to->typeofElement, fr->typeofElement, VA_TYPENAMELEN);
  } else {
    ResetToRange_VA(to, 0);
  }

  MakeRoom_VA(to, fr->numElements);
  EnableRange_VA(to, fr->numElements);

  if (fr->numElements > 0)
    memcpy(to->Elements, fr->Elements, fr->sizeofElement * fr->numElements);
}

static
void
ReadVA(FILE *fp, VarArrayType *va, FileVarArrayType *vat) {
  MakeRoom_VA(va, vat->numElements);
  EnableRange_VA(va, vat->numElements);

  if (vat->numElements > 0) {
    assert(va->Elements != NULL);
    assert(vat->sizeofElement == va->sizeofElement);

    size_t numRead = AS_UTL_safeRead(fp, va->Elements, "LoadFromFile_VA", va->sizeofElement, va->numElements);

    if (va->numElements != numRead)
      fprintf(stderr, "ReadVA()-- Short read from va <%s>; expected "F_SIZE_T" elements, read "F_SIZE_T" elements.\n",
              va->typeofElement, va->numElements, numRead), exit(1);
  }
}

void
LoadFromFile_VA(FILE *fp,
                VarArrayType *va) {

  FileVarArrayType    vat = {0, 0, 0, 0, {0}};

  if (1 != AS_UTL_safeRead(fp, &vat, "LoadFromFile_VA (vat)", sizeof(FileVarArrayType), 1))
    fprintf(stderr, "LoadFromFile_VA()-- Failed to read vat\n"), exit(1);

  assert(vat.numElements <= vat.allocatedElements);

  if(strncmp(va->typeofElement, vat.typeofElement, VA_TYPENAMELEN))
    fprintf(stderr,"* Expecting array of type <%s> but read array of type <%s>\n",
            va->typeofElement, vat.typeofElement), exit(1);

  ReadVA(fp, va, &vat);
}


VarArrayType *
CreateFromFile_VA(FILE *fp,
                  const char *thetype) {

  FileVarArrayType    vat = {0, 0, 0, 0, {0}};
  VarArrayType       *va  = (VarArrayType *)safe_calloc(1, sizeof(VarArrayType));

  if (1 != AS_UTL_safeRead(fp, &vat, "CreateFromFile_VA (vat)", sizeof(FileVarArrayType), 1))
    fprintf(stderr, "LoadFromFile_VA()-- Failed to read vat\n"), exit(1);

  assert(vat.numElements <= vat.allocatedElements);

  if(strncmp(vat.typeofElement,thetype,VA_TYPENAMELEN))
    fprintf(stderr,"* Expecting array of type <%s> but read array of type <%s>\n",
	    thetype, vat.typeofElement), exit(1);

  // We construct a VA just big enough to hold all the elements on disk.

  va->Elements          = NULL;
  va->sizeofElement     = vat.sizeofElement;
  va->numElements       = 0;
  va->allocatedElements = 0;

  strncpy(va->typeofElement, vat.typeofElement, VA_TYPENAMELEN);
  va->typeofElement[VA_TYPENAMELEN-1] = 0;

  ReadVA(fp, va, &vat);

  return(va);
}



size_t CopyToFile_VA(VarArrayType *va,FILE *fp){
  FileVarArrayType vat = {0, 0, 0, 0, {0}};

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













void
LoadFromMemory_VA(char *&memory,
                  VarArrayType *va) {

  assert(memory != NULL);

  FileVarArrayType    vat = {0, 0, 0, 0, {0}};

  memcpy(&vat, memory, sizeof(FileVarArrayType));
  memory += sizeof(FileVarArrayType);

  assert(vat.numElements <= vat.allocatedElements);

  if(strncmp(va->typeofElement, vat.typeofElement, VA_TYPENAMELEN))
    fprintf(stderr,"* Expecting array of type <%s> but read array of type <%s>\n",
            va->typeofElement, vat.typeofElement), exit(1);

  MakeRoom_VA(va, vat.numElements);
  EnableRange_VA(va, vat.numElements);

  if (vat.numElements > 0) {
    assert(va->Elements != NULL);
    assert(vat.sizeofElement == va->sizeofElement);

    memcpy(va->Elements, memory, va->sizeofElement * va->numElements);
    memory += va->sizeofElement * va->numElements;
  }
}


VarArrayType *
CreateFromMemory_VA(char *&memory,
                    const char  *thetype) {

  assert(memory != NULL);

  FileVarArrayType    vat = {0, 0, 0, 0, {0}};
  VarArrayType       *va  = (VarArrayType *)safe_calloc(1, sizeof(VarArrayType));

  memcpy(&vat, memory, sizeof(FileVarArrayType));
  memory += sizeof(FileVarArrayType);

  assert(vat.numElements <= vat.allocatedElements);

  if(strncmp(vat.typeofElement,thetype,VA_TYPENAMELEN))
    fprintf(stderr,"* Expecting array of type <%s> but read array of type <%s>\n",
	    thetype, vat.typeofElement), exit(1);

  // We construct a VA just big enough to hold all the elements on disk.

  va->Elements          = NULL;
  va->sizeofElement     = vat.sizeofElement;
  va->numElements       = 0;
  va->allocatedElements = 0;

  strncpy(va->typeofElement, vat.typeofElement, VA_TYPENAMELEN);
  va->typeofElement[VA_TYPENAMELEN-1] = 0;

  MakeRoom_VA(va, vat.numElements);
  EnableRange_VA(va, vat.numElements);

  if (vat.numElements > 0) {
    assert(va->Elements != NULL);
    assert(vat.sizeofElement == va->sizeofElement);

    memcpy(va->Elements, memory, va->sizeofElement * va->numElements);
    memory += va->sizeofElement * va->numElements;
  }

  return(va);
}


size_t
CopyToMemory_VA(VarArrayType *va,
                char         *&memory) {

  if (memory == NULL)
    return(sizeof(FileVarArrayType) + va->sizeofElement * va->numElements);

  assert(va != NULL);
  assert(va->numElements == 0 || va->Elements != NULL);
  assert(va->sizeofElement > 0);

  FileVarArrayType vat = {0, 0, 0, 0, {0}};

  vat.Elements           = 0;
  vat.sizeofElement      = va->sizeofElement;
  vat.numElements        = va->numElements;
  vat.allocatedElements  = va->allocatedElements;

  strncpy(vat.typeofElement, va->typeofElement, VA_TYPENAMELEN);

  memcpy(memory, &vat, sizeof(FileVarArrayType));
  memory += sizeof(FileVarArrayType);

  memcpy(memory,  va->Elements, va->numElements * va->sizeofElement);
  memory += va->numElements * va->sizeofElement;

  return(sizeof(FileVarArrayType) + va->sizeofElement * va->numElements);
}

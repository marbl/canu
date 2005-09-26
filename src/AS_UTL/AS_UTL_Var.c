
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
static char CM_ID[] = "$Id: AS_UTL_Var.c,v 1.9 2005-09-26 20:00:02 brianwalenz Exp $";
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
 * Revision: $Revision: 1.9 $
 * Date:     $Date: 2005-09-26 20:00:02 $
 * CMM, 1999/03/29:  Ported to large arrays on the Digital systems by declaring
 * array sizes using size_t, rather than unit32.
 *
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


//  We write a VarArrayType as a known-size structure -- The size of a
//  VarArrayType changes between 32-bit and 64-bit platforms.
//
typedef struct {
  uint64 Elements;
  uint64 sizeofElement;
  uint64 numElements;
  uint64 allocatedElements;
  char typeofElement[VA_TYPENAMELEN];
} FileVarArrayType;


  
// Re-allocations are rounded up to a power of two number of bytes
// when the EnableRange functionality is used.  This can be avoided by
// using Initialize_VA with the number of elements to pre-allocate.

#define DEBUG
#undef DEBUG

// *******************************************************************
//
// The routines that directly allocate or free array memory are
// MakeRoom_VA and Clear_VA.
//

int MakeRoom_VA
(VarArrayType * const va,
 const size_t         maxElements,
 const int            pad_to_a_power_of_two
)
{
  /* CMM 1999/03/29
     The return value is a flag that indicates whether
     a reallocation occurred. 

     PREVIOUSLY, maxElements was the index of the maximum
     element.  Now, it is the maximum number of elements.
  */

  size_t newElements, newSize, tentativeNewSize, oldSize;
  char *mem = NULL;

#ifdef DEBUG
  fprintf(stderr,"* MakeRoom_VA for handle %p\n",va);
  fprintf(stderr,"* Elements = %p\n", va->Elements);
  fprintf(stderr,"* numElements = " F_SIZE_T "\n", va->numElements);
  fprintf(stderr,"* allocatedElements = " F_SIZE_T "\n",
          va->allocatedElements);
  fprintf(stderr,"* sizeofElement = " F_SIZE_T "\n", va->sizeofElement);
  fprintf(stderr,"* typeofElement = %s\n", va->typeofElement);
  fprintf(stderr,"* requested maxElements = " F_SIZE_T "\n", maxElements);
#endif

  /* If we have enough space, then return now. */
  if(maxElements <= (va->allocatedElements)){
#ifdef DEBUG
    fprintf(stderr,"* MakeRoom_VA called with "
	    "maxElements " F_SIZE_T " < allocElements " F_SIZE_T "...returning\n",
	    maxElements, va->allocatedElements);
#endif
    // assert(va->Elements != NULL);
    return FALSE;
  }

  oldSize = (va->allocatedElements)*(va->sizeofElement);
  
#ifdef DEBUG
  fprintf(stderr,"* MakeRoom_VA oldSize=" F_SIZE_T "\n",oldSize);
#endif

  // Minimimum allocation is one element;
  newSize = max(maxElements, 1)*(va->sizeofElement);
    
  if(pad_to_a_power_of_two) {
    /* Compute a power-of-two allocation size for the va */
    
    // Only allocate a power of 2 number of bytes.
    tentativeNewSize = (((size_t)1) << ceil_log2(newSize));

    // Cap malloc size at 512MB to decrease failure rate
    if ( tentativeNewSize - newSize > (2 << 28) )
        tentativeNewSize = oldSize + (2 << 28);
    
    // If we need to use the end of the block, do it
    newSize = max(newSize, tentativeNewSize);
  }
  
  assert(oldSize <= newSize);
  newElements = (newSize)/va->sizeofElement;
  assert( va->allocatedElements <= newElements);
  assert( maxElements <= newElements);
  
  if(newSize <= oldSize ) return FALSE; /* Do not need to do anything. */
  assert(newSize > oldSize);

  if( NULL == va->Elements ) {
    mem = (char *)malloc(newSize);
  } else {
    mem = (char *)realloc(va->Elements, newSize);
  }

  assert(mem != NULL);
  va->Elements = mem;
  va->allocatedElements = newElements;

  // Initialize the new memory space to zero.
  mem = va->Elements + oldSize;
  memset(mem,0,newSize - oldSize);

  assert(va->Elements != NULL);
  return TRUE;
}

void Clear_VA(VarArrayType * const va){
  /* Handle NULL va's gracefully */
  if(NULL == va)
    return;
  if( NULL != va->Elements ) {
    free(va->Elements);
  }
  va->allocatedElements = 0;
  va->numElements = 0;
  va->sizeofElement = 0;
  va->Elements = NULL;
  va->typeofElement[0] = 0; // A zero length string.
}


// *******************************************************************
//
// Utility functions
//

void Concat_VA(VarArrayType * const va,const VarArrayType * const vb){
  char *aoffset;
  size_t asize, bsize;

  assert( NULL != va );
  assert( NULL != vb );
  assert( va != vb );
  
  /* Compute the size of a */
  asize = va->numElements * (size_t) va->sizeofElement;
  /* Compute the size of b */
  bsize = vb->numElements * (size_t) vb->sizeofElement;

  /* Make sure we have enough space */
  MakeRoom_VA(va, va->numElements + vb->numElements, FALSE);

  /* Update the number of elements in va */
  va->numElements += vb->numElements;

  if( asize + bsize > 0) {
    /* We also assume that the elements of va and vb are not aliased.
       Since memcpy would fail.*/
    assert(va->Elements != vb->Elements);
    
    /* Append vb's data onto va */
    if( bsize > 0 ) {
      memcpy(va->Elements + asize, vb->Elements, bsize);
    }
  }

  return;

}

void ResetToRange_VA(VarArrayType * const va, const size_t indx){
  size_t oldSize;
  char *oldSuffix;
  assert(va != NULL);

  /* Resetting to a larger array is equivalent to EnableRange */
  if(indx > va->numElements){
    EnableRange_VA(va,indx);
    return;
  }

  if(indx == va->numElements) {
    // va->sizeofElement could be NULL if va->numElements == 0.
    return;
  }

  assert( indx < va->numElements );
  assert( NULL != va->Elements );
  
  /* Resetting to a smaller array, zeros out the old suffix and
     resets numElements.
  */

  oldSize = (size_t) va->sizeofElement * (va->numElements - indx);
  oldSuffix = va->Elements + ((size_t)va->sizeofElement * indx);

  /* Zero out previously used elements of va */
  memset(oldSuffix,0,oldSize);

  va->numElements = indx;

  return;
}

void EnableRange_VA(VarArrayType * const va, const size_t maxElements){
  assert(va != NULL);

#ifdef DEBUG
  fprintf(stderr,"* EnableRange " F_SIZE_T " before (" F_SIZE_T ") %p\n",
	  maxElements, va->numElements, va->Elements);
#endif
  if(maxElements >  va->allocatedElements){  /* Allocate more space */
    MakeRoom_VA(va,maxElements,TRUE);
  } 
#ifdef DEBUG 
    fprintf(stderr,"* EnableRange " F_SIZE_T " after (" F_SIZE_T ") %p\n",
	    maxElements, va->allocatedElements, va->Elements);
#endif
  assert(maxElements == 0 || va->Elements != NULL);

  if(maxElements >= va->numElements) { va->numElements = maxElements;}
}

void SetElement_VA(VarArrayType * const va,
		   const size_t indx, 
		   const void * const data){
  size_t indexOffset;

  assert(va != NULL);
  EnableRange_VA(va, (indx+1));

  /* Copy the data, overwriting what was there before */
  indexOffset = indx * (size_t) va->sizeofElement;
  assert( NULL != va->Elements );
  memcpy(va->Elements + indexOffset, (char *)data, (size_t)va->sizeofElement);

  return ;
}

void SetRange_VA(VarArrayType * const va,
		   const size_t indx, 
		   const size_t num_elements, 
		   const void * const data){
  size_t indexOffset;

  assert(va != NULL);
  EnableRange_VA(va, (indx+num_elements));

  /* Copy the data, overwriting what was there before */
  indexOffset = indx * (size_t) va->sizeofElement;
  assert( NULL != va->Elements );
  memcpy(va->Elements + indexOffset, (char *)data,
         (size_t)va->sizeofElement*num_elements);

  return ;
}

void Initialize_VA
( VarArrayType * const va,
  const size_t arraySize,
  const size_t sizeofElement,
  const char * const thetype)
{
  /* 
     Initialize new variable length array handle.
     
     arraySize: The initial number of elements. This can be zero.
     sizeofElement: The size in bytes of the elements.
     thetype: A character string identifying the type.
  */
     
  assert(va != NULL);
  va->Elements = NULL;
  va->sizeofElement = sizeofElement;
  va->numElements = 0;
  va->allocatedElements = 0;
  strncpy(va->typeofElement,thetype,VA_TYPENAMELEN);
  /* Now make sure that the character string is zero terminated. */
  assert(VA_TYPENAMELEN > 0);
  va->typeofElement[VA_TYPENAMELEN-1] = (char)0;

  MakeRoom_VA(va,arraySize,FALSE);
}

/*************************************************/
void ReInitialize_VA
( VarArrayType * const va,
  const size_t arraySize,
  const size_t sizeofElement,
  const char * const thetype)
{
  /* 
     Initialize new variable length array handle.
     
     arraySize: The initial number of elements. This can be zero.
     sizeofElement: The size in bytes of the elements.
     thetype: A character string identifying the type.
  */
     
  assert(va != NULL);
  if(va->sizeofElement != sizeofElement ||
     strcmp(va->typeofElement, thetype)){
    free(va->Elements);
    va->Elements = NULL;
    va->sizeofElement = sizeofElement;
    va->numElements = 0;
    va->allocatedElements = 0;
    strncpy(va->typeofElement,thetype,VA_TYPENAMELEN);
    /* Now make sure that the character string is zero terminated. */
    assert(VA_TYPENAMELEN > 0);
    va->typeofElement[VA_TYPENAMELEN-1] = (char)0;
  }else{
    ResetToRange_VA(va,0);
  }
  MakeRoom_VA(va,arraySize,FALSE);
}

void InitializeFromArray_VA
( VarArrayType * const va,
  const size_t arraySize,
  const size_t sizeofElement,
  const char * const thetype,
  const void * const data
  )
{
  /* 
     Take a variable length array and initialize
     with a copy of the input array 'data'.
     
     arraySize: The initial number of elements. This can be zero.
     sizeofElement: The size in bytes of the elements.
     thetype: A character string identifying the type.
  */
     
  Initialize_VA( va, arraySize, sizeofElement, thetype);
  if(arraySize > 0){
    EnableRange_VA(va,arraySize);
    assert( NULL != va->Elements );
    memcpy(va->Elements,data,sizeofElement*arraySize);
  }
}

void ReInitializeFromArray_VA
( VarArrayType * const va,
  const size_t arraySize,
  const size_t sizeofElement,
  const char * const thetype,
  const void * const data
  )
{
  /* 
     Take a variable length array and initialize
     with a copy of the input array 'data'.
     
     arraySize: The initial number of elements. This can be zero.
     sizeofElement: The size in bytes of the elements.
     thetype: A character string identifying the type.
  */
     
  ReInitialize_VA( va, arraySize, sizeofElement, thetype);
  if(arraySize > 0){
    EnableRange_VA(va,arraySize);
    assert( NULL != va->Elements );
    memcpy(va->Elements,data,sizeofElement*arraySize);
  }
}

void ReInitializeFromFile_VA
( FILE *fp,
  VarArrayType * const va,
  const char * const thetype,
  const size_t allocate_additional_elements
  ){
  FileVarArrayType vat;
  size_t nitems=0;
  assert(fp != NULL);
  /* we should use the safe read and write routines here. */
  nitems = fread(&vat,sizeof(FileVarArrayType),1,fp);
  assert(1 == nitems);
  assert(vat.numElements <= vat.allocatedElements);
  if(strncmp(vat.typeofElement,thetype,VA_TYPENAMELEN)){
    fprintf(stderr,"* Expecting array of type <%s> read array of type <%s>\n",
	    thetype, vat.typeofElement);
    assert(0);
  }
  // We construct a VA big enough to hold numElements.  This will
  // 'compress' the VA, if it was initially allocated much larger than
  // necessary.  Initialize_VA will actually allocate the next power of
  // two size necessary to hold numElements.
  //
  ReInitialize_VA( va, vat.numElements + allocate_additional_elements, 
		 vat.sizeofElement, vat.typeofElement);
  assert(va != NULL);
  EnableRange_VA(va,vat.numElements);
  {  
    char * elems = va->Elements;
    const size_t nsize = va->sizeofElement;
    const size_t nelem = va->numElements;

    assert(vat.numElements == 0 || (NULL != elems));
    assert(vat.sizeofElement == nsize);
    assert(vat.numElements == nelem);
    assert(nsize > 0);

    AS_UTL_safeRead(fp, elems, "ReInitializeFromFile_VA", nelem * nsize);
  }
}

void LoadFromFile_VA (FILE * const fp, 
		      VarArrayType *va, 
		      const char * const typetype, 
		      size_t growth_space){
    ReInitializeFromFile_VA( fp, va,typetype, growth_space);
}



void InitializeFromFile_VA
( FILE *fp,
  VarArrayType * const va,
  const char * const thetype,
  const size_t allocate_additional_elements
  ){
  FileVarArrayType vat;
  size_t nitems=0;
  assert(fp != NULL);
#ifdef DEBUG3
  fprintf(stderr,"sizeof(VarArrayType)=" F_SIZE_T " sizeof(FileVarArrayType)=" F_SIZE_T "\n",
          sizeof(VarArrayType), sizeof(FileVarArrayType));
#endif  
  /* we should use the safe read and write routines here. */
#if 0
  nitems = fread(&vat,sizeof(FileVarArrayType),1,fp);
  assert(1 == nitems);
#else
  nitems = fread(&vat.Elements,sizeof(uint64),1,fp);
  assert(1 == nitems);
  nitems = fread(&vat.sizeofElement,sizeof(uint64),1,fp);
  assert(1 == nitems);
  nitems = fread(&vat.numElements,sizeof(uint64),1,fp);
  assert(1 == nitems);
  nitems = fread(&vat.allocatedElements,sizeof(uint64),1,fp);
  assert(1 == nitems);
  nitems = fread(&vat.typeofElement,sizeof(char),VA_TYPENAMELEN,fp);
  assert(VA_TYPENAMELEN == nitems);
#ifdef DEBUG  
  fprintf(stderr,"vat.Elements=%p vat.sizeofElement=" F_SIZE_T " "
          "vat.numElements=" F_SIZE_T " vat.allocatedElements=" F_SIZE_T " "
          "vat.typeofElement=<%s>\n",
          vat.Elements, vat.sizeofElement,
          vat.numElements, vat.allocatedElements,
          vat.typeofElement );
#endif  
#endif
  assert(vat.numElements <= vat.allocatedElements);
  if(strncmp(vat.typeofElement,thetype,VA_TYPENAMELEN)){
    fprintf(stderr,"* Expecting array of type <%s> but read array of type <%s>\n",
	    thetype, vat.typeofElement);
    assert(0);
  }

  // We construct a VA big enough to hold numElements.  This will
  // 'compress' the VA, if it was initially allocated much larger than
  // necessary.  Initialize_VA will actually allocate the next power of
  // two size necessary to hold numElements.
  //
  Initialize_VA( va, vat.numElements + allocate_additional_elements, 
		 vat.sizeofElement, vat.typeofElement);
  assert(va != NULL);
  EnableRange_VA(va,vat.numElements);
  {  
    char * elems = va->Elements;
    const size_t nsize = va->sizeofElement;
    const size_t nelem = va->numElements;

    assert(vat.numElements == 0 || (NULL != elems));
    assert(vat.numElements == va->numElements);
    assert(vat.sizeofElement == va->sizeofElement);
    assert(nsize > 0);

    AS_UTL_safeRead(fp, elems, "InitializeFromFile_VA", nelem * nsize);
  }
}

size_t CopyToFile_VA(const VarArrayType * const va,FILE *fp){
  FileVarArrayType vat;
  size_t totalSize = 0;
  size_t nitems = 0;
  assert(fp != NULL);
  assert(va != NULL);

  /* Make a handle suitable for regression testing. */
  vat.Elements = 0;
  vat.sizeofElement = va->sizeofElement;
  vat.numElements = va->numElements;
  vat.allocatedElements = va->allocatedElements;
  strncpy(vat.typeofElement,va->typeofElement,VA_TYPENAMELEN);

  /* we should use the safe read and write routines here. */
  nitems = fwrite(&vat,sizeof(FileVarArrayType),1,fp);
  totalSize += sizeof(FileVarArrayType);
  assert(1 == nitems);
  { 
    const char * elems = va->Elements;
    const size_t nsize = va->sizeofElement;
    const size_t nelem = va->numElements;
    
    assert(nelem == 0 || NULL != elems);
    assert(nsize > 0);
    
    totalSize += nsize * nelem;

    AS_UTL_safeWrite(fp, elems, "CopyToFile_VA", nelem * nsize);
  }
  //  fprintf(stderr,"* CopyFileToVa appended " F_SIZE_T " bytes\n", totalSize);
  return totalSize;
}

void CheckFile_VA
(
 const VarArrayType * const va,
 FILE *fp
 ){
  // Check the contents of the variable length array with the contents
  // of the file based check point.
  FileVarArrayType vat;
  size_t nitems=0;
  assert(fp != NULL);
  assert(va != NULL);
  /* We should use the safe read and write routines here. */
  nitems = fread(&vat,sizeof(FileVarArrayType),1,fp);
  assert(1 == nitems);
  assert(vat.numElements <= vat.allocatedElements);
  nitems = strncmp(vat.typeofElement,va->typeofElement,VA_TYPENAMELEN);
  assert(nitems == 0);

  { 
    const char * const elems = va->Elements;
    const size_t nsize = va->sizeofElement;
    const size_t nelem = va->numElements;
    
    assert((nelem == 0) || (NULL != elems));
    assert(vat.sizeofElement == nsize);
    assert(vat.numElements == nelem);
    assert(nelem <= va->allocatedElements);
    assert(nsize > 0);

    //  XXX: This should read in pieces at a time, rather than the
    //  whole thing!
    //
    //  XXX: Better, we should convince ourselves that this does
    //  nothing when the writer properly catches errors.

    {
      char *tmp = (char *)malloc(nelem * nsize);
      size_t ii; 

      AS_UTL_safeRead(fp, tmp, "CheckFile_VA", nelem * nsize);

      for(ii=0;ii<nelem*nsize;ii++)
        assert(elems[ii] == tmp[ii]);

      free(tmp);
    }
  }
  return;
}

void ScatterInPlace_VA
(
 VarArrayType * const va,
 const size_t         nrange, // The length of the image of rank[].
 const size_t * const rank    // The rank of each item
 ){
  // Written by Clark Mobarry, 1999 Oct 14. 

  // This routine has the problem that it assumes that rank[] is a
  // correct permutation!!
  const size_t num  = va->numElements;
  const size_t mitems = max(nrange,(va->allocatedElements));
  const size_t size = va->sizeofElement;
  char * const old_array  = va->Elements;
  char * const new_array  = (char *) calloc(mitems,size);
  size_t io,in;
  assert(NULL != old_array);
  assert(NULL != new_array);
  // assert(nrange == num); // For a permutation 
  assert(nrange >= num); // Allow an unpack operation.
  for(io=0;io<num;io++) {
    in = rank[io];
    //assert(in >= 0);
    assert(in < nrange);
    memmove(new_array+in*size,old_array+io*size,size);
  }
  va->Elements = new_array;
  va->numElements = nrange;
  free(old_array);
}

void GatherInPlace_VA
(
 VarArrayType * const va,
 const size_t         nitems, // the length of indx[]
 const size_t * const indx    // the list of the new order
 ){
  // Written by Clark Mobarry, 1999 Oct 14. 

  // This routine has the problem that it assumes that indx[] is a
  // correct permutation!!
  const size_t num  = va->numElements;
  const size_t size = va->sizeofElement;
  char * const old_array = va->Elements;
  char * const new_array = (char *) calloc((va->allocatedElements),size);
  size_t io,in;
  assert(NULL != old_array);
  assert(NULL != new_array);
  assert(num >= nitems);
  for(in=0;in<nitems;in++) {
    io = indx[in];
    //assert(io >= 0);
    assert(io < num);
    memmove(new_array+in*size,old_array+io*size,size);
  }
  va->numElements = nitems;
  va->Elements = new_array;
  free(old_array);
}

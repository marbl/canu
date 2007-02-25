
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
 *
 * Once a VA (variable array) has been defined, using the VA_DEF(<Type>)
 * macro, it may be manipulated as follows:
 *
 * A typedef is defined VarArray<Type> for the VA, and is accessible via the
 * VA_TYPE(<Type>) macro.
 *
 * Allocate a new VA_Type(<Type>)
 * VarArray<Type> * CreateVA_<Type>(size_t numElements, size_t sizeofElement)
 *
 * Delete a VA_Type(<Type>)
 * void DeleteVA_<Type>(VarArray<Type> *va);
 *
 * Get a pointer to an element of an array, returning NULL if
 * index is out of bounds for the array.
 * <Type> *GetVA_<Type>(VarArray<Type> *va, size_t index); 
 *
 * Copy the data into the array, enlarging the array as necessary.
 * This operation also expands the number of elements in the array.
 * void SetVA_<Type>(VarArray<Type> *va, size_t index, <Type> *val); 
 *
 *
 * Get the number of elements in the array
 * size_t GetNumVA_<Type>(VarArray<Type> *va); 
 *
 * Get the number of elements for which space is allocated.
 * size_t GetAllocatedVA_<Type>(VarArray<Type> *va);
 *
 * When space is allocated, it is initialized to zero. 
 *
 * CMM 1999/04/12: Added CopyToFile_VA and CreateFromFile_VA.
 *
 * KAR 1999/06/30: Added CreateFromArray_VA.
 */

#ifndef AS_UTL_VAR_H
#define AS_UTL_VAR_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"

#define VA_TYPENAMELEN 32 /* The number of significant characters used to 
			     distinguish a array element type. */
typedef struct {
  char *Elements; /* The Data pointer. Must be cast to the appropriate type */
  size_t sizeofElement;   /* The size in bytes of the appropriate type. */
  size_t numElements;     
  size_t allocatedElements;
  char typeofElement[VA_TYPENAMELEN]; /* The name of the data type of 
					  each element. */
} VarArrayType;
  

void Initialize_VA(VarArrayType * const va,
                   const size_t arraySize,
                   const size_t sizeofElement,
                   const char * const thetype);

void InitializeFromArray_VA(VarArrayType * const va,
                            const size_t arraySize,
                            const size_t sizeofElement,
                            const char * const thetype,
                            const void * const data);

void ReInitializeFromArray_VA(VarArrayType * const va,
                              const size_t arraySize,
                              const size_t sizeofElement,
                              const char * const thetype,
                              const void * const data);

void InitializeFromFile_VA(FILE *fp,
                           VarArrayType * const va,
                           const char * const thetype,
                           const size_t allocate_additional_elements);

void ReInitializeFromFile_VA(FILE *fp,
                             VarArrayType * const va,
                             const char * const thetype,
                             const size_t allocate_additional_elements);

size_t CopyToFile_VA(const VarArrayType * const va,
                     FILE *fp);

void CheckFile_VA(const VarArrayType * const va,
                  FILE *fp);

int MakeRoom_VA(VarArrayType * const va,
                const size_t         maxElements,
                const int            pad_to_a_power_of_two);


// ResetToRange_VA grows or shrinks the size of the VA as a function of index
//
void ResetToRange_VA(VarArrayType * const va, const size_t index);

#define Reset_VA(V)   ResetToRange_VA((V), 0)

void Clear_VA(VarArrayType * const va);

#define GetElement_VA(V, I)  ((I) < (V)->numElements ? ((V)->Elements + ((size_t)(I) * (size_t)((V)->sizeofElement))) : NULL)

// EnableRange grows size of VA, as necessary, to size index
void EnableRange_VA ( VarArrayType * const va, const size_t index);
void SetElement_VA ( VarArrayType * const va, const size_t index, const void * const data);
void SetRange_VA ( VarArrayType * const va, const size_t index, const size_t num_elements, const void * const data);
void Concat_VA ( VarArrayType * const va, const VarArrayType * const vb);

#define GetNumElements_VA(V)        (size_t)((V)->numElements)
#define GetAllocatedElements_VA(V)  (size_t)((V)->allocatedElements)
#define GetsizeofElement_VA(V)      (size_t)((V)->sizeofElement)


// *************************************************************
//
// Functions that allocate from the heap and initialize a VarArrayType.
//

static VarArrayType *Create_VA(const size_t arraySize,
			       const size_t sizeofElement,
			       const char * const thetype)
{
  /* 
     Return a handle to a new variable length array.
     
     arraySize: The initial number of elements. This can be zero.
     sizeofElement: The size in bytes of the elements.
     thetype: A character string identifying the type.
  */
     
  VarArrayType * const va = (VarArrayType *)safe_malloc(sizeof(VarArrayType));
  Initialize_VA( va, arraySize, sizeofElement, thetype);
  return va;
}

static VarArrayType *CreateFromArray_VA(const void * const data,
                                        const size_t arraySize,
                                        const size_t sizeofElement,
                                        const char * const thetype) {
  /*
     Return a handle to a new variable length array, initialized
     with a copy of the input array 'data'.
     
     arraySize: The initial number of elements. This can be zero.
     sizeofElement: The size in bytes of the elements.
     thetype: A character string identifying the type.
  */
     
  VarArrayType *va;
  va = (VarArrayType *)safe_malloc(sizeof(VarArrayType));
  InitializeFromArray_VA( va, arraySize, sizeofElement, thetype, data);
  return va;
}

static VarArrayType *CreateFromFile_VA(FILE *fp,const char * const thetype,size_t space_in_elements_for_growth){
  VarArrayType *va;
  va = (VarArrayType *)safe_malloc(sizeof(VarArrayType));
  InitializeFromFile_VA( fp, va, thetype, space_in_elements_for_growth);
  return va;
}

void LoadFromFile_VA (FILE * const fp, 
		      VarArrayType *va, 
		      const char * const typetype, 
		      size_t growth_space);


static VarArrayType *Clone_VA(const VarArrayType *va){
  VarArrayType *newva;
  newva = (VarArrayType *)safe_malloc(sizeof(VarArrayType));
  InitializeFromArray_VA(newva, va->numElements, va->sizeofElement, va->typeofElement, va->Elements);
  return newva;
}

// Like Clone, except the vato VA is recycled, rather than allocated
#define ReuseClone_VA(T, F)  ReInitializeFromArray_VA(T, (F)->numElements, (F)->sizeofElement, (F)->typeofElement, (F)->Elements)

#define Delete_VA(V)         { Clear_VA(V); safe_free(V); }

#define GetMemorySize_VA(V)   (size_t)(((V) ? (V)->allocatedElements * (V)->sizeofElement : 0))


static size_t ReportMemorySize_VA(VarArrayType * const va,
                                  const char * const name,
                                  FILE * stream ) {
  size_t numElements       = (NULL == va ? 0 : va->numElements);
  size_t allocatedElements = (NULL == va ? 0 : va->allocatedElements);
  size_t sizeofElement     = (NULL == va ? 0 : va->sizeofElement);
  size_t memorySize        = allocatedElements * sizeofElement;

  assert(NULL != name);
  assert(NULL != stream);

  fprintf(stream,
          "VA"
          " %10" F_SIZE_TP " bytes "
          " %10" F_SIZE_TP " elements active"
          " %10" F_SIZE_TP " elements allocated"
          " %5" F_SIZE_TP " bytes per element"
          " for %s\n"
          ,memorySize
          ,numElements
          ,allocatedElements
          ,sizeofElement
          ,name
          );

  return memorySize;
}



// *************************************************************
//
// The user interface:
//

#define VA_TYPE(Type) VarArray ## Type


#define VA_DEF(Type)\
typedef VarArrayType VarArray ## Type ;\
static VA_TYPE(Type) InitVA_ ## Type ( const size_t numElements) { \
     VA_TYPE(Type) va ; \
     Initialize_VA( &va, numElements, sizeof(Type), #Type);\
     return (va); }\
static void ClearVA_ ## Type (VA_TYPE(Type) * const va){\
     Clear_VA(va);}\
static VA_TYPE(Type) * CreateVA_ ## Type (const size_t numElements){\
     return ( (VA_TYPE(Type) *)Create_VA(numElements, sizeof(Type), #Type)); }\
static void DeleteVA_ ## Type (VA_TYPE(Type) *va){\
     Delete_VA(va); }\
static void ConcatVA_ ## Type ( VA_TYPE(Type) *va, const VA_TYPE(Type) * const vb){\
     Concat_VA(va,vb); }\
static Type *GetVA_ ## Type (const VA_TYPE(Type) * const va, size_t index){\
     return ( (Type *)GetElement_VA(va,index));\
}\
static size_t GetVAIndex_ ## Type (const VA_TYPE(Type) * const va, Type *elem){\
     size_t index = (size_t)(elem - GetVA_ ## Type (va, 0));\
     assert((size_t)elem >= (size_t)GetVA_##Type (va,0));\
     assert(index <= va->numElements);\
     return index;\
}\
static void ResetVA_ ## Type (VA_TYPE(Type) * const va){\
      Reset_VA(va);\
}\
static void ResetToRangeVA_ ## Type (VA_TYPE(Type) * const va, size_t index){\
      ResetToRange_VA(va,index);\
}\
static void EnableRangeVA_ ## Type (VA_TYPE(Type) * const va, size_t index){\
      EnableRange_VA(va,index);\
}\
static void SetVA_ ## Type (VA_TYPE(Type) * const va, \
			 const size_t index, \
			 const Type * const data){\
      SetElement_VA(va,index,data);\
}\
static void SetRangeVA_ ## Type (VA_TYPE(Type) * const va, \
			 const size_t index, \
			 const size_t num_elements, \
			 const Type * const data){\
      SetRange_VA(va,index,num_elements,data);\
}\
static void AppendVA_ ## Type (VA_TYPE(Type) * const va, \
			       const Type * const data){ \
      SetElement_VA(va,GetNumElements_VA(va),data);\
}\
static void AppendRangeVA_ ## Type (VA_TYPE(Type) * const va, \
			       const size_t num_elements, const Type * const data){ \
      SetRange_VA(va,GetNumElements_VA(va),num_elements,data);\
}\
static VA_TYPE(Type) * CreateFromArrayVA_ ## Type (const void * const data, size_t numElements){\
     return ( (VA_TYPE(Type) *)CreateFromArray_VA(data, numElements, sizeof(Type), #Type));\
}\
static size_t GetNumVA_ ## Type (const VA_TYPE(Type) * const va){\
  return GetNumElements_VA(va);\
}\
static size_t GetAllocatedVA_ ## Type (const VA_TYPE(Type) * const va){\
  return GetAllocatedElements_VA(va);\
}\
static VA_TYPE(Type) * CreateFromFileVA_ ## Type (FILE * const fp,size_t growth_space){\
 return (VA_TYPE(Type) *)CreateFromFile_VA(fp, #Type, growth_space);\
}\
static size_t CopyToFileVA_ ## Type (const VA_TYPE(Type) * const va,FILE *fp){\
 return CopyToFile_VA(va,fp);\
}\
static Type *Get ## Type (const VA_TYPE(Type) * const va, size_t index){\
     return ( (Type *)GetElement_VA(va,index));\
}\
static void Reset ## Type (VA_TYPE(Type) * const va){\
      Reset_VA(va);\
}\
static void ResetToRange_ ## Type (VA_TYPE(Type) * const va, size_t index){\
      ResetToRange_VA(va,index);\
}\
static void Set ## Type (VA_TYPE(Type) * const va, \
			 const size_t index, \
			 const Type * const data){\
      SetElement_VA(va,index,data);\
}\
static void Append ## Type (VA_TYPE(Type) * const va, const Type * const data){\
      SetElement_VA(va,GetNumElements_VA(va),data);\
}\
static void AppendRange ## Type (VA_TYPE(Type) * const va, \
                          const size_t num_elements, const Type * const data){\
      SetRange_VA(va,GetNumElements_VA(va),num_elements,data);\
}\
static size_t GetNum ## Type ##s(const VA_TYPE(Type) * const va){\
  return GetNumElements_VA(va);\
}\
static size_t GetAllocated ## Type ##s(const VA_TYPE(Type) * const va){\
  return GetAllocatedElements_VA(va);\
}\
static void LoadFromFileVA_ ## Type (FILE * const fp, VA_TYPE(Type) *va, size_t growth_space){\
 LoadFromFile_VA(fp, va, #Type, growth_space);\
}


/*** Parametrized Stack Type Implemented with VAs ***/

#define STACK_DEF(type)\
typedef struct{\
  VA_TYPE(type) *stack;\
  int top;\
}Stack_##type;\
\
static Stack_##type *CreateStack_##type (int size){\
  Stack_##type *new = (Stack_##type *)safe_malloc(sizeof(Stack_##type ));\
  new->stack = CreateVA_##type (size);\
  new->top = NULLINDEX; \
  return new;\
}\
\
static void DeleteStack_##type(Stack_##type *stack){\
  DeleteVA_##type (stack->stack);\
  safe_free(stack);\
}\
static void ResetStack_##type(Stack_##type *stack){\
  ResetVA_##type (stack->stack);\
  stack->top = NULLINDEX;\
}\
static void PushStack_##type (Stack_##type *stack, void *item){\
  Set##type (stack->stack, ++(stack->top), &item);\
}\
\
static void *PopStack_##type (Stack_##type *stack){\
  void **item;\
\
  if(stack->top < 0)\
    return NULL;\
  item = Get##type (stack->stack, (stack->top)--);\
  return *item;\
}


#endif // AS_UTL_VAR_H

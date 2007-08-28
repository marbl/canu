
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

#ifndef AS_UTL_VAR_H
#define AS_UTL_VAR_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_UTL_alloc.h"

//  The number of significant characters used to distinguish a array
//  element type.
#define VA_TYPENAMELEN 32

typedef struct {
  char      *Elements;        // The Data pointer. Must be cast to the appropriate type
  size_t     sizeofElement;   // The size in bytes of the appropriate type
  size_t     numElements;     
  size_t     allocatedElements;
  char typeofElement[VA_TYPENAMELEN]; // The name of the data type of each element
} VarArrayType;


int
MakeRoom_VA(VarArrayType * const va,
            const size_t         maxElements,
            const int            pad_to_a_power_of_two);


VarArrayType *
Create_VA(const size_t arraySize,
          const size_t sizeofElement,
          const char * const thetype);

void
Clear_VA(VarArrayType * const va);

void
Concat_VA(VarArrayType * const va, const VarArrayType * const vb);

void
ResetToRange_VA(VarArrayType * const va, const size_t index);

void
EnableRange_VA(VarArrayType * const va, const size_t index);

void
SetElement_VA(VarArrayType * const va, const size_t index, const void * const data);

void
SetRange_VA(VarArrayType * const va, const size_t index, const size_t num_elements, const void * const data);

VarArrayType *
Clone_VA(const VarArrayType *fr);

void
ReuseClone_VA(VarArrayType *to, const VarArrayType *fr);

void
LoadFromFile_VA(FILE * const fp, VarArrayType *va);

VarArrayType *
CreateFromFile_VA(FILE *fp,
                  const char * const thetype,
                  size_t space_in_elements_for_growth);

size_t
CopyToFile_VA(const VarArrayType * const va, FILE *fp);


#define Delete_VA(V)         { Clear_VA(V); safe_free(V); }

#define GetMemorySize_VA(V)   (size_t)(((V) ? (V)->allocatedElements * (V)->sizeofElement : 0))

#define GetElement_VA(V, I)         ((I) < (V)->numElements ? ((V)->Elements + ((size_t)(I) * (size_t)((V)->sizeofElement))) : NULL)
#define GetNumElements_VA(V)        ((V)->numElements)
#define GetAllocatedElements_VA(V)  ((V)->allocatedElements)
#define GetsizeofElement_VA(V)      ((V)->sizeofElement)


static
size_t
ReportMemorySize_VA(VarArrayType * const va,
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
      ResetToRange_VA(va, 0);\
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
static size_t GetNumVA_ ## Type (const VA_TYPE(Type) * const va){\
  return GetNumElements_VA(va);\
}\
static size_t GetAllocatedVA_ ## Type (const VA_TYPE(Type) * const va){\
  return GetAllocatedElements_VA(va);\
}\
static VA_TYPE(Type) * CreateFromFileVA_ ## Type (FILE * const fp,size_t growth_space){\
 return (VA_TYPE(Type) *)CreateFromFile_VA(fp, #Type, growth_space);\
}\
static void LoadFromFileVA_ ## Type (FILE * const fp,VA_TYPE(Type) *va){\
 LoadFromFile_VA(fp, va);\
}\
static size_t CopyToFileVA_ ## Type (const VA_TYPE(Type) * const va,FILE *fp){\
 return CopyToFile_VA(va,fp);\
}\
static Type *Get ## Type (const VA_TYPE(Type) * const va, size_t index){\
     return ( (Type *)GetElement_VA(va,index));\
}\
static void Reset ## Type (VA_TYPE(Type) * const va){\
      ResetToRange_VA(va, 0);\
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
}


/*** Parametrized Stack Type Implemented with VAs ***/

#define STACK_DEF(type)\
typedef struct{\
  VA_TYPE(type) *stack;\
  int top;\
}Stack_##type;\
\
static Stack_##type *CreateStack_##type (int size){\
  Stack_##type *new = (Stack_##type *)safe_calloc(1, sizeof(Stack_##type ));\
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

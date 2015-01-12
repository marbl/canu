
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

#ifndef AS_UTL_ALLOC_H
#define AS_UTL_ALLOC_H

static const char *rcsid_AS_UTL_ALLOC_H = "$Id$";

//
//  The safe_*alloc routines are the same as the normal routines,
//  except they print a standard message and assert if memory cannot
//  be allocated.
//
//  safe_free is somewhat useless.  It frees the memory, then sets the
//  pointer to NULL.  It's useless because many times the pointer we
//  set to NULL is a local copy, e.g.:
//
//     void freeSomeStructure(SomeStructure *p) {
//       safe_free(p->data);
//       safe_free(p);
//     }
//

void *safe_calloc(size_t num, size_t len);
void *safe_malloc(size_t len);
void *safe_realloc(void *q, size_t len);
void  safe_free2(void *);

#define safe_free(Q) { safe_free2(Q); Q = NULL; }

//  And, thanks, GNU.  strdup() (and, sigh, probably lots others) are
//  implemented as a macro that calls free(), which then gets expanded
//  into our bogus function.
//
//  And, thanks, Apple.  Your iostream implementation uses malloc. (!)
//
#ifndef BROKEN_CLANG_OpenMP
#ifndef X86_GCC_LINUX

#define malloc(X)     use_safe_malloc_instead(X)
#define calloc(X,Y)   use_safe_calloc_instead(X,Y)
#define realloc(X,Y)  use_safe_realloc_instead(X,Y)
#define free(X)       use_safe_free_instead(X)

#endif
#endif

#if 0
template<class T>
inline void zero(T* ptr, size_t num)
{
	for (size_t a = 0; a < num; a++)
	{
		ptr[a] = NULL;
	}
}
#endif




template<typename TT, typename LL>
void
duplicateArray(TT*& to, LL &toLen, LL &toMax, TT *fr, LL frLen, LL frMax, bool forceAlloc=false) {

  if ((toMax < frLen) || (forceAlloc)) {
    delete [] to;

    toLen = frLen;
    toMax = toLen;

    to = new TT [toMax];
  }

  memcpy(to, fr, sizeof(TT) * toLen);
}

const uint32  resizeArray_doNothing = 0x00;
const uint32  resizeArray_copyData  = 0x01;
const uint32  resizeArray_clearNew  = 0x02;


template<typename TT, typename LL>
void
increaseArray(TT*& array, uint64 arrayLen, LL &arrayMax, uint64 increment) {

  if (arrayLen + increment <= arrayMax)
    return;

  resizeArray(array, arrayLen, arrayMax, 2 * arrayMax, resizeArray_copyData);
}


template<typename T1, typename T2, typename LL>
void
increaseArrayPair(T1*& array1, T2*& array2, uint64 arrayLen, LL &arrayMax, uint64 increment) {

  if (arrayLen + increment <= arrayMax)
    return;

  resizeArrayPair(array1, array2, arrayLen, arrayMax, 2 * arrayMax, resizeArray_copyData);
}


template<typename TT, typename LL>
void
resizeArray(TT*& array, uint64 arrayLen, LL &arrayMax, uint64 newMax, uint32 op=resizeArray_copyData) {

  if (newMax <= arrayMax)
    return;

  arrayMax = newMax;

  TT *copy = new TT [arrayMax];

  if (op & resizeArray_copyData)
    memcpy(copy, array, sizeof(TT) * arrayLen);

  delete [] array;
  array = copy;

  if (op & resizeArray_clearNew)
    memset(array + sizeof(TT) * arrayLen, 0, sizeof(TT) * (newMax - arrayLen));
}


template<typename T1, typename T2, typename LL>
void
resizeArrayPair(T1*& array1, T2*& array2, uint64 arrayLen, LL &arrayMax, LL newMax, uint32 op=resizeArray_copyData) {

  if (newMax <= arrayMax)
    return;

  arrayMax = newMax;

  T1 *copy1 = new T1 [arrayMax];
  T2 *copy2 = new T2 [arrayMax];

  if (op & resizeArray_copyData) {
    memcpy(copy1, array1, sizeof(T1) * arrayLen);
    memcpy(copy2, array2, sizeof(T2) * arrayLen);
  }

  delete [] array1;
  delete [] array2;

  array1 = copy1;
  array2 = copy2;

  if (op & resizeArray_clearNew) {
    memset(array1 + sizeof(T1) * arrayLen, 0, sizeof(T1) * (newMax - arrayLen));
    memset(array2 + sizeof(T2) * arrayLen, 0, sizeof(T2) * (newMax - arrayLen));
  }
}



#endif // AS_UTL_ALLOC_H

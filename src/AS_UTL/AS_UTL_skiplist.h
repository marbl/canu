
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
   CVS_ID:  $Id: AS_UTL_skiplist.h,v 1.6 2008-06-27 06:29:21 brianwalenz Exp $
 *********************************************************************/

/********************************************************************/
/* dynamic sorted sequences using skiplist (Pugh)
 * the insert, delete and lookup operations take exspected time O(log n)
 * with an exponential tail estimate.
 *
 *     Knut Reinert
 *     Juli 7 1999
 *
 *
 * Once a SL (skiplist) has been defined, using the SL_DEF(<Type>)
 * macro, it may be manipulated as follows:
 *
 * A typedef is defined SkipListType for the SL, and is accessible via the
 * SL_TYPE(Type) macro.
 *
 * Allocate a new SL_TYPE(Type)
 * SL_TYPE(Type) *CreateSL_Type(int freeData, SLF free)
 *
 * Free a SL_Type(Type)
 * void FreeSL_Type(SL_TYPE(Type) *sl);
 *
 * Lookup a key in the skiplist returning the skiplist item
 * element with the biggest key less than or equal to key
 * sl_item LookupSL_Type(double key,SL_TYPE(Type)* sl);
 *
 * Insert a key,value pair in the skiplist returning the inserted
 * skiplist item or the item that was already inserted with that key
 * sl_item InsertSL_Type(double key,Type,SL_TYPE(Type),SL_TYPE(Type)* sl);
 *
 * Delete a key,value pair in the skiplist. The function does nothing
 * if the element is not present.
 * sl_item DeleteSL_Type(double key, SL_TYPE(Type) *sl);
 *
 * Return the minimum resp. maximum item
 * sl_item MinSL_Type(SL_TYPE(Type) *sl);
 * sl_item MaxSL_Type(SL_TYPE(Type) *sl);
 *
 * Get the number of elements in the array
 * size_t GetNumSL_Type(SL_TYPE(Type) *sl);
 *
 * Look in the file AS_UTL_skiplist_test for examples of how to use
 * Skiplist
 */

#ifndef AS_UTL_SKIPLIST_H
#define AS_UTL_SKIPLIST_H


#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<float.h>

#include "AS_global.h"
#include <assert.h>
#include "AS_UTL_rand.h"



#define minf -FLT_MAX
#define pinf  FLT_MAX

typedef double keyType;
typedef void*  valueType;
typedef void (*SLF)(void* v);

typedef struct element *sl_item;
struct element{
  	keyType   key;
        valueType value;
	sl_item pred;
	sl_item succ;
	sl_item down;
	sl_item up;
};


typedef struct sl{
  SLF free_value;
  int freeData;
  sl_item head;
  sl_item tail;
  int no_of_levels;
  int no_of_elements;
} SkipList;


/* interface functions */

void	 Delete_SL(keyType key, SkipList *sl);
sl_item  Insert_SL(keyType key, valueType val, SkipList *sl);
sl_item	 Lookup_SL(keyType key, SkipList *sl);
SkipList *Create_SL(int fd, SLF free_value);
void	 Free_SL(SkipList *sl);
int      GetNum_SL(SkipList *sl);
void	 Print_SL(SkipList *sl);
sl_item  Max_SL(SkipList*sl);
sl_item  Min_SL(SkipList*sl);

/* Macros */
#define new_item (sl_item) safe_malloc(sizeof(struct element))

#define SL_TYPE(Type) SkipList ## Type

#define SL_DEF(Type)\
typedef SkipList SkipList ## Type;\
static SL_TYPE(Type) * CreateSL_ ## Type (int freeData, SLF free_value){\
     return ((SL_TYPE(Type)*) Create_SL(freeData,free_value));\
}\
static void FreeSL_ ## Type (SL_TYPE(Type) *sl){\
     Free_SL(sl);\
}\
static sl_item LookupSL_ ## Type (keyType key, SL_TYPE(Type) *sl){\
     return Lookup_SL(key,sl);\
}\
static sl_item InsertSL_ ## Type (keyType key, Type* val, SL_TYPE(Type) *sl){\
     return Insert_SL(key,val,sl);\
}\
static void DeleteSL_ ## Type (keyType key, SL_TYPE(Type) *sl){\
     Delete_SL(key,sl);\
}\
static size_t GetNumSL_ ## Type (SL_TYPE(Type) *sl){\
  return GetNum_SL(sl);\
}\
static void PrintSL_ ## Type (SL_TYPE(Type) *sl){\
  Print_SL(sl);\
}\
static sl_item MinSL_ ## Type (SL_TYPE(Type) *sl){\
     return Min_SL(sl);\
}\
static sl_item MaxSL_ ## Type (SL_TYPE(Type) *sl){\
     return Max_SL(sl);\
}

#endif









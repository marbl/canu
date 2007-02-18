
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

//  We explicitly do not include AS_UTL_alloc.h here, because it
//  redefines malloc(), calloc(), realloc() and free() to be errors.
//  We want everyone to use the safe_*() versions supplied here.
//
//#include "AS_UTL_alloc.h"

//  We want to include AS_global.h to get the F_SIZE_T definition, but that
//  includes AS_UTL_alloc.h, so we pretend we've already included it.
//
#define AS_UTL_ALLOC_H
#include "AS_global.h"


// Allocate and return a pointer to an array of  num  elements of
// len  bytes each.  All are set to 0.  Exit if fai.
//
void *
safe_calloc(size_t num, size_t len) {
  void  *p;

   p = calloc(num, len);
   if (p == NULL) {
     fprintf(stderr, "Could not calloc memory ("F_SIZE_T" * "F_SIZE_T" bytes = "F_SIZE_T")\n",
             num, len, num*len);
     assert(p != NULL);
   }
   return(p);
}



// Allocate and return a pointer to len bytes of memory.
// Len  bytes each.  Exit if fail.
//
void *
safe_malloc(size_t len) {
  void  *p;

  p = malloc(len);
  if (p == NULL) {
    fprintf(stderr, "Could not malloc memory ("F_SIZE_T" bytes)\n", len);
    assert(p != NULL);
  }

#undef TRASH_MEMORY_FIRST
#ifdef TRASH_MEMORY_FIRST
  memset(p, 0xff, len);
#endif
  
  return(p);
}




// Reallocate memory for q to len  bytes and return a pointer
// to the new memory.  Exit if fail.
//
void *
safe_realloc(void *q, size_t len) {
  void  *p;

  p = realloc(q, len);
  if (p == NULL) {
    fprintf(stderr, "Could not realloc memory ("F_SIZE_T" bytes)\n", len);
    assert(p != NULL);
  }
  
  return(p);
}



void
safe_free2(void *q) {
  free(q);
}

#include <stdio.h>
#include <stdlib.h>
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
  memset(p, 0, len);

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

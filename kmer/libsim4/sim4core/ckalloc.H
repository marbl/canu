#ifndef CKALLOC_H
#define CKALLOC_H

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#include "libbritypes.h"

//#define CHECK_MEMORY

inline
void *
ckalloc(size_t amount) {

  if (amount == 0)
    amount = 8;

  void *p = malloc(amount);

  if (p == NULL) {
    fprintf(stderr, "Can't allocate "u64bitFMT" bytes.\n", (u64bit)amount);
    kill(getpid(), SIGKILL);
  }

#ifdef CHECK_MEMORY
  fprintf(stderr, "ckalloc "u64bitHEX" "u64bitFMT"\n", (u64bit)p, (u64bit)amount);
#endif

  return(p);
}



#ifdef CHECK_MEMORY
#define ckfree(P) { fprintf(stderr, "ckfree "u64bitHEX" at %s,%d\n", (u64bit)(P), __FILE__, __LINE__); free(P); }
#else
#define ckfree(P) free(P)
#endif


#endif  //  CKALLOC_H

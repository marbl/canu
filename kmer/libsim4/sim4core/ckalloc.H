#ifndef CKALLOC_H
#define CKALLOC_H

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#include "bri.h"

//#define CHECK_MEMORY

//#define WITH_PALLOC



#ifdef WITH_PALLOC

inline
void*
ckalloc(size_t amount) {
  return(palloc(amount));
}


#define ckfree(P)

#else

inline
void *
ckalloc(size_t amount) {

#ifdef __APPLE__
  if (amount == 0)
    amount = 16;
#else
  if (amount == 0)
    amount = 8;
#endif

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
#endif  //  CHECK_MEMORY

#endif  //  WITH_PALLOC

#endif  //  CKALLOC_H

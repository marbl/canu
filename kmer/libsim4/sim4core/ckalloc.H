#ifndef CKALLOC_H
#define CKALLOC_H

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

//#define CHECK_MEMORY

inline
void *
ckalloc(size_t amount) {

  if (amount == 0)
    amount = 8;

  void *p = malloc(amount);

  if (p == NULL) {
#ifdef TRUE64BIT
    fprintf(stderr, "Can't allocate %lu bytes.\n", amount);
#else
    fprintf(stderr, "Can't allocate %d bytes.\n", amount);
#endif
    kill(getpid(), SIGKILL);
  }

#ifdef CHECK_MEMORY
#ifdef TRUE64BIT
  fprintf(stderr, "ckalloc 0x%016lx %lu\n", p, amount);
#else
  fprintf(stderr, "ckalloc 0x%08lx %ld\n", p, amount);
#endif
#endif

  return(p);
}


#ifdef CHECK_MEMORY

#define ckfree(P) { fprintf(stderr, "ckfree 0x%08lx at %s,%d\n", (P), __FILE__, __LINE__); free(P); }

#else

#define ckfree(P) free(P)

#endif


#endif  //  CKALLOC_H

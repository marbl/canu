#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "posix.H"
#include "searchGENOME.H"

void*
loaderThread(void *) {
  encodedQuery         *Q = 0L;
  seqInCore            *B = 0L;

  try {
    B = config._qsFASTA->getSequenceInCore();
  } catch (std::bad_alloc) {
    fprintf(stderr, "loaderThread()-- Failed to load next query sequence\ncaught bad_alloc in %s at line %d\n", __FILE__, __LINE__);
    exit(1);
  }

  if (B) {
    Q = new encodedQuery(B, config._merSize);
    delete B;
  }

  return(Q);
}

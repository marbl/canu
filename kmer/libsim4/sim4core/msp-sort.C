#include <stdio.h>
#include <stdlib.h>
#include "mspManager.H"

static
int
mspManager_msp_compare(const void *A, const void *B) {
  msp const  *a = (msp const *)A;
  msp const  *b = (msp const *)B;

  if (a->pos2 < b->pos2)
    return(-1);

  if (a->pos2 > b->pos2)
    return(1);

  if (a->pos1 < b->pos1)
    return(-1);

  if (a->pos1 > b->pos1)
    return(1);

  return(0);
}

void
mspManager::sort(void) {
  qsort(_allMSPs, _numMSPs, sizeof(struct msp), mspManager_msp_compare);
}

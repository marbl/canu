#include <stdio.h>
#include <stdlib.h>
#include "mspManager.H"

//
//  XXX:  This might be negated.
//

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



#if 0

//  sort_msps();    /* sort in order of mp->pos2, in the shorter seq */


/* smaller - determine ordering relationship between two MSPs */
int
Sim4::smaller(int i, int j) {
  Msp_ptr ki = msp[i], kj = msp[j];
  
  if (ki->pos2 > kj->pos2) return i;
  if (ki->pos2 < kj->pos2) return j;
  
  return ((ki->pos1 >= kj->pos1) ? i : j);
}


/* heapify - core procedure for heapsort */
void
Sim4::heapify(int i, int last) {
  int lim = (last-1)/2, left_son, small_son;
  Msp_ptr mp;

  while (i <= lim) {
    left_son = 2*i + 1;
    if (left_son == last)
      small_son = left_son;
    else
      small_son = smaller(left_son, left_son+1);
    if (smaller(i, small_son) == small_son) {
      mp = msp[i];
      msp[i] = msp[small_son];
      msp[small_son] = mp;
      i = small_son;
    } else
      break;
  }
}

/* sort_msps - order database sequence for printing */
void
Sim4::sort_msps(void) {
  int i;
  Msp_ptr mp;

  for (i = (numMSPs/2) - 1; i >= 0; --i)
    heapify(i, (int) numMSPs-1);
  for (i = numMSPs-1; i > 0; --i) {
    mp = msp[0];
    msp[0] = msp[i];
    msp[i] = mp;
    if (i > 1)
      heapify(0, i-1);
  }
}

#endif

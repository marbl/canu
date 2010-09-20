#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "bio++.H"

#include "sim4polishList.H"
#include "sim4polishReader.H"

sim4polishList::sim4polishList() {
  len  = 0;
  max  = 4;
  list = new sim4polish* [max];
}

sim4polishList::sim4polishList(char const *filename) {
  len  = 0;
  max  = 4;
  list = new sim4polish* [max];

  sim4polishReader *R = new sim4polishReader(filename);
  sim4polish       *p = 0L;

  while (R->nextAlignment(p))
    push(p);

  delete R;
}

sim4polishList::~sim4polishList() {
  for (u32bit i=0; i<len; i++)
    delete list[i];
  delete [] list;
}

void
sim4polishList::push(sim4polish *p) {

  if (p == 0L)
    return;

  if (len >= max) {
    max *= 2;
    sim4polish **l = new sim4polish* [max];
    memcpy(l, list, len * sizeof(sim4polish*));
    delete [] list;
    list = l;
  }

  list[len++] = p;
}

void
sim4polishList::remove(u32bit i) {

  if (i >= len)
    return;

  delete list[i];

  len--;
  for ( ; i < len; i++)
    list[i] = list[i+1];
}


void
sim4polishList::sortBycDNAIID(void) {
  qsort(list, len, sizeof(sim4polish *), s4p_estIDcompare);
}

void
sim4polishList::sortByGenomicIID(void) {
  qsort(list, len, sizeof(sim4polish *), s4p_genIDcompare);
}


void
sim4polishList::filterByQuality(u32bit minI, u32bit minC) {
  u32bit save = 0;
  u32bit next = 0;

  while (next < len) {
    if ((list[next]->_percentIdentity  >= minI) &&
        (list[next]->_querySeqIdentity >= minC)) {
      list[save++] = list[next++];
    } else {
      delete list[next];
      list[next++] = 0L;
    }
  }

  len = save;
}

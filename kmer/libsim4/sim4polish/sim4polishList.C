#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "bri++.H"
#include "sim4polish.h"
#include "sim4polishList.H"


sim4polishList::sim4polishList() {
  len  = 0;
  max  = 4;
  list = new sim4polish* [max];
}

sim4polishList::sim4polishList(char const *filename) {
  len  = 0;
  max  = 4;
  list = new sim4polish* [max];

  errno=0;
  FILE *F = fopen(filename, "r");
  if (errno) {
    fprintf(stderr, "sim4polishList()--  Can't open '%s' for reading\n%s\n", filename, strerror(errno));
    exit(1);
  }

  while (!feof(F))
    push(s4p_readPolish(F));

  fclose(F);
}

sim4polishList::~sim4polishList() {
  for (u32bit i=0; i<len; i++)
    s4p_destroyPolish(list[i]);
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
    if ((list[next]->percentIdentity  >= minI) &&
        (list[next]->querySeqIdentity >= minC)) {
      list[save++] = list[next++];
    } else {
      s4p_destroyPolish(list[next]);
      list[next++] = 0L;
    }
  }

  len = save;
}

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "util++.H"

//
//  mbri && CC -g -o test-bigQueue test-bigQueue.C -L. -lutil && ./bigQueue-test | & more
//

struct thing_s {
  int    a;
  int    b;
  double c;
  int    d;
};


int
sortthing(const void *a, const void *b) {
  thing_s *A = *((thing_s **)a);
  thing_s *B = *((thing_s **)b);

  if (A->a < B->a)
    return(-1);
  if (A->a > B->a)
    return(1);
  if (A->b < B->b)
    return(-1);
  if (A->b > B->b)
    return(1);
  return(0);
}


int
main(int argc, char **argv) {
  bigQueue *T = new bigQueue(sortthing, 0L, 0L, 0L, sizeof(thing_s), 1, 0L);

  mt_s  *mtctx = mtInit(3);

  int   testSize = 2000000;

  FILE *out = fopen("junk-bigQueue-out-1", "w");

  for (int i=0; i<testSize; i++) {
    thing_s *t = new thing_s;
    t->a =  mtRandom32(mtctx) / 4;
    t->b =  i;
    t->c =  (double)i;
    t->d = -i;

    fprintf(out, "%012d %08d %12.3f %08d\n", t->a, t->b, t->c, t->d);

    T->add(t);
  }

  fclose(out);
  out = fopen("junk-bigQueue-out-2", "w");

  T->sort();

  while (T->next()) {
    thing_s *t = (thing_s *)T->get();

    fprintf(out, "%012d %08d %12.3f %08d\n", t->a, t->b, t->c, t->d);
  }

  delete T;

  fclose(out);
}

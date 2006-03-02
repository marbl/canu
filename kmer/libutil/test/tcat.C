#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "sweatShop.H"

//  Reads stdin, writes stdout.  Uses threads.

typedef struct {
  int     frameNumber;
  int     dataLen;
  char   *data;
} tcat_s;

typedef struct {
  int     globalFrameNumber;
} tcatGlobal_s;


void*
tcatReader(void *G) {
  tcatGlobal_s  *g = (tcatGlobal_s *)G;
  tcat_s        *s = new tcat_s;

  //  XXX:  Should use safeRead()

  s->frameNumber = g->globalFrameNumber++;
  s->data        = new char [1024];
  s->dataLen     = fread(s->data, sizeof(char), 1024, stdin);

  if (s->dataLen == 0L) {
    delete [] s->data;
    delete    s;
    return(0L);
  }

  return(s);
}

void
tcatWorker(void *G, void *S) {
  //tcatGlobal_s  *g = (tcatGlobal_s *)G;
  //tcat_s        *s = (tcat_s       *)S;

  //  Noop!

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = lrand48() % 10000000;

  nanosleep(&naptime, 0L);
}

void
tcatWriter(void *G, void *S) {
  //tcatGlobal_s  *g = (tcatGlobal_s *)G;
  tcat_s        *s = (tcat_s       *)S;

  //  XXX:  Should use safeWrite()

  fwrite(s->data, sizeof(char), s->dataLen, stdout);

  delete [] s->data;
  delete    s;
}


int
main(int argc, char **argv) {
  sweatShop   *ss = new sweatShop(tcatReader, tcatWorker, tcatWriter);

  tcatGlobal_s  *G = new tcatGlobal_s;
  G->globalFrameNumber = 0;

  srand48(time(NULL));

  ss->run(G);

  exit(0);
}

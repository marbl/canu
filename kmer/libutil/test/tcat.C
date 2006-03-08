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


#define BLOCKSIZE 8192

void*
tcatReader(void *G) {
  tcatGlobal_s  *g = (tcatGlobal_s *)G;
  tcat_s        *s = new tcat_s;

  s->frameNumber = g->globalFrameNumber++;
  s->data        = new char [BLOCKSIZE];
  s->dataLen     = safeRead(STDIN_FILENO, s->data, "tcatReader", sizeof(char) * BLOCKSIZE);

  if (s->dataLen == 0) {
    delete [] s->data;
    delete    s;
    return(0L);
  }

  return(s);
}

void
tcatWorker(void *G, void *T, void *S) {
  //tcatGlobal_s  *g = (tcatGlobal_s *)G;
  //tcat_s        *s = (tcat_s       *)S;

  //  Noop!
}

void
tcatWriter(void *G, void *S) {
  //tcatGlobal_s  *g = (tcatGlobal_s *)G;
  tcat_s        *s = (tcat_s       *)S;

  safeWrite(STDOUT_FILENO, s->data, "tcatWriter", sizeof(char) * s->dataLen);

  delete [] s->data;
  delete    s;
}


int
main(int argc, char **argv) {
  sweatShop   *ss = new sweatShop(tcatReader, tcatWorker, tcatWriter);

  ss->loaderQueueSize(64 * 1024 * 1024 / BLOCKSIZE);
  ss->writerQueueSize(64 * 1024 * 1024 / BLOCKSIZE);

  tcatGlobal_s  *G = new tcatGlobal_s;
  G->globalFrameNumber = 0;

  ss->run(G);

  exit(0);
}

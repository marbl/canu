#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "sweatShop.H"

//  Reads stdin, writes stdout.  Uses threads.

int  blockSize  = 8192;
int  readBuf    = 64;
int  writBuf    = 64;

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

  s->frameNumber = g->globalFrameNumber++;
  s->data        = new char [blockSize];
  s->dataLen     = safeRead(STDIN_FILENO, s->data, "tcatReader", sizeof(char) * blockSize);

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

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-r") == 0) {
      readBuf = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-w") == 0) {
      writBuf = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-b") == 0) {
      blockSize = atoi(argv[++arg]);
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }

  if (err) {
    fprintf(stderr, "usage: %s [-b blockSizeBytes] [-r readBufferSizeMB] [-w writeBufferSizeMB]\n", argv[0]);
    exit(1);
  }

  ss->loaderQueueSize(readBuf * 1024 * 1024 / blockSize);
  ss->writerQueueSize(writBuf * 1024 * 1024 / blockSize);

  tcatGlobal_s  *G = new tcatGlobal_s;
  G->globalFrameNumber = 0;

  ss->run(G);

  exit(0);
}

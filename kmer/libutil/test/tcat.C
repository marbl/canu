#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "sweatShop.H"

//  Reads stdin, writes stdout.  Uses threads.

int  blockSize  = 8192;

struct tcat_s {
  int     dataLen;
  char   *data;
};

void*
tcatReader(void *) {
  tcat_s        *s = new tcat_s;

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
tcatWorker(void *, void *, void *) {
  //  Noop!
}

void
tcatWriter(void *, void *S) {
  tcat_s        *s = (tcat_s *)S;

  safeWrite(STDOUT_FILENO, s->data, "tcatWriter", sizeof(char) * s->dataLen);

  delete [] s->data;
  delete    s;
}


int
main(int argc, char **argv) {
  int  readBuf    = 64;
  int  writBuf    = 64;

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

  sweatShop   *ss = new sweatShop(tcatReader, tcatWorker, tcatWriter);

  ss->setLoaderQueueSize(readBuf * 1024 * 1024 / blockSize);
  ss->setWriterQueueSize(writBuf * 1024 * 1024 / blockSize);

  ss->run();

  exit(0);
}

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

int
main(int argc, char **argv) {

  if (argc == 1) {
    fprintf(stderr, "usage: %s <count-1> <count-2> <....>\n", argv[0]);
    fprintf(stderr, "       This is part of the ESTmapper; you don't want to run it by hand.\n");
    exit(1);
  }

  int numFiles = argc-1;

  FILE **Fs = new FILE * [numFiles];
  for (int i=1; i<argc; i++) {
    errno = 0;
    Fs[i-1] = fopen(argv[i], "r");
    if (errno) {
      fprintf(stderr, "%s: ERROR: couldn't open %s: %s\n", argv[0], argv[i], strerror(errno));
      exit(1);
    }
  }

  char  buf[256];
  int   eof = 0;

  while (eof == 0) {
    int count = 0;
    for (int i=0; i<numFiles; i++) {
      fgets(buf, 256, Fs[i]);
      if (feof(Fs[i])) {
        eof++;
      } else {
        count += atoi(buf);
      }
    }
    if (eof == 0)
      fprintf(stdout, "%d\n", count);
  }

  if (eof != numFiles) {
    fprintf(stderr, "%s: ERROR: Short read on a count file.\n", argv[0]);
    exit(1);
  }

  //  Let the OS cleanup for us...

  exit(0);
}

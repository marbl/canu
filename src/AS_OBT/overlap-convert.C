#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>

#include "util++.H"
#include "overlap.H"
#include "constants.H"

//  Reads parital overlaps, rewrites as either binary or ascii.

int
main(int argc, char **argv) {
  char  *inPath = 0L;
  char  *otPath = 0L;
  FILE  *inFile = 0L;
  FILE  *otFile = 0L;
  bool   isBinary = false;
  bool   isASCII  = false;

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-b") == 0) {
      inPath = argv[++arg];
      otPath = argv[++arg];
      inFile = fopen(inPath, "r");
      otFile = fopen(otPath, "wb");
      isBinary = true;
    } else if (strcmp(argv[arg], "-a") == 0) {
      inPath = argv[++arg];
      inFile = fopen(inPath, "r");
      isASCII = true;
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (inPath == 0L) {
    fprintf(stderr, "usage: %s -b <input.ovl> <output.ovb>\n", argv[0]);
    fprintf(stderr, "       %s -a <input.ovb>\n", argv[0]);
    fprintf(stderr, " -b -- convert from ascii to binary\n");
    fprintf(stderr, " -a -- convert from binary to ascii, writing to stdout\n");
    exit(1);
  }

  if (inFile == 0L)
    fprintf(stderr, "Can't open '%s' for input: %s\n", inPath, strerror(errno));
  if ((otFile == 0L) && isBinary)
    fprintf(stderr, "Can't open '%s' for output: %s\n", otPath, strerror(errno));

  overlap_t    overlap;

#ifdef SPEEDCOUNTER_H
  speedCounter  *C = new speedCounter("%7.2f Moverlaps -- %5.2f Moverlaps/second\r",
                                      1000000.0, 0x7fff, true);
  C->enableLiner();
#endif


  if (isBinary) {
    char line[1024];

    fgets(line, 1024, inFile);
    while (!feof(inFile)) {
      overlap.decode(line, false);
      if (overlap.acceptable())
        overlap.dump(otFile);
      fgets(line, 1024, inFile);
#ifdef SPEEDCOUNTER_H
      C->tick();
#endif
    }

    fclose(inFile);
    fclose(otFile);
  }


  if (isASCII) {
    overlap.load(inFile);
    while (!feof(inFile)) {
      if (overlap.acceptable())
        overlap.print(stdout);
      overlap.load(inFile);
#ifdef SPEEDCOUNTER_H
      C->tick();
#endif
    }

    fclose(inFile);
  }

#ifdef SPEEDCOUNTER_H
  delete C;
#endif

  exit(0);
}


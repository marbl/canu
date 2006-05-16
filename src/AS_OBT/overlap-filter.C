#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>

#include "util++.H"
#include "overlap.H"
#include "constants.H"

//  Reads parital overlaps (as either binary or ascii), applies a
//  simple filter, and writes (either binary or ascii).

int
main(int argc, char **argv) {
  char  *inPath = 0L;
  FILE  *inFile = 0L;
  bool   inBinary = false;
  bool   inASCII  = false;

  char  *otPath = 0L;
  FILE  *otFile = 0L;
  bool   otBinary = false;
  bool   otASCII  = false;

  int    minLength = 0;
  int    maxError  = (int)(100.0 * 100);  //  100% error

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-ia") == 0) {
      inPath = argv[++arg];
      inFile = fopen(inPath, "r");
      inASCII = true;
    } else if (strcmp(argv[arg], "-ib") == 0) {
      inPath = argv[++arg];
      inFile = fopen(inPath, "r");
      inBinary = true;

    } else if (strcmp(argv[arg], "-oa") == 0) {
      otPath = argv[++arg];
      otFile = fopen(otPath, "w");
      otASCII = true;
    } else if (strcmp(argv[arg], "-ob") == 0) {
      otPath = argv[++arg];
      otFile = fopen(otPath, "w");
      otBinary = true;

    } else if (strcmp(argv[arg], "-l") == 0) {
      minLength = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-e") == 0) {
      maxError = (int)(atof(argv[++arg]) * 100);

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if ((inPath == 0L) ||
      (otPath == 0L) ||
      ((minLength == 0) && (maxError == 100))) {
    fprintf(stderr, "usage: %s [-ia | -ib infile] [-oa | -ob outfile] [-l minlen] [-e maxerror (float, percent)]\n", argv[0]);
    exit(1);
  }

  if (inFile == 0L)
    fprintf(stderr, "Can't open '%s' for input: %s\n", inPath, strerror(errno));
  if (otFile == 0L)
    fprintf(stderr, "Can't open '%s' for output: %s\n", otPath, strerror(errno));

  fprintf(stderr, "filtering out overlaps shorter than %d and with more than %f%% error.\n", minLength, maxError / 100.0);

  char         line[1024];
  overlap_t    overlap;

#ifdef SPEEDCOUNTER_H
  speedCounter  *C = new speedCounter("%7.2f Moverlaps -- %5.2f Moverlaps/second\r",
                                      1000000.0, 0x7ffff, true);
  C->enableLiner();
#endif

  while (!feof(inFile)) {
    if (inASCII) {
      fgets(line, 1024, inFile);
      overlap.decode(line, false);
    }

    if (inBinary) {
      overlap.load(inFile);
    }

    if (!feof(inFile)) {

      if ((overlap.erate <= maxError) &&
          ((overlap.Aend - overlap.Abeg) >= minLength)) {

        if (otASCII)
          overlap.print(otFile);
        if (otBinary)
          overlap.dump(otFile);
      }
    }

#ifdef SPEEDCOUNTER_H
    C->tick();
#endif
  }

  fclose(inFile);
  fclose(otFile);

#ifdef SPEEDCOUNTER_H
  delete C;
#endif

  exit(0);
}

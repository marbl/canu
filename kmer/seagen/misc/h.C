#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"

//  Generates a histogram of a hit file

typedef struct {
  u32bit dir;
  u32bit estID;
  u32bit scfID;
  u32bit scfLo;
  u32bit scfHi;
} hit_s;


int
main(int argc, char **argv) {

  if (argc < 2) {
    fprintf(stderr, "usage: %s ....\n", argv[0]);
    exit(1);
  }

  int arg = 1;

  //  Things for reading hits
  //
  FILE    *file;
  char     b[1024];
  aHit     a;
  bool     isBINARY;
  int      histogram[10] = {0};

  while (arg < argc) {

    //  Open the file, fatally failing if we cannot do it.
    //
    errno = 0;
    file = fopen(argv[arg], "r");
    if (file == 0L) {
      fprintf(stderr, "ESTmapper/filterEST-- ERROR opening '%s'\n%s\n", argv[arg], strerror(errno));
      exit(1);
    }

    //  Binary or ASCII input?
    //
    char x = (char)fgetc(file);
    ungetc(x, file);

    isBINARY = (x != '-');

    if (isBINARY)
      fprintf(stderr, "reading BINARY hits from '%s'\n", argv[arg]);
    else
      fprintf(stderr, "reading ASCII hits from '%s'\n", argv[arg]);

    //  Read hits until we run out of space
    //
    while (!feof(file)) {
      if (isBINARY) {
        ahit_readBinary(&a, file);
      } else {
        fgets(b, 1024, file);
        ahit_parseString(&a, b);
      }

      //  Fill the histogram
      //
      int len = a._dsHi - a._dsLo;
      if      (len < 25000)
        histogram[0]++;
      else if (len < 50000)
        histogram[1]++;
      else if (len < 100000)
        histogram[2]++;
      else if (len < 200000)
        histogram[3]++;
      else if (len < 400000)
        histogram[4]++;
      else if (len < 800000)
        histogram[5]++;
      else if (len < 1600000)
        histogram[6]++;
      else if (len < 3200000)
        histogram[7]++;
      else if (len < 6400000)
        histogram[8]++;
      else
        histogram[9]++;
    }

    fclose(file);

    arg++;
  }

  for (int i=0; i<10; i++)
    fprintf(stderr, "%2d] %d\n", i, histogram[i]);

  return(0);
}

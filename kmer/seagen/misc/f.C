#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"

//  Filters a hit file based on the length of the genomic region

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
  int      histogram[12] = {0};
  FILE    *outf[12];

  outf[0] = fopen("filteredHits.0", "w");
  outf[1] = fopen("filteredHits.1", "w");
  outf[2] = fopen("filteredHits.2", "w");
  outf[3] = fopen("filteredHits.3", "w");
  outf[4] = fopen("filteredHits.4", "w");
  outf[5] = fopen("filteredHits.5", "w");
  outf[6] = fopen("filteredHits.6", "w");
  outf[7] = fopen("filteredHits.7", "w");
  outf[8] = fopen("filteredHits.8", "w");
  outf[9] = fopen("filteredHits.9", "w");
  outf[10] = fopen("filteredHits.a", "w");
  outf[11] = fopen("filteredHits.b", "w");

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
      if      (len < 25000) {
        fprintf(outf[0], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[0]++;
      } else if (len < 50000) {
        fprintf(outf[1], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[1]++;
      } else if (len < 100000) {
        fprintf(outf[2], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[2]++;
      } else if (len < 200000) {
        fprintf(outf[3], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[3]++;
      } else if (len < 400000) {
        fprintf(outf[4], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[4]++;
      } else if (len < 800000) {
        fprintf(outf[5], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[5]++;
      } else if (len < 1600000) {
        fprintf(outf[6], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[6]++;
      } else if (len < 3200000) {
        fprintf(outf[7], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[7]++;
      } else if (len < 6400000) {
        fprintf(outf[8], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[8]++;
      } else if (len < 12800000) {
        fprintf(outf[9], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[9]++;
      } else if (len < 25600000) {
        fprintf(outf[10], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[10]++;
      } else {
        fprintf(outf[11], "-%c -e %u -D %u %u %u\n", a._direction ? 'f' : 'r', a._qsIdx, a._dsIdx, a._dsLo, a._dsHi);
        histogram[11]++;
      }
    }

    fclose(file);

    arg++;
  }

  for (int i=0; i<12; i++)
    fprintf(stderr, "%2d] %d\n", i, histogram[i]);

  return(0);
}

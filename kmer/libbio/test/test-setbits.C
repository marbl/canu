#include "bio++.H"

//g++ -o test-setbits test-setbits.C -I../libutil -I. -L../libutil -L. -lbio -lutil

int
main(int argc, char **argv) {
  
  kMer   x(96);
  char   str[256];

  if (KMER_WORDS < 3) {
    fprintf(stderr, "I need at least KMER_WORDS == 3; test not run.\n");
    exit(0);
  }

  for (u32bit i=0; i<168; i++) {
    x.clear();
    x.setBits(i, 24, 0x535);
    fprintf(stderr, u32bitFMTW(3)" -- %s -- "u64bitHEX"\n", i, x.merToString(str), x.getBits(i, 16));

    if (x.getBits(i, 16) != 0x535) {
      fprintf(stderr, "decode error.\n");
      exit(1);
    }
  }

  exit(0);
}

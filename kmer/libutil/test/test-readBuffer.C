#include <stdio.h>

#include "util++.H"

char            *filename = 0L;
md5_s           *full     = 0L;
md5_s           *part     = 0L;


int
doTest(readBuffer *B, md5_s *correct, const char *description) {
  int              error    = 0;
  md5_increment_s *testing  = 0L;
  int              bread    = 0;

  fprintf(stderr, "readBuffer test %s.\n", description);

  for (char x = B->read(); !B->eof(); x = B->read()) {
    testing = md5_increment_char(testing, x);
    bread++;
  }

  md5_increment_finalize(testing);

  if ((testing->a != correct->a) || (testing->b != correct->b)) {
    fprintf(stderr, "readBuffer test %s failed (read %d bytes).\n", description, bread);
    fprintf(stderr, "Got correct md5 of "u64bitHEX" "u64bitHEX"\n", correct->a, correct->b);
    fprintf(stderr, "Got testing md5 of "u64bitHEX" "u64bitHEX"\n", testing->a, testing->b);
    error = 1;
  }

  md5_increment_destroy(testing);

  return(error);
}


int
doTestRead(readBuffer *B, md5_s *correct, size_t bufferSize, const char *description) {
  int      error         = 0;
  char    *buffer        = new char [bufferSize];
  size_t   bufferLen     = 0;

  md5_increment_s *testing  = 0L;

  fprintf(stderr, "readBuffer test %s.\n", description);

  while (!B->eof()) {
    bufferLen = B->read(buffer, bufferSize);
    //fprintf(stderr, "Read bufferLen=%d bufferSize=%d\n", bufferLen, bufferSize);
    testing = md5_increment_block(testing, buffer, bufferLen);
  }

  md5_increment_finalize(testing);

  if ((testing->a != correct->a) || (testing->b != correct->b)) {
    fprintf(stderr, "readBuffer test %s failed.\n", description);
    fprintf(stderr, "Got correct md5 of "u64bitHEX" "u64bitHEX"\n", correct->a, correct->b);
    fprintf(stderr, "Got testing md5 of "u64bitHEX" "u64bitHEX"\n", testing->a, testing->b);
    error = 1;
  }

  md5_increment_destroy(testing);

  return(error);
}


int
main(int argc, char **argv) {
  int         error = 0;
  readBuffer *B = 0L;

  size_t      L = 0;
  size_t      H = 0;
  size_t      R = 0;

  //  If we are given a file, use that, otherwise, use ourself.
  //
  filename = argv[argc-1];

  L = sizeOfFile(filename);
  H = L/2;
  R = L - H;

  fprintf(stderr, "L=%d H=%d R=%d\n", L, H, R);

  //  Suck in the whole file, compute the correct md5 checksum on it
  //
  char *c = new char [L];

  FILE *F = fopen(filename, "r");
  fread(c, sizeof(char), L, F);
  fclose(F);
  full = md5_string(0L, c,   L);
  part = md5_string(0L, c+H, R);

  delete [] c;


  B = new readBuffer(filename, 999);
  error += doTest(B, full, "#1 (read)");
  B->seek(0);
  error += doTest(B, full, "#2 (seek)");
  B->seek(H);
  error += doTest(B, part, "#2 (seek half)");
  delete B;

  B = new readBuffer(filename, 0);
  error += doTest(B, full, "#3 (mmap)");
  B->seek(0);
  error += doTest(B, full, "#2 (mmap seek)");
  B->seek(H);
  error += doTest(B, part, "#2 (mmap seek half)");
  delete B;

  B = new readBuffer(filename, 0);
  error += doTestRead(B, full, 10000, "#4 (read buffer=mmap readsize=10000)");
  delete B;

  B = new readBuffer(filename, 100);
  error += doTestRead(B, full, 10000, "#4 (read buffer=100 readsize=10000)");
  delete B;

  B = new readBuffer(filename, 2000);
  error += doTestRead(B, full, 1000, "#4 (read buffer=2000 readsize=1000)");
  delete B;

  B = new readBuffer(filename, L);
  error += doTestRead(B, full, L+1000, "#5 (read buffer=filesize readsize=filesize+1000)");
  delete B;

  return(error);
}


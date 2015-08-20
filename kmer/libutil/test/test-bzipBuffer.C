
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2006-JUN-23 to 2014-APR-11
 *      are Copyright 2006,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>

#include "util++.H"

char            *filename = 0L;
md5_s           *correct  = 0L;

int
doTest(bzipBuffer *B, char *description) {
  int              error    = 0;
  md5_increment_s *testing  = 0L;

  while (!B->eof())
    testing = md5_increment_char(testing, B->getnext());

  md5_increment_finalize(testing);

  if ((testing->a != correct->a) || (testing->b != correct->b)) {
    fprintf(stderr, "bzipBuffer test %s failed.\n", description);
    fprintf(stderr, "Got correct md5 of "uint64HEX" "uint64HEX"\n", correct->a, correct->b);
    fprintf(stderr, "Got testing md5 of "uint64HEX" "uint64HEX"\n", testing->a, testing->b);
    error = 1;
  }

  md5_increment_destroy(testing);

  return(error);
}

int
doTestRead(bzipBuffer *B, size_t bufferSize, char *description) {
  int      error         = 0;
  char    *buffer        = new char [bufferSize];
  size_t   bufferLen     = 0;

  md5_increment_s *testing  = 0L;

  while (!B->eof()) {
    bufferLen = B->read(buffer, bufferSize);
    testing = md5_increment_block(testing, buffer, bufferLen);
  }

  md5_increment_finalize(testing);

  if ((testing->a != correct->a) || (testing->b != correct->b)) {
    fprintf(stderr, "bzipBuffer test %s failed.\n", description);
    fprintf(stderr, "Got correct md5 of "uint64HEX" "uint64HEX"\n", correct->a, correct->b);
    fprintf(stderr, "Got testing md5 of "uint64HEX" "uint64HEX"\n", testing->a, testing->b);
    error = 1;
  }

  md5_increment_destroy(testing);

  return(error);
}


int
main(int argc, char **argv) {
  int         error = 0;
  bzipBuffer *B = 0L;

  //  If we are given a file, use that, otherwise, use ourself.
  //
  filename = argv[argc-1];


  //  Suck in the whole file, compute the correct md5 checksum on it
  //
  char *c = new char [sizeOfFile(filename)];
  FILE *F = fopen(filename, "r");
  fread(c, sizeof(char), sizeOfFile(filename), F);
  fclose(F);
  correct = md5_string(0L, c, sizeOfFile(filename));
  delete [] c;


  //  Test just reading, with a small buffer
  //
  B = new bzipBuffer(filename, 999);
  error += doTest(B, "#1 (read)");


  exit(1);


  //  Test read() with a small buffer, reading large chunks
  //
  B = new bzipBuffer(filename, 100);
  error += doTestRead(B, 10000, "#4 (read)");
  delete B;


  //  Test read() with a small buffer, reading small chunks that are a
  //  factor of the buffersize.
  //
  B = new bzipBuffer(filename, 2000);
  error += doTestRead(B, 1000, "#4 (read)");
  delete B;


  //  Test read() with a large buffer, reading even larger pieces
  //
  B = new bzipBuffer(filename, sizeOfFile(filename));
  error += doTestRead(B, sizeOfFile(filename) + 100000, "#5 (read)");
  delete B;

  return(error);
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "meryl.H"

const char *err1 = "ERROR: %s._merSizeInBases (%u) != %s._merSizeInBases (%u)!\n";
const char *err2 = "ERROR: %s._tableSizeInBits (%u) != %s._tableSizeInBits (%u)!\n";
const char *err3 = "ERROR: %s._chckBits (%u) != %s._chckBits (%u)!\n";
const char *err4 = "ERROR: %s._hashWidth (%u) != %s._hashWidth (%u)!\n";


bool
checkSingleDescription(mcDescription *a, char *A,
                       mcDescription *b, char *B) {
  bool fail = false;

  if (a->_merSizeInBases != b->_merSizeInBases) {
    fprintf(stderr, err1, A, a->_merSizeInBases, B, b->_merSizeInBases);
    fail = true;
  }

  if (a->_tableSizeInBits != b->_tableSizeInBits) {
    fprintf(stderr, err2, A, a->_tableSizeInBits, B, b->_tableSizeInBits);
    fail = true;
  }

  if (a->_chckBits != b->_chckBits) {
    fprintf(stderr, err3, A, a->_chckBits, B, b->_chckBits);
    fail = true;
  }

  if (a->_hashWidth != b->_hashWidth) {
    fprintf(stderr, err4, A, a->_hashWidth, B, b->_hashWidth);
    fail = true;
  }

  return(fail);
}

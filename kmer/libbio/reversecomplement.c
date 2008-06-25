#include "bio.h"

#include <string.h>

//  Inplace reverse-complement an ACGT sequence.  A pointer the the
//  string is returned.
//
char *
reverseComplementSequence(char *seq, u32bit seqlen) {
  char   *s = seq;
  char   *e = seq + seqlen - 1;
  char    t;
  u32bit  c = seqlen / 2;

  while (c--) {
    t = complementSymbol[*s];
    *(s++) = complementSymbol[*e];
    *(e--) = t;
  }

  if (s == e)
    *s = complementSymbol[*s];

  return(seq);
}


//  Inplace reverse a string.  A pointer the the string is returned.
//
char *
reverseString(char *seq, u32bit seqlen) {
  char   *s = seq;
  char   *e = seq + seqlen - 1;
  char    t;
  u32bit  c = seqlen / 2;

  while (c--) {
    t = *s;
    *(s++) = *e;
    *(e--) = t;
  }

  return(seq);
}

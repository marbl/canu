#include "bri.h"

//  Inplace reverse-complement the sequence.  A pointer
//  the the character string is returned.
//
char *
reverseComplementSequence(char *seq, u32bit seqlen) {
  char   *s = seq;
  char   *e = seq + seqlen - 1;
  char    t;
  u32bit  c = seqlen / 2;

  while (c--) {
    t = complementSymbol[(int)*s];
    *(s++) = complementSymbol[(int)*e];
    *(e--) = t;
  }

  if (s == e)
    *s = complementSymbol[(int)*s];

  return(seq);
}

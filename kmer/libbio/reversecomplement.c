#include "libbritypes.h"
#include "alphabet.h"

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
    t = complementSymbol[*s];
    *(s++) = complementSymbol[*e];
    *(e--) = t;
  }

  if (s == e)
    *s = complementSymbol[*s];

  return(seq);
}

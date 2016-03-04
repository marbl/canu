
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
 *    Brian P. Walenz from 2008-MAY-15 to 2013-AUG-01
 *      are Copyright 2008-2009,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

static char inv[256] = {0};


static
void
initRC(void) {
  if (inv['a'] == 't')
    return;

  inv['a'] = 't';
  inv['c'] = 'g';
  inv['g'] = 'c';
  inv['t'] = 'a';
  inv['n'] = 'n';
  inv['A'] = 'T';
  inv['C'] = 'G';
  inv['G'] = 'C';
  inv['T'] = 'A';
  inv['N'] = 'N';
  inv['-'] = '-';
}


void
reverseComplementSequence(char *seq, int len) {
  char   c=0;
  char  *s=seq,  *S=seq+len-1;

  initRC();

  if (len == 0) {
    len = strlen(seq);
    S = seq + len - 1;
  }

  while (s < S) {
    c    = *s;
    *s++ =  inv[*S];
    *S-- =  inv[c];
  }

  if (s == S)
    *s = inv[*s];
}

//  Inplace reverse-complement an ACGT sequence.  A pointer the the
//  string is returned.
//
#if 0
//  From kmer
char *
reverseComplementSequence(char *seq, uint32 seqlen) {
  char   *s = seq;
  char   *e = seq + seqlen - 1;
  char    t;
  uint32  c = seqlen / 2;

  while (c--) {
    t = complementSymbol[*s];
    *(s++) = complementSymbol[*e];
    *(e--) = t;
  }

  if (s == e)
    *s = complementSymbol[*s];

  return(seq);
}
#endif

void
reverseComplement(char *seq, char *qlt, int len) {
  char   c=0;
  char  *s=seq,  *S=seq+len-1;
  char  *q=qlt,  *Q=qlt+len-1;

  if (qlt == NULL) {
    reverseComplementSequence(seq, len);
    return;
  }

  initRC();

  if (len == 0) {
    len = strlen(seq);
    S = seq + len - 1;
    Q = qlt + len - 1;
  }

  while (s < S) {
    c    = *s;
    *s++ =  inv[*S];
    *S-- =  inv[c];

    c    = *q;
    *q++ = *Q;
    *Q-- =  c;
  }

  if (s == S)
    *s = inv[*s];
}


void
reverse(char *a, char *b, int len) {
  char   c=0;
  char  *s=a,  *S=a+len-1;
  char  *q=b,  *Q=b+len-1;

  while (s < S) {
    c    = *s;
    *s++ =  *S;
    *S-- =  c;

    c    = *q;
    *q++ = *Q;
    *Q-- =  c;
  }
}


//  Inplace reverse a string.  A pointer the the string is returned.
//
#if 0
//  From kmer
char *
reverseString(char *seq, uint32 seqlen) {
  char   *s = seq;
  char   *e = seq + seqlen - 1;
  char    t;
  uint32  c = seqlen / 2;

  while (c--) {
    t = *s;
    *(s++) = *e;
    *(e--) = t;
  }

  return(seq);
}
#endif

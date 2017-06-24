
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


static
char
inv[256] = {
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x00 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x08 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x10 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x18 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x20 -  !"#$%&'
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x28 - ()*+,-./
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x30 - 01234567
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x38 - 89:;<=>?
   0,'T',  0,'G',  0,  0,  0,'C',  //  0x40 - @ABCDEFG
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x48 - HIJKLMNO
   0,  0,  0,  0,'A',  0,  0,  0,  //  0x50 - PQRSTUVW
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x58 - XYZ[\]^_
   0,'t',  0,'g',  0,  0,  0,'c',  //  0x60 - `abcdefg
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x68 - hijklmno
   0,  0,  0,  0,'a',  0,  0,  0,  //  0x70 - pqrstuvw
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x78 - xyz{|}~
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x80 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x88 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x90 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x98 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xa0 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xa8 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xb0 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xb8 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xc0 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xc8 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xd0 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xd8 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xe0 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xe8 - 
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xf0 - 
   0,  0,  0,  0,  0,  0,  0,  0   //  0xf8 - 
};



void
reverseComplementSequence(char *seq, int len) {
  char   c=0;
  char  *s=seq,  *S=seq+len-1;

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



char *
reverseComplementCopy(char *seq, int len) {
  char  *rev = new char [len+1];

  assert(len > 0);

  for (int32 p=len, q=0; p>0; )
    rev[q++] = inv[seq[--p]];

  rev[len] = 0;

  return(rev);
}



void
reverseComplement(char *seq, char *qlt, int len) {
  char   c=0;
  char  *s=seq,  *S=seq+len-1;
  char  *q=qlt,  *Q=qlt+len-1;

  if (qlt == NULL) {
    reverseComplementSequence(seq, len);
    return;
  }

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


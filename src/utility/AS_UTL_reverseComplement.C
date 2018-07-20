
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
   0,  0,  0,  0,  0,  0,  'N',  0,  //  0x48 - HIJKLMNO
   0,  0,  0,  0,'A',  0,  0,  0,  //  0x50 - PQRSTUVW
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x58 - XYZ[\]^_
   0,'t',  0,'g',  0,  0,  0,'c',  //  0x60 - `abcdefg
   0,  0,  0,  0,  0,  0,  'n',  0,  //  0x68 - hijklmno
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



template<typename qvType>
void
reverseComplement(char *seq, qvType *qlt, int len) {
  char    c=0;
  char   *s=seq,  *S=seq+len-1;
  qvType *q=qlt,  *Q=qlt+len-1;

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

template void reverseComplement<char> (char *seq, char  *qlt, int len);   //  Give the linker
template void reverseComplement<uint8>(char *seq, uint8 *qlt, int len);   //  something to link

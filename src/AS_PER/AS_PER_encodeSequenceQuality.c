
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static char *rcsid = "$Id: AS_PER_encodeSequenceQuality.c,v 1.9 2008-10-09 00:48:12 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "AS_PER_encodeSequenceQuality.h"


//  This module encodes/decodes sequence/quality strings into a string
//  of one char per seq/quality value pair
//
// Quality values are encoded as follows:
//
// lowest 2 bits for sequence value:
//    a = 0
//    c = 1
//    t = 2
//    g = 3
//
// highest 6 bits for quality value:
//    0 is unused (qv=0 + 'a' == 0, which terminates the string)
//    1 - 61 are valid quality values
//    62 is unused
//    63 encodes an n sequence value
//
//
// If no quality values are supplied (the encodeSequence() and
// decodeSequence() interfaces) sequence is 2-bit encoded if only ACGT
// is present, otherwise, 4-bit encoded (ACGTN; four bits for
// simplicity).

#define SEQ_A 0x00
#define SEQ_C 0x01
#define SEQ_G 0x02
#define SEQ_T 0x03
#define SEQ_N 0xff

#define QUALITY_MAX  60

void
encodeSequenceQuality(char *enc,
                      char *seq,
                      char *qlt) {

  while ((*seq != 0) && (*qlt != 0)) {
    unsigned char qv;
    unsigned char sv;

    switch (*seq) {
      case 'a':
      case 'A':
        sv = SEQ_A;
        break;
      case 'c':
      case 'C':
        sv = SEQ_C;
        break;
      case 'g':
      case 'G':
        sv = SEQ_G;
        break;
      case 't':
      case 'T':
        sv = SEQ_T;
        break;
      case 'n':
      case 'N':
        sv = SEQ_N;
        break;
      default:
        fprintf(stderr,"encodeSequenceQuality()-- Illegal char %c detected!  Kaboom!\n", *seq);
        assert(0);
        break;
    }

    assert(*qlt >= '0');
    assert(*qlt <= QUALITY_MAX + '0');

    qv   = *qlt - '0' + 1;

    *enc = (qv << 2) | sv;

    seq++;
    qlt++;
    enc++;
  }

  *enc = 0;

  assert(*seq == 0);
  assert(*qlt == 0);
}

void
decodeSequenceQuality(char *enc,
                      char *seq,
                      char *qlt) {
  const char sm[5] = {'A', 'C', 'G', 'T', 'N'};

  while (*enc) {
    *seq = sm[*enc & 0x03];
    *qlt = ((*enc >> 2) & 0x3f) - 1;

    if (*qlt > QUALITY_MAX) {
      *seq = 'N';
      *qlt = 0;
    }

    *qlt += '0';

    seq++;
    qlt++;
    enc++;
  }
  *seq = 0;
  *qlt = 0;
}





int
encodeSequence(char *enc,
               char *seq) {
  int   len    = strlen(seq);

  int   encLen = 0;

  char  eee    = 0;
  int   eel    = 0;

  int   twob   = 0;
  int   four   = 0;

  int   i;

  for (i=0; i<len; i++) {
    switch (seq[i]) {
      case 'a':
      case 'A':
      case 'c':
      case 'C':
      case 'g':
      case 'G':
      case 't':
      case 'T':
        twob++;
        break;
      default:
        four++;
        break;
    }
  }

  if (four == 0) {
    enc[encLen++] = 't';

    for (i=0; i<len; i++) {
      if (eel == 4) {
        enc[encLen++] = eee;
        eel = 0;
      }

      switch (seq[i]) {
        case 'a':
        case 'A':
          eee <<= 2;
          eee  |= 0x00;  //  %00000000
          eel++;
          break;
        case 'c':
        case 'C':
          eee <<= 2;
          eee  |= 0x01;  //  %00000001
          eel++;
          break;
        case 'g':
        case 'G':
          eee <<= 2;
          eee  |= 0x02;  //  %00000010
          eel++;
          break;
        case 't':
        case 'T':
          eee <<= 2;
          eee  |= 0x03;  //  %00000011
          eel++;
          break;
      }
    }

    eee <<= 2 * (4 - eel);

  } else {
    enc[encLen++] = 'f';

    for (i=0; i<len; i++) {
      if (eel == 2) {
        enc[encLen++] = eee;
        eel = 0;
      }

      switch (seq[i]) {
        case 'a':
        case 'A':
          eee <<= 4;
          eee  |= 0x00;  //  %00000001
          eel++;
          break;
        case 'c':
        case 'C':
          eee <<= 4;
          eee  |= 0x01;  //  %00000001
          eel++;
          break;
        case 'g':
        case 'G':
          eee <<= 4;
          eee  |= 0x02;  //  %00000001
          eel++;
          break;
        case 't':
        case 'T':
          eee <<= 4;
          eee  |= 0x03;  //  %00000001
          eel++;
          break;
        default:
          eee <<= 4;
          eee  |= 0x0f;
          eel++;
          break;
      }
    }

    eee <<= 4 * (2 - eel);
  }

  enc[encLen++] = eee;
  enc[encLen]   = 0;

  return(encLen);
}


void
decodeSequence(char *enc,
               char *seq,
               int   seqLen) {

  int   wrdpos = 0;
  int   encpos = 1;
  char  wrd    = enc[encpos++];

  int   len;

  if        (enc[0] == 't') {
    //  Sequence is two-bit encoded

    for (len=0; len < seqLen; len++) {
      switch (wrd & 0xc0) {
        case 0x00:  //  %00000000
          seq[len] = 'A';
          break;
        case 0x40:  //  %01000000
          seq[len] = 'C';
          break;
        case 0x80:  //  %10000000
          seq[len] = 'G';
          break;
        case 0xc0:  //  %11000000
          seq[len] = 'T';
          break;
      }

      wrd <<= 2;
      wrdpos++;

      if (wrdpos == 4) {
        wrdpos = 0;
        wrd = enc[encpos++];
      }
    }

  } else if (enc[0] == 'f') {
    //  Sequence is four-bit encoded

    for (len=0; len < seqLen; len++) {
      switch (wrd & 0xf0) {
        case 0x00:  //  %00000000
          seq[len] = 'A';
          break;
        case 0x10:  //  %00010000
          seq[len] = 'C';
          break;
        case 0x20:  //  %00100000
          seq[len] = 'G';
          break;
        case 0x30:  //  %00110000
          seq[len] = 'T';
          break;
        default:
          seq[len] = 'N';
          break;
      }

      wrd <<= 4;
      wrdpos++;

      if (wrdpos == 2) {
        wrdpos = 0;
        wrd = enc[encpos++];
      }
    }

  } else {
    //  Sequence is not encoded.

    strncpy(seq, enc, seqLen);
  }

  seq[seqLen] = 0;
}

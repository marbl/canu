
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

static char CM_ID[] = "$Id: AS_PER_encodeSequenceQuality.c,v 1.6 2007-02-27 07:11:03 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
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

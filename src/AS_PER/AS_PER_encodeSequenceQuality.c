
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
static char CM_ID[] = "$Id: AS_PER_encodeSequenceQuality.c,v 1.2 2004-09-23 20:25:26 mcschatz Exp $";
/*************************************************************************
 Module:  AS_PER_encodeSequenceQuality
 Description:
     This module encodes/decodes sequence/quality strings into a string of one char per seq/quality value pair

 Assumptions:

 Document:

 *************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_PER_encodeSequenceQuality.h"



/*************************************************************************************************/
/* Quality values are encoded as follows:
   2 bits for sequence value:
      a = 0
      c = 1
      t = 2
      g = 3
   6 bits for quality value:
      0 - 60 are valid quality values
      63 encodes an n sequence value
      62 encodes a 99 quality value
*/

#define SEQUENCE_MASK 0xC0
#define QUALITY_MASK  (~ SEQUENCE_MASK)

#define GET_SEQUENCE(ch) (((ch) & SEQUENCE_MASK))
#define GET_QUALITY(ch) ((ch) & QUALITY_MASK)

#define SEQ_A 00
#define SEQ_C 64
#define SEQ_T 128
#define SEQ_G 192
#define SEQ_N 4
#define QUALITY_MAX 60
#define QUALITY_N 63
#define QUALITY_99 62
#define EOS '\0'        

const char SeqChars[] = {
  'A',
  'C',
  'T',
  'G',
  'N'
};

  /*** NOTE -- encoded value is NOT a null terminated string!!!! */
int encodeSequenceQuality
( char *encoded, char *sequence, char *quality, uint hasQuality){
  char *s = sequence;
  char *q = quality;
  char *e = encoded;
  while(*s != EOS && (!hasQuality || *q != EOS)){
    unsigned char qv = 0;
    unsigned char sv;


    *e = 0;

    if(hasQuality){
      qv = *q - '0';
      assert(qv <= QUALITY_MAX);
      q++;
    }
    switch(*s){
    case 'a':
    case 'A':
      sv = SEQ_A;
      break;
    case 't':
    case 'T':
      sv = SEQ_T;
      break;
    case 'c':
    case 'C':
      sv = SEQ_C;
      break;
    case 'g':
    case 'G':
      sv = SEQ_G;
      break;
    case 'n':
    case 'N':
      sv = SEQ_N;
      break;
    default:
      fprintf(stderr,"*** Illegal sequence char %c detected in sequence:%s\n",
	      *s, sequence);
      assert(0);
      break;
    }
    if(sv == SEQ_N){
      sv = SEQ_A;
      qv = QUALITY_N;
    }

    *e = sv | qv;
    
    s++;
    e++;
  }
  /*** NOTE -- encoded value is NOT a null terminated string!!!! */

  assert(*s == *q); // Should both be EOS
  return 0;
}

int decodeSequenceQuality
( char *encoded, int encodedLength, char *sequence, char *quality, 
  uint hasQuality){
  char *s = sequence;
  char *q = quality;
  char *e = encoded;   /*** NOTE -- encoded value is NOT a null terminated string!!!! */

  int i;

  *q = '\0';

  for(i = 0; i < encodedLength; i++){
      *s = SeqChars[ GET_SEQUENCE(*e)/ SEQ_C ]; // 0-3
    if(hasQuality){
      *q = GET_QUALITY((*e));
      if(*q == QUALITY_N){
	*s = SeqChars[SEQ_N];
	*q = 0;
      }
      *q += '0';
      q++;
    }

    s++;
    e++;
  }
  if(hasQuality)
    *q = EOS;

  *s = EOS;

  assert(strlen(sequence) == encodedLength);
#ifdef DEBUG
  fprintf(stderr,"decode seq = %s\nqual = %s\n",
	  sequence, (hasQuality?quality:""));
#endif
  return 0;
}



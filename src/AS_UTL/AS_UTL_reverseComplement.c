
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

#include "AS_global.h"

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


void
reverseComplement(char *seq, char *qlt, int len) {
  char   c=0;
  char  *s=seq,  *S=seq+len-1;
  char  *q=qlt,  *Q=qlt+len-1;

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

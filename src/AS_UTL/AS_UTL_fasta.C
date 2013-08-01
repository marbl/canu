
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

static const char *rcsid = "$Id: AS_UTL_fasta.c,v 1.11 2013-03-20 17:25:27 skoren Exp $";

#include "AS_UTL_fasta.H"
#include "AS_UTL_fileIO.H"

#include <stdarg.h>

int 
AS_UTL_isValidSequence(char *s, int sl) {
  AS_UTL_initValidSequence();
  int p = 0;
   
  for (p = 0; s[p] && p < sl; p++) {
    if ((AS_UTL_isspacearray[s[p]]) || (AS_UTL_isvalidACGTN[s[p]])) {
    } else {
      return FALSE;
    }
  }
    
  return TRUE;
}

void
AS_UTL_writeFastA(FILE *f,
                  char *s, int sl, int bl,
                  char *h, ...) {
  va_list ap;
  char   *o  = (char *)safe_malloc(sizeof(char) * (sl + sl / 60 + 2));
  int     si = 0;
  int     oi = 0;

  while (si < sl) {
    o[oi++] = s[si++];

    if (bl != 0 && (si % bl) == 0)
      o[oi++] = '\n';
  }
  if (o[oi-1] != '\n')
    o[oi++] = '\n';
  o[oi] = 0;

  va_start(ap, h);
  vfprintf(f, h, ap);
  va_end(ap);

  AS_UTL_safeWrite(f, o, "AS_UTL_writeFastA", sizeof(char), oi);

  safe_free(o);
}


void
AS_UTL_writeQVFastA(FILE *f,
                    char *q, int ql, int bl,
                    char *h, ...) {
  va_list ap;
  char   *o  = (char *)safe_malloc(sizeof(char) * (3*ql + 3*ql / 60 + 2));
  int     qi = 0;
  int     oi = 0;

  //
  //  20 values per line -> 60 letters per line.
  //  |xx xx xx xx xx ..... xx|
  //

  while (qi < ql) {
    // decode the quality value
    // we convert the qlt character to the integer value by subtracting '0' and take the significant digit by dividing by ten. Back to character by adding '0'
    o[oi++] = (((((int)q[qi])-'0') / 10) + '0');
    // same thing except now use mod
    o[oi++] = (((((int)q[qi])-'0') % 10) + '0');
    o[oi++] = ' ';

    qi++;

    if (bl != 0 && (qi % bl) == 0)
      o[oi-1] = '\n';
  }
  if (o[oi-1] != '\n')
    o[oi++] = '\n';
  o[oi] = 0;

  va_start(ap, h);
  vfprintf(f, h, ap);
  va_end(ap);

  AS_UTL_safeWrite(f, o, "AS_UTL_writeQVFastA", sizeof(char), oi);

  safe_free(o);
}


void
AS_UTL_writeFastQ(FILE *f,
                  char *s, int sl,
                  char *q, int ql,
                  char *h, ...) {
  va_list ap;
  char   *o  = (char *)safe_malloc(sizeof(char) * (ql + 1));
  int     qi = 0;
  int     oi = 0;

  assert(sl == ql);

  //  Reencode the QV to the Sanger spec.
  while (qi < ql)
    o[oi++] = q[qi++] - '0' + '!';
  o[oi] = 0;

  va_start(ap, h);
  vfprintf(f, h, ap);
  va_end(ap);

  AS_UTL_safeWrite(f, s, "AS_UTL_writeFastQ", sizeof(char), sl);
  fprintf(f, "\n");

  fprintf(f, "+\n");
  AS_UTL_safeWrite(f, o, "AS_UTL_writeFastQ", sizeof(char), ql);
  fprintf(f, "\n");

  safe_free(o);
}

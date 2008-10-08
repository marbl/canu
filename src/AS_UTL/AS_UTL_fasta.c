
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

static const char *rcsid = "$Id: AS_UTL_fasta.c,v 1.5 2008-10-08 22:03:00 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"

void
AS_UTL_writeFastA(FILE *f,
                  char *s, int sl,
                  char *h, ...) {
  va_list ap;
  char   *o  = (char *)safe_malloc(sizeof(char) * (sl + sl / 70 + 2));
  int     si = 0;
  int     oi = 0;

  while (si < sl) {
    o[oi++] = s[si++];

    if ((si % 70) == 0)
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
                    char *q, int ql,
                    char *h, ...) {
  va_list ap;
  char   *o  = (char *)safe_malloc(sizeof(char) * (3*ql + 3*ql / 70 + 2));
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

    if ((qi % 20) == 0)
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


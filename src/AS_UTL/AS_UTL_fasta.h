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

#ifndef AS_UTL_FASTA_H
#define AS_UTL_FASTA_H

static const char *rcsid_AS_UTL_FASTA_H = "$Id: AS_UTL_fasta.h,v 1.6 2009-03-31 20:32:31 skoren Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "AS_global.h"

static int    AS_UTL_isspacearray[256] = {0};
static int    AS_UTL_isvalidACGTN[256] = {0};

static inline void
AS_UTL_initValidSequence(void) {
  if (AS_UTL_isvalidACGTN['a'] == 0) {
    int i;

    for (i=0; i<256; i++)
      AS_UTL_isspacearray[i] = isspace(i);

    AS_UTL_isvalidACGTN['a'] = 'A';
    AS_UTL_isvalidACGTN['c'] = 'C';
    AS_UTL_isvalidACGTN['g'] = 'G';
    AS_UTL_isvalidACGTN['t'] = 'T';
    AS_UTL_isvalidACGTN['n'] = 'N';
    AS_UTL_isvalidACGTN['A'] = 'A';
    AS_UTL_isvalidACGTN['C'] = 'C';
    AS_UTL_isvalidACGTN['G'] = 'G';
    AS_UTL_isvalidACGTN['T'] = 'T';
    AS_UTL_isvalidACGTN['N'] = 'N';
  }
}

static inline int*
AS_UTL_getValidACGTN() {
   AS_UTL_initValidSequence();
   return AS_UTL_isvalidACGTN;
}

static inline int*
AS_UTL_getSpaceArray() {
   AS_UTL_initValidSequence();
   return AS_UTL_isspacearray;
}

int 
AS_UTL_isValidSequence(char *s, int sl);

void 
AS_UTL_writeFastAWithBreaks(FILE *f,
                  char *s, int sl,
                  int lineBreaks,                  
                  char *h, ...);

static inline void
AS_UTL_writeFastA(FILE *f,
                  char *s, int sl,                  
                  char *h, ...) {
   va_list ap;
   va_start(ap, h);
   vfprintf(f, h, ap);
   va_end(ap);
   
   AS_UTL_writeFastAWithBreaks(f, s, sl, 70, "", NULL);
}

//  Assumes a CA encoded QV string, and writes it to a fasta formated
//  as 2-digit white space separated.
//
void
AS_UTL_writeQVFastAWithBreaks(FILE *f,
                    char *q, int ql,
                    int lineBreaks,
                    char *h, ...);

static inline void
AS_UTL_writeQVFastA(FILE *f,
                    char *q, int ql,
                    char *h, ...) {
   va_list ap;
   va_start(ap, h);
   vfprintf(f, h, ap);
   va_end(ap);

   AS_UTL_writeQVFastAWithBreaks(f, q, ql, 20, "", NULL);
                    }

void
AS_UTL_writeQVFastQWithBreaks(FILE *f,
                    char *q, int ql,
                    int lineBreaks, char *h, ...);

#endif

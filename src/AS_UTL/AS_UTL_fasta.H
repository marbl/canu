
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
 *    Brian P. Walenz from 2007-NOV-02 to 2013-AUG-01
 *      are Copyright 2007-2008,2010,2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-OCT-29 to 2009-MAR-31
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef AS_UTL_FASTA_H
#define AS_UTL_FASTA_H

static const char *rcsid_AS_UTL_FASTA_H = "$Id$";

#include "AS_global.H"

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

//  Writes sequence as fasta, with at most 'bl' letters per line (unlimited if 0).
void
AS_UTL_writeFastA(FILE *f,
                  char *s, int sl, int bl,
                  char *h, ...);

//  Writes QVs as decimal 'fasta' ("00 00 00 00 ...") with up to 'bl' QVs per line.
void
AS_UTL_writeQVFastA(FILE *f,
                    char *q, int ql, int bl,
                    char *h, ...);

//  Writes FastQ, converting CA QVs into Sanger QVs.
void
AS_UTL_writeFastQ(FILE *f,
                  char *s, int sl,
                  char *q, int ql,
                  char *h, ...);

#endif


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
 *    Sergey Koren from 2008-JUN-04 to 2009-MAR-31
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2013-MAR-20
 *      are Copyright 2013 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz on 2015-FEB-27
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-23
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_UTL_fasta.H"
#include "AS_UTL_fileIO.H"

#include <stdarg.h>


void
AS_UTL_writeFastA(FILE *f,
                  char *s, int sl, int bl,
                  char *h, ...) {
  va_list ap;
  char   *o  = new char [sl + sl / 60 + 2];
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

  delete [] o;
}


void
AS_UTL_writeFastQ(FILE *f,
                  char *s, int sl,
                  char *q, int ql,
                  char *h, ...) {
  va_list ap;
  char   *o  = new char [ql + 1];
  int     qi = 0;
  int     oi = 0;

  assert(sl == ql);

  //  Reencode the QV to the Sanger spec.  This is a copy, so we don't need to undo it.
  while (qi < ql)
    o[oi++] = q[qi++] + '!';
  o[oi] = 0;

  va_start(ap, h);
  vfprintf(f, h, ap);
  va_end(ap);

  AS_UTL_safeWrite(f, s, "AS_UTL_writeFastQ", sizeof(char), sl);
  fprintf(f, "\n");

  fprintf(f, "+\n");
  AS_UTL_safeWrite(f, o, "AS_UTL_writeFastQ", sizeof(char), ql);
  fprintf(f, "\n");

  delete [] o;
}


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
 *    Brian P. Walenz from 2005-AUG-24 to 2013-AUG-01
 *      are Copyright 2005-2008,2011-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter on 2007-AUG-30
 *      are Copyright 2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static const char *rcsid = "$Id$";

//  We explicitly do not include AS_UTL_alloc.h here, because it
//  redefines malloc(), calloc(), realloc() and free() to be errors.
//  We want everyone to use the safe_*() versions supplied here.
//
//#include "AS_UTL_alloc.H"

//  We want to include AS_global.h to get the F_SIZE_T definition, but that
//  includes AS_UTL_alloc.h, so we pretend we've already included it.
//
#define AS_UTL_ALLOC_H

#include "AS_global.H"


void *
safe_calloc(size_t num, size_t len) {

  if ((num == 0) || (len == 0))
    return(NULL);    //  Bail, user didn't request anything.

  void  *p = calloc(num, len);

  if (p == NULL)
    fprintf(stderr, "Could not calloc memory ("F_SIZE_T" * "F_SIZE_T" bytes = "F_SIZE_T")\n",
            num, len, num*len);
  assert(p != NULL);

  return(p);
}


void *
safe_malloc(size_t len) {

  if (len == 0)
    return(NULL);    //  Bail, user didn't request anything.

  void  *p = malloc(len);

  if (p == NULL)
    fprintf(stderr, "Could not malloc memory ("F_SIZE_T" bytes)\n", len);
  assert(p != NULL);

  return(p);
}


void *
safe_realloc(void *q, size_t len) {

  if (len == 0)
    return(NULL);    //  Bail, user didn't request anything.

  void  *p = realloc(q, len);

  if (p == NULL)
    fprintf(stderr, "Could not realloc memory ("F_SIZE_T" bytes)\n", len);
  assert(p != NULL);

  return(p);
}


void
safe_free2(void *q) {
  free(q);
}

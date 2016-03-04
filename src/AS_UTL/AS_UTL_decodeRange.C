
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
 *    Brian P. Walenz from 2012-FEB-12 to 2013-OCT-11
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2015-MAY-28
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_UTL_decodeRange.H"


void
AS_UTL_decodeRange(char *range, set<uint64> &ranges) {
  char    *ap = range;
  uint32   av = 0,      bv = 0;

  while (*ap != 0) {
    av = strtoull(ap, &ap, 10);

    if (*ap == ',') {
      ap++;
      ranges.insert(av);

    } else if (*ap == 0) {
      ranges.insert(av);

    } else if (*ap == '-') {
      ap++;
      bv = strtoull(ap, &ap, 10);

      for (uint32 xx=av; xx<=bv; xx++)
        ranges.insert(xx);

      if (*ap == ',')
        ap++;

    } else if (*ap != 0) {
      fprintf(stderr, "ERROR: invalid range '%s'\n", range);
      exit(1);
    }
  }
}


void
AS_UTL_decodeRange(char *range, set<uint32> &ranges) {
  char    *ap = range;
  uint32   av = 0,      bv = 0;

  while (*ap != 0) {
    av = strtoul(ap, &ap, 10);

    if (*ap == ',') {
      ap++;
      ranges.insert(av);

    } else if (*ap == 0) {
      ranges.insert(av);

    } else if (*ap == '-') {
      ap++;
      bv = strtoull(ap, &ap, 10);

      for (uint32 xx=av; xx<=bv; xx++)
        ranges.insert(xx);

      if (*ap == ',')
        ap++;

    } else if (*ap != 0) {
      fprintf(stderr, "ERROR: invalid range '%s'\n", range);
      exit(1);
    }
  }
}


void
AS_UTL_decodeRange(char *range, uint64 &lo, uint64 &hi) {
  char    *ap = range;

  lo = hi = strtoull(ap, &ap, 10);

  if (*ap == '-') {
    ap++;
    hi = strtoull(ap, &ap, 10);

  } else if (*ap != 0) {
    fprintf(stderr, "ERROR: invalid range '%s'\n", range);
    exit(1);
  }
}


void
AS_UTL_decodeRange(char *range, int64 &lo, int64 &hi) {
  char    *ap = range;

  lo = hi = strtoll(ap, &ap, 10);

  if (*ap == '-') {
    ap++;
    hi = strtoll(ap, &ap, 10);

  } else if (*ap != 0) {
    fprintf(stderr, "ERROR: invalid range '%s'\n", range);
    exit(1);
  }
}


void
AS_UTL_decodeRange(char *range, uint32 &lo, uint32 &hi) {
  char    *ap = range;

  lo = hi = strtoul(ap, &ap, 10);

  if (*ap == '-') {
    ap++;
    hi = strtoul(ap, &ap, 10);

  } else if (*ap != 0) {
    fprintf(stderr, "ERROR: invalid range '%s'\n", range);
    exit(1);
  }
}


void
AS_UTL_decodeRange(char *range, int32 &lo, int32 &hi) {
  char    *ap = range;

  lo = hi = strtol(ap, &ap, 10);

  if (*ap == '-') {
    ap++;
    hi = strtol(ap, &ap, 10);

  } else if (*ap != 0) {
    fprintf(stderr, "ERROR: invalid range '%s'\n", range);
    exit(1);
  }
}


void
AS_UTL_decodeRange(char *range, double &lo, double &hi) {
  char    *ap = range;

  lo = hi = strtod(ap, &ap);

  if (*ap == '-') {
    ap++;
    hi = strtod(ap, &ap);

  } else if (*ap != 0) {
    fprintf(stderr, "ERROR: invalid range '%s'\n", range);
    exit(1);
  }
}


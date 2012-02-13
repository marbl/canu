
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2012, J. Craig Venter Institute.
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

static const char *rcsid = "$Id: AS_UTL_decodeRange.C,v 1.2 2012-02-13 08:30:41 brianwalenz Exp $";

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


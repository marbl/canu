
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, 2008, J. Craig Venter Institute.
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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include "AS_global.h"

//  Examines all VAR records, counts the number of fragments that
//  contain a var (partially or completely), and the number that are
//  contained within a single VAR.

typedef struct {
  uint64   fid;
  uint64   sid;
  uint64   ori:1;
  uint64   beg:31;
  uint64   end:31;
} frag;

typedef struct {
  uint64   vid;
  uint64   sid;
  uint32   beg;
  uint32   end;
  uint32   var;  //  Offset into global var string
} var;

int
fragCompare(const void *a, const void *b) {
  const frag *A = (const frag *)a;
  const frag *B = (const frag *)b;
  if (A->sid < B->sid)  return(-1);
  if (A->sid > B->sid)  return( 1);
  if (A->beg < B->beg)  return(-1);
  if (A->beg > B->beg)  return( 1);
  if (A->end < B->end)  return(-1);
  if (A->end > B->end)  return( 1);
  return(0);
}

int
varCompare(const void *a, const void *b) {
  const var *A = (const var *)a;
  const var *B = (const var *)b;
  if (A->sid < B->sid)  return(-1);
  if (A->sid > B->sid)  return( 1);
  if (A->beg < B->beg)  return(-1);
  if (A->beg > B->beg)  return( 1);
  if (A->end < B->end)  return(-1);
  if (A->end > B->end)  return( 1);
  return(0);
}


void
readFrags(FILE *frgscf, frag **frags, uint32 *fragsLen, uint32 *fragsMax) {
  char          L[1024];
  char         *l;

  fgets(L, 1024, frgscf);

  while (!feof(frgscf)) {

    if (*fragsLen >= *fragsMax) {
      *fragsMax *= 2;
      *frags     = (frag *)safe_realloc(*frags, sizeof(frag) * *fragsMax);
    }

    l = L;

    (*frags)[*fragsLen].fid = strtoull(l, &l, 10);
    (*frags)[*fragsLen].sid = strtoull(l, &l, 10);
    (*frags)[*fragsLen].beg = strtoull(l, &l, 10);
    (*frags)[*fragsLen].end = strtoull(l, &l, 10);
    (*frags)[*fragsLen].ori = strtoull(l, &l, 10);

    (*fragsLen)++;

    fgets(L, 1024, frgscf);
  }

  fprintf(stderr, "Read "F_U32" fragments.\n", *fragsLen);
}


void
readVars(FILE    *varscf,
         var    **vars, uint32 *varsLen, uint32 *varsMax,
         char   **vardata,
         uint32  *vardataLen,
         uint32  *vardataMax) {
  char          L[10240];
  char         *l;

  uint64        vid = 0;
  L[10240 - 1] = 0;
  fgets(L, 10240, varscf);

  while (!feof(varscf)) {
    if (L[10240 - 1] != 0) {
      fprintf(stderr, "WARNING:  Potentially long line at vid="F_U64": '%s'\n", vid, L);
      exit(1);
    }

#if 0
    //  Check if this VAR looks like it has more than 2 alleles.  We
    //  don't usually want to limit our reads to this.  If enabled,
    //  you need to also include "&& suspicious) in the "ignore SNP"
    //  test below.
    //
    bool  suspicious = false;
    if (0) {
      char  *vt = L;

      uint32  beforeSlash = 0;
      uint32  afterSlash  = 0;

      while (*vt != '/') {
        if (*vt == '-')
          beforeSlash++;
        vt++;
      }
      vt++;
      while (*vt) {
        if (*vt == '-')
          afterSlash++;
        vt++;
      }

      if ((beforeSlash > 0) && (afterSlash > 0))
        suspicious++;
    }
#endif

    //  Ignore SNP vars.
    if (L[1] != '/') {
      if (*varsLen >= *varsMax) {
        *varsMax *= 2;
        *vars     = (var *)safe_realloc(*vars, sizeof(var) * *varsMax);
      }

      fprintf(stderr, "L=%s", L);

      l = L;
      while (!isspace(*l))
        l++;

      (*vars)[*varsLen].var = *vardataLen;
      (*vars)[*varsLen].vid = vid++;
      (*vars)[*varsLen].sid = strtoull(l, &l, 10);
      (*vars)[*varsLen].beg = strtoull(l, &l, 10);
      (*vars)[*varsLen].end = strtoull(l, &l, 10);

      (*varsLen)++;

      if (*vardataLen + strlen(L) + 256 > *vardataMax) {
        *vardataMax *= 2;
        *vardata = (char *)safe_realloc(*vardata, sizeof(char) * *vardataMax);
      }

      char  *vt = L;
      while (!isspace(*vt))
        (*vardata)[(*vardataLen)++] = *(vt++);
      (*vardata)[(*vardataLen)++] = 0;
    }

    fgets(L, 10240, varscf);
  }

  fprintf(stderr, "Read "F_U32" var records.\n", *varsLen);
  fprintf(stderr, "Using "F_U32" letters to store variants.\n", *vardataLen);
}


int
main(int argc, char **argv) {
  FILE   *varscf   = 0L;
  FILE   *frgscf   = 0L;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-varscf") == 0) {
      errno = 0;
      varscf = fopen(argv[++arg], "r");
      if (errno)
        fprintf(stderr, "Failed to open varscf '%s': %s\n", argv[arg], strerror(errno)), exit(1);
    } else if (strcmp(argv[arg], "-frgscf") == 0) {
      errno = 0;
      frgscf = fopen(argv[++arg], "r");
      if (errno)
        fprintf(stderr, "Failed to open varscf '%s': %s\n", argv[arg], strerror(errno)), exit(1);
    } else {
      fprintf(stderr, "error: unknown arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err > 0) || (varscf == 0L) || (frgscf == 0L)) {
    fprintf(stderr, "usage: %s -frgscf X -varscf X\n", argv[0]);
    fprintf(stderr, "  Examines all VAR records, counts the number of fragments that\n");
    fprintf(stderr, "  contain a var (partially or completely), and the number that are\n");
    fprintf(stderr, "  contained within a single VAR.\n");
    exit(1);
  }

  //  Read all the frags.

  uint32  fragsLen = 0;
  uint32  fragsMax = 32;
  frag   *frags = (frag *)safe_malloc(sizeof(frag) * fragsMax);

  readFrags(frgscf, &frags, &fragsLen, &fragsMax);

  uint32  vardataLen = 0;
  uint32  vardataMax = 1024 * 1024;
  char   *vardata = (char *)safe_malloc(sizeof(char) * vardataMax);

  uint32  varsLen = 0;
  uint32  varsMax = 2;
  var    *vars = (var *)safe_malloc(sizeof(var) * varsMax);

  readVars(varscf, &vars, &varsLen, &varsMax, &vardata, &vardataLen, &vardataMax);

  fprintf(stderr, "Sorting frags.\n");
  qsort(frags, fragsLen, sizeof(frag), fragCompare);

  fprintf(stderr, "Sorting vars.\n");
  qsort(vars,  varsLen,  sizeof(var),  varCompare);

  fprintf(stderr, "Scanning.\n");

  uint32 fragsContainV  = 0;
  uint32 fragsRight     = 0;
  uint32 fragsLeft      = 0;
  uint32 fragsContained = 0;

  uint32 varsNoFrags    = 0;

  uint32  fp = 0;
  uint32  vp = 0;

  while ((fp < fragsLen) &&
         (vp < varsLen)) {

    fprintf(stdout, "VAR\t"F_U64"\t%s\t"F_U64"\t"F_U32"\t"F_U32"\n",
            vars[vp].vid, vardata + vars[vp].var, vars[vp].sid, vars[vp].beg, vars[vp].end);

    //  Advance the frag to the scaffold this var is in.  We can
    //  assume a fragment is found.
    //
    while (frags[fp].sid > vars[vp].sid)
      fp--;
    while (frags[fp].sid < vars[vp].sid)
      fp++;

    //  We need to backup to the first fragment that intersects this
    //  region.  We'll later iterate until we get to the first
    //  fragment whose begin is after the end.
    //
    //  Unfortunately, the frags are sorted by begin position, and so
    //  we really don't know when we find the first fragment that
    //  intersects this region.
    //
    //  Example: if we backup until the first fragment that doesn't
    //  intersect, we can miss fragments due to short guys:
    //
    //     ---------------------|------|--------------
    //                  -------
    //                     --------
    //                       -- <-------- short guy that stops backtracking.
    //                        -------
    //                          -----------
    //
    //  We can't simply sort by the end coordinate, since we'll have
    //  the exact same problem when we iterate forward to find the
    //  intersecting reads.
    //
    //  It runs too (way too) slow to backup to the start of the
    //  scaffold, so we backup to the first fragment that is 2KB
    //  before the region.  2KB is, of course, the maximum size of a
    //  read.
    //
    while ((frags[fp].sid == vars[vp].sid) &&
           (frags[fp].end + 2048 > vars[vp].beg) &&
           (fp > 0))
      fp--;

    //  Stupid special case
    if (fp == 1)
      fp = 0;

    while (frags[fp].sid < vars[vp].sid)
      fp++;

    //  Scan forward in the frags until we get past the var.  Remember
    //  any frag that intersects.  We don't care about frags that
    //  cover half the variant, since we want to find out if any frag
    //  has the whole variant!

    uint32 fV = 0;
    uint32 fR = 0;
    uint32 fL = 0;
    uint32 fC = 0;

    while ((frags[fp].sid == vars[vp].sid) &&
           (frags[fp].beg <= vars[vp].end) &&
           (fp < fragsLen)) {

      //  We could be a little more aggressive on the above while
      //  (frag.beg < var.beg), but we also want to count the number
      //  of frags with ends in the var.

      if ((frags[fp].beg < vars[vp].beg) &&
          (frags[fp].end > vars[vp].end)) {
        //  var is contained in frag.
        fprintf(stdout, "FRAG\t"F_U64"\t"F_U64"\t"F_U32"\t"F_U32" COMPLETE\n",
                frags[fp].fid, frags[fp].sid, frags[fp].beg, frags[fp].end);
        fV++;

      } else if ((frags[fp].beg < vars[vp].beg) &&
                 (frags[fp].end > vars[vp].beg)) {
        //  var hangs off right end of frag
        fprintf(stdout, "FRAG\t"F_U64"\t"F_U64"\t"F_U32"\t"F_U32" PARTIAL_RIGHT\n",
                frags[fp].fid, frags[fp].sid, frags[fp].beg, frags[fp].end);
        fR++;

      } else if ((frags[fp].beg < vars[vp].end) &&
                 (frags[fp].end > vars[vp].end)) {
        //  var hangs off left end of frag
        fprintf(stdout, "FRAG\t"F_U64"\t"F_U64"\t"F_U32"\t"F_U32" PARTIAL_LEFT\n",
                frags[fp].fid, frags[fp].sid, frags[fp].beg, frags[fp].end);
        fL++;

      } else if ((frags[fp].end > vars[vp].beg) &&
                 (frags[fp].beg < vars[vp].end)) {
        //  frag is contained in var
        fprintf(stdout, "FRAG\t"F_U64"\t"F_U64"\t"F_U32"\t"F_U32" CONTAINED\n",
                frags[fp].fid, frags[fp].sid, frags[fp].beg, frags[fp].end);
        fC++;
      }

      fp++;
    }

    if (fV + fR + fL + fC == 0) {
      fprintf(stderr, "WARNING:  No fragments found for vid "F_U64"\n", vars[vp].vid);
      varsNoFrags++;
    }

    fragsContainV  += fV;
    fragsRight     += fR;
    fragsLeft      += fL;
    fragsContained += fC;

    vp++;
  }

  safe_free(frags);
  safe_free(vardata);
  safe_free(vars);

  fprintf(stderr, "fragsContain VAR: "F_U32"\n", fragsContainV);
  fprintf(stderr, "fragsRight:       "F_U32"\n", fragsRight);
  fprintf(stderr, "fragsLeft:        "F_U32"\n", fragsLeft);
  fprintf(stderr, "fragsContained:   "F_U32"\n", fragsContained);
  fprintf(stderr, "varsNoFrags:      "F_U32"\n", varsNoFrags);
}

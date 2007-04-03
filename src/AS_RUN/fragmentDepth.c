
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

#include "AS_global.h"  //  only for CDS_UID_t, sigh.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//  A quick hack to compute a histogram of coverage depth using
//  the runCA-OBT posmap files.

#define HISTMAX     (8192)
#define DEPTHSIZE   (128 * 1024 * 1024)
#define FRAGMAX     (1024 * 1024)


typedef struct {
  uint32    lo;
  uint32    hi;
  uint32    de;
} intDep;


static
int
intDep_sort(const void *a, const void *b) {
  intDep *A = (intDep *)a;
  intDep *B = (intDep *)b;

  if (A->lo < B->lo) return(-1);
  if (A->lo > B->lo) return(1);
  return(0);
}


int
main(int argc, char **argv) {
  int              i = 0;

  CDS_UID_t        uidjunk = 0;
  CDS_UID_t        uid = 0;
  int              beg = 0;
  int              end = 0;

  CDS_UID_t        lastuid = 0;
  int              lastend = 0;

  int              histogram[HISTMAX] = { 0 };
  int              histmax = 0;

  int              minSize = 0;
  int              maxSize = DEPTHSIZE;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-min") == 0) {
      minSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-max") == 0) {
      maxSize = atoi(argv[++arg]);
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if (err) {
    fprintf(stderr, "usage: %s [-min N] [-max N] < x.posmap.frgscf", argv[0]);
    fprintf(stderr, "  -min N     count scaffolds at least N bases long.\n");
    fprintf(stderr, "  -max N     count scaffolds at most N bases long.\n");
  }

  uint32   inlen = 0;
  uint32   inmax = 1048576;
  intDep  *in    = (intDep *)safe_malloc(sizeof(intDep) * inmax);

  while (4 == fscanf(stdin, " "F_UID" "F_UID" %d %d %*d ", &uidjunk, &uid, &beg, &end)) {
    //fprintf(stderr, "read "F_UID" "F_UID" %d %d\n", uidjunk, uid, beg, end);

    //  Did we switch to a new scaffold?  Process this set of intervals.
    //
    if ((uid != lastuid) &&
        (inlen > 0)) {

      //  This scaffold is the correct size
      //
      if ((minSize <= lastend) &&
          (lastend <= maxSize)) {
        uint32   idlen = 0;
        intDep  *id    = NULL;

        //  Convert the list of overlapping intervals into a list
        //  of non-intersecting intervals annotated with depth
        uint32   islen = inlen * 2;
        intDep  *is    = (intDep *)safe_malloc(sizeof(intDep) * islen);

        for (i=0; i<inlen; i++) {
          is[2*i  ].lo = in[i].lo;
          is[2*i  ].hi = 0;
          is[2*i  ].de = 1;
          is[2*i+1].lo = in[i].hi;
          is[2*i+1].hi = 0;
          is[2*i+1].de = 0;
        }

        qsort(is, islen, sizeof(intDep), intDep_sort);

        //  Scan the list, counting how many times we change depth.
        //
        idlen = 1;
        for (i=1; i<islen; i++) {
          if (is[i-1].lo != is[i].lo)
            idlen++;
        }

        //  Allocate the real depth of coverage intervals
        //
        id    = (intDep *)safe_malloc(sizeof(intDep) * idlen);
        idlen = 0;

        //  Build new intervals
        //
        //  Initialize the first interval
        //
        id[idlen].lo = is[0].lo;
        id[idlen].hi = is[0].lo;
        id[idlen].de = 1;

        for (i=1; i<islen; i++) {

          if (id[idlen].de == 0) {
            //  Update the start position if the current interval is at zero
            //  depth.
            //
            id[idlen].lo = is[i].lo;
          } else {

            //  If we are at a position different from the start, we need to
            //  close out the current interval and make a new one.
            //
            if (is[i-1].lo != is[i].lo) {
              id[idlen].hi = is[i].lo;

              idlen++;

              id[idlen].lo = is[i].lo;
              id[idlen].hi = is[i].lo;
              id[idlen].de = id[idlen-1].de;
            }
          }

          //  Finally, update the depth of the current interval
          //
          if (is[i].de)
            id[idlen].de++;
          else
            id[idlen].de--;
        }

        //  Toss out the last one if it's zero length -- I think it's always
        //  zero length, just can't convince myself.
        //
        if (id[idlen].lo == id[idlen].hi)
          idlen--;

        safe_free(is);

        //  Update the histogram
        //
        for (i=0; i<idlen; i++) {
          if (id[i].de < HISTMAX) {
            histogram[id[i].de] += id[i].hi - id[i].lo;
            if (histmax < id[i].de)
              histmax = id[i].de;
          }
        }

        safe_free(id);
      }  //  scaffold is correct size

      //  Setup for the next scaffold
      //
      inlen = 0;

      lastuid = uid;
      lastend = 0;
    }  //  got a new scaffold

    //  Save this fragment.
    //
    in[inlen].lo = beg;
    in[inlen].hi = end;
    inlen++;

    if (lastend < end)
      lastend = end;
  }

  for (i=0; i<=histmax; i++)
    fprintf(stdout, "%d\t%d\n", i, histogram[i]);

  safe_free(in);

  exit(0);
}

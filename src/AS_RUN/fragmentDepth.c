
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute.
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
#include <assert.h>

//  A quick hack to compute a histogram of coverage depth using
//  the runCA-OBT posmap files.

#define HISTMAX     (8192)
#define DEPTHSIZE   (128 * 1024 * 1024)
#define FRAGMAX     (1024 * 1024)
#define DEFAULT_MEMORY_LIMIT 128


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


void
computeStuff(uint32 *V, uint32 N,
             uint32  B,
             uint32  E,
             uint32  *mode,
             double  *mean,
             uint32  *median) {

  uint32  histogramMax = 128 * 1024;
  uint32 *histogram    = (uint32 *)safe_malloc(sizeof(uint32) * histogramMax);
  uint32  histogramBig = 0;
  uint32  meanCount    = 0;
  uint32  i;

  if (E > N)
    E = N;

  *mean = 0;
  for (i=B; i<E; i++) {
    if (V[i] > 0) {
      *mean += V[i];
      meanCount++;
    }

    if (V[i] < histogramMax)
      histogram[V[i]]++;
    else 
      histogramBig++;
  }

  if (histogramBig) {
    fprintf(stderr, "histogramBig: "F_U32"\n", histogramBig);
    exit(1);
  }

  //  Find the mode -- except for 0.
  //
  *mode = 1;
  for (i=1; i<histogramMax; i++) {
    if (histogram[*mode] < histogram[i])
      *mode = i;
  }

  //  Find the mean
  //
  if (meanCount == 0) {
   *mean = 0;
  }
  else {
   *mean = *mean / meanCount;
  }


  //  Find the median
  //
  meanCount /= 2;
  *median    = 1;

  for (i=1; i<histogramMax; i++)
    if (meanCount >= histogram[i]) {
      meanCount -= histogram[i];
    } else {
      *median = i;
      break;
    }

  safe_free(histogram);
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

  int              doScaffold = 0;
  int              stepSize = 0;
  uint32           memLimit = DEFAULT_MEMORY_LIMIT;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-min") == 0) {
      minSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-max") == 0) {
      maxSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-scaffold") == 0) {
      doScaffold = 1;
    } else if (strcmp(argv[arg], "-stepSize") == 0) {
      stepSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-memLimit") == 0) {
      memLimit = atoi(argv[++arg]);
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if (err) {
    fprintf(stderr, "usage: %s [-min N] [-max N] [-scaffold] [-stepSize N] [-memLimit N] < x.posmap.frgscf", argv[0]);
    fprintf(stderr, "  Default is to compute a histogram of the number of bases at some\n");
    fprintf(stderr, "  depth of coverage.  Options for this are:\n");
    fprintf(stderr, "    -min N     use scaffolds at least N bases long.\n");
    fprintf(stderr, "    -max N     use scaffolds at most N bases long.\n");
    fprintf(stderr, "    -stepSize N     Used together with -scaffold, compute stats per stepSize per scaffold.\n");
    fprintf(stderr, "    -memLimit N     Used to limit the maximum memory to be used while running. Default = %dMB.\n", DEFAULT_MEMORY_LIMIT);
    fprintf(stderr, "\n");
    fprintf(stderr, "  The -scaffold switch disables the default output, and instead reports\n");
    fprintf(stderr, "  the mode, mean, median per scaffold.\n");
    exit(1);
  }

  uint32   inlen = 0;
  uint32     inmax = (((uint64) memLimit) * 1024 * 1024) / (sizeof(intDep));
  intDep    *in    = (intDep *)safe_malloc(sizeof(intDep) * inmax);
  if (doScaffold)
    fprintf(stdout, "uid\tstart\tend\tmode\tmean\tmedian\n");

  while (4 == fscanf(stdin, " "F_UID" "F_UID" %d %d %*d ", &uidjunk, &uid, &beg, &end)) {
    if (lastuid == 0)
      lastuid = uid;

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

        //  The way the loop is constructed above, the number of id
        //  intervals is idlen+1.  The last interval is always zero
        //  (thats id[idlen]) and so our later loops are supposed to
        //  be i<idlen.
        assert(id[idlen].lo == id[idlen].hi);

        safe_free(is);

        //  Update the histogram
        //
        if (doScaffold == 0) {
          for (i=0; i<idlen; i++) {
            if (id[i].de < HISTMAX) {
              histogram[id[i].de] += id[i].hi - id[i].lo;
              if (histmax < id[i].de)
                histmax = id[i].de;
            }
          }
        }

        //  Report mode, mean and median for thsi scaffold
        //
        if (doScaffold) {
          uint32  N      = id[idlen-1].hi;
          uint32 *V      = (uint32 *)safe_calloc(N, sizeof(uint32));
          uint32  mode   = 0;
          double  mean   = 0.0;
          uint32  median = 0;
          uint32  currStep = 0;

          for (i=0; i<idlen; i++) {
            int j;
            for (j=id[i].lo; j<id[i].hi; j++) {
              V[j] = id[i].de;
            }
          }
         
          if (stepSize == 0) { 
            currStep = N; 
          } 
          else {
            currStep = stepSize;
          }
          
          for (i = 0; i < N; i+=currStep) {                    
            int E  = i+currStep;            
            if (E > N) { E = N; }
            
            computeStuff(V, N, i, E, &mode, &mean, &median);
            fprintf(stdout, F_UID"\t"F_U32"\t"F_U32"\t"F_U32"\t%f\t"F_U32"\n", lastuid, i, E, mode, mean, median);
          }
          safe_free(V);
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

    if (inlen >= inmax) {
      fprintf(stderr, "The maxmimum limit of fragments per scaffold is %d and have already read %d\n", inmax, inlen);
      exit(1);
    }

    if (lastend < end)
      lastend = end;
  }

  if (doScaffold == 0)
    for (i=0; i<=histmax; i++)
      fprintf(stdout, "%d\t%d\n", i, histogram[i]);
  safe_free(in);

  exit(0);
}

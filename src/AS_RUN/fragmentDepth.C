
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
 *    Brian P. Walenz from 2007-MAR-28 to 2013-AUG-01
 *      are Copyright 2007-2008,2010,2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2007-SEP-04 to 2009-AUG-14
 *      are Copyright 2007-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-NOV-21
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "AS_UTL_fasta.H"

//  A quick hack to compute a histogram of coverage depth using
//  the runCA-OBT posmap files.

#define HISTMAX     (8192)
#define DEPTHSIZE   (128 * 1024 * 1024)
#define FRAGMAX     (1024 * 1024)

#define MODE_HISTOGRAM 0
#define MODE_SCAFFOLD  1
#define MODE_DEPTH     2

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
  uint32 *histogram    = (uint32 *)safe_calloc(histogramMax, sizeof(uint32));
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

void outputResult(AS_IID lastiid,
                  intDep *id,
                  uint32 idlen,
                  int mode,
                  uint32 *histogram,
                  uint32 *histmax,
                  int stepSize) {
   uint32 i = 0;
   uint32 lastpos = 0;

   switch (mode) {
      case MODE_HISTOGRAM:
      //  Update the histogram
      //
      for (i=0; i<idlen; i++) {
         if (id[i].de < HISTMAX) {
            // if there is a gap between the previous interval and the current one, add to our 0 coverage count
            if ((id[i].lo - lastpos) > 0) {
               histogram[0] += id[i].lo - lastpos;
            }

            histogram[id[i].de] += id[i].hi - id[i].lo;
            if ((*histmax) < id[i].de)
               (*histmax) = id[i].de;
         }
         lastpos = id[i].hi;
      }
      break;


      case MODE_SCAFFOLD:
         //  Report mode, mean and median for this scaffold
         //
         {
            uint32  N      = id[idlen-1].hi;
            uint32 *V      = (uint32 *)safe_calloc(N, sizeof(uint32));
            uint32  mode   = 0;
            double  mean   = 0.0;
            uint32  median = 0;
            uint32  currStep = 0;

            for (i=0; i<idlen; i++) {
              uint32 j;
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
              uint32 E  = i+currStep;
              if (E > N) { E = N; }

              computeStuff(V, N, i, E, &mode, &mean, &median);
              fprintf(stdout, "%s\t"F_U32"\t"F_U32"\t"F_U32"\t%f\t"F_U32"\n", lastiid, i, E, mode, mean, median);
            }
            safe_free(V);
          }
       break;

       case MODE_DEPTH:
          {
            char   *seq = (char *)safe_malloc((id[idlen-1].hi + 1) * sizeof(char));
            uint32     j;

            memset(seq, '0', id[idlen-1].hi);

            for (i=0; i<idlen; i++) {
              for (j=id[i].lo; j<id[i].hi; j++) {
                if (id[i].de < 10)
                  seq[j] = '0' + id[i].de;
                else if (id[i].de < 68)
                  seq[j] = 'A' + id[i].de - 10;
                else
                  seq[j] = '~';
               }
            }

            seq[id[idlen-1].hi] = 0;

            AS_UTL_writeFastA(stdout, seq, id[idlen-1].hi, 0, ">%d\n", lastiid);
            }
         break;
   }
}

void processScaffold(AS_IID lastiid,
                     intDep *in,
                     uint32 inlen,
                     int mode,
                     uint32 *histogram,
                     uint32 *histmax,
                     int stepSize) {
   uint32           i      = 0;
   uint32           idlen  = 0;
   intDep          *id     = NULL;

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
   outputResult(lastiid, id, idlen, mode, histogram, histmax, stepSize);
   safe_free(id);
}


int
main(int argc, char **argv) {
  uint32           i = 0;

  AS_IID           lastiid = NO_IID;
  int              lastend = 0;

  uint32           histogram[HISTMAX] = { 0 };
  uint32           histmax = 0;

  int              minSize = 0;
  int              maxSize = DEPTHSIZE;

  int              mode     = MODE_HISTOGRAM;
  int              stepSize = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-min") == 0) {
      minSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-max") == 0) {
      maxSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-stepSize") == 0) {
      stepSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-histogram") == 0) {
      mode = MODE_HISTOGRAM;
    } else if (strcmp(argv[arg], "-scaffold") == 0) {
      mode = MODE_SCAFFOLD;
    } else if (strcmp(argv[arg], "-depth") == 0) {
      mode = MODE_DEPTH;
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if (err || isatty(fileno(stdin))) {
    fprintf(stderr, "usage: %s MODE [-min N] [-max N] [-stepSize N] < x.posmap.frgscf\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -min N     use scaffolds at least N bases long.\n");
    fprintf(stderr, "  -max N     use scaffolds at most N bases long.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "MODES:  -histogram, -scaffold or -depth\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The default mode is to compute a histogram of the number of bases at some\n");
    fprintf(stderr, "depth of coverage.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The -scaffold mode reports the mode, mean, median depth per scaffold.  The\n");
    fprintf(stderr, "-stepSize option will compute those stats, in blocks of N bases (e.g., for bases\n");
    fprintf(stderr, "0 through N, then N through 2N, then 2N through 3N, etc.)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The -depth mode writes a multi-fasta file with the actual depth at each base\n");
    fprintf(stderr, "encoded.  The encoding is somewhat complicated to avoid using the '>' letter.\n");
    fprintf(stderr, "Depth 0 through 9 is encoded as '0' through '9'.  Depth 10 through 68 is\n");
    fprintf(stderr, "encoded as A-Z[\\]^_`a-z{|}, and depth more than 68 is encoded as ~.  Decode as:\n");
    fprintf(stderr, "  depth = letter - '0';\n");
    fprintf(stderr, "  if (depth > 9)\n");
    fprintf(stderr, "    depth -= 7;\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "!!WARNING -- The input frgscf MUST be sorted by scaffold ID -- WARNING!!\n");
    fprintf(stderr, "\n");
    exit(1);
  }






  uint32     inlen = 0;
  uint32     inmax = 4194304;  //  4 million fragments per scaffold should be enough (but we'll realloc later if not)
  intDep    *in    = (intDep *)safe_malloc(sizeof(intDep) * inmax);

  char       line[1024] = {0};
  char      *cont       = NULL;

  if (mode == MODE_SCAFFOLD)
    fprintf(stdout, "iid\tstart\tend\tmode\tmean\tmedian\n");

  while (fgets(line, 1024, stdin) != NULL) {
    AS_IID iidjunk = strtol(cont, &cont, 10);  //  read id
    AS_IID iid     = strtol(cont, &cont, 10);  //  scaffold id
    int32  beg     = strtol(cont, &cont, 10);
    int32  end     = strtol(cont, &cont, 10);

    if (lastiid == NO_IID)
      lastiid = iid;

    //  Did we switch to a new scaffold?  Process this set of intervals.
    //
    if ((iid != lastiid) && (inlen > 0)) {

      //  This scaffold is the correct size
      //
      if ((minSize <= lastend) &&
          (lastend <= maxSize)) {
         processScaffold(lastiid, in, inlen, mode, histogram, &histmax, stepSize);
      }

      //  Setup for the next scaffold
      //
      inlen = 0;

      lastiid = iid;
      lastend = 0;
    }  //  got a new scaffold

    //  Save this fragment.
    //
    in[inlen].lo = beg;
    in[inlen].hi = end;
    inlen++;

    if (inlen >= inmax) {
      inmax *= 2;
      in     = (intDep *)safe_realloc(in, inmax * sizeof(intDep));
    }

    if (lastend < end)
      lastend = end;
  }

  // process last scaffold
  if ((lastiid != NO_IID) && (inlen > 0)) {

     //  This scaffold is the correct size
     //
     if ((minSize <= lastend) &&
       (lastend <= maxSize)) {
       processScaffold(lastiid, in, inlen, mode, histogram, &histmax, stepSize);
    }
  }


  if (mode == MODE_HISTOGRAM)
    for (i=0; i<=histmax; i++)
      fprintf(stdout, "%d\t%d\n", i, histogram[i]);

  safe_free(in);

  exit(0);
}

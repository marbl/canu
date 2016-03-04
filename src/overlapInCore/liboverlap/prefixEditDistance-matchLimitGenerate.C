
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
 *    Brian P. Walenz from 2015-JUN-02 to 2015-AUG-14
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-FEB-05
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-15
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"

#include "Binomial_Bound.H"

#ifndef BROKEN_CLANG_OpenMP
#include <omp.h>
#endif

//  To use:
//
//  Set read length in gkStore.H.  Run this "50 5000 50" to compute data from 0.50% to 50.00%
//  in steps of 0.50%.  It will write three output files with the data:
//    *.C    - for inclusion in programs
//    *.bin  - binary dump of the array
//    *.dat  - ascii integer dump of the array
//
//  WARNING!  BITS=21 needs about 40 CPU hours.  (5000 took 2:40; 3800 took 1:28
//
//  The *.dat also contains the slope of the line from [0] to [i] and from [i] to [max] for each i.
//
//  To generate the slope parameters used in Binomial_Bound.C (BITS=20 and BITS=21 are expensive):
//
//  prefixEditDistance-matchLimitGenerate 100 5000 100
//
//  grep '^2000 ' *21/*dat | sed 's/.dat:/ /' | sed 's/prefixEditDistance-matchLimitData-BITS=[0-9][0-9]\/prefixEditDistance-matchLimit-//' > slopes
//
//  gnuplot:
//    f(x) = a/x+b
//    fit f(x) 'slopes' using 1:5 via a,b
//    plot 'slopes' using 1:5 with lines, f(x)
//    show var
//  plug a and b into Binomial_Bound.C
//

int
main(int argc, char **argv) {
  int32 minEvalue = 0;
  int32 maxEvalue = 0;
  int32 step      = 1;

  char D[FILENAME_MAX];
  char O[FILENAME_MAX];

  if        (argc == 2) {
    minEvalue  = atoi(argv[1]);
    maxEvalue  = minEvalue;

  } else if (argc == 3) {
    minEvalue  = atoi(argv[1]);
    maxEvalue  = atoi(argv[2]);

  } else if (argc == 4) {
    minEvalue  = atoi(argv[1]);
    maxEvalue  = atoi(argv[2]);
    step       = atoi(argv[3]);

  } else {
    fprintf(stderr, "usage: %s minEvalue [maxEvalue [step]]\n", argv[0]);
    fprintf(stderr, "  computes overlapper probabilities for minEvalue <= eValue <= maxEvalue'\n");
    fprintf(stderr, "    eValue 100 == 0.01 fraction error == 1%% error\n");
    exit(1);
  }

  fprintf(stderr, "Computing Edit_Match_Limit data for reads of length %ubp (bits = %u).\n", AS_MAX_READLEN, AS_MAX_READLEN_BITS);

  sprintf(D, "prefixEditDistance-matchLimitData-BITS=%01d", AS_MAX_READLEN_BITS);
  AS_UTL_mkdir(D);

#pragma omp parallel for schedule(dynamic, 1)
  for (int32 evalue=maxEvalue; evalue>=minEvalue; evalue -= step) {
    char    N[FILENAME_MAX];  //  Local to this thread!

    double  erate             = evalue / 10000.0;
    int32   start             = 1;

    int32   MAX_ERRORS        = (1 + (int) (erate * AS_MAX_READLEN));
    int32   ERRORS_FOR_FREE   = 1;

    int32  *starts            = new int32 [MAX_ERRORS + 1];

    memset(starts, 0, sizeof(int32) * (MAX_ERRORS + 1));

    sprintf(N, "%s/prefixEditDistance-matchLimit-%04d.bin", D, evalue);

    if (AS_UTL_fileExists(N)) {
      fprintf(stderr, "eValue %04d -- eRate %6.4f -- %7.4f%% error -- %8d values -- thread %2d - LOAD\n",
              evalue, erate, erate * 100.0, MAX_ERRORS, omp_get_thread_num());

      errno = 0;
      FILE *F = fopen(N, "r");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for reading: %s\n", N, strerror(errno)), exit(1);

      int32  me = 0;
      double er = 0.0;

      fread(&me,     sizeof(int32),  1,          F);
      fread(&er,     sizeof(double), 1,          F);
      fread( starts, sizeof(int32),  MAX_ERRORS, F);

      assert(me == MAX_ERRORS);
      assert(er == erate);

      fclose(F);

    } else {
      fprintf(stderr, "eValue %04d -- eRate %6.4f -- %7.4f%% error -- %8d values -- thread %2d - COMPUTE\n",
              evalue, erate, erate * 100.0, MAX_ERRORS, omp_get_thread_num());

      for (int32 e=ERRORS_FOR_FREE + 1; e<MAX_ERRORS; e++) {
        start = Binomial_Bound(e - ERRORS_FOR_FREE, erate, start);
        starts[e] = start - 1;
      }
    }



    {
      sprintf(O, "%s/prefixEditDistance-matchLimit-%04d.bin", D, evalue);

      errno = 0;
      FILE *F = fopen(O, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

      fwrite(&MAX_ERRORS, sizeof(int32),  1,          F);
      fwrite(&erate,      sizeof(double), 1,          F);
      fwrite( starts,     sizeof(int32),  MAX_ERRORS, F);

      fclose(F);
    }



    {
      sprintf(O, "%s/prefixEditDistance-matchLimit-%04d.dat", D, evalue);

      errno = 0;
      FILE *F = fopen(O, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

      fprintf(F, "#length     limit   slope0toX slopeXtoMAX for erate=%0.4f MAX_ERRORS=%d\n", erate, MAX_ERRORS);

      for (uint32 mm=MAX_ERRORS-1, ii=1; ii<MAX_ERRORS; ii++)
        fprintf(F, "%-8d %8d %11.6f %11.6f\n",
                ii,
                starts[ii],
                (double)(starts[ii] - starts[1])  / (ii -  1 + 1),
                (double)(starts[mm] - starts[ii]) / (mm - ii + 1));

      fclose(F);
    }



    {
      sprintf(O, "%s/prefixEditDistance-matchLimit-%04d.C", D, evalue);

      errno = 0;
      FILE *F = fopen(O, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

      fprintf(F, "//\n");
      fprintf(F, "//  Automagically generated.  Do not edit.\n");
      fprintf(F, "//\n");
      fprintf(F, "\n");
      fprintf(F, "#include \"gkStore.H\"\n");
      fprintf(F, "\n");
      fprintf(F, "#if (AS_MAX_READLEN_BITS == %d)\n", AS_MAX_READLEN_BITS);
      fprintf(F, "\n");
      fprintf(F, "extern\n");
      fprintf(F, "const\n");
      fprintf(F, "int32\n");
      fprintf(F, "Edit_Match_Limit_%04d[%d] = {\n", evalue, MAX_ERRORS + 1);

      uint32  i=0;

      while (i < MAX_ERRORS) {
        uint32  j=0;

        fprintf(F, "  ");

        while ((j < 16) && (i < MAX_ERRORS)) {
          if (i < MAX_ERRORS-1)
            fprintf(F, "0x%08x,", starts[i]);
          else
            fprintf(F, "0x%08x", starts[i]);

          i++;
          j++;
        }

        fprintf(F, "\n");
      }

      fprintf(F, "};\n");
      fprintf(F, "\n");
      fprintf(F, "#endif\n");

      fclose(F);
    }
  }
}

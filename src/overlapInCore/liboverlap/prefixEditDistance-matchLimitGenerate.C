#include "AS_global.H"
#include "gkStore.H"

#include "Binomial_Bound.H"

#include <omp.h>

//  For testing:
//    ca83 - 18
//    ca3g - 21
//
//#undef AS_MAX_READLEN_BITS
//#undef AS_MAX_READLEN
//
//#define AS_MAX_READLEN_BITS        18
//#define AS_MAX_READLEN             (((uint32)1 << AS_MAX_READLEN_BITS) - 1)
//

int
main(int argc, char **argv) {
  int32 minEvalue = 0;
  int32 maxEvalue = 0;
  int32 step      = 0;

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

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 evalue=minEvalue; evalue<=maxEvalue; evalue += step) {
    double  erate             = evalue / 10000.0;
    int32   start             = 1;

    int32   MAX_ERRORS        = (1 + (int) (erate * AS_MAX_READLEN));
    int32   ERRORS_FOR_FREE   = 1;

    int32  *starts            = new int32 [MAX_ERRORS + 1];

    memset(starts, 0, sizeof(int32) * (MAX_ERRORS + 1));

    char N[FILENAME_MAX];

    sprintf(N, "prefixEditDistance-matchLimitData/prefixEditDistance-matchLimit-%04d.dat", evalue);


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
      sprintf(N, "prefixEditDistance-matchLimitData/prefixEditDistance-matchLimit-%04d.dat", evalue);

      errno = 0;
      FILE *F = fopen(N, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

      fwrite(&MAX_ERRORS, sizeof(int32),  1,          F);
      fwrite(&erate,      sizeof(double), 1,          F);
      fwrite( starts,     sizeof(int32),  MAX_ERRORS, F);
      
      fclose(F);
    }



    {
      sprintf(N, "prefixEditDistance-matchLimit-%04d.C", evalue);

      errno = 0;
      FILE *F = fopen(N, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

      fprintf(F, "//\n");
      fprintf(F, "//  Automagically generated.  Do not edit.\n");
      fprintf(F, "//\n");
      fprintf(F, "\n");
      fprintf(F, "#include \"AS_global.H\"\n");
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

      fclose(F);
    }
  }
}
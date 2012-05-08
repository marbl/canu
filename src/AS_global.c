
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

static const char *rcsid = "$Id: AS_global.c,v 1.22 2012-05-08 23:31:32 brianwalenz Exp $";

#include "AS_global.h"

#ifdef X86_GCC_LINUX
#include <fpu_control.h>
#endif

#ifdef _GLIBCXX_PARALLEL
#include <parallel/algorithm>
#include <parallel/settings.h>
#endif

//  Nonsense values, mostly for making sure everybody that uses an
//  error rate calls AS_configure() at startup.

double AS_OVL_ERROR_RATE   = -90.0;
double AS_CGW_ERROR_RATE   = -90.0;
double AS_CNS_ERROR_RATE   = -90.0;
double AS_MAX_ERROR_RATE   =   0.25;

//  Historical Sanger defaults.  Seem to work well with 454 and Illumina,
//  but occasionally they need to be reduced (for, say, 50bp Illumina reads).

uint32 AS_READ_MIN_LEN     = 64;
uint32 AS_OVERLAP_MIN_LEN  = 40;

//  EVERY main program should define mainid.  The release manager
//  should fill in releaseid with the release name.

extern
const char *mainid;
const char *releaseid = "CVS TIP";

//  We take argc and argv, so, maybe, eventually, we'll want to parse
//  something out of there.  We return argc in case what we parse we
//  want to remove.
//
int
AS_configure(int argc, char **argv) {
  char *p = NULL;
  int   i, j;

#ifdef X86_GCC_LINUX
  //  Set the x86 FPU control word to force double precision rounding
  //  rather than `extended' precision rounding. This causes base
  //  calls and quality values on x86 GCC-Linux (tested on RedHat
  //  Linux) machines to be identical to those on IEEE conforming UNIX
  //  machines.
  //
  fpu_control_t fpu_cw = ( _FPU_DEFAULT & ~_FPU_EXTENDED ) | _FPU_DOUBLE;

  _FPU_SETCW( fpu_cw );
#endif

#ifdef _GLIBCXX_PARALLEL_SETTINGS_H
  __gnu_parallel::_Settings s;

  //  Force all algorithms to be parallel.
  //  Force some algs to be sequential by using a tag, eg:
  //    sort(a, a+end, __gnu_parallel::sequential_tag());
  //
  //s.algorithm_strategy = __gnu_parallel::force_parallel;

  //  The default seems to be 1000, way too small for us.
  s.sort_minimal_n = 128 * 1024;

  __gnu_parallel::_Settings::set(s);
#endif

  //
  //  Default values
  //

  AS_OVL_ERROR_RATE  = 0.06;
  AS_CGW_ERROR_RATE  = 0.10;
  AS_CNS_ERROR_RATE  = 0.06;

  //
  //  Environment
  //

  p = getenv("AS_OVL_ERROR_RATE");
  if (p)
    AS_OVL_ERROR_RATE = atof(p);

  p = getenv("AS_CGW_ERROR_RATE");
  if (p)
    AS_CGW_ERROR_RATE = atof(p);

  p = getenv("AS_CNS_ERROR_RATE");
  if (p)
    AS_CNS_ERROR_RATE = atof(p);

  p = getenv("AS_READ_MIN_LEN");
  if (p)
    AS_READ_MIN_LEN = atoi(p);

  p = getenv("AS_OVERLAP_MIN_LEN");
  if (p)
    AS_OVERLAP_MIN_LEN = atoi(p);

  //
  //  Command line
  //

  for (i=0; i<argc; i++) {
    if        (strcasecmp(argv[i], "--ovlErrorRate") == 0) {
      AS_OVL_ERROR_RATE = atof(argv[i+1]);
      for (j=i+2; j<argc; j++)
        argv[j-2] = argv[j];
      argv[--argc] = NULL;
      argv[--argc] = NULL;
      i--;

    } else if (strcasecmp(argv[i], "--cgwErrorRate") == 0) {
      AS_CGW_ERROR_RATE = atof(argv[i+1]);
      for (j=i+2; j<argc; j++)
        argv[j-2] = argv[j];
      argv[--argc] = NULL;
      argv[--argc] = NULL;
      i--;

    } else if (strcasecmp(argv[i], "--cnsErrorRate") == 0) {
      AS_CNS_ERROR_RATE = atof(argv[i+1]);
      for (j=i+2; j<argc; j++)
        argv[j-2] = argv[j];
      argv[--argc] = NULL;
      argv[--argc] = NULL;
      i--;

    } else if (strcasecmp(argv[i], "--frgMinLen") == 0) {
      AS_READ_MIN_LEN = atoi(argv[i+1]);
      for (j=i+2; j<argc; j++)
        argv[j-2] = argv[j];
      argv[--argc] = NULL;
      argv[--argc] = NULL;
      i--;

    } else if (strcasecmp(argv[i], "--ovlMinLen") == 0) {
      AS_OVERLAP_MIN_LEN = atoi(argv[i+1]);
      for (j=i+2; j<argc; j++)
        argv[j-2] = argv[j];
      argv[--argc] = NULL;
      argv[--argc] = NULL;
      i--;
    }
  }

  //
  //  Et cetera.
  //

  for (i=0; i<argc; i++) {
    if (strcmp(argv[i], "--version") == 0) {
      fprintf(stderr, "CA version %s (%s).\n", releaseid, mainid);
      exit(0);
    }
  }

  //
  //  Checking.
  //

  if ((AS_OVL_ERROR_RATE < 0.0) || (AS_MAX_ERROR_RATE < AS_OVL_ERROR_RATE))
    fprintf(stderr, "%s: ERROR:  Invalid AS_OVL_ERROR_RATE (%0.2f); should be between 0.0 and %0.2f\n",
            argv[0], AS_OVL_ERROR_RATE, AS_MAX_ERROR_RATE), exit(1);

  if ((AS_CGW_ERROR_RATE < 0.0) || (AS_MAX_ERROR_RATE < AS_CGW_ERROR_RATE))
    fprintf(stderr, "%s: ERROR:  Invalid AS_CGW_ERROR_RATE (%0.2f); should be between 0.0 and %0.2f\n",
            argv[0], AS_CGW_ERROR_RATE, AS_MAX_ERROR_RATE), exit(1);

  if ((AS_CNS_ERROR_RATE < 0.0) || (AS_MAX_ERROR_RATE < AS_CNS_ERROR_RATE))
    fprintf(stderr, "%s: ERROR:  Invalid AS_CNS_ERROR_RATE (%0.2f); should be between 0.0 and %0.2f\n",
            argv[0], AS_CNS_ERROR_RATE, AS_MAX_ERROR_RATE), exit(1);

  if (AS_READ_MIN_LEN < 1)
    AS_READ_MIN_LEN = 1;

  if (AS_OVERLAP_MIN_LEN < 1)
    AS_OVERLAP_MIN_LEN = 1;

#if 0
  if (AS_OVL_ERROR_RATE != 0.06)
    fprintf(stderr, "%s: AS_configure()-- AS_OVL_ERROR_RATE set to %0.2f\n", argv[0], AS_OVL_ERROR_RATE);

  if (AS_CGW_ERROR_RATE != 0.10)
    fprintf(stderr, "%s: AS_configure()-- AS_CGW_ERROR_RATE set to %0.2f\n", argv[0], AS_CGW_ERROR_RATE);

  if (AS_CNS_ERROR_RATE != 0.06)
    fprintf(stderr, "%s: AS_configure()-- AS_CNS_ERROR_RATE set to %0.2f\n", argv[0], AS_CNS_ERROR_RATE);
#endif

  //
  //  Logging.
  //

  p = getenv("AS_RUNCA_DIRECTORY");
  if (p) {
    char  D[FILENAME_MAX] = {0};
    char  N[FILENAME_MAX] = {0};
    char  H[1024]         = {0};  //  HOST_NAME_MAX?  Undefined.
    char *E;
    FILE *F;

    //  Make a directory for logs.  Allow all errors, in particular,
    //  the error of "directory already exists".
    //
    sprintf(D, "%s/runCA-logs", p);
    mkdir(D, S_IRWXU | S_IRWXG | S_IRWXO);

    //  Our hostname is part of our unique filename.
    //
    gethostname(H, 1024);

    //  Our executable name is part of our unique filename too.
    //
    E = argv[0] + strlen(argv[0]) - 1;
    while ((E != argv[0]) && (*E != '/'))
      E--;
    if (*E == '/')
      E++;

    //  Construct a name for this log, and open it.  If we can't open
    //  it, just skip the log.
    //
    sprintf(N, "%s/"F_U64"_%s_"F_U64"_%s",
            D,
            (uint64)time(NULL),
            H,
            (uint64)getpid(),
            E);

    //fprintf(stderr, "%s\n", N);

    errno = 0;
    F = fopen(N, "w");
    if ((errno == 0) && (F)) {
      fprintf(F, "CA version %s (%s).\n", releaseid, mainid);
      fprintf(F, "\n");

      fprintf(F, "Error Rates:\n");
      fprintf(F, "AS_OVL_ERROR_RATE %f\n", AS_OVL_ERROR_RATE);
      //fprintf(F, "AS_UTG_ERROR_RATE %f\n", AS_UTG_ERROR_RATE);
      fprintf(F, "AS_CNS_ERROR_RATE %f\n", AS_CNS_ERROR_RATE);
      fprintf(F, "AS_CGW_ERROR_RATE %f\n", AS_CGW_ERROR_RATE);
      fprintf(F, "AS_MAX_ERROR_RATE %f\n", AS_MAX_ERROR_RATE);
      fprintf(F, "\n");

      fprintf(F, "Current Working Directory:\n");
      fprintf(F, "%s\n", getcwd(N, FILENAME_MAX));
      fprintf(F, "\n");
      fprintf(F, "Command:\n");
      fprintf(F, "%s", argv[0]);

      for (i=1; i<argc; i++)
        if (argv[i][0] == '-')
          fprintf(F, " \\\n  %s", argv[i]);
        else
          fprintf(F, " %s", argv[i]);

      fprintf(F, "\n");

      fclose(F);
    }
  }

  return(argc);
}

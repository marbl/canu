
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

static const char *rcsid = "$Id$";

#include "AS_global.H"

#include "AS_UTL_stackTrace.H"

#ifdef X86_GCC_LINUX
#include <fpu_control.h>
#endif

#ifdef _GLIBCXX_PARALLEL
#include <parallel/algorithm>
#include <parallel/settings.h>
#endif

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
  __gnu_parallel::_Settings s = __gnu_parallel::_Settings::get();

  //  Force all algorithms to be parallel.
  //  Force some algs to be sequential by using a tag, eg:
  //    sort(a, a+end, __gnu_parallel::sequential_tag());
  //
  //s.algorithm_strategy = __gnu_parallel::force_parallel;

  //  The default seems to be 1000, way too small for us.
  s.sort_minimal_n = 128 * 1024;

  //  The default is MWMS, which, at least on FreeBSD 8.2 w/gcc46, is NOT inplace.
  //  Then again, the others also appear to be NOT inplace as well.
  //s.sort_algorithm = __gnu_parallel::MWMS;
  //s.sort_algorithm = __gnu_parallel::QS_BALANCED;
  //s.sort_algorithm = __gnu_parallel::QS;

  __gnu_parallel::_Settings::set(s);
#endif

  //  Install a signal handler to catch seg faults and errors.

  AS_UTL_installCrashCatcher();

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

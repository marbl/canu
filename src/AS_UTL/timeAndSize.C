
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
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>



double
getTime(void) {
  struct timeval  tp;
  gettimeofday(&tp, NULL);
  return(tp.tv_sec + (double)tp.tv_usec / 1000000.0);
}



static
bool
getrusage(struct rusage &ru) {

  errno = 0;

  if (getrusage(RUSAGE_SELF, &ru) == -1) {
    fprintf(stderr, "getrusage(RUSAGE_SELF, ...) failed: %s\n",
            strerror(errno));
    return(false);
  }

  return(true);
}



static
bool
getrlimit(struct rlimit &rl) {

  errno = 0;

  if (getrlimit(RLIMIT_DATA, &rl) == -1) {
    fprintf(stderr, "getrlimit(RLIMIT_DATA, ...) failed: %s\n",
            strerror(errno));
    return(false);
  }

  return(true);
}



double
getCPUTime(void) {
  struct rusage  ru;
  double         tm = 0;

  if (getrusage(ru) == true)
    tm  = ((ru.ru_utime.tv_sec + ru.ru_utime.tv_usec / 1000000.0) + 
           (ru.ru_stime.tv_sec + ru.ru_stime.tv_usec / 1000000.0));

  return(tm);
}



double
getProcessTime(void) {
  struct timeval tp;
  static double  st = 0.0;
  double         tm = 0;

  if (gettimeofday(&tp, NULL) == 0)
    tm  = tp.tv_sec + tp.tv_usec / 100000.0;

  if (st == 0.0)
    st = tm;

  return(tm - st);
}



uint64
getProcessSize(void) {
  struct rusage  ru;
  uint64         sz = 0;

  if (getrusage(ru) == true) {
    sz  = ru.ru_maxrss;
    sz *= 1024;
  }

  return(sz);
}



uint64
getProcessSizeLimit(void) {
  struct rlimit rl;
  uint64        sz = ~uint64ZERO;

  if (getrlimit(rl) == true)
    sz = rl.rlim_cur;

  return(sz);
}


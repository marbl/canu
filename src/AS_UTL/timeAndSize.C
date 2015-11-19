
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


uint64
getProcessSizeCurrent(void) {
  struct rusage  ru;
  uint64         sz = 0;

  errno = 0;
  if (getrusage(RUSAGE_SELF, &ru) == -1) {
    fprintf(stderr, "getProcessSizeCurrent()-- getrusage(RUSAGE_SELF, ...) failed: %s\n",
            strerror(errno));
  } else {
    sz  = ru.ru_maxrss;
    sz *= 1024;
  }

  return(sz);
}


uint64
getProcessSizeLimit(void) {
  struct rlimit rlp;
  uint64        sz = ~uint64ZERO;

  errno = 0;
  if (getrlimit(RLIMIT_DATA, &rlp) == -1) {
    fprintf(stderr, "getProcessSizeLimit()-- getrlimit(RLIMIT_DATA, ...) failed: %s\n",
            strerror(errno));
  } else {
    sz = rlp.rlim_cur;
  }

  return(sz);
}


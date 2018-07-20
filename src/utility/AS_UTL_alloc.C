
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_UTL_alloc.H"




#if !defined(__CYGWIN__) && !defined(_WIN32)
#include <sys/sysctl.h>
#endif

#ifdef HW_PHYSMEM

uint64
getPhysicalMemorySize(void) {
  uint64  physMemory = 0;

  int     mib[2] = { CTL_HW, HW_PHYSMEM };
  size_t  len    = sizeof(uint64);

  errno = 0;

  if (sysctl(mib, 2, &physMemory, &len, NULL, 0) != 0)
    fprintf(stderr, "getPhysicalMemorySize()-- sysctl() failed to return CTL_HW, HW_PHYSMEM: %s\n", strerror(errno)), exit(1);

  if (len != sizeof(uint64)) {
#ifdef HW_MEMSIZE
    mib[1] = HW_MEMSIZE;
    len = sizeof(uint64);
    if (sysctl(mib, 2, &physMemory, &len, NULL, 0) != 0 || len != sizeof(uint64))
#endif
      fprintf(stderr, "getPhysicalMemorySize()-- sysctl() failed to return CTL_HW, HW_PHYSMEM: %s\n", strerror(errno)), exit(1);
  }

  return(physMemory);
}

#else

uint64
getPhysicalMemorySize(void) {
  uint64  physPages  = sysconf(_SC_PHYS_PAGES);
  uint64  pageSize   = sysconf(_SC_PAGESIZE);
  uint64  physMemory = physPages * pageSize;

  return(physMemory);
}

#endif

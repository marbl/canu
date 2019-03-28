
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
 *  This file is derived from:
 *
 *    src/utility/timeAndSize.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef SYSTEM_H
#define SYSTEM_H

#include "AS_global.H"
#include <signal.h>



double   getTime(void);

double   getCPUTime(void);
double   getProcessTime(void);

uint64   getProcessSize(void);
uint64   getProcessSizeLimit(void);

uint64   getBytesAllocated(void);

uint64   getPhysicalMemorySize(void);

uint64   getPageSize(void);


void  AS_UTL_catchCrash(int sig_num, siginfo_t *info, void *ctx);

void  AS_UTL_installCrashCatcher(const char *filename);



#endif  //  SYSTEM_H

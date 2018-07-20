
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

#ifndef AS_UTL_STACKTRACE_H
#define AS_UTL_STACKTRACE_H

#include <signal.h>

void
AS_UTL_catchCrash(int sig_num, siginfo_t *info, void *ctx);

void
AS_UTL_installCrashCatcher(const char *filename);

#endif  //  AS_UTL_STACKTRACE_H

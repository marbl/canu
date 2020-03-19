
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
 *    src/AS_BAT/AS_BAT_Logging.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-JUL-29 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-AUG-11 to 2014-DEC-19
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

#ifndef INCLUDE_AS_BAT_LOGGING
#define INCLUDE_AS_BAT_LOGGING

#include "runtime.H"
#include "files.H"

void    setLogFile(char const *prefix, char const *name);
char   *getLogFilePrefix(void);

void    writeStatus(char const *fmt, ...);
void    writeLog(char const *fmt, ...);

void    flushLog(void);

#define logFileFlagSet(L) ((logFileFlags & L) == L)

extern uint64  logFileFlags;
extern uint32  logFileOrder;  //  Used debug tigStore dumps, etc

extern uint64 LOG_OVERLAP_SCORING;
extern uint64 LOG_BEST_EDGES;
extern uint64 LOG_BEST_OVERLAPS;
extern uint64 LOG_ERROR_PROFILES;
extern uint64 LOG_OPTIMIZE_POSITIONS;
extern uint64 LOG_CHUNK_GRAPH;
extern uint64 LOG_BUILD_UNITIG;
extern uint64 LOG_PLACE_UNPLACED;
extern uint64 LOG_ORPHAN_DETAIL;
extern uint64 LOG_SPLIT_DISCONTINUOUS;
extern uint64 LOG_INTERMEDIATE_TIGS;
extern uint64 LOG_SET_PARENT_AND_HANG;
extern uint64 LOG_STDERR;

extern uint64 LOG_PLACE_READ;

extern char const *logFileFlagNames[64];

#endif  //  INCLUDE_AS_BAT_LOGGING

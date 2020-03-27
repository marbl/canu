
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
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

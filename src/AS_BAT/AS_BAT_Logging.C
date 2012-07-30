
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id: AS_BAT_Logging.C,v 1.1 2012-07-30 01:21:01 brianwalenz Exp $";

#include "AS_BAT_Logging.H"

class logFileInstance {
public:
  logFileInstance() {
    file    = NULL;
    name[0] = 0;
  };
  ~logFileInstance() {
    if (file) {
      fprintf(stderr, "WARNING: open file\n");
      fclose(file);
    }
  };

  FILE *file;
  char  name[FILENAME_MAX];
};


logFileInstance   *logFile      = NULL;
uint32             logFileOrder = 0;
uint64             logFileFlags = 0;

uint64 LOG_OVERLAP_QUALITY             = 0x0000000000000001;  //  Debug, scoring of overlaps
uint64 LOG_OVERLAPS_USED               = 0x0000000000000002;  //  Report overlaps used/not used
uint64 LOG_CHUNK_GRAPH                 = 0x0000000000000004;  //  Report the chunk graph as we build it
uint64 LOG_INTERSECTIONS               = 0x0000000000000008;  //  Report intersections found when building initial unitigs
uint64 LOG_POPULATE_UNITIG             = 0x0000000000000010;  //  Report building of initial unitigs (both unitig creation and fragment placement)
uint64 LOG_INTERSECTION_BREAKING       = 0x0000000000000020;  //  
uint64 LOG_INTERSECTION_BUBBLES        = 0x0000000000000040;  //  
uint64 LOG_INTERSECTION_BUBBLES_DEBUG  = 0x0000000000000080;  //  
uint64 LOG_INTERSECTION_JOINING        = 0x0000000000000100;  //
uint64 LOG_INTERSECTION_JOINING_DEBUG  = 0x0000000000000200;  //
uint64 LOG_INITIAL_CONTAINED_PLACEMENT = 0x0000000000000400;  //
uint64 LOG_HAPPINESS                   = 0x0000000000000800;  //
uint64 LOG_INTERMEDIATE_UNITIGS        = 0x0000000000001000;  //  At various spots, dump the current unitigs
uint64 LOG_MATE_SPLIT_ANALYSIS         = 0x0000000000002000;  //
uint64 LOG_MATE_SPLIT_DISCONTINUOUS    = 0x0000000000004000;  //
uint64 LOG_MATE_SPLIT_UNHAPPY_CONTAINS = 0x0000000000008000;  //
uint64 LOG_MATE_SPLIT_COVERAGE_PLOT    = 0x0000000000010000;  //
uint64 LOG_SET_PARENT_AND_HANG         = 0x0000000000020000;  //
uint64 LOG_STDERR                      = 0x0000000000040000;  //  Write ALL logging to stderr, not the files.

uint64 LOG_PLACE_FRAG                  = 0x8000000000000000;  //  Internal use only.

char const *logFileFlagNames[64] = { "overlapQuality",
                                     "overlapsUsed",
                                     "chunkGraph",
                                     "intersections",
                                     "populate",
                                     "intersectionBreaking",
                                     "intersectionBubbles",
                                     "intersectionBubblesDebug",
                                     "intersectionJoining",
                                     "intersectionJoiningDebug",
                                     "containedPlacement",
                                     "happiness",
                                     "intermediateUnitigs",
                                     "mateSplitAnalysis",
                                     "mateSplitDiscontinuous",
                                     "mateSplitUnhappyContains",
                                     "mateSplitCoveragePlot",
                                     "setParentAndHang",
                                     "stderr",
                                     NULL
};

//  Closes the current logFile, opens a new one called 'prefix.logFileOrder.name'.  If 'name' is
//  NULL, the logFile is reset to stderr.
void
setLogFile(char const *prefix, char const *name) {

  if (logFileFlagSet(LOG_STDERR))
    //  Write everything to stderr
    return;

  int32  nt = omp_get_num_threads();
  int32  tn = omp_get_thread_num();

  if (logFile == NULL)
    logFile = new logFileInstance [omp_get_max_threads()];

  for (int32 i=0; i<omp_get_max_threads(); i++)
    if ((logFile[tn].file != NULL) && (logFile[tn].file != stderr)) {
      fclose(logFile[tn].file);
      logFile[tn].file    = NULL;
    }

  if (name == NULL) {
    logFile[tn].file    = stderr;
    logFile[tn].name[0] = 0;
    return;
  }

#pragma omp critical
  if (nt == 1)
    sprintf(logFile[tn].name, "%s.%03u.%s.log", prefix, ++logFileOrder, name);
  else
    sprintf(logFile[tn].name, "%s.%03u.%s.thread%03d.log", prefix, ++logFileOrder, name, tn);

  errno = 0;
  logFile[tn].file = fopen(logFile[tn].name, "w");
  if (errno) {
    fprintf(stderr, "setLogFile()-- Failed to open logFile '%s': %s.\n", logFile[tn].name, strerror(errno));
    fprintf(stderr, "setLogFile()-- Will now log to stderr instead.\n");
    logFile[tn].file = stderr;
  }

  fprintf(stderr,  "setLogFile()-- Now logging to '%s'\n", logFile[tn].name);
}


void
resetLogFile(char const *prefix, char const *name) {

  if (logFileFlagSet(LOG_STDERR))
    //  Write everything to stderr
    return;

  int32  tn = omp_get_thread_num();

  if ((logFile[tn].name[0] == 0) ||
      (AS_UTL_sizeOfFile(logFile[tn].name) > 512 * 1024 * 1024))
    setLogFile(prefix, name);
}


void
writeLog(char const *fmt, ...) {
  va_list   ap;

  va_start(ap, fmt);

  vfprintf(logFile[omp_get_thread_num()].file, fmt, ap);
}


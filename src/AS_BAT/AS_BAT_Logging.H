
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

#ifndef INCLUDE_AS_BAT_LOGGING
#define INCLUDE_AS_BAT_LOGGING

static const char *rcsid_INCLUDE_AS_BAT_LOGGING = "$Id$";

#include "AS_global.H"
#include "AS_UTL_fileIO.H"

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>

#include <omp.h>

void  setLogFile(char const *prefix, char const *name);
void  writeLog(char const *fmt, ...);

#define logFileFlagSet(L) ((logFileFlags & L) == L)

extern uint64  logFileFlags;
extern uint32  logFileOrder;  //  Used debug tigStore dumps, etc

extern uint64 LOG_OVERLAP_QUALITY;
extern uint64 LOG_OVERLAPS_USED;
extern uint64 LOG_CHUNK_GRAPH;
extern uint64 LOG_INTERSECTIONS;
extern uint64 LOG_POPULATE_UNITIG;
extern uint64 LOG_INTERSECTION_BREAKING;
extern uint64 LOG_INTERSECTION_BUBBLES;
extern uint64 LOG_INTERSECTION_BUBBLES_DEBUG;
extern uint64 LOG_INTERSECTION_JOINING;
extern uint64 LOG_INTERSECTION_JOINING_DEBUG;
extern uint64 LOG_INITIAL_CONTAINED_PLACEMENT;
extern uint64 LOG_HAPPINESS;
extern uint64 LOG_INTERMEDIATE_UNITIGS;
extern uint64 LOG_MATE_SPLIT_ANALYSIS;
extern uint64 LOG_MATE_SPLIT_DISCONTINUOUS;
extern uint64 LOG_MATE_SPLIT_UNHAPPY_CONTAINS;
extern uint64 LOG_MATE_SPLIT_COVERAGE_PLOT;
extern uint64 LOG_SET_PARENT_AND_HANG;
extern uint64 LOG_STDERR;

extern uint64 LOG_PLACE_FRAG;

extern char const *logFileFlagNames[64];

#endif  //  INCLUDE_AS_BAT_LOGGING

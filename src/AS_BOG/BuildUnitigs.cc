
/**************************************************************************
 * Copyright (C) 2005, J Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id: BuildUnitigs.cc,v 1.74 2010-09-28 09:17:54 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_MateChecker.hh"

//  Used in AS_BOG_Unitig.cc
FragmentInfo     *debugfi = 0L;
BestOverlapGraph *bog     = 0L;

FILE             *logFile      = stderr;
uint32            logFileOrder = 1;
uint64            logFileFlags = 0;  //  defined in AS_BOG_Datatypes.hh

uint64 LOG_OVERLAP_QUALITY             = 0x0000000000000001;  //  Debug, scoring of overlaps
uint64 LOG_CHUNK_GRAPH                 = 0x0000000000000002;  //  Report the chunk graph as we build it
uint64 LOG_INTERSECTIONS               = 0x0000000000000004;  //  Report intersections found when building initial unitigs
uint64 LOG_POPULATE_UNITIG             = 0x0000000000000008;  //  Report building of initial unitigs (both unitig creation and fragment placement)
uint64 LOG_INTERSECTION_BREAKING       = 0x0000000000000010;  //  
uint64 LOG_INTERSECTION_BUBBLES        = 0x0000000000000020;  //  
uint64 LOG_INTERSECTION_BUBBLES_DEBUG  = 0x0000000000000040;  //  
uint64 LOG_INTERSECTION_JOINING        = 0x0000000000000080;  //
uint64 LOG_INTERSECTION_JOINING_DEBUG  = 0x0000000000000100;  //
uint64 LOG_INITIAL_CONTAINED_PLACEMENT = 0x0000000000000200;  //
uint64 LOG_HAPPINESS                   = 0x0000000000000400;  //

const char *logFileFlagNames[64] = { "overlapQuality",
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
                                     NULL
};

//  Closes the current logFile, opens a new one called 'prefix.logFileOrder.name'.  If 'name' is
//  NULL, the logFile is reset to stderr.
void
setLogFile(char *prefix, char *name) {
  char  logFileName[FILENAME_MAX];

  if (logFile != stderr)
    fclose(logFile);

  if (name == NULL) {
    logFile = stderr;
    return;
  }

  sprintf(logFileName, "%s.%03u.%s.log", prefix, logFileOrder++, name);

  errno = 0;
  logFile = fopen(logFileName, "w");
  if (errno) {
    fprintf(stderr, "setLogFile()-- Failed to open logFile '%s': %s.\n", logFileName, strerror(errno));
    fprintf(stderr, "setLogFile()-- Will now log to stderr instead.\n");
    logFile = stderr;
  }

  fprintf(stderr,  "setLogFile()-- Now logging to '%s'\n", logFileName);
  //fprintf(logFile, "setLogFile()-- Search for '%s' in unitigger.err for messages before/after this file.\n", logFileName);
}



int
main (int argc, char * argv []) {
  char      *gkpStorePath           = NULL;
  char      *ovlStoreUniqPath       = NULL;
  char      *ovlStoreReptPath       = NULL;
  char      *tigStorePath           = NULL;

  double    erate                   = 0.015;
  double    elimit                  = 0.0;
  long      genome_size             = 0;

  int       fragment_count_target   = 0;
  char     *output_prefix           = NULL;

  bool      popBubbles              = false;
  bool      breakIntersections      = false;
  bool      joinUnitigs             = false;
  int       badMateBreakThreshold   = -7;

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-B") == 0) {
      fragment_count_target = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      output_prefix = argv[++arg];

    } else if (strcmp(argv[arg], "-G") == 0) {
      gkpStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      if      (ovlStoreUniqPath == NULL)
        ovlStoreUniqPath = argv[++arg];
      else if (ovlStoreReptPath == NULL)
        ovlStoreReptPath = argv[++arg];
      else
        err++;

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-U") == 0) {
      popBubbles = true;

    } else if (strcmp(argv[arg], "-b") == 0) {
      breakIntersections = true;

    } else if (strcmp(argv[arg], "-J") == 0) {
      joinUnitigs = true;

    } else if (strcmp(argv[arg], "-e") == 0) {
      erate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      elimit = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-m") == 0) {
      badMateBreakThreshold = -atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-s") == 0) {
      genome_size = atol(argv[++arg]);

    } else if (strcmp(argv[arg], "-D") == 0) {
      uint32  opt = 0;
      uint64  flg = 1;
      bool    fnd = false;
      for (arg++; logFileFlagNames[opt]; flg <<= 1, opt++) {
        if (strcasecmp(logFileFlagNames[opt], argv[arg]) == 0) {
          fprintf(stderr, "SET flag to "F_U64"\n", flg);
          logFileFlags |= flg;
          fnd = true;
        }
      }
      if (strcasecmp("all", argv[arg]) == 0) {
        logFileFlags = 0xffffffffffffffff;
        fnd = true;
      }
      if (fnd == false) {
        fprintf(stderr, "ERROR:  Unknown '-D' option '%s'.\n", argv[arg]);
        err++;
      }

    } else {
      err++;
    }

    arg++;
  }

  if ((erate < 0.0) || (AS_MAX_ERROR_RATE < erate))
    err++;
  if (elimit < 0.0)
    err++;
  if (output_prefix == NULL)
    err++;
  if (gkpStorePath == NULL)
    err++;
  if (ovlStoreUniqPath == NULL)
    err++;
  if (tigStorePath == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -o outputName -O ovlStore -G gkpStore -T tigStore\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -O         Mandatory path to an ovlStore.\n");
    fprintf(stderr, "  -G         Mandatory path to a gkpStore.\n");
    fprintf(stderr, "  -T         Mandatory path to a tigStore (can exist or not).\n");
    fprintf(stderr, "  -o prefix  Mandatory name for the output files\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -B b       Target number of fragments per tigStore (consensus) partition\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -s size    If the genome size is set to 0, this will cause the unitigger\n");
    fprintf(stderr, "             to try to estimate the genome size based on the constructed\n");
    fprintf(stderr, "             unitig lengths.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -U         Enable EXPERIMENTAL short unitig merging (aka bubble popping).\n");
    fprintf(stderr, "  -J         Enable EXPERIMENTAL long unitig joining.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -b         Break promisciuous unitigs at unitig intersection points\n");
    fprintf(stderr, "  -m 7       Break a unitig if a region has more than 7 bad mates\n");
    fprintf(stderr, " \n");
    fprintf(stderr, "Overlap Selection - an overlap will be considered for use in a unitig if either of\n");
    fprintf(stderr, "                    the following conditions hold:\n");
    fprintf(stderr, "  -e 0.015   no more than 0.015 fraction (1.5%%) error\n");
    fprintf(stderr, "  -E 0       no more than 0 errors\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Debugging and Logging\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -D <name>  enable logging/debugging for a specific component.\n");
    for (uint32 l=0; logFileFlagNames[l]; l++)
      fprintf(stderr, "               %s\n", logFileFlagNames[l]);
    fprintf(stderr, "\n");

    if ((erate < 0.0) || (AS_MAX_ERROR_RATE < erate))
      fprintf(stderr, "Invalid overlap error threshold (-e option); must be between 0.00 and %.2f.\n",
              AS_MAX_ERROR_RATE);

    if (elimit < 0.0)
      fprintf(stderr, "Invalid overlap error limit (-E option); must be above 0.\n");

    if (output_prefix == NULL)
      fprintf(stderr, "No output prefix name (-o option) supplied.\n");

    if (gkpStorePath == NULL)
      fprintf(stderr, "No gatekeeper store (-G option) supplied.\n");

    if (ovlStoreUniqPath == NULL)
      fprintf(stderr, "No overlap store (-O option) supplied.\n");

    if ((ovlStoreUniqPath != NULL) && (ovlStoreUniqPath == ovlStoreReptPath))
      fprintf(stderr, "Too many overlap stores (-O option) supplied.\n");

    if (tigStorePath == NULL)
      fprintf(stderr, "No output tigStore (-T option) supplied.\n");

    exit(1);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Bubble popping        = %s\n", (popBubbles) ? "on" : "off");
  fprintf(stderr, "Intersection breaking = %s\n", (breakIntersections) ? "on" : "off");
  fprintf(stderr, "Bad mate threshold    = %d\n", badMateBreakThreshold);
  fprintf(stderr, "Error threshold       = %.3f (%.3f%%)\n", erate, erate * 100);
  fprintf(stderr, "Error limit           = %.3f errors\n", elimit);
  fprintf(stderr, "Genome Size           = "F_S64"\n", genome_size);
  fprintf(stderr, "\n");
  fprintf(stderr, "sizeof(DoveTailNode)  = %d\n", (int)sizeof(DoveTailNode));
  fprintf(stderr, "\n");

  gkStore          *gkpStore     = new gkStore(gkpStorePath, FALSE, FALSE);
  OverlapStore     *ovlStoreUniq = AS_OVS_openOverlapStore(ovlStoreUniqPath);
  OverlapStore     *ovlStoreRept = ovlStoreReptPath ? AS_OVS_openOverlapStore(ovlStoreReptPath) : NULL;

  FragmentInfo     *fragInfo = new FragmentInfo(gkpStore);

  debugfi = fragInfo;

  BestOverlapGraph      *BOG = new BestOverlapGraph(fragInfo, ovlStoreUniq, ovlStoreRept, erate, elimit);

  bog = BOG;

  ChunkGraph *cg = new ChunkGraph(fragInfo, BOG);
  UnitigGraph utg(fragInfo, BOG);
  utg.build(cg, ovlStoreUniq, ovlStoreRept, breakIntersections, joinUnitigs, popBubbles, output_prefix);

  MateChecker  mateChecker(fragInfo);
  mateChecker.checkUnitigGraph(utg, badMateBreakThreshold);

  setLogFile("unitigger", "output");

  utg.setParentAndHang(cg);

  float globalARate = utg.getGlobalArrivalRate(gkpStore->gkStore_getNumRandomFragments(), genome_size);
  Unitig::setGlobalArrivalRate(globalARate);

  utg.writeIUMtoFile(output_prefix, tigStorePath, fragment_count_target);
  utg.writeOVLtoFile(output_prefix);
  utg.writeCGAtoFile(output_prefix, globalARate);

  delete BOG;
  delete cg;
  delete fragInfo;

  delete gkpStore;
  AS_OVS_closeOverlapStore(ovlStoreUniq);
  AS_OVS_closeOverlapStore(ovlStoreRept);

  setLogFile(NULL, NULL);
  fprintf(stderr, "Bye.\n");

  return(0);
}

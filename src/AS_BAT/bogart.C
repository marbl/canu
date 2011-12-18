
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

const char *mainid = "$Id: bogart.C,v 1.12 2011-12-18 08:14:34 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_ChunkGraph.H"
#include "AS_BAT_Unitig.H"

#include "AS_BAT_OverlapCache.H"

#include "AS_BAT_InsertSizes.H"

#include "AS_BAT_PopulateUnitig.H"
#include "AS_BAT_Instrumentation.H"
#include "AS_BAT_EvaluateMates.H"
#include "AS_BAT_PlaceContains.H"
#include "AS_BAT_PlaceZombies.H"

#include "AS_BAT_MergeSplitJoin.H"
#include "AS_BAT_SplitDiscontinuous.H"

#include "AS_BAT_SetParentAndHang.H"
#include "AS_BAT_ArrivalRate.H"
#include "AS_BAT_Outputs.H"


FragmentInfo     *FI  = 0L;
OverlapCache     *OC  = 0L;
BestOverlapGraph *OG  = 0L;
ChunkGraph       *CG  = 0L;
InsertSizes      *IS  = 0L;


FILE             *logFile      = stderr;
uint32            logFileOrder = 0;
uint64            logFileFlags = 0;  //  defined in AS_BAT_Datatypes.hh

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

const char *logFileFlagNames[64] = { "overlapQuality",
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
setLogFile(const char *prefix, const char *name) {
  char  logFileName[FILENAME_MAX];

  if (logFileFlagSet(LOG_STDERR))
    //  Write everything to stderr
    return;

  if (logFile != stderr)
    fclose(logFile);

  if (name == NULL) {
    logFile = stderr;
    return;
  }

  sprintf(logFileName, "%s.%03u.%s.log", prefix, ++logFileOrder, name);

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

  double    erateGraph              = 0.020;
  double    elimitGraph             = 2.0;
  double    erateMerge              = 0.045;
  double    elimitMerge             = 4.0;

  uint64    ovlCacheMemory          = UINT64_MAX;
  uint32    ovlCacheLimit           = UINT32_MAX;
  uint64    genomeSize              = 0;

  int       fragment_count_target   = 0;
  char     *output_prefix           = NULL;

  bool      shatterRepeats          = false;

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

    } else if (strcmp(argv[arg], "-s") == 0) {
      genomeSize = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-R") == 0) {
      shatterRepeats = true;

    } else if (strcmp(argv[arg], "-eg") == 0) {
      erateGraph = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-Eg") == 0) {
      elimitGraph = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-em") == 0) {
      erateMerge = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-Em") == 0) {
      elimitMerge = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-M") == 0) {
      ovlCacheMemory  = (uint64)(atof(argv[++arg]) * 1024 * 1024 * 1024);

    } else if (strcmp(argv[arg], "-N") == 0) {
      ovlCacheLimit   = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-D") == 0) {
      uint32  opt = 0;
      uint64  flg = 1;
      bool    fnd = false;
      for (arg++; logFileFlagNames[opt]; flg <<= 1, opt++) {
        if (strcasecmp(logFileFlagNames[opt], argv[arg]) == 0) {
          logFileFlags |= flg;
          fnd = true;
        }
      }
      if (strcasecmp("all", argv[arg]) == 0) {
        for (flg=1, opt=0; logFileFlagNames[opt]; flg <<= 1, opt++)
          if (strcasecmp(logFileFlagNames[opt], "stderr") != 0)
            logFileFlags |= flg;
        fnd = true;
      }
      if (strcasecmp("most", argv[arg]) == 0) {
        for (flg=1, opt=0; logFileFlagNames[opt]; flg <<= 1, opt++)
          if ((strcasecmp(logFileFlagNames[opt], "stderr") != 0) &&
              (strcasecmp(logFileFlagNames[opt], "mateSplitCoveragePlot") != 0) &&
              (strcasecmp(logFileFlagNames[opt], "overlapQuality") != 0) &&
              (strcasecmp(logFileFlagNames[opt], "happiness") != 0) &&
              (strcasecmp(logFileFlagNames[opt], "setParentAndHang") != 0))
            logFileFlags |= flg;
        fnd = true;
      }
      if (fnd == false) {
        fprintf(stderr, "ERROR:  Unknown '-D' option '%s'.\n", argv[arg]);
        err++;
      }

    } else if (strcmp(argv[arg], "-d") == 0) {
      uint32  opt = 0;
      uint64  flg = 1;
      bool    fnd = false;
      for (arg++; logFileFlagNames[opt]; flg <<= 1, opt++) {
        if (strcasecmp(logFileFlagNames[opt], argv[arg]) == 0) {
          logFileFlags &= ~flg;
          fnd = true;
        }
      }
      if (fnd == false) {
        fprintf(stderr, "ERROR:  Unknown '-d' option '%s'.\n", argv[arg]);
        err++;
      }

    } else {
      err++;
    }

    arg++;
  }

  if ((erateGraph < 0.0) || (AS_MAX_ERROR_RATE < erateGraph))
    err++;
  if (elimitGraph < 0.0)
    err++;
  if ((erateMerge < 0.0) || (AS_MAX_ERROR_RATE < erateMerge))
    err++;
  if (elimitMerge < 0.0)
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
    fprintf(stderr, "Algorithm Options\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -R         Shatter and rebuild repeat unitigs\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Overlap Selection - an overlap will be considered for use in a unitig if either of\n");
    fprintf(stderr, "                    the following conditions hold:\n");
    fprintf(stderr, "  When constructing the Best Overlap Graph and Promiscuous Unitigs ('g'raph):\n");
    fprintf(stderr, "    -eg 0.020   no more than 0.020 fraction (2.0%%) error\n");
    fprintf(stderr, "    -Eg 2       no more than 2 errors (useful with short reads)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  When popping bubbles and splitting repeat/unique junctions ('m'erging):\n");
    fprintf(stderr, "    -em 0.045   no more than 0.045 fraction (4.5%%) error when bubble popping and repeat splitting\n");
    fprintf(stderr, "    -Em 4       no more than r errors (useful with short reads)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Overlap Storage\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -M gb    Use at most 'gb' gigabytes of memory for storing overlaps.\n");
    fprintf(stderr, "    -N num   Load at most 'num' overlaps per read.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Debugging and Logging\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -D <name>  enable logging/debugging for a specific component.\n");
    for (uint32 l=0; logFileFlagNames[l]; l++)
      fprintf(stderr, "               %s\n", logFileFlagNames[l]);
    fprintf(stderr, "\n");

    if ((erateGraph < 0.0) || (AS_MAX_ERROR_RATE < erateGraph))
      fprintf(stderr, "Invalid overlap error threshold (-eg option); must be between 0.00 and %.2f.\n",
              AS_MAX_ERROR_RATE);

    if (elimitGraph < 0.0)
      fprintf(stderr, "Invalid overlap error limit (-Eg option); must be above 0.\n");

    if ((erateMerge < 0.0) || (AS_MAX_ERROR_RATE < erateMerge))
      fprintf(stderr, "Invalid overlap error threshold (-em option); must be between 0.00 and %.2f.\n",
              AS_MAX_ERROR_RATE);

    if (elimitMerge < 0.0)
      fprintf(stderr, "Invalid overlap error limit (-Em option); must be above 0.\n");

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
  fprintf(stderr, "Graph error threshold = %.3f (%.3f%%)\n", erateGraph, erateGraph * 100);
  fprintf(stderr, "Graph error limit     = %.3f errors\n", elimitGraph);
  fprintf(stderr, "Merge error threshold = %.3f (%.3f%%)\n", erateMerge, erateMerge * 100);
  fprintf(stderr, "Merge error limit     = %.3f errors\n", elimitMerge);
  fprintf(stderr, "Genome Size           = "F_U64"\n", genomeSize);
  fprintf(stderr, "\n");

  for (uint64 i=0, j=1; i<64; i++, j<<=1)
    if (logFileFlagSet(j))
      fprintf(stderr, "DEBUG                 = %s\n", logFileFlagNames[i]);

  gkStore          *gkpStore     = new gkStore(gkpStorePath, FALSE, FALSE);
  OverlapStore     *ovlStoreUniq = AS_OVS_openOverlapStore(ovlStoreUniqPath);
  OverlapStore     *ovlStoreRept = ovlStoreReptPath ? AS_OVS_openOverlapStore(ovlStoreReptPath) : NULL;

  UnitigVector      unitigs;

  FI = new FragmentInfo(gkpStore, output_prefix);

  // Initialize where we've been to nowhere, and create the non-existent 0th unitig.
  Unitig::resetFragUnitigMap(FI->numFragments());
  unitigs.push_back(NULL);

  OC = new OverlapCache(ovlStoreUniq, ovlStoreRept, MAX(erateGraph, erateMerge), MAX(elimitGraph, elimitMerge), ovlCacheMemory, ovlCacheLimit);
  OG = new BestOverlapGraph(erateGraph, elimitGraph, output_prefix);
  CG = new ChunkGraph(output_prefix);
  IS = NULL;

  AS_OVS_closeOverlapStore(ovlStoreUniq);  ovlStoreUniq = NULL;
  AS_OVS_closeOverlapStore(ovlStoreRept);  ovlStoreRept = NULL;




  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Build the initial unitig path from non-contained fragments.  The first pass is usually the
  //  only one needed, but occasionally (maybe) we miss fragments, so we make an explicit pass
  //  through all fragments and place whatever isn't already placed.
  //

  setLogFile(output_prefix, "buildUnitigs");
  fprintf(logFile, "==> BUILDING UNITIGS from %d fragments.\n", FI->numFragments());

  for (uint32 fi=CG->nextFragByChunkLength(); fi>0; fi=CG->nextFragByChunkLength())
    populateUnitig(unitigs, fi);

  //setLogFile(output_prefix, "buildUnitigs-MissedFragments");
  fprintf(logFile, "==> BUILDING UNITIGS catching missed fragments.\n");

  for (uint32 fi=1; fi <= FI->numFragments(); fi++)
    populateUnitig(unitigs, fi);

  reportOverlapsUsed(unitigs, output_prefix, "buildUnitigs");
  reportUnitigs(unitigs, output_prefix, "buildUnitigs");
  evaluateMates(unitigs, output_prefix, "buildUnitigs");


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Place contained fragments.  There is a variety of evidence we can use to determine
  //  where to place these:
  //
  //    1) By using the best containment overlap (no way to measure if this is an uncontested
  //    placement)
  //
  //    2) By using all overlaps (we can measure if the placement is unique)
  //
  //    3) By using any mate relationship to a non-contained fragment
  //
  //    4) By using any mate relationship to an unambiguously placed (#2) contained fragment
  //
  //  For now, we're using the standard method #1.
  //
  //
  //  After placing contains, every fragment should be in a unitig.  Sometimes this
  //  is not true and we have zombie fragments (not dead, but not in a unitig).  We
  //  place any of those into their own unitig.
  //

  setLogFile(output_prefix, "placeContainsZombies");

  placeContainsUsingBestOverlaps(unitigs);
  //placeContainsUsingAllOverlaps(bool withMatesToNonContained,
  //                              bool withMatesToUnambiguousContain);

  placeZombies(unitigs, erateMerge, elimitMerge);

  checkUnitigMembership(unitigs);
  reportOverlapsUsed(unitigs, output_prefix, "placeContainsZombies");
  reportUnitigs(unitigs, output_prefix, "placeContainsZombies");
  evaluateMates(unitigs, output_prefix, "placeContainsZombies");

  setLogFile(output_prefix, "mergeSplitJoin");

  mergeSplitJoin(unitigs, shatterRepeats);

  //  rebuildRepeatUnitigs()
  if (shatterRepeats) {
    BestOverlapGraph  *OGsave = OG;
    ChunkGraph        *CGsave = CG;

    OG = new BestOverlapGraph(erateGraph, elimitGraph, output_prefix);
    CG = new ChunkGraph(output_prefix);

    setLogFile(output_prefix, "buildRepeatUnitigs");
    fprintf(logFile, "==> BUILDING REPEAT UNITIGS from %d fragments.\n", FI->numFragments());

    for (uint32 fi=CG->nextFragByChunkLength(); fi>0; fi=CG->nextFragByChunkLength())
      populateUnitig(unitigs, fi);

    fprintf(logFile, "==> BUILDING REPEAT UNITIGS catching missed fragments.\n");

    for (uint32 fi=1; fi <= FI->numFragments(); fi++)
      populateUnitig(unitigs, fi);

    delete OG;
    delete CG;

    OG = OGsave;
    CG = CGsave;
  }

  placeContainsUsingBestOverlaps(unitigs);

  splitDiscontinuousUnitigs(unitigs);       //  Clean up splitting problems.
  placeContainsUsingBestOverlaps(unitigs);

  checkUnitigMembership(unitigs);
  reportOverlapsUsed(unitigs, output_prefix, "mergeSplitJoin");
  reportUnitigs(unitigs, output_prefix, "mergeSplitJoin");
  evaluateMates(unitigs, output_prefix, "mergeSplitJoin");

  //  OUTPUT

  setLogFile(output_prefix, "setParentAndHang");
  setParentAndHang(unitigs);

  setLogFile(output_prefix, "output");

  float globalARate = getGlobalArrivalRate(unitigs, gkpStore->gkStore_getNumRandomFragments(), genomeSize);
  Unitig::setGlobalArrivalRate(globalARate);

  writeIUMtoFile(unitigs, output_prefix, tigStorePath, fragment_count_target);
  writeOVLtoFile(unitigs, output_prefix);
  writeCGAtoFile(unitigs, output_prefix, globalARate);

  delete IS;
  delete CG;
  delete OG;
  delete OC;
  delete FI;

  delete gkpStore;

  for (uint32  ti=0; ti<unitigs.size(); ti++)
    delete unitigs[ti];

  setLogFile(NULL, NULL);
  fprintf(logFile, "Bye.\n");

  return(0);
}

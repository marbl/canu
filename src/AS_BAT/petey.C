
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

const char *mainid = "$Id: petey.C,v 1.11 2012-07-30 01:21:01 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_ChunkGraph.H"
#include "AS_BAT_Unitig.H"

#include "AS_BAT_InsertSizes.H"
#include "AS_BAT_MateLocation.H"
#include "AS_BAT_EvaluateMates.H"
#include "AS_BAT_Instrumentation.H"
#include "AS_BAT_Breaking.H"
#include "AS_BAT_PlaceContains.H"
#include "AS_BAT_MoveContains.H"
#include "AS_BAT_SplitDiscontinuous.H"
#include "AS_BAT_MateChecker.H"
#include "AS_BAT_SetParentAndHang.H"
#include "AS_BAT_Outputs.H"


FragmentInfo     *FI  = 0L;
OverlapCache     *OC  = 0L;
BestOverlapGraph *OG  = 0L;
ChunkGraph       *CG  = 0L;
InsertSizes      *IS  = 0L;


int
main (int argc, char * argv []) {
  char      *gkpStorePath           = NULL;
  char      *ovlStoreUniqPath       = NULL;
  char      *ovlStoreReptPath       = NULL;
  char      *tigStorePath           = NULL;

  double    erate                   = 0.015;
  double    elimit                  = 0.0;

  int       fragment_count_target   = 0;
  char     *output_prefix           = NULL;

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

    } else if (strcmp(argv[arg], "-e") == 0) {
      erate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      elimit = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-m") == 0) {
      badMateBreakThreshold = -atoi(argv[++arg]);

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
              (strcasecmp(logFileFlagNames[opt], "overlapQuality") != 0))
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
  fprintf(stderr, "Bad mate threshold    = %d\n", badMateBreakThreshold);
  fprintf(stderr, "Error threshold       = %.3f (%.3f%%)\n", erate, erate * 100);
  fprintf(stderr, "Error limit           = %.3f errors\n", elimit);
  fprintf(stderr, "\n");

  for (uint64 i=0, j=1; i<64; i++, j<<=1)
    if (logFileFlagSet(j))
      fprintf(stderr, "DEBUG                 = %s\n", logFileFlagNames[i]);

  gkStore          *gkpStore     = new gkStore(gkpStorePath, FALSE, FALSE);
  OverlapStore     *ovlStoreUniq = AS_OVS_openOverlapStore(ovlStoreUniqPath);
  OverlapStore     *ovlStoreRept = ovlStoreReptPath ? AS_OVS_openOverlapStore(ovlStoreReptPath) : NULL;

  FI = new FragmentInfo(gkpStore, output_prefix);
  OC = new OverlapCache(ovlStoreUniq, ovlStoreRept, output_prefix, 6.0, 1.5, 0, 0, false, false);
  OG = new BestOverlapGraph(erate, elimit, output_prefix, false, false);
  CG = new ChunkGraph(output_prefix);
  IS = NULL;

  UnitigVector   unitigs;







  ////////////////////////////////////////////////////////////////////////////////
  //
  //  This is the MateChecker
  //
  ////////////////////////////////////////////////////////////////////////////////

  //  Move unhappy contained fragments before attempting mate based
  //  splitting.  We know they're broken already, and if scaffolder
  //  can put them back, we'll let it.

  setLogFile(output_prefix, "moveContains1");
  writeLog("==> MOVE CONTAINS #1\n");
  moveContains(unitigs);

  reportOverlapsUsed(unitigs, output_prefix, "moveContains1");
  checkUnitigMembership(unitigs);
  evaluateMates(unitigs, output_prefix, "moveContains1");

  //  This should do absolutely nothing.  If it does, something is
  //  broken.  By just ejecting unhappy contains, nothing should be
  //  disconnected.

  setLogFile(output_prefix, "splitDiscontinuous1");
  writeLog("==> SPLIT DISCONTINUOUS #1\n");
  splitDiscontinuousUnitigs(unitigs);

  reportOverlapsUsed(unitigs, output_prefix, "splitDiscontinuous1");
  checkUnitigMembership(unitigs);
  evaluateMates(unitigs, output_prefix, "splitDiscontinuous1");

  setLogFile(output_prefix, "splitBadMates");
  writeLog("==> SPLIT BAD MATES\n");

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];

    if ((tig == NULL) || (tig->getNumFrags() < 2))
      continue;

    MateLocation       *mates  = new MateLocation(unitigs, tig);
    vector<breakPoint> *breaks = computeMateCoverage(tig, mates, badMateBreakThreshold);

    if (breakUnitigAt(unitigs, tig, *breaks, true) == true)
      if (logFileFlagSet(LOG_MATE_SPLIT_COVERAGE_PLOT))
        mates->dumpHappiness(output_prefix, "splitBadMates");

    delete breaks;
    delete mates;
  }

  placeContainsUsingBestOverlaps(unitigs);

  reportOverlapsUsed(unitigs, output_prefix, "splitBadMates");
  checkUnitigMembership(unitigs);
  evaluateMates(unitigs, output_prefix, "splitBadMates");

  //  The splitting code above is not smart enough to move contained
  //  fragments with containers.  This leaves unitigs disconnected.
  //  We break those unitigs here.

  setLogFile(output_prefix, "splitDiscontinuous2");
  writeLog("==> SPLIT DISCONTINUOUS #2\n");
  splitDiscontinuousUnitigs(unitigs);

  reportOverlapsUsed(unitigs, output_prefix, "splitDiscontinuous2");
  checkUnitigMembership(unitigs);
  evaluateMates(unitigs, output_prefix, "splitDiscontinuous2");

  //  But now, all the splitting probably screwed up happiness of
  //  contained fragments, of left some unhappy fragments in a unitig
  //  that just lost the container.

  setLogFile(output_prefix, "moveContains2");
  writeLog("==> MOVE CONTAINS #2\n");

  moveContains(unitigs);

  reportOverlapsUsed(unitigs, output_prefix, "moveContains2");
  checkUnitigMembership(unitigs);
  evaluateMates(unitigs, output_prefix, "moveContains2");

  //  Do one last check for disconnected unitigs.

  setLogFile(output_prefix, "splitDiscontinuous3");
  writeLog("==> SPLIT DISCONTINUOUS #3\n");
  splitDiscontinuousUnitigs(unitigs);

  reportOverlapsUsed(unitigs, output_prefix, "splitDiscontinuous3");
  checkUnitigMembership(unitigs);
  evaluateMates(unitigs, output_prefix, "splitDiscontinuous3");

  //  OUTPUT

  setLogFile(output_prefix, "setParentAndHang");
  setParentAndHang(unitigs);

  setLogFile(output_prefix, "output");

  writeIUMtoFile(unitigs, output_prefix, tigStorePath, fragment_count_target);
  writeOVLtoFile(unitigs, output_prefix);

  delete IS;
  delete CG;
  delete OG;
  delete FI;

  delete gkpStore;
  AS_OVS_closeOverlapStore(ovlStoreUniq);
  AS_OVS_closeOverlapStore(ovlStoreRept);

  for (uint32  ti=0; ti<unitigs.size(); ti++)
    delete unitigs[ti];

  setLogFile(NULL, NULL);
  writeLog("Bye.\n");

  return(0);
}

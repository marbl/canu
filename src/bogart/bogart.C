
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
 *    src/AS_BAT/bogart.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2014-JAN-29
 *      are Copyright 2010-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-21 to 2015-AUG-07
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-MAR-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_ChunkGraph.H"

#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_PopulateUnitig.H"
#include "AS_BAT_Instrumentation.H"
#include "AS_BAT_PlaceContains.H"

#include "AS_BAT_MergeUnitigs.H"
#include "AS_BAT_PopBubbles.H"
#include "AS_BAT_MarkRepeatReads.H"

#include "AS_BAT_SplitDiscontinuous.H"

#include "AS_BAT_PromoteToSingleton.H"

#include "AS_BAT_SetParentAndHang.H"
#include "AS_BAT_Outputs.H"


FragmentInfo     *FI  = 0L;
OverlapCache     *OC  = 0L;
BestOverlapGraph *OG  = 0L;
ChunkGraph       *CG  = 0L;

int
main (int argc, char * argv []) {
  char      *gkpStorePath            = NULL;
  char      *ovlStoreUniqPath        = NULL;
  char      *ovlStoreReptPath        = NULL;
  char      *tigStorePath            = NULL;

  double    erateGraph               = 0.030;
  double    erateMax                 = 0.050;

  uint64    genomeSize               = 0;

  uint32    fewReadsNumber           = 2;      //  Parameters for labeling of unassembled; also set in pipelines/canu/Defaults.pm
  uint32    tooShortLength           = 1000;
  double    spanFraction             = 0.75;
  double    lowcovFraction           = 0.75;
  uint32    lowcovDepth              = 2;

  double    deviationGraph           = 5.0;
  double    deviationBubble          = 5.0;
  double    deviationRepeat          = 3.0;

  uint32    confusedAbsolute         = 250;
  double    confusedPercent          = 100.0;

  int32     numThreads               = 0;

  uint64    ovlCacheMemory           = UINT64_MAX;
  uint32    ovlCacheLimit            = UINT32_MAX;

  bool      onlySave                 = false;
  bool      doSave                   = false;

  int       fragment_count_target    = 0;
  char     *prefix                   = NULL;

  bool      enableJoining            = false;

  bool      enableShatterRepeats     = false;
  bool      enableReconstructRepeats = false;

  uint32    minReadLen               = 0;
  uint32    minOverlap               = 500;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-B") == 0) {
      fragment_count_target = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      prefix = argv[++arg];

    } else if (strcmp(argv[arg], "-G") == 0) {
      gkpStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      if      (ovlStoreUniqPath == NULL)
        ovlStoreUniqPath = argv[++arg];
      else if (ovlStoreReptPath == NULL)
        ovlStoreReptPath = argv[++arg];
      else
        err.push_back(NULL);

    } else if (strcmp(argv[arg], "-gs") == 0) {
      genomeSize = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-J") == 0) {
      enableJoining = true;

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-SR") == 0) {
      enableShatterRepeats     = true;

    } else if (strcmp(argv[arg], "-R") == 0) {
      enableShatterRepeats     = true;
      enableReconstructRepeats = true;

    } else if (strcmp(argv[arg], "-unassembled") == 0) {
      fewReadsNumber  = atoi(argv[++arg]);
      tooShortLength  = atoi(argv[++arg]);
      spanFraction    = atof(argv[++arg]);
      lowcovFraction  = atof(argv[++arg]);
      lowcovDepth     = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-RL") == 0) {
      minReadLen = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-threads") == 0) {
      numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-eg") == 0) {
      erateGraph = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-eM") == 0) {
      erateMax = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-el") == 0) {
      minOverlap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-ca") == 0) {  //  Edge confused, based on absolute difference
      confusedAbsolute = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-cp") == 0) {  //  Edge confused, based on percent difference
      confusedPercent = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-dg") == 0) {  //  Deviations, graph
      deviationGraph = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-db") == 0) {  //  Deviations, bubble
      deviationBubble = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-dr") == 0) {  //  Deviations, repeat
      deviationRepeat = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-M") == 0) {
      ovlCacheMemory  = (uint64)(atof(argv[++arg]) * 1024 * 1024 * 1024);

    } else if (strcmp(argv[arg], "-N") == 0) {
      ovlCacheLimit   = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-create") == 0) {
      onlySave = true;
      doSave   = true;

    } else if (strcmp(argv[arg], "-save") == 0) {
      doSave = true;

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
              (strcasecmp(logFileFlagNames[opt], "overlapScoring") != 0) &&
              (strcasecmp(logFileFlagNames[opt], "chunkGraph") != 0) &&
              (strcasecmp(logFileFlagNames[opt], "setParentAndHang") != 0))
            logFileFlags |= flg;
        fnd = true;
      }
      if (fnd == false) {
        char *s = new char [1024];
        sprintf(s, "Unknown '-D' option '%s'.\n", argv[arg]);
        err.push_back(s);
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
        char *s = new char [1024];
        sprintf(s, "Unknown '-d' option '%s'.\n", argv[arg]);
        err.push_back(s);
      }

    } else {
      char *s = new char [1024];
      sprintf(s, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (erateGraph < 0.0)
    err.push_back(NULL);
  if (erateMax < 0.0)
    err.push_back(NULL);

  if (prefix == NULL)
    err.push_back(NULL);
  if (gkpStorePath == NULL)
    err.push_back(NULL);
  if (ovlStoreUniqPath == NULL)
    err.push_back(NULL);
  if (tigStorePath == NULL)
    err.push_back(NULL);

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -o outputName -O ovlStore -G gkpStore -T tigStore\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -O         Mandatory path to an ovlStore.\n");
    fprintf(stderr, "  -G         Mandatory path to a gkpStore.\n");
    fprintf(stderr, "  -T         Mandatory path to a tigStore (can exist or not).\n");
    fprintf(stderr, "  -o prefix  Mandatory name for the output files\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -B b       Target number of fragments per tigStore (consensus) partition\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Algorithm Options\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -gs        Genome size in bases.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -J         Join promiscuous unitigs using unused best edges.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -SR        Shatter repeats, don't rebuild.\n");
    fprintf(stderr, "  -R         Shatter repeats (-SR), then rebuild them\n");
    fprintf(stderr, "  -RL len    Force reads below 'len' bases to be singletons.\n");
    fprintf(stderr, "               This WILL cause CGW to fail; diagnostic only.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -threads N Use N compute threads during repeat detection.\n");
    fprintf(stderr, "               0 - use OpenMP default (default)\n");
    fprintf(stderr, "               1 - use one thread\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Overlap Selection - an overlap will be considered for use in a unitig under\n");
    fprintf(stderr, "                    the following conditions:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  When constructing the Best Overlap Graph and Promiscuous Unitigs ('g'raph):\n");
    fprintf(stderr, "    -eg 0.020   no more than 0.020 fraction (2.0%%) error   ** DEPRECATED **\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  When loading overlaps, an inflated maximum (to allow reruns with different error rates):\n");
    fprintf(stderr, "    -eM 0.05   no more than 0.05 fraction (5.0%%) error in any overlap loaded into bogart\n");
    fprintf(stderr, "               the maximum used will ALWAYS be at leeast the maximum of the four error rates\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  For all, the lower limit on overlap length\n");
    fprintf(stderr, "    -el 500     no shorter than 40 bases\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Overlap Storage\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -M gb    Use at most 'gb' gigabytes of memory for storing overlaps.\n");
    fprintf(stderr, "    -N num   Load at most 'num' overlaps per read.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -create  Only create the overlap graph, save to disk and quit.\n");
    fprintf(stderr, "    -save    Save the overlap graph to disk, and continue.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Debugging and Logging\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -D <name>  enable logging/debugging for a specific component.\n");
    fprintf(stderr, "  -d <name>  disable logging/debugging for a specific component.\n");
    for (uint32 l=0; logFileFlagNames[l]; l++)
      fprintf(stderr, "               %s\n", logFileFlagNames[l]);
    fprintf(stderr, "\n");

    if (erateGraph < 0.0)
      fprintf(stderr, "Invalid overlap error threshold (-eg option); must be at least 0.0.\n");
    if (erateMax < 0.0)
      fprintf(stderr, "Invalid overlap error threshold (-eM option); must be at least 0.0.\n");

    if (prefix == NULL)
      fprintf(stderr, "No output prefix name (-o option) supplied.\n");

    if (gkpStorePath == NULL)
      fprintf(stderr, "No gatekeeper store (-G option) supplied.\n");

    if (ovlStoreUniqPath == NULL)
      fprintf(stderr, "No overlap store (-O option) supplied.\n");

    if ((ovlStoreUniqPath != NULL) && (ovlStoreUniqPath == ovlStoreReptPath))
      fprintf(stderr, "Too many overlap stores (-O option) supplied.\n");

    if (tigStorePath == NULL)
      fprintf(stderr, "No output tigStore (-T option) supplied.\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Graph  error threshold  = %.3f (%.3f%%)\n", erateGraph,  erateGraph  * 100);
  fprintf(stderr, "Max    error threshold  = %.3f (%.3f%%)\n", erateMax, erateMax * 100);
  fprintf(stderr, "\n");
  fprintf(stderr, "Minimum overlap length = %u bases\n", minOverlap);
  fprintf(stderr, "\n");

  if (numThreads > 0) {
    omp_set_num_threads(numThreads);
    fprintf(stderr, "number of threads     = %d (command line)\n", numThreads);
    fprintf(stderr, "\n");
  } else {
    fprintf(stderr, "number of threads     = %d (OpenMP default)\n", omp_get_max_threads());
    fprintf(stderr, "\n");
  }

  for (uint64 i=0, j=1; i<64; i++, j<<=1)
    if (logFileFlagSet(j))
      fprintf(stderr, "DEBUG                 = %s\n", logFileFlagNames[i]);

  gkStore          *gkpStore     = gkStore::gkStore_open(gkpStorePath);
  ovStore          *ovlStoreUniq = new ovStore(ovlStoreUniqPath, gkpStore);
  ovStore          *ovlStoreRept = ovlStoreReptPath ? new ovStore(ovlStoreReptPath, gkpStore) : NULL;

  UnitigVector      unitigs;

  setLogFile(prefix, NULL);

  FI = new FragmentInfo(gkpStore, prefix, minReadLen);

  // Initialize where we've been to nowhere
  Unitig::resetFragUnitigMap(FI->numFragments());

  OC = new OverlapCache(ovlStoreUniq, ovlStoreRept, prefix, MAX(erateMax, erateGraph), minOverlap, ovlCacheMemory, ovlCacheLimit, onlySave, doSave);
  OG = new BestOverlapGraph(erateGraph, deviationGraph, prefix);
  CG = new ChunkGraph(prefix);

  delete ovlStoreUniq;  ovlStoreUniq = NULL;
  delete ovlStoreRept;  ovlStoreRept = NULL;

  gkpStore->gkStore_close();
  gkpStore = NULL;

  //
  //  Build the initial unitig path from non-contained fragments.  The first pass is usually the
  //  only one needed, but occasionally (maybe) we miss fragments, so we make an explicit pass
  //  through all fragments and place whatever isn't already placed.
  //

  setLogFile(prefix, "buildUnitigs");
  writeLog("==> BUILDING UNITIGS from %d fragments.\n", FI->numFragments());

  for (uint32 fi=CG->nextFragByChunkLength(); fi>0; fi=CG->nextFragByChunkLength())
    populateUnitig(unitigs, fi);

  delete CG;
  CG = NULL;

  breakSingletonTigs(unitigs);

  reportOverlaps(unitigs, prefix, "buildUnitigs");
  reportUnitigs(unitigs, prefix, "buildUnitigs", genomeSize);

  //
  //  Place contained reads.
  //

#if 1
  setLogFile(prefix, "placeContains");

  unitigs.computeArrivalRate(prefix, "initial");
  unitigs.computeErrorProfiles(prefix, "initial");
  //unitigs.reportErrorProfiles(prefix, "initial");

  placeUnplacedUsingAllOverlaps(unitigs, prefix);

  reportOverlaps(unitigs, prefix, "placeContains");
  reportUnitigs(unitigs, prefix, "placeContains", genomeSize);
#endif

  //
  //  Merge tigs (and detect ciruclar ones too).  Contained reads need to be placed to 'clean up'
  //  the error rate.  Dovetail alone is too 'clean' for circular to be detected (in ecoli).
  //

#if 0
  setLogFile(prefix, "merge");

  computeErrorProfiles(unitigs, prefix, "merge");
  //reportErrorProfiles(unitigs, prefix, "merge");

  mergeUnitigs(unitigs, deviationGraph, false);

  reportOverlaps(unitigs, prefix, "merge");
  reportUnitigs(unitigs, prefix, "merge", genomeSize);
#endif


  //
  //  Pop bubbles
  //

#if 1
  setLogFile(prefix, "popBubbles");

  unitigs.computeErrorProfiles(prefix, "unplaced");
  //unitigs.reportErrorProfiles(prefix, "unplaced");

  popBubbles(unitigs,
             deviationBubble);

  //checkUnitigMembership(unitigs);
  reportOverlaps(unitigs, prefix, "popBubbles");
  reportUnitigs(unitigs, prefix, "popBubbles", genomeSize);
#endif

  //
  //  Detect and break repeats.  Annotate each read with overlaps to reads not overlapping in the tig,
  //  project these regions back to the tig, and break unless there is a read spanning the region.
  //

  setLogFile(prefix, "markRepeatReads");

  unitigs.computeErrorProfiles(prefix, "repeats");

  markRepeatReads(unitigs, deviationRepeat, confusedAbsolute, confusedPercent);

  //checkUnitigMembership(unitigs);
  reportOverlaps(unitigs, prefix, "markRepeatReads");
  reportUnitigs(unitigs, prefix, "markRepeatReads", genomeSize);

  //
  //  Try to reassemble just the split repeats.
  //

#if 0
  if (enableReconstructRepeats) {
    assert(enableShatterRepeats);
    setLogFile(prefix, "reconstructRepeats");

    reconstructRepeats(unitigs, erateGraph);

    //checkUnitigMembership(unitigs);
    reportOverlaps(unitigs, prefix, "reconstructRepeats");
    reportUnitigs(unitigs, prefix, "reconstructRepeats", genomeSize);
  }
#endif

  //
  //  Cleanup unitigs.  Break those that have gaps in them.  Place contains again.  For any read
  //  still unplaced, make it a singleton unitig.
  //

  setLogFile(prefix, "cleanup");

  splitDiscontinuousUnitigs(unitigs, minOverlap);

  breakSingletonTigs(unitigs);

  //unitigs.computeErrorProfiles(prefix, "final");
  //unitigs.reportErrorProfiles(prefix, "final");

  //placeUnplacedUsingAllOverlaps(unitigs, prefix);

  promoteToSingleton(unitigs);

  classifyUnitigsAsUnassembled(unitigs,
                               fewReadsNumber,
                               tooShortLength,
                               spanFraction,
                               lowcovFraction, lowcovDepth);

  //checkUnitigMembership(unitigs);
  reportOverlaps(unitigs, prefix, "final");
  reportUnitigs(unitigs, prefix, "final", genomeSize);

  //
  //  Generate outputs.  The graph MUST come after output, because it needs
  //  the tigStore tigID.
  //

  setLogFile(prefix, "output");

  setParentAndHang(unitigs);
  writeUnitigsToStore(unitigs, prefix, tigStorePath, fragment_count_target, true);

  setLogFile(prefix, "graph");

  writeUnusedEdges(unitigs, prefix);

  //
  //  Tear down bogart.
  //

  delete CG;
  delete OG;
  delete OC;
  delete FI;

  setLogFile(prefix, NULL);

  writeLog("Bye.\n");

  return(0);
}

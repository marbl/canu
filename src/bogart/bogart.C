
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

#include "system.H"

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_ChunkGraph.H"
#include "AS_BAT_AssemblyGraph.H"

#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_PopulateUnitig.H"
#include "AS_BAT_Instrumentation.H"
#include "AS_BAT_PlaceContains.H"

#include "AS_BAT_DetectSpurs.H"

#include "AS_BAT_MergeOrphans.H"
#include "AS_BAT_MarkRepeatReads.H"

#include "AS_BAT_SplitDiscontinuous.H"

#include "AS_BAT_FindCircular.H"

#include "AS_BAT_PromoteToSingleton.H"

#include "AS_BAT_SetParentAndHang.H"
#include "AS_BAT_Outputs.H"


ReadInfo         *RI  = 0L;
OverlapCache     *OC  = 0L;
BestOverlapGraph *OG  = 0L;
ChunkGraph       *CG  = 0L;



//  Check if we should terminate bogart early.  If so, close the log files,
//  stop threads, delete the global data, write a nice message and return
//  true.  Then, in main(), we can call return(0) to actually exit cleanly
//  (since exit(0) doesn't exit cleanly).
bool
terminateBogart(uint64 stop, char const *message) {

  if (stopFlagSet(stop) == false)
    return(false);

  setLogFile(nullptr, nullptr);   //  Close files.
  setNumThreads(1);               //  Hopefully kills off other threads.

  delete CG;
  delete OG;
  delete OC;
  delete RI;

  writeStatus("\n");
  writeStatus(message);

  return(true);
}



int
main (int argc, char * argv []) {
  char const  *seqStorePath            = NULL;
  char const  *ovlStorePath            = NULL;

  double       erateGraph               = 0.075;
  double       erateMax                 = 0.100;
  double       erateForced              = 1.0;     //  Used only if erateForced < 1.0

  bool         filterHighError          = true;
  bool         filterLopsided           = true;
  bool         filterSpur               = true;
  uint32       spurDepth                = 3;
  bool         filterDeadEnds           = true;

  covgapType   covGapType               = covgapUncovered;
  uint32       covGapOlap               = 500;     //  Require overlap of x bp when detecting coverage gaps.
  double       lopsidedDiff             = 25.0;    //  Call reads lopsided if diff between is more than x percent.
  double       minOlapPercent           =  0.0;
  double       minReadsBest             =  0.8;    //  80% of reads must have best edges.
  double       percentileError          =  0.9;

  uint64       genomeSize               = 0;

  uint32       fewReadsNumber           = 2;      //  Parameters for labeling of unassembled; also set in pipelines/canu/Defaults.pm
  uint32       tooShortLength           = 0;
  double       spanFraction             = 1.0;
  double       lowcovFraction           = 0.5;
  uint32       lowcovDepth              = 3;

  double       deviationGraph           = 6.0,    similarityGraph  = 0.00;
  double       deviationBubble          = 6.0,    similarityBubble = 0.01;
  double       deviationRepeat          = 3.0,    similarityRepeat = 0.01;

  uint32       confusedAbsolute         = 2500;
  double       confusedPercent          = 15.0;

  uint64       ovlCacheMemory           = UINT64_MAX;

  char const  *prefix                   = NULL;

  uint32       minReadLen               = 0;
  uint32       maxReadLen               = UINT32_MAX;

  uint32       minOverlapLen            = 500;
  uint32       minIntersectLen          = 500;
  uint32       maxPlacements            = 2;

  argc = AS_configure(argc, argv);

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      prefix = argv[++arg];


    } else if (strcmp(argv[arg], "-threads") == 0) {
      setNumThreads(argv[++arg]);

    } else if (strcmp(argv[arg], "-M") == 0) {
      ovlCacheMemory  = (uint64)(atof(argv[++arg]) * 1024 * 1024 * 1024);


    } else if (strcmp(argv[arg], "-gs") == 0) {
      genomeSize = strtoull(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-unassembled") == 0) {
      uint32  invalid = 0;

      if ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        fewReadsNumber  = atoi(argv[++arg]);
      else
        invalid++;

      if ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        tooShortLength  = atoi(argv[++arg]);
      else
        invalid++;

      if ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        spanFraction    = atof(argv[++arg]);
      else
        invalid++;

      if ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        lowcovFraction  = atof(argv[++arg]);
      else
        invalid++;

      if ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        lowcovDepth     = atoi(argv[++arg]);
      else
        invalid++;

      if (invalid) {
        char *s = new char [1024];
        snprintf(s, 1024, "Too few parameters to -unassembled option.\n");
        err.push_back(s);
      }

    } else if (strcmp(argv[arg], "-readlen") == 0) {
      decodeRange(argv[++arg], minReadLen, maxReadLen);

    } else if (strcmp(argv[arg], "-mr") == 0) {
      minReadLen = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-mo") == 0) {
      minOverlapLen = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-mi") == 0) {
      minIntersectLen = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-mp") == 0) {
      maxPlacements = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-eg") == 0) {
      erateGraph = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-eM") == 0) {
      erateMax = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-ef") == 0) {
      erateForced = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-ep") == 0) {
      percentileError = atof(argv[++arg]);

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

    } else if (strcmp(argv[arg], "-sg") == 0) {  //  Similarity threshold, graph, UNUSED
      similarityGraph = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-sb") == 0) {  //  Similarity threshold, bubble
      similarityBubble = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-sr") == 0) {  //  Similarity threshold, repeat, UNUSED
      similarityRepeat = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-sd") == 0) {  //  Depth to look for spurs
      spurDepth = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-nofilter") == 0) {
      ++arg;

      if ((arg >= argc) || (argv[arg][0] == '-')) {
        char *s = new char[1024];
        snprintf(s, 1024, "No filters specified for '-nofilter' option.\n");
        err.push_back(s);
      }

      else if (strcasecmp(argv[arg], "higherror") == 0)   { filterHighError = false; }
      else if (strcasecmp(argv[arg], "lopsided")  == 0)   { filterLopsided  = false; }
      else if (strcasecmp(argv[arg], "spur")      == 0)   { filterSpur      = false; }
      else if (strcasecmp(argv[arg], "deadends")  == 0)   { filterDeadEnds  = false; }
      else {
        char *s = new char[1024];
        snprintf(s, 1024, "Invalid filter '%s' specified for '-nofilter' option.\n", argv[arg]);
        err.push_back(s);
      }

    } else if (strcmp(argv[arg], "-minolappercent") == 0) {
      minOlapPercent = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-minreadsbest") == 0) {
      minReadsBest = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-covgaptype") == 0) {
      ++arg;
      if      (strcasecmp(argv[arg], "none")      == 0)   covGapType = covgapNone;
      else if (strcasecmp(argv[arg], "chimer")    == 0)   covGapType = covgapChimer;
      else if (strcasecmp(argv[arg], "uncovered") == 0)   covGapType = covgapUncovered;
      else if (strcasecmp(argv[arg], "deadend")   == 0)   covGapType = covgapDeadend;
      else {
        char *s = new char [1024];
        snprintf(s, 1024, "Unknown '-covgaptype' option '%s'.\n", argv[arg]);
        err.push_back(s);
      }

    } else if (strcmp(argv[arg], "-covgapolap") == 0) {
      covGapOlap = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-lopsided") == 0) {
      lopsidedDiff = strtodouble(argv[++arg]);

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
              (strcasecmp(logFileFlagNames[opt], "bestOverlaps") != 0) &&
              (strcasecmp(logFileFlagNames[opt], "overlapScoring") != 0) &&
              (strcasecmp(logFileFlagNames[opt], "errorProfiles") != 0) &&
              (strcasecmp(logFileFlagNames[opt], "optimizePositions") != 0) &&
              (strcasecmp(logFileFlagNames[opt], "chunkGraph") != 0) &&
              (strcasecmp(logFileFlagNames[opt], "setParentAndHang") != 0))
            logFileFlags |= flg;
        fnd = true;
      }
      if (fnd == false) {
        char *s = new char [1024];
        snprintf(s, 1024, "Unknown '-D' option '%s'.\n", argv[arg]);
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
        snprintf(s, 1024, "Unknown '-d' option '%s'.\n", argv[arg]);
        err.push_back(s);
      }

    } else if (strcmp(argv[arg], "-stop") == 0) {
      uint32  opt = 0;
      uint64  flg = 1;
      bool    fnd = false;
      for (arg++; stopFlagNames[opt]; flg <<= 1, opt++) {
        if (strcasecmp(stopFlagNames[opt], argv[arg]) == 0) {
          stopFlags |= flg;
          fnd = true;
        }
      }
      if (fnd == false) {
        char *s = new char [1024];
        snprintf(s, 1024, "Unknown '-stop' option '%s'.\n", argv[arg]);
        err.push_back(s);
      }


    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  if (genomeSize   == 0)       err.push_back("Genome size (-gs option) must be supplied\n");
  if (erateGraph    < 0.0)     err.push_back("Invalid overlap error threshold (-eg option); must be at least 0.0.\n");
  if (erateMax      < 0.0)     err.push_back("Invalid overlap error threshold (-eM option); must be at least 0.0.\n");
  if (prefix       == NULL)    err.push_back("No output prefix name (-o option) supplied.\n");
  if (seqStorePath == NULL)    err.push_back("No sequence store (-S option) supplied.\n");
  if (ovlStorePath == NULL)    err.push_back("No overlap store (-O option) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S seqPath -O ovlPath -T tigPath -o outPrefix ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Mandatory Parameters:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S seqPath     Mandatory path to an existing seqStore.\n");
    fprintf(stderr, "  -O ovlPath     Mandatory path to an existing ovlStore.\n");
    fprintf(stderr, "  -T tigPath     Mandatory path to an output tigStore (can exist or not).\n");
    fprintf(stderr, "  -o outPrefix   Mandatory prefix for the output files.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Process Options:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -threads T     Use at most T compute threads.\n");
    fprintf(stderr, "  -M gb          Use at most 'gb' gigabytes of memory.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Algorithm Options:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -gs            Genome size in bases.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -mr len        Force reads below 'len' bases to be singletons.\n");
    fprintf(stderr, "  -mo len        Ignore overlaps shorter than 'len' bases.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -mi len        Create unitigs from contig intersections of at least 'len' bases.\n");
    fprintf(stderr, "  -mp num        Create unitigs from contig intersections with at most 'num' placements.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nofilter [coverageGap | highError | lopsided | spur]\n");
    fprintf(stderr, "                 Disable filtering of:\n");
    fprintf(stderr, "                   coverageGap - reads that have a suspicious lack of overlaps in the middle\n");
    fprintf(stderr, "                   highError   - overlaps that have error rates well outside the observed\n");
    fprintf(stderr, "                   lopsided    - reads that have unusually asymmetric best overlaps\n");
    fprintf(stderr, "                   spur        - reads that have no overlaps on one end\n");
    fprintf(stderr, "                 Exactly one filter type must be supplied.  Use multiple -nofilter options\n");
    fprintf(stderr, "                 to disable multiple filters.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -eg F          Do not use overlaps more than F fraction error when when finding initial best edges.\n");
    fprintf(stderr, "  -eM F          Do not load overlaps more then F fraction error (useful only for -save).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -ef F          Force the final best edge threshold to be F fraction error.  By default,\n");
    fprintf(stderr, "                 bogart will analyze the overlaps loaded (-eM) that are also below the initial\n");
    fprintf(stderr, "                 edge threshold (-eg) and pick a final edge threshold; this option forces the \n");
    fprintf(stderr, "                 final edge threshold to exactly F.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -ep P          When deciding which overlaps to use, fall back to percentile P (0.0-1.0) if\n");
    fprintf(stderr, "                 the median error is 0.0, as commonly found in PacBio HiFi reads.  Default: 0.9\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -ca L          Split a contig if there is an alternate path from an overlap of at least L bases.\n");
    fprintf(stderr, "                 Default: 2100.\n");
    fprintf(stderr, "  -cp P          Split a contig if there is an alternate path from an overlap at most P percent\n");
    fprintf(stderr, "                 different from the length of the best overlap.  Default: 200.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -dg D          Use overlaps up to D standard deviations from the mean when building the best\n");
    fprintf(stderr, "                 overlap graph.  Default 6.0.\n");
    fprintf(stderr, "  -db D          Like -dg, but for merging bubbles into primary contigs.  Default 6.0.\n");
    fprintf(stderr, "  -dr D          Like -dg, but for breaking repeats.  Default 3.0.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Secret Algorithmic Options:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -covgaptype n      Set covgap pattern: none, chimer, uncovered, deadend.\n");
    fprintf(stderr, "  -covgapolap n      Require overlaps to overlap by at least n bases.\n");
    fprintf(stderr, "  -lopsided m n      Set how lopsided reads are detected and/or treated.\n");
    fprintf(stderr, "                       m = off        - don't detect at all (omit n parameter)\n");
    fprintf(stderr, "                       m = noseed n   - detect, n%% difference, allow edges to but don't seed overlaps with them\n");
    fprintf(stderr, "                       m = nobest n   - detect, n%% difference, exclude from bog graph completely\n");
    fprintf(stderr, "  -minolappercent f  Set a minimum overlap length, per overlap, as f*min(readAlen, readBlen)\n");
    fprintf(stderr, "  -minreadsbest f    Automatically relax overlap quality requirements if there are fewer\n");
    fprintf(stderr, "                     than fraction f (0.0-1.0) reads that have two best edges.  Default: 0.9\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Debugging and Logging\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -stop when     Stop after computing step 'when'..  Useful for debugging and testing.\n");
    fprintf(stderr, "                   'edges'      - after computing best edges; output asm.best.edges will exist.\n");
    fprintf(stderr, "                   'chunkgraph' - after building the ChunkGraph and logging.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -D <name>  enable logging/debugging for a specific component.\n");
    fprintf(stderr, "  -d <name>  disable logging/debugging for a specific component.\n");
    for (uint32 l=0; logFileFlagNames[l]; l++)
      fprintf(stderr, "               %s\n", logFileFlagNames[l]);
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "==> PARAMETERS.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Resources:\n");
  fprintf(stderr, "  Memory                " F_U64 " GB\n", ovlCacheMemory >> 30);
  fprintf(stderr, "  Compute Threads       %d\n", getNumThreads());
  fprintf(stderr, "\n");
  fprintf(stderr, "Lengths:\n");
  fprintf(stderr, "  Minimum read          %u bases\n",     minReadLen);
  fprintf(stderr, "  Maximum read          %u bases\n",     maxReadLen);
  fprintf(stderr, "  Minimum overlap       %u bases\n",     minOverlapLen);
  fprintf(stderr, "\n");
  fprintf(stderr, "Overlap Error Rates:\n");
  fprintf(stderr, "  Graph                 %.3f (%.3f%%)\n", erateGraph, erateGraph  * 100);
  fprintf(stderr, "  Max                   %.3f (%.3f%%)\n", erateMax,   erateMax    * 100);
  if (erateForced < 1.0)
  fprintf(stderr, "  Forced                %.3f (%.3f%%)\n", erateForced, erateForced  * 100);
  else
  fprintf(stderr, "  Forced                -.--- (-.---%%)   (not used)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Deviations:\n");
  fprintf(stderr, "  Graph                 %.3f\n", deviationGraph);
  fprintf(stderr, "  Bubble                %.3f\n", deviationBubble);
  fprintf(stderr, "  Repeat                %.3f\n", deviationRepeat);
  fprintf(stderr, "\n");
  fprintf(stderr, "Similarity Thresholds:\n");
  fprintf(stderr, "  Graph                 %.3f\n", similarityGraph);
  fprintf(stderr, "  Bubble                %.3f\n", similarityBubble);
  fprintf(stderr, "  Repeat                %.3f\n", similarityRepeat);
  fprintf(stderr, "\n");
  fprintf(stderr, "Edge Confusion:\n");
  fprintf(stderr, "  Absolute              %d\n",   confusedAbsolute);
  fprintf(stderr, "  Percent               %.4f\n", confusedPercent);
  fprintf(stderr, "\n");
  fprintf(stderr, "Unitig Construction:\n");
  fprintf(stderr, "  Minimum intersection  %u bases\n",     minIntersectLen);
  fprintf(stderr, "  Maxiumum placements   %u positions\n", maxPlacements);
  fprintf(stderr, "\n");
  fprintf(stderr, "Debugging Enabled:\n");

  if (logFileFlags == 0)
    fprintf(stderr, "  (none)\n");

  for (uint64 i=0, j=1; i<64; i++, j<<=1)
    if (logFileFlagSet(j))
      fprintf(stderr, "  %s\n", logFileFlagNames[i]);

  writeStatus("\n");
  writeStatus("==> LOADING AND FILTERING OVERLAPS.\n");
  writeStatus("\n");

  setLogFile(prefix, "filterOverlaps");

  RI = new ReadInfo(seqStorePath, prefix, minReadLen, maxReadLen);
  OC = new OverlapCache(ovlStorePath, prefix, std::max(erateMax, erateGraph), minOverlapLen, ovlCacheMemory, genomeSize);
  OG = new BestOverlapGraph(erateGraph,
                            std::max(erateMax, erateGraph),
                            erateForced,
                            percentileError,
                            deviationGraph,
                            minOlapPercent,
                            minReadsBest,
                            prefix,
                            covGapType, covGapOlap,
                            filterHighError,
                            filterLopsided, lopsidedDiff,
                            filterSpur, spurDepth);

  if (terminateBogart(STOP_BEST_EDGES, "Stopping after BestOverlapGraph() construction.\n"))
    return(0);

  CG = new ChunkGraph(prefix);

  if (terminateBogart(STOP_CHUNK_GRAPH, "Stopping after ChunkGraph() construction.\n"))
    return(0);

  //
  //  OG is used:
  //    in AssemblyGraph.C to decide if contained
  //    in DetectSpurs.C to find spurs
  //    in Instrumentation
  //    in MarkRepeatReads - contained, coverageGap, bestedges too
  //    in MergeOrphans to find contains
  //    in PlaceContains to find contains
  //    in PopulateUnitig
  //


  //
  //  Build the initial unitig path from non-contained reads.  The first pass is usually the
  //  only one needed, but occasionally (maybe) we miss reads, so we make an explicit pass
  //  through all reads and place whatever isn't already placed.
  //

  TigVector         contigs(RI->numReads());

  writeStatus("\n");
  writeStatus("==> BUILDING GREEDY TIGS.\n");
  writeStatus("\n");

  setLogFile(prefix, "buildGreedy");

  for (uint32 fi=CG->nextReadByChunkLength(); fi>0; fi=CG->nextReadByChunkLength())
    populateUnitig(contigs, fi);

  delete CG;
  CG = NULL;

  breakSingletonTigs(contigs);

  reportTigs(contigs, prefix, "buildGreedy", genomeSize);

  //  populateUnitig() uses only one hang from one overlap to compute the
  //  positions of reads.  Once all reads are (approximately) placed, compute
  //  positions using all overlaps.

  setLogFile(prefix, "buildGreedyOpt");
  contigs.optimizePositions(prefix, "buildGreedyOpt");
  reportTigs(contigs, prefix, "buildGreedyOpt", genomeSize);

  //  Break any tigs that aren't contiguous.

  setLogFile(prefix, "splitDiscontinuous");
  splitDiscontinuous(contigs, minOverlapLen);
  //reportOverlaps(contigs, prefix, "splitDiscontinuous");
  reportTigs(contigs, prefix, "splitDiscontinuous", genomeSize);

  //  Detect and fix spurs.

  setLogFile(prefix, "detectSpurs");
  detectSpurs(contigs);
  reportTigs(contigs, prefix, "detectSpurs", genomeSize);

  //
  //  For future use, remember the reads in contigs.
  //

  for (uint32 fid=1; fid<RI->numReads()+1; fid++)    //  This really should be incorporated
    if (contigs.inUnitig(fid) != 0)                  //  into populateUnitig()
      OG->setBackbone(fid);

  //
  //  Place contained reads.
  //

  writeStatus("\n");
  writeStatus("==> PLACE CONTAINED READS.\n");
  writeStatus("\n");

  setLogFile(prefix, "placeContains");

  //contigs.computeArrivalRate(prefix, "initial");
  contigs.computeErrorProfiles(prefix, "initial");
  contigs.reportErrorProfiles(prefix, "initial");

  std::set<uint32>   placedReads;

  placeUnplacedUsingAllOverlaps(contigs, deviationGraph, OG->reportErrorLimit(), prefix, placedReads);

  //  Compute positions again.  This fixes issues with contains-in-contains that
  //  tend to excessively shrink reads.  The one case debugged placed contains in
  //  a three read nanopore contig, where one of the contained reads shrank by 10%,
  //  which was enough to swap bgn/end coords when they were computed using hangs
  //  (that is, sum of the hangs was bigger than the placed read length).

  reportTigs(contigs, prefix, "placeContains", genomeSize);

  setLogFile(prefix, "placeContainsOpt");
  contigs.optimizePositions(prefix, "placeContainsOpt");
  reportTigs(contigs, prefix, "placeContainsOpt", genomeSize);

  setLogFile(prefix, "splitDiscontinuous");
  splitDiscontinuous(contigs, minOverlapLen);
  //reportOverlaps(contigs, prefix, "placeContains");
  reportTigs(contigs, prefix, "splitDiscontinuous", genomeSize);

  //
  //  Merge orphans.
  //

  writeStatus("\n");
  writeStatus("==> MERGE ORPHANS.\n");
  writeStatus("\n");

  setLogFile(prefix, "mergeOrphans");

  contigs.computeErrorProfiles(prefix, "unplaced");
  contigs.reportErrorProfiles(prefix, "unplaced");

  // we call this twice, once to merge in orphans, a second for bubbles
  mergeOrphans(contigs, deviationGraph, OG->reportErrorLimit(), false);

  writeStatus("\n");
  writeStatus("==> MARK SIMPLE BUBBLES.\n");
  writeStatus("    using %f user-specified threshold\n",  similarityBubble);
  writeStatus("\n");
  mergeOrphans(contigs, deviationBubble, similarityBubble, true);

  //checkUnitigMembership(contigs);
  //reportOverlaps(contigs, prefix, "mergeOrphans");
  reportTigs(contigs, prefix, "mergeOrphans", genomeSize);

#if 0
  {
    setLogFile(prefix, "reducedGraph");

    //  Build a new BestOverlapGraph, let it dump logs to 'reduced',
    //  then destroy the graph.
    //
    //  Note that this will fail an assert in BestOverlapGraph::removeLopsidedEdges(),
    //  and you'll need to disable it manually.

    fprintf(stderr, "\n");
    fprintf(stderr, "----------------------------------------\n");
    fprintf(stderr, "Building new graph after removing %u placed reads and %u bubble reads.\n",
            OG->numOrphan(),
            OG->numBubble());

    BestOverlapGraph *OGbf = new BestOverlapGraph(erateGraph,
                                                  deviationGraph, minOlapPercent,
                                                  "reduced",
                                                  filterCoverageGap, covGapOlap,
                                                  filterHighError,
                                                  filterLopsided, lopsidedDiff,
                                                  filterSpur, spurDepth,
                                                  OG);
    delete OGbf;

    //fprintf(stderr, "STOP after emitting OGbf.\n");
    //return(1);
    //exit(1);
  }
#endif

  //
  //  Initial construction done.  Classify what we have as assembled or unassembled.
  //

  classifyTigsAsUnassembled(contigs,
                            fewReadsNumber,
                            tooShortLength,
                            spanFraction,
                            lowcovFraction, lowcovDepth);

  //
  //  Generate a new graph using only edges that are compatible with existing tigs.
  //

  writeStatus("\n");
  writeStatus("==> GENERATING ASSEMBLY GRAPH.\n");
  writeStatus("\n");

  setLogFile(prefix, "assemblyGraph");

  contigs.computeErrorProfiles(prefix, "assemblyGraph");
  contigs.reportErrorProfiles(prefix, "assemblyGraph");

  AssemblyGraph *AG = new AssemblyGraph(prefix,
                                        deviationRepeat,
                                        OG->reportErrorLimit(),
                                        contigs);

  //AG->reportReadGraph(contigs, prefix, "initial");

  //
  //  Detect and break repeats.  Annotate each read with overlaps to reads not overlapping in the tig,
  //  project these regions back to the tig, and break unless there is a read spanning the region.
  //

  writeStatus("\n");
  writeStatus("==> BREAK REPEATS.\n");
  writeStatus("\n");

  setLogFile(prefix, "breakRepeats");

  contigs.computeErrorProfiles(prefix, "repeats");
  contigs.reportErrorProfiles(prefix, "repeats");

  markRepeatReads(AG, contigs, deviationRepeat, confusedAbsolute, confusedPercent);

  delete AG;
  AG = NULL;

  //checkUnitigMembership(contigs);
  //reportOverlaps(contigs, prefix, "markRepeatReads");
  reportTigs(contigs, prefix, "markRepeatReads", genomeSize);

  //
  //  Cleanup tigs.  Break those that have gaps in them.  Place contains again.  For any read
  //  still unplaced, make it a singleton unitig.
  //

  writeStatus("\n");
  writeStatus("==> CLEANUP MISTAKES.\n");
  writeStatus("\n");

  setLogFile(prefix, "cleanupMistakes");

  splitDiscontinuous(contigs, minOverlapLen);
  promoteToSingleton(contigs);

  if (filterDeadEnds) {
    splitDiscontinuous(contigs, minOverlapLen);
    promoteToSingleton(contigs);
  }

  writeStatus("\n");
  writeStatus("==> CLEANUP GRAPH.\n");
  writeStatus("\n");

  writeStatus("\n");
  writeStatus("==> GENERATE OUTPUTS.\n");
  writeStatus("\n");

  setLogFile(prefix, "generateOutputs");

  findCircularContigs(contigs, prefix);

  //checkUnitigMembership(contigs);
  reportOverlaps(contigs, prefix, "final");
  reportTigs(contigs, prefix, "final", genomeSize);

  setParentAndHang(contigs);
  writeTigsToStore(contigs, prefix, "ctg", true);

  //
  //  Tear down bogart.
  //

  terminateBogart(STOP_AT_END, "Bye.\n");
  return(0);
}


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

#include "AS_BAT_DropDeadEnds.H"

#include "AS_BAT_PromoteToSingleton.H"

#include "AS_BAT_CreateUnitigs.H"

#include "AS_BAT_SetParentAndHang.H"
#include "AS_BAT_Outputs.H"

#include "AS_BAT_TigGraph.H"


ReadInfo         *RI  = 0L;
OverlapCache     *OC  = 0L;
BestOverlapGraph *OG  = 0L;
ChunkGraph       *CG  = 0L;


int
main (int argc, char * argv []) {
  char const  *seqStorePath            = NULL;
  char const  *ovlStorePath            = NULL;

  double       erateGraph               = 0.075;
  double       erateMax                 = 0.100;

  bool         filterCoverageGap        = true;
  bool         filterHighError          = true;
  bool         filterLopsided           = true;
  bool         filterSpur               = true;
  uint32       spurDepth                = 3;
  bool         filterDeadEnds           = true;

  uint32       covGapOlap               = 500;     //  Require overlap of x bp when detecting coverage gaps.
  double       lopsidedDiff             = 25.0;    //  Call reads lopsided if diff between is more than x percent.
  double       minOlapPercent           =  0.0;

  uint64       genomeSize               = 0;

  uint32       fewReadsNumber           = 2;      //  Parameters for labeling of unassembled; also set in pipelines/canu/Defaults.pm
  uint32       tooShortLength           = 0;
  double       spanFraction             = 1.0;
  double       lowcovFraction           = 0.5;
  uint32       lowcovDepth              = 3;

  double       deviationGraph           = 6.0,    similarityGraph  = 0.0;
  double       deviationBubble          = 6.0,    similarityBubble = 0.1;
  double       deviationRepeat          = 3.0,    similarityRepeat = 0.1;

  uint32       confusedAbsolute         = 2100;
  double       confusedPercent          = 200.0;

  int32        numThreads               = 0;

  uint64       ovlCacheMemory           = UINT64_MAX;

  char const  *prefix                   = NULL;

  uint32       minReadLen               = 0;
  uint32       maxReadLen               = UINT32_MAX;

  uint32       minOverlapLen            = 500;
  uint32       minIntersectLen          = 500;
  uint32       maxPlacements            = 2;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      prefix = argv[++arg];


    } else if (strcmp(argv[arg], "-threads") == 0) {
      if ((numThreads = atoi(argv[++arg])) > 0)
        omp_set_num_threads(numThreads);

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
      filterCoverageGap = ((arg >= argc) || (strcasestr(argv[arg], "coverageGap") == NULL));
      filterCoverageGap = ((arg >= argc) || (strcasestr(argv[arg], "suspicious")  == NULL));   //  Deprecated!
      filterHighError   = ((arg >= argc) || (strcasestr(argv[arg], "higherror")   == NULL));
      filterLopsided    = ((arg >= argc) || (strcasestr(argv[arg], "lopsided")    == NULL));
      filterSpur        = ((arg >= argc) || (strcasestr(argv[arg], "spur")        == NULL));
      filterDeadEnds    = ((arg >= argc) || (strcasestr(argv[arg], "deadends")    == NULL));

    } else if (strcmp(argv[arg], "-minolappercent") == 0) {
      minOlapPercent = strtodouble(argv[++arg]);

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

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

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
    fprintf(stderr, "  -nofilter [coverageGap],[highError],[lopsided],[spur]\n");
    fprintf(stderr, "                 Disable filtering of:\n");
    fprintf(stderr, "                   coverageGap - reads that have a suspicious lack of overlaps in the middle\n");
    fprintf(stderr, "                   highError   - overlaps that have error rates well outside the observed\n");
    fprintf(stderr, "                   lopsided    - reads that have unusually asymmetric best overlaps\n");
    fprintf(stderr, "                   spur        - reads that have no overlaps on one end\n");
    fprintf(stderr, "                 The value supplied to -nofilter must be one word, case, order and punctuation\n");
    fprintf(stderr, "                 do not matter.  The following examples behave the same:\n");
    fprintf(stderr, "                    '-nofilter coverageGap,higherror'\n");
    fprintf(stderr, "                    '-nofilter coveragegap-and-HIGHERROR'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -eg F          Do not use overlaps more than F fraction error when when finding initial best edges.\n");
    fprintf(stderr, "  -eM F          Do not load overlaps more then F fraction error (useful only for -save).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -ca L          Split a contig if there is an alternate path from an overlap of at least L bases.\n");
    fprintf(stderr, "                 Default: 2100.\n");
    fprintf(stderr, "  -cp P          Split a contig if there is an alternate path from an overlap at most P percent\n");
    fprintf(stderr, "                 different from the length of the best overlap.  Default: 200.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -dg D          Use overlaps upto D standard deviations from the mean when building the best\n");
    fprintf(stderr, "                 overlap graph.  Default 6.0.\n");
    fprintf(stderr, "  -db D          Like -dg, but for merging bubbles into primary contigs.  Default 6.0.\n");
    fprintf(stderr, "  -dr D          Like -dg, but for breaking repeats.  Default 3.0.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Secret Algorithmic Options:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -covgapolap n      Require overlaps to overlap by at least n bases.\n");
    fprintf(stderr, "  -lopsided m n      Set how lopsided reads are detected and/or treated.\n");
    fprintf(stderr, "                       m = off        - don't detect at all (omit n parameter)\n");
    fprintf(stderr, "                       m = noseed n   - detect, n%% difference, allow edges to but don't seed overlaps with them\n");
    fprintf(stderr, "                       m = nobest n   - detect, n%% difference, exclude from bog graph completely\n");
    fprintf(stderr, "  -minolappercent f  Set a minimum overlap length, per overlap, as f*min(readAlen, readBlen)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Debugging and Logging\n");
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
  fprintf(stderr, "  Compute Threads       %d (%s)\n", omp_get_max_threads(), (numThreads > 0) ? "command line" : "OpenMP default");
  fprintf(stderr, "\n");
  fprintf(stderr, "Lengths:\n");
  fprintf(stderr, "  Minimum read          %u bases\n",     minReadLen);
  fprintf(stderr, "  Maximum read          %u bases\n",     maxReadLen);
  fprintf(stderr, "  Minimum overlap       %u bases\n",     minOverlapLen);
  fprintf(stderr, "\n");
  fprintf(stderr, "Overlap Error Rates:\n");
  fprintf(stderr, "  Graph                 %.3f (%.3f%%)\n", erateGraph, erateGraph  * 100);
  fprintf(stderr, "  Max                   %.3f (%.3f%%)\n", erateMax,   erateMax    * 100);
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
  OC = new OverlapCache(ovlStorePath, prefix, max(erateMax, erateGraph), minOverlapLen, ovlCacheMemory, genomeSize);
  OG = new BestOverlapGraph(erateGraph,
                            deviationGraph, minOlapPercent,
                            prefix,
                            filterCoverageGap, covGapOlap,
                            filterHighError,
                            filterLopsided, lopsidedDiff,
                            filterSpur, spurDepth);
  CG = new ChunkGraph(prefix);

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

  TigVector         contigs(RI->numReads());  //  Both initial greedy tigs and final contigs
  TigVector         unitigs(RI->numReads());  //  The 'final' contigs, split at every intersection in the graph

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
  //  For future use, remember the reads in contigs.  When we make unitigs, we'll
  //  require that every unitig end with one of these reads -- this will let
  //  us reconstruct contigs from the unitigs.
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

  set<uint32>   placedReads;

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

  vector<confusedEdge>  confusedEdges;

  markRepeatReads(AG, contigs, deviationRepeat, confusedAbsolute, confusedPercent, confusedEdges);

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

  //checkUnitigMembership(contigs);
  reportOverlaps(contigs, prefix, "final");
  reportTigs(contigs, prefix, "final", genomeSize);

  //
  //  unitigSource:
  //
  //  We want some way of tracking unitigs that came from the same contig.  Ideally,
  //  we'd be able to emit only the edges that would join unitigs into the original
  //  contig, but it's complicated by containments.  For example:
  //
  //    [----------------------------------]   CONTIG
  //    -------------                          UNITIG
  //              --------------------------   UNITIG
  //                         -------           UNITIG
  //
  //  So, instead, we just remember the set of unitigs that were created from each
  //  contig, and assume that any edge between those unitigs represents the contig.
  //  Which it totally doesn't -- any repeat in the contig collapses -- but is a
  //  good first attempt.
  //

  vector<tigLoc>  unitigSource;

  //  The graph must come first, to find circular contigs.

  reportTigGraph(contigs, unitigSource, prefix, "contigs");

  setParentAndHang(contigs);
  writeTigsToStore(contigs, prefix, "ctg", true);

  setLogFile(prefix, "tigGraph");

  writeStatus("\n");
  writeStatus("==> GENERATE UNITIGS.\n");
  writeStatus("\n");

  setLogFile(prefix, "generateUnitigs");

  contigs.computeErrorProfiles(prefix, "generateUnitigs");
  contigs.reportErrorProfiles(prefix, "generateUnitigs");

  createUnitigs(contigs, unitigs, minIntersectLen, maxPlacements, confusedEdges, unitigSource);

  splitDiscontinuous(unitigs, minOverlapLen, unitigSource);

  reportTigGraph(unitigs, unitigSource, prefix, "unitigs");

  setParentAndHang(unitigs);
  writeTigsToStore(unitigs, prefix, "utg", true);

  //
  //  Tear down bogart.
  //

  //  How bizarre.  Human regression of 2017-07-28-2128 deadlocked (apparently) when deleting OC.
  //  It had 31 threads in futex_wait, thread 1 was in delete of the second block of data.  CPU
  //  usage was 100% IIRC.  Reproducable, at least twice, possibly three times.  setLogFilePrefix
  //  was moved before the deletes in hope that it'll close down threads.  Certainly, it should
  //  close thread output files from createUnitigs.

  setLogFile(prefix, NULL);    //  Close files.
  omp_set_num_threads(1);      //  Hopefully kills off other threads.

  delete CG;
  delete OG;
  delete OC;
  delete RI;

  writeStatus("\n");
  writeStatus("Bye.\n");

  return(0);
}

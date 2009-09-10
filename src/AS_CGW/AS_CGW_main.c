
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

const char *mainid = "$Id: AS_CGW_main.c,v 1.76 2009-09-10 14:58:11 skoren Exp $";


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include <math.h>
#include <time.h>

#undef INSTRUMENT_CGW
#undef CHECK_CONTIG_ORDERS
#undef CHECK_CONTIG_ORDERS_INCREMENTAL

#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
#define POPULATE_COC_HASHTABLE 1
#else
#define POPULATE_COC_HASHTABLE 0
#endif

//  If -1, do not test or fix edges.  If 0, test but do not fix.  If
//  1, test and fix edges.
//
#define FIX_CONTIG_EDGES -1

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Input_CGW.h"
#include "Output_CGW.h"
#include "SplitChunks_CGW.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "Stats_CGW.h"
#include "AS_ALN_aligners.h"
#include "Instrument_CGW.h"
#include "AS_CGW_EdgeDiagnostics.h"
#include "fragmentPlacement.h"  //  for resolveSurrogates()


//  Defines the logical checkpoints

#define CHECKPOINT_AFTER_BUILDING_SCAFFOLDS         "ckp01-ABS"
#define CHECKPOINT_DURING_CLEANING_SCAFFOLDS        "ckp02-DCS"
#define CHECKPOINT_AFTER_CLEANING_SCAFFOLDS         "ckp03-ACD"
#define CHECKPOINT_DURING_1ST_SCAFF_MERGE           "ckp04-1SM-partial"
#define CHECKPOINT_AFTER_1ST_SCAFF_MERGE            "ckp05-1SM"
#define CHECKPOINT_AFTER_STONES                     "ckp06-AS"
#define CHECKPOINT_DURING_2ND_SCAFF_MERGE           "ckp07-2SM-partial"
#define CHECKPOINT_AFTER_2ND_SCAFF_MERGE            "ckp08-2SM"
#define CHECKPOINT_AFTER_FINAL_ROCKS                "ckp09-FR"
#define CHECKPOINT_AFTER_PARTIAL_STONES             "ckp10-PS"
#define CHECKPOINT_AFTER_FINAL_CONTAINED_STONES     "ckp11-FCS"
#define CHECKPOINT_AFTER_FINAL_CLEANUP              "ckp12-FC"
#define CHECKPOINT_AFTER_RESOLVE_SURROGATES         "ckp13-RS"



void AddUnitigOverlaps(GraphCGW_T *graph, char       *ovlFileName);
void DemoteUnitigsWithRBP(FILE *stream, GraphCGW_T *graph);
void CheckCITypes(ScaffoldGraphT *sgraph);
void RemoveSurrogateDuplicates(void);

int
main(int argc, char **argv) {

  //  Options controlling main

  int    generateOutput = 1;

  int    preMergeRezLevel = -1;
  int    repeatRezLevel   = 0;

  int    minSamplesForOverride = 100;

  int    outputOverlapOnlyContigEdges = 0;
  int    demoteSingletonScaffolds = TRUE;
  int    checkRepeatBranchPattern = FALSE;

  int    restartFromCheckpoint = -1;
  char  *restartFromLogical    = "ckp00";

  int    doResolveSurrogates               = 1;      //  resolveSurrogates
  int    placeAllFragsInSinglePlacedSurros = 0;      //  resolveSurrogates
  double cutoffToInferSingleCopyStatus     = 0.666;  //  resolveSurrogates

  int    firstFileArg = 0;

  char  filepath[2048];

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
  ContigOrientChecker * coc;
  coc = CreateContigOrientChecker();
  assert(coc != NULL);
#endif

  GlobalData = CreateGlobal_CGW();  

  GlobalData->maxSequencedbSize            = MAX_SEQUENCEDB_SIZE;
  GlobalData->repeatRezLevel               = repeatRezLevel;
  GlobalData->stoneLevel                   = 0;
  GlobalData->ignoreChaffUnitigs           = 0;
  GlobalData->performCleanupScaffolds      = 1;
  GlobalData->debugLevel                   = 0;
  GlobalData->cgbUniqueCutoff              = CGB_UNIQUE_CUTOFF;
  GlobalData->cgbDefinitelyUniqueCutoff    = CGB_UNIQUE_CUTOFF;
  GlobalData->cgbApplyMicrohetCutoff       = -1;         // This basically turns it off, unless enabled
  GlobalData->cgbMicrohetProb              = 1.e-05;     // scores less than this are considered repeats
  GlobalData->doInterleavedScaffoldMerging = 1;
  GlobalData->allowDemoteMarkedUnitigs     = TRUE;       // allow toggled unitigs to be demoted to be repeat if they were marked unique

  argc = AS_configure(argc, argv);

  int arg     = 1;
  int err     = 0;
  int unk[64] = {0};
  int unl     = 0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-C") == 0) {
      GlobalData->performCleanupScaffolds = 0;

    } else if (strcmp(argv[arg], "-p") == 0) {
      GlobalData->closurePlacement = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-D") == 0) {
      GlobalData->debugLevel = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      outputOverlapOnlyContigEdges = 1;

    } else if (strcmp(argv[arg], "-e") == 0) {
      GlobalData->cgbMicrohetProb = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-F") == 0) {
      GlobalData->allowDemoteMarkedUnitigs = FALSE;

    } else if (strcmp(argv[arg], "-G") == 0) {
      generateOutput = 0;

    } else if (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->Gatekeeper_Store_Name, argv[++arg]);

    } else if (strcmp(argv[arg], "-I") == 0) {
      GlobalData->ignoreChaffUnitigs = 1;

    } else if (strcmp(argv[arg], "-i") == 0) {
      GlobalData->cgbApplyMicrohetCutoff = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-j") == 0) {
      GlobalData->cgbUniqueCutoff = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-k") == 0) {
      GlobalData->cgbDefinitelyUniqueCutoff = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-M") == 0) {
      GlobalData->doInterleavedScaffoldMerging = 0;

    } else if (strcmp(argv[arg], "-m") == 0) {
      minSamplesForOverride = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-N") == 0) {
      restartFromLogical = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      strcpy(GlobalData->File_Name_Prefix, argv[++arg]);

    } else if (strcmp(argv[arg], "-p") == 0) {
      preMergeRezLevel = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-R") == 0) {
      restartFromCheckpoint = atoi(argv[++arg]);
    
    } else if (strcmp(argv[arg], "-r") == 0) {
      repeatRezLevel = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-S") == 0) {
      doResolveSurrogates               = 1;
      cutoffToInferSingleCopyStatus     = atof(argv[++arg]);
      placeAllFragsInSinglePlacedSurros = 0;

      if (cutoffToInferSingleCopyStatus == 0.0)
        doResolveSurrogates               = 0;

      if (cutoffToInferSingleCopyStatus < 0) {
        cutoffToInferSingleCopyStatus     = 0.0;
        placeAllFragsInSinglePlacedSurros = 1;
      }

    } else if (strcmp(argv[arg], "-s") == 0) {
      GlobalData->stoneLevel = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-u") == 0) {
      strcpy(GlobalData->unitigOverlaps, argv[++arg]);

    } else if (strcmp(argv[arg], "-Z") == 0) {
      demoteSingletonScaffolds = FALSE;

    } else if (strcmp(argv[arg], "-z") == 0) {
      checkRepeatBranchPattern = TRUE;

    } else if ((argv[arg][0] != '-') && (firstFileArg == 0)) {
      firstFileArg = arg;
      arg = argc;

    } else {
      unk[unl++] = arg;
      err++;
    }

    arg++;
  }

  if (GlobalData->Gatekeeper_Store_Name[0] == 0)
    err++;

  if (GlobalData->File_Name_Prefix[0] == 0)
    err++;

  if (cutoffToInferSingleCopyStatus > 1.0)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [options] -g <GatekeeperStoreName> -o <OutputPath> <unitigs*.cgb>\n", argv[0]);
    fprintf(stderr, "   -C           Don't cleanup scaffolds\n");
    fprintf(stderr, "   -p <int>     how to place closure reads. 0 - place at first location found, 1 - place at best gap, 2 - allow to be placed in multiple gaps\n");
    fprintf(stderr, "   -D <lvl>     Debug\n");
    fprintf(stderr, "   -E           output overlap only contig edges\n");
    fprintf(stderr, "   -e <thresh>  Microhet score probability cutoff\n");
    fprintf(stderr, "   -F           strongly enforce unique/repeat flag set in unitig, default if not set is to still\n");
    fprintf(stderr, "                  allow those marked unique to be demoted due to Repeat Branch Pattern or being\n");
    fprintf(stderr, "                  too small\n");
    fprintf(stderr, "   -g           gkp Store path (required)\n");
    fprintf(stderr, "   -G           Don't generate output (cgw or cam)\n");
    fprintf(stderr, "   -I           ignore chaff unitigs\n");
    fprintf(stderr, "   -i <thresh>  Set max coverage stat for microhet determination of non-uniqueness (default -1)\n");
    fprintf(stderr, "   -j <thresh>  Set min coverage stat for definite uniqueness\n");
    fprintf(stderr, "   -k <thresh>  Set max coverage stat for possible uniqueness\n");
    fprintf(stderr, "   -M           don't do interleaved scaffold merging\n");
    fprintf(stderr, "   -m <min>     Number of mate samples to recompute an insert size, default is 100\n");
    fprintf(stderr, "   -N <ckp>     restart from checkpoint location 'ckp' (see the timing file)\n");
    fprintf(stderr, "   -o           Output Name (required)\n");
    fprintf(stderr, "   -R <ckp>     restart from checkpoint file number 'ckp'\n");
    fprintf(stderr, "   -r <lvl>     repeat resolution level\n");
    fprintf(stderr, "   -S <t>       place all frags in singly-placed surrogates if at least fraction <x> can be placed\n");
    fprintf(stderr, "                  two special cases:\n");
    fprintf(stderr, "                    if <t> = -1, place all frags in singly-placed surrogates aggressively\n");
    fprintf(stderr, "                                 (which really mean t = 0.0, but triggers a better algorithm)\n");
    fprintf(stderr, "                    if <t> =  0, do not resolve surrogate fragments\n");
    fprintf(stderr, "   -s <lvl>     stone throwing level\n");
    fprintf(stderr, "   -u <file>    load these overlaps (from BOG) into the scaffold graph\n");
    fprintf(stderr, "   -v           verbose\n");
    fprintf(stderr, "   -Z           Don't demote singleton scaffolds\n");
    fprintf(stderr, "   -z           Turn on Check for Repeat Branch Pattern (demotes some unique unitigs to repeat)\n");

    fprintf(stderr, "\n");

    if (GlobalData->Gatekeeper_Store_Name[0] == 0)
      fprintf(stderr, "ERROR:  No gatekeeper (-g) supplied.\n");

    if (GlobalData->File_Name_Prefix[0] == 0)
      fprintf(stderr, "ERROR:  No output prefix (-o) supplied.\n");

    if (cutoffToInferSingleCopyStatus > 1.0)
      fprintf(stderr, "ERROR:  surrogate fraction cutoff (-S) must be between 0.0 and 1.0.\n");

    if (unl) {
      for (arg=0; arg<unl; arg++)
        fprintf(stderr, "ERROR:  Unknown option '%s'\n", argv[unk[arg]]);
    }

    exit(1);
  }


  if(GlobalData->cgbDefinitelyUniqueCutoff < GlobalData->cgbUniqueCutoff)
    GlobalData->cgbDefinitelyUniqueCutoff = GlobalData->cgbUniqueCutoff;


  if (preMergeRezLevel >= 0)
    GlobalData->repeatRezLevel = preMergeRezLevel;
  else
    GlobalData->repeatRezLevel = repeatRezLevel;


  if (generateOutput) {
    sprintf(filepath, "%s.cgw", GlobalData->File_Name_Prefix);
    GlobalData->cgwfp = File_Open(filepath, "w", TRUE);

    sprintf(filepath, "%s.cgw_contigs", GlobalData->File_Name_Prefix);
    GlobalData->ctgfp = File_Open(filepath, "w", TRUE);

    sprintf(filepath, "%s.cgw_scaffolds", GlobalData->File_Name_Prefix);
    GlobalData->scffp = File_Open(filepath, "w", TRUE);
  }

  if (strcasecmp(restartFromLogical, CHECKPOINT_AFTER_BUILDING_SCAFFOLDS) < 0) {
    fprintf(stderr, "Beginning CHECKPOINT_AFTER_BUILDING_SCAFFOLDS\n");

    //  Create the checkpoint from scratch
    //  TRUE -- doRezOnContigs
    ScaffoldGraph = CreateScaffoldGraph(TRUE, GlobalData->File_Name_Prefix);

    ProcessInput(GlobalData, firstFileArg, argc, argv);    // Handle the rest of the first file

    LoadDistData();

    fprintf(stderr,"* Splitting chimeric input unitigs\n");

    ComputeMatePairStatisticsRestricted(UNITIG_OPERATIONS, minSamplesForOverride, "unitig_preinitial");
    SplitInputUnitigs(ScaffoldGraph);
    ComputeMatePairStatisticsRestricted(UNITIG_OPERATIONS, minSamplesForOverride, "unitig_initial");

    BuildGraphEdgesDirectly(ScaffoldGraph->CIGraph);

    if (GlobalData->unitigOverlaps[0])
      AddUnitigOverlaps(ScaffoldGraph->CIGraph, GlobalData->unitigOverlaps);

    // Compute all overlaps implied by mate links between pairs of unique unitigs
    ComputeOverlaps(ScaffoldGraph->CIGraph, TRUE, TRUE);

    CheckEdgesAgainstOverlapper(ScaffoldGraph->CIGraph);
    CheckSurrogateUnitigs();

    MergeAllGraphEdges(ScaffoldGraph->CIGraph, FALSE, FALSE);

    CheckEdgesAgainstOverlapper(ScaffoldGraph->CIGraph);
    CheckSurrogateUnitigs();

    //  Mark some Unitigs/Chunks/CIs as repeats based on overlaps GRANGER 2/2/07
    //
    if (checkRepeatBranchPattern)
      DemoteUnitigsWithRBP(stderr, ScaffoldGraph->CIGraph);

    //  At this Point we've constructed the CIGraph

    {
      /* Looks like we never keep BOTH the CIGraph overlaps and the ContigGraph overlaps */
      ChunkOverlapperT *tmp = ScaffoldGraph->ContigGraph->overlapper;
      ScaffoldGraph->ContigGraph->overlapper = ScaffoldGraph->CIGraph->overlapper;
      ScaffoldGraph->CIGraph->overlapper = tmp;

      ScaffoldGraph->RezGraph = ScaffoldGraph->ContigGraph;
    }

    BuildInitialContigs(ScaffoldGraph);

    CheckCIScaffoldTs(ScaffoldGraph);

    if(GlobalData->debugLevel > 0){
      CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
      CheckSurrogateUnitigs();
    }

    CheckpointScaffoldGraph(CHECKPOINT_AFTER_BUILDING_SCAFFOLDS, "after building scaffolds");
  } else {
    LoadScaffoldGraphFromCheckpoint(GlobalData->File_Name_Prefix,restartFromCheckpoint, TRUE);

    //  Dump stats on the loaded checkpoint
    //GeneratePlacedContigGraphStats(tmpBuffer,0);
    //GenerateScaffoldGraphStats(tmpBuffer,0);
  }


  clearCacheSequenceDB(ScaffoldGraph->sequenceDB);



  if ((strcasecmp(restartFromLogical, CHECKPOINT_AFTER_CLEANING_SCAFFOLDS) < 0) &&
      (GlobalData->repeatRezLevel > 0)) {
    fprintf(stderr, "Beginning CHECKPOINT_AFTER_CLEANING_SCAFFOLDS\n");

    if(GlobalData->debugLevel > 0)
      DumpContigs(stderr,ScaffoldGraph, FALSE);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);

    // Transitive reduction of RezGraph followed by construction of SEdges
    RebuildScaffolds(ScaffoldGraph, TRUE);

    //  rocks is called inside of here
    //  checkpoints are written inside of here

    int iter     = 0;
    int iterMax  = 10;  //  MAX_OUTPUT_REZ_ITERATIONS
    int ctme     = time(0);
    int changed  = TRUE;

    fprintf(stderr,"** Running Level 1 Repeat Rez **\n");

    while ((changed) && (iter < iterMax)) {
      CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
      CheckCITypes(ScaffoldGraph);

      changed = RepeatRez(GlobalData->repeatRezLevel, GlobalData->File_Name_Prefix);

      if (changed){
        CheckCIScaffoldTs(ScaffoldGraph);

        // merge in stuff placed by rocks, assuming its position is correct!
        CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

        // Transitive reduction of RezGraph followed by construction of SEdges
        RebuildScaffolds(ScaffoldGraph, FALSE);

        //  If we've been running for 2 hours, AND we've not just
        //  completed the last iteration, checkpoint.
        //
        if ((time(0) - ctme > 120 * 60) && (iter+1 < iterMax)) {
          ctme = time(0);
          CheckpointScaffoldGraph(CHECKPOINT_DURING_CLEANING_SCAFFOLDS, "during scaffold cleaning");
        }

        iter++;
      }
    }

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    if(GlobalData->debugLevel > 0)
      DumpCIScaffolds(stderr,ScaffoldGraph, FALSE);

    CheckpointScaffoldGraph(CHECKPOINT_AFTER_CLEANING_SCAFFOLDS, "after scaffold cleaning");
  }
  //  else TidyUpScaffolds (ScaffoldGraph);


  if (strcasecmp(restartFromLogical, CHECKPOINT_AFTER_1ST_SCAFF_MERGE) < 0) {
    fprintf(stderr, "Beginning CHECKPOINT_AFTER_1ST_SCAFF_MERGE\n");

    CleanupScaffolds(ScaffoldGraph,FALSE, NULLINDEX, FALSE);

    CheckCIScaffoldTs(ScaffoldGraph);

    /* First we try to merge Scaffolds agressively */
    MergeScaffoldsAggressive(ScaffoldGraph, CHECKPOINT_DURING_1ST_SCAFF_MERGE, FALSE);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
    CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr, "---Checking contig orders after MergeScaffoldsAggressive (1)\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif

#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    CheckpointScaffoldGraph(CHECKPOINT_AFTER_1ST_SCAFF_MERGE, "after 1st scaffold merge");
  }


  clearCacheSequenceDB(ScaffoldGraph->sequenceDB);


  /*
    now that we are done with initial scaffold merge, we want to use the
    standard/default repeatRezLevel. Up to now, the value of preMergeRezLevel
    was in use if set on the command line
  */
  GlobalData->repeatRezLevel = repeatRezLevel;


  // Convert single-contig scaffolds that are marginally unique back
  // to unplaced contigs so they might be placed as stones
  //
  if (demoteSingletonScaffolds)
    DemoteSmallSingletonScaffolds();


  /* Now we throw stones */
  if ((strcasecmp(restartFromLogical, CHECKPOINT_AFTER_STONES) < 0) &&
      (GlobalData->stoneLevel > 0)) {
    fprintf(stderr, "Beginning CHECKPOINT_AFTER_STONES\n");

    CheckCIScaffoldTs(ScaffoldGraph);
    Throw_Stones(GlobalData->File_Name_Prefix, GlobalData->stoneLevel, FALSE);
    CheckCIScaffoldTs(ScaffoldGraph);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr, "---Checking contig orders after Throw_Stones\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif

#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    CheckpointScaffoldGraph(CHECKPOINT_AFTER_STONES, "after stone throwing");

    //GenerateLinkStats(ScaffoldGraph->CIGraph, "Stones", 0);
    //GeneratePlacedContigGraphStats("Stones", 0);
    //GenerateLinkStats(ScaffoldGraph->ContigGraph, "Stones", 0);
    //GenerateScaffoldGraphStats("Stones", 0);
  }



  if (strcasecmp(restartFromLogical, CHECKPOINT_AFTER_2ND_SCAFF_MERGE) < 0) {
    fprintf(stderr, "Beginning CHECKPOINT_AFTER_2ND_SCAFF_MERGE\n");

    CheckCIScaffoldTs(ScaffoldGraph);

    MergeScaffoldsAggressive(ScaffoldGraph, CHECKPOINT_DURING_2ND_SCAFF_MERGE, FALSE);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
    CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr, "---Checking contig orders after MergeScaffoldsAggressive (2)\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif

#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    CheckpointScaffoldGraph(CHECKPOINT_AFTER_2ND_SCAFF_MERGE, "after 2nd scaffold merge");
  }


  clearCacheSequenceDB(ScaffoldGraph->sequenceDB);

  if ((strcasecmp(restartFromLogical, CHECKPOINT_AFTER_FINAL_ROCKS) < 0) &&
      (GlobalData->repeatRezLevel > 0)) {
    const int  MAX_EXTRA_ROCKS_ITERS = 5;
    int  iter = 0, extra_rocks;

    fprintf(stderr, "Beginning CHECKPOINT_AFTER_FINAL_ROCKS\n");

    do {
      extra_rocks = Fill_Gaps(GlobalData, GlobalData->File_Name_Prefix, GlobalData->repeatRezLevel, iter);
      fprintf(stderr, "Threw additional %d rocks on iter %d\n", extra_rocks, iter);

      clearCacheSequenceDB(ScaffoldGraph->sequenceDB);
    } while (extra_rocks > 1 && iter < MAX_EXTRA_ROCKS_ITERS);

    CheckpointScaffoldGraph(CHECKPOINT_AFTER_FINAL_ROCKS, "after final rocks");
  }

  if ((strcasecmp(restartFromLogical, CHECKPOINT_AFTER_PARTIAL_STONES) < 0) &&
      (GlobalData->stoneLevel > 0)) {
    fprintf(stderr, "Beginning CHECKPOINT_AFTER_PARTIAL_STONES\n");

    CheckCIScaffoldTs (ScaffoldGraph);

    int partial_stones = Throw_Stones(GlobalData->File_Name_Prefix, GlobalData->stoneLevel, TRUE);

    CheckCIScaffoldTs (ScaffoldGraph);
    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);

    clearCacheSequenceDB(ScaffoldGraph->sequenceDB);

    fprintf (stderr, "Threw %d partial stones\n", partial_stones);
#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr,
            "---Checking contig orders after partial_stones\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    CheckpointScaffoldGraph(CHECKPOINT_AFTER_PARTIAL_STONES, "after partial stones");

    //GenerateLinkStats (ScaffoldGraph->CIGraph, "PStones", 0);
    //GeneratePlacedContigGraphStats ("PStones", 0);
    //GenerateLinkStats(ScaffoldGraph->ContigGraph, "PStones", 0);
    //GenerateScaffoldGraphStats ("PStones", 0);
  }

  if ((strcasecmp(restartFromLogical, CHECKPOINT_AFTER_FINAL_CONTAINED_STONES) < 0) &&
      (GlobalData->stoneLevel > 0)) {
    fprintf(stderr, "Beginning CHECKPOINT_AFTER_FINAL_CONTAINED_STONES\n");

    CheckCIScaffoldTs (ScaffoldGraph);
    int contained_stones = Toss_Contained_Stones (GlobalData->File_Name_Prefix, GlobalData->stoneLevel, 0);
    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
    CheckCIScaffoldTs (ScaffoldGraph);
    fprintf (stderr, "**** Finished Final Contained Stones level %d ****\n", GlobalData->stoneLevel);

    CleanupScaffolds (ScaffoldGraph, FALSE, NULLINDEX, FALSE);
    clearCacheSequenceDB(ScaffoldGraph->sequenceDB);

    fprintf(stderr, "Threw %d contained stones\n", contained_stones);

    // Remove copies of surrogates which are placed multiple times in the same place in a contig

    RemoveSurrogateDuplicates();

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr, "---Checking contig orders after contained_stones\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    CheckpointScaffoldGraph(CHECKPOINT_AFTER_FINAL_CONTAINED_STONES, "after final contained stones");

    //GenerateLinkStats (ScaffoldGraph->CIGraph, "CStones", 0);
    //GeneratePlacedContigGraphStats ("CStones", 0);
    //GenerateLinkStats(ScaffoldGraph->ContigGraph, "CStones", 0);
    //GenerateScaffoldGraphStats ("CStones", 0);
  }


  if (strcasecmp(restartFromLogical, CHECKPOINT_AFTER_FINAL_CLEANUP) < 0) {
    fprintf(stderr, "Beginning CHECKPOINT_AFTER_FINAL_CLEANUP\n");

    // Try to cleanup failed merges, and if we do, generate a checkpoint
    if(CleanupFailedMergesInScaffolds(ScaffoldGraph)){
      ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);

      // This call deletes surrogate-only contigs that failed to merge
      if(CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, TRUE)){

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
        fprintf(stderr, "---Checking contig orders after final cleanup\n\n");
        CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
        ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
      }
      CheckpointScaffoldGraph(CHECKPOINT_AFTER_FINAL_CLEANUP, "after final cleanup");
    }
  }


  if ((strcasecmp(restartFromLogical, CHECKPOINT_AFTER_RESOLVE_SURROGATES) < 0) &&
      (doResolveSurrogates > 0)) {
    fprintf(stderr, "Beginning CHECKPOINT_AFTER_RESOLVE_SURROGATES\n");

    resolveSurrogates(placeAllFragsInSinglePlacedSurros, cutoffToInferSingleCopyStatus);
    // Call resolve surrogate twice, this is necessary for finishing (closure) reads.
    // Consider a closure read and its two bounding reads, named left and right:
    //    If one (right) is placed in a unique region while the other (left) is in a surrogate itself, the closure read cannot be placed
    //    However, once the surrogate bounding read is placed (and fully incorporated which happens at the very end of resolveSurrogates)
    //    the closure read can be placed. 
    //    Therefore, we run resolve surrogates twice. 
    // Note that is closure reads are themselves mated, it may be necessary to do a third round of placement.  
    resolveSurrogates(placeAllFragsInSinglePlacedSurros, cutoffToInferSingleCopyStatus);
    
    CheckpointScaffoldGraph(CHECKPOINT_AFTER_RESOLVE_SURROGATES, "after resolve surrogates");
  }

  //  This generates the 'rezlog/gapreads' file.  It's hugely
  //  expensive, usually dies on a negative variance assert, and as
  //  far as BPW knows, unused.
  //
  //Show_Reads_In_Gaps (GlobalData->File_Name_Prefix);

  ComputeMatePairStatisticsRestricted(SCAFFOLD_OPERATIONS, minSamplesForOverride, "scaffold_final");
  ComputeMatePairStatisticsRestricted(CONTIG_OPERATIONS, minSamplesForOverride, "contig_final");

  GenerateCIGraph_U_Stats();
  GenerateLinkStats(ScaffoldGraph->CIGraph,"final",0);
  GeneratePlacedContigGraphStats("final",0);
  GenerateLinkStats(ScaffoldGraph->ContigGraph,"final",0);
  GenerateScaffoldGraphStats("final",0);
  GenerateSurrogateStats("final");
  
#ifdef DEBUG
  int j = 0;
  for (j = 0; j < GetNumVA_CIFragT(ScaffoldGraph->CIFrags); j++) {
    CIFragT * frag = GetCIFragT(ScaffoldGraph->CIFrags, j);
         
    if (ScaffoldGraph->gkpStore->gkStore_getFRGtoPLC(frag->iid) != 0) {
      AS_UID uid = getGatekeeperIIDtoUID(ScaffoldGraph->gkpStore, frag->iid, AS_IID_FRG);
      if (frag->contigID != -1) {
        ChunkInstanceT * ctg = GetGraphNode(ScaffoldGraph->RezGraph, frag->contigID);            
        fprintf(stderr, "CLOSURE_READS: CLOSURE READ %s PLACED=%d CHAFF=%d SINGLETON=%d IN ASM type %c in SCF %d\n", AS_UID_toString(uid), frag->flags.bits.isPlaced, frag->flags.bits.isChaff, frag->flags.bits.isSingleton, frag->type, ctg->scaffoldID);
      }
    }
  }
#endif

  clearCacheSequenceDB(ScaffoldGraph->sequenceDB);

  FixupLengthsScaffoldTs(ScaffoldGraph);

  if(generateOutput){
    CelamyAssembly(GlobalData->File_Name_Prefix);

    MarkContigEdges();
    OutputMateDists(ScaffoldGraph);

    ComputeMatePairDetailedStatus();
    OutputFrags(ScaffoldGraph);

    // We always have multiAlignments for Unitigs
    OutputUnitigsFromMultiAligns();
    OutputUnitigLinksFromMultiAligns();

    if(GlobalData->debugLevel > 0){
      DumpContigs(stderr,ScaffoldGraph, FALSE);
      DumpCIScaffolds(stderr,ScaffoldGraph, FALSE);
    }

    OutputContigsFromMultiAligns();
    OutputContigLinks(ScaffoldGraph, outputOverlapOnlyContigEdges);

    OutputScaffolds(ScaffoldGraph);
    OutputScaffoldLinks(ScaffoldGraph);
  }

  DestroyScaffoldGraph(ScaffoldGraph);

  DeleteGlobal_CGW(GlobalData);

  fprintf(stderr,"* Bye *\n");

  exit(0);
}

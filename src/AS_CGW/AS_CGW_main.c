
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

const char *mainid = "$Id: AS_CGW_main.c,v 1.111 2012-09-10 10:55:44 brianwalenz Exp $";

#undef CHECK_CONTIG_ORDERS
#undef CHECK_CONTIG_ORDERS_INCREMENTAL

#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
#define POPULATE_COC_HASHTABLE 1
#else
#define POPULATE_COC_HASHTABLE 0
#endif

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
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "Stats_CGW.h"
#include "ChunkOverlap_CGW.h"

#include "Instrument_CGW.h"
#include "fragmentPlacement.h"  //  for resolveSurrogates()

#include "CIScaffoldT_Analysis.H"  //  For checking mates on load

#include <omp.h>

//  Defines the logical checkpoints

char *ckpNames[16] = { "ckp00-NUL",
                       "ckp01-INI",
                       "ckp02-EDG",
                       "ckp03-SCF-partial",
                       "ckp04-SCF",
                       "ckp05-1SM-partial",
                       "ckp06-1SM",
                       "ckp07-AS",
                       "ckp08-2SM-partial",
                       "ckp09-2SM",
                       "ckp10-FR",
                       "ckp11-PS",
                       "ckp12-FCS",
                       "ckp13-FC",
                       "ckp14-RS",
                       "ckp15-FIN" };

#define CHECKPOINT_AFTER_LOADING                    1
#define CHECKPOINT_AFTER_EDGE_BUILDING              2
#define CHECKPOINT_DURING_INITIAL_SCAFFOLDING       3
#define CHECKPOINT_AFTER_INITIAL_SCAFFOLDING        4
#define CHECKPOINT_DURING_1ST_SCAFF_MERGE           5
#define CHECKPOINT_AFTER_1ST_SCAFF_MERGE            6
#define CHECKPOINT_AFTER_STONES                     7
#define CHECKPOINT_DURING_2ND_SCAFF_MERGE           8
#define CHECKPOINT_AFTER_2ND_SCAFF_MERGE            9
#define CHECKPOINT_AFTER_FINAL_ROCKS                10
#define CHECKPOINT_AFTER_PARTIAL_STONES             11
#define CHECKPOINT_AFTER_FINAL_CONTAINED_STONES     12
#define CHECKPOINT_AFTER_FINAL_CLEANUP              13
#define CHECKPOINT_AFTER_RESOLVE_SURROGATES         14
#define CHECKPOINT_AFTER_OUTPUT                     15


void
isValidCheckpointName(char *ckpName) {
  uint32 i=0;

  for (; i<16; i++)
    if (strcasecmp(ckpName, ckpNames[i]) == 0)
      break;

  if (i < 16)
    return;

  fprintf(stderr, "Invalid checkpoint name (-N) '%s'\n", ckpName);
  fprintf(stderr, "Valid names are: \n");

  for (i=1; i<16; i++)
    fprintf(stderr, "  %s\n", ckpNames[i]);

  exit(1);
}


bool
isThisCheckpoint(char *ckpName, uint32 level) {
  return(strcasecmp(ckpName, ckpNames[level]) == 0);
}

//  True if the checkpoint we loaded (ckpName) is earlier than the checkpoint we are at (level).
bool
runThisCheckpoint(char *ckpName, uint32 level) {

  if (strcasecmp(ckpName, ckpNames[level]) < 0) {
    fprintf(stderr, "Beginning %s (loaded %s).\n", ckpNames[level], ckpName);
    return(true);
  }

  fprintf(stderr, "Skipping %s (loaded %s).\n", ckpNames[level], ckpName);

  return(false);
}





void DemoteUnitigsWithRBP(FILE *stream, GraphCGW_T *graph);
void CheckCITypes(ScaffoldGraphT *sgraph);
void RemoveSurrogateDuplicates(void);

int
main(int argc, char **argv) {

  //  Options controlling main

  int    generateOutput = 1;

  int    preMergeRezLevel = -1;
  int    repeatRezLevel   = 0;

  int    restartFromCheckpoint = -1;
  char  *restartFromLogical    = "ckp00-NUL";

  bool   recomputeLeastSquaresOnLoad = false;

  int    doResolveSurrogates               = 1;      //  resolveSurrogates
  int    placeAllFragsInSinglePlacedSurros = 0;      //  resolveSurrogates
  double cutoffToInferSingleCopyStatus     = 0.666;  //  resolveSurrogates

  int    firstFileArg = 0;

  int32  outputFragsPerPartition = 0;

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
  ContigOrientChecker * coc;
  coc = CreateContigOrientChecker();
  assert(coc != NULL);
#endif

  //  temporary!
  fprintf(stderr, "Using up to %d OpenMP threads.\n", omp_get_max_threads());

  GlobalData = new Globals_CGW();

  argc = AS_configure(argc, argv);

  int arg     = 1;
  int err     = 0;
  int unk[64] = {0};
  int unl     = 0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-C") == 0) {
      GlobalData->performCleanupScaffolds = 0;

    } else if (strcmp(argv[arg], "-D") == 0) {
      GlobalData->debugLevel = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      GlobalData->outputOverlapOnlyContigEdges = 1;

    } else if (strcmp(argv[arg], "-e") == 0) {
      GlobalData->cgbMicrohetProb = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-F") == 0) {
      GlobalData->allowDemoteMarkedUnitigs = FALSE;

    } else if (strcmp(argv[arg], "-G") == 0) {
      generateOutput = 0;

    } else if (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->gkpStoreName, argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {
      strcpy(GlobalData->tigStoreName, argv[++arg]);

    } else if (strcmp(argv[arg], "-I") == 0) {
      GlobalData->ignoreChaffUnitigs = 1;

    } else if (strcmp(argv[arg], "-i") == 0) {
      GlobalData->cgbApplyMicrohetCutoff = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-j") == 0) {
      GlobalData->cgbUniqueCutoff = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-K") == 0) {
      GlobalData->removeNonOverlapingContigsFromScaffold = 1;

    } else if (strcmp(argv[arg], "-k") == 0) {
      GlobalData->cgbDefinitelyUniqueCutoff = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-m") == 0) {
      GlobalData->minSamplesForOverride = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-N") == 0) {
      restartFromLogical = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      strcpy(GlobalData->outputPrefix, argv[++arg]);

    } else if (strcmp(argv[arg], "-B") == 0) {
      outputFragsPerPartition = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-P") == 0) {
      GlobalData->closurePlacement = atoi(argv[++arg]);

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

    } else if (strcmp(argv[arg], "-filter") == 0) {
      GlobalData->mergeFilterLevel = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-shatter") == 0) {
      GlobalData->shatterLevel = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-missingMate") == 0) {
      GlobalData->mergeScaffoldMissingMates = atof(argv[++arg]);

      // the value is a percentage between 0 and 1 so make sure it never goes out of those bounds
      if (GlobalData->mergeScaffoldMissingMates < 0) {
    	  GlobalData->mergeScaffoldMissingMates = -1;
      } else if (GlobalData->mergeScaffoldMissingMates > 1) {
    	  GlobalData->mergeScaffoldMissingMates = 1;
      }

    } else if (strcmp(argv[arg], "-U") == 0) {
      GlobalData->doUnjiggleWhenMerging = 1;

    } else if (strcmp(argv[arg], "-u") == 0) {
      fprintf(stderr, "Option -u is broken.\n");
      exit(1);
      strcpy(GlobalData->unitigOverlaps, argv[++arg]);

    } else if (strcmp(argv[arg], "-Z") == 0) {
      GlobalData->demoteSingletonScaffolds = FALSE;

    } else if (strcmp(argv[arg], "-z") == 0) {
      GlobalData->checkRepeatBranchPattern = TRUE;

    } else if (strcmp(argv[arg], "-minmergeweight") == 0) {
      GlobalData->minWeightToMerge = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-recomputegaps") == 0) {
      recomputeLeastSquaresOnLoad = true;

    } else if ((argv[arg][0] != '-') && (firstFileArg == 0)) {
      firstFileArg = arg;
      arg = argc;

    } else {
      unk[unl++] = arg;
      err++;
    }

    arg++;
  }

  if (GlobalData->gkpStoreName[0] == 0)
    err++;

  if (GlobalData->outputPrefix[0] == 0)
    err++;

  if (cutoffToInferSingleCopyStatus > 1.0)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [options] -g <GatekeeperStoreName> -o <OutputPath> <unitigs*.cgb>\n", argv[0]);
    fprintf(stderr, "   -C                     Don't cleanup scaffolds\n");    
    fprintf(stderr, "   -D <lvl>               Debug\n");
    fprintf(stderr, "   -E                     output overlap only contig edges\n");
    fprintf(stderr, "   -e <thresh>            Microhet score probability cutoff\n");
    fprintf(stderr, "   -F                     strongly enforce unique/repeat flag set in unitig, default if not set is to still\n");
    fprintf(stderr, "                              allow those marked unique to be demoted due to Repeat Branch Pattern or being\n");
    fprintf(stderr, "                              too small\n");
    fprintf(stderr, "   -g                     gkp Store path (required)\n");
    fprintf(stderr, "   -G                     Don't generate output (cgw or cam)\n");
    fprintf(stderr, "   -I                     ignore chaff unitigs\n");
    fprintf(stderr, "   -i <thresh>            Set max coverage stat for microhet determination of non-uniqueness (default -1)\n");
    fprintf(stderr, "   -j <thresh>            Set min coverage stat for definite uniqueness\n");
    fprintf(stderr, "   -K                     Allow kicking out a contig placed in a scaffold by mate pairs that has no overlaps\n");
    fprintf(stderr, "                            to both its left and right neighbor contigs.\n");
    fprintf(stderr, "   -k <thresh>            Set max coverage stat for possible uniqueness\n");
    fprintf(stderr, "   -M                     don't do interleaved scaffold merging\n");
    fprintf(stderr, "   -m <min>               Number of mate samples to recompute an insert size, default is 100\n");
    fprintf(stderr, "   -N <ckp>               restart from checkpoint location 'ckp' (see the timing file)\n");
    fprintf(stderr, "   -o                     Output Name (required)\n");
    fprintf(stderr, "   -P <int>               how to place closure reads.\n");
    fprintf(stderr, "                              0 - place at first location found\n");
    fprintf(stderr, "                              1 - place at best gap\n");
    fprintf(stderr, "                              2 - allow to be placed in multiple gaps\n");
    fprintf(stderr, "   -R <ckp>               restart from checkpoint file number 'ckp'\n");
    fprintf(stderr, "   -r <lvl>               repeat resolution level\n");
    fprintf(stderr, "   -S <t>                 place all frags in singly-placed surrogates if at least fraction <x> can be placed\n");
    fprintf(stderr, "                          two special cases:\n");
    fprintf(stderr, "                              if <t> = -1, place all frags in singly-placed surrogates aggressively\n");
    fprintf(stderr, "                                           (which really mean t = 0.0, but triggers a better algorithm)\n");
    fprintf(stderr, "                              if <t> =  0, do not resolve surrogate fragments\n");
    fprintf(stderr, "   -s <lvl>               stone throwing level\n");
    fprintf(stderr, "   -shatter <thresh>      Set threshold for shattering scaffolds when loading from checkpoint. Any contigs\n");
    fprintf(stderr, "                            connected to a scaffold only by edges with less weight than the threshold will be\n");
    fprintf(stderr, "                            split into a new scaffold (default OFF)\n");
    fprintf(stderr, "   -missingMate <thresh>  Set threshold (0-1) for the percentage of mates (out of total) that are allowed to be\n");
    fprintf(stderr, "                            missing when attempting a scaffold merge (default 0). A value of -1 will ignore all\n");
    fprintf(stderr, "                            missing mates\n");
    fprintf(stderr, "   -minmergeweight <w>    Only use weight w or better edges for merging scaffolds.\n");
    fprintf(stderr, "   -recomputegaps         if loading a checkpoint, recompute gaps, merging contigs and splitting low weight scaffolds.\n");
    fprintf(stderr, "   -U                     after inserting rocks/stones try shifting contig positions back to their original location\n");
    fprintf(stderr, "                            when computing overlaps to see if they overlap with the rock/stone and allow them to merge\n");
    fprintf(stderr, "                            if they do\n");
    fprintf(stderr, "   -u <file>              load these overlaps (from BOG) into the scaffold graph\n");
    fprintf(stderr, "   -v                     verbose\n");
    fprintf(stderr, "   -Z                     Don't demote singleton scaffolds\n");
    fprintf(stderr, "   -z                     Turn on Check for Repeat Branch Pattern (demotes some unique unitigs to repeat)\n");

    fprintf(stderr, "\n");

    if (GlobalData->gkpStoreName[0] == 0)
      fprintf(stderr, "ERROR:  No gatekeeper (-g) supplied.\n");

    if (GlobalData->outputPrefix[0] == 0)
      fprintf(stderr, "ERROR:  No output prefix (-o) supplied.\n");

    if (cutoffToInferSingleCopyStatus > 1.0)
      fprintf(stderr, "ERROR:  surrogate fraction cutoff (-S) must be between 0.0 and 1.0.\n");

    if (unl) {
      for (arg=0; arg<unl; arg++)
        fprintf(stderr, "ERROR:  Unknown option '%s'\n", argv[unk[arg]]);
    }

    exit(1);
  }

  isValidCheckpointName(restartFromLogical);

  if(GlobalData->cgbDefinitelyUniqueCutoff < GlobalData->cgbUniqueCutoff)
    GlobalData->cgbDefinitelyUniqueCutoff = GlobalData->cgbUniqueCutoff;


  if (preMergeRezLevel >= 0)
    GlobalData->repeatRezLevel = preMergeRezLevel;
  else
    GlobalData->repeatRezLevel = repeatRezLevel;


  if (runThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_LOADING) == true) {
    int ctme     = time(0);

    //  Create the checkpoint from scratch
    ScaffoldGraph = CreateScaffoldGraph(GlobalData->outputPrefix);

    ProcessInput(firstFileArg, argc, argv);

    LoadDistData();

    //  This also labels unitigs as potential rocks / stones.
    //  That is the only real reason to call it here.  Insert sizes are set already.
    ComputeMatePairStatisticsRestricted(UNITIG_OPERATIONS, GlobalData->minSamplesForOverride, "unitig_initial");

    if (time(0) - ctme > 60 * 60)
      CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_LOADING], "after loading");

  } else if (isThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_LOADING) == true) {
    //  Load the checkpoint if we are exactly after loading, otherwise, fall through to the
    //  real load.
    LoadScaffoldGraphFromCheckpoint(GlobalData->outputPrefix,restartFromCheckpoint, TRUE);
  }


  if (runThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_EDGE_BUILDING) == true) {
    vector<CDS_CID_t>  rawEdges;

    BuildGraphEdgesDirectly(ScaffoldGraph->CIGraph, rawEdges);

    //  Broken, see comments in ChunkOverlap_CGW.c
    //
    //if (GlobalData->unitigOverlaps[0])
    //  AddUnitigOverlaps(ScaffoldGraph->CIGraph, GlobalData->unitigOverlaps, rawEdges);

    // Compute all overlaps implied by mate links between pairs of unique unitigs
    ComputeOverlaps(ScaffoldGraph->CIGraph, rawEdges);

    MergeAllGraphEdges(ScaffoldGraph->CIGraph, rawEdges, FALSE, FALSE);

    CheckEdgesAgainstOverlapper(ScaffoldGraph->CIGraph);
    CheckSurrogateUnitigs();

    //  Mark some Unitigs/Chunks/CIs as repeats based on overlaps GRANGER 2/2/07
    //
    if (GlobalData->checkRepeatBranchPattern)
      DemoteUnitigsWithRBP(stderr, ScaffoldGraph->CIGraph);

    //  At this Point we've constructed the CIGraph

    BuildInitialContigs(ScaffoldGraph);

    if(GlobalData->debugLevel > 0){
      CheckEdgesAgainstOverlapper(ScaffoldGraph->ContigGraph);
      CheckSurrogateUnitigs();
    }

    CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_EDGE_BUILDING], "after building edges");
  } else {
    LoadScaffoldGraphFromCheckpoint(GlobalData->outputPrefix,restartFromCheckpoint, TRUE);

    //  Dump stats on the loaded checkpoint
    //GeneratePlacedContigGraphStats(tmpBuffer,0);
    //GenerateScaffoldGraphStats(tmpBuffer,0);

    // shatter scaffolds if requested
    if (GlobalData->shatterLevel > 0) {
    	ShatterScaffoldsConnectedByLowWeight(stderr, ScaffoldGraph, GlobalData->shatterLevel, TRUE);
    }

    //  Useful for checking mate happiness on loading.  Currently only checks one scaffold.
    if (0) {
      vector<instrumentLIB>   libs;

      for (int32 i=0; i<GetNumDistTs(ScaffoldGraph->Dists); i++) {
        DistT *dptr = GetDistT(ScaffoldGraph->Dists, i);

        libs.push_back(instrumentLIB(i, dptr->mu, dptr->sigma, true));
      }

      for (int32 sID=287340; sID < GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds); sID++) {
        CIScaffoldT *scaffold = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, sID);

        fprintf(stderr, "ANALYZING SCAFFOLD %d\n", sID);

        if (scaffold->flags.bits.isDead == true)
          continue;

        instrumentSCF   A(scaffold);
        A.analyze(libs);
        A.report();

        exit(0);
      }
    }

    if (recomputeLeastSquaresOnLoad) {
      for (int32 sID=0; sID < GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds); sID++) {
        CIScaffoldT *scaffold = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, sID);

        if (scaffold->flags.bits.isDead == true)
          continue;

        if (true == LeastSquaresGapEstimates(ScaffoldGraph, GetCIScaffoldT(ScaffoldGraph->CIScaffolds, sID), LeastSquares_Cleanup | LeastSquares_Split))
          ScaffoldSanity(ScaffoldGraph, scaffold);
      }
    }
  }


  //  We DO want to flush unused unitigs/contigs at this point.  They're not in
  //  a scaffold, and possibly will never be used again (except as rocks/stones).
  //
  ScaffoldGraph->tigStore->flushCache();


  if ((runThisCheckpoint(restartFromLogical, CHECKPOINT_DURING_INITIAL_SCAFFOLDING) == true) &&
      (GlobalData->repeatRezLevel > 0)) {
    int ctme     = time(0);

    if(GlobalData->debugLevel > 0)
      DumpContigs(stderr,ScaffoldGraph, FALSE);

    // Transitive reduction of ContigGraph followed by construction of SEdges

    //  With markShakyBifurcations enabled.
    BuildUniqueCIScaffolds(ScaffoldGraph, TRUE, FALSE);

    CheckEdgesAgainstOverlapper(ScaffoldGraph->ContigGraph);

    //  Equivalent to TidyUpScaffolds().
    //
    for (int32 sID=0; sID < GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds); sID++) {
      CIScaffoldT *scaffold = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, sID);

      if (true == LeastSquaresGapEstimates(ScaffoldGraph, scaffold, LeastSquares_Cleanup | LeastSquares_Split))
        ScaffoldSanity(ScaffoldGraph, scaffold);
    }

    if (time(0) - ctme > 60 * 60)
      CheckpointScaffoldGraph(ckpNames[CHECKPOINT_DURING_INITIAL_SCAFFOLDING], "during initial scaffolding");
  }


  if ((runThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_INITIAL_SCAFFOLDING) == true) &&
      (GlobalData->repeatRezLevel > 0)) {

    //CheckAllTrustedEdges(ScaffoldGraph);

    {
      vector<CDS_CID_t>  rawEdges;

      BuildSEdges(rawEdges, FALSE);
      MergeAllGraphEdges(ScaffoldGraph->ScaffoldGraph, rawEdges, TRUE, FALSE);
    }

    //ScaffoldSanity(ScaffoldGraph);

    //  rocks is called inside of here
    //  checkpoints are written inside of here

    int iter     = 0;
    int iterMax  = 10;  //  MAX_OUTPUT_REZ_ITERATIONS
    int ctme     = time(0);
    int changed  = TRUE;

    fprintf(stderr,"** Running Level 1 Repeat Rez **\n");

    while ((changed) && (iter < iterMax)) {
      CheckEdgesAgainstOverlapper(ScaffoldGraph->ContigGraph);
      CheckCITypes(ScaffoldGraph);

      changed = RepeatRez(GlobalData->repeatRezLevel, GlobalData->outputPrefix);

      if (changed){
        CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);
        ScaffoldSanity(ScaffoldGraph);

        //  With markShakyBifurcations disabled.
        BuildUniqueCIScaffolds(ScaffoldGraph, FALSE, FALSE);

        CheckEdgesAgainstOverlapper(ScaffoldGraph->ContigGraph);

        for (int32 sID=0; sID < GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds); sID++) {
          CIScaffoldT *scaffold = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, sID);

          if (true == LeastSquaresGapEstimates(ScaffoldGraph, scaffold, LeastSquares_Cleanup | LeastSquares_Split))
            ScaffoldSanity(ScaffoldGraph, scaffold);
        }

        //CheckAllTrustedEdges(ScaffoldGraph);

        //  This shouldn't be necessary (RepeatRez() calling TidyUpScaffolds() should be doing it),
        //  but it is infrequent (at most iterMax=10 times).
        {
          vector<CDS_CID_t>  rawEdges;

          BuildSEdges(rawEdges, FALSE);
          MergeAllGraphEdges(ScaffoldGraph->ScaffoldGraph, rawEdges, TRUE, FALSE);
        }

        //  If we've been running for 2 hours, AND we've not just
        //  completed the last iteration, checkpoint.
        //
        if ((time(0) - ctme > 120 * 60) && (changed) && (iter+1 < iterMax)) {
          ctme = time(0);
          CheckpointScaffoldGraph(ckpNames[CHECKPOINT_DURING_INITIAL_SCAFFOLDING], "during initial scaffolding");
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

    CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_INITIAL_SCAFFOLDING], "after initial scaffolding");
  }
  //  else TidyUpScaffolds (ScaffoldGraph);


  //  We DO want to flush unused unitigs/contigs at this point.  They're not in
  //  a scaffold, and possibly will never be used again (except as rocks/stones).
  //
  ScaffoldGraph->tigStore->flushCache();


  if (runThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_1ST_SCAFF_MERGE) == true) {
    CleanupScaffolds(ScaffoldGraph,FALSE, NULLINDEX, FALSE);

    ScaffoldSanity(ScaffoldGraph);

    /* First we try to merge Scaffolds agressively */
    MergeScaffoldsAggressive(ScaffoldGraph, ckpNames[CHECKPOINT_DURING_1ST_SCAFF_MERGE], FALSE);
    CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr, "---Checking contig orders after MergeScaffoldsAggressive (1)\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif

#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_1ST_SCAFF_MERGE], "after 1st scaffold merge");
  }


  //  We DO want to flush unused unitigs/contigs at this point.  They're not in
  //  a scaffold, and possibly will never be used again (except as rocks/stones).
  //
  ScaffoldGraph->tigStore->flushCache();


  /*
    now that we are done with initial scaffold merge, we want to use the
    standard/default repeatRezLevel. Up to now, the value of preMergeRezLevel
    was in use if set on the command line
  */
  GlobalData->repeatRezLevel = repeatRezLevel;



  /* Now we throw stones */
  if ((runThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_STONES) == true) &&
      (GlobalData->stoneLevel > 0)) {

    // Convert single-contig scaffolds that are marginally unique back
    // to unplaced contigs so they might be placed as stones
    //
    //  If we removed any scaffolds, rebuild all the edges.
    //
    if ((GlobalData->demoteSingletonScaffolds == true) &&
        (DemoteSmallSingletonScaffolds() == true)) {
      vector<CDS_CID_t>  rawEdges;

      BuildSEdges(rawEdges, TRUE);
      MergeAllGraphEdges(ScaffoldGraph->ScaffoldGraph, rawEdges, TRUE, TRUE);
    }

    ScaffoldSanity(ScaffoldGraph);
    Throw_Stones(GlobalData->outputPrefix, GlobalData->stoneLevel, FALSE);

    //  If Throw_Stones splits scaffolds, rebuild edges.
    //  Throw_Stones should return if 'splitscaffolds > 0' so we can avoid
    //  rebuilding.
    {
      vector<CDS_CID_t>  rawEdges;

      BuildSEdges(rawEdges, TRUE);
      MergeAllGraphEdges(ScaffoldGraph->ScaffoldGraph, rawEdges, TRUE, TRUE);
    }

    ScaffoldSanity(ScaffoldGraph);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr, "---Checking contig orders after Throw_Stones\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif

#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_STONES], "after stone throwing");

    //GenerateLinkStats(ScaffoldGraph->CIGraph, "Stones", 0);
    //GeneratePlacedContigGraphStats("Stones", 0);
    //GenerateLinkStats(ScaffoldGraph->ContigGraph, "Stones", 0);
    //GenerateScaffoldGraphStats("Stones", 0);
  }


  if ((runThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_2ND_SCAFF_MERGE) == true) &&
      (GlobalData->stoneLevel > 0)) {

    ScaffoldSanity(ScaffoldGraph);

    MergeScaffoldsAggressive(ScaffoldGraph, ckpNames[CHECKPOINT_DURING_2ND_SCAFF_MERGE], FALSE);

    CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr, "---Checking contig orders after MergeScaffoldsAggressive (2)\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif

#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_2ND_SCAFF_MERGE], "after 2nd scaffold merge");
  }

  //  We DO want to flush unused unitigs/contigs at this point.  They're not in
  //  a scaffold, and possibly will never be used again (except as rocks/stones).
  //
  ScaffoldGraph->tigStore->flushCache();

  if ((runThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_FINAL_ROCKS) == true) &&
      (GlobalData->repeatRezLevel > 0)) {
    const int  MAX_EXTRA_ROCKS_ITERS = 5;
    int  iter = 0, extra_rocks;

    do {
      extra_rocks = Fill_Gaps(GlobalData->outputPrefix, GlobalData->repeatRezLevel, iter);
      fprintf(stderr, "Threw additional %d rocks on iter %d\n", extra_rocks, iter);

      //ScaffoldGraph->tigStore->flushCache();
    } while (extra_rocks > 1 && iter < MAX_EXTRA_ROCKS_ITERS);

    CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_FINAL_ROCKS], "after final rocks");
  }

  if ((runThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_PARTIAL_STONES) == true) &&
      (GlobalData->stoneLevel > 0)) {

    ScaffoldSanity (ScaffoldGraph);

    int partial_stones = Throw_Stones(GlobalData->outputPrefix, GlobalData->stoneLevel, TRUE);

    //  If throw_stones splits scaffolds, rebuild edges
    {
      vector<CDS_CID_t>  rawEdges;

      BuildSEdges(rawEdges, TRUE);
      MergeAllGraphEdges(ScaffoldGraph->ScaffoldGraph, rawEdges, TRUE, TRUE);
    }

    ScaffoldSanity (ScaffoldGraph);

    //ScaffoldGraph->tigStore->flushCache();

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

    CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_PARTIAL_STONES], "after partial stones");

    //GenerateLinkStats (ScaffoldGraph->CIGraph, "PStones", 0);
    //GeneratePlacedContigGraphStats ("PStones", 0);
    //GenerateLinkStats(ScaffoldGraph->ContigGraph, "PStones", 0);
    //GenerateScaffoldGraphStats ("PStones", 0);
  }

  if ((runThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_FINAL_CONTAINED_STONES) == true) &&
      (GlobalData->stoneLevel > 0)) {

    ScaffoldSanity (ScaffoldGraph);

    int contained_stones = Toss_Contained_Stones (GlobalData->outputPrefix, GlobalData->stoneLevel, 0);
    fprintf(stderr, "Threw %d contained stones\n", contained_stones);
    fprintf (stderr, "**** Finished Final Contained Stones level %d ****\n", GlobalData->stoneLevel);

    CleanupScaffolds (ScaffoldGraph, FALSE, NULLINDEX, FALSE);
    ScaffoldSanity (ScaffoldGraph);

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

    CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_FINAL_CONTAINED_STONES], "after final contained stones");

    //GenerateLinkStats (ScaffoldGraph->CIGraph, "CStones", 0);
    //GeneratePlacedContigGraphStats ("CStones", 0);
    //GenerateLinkStats(ScaffoldGraph->ContigGraph, "CStones", 0);
    //GenerateScaffoldGraphStats ("CStones", 0);
  }

  //  We DO want to flush unused unitigs/contigs at this point.  They're not in
  //  a scaffold, and possibly will never be used again (except as rocks/stones).
  //
  ScaffoldGraph->tigStore->flushCache();


  if (runThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_FINAL_CLEANUP) == true) {

    // Try to cleanup failed merges, and if we do, generate a checkpoint
    if(CleanupFailedMergesInScaffolds(ScaffoldGraph)){
      // This call deletes surrogate-only contigs that failed to merge
      if(CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, TRUE)){

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
        fprintf(stderr, "---Checking contig orders after final cleanup\n\n");
        CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
      }
      CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_FINAL_CLEANUP], "after final cleanup");
    }
  }


  if ((runThisCheckpoint(restartFromLogical, CHECKPOINT_AFTER_RESOLVE_SURROGATES) == true) &&
      (doResolveSurrogates > 0)) {

    resolveSurrogates(placeAllFragsInSinglePlacedSurros, cutoffToInferSingleCopyStatus);
    // Call resolve surrogate twice, this is necessary for finishing (closure) reads.
    // Consider a closure read and its two bounding reads, named left and right:
    //    If one (right) is placed in a unique region while the other (left) is in a surrogate itself, the closure read cannot be placed
    //    However, once the surrogate bounding read is placed (and fully incorporated which happens at the very end of resolveSurrogates)
    //    the closure read can be placed. 
    //    Therefore, we run resolve surrogates twice. 
    // Note that is closure reads are themselves mated, it may be necessary to do a third round of placement.  
    resolveSurrogates(placeAllFragsInSinglePlacedSurros, cutoffToInferSingleCopyStatus);
    
    CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_RESOLVE_SURROGATES], "after resolve surrogates");
  }

  //  This generates the 'rezlog/gapreads' file.  It's hugely
  //  expensive, usually dies on a negative variance assert, and as
  //  far as BPW knows, unused.
  //
  //Show_Reads_In_Gaps (GlobalData->outputPrefix);

  ComputeMatePairStatisticsRestricted(SCAFFOLD_OPERATIONS, GlobalData->minSamplesForOverride, "scaffold_final");
  ComputeMatePairStatisticsRestricted(CONTIG_OPERATIONS, GlobalData->minSamplesForOverride, "contig_final");

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
         
    if (ScaffoldGraph->gkpStore->gkStore_getFRGtoPLC(frag->read_iid) != 0) {
      AS_UID uid = getGatekeeperIIDtoUID(ScaffoldGraph->gkpStore, frag->read_iid, AS_IID_FRG);
      if (frag->contigID != -1) {
        ChunkInstanceT * ctg = GetGraphNode(ScaffoldGraph->ContigGraph, frag->contigID);            
        fprintf(stderr, "CLOSURE_READS: CLOSURE READ %s PLACED=%d CHAFF=%d SINGLETON=%d IN ASM type %c in SCF %d\n", AS_UID_toString(uid), frag->flags.bits.isPlaced, frag->flags.bits.isChaff, frag->flags.bits.isSingleton, frag->type, ctg->scaffoldID);
      }
    }
  }
#endif

  //  We DO want to flush unused unitigs/contigs at this point.  They're not in
  //  a scaffold, and possibly will never be used again (except as rocks/stones).
  //
  //  (This assumes that output doesn't load unitigs/contigs again)
  //
  ScaffoldGraph->tigStore->flushCache();

  SetCIScaffoldTLengths(ScaffoldGraph);

  if(generateOutput){
    CelamyAssembly(GlobalData->outputPrefix);

    MarkContigEdges();
    ComputeMatePairDetailedStatus();

    //  Note that OutputContigs partitions the tigStore, and closes ScaffoldGraph->tigStore.  The
    //  only operation valid after this function is CheckpointScaffoldGraph().

    OutputUnitigsFromMultiAligns();
    OutputContigsFromMultiAligns(outputFragsPerPartition);

    CheckpointScaffoldGraph(ckpNames[CHECKPOINT_AFTER_OUTPUT], "after output");
  }

  DestroyScaffoldGraph(ScaffoldGraph);

  delete GlobalData;

  fprintf(stderr,"* Bye *\n");

  exit(0);
}

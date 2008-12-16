
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

const char *mainid = "$Id: AS_CGW_main.c,v 1.65 2008-12-16 22:32:36 skoren Exp $";

static const char *usage =
"usage: %s [options] -g <GatekeeperStoreName> -o <OutputPath> <InputCGB.ext>\n"
"\n"
"   [-a]           align overlaps  (default)\n"
"   [-d]           DumpGraphs\n"
"   [-e <thresh>]  Microhet score probability cutoff\n"
"   [-g]           gkp Store path (required)\n"
"   [-h]           Don't fail on merge alignment failure\n"
"   [-i <thresh>]  Set max coverage stat for microhet determination of non-uniqueness (default -1)\n"
"   [-j <thresh>]  Set min coverage stat for definite uniqueness\n"
"   [-k <thresh>]  Set max coverage stat for possible uniqueness\n"
"   [-l <maxdegree> ]\n"
"   [-m <min>]     Number of mate samples to recompute an insert size, default is 100\n"
"   [-n] <scaffnum>  Starting scaffold number for stones\n"
"   [-o]           Output Name (required)\n"
"   [-q <cutoff>]  Transquality cutoff\n"
"   [-r <repeatRezlevel>]\n"
"   [-s <stoneLevel> ]\n"
"   [-v]           verbose\n"
"   [-x]           Dump stats on checkpoint load\n"
"   [-y]           Turn off Check for -R < -N, for restarting small assemblies\n"
"   [-z]           Turn on Check for Repeat Branch Pattern (demotes some unique unitigs to repeat)\n"
"\n"
"   [-A]           don't align overlaps (quality values reflect bayesian)\n"
"   [-C]           Don't cleanup scaffolds\n"
"   [-D <debugLevel>]\n"
"   [-E]           output overlap only contig edges\n"
"   [-F]           strongly enforce unique/repeat flag set in unitig, default if not set is to still allow those marked unique to be demoted due to Repeat Branch Pattern or being too small\n"
"   [-G]           Don't generate output (cgw or cam)\n"
"   [-H]           fail on merge alignment failure   (default)\n"
"   [-I]           ignore chaff unitigs\n"
"   [-J]           annotate celamy output with Bactig content (default = false)\n"
"   [-M]           don't do interleaved scaffold merging\n"
"   [-N <ckp>]     restart from checkpoint location ckp\n"
"   [-O [icC]]     i = immediate, c = continue, C = celamy output only\n"
"   [-P]                            Proto Output\n"
"   [-R <ckp>]     restart from checkpoint file ckp\n"
"   [-S <t>]       place all frags in singly-placed surrogates if at least fraction\n"
"                  <t> can be placed; two special cases:\n"
"                    if <t> = -1, place all frags in singly-placed surrogates aggressively\n"
"                                 (which really mean t = 0.0, but triggers a better algorithm)\n"
"                    if <t> =  0, do not resolve surrogate fragments\n"
"   [-X <Estimated number of nodes>]\n"
"   [-Y <Estimated number of edges>]\n"
"   [-Z]           Don't demote singleton scaffolds\n";


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include <math.h>
#include <time.h>

#define DEBUG_MERGE_SCAF FALSE

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
#include "AS_ALN_forcns.h"
#include "Instrument_CGW.h"
#include "AS_CGW_EdgeDiagnostics.h"
#include "Checkpoints_CGW.h"
#include "fragmentPlacement.h"  //  for resolveSurrogates()

void DemoteUnitigsWithRBP(FILE *stream, GraphCGW_T *graph);
void RemoveSurrogateDuplicates(void);

FILE *  File_Open (const char * Filename, const char * Mode, int exitOnFailure);


int main(int argc, char *argv[]){
  int checkRepeatBranchPattern = FALSE;
  int preMergeRezLevel = -1;
  int generateOutput = 1;
  int annotateUnitigs = 0;
  int doInterleavedScaffoldMerging = 1;
  int debugLevel = 0;
  int repeatRezLevel = 0;
  int stoneLevel = 0;
  int starting_stone_scaffold = 0;
  int numNodes = 1024;
  int numEdges = 1024;
  int minSamplesForOverride = 100;
  int ignoreChaffUnitigs = 0;
  int failOn_NoOverlapFound = 0;
  int geneOutput = 1;  // output simulated coordinates
  int dumpGraphs = 0;
  int outputOverlapOnlyContigEdges = 0;
  int doRezOnContigs = 1;
  int immediateOutput = 0;
  int camFileOnly = 0;
  int performCleanupScaffolds = 1;
  int alignOverlaps = 1; // if true (set to true by -a, false -A) recompute all overlaps supplied by cgb
  int dumpStatsOnRestart = FALSE;
  int32 restartFromCheckpoint = NULLINDEX;
  int32 restartFromLogicalCheckpoint = NULLINDEX;
  int demoteSingletonScaffolds = TRUE;
  Global_CGW *data;
  FILE *infp = NULL;
  float transQualityCutoff = 0.1; // quality cutoff for TransChunkEdges
  float cgbMicrohetProb = 1.e-05;      // scores less than this are considered repeats
  float cgbApplyMicrohetCutoff = -1; // This basically turns it off, unless enabled
  float cgbUniqueCutoff = CGB_UNIQUE_CUTOFF;
  float cgbDefinitelyUniqueCutoff = CGB_UNIQUE_CUTOFF;
  int maxDegree = 3; // maximum edges to keep for 'nonUnique' nodes
  int maxDegreeUnique = 30; // maximum edges to keep for 'Unique' nodes
  int setGatekeeperStore = 0;
  char *outputPath = NULL;
  int checkpointChecker = 1;
  int    doResolveSurrogates               = 1;      //  resolveSurrogates
  int    placeAllFragsInSinglePlacedSurros = 0;      //  resolveSurrogates
  double cutoffToInferSingleCopyStatus     = 0.666;  //  resolveSurrogates
  int    allowDemoteMarkedUnitigs          = TRUE;      // allow toggled unitigs to be demoted to be repeat if they were marked unique

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
  ContigOrientChecker * coc;
  coc = CreateContigOrientChecker();
  assert(coc != NULL);
#endif

  GlobalData  = data = CreateGlobal_CGW();  
  GlobalData->stderrc = stderr;
  GlobalData->aligner=DP_Compare;
  GlobalData->closureReads = NULL;

  argc = AS_configure(argc, argv);

  { /* Parse the argument list using "man 3 getopt". */
    int ch,errflg=0;
    char  * p;

    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
                                    "ade:f:g:hi:j:k:l:m:n:o:p:q:r:s:vxyzACD:EFGHIJK:L:N:MO:PQR:S:V:X:Y:Z")) != EOF)){
      switch(ch) {
        case 'O':
          fprintf(GlobalData->stderrc,"* Immediate output specified: arg %s\n", optarg);
          switch(*optarg){
            case 'C':
              immediateOutput = 3;
              camFileOnly = 1;
              generateOutput = 0;
              break;
            case 'i':
              immediateOutput = 2;
              break;
            case 'c':
              immediateOutput = 1;
              break;
            case 'm': //MergeScaffolds no cgw file
              immediateOutput = -1;
              camFileOnly = 1;
              generateOutput = 0;
              break;
            default:
              break;
          }
        case 'a':
          alignOverlaps = 1;
          break;
        case 'A':
          alignOverlaps = 0;
          break;
        case 'C':
          performCleanupScaffolds = 0;
          fprintf(GlobalData->stderrc,"* -C ==> No Cleanup Scaffolds!\n");
          fflush(stderr);
          break;
        case 'g':
          strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
          setGatekeeperStore = 1;
          break;
        case 'h':
          failOn_NoOverlapFound = 0;
          break;
        case 'G':
          generateOutput = 0;
          break;
        case 'H':
          failOn_NoOverlapFound = 1;
          break;
        case 'I':
          ignoreChaffUnitigs = 1;
          break;
        case 'J':
          annotateUnitigs = 1;
          break;
        case 'M':
          doInterleavedScaffoldMerging = 0;
          break;
        case 'K':
          switch(optarg[0]) {
            case 'D':
              GlobalData->aligner=DP_Compare;
              break;
            case 'L':
              GlobalData->aligner=Local_Overlap_AS_forCNS;
              break;
            case 'A':
              GlobalData->aligner=Affine_Overlap_AS_forCNS;
              break;
            default:
              GlobalData->aligner=DP_Compare;
              fprintf(stderr,"Unrecognized value for option -%c (%s)",optopt,optarg);
          }
          break;
        case 'q':
          transQualityCutoff = (float)atof(optarg);
          fprintf(GlobalData->stderrc,"* transQualityCutoff set to %f\n", transQualityCutoff);
          break;
        case 'i':
          cgbApplyMicrohetCutoff = (float)atof(optarg);
          fprintf(GlobalData->stderrc,"* cgbApplyMicrohetCutoff set to %f\n", cgbApplyMicrohetCutoff);
          break;
        case 'j':
          cgbUniqueCutoff = (float)atof(optarg);
          fprintf(GlobalData->stderrc,"* cgbUniqueCutoff set to %f\n", cgbUniqueCutoff);
          break;
        case 'k':
          cgbDefinitelyUniqueCutoff = (float)atof(optarg);
          fprintf(GlobalData->stderrc,"* cgbDefinitelyUniqueCutoff set to %f\n", cgbDefinitelyUniqueCutoff);
          break;
        case 'l':
          maxDegree = atoi(optarg);
          fprintf(GlobalData->stderrc,"* maxDegree set to %d\n", maxDegree);
          break;
        case 'm':
          minSamplesForOverride = atoi(optarg);
          fprintf(GlobalData->stderrc,"* minSamplesForOverride set to %d\n", minSamplesForOverride);
          break;
        case 'n':
          starting_stone_scaffold = atoi(optarg);
          fprintf (GlobalData -> stderrc,"* Starting Scaffold for Stones set to %d\n",
                   starting_stone_scaffold);
          break;
        case 'L':
          maxDegreeUnique = atoi(optarg);
          fprintf(GlobalData->stderrc,"* maxDegreeUnique set to %d\n", maxDegreeUnique);
          break;
        case 'X':
          numNodes = atoi(optarg);
          fprintf(GlobalData->stderrc,"* numNodes set to %d\n", numNodes);
          break;
        case 'y':
          checkpointChecker = 0;
          break;
        case 'Y':
          numEdges = atoi(optarg);
          fprintf(GlobalData->stderrc,"* numEdges set to %d\n", numEdges);
          break;

        case 'r':
          repeatRezLevel = atoi(optarg);
          fprintf(GlobalData->stderrc,"* repeatRezLevel set to %d\n", repeatRezLevel);
          break;
        case 's':
          stoneLevel = atoi(optarg);
          fprintf(GlobalData->stderrc,"* stoneLevel set to %d\n", stoneLevel);
          break;
        case 'S':
          doResolveSurrogates               = 1;
          cutoffToInferSingleCopyStatus     = atof(optarg);
          placeAllFragsInSinglePlacedSurros = 0;

          if (cutoffToInferSingleCopyStatus == 0.0)
            doResolveSurrogates               = 0;

          if (cutoffToInferSingleCopyStatus < 0) {
            cutoffToInferSingleCopyStatus     = 0.0;
            placeAllFragsInSinglePlacedSurros = 1;
          }
          if (doResolveSurrogates)
            fprintf(GlobalData->stderrc,"* resolveSurrogates: -S %f%s\n",
                    cutoffToInferSingleCopyStatus, placeAllFragsInSinglePlacedSurros ? " -1" : "");
          else
            fprintf(GlobalData->stderrc,"* resolveSurrogates: DISABLED\n");
          break;
        case 'D':
          debugLevel = atoi(optarg);
          fprintf(GlobalData->stderrc,"* debugLevel set to %d\n", debugLevel);
          break;
        case 'o':
          outputPath = strdup(optarg);
          break;
        case 'p':
          preMergeRezLevel = atoi(optarg);
          if(preMergeRezLevel < 0){
            fprintf(GlobalData->stderrc,"* pre-first-scaffold-merge repeatRezLevel must be >= 0\n");
            fprintf(stderr,"* pre-first-scaffold-merge repeatRezLevel must be >= 0\n");
            exit(1);
          }
          fprintf(GlobalData->stderrc,"* pre-first-scaffold-merge repeatRezLevel set to %d\n", preMergeRezLevel);
          break;
        case 'd':
          dumpGraphs = 1;
          fprintf(GlobalData->stderrc,"* dumpGraphs set to 1\n");
          break;
        case 'e':
          cgbMicrohetProb = atof(optarg);
          fprintf(GlobalData->stderrc,"* Microhet probability threshhold set to %g\n",cgbMicrohetProb);
          break;
        case 'E':
          outputOverlapOnlyContigEdges = 1;
          fprintf(GlobalData->stderrc,"* Output overlap only  contig edges set to TRUE\n");
          break;
        case 'P':
          fprintf(stderr, "-P default.\n");
          break;
        case 'R':
          restartFromCheckpoint = atoi(optarg);
          fprintf(GlobalData->stderrc,"* Restarting from checkpoint %d\n", restartFromCheckpoint);
          break;
        case 'N':
          restartFromLogicalCheckpoint = atoi(optarg);
          fprintf(GlobalData->stderrc,"* Restarting from logical checkpoint %d\n", restartFromLogicalCheckpoint);
          break;
        case 'x':
          dumpStatsOnRestart = TRUE;
          fprintf(GlobalData->stderrc,"* Dumping Stats on restart\n");
          break;
        case 'Z':
          fprintf(GlobalData->stderrc,"* DON'T demote singleton scaffolds\n");
          demoteSingletonScaffolds = FALSE;
          break;
        case 'z':
          fprintf(GlobalData->stderrc,"* Check for Repeat Branch Pattern\n");
          checkRepeatBranchPattern = TRUE;
          break;
        case 'F':
          fprintf(GlobalData->stderrc, "* Allow Demote Unique Unitigs\n");
          allowDemoteMarkedUnitigs = FALSE;
          break;

        case '?':
          fprintf(GlobalData->stderrc,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if ((argc - optind < 1 ) || (setGatekeeperStore == 0) || (outputPath == NULL)) {
      fprintf(GlobalData->stderrc,"* argc = %d optind = %d setGatekeeperStore = %d outputPath = %s\n",
              argc, optind,setGatekeeperStore, outputPath);
      fprintf(GlobalData->stderrc, usage, argv[0]);
      exit (EXIT_FAILURE);
    }


    if(cgbDefinitelyUniqueCutoff < cgbUniqueCutoff){
      cgbDefinitelyUniqueCutoff = cgbUniqueCutoff;
      fprintf(GlobalData->stderrc,"* cgbDefinitelyUniqueCutoff set to %f\n", cgbUniqueCutoff);
    }
    data->annotateUnitigs = annotateUnitigs;
    data->doInterleavedScaffoldMerging = doInterleavedScaffoldMerging;
    data->maxSequencedbSize = MAX_SEQUENCEDB_SIZE;
    data->maxDegreeUnique = maxDegreeUnique;
    data->maxDegree = maxDegree;
    data->transQualityCutoff = transQualityCutoff;
    data->cgbUniqueCutoff = cgbUniqueCutoff;
    if(preMergeRezLevel >= 0){
      data->repeatRezLevel = preMergeRezLevel;
    } else {
      data->repeatRezLevel = repeatRezLevel;
    }
    data->stoneLevel     = stoneLevel;
    data->debugLevel     = debugLevel;
    data->failOn_NoOverlapFound = failOn_NoOverlapFound;
    data->ignoreChaffUnitigs = ignoreChaffUnitigs;
    data->performCleanupScaffolds = performCleanupScaffolds;
    data->cgbUniqueCutoff = cgbUniqueCutoff;
    data->cgbApplyMicrohetCutoff = cgbApplyMicrohetCutoff;
    data->cgbDefinitelyUniqueCutoff = cgbDefinitelyUniqueCutoff;
    data->cgbMicrohetProb = cgbMicrohetProb;
    data -> starting_stone_scaffold = starting_stone_scaffold;
    data->allowDemoteMarkedUnitigs = allowDemoteMarkedUnitigs;

    if(optind < argc)
      {
        char  filepath[2048];

	strcpy(data->Input_File_Name, argv[optind]);
	infp = File_Open (data->Input_File_Name, "r", TRUE);     // frg file
	AssertPtr(infp);
	strcpy(data->File_Name_Prefix, outputPath);
        if (generateOutput) {
          sprintf(filepath, "%s.cgw", outputPath);
          data->cgwfp = File_Open(filepath, "w", TRUE);

          sprintf(filepath, "%s.cgw_contigs", outputPath);
          data->ctgfp = File_Open(filepath, "w", TRUE);

          sprintf(filepath, "%s.cgw_scaffolds", outputPath);
          data->scffp = File_Open(filepath, "w", TRUE);
        }

        sprintf(filepath, "%s.timing", outputPath);
        data->timefp = File_Open(filepath, "a", TRUE); // timing file

	{
	  time_t t;
	  t = time(0);
	  fprintf(data->timefp,"\n\n>>>>*************************************************************************<<<<\n");
	  fprintf(data->timefp,"* Restarting from checkpoint %d at time %s\n", restartFromCheckpoint,ctime(&t));
	}

	//	optind++;
      }
    /* End of command line parsing */
  }

  if(!doRezOnContigs && (stoneLevel > 0)){
    fprintf(GlobalData->stderrc,"#### Stone Throwing and Gap Walking can only be run in conjunction with the -C option ####\n");
    fprintf(GlobalData->stderrc,".....exiting.....\n");
    exit(1);
  }

  if (cutoffToInferSingleCopyStatus > 1.0) {
    fprintf(GlobalData->stderrc,"#### -S must be between 0.0 and 1.0.\n");
    exit(1);
  }

  if(immediateOutput > 0 &&
     restartFromCheckpoint <= CHECKPOINT_AFTER_BUILDING_CIGRAPH){
    fprintf(GlobalData->stderrc,"* Can't generate output unless -R n, n > %d is specified...ignored\n",
	    CHECKPOINT_AFTER_BUILDING_CIGRAPH);
    immediateOutput = 0;
  }

  fprintf(GlobalData->stderrc,"* Scaffolding will be performed on %s\n",
	  (doRezOnContigs?"Contigs":"CIs"));

  data->outputCalculatedOffsets = (geneOutput == 0);

  if(checkpointChecker && restartFromCheckpoint < restartFromLogicalCheckpoint &&
     restartFromLogicalCheckpoint !=  CHECKPOINT_BEFORE_FINAL_CLEANUP){
    fprintf(GlobalData->stderrc,"* Logical Checkpoint (%d)  MUST be no greater than restart checkpoint (%d)\n", restartFromLogicalCheckpoint, restartFromCheckpoint);
    exit(1);
  }

  if(restartFromLogicalCheckpoint == NULLINDEX){
    restartFromLogicalCheckpoint = restartFromCheckpoint;
  }


  if(restartFromLogicalCheckpoint < CHECKPOINT_AFTER_READING_INPUT){
    // create the checkpoint from scratch
    ScaffoldGraph = CreateScaffoldGraph(doRezOnContigs, data->File_Name_Prefix, numNodes, numEdges);

    ScaffoldGraph->alignOverlaps = alignOverlaps;
    ScaffoldGraph->doRezOnContigs = doRezOnContigs;
    ScaffoldGraph->checkPointIteration = 0;

    ProcessInput(data, optind, argc, argv);    // Handle the rest of the first file

    LoadDistData();
  } else {
    // Load the checkpoint from a file
    ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(data->File_Name_Prefix,restartFromCheckpoint, !immediateOutput);

    if(dumpStatsOnRestart &&
       restartFromCheckpoint >= CHECKPOINT_AFTER_BUILDING_AND_CLEANING_SCAFFOLDS){
      char tmpBuffer[2048];

      fprintf(GlobalData->stderrc," Dumping Stats from Checkpoint ");
      sprintf(tmpBuffer,"LoadCKP%d",restartFromCheckpoint);
      GeneratePlacedContigGraphStats(tmpBuffer,0);
      GenerateScaffoldGraphStats(tmpBuffer,0);
    }
  }
  #warning SK - adding closure info from file now, it should be in GKP store
  if (GlobalData->closureReads == NULL) {
   GlobalData->closureReads = CreateScalarHashTable_AS();
   int line_len = ( 16 * 1024 * 1024);
   char *currLine = safe_malloc(sizeof(char)*line_len);
   char fileName[1024];
   sprintf(fileName, "%s.closureEdges", GlobalData->Gatekeeper_Store_Name);
   errno = 0;
   FILE *file = fopen(fileName, "r");
   if (errno) {
      errno = 0;
   } else {
      while (fgets(currLine, line_len-1, file) != NULL) {
         AS_UID read = AS_UID_lookup(currLine, NULL);
         InsertInHashTable_AS(GlobalData->closureReads, getGatekeeperUIDtoIID(ScaffoldGraph->gkpStore, read, NULL), 0, 1, 0);
      }
      fclose(file);
   }
   safe_free(currLine);
  }

  clearCacheSequenceDB(ScaffoldGraph->sequenceDB);

  if(restartFromLogicalCheckpoint < CHECKPOINT_AFTER_BUILDING_CIGRAPH) {

    // Split chimeric unitigs
    if( restartFromLogicalCheckpoint < CHECKPOINT_AFTER_UNITIG_SPLITTING ){
      fprintf(GlobalData->stderrc,"* Splitting chimeric input unitigs\n");
      fflush(GlobalData->stderrc);
      SplitInputUnitigs(ScaffoldGraph);
      clearCacheSequenceDB(ScaffoldGraph->sequenceDB);
    }

    if(restartFromLogicalCheckpoint < CHECKPOINT_AFTER_BUILDING_EDGES){
      ComputeMatePairStatisticsRestricted(UNITIG_OPERATIONS,
                                          minSamplesForOverride,
                                          "unitig_initial");

      fprintf(data->stderrc,"** Before BUILDCIEDGES **\n");
      if(GlobalData->debugLevel > 0){
        DumpChunkInstances(data->stderrc, ScaffoldGraph, FALSE, FALSE, FALSE, FALSE);
      }

      // Allocate space for edges
      ReallocGraphEdges(ScaffoldGraph->CIGraph, numEdges);

      //    BuildCIEdges(ScaffoldGraph);
      BuildGraphEdgesDirectly(ScaffoldGraph->CIGraph);
    }  //  checkpoint < CHECKPOINT_AFTER_BUILDING_EDGES

    // Compute all overlaps implied by mate links between pairs of unique unitigs
    ComputeOverlaps( ScaffoldGraph->CIGraph, TRUE, alignOverlaps);

    clearCacheSequenceDB(ScaffoldGraph->sequenceDB);

    CheckEdgesAgainstOverlapper(ScaffoldGraph->CIGraph);
    CheckSurrogateUnitigs();

    MergeAllGraphEdges(ScaffoldGraph->CIGraph, FALSE);

    CheckEdgesAgainstOverlapper(ScaffoldGraph->CIGraph);
    CheckSurrogateUnitigs();

    /* Mark some Unitigs/Chunks/CIs as repeats based on overlaps GRANGER 2/2/07 */

    if(checkRepeatBranchPattern){
      DemoteUnitigsWithRBP(GlobalData->stderrc, ScaffoldGraph->CIGraph);
    }

    /* At this Point we've constructed the CIGraph */

    /*    if(dumpGraphs)
	  DumpScaffoldGraph(ScaffoldGraph);*/

    if(GlobalData->debugLevel > 0){
      DumpChunkInstances(data->stderrc, ScaffoldGraph, FALSE, FALSE, FALSE, FALSE);
      CheckEdgesAgainstOverlapper(ScaffoldGraph->CIGraph);
    }

    /* Insert the guide edges into the graph AFTER merging */
    //    InsertGuideEdges(ScaffoldGraph->CIGraph);
    {
      /* Looks like we never keep BOTH the CIGraph overlaps and the ContigGraph overlaps */
      ChunkOverlapperT *tmp = ScaffoldGraph->ContigGraph->overlapper;
      ScaffoldGraph->ContigGraph->overlapper = ScaffoldGraph->CIGraph->overlapper;
      ScaffoldGraph->CIGraph->overlapper = tmp;
      ScaffoldGraph->RezGraph = ScaffoldGraph->ContigGraph;
    }

    clearCacheSequenceDB(ScaffoldGraph->sequenceDB);

    BuildInitialContigs(ScaffoldGraph);
    if(GlobalData->debugLevel){
      CheckEdgesAgainstOverlapper(ScaffoldGraph->ContigGraph);
      CheckSurrogateUnitigs();
    }

    clearCacheSequenceDB(ScaffoldGraph->sequenceDB);

    CheckCIScaffoldTs(ScaffoldGraph);
    if(GlobalData->debugLevel > 0){
      CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
      CheckSurrogateUnitigs();
    }
  }  //  checkpoint < CHECKPOINT_AFTER_BUILDING_CIGRAPH


  // Build scaffolds and do rocks

  if(immediateOutput == 0 &&
     ((restartFromLogicalCheckpoint < CHECKPOINT_AFTER_BUILDING_AND_CLEANING_SCAFFOLDS) &&
      data->repeatRezLevel > 0)){
    int skipInitialScaffolds = 0;

    //#include obsolete/rat_lbac_reactivation

    if(GlobalData->debugLevel > 0)
      DumpContigs(data->stderrc,ScaffoldGraph, FALSE);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);

    skipInitialScaffolds = (restartFromLogicalCheckpoint >= CHECKPOINT_AFTER_BUILDING_SCAFFOLDS);

    fprintf(GlobalData->stderrc,"**** Running BuildScaffoldsFromFirstPriniciples ****\n");
    fprintf(GlobalData->stderrc,"**** with skipInitialScaffolds = %d\n", skipInitialScaffolds);

    BuildScaffoldsFromFirstPriniciples(ScaffoldGraph, skipInitialScaffolds);  // rocks is called inside of here

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    if(GlobalData->debugLevel > 0)
      DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph, FALSE);

    //  The checkpoints are written in BuildScaffoldsFromFirstPriniciples(), either for
    //  AFTER_BUILDING_SCAFFOLDS or AFTER_BUILDING_AND_CLEANING_SCAFFOLDS, depending on
    //  (GlobalData->repeatRezLevel > 0).

  }


  if((immediateOutput == -1) || (immediateOutput <= 1 &&
                                 (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_1ST_SCAFF_MERGE))){
    fprintf(GlobalData->stderrc,"*** CS before SM ***\n");
    fflush(GlobalData->stderrc);

    if(GlobalData->debugLevel > 0){
      fprintf(GlobalData->stderrc,"**** BEFORE cleanupScaffolds ***\n");
      DumpCIScaffolds(GlobalData->stderrc, ScaffoldGraph, FALSE);
    }
    CleanupScaffolds(ScaffoldGraph,FALSE, NULLINDEX, FALSE);

    if(GlobalData->debugLevel > 0){
      fprintf(GlobalData->stderrc,"**** AFTER cleanupScaffolds ***\n");
      DumpCIScaffolds(GlobalData->stderrc, ScaffoldGraph, FALSE);
    }

    CheckCIScaffoldTs(ScaffoldGraph);

    /* First we try to merge Scaffolds agressively */
    MergeScaffoldsAggressive(ScaffoldGraph, CHECKPOINT_BEFORE_1ST_SCAFF_MERGE, DEBUG_MERGE_SCAF);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
    CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr,
            "---Checking contig orders after MergeScaffoldsAggressive (1)\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    fprintf(data->stderrc, "Checkpoint %d written after 1st Scaffold Merge\n",ScaffoldGraph->checkPointIteration);
    fprintf(data->timefp,  "Checkpoint %d written after 1st Scaffold Merge\n",ScaffoldGraph->checkPointIteration);
    CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_1ST_SCAFF_MERGE+1);
  } // No immediate output

  clearCacheSequenceDB(ScaffoldGraph->sequenceDB);

  /*
    now that we are done with initial scaffold merge, we want to use the
    standard/default repeatRezLevel. Up to now, the value of preMergeRezLevel
    was in use if set on the command line
  */
  data->repeatRezLevel = repeatRezLevel;

  // Convert single-contig scaffolds that are marginally unique back
  // to unplaced contigs so they might be placed as stones
  //
  if (demoteSingletonScaffolds)
    DemoteSmallSingletonScaffolds();

  /* Now we throw stones */
  if  (((immediateOutput == -1) ||
        ((immediateOutput <= 1) &&
         (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_STONES))) &&
       (GlobalData->stoneLevel > 0)) {
    fprintf(GlobalData->stderrc, "**** Running Stone Throwing level %d ****\n",
            GlobalData->stoneLevel);

    CheckCIScaffoldTs(ScaffoldGraph);

    Throw_Stones(GlobalData->File_Name_Prefix, GlobalData->stoneLevel, FALSE);

    CheckCIScaffoldTs(ScaffoldGraph);

    //  Throw_Stones does this now.
    //Force_Increasing_Variances();

    fprintf(GlobalData->stderrc, "**** Finished Stone Throwing level %d ****\n",
            GlobalData->stoneLevel);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);

    //  Throw_Stones does this now.
    //CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr, "---Checking contig orders after Throw_Stones\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif

#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    fprintf(data->stderrc, "Checkpoint %d written after Stone Throwing and CleanupScaffolds\n", ScaffoldGraph->checkPointIteration);
    fprintf(data->timefp,  "Checkpoint %d written after Stone Throwing and CleanupScaffolds\n", ScaffoldGraph->checkPointIteration);
    CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_STONES+1);

    GenerateLinkStats(ScaffoldGraph->CIGraph, "Stones", 0);
    GeneratePlacedContigGraphStats("Stones", 0);
    GenerateLinkStats(ScaffoldGraph->ContigGraph, "Stones", 0);
    GenerateScaffoldGraphStats("Stones", 0);
  }



  if(immediateOutput == 0 &&
     (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_2ND_SCAFF_MERGE)){
    fprintf(GlobalData->stderrc,"* Before Final MergeScaffoldsAggressive\n");
    CheckCIScaffoldTs(ScaffoldGraph);

    MergeScaffoldsAggressive(ScaffoldGraph, CHECKPOINT_BEFORE_2ND_SCAFF_MERGE, DEBUG_MERGE_SCAF);

    ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
    CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
    fprintf(stderr,
            "---Checking contig orders after MergeScaffoldsAggressive (2)\n\n");
    CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
    ResetContigOrientChecker(coc);
    AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

    fprintf(data->stderrc, "Checkpoint %d written after 2nd Aggressive Scaffold Merge\n",ScaffoldGraph->checkPointIteration);
    fprintf(data->timefp,  "Checkpoint %d written after 2nd Aggressive Scaffold Merge\n",ScaffoldGraph->checkPointIteration);
    CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_2ND_SCAFF_MERGE+1);
  }  // No immediate output


  clearCacheSequenceDB(ScaffoldGraph->sequenceDB);

  if  (immediateOutput == 0
       && (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_FINAL_ROCKS)
       && GlobalData -> repeatRezLevel > 0)
    {
      const int  MAX_EXTRA_ROCKS_ITERS = 5;
      int  iter = 0, extra_rocks;

      fprintf (GlobalData -> stderrc,
               "** Running Final Round of RepeatRez (Rocks) **\n");
      do
        {
          extra_rocks = Fill_Gaps (GlobalData, GlobalData -> File_Name_Prefix,
                                   GlobalData -> repeatRezLevel, iter);
          fprintf (GlobalData -> stderrc, "Threw additional %d rocks on iter %d\n",
                   extra_rocks, iter);

          clearCacheSequenceDB(ScaffoldGraph->sequenceDB);
        }  while  (extra_rocks > 1 && iter < MAX_EXTRA_ROCKS_ITERS);

      fprintf(data->stderrc, "Checkpoint %d written after Final Rocks\n", ScaffoldGraph -> checkPointIteration);
      fprintf(data->timefp,  "Checkpoint %d written after Final Rocks\n", ScaffoldGraph -> checkPointIteration);
      CheckpointScaffoldGraph (ScaffoldGraph, CHECKPOINT_BEFORE_FINAL_ROCKS+1);
    }


  if  (immediateOutput == 0
       && (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_PARTIAL_STONES)
       && GlobalData -> stoneLevel > 0)
    {
      int  partial_stones;

      fprintf (GlobalData -> stderrc, "* Starting Partial_Stones\n");

      CheckCIScaffoldTs (ScaffoldGraph);

      partial_stones = Throw_Stones (GlobalData -> File_Name_Prefix,
                                     GlobalData -> stoneLevel,
                                     TRUE);

      CheckCIScaffoldTs (ScaffoldGraph);
      ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);

      //  We used to CleanupScaffolds() here, but Throw_Stones now
      //  does that inline.

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
      fprintf(data->stderrc, "Checkpoint %d written after Partial Stones\n", ScaffoldGraph->checkPointIteration);
      fprintf(data->timefp,  "Checkpoint %d written after Partial Stones\n", ScaffoldGraph->checkPointIteration);
      CheckpointScaffoldGraph (ScaffoldGraph, CHECKPOINT_BEFORE_PARTIAL_STONES+1);

      GenerateLinkStats (ScaffoldGraph -> CIGraph, "PStones", 0);
      GeneratePlacedContigGraphStats ("PStones", 0);
      GenerateLinkStats(ScaffoldGraph -> ContigGraph, "PStones", 0);
      GenerateScaffoldGraphStats ("PStones", 0);
    }


  if  (immediateOutput == 0
       && (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_FINAL_CONTAINED_STONES)
       && GlobalData -> stoneLevel > 0)
    {
      int  contained_stones;

      fprintf (GlobalData -> stderrc, "* Starting Final Contained Stones\n");

      CheckCIScaffoldTs (ScaffoldGraph);
      contained_stones
        = Toss_Contained_Stones (GlobalData -> File_Name_Prefix,
                                 GlobalData -> stoneLevel, 0);
      ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
      CheckCIScaffoldTs (ScaffoldGraph);
      fprintf (GlobalData -> stderrc,
               "**** Finished Final Contained Stones level %d ****\n",
               GlobalData -> stoneLevel);

      CleanupScaffolds (ScaffoldGraph, FALSE, NULLINDEX, FALSE);
      clearCacheSequenceDB(ScaffoldGraph->sequenceDB);

      fprintf (stderr, "Threw %d contained stones\n", contained_stones);

      // Remove copies of surrogates which are placed multiple times in the same place in a contig

      RemoveSurrogateDuplicates();

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
      fprintf(stderr,
              "---Checking contig orders after contained_stones\n\n");
      CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
#ifdef CHECK_CONTIG_ORDERS_INCREMENTAL
      ResetContigOrientChecker(coc);
      AddAllScaffoldsToContigOrientChecker(ScaffoldGraph, coc);
#endif

      fprintf(data->stderrc, "Checkpoint %d written after Final Contained Stones\n", ScaffoldGraph->checkPointIteration);
      fprintf(data->timefp,  "Checkpoint %d written after Final Contained Stones\n", ScaffoldGraph->checkPointIteration);
      CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_FINAL_CONTAINED_STONES+1);

      GenerateLinkStats (ScaffoldGraph -> CIGraph, "CStones", 0);
      GeneratePlacedContigGraphStats ("CStones", 0);
      GenerateLinkStats(ScaffoldGraph -> ContigGraph, "CStones", 0);
      GenerateScaffoldGraphStats ("CStones", 0);
    }


  if(immediateOutput == 0 &&
     (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_FINAL_CLEANUP)){
    fprintf(GlobalData->stderrc,"* Before CleanupFailedMerges\n");


    // Try to cleanup failed merges, and if we do, generate a checkpoint
    if(CleanupFailedMergesInScaffolds(ScaffoldGraph)){

      fprintf(data->stderrc, "Checkpoint %d written after CleanupFailedMergesInScaffolds\n",ScaffoldGraph->checkPointIteration);
      fprintf(data->timefp,  "Checkpoint %d written after CleanupFailedMergesInScaffolds\n",ScaffoldGraph->checkPointIteration);
      CheckpointScaffoldGraph(ScaffoldGraph, -1);

      ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);

      // This call deletes surrogate-only contigs that failed to merge
      if( CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, TRUE)){

        fprintf(data->stderrc, "Checkpoint %d written after DeleteAllSurrogateContigsFromFailedMerges\n",ScaffoldGraph->checkPointIteration);
        fprintf(data->timefp,  "Checkpoint %d written after DeleteAllSurrogateContigsFromFailedMerges\n",ScaffoldGraph->checkPointIteration);
        CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_FINAL_CLEANUP+1);

#if defined(CHECK_CONTIG_ORDERS) || defined(CHECK_CONTIG_ORDERS_INCREMENTAL)
        fprintf(stderr, "---Checking contig orders after final cleanup\n\n");
        CheckAllContigOrientationsInAllScaffolds(ScaffoldGraph, coc, POPULATE_COC_HASHTABLE);
#endif
        ValidateAllContigEdges(ScaffoldGraph, FIX_CONTIG_EDGES);
      }
    }
  }


  if(immediateOutput == 0 &&
     (restartFromLogicalCheckpoint <= CHECKPOINT_BEFORE_RESOLVE_SURROGATES) &&
     (doResolveSurrogates > 0)) {
    fprintf(GlobalData->stderrc,"* Before resolveSurrogates (-S=%f -1=%d)\n",
            cutoffToInferSingleCopyStatus, placeAllFragsInSinglePlacedSurros);

    resolveSurrogates(placeAllFragsInSinglePlacedSurros, cutoffToInferSingleCopyStatus);

    fprintf(data->stderrc, "Checkpoint %d written after resolveSurrogates\n", ScaffoldGraph->checkPointIteration);
    fprintf(data->timefp,  "Checkpoint %d written after resolveSurrogates\n", ScaffoldGraph->checkPointIteration);
    CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_BEFORE_RESOLVE_SURROGATES+1);
  }

  //  This generates the 'rezlog/gapreads' file.  It's hugely
  //  expensive, usually dies on a negative variance assert, and as
  //  far as BPW knows, unused.
  //
  //Show_Reads_In_Gaps (GlobalData -> File_Name_Prefix);

  // now recompute mate pair statistics, once on scaffolds, once on
  // contigs, with the results on contigs being the ones that are
  // output in OutputMateDists

  if (immediateOutput == 0) {
    ComputeMatePairStatisticsRestricted(SCAFFOLD_OPERATIONS,
                                        minSamplesForOverride,
                                        "scaffold_final");
    ComputeMatePairStatisticsRestricted(CONTIG_OPERATIONS,
                                        minSamplesForOverride,
                                        "contig_final");
  }

  fprintf(GlobalData->stderrc,"* Output cam files *\n");
  GenerateCIGraph_U_Stats();
  GenerateLinkStats(ScaffoldGraph->CIGraph,"final",0);
  GeneratePlacedContigGraphStats("final",0);
  GenerateLinkStats(ScaffoldGraph->ContigGraph,"final",0);
  GenerateScaffoldGraphStats("final",0);
  GenerateSurrogateStats("final");
  //  GenerateContigAlignmentStats("final");
  
  #ifdef DEBUG
      int j = 0;
      for (j = 0; j < GetNumVA_CIFragT(ScaffoldGraph->CIFrags); j++) {
         CIFragT * frag = GetCIFragT(ScaffoldGraph->CIFrags, j);
         
         if (LookupValueInHashTable_AS(GlobalData->closureReads, frag->iid, 0)) {
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

  if(camFileOnly || generateOutput){
    CelamyAssembly(data->File_Name_Prefix);
  }

  /************* Output ***************/

  if(generateOutput){
    MarkContigEdges();
    OutputMateDists(ScaffoldGraph);

    ComputeMatePairDetailedStatus();
    OutputFrags(ScaffoldGraph);

    // We always have multiAlignments for Unitigs
    OutputUnitigsFromMultiAligns();
    OutputUnitigLinksFromMultiAligns();

    if(GlobalData->debugLevel > 0){
      DumpContigs(data->stderrc,ScaffoldGraph, FALSE);
      DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph, FALSE);
    }

    if (!ScaffoldGraph->doRezOnContigs)
      assert(0);

    OutputContigsFromMultiAligns();
    OutputContigLinks(ScaffoldGraph, outputOverlapOnlyContigEdges);

    OutputScaffolds(ScaffoldGraph);
    OutputScaffoldLinks(ScaffoldGraph);
  }

  DestroyScaffoldGraph(ScaffoldGraph);

  if(restartFromCheckpoint == NULLINDEX)
    DeleteGlobal_CGW(data);

  fprintf(stderr,"* Bye *\n");

  exit(0);
}

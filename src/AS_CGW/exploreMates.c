
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
static char CM_ID[] = "$Id: exploreMates.c,v 1.4 2005-03-22 19:48:37 jason_miller Exp $";


/*********************************************************************
 * Module:  based on AS_CGW_LoadCheckpoint
 * Description:
 *    For use with debugger to query values in a checkpoint
 * 
 *    Reference: 
 *
 *    Command Line Interface:
 *        $ loadcgw checkpointPath
 *
 *       CGBInputFiles: The file with new IUM,OUM, etc records to process. 
 *
 *       Checkpoint File: File named <outputName>.ckp.n
 *
 * 
 *********************************************************************/
//#define DEBUG 1
//#define DEBUG_BUCIS 1
//#define DEBUG_MERGE_SCAF 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "Stats_CGW.h"
#include "Instrument_CGW.h"

void DumpScaffoldContigUnitigPositions(CIScaffoldT * scaffold, FILE * fp)
{
  CIScaffoldTIterator CIs;
  ContigT * contig;
    
  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);
  while((contig = NextCIScaffoldTIterator(&CIs)) != NULL)
  {
    ChunkInstanceT * unitig;
    ContigTIterator unitigIterator;
    LengthT length;
    
    length = (contig->offsetAEnd.mean < contig->offsetBEnd.mean) ?
      contig->offsetAEnd : contig->offsetBEnd;
    
    fprintf(fp, "Contig " F_CID ": AEnd (%.3f,%.3f), BEnd (%.3f,%.3f)\n",
            contig->id,
            contig->offsetAEnd.mean, contig->offsetAEnd.variance,
            contig->offsetBEnd.mean, contig->offsetBEnd.variance);
    // Iterate over unitigs in contig & add data to contig instrumenter
    InitContigTIterator(ScaffoldGraph, contig->id,
                        TRUE, FALSE, &unitigIterator);
    while((unitig = NextContigTIterator(&unitigIterator)) != NULL)
    {
      fprintf(fp,
              "\tUnitig " F_CID ": AEnd (%.3f,%.3f), BEnd (%.3f,%.3f)\n",
              unitig->id,
              length.mean + unitig->offsetAEnd.mean,
              length.variance + unitig->offsetAEnd.variance,
              length.mean + unitig->offsetBEnd.mean,
              length.variance + unitig->offsetBEnd.variance);
    }
  }
}


void DumpContigUnitigPositions(FILE * fp)
{
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold = NULL;
  
  InitGraphNodeIterator(&scaffolds,
                        ScaffoldGraph->ScaffoldGraph,
                        GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL)
  {
    if(isDeadCIScaffoldT(scaffold) || scaffold->type != REAL_SCAFFOLD ||
       scaffold->info.Scaffold.numElements < 2)
      continue;

    fprintf(fp,
            "Scaffold " F_CID " - contig positions (mean,variance)\n",
            scaffold->id);
    DumpScaffoldContigUnitigPositions(scaffold, fp);
  }
}

void Usage(char * progName)
{
  fprintf (stderr,
           "USAGE:  %s -f <FragStoreName> -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum>\n",
           progName);
  exit (EXIT_FAILURE);
}


void DumpScaffoldGapSizes(CIScaffoldT * scaffold, FILE * fp)
{
  CIScaffoldTIterator CIs;
  ContigT * contig;
  float64 lastMax = -1;
  float64  thisMin = -1;
  CDS_CID_t lastID = NULLINDEX;
    
  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);
  while((contig = NextCIScaffoldTIterator(&CIs)) != NULL)
  {
    thisMin = min(contig->offsetAEnd.mean, contig->offsetBEnd.mean);
    if(lastID != -1)
      fprintf(fp, F_CID " " F_CID " " F_CID " %.f\n",
              scaffold->id, lastID, contig->id, thisMin - lastMax);
    lastMax = max(contig->offsetAEnd.mean, contig->offsetBEnd.mean);
    lastID = contig->id;
  }
}

void DumpAllGapSizes(FILE * fp)
{
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold = NULL;
  
  InitGraphNodeIterator(&scaffolds,
                        ScaffoldGraph->ScaffoldGraph,
                        GRAPH_NODE_DEFAULT);

  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL)
  {
    if(isDeadCIScaffoldT(scaffold) || scaffold->type != REAL_SCAFFOLD ||
       scaffold->info.Scaffold.numElements < 2)
      continue;
    DumpScaffoldGapSizes(scaffold, fp);
  }  
}

int main(int argc, char *argv[]){
  Global_CGW *data;
  char *outputPath = NULL;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int ckptNum = NULLINDEX;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->stderro = stderr;
  data->stderrfp = fopen("loadcgw.stderr","w");
  data->logfp = fopen("loadcgw.log","w");

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "f:g:n:c:h")) != EOF)){
      switch(ch) {
      case 'n':
	ckptNum = atoi(argv[optind - 1]);
	break;
      case 'c':
	{
	  strcpy( data->File_Name_Prefix, argv[optind - 1]);
	  setPrefixName = 1;

	}
	break;
      case 'f':
	{
	  strcpy( data->Frag_Store_Name, argv[optind - 1]);
	  setFragStore = 1;
	}
	break;
      case 'g':
	{
	  strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
	  setGatekeeperStore = 1;
	}
	break;
        case 'h':
          Usage(argv[0]);
          break;
      case '?':
	fprintf(stderr,"Unrecognized option -%c\n",optopt);
        Usage(argv[0]);
        break;
      default :
	errflg++;
      }
    }
    if((setPrefixName == FALSE) || (setFragStore == 0) || (setGatekeeperStore == 0))
      {
	fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setGatekeeperStore = %d outputPath = %s\n",
		argc, optind, setFragStore,setGatekeeperStore, outputPath);
        Usage(argv[0]);
      }
  }

    ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(data->File_Name_Prefix,ckptNum, FALSE);


    {
      /*
      FILE * fp;
      char filename[1024];
      
      ScaffoldGraphInstrumenter * sgi;

      sgi = CreateScaffoldGraphInstrumenter(ScaffoldGraph, INST_OPT_ALL);
      InstrumentScaffoldGraph(ScaffoldGraph, sgi, 0, CDS_COORD_MAX,
                              InstrumenterSilent, GlobalData->stderrc);
      */

      /*
      sprintf(filename, "contigUnitigPositions_%d.txt", ckptNum);
      fp = fopen(filename, "w");
      DumpContigUnitigPositions(fp);
      fclose(fp);
      */
      DumpAllGapSizes(stdout);
      /*
        This function iterates over all scaffolds in the graph & calls
        RecomputeOffsetsInScaffold() for each one
        Both functions are in LeastSquaresGaps_CGW.c
      */
      /*
      LeastSquaresGapEstimates(ScaffoldGraph, TRUE, FALSE, FALSE, FALSE, FALSE);
    
      sprintf(filename, "posAfterGapReestimation_%d.txt", ckptNum);
      fp = fopen(filename, "w");
      DumpContigUnitigPositions(fp);
      fclose(fp);
      
      sgi = CreateScaffoldGraphInstrumenter(ScaffoldGraph, INST_OPT_ALL);
      InstrumentScaffoldGraph(ScaffoldGraph, sgi, 0, CDS_COORD_MAX,
                              InstrumenterSilent, GlobalData->stderrc);
      */
    }
      

    /*
    fprintf(stderr,"* Hit Interrupt in Debugger to Proceed!!!! *\n");
    fflush(stderr);
    while(1){
      // wait for interrupt in debugger
    }
    */
    return 0;
}


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
static char CM_ID[] = "$Id: smallLargeScaffolds.c,v 1.6 2006-09-21 21:34:01 brianwalenz Exp $";

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
#include "AS_UTL_Hash.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "Instrument_CGW.h"
#include "InterleavedMerging.h"

#define SLOP 20000


extern void BuildMergedScaffoldEdges(ScaffoldGraphT *);

void Usage(char * progName)
{
  fprintf(stderr, "Usage: %s\n"
          "\t<-a  assemblyName>\n"
          "\t<-n  checkpointNumber>\n"
          "\t<-g  gkpStore>\n"
          "\t<-s  size>               - separating small from large scaffolds\n"
          "\t                           default = 30kb\n",
          progName);
  exit(1);
}

void InterleaveScaffolds(CIScaffoldT * scaffold,
                         CIScaffoldT * otherScaffold,
                         SEdgeT * sEdge,
                         ScaffoldAlignmentInterface * sai)
{
  SeqInterval interval;
  ChunkOrientationType orient =
    sEdge->idA == scaffold->id ? sEdge->orient :
    FlipEdgeOrient(sEdge->orient);
  
            
  PopulateScaffoldAlignmentInterface(scaffold, otherScaffold,
                                     sEdge, sai);
  fprintf(stdout, F_CID "(%d,%d) - " F_CID "(%d,%d) ",
          scaffold->id,
          (int) scaffold->bpLength.mean,
          (int) sqrt((double) scaffold->bpLength.variance),
          otherScaffold->id,
          (int) otherScaffold->bpLength.mean,
          (int) sqrt((double) otherScaffold->bpLength.variance));
  sai->segmentList = Align_Scaffold(sai->segmentList,
                                    sai->numSegs,
                                    sai->varWin,
                                    sai->scaffoldA->scaffold,
                                    sai->scaffoldB->scaffold,
                                    &(sai->best),
                                    sai->scaffoldA->bandBeg,
                                    sai->scaffoldA->bandEnd);
  interval.bgn =
    max(0, (-sEdge->distance.mean - 3 *
            sqrt((double) sEdge->distance.variance)) -
        (scaffold->bpLength.mean + 3 *
         sqrt((double) scaffold->bpLength.variance)));
  interval.end =
    min(otherScaffold->bpLength.mean + 3 *
        sqrt((double) otherScaffold->bpLength.variance),
        -sEdge->distance.mean + 3 *
        sqrt((double) sEdge->distance.variance));
  
  if(orient == AB_BA ||
     (scaffold->id == sEdge->idA && orient == BA_BA) ||
     (scaffold->id == sEdge->idB && orient == AB_AB))
    {
      CDS_COORD_t temp = interval.bgn;
      interval.bgn =
        max(0, otherScaffold->bpLength.mean -
            3 * sqrt((double) otherScaffold->bpLength.variance) -
            interval.end);
      interval.end = otherScaffold->bpLength.mean +
        3 * sqrt((double) otherScaffold->bpLength.variance) -
        temp;
    }
  
  if(sai->segmentList != NULL)
    {
      fprintf(stdout, "MERGEABLE ");
    }
  else if(sai->best == 0)
    {
      fprintf(stdout, "may be merged by interleaving ");
    }
  else
    {
      fprintf(stdout, "not mergeable ");
    }
  fprintf(stdout, "by edge: (%.f,%.f) %s weight " F_CID " - ",
          sEdge->distance.mean,
          sqrt((double) sEdge->distance.variance),
          orient == AB_AB ? "AB_AB" :
          (orient == AB_BA ? "AB_BA" :
           (orient == BA_AB ? "BA_AB" : "BA_BA")),
          sEdge->edgesContributing);
  interval.bgn = (interval.bgn + interval.end - scaffold->bpLength.mean) / 2;
  interval.end = interval.bgn + scaffold->bpLength.mean;
  fprintf(stdout, "interval in larger scaffold: (" F_COORD "," F_COORD ")\n",
          interval.bgn, interval.end);
}


void ChangeToGappedCoordinatesAndLengths(ScaffoldGraphT * graph)
{
  // Change contig lengths to be gapped lengths
  fprintf(stderr, "Changing contig lengths to be gapped\n");
  {
    GraphNodeIterator nodes;
    ChunkInstanceT * ci;
    CIScaffoldT * scaffold;
    CDS_COORD_t scaffoldLength;
    InitGraphNodeIterator(&nodes, graph->RezGraph, GRAPH_NODE_UNIQUE_ONLY);
    while((ci = NextGraphNodeIterator(&nodes)) != NULL)
      {
        MultiAlignT * ma;
        ma = LoadMultiAlignTFromSequenceDB(graph->sequenceDB, ci->id, FALSE);
        ci->bpLength.mean = GetMultiAlignLength(ma);
        UnloadMultiAlignTFromSequenceDB(graph->sequenceDB, ci->id, FALSE);
      }
  }

  // change contig coordinates & scaffold lengths to be gapped
  fprintf(stderr, "Changing contig positions in scaffolds & scaffold lengths to be gapped\n");
  {
    GraphNodeIterator scaffolds;
    CIScaffoldT *scaffold;
    
    InitGraphNodeIterator(&scaffolds,
                          ScaffoldGraph->ScaffoldGraph,
                          GRAPH_NODE_DEFAULT);
    while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL)
      {
        CDS_COORD_t lastEndOld = 0;
        CDS_COORD_t lastEndNew = 0;
        CIScaffoldTIterator cis;
        ChunkInstanceT * ci;

        /*
          -------------                     ------------------------
          0       c1EndOld             c2StartOld              c2EndOld
          |                     |
          c2StartOld - c1EndOld

          0       c1LenNew
        */
        InitCIScaffoldTIterator(ScaffoldGraph,
                                scaffold, TRUE, FALSE, &cis);
        while((ci = NextCIScaffoldTIterator(&cis)) != NULL)
          {
            if(ci->offsetAEnd.mean < ci->offsetBEnd.mean)
              {
                ci->offsetAEnd.mean = lastEndNew +
                  (ci->offsetAEnd.mean - lastEndOld);
                lastEndOld = ci->offsetBEnd.mean;
                ci->offsetBEnd.mean = ci->offsetAEnd.mean + ci->bpLength.mean;
                lastEndNew = ci->offsetBEnd.mean;
              }
            else
              {
                ci->offsetBEnd.mean = lastEndNew +
                  (ci->offsetBEnd.mean - lastEndOld);
                lastEndOld = ci->offsetAEnd.mean;
                ci->offsetAEnd.mean = ci->offsetBEnd.mean + ci->bpLength.mean;
                lastEndNew = ci->offsetAEnd.mean;
              }
          }
        scaffold->bpLength.mean = lastEndNew;
      }
  }
}

    
int main(int argc, char *argv[])
{
  int32  checkPointNumber = -1;
  ScaffoldGraphT * graph = NULL;
  int sizeLimit = 30000;

  GlobalData = CreateGlobal_CGW();
  GlobalData->stderrc = stderr;
  GlobalData->File_Name_Prefix[0] = '\0';
  GlobalData->Gatekeeper_Store_Name[0] = '\0';

  // parse the command line
  {
    int ch, errflg = 0;
    char * temp_char;
    
    optarg = NULL;
    /*
      a: assembly name
      n: checkpoint number
      g: gatekeeper store name
      s: size distinguishing small from large scaffolds
    */
    while( !errflg &&
           ((ch = getopt( argc,
                          argv,
                          "a:n:g:s:" )) != EOF) )
      {
        switch(ch)
          {
            case 'a':
              strcpy(GlobalData->File_Name_Prefix, optarg);
              break;
            case 'n':
              checkPointNumber = atoi(optarg);
              break;
            case 'g':
              strcpy(GlobalData->Gatekeeper_Store_Name, optarg);
              break;
            case 's':
              sizeLimit = atoi(optarg);
              break;
            default:
              Usage(argv[0]);
              break;
          }
      }
  }

  if(checkPointNumber < 0 ||
     strlen(GlobalData->File_Name_Prefix) == 0 ||
     strlen(GlobalData->Gatekeeper_Store_Name) == 0)
    Usage(argv[0]);

  // load the checkpoint
  fprintf( stderr, "Loading scaffold graph %s, checkpoint %d.\n",
           GlobalData->File_Name_Prefix, checkPointNumber);
  graph = LoadScaffoldGraphFromCheckpoint(GlobalData->File_Name_Prefix,
                                          checkPointNumber,
                                          FALSE );
  if(graph == NULL)
    {
      fprintf(stderr, "Failed to load graph from checkpoint!\n");
      return 1;
    }

  // open the gatekeeper store
  InitGateKeeperStore(&(graph->gkpStore), GlobalData->Gatekeeper_Store_Name);
  OpenReadOnlyGateKeeperStore(&(graph->gkpStore));

  ChangeToGappedCoordinatesAndLengths(graph);
  
  BuildMergedScaffoldEdges(graph);

  // now try interleaving as appropriate
  {
    GraphNodeIterator scaffolds;
    CIScaffoldT *scaffold;
    int numScaffolds = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
    ScaffoldAlignmentInterface * sai = CreateScaffoldAlignmentInterface();

    InitGraphNodeIterator(&scaffolds,
                          ScaffoldGraph->ScaffoldGraph,
                          GRAPH_NODE_DEFAULT);
    while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL)
      {
        SEdgeTIterator SEdges;
        SEdgeT * sEdge;
        SEdgeT * lastEdge = NULL;
        int numEdges = 0;
        SeqInterval distanceRange;
        SeqInterval sigmaRange;
        int weight = 0;

        if(isDeadCIScaffoldT(scaffold) ||
           scaffold->type != REAL_SCAFFOLD ||
           scaffold->bpLength.mean > sizeLimit)
          continue;

        InitSEdgeTIterator(ScaffoldGraph, scaffold->id,
                           FALSE, FALSE, ALL_END, FALSE, &SEdges);
        while((sEdge = NextSEdgeTIterator(&SEdges)) != NULL)
          {
            CIScaffoldT * otherScaffold =
              GetCIScaffoldT(ScaffoldGraph->CIScaffolds,
                             sEdge->idA == scaffold->id ?
                             sEdge->idB : sEdge->idA);

            if(!isDeadCIScaffoldT(otherScaffold) &&
               otherScaffold->type == REAL_SCAFFOLD &&
               otherScaffold->bpLength.mean > sizeLimit &&
               sEdge->edgesContributing >= 2)
              {
                if(sEdge->distance.mean < 0)
                  {
                    if(numEdges == 0)
                      {
                        distanceRange.bgn = sEdge->distance.mean;
                        distanceRange.end = sEdge->distance.mean;
                        sigmaRange.bgn = sEdge->distance.mean -
                          3 * sqrt((double) sEdge->distance.variance);
                        sigmaRange.end = sEdge->distance.mean +
                          3 * sqrt((double) sEdge->distance.variance);
                        weight = sEdge->edgesContributing;
                      }
                    numEdges++;
                    if(lastEdge == NULL)
                      lastEdge = sEdge;
                    else if(sEdge->idA == lastEdge->idA &&
                            sEdge->idB == lastEdge->idB &&
                            sEdge->orient == lastEdge->orient &&
                            sEdge->distance.mean + SLOP >
                            distanceRange.bgn &&
                            sEdge->distance.mean - SLOP <
                            distanceRange.end)
                      {
                        distanceRange.bgn = min(distanceRange.bgn, sEdge->distance.mean);
                        distanceRange.end = max(distanceRange.end, sEdge->distance.mean);
                        sigmaRange.bgn = min(sigmaRange.bgn,
                                             sEdge->distance.mean -
                                             3 * sqrt((double) sEdge->distance.variance));
                        sigmaRange.end = max(sigmaRange.end,
                                             sEdge->distance.mean +
                                             3 * sqrt((double) sEdge->distance.variance));
                        weight += sEdge->edgesContributing;
                        continue;
                      }
                    else
                      {
                        lastEdge = NULL;
                        break;
                      }
                  }
              }
          }
        if(lastEdge != NULL)
          {
            CIScaffoldT * otherScaffold =
              GetCIScaffoldT(ScaffoldGraph->CIScaffolds,
                             lastEdge->idA == scaffold->id ?
                             lastEdge->idB : lastEdge->idA);
            if(numEdges > 1)
              {
                lastEdge->distance.mean =
                  (distanceRange.bgn + distanceRange.end) / 2;
                lastEdge->distance.variance =
                  max(sigmaRange.end - lastEdge->distance.mean,
                      lastEdge->distance.mean - sigmaRange.bgn);
                lastEdge->distance.variance =
                  lastEdge->distance.variance * lastEdge->distance.variance / 9;
                lastEdge->edgesContributing = weight;
              }
            InterleaveScaffolds(scaffold, otherScaffold, lastEdge, sai);
          }
        else if(numEdges == 0)
          {
            fprintf(stdout,
                    F_CID " - no negative distance edges to large scaffolds.\n",
                    scaffold->id);
          }
        else
          {
            fprintf(stdout,
                    F_CID " - %d conflicting negative distance edges to large scaffolds.\n",
                    scaffold->id, numEdges);
          }
      }
  }

  return 0;
}

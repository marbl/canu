
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
static char CM_ID[] = "$Id: GreedyOverlapCIREZ.c,v 1.2 2004-09-23 20:25:28 mcschatz Exp $";

#include <math.h>
#include <string.h>
#include <assert.h>

#include "AS_global.h"
#include "GreedyOverlapREZ.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_ALN_aligners.h"
#include "UtilsREZ.h"


#define REZ_DP_ERATE  .1
#define REZ_DP_THRESH 1e-6
#define REZ_DP_MINLEN 30



FILE *logFile=NULL;


int isRepeatOverlap(EdgeCGW_T *cie, GraphCGW_T *graph)
{
  uint64 start1,end1,start2,end2;
  uint min1,max1,min2,max2;
  
  start1 = GetGraphNode(graph, cie->idA)->aEndCoord;
  end1   = GetGraphNode(graph, cie->idA)->bEndCoord;
  start2 = GetGraphNode(graph, cie->idB)->aEndCoord;
  end2   = GetGraphNode(graph, cie->idB)->bEndCoord;

  if(logFile == NULL)
    logFile = stderr;

#if GREEDYDEBUG > 1
  fprintf(logFile,"+++ SIMULATOR : IdA %ld=%ld,%ld  : IdB %ld=%ld,%ld\n",cie->idA,start1,end1,cie->idB,start2,end2);
#endif

  if( start1 < end1 )
    {
      min1 = start1;
      max1 = end1;
    }
  else
    {
      min1 = end1;
      max1 = start1;
    }

  if( start2 < end2 )
    {
      min2 = start2;
      max2 = end2;
    }
  else
    {
      min2 = end2;
      max2 = start2;
    }

  if( min1 >= min2 && min1 <= max2)
    return FALSE;
  if( min2 >= min1 && min2 <= max1)
    return FALSE;
  
  return TRUE;
}



extern char *Filename_Prefix;


void checkEdgeQuality(GraphCGW_T* graph,int minLength, int qt, char* fp)
{
  int cid;
  GraphEdgeIterator iterator;
  EdgeCGW_T         *edge;
  ChunkOverlapCheckT olap;
  float quality, best = 1.0;
  EdgeCGW_T    *bestEdge;
  int fpcount = 0;
  int count   = 0;
  int repeatNodes = 0; // number of nodes that have an edge to a repeat
  int hasRepeat;

  if(logFile == NULL)
    logFile = stderr;
#if GREEDYDEBUG > 0
  char fileName[100];
  Filename_Prefix = fp;
  sprintf(fileName,"%s.greedy.%d.%d.log",Filename_Prefix,minLength,qt);
  logFile = file_open(fileName,"w");
  fprintf(logFile,"+++ Entering check edge quality for %d nodes +++\n",GetNumGraphNodes(graph));
#endif

  if(MatchTableREZ == NULL)
    MatchTableREZ    = fill_QV_Match_table(prob_match);
  if(MismatchTableREZ == NULL)
    MismatchTableREZ = fill_QV_Mismatch_table(prob_mismatch);
  
  for(cid=0; cid<GetNumGraphNodes(graph); cid++)
    {
      int nEdges = 0;
      int nRep   = 0;
      int nNorm  = 0;
      hasRepeat = FALSE; 
      // none of the outgoing overlap edges goes to a repeat

      if( GetGraphNode(graph, cid) == NULL ) 
	continue;
      
      InitGraphEdgeIterator(graph,cid,ALL_END,ALL_EDGES,GRAPH_EDGE_DEFAULT,&iterator);

      // we iterate over all edges and keep record of the best for each node
      best = 1.0;
      bestEdge = NULL;
      while( (edge = NextGraphEdgeIterator(&iterator)) != NULL )
	{
	  nEdges++;
	  if( TRUE == LookupQualityOverlap(graph,edge,
					   edge->orient,&olap,BAYESIAN,
					   &quality,logFile) )
	    {
	      // we regard an edge as better if its quality is better and
	      // its Minimum length is at least MinLength
	      if(quality < best && olap.overlap > minLength )
		{
		  best = quality;
		  bestEdge = edge;
		}
	      if(isRepeatOverlap(edge,graph))
		{
		  nRep++; 
		  // We count all repetitiv edges 
#if GREEDYDEBUG > 0 
 		  fprintf(logFile,"+++ REPEAT overlap from %ld with quality %f and length %ld\n",cid,quality,olap.overlap);
#endif
		  if( hasRepeat == FALSE)
		    {
		      repeatNodes++;
		      hasRepeat = TRUE;	// This node has a repetitive overlap
		    }	  
		}
	      else
		{
		  nNorm++;
#if GREEDYDEBUG > 0 
		  fprintf(logFile,"+++ NORMAL overlap from %ld with quality %f and length %ld\n",cid,quality,olap.overlap);
#endif
		}
	    }
	}
#if GREEDYDEBUG > 0   
      fprintf(logFile,"There were %d edges overall. %d repetitive and %d normal\n",nEdges,nRep,nNorm);
#endif
      if( bestEdge != NULL )
	if( hasRepeat == TRUE )
	  if( isRepeatOverlap(bestEdge,graph) == TRUE )
	    {
#if GREEDYDEBUG > 0    
	      fprintf(logFile,"=== Best edge has quality %f and is repeat\n",best);
#endif
	      fpcount++;
	    }
	  else
	    {
#if GREEDYDEBUG > 0 
	      fprintf(logFile,"+++ Best edge has quality %f and is normal\n",best);
#endif
	      count++;
	    }
    }

#if GREEDYDEBUG > 0
  fprintf(logFile,"+++ STAT : Total number of nodes with repeat edges =%d\nNo of best repeats = %d, No of best normals %d\n",fpcount+count,fpcount,count);
#endif

    if(MatchTableREZ != NULL)
      {
	free(MatchTableREZ);
	MatchTableREZ = NULL;
      }
    if(MismatchTableREZ != NULL)
      {
	free(MismatchTableREZ);
	MismatchTableREZ = NULL;
      }
#if GREEDYDEBUG > 0
    fclose(logFile);
#endif
}


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
/**********************************************************************
 Module:       ConfirmStones.c
 Description:  File contains function that takes a Scaffold_Fill_t
               and tries to confirm single linked "stones" with overlap paths
  
    Programmer:  K. Remington
       Written:  18 August 99
   Assumptions:  The chunk and scaffold ids are consistent with the global
                 variable Global_CGW in file ../AS_CGW/Globals_CGW.h

**********************************************************************/

static char CM_ID[] = "$Id: ConfirmStones.c,v 1.1.1.1 2004-04-14 13:52:57 catmandew Exp $";


/* ---------------------------------------------------- */
/* standard includes */
/* ---------------------------------------------------- */
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


/* ---------------------------------------------------- */
/* REZ includes */
/* ---------------------------------------------------- */
#include "ConsistencyChecksREZ.h"
#include "DataTypesREZ.h"
#include "CommonREZ.h"
#include "UtilsREZ.h"
#include "GapWalkerREZ.h"
#include "BccREZ.h"
#include "SubgraphREZ.h"
#include "GWDriversREZ.h"

/* ---------------------------------------------------- */
/* AS includes */
/* ---------------------------------------------------- */
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UnionFind_AS.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ChiSquareTest_CGW.h"


/* ---------------------------------------------------- */
/* local #defines and typedefs */
/* ---------------------------------------------------- */


#define DEBUG 1
// Define Level of debugging. 0 Means no debug information

/* ---------------------------------------------------- */
/* variable declarations */
/* ---------------------------------------------------- */

extern char * Filename_Prefix;

#if DEBUG > 0
static  FILE* logFile;
#endif

int ConfirmStones(Scaffold_Fill_t *tentative, int noScaff, int iteration){
  /******************************************************************/
  // Precondition: The chunk and scaffold ids are consistent with the
  //global variable Global_CGW in file ../AS_CGW/Globals_CGW.h
  /*****************************************************************/
  GraphEdgeIterator iterator;
  CIEdgeT         *edge;
  chunk_subgraph *g;

  int i,j,k;
  int numGaps;
  int numChunks;

  char filename[100];

  LengthT covered;
  covered.mean     = 0.0;
  covered.variance = 0.0;

  sprintf(filename,"%s.rezsc.i%d.log",Filename_Prefix,iteration);
  logFile = file_open(filename,"w");


#if DEBUG > 1
  fprintf(logFile,"Stone confirmation was called for %d Scaffolds.\n",noScaff);
#endif

  for(i=0; i<noScaff; i++){ 
    /* This is the main loop that iterates over all ChunkOfScaffolds */
    numGaps = tentative[i].num_gaps;
    
#if DEBUG > 1
    fprintf(logFile,"********* Scaffold %d *********\n",i);
    fprintf(logFile,"Number of gaps: %d.\n\n",numGaps);
#endif

    for(j=0; j<numGaps; j++){ 
      /* This loop iterates over all gaps in a particular Scaffold 
	 Gap 0 is before the first chunk and the last gap is after
	 last chunk */
      numChunks = tentative[i].gap[j].num_chunks;

#if DEBUG > 1
      fprintf(logFile,"\nGAP %d\n",j);
      fprintf(logFile,"Number of unitigs assigned to the gap: %d\n",numChunks);
#endif
      g = Build_Subgraph_Gap(tentative, i, j, Is_Not_Bogus_Not_Guide);

      for(k=0; k<numChunks; k++){  
	/* for each gap we iterate over the chunks assigned and 
	   determine whether we keep them or not 
           Here, the chunks are assumed to be the 'stones' thrown in
           the gap, that is, the chunks which have a single mate link
           placing them in the gap.
        */
	int cidk = tentative[i].gap[j].chunk[k].chunk_id;
	ChunkInstanceT  *chunkk;
	chunkk = GetChunkInstanceT(ScaffoldGraph->ChunkInstances, cidk);
	
#if DEBUG > 1
	fprintf (logFile, "\nArt's ID %6ld\n",cidk);
	fprintf (logFile, "\n unitig %ld\n",
		 chunkk->id
		 //		 chunk1->offsetAEnd.mean,
		 //		 chunk1->offsetAEnd.stdDev,
		 //		 chunk1->offsetBEnd.mean,
		 //		 chunk1->offsetBEnd.stdDev
		 );
#endif

	// First we do not want to keep the chunk
	tentative[i].gap[j].chunk[k].keep = FALSE;


	// Initialize the iterator for the CIEdges
	InitGraphEdgeIterator(ScaffoldGraph->RezGraph,cidk,ALL_END,ALL_EDGES,GRAPH_EDGE_DEFAULT,&iterator);


	while( (edge = NextGraphEdgeIterator(&iterator)) != NULL )
	  {
	    /* In this loop we iterate over all merged Chunk Instance edges
	       of the assigend chunks */
	    ChunkInstanceT *chunk2;
	    int citer;
	    
	    /* which is the other unitig ? */
	    if( edge->idA == cidk )
	      citer = edge->idB;
	    else
	      citer = edge->idA;

	    /* Now we check whether the CI edge points
	       to a unique chunk. We only test CI edges
	       between assigned chunks */

	    chunk2 = GetGraphNode(ScaffoldGraph->RezGraph, citer);  
	    assert( citer == chunk2->id );

#if DEBUG > 1
	    fprintf(logFile,"CI edge goes from %d to %d\n",
		    edge->idA,edge->idB);
	    fprintf(logFile,"Length of gap (%7lf, %5.3lf)\n",
		    edge->distance.mean,sqrt(edge->distance.variance));		      
#endif
	  
         }
       }
     }
  }
  
  fclose(logFile);
  return(0);
}

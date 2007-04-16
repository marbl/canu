
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
 Module:       ConsistencyChecksRez.c
 Description:  File contains functions that take a Scaffold_Fill_t
               and checks the contig->gap assignment given by it.
  
    Programmer:  K. Reinert
       Written:  17 May 99
   Assumptions:  The chunk and scaffold ids are consistent with the global
                 variable Global_CGW in file ../AS_CGW/Globals_CGW.h

**********************************************************************/

static char CM_ID[] = "$Id: ConsistencyChecksREZ.c,v 1.8 2007-04-16 17:34:15 brianwalenz Exp $";


/* ---------------------------------------------------- */
/* standard includes */
/* ---------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>


/* ---------------------------------------------------- */
/* REZ includes */
/* ---------------------------------------------------- */
#include "ConsistencyChecksREZ.h"
#include "DataTypesREZ.h"
#include "CommonREZ.h"
#include "UtilsREZ.h"
#include "GapWalkerREZ.h"
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
#include "ChiSquareTest_CGW.h"


extern int  Global_Debug_Flag;

/* ---------------------------------------------------- */
/* static function and variable declaration */
/* ---------------------------------------------------- */

static void combine_two_distrib(LengthT dist1, LengthT dist2, LengthT *cDist);

static OverlapStatusREZ check_overlap(Gap_Chunk_t cidA, Gap_Chunk_t cidB, 
                                      int addtoCG, ChunkOverlapCheckT *olap);

static void  Canonicalize
    (CIEdgeT * edge, const Gap_Chunk_t * c1, const Gap_Chunk_t * c2);

static bool check_distribs(LengthT *, 
			   LengthT *, 
			   LengthT *);

static bool left_of(const Gap_Chunk_t *,
		    const Gap_Chunk_t *);

static void estimate_gap_distrib(const Gap_Chunk_t *,
				 const Gap_Chunk_t *,
				 LengthT *);



#if  CHECK_CELSIM_COORDS
static bool test_coordinates(const Gap_Chunk_t *);
#endif

/* ---------------------------------------------------- */
/* local #defines and typedefs */
/* ---------------------------------------------------- */


#define DEBUG 0
// Define Level of debugging. 0 Means no debug information

/* ---------------------------------------------------- */
/* variable declarations */
/* ---------------------------------------------------- */

#if DEBUG > 1
extern char * Filename_Prefix;
static  FILE* logFile;
#endif

CIOrient Get_Gap_Chunk_Orient(const Gap_Chunk_t * chunk) {
  if (chunk->end.mean > chunk->start.mean) {
    return A_B;
  }else{
    return B_A;
  }
}

int Is_Edge_Orientation_Consistent(CIEdgeT * edge,
				   const Gap_Chunk_t * left,
				   const Gap_Chunk_t * right) {
  //
  // checks if the orientation of the edge is consistent with the orientation
  // of the chunks
  //
  CIOrient leftCIorient = X_X, rightCIorient = X_X;
  CIEdgeT  tmp_edge = * edge;
  
  // Contained edges can have two different orientations
  // Make  tmp_edge  be the one that matches how
  // this test is done.
  if  (tmp_edge . distance . mean < 0.0)
      Canonicalize (& tmp_edge, left, right);

  switch (GetEdgeOrientationWRT(& tmp_edge, left->chunk_id)) {
  case AB_AB:
    //      leftCI                                        rightCI
    //  A --------------------- B               A --------------------- B
    //    5'----->                                           <------5'
    leftCIorient = A_B;
    rightCIorient = A_B;
    break;
  case AB_BA:
    //      leftCI                                        rightCI
    //  A --------------------- B               B --------------------- A
    //    5'----->                                           <------5'
    leftCIorient = A_B;
    rightCIorient = B_A;
    break;
  case BA_BA:
    //      leftCI                                        rightCI
    //  B --------------------- A               B --------------------- A
    //    5'----->                                           <------5'
    leftCIorient = B_A;
    rightCIorient = B_A;
    break;
  case BA_AB:
    //      thisCI                                        otherCI
    //  B --------------------- A               A --------------------- B
    //    5'----->                                           <------5'
    leftCIorient = B_A;
    rightCIorient = A_B;
    break;
  default:
    assert(0);
    break;
  }
  
  //
  // now we take care of which is to left of which
  //
  if(left_of(left,right)) 
    return ((Get_Gap_Chunk_Orient(left) == leftCIorient) &&
	    (Get_Gap_Chunk_Orient(right) == rightCIorient));
  else 
    return ((Get_Gap_Chunk_Orient(left) != leftCIorient) &&
	    (Get_Gap_Chunk_Orient(right) != rightCIorient));
  
  //return ((left->flipped == right->flipped) && ((orient == AB_AB) || (orient = BA_BA)) ||
  //  (left->flipped != right->flipped) && ((orient == AB_BA) || (orient = BA_AB)));
}


int Is_Edge_Consistent(CIEdgeT *edge,
		       const Gap_Chunk_t *left,
		       const Gap_Chunk_t *right) {
  //
  // consistency check for an edge and the two flanking chunks positions
  //
  LengthT
    estimated_distribution,
    combined_distribution;
  CIEdgeT  tmp_edge;
  int retVal;
  assert(edge != NULL);
  assert(left != NULL);
  assert(right != NULL);

  // Contained edges can have two different orientations
  // Make  tmp_edge  be the one that matches how
  // this test is done.
  tmp_edge = * edge;
  if  (tmp_edge . distance . mean < 0.0)
      Canonicalize (& tmp_edge, left, right);
      
  estimate_gap_distrib(left, right, &estimated_distribution);
  combine_two_distrib(estimated_distribution, tmp_edge . distance, &combined_distribution);
  retVal = check_distribs(&estimated_distribution, &(tmp_edge . distance),
       &combined_distribution);

  /*  if(!retVal){
    fprintf(stderr,"* Edge (%d,%d,%c) is inconsistent\n",
	    edge->idA, edge->idB, edge->orient);
  }
  */

  return retVal;
}

/* ---------------------------------------------------- */
/* non static functions */
/* -----------------------------------------------------*/


/* declaration of two function in ChunkOverlap_CGW.c */
void ComputeCanonicalOverlap(Global_CGW *data, 
			     ChunkOverlapperT *chunkOverlapper, 
			     ChunkOverlapCheckT *canOlap);

int InitCanonicalOverlapSpec(int cidA, int cidB, 
			     ChunkOrientationType orientation, ChunkOverlapSpecT *spec);



int check_consistency(Scaffold_Fill_t *gapAssignment, int noScaff, int iteration){
  /******************************************************************/
  // Precondition: The chunk and scaffold ids are consistent with the
  // global variable Global_CGW in file ../AS_CGW/Globals_CGW.h
  //               
  /*****************************************************************/
  GraphEdgeIterator iterator;
  CIEdgeT         *edge;

  int i,j,k,l;
  int numGaps;
  int numChunks;

  int round1Rejected = 0;
  int round1Accepted = 0;
  int round1All = 0;
  int round2Rejected = 0;
  int round2Accepted = 0;
  int round2All = 0;
#if DEBUG > 1
  char filename[100];
#endif

  LengthT covered;
  covered.mean     = 0.0;
  covered.variance = 0.0;

#if DEBUG > 1
      sprintf(filename,"%s.rezcc.i%d.log",Filename_Prefix,iteration);
      logFile = file_open(filename,"w");
#endif  

#if DEBUG > 1
  fprintf(logFile,"check_consistency was called for %d Scaffolds.\n",noScaff);
#endif

  for(i=0; i<noScaff; i++){ 
    /* This is the main loop that iterates over all ChunkOfScaffolds */
    numGaps = gapAssignment[i].num_gaps;
    
#if DEBUG > 1
    fprintf(logFile,"********* Scaffold %d *********\n",i);
    fprintf(logFile,"Number of gaps: %d.\n\n",numGaps);
#endif

    for(j=0; j<numGaps; j++){ 
      /* This loop iterates over all gaps in a particular Scaffold 
	 Gap 0 is before the first chunk and the last gap is after
	 last chunk */
      numChunks = gapAssignment[i].gap[j].num_chunks;
     
#if DEBUG > 1
      fprintf(logFile,"\nGAP %d\n",j);
      fprintf(logFile,"Number of unitigs assigned to the gap: %d\n",numChunks);
#endif

      for(k=0; k<numChunks; k++){  
	/* for each gap we iterate over the chunks assigned and 
	   determine whether we keep them or not */
	int cid1 = gapAssignment[i].gap[j].chunk[k].chunk_id;
	bool allSucceeded = TRUE; 
	// if the distribution test fails for one edge mate
	// this variable is set to FALSE
	bool tested1 = FALSE;
	// If the chunk is tested in round 1 this is set to TRUE

	ChunkInstanceT  *chunk1;
	chunk1 = GetGraphNode(ScaffoldGraph->RezGraph, cid1);
	round1All++;
	
#if DEBUG > 1
	fprintf (logFile, "\nArt's ID %6ld\n",cid1);
	fprintf (logFile, "\n unitig %ld\n",
		 chunk1->id
		 //		 chunk1->offsetAEnd.mean,
		 //		 chunk1->offsetAEnd.stdDev,
		 //		 chunk1->offsetBEnd.mean,
		 //		 chunk1->offsetBEnd.stdDev
		 );
#endif

	// First we do not want to keep the chunk
	gapAssignment[i].gap[j].chunk[k].keep = FALSE;


	// Initialize the iterator for the CIEdges
	InitGraphEdgeIterator(ScaffoldGraph->RezGraph,cid1,ALL_END,ALL_EDGES,GRAPH_EDGE_DEFAULT,&iterator);

	while( (edge = NextGraphEdgeIterator(&iterator)) != NULL )
	  {
	    /* In this loop we iterate over all merged Chunk Instance edges
	       of the assigend chunks */
	    ChunkInstanceT *chunk2;
	    int citer;
	    
	    /* which is the other unitig ? */
	    if( edge->idA == cid1 )
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
	  
	    
	    if( ! chunk2->flags.bits.isUnique ){ 
	      /* Now we check whether the link goes to another
		 assigned chunk. This is done by scanning over all
		 other chunks and hence takes overall quadratic time
		 in the number of assigned chunks. This should be
		 changed to set or hash operations */
	      
	      for(l=0; l<numChunks; l++){// l loop
		int cid2 = gapAssignment[i].gap[j].chunk[l].chunk_id;
		/* do not check against the same chunk */
		if( l == k )
		  continue;
		
		if( citer == cid2 ){ 
		  LengthT estDist,combDist;
		
#if DEBUG > 5
		  // overlapper test
		  {//start
		    int isCanonical;
		    ChunkOrientationType orientation;
		    ChunkOrientType orientA, orientB;
		    ChunkOverlapCheckT *lookup;
		    Gap_Chunk_t cidA = gapAssignment[i].gap[j].chunk[k];
		    Gap_Chunk_t cidB = gapAssignment[i].gap[j].chunk[l];

		    fprintf(logFile,"checking overlap precomputed overlap\n");

		    if( left_of(&cidA,&cidB) ){
		      if( cidA.start.mean < cidA.end.mean )
			orientA = A_B;
		      else
			orientA = B_A;
		      if( cidB.start.mean < cidB.end.mean )
			orientB = A_B;
		      else
			orientB = B_A;
		    }
		    else{
		      if( cidA.start.mean < cidA.end.mean )
			orientA = B_A;
		      else
			orientA = A_B;
		      if( cidB.start.mean < cidB.end.mean )
			orientB = B_A;
		      else
			orientB = A_B;
		    }
		  
		    orientation = GetChunkPairOrientation(orientA,orientB);
		    fprintf(logFile,"chunkA %d chunkB %d computed orient %c\n",
			    edge->idA,edge->idB,orientation);
		
		    fprintf(logFile,"chunkA %d chunkB %d orient %c\n",edge->idA,edge->idB, edge->orient);

		    if( isOverlapCIEdgeT(edge)){
		    
		      ChunkOverlapCheckT olap;
		      LengthT dist;
		      estimate_gap_distrib(&cidA,&cidB,&dist);
		      isCanonical = InitCanonicalOverlapSpec(cidA.chunk_id,
							     cidB.chunk_id,
							     orientation,
							     &olap.spec);
		  
		      lookup = LookupCanonicalOverlap(ScaffoldGraph->overlapper,&olap.spec);
		      if( lookup != NULL){
			int oldlook;
			ChunkOverlapCheckT olap;
			oldlook=lookup->overlap;

			olap = OverlapChunks(ScaffoldGraph->overlapper,      // doesn't handle suspicious -- debug code
					       cidA.chunk_id,cidB.chunk_id,
					       orientation, 
					       lookup->minOverlap,
					       lookup->maxOverlap,
					       0.1);
			/*	lookup = OverlapChunks(ScaffoldGraph->overlapper, 
					       cidA.chunk_id,cidB.chunk_id,
					       orientation, 
					       lookup->minOverlap,
					       lookup->maxOverlap,
					       0.1);*/
			
			fprintf(logFile,"computed overlap %d\n",olap.overlap);
			//	CheckCIEdgesAgainstChunkOverlapper(ScaffoldGraph);
		      }
		    }
		  }
#endif


		  /* we found an assigned unitig */
		
#if DEBUG > 2
		  fprintf(logFile,"*** %d is an assigned unitig ***\n",cid2);
#endif

		  /* Now we estimate the distance between the two unitigs
		     based on their computed coordinates */
		  estimate_gap_distrib(&gapAssignment[i].gap[j].chunk[k],
				       &gapAssignment[i].gap[j].chunk[l],
				       &estDist);

#if DEBUG > 2
		  fprintf(logFile,"estimated mean %5.3lf and var %5.3lf\n",
			  estDist.mean,estDist.variance);
		  fprintf(logFile,"edge      mean %5.3lf and var %5.3lf\n",
			  edge->distance.mean,edge->distance.variance);
#endif
		
		  /* Now we combine the estimated distance with the distance
		     in the edge */
		  combine_two_distrib(estDist,edge->distance,
				      &combDist);
#if DEBUG > 2
		  fprintf(logFile,"combined  mean %5.3lf and var %5.3lf\n",
			  combDist.mean,combDist.variance);
#endif
		  tested1 = TRUE;
		  /* finally we check whether the estimated and the combined
		     and the edge and the combined (mean,3*stdDev) intervals
		     intersect */
		  if( check_distribs(&estDist,&edge->distance,&combDist) == TRUE ){
#if DEBUG > 2
		    fprintf(logFile,"+++ check_distrib succeeded for chunk %d +++\n",cid1);
		  //	  if( Intersection(&estDist,&edge->distance) <= 0 ){
		  //		    test2Failed++;
		  //		    fprintf(logFile,"+++ check_distrib2 DID NOT succeeded +++\n");
		  //		  }
		  
#endif
		  
		    //		    gapAssignment[i].gap[j].chunk[k].keep = TRUE;

		    // we force the overlap test !!!!
		    gapAssignment[i].gap[j].chunk[k].keep = FALSE;

		  }
		  else{
		    allSucceeded = FALSE;
#if DEBUG > 1
		    fprintf(logFile,"--- check_distrib failed ---\n");
#endif
		    /* The first test failed. Before rejecting anything we
		       check the quality of the three CI edges C,m,O
		       where m is the current CI edge between the two unitigs and 
		       C and O are the sets of CI edges joining the current
		       (k loop) resp. other (l loop) unitig to the scaffold.
		       If the average quality of C is higher than that of m, we
		       discard m and accept the positioning of the current chunk
		       for the moment. Otherwise we mark it and proceed to the
		       next current chunk.

		       THIS IS NOT DONE. SHOULD WE DO IT ? 
		    */
		  
		    gapAssignment[i].gap[j].chunk[k].keep = FALSE;
		    break;
		  }	      
		}
	      }
	    }	  
	    if( ! allSucceeded )
	      break;
	  }
	
	if( allSucceeded && tested1 ){
#if DEBUG > 2 && CHECK_CELSIM_COORDS
	  if( test_coordinates(&gapAssignment[i].gap[j].chunk[k]) == FALSE )
	    fprintf(logFile,"!!! Should have not kept chunk %d !!!\n",
		    gapAssignment[i].gap[j].chunk[k].chunk_id); 
#endif
	  round1Accepted++;
	}
	else{
	  round1Rejected++;
	}
	
	if( gapAssignment[i].gap[j].chunk[k].keep == FALSE ){
	  int overlaps   = 0;
	  int noOverlaps = 0;
	  int wrongOverlapChunk = 0;
	  int goodChunks[numChunks]; 
	  int goodChunksTop=0;
	  // in this stack we keep track
	  // of the chunks that were tested for overlap with chunk k and passed

	  /*   If the above tests fail or if a unitig has not been tested
	       because he has no CI edges to other assigned unitigs, 
	       we finally check whether
	       the assigned positions indicate an overlap. If so,
	       we check specifically for this overlap. If it is present
	       we set the appropriate overlap flag.
	       We accept the positioning of a chunk if two conditions
	       are fullfilled
	       1) all or all but one
	       of the indicated overlaps is true. 
	       This takes care of
	       the following situation: ------     1
	                                  -------   2
 	                                     ------  3
	       If chunk 1 is positioned wrongly it will miss two overlaps
	       whereas chunk 2 and 3 have a confirmed overlap with each other.
	       In the case that one overlap is false we check, 
	       whether it is consistently
	       false with all chunks that form good overlaps
	  */ 
	  round2All++;
	  
#if DEBUG > 1
	  fprintf(logFile,"ooo Testing overlaps of rejected chunks ooo\n");
#endif
	  for(l=0; l<numChunks; l++){
	    // l loop
	    int first,second;
	    OverlapStatusREZ oStatus;
	    ChunkOverlapCheckT olap;
	    
	    if( l == k )
	      continue;
	    
	    if( left_of(&gapAssignment[i].gap[j].chunk[k],
			&gapAssignment[i].gap[j].chunk[l]) ){
	      first  = k;
	      second = l;
	    }
	    else{
	      first  = l;
	      second = k;
	    }
	    
	    oStatus = check_overlap(gapAssignment[i].gap[j].chunk[first],
				    gapAssignment[i].gap[j].chunk[second],
				    TRUE,&olap);
	    
	    if( oStatus == ASS_OVLP_TRUE ){
	      overlaps++;
	      goodChunks[goodChunksTop++] = l;
	    }
	    if( oStatus == ASS_OVLP_FALSE ){
	      noOverlaps++;
	      wrongOverlapChunk = l;
	    }	    
	  }
	  
	  
	  if( noOverlaps == 0 ){
	    gapAssignment[i].gap[j].chunk[k].keep = TRUE;
#if DEBUG > 1
	    fprintf(logFile,"+++ check_overlap succeeded  for unitig %d +++\n",cid1);
#endif


	    round2Accepted++;
#if DEBUG > 1 && CHECK_CELSIM_COORDS
	    if( test_coordinates(&gapAssignment[i].gap[j].chunk[k]) == FALSE )
	      fprintf(logFile,"!!! Coordinates are bad for unitig %d !!!\n",
		      gapAssignment[i].gap[j].chunk[k].chunk_id); 
#endif
	  }
	  else
	    if( noOverlaps == 1 ){
	      bool keep=TRUE;
	      /* If there is no good chunk, we do not keep the chunk */
	      if( goodChunksTop == 0 )
		keep = FALSE;
	      
	      while( --goodChunksTop > 0 ){
		int cidA = goodChunks[goodChunksTop];
		int cidB = wrongOverlapChunk;
		int first,second;
		OverlapStatusREZ oStatus;
		ChunkOverlapCheckT olap;
		assert( cidA != cidB );
		
		if( left_of(&gapAssignment[i].gap[j].chunk[cidA],
			    &gapAssignment[i].gap[j].chunk[cidB]) ){
		  first  = cidA;
		  second = cidB;
		}
		else{
		  first  = cidB;
		  second = cidA;
		}
		
		oStatus = check_overlap(gapAssignment[i].gap[j].chunk[first],
					gapAssignment[i].gap[j].chunk[second],
					TRUE,&olap);
		
#if DEBUG > 1
		fprintf(logFile,"+++ check_overlap for unitigs %d and %d +++\n",
			gapAssignment[i].gap[j].chunk[first].chunk_id,
			gapAssignment[i].gap[j].chunk[second].chunk_id);
#endif
		if( oStatus == ASS_OVLP_TRUE ){
#if DEBUG > 1
		  fprintf(logFile,"Oops. It succeeded\n");
#endif
		  keep = FALSE;
		  break;
		}
		// All assumed overlaps with the good overlap should be false
	      }

	      
	      if( keep == TRUE ){
#if DEBUG > 2
		fprintf(logFile,"+++ check_overlap (with one failed) succeeded  for unitig %d +++\n",cid1);
#endif

		gapAssignment[i].gap[j].chunk[k].keep = TRUE;
		round2Accepted++;
#if DEBUG > 2 && CHECK_CELSIM_COORDS
		if( test_coordinates(&gapAssignment[i].gap[j].chunk[k]) == FALSE )
		  fprintf(logFile,"!!! Coordinates are bad for unitig %d !!!\n",
			  gapAssignment[i].gap[j].chunk[k].chunk_id); 
#endif
	      }
	      else{
#if DEBUG > 2
		fprintf(logFile,"--- check_overlap (with one failed) failed for unitig %d ---\n",cid1);
#endif
		round2Rejected++;
	      }
	    }
	    else{
#if DEBUG > 2
	      fprintf(logFile,"--- check_overlap failed for unitig %d ---\n",cid1);
#endif
	      round2Rejected++;
	    } 
	}
      }
    }
  }
#if DEBUG > 2
  fprintf(logFile,
	  "================\n"
	  "Check Statistics\n"
	  "================\n"
	  "No. of round 1 chunks          : %d\n"
	  "No. of round 1 chunks accepted : %d \n"
	  //	  "No. of round 1 chunks with test2 failed : %d \n"
	  "No. of round 1 chunks rejected : %d \n"
	  "No. of round 2 chunks          : %d\n"
	  "No. of round 2 chunks accepted : %d \n"
	  "No. of round 2 chunks rejected : %d \n",
	  round1All,round1Accepted,round1Rejected,
	  round2All,round2Accepted,round2Rejected);
#endif
#if DEBUG > 1
  fclose(logFile);
#endif
  //  return(round1Accepted+round2Accepted);
  return(round2Accepted);
  // since we forced the overlap check the final number of accepted
  // chunks is now round2Accepted
}


/* ----------------------------------------------*/

static void combine_two_distrib(LengthT dist1, LengthT dist2, LengthT *cDist){
  // note that we take the stdDev field instead of directly the variance
  // field, as some LengthT structs have it not set.

  //  float v1 = dist1.stdDev*dist1.stdDev;
  //  float v2 = dist2.stdDev*dist2.stdDev;  // Bug was fixed by Saul

  float64 v1 = dist1.variance;
  float64 v2 = dist2.variance;

  if( v1 == 0.0)
    {
      cDist->variance = v2;
      cDist->mean     =dist2.mean;
    }
  else
    if( v2 == 0.0 )
      {
	cDist->variance = v1;
	cDist->mean     =dist1.mean;
      }
    else
      {
	cDist->variance = 1.0 / ((1.0 / v1 ) + (1.0 / v2));
	cDist->mean = (int) (((dist1.mean / v1) + (dist2.mean / v2)) * cDist->variance);
      }
  return;
}


/* ----------------------------------------------*/

static OverlapStatusREZ check_overlap(Gap_Chunk_t cidA, Gap_Chunk_t cidB, 
			       int addtoCG, ChunkOverlapCheckT *olap){
  /* We pass as arguments the two Gap_Chunk_t structs that contain
     the assumed positions. 

     If addToCG is TRUE, we add the overlap !!! NOT DONE !!
     information to the chunk graph. 

     The functions returns NO_ASS_OVLP
     if the assumed positions do not indicate an possible overlap.
     It returns ASS_OVLP_TRUE if the positions indicate an overlap
     and their is one and ASS_OVLP_FALSE otherwise.
     In addition a more detailed description of the overlap is returned 
     in olap. 

     Precondition : The function assumes that the minimal coordinate
                    of cidA is less than the minimal coordinate
		    of cidB.
  */

  int minOlap;
  int maxOlap;
  // these value indicate the range for an overlap within these bounds

  LengthT dist;
  // this struct contains the presumed length of the gap between
  // the two chunks.

  ChunkOrientationType orientation;
  // Holds the orientation (AB_AB, AB_BA...) we compute according to the
  // assumed positions

  ChunkOrientType orientA, orientB;
  // The orientation of the two chunks
  
  float relError=0.0;
  // the relative error of the assumed overlap and the computed

  estimate_gap_distrib(&cidA,&cidB,&dist);


  if( dist.mean+AS_REZ_MIN_OVERLAP > 0 ){
#if DEBUG > 2
    fprintf(logFile,"do not assume overlap %lf,%lf\n",dist.mean,sqrt(dist.variance));
#endif
    return NO_ASS_OVLP;
  }
  else{
#if DEBUG > 2
    fprintf(logFile,"assuming overlap %lf,%lf\n",dist.mean,sqrt(dist.variance));
#endif
    if( cidA.start.mean < cidA.end.mean )
      orientA = A_B;
    else
      orientA = B_A;
    if( cidB.start.mean < cidB.end.mean )
      orientB = A_B;
    else
      orientB = B_A;
    
    orientation = GetChunkPairOrientation(orientA,orientB);
    // Compute the orientation of the two chunks
    //    printf("Orient = %c cidA %d, cidB %d \n",orientation,cidA.chunk_id,cidB.chunk_id);
   
    // Old test
    //   minOlap = -dist.mean-3*sqrt(dist.variance);
    //   maxOlap = -dist.mean+3*sqrt(dist.variance);
    // we changed it so that we assume the biggest possible overlap

   minOlap = AS_REZ_MIN_OVERLAP;
   if( abs(cidA.end.mean-cidA.start.mean) > abs(cidB.end.mean-cidB.start.mean))
     maxOlap = abs(cidB.end.mean-cidB.start.mean);
   else
     maxOlap = abs(cidA.end.mean-cidA.start.mean);

   *olap = OverlapChunks(ScaffoldGraph->RezGraph,    // handles suspicious
			 cidA.chunk_id, cidB.chunk_id,
			 orientation, 
			 minOlap,maxOlap,
			 AS_REZ_ERROR_RATE, TRUE);
 
   if(olap->suspicious){
	      fprintf(stderr,"* SUSPICIOUS Overlap found! Looked for (%d,%d,%c) found (%d,%d,%c)\n",
		      cidA.chunk_id, cidB.chunk_id, orientation,
		      olap->spec.cidA, olap->spec.cidB, olap->spec.orientation);
   }else if(olap->overlap > 0){
     relError = fabs((float)(olap->overlap+dist.mean)/(float)olap->overlap);
   }
   else
     relError = 2*AS_REZ_MAX_REL_ERROR;

    if( relError < AS_REZ_MAX_REL_ERROR )
      {
	// here we might add an edge mate to the chunk graph
#if DEBUG > 2
	fprintf(logFile,"overlap of length %d\n",olap->overlap);
#endif
	return ASS_OVLP_TRUE;
      }
    else{
#if DEBUG > 2
      fprintf(logFile,"overlap between %d and %d failed with relative error of %f\n",cidA.chunk_id,cidB.chunk_id,relError);
      fprintf(logFile,"overlap of length %d\n",olap->overlap);
#endif
      return ASS_OVLP_FALSE; 
    }
  }
}




/* ----------------------------------------------*/
/* STATIC functions */
/* ----------------------------------------------*/



static void estimate_gap_distrib(const Gap_Chunk_t *cT1,
				 const Gap_Chunk_t *cT2,
				 LengthT *gDist){
  /* the function takes two Gap_chunk_t structs
     and computes an estimate of the mean and variance of the gap between
     the two chunks. In order to do this one has to make a case distinction
     depending of the orientation of both chunks.
     
     The function returns the mean and variance in the struct gDist

     Precondition : Each of the CIEdges must have an orientation
  */

  const  Gap_Chunk_t *chunkT1, *chunkT2;
  ChunkOrientType orient1, orient2;


  if( left_of(cT1,cT2) ){
    chunkT1 = cT1;
    chunkT2 = cT2;
  }
  else{
    chunkT1 = cT2;
    chunkT2 = cT1;
  }  	    


  if( chunkT1->start.mean < chunkT1->end.mean )
    orient1 = A_B;
  else
    orient1 = B_A;
  if( chunkT2->start.mean < chunkT2->end.mean )
    orient2 = A_B;
  else
    orient2 = B_A;


  
#if 0
  fprintf(stderr,"Unitig1 %d Unitig2 %d (%lf,%lf,%lf,%lf) (%lf,%lf,%lf,%lf)\n",
	  cT1->chunk_id,
	  cT2->chunk_id,
	  cT1->start.mean,
	  cT1->start.variance,
	  cT1->end.mean,
	  cT1->end.variance,
	  cT2->start.mean,
	  cT2->start.variance,
	  cT2->end.mean,
	  cT2->end.variance
	  );  
#endif

  switch(orient1){
  case  A_B :
#if 0
    fprintf(logFile,"orient1=A_B\n");
#endif
    switch(orient2){
    case A_B :
#if 0
      fprintf(logFile,"orient2=A_B\n");
#endif
      gDist->mean     = chunkT2->start.mean - chunkT1->end.mean;
      gDist->variance = chunkT2->start.variance + chunkT1->end.variance;
      break;
    case B_A :
#if 0
      fprintf(logFile,"orient2=B_A\n");
#endif
      gDist->mean     = chunkT2->end.mean - chunkT1->end.mean;
      gDist->variance = chunkT2->end.variance + chunkT1->end.variance;
      break;
    default:
      assert(0);
      break;
    }
    break;
  case B_A :
#if 0
    fprintf(logFile,"orient1=B_A\n");
#endif
    switch(orient2){
    case A_B :  
#if 0
      fprintf(logFile,"orient2=A_B\n");
#endif
      gDist->mean     = chunkT2->start.mean - chunkT1->start.mean;
      gDist->variance = chunkT2->start.variance + chunkT1->start.variance;
      break;
    case B_A :	    
#if 0
      fprintf(logFile,"orient2=B_A\n");
#endif
      gDist->mean     = chunkT2->end.mean - chunkT1->start.mean;
      gDist->variance = chunkT2->end.variance + chunkT1->start.variance;
      break;
    default:
      assert(0);
      break;
    }
    break;
  default:
    assert(0);
    break;
  }
  return;
}



static void  Canonicalize
    (CIEdgeT * edge, const Gap_Chunk_t * c1, const Gap_Chunk_t * c2)

//  Flip  edge  if necessary so that it represents the gap between
//   c1  and  c2 , depending on which begins to the left of the other.

  {
   bool  c1_flipped, c2_flipped, need_change;
   double  left_c1, left_c2;

   // Make edge go from c1 to c2 if necessary
   if  (edge -> idA != c1 -> chunk_id)
       {
        const Gap_Chunk_t  * tmp;

        tmp = c1;
        c1 = c2;
        c2 = tmp;
       }
assert (edge -> idA == c1 -> chunk_id && edge -> idB == c2 -> chunk_id);

   if  (c1 -> start . mean < c1 -> end . mean)
       {
        c1_flipped = FALSE;
        left_c1 = c1 -> start . mean;
       }
     else
       {
        c1_flipped = TRUE;
        left_c1 = c1 -> end . mean;
       }
   if  (c2 -> start . mean < c2 -> end . mean)
       {
        c2_flipped = FALSE;
        left_c2 = c2 -> start . mean;
       }
     else
       {
        c2_flipped = TRUE;
        left_c2 = c2 -> end . mean;
       }

   need_change = FALSE;
   if  (left_c1 < left_c2)
       {
        switch (edge -> orient)
          {
           case  AB_AB :
           case  AB_BA :
             if  (c1_flipped)
                 need_change = TRUE;
             break;
           case  BA_AB :
           case  BA_BA :
             if  (! c1_flipped)
                 need_change = TRUE;
             break;
           default :
             fprintf (stderr, "YIKES:  Bad orientation = %d at line %d file %s\n",
                  (int) edge -> orient, __LINE__, __FILE__);
             assert (FALSE);
          }
       }
     else
       {
        switch (edge -> orient)
          {
           case  AB_AB :
           case  AB_BA :
             if  (! c1_flipped)
                 need_change = TRUE;
             break;
           case  BA_AB :
           case  BA_BA :
             if  (c1_flipped)
                 need_change = TRUE;
             break;
           default :
             fprintf (stderr, "YIKES:  Bad orientation = %d at line %d file %s\n",
                  (int) edge -> orient, __LINE__, __FILE__);
             assert (FALSE);
          }
       }
       
   if  (need_change)
       {
        double  c1_len, c2_len;
        double  neg_half;

        c1_len = fabs (c1 -> end . mean - c1 -> start . mean);
        c2_len = fabs (c2 -> end . mean - c2 -> start . mean);
        neg_half = (c1_len + c2_len) / -2.0;
        edge -> distance . mean = - c1_len - c2_len - edge -> distance . mean;

        switch (edge -> orient)
          {
           case  AB_AB :
             edge -> orient = BA_BA;
             break;
           case  AB_BA :
             edge -> orient = BA_AB;
             break;
           case  BA_AB :
             edge -> orient = AB_BA;
             break;
           case  BA_BA :
             edge -> orient = AB_AB;
             break;
           default :
             fprintf (stderr, "YIKES:  Bad orientation = %d at line %d file %s\n",
                  (int) edge -> orient, __LINE__, __FILE__);
             assert (FALSE);
          }
       }

   return;
  }



static bool check_distribs(LengthT *est, 
			   LengthT *given, 
			   LengthT *comb){
  
  /* this function returns TRUE if the 3*sdtDev intervals around the
     means of est and given intersect both with comb. This function
     could be refined to yield TRUE only if a certain percentage of
     overlap is achieved. */

  if( (Intersection(est,comb) > 0) && (Intersection(given,comb) > 0) )
    return TRUE;
  else
    return FALSE;

}




static bool left_of(const Gap_Chunk_t *cT1,
		    const Gap_Chunk_t *cT2){
  /* the function takes two Gap_chunk_t structs
     and returns true if the minimum coordinate of the first is
     less than the minimum coordinate of the second.
  */

  if(MIN(cT1->start.mean,cT1->end.mean) < MIN(cT2->start.mean,cT2->end.mean) )
    return TRUE;
  else
    return FALSE;

}  	    




#if  CHECK_CELSIM_COORDS
static bool test_coordinates(const Gap_Chunk_t *ct1){
#if DEBUG > 2
  fprintf(logFile,"simulated (%d,%d), computed (%lf,%lf)\n",ct1->sim_start,ct1->sim_end,ct1->start.mean,ct1->end.mean);
#endif

  if( ct1->sim_start < ct1->start.mean-AS_REZ_SIMTEST)
    return FALSE;
  if( ct1->sim_start > ct1->start.mean+AS_REZ_SIMTEST)
    return FALSE;
  if( ct1->sim_end < ct1->end.mean-AS_REZ_SIMTEST)
    return FALSE;
  if( ct1->sim_end > ct1->end.mean+AS_REZ_SIMTEST)
    return FALSE;

  return TRUE;
}
#endif

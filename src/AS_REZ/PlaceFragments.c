
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

        Module:  PlaceFragments.c

   Description:  Place fragments into the consensus after everybody else is done

    Programmer:  M. Flanigan

       Written:  28 October 99

 **********************************************************************/

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_ALN_aligners.h"
#include "AS_CNS_multialign.h"
//
// AS_CGW
//
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"

//
// AS_REZ
//
#include "DataTypesREZ.h"
#include "UtilsREZ.h"
#include "CommonREZ.h"
#include "GapWalkerREZ.h"
#include "BccREZ.h"
#include "SubgraphREZ.h"
#include "ConsistencyChecksREZ.h"
#include "GWDriversREZ.h"

//
// AS_CNS
//
#include "PlaceInAlignment.h"

//
// globals
//
char
* GW_Filename_Prefix;

int Print_Scaffold_Cam(CDS_CID_t);
int Print_Mates_Cam(CIFragT *frag, CDS_CID_t fragSurrogate,
                    CIFragT *mateFrag, CDS_CID_t mateFragChunkID);
int Print_Fragment_Cam(CIFragT *frag, CDS_CID_t fragSurrogate);
int GetChunkPositionInScaffold(ChunkInstanceT *chunk,
                               CDS_COORD_t *left_end,
                               CDS_COORD_t *right_end,
                               int *chunkScaffoldOrientation);
// next one should go in GraphCGW_t.h
int SurrogateChunkDistance(GraphCGW_T *graph,  
                           CIFragT *frag, 
                           CIFragT *mfrag, 
                           DistT *dist, 
                           int type, 
                           int orientIn, 
                           GraphEdgeStatT *stat,
                           LengthT *distance);
int PlaceUnscaffoldedContigs(float (* quality)(CIEdgeT *, CDS_CID_t),
                             float qualityThresh);
int CountFragments(float (* quality)(CIEdgeT *, CDS_CID_t),
                   float qualityThresh);
int CountFragments2(float (* quality)(CIEdgeT *, CDS_CID_t),
                    float qualityThresh);

int attemptedBasedOnMate = 0,
  placedBasedOnMate = 0,
  attemptedBasedOnOverlaps = 0,
  attemptedBasedOnOverlapsOnly = 0,
  placedBasedOnOverlaps = 0,
  placedBasedOnOverlapsOnly = 0,
  attemptedMatePairs = 0,
  placedMatePairs = 0,
  totalFragmentsPlaced,
  totalFragmentsAttempted;

int PlaceFragments()
{

  CountFragments(Bayesian_Quality, 0.5);
  CountFragments2(Bayesian_Quality, 0.5);
  return(0);

  fprintf(stderr,"------------------------------: fragment placing running\n");
  Print_Scaffold_Cam(0);

  PlaceFragmentsInSingletonSurrogates();

  PlaceFragmentsWithPlacedMates();
  // we need to force re-contigging before proceeding

  PlaceFragmentsWithUnplacedMates();

  PlaceRemainingFragments();


  // performance calculations
  totalFragmentsAttempted = attemptedBasedOnMate + 2 * attemptedMatePairs + attemptedBasedOnOverlapsOnly;
  totalFragmentsPlaced = placedBasedOnMate + 2 * placedMatePairs + placedBasedOnOverlaps;
  
  fprintf(stderr, "totalFragementsAttempted = %d\n", totalFragmentsAttempted);
  if (totalFragmentsAttempted > 0)
	fprintf(stderr,
                "numFragmentsPlaced: %d (%.2f %%)\n", totalFragmentsPlaced, 
                100.0 * (float) totalFragmentsPlaced / totalFragmentsAttempted);
  else
	fprintf(stderr, "numFragmentsPlaced: %d (%.2f %%)\n",
                totalFragmentsPlaced, 0.0);
  fprintf(stderr, "placed / attempted BasedOnMate: %d / %d\n",
          placedBasedOnMate, attemptedBasedOnMate);
  fprintf(stderr, "placed / attempted MatePairs: %d / %d\n",
          placedMatePairs, attemptedMatePairs);
  fprintf(stderr, "placed / attempted BasedOnOverlaps (purely): %d / %d\n",
          placedBasedOnOverlapsOnly, 
          attemptedBasedOnOverlapsOnly);
  fprintf(stderr, "placed / attempted BasedOnOverlaps: %d / %d\n",
          placedBasedOnOverlaps, attemptedBasedOnOverlaps);


  Print_Scaffold_Cam(1);
  fprintf(stderr, "------------------------------: fragment placing done\n");

  PlaceUnscaffoldedContigs(Bayesian_Quality, 0.5);

}


typedef struct 
{
	  CDS_CID_t fragID;
	  int start;                   // where this fragment statrts in the fragmentInfo array
	  int count;                   // the number of possible positions of the frag (== # entries in fragmentInfo array)
	  double score;                // presumably some combo of coverage and bayesianQuality
	  CDS_CID_t parentChunkID;           // the chunk that this frag belongs to initially
	  CDS_CID_t candidateChunkID;        // so we know where to place a fragment when its time comes
	  char placed;
} fragmentSummary;

// fragment contig contig_start contig_end coverage quality 
typedef struct 
{
	  CDS_CID_t fragID;
	  SeqInterval contigPosition;  // so we can see when fragments interact
	  int coverage;
	  double bayesianQuality;
	  CDS_CID_t chunkID;                 // so we know where to place a fragment when its time comes
	  CDS_CID_t contigID;
	  int fragContigOrientation;   // needed in DP_Compare_AS
} fragmentInformation;

// When a frag has changed coverage, we need to take the candidateChunkID, find the parent and recompute
// the score.  We really should figure out how to do it just based on the candidateChunkID, but perhaps
// this is okay if the score is going to involve the various possible places for this chunk.

static int compareFragmentScores(const void *fragSumIn1,
                                 const void *fragSumIn2)
{
  //
  // The function that compares two fragmentSummary structs
  //
  fragmentSummary *fragSum1 = (fragmentSummary *) fragSumIn1;
  fragmentSummary *fragSum2 = (fragmentSummary *) fragSumIn2;
  double fragScore1, fragScore2;
  
  assert((fragSum1 != NULL) && (fragSum2 != NULL));

  if (fragSum1->placed) 
	fragScore1 = 10.0;
  else
	fragScore1 = fragSum1->score;

  if (fragSum2->placed) 
	fragScore2 = 10.0;
  else
	fragScore2 = fragSum2->score;  
  
  if (fragScore1 > fragScore2)
    return 1;
  else if (fragScore1 < fragScore2)
    return -1;
  else
    return 0;
}


int PlaceRemainingFragments() {
  int
    i,
    canPlace;
  chunk_subgraph
    * f;  
  Gap_Fill_t
    * fc;
  ChunkInstanceT
    * parentChunk;
  LengthT
    gap;
  GraphNodeIterator
	CIGraphIterator;
  int
	begin = 0,
	end = 100,
	scaffoldID = 0,
	numFragmentPositions = 0,
	numFragmentInfos = 0,
	numFragmentsToPlace = 0,
	fragmentCount = 0,
	fragmentSumIndex, 
	fragmentInfoIndex,
	placed;
  fragmentSummary   // fragmentSummary holds where the fragments members are in fragmentInfo array and their best score
	*fragmentSum;   // its purpose is to allow us to pick the best fragment
  fragmentInformation 
	*fragmentInfo;		// fragmentInfo holds data on the possible placements of each fragment
  
  // count fragments we're interested in * number of surrogates that fragment is in
  // scan all the scaffolds
  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while (parentChunk = NextGraphNodeIterator(&CIGraphIterator))
  {
    // scan all the nodes (chunks)
	// if numInstances >= 2 then it has at least two surrogates - all parents with a single surrogate have already been handled
	// using 1 here just for testing purposes, and below
	if (parentChunk->info.CI.numInstances >= 1) // || parentChunk->scaffoldID == -1
	{
	  // get the MultiAlign for this parentChunk
	  MultiAlignT *maParent = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, parentChunk->id);
	  assert (maParent);
	 
	  // count fragments and positions
	  numFragmentsToPlace += GetNumIntMultiPoss(maParent->f_list);
	  numFragmentPositions += GetNumIntMultiPoss(maParent->f_list) * parentChunk->info.CI.numInstances;
	}
  }
  
  fprintf( stderr, "numFragmentPositions = %d\n", numFragmentPositions);
  if (numFragmentsToPlace > 0)
  {
	fragmentSum = (fragmentSummary *) malloc (sizeof(fragmentSummary) * numFragmentsToPlace);
	assert (fragmentSum);
	fragmentInfo = (fragmentInformation *) malloc (sizeof(fragmentInformation) * numFragmentPositions);
	assert (fragmentInfo);

	// scan all the scaffolds
	InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
	while (parentChunk = NextGraphNodeIterator(&CIGraphIterator))
	{
	  // scan all the nodes (chunks)
	  // if numInstances >= 2 then it has at least two surrogates - all parents with a single surrogate have already been handled
	  // using 1 here just for testing purposes, as above
	  if (parentChunk->info.CI.numInstances >= 1) // || parentChunk->scaffoldID == -1
	  {
		MultiAlignT *maParent;
		
		fprintf(stderr, "\nparentChunk " F_CID " has %d instances\n",
                        parentChunk->id, parentChunk->info.CI.numInstances);
		fprintf(stderr, "parentChunk " F_CID " has %d fragments\n",
                        parentChunk->id, parentChunk->info.CI.numFragments);
		fprintf(stderr, "parentChunk " F_CID " scaffoldID: " F_CID "\n",
                        parentChunk->id, parentChunk->scaffoldID);
		
		// Now cycle through the fragments in a parentChunk.  Check its relative fitness to go with 
		// each instance based on matepair consistency.
		
		// get the MultiAlign for this parentChunk
		maParent = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, parentChunk->id);
		assert (maParent);
		
		// cycle through fragments 
		for (i = 0; i < GetNumIntMultiPoss(maParent->f_list); i++)
		{
		  IntMultiPos *mpParent = GetIntMultiPos(maParent->f_list, i);
		  CIFragT 
			*frag = GetCIFragT(ScaffoldGraph->CIFrags, (CDS_CID_t)mpParent->source);
		  CDS_CID_t fragID = (CDS_CID_t)mpParent->source;
		  
		  totalFragmentsAttempted++;
		  
		  fprintf( stderr, "parentChunk " F_CID " fragID: " F_CID "\n", parentChunk->id, fragID);
		  fprintf( stderr, "fragID: " F_CID " frag->iid: " F_CID "\n", fragID, frag->iid);
		  fprintf( stderr, "       frag " F_CID " mateOf: " F_CID "\n", fragID, frag->mateOf);
		  
		  fragmentSum[fragmentCount].fragID = frag->iid;
		  fragmentSum[fragmentCount].placed = FALSE;
		  fragmentSum[fragmentCount].start = numFragmentInfos;
		  fragmentSum[fragmentCount].count = parentChunk->info.CI.numInstances;
		  fragmentSum[fragmentCount].parentChunkID = parentChunk->id;
		  fragmentSum[fragmentCount].score = 2.0;

		  canPlace = ScoreFragmentBasedOnOverlaps(frag, parentChunk, fragmentInfo, &numFragmentInfos);
		  fragmentCount++;
		  
		  attemptedBasedOnOverlaps++;
		  attemptedBasedOnOverlapsOnly++;
		  if (canPlace == 1) 
		  {
			placedBasedOnOverlaps++;
			placedBasedOnOverlapsOnly++;
			fprintf(stderr, "pointpoint: frag->iid was just placed 6: " F_CID "\n", frag->iid);
		  }
		  else
			fprintf(stderr, "pointpoint: frag->iid was just placed 6: " F_CID " not! pure overlap-based placement failure \n", 
					frag->iid);
		}
	  }
	}
  }

  placed = 1;
  while (placed == 1)
  {
	placed = 0;
	PickBestFragment(fragmentSum, numFragmentsToPlace, fragmentInfo, numFragmentInfos, 
					 &fragmentSumIndex, &fragmentInfoIndex);
	if (fragmentSumIndex != -1)
	{
	  placed = PlaceFragment(fragmentSumIndex, fragmentSum);

	  fprintf(stderr, "at qsort, numFragmentsToPlace = %d\n",
                  numFragmentsToPlace);
	  
	  if (placed == 1)
	  {
		numFragmentsToPlace--;
		qsort(fragmentSum, numFragmentsToPlace, sizeof(fragmentSummary), compareFragmentScores);
		rescoreOverlappingFragments(fragmentSum, numFragmentsToPlace, fragmentInfo, numFragmentInfos, fragmentInfoIndex);
	  }
	}
  }  

  for (i = 0; i < numFragmentsToPlace; i++)
  {
	fprintf (stderr, "          fragmentSum[%d].fragID = " F_CID "\n",
                 i, fragmentSum[i].fragID);
	fprintf (stderr, "          fragmentSum[%d].placed = %d\n",
                 i, fragmentSum[i].placed);
	fprintf (stderr, "          fragmentSum[%d].start  = %d\n",
                 i, fragmentSum[i].start);
	fprintf (stderr, "          fragmentSum[%d].count  = %d\n", i
                 , fragmentSum[i].count);
	fprintf (stderr, "          fragmentSum[%d].score  = %f\n",
                 i, fragmentSum[i].score);
	fprintf (stderr, "   fragmentSum[%d].parentChunkID = " F_CID "\n",
                 i, fragmentSum[i].parentChunkID);
	fprintf (stderr, "fragmentSum[%d].candidateChunkID = " F_CID "\n\n",
                 i, fragmentSum[i].candidateChunkID);
  }
  free(fragmentSum);

  for (i = 0; i < numFragmentInfos; i++)
  {
	fprintf (stderr, "            fragmentInfo[%d].fragID = " F_CID "\n",
                 i, fragmentInfo[i].fragID);
	fprintf (stderr, "          fragmentInfo[%d].coverage = %d\n",
                 i, fragmentInfo[i].coverage);
	fprintf (stderr, "   fragmentInfo[%d].bayesianQuality = %f\n",
                 i, fragmentInfo[i].bayesianQuality);
	fprintf (stderr, "  fragmentInfo[%d].chunkID = " F_CID "\n",
                 i, fragmentInfo[i].chunkID);
	fprintf (stderr, "fragmentInfo[%d].contigPosition.bgn = " F_COORD "\n",
                 i, fragmentInfo[i].contigPosition.bgn);
	fprintf (stderr, "fragmentInfo[%d].contigPosition.end = " F_COORD "\n\n",
                 i, fragmentInfo[i].contigPosition.end);
  }
  free(fragmentInfo);
  
  return 0;
}

int rescoreOverlappingFragments(fragmentSummary *fragmentSum,
                                int numFragmentsToPlace,
                                fragmentInformation *fragmentInfo,
                                int numFragmentInfos,
                                int fragmentInfoIndex)
{
  int 
	i,
	j,
	leftPosition, 
	rightPosition,
	contigOfOverlappingFragment;
  ChunkInstanceT
	*chunkOfPlacedFragment,
	*chunkOfOverlappingFragment;
  
  leftPosition = fragmentInfo[fragmentInfoIndex].contigPosition.bgn;
  rightPosition = fragmentInfo[fragmentInfoIndex].contigPosition.end;

  // now step through all fragments
  for (j = 0; j < numFragmentsToPlace; j++)
  {
	if (fragmentSum[j].placed == FALSE)
	{
	  for (i = fragmentSum[j].start; i < fragmentSum[j].start + fragmentSum[j].count; i++)
	  {
		if (i != fragmentInfoIndex)  // don't check the fragment we just placed
		{
		  if ((fragmentInfo[fragmentInfoIndex].contigID != fragmentInfo[i].contigID) &&
			  ((fragmentInfo[i].contigPosition.bgn > leftPosition && fragmentInfo[i].contigPosition.bgn < rightPosition) ||
			   (fragmentInfo[i].contigPosition.end > leftPosition && fragmentInfo[i].contigPosition.end < rightPosition)))
			   
		  {
			char 
			  fragSequence[CNS_MAX_SEQUENCE_LENGTH], 
			  fragQuality[CNS_MAX_SEQUENCE_LENGTH];
			MultiAlignT 
			  *maCandidate;
			ChunkInstanceT 
			  *candidateChunk;
		  
			// fprintf(stderr, "in rescoreOverlappingFragments, placed frag from ");
			
			candidateChunk = GetGraphNode(ScaffoldGraph->CIGraph, fragmentInfo[i].chunkID);
			assert (candidateChunk);
			
			maCandidate = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, candidateChunk->info.CI.contigID);
			assert (maCandidate);
			
			// grab the frags clear sequence
			getFragmentSequenceAndQuality( fragmentInfo[i].fragID, fragSequence, fragQuality);
			
			// rescore surrogate
			ScoreFragmentBasedOnOverlapsInSurrogate(maCandidate, &fragmentInfo[i], fragSequence, fragQuality);
		  }
		}
	  }
	}
  }
}



// idea from Art: maybe we shouldn't place any fragment that has a link to a placed fragment at this 
// point - it must not be consistent with any surrogate since we placed all those already

int PickBestFragment(fragmentSummary *fragmentSum, int numFragmentsToPlace, fragmentInformation *fragmentInfo, int numFragmentInfos,
					 int *indexOfBestFragmentSum, int *indexOfBestFragmentInfo)
{
  int 
	i, 
	j,
	tmpIndexOfBestFragmentInfo = -1;
  double 
	bestScore = 2.0;

  fprintf(stderr, "in PickBestFragment\n");
  
  *indexOfBestFragmentInfo = -1;
  *indexOfBestFragmentSum = -1;

  // we're doing the scoring right here, need something more complicated

  // look through all the fragments
  for (i = 0; i < numFragmentsToPlace; i++)
  {
	// look at all the places each frag can go
	for (j = 0; j < fragmentSum[i].count; j++)
	{
	  // if the score of the frag at this place is the best so far for this frag, mark it
	  if (fragmentSum[i].score >= fragmentInfo[fragmentSum[i].start + j].bayesianQuality)
	  {
		fragmentSum[i].score = fragmentInfo[fragmentSum[i].start + j].bayesianQuality;
		fragmentSum[i].candidateChunkID = fragmentInfo[fragmentSum[i].start + j].chunkID;
		tmpIndexOfBestFragmentInfo = j;
	  }
	}
	// if this is the best fragment of all so far, mark its position in both arrays
	if (fragmentSum[i].score < bestScore)
	{
	  bestScore = fragmentSum[i].score;
	  *indexOfBestFragmentSum = i;
	  *indexOfBestFragmentInfo = tmpIndexOfBestFragmentInfo;
	}
  }
  fprintf (stderr, "indexOfBestFragmentSum = %d\n", *indexOfBestFragmentSum);
  fprintf (stderr, "indexOfBestFragmentInfo = %d\n", *indexOfBestFragmentInfo);
  return (1);
}


int PlaceFragment(int fragmentIndex, fragmentSummary *fragmentSum)
{

  fprintf(stderr, "in PlaceFragment, fragmentSum[fragmentIndex].score = %f\n",
          fragmentSum[fragmentIndex].score);

  if (fragmentSum[fragmentIndex].score < 0.5)
  {
	// AssignFragsToResolvedCI(ScaffoldGraph->CIGraph, parentChunk->id, bestPlacementChunk->id, frag);
	fragmentSum[fragmentIndex].placed = TRUE;
	return 1;
  }
  else
	return 0;
}


// Saul says: perhaps we can do this by deleting the surrogate/inserting the parent in the u_list of the contig
// this avoids the construction of a multi-alignment for the surrogate after all the chunks have been transferred
// that just essentially reconstructs the parent

int PlaceFragmentsInSingletonSurrogates()
{
  //
  // if just one surrogate of a parent, place all the fragments in that surrogate
  //
  int
    i,
	placed;
  chunk_subgraph
    * f;  
  Gap_Fill_t
    * fc;
  ChunkInstanceT
    * parentChunk;
  GraphNodeIterator
	CIGraphIterator;

  fprintf(stderr, "in PlaceFragmentsInSingletonSurrogates\n");

  // scan all the scaffolds
  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while (parentChunk = NextGraphNodeIterator(&CIGraphIterator))
  {
    // scan all the nodes (chunks)
	// find those with exactly one surrogate
	fprintf(stderr, "\nparentChunk " F_CID " has %d instances\n",
                parentChunk->id, parentChunk->info.CI.numInstances);
	if (parentChunk->info.CI.numInstances == 1) 
	{
	  MultiAlignT *maParent;
	  ChunkInstanceT *candidateChunk;
	  int index = getChunkInstanceID(parentChunk, 0);
	  
	  candidateChunk = GetGraphNode(ScaffoldGraph->CIGraph, index);
	  assert (candidateChunk);

	  fprintf(stderr, "parentChunk " F_CID " has %d fragments\n", 
			  parentChunk->id, parentChunk->info.CI.numFragments);
	  fprintf(stderr, "parentChunk " F_CID " scaffoldID: " F_CID "\n",
                  parentChunk->id, parentChunk->scaffoldID);
	  fprintf(stderr, "placing all\n");
	  
	  // Now cycle through the fragments in a parentChunk.  Check its relative fitness to go with 
      // each instance based on matepair consistency.
	  
	  // get the MultiAlign for this parentChunk
	  maParent = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, parentChunk->id);
	  assert (maParent);
	 
	  // cycle through fragments 
	  for (i = 0; i < GetNumIntMultiPoss(maParent->f_list); i++)
	  {
		IntMultiPos *mpParent = GetIntMultiPos(maParent->f_list, i);
		CIFragT 
		  *frag = GetCIFragT(ScaffoldGraph->CIFrags, (CDS_CID_t)mpParent->source);
		CDS_CID_t fragID = (CDS_CID_t)mpParent->source;

		// AssignFragsToResolvedCI(ScaffoldGraph->CIGraph, parentChunk->id, candidateChunk->id, frag);


	  }
	}
  }
  fprintf(stderr, "leaving PlaceFragmentsInSingletonSurrogates\n");
  return 0;
}


int PlaceFragmentsWithPlacedMates()
{
  //
  // place all fragments who have placed mates and an alignment with a quality above the threshold
  //
  int
    i,
	placed;
  chunk_subgraph
    * f;  
  Gap_Fill_t
    * fc;
  ChunkInstanceT
    * parentChunk;
  LengthT
    gap;
  GraphNodeIterator
	CIGraphIterator;
  int
	begin = 0,
	end = 100,
	scaffoldID = 0;

  // scan all the scaffolds
  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while (parentChunk = NextGraphNodeIterator(&CIGraphIterator))
  {
    // scan all the nodes (chunks)
	// if numInstances > 0 then it has at least one surrogate
	if (parentChunk->info.CI.numInstances > 0) // can't use parentChunk->scaffoldID == -1, we might place parents later
	{
	  MultiAlignT *maParent;
	  
	  fprintf(stderr, "\nparentChunk " F_CID " has %d instances\n",
                  parentChunk->id, parentChunk->info.CI.numInstances);
	  fprintf(stderr, "parentChunk " F_CID " has %d fragments\n",
                  parentChunk->id, parentChunk->info.CI.numFragments);
	  fprintf(stderr, "parentChunk " F_CID " scaffoldID: " F_CID "\n",
                  parentChunk->id, parentChunk->scaffoldID);
	  
	  // Now cycle through the fragments in a parentChunk.  Check its relative fitness to go with 
      // each instance based on matepair consistency.
	  
	  // get the MultiAlign for this parentChunk
	  maParent = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, parentChunk->id);
	  assert (maParent);
	 
	  // cycle through fragments 
	  for (i = 0; i < GetNumIntMultiPoss(maParent->f_list); i++)
	  {
		IntMultiPos *mpParent = GetIntMultiPos(maParent->f_list, i);
		CIFragT 
		  *mateFrag, 
		  *frag = GetCIFragT(ScaffoldGraph->CIFrags, (CDS_CID_t)mpParent->source);
		CDS_CID_t fragID = (CDS_CID_t)mpParent->source;
		ChunkInstanceT * mateChunk;

		totalFragmentsAttempted++;
		
		fprintf( stderr, "parentChunk " F_CID " fragID: " F_CID "\n",
                         parentChunk->id, fragID);
		fprintf( stderr, "fragID: " F_CID " frag->iid: " F_CID "\n",
                         fragID, frag->iid);
		fprintf( stderr, "       frag " F_CID " mateOf: " F_CID "\n",
                         fragID, frag->mateOf);
		
		if (frag->numLinks == 1 && frag->mateOf != NULLINDEX)
		{
		  char orientation[STR_LEN];

		  mateFrag = GetCIFragT(ScaffoldGraph->CIFrags, frag->mateOf);
		  mateChunk = GetGraphNode(ScaffoldGraph->CIGraph, mateFrag->CIid);
		  assert(mateChunk != NULL);
		  fprintf( stderr, "       frag " F_CID " placed: %d\n",
                           frag->mateOf,
                           mateFrag->flags.bits.isPlaced == TRUE);
		  fprintf( stderr, "       fragment " F_CID " is in chunk: " F_CID "\n",
                           frag->mateOf, mateChunk->id);
		  fprintf( stderr, "       mateChunk: " F_CID " is in scaffoldID " F_CID "\n",
                           mateChunk->id, mateChunk->scaffoldID);
		  fprintf( stderr,
                           "       mateChunk: " F_CID " has %d instances\n",
                           mateChunk->id, mateChunk->info.CI.numInstances);

		  if (mateFrag->flags.bits.isPlaced)
		  {
			// now cycle through all surrogates of the fragment we are trying to place, picking
			// the best of the acceptable (if any) placements

			placed = PlaceFragmentBasedOnMate(frag, parentChunk, mateFrag, mateChunk);
			attemptedBasedOnMate++;
			if (placed == 1) 
			{
			  placedBasedOnMate++;
			  fprintf(stderr, "pointpoint: frag->iid was just placed 1: " F_CID "\n", frag->iid);
			}
		  }
		}
	  }
	}
  }
  return 0;
}


int PlaceFragmentsWithUnplacedMates()
{
  //
  // place all fragments who have placed mates and an alignment with a quality above the threshold
  //
  int
    i,
	placed;
  chunk_subgraph
    * f;  
  Gap_Fill_t
    * fc;
  ChunkInstanceT
    * parentChunk;
  LengthT
    gap;
  GraphNodeIterator
	CIGraphIterator;
  int
	begin = 0,
	end = 100,
	scaffoldID = 0;

  // scan all the scaffolds
  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while (parentChunk = NextGraphNodeIterator(&CIGraphIterator))
  {
    // scan all the nodes (chunks)
	// if numInstances > 0 then it has at least one surrogate
	if (parentChunk->info.CI.numInstances > 0) // || parentChunk->scaffoldID == -1
	{
	  MultiAlignT *maParent;
	  
	  fprintf(stderr, "\nparentChunk " F_CID " has %d instances\n",
                  parentChunk->id, parentChunk->info.CI.numInstances);
	  fprintf(stderr, "parentChunk " F_CID " has %d fragments\n",
                  parentChunk->id, parentChunk->info.CI.numFragments);
	  fprintf(stderr, "parentChunk " F_CID " scaffoldID: " F_CID "\n",
                  parentChunk->id, parentChunk->scaffoldID);
	  
	  // Now cycle through the fragments in a parentChunk.  Check its relative fitness to go with 
      // each instance based on matepair consistency.
	  
	  // get the MultiAlign for this parentChunk
	  maParent = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, parentChunk->id);
	  assert (maParent);
	 
	  // cycle through fragments 
	  for (i = 0; i < GetNumIntMultiPoss(maParent->f_list); i++)
	  {
		IntMultiPos *mpParent = GetIntMultiPos(maParent->f_list, i);
		CIFragT 
		  *mateFrag, 
		  *frag = GetCIFragT(ScaffoldGraph->CIFrags, (CDS_CID_t)mpParent->source);
		CDS_CID_t fragID = (CDS_CID_t)mpParent->source;
		ChunkInstanceT * mateChunk;

		totalFragmentsAttempted++;
		
		fprintf( stderr, "parentChunk " F_CID " fragID: " F_CID "\n",
                         parentChunk->id, fragID);
		fprintf( stderr, "fragID: " F_CID " frag->iid: " F_CID "\n",
                         fragID, frag->iid);
		fprintf( stderr, "       frag " F_CID " mateOf: " F_CID "\n",
                         fragID, frag->mateOf);
		
		if (frag->numLinks == 1 && frag->mateOf != NULLINDEX)
		{
		  char orientation[STR_LEN];

		  mateFrag = GetCIFragT(ScaffoldGraph->CIFrags, frag->mateOf);
		  mateChunk = GetGraphNode(ScaffoldGraph->CIGraph, mateFrag->CIid);
		  assert(mateChunk != NULL);
		  fprintf( stderr, "       frag " F_CID " placed: %d\n",
                           frag->mateOf,
                           mateFrag->flags.bits.isPlaced == TRUE);
		  fprintf( stderr, "       fragment " F_CID " is in chunk: " F_CID "\n",
                           frag->mateOf, mateChunk->id);
		  fprintf( stderr, "       mateChunk: " F_CID " is in scaffoldID " F_CID "\n",
                           mateChunk->id, mateChunk->scaffoldID);
		  fprintf( stderr, "       mateChunk: " F_CID " has %d instances\n",
                           mateChunk->id, mateChunk->info.CI.numInstances);

		  if (!mateFrag->flags.bits.isPlaced)
		  {
			fprintf(stderr, "matepair: mateChunk->scaffoldID: " F_CID ", parentChunk->scaffoldID: " F_CID "\n", 
                                mateChunk->scaffoldID,
                                parentChunk->scaffoldID);
			
			if (mateChunk->scaffoldID == parentChunk->scaffoldID  // these should both be -1
				&& mateChunk->info.CI.numInstances > 0)  // don't place right now except in surrogates
			{
			  placed = PlaceMatePair(frag, parentChunk, mateFrag, mateChunk);
			  attemptedMatePairs++;

			  if (placed == 1)
			  {
				placedMatePairs++;
				fprintf(stderr, "pointpoint: frag->iid was just placed 3: " F_CID " and " F_CID "\n",
                                        frag->iid, mateFrag->iid);
			  }
			}
		  }
		}
	  }
	}
  }
  return 0;
}


int PlaceMatePair(CIFragT *frag,
                  ChunkInstanceT *parentChunk,
                  CIFragT *mateFrag,
                  ChunkInstanceT *mateChunk)
{
  CDS_CID_t
	mateChunkID,
	parentChunkID;
  int 
	i, 
	j, 
        placed = 0;
  ChunkInstanceT 
	* candidateMateChunk,
	* candidateParentChunk;
  
  fprintf( stderr,
           "in PlaceMatePair, attempting with frags " F_CID " and " F_CID "\n",
           frag->iid, mateFrag->iid);
  
  // we need to account for the fact that mate pairs can be in surrogates and non-surrogates (all chunks with scaffoldID -1)
  
  i = 0;
  while ((mateChunkID = getChunkInstanceID(mateChunk, i)) != -1 && placed == 0)
  {
	i++;
	candidateMateChunk = GetGraphNode(ScaffoldGraph->CIGraph, mateChunkID);
	assert (candidateMateChunk);
	
	fprintf( stderr, "in PlaceMatePair, attempting PlaceFragmentBasedOnMate with chunks " F_CID " and " F_CID "\n", 
			 parentChunk->id, candidateMateChunk->id);

	placed = PlaceFragmentBasedOnMate(frag, parentChunk, mateFrag, candidateMateChunk);
	
	fprintf( stderr, "in PlaceMatePair, result of PlaceFragmentBasedOnMate with chunks " F_CID " and " F_CID ": %d\n", 
			 parentChunk->id, candidateMateChunk->id, placed);
  }
  
  if (placed == 1)  // frag has been placed in the appropriate surrogate by PlaceFragmentBasedOnMate
  {
	fprintf( stderr, "just placed mate pair " F_CID ", " F_CID "\n",
                 frag->id, mateFrag->id);
	// AssignFragsToResolvedCI(ScaffoldGraph->CIGraph, mateChunk->id, candidateChunk->id, frag);
  }

  // If we are here, no pair-wise placement is available based on where the two chunks have surrogates (if anywhere).
  // The only hope is that mateFrag has no surrogates, then we are allowed to place mateFrag based on frag's possible 
  // positions, i.e., see if there is an alignment with the consensus for mateFrag in locations based on frag.
  // Then we assign frag to its surrogate, mateFrag to whichever chunk is in the consensus at that point.
  // But do we need to account for all of the fragments in mateChunk, i.e., see if the whole chunk/contig aligns 
  // with the consensus where we want to place mateFrag?  I think so, so therefore we can count on the chunk being
  // placed as a whole when we assign all chunks with a scaffoldID == -1 in PlaceUnscaffoldedContigs().

  // Of course, it might not be possible to assign a frag to a contig that does not at least contain
  // a surrogate of the frag's parent.  Ask Saul/Karin.

  return placed;
}


// see if frag placement into a surrogate of parentChunk is supported by placed mate in mateChunk
int PlaceFragmentBasedOnMate(CIFragT *frag, ChunkInstanceT *parentChunk, CIFragT *mateFrag, ChunkInstanceT *mateChunk)
{
  CDS_COORD_t
	candidate_chunk_left_end,
	candidate_chunk_right_end,
	mate_chunk_left_end,
        mate_chunk_right_end;
  int 
	chunkScaffoldOrientation,
	index,
	j;

  ChunkInstanceT *candidateChunk,
	*bestCandidateChunk = NULL;
  
  // cycle through all the instances of the parent chunk of the frag we are trying to place
  // we are trying to find a chunk in which it matches well with its mate

  j = 0;
  while ((index = getChunkInstanceID(parentChunk, j)) != -1)
  {
	MultiAlignT *maCandidate;
	DistT scaffoldGap;
	NodeOrient lorient, rorient;

	j++;
	candidateChunk = GetGraphNode(ScaffoldGraph->CIGraph, index);
	assert (candidateChunk);
	
	fprintf( stderr, "       candidateChunk: " F_CID " is in scaffoldID " F_CID "\n", 
			 candidateChunk->id, candidateChunk->scaffoldID);
	
	// for now, we want the surrogate and the placed chunk to be in the same scaffold
	if (mateChunk->scaffoldID == candidateChunk->scaffoldID && mateChunk->scaffoldID != -1)   
	{
	  GraphEdgeStatT stat;
	  DistT *fragDist;
	  LengthT distance;
	  
	  lorient = GetNodeOrient(candidateChunk);
	  rorient = GetNodeOrient(mateChunk);

	  // need to add in the offsets of their containing contigs to figure out scaffold coords
	  GetChunkPositionInScaffold(candidateChunk, &candidate_chunk_left_end, &candidate_chunk_right_end, 
								 &chunkScaffoldOrientation);

	  GetChunkPositionInScaffold(mateChunk, &mate_chunk_left_end, &mate_chunk_right_end, 
								 &chunkScaffoldOrientation);

	  fprintf(stderr,
			  "  positions in scaffold: candidateChunk (" F_CID ") [" F_COORD "," F_COORD "] o:%s - mateChunk (" F_CID ") [" F_COORD "," F_COORD "] o:%s\n",
			  candidateChunk->id,
			  candidate_chunk_left_end,
			  candidate_chunk_right_end,
			  (lorient == A_B) ? "A_B" : "B_A",
			  mateChunk->id,
			  mate_chunk_left_end,
			  mate_chunk_right_end,
			  (rorient == A_B) ? "A_B" : "B_A");

	  if (candidateChunk->id == mateChunk->id)
		scaffoldGap.mean = 0.0;
	  else
	  {
		if (candidate_chunk_left_end < mate_chunk_left_end)
		{
		  // fprintf(stderr, "lchunk is candidate, rchunk is mate\n");
		  scaffoldGap.mean = mate_chunk_left_end - candidate_chunk_right_end;
		}
		else 
		{
		  // fprintf(stderr, "lchunk is mate, rchunk is candidate\n");
		  scaffoldGap.mean = candidate_chunk_left_end - mate_chunk_right_end;
		}
	  }
	  scaffoldGap.stddev = 0.0;  // really don't need variance, since we are going to judge valid
                                 // links by the variance in the library for that type of link

	  fprintf(stderr, "According to us, the scaffoldGap.mean = %f and scaffoldGap.stddev = %f\n",
			  scaffoldGap.mean, scaffoldGap.stddev);

	  // we set AS_GKP_INNIE since we are dealing with a mate pair, and always want the offset
	  // to be from the 5' end when calling FragOffsetAndOrientation in SurrogateChunkDistance
	  SurrogateChunkDistance(ScaffoldGraph->CIGraph, mateFrag, frag, &scaffoldGap, mateFrag->type,
							 AS_GKP_INNIE, &stat, &distance);

	  fprintf(stderr, "According to us, the distance.mean = %f and distance.variance = %f\n",
			  distance.mean, distance.variance);

	  fragDist = GetDistT(ScaffoldGraph->Dists, frag->dist);

	  fprintf(stderr, "According to Dists, the fragDist->mean = %f and fragDist->stddev = %f\n",
			  fragDist->mean, fragDist->stddev);
	  fprintf(stderr, "The difference in the two estimates is %f, and 3.0 * fragDist->stddev is %f\n",
			  abs (distance.mean - fragDist->mean), 3.0 * fragDist->stddev);

	  if (abs (distance.mean - fragDist->mean) < 3.0 * fragDist->stddev)
	  {
		// what if more fragment is consistent with mate in more than one surrogate?
		if (1)
		  bestCandidateChunk = candidateChunk;
		
	  }
	  else
	  {
		fprintf(stderr, "we CAN NOT place fragment " F_CID " in surrogate " F_CID " based on mate\n", frag->iid, candidateChunk->id);	
	  }
	}
  }

  if (bestCandidateChunk != NULL)
  {
	fprintf(stderr, "we CAN place fragment with iid " F_CID " in surrogate " F_CID " based on mate\n", 
			frag->iid, candidateChunk->id);	
	Print_Mates_Cam(frag, candidateChunk->id, mateFrag, mateChunk->id);
	



	// DO WE NEED TO CHECK QUALITY OF PLACEMENT???????
	
	// AssignFragsToResolvedCI(ScaffoldGraph->CIGraph, parentChunk->id, candidateChunk->id, frag);
	// only try to place once, and taking the last consistent place
	return(1);
  }
  else
	return(0);
}

// this needs to have a subroutine that calculates based on surrogates
// then when this function cycles through a parent, it calls this routine
// the other functionality is needed elsewhere in the code when we want to update individual
// surrogates who have new fragment overlaps
// old
int ScoreFragmentBasedOnOverlapsOld(CIFragT *frag, ChunkInstanceT *parentChunk, 
								 fragmentInformation *fragmentInfo, int *numFragmentInfos)
{
  CDS_CID_t
        chunkIndex;
  CDS_COORD_t
	candidate_chunk_left_end,
	candidate_chunk_right_end,
	bestPlacementAhg,
	bestPlacementBhg,
	bgn, 
	end,
	frag_left_end,
	frag_right_end,
	chunk_left_end,
        chunk_right_end;
  int 
	j,
	fragContigOrientation,
	chunkContigOrientation;
  float 
	bestPlacementQuality = 1.1;  // 1.0 is max quality value, the lower the better
  ChunkInstanceT 
	*bestPlacementChunk;
  SeqInterval 
	range;
  char *untrimmedSequence,
	*untrimmedQuality;
  VA_TYPE(char) *consensusSequence = NULL;
  VA_TYPE(char) *consensusQuality = NULL;
  VA_TYPE(int) *consensusCoverage = NULL;

  fprintf( stderr, "in ScoreFragmentBasedOnOverlaps\n");
  
  fprintf( stderr, "frag " F_CID ", numInstances of parentChunk " F_CID ": %d\n", 
		   frag->iid, parentChunk->id, parentChunk->info.CI.numInstances);
  fprintf( stderr, "parentChunk " F_CID " is in scaffold: " F_CID "\n",
           parentChunk->id, parentChunk->scaffoldID);

  j = 0;
  while ((chunkIndex = getChunkInstanceID(parentChunk, j++)) != -1)
  {
	MultiAlignT *maCandidate, *maParent;
	ChunkInstanceT *candidateChunk;
	char fragSequence[CNS_MAX_SEQUENCE_LENGTH], fragQuality[CNS_MAX_SEQUENCE_LENGTH];
	
	candidateChunk = GetGraphNode(ScaffoldGraph->CIGraph, chunkIndex);
	assert (candidateChunk);

	fprintf( stderr, "       candidateChunk: " F_CID " is in scaffoldID " F_CID "\n", 
			 candidateChunk->id, candidateChunk->scaffoldID);

	fragmentInfo[*numFragmentInfos].chunkID = candidateChunk->id;
	fragmentInfo[*numFragmentInfos].fragID = frag->iid;
	
	// grab the frags clear sequence
	getFragmentSequenceAndQuality( frag->iid, fragSequence, fragQuality);
	
	if(consensusSequence == NULL)
	{
	  consensusSequence = CreateVA_char(2048);
	  consensusQuality = CreateVA_char(2048);
	  consensusCoverage = CreateVA_int(2048);
	}
	
	// get the multiAlign for the contig
	maCandidate = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, candidateChunk->info.CI.contigID);
	assert (maCandidate);
	
	// find out where the surrogate is in the contig
	GetChunkGappedPositionInContig(maCandidate, candidateChunk, &chunk_left_end, &chunk_right_end, 
								   &chunkContigOrientation);
	fprintf( stderr, "chunk position in contig: chunk_left_end: " F_COORD ", chunk_right_end: " F_COORD ", orientation: %d\n", 
                 chunk_left_end, chunk_right_end,
                 chunkContigOrientation);
	
	// get the multiAlign for the parent
	maParent = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, parentChunk->id);
	assert (maParent);
	
	// now we can get frags position in parent to deduce location in surrogate which is in the contig
	GetFragmentGappedPositionInContigFromChunk(maParent, frag, &frag_left_end, &frag_right_end, 
											   &fragContigOrientation, 
											   chunkContigOrientation, chunk_left_end, chunk_right_end);
	fprintf( stderr, "frag position in contig: frag_left_end: " F_COORD ", frag_right_end: " F_COORD ", orientation: %d\n", 
                 frag_left_end, frag_right_end, fragContigOrientation);
	
	fragmentInfo[*numFragmentInfos].contigPosition.bgn = frag_left_end;
	fragmentInfo[*numFragmentInfos].contigPosition.end = frag_right_end;

	range.bgn = frag_left_end;
	if (range.bgn < 0) range.bgn = 0;
	
	range.end = frag_right_end;
	// subtract 1 because numElements counts a null character at the end of the consensus string
	if (range.end > maCandidate->consensus->numElements - 1) range.end = maCandidate->consensus->numElements - 1;
	
	// these are used in DP_Compare_AS below
	bgn = 0;
	end = range.end - range.bgn;
	
	fprintf( stderr,
                 "the new way2: maCandidate->consensus->numElements: %d\n",
                 maCandidate->consensus->numElements);
	fprintf( stderr, "the new way2: range.bgn: " F_COORD ", range.end: " F_COORD ", |range| = " F_COORD "\n", range.bgn, range.end, range.end - range.bgn);
	
	GetCoverageInMultiAlignTFast(maCandidate, range, consensusCoverage);
	
	// overlap 'em
	{
	  CDS_COORD_t A_hang_estimate,A_hang_lower,A_hang_upper;
	  InternalFragMesg AFR, BFR;
	  OverlapMesg *O = NULL;
	  int where, opposite = FALSE;
	  FILE* log = NULL;
	  CDS_COORD_t length;
	  char *tempAFRsequence, *tempAFRquality;
	  
	  if (fragContigOrientation == 1) opposite = TRUE;
	  
	  tempAFRsequence = Getchar(maCandidate->consensus, range.bgn);
	  tempAFRquality = Getchar(maCandidate->quality, range.bgn);
	  AFR.sequence = tempAFRsequence;
	  AFR.sequence[range.end + 1] = '\0';
	  AFR.quality = tempAFRquality;
	  AFR.quality[range.end + 1] = '\0';
	  
	  AFR.eaccession = 0;
	  AFR.iaccession = 0;
	  BFR.sequence   = fragSequence;
	  BFR.quality    = fragQuality;
	  BFR.eaccession = 0;
	  BFR.iaccession = 1;		
	  
	  // now we need to mask out quality where coverage is zero
	  fragmentInfo[*numFragmentInfos].coverage = maskQuality(AFR.quality, consensusCoverage, range.end - range.bgn);
	  fprintf(stderr, "coverage for frag " F_CID " = %d\n",
                  frag->iid, fragmentInfo[*numFragmentInfos].coverage);
	  
#if 1
	  fprintf(stderr,"* ScoreFragmentBasedOnOverlaps, Calling DP_Compare_AS with %s sequences lengths " F_SIZE_T "," F_SIZE_T " bgn,end =[" F_COORD "," F_COORD "]\n sequence: %s\n fragment: %s\n",
			  (!strcmp(AFR.sequence, BFR.sequence)?"EQUAL":"DIFFERENT"),
			  strlen(AFR.sequence), strlen(BFR.sequence),
			  bgn, end, AFR.sequence, BFR.sequence);
#endif

	  if (fragmentInfo[*numFragmentInfos].coverage != 0)
	  {
		//CGW_DP_ERATE=0.1 changed parameter to canOlap->errRate
		// a new field in ChunkOverlapCheckT
		O = DP_Compare_AS(&AFR, &BFR, 
						  bgn, end, opposite,
						  CGW_DP_ERATE, CGW_DP_THRESH, CGW_DP_MINLEN,
						  AS_FIND_ALIGN, //_NO_TRACE, // Slower, computes an alignment, throw away delta encoding
						  &where);
	  }

	  if (O != NULL)
	  {
		myPrint_Overlap_AS(stderr, &AFR, &BFR, O);
		if(log == NULL)
		  log = stderr;			
		
		// we compute quality with a quality threshold of 1 - this masks those positions that have quality = 0, eg, 
		// base pairs with no fragment coverage
		compute_bayesian_quality(&AFR, &BFR, O, 1, &length, log);
		fprintf(stderr, "In ScoreFragmentBasedOnOverlaps, Quality between WITHOUT QV Realigning = %lf\n", O->quality);
		fragmentInfo[*numFragmentInfos].bayesianQuality = O->quality;
		if (O->quality < bestPlacementQuality)
		{
		  bestPlacementQuality = O->quality;
		  bestPlacementChunk = candidateChunk;
		  bestPlacementAhg = O->ahg;
		  bestPlacementBhg = O->bhg;			
		}
	  }
	  else  // there was no overlap or insufficient coverage 
		fragmentInfo[*numFragmentInfos].bayesianQuality = 2.0;  // this is safely off the quality scale of [0.0, 1.0]
	}
	DeleteVA_char(consensusSequence); 
	DeleteVA_char(consensusQuality); 
	DeleteVA_int(consensusCoverage); 
	(*numFragmentInfos)++;
  }

  // if best quality is less than threshold place
  if (bestPlacementQuality < 0.5) // does Knut have a #define we can borrow?
  {
	fprintf(stderr, "we CAN place fragment %d in surrogate %d based on overlaps\n", 
			frag->iid, bestPlacementChunk->id);
	fprintf(stderr, "we are placing fragment at %d, %d\n", bestPlacementAhg, range.end - bestPlacementBhg);
	
	// AssignFragsToResolvedCI(ScaffoldGraph->CIGraph, parentChunk->id, bestPlacementChunk->id, frag);
	Print_Fragment_Cam(frag, bestPlacementChunk->id);
	return(1);
  }  
  else
  {
	fprintf(stderr, "we CAN NOT place fragment %d in any surrogate of %d\n", frag->iid, parentChunk->id);	
	return(0);
  }
}

// the new way
// this function also initializes the fragmentInfo data structure
int ScoreFragmentBasedOnOverlaps(CIFragT *frag, ChunkInstanceT *parentChunk,
								 fragmentInformation *fragmentInfo, int *numFragmentInfos)
{
  int 
	candidate_chunk_left_end,
	candidate_chunk_right_end,
	j,
	chunkIndex,
	frag_left_end,
	frag_right_end,
	fragContigOrientation,
	chunk_left_end,
	chunk_right_end,
	chunkContigOrientation,
	ableToPlace = 0;
  char 
	fragSequence[CNS_MAX_SEQUENCE_LENGTH], 
	fragQuality[CNS_MAX_SEQUENCE_LENGTH];
  MultiAlignT 
	*maParent;

  fprintf( stderr, "in ScoreFragmentBasedOnOverlaps\n");
  
  fprintf( stderr, "frag %d, numInstances of parentChunk %d: %d\n", 
		   frag->iid, parentChunk->id, parentChunk->info.CI.numInstances);
  fprintf( stderr, "parentChunk %d is in scaffold: %d\n", parentChunk->id, parentChunk->scaffoldID);

  // grab the frags clear sequence
  getFragmentSequenceAndQuality( frag->iid, fragSequence, fragQuality);
	
  // get the multiAlign for the parent
  maParent = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, parentChunk->id);
  assert (maParent);
	
  j = 0;
  while ((chunkIndex = getChunkInstanceID(parentChunk, j++)) != -1)
  {
	MultiAlignT *maCandidate;
	ChunkInstanceT *candidateChunk;
	
	candidateChunk = GetGraphNode(ScaffoldGraph->CIGraph, chunkIndex);
	assert (candidateChunk);

	fprintf( stderr, "       candidateChunk: %ld is in scaffoldID %d\n", 
			 candidateChunk->id, candidateChunk->scaffoldID);

	// get the multiAlign for the contig
	maCandidate = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, candidateChunk->info.CI.contigID);
	assert (maCandidate);
	
	// find out where the surrogate is in the contig
	GetChunkGappedPositionInContig(maCandidate, candidateChunk, &chunk_left_end, &chunk_right_end, 
								   &chunkContigOrientation);
	fprintf( stderr, "chunk position in contig: chunk_left_end: %d, chunk_right_end: %d, orientation: %d\n", 
			 chunk_left_end, chunk_right_end, chunkContigOrientation);
	
	// now we can get frags position in parent to deduce location in surrogate which is in the contig
	GetFragmentGappedPositionInContigFromChunk(maParent, frag, &frag_left_end, &frag_right_end, 
											   &fragContigOrientation, 
											   chunkContigOrientation, chunk_left_end, chunk_right_end);
	fprintf( stderr, "frag position in contig: frag_left_end: %d, frag_right_end: %d, orientation: %d\n", 
			 frag_left_end, frag_right_end, fragContigOrientation);


	// now initialize fragmentInfo[*numFragmentInfos] with the information we have gathered
	fragmentInfo[*numFragmentInfos].chunkID = candidateChunk->id;
	fragmentInfo[*numFragmentInfos].fragID = frag->iid;
	fragmentInfo[*numFragmentInfos].fragContigOrientation = fragContigOrientation;
	fragmentInfo[*numFragmentInfos].contigID = candidateChunk->info.CI.contigID;
	
	assert (frag_left_end >= 0);
	fragmentInfo[*numFragmentInfos].contigPosition.bgn = frag_left_end;
	
	// subtract 1 because numElements counts a null character at the end of the consensus string
	assert ( frag_right_end <= maCandidate->consensus->numElements - 1);
	fragmentInfo[*numFragmentInfos].contigPosition.end = frag_right_end;
	
	fprintf( stderr, "the new way2: maCandidate->consensus->numElements: %d\n", maCandidate->consensus->numElements);
	fprintf( stderr, "the new way2: frag_left_end: %d, frag_right_end: %d, |range| = %d\n", 
			 frag_left_end, frag_right_end, frag_right_end - frag_left_end);
	
	// find out if we can place this frag in this surrogate
	ableToPlace = ScoreFragmentBasedOnOverlapsInSurrogate(maCandidate, &fragmentInfo[*numFragmentInfos],
														  fragSequence, fragQuality);
	(*numFragmentInfos)++;
  }

  if (ableToPlace == 0)
	fprintf(stderr, "we CAN NOT place fragment %d in any surrogate of %d\n", frag->iid, parentChunk->id);	
  return(ableToPlace);
}


int ScoreFragmentBasedOnOverlapsInSurrogate(MultiAlignT *maCandidate, fragmentInformation *fragmentInfo,
											char *fragSequence, char *fragQuality)
{
  int ableToPlace = 0;
  VA_TYPE(int) *consensusCoverage = CreateVA_int(2048);
  
  GetCoverageInMultiAlignTFast(maCandidate, fragmentInfo->contigPosition, consensusCoverage);
	  
  // overlap 'em
  {
	int A_hang_estimate,A_hang_lower,A_hang_upper;
	InternalFragMesg AFR, BFR;
	OverlapMesg *O = NULL;
	int where, opposite = FALSE;
	FILE* log = stderr;
	int length;
	
	if (fragmentInfo->fragContigOrientation == 1) opposite = TRUE;
	
	AFR.sequence = Getchar(maCandidate->consensus, fragmentInfo->contigPosition.bgn);
	AFR.sequence[fragmentInfo->contigPosition.end + 1] = '\0';
	AFR.quality = Getchar(maCandidate->quality, fragmentInfo->contigPosition.bgn);
	AFR.quality[fragmentInfo->contigPosition.end + 1] = '\0';
	AFR.eaccession = 0;
	AFR.iaccession = 0;
	
	BFR.sequence   = fragSequence;
	BFR.quality    = fragQuality;
	BFR.eaccession = 0;
	BFR.iaccession = 1;		
		
	// now we need to mask out quality where coverage is zero
	fragmentInfo->coverage = 
	  maskQuality(AFR.quality, consensusCoverage, 
				  fragmentInfo->contigPosition.end - fragmentInfo->contigPosition.bgn);
	fprintf(stderr, "coverage for frag %d = %d\n", fragmentInfo->fragID, fragmentInfo->coverage);
		
#if 0
	fprintf(stderr, "* ScoreFragmentBasedOnOverlaps, Calling DP_Compare_AS with %s sequences lengths %ld,%ld bgn,end =[%d,%d]\n",
			(!strcmp(AFR.sequence, BFR.sequence)?"EQUAL":"DIFFERENT"),
			strlen(AFR.sequence), strlen(BFR.sequence),
			0, fragmentInfo->contigPosition.end - fragmentInfo->contigPosition.bgn);
	fprintf(stderr, "sequence: %s\n fragment: %s\n", AFR.sequence, BFR.sequence);
#endif
	
	if (fragmentInfo->coverage != 0)
	{
	  //CGW_DP_ERATE=0.1 changed parameter to canOlap->errRate
	  // a new field in ChunkOverlapCheckT
	  O = DP_Compare_AS(&AFR, &BFR, 
						0, fragmentInfo->contigPosition.end - fragmentInfo->contigPosition.bgn, 
						opposite,
						CGW_DP_ERATE, CGW_DP_THRESH, CGW_DP_MINLEN,
						AS_FIND_ALIGN, //_NO_TRACE, // Slower, computes an alignment, throw away delta encoding
						&where);
	}
	
	if (O != NULL)
	{
	  // myPrint_Overlap_AS(stderr, &AFR, &BFR, O);
		  
	  // we compute quality with a quality threshold of 1 - this masks those positions that have quality = 0, eg, 
	  // base pairs with no fragment coverage
	  compute_bayesian_quality(&AFR, &BFR, O, 1, &length, log);
	  fprintf(stderr, "In ScoreFragmentBasedOnOverlaps, Quality between WITHOUT QV Realigning = %lf\n", O->quality);
	  fragmentInfo->bayesianQuality = O->quality;
	  
	  if (O->quality < 0.5) // does Knut have a #define we can borrow?
	  {
		fprintf(stderr, "we CAN place fragment %d in surrogate %d based on overlaps\n", 
				fragmentInfo->fragID, fragmentInfo->chunkID);
		ableToPlace = 1;
	  }
	}
	else  // there was no overlap or insufficient coverage 
	  fragmentInfo->bayesianQuality = 2.0;  // this is safely off the quality scale of [0.0, 1.0]
  }
  DeleteVA_int(consensusCoverage); 

  return(ableToPlace);
}


int maskQuality(char *quality, VA_TYPE(int) *coverage, int length)
{
  int 
	i,
	numBasePairsCovered = 0;
  
  for (i = 0; i < length; i++)
  {
	if (*GetVA_int(coverage, i) != 0)   // if there are real fragments under this base pair
	{
	  quality++;                        // just leave quality alone
	  numBasePairsCovered++;            // we would like this info returned
	}
	else                                // if not
	  *quality++ = '0';                 // set quality to "zero"
  }
  return numBasePairsCovered;
}



#define NUM_COLORS 4

int Print_Mates_Cam(CIFragT *frag, int fragSurrogateID, CIFragT *mateFrag, int mateFragChunkID)
{
  int
    color,
    k,
    o,
    id,
    i;
  CIEdgeT
    * edge;
  int32
    low,
    high,
    left_end,
    right_end,
	frag_left_end, 
	frag_right_end,
	fragScaffoldOrientation,
	mateFrag_left_end, 
	mateFrag_right_end,
    cid,
    other_cid;
  ChunkInstanceT
    * other_chunk,
    * chunk;
  char
    unique[STR_LEN],
    orientation[STR_LEN],
    filename[STR_LEN],
    * Colour[NUM_COLORS] = {
      "CFFFF00 T2 S # real node",
      "C0000F0 T2 S # surrogate_node",
      "CF0F0AA T1 S # fragment",
      "C00AAAA T1 S # CIEdge"};
  FILE
    * cam_file;
  GraphNodeIterator
	CIGraphIterator;
  ContigT 
	* containingContig;

  //
  // open the cam file
  //
  sprintf(filename, "./cam/mates%d_%d.cam", frag->iid, mateFrag->iid);
  cam_file = file_open (filename, "w");
  assert(cam_file != NULL);
  
  //
  // output the colors
  //
  for (i = 0; i < NUM_COLORS; i++)
    fprintf(cam_file, "%d: %s\n",
			i,
			Colour[i]);

  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while (chunk = NextGraphNodeIterator(&CIGraphIterator))
  {
	int chunkScaffoldOrientation;
	
    assert(chunk != NULL);

	if (chunk->scaffoldID != -1)
	{
	  if (chunk->type == DISCRIMINATORUNIQUECHUNK_CGW) 
	  {
		strcpy(unique, "discriminator unique chunk");
		color = 0;
	  } 
	  else if (chunk->type == RESOLVEDREPEATCHUNK_CGW)
	  {
		strcpy(unique, "surrogate chunk");
		color = 1;
	  }
	  
	  GetChunkPositionInScaffold(chunk, &left_end, &right_end, &chunkScaffoldOrientation);	  

	  if (chunkScaffoldOrientation == 0)
		strcpy(orientation, "direct in scaffold");
	  else
		strcpy(orientation, "reverse in scaffold");

	  // fprintf(cam_file, "%f %f %f %f\n", containingContig->offsetAEnd.mean, containingContig->offsetBEnd.mean,
	  //	chunk->offsetAEnd.mean, chunk->offsetBEnd.mean);
		
	  fprintf(cam_file,
			  "%d: %d A%d %d # %s %s chunk %d scaff_id %d\n",
			  chunk->id,
			  left_end,
			  color,
			  right_end,
			  unique,
			  orientation,
			  chunk->id,
			  chunk->scaffoldID);

	  // fprintf(cam_file, "chunk->id, fragSurrogateID: %d, %d\n", chunk->id, fragSurrogateID);

	  if (chunk->id == fragSurrogateID)
	  {
		GetFragmentPositionInScaffoldFromChunk(frag, &frag_left_end, &frag_right_end, 
											   &fragScaffoldOrientation, 
											   chunkScaffoldOrientation, left_end, right_end);
		if (fragScaffoldOrientation == 0)
		  strcpy(orientation, "direct in scaffold");
		else
		  strcpy(orientation, "reverse in scaffold");

		fprintf(cam_file,
				"%d: %d A%d %d # %s %s frag iid %d chunk_id %d surrogate_id %d\n",
				frag->iid,
				frag_left_end,
				2,
				frag_right_end,
				"fragment",
				orientation,
				frag->iid,
				frag->cid,
				fragSurrogateID);
	  }

	  if (chunk->id == mateFragChunkID)
	  {
		GetFragmentPositionInScaffoldFromChunk(mateFrag, &mateFrag_left_end, &mateFrag_right_end, 
											   &fragScaffoldOrientation, 
											   chunkScaffoldOrientation, left_end, right_end);
		if (fragScaffoldOrientation == 0)
		  strcpy(orientation, "direct in scaffold");
		else
		  strcpy(orientation, "reverse in scaffold");

		fprintf(cam_file,
				"%d: %d A%d %d # %s %s mateFrag iid %d chunk_id %d\n",
				mateFrag->iid,
				mateFrag_left_end,
				2,
				mateFrag_right_end,
				"mateFragment",
				orientation,
				mateFrag->iid,
				mateFrag->cid);
	  }
	}
  }
  
  
#if 1
  //
  // draw the edges
  //
  
  {
	int distance;
	
	if (mateFrag_left_end < frag_right_end)
	  distance = frag_right_end - mateFrag_left_end;
  	else
	  distance = mateFrag_right_end - frag_left_end;

	fprintf(cam_file, "LNK: %d %d A2 # edge distance %d\n",
			frag->iid,
			mateFrag->iid,
			distance);
  }


#endif

  fclose (cam_file);
  return(0);
}

int Print_Fragment_Cam(CIFragT *frag, int fragSurrogateID)
{
  int
    color,
    k,
    o,
    id,
    i;
  CIEdgeT
    * edge;
  int32
    low,
    high,
    left_end,
    right_end,
	frag_left_end, 
	frag_right_end,
	fragScaffoldOrientation,
    cid,
    other_cid;
  ChunkInstanceT
    * other_chunk,
    * chunk;
  char
    unique[STR_LEN],
    orientation[STR_LEN],
    filename[STR_LEN],
    * Colour[NUM_COLORS] = {
      "CFFFF00 T2 S # real node",
      "C0000F0 T2 S # surrogate_node",
      "CF0F0AA T1 S # fragment",
      "C00AAAA T1 S # CIEdge"};
  FILE
    * cam_file;
  GraphNodeIterator
	CIGraphIterator;
  ContigT 
	* containingContig;

  //
  // open the cam file
  //
  sprintf(filename, "./cam/fragment%d.cam", frag->iid);
  cam_file = file_open (filename, "w");
  assert(cam_file != NULL);
  
  //
  // output the colors
  //
  for (i = 0; i < NUM_COLORS; i++)
    fprintf(cam_file, "%d: %s\n",
			i,
			Colour[i]);

  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while (chunk = NextGraphNodeIterator(&CIGraphIterator))
  {
	int chunkScaffoldOrientation;
	
    assert(chunk != NULL);

	if (chunk->scaffoldID != -1)
	{
	  if (chunk->type == DISCRIMINATORUNIQUECHUNK_CGW) 
	  {
		strcpy(unique, "discriminator unique chunk");
		color = 0;
	  } 
	  else if (chunk->type == RESOLVEDREPEATCHUNK_CGW)
	  {
		strcpy(unique, "surrogate chunk");
		color = 1;
	  }
	  
	  GetChunkPositionInScaffold(chunk, &left_end, &right_end, &chunkScaffoldOrientation);	  

	  if (chunkScaffoldOrientation == 0)
		strcpy(orientation, "direct in scaffold");
	  else
		strcpy(orientation, "reverse in scaffold");

	  // fprintf(cam_file, "%f %f %f %f\n", containingContig->offsetAEnd.mean, containingContig->offsetBEnd.mean,
	  //	chunk->offsetAEnd.mean, chunk->offsetBEnd.mean);
		
	  fprintf(cam_file,
			  "%d: %d A%d %d # %s %s chunk %d scaff_id %d\n",
			  chunk->id,
			  left_end,
			  color,
			  right_end,
			  unique,
			  orientation,
			  chunk->id,
			  chunk->scaffoldID);

	  // fprintf(cam_file, "chunk->id, fragSurrogateID: %d, %d\n", chunk->id, fragSurrogateID);

	  if (chunk->id == fragSurrogateID)
	  {
		GetFragmentPositionInScaffoldFromChunk(frag, &frag_left_end, &frag_right_end, 
											   &fragScaffoldOrientation, 
											   chunkScaffoldOrientation, left_end, right_end);
		if (fragScaffoldOrientation == 0)
		  strcpy(orientation, "direct in scaffold");
		else
		  strcpy(orientation, "reverse in scaffold");

		fprintf(cam_file,
				"%d: %d A%d %d # %s %s frag iid %d chunk_id %d surrogate_id %d\n",
				frag->iid,
				frag_left_end,
				2,
				frag_right_end,
				"fragment",
				orientation,
				frag->iid,
				frag->cid,
				fragSurrogateID);
	  }
	}
  }

  fclose (cam_file);
  return(0);
}

int Print_Scaffold_Cam(int index)
{
  int
    color,
    k,
    o,
    id,
    i;
  CIEdgeT
    * edge;
  int32
    low,
    high,
    left_end,
    right_end,
    cid,
    other_cid;
  ChunkInstanceT
    * other_chunk,
    * chunk;
  char
    unique[STR_LEN],
    orientation[STR_LEN],
    filename[STR_LEN],
    * Colour[NUM_COLORS] = {
      "CFFFF00 T2 S # real node",
      "C0000F0 T2 S # surrogate_node",
      "CF0F0AA T1 S # fragment",
      "C00AAAA T1 S # CIEdge"};
  FILE
    * cam_file;
  GraphNodeIterator
	ContigGraphIterator;
  ContigT
	* containingContig;
  
  //
  // open the cam file
  //
  sprintf(filename, "./cam/scaffold%d.cam", index);
  cam_file = file_open (filename, "w");
  assert(cam_file != NULL);
  
  //
  // output the colors
  //
  for (i = 0; i < NUM_COLORS; i++)
    fprintf(cam_file, "%d: %s\n",
			i,
			Colour[i]);

  InitGraphNodeIterator(&ContigGraphIterator, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while (chunk = NextGraphNodeIterator(&ContigGraphIterator))
  {
    assert(chunk != NULL);

	strcpy(unique, "discriminator unique chunk");
	color = 0;
	
	if (chunk->offsetAEnd.mean <= chunk->offsetBEnd.mean) 
	{
	  left_end = chunk->offsetAEnd.mean;
	  right_end = chunk->offsetBEnd.mean;
	  strcpy(orientation, "direct in scaffold");
	} 
	else 
	{
	  left_end = chunk->offsetBEnd.mean;
	  right_end = chunk->offsetAEnd.mean;
	  strcpy(orientation, "reverse in scaffold");
	}

	if (chunk->scaffoldID == 0)
	  fprintf(cam_file,
			  "%d: %d A%d %d # %s %s chunk %d scaff_id %d\n",
			  chunk->id,
			  left_end,
			  0,
			  right_end,
			  "contig",
			  orientation,
			  chunk->id,
			  chunk->scaffoldID);
  }
  
  fclose (cam_file);
  return(0);
};


int Print_Scaffold_Cam_CIs(int index)
{
  int
    color,
    k,
    o,
    id,
    i;
  CIEdgeT
    * edge;
  int32
    low,
    high,
    left_end,
    right_end,
    cid,
    other_cid;
  ChunkInstanceT
    * other_chunk,
    * chunk;
  char
    unique[STR_LEN],
    orientation[STR_LEN],
    filename[STR_LEN],
    * Colour[NUM_COLORS] = {
      "CFFFF00 T2 S # real node",
      "C0000F0 T2 S # surrogate_node",
      "CF0F0AA T1 S # fragment",
      "C00AAAA T1 S # CIEdge"};
  FILE
    * cam_file;
  GraphNodeIterator
	CIGraphIterator;
  ContigT
	* containingContig;
  
  //
  // open the cam file
  //
  sprintf(filename, "./cam/scaffold%d.cam", index);
  cam_file = file_open (filename, "w");
  assert(cam_file != NULL);
  
  //
  // output the colors
  //
  for (i = 0; i < NUM_COLORS; i++)
    fprintf(cam_file, "%d: %s\n",
			i,
			Colour[i]);

  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while (chunk = NextGraphNodeIterator(&CIGraphIterator))
  {
    assert(chunk != NULL);

	if (chunk->scaffoldID != -1)
	{
	  if (chunk->type == DISCRIMINATORUNIQUECHUNK_CGW) 
	  {
		strcpy(unique, "discriminator unique chunk");
		color = 0;
	  } 
	  else if (chunk->type == RESOLVEDREPEATCHUNK_CGW)
	  {
		strcpy(unique, "surrogate chunk");
		color = 1;
	  }

	  containingContig = GetGraphNode(ScaffoldGraph->ContigGraph, chunk->info.CI.contigID);
	  
	  if (containingContig->offsetAEnd.mean < containingContig->offsetBEnd.mean)  // contig is normal
	  {
		if (chunk->offsetAEnd.mean <= chunk->offsetBEnd.mean) 
		{
		  left_end = containingContig->offsetAEnd.mean + chunk->offsetAEnd.mean;
		  right_end = containingContig->offsetAEnd.mean + chunk->offsetBEnd.mean;
		  strcpy(orientation, "direct in scaffold");
		} 
		else 
		{
		  left_end = containingContig->offsetAEnd.mean + chunk->offsetBEnd.mean;
		  right_end = containingContig->offsetAEnd.mean + chunk->offsetAEnd.mean;
		  strcpy(orientation, "reverse in scaffold");
		}
	  }
	  else  // contig is reversed in scaffold
	  {
		if (chunk->offsetAEnd.mean <= chunk->offsetBEnd.mean) 
		{
		  left_end = containingContig->offsetAEnd.mean - chunk->offsetBEnd.mean;
		  right_end = containingContig->offsetAEnd.mean - chunk->offsetAEnd.mean;
		  strcpy(orientation, "reverse in scaffold");
		} 
		else 
		{
		  left_end = containingContig->offsetAEnd.mean - chunk->offsetAEnd.mean;
		  right_end = containingContig->offsetAEnd.mean - chunk->offsetBEnd.mean;
		  strcpy(orientation, "direct in scaffold");
		}
	  }
	  
	  // fprintf(cam_file, "%f %f %f %f\n", containingContig->offsetAEnd.mean, containingContig->offsetBEnd.mean,
	  //	chunk->offsetAEnd.mean, chunk->offsetBEnd.mean);
	  
	  fprintf(cam_file,
			  "%d: %d A%d %d # %s %s chunk %d scaff_id %d\n",
			  chunk->id,
			  left_end,
			  color,
			  right_end,
			  unique,
			  orientation,
			  chunk->id,
			  chunk->scaffoldID);
	}
  }
  fclose (cam_file);
  return(0);
};

int GetChunkPositionInScaffold(ChunkInstanceT *chunk, int *left_end, int *right_end, 
							   int *chunkScaffoldOrientation)
{  
  NodeCGW_T *containingContig = GetGraphNode(ScaffoldGraph->ContigGraph, chunk->info.CI.contigID);

  if (containingContig->offsetAEnd.mean < containingContig->offsetBEnd.mean)  // contig is direct in scaffold
  {
	if (chunk->offsetAEnd.mean <= chunk->offsetBEnd.mean) 
	{
	  *left_end = containingContig->offsetAEnd.mean + chunk->offsetAEnd.mean;
	  *right_end = containingContig->offsetAEnd.mean + chunk->offsetBEnd.mean;
	  *chunkScaffoldOrientation = 0;
	} 
	else 
	{
	  *left_end = containingContig->offsetAEnd.mean + chunk->offsetBEnd.mean;
	  *right_end = containingContig->offsetAEnd.mean + chunk->offsetAEnd.mean;
	  *chunkScaffoldOrientation = 1;
	}
  }
  else  // contig is reversed in scaffold
  {
	if (chunk->offsetAEnd.mean <= chunk->offsetBEnd.mean) 
	{
	  *left_end = containingContig->offsetAEnd.mean - chunk->offsetBEnd.mean;
	  *right_end = containingContig->offsetAEnd.mean - chunk->offsetAEnd.mean;
	  *chunkScaffoldOrientation = 1;
	} 
	else 
	{
	  *left_end = containingContig->offsetAEnd.mean - chunk->offsetAEnd.mean;
	  *right_end = containingContig->offsetAEnd.mean - chunk->offsetBEnd.mean;
	  *chunkScaffoldOrientation = 0;
	}
  }
  return 0;
}

int GetContigPositionInScaffold(ChunkInstanceT *contig, int *left_end, int *right_end, 
								int *contigScaffoldOrientation)
{  
  if (contig->offsetAEnd.mean <= contig->offsetBEnd.mean) 
  {
	*left_end = contig->offsetAEnd.mean;
	*right_end = contig->offsetBEnd.mean;
	*contigScaffoldOrientation = 0;
  } 
  else 
  {
	*left_end = contig->offsetBEnd.mean;
	*right_end = contig->offsetAEnd.mean;
	*contigScaffoldOrientation = 1;
  }
  return 0;
}

int GetChunkPositionInContig(ChunkInstanceT *chunk, int *left_end, int *right_end, 
							 int *chunkScaffoldOrientation)
{  
  if (chunk->offsetAEnd.mean <= chunk->offsetBEnd.mean) 
  {
	*left_end = chunk->offsetAEnd.mean;
	*right_end = chunk->offsetBEnd.mean;
	*chunkScaffoldOrientation = 0;
  } 
  else 
  {
	*left_end = chunk->offsetBEnd.mean;
	*right_end = chunk->offsetAEnd.mean;
	*chunkScaffoldOrientation = 1;
  }
  return 0;
}


int GetFragmentPositionInContigFromChunk(CIFragT *frag, int *left_end, int *right_end, 
										 int *fragmentContigOrientation, 
										 int chunkContigOrientation, int chunkLeftEnd, int chunkRightEnd)
{
  fprintf(stderr, "for frag %d, chunkLeftEnd: %d, chunkRightEnd: %d\n", frag->iid, chunkLeftEnd, chunkRightEnd);
  fprintf(stderr, "for frag %d, frag->offset5p.mean: %f, frag->offset3p.mean : %f\n", 
		  frag->iid, frag->offset5p.mean, frag->offset3p.mean);
  
  if (chunkContigOrientation == 0)  // chunk is direct in contig
  {
	if (frag->offset5p.mean < frag->offset3p.mean)  // frag is direct in chunk
	{
	  *left_end = chunkLeftEnd + frag->offset5p.mean;
	  *right_end = chunkLeftEnd + frag->offset3p.mean;
	  *fragmentContigOrientation = 0;
	}
	else  // frag is reversed in chunk
	{
	  *left_end = chunkLeftEnd + frag->offset3p.mean;
	  *right_end = chunkLeftEnd + frag->offset5p.mean;
	  *fragmentContigOrientation = 1;
	}
  }
  else   // chunk is reversed in contig
  {
	if (frag->offset5p.mean < frag->offset3p.mean)  // frag is direct in chunk
	{
	  *left_end = chunkRightEnd - frag->offset3p.mean;
	  *right_end = chunkRightEnd - frag->offset5p.mean;
	  *fragmentContigOrientation = 1;
	}
	else  // frag is reversed in chunk
	{
	  *left_end = chunkRightEnd - frag->offset5p.mean;
	  *right_end = chunkRightEnd - frag->offset3p.mean;
	  *fragmentContigOrientation = 0;
	}
  }
  return 0;
}

int GetChunkGappedPositionInContig(MultiAlignT *maContig, ChunkInstanceT *chunk, int *left_end, int *right_end, 
								   int *chunkScaffoldOrientation)
{  
  int i;
  
  // cycle through unitigs
  for (i = 0; i < GetNumIntUnitigPoss(maContig->u_list); i++)
  {
	IntUnitigPos *upContig = GetIntUnitigPos(maContig->u_list, i);

	if (upContig->ident != chunk->id)
	  continue;

	if ( upContig->position.bgn <= upContig->position.end) 
	{
	  *left_end = upContig->position.bgn;
	  *right_end = upContig->position.end;
	  *chunkScaffoldOrientation = 0;
	} 
	else 
	{
	  *left_end = upContig->position.end;
	  *right_end = upContig->position.bgn;
	  *chunkScaffoldOrientation = 1;
	}
	return 0;
  }
}


int GetFragmentGappedPositionInContigFromChunk(MultiAlignT *maChunk, CIFragT *frag, int *left_end, int *right_end, 
											   int *fragmentContigOrientation, 
											   int chunkContigOrientation, int chunkLeftEnd, int chunkRightEnd)
{
  int i;
  
  // cycle through fragments
  for (i = 0; i < GetNumIntMultiPoss(maChunk->f_list); i++)
  {
	IntMultiPos *mpChunk = GetIntMultiPos(maChunk->f_list, i);
	
	if (mpChunk->ident != frag->iid)
	  continue;
	
	fprintf(stderr, "for frag %d, chunkLeftEnd: %d, chunkRightEnd: %d\n", frag->iid, chunkLeftEnd, chunkRightEnd);
	fprintf(stderr, "for frag %d, frag->offset5p.mean: %f, frag->offset3p.mean : %f\n", 
			frag->iid, frag->offset5p.mean, frag->offset3p.mean);
	
	//*left_end = mpChunk->position.bgn;
	//*right_end =  mpChunk->position.end;

	if (chunkContigOrientation == 0)  // chunk is direct in contig
	{
	  if (mpChunk->position.bgn < mpChunk->position.end)  // frag is direct in chunk
	  {
		*left_end = chunkLeftEnd + mpChunk->position.bgn;
		*right_end = chunkLeftEnd + mpChunk->position.end;
		*fragmentContigOrientation = 0;
	  }
	  else  // frag is reversed in chunk
	  {
		*left_end = chunkLeftEnd + mpChunk->position.end;
		*right_end = chunkLeftEnd + mpChunk->position.bgn;
		*fragmentContigOrientation = 1;
	  }
	}
	else   // chunk is reversed in contig
	{
	  if (mpChunk->position.bgn < mpChunk->position.end)  // frag is direct in chunk
	  {
		*left_end = chunkRightEnd - mpChunk->position.end;
		*right_end = chunkRightEnd - mpChunk->position.bgn;
		*fragmentContigOrientation = 1;
	  }
	  else  // frag is reversed in chunk
	  {
		*left_end = chunkRightEnd - mpChunk->position.bgn;
		*right_end = chunkRightEnd - mpChunk->position.end;
		*fragmentContigOrientation = 0;
	  }
	}
	return 0;
  }  
}

// might need another version of this that goes from scratch, assuming frag is in a contig in a scaffold
// just use a flag and GetChunkPositionInScaffold

int GetFragmentPositionInScaffoldFromChunk(CIFragT *frag, int *left_end, int *right_end, 
										   int *fragmentScaffoldOrientation, 
										   int chunkScaffoldOrientation, int chunkLeftEnd, int chunkRightEnd)
{

  fprintf(stderr, "for frag %d, chunkLeftEnd: %d, chunkRightEnd: %d\n", frag->iid, chunkLeftEnd, chunkRightEnd);
  fprintf(stderr, "for frag %d, frag->offset5p.mean: %f, frag->offset3p.mean : %f\n", 
		  frag->iid, frag->offset5p.mean, frag->offset3p.mean);
  
  if (chunkScaffoldOrientation == 0)  // chunk is direct in scaffold
  {
	if (frag->offset5p.mean < frag->offset3p.mean)  // frag is direct in chunk
	{
	  *left_end = chunkLeftEnd + frag->offset5p.mean;
	  *right_end = chunkLeftEnd + frag->offset3p.mean;
	  *fragmentScaffoldOrientation = 0;
	}
	else  // frag is reversed in chunk
	{
	  *left_end = chunkLeftEnd + frag->offset3p.mean;
	  *right_end = chunkLeftEnd + frag->offset5p.mean;
	  *fragmentScaffoldOrientation = 1;
	}
  }
  else   // chunk is reversed in scaffold
  {
	if (frag->offset5p.mean < frag->offset3p.mean)  // frag is direct in chunk
	{
	  *left_end = chunkRightEnd - frag->offset3p.mean;
	  *right_end = chunkRightEnd - frag->offset5p.mean;
	  *fragmentScaffoldOrientation = 1;
	}
	else  // frag is reversed in chunk
	{
	  *left_end = chunkRightEnd - frag->offset5p.mean;
	  *right_end = chunkRightEnd - frag->offset3p.mean;
	  *fragmentScaffoldOrientation = 0;
	}
  }
  return 0;
}

int myPrint_Overlap_AS(FILE *file, InternalFragMesg *a,
					   InternalFragMesg *b, OverlapMesg *align)
{
//  if (a->iaccession == align->bifrag)
  //  {
  //  InternalFragMesg *c;
  //  c = a;
  //  a = b;
  //  b = c;
  //}
  
  fprintf(file,"\nOVERLAP BETWEEN");
  fprintf(file," A = (%ld,%d)",a->eaccession,a->iaccession);
  fprintf(file," and");
  fprintf(file," B = (%ld,%d)",b->eaccession,b->iaccession);
  fprintf(file,"\n\n");

  switch (align->orientation)
  { case AS_NORMAL:
	if (align->bhg <= 0)
	{ fprintf(file,"  A -----+------+----> %-4d\n",align->bhg);
	fprintf(file,"    %4d -------> B\n",align->ahg);
	}
	else
	{ fprintf(file,"  A -----+------> %-4d\n",align->bhg);
	fprintf(file,"    %4d -------+----> B\n",align->ahg);
	}
	break;
    case AS_INNIE:
      if (align->bhg <= 0)
	  { fprintf(file,"  A -----+------+----> %-4d\n",align->bhg);
	  fprintf(file,"    %4d <------- B\n",align->ahg);
	  }
      else
	  { fprintf(file,"  A -----+------> %-4d\n",align->bhg);
	  fprintf(file,"    %4d <------+----- B\n",align->ahg);
	  }
      break;
    case AS_OUTTIE:
      if (align->bhg <= 0)
	  { fprintf(file,"  A <----+------+----- %-4d\n",align->bhg);
	  fprintf(file,"    %4d -------> B\n",align->ahg);
	  }
      else
	  { fprintf(file,"  A <----+------- %-4d\n",align->bhg);
	  fprintf(file,"    %4d -------+----> B\n",align->ahg);
	  }
      break;
  }
  fprintf(file,"\n");

  if (align->delta != NULL)
  { int *trace;

  fprintf(stderr, "align->delta != NULL in myPrint_Overlap_AS\n");
	
  trace = Unpack_Alignment_AS(align);

#ifdef DEBUG
  { int i;

  fprintf(file,"\nUncompressed trace:\n");
  for (i = 0; trace[i] != 0; i++)
	fprintf(file,"  %3d\n",trace[i]);
  }
#endif

  if (align->orientation == AS_INNIE)
	Complement_Fragment_AS(b);
  else if (align->orientation == AS_OUTTIE)
	Complement_Fragment_AS(a);

  myPrintAlign(file,align->ahg,align->bhg,a->sequence,b->sequence,trace);

  if (align->orientation == AS_INNIE)
	Complement_Fragment_AS(b);
  else if (align->orientation == AS_OUTTIE)
	Complement_Fragment_AS(a);
  }
} 

/*** OVERLAP PRINT ROUTINE ***/

#define PRINT_WIDTH  50   /* Width of each line of a printed alignment */
#define BAND_WIDTH    3   /* Width of band about original for realignment */

/* Print an alignment to file between a and b given in trace (unpacked).
   Prefix gives the length of the initial prefix of a that is unaligned.  */

static int myPrintAlign(FILE *file, int prefix, int suffix,
						char *a, char *b, int *trace)
{ int i, j, o;
 static char Abuf[PRINT_WIDTH+1], Bbuf[PRINT_WIDTH+1];
 static int  Firstime = 1;

 if (Firstime)
 { Firstime = 0;
 Abuf[PRINT_WIDTH] = Bbuf[PRINT_WIDTH] = '\0';
 }
 /* buffer/output next column */
#define COLUMN(x,y)			\
 { if (o >= PRINT_WIDTH)			\
   { fprintf(file,"\n\t%s\n",Abuf);	\
									  fprintf(file,"\t%s\n",Bbuf);	\
																	  o = 0;				\
																							  }					\
																												  Abuf[o] = (x);			\
																																			  Bbuf[o] = (y);			\
																																										  o += 1;				\
																																																  }

 a -= 1;
 b -= 1;
 o  = 0;
 i = j = 1;

 if (prefix > AS_READ_MAX_LEN)
 { i = prefix-24;
 prefix = 25;
 }

 while (prefix-- > 0)     /* Output unaligned prefix */
   COLUMN(a[i++],' ')

   { int p, c;      /* Output columns of alignment til reach trace end */

   p = 0;
   while ((c = trace[p++]) != 0)
	 if (c < 0)
	 { c = -c;
	 while (i != c)
	   COLUMN(a[i++],b[j++])
		 COLUMN('-',b[j++])
		 }
	 else
	 { while (j != c)
	   COLUMN(a[i++],b[j++])
		 COLUMN(a[i++],'-')
		 }
   }

 if (suffix < 0) suffix = -suffix;
 if (suffix > AS_READ_MAX_LEN)
   suffix = 25;

 { int x, y, s;     /* Output remaining column including unaligned suffix */

 s = 0;
 y = 1;
 while ((x = a[i++]) != 0)
 { if ((y = b[j++]) != 0)
   COLUMN(x,y)
	 else
	 { do 
	 { COLUMN(x,' ')
		 s += 1;
	 }
	   while ((x = a[i++]) != 0 && s < suffix);
	 break;
	 }
 }
 if (y)
   while ((y = b[j++]) != 0 && s < suffix)
   { COLUMN(' ',y)
	   s += 1;
   }
 }

 fprintf(file,"\n\t%.*s\n",o,Abuf);   /* Print remainder of buffered col.s */
 fprintf(file,"\t%.*s\n",o,Bbuf);
}

int getFragmentSequenceAndQuality(int iid, char *fragSequence, char *fragQuality)
{
  ScreenedFragMesg fmesg;  
  MesgWriter   writer;
  static char fseq[1025];
  static char fqual[1025];
  ReadStructp fsread = new_ReadStruct();

  fmesg.sequence = fseq;
  fmesg.quality = fqual;

  // fprintf(stderr, "calling getFragStore with iid: %d, fragStore: %p\n", iid, ScaffoldGraph->fragStore);
  
  getFragStore( ScaffoldGraph->fragStore, (IntFragment_ID) iid, FRAG_S_ALL, fsread);
  getSequence_ReadStruct( fsread, fmesg.sequence, fmesg.quality, CNS_MAX_SEQUENCE_LENGTH);
  // fprintf(stderr, "%s\n", fmesg.sequence);
  // fprintf(stderr, "%s\n", fmesg.quality);
  
  //getLocID_ReadStruct( fsread, &fmesg.elocale);
  //getLocalePos_ReadStruct( fsread, (uint32 *)&fmesg.locale_pos.bgn, (uint32 *)&fmesg.locale_pos.end);

  getClearRegion_ReadStruct( fsread, (uint32 *)&fmesg.clear_rng.bgn, (uint32 *)&fmesg.clear_rng.end);

  delete_ReadStruct(fsread);

  // fprintf(stderr, "fmesg.clear_rng.bgn: %d, fmesg.clear_rng.end: %d\n", fmesg.clear_rng.bgn, fmesg.clear_rng.end);

  if (fmesg.clear_rng.end - fmesg.clear_rng.bgn > CNS_MAX_SEQUENCE_LENGTH)
	assert(0);

  strncpy(fragSequence, &fmesg.sequence[fmesg.clear_rng.bgn], fmesg.clear_rng.end - fmesg.clear_rng.bgn);
  fragSequence[fmesg.clear_rng.end] = '\0';

  strncpy(fragQuality, &fmesg.quality[fmesg.clear_rng.bgn], fmesg.clear_rng.end - fmesg.clear_rng.bgn);
  fragQuality[fmesg.clear_rng.end] = '\0';
}

// types go from 1 to 10
#define NUM_REPEAT_TYPES 10

int getFragmentScreenInfo(int iid, int *repeatIDCount)
{
  ScreenedFragMesg 
	fmesg;
  ReadStructp 
	fsread = new_ReadStruct();
  int 
	nm,
	numMatches;
  
  // fprintf(stderr, "calling getFragStore with iid: %d, fragStore: %p\n", iid, ScaffoldGraph->fragStore);
  
  getFragStore( ScaffoldGraph->fragStore, (IntFragment_ID) iid, FRAG_S_ALL, fsread);

  fmesg.screened = NULL;
  numMatches = getScreenMatches_ReadStruct(fsread, fmesg.screened, 0);
  if (numMatches > 0) 
  {
	fmesg.screened = (IntScreenMatch *) malloc(numMatches * sizeof(IntScreenMatch));
	getScreenMatches_ReadStruct(fsread, fmesg.screened, numMatches);

	repeatIDCount[0]++;  // total fragments analyzed
	for (nm = 0; nm < numMatches; nm++)
	{
	  repeatIDCount[fmesg.screened[nm].repeat_id]++;
	}
  }
  
  //for (nm = 0; nm < numMatches; nm++) 
  //fprintf(stderr, " (%d [%d, %d])", fmesg.screened[nm].iwhat, fmesg.screened[nm].where.bgn,
  //		fmesg.screened[nm].where.end);

  delete_ReadStruct(fsread);
  if (fmesg.screened) 
	free(fmesg.screened);
}

int getChunkInstanceID(ChunkInstanceT *chunk, int index)
{
  if (chunk->info.CI.numInstances == 0)  // chunk is not a surrogate
  {
	if (index == 0)  // just return chunk's id
	  return(chunk->id);
	else
	  return(-1);
  }
  else  // chunk is a surrogate
  {
	if (chunk->info.CI.numInstances == 1  && index == 0)
	  return( chunk->info.CI.instances.in_line.instance1 );
	else if (chunk->info.CI.numInstances == 2 && (index == 0 || index == 1))
	{
	  if (index == 0)
		return( chunk->info.CI.instances.in_line.instance1);
	  else if (index == 1)
		return( chunk->info.CI.instances.in_line.instance2);
	}
	else if (index < chunk->info.CI.numInstances)
	  return( * (int32 *) Getint32(chunk->info.CI.instances.va, index));
	else
	  return(-1);
  }
}

int GetContigPlacementLimits(CIEdgeT *e, int overlappingContigScaffoldOrientation, 
							 int overlappingContigLeftEnd, int overlappingContigRightEnd,
							 ChunkInstanceT *contig, int *contigLeftLimit, int *contigRightLimit, int *contigOrientation,
							 int tight)
{
  if (overlappingContigScaffoldOrientation == 0) // direct in scaffold
  {
	if ((e->orient == AB_AB && contig->id == e->idA) || (e->orient == BA_BA && contig->id == e->idB))
	{
	  *contigLeftLimit = overlappingContigRightEnd + e->distance.mean;
	  *contigRightLimit = *contigLeftLimit + contig->bpLength.mean;
	  *contigOrientation = 0;
	}
	else if ((e->orient == AB_BA && contig->id == e->idA) || (e->orient == AB_BA && contig->id == e->idB))
	{
	  *contigLeftLimit = overlappingContigRightEnd + e->distance.mean;
	  *contigRightLimit = *contigLeftLimit + contig->bpLength.mean;
	  *contigOrientation = 1;
	}
	else if ((e->orient == BA_AB && contig->id == e->idA) || (e->orient == BA_AB && contig->id == e->idB))
	{
	  *contigRightLimit = overlappingContigLeftEnd - e->distance.mean;
	  *contigLeftLimit = *contigRightLimit - contig->bpLength.mean;
	  *contigOrientation = 1;
	}
	else if ((e->orient == BA_BA && contig->id == e->idA) || (e->orient == AB_AB && contig->id == e->idB))
	{
	  *contigRightLimit = overlappingContigLeftEnd - e->distance.mean;
	  *contigLeftLimit = *contigRightLimit - contig->bpLength.mean;
	  *contigOrientation = 0;
	}
  }
  else if (overlappingContigScaffoldOrientation == 1) // reversed in scaffold
  {
	if ((e->orient == AB_AB && contig->id == e->idA) || (e->orient == BA_BA && contig->id == e->idB))
	{
	  *contigRightLimit = overlappingContigLeftEnd - e->distance.mean;
	  *contigLeftLimit = *contigRightLimit - contig->bpLength.mean;
	  *contigOrientation = 1;
	}
	else if ((e->orient == AB_BA && contig->id == e->idA) || (e->orient == AB_BA && contig->id == e->idB))
	{
	  *contigRightLimit = overlappingContigLeftEnd - e->distance.mean;
	  *contigLeftLimit = *contigRightLimit - contig->bpLength.mean;
	  *contigOrientation = 0;
	}
	else if ((e->orient == BA_AB && contig->id == e->idA) || (e->orient == BA_AB && contig->id == e->idB))
	{
	  *contigLeftLimit = overlappingContigRightEnd + e->distance.mean;
	  *contigRightLimit = *contigLeftLimit + contig->bpLength.mean;
	  *contigOrientation = 0;
	}
	else if ((e->orient == BA_BA && contig->id == e->idA) || (e->orient == AB_AB && contig->id == e->idB))
	{
	  *contigLeftLimit = overlappingContigRightEnd + e->distance.mean;
	  *contigRightLimit = *contigLeftLimit + contig->bpLength.mean;
	  *contigOrientation = 1;
	}
  }
  if (tight == FALSE)
  {
	*contigLeftLimit -= 3.0 * sqrt(e->distance.variance);
	*contigRightLimit += 3.0 * sqrt(e->distance.variance);
  }
}

// look at all contigs w/scaffoldID == -1
// place them based on quality values and existing edges of the contig graph
int PlaceUnscaffoldedContigs(float (* quality)(CIEdgeT *, int32), float qualityThresh)
{
  ChunkInstanceT 
	* contig;
  GraphNodeIterator
	ContigGraphIterator;
  int
	unscaffoldedContigs = 0,
	contigsWithValidEdges = 0,
	previousContigId = -1,
	containedContigsWithValidEdges = 0,
	previousOverlappedContigId = 0,
	totalEdges = 0,
	validEdges = 0,
	validEdgesToPlacedContigs = 0,
	placedContigs = 0,
	overlapsAttempted = 0;

    fprintf(stderr,
		  "------------------------------: fragment placing (unscaffolded contigs) starting\n");

  InitGraphNodeIterator(&ContigGraphIterator, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while (contig = NextGraphNodeIterator(&ContigGraphIterator))
  {
	CIEdgeT
	  * e;
	GraphEdgeIterator
	  edges;
	int32
	  overlapping_contig_cid,
	  overlapping_contig_left_end, overlapping_contig_right_end,
	  overlappingContigScaffoldOrientation,
	  contigLeftLimit, contigRightLimit, contigOrientation;
	ChunkInstanceT
	  * overlapping_contig,
	  * bestContig;
	double bestQuality = 1.0;
	
    //
    // scan all the nodes (contigs)
    //
	
	// if numInstances > 0 then it has a surrogate, and we don't want it
	// but does numInstances apply only to chunks or to contigs as well? Ask Saul.
	if (!(contig->scaffoldID == -1 && contig->info.CI.numInstances == 0))
	  continue;
	
	unscaffoldedContigs++;
	fprintf( stderr, "\nExamining contig: %d\n", contig->id);
	
	// loop over all the edges from this contig
	InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, contig->id, 
						  ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
	while (e = NextGraphEdgeIterator(&edges)) 
	{
	  MultiAlignT *maOverlappingContig, *maContig;
	  static VA_TYPE(char) *overlappingContigConsensusSequence = NULL;
	  static VA_TYPE(char) *overlappingContigConsensusQuality = NULL;
	  static VA_TYPE(char) *contigConsensusSequence = NULL;
	  static VA_TYPE(char) *contigConsensusQuality = NULL;
	  
	  assert(e != NULL);
	  if (contig->id == e->idA)
		overlapping_contig_cid = e->idB;
	  else
		overlapping_contig_cid = e->idA;
	  overlapping_contig = GetGraphNode(ScaffoldGraph->ContigGraph, overlapping_contig_cid);
	  assert(overlapping_contig != NULL);

	  totalEdges++;

	  if (!Is_Not_Bogus_Not_Guide(e))  // don't play with funky edges
		continue;
	  validEdges++;

	  if (overlapping_contig->scaffoldID != -1)  // we only want to place off of contigs that are themselves placed
	  {
		validEdgesToPlacedContigs++;
		if (contig->id != previousContigId)
		{
		  contigsWithValidEdges++;
		  previousContigId = contig->id;		  
		}

		// find out where the contig with the edge to it is placed
		GetContigPositionInScaffold(overlapping_contig, &overlapping_contig_left_end, &overlapping_contig_right_end, 
									&overlappingContigScaffoldOrientation);
		
		// now determine the limits of where our contig can be
		GetContigPlacementLimits(e, overlappingContigScaffoldOrientation, overlapping_contig_left_end, overlapping_contig_right_end,
								 contig, &contigLeftLimit, &contigRightLimit, &contigOrientation, FALSE);

		// experiment with Karin's stuff
		if (0)
		{
		  MultiAlignT  *newMultiAlign;
		  static VA_TYPE(IntElementPos) *ContigPositions = NULL;
		  IntMultiPos contigPos;
		  
		  if(ContigPositions == NULL)
			ContigPositions = CreateVA_IntElementPos(2);

		  // contigPos.type = 'C';
		  contigPos.type &= AS_CONTIG;
		  contigPos.ident = contig->id;
		  contigPos.position.bgn = overlapping_contig_left_end;
		  contigPos.position.end = overlapping_contig_right_end;
		  AppendIntElementPos(ContigPositions, &contigPos);

		  contigPos.type &= AS_CONTIG;
		  contigPos.ident = overlapping_contig->id;
		  contigPos.position.bgn = contigLeftLimit;
		  contigPos.position.end = contigRightLimit;
		  AppendIntElementPos(ContigPositions, &contigPos);
		  
		  newMultiAlign = MergeMultiAligns(ScaffoldGraph->CIGraph->maStore, 
										   ScaffoldGraph->ContigGraph->maStore, 
										   ScaffoldGraph->fragStore,
										   ContigPositions, FALSE);
		  
		  if(!newMultiAlign)
			fprintf (stderr, "MergeMultiAligns for contig %d and overlapping_contig %d failed!\n", 
					 contig->id, overlapping_contig->id);
		  else
			fprintf (stderr, "MergeMultiAligns for contig %d and overlapping_contig %d succeeded!\n", 
					 contig->id, overlapping_contig->id);
#if 0
			if(!newMultiAlign)
		  {
			fprintf(stderr,"* Calling MergeMultiAligns with %ld elements*\n",
					GetNumIntElementPoss(ContigPositions));
			{
			  int i;
			  
			  fprintf(stderr,"*** Contigging ContigPositions\n");
			  for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++)
			  {
				IntElementPos *pos = GetIntElementPos(ContigPositions,i);
				ContigT *contig = GetGraphNode(ScaffoldGraph->ContigGraph, pos->ident);
				int32 originalID = GetOriginalContigID(contig->id);
				
				fprintf(stderr,"* Contig %d (%d)   bgn:%d end:%d [%d,%d]\n",
						pos->ident, originalID, pos->position.bgn, pos->position.end, contig->aEndCoord, contig->bEndCoord);
			  }
			}
			DumpCIScaffold(stderr, ScaffoldGraph, scaffold, FALSE);
			fprintf(stderr,"* MergeMultiAligns failed....bye\n");
			fflush(NULL);
			assert(0);
		  }
#endif
		} // end of Karin's stuff

		
		// then overlap them

		// see if contig is contained in placed contig - is this limitation necessary?
		if (overlapping_contig_left_end < contigLeftLimit && contigRightLimit < overlapping_contig_right_end)
		{
		  if (contig->id != previousOverlappedContigId)
		  {
			containedContigsWithValidEdges++;
			previousOverlappedContigId = contig->id;			
		  }
		  overlapsAttempted++;

		  maOverlappingContig = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, overlapping_contig->id);
		  assert (maOverlappingContig);

		  maContig = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id);
		  assert (maContig);
		  
		  fprintf( stderr, "trying to place: %d between (%d, %d) edge (%f, %c) \nplaced contig: %d, position (%d, %d), scaffold: %d, previous: %d, next: %d\n",
				   contig->id, contigLeftLimit, contigRightLimit,
				   e->distance.mean, e->orient,
				   overlapping_contig_cid, overlapping_contig_left_end, overlapping_contig_right_end,
				   overlapping_contig->scaffoldID, overlapping_contig->AEndNext, overlapping_contig->BEndNext);
		  
		  if(overlappingContigConsensusSequence == NULL)
		  {
			overlappingContigConsensusSequence = CreateVA_char(2048);
			overlappingContigConsensusQuality = CreateVA_char(2048);
			contigConsensusSequence = CreateVA_char(2048);
			contigConsensusQuality = CreateVA_char(2048);
		  }
		  
		  // grab the corresponding consensus plus slop
		  GetMultiAlignUngappedConsensus(maOverlappingContig, overlappingContigConsensusSequence, overlappingContigConsensusQuality);

		  // grab the corresponding consensus plus slop
		  GetMultiAlignUngappedConsensus(maContig, contigConsensusSequence, contigConsensusQuality);		  
	  
#if 1
		  // overlap 'em
		  {
			int A_hang_estimate,A_hang_lower,A_hang_upper;
			InternalFragMesg AFR, BFR;
			OverlapMesg     *O;
			int where, opposite = FALSE;
			FILE* log = NULL;
			int length;
			
			AFR.sequence   = Getchar(overlappingContigConsensusSequence, 0);
			AFR.quality    = Getchar(overlappingContigConsensusQuality, 0);
			AFR.eaccession = 0;
			AFR.iaccession = 0;
			BFR.sequence   = Getchar(contigConsensusSequence, 0);
			BFR.quality    = Getchar(contigConsensusQuality, 0);
			BFR.eaccession = 0;			
			BFR.iaccession = 1;

			if (e->orient == 'I' || e->orient == 'O')
			  opposite = TRUE;
			
			// add some slop
			// contigLeftLimit -= 500;
			// contigRightLimit += 500;
			contigLeftLimit = 0;
			contigRightLimit = strlen(AFR.sequence);
			
			fprintf(stderr,"* Calling DP_Compare with %s sequences lengths %ld,%ld bgn,end =[%d,%d] %c %c\n",
					(!strcmp(AFR.sequence, BFR.sequence)?"EQUAL":"DIFFERENT"),
					strlen(AFR.sequence), strlen(BFR.sequence),
					contigLeftLimit, contigRightLimit, *AFR.sequence, *BFR.sequence);

			//CGW_DP_ERATE=0.1 changed parameter to canOlap->errRate
			// a new field in ContigOverlapCheckT
			O = DP_Compare_AS(&AFR, &BFR, 
							  contigLeftLimit, contigRightLimit, opposite,
							  CGW_DP_ERATE, CGW_DP_THRESH, CGW_DP_MINLEN,
							  AS_FIND_ALIGN, //_NO_TRACE, // Slower, computes an alignment, throw away delta encoding
							  &where);
			if (O != NULL)
			{
			  // myPrint_Overlap_AS(stderr, &AFR, &BFR, O);
			  if(log == NULL)
				log = stderr;
			  
			  compute_bayesian_quality(&AFR, &BFR, O, 0, &length, log);
			  fprintf(stderr, "dorf, Quality between contig %4d and overlapping contig %4d WITHOUT QV Realigning = %lf\n", 
					  contig->id, overlapping_contig->id, O->quality);
			  if (O->quality < bestQuality)
			  {
				bestQuality = O->quality;
				bestContig = overlapping_contig;
			  }
			}		
		  }
#endif
		}
	  }
	}
	// pick the maximum, place it there
	if (bestQuality < qualityThresh)
	{
	  fprintf(stderr, "placing contig %d based on a quality value of %f with contig %d\n", contig->id, bestQuality, bestContig->id);
	  placedContigs++;





	  // wait a second!!! we need to place this guy based on the coords calculated by DP_Align, not where we started
	  // what about this?





	  //InsertCIInScaffold(ScaffoldGraphT *sgraph, int32 bestContig, int32 sid, 
	  //				   LengthT aEndOffset, LengthT bEndOffset, int AEndToBend, int contigNow);
	}
  }
  fprintf(stderr, "unscaffoldedContigs: %d\n", unscaffoldedContigs);
  fprintf(stderr, "contigsWithValidEdges: %d\n", contigsWithValidEdges);
  fprintf(stderr, "totalEdges: %d\n", totalEdges);
  fprintf(stderr, "validEdges: %d\n", validEdges);
  fprintf(stderr, "validEdgesToPlacedContigs: %d\n", validEdgesToPlacedContigs);
  fprintf(stderr, "overlapsAttempted: %d\n", overlapsAttempted);  
  fprintf(stderr, "containedContigsWithValidEdges: %d\n", containedContigsWithValidEdges);
  fprintf(stderr, "placedContigs: %d\n", placedContigs);
  if (unscaffoldedContigs > 0)
	fprintf(stderr, "placedContigs/unscaffoldedContigs: %d/%d (%.2f %%)\n", placedContigs, unscaffoldedContigs, 
			(100.0 * placedContigs) / unscaffoldedContigs);  
  else
	fprintf(stderr, "placedContigs/unscaffoldedContigs: %d/%d (%.2f %%)\n", placedContigs, unscaffoldedContigs, 
			0.0);

  if (containedContigsWithValidEdges > 0)
	fprintf(stderr, "placedContigs/containedContigsWithValidEdges: %d/%d (%.2f %%)\n", placedContigs, containedContigsWithValidEdges, 
			(100.0 * placedContigs) / containedContigsWithValidEdges);  
  else
	fprintf(stderr, "placedContigs/containedContigsWithValidEdges: %d/%d (%.2f %%)\n", placedContigs, containedContigsWithValidEdges, 
			0.0);  

  fprintf(stderr,
		  "------------------------------: fragment placing (unscaffolded contigs) ending\n");
}


int CountFragments(float (* quality)(CIEdgeT *, int32), float qualityThresh)
{
  ChunkInstanceT 
	* chunk;
  GraphNodeIterator
	CIGraphIterator;
  int
	i,
	contigsWithValidEdges = 0,
	previousContigId = -1,
	containedContigsWithValidEdges = 0,
	previousOverlappedContigId = 0,
	totalEdges = 0,
	validEdges = 0,
	validEdgesToPlacedContigs = 0,
	placedContigs = 0,
	fragmentsInScaffoldedChunks = 0,
	* fragmentsInUnscaffoldedChunks,
	fragmentsInUnscaffoldedChunksWithValidEdges = 0,
	fragmentsInParentChunks = 0,
	scaffoldedChunks = 0,
	* unscaffoldedChunks,
	* timesContained,
	unscaffoldedChunksWithValidEdges = 0,
	parentChunks = 0,
	fragmentCount = 0,
	overlapsAttempted = 0,
	numChunks = 0;
  int 
	chunkDepth,
	* chunkEdgeDepth,
	chunkLabelled,
	fragmentsInChunkDepth1Contains = 0,
	chunkDepth1Contains = 0;
  fprintf(stderr,
		  "------------------------------: fragment counting starting\n");

#if 0
  // is there a data struct that already holds the number of chunks?
  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while (chunk = NextGraphNodeIterator(&CIGraphIterator))
	numChunks++;
  fprintf(stderr, "\n");
  fprintf(stderr, "numChunks = %d\n", numChunks);
#endif

  numChunks = GetNumGraphNodes(ScaffoldGraph->CIGraph);
  fprintf(stderr, "\n");
  fprintf(stderr, "numChunks = %d\n", numChunks);

  chunkEdgeDepth = (int *) malloc(numChunks * sizeof(int));
  assert (chunkEdgeDepth);
  
  fragmentsInUnscaffoldedChunks = (int *) malloc(numChunks * sizeof(int));
  assert (fragmentsInUnscaffoldedChunks);

  unscaffoldedChunks = (int *) malloc(numChunks * sizeof(int));
  assert (unscaffoldedChunks);

  timesContained = (int *) malloc(numChunks * sizeof(int));
  assert (timesContained);

  // scan all the nodes (chunks)
  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while (chunk = NextGraphNodeIterator(&CIGraphIterator))
  {
	chunkEdgeDepth[chunk->id] = -2; // initialize
	timesContained[chunk->id] = 0;
	
	// count all the chunks that are scaffolded
	// those contigs that have a valid scaffold number and are not surrogates
	if (chunk->scaffoldID != -1 && chunk->type != RESOLVEDREPEATCHUNK_CGW)
	{
	  MultiAlignT *maChunk = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, chunk->id);
	  assert (maChunk);

	  // fprintf(stderr, "assigning chunk %d as a scaffolded chunk\n", chunk->id);
	  
	  chunkEdgeDepth[chunk->id] = 0; // 0 => chunk in scaffold
	  scaffoldedChunks++;
	  fragmentsInScaffoldedChunks += GetNumIntMultiPoss(maChunk->f_list);
	}
	
	// those contigs that have a valid scaffold number and are surrogates
	if (chunk->scaffoldID != -1 && chunk->type == RESOLVEDREPEATCHUNK_CGW)
	{
	  MultiAlignT *maChunk = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, chunk->id);
	  assert (maChunk);
	  
	  // fprintf(stderr, "assigning chunk %d as a surrogate chunk\n", chunk->id);

	  chunkEdgeDepth[chunk->id] = 0; // 0 => chunk in scaffold
	  scaffoldedChunks++;
	  fragmentsInScaffoldedChunks += GetNumIntMultiPoss(maChunk->f_list);
	}
	
	// those contigs that are parents for surrogates
	if (chunk->scaffoldID == -1 && chunk->info.CI.numInstances >= 1)
	{
	  MultiAlignT *maChunk = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, chunk->id);
	  assert (maChunk);
	  
	  // fprintf(stderr, "assigning chunk %d as a parent chunk\n", chunk->id);

	  chunkEdgeDepth[chunk->id] = -1; // -1 => an unscaffolded chunk that has surrogates 
	  parentChunks++;
	  fragmentsInParentChunks += GetNumIntMultiPoss(maChunk->f_list);
	}
  }
  
  chunkDepth = 0;
  chunkLabelled = 1;
  while (chunkLabelled)
  {
	chunkLabelled = 0;
	chunkDepth++;
	unscaffoldedChunks[chunkDepth] = 0;
	fragmentsInUnscaffoldedChunks[chunkDepth] = 0;
	
	// fprintf( stderr, "\nAt chunkDepth: %d\n", chunkDepth);

	// scan all the nodes (chunks)
	for (i = 0; i < numChunks; i++)
	{
	  CIEdgeT
		* e;
	  GraphEdgeIterator
		edges;
	  int32
		paired_chunk_id;
	  ChunkInstanceT
		* paired_chunk;
	  
	  chunk = GetGraphNode(ScaffoldGraph->CIGraph, i);
	  
	  // deal only with those chunks that are not in a scaffold and are not parents and not already labelled
	  if (chunk->scaffoldID == -1 && chunk->info.CI.numInstances == 0 && chunkEdgeDepth[i] == -2)
	  {		
		// fprintf( stderr, "\nExamining chunk: %d\n", chunk->id);
		
		// loop over all the edges from this chunk
		InitGraphEdgeIterator(ScaffoldGraph->CIGraph, chunk->id, 
							  ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
		while (e = NextGraphEdgeIterator(&edges)) 
		{
		  assert(e != NULL);
		  if (chunk->id == e->idA)
			paired_chunk_id = e->idB;
		  else
			paired_chunk_id = e->idA;
		  paired_chunk = GetGraphNode(ScaffoldGraph->CIGraph, paired_chunk_id);
		  assert(paired_chunk != NULL);
		  
		  // fprintf( stderr, "\nExamining edge from chunk %d to chunk %d\n", chunk->id, paired_chunk->id);
		  
		  if (!Is_Not_Bogus_Not_Guide(e))  // don't play with funky edges
			continue;
		  
          // we only want to link off of chunks that are themselves assigned previously (including parents)
		  if (chunkEdgeDepth[paired_chunk->id] < chunkDepth && chunkEdgeDepth[paired_chunk->id] >= -1)  
		  {
			// fprintf(stderr, "chunk->info.CI.numFragments = %d\n", chunk->info.CI.numFragments);

			// fprintf(stderr, "assigning chunk %d as a linked chunk\n", chunk->id);
			
			if (chunkEdgeDepth[i] == -2)  // don't add in frags for this chunk more than once on different edges
			{
			  fragmentsInUnscaffoldedChunks[chunkDepth] += chunk->info.CI.numFragments;
			  unscaffoldedChunks[chunkDepth]++;
			  chunkEdgeDepth[i] = chunkDepth;
			  chunkLabelled = 1;
			}
			
			// fprintf(stderr, "linked chunk with %d fragments\n", chunk->info.CI.numFragments);

			// determine if chunks with a single link are completely contained in a contig
			if (chunkDepth == 1)
			{
			  int32
				paired_contig_left_end, paired_contig_right_end,
				pairedContigScaffoldOrientation,
				contigLeftLimit, contigRightLimit, contigOrientation;
			  ChunkInstanceT 
				*contig,
				*paired_contig;

			  contig = GetGraphNode(ScaffoldGraph->ContigGraph, chunk->info.CI.contigID);
			  paired_contig = GetGraphNode(ScaffoldGraph->ContigGraph, paired_chunk->info.CI.contigID);

			  // find out where the contig with the edge to it is placed
			  GetContigPositionInScaffold(paired_contig, &paired_contig_left_end, &paired_contig_right_end, 
										  &pairedContigScaffoldOrientation);
		
			  // now determine the limits of where our contig can be
			  GetContigPlacementLimits(e, pairedContigScaffoldOrientation, paired_contig_left_end, paired_contig_right_end,
									   contig, &contigLeftLimit, &contigRightLimit, &contigOrientation, TRUE);

			  // see if contig is contained in placed contig - is this limitation necessary?
			  if (paired_contig_left_end < contigLeftLimit && contigRightLimit < paired_contig_right_end)
			  {
				timesContained[i]++;
			  }
			  else if (0)
			  {
				fprintf( stderr, "contig: %d, %d\t paired_contig: %d %d\n", contigLeftLimit, contigRightLimit, 
						 paired_contig_left_end, paired_contig_right_end);
			  }
			}
			// done with this chunk
			// break;
		  }
		}
	  }
	}
  }
  
  fprintf(stderr, "\n");
  {
	int totalUnscaffoldedChunks = 0,
	  unlinkedChunks = 0;

	fprintf(stderr, "scaffoldedChunks: %d\n", scaffoldedChunks);
	fprintf(stderr, "parentChunks: %d\n", parentChunks);
	for (i = 1; i < chunkDepth; i++)
	{
	  totalUnscaffoldedChunks += unscaffoldedChunks[i];
	  fprintf(stderr, "unscaffoldedChunks[%d]: %d\n", i, unscaffoldedChunks[i]);
	}
	fprintf(stderr, "totalUnscaffoldedChunks = %d\n", totalUnscaffoldedChunks);
	for (i = 0; i < numChunks; i++)
	{
	  if (chunkEdgeDepth[i] == -2)
	  {
		// fprintf(stderr, "assigning chunk %d as an unlinked chunk\n", i);
		unlinkedChunks++;
	  }
	}
	fprintf(stderr, "unlinkedChunks: %d\n", unlinkedChunks);
	fprintf(stderr, "Total Chunks: %d\n", scaffoldedChunks + parentChunks + totalUnscaffoldedChunks + unlinkedChunks);
	fprintf(stderr, "\n");
  }
  
  {
	int totalFragmentsInUnscaffoldedChunks = 0,
	  fragmentsInUnlinkedChunks = 0;

	fprintf(stderr, "fragmentsInScaffoldedChunks: %d\n", fragmentsInScaffoldedChunks);
	fprintf(stderr, "fragmentsInParentChunks: %d\n", fragmentsInParentChunks);
	for (i = 1; i < chunkDepth; i++)
	{
	  totalFragmentsInUnscaffoldedChunks += fragmentsInUnscaffoldedChunks[i];
	  fprintf(stderr, "fragmentsInUnscaffoldedChunks[%d]: %d\n", i, fragmentsInUnscaffoldedChunks[i]);
	}
	fprintf(stderr, "totalFragmentsInUnscaffoldedChunks = %d\n", totalFragmentsInUnscaffoldedChunks);

	for (i = 0; i < numChunks; i++)
	{
	  int j;

	  // if (i < 10) chunkEdgeDepth[i] = -2; // testing printf

	  if (timesContained[i] > 0)
	  {
		chunkDepth1Contains++;
		chunk = GetGraphNode(ScaffoldGraph->CIGraph, i);
		fragmentsInChunkDepth1Contains += chunk->info.CI.numFragments;
	  }
	  
	  if (chunkEdgeDepth[i] == -2)
	  {
		chunk = GetGraphNode(ScaffoldGraph->CIGraph, i);
		fragmentsInUnlinkedChunks += chunk->info.CI.numFragments;

		// fprintf(stderr, "unlinked chunk %d has %d fragments\n", chunk->id, chunk->info.CI.numFragments);
		// fprintf(stderr, "unlinked chunk %d baseID: %d\n", chunk->id, chunk->info.CI.baseID);
		// fprintf(stderr, "unlinked chunk %d numInstances: %d\n\n", chunk->id, chunk->info.CI.numInstances);

		if (chunk->info.CI.numFragments >= 5)
		{
		  MultiAlignT *maChunk = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, chunk->id);
		  assert (maChunk);

		  for (j = 0; j < GetNumIntMultiPoss(maChunk->f_list); j++)
		  {
			IntMultiPos *mpChunk = GetIntMultiPos(maChunk->f_list, j);
			fprintf(stderr, "unlinked chunk %d contains fragment: %d\n", chunk->id, mpChunk->ident);
		  }
		  printEdgeInfoOnChunk(chunk->id, chunkEdgeDepth);
		  fprintf(stderr, "unlinked chunk %d consensus: %s\n\n", chunk->id, Getchar(maChunk->consensus, 0));
		}
	  }
	}
	fprintf(stderr, "fragmentsInUnlinkedChunks = %d\n", fragmentsInUnlinkedChunks);
	fprintf(stderr, "Total Fragments: %d\n", fragmentsInScaffoldedChunks + fragmentsInParentChunks + 
			totalFragmentsInUnscaffoldedChunks + fragmentsInUnlinkedChunks);
	fprintf(stderr, "\n");
	fprintf(stderr, "chunkDepth1Contains: %d\n", chunkDepth1Contains);
	fprintf(stderr, "fragmentsInChunkDepth1Contains: %d\n", fragmentsInChunkDepth1Contains);
	fprintf(stderr, "\n");
  }
  fprintf(stderr,
		  "------------------------------: fragment counting ending\n");
}

int CountFragments2(float (* quality)(CIEdgeT *, int32), float qualityThresh)
{
  ChunkInstanceT 
	* contig;
  GraphNodeIterator
	ContigGraphIterator;
  int
	iloop,
	contigsWithValidEdges = 0,
	previousContigId = -1,
	containedContigsWithValidEdges = 0,
	previousOverlappedContigId = 0,
	totalEdges = 0,
	validEdges = 0,
	validEdgesToPlacedContigs = 0,
	placedContigs = 0,
	fragmentsInScaffoldedContigs = 0,
	* fragmentsInUnscaffoldedContigs,
	fragmentsInUnscaffoldedContigsWithValidEdges = 0,
	fragmentsInParentChunks = 0,
	scaffoldedContigs = 0,
	* unscaffoldedContigs,
	* timesContained,
	unscaffoldedContigsWithValidEdges = 0,
	parentChunks = 0,
	fragmentCount = 0,
	overlapsAttempted = 0,
	numContigs = 0;
  int 
	contigDepth,
	* contigEdgeDepth,
	contigLabelled,
	fragmentsInContigDepth1Contains = 0,
	contigDepth1Contains = 0,
	repeatIDCount[NUM_REPEAT_TYPES + 1];  // since they are numbered starting at 1
  char 
	*repeatIDType[NUM_REPEAT_TYPES + 1] = { "", "rDNA(5S)", "Histone", "Simple Repeat", "Satellite", "DROS Heterochromatin", 
											"Moderate Repeat", "Transposon", "", "", "On the fly"};
  

#define FRAG_DIST_MAX 2048
  int unlinkedContigFragmentDist[FRAG_DIST_MAX];

  fprintf(stderr,
		  "------------------------------: fragment counting2 starting\n");

  numContigs = GetNumGraphNodes(ScaffoldGraph->ContigGraph);
  fprintf(stderr, "\n");
  fprintf(stderr, "numContigs = %d\n", numContigs);

  contigEdgeDepth = (int *) malloc(numContigs * sizeof(int));
  assert (contigEdgeDepth);
  
  fragmentsInUnscaffoldedContigs = (int *) malloc(numContigs * sizeof(int));
  assert (fragmentsInUnscaffoldedContigs);

  unscaffoldedContigs = (int *) malloc(numContigs * sizeof(int));
  assert (unscaffoldedContigs);

  timesContained = (int *) malloc(numContigs * sizeof(int));
  assert (timesContained);

  // scan all the nodes (contigs)
  InitGraphNodeIterator(&ContigGraphIterator, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while (contig = NextGraphNodeIterator(&ContigGraphIterator))
  {
	contigEdgeDepth[contig->id] = -2; // initialize
	timesContained[contig->id] = 0;
	
	// count all the contigs that are scaffolded
	if (contig->scaffoldID != -1)
	{
	  MultiAlignT *maContig = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id);
	  assert (maContig);

	  // fprintf(stderr, "assigning contig %d as a scaffolded contig\n", contig->id);
	  
	  contigEdgeDepth[contig->id] = 0; // 0 => contig in scaffold
	  scaffoldedContigs++;
	  fragmentsInScaffoldedContigs += GetNumIntMultiPoss(maContig->f_list);
	}	
  }
  
  contigDepth = 0;
  contigLabelled = 1;
  while (contigLabelled)
  {
	contigLabelled = 0;
	contigDepth++;
	unscaffoldedContigs[contigDepth] = 0;
	fragmentsInUnscaffoldedContigs[contigDepth] = 0;
	
	for (iloop = 0; iloop < FRAG_DIST_MAX; iloop++)
	  unlinkedContigFragmentDist[iloop] = 0;
	
	// fprintf( stderr, "\nAt contigDepth: %d\n", contigDepth);

	// scan all the nodes (contigs)
	InitGraphNodeIterator(&ContigGraphIterator, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
	while (contig = NextGraphNodeIterator(&ContigGraphIterator))
	{
	  CIEdgeT
		* e;
	  GraphEdgeIterator
		edges;
	  int32
		paired_contig_id;
	  ChunkInstanceT
		* paired_contig;
	  
	  // deal only with those contigs that are not in a scaffold and are not parents and not already labelled
	  if (contig->scaffoldID == -1)
	  {		
		// fprintf( stderr, "\nExamining contig: %d\n", contig->id);
		
		// loop over all the edges from this contig
		InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, contig->id, 
							  ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
		while (e = NextGraphEdgeIterator(&edges)) 
		{
		  assert(e != NULL);
		  if (contig->id == e->idA)
			paired_contig_id = e->idB;
		  else
			paired_contig_id = e->idA;
		  paired_contig = GetGraphNode(ScaffoldGraph->ContigGraph, paired_contig_id);
		  assert(paired_contig != NULL);
		  
		  // fprintf( stderr, "\nExamining edge from contig %d to contig %d\n", contig->id, paired_contig->id);
		  
		  if (!Is_Not_Bogus_Not_Guide(e))  // don't play with funky edges
			continue;
		  
          // we only want to link off of contigs that are themselves assigned previously
		  if (contigEdgeDepth[paired_contig->id] < contigDepth && contigEdgeDepth[paired_contig->id] >= 0)  
		  {

			// fprintf(stderr, "assigning contig %d as a linked contig\n", contig->id);
			
			if (contigEdgeDepth[contig->id] == -2)  // don't add in frags for this contig more than once on different edges
			{
			  MultiAlignT *maContig = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id);
			  assert (maContig);

			  fragmentsInUnscaffoldedContigs[contigDepth] += GetNumIntMultiPoss(maContig->f_list);;
			  unscaffoldedContigs[contigDepth]++;
			  contigEdgeDepth[contig->id] = contigDepth;
			  contigLabelled = 1;
			}			

#if 0
			// determine if contigs with a single link are completely contained in a contig
			if (contigDepth == 1)
			{
			  int32
				paired_contig_left_end, paired_contig_right_end,
				pairedContigScaffoldOrientation,
				contigLeftLimit, contigRightLimit, contigOrientation;
			  ContigInstanceT 
				*contig,
				*paired_contig;

			  contig = GetGraphNode(ScaffoldGraph->ContigGraph, contig->info.CI.contigID);
			  paired_contig = GetGraphNode(ScaffoldGraph->ContigGraph, paired_contig->info.CI.contigID);

			  // find out where the contig with the edge to it is placed
			  GetContigPositionInScaffold(paired_contig, &paired_contig_left_end, &paired_contig_right_end, 
										  &pairedContigScaffoldOrientation);
		
			  // now determine the limits of where our contig can be
			  GetContigPlacementLimits(e, pairedContigScaffoldOrientation, paired_contig_left_end, paired_contig_right_end,
									   contig, &contigLeftLimit, &contigRightLimit, &contigOrientation, TRUE);

			  // see if contig is contained in placed contig - is this limitation necessary?
			  if (paired_contig_left_end < contigLeftLimit && contigRightLimit < paired_contig_right_end)
			  {
				timesContained[i]++;
			  }
			  else if (0)
			  {
				fprintf( stderr, "contig: %d, %d\t paired_contig: %d %d\n", contigLeftLimit, contigRightLimit, 
						 paired_contig_left_end, paired_contig_right_end);
			  }
			}
#endif
			// done with this contig
			// break;
		  }
		}
	  }
	}
  }


  // we are looping throught the contig graph twice to do the unlinked contigs
  // need only do it once
  fprintf(stderr, "\n");
  {
	int totalUnscaffoldedContigs = 0,
	  unlinkedContigs = 0;

	fprintf(stderr, "scaffoldedContigs: %d\n", scaffoldedContigs);
	for (iloop = 1; iloop < contigDepth; iloop++)
	{
	  totalUnscaffoldedContigs += unscaffoldedContigs[iloop];
	  fprintf(stderr, "unscaffoldedContigs[%d]: %d\n", iloop, unscaffoldedContigs[iloop]);
	}
	fprintf(stderr, "totalUnscaffoldedContigs = %d\n", totalUnscaffoldedContigs);

	InitGraphNodeIterator(&ContigGraphIterator, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
	while (contig = NextGraphNodeIterator(&ContigGraphIterator))
	{
	  if (contigEdgeDepth[contig->id] == -2)
	  {
		// fprintf(stderr, "assigning contig %d as an unlinked contig\n", contig->id);
		unlinkedContigs++;
	  }
	}
	fprintf(stderr, "unlinkedContigs: %d\n", unlinkedContigs);
	fprintf(stderr, "Total Contigs: %d\n", scaffoldedContigs + totalUnscaffoldedContigs + unlinkedContigs);
	fprintf(stderr, "\n");
  }
  
  {
	int totalFragmentsInUnscaffoldedContigs = 0,
	  fragmentsInUnlinkedContigs = 0;

	fprintf(stderr, "fragmentsInScaffoldedContigs: %d\n", fragmentsInScaffoldedContigs);
	for (iloop = 1; iloop < contigDepth; iloop++)
	{
	  totalFragmentsInUnscaffoldedContigs += fragmentsInUnscaffoldedContigs[iloop];
	  fprintf(stderr, "fragmentsInUnscaffoldedContigs[%d]: %d\n", iloop, fragmentsInUnscaffoldedContigs[iloop]);
	}
	fprintf(stderr, "totalFragmentsInUnscaffoldedContigs = %d\n", totalFragmentsInUnscaffoldedContigs);

	for (iloop = 0; iloop < NUM_REPEAT_TYPES + 1; iloop++)
	  repeatIDCount[iloop] = 0;

	InitGraphNodeIterator(&ContigGraphIterator, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
	while (contig = NextGraphNodeIterator(&ContigGraphIterator))
	{
	  int j;

	  // if (contig->id < 10) contigEdgeDepth[contig->id] = -2; // testing printf

#if 0
	  if (timesContained[i] > 0)
	  {
		MultiAlignT *maContig;
		assert (maContig);

		contigDepth1Contains++;
		contig = GetGraphNode(ScaffoldGraph->ContigGraph, i);
		maContig = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id);
		fragmentsInContigDepth1Contains += GetNumIntMultiPoss(maContig->f_list);
	  }
#endif	  

	  if (contigEdgeDepth[contig->id] == -2)
	  {
		MultiAlignT *maContig;

		// contig = GetGraphNode(ScaffoldGraph->ContigGraph, i);
		maContig = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id);
		assert (maContig);

		fragmentsInUnlinkedContigs += GetNumIntMultiPoss(maContig->f_list);

		// fprintf(stderr, "unlinked contig %d has %d fragments\n", contig->id, GetNumIntMultiPoss(maContig->f_list));
		// fprintf(stderr, "unlinked contig %d baseID: %d\n", contig->id, contig->info.CI.baseID);
		// fprintf(stderr, "unlinked contig %d numInstances: %d\n\n", contig->id, contig->info.CI.numInstances);

		unlinkedContigFragmentDist[GetNumIntMultiPoss(maContig->f_list)]++;

		// if (GetNumIntMultiPoss(maContig->f_list) >= 5)
		{
		  MultiAlignT *maContig = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id);
		  assert (maContig);

		  for (j = 0; j < GetNumIntMultiPoss(maContig->f_list); j++)
		  {
			IntMultiPos *mpContig = GetIntMultiPos(maContig->f_list, j);
			// fprintf(stderr, "unlinked contig %d contains fragment: %d\n", contig->id, mpContig->ident);
			getFragmentScreenInfo(mpContig->ident, repeatIDCount);
		  }
		  if (GetNumIntMultiPoss(maContig->f_list) >= 5)
			printEdgeInfoOnContig(contig, contigEdgeDepth);
		  // fprintf(stderr, "unlinked contig %d consensus: %s\n\n", contig->id, Getchar(maContig->consensus, 0));
		}
	  }
	}
	fprintf(stderr, "fragmentsInUnlinkedContigs = %d\n", fragmentsInUnlinkedContigs);
	fprintf(stderr, "Total Fragments: %d\n", fragmentsInScaffoldedContigs + 
			totalFragmentsInUnscaffoldedContigs + fragmentsInUnlinkedContigs);
	fprintf(stderr, "\n");
	for (iloop = 1; iloop < FRAG_DIST_MAX; iloop++)
	  if (unlinkedContigFragmentDist[iloop] > 0)
		fprintf( stderr, "unlinkedContigFragmentDist[%d] = %d\n", iloop, unlinkedContigFragmentDist[iloop]);
	fprintf(stderr, "\n");
	fprintf( stderr, "Total fragments in unlinked contigs with screen matches: = %d\n", repeatIDCount[0]);
	for (iloop = 1; iloop < NUM_REPEAT_TYPES + 1; iloop++)
	  fprintf(stderr, "repeatIDCount[%d] = %d (%s)\n", iloop, repeatIDCount[iloop], repeatIDType[iloop]);
	fprintf(stderr, "\n");
	fprintf(stderr, "contigDepth1Contains: %d\n", contigDepth1Contains);
	fprintf(stderr, "whereis ****************** fragmentsInContigDepth1Contains: %d\n", fragmentsInContigDepth1Contains);
	fprintf(stderr, "\n");
  }
  fprintf(stderr,
		  "------------------------------: fragment counting2 ending\n");
}

int printEdgeInfoOnChunk(int chunkID, int *chunkEdgeDepth)
{
  CIEdgeT
	* e;
  GraphEdgeIterator
	edges;
  int32
	paired_chunk_id;
  ChunkInstanceT
	* chunk,
	* paired_chunk;
  
  chunk = GetGraphNode(ScaffoldGraph->CIGraph, chunkID);
  
  // loop over all the edges from this chunk
  InitGraphEdgeIterator(ScaffoldGraph->CIGraph, chunk->id, 
						ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
  while (e = NextGraphEdgeIterator(&edges)) 
  {
	assert(e != NULL);
	if (chunk->id == e->idA)
	  paired_chunk_id = e->idB;
	else
	  paired_chunk_id = e->idA;
	paired_chunk = GetGraphNode(ScaffoldGraph->CIGraph, paired_chunk_id);
	assert(paired_chunk != NULL);
	
	fprintf( stderr, "Edge from chunk %d to chunk %d, chunkEdgeDepth[%d]: %d, Is_Not_Bogus_Not_Guide: %d\n", 
			 chunk->id, paired_chunk->id, paired_chunk->id, chunkEdgeDepth[paired_chunk->id], Is_Not_Bogus_Not_Guide(e));	
  }
}

int printEdgeInfoOnContig(ChunkInstanceT *contig, int *contigEdgeDepth)
{
  CIEdgeT
	* e;
  GraphEdgeIterator
	edges;
  int32
	paired_contig_id;
  ChunkInstanceT
	* paired_contig;
  
  // loop over all the edges from this contig
  InitGraphEdgeIterator(ScaffoldGraph->CIGraph, contig->id, 
						ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
  while (e = NextGraphEdgeIterator(&edges)) 
  {
	assert(e != NULL);
	if (contig->id == e->idA)
	  paired_contig_id = e->idB;
	else
	  paired_contig_id = e->idA;
	paired_contig = GetGraphNode(ScaffoldGraph->CIGraph, paired_contig_id);
	assert(paired_contig != NULL);
	
	fprintf( stderr, "Edge from contig %d to contig %d,Is_Not_Bogus_Not_Guide: %d\n", 
			 contig->id, paired_contig->id, Is_Not_Bogus_Not_Guide(e));	
	fprintf( stderr, "contigEdgeDepth[%d]: %d, contigEdgeDepth[%d]: %d\n", 
			 contig->id, contigEdgeDepth[contig->id], paired_contig->id, contigEdgeDepth[paired_contig->id]);	
	fprintf( stderr, "\n");
  }
}

	
#if 0
int IsChunkContainedInContig(float (* quality)(CIEdgeT *, int32), float qualityThresh)
{
  ChunkInstanceT 
	* contig;
  GraphNodeIterator
	ContigGraphIterator;
  int
	unscaffoldedContigs = 0,
	contigsWithValidEdges = 0,
	previousContigId = -1,
	containedContigsWithValidEdges = 0,
	previousOverlappedContigId = 0,
	totalEdges = 0,
	validEdges = 0,
	validEdgesToPlacedContigs = 0,
	placedContigs = 0,
	overlapsAttempted = 0;

    fprintf(stderr,
		  "------------------------------: fragment placing (unscaffolded contigs) starting\n");

  InitGraphNodeIterator(&ContigGraphIterator, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while (contig = NextGraphNodeIterator(&ContigGraphIterator))
  {
	CIEdgeT
	  * e;
	GraphEdgeIterator
	  edges;
	int32
	  overlapping_contig_cid,
	  overlapping_contig_left_end, overlapping_contig_right_end,
	  overlappingContigScaffoldOrientation,
	  contigLeftLimit, contigRightLimit, contigOrientation;
	ChunkInstanceT
	  * overlapping_contig,
	  * bestContig;
	double bestQuality = 1.0;
	
    //
    // scan all the nodes (contigs)
    //
	
	// if numInstances > 0 then it has a surrogate, and we don't want it
	// but does numInstances apply only to chunks or to contigs as well? Ask Saul.
	if (!(contig->scaffoldID == -1 && contig->info.CI.numInstances == 0))
	  continue;
	
	unscaffoldedContigs++;
	fprintf( stderr, "\nExamining contig: %d\n", contig->id);
	
	// loop over all the edges from this contig
	InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, contig->id, 
						  ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
	while (e = NextGraphEdgeIterator(&edges)) 
	{
	  MultiAlignT *maOverlappingContig, *maContig;
	  static VA_TYPE(char) *overlappingContigConsensusSequence = NULL;
	  static VA_TYPE(char) *overlappingContigConsensusQuality = NULL;
	  static VA_TYPE(char) *contigConsensusSequence = NULL;
	  static VA_TYPE(char) *contigConsensusQuality = NULL;
	  
	  assert(e != NULL);
	  if (contig->id == e->idA)
		overlapping_contig_cid = e->idB;
	  else
		overlapping_contig_cid = e->idA;
	  overlapping_contig = GetGraphNode(ScaffoldGraph->ContigGraph, overlapping_contig_cid);
	  assert(overlapping_contig != NULL);

	  totalEdges++;

	  if (!Is_Not_Bogus_Not_Guide(e))  // don't play with funky edges
		continue;
	  validEdges++;

	  if (overlapping_contig->scaffoldID != -1)  // we only want to place off of contigs that are themselves placed
	  {
		validEdgesToPlacedContigs++;
		if (contig->id != previousContigId)
		{
		  contigsWithValidEdges++;
		  previousContigId = contig->id;		  
		}

		// find out where the contig with the edge to it is placed
		GetContigPositionInScaffold(overlapping_contig, &overlapping_contig_left_end, &overlapping_contig_right_end, 
									&overlappingContigScaffoldOrientation);
		
		// now determine the limits of where our contig can be
		GetContigPlacementLimits(e, overlappingContigScaffoldOrientation, overlapping_contig_left_end, overlapping_contig_right_end,
								 contig, &contigLeftLimit, &contigRightLimit, &contigOrientation, FALSE);

		// experiment with Karin's stuff
		if (0)
		{
		  MultiAlignT  *newMultiAlign;
		  static VA_TYPE(IntElementPos) *ContigPositions = NULL;
		  IntMultiPos contigPos;
		  
		  if(ContigPositions == NULL)
			ContigPositions = CreateVA_IntElementPos(2);

		  // contigPos.type = 'C';
		  contigPos.type &= AS_CONTIG;
		  contigPos.ident = contig->id;
		  contigPos.position.bgn = overlapping_contig_left_end;
		  contigPos.position.end = overlapping_contig_right_end;
		  AppendIntElementPos(ContigPositions, &contigPos);

		  contigPos.type &= AS_CONTIG;
		  contigPos.ident = overlapping_contig->id;
		  contigPos.position.bgn = contigLeftLimit;
		  contigPos.position.end = contigRightLimit;
		  AppendIntElementPos(ContigPositions, &contigPos);
		  
		  newMultiAlign = MergeMultiAligns(ScaffoldGraph->CIGraph->maStore, 
										   ScaffoldGraph->ContigGraph->maStore, 
										   ScaffoldGraph->fragStore,
										   ContigPositions, FALSE);
		  
		  if(!newMultiAlign)
			fprintf (stderr, "MergeMultiAligns for contig %d and overlapping_contig %d failed!\n", 
					 contig->id, overlapping_contig->id);
		  else
			fprintf (stderr, "MergeMultiAligns for contig %d and overlapping_contig %d succeeded!\n", 
					 contig->id, overlapping_contig->id);
#if 0
			if(!newMultiAlign)
		  {
			fprintf(stderr,"* Calling MergeMultiAligns with %ld elements*\n",
					GetNumIntElementPoss(ContigPositions));
			{
			  int i;
			  
			  fprintf(stderr,"*** Contigging ContigPositions\n");
			  for(i = 0; i < GetNumIntElementPoss(ContigPositions); i++)
			  {
				IntElementPos *pos = GetIntElementPos(ContigPositions,i);
				ContigT *contig = GetGraphNode(ScaffoldGraph->ContigGraph, pos->ident);
				int32 originalID = GetOriginalContigID(contig->id);
				
				fprintf(stderr,"* Contig %d (%d)   bgn:%d end:%d [%d,%d]\n",
						pos->ident, originalID, pos->position.bgn, pos->position.end, contig->aEndCoord, contig->bEndCoord);
			  }
			}
			DumpCIScaffold(stderr, ScaffoldGraph, scaffold, FALSE);
			fprintf(stderr,"* MergeMultiAligns failed....bye\n");
			fflush(NULL);
			assert(0);
		  }
#endif
		} // end of Karin's stuff

		
		// then overlap them

		// see if contig is contained in placed contig - is this limitation necessary?
		if (overlapping_contig_left_end < contigLeftLimit && contigRightLimit < overlapping_contig_right_end)
		{
		  if (contig->id != previousOverlappedContigId)
		  {
			containedContigsWithValidEdges++;
			previousOverlappedContigId = contig->id;			
		  }
		  overlapsAttempted++;

		  maOverlappingContig = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, overlapping_contig->id);
		  assert (maOverlappingContig);

		  maContig = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id);
		  assert (maContig);
		  
		  fprintf( stderr, "trying to place: %d between (%d, %d) edge (%f, %c) \nplaced contig: %d, position (%d, %d), scaffold: %d, previous: %d, next: %d\n",
				   contig->id, contigLeftLimit, contigRightLimit,
				   e->distance.mean, e->orient,
				   overlapping_contig_cid, overlapping_contig_left_end, overlapping_contig_right_end,
				   overlapping_contig->scaffoldID, overlapping_contig->AEndNext, overlapping_contig->BEndNext);
		  
		  if(overlappingContigConsensusSequence == NULL)
		  {
			overlappingContigConsensusSequence = CreateVA_char(2048);
			overlappingContigConsensusQuality = CreateVA_char(2048);
			contigConsensusSequence = CreateVA_char(2048);
			contigConsensusQuality = CreateVA_char(2048);
		  }
		  
		  // grab the corresponding consensus plus slop
		  GetMultiAlignUngappedConsensus(maOverlappingContig, overlappingContigConsensusSequence, overlappingContigConsensusQuality);

		  // grab the corresponding consensus plus slop
		  GetMultiAlignUngappedConsensus(maContig, contigConsensusSequence, contigConsensusQuality);		  
	  
#if 1
		  // overlap 'em
		  {
			int A_hang_estimate,A_hang_lower,A_hang_upper;
			InternalFragMesg AFR, BFR;
			OverlapMesg     *O;
			int where, opposite = FALSE;
			FILE* log = NULL;
			int length;
			
			AFR.sequence   = Getchar(overlappingContigConsensusSequence, 0);
			AFR.quality    = Getchar(overlappingContigConsensusQuality, 0);
			AFR.eaccession = 0;
			AFR.iaccession = 0;
			BFR.sequence   = Getchar(contigConsensusSequence, 0);
			BFR.quality    = Getchar(contigConsensusQuality, 0);
			BFR.eaccession = 0;			
			BFR.iaccession = 1;

			if (e->orient == 'I' || e->orient == 'O')
			  opposite = TRUE;
			
			// add some slop
			// contigLeftLimit -= 500;
			// contigRightLimit += 500;
			contigLeftLimit = 0;
			contigRightLimit = strlen(AFR.sequence);
			
			fprintf(stderr,"* Calling DP_Compare with %s sequences lengths %ld,%ld bgn,end =[%d,%d] %c %c\n",
					(!strcmp(AFR.sequence, BFR.sequence)?"EQUAL":"DIFFERENT"),
					strlen(AFR.sequence), strlen(BFR.sequence),
					contigLeftLimit, contigRightLimit, *AFR.sequence, *BFR.sequence);

			//CGW_DP_ERATE=0.1 changed parameter to canOlap->errRate
			// a new field in ContigOverlapCheckT
			O = DP_Compare_AS(&AFR, &BFR, 
							  contigLeftLimit, contigRightLimit, opposite,
							  CGW_DP_ERATE, CGW_DP_THRESH, CGW_DP_MINLEN,
							  AS_FIND_ALIGN, //_NO_TRACE, // Slower, computes an alignment, throw away delta encoding
							  &where);
			if (O != NULL)
			{
			  // myPrint_Overlap_AS(stderr, &AFR, &BFR, O);
			  if(log == NULL)
				log = stderr;
			  
			  compute_bayesian_quality(&AFR, &BFR, O, 0, &length, log);
			  fprintf(stderr, "dorf, Quality between contig %4d and overlapping contig %4d WITHOUT QV Realigning = %lf\n", 
					  contig->id, overlapping_contig->id, O->quality);
			  if (O->quality < bestQuality)
			  {
				bestQuality = O->quality;
				bestContig = overlapping_contig;
			  }
			}		
		  }
#endif
		}
	  }
	}
	// pick the maximum, place it there
	if (bestQuality < qualityThresh)
	{
	  fprintf(stderr, "placing contig %d based on a quality value of %f with contig %d\n", contig->id, bestQuality, bestContig->id);
	  placedContigs++;





	  // wait a second!!! we need to place this guy based on the coords calculated by DP_Align, not where we started
	  // what about this?





	  //InsertCIInScaffold(ScaffoldGraphT *sgraph, int32 bestContig, int32 sid, 
	  //				   LengthT aEndOffset, LengthT bEndOffset, int AEndToBend, int contigNow);
	}
  }
  fprintf(stderr, "unscaffoldedContigs: %d\n", unscaffoldedContigs);
  fprintf(stderr, "contigsWithValidEdges: %d\n", contigsWithValidEdges);
  fprintf(stderr, "totalEdges: %d\n", totalEdges);
  fprintf(stderr, "validEdges: %d\n", validEdges);
  fprintf(stderr, "validEdgesToPlacedContigs: %d\n", validEdgesToPlacedContigs);
  fprintf(stderr, "overlapsAttempted: %d\n", overlapsAttempted);  
  fprintf(stderr, "containedContigsWithValidEdges: %d\n", containedContigsWithValidEdges);
  fprintf(stderr, "placedContigs: %d\n", placedContigs);
  if (unscaffoldedContigs > 0)
	fprintf(stderr, "placedContigs/unscaffoldedContigs: %d/%d (%.2f %%)\n", placedContigs, unscaffoldedContigs, 
			(100.0 * placedContigs) / unscaffoldedContigs);  
  else
	fprintf(stderr, "placedContigs/unscaffoldedContigs: %d/%d (%.2f %%)\n", placedContigs, unscaffoldedContigs, 
			0.0);

  if (containedContigsWithValidEdges > 0)
	fprintf(stderr, "placedContigs/containedContigsWithValidEdges: %d/%d (%.2f %%)\n", placedContigs, containedContigsWithValidEdges, 
			(100.0 * placedContigs) / containedContigsWithValidEdges);  
  else
	fprintf(stderr, "placedContigs/containedContigsWithValidEdges: %d/%d (%.2f %%)\n", placedContigs, containedContigsWithValidEdges, 
			0.0);  

  fprintf(stderr,
		  "------------------------------: fragment placing (unscaffolded contigs) ending\n");
}
#endif



// how to call Karin's stuff (from MergeAlignments.c)
// for us, unitig_store == ScaffoldGraph->ContigGraph->maStore

#if 0
    contig = GetMultiAlignInStore(contig_store,iums[0].ident);
    if (0) { // test for GetCoverage
      SeqInterval range;
      VA_TYPE(char) *cns_window;
      VA_TYPE(char) *qlt_window;
      VA_TYPE(int)  *cov_window;
      cns_window = CreateVA_char(10000);
      qlt_window = CreateVA_char(10000);
      cov_window = CreateVA_int(10000);
      range.bgn = 10;
      range.end = 30;
      GetCoverageInMultiAlignT(unitig_store, contig, range, cns_window,qlt_window, 
                cov_window, frag_store);
      fprintf(stderr,"Successful run of GetCoverage\n");
      DeleteVA_char(cns_window); 
      DeleteVA_char(qlt_window); 
      DeleteVA_int(cov_window); 
    }
#endif

#if 0

// the func from GraphCGW_T.c
int SurrogateChunkDistance(GraphCGW_T *graph,  
		  CIFragT *frag, 
		  CIFragT *mfrag, 
		  DistT *dist, 
		  int type, 
		  int orientIn, 
	      GraphEdgeStatT *stat,
          LengthT *distance){
  int num = 0;
  int mate;
  int isCI;
  NodeCGW_T *node, *mnode;
  CIEdgeT ciEdge;
  int32 fragID = GetVAIndex_CIFragT(ScaffoldGraph->CIFrags, frag);
  int32 mfragID = GetVAIndex_CIFragT(ScaffoldGraph->CIFrags, mfrag);
  LengthT ciOffset, mciOffset;
  int32 extremalA, extremalB;
  int32 minDistance, maxDistance;
  ChunkOrientationType ciOrient, mciOrient, ciEdgeOrient;
  LengthT chunkOffset;
  FragOrient chunkOrient;
  EdgeStatus status;
  int insert = TRUE ; // (type == AS_MATE); // BOGUS

  fprintf(stderr,"* SurrogateChunkDistance frags (%d,%d) dist:%g orientIn %d\n",
		  fragID, mfragID, dist->mean, orientIn);

  switch(graph->type){
  case CI_GRAPH:
#if 0
    if(frag->cid == mfrag->cid)
      return FALSE;
#endif
    isCI = TRUE;
    node = GetGraphNode(graph, frag->cid);
    mnode = GetGraphNode(graph, mfrag->cid);
    break;
  case CONTIG_GRAPH:
    if(frag->contigID == mfrag->contigID)
      return FALSE;
    isCI = FALSE;
    node = GetGraphNode(graph, frag->contigID);
    mnode = GetGraphNode(graph, mfrag->contigID);
    break;
  default:
    assert(0);
  }

  fprintf(stderr,"* SurrogateChunkDistance; get graph nodes\n");

  if(0)
    fprintf(stderr,
	    "* Found mate of frag %d (%d) in chunk:%d mate:%d (%d) chunk:%d\n",
	    fragID, frag->iid,
	    frag->cid, 
	    mfragID, mfrag->iid,
		mfrag->cid);
  


  if(!FragOffsetAndOrientation(frag,
				 node,
				 &ciOffset,
				 &ciOrient,
				 &extremalA,
				 (orientIn == AS_GKP_INNIE ||
				  orientIn == AS_GKP_NORMAL))){
    fprintf(stderr,"* CIFragOffset returned false for frag %d (%d)\n",
	    fragID, frag->iid);
    return FALSE;
  }

  fprintf(stderr,"* CIFragOffset returned true for frag %d (%d)\n",
		  fragID, frag->iid);

  if(!FragOffsetAndOrientation(mfrag,
				 mnode,
				 &mciOffset,
				 &mciOrient,
				 &extremalB,
				 (orientIn == AS_GKP_INNIE ||
				  orientIn == AS_GKP_ANTINORMAL))){
    fprintf(stderr,"* CIFragOffset returned false for frag %d (%d)\n",
	    mfragID, mfrag->iid);
    return FALSE;
  }

  fprintf(stderr,"* CIFragOffset returned true for frag %d (%d)\n",
		  mfragID, mfrag->iid);

  /* The following triply nested case statement captures all of the cases that arise from different
     relative alignments of the fragments in the LKG relationship, and their alignment with their 
     respective chunks.
     There is probably a better way to do this, but I think this is the clearest way to codify the relationships,
     complete with 'drawings'
  */

  switch(orientIn){
	  
  case AS_GKP_INNIE: /******************************* AB_BA *******************************/
    switch(ciOrient){
      //
    case A_B:

      switch(mciOrient){
      case A_B:
	//           length - 5'             gap            length - 5'
	//      |------------------------||---------------||-----------|
	//  A --------------------------- B               B --------------------------- A
	//    5'----->                                           <------5'
	//      |-------------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = AB_BA;  // A
	break;
      case B_A:
	//           length - 5'             gap                5'
	//      |------------------------||---------------||-----------|
	//  A --------------------------- B               A --------------------------- B
	//    5'----->                                           <------5'
	//      |-------------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = AB_AB; // N
	break;
      default:
	assert(0);
	break;
      }
      break;
    case B_A:

      switch(mciOrient){
      case A_B:
	//                     5'             gap            length - 5'
	//      |------------------------||---------------||-----------|
	//  B --------------------------- A               B --------------------------- A
	//    5'----->                                           <------5'
	//      |-------------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = BA_BA; // I
	break;
      case B_A:
	//                     5'             gap                5'
	//      |------------------------||---------------||-----------|
	//  B --------------------------- A               A --------------------------- B
	//    5'----->                                           <------5'
	//      |-------------------------------------------------------|
	//                             mate/guide distance
	//
	ciEdgeOrient = BA_AB; // O
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
    break;

  case AS_GKP_NORMAL: /******************************* AB_AB *******************************/
    switch(ciOrient){
      //
    case A_B:

      switch(mciOrient){
      case B_A:
	//           length - 5'             gap              Length - 5'
	//      |------------------------||---------------||-----------|
	//  A --------------------------- B               B --------------------------- A
	//    5'----->                                                  5'------>
	//      |-------------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = AB_BA;  // A
	break;
      case A_B:
	//           length - 5'             gap                5'
	//      |------------------------||---------------||------|
	//  A --------------------------- B               A --------------------------- B
	//    5'----->                                           5'------>
	//      |-------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = AB_AB; // N
	break;
      default:
	assert(0);
	break;
      }
      break;
    case B_A:

      switch(mciOrient){
      case A_B:
	//                     5'             gap            5'
	//      |------------------------||---------------||----|
	//  B --------------------------- A               A --------------------------- B
	//    5'----->                                          5'------>
	//      |-----------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = BA_AB; // O
	break;
      case B_A:
	//                     5'             gap                Length - 5'
	//      |------------------------||---------------||----|
	//  B --------------------------- A               B --------------------------- A
	//    5'----->                                          5'------>
	//      |-----------------------------------------------|
	//                             mate/guide distance
	//
	ciEdgeOrient = BA_BA; // A
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
    break;
  case AS_GKP_ANTINORMAL: /******************************* BA_BA *******************************/
    switch(ciOrient){
      //
    case A_B:

      switch(mciOrient){
      case B_A:
	//                 5'             gap                     5'
	//          |---------------------||---------------||------------------|
	//  B --------------------------- A               A --------------------------- B
	//    <-----5'                                                  <------5'
	//          |-----------------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = BA_AB; 
	break;
      case A_B:
	//                 5'             gap                     Length - 5'
	//          |---------------------||---------------||------------------|
	//  B --------------------------- A               B --------------------------- A
	//    <-----5'                                                  <------5'
	//          |-----------------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = BA_BA; // N
	break;
      default:
	assert(0);
	break;
      }
      break;
    case B_A:

      switch(mciOrient){
      case A_B:
	//                 Length - 5'             gap               Length - 5'
	//          |---------------------||---------------||------------------|
	//  A --------------------------- B               B --------------------------- A
	//    <-----5'                                                  <------5'
	//          |-----------------------------------------------------------|
	//                             mate distance
	//

	ciEdgeOrient = AB_BA; 
	break;
      case B_A:
	//                 length - 5'           gap                     5'
	//          |---------------------||---------------||------------------|
	//  A --------------------------- B               A --------------------------- B
	//    <-----5'                                                  <------5'
	//          |-----------------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = BA_BA; // A
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
    break;
  case AS_GKP_OUTTIE: /******************************* BA_AB *******************************/
    switch(ciOrient){
      //
    case A_B:

      switch(mciOrient){
      case B_A:
	//                 5'             gap                     Length - 5'
	//          |---------------------||---------------||------------ |
	//  B --------------------------- A               B --------------------------- A
	//    <-----5'                                                  5'------>
	//          |-----------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = BA_BA; 
	break;
      case A_B:
	//                 5'             gap                      5'
	//          |---------------------||---------------||-----------|
	//  B --------------------------- A               A --------------------------- B
	//    <-----5'                                                  5' ------>
	//          |---------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = BA_AB; // N
	break;
      default:
	assert(0);
	break;
      }
      break;
    case B_A:

      switch(mciOrient){
      case B_A:
	//                 Length - 5'          gap                     Length - 5'
	//          |---------------------||---------------||------------ |
	//  A --------------------------- B               B --------------------------- A
	//    <-----5'                                                  5'------>
	//          |-----------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = AB_BA; 
	break;
      case A_B:
	//                 Length - 5'            gap                      5'
	//          |---------------------||---------------||-----------|
	//  A --------------------------- B               A --------------------------- B
	//    <-----5'                                                  5' ------>
	//          |---------------------------------------------------|
	//                             mate distance
	//
	ciEdgeOrient = AB_AB; // N
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
    break;
  case XX_XX:
  default:
    assert(0);
  }


  fprintf(stderr, "ciOffset.mean = %f\n", ciOffset.mean);
  fprintf(stderr, "ciOffset.variance = %f\n", ciOffset.variance);
  fprintf(stderr, "mciOffset.mean = %f\n", mciOffset.mean);
  fprintf(stderr, "mciOffset.variance = %f\n", mciOffset.variance);  

  fprintf(stderr, "distance->mean = %f\n", distance->mean);  
  fprintf(stderr, "dist->mean = %f\n", dist->mean);  

  distance->mean = (int32)dist->mean + ciOffset.mean + mciOffset.mean;

  fprintf(stderr, "dist->stddev = %f\n", dist->stddev);  

  // Since the two offsets and the dist are independent we SUM their variances
  distance->variance = dist->stddev * dist->stddev  + ciOffset.variance + mciOffset.variance;
  
  return TRUE;
}

#endif

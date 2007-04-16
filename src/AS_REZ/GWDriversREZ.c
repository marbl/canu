
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

        Module:  GWDriversREZ.c

   Description:  main functions for the gap walker and the biconnected
                 components

		 log info is sent to <inputFile>.gwlog

    Programmer:  S. Lonardi (stelo@cs.purdue.edu)

       Written:  7 July 99

 **********************************************************************/


static char fileID[] = "$Id: GWDriversREZ.c,v 1.10 2007-04-16 17:34:15 brianwalenz Exp $";

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_timer.h"
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
#include "SubgraphREZ.h"
#include "ConsistencyChecksREZ.h"
#include "GWDriversREZ.h"
#include "StatisticsREZ.h"
#include "FbacREZ.h"
#include "RepeatRez.h"

//
// globals
//
char
* GW_Filename_Prefix;
FILE
* gwlogfp = NULL;








//
// constants for testing the GapWalker (to be removed at some point)
//
#define WALK_BEGIN_POS     1090000   
// position for testing
#define WALK_END_POS       1170000   
// "
#define WALK_BEGIN_S_ID          8   
// left scaffold (to test the walker)
#define WALK_END_S_ID            6   
// right scaffold (to test the walker)
#define WALK_BEGIN_CID          50   
// begin chunk id
#define WALK_END_CID            60   
// end chunk id
#define ANNOTATE_CONFIRMING_PATH 0   
// set this to 1 to run the path confirmation on CIedges
#define CHECK_MULTIPLE_USE       0   
// set this to 1 to check if chunk happen to be used in multiple paths

#define LOW_BAYESIAN_THRESHOLD  0.3
#define HIGH_BAYESIAN_THRESHOLD 1.0

#define BAC_WALKING_THRESHOLD  0.3

#define CHECKPOINT_DISTANCE 10000000    
// this is the number of basepairs we have inspected in scaffolds before
// gapwalking issues another checkpoint


typedef  struct
{
  int  scaff_id;
  int  scaff_length;
}  Scaffold_Size_t;















int Walk_Gaps(Global_CGW * data,
			  char * prefix,
			  int level,
			  int startWalkFrom,
			  double gapSizeStdDevs) 
{
  //
  // main Walker function. it calls the appropriate function given the
  // <level> as input.  The debug information are stored in "*.gwlog"
  //
  //
  char
    filename[STR_LEN];
  int doSmallContigRemoval;

  GW_Filename_Prefix = prefix;

  sprintf(filename, "%s.gwlog", GW_Filename_Prefix);
  gwlogfp = file_open(filename,"w");
  assert(gwlogfp != NULL);
  
  if (gapSizeStdDevs == AGGRESSIVE_WALKING_STD_DEVS)
	doSmallContigRemoval = TRUE;
  else
	doSmallContigRemoval = FALSE;

  switch (level) {
	case 1: // intra only
	{
	  Intra_Scaffold_Path_Finding( startWalkFrom, gapSizeStdDevs, doSmallContigRemoval, FALSE );
	  // Inter_Scaffold_Walking();
	  // Inter_Scaffold_Analysis();
	  // localeCam();
	  break;
	}
	case 2: // inter only
	  // Inter_Scaffold_Analysis();
	  break;
	case 3: // intra then inter
	  // Intra_Scaffold_Path_Finding( startWalkFrom , gapSizeStdDevs );
	  // Inter_Scaffold_Analysis();
	  break;
	default:
	  fprintf(stderr, "*** Error: walk level %d not defined\n", level);

	  fclose(gwlogfp);
	  exit(-1);
	  break;
  }

  fclose(gwlogfp);

  return 0;
}

int Compare_Scaff_Size(const void *e1,
					   const void *e2) {
  //
  // The function that compares two scaffold sizes
  //
  Scaffold_Size_t *s1 = (Scaffold_Size_t *)e1;
  Scaffold_Size_t *s2 = (Scaffold_Size_t *)e2;

  assert((s1 != NULL) && (s2 != NULL));

  if(GlobalData->walkScaffoldsBiggestFirst == TRUE)  
    // the default, biggest scaffolds first
  {
	if (s1->scaff_length > s2->scaff_length)
	  return -1;
	else if (s1->scaff_length < s2->scaff_length)
	  return 1;
	else
	  return 0;
  }
  else // smallest scaffolds first
  {
	if (s1->scaff_length > s2->scaff_length)
	  return 1;
	else if (s1->scaff_length < s2->scaff_length)
	  return -1;
	else
	  return 0;
  }
}

int Sort_ScaffoldEnds_By_Locale(const void *e1,	const void *e2) 
{
  //
  // The function that compares two scaffold ends by locale
  //
  ScaffoldEndT *s1 = (ScaffoldEndT *) e1;
  ScaffoldEndT *s2 = (ScaffoldEndT *) e2;
  
  assert((s1 != NULL) && (s2 != NULL));
  
  // smallest localeNumber first
  if (s1->localeNumber > s2->localeNumber)
	return 1;
  else if (s1->localeNumber < s2->localeNumber)
	return -1;
  else  // then smallest leftExtremalFragID
  {
	if (s1->leftBacOffset > s2->leftBacOffset)
	  return 1;
	else
	  return -1;
  }
}

int Sort_ScaffoldEnds_By_Gap(const void *e1, const void *e2) 
{
  //
  // The function that compares two scaffoldEndTs by contig pair
  //
  ScaffoldEndT *s1 = (ScaffoldEndT *) e1;
  ScaffoldEndT *s2 = (ScaffoldEndT *) e2;
  
  assert((s1 != NULL) && (s2 != NULL));
  
  // smallest leftContigID first
  if (s1->leftContigID > s2->leftContigID)
	return 1;
  else if (s1->leftContigID < s2->leftContigID)
	return -1;
  else  // then smallest rightContigID
  {
	if (s1->rightContigID > s2->rightContigID)
	  return 1;
	else if (s1->rightContigID < s2->rightContigID)
	  return -1;
	else
	  return 0;
  }
}

int Sort_Gaps_By_BAC_Count(const void *e1, const void *e2) 
{
  //
  // The function that compares two scaffold ends by BAC count, then by frag separation
  //
  ScaffoldEndT *s1 = (ScaffoldEndT *) e1;
  ScaffoldEndT *s2 = (ScaffoldEndT *) e2;
  
  assert((s1 != NULL) && (s2 != NULL));
  
  if (s1->spanningBACs > s2->spanningBACs)
	return -1;
  else if (s1->spanningBACs < s2->spanningBACs)
	return 1;
  else
  {
	// smallest leftContigID first
	if (s1->leftContigID > s2->leftContigID)
	  return 1;
	else if (s1->leftContigID < s2->leftContigID)
	  return -1;
	else  // then smallest rightContigID
	{
	  if (s1->rightContigID > s2->rightContigID)
		return 1;
	  else if (s1->rightContigID < s2->rightContigID)
		return -1;
	  else  // then frag separation
	  {
		int s1_frag_separation = s1->rightExtremalFragID - s1->leftExtremalFragID;
		int s2_frag_separation = s2->rightExtremalFragID - s2->leftExtremalFragID;
		
		if (s1_frag_separation > s2_frag_separation)
		  return -1;
		else if (s1_frag_separation < s2_frag_separation)
		  return 1;
		else
		  return 0;
	  }
	}
  }
}

void Intra_Scaffold_Path_Finding( int startWalkFrom, double gapSizeStdDevs, int doSmallContigRemoval, int trimScaffoldEnds) 
{
  //
  // the main intrascaffold gap walking
  //

  int
    sc;
  
  WalkStatisticsT stats;
  Scaffold_Size_t* sortedScaffolds;
  // ScaffoldWalkStatisticsT *scaffStatp;
  
  int walksAttempted = 0,
	walksSuccessful = 0,
	walksUnsuccessful = 0,
        // walksUnsuccessfulGapTooSmall = 0,
	walksUnsuccessfulFailedOverlap = 0,
	walksUnsuccessfulFragOutOfBounds = 0,
	walksUnsuccessfulContigInertiaInScaffold = 0,
	walksUnsuccessfulFailedOverlapStack = 0,
	walksUnsuccessfulRchunkContained = 0,
	walksUnsuccessfulFragNotSet = 0,
	walksUnsuccessfulContigInertiaOtherScaffold = 0,
	walksUnsuccessfulFragsDoNotOverlapInLocale = 0,
	walksUnsuccessfulDistanceDifference = 0,
	walksUnsuccessfulNoCommonLocales = 0;

  int numGapsBeforeWalking = 0;
  double basesInThisGap, 
	basesInGapsBeforeWalking = 0.0;

  long cycles = 0;
  // for timing purposes

  ChunkInstanceT
    * lchunk = NULL,
    * rchunk = NULL;
  LengthT
    gapEstimate;
  int numScaffoldsWalked = 0;
  int checkScaffold[MAX_SCAFFOLDS_CHECKED], checkScaffoldCount = 0;
  int initialNumGraphNodes;
  
  FILE *walkingStatsFile;
  
  fprintf(stderr,"------------------------------: gap walking running\n");
  
  // check for the existence of a walk statistics file
  walkingStatsFile = fopen("walkingStats", "r+");
  if (1 || walkingStatsFile == NULL)
  {
	fprintf( stderr, "File walkingStats not found.  Creating.\n");
	walkingStatsFile = fopen("walkingStats", "w");
  }
  else
  {
	fprintf( stderr, "File walkingStats found.  Reading last set of stats.\n");

	fscanf( walkingStatsFile, "%d", &numScaffoldsWalked);
	fscanf( walkingStatsFile, "%d", &numGapsBeforeWalking);
	fscanf( walkingStatsFile, "%lf", &basesInGapsBeforeWalking);
	// fscanf( walkingStatsFile, "%d", &lastWalkedScaffoldID);
	fscanf( walkingStatsFile, "%d", &walksSuccessful);
	fscanf( walkingStatsFile, "%d", &walksUnsuccessful);
	fscanf( walkingStatsFile, "%d", &walksUnsuccessfulFailedOverlap);
	fscanf( walkingStatsFile, "%d", &walksUnsuccessfulFragOutOfBounds);
	fscanf( walkingStatsFile, "%d", &walksUnsuccessfulContigInertiaInScaffold);
	fscanf( walkingStatsFile, "%d", &walksUnsuccessfulFailedOverlapStack);
	fscanf( walkingStatsFile, "%d", &walksUnsuccessfulRchunkContained);
	fscanf( walkingStatsFile, "%d", &walksUnsuccessfulFragNotSet);
	fscanf( walkingStatsFile, "%d", &walksUnsuccessfulContigInertiaOtherScaffold);
	fscanf( walkingStatsFile, "%d", &walksUnsuccessfulFragsDoNotOverlapInLocale);
  }

  fprintf( stderr, "numScaffoldsWalked: %d\n", numScaffoldsWalked);
  fprintf( stderr, "numGapsBeforeWalking: %d\n", numGapsBeforeWalking);
  fprintf( stderr, "basesInGapsBeforeWalking: %f\n", basesInGapsBeforeWalking);
  // fprintf( stderr, "lastWalkedScaffoldID: %d\n", lastWalkedScaffoldID);
  fprintf( stderr, "walksSuccessful: %d\n", walksSuccessful);
  fprintf( stderr, "walksUnsuccessful: %d\n", walksUnsuccessful);
  fprintf( stderr, "walksUnsuccessfulFailedOverlap: %d\n", walksUnsuccessfulFailedOverlap);
  fprintf( stderr, "walksUnsuccessfulFragOutOfBounds: %d\n", walksUnsuccessfulFragOutOfBounds);
  fprintf( stderr, "walksUnsuccessfulContigInertiaInScaffold: %d\n", walksUnsuccessfulContigInertiaInScaffold);
  fprintf( stderr, "walksUnsuccessfulFailedOverlapStack: %d\n", walksUnsuccessfulFailedOverlapStack);
  fprintf( stderr, "walksUnsuccessfulRchunkContained: %d\n", walksUnsuccessfulRchunkContained);
  fprintf( stderr, "walksUnsuccessfulFragNotSet: %d\n", walksUnsuccessfulFragNotSet);
  fprintf( stderr, "walksUnsuccessfulContigInertiaOtherScaffold: %d\n", walksUnsuccessfulContigInertiaOtherScaffold);
  fprintf( stderr, "walksUnsuccessfulFragsDoNotOverlapInLocale: %d\n", walksUnsuccessfulFragsDoNotOverlapInLocale);

  // sort all the scaffolds according to their size
  // then fill gaps starting in the biggest one.
  // make sure that a scaffold has not been declared dead.
  //
  sortedScaffolds = (Scaffold_Size_t*) safe_calloc(GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph), sizeof(Scaffold_Size_t));
  
  for (sc = 0; sc < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sc++)
  {
	CIScaffoldT  * scaff = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, sc);
	sortedScaffolds[sc].scaff_id = sc;
	sortedScaffolds[sc].scaff_length = (int) scaff->bpLength.mean;
  }

  qsort(sortedScaffolds,GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds), 
		sizeof(Scaffold_Size_t),Compare_Scaff_Size); 


  // allocate the statistic structs
  allocate_walk_statistics(&stats,GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));

  // we read in the statistics previously computed
  // for this computation. Note that anything read here is overriden
  // by the checkpoints, i.e. if a scaffold was not walked in a checkpoint
  // but its statistics are present from a previous run
  // it is overwritten.

#if DEBUG_GAP_WALKER > -1
  // print the scaffolds in sorted order with information whether they are dead
  for (sc = 0; sc < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sc++)
  {
	int sID = sortedScaffolds[sc].scaff_id;
	CIScaffoldT  * scaff = GetCIScaffoldT(ScaffoldGraph->CIScaffolds,sID);

	// all scaffolds that are not of relevance for walking
	// get an ID NO_SCAFFOLD. They are disregraded in the accumulative
	// statistics
	if(isDeadCIScaffoldT(scaff))
	{
	  continue;
	  fprintf(stderr,"=== Scaffold dead\n");
	}
	if( scaff->type != REAL_SCAFFOLD )
	{
	  continue;
	  fprintf(stderr,"=== Scaffold not real\n");
	}
	if( scaff->info.Scaffold.numElements < 2)
	{
	  continue;
	  fprintf(stderr,"=== Scaffold has less than two elements\n");	  
	}      
	fprintf(stderr,"=== Scaffold %d, size %d\n",sortedScaffolds[sc].scaff_id,sortedScaffolds[sc].scaff_length);
  }  
#endif

  //
  // scan all the scaffolds
  //

  initialNumGraphNodes = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
  for (sc = 0; sc < initialNumGraphNodes; sc++) 
  {
    CIScaffoldTIterator
      CIsTemp;
    CIScaffoldT
      * scaff;
	int i, currentLeftNodeId, currentRightNodeId;
	
    // Scaffold_Fill_t *fill_chunks; 
    int icnt;
    int *chunkArray;
    TimerT singleWalkAndUpdateTime;
    int sId = sortedScaffolds[sc].scaff_id;

    /* Initialize and start the timer for that scaffold */
    InitTimerT(&singleWalkAndUpdateTime); 
    StartTimerT(&singleWalkAndUpdateTime);    
    
    // if the argument startWalkFrom is not NULLINDEX we skip all scaffolds
    // except the one specified by that argument
    if( startWalkFrom != NULLINDEX )
      if( startWalkFrom != sId )
		continue;
    
    scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sId);
    // make sure the scaffold is there
    assert(scaff != NULL);   

    // if we have already walked that scaffold we skip it
    // this is when we start from a checkpoint
    //if( scaff->flags.bits.walkedAlready )
    //  continue;    
    
    // print some stats of the scaffold
#if DEBUG_GAP_WALKER > -1
    if ((isDeadCIScaffoldT(scaff)))
      fprintf(stderr,"=== Scaffold dead\n");
    if( scaff->type != REAL_SCAFFOLD )
      fprintf(stderr,"=== Scaffold not real\n");
    if(scaff->info.Scaffold.numElements < 2)
      fprintf(stderr,"=== Scaffold has less than two elements\n");
#endif
    
    //
    // not interested in dead scaffold, not real scaffolds, or singleton
    // scaffolds
    
    if ((isDeadCIScaffoldT(scaff)) ||
       	(scaff->type != REAL_SCAFFOLD) ||
		(scaff->info.Scaffold.numElements < 2))
	{
	  continue;
	}
#if DEBUG_GAP_WALKER > -1
    fprintf(stderr,"\n=====================================================================\n");
    fprintf(stderr,"=== Walking scaffold %d, size %d\n",sId,sortedScaffolds[sc].scaff_length);
#endif

	// only remove small contigs if we are doing aggressive walking
	if ( doSmallContigRemoval == TRUE)
          removeSmallContigs( scaff );
        //removeSmallContigs( scaff, trimScaffoldEnds );

	// removeSmallContigs can potentially add a scaffold, so reinitialize pointer
    scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sId);
    // make sure the scaffold is there
    assert(scaff != NULL);   

    // Scan the gaps. Here is the same
    // Scan_Gaps should also only scan gaps
    // in the current scaffold
    // fill_chunks = Scan_Gaps_In_Scaffold (sId);     

    chunkArray = (int *) safe_malloc( scaff->info.Scaffold.numElements * sizeof(int));
	
    //
    // print all the gaps of this scaffold
	// we should replace this routine and its sister (AfterWalking) below as
	// a routine called CountGapsInScaffolds or such
    //
	icnt = 0;
    InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE,
							FALSE, &CIsTemp);
    while (NextCIScaffoldTIterator(&CIsTemp))
	{
	  ChunkInstanceT
		* lchunkTemp,
		* rchunkTemp;

	  // not walking off of scaffolds currently
	  if (CIsTemp.next == -1)
		break;
	  
      //
      // find the chunks in the gap by walking between the chunk <CIs.curr>
      // and the chunk <CIs.next>
      //
            
      lchunkTemp = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.curr);
	  rchunkTemp = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.next);
	  chunkArray[icnt++] = lchunkTemp->id;
	  chunkArray[icnt] = rchunkTemp->id;
      
      assert(lchunkTemp != NULL);
      assert(rchunkTemp != NULL);

	  // fprintf( stderr, "tail iterator, gap to walk: lchunk: %d, rchunk: %d\n", lchunkTemp->id, rchunkTemp->id);

	  numGapsBeforeWalking++;
	  basesInThisGap = MIN (rchunkTemp->offsetAEnd.mean, rchunkTemp->offsetBEnd.mean) - 
		MAX (lchunkTemp->offsetAEnd.mean, lchunkTemp->offsetBEnd.mean);
	  if (basesInThisGap > 0)
		basesInGapsBeforeWalking += basesInThisGap;
	}	


#if DEBUG_GAP_WALKER > -1
	fprintf( stderr, "gap to walk: lchunk: %d, rchunk: %d\n", -1, -1);
	fprintf( stderr, "for scaffold %d, scaff->info.Scaffold.numElements: %d\n", 
			 scaff->id, scaff->info.Scaffold.numElements);
	for (icnt = 0; icnt < scaff->info.Scaffold.numElements; icnt++)
	  fprintf( stderr, "chunkArray[%d] = %d\n", icnt, chunkArray[icnt]);
#endif

    //
    // scan all the gaps of this scaffold
    //

	currentLeftNodeId = scaff->info.Scaffold.AEndCI;
	while (currentLeftNodeId != scaff->info.Scaffold.BEndCI)
	{
	  int atScaffoldEnd, completeOverlapPath = 0;
	  
	  walksAttempted++;
	  
      //
      // find the chunks in the gap by walking between the chunk <CIs.curr>
      // and the chunk <CIs.next>
      //
            
	  lchunk = GetGraphNode(ScaffoldGraph->RezGraph, currentLeftNodeId);
      assert(lchunk != NULL);
	  rchunk = GetGraphNode(ScaffoldGraph->RezGraph, lchunk->BEndNext);
      assert(rchunk != NULL);

	  // if the rchunk is contained by the lchunk keep going until we are past end of lchunk
	  // otherwise we end trying to walk a gap up "inside of" rchunk
	  atScaffoldEnd = FALSE;
	  while ( MAX( rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean) < 
			  MAX( lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean))
	  {
		fprintf( stderr, "skipping contig %d, setting rchunk to %d\n",
				 rchunk->id, rchunk->BEndNext);
		if (rchunk->BEndNext != NULLINDEX)
		  rchunk = GetGraphNode(ScaffoldGraph->RezGraph, rchunk->BEndNext);
		else  // we are at end of scaffold, the last contig is inside of rchunk
		{
		  atScaffoldEnd = TRUE;
		  break;
		}
	  }

	  // no more gaps to walk if atScaffoldEnd == TRUE
	  if (atScaffoldEnd)
		break;
	  
	  currentRightNodeId = rchunk->id;  // we need to save this because the rchunk pointer may be invalidated
	                                    // if contig array is reallocated, eg, in SplitUnresolvedContig
	  
      // We set the following bits to FALSE in the the left chunk to
      // gather statistics on the walk
      lchunk->flags.bits.walkMaxedOut   = FALSE;
      lchunk->flags.bits.walkedTooLong  = FALSE;
      lchunk->flags.bits.walkedTooShort = FALSE;
      lchunk->flags.bits.walkedTrivial  = FALSE;
      
      
#if DEBUG_GAP_WALKER > -1
	  fprintf(stderr, "\n\n\n-----------------------------------------------------\n");
	  fprintf(stderr, "Examining gap between lchunk: %d and rchunk: %d\n", lchunk->id, rchunk->id);
#endif 

      // compute the estimated gap length
      // gapEstimate = Compute_Gap_Length(lchunk, rchunk, TRUE);
      gapEstimate = FindGapLength(lchunk, rchunk, TRUE);

      // possible that the variance is set to 0.0 - makes walking impossible
      // give us some slack if necessary
      if (gapEstimate.variance < 100.0)
		gapEstimate.variance = 100.0;


#if DEBUG_GAP_WALKER > -1
      // fprintf(gwlogfp, "\nlchunk: %d, is surrogate %d  rchunk: %d, is surrogate %d gap length: %f\n", CIs.curr,IsSurrogateNode(lchunk), CIs.next, IsSurrogateNode(rchunk), gapEstimate.mean);
#endif      

	  {
		GapInfoT *gapInfoArray;
		int ilocale;
		LengthT rchunk_delta;
		ChunkOrientationType olapOrientation;
		Overlap rchunkOverlap;
		double rchunkNewAEnd, rchunkNewBEnd;
		int lchunkGapEnd, rchunkGapEnd;
		int numCommonLocales;

#if DEBUG_GAP_WALKER > -1
		fprintf(gwlogfp,"Entering Intra_Scaffold_Gap_Walker for scaffold number %d for gap of length %f (stddev : %f)\n",sId, gapEstimate.mean,sqrt(gapEstimate.variance));
#endif	
		// ****** Call the main walking routine for that gap *****		
		
		// find the locales that are at the apropos end of the left contig
		if (lchunk->offsetAEnd.mean < lchunk->offsetBEnd.mean)
		  lchunkGapEnd = B_END;
		else
		  lchunkGapEnd = A_END;
		
		// find the locales that are at the apropos end of the right contig
		if (rchunk->offsetAEnd.mean < rchunk->offsetBEnd.mean)
		  rchunkGapEnd = A_END;
		else
		  rchunkGapEnd = B_END;
		
		// first get the locales common to both sides of the gap (and check and rank by sizes)
		numCommonLocales = FindCommonLocales( lchunk, lchunkGapEnd, 
						      rchunk, rchunkGapEnd, &gapInfoArray);
		
		// temporarily set the contig positions of frags that are in surrogates
		SetFragPositions( numCommonLocales, &gapInfoArray);
		
		// dumpContigInfo(lchunk);
		// dumpContigInfo(rchunk);		  
		
		// for each locale compute how well its gap estimate matches our gap estimate 
		// and sort based on that measure
		if (numCommonLocales > 0)
		  SortCommonLocales( numCommonLocales, &gapInfoArray, 1, gapEstimate);
		else
		  walksUnsuccessfulNoCommonLocales++;
		
		// loop over locales until and if gap is walked
		for (ilocale = 0; ilocale < numCommonLocales; ilocale++)
		{
		  ChunkInsertInfoT *chunksWalked;
		  
#if DEBUG_GAP_WALKER > -1		  
		  {
			fprintf( stderr, "attempting walk using locale: %d, distanceDifference: %f\n", 
					 gapInfoArray[ilocale].localeNumber,
					 gapInfoArray[ilocale].distanceDifference);
			fprintf( stderr, "lchunk->offsetAEnd.mean: %f, lchunk->offsetBEnd.mean: %f\n",
					 lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean);
			fprintf( stderr, "rchunk->offsetAEnd.mean: %f, rchunk->offsetBEnd.mean: %f\n",
					 rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean);
			fprintf( stderr, "gapInfoArray[%d].fragIDLeftContig: %d\n",
					 gapInfoArray[ilocale].localeNumber,
					 gapInfoArray[ilocale].fragIDLeftContig);
			fprintf( stderr, "gapInfoArray[%d].fragIDRightContig: %d\n",
					 gapInfoArray[ilocale].localeNumber,
					 gapInfoArray[ilocale].fragIDRightContig);
			fprintf( stderr, "gapInfoArray[%d].unitigIDLeftContig: %d\n",
					 gapInfoArray[ilocale].localeNumber,
					 gapInfoArray[ilocale].unitigIDLeftContig);
			fprintf( stderr, "gapInfoArray[%d].unitigIDRightContig: %d\n",
					 gapInfoArray[ilocale].localeNumber,
					 gapInfoArray[ilocale].unitigIDRightContig);
			fprintf( stderr, "gapInfoArray[%d].fragInSurrogateLeftContig: %d\n",
					 gapInfoArray[ilocale].localeNumber,
					 gapInfoArray[ilocale].fragInSurrogateLeftContig);
			fprintf( stderr, "gapInfoArray[%d].fragInSurrogateRightContig: %d\n",
					 gapInfoArray[ilocale].localeNumber,
					 gapInfoArray[ilocale].fragInSurrogateRightContig);			

			// we filter based on how well the locale's measure of the gap
			// agrees with ours
			fprintf(stderr, "gapEstimate.mean: %f\n", gapEstimate.mean);
			fprintf(stderr, "gapEstimate.variance: %f\n", gapEstimate.variance);
			fprintf(stderr, "sqrt(gapEstimate.variance): %f\n", sqrt(gapEstimate.variance));
			fprintf(stderr, "%f * sqrt(gapEstimate.variance): %f\n", gapSizeStdDevs,
					gapSizeStdDevs * sqrt(gapEstimate.variance));
			fprintf(stderr, "gapInfoArray[ilocale].distanceDifference: %f\n", 
					gapInfoArray[ilocale].distanceDifference);
		  }
#endif		  
		  if (gapInfoArray[ilocale].distanceDifference > gapSizeStdDevs * sqrt(gapEstimate.variance))
		  {
		    int startFragIid, endFragIid;

		    walksUnsuccessfulDistanceDifference++;
		    startFragIid = GetCIFragT(ScaffoldGraph->CIFrags, (int32) gapInfoArray[ilocale].fragIDLeftContig)->iid;
		    endFragIid = GetCIFragT(ScaffoldGraph->CIFrags, (int32) gapInfoArray[ilocale].fragIDRightContig)->iid;
		    fprintf( stderr, "walking from startFragIid: %d (contig: %d) to endFragIid: %d (contig: %d)\n",
			     startFragIid, lchunk->id, endFragIid, rchunk->id);
			fprintf( stderr, " fragSeparationScaffold: %f\n", gapInfoArray[ilocale].fragSeparationScaffold);
			fprintf( stderr, "   fragSeparationLocale: %f\n", gapInfoArray[ilocale].fragSeparationLocale);
			fprintf( stderr, " basesFromEndLeftContig: %d\n", gapInfoArray[ilocale].basesFromEndLeftContig);
			fprintf( stderr, "basesFromEndRightContig: %d\n", gapInfoArray[ilocale].basesFromEndRightContig);
		    fprintf( stderr, "distanceDifference too large for walking (%s separation is larger)\n\n",
			  gapInfoArray[ilocale].fragSeparationScaffold > gapInfoArray[ilocale].fragSeparationLocale ? "scaff" : "locale");
		    dumpContigInfo( lchunk );
		    dumpContigInfo( rchunk );
			dumpGapInfo( lchunk, rchunk);
		    break;  // don't look at remaining locales, we've already sorted by distanceDifference
		  }
		  
		  // fprintf( stderr, "calling dumpWalkContigs\n");

		  // dumpWalkContigs( gapInfoArray[ilocale] );
		  // then compute the BAC-induced overlaps among all the contigs in current locale
		  chunksWalked = BuildLocaleOverlaps(&gapInfoArray[ilocale],
						     lchunk, rchunk, &completeOverlapPath, 
						     &rchunkNewAEnd, &rchunkNewBEnd, 
						     &olapOrientation,
						     &rchunkOverlap, gapSizeStdDevs != CONSERVATIVE_WALKING_STD_DEVS);
		  
		  if (completeOverlapPath == 1)
		  {
			rchunk_delta.mean = MIN( rchunkNewAEnd, rchunkNewBEnd) - 
			  MIN( rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean);
			if (chunksWalked->contigID == rchunk->id)
			{
			  //fprintf( stderr, "before CheckOrientation, rchunk %d: AEndOffset.mean: %f, BEndOffset.mean: %f\n",
			  //	   rchunk->id, rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean);
			  CheckOrientation( lchunk, rchunk, &rchunk_delta, olapOrientation, &rchunkOverlap);
			}
#if DEBUG_GAP_WALKER > -1		  
			fprintf( stderr, "before adjustment, rchunk %d: AEndOffset.mean: %f, BEndOffset.mean: %f\n",
					 rchunk->id, rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean);
#endif
			AdjustRightContigPos( lchunk, rchunk, rchunk_delta, scaff);

#if DEBUG_GAP_WALKER > -1		  
			fprintf( stderr, "after adjustment, rchunk %d: AEndOffset.mean: %f, BEndOffset.mean: %f\n",
					 rchunk->id, rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean);
#endif
			InsertWalkedChunks( chunksWalked, lchunk, rchunk, completeOverlapPath, 
								checkScaffold, &checkScaffoldCount );

			Force_Increasing_Variances_One_Scaffold( rchunk->scaffoldID );

			CheckCIScaffoldTLength( ScaffoldGraph, scaff);
#if DEBUG_GAP_WALKER > -1		  
			fprintf( stderr, "after insertion, rchunk %d: AEndOffset.mean: %f, BEndOffset.mean: %f\n",
					 rchunk->id, rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean);
#endif			
			CheckScaffoldOrder(scaff, ScaffoldGraph);
			
			// no need to check the other common locales
			break;
		  }
		  // to do (and below): FreeChunksWalked( chunksWalked );
		}

		// undo the temporary setting of contig positions of frags that are in surrogates
		ResetFragPositions( numCommonLocales, &gapInfoArray);

		if (completeOverlapPath == 1)
		  walksSuccessful++;
		else
		{
		  walksUnsuccessful++;
		  if (completeOverlapPath == -1)
			walksUnsuccessfulFailedOverlap++;
		  else if (completeOverlapPath == -2)
			walksUnsuccessfulFragOutOfBounds++;
		  else if (completeOverlapPath == -3)
			walksUnsuccessfulContigInertiaInScaffold++;
		  else if (completeOverlapPath == -4)
			assert(0);
		  else if (completeOverlapPath == -5)
			walksUnsuccessfulFailedOverlapStack++;
		  else if (completeOverlapPath == -6)
			walksUnsuccessfulRchunkContained++;
		  else if (completeOverlapPath == -7)
			walksUnsuccessfulFragNotSet++;
		  else if (completeOverlapPath == -8)
			walksUnsuccessfulContigInertiaOtherScaffold++;
		  else if (completeOverlapPath == -9)
			walksUnsuccessfulFragsDoNotOverlapInLocale++;
		}
		safe_free(gapInfoArray);  // this is malloced in FindCommonLocales
      }
	  currentLeftNodeId = currentRightNodeId;
	}

	StopTimerT(&singleWalkAndUpdateTime);
	
	StartTimerT(&GlobalData->WalkUpdateTimer);
	// scaffStatp->insertedChunks = Bac_Update_Scaffold_Graph(ScaffoldGraph, fill_chunks, FALSE, TRUE, FALSE, FALSE, sId);
	
	/* Make sure all variances are increasing */
	// fprintf(stderr,"**** FORCE variances after updating scaffold graph   ****\n");
	Force_Increasing_Variances_One_Scaffold( rchunk->scaffoldID );
	StopTimerT(&GlobalData->WalkUpdateTimer);
	
	// fprintf(stderr,"=== Inserted: %d chunks in Scaffold %d \n", scaffStatp->insertedChunks,sId);
	
	CleanupAScaffold(ScaffoldGraph,scaff,FALSE, NULLINDEX, FALSE);
	
	// CheckAlteredScaffolds( checkScaffold, checkScaffoldCount );

	// mark this scaffold as walked
	scaff->flags.bits.walkedAlready = TRUE;
	numScaffoldsWalked++;
	
	fprintf(gwlogfp,"*** Time for walking and updating scaffold %d = %g seconds\n",
			sId, TotalTimerT(&singleWalkAndUpdateTime, &cycles));
  }

  fprintf( stderr, "numScaffoldsWalked: %d\n", numScaffoldsWalked);

  fprintf( stderr, "Before walking, number of gaps: %d\n", numGapsBeforeWalking);
  fprintf( stderr, "Before walking, bases in gaps: %f\n", basesInGapsBeforeWalking);

  fprintf( stderr, "Walks attempted: %d\n", walksAttempted);
  fprintf( stderr, "Walks successful: %d\n", walksSuccessful);
  fprintf( stderr, "Walks unsuccessful: %d\n", walksUnsuccessful);
  // fprintf( stderr, "  Walks unsuccessful (gap too small): %d\n", walksUnsuccessfulGapTooSmall);
  fprintf( stderr, "  Walks unsuccessful (distance estimate too large): %d\n", walksUnsuccessfulDistanceDifference);
  fprintf( stderr, "  Walks unsuccessful (no common locale): %d\n", walksUnsuccessfulNoCommonLocales);
  fprintf( stderr, "  Walks unsuccessful (failed overlap - pairwise) (c1): %d\n", walksUnsuccessfulFailedOverlap);
  fprintf( stderr, "  Walks unsuccessful (fragment out of bounds) (c2): %d\n", walksUnsuccessfulFragOutOfBounds);
  fprintf( stderr, "  Walks unsuccessful (contig in this scaffold) (c3): %d\n", walksUnsuccessfulContigInertiaInScaffold);
  fprintf( stderr, "  Walks unsuccessful (failed overlap - stack) (c5): %d\n", walksUnsuccessfulFailedOverlapStack);
  fprintf( stderr, "  Walks unsuccessful (rchunk contained) (c6): %d\n", walksUnsuccessfulRchunkContained);
  fprintf( stderr, "  Walks unsuccessful (frag not set) (c7): %d\n", walksUnsuccessfulFragNotSet);
  fprintf( stderr, "  Walks unsuccessful (contig in other scaffold) (c8): %d\n", walksUnsuccessfulContigInertiaOtherScaffold);
  fprintf( stderr, "  Walks unsuccessful (frags do not overlap) (c9): %d\n", walksUnsuccessfulFragsDoNotOverlapInLocale);
  fprintf( stderr, "Writing file walkingStats.\n");
  
  fprintf( walkingStatsFile, "%d\n", numScaffoldsWalked);
  fprintf( walkingStatsFile, "%d\n", numGapsBeforeWalking);
  fprintf( walkingStatsFile, "%f\n", basesInGapsBeforeWalking);
  fprintf( walkingStatsFile, "%d\n", walksSuccessful);
  fprintf( walkingStatsFile, "%d\n", walksUnsuccessful);
  fprintf( walkingStatsFile, "%d\n", walksUnsuccessfulFailedOverlap);
  fprintf( walkingStatsFile, "%d\n", walksUnsuccessfulFragOutOfBounds);
  fprintf( walkingStatsFile, "%d\n", walksUnsuccessfulContigInertiaInScaffold);
  fprintf( walkingStatsFile, "%d\n", walksUnsuccessfulFailedOverlapStack);
  fprintf( walkingStatsFile, "%d\n", walksUnsuccessfulRchunkContained);
  fprintf( walkingStatsFile, "%d\n", walksUnsuccessfulFragNotSet);
  fprintf( walkingStatsFile, "%d\n", walksUnsuccessfulContigInertiaOtherScaffold);
  fprintf( walkingStatsFile, "%d\n", walksUnsuccessfulFragsDoNotOverlapInLocale);
  fclose( walkingStatsFile );
  
  // ComputeGapStatisitics();

  // After walking the scaffolds compute accumulative statistics
  // write them to a file and to stderr
  // compute the accumulative statistics and print them
  
  fprintf(stderr,
		  "\n------------------------------: gap walking done\n");
  
  return;
}

void Inter_Scaffold_Analysis(void)
{
  CDS_CID_t sid1, sid2;
  GapInfoT *scaffoldGapInfoArray;
  int numCommonLocales;
  NodeCGW_T *contig1, *contig2;
  int i, contig1GapEnd, contig2GapEnd;
  
  fprintf(stderr,"------------------------------: interscaffold gap walking running\n");


  for (sid1 = 0; sid1 < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid1++) 
  {
    CIScaffoldT *scaff1;
	
    scaff1 = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid1);
    // make sure the scaffold is there
    assert(scaff1 != NULL);   

    //
    // not interested in dead scaffold, not real scaffolds, or singleton
    // scaffolds
    
    if ((isDeadCIScaffoldT(scaff1)) ||
       	(scaff1->type != REAL_SCAFFOLD))
	{
	  continue;
	}

	// scaff1, AEnd
	contig1 = GetGraphNode( ScaffoldGraph->ContigGraph, scaff1->info.Scaffold.AEndCI);
	assert(contig1 != NULL);

	// find the locales that are at the apropos end of contig1
	contig1GapEnd = ( contig1->offsetAEnd.mean < contig1->offsetBEnd.mean ? A_END : B_END);
	
	for (sid2 = sid1 + 1; sid2 < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid2++) 
	{
	  CIScaffoldT *scaff2;

	  if (scaff1->flags.bits.isDead)
		continue;

	  scaff2 = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid2);
	  // make sure the scaffold is there
	  assert(scaff2 != NULL);   
	  
	  //
	  // not interested in dead scaffold, not real scaffolds, or singleton
	  // scaffolds
	  
	  if ((isDeadCIScaffoldT(scaff2)) ||
		  (scaff2->type != REAL_SCAFFOLD) ||
		  (scaff2->info.Scaffold.numElements < 2))
	  {
		continue;
	  }
	  
	  // get the contig at the AEnd
	  contig2 = GetGraphNode( ScaffoldGraph->ContigGraph, scaff2->info.Scaffold.AEndCI);
	  assert(contig2 != NULL);

	  // find the locales that are at the apropos end of contig2
	  contig2GapEnd = ( contig2->offsetAEnd.mean < contig2->offsetBEnd.mean ? A_END : B_END);
	  
	  numCommonLocales = FindCommonLocales( contig1, contig1GapEnd, contig2, contig2GapEnd, &scaffoldGapInfoArray);

	  fprintf( stderr, "The A_END of scaff %d (contig %d) and the A_END end of scaff %d (contig %d) have %d common locales\n",
			   scaff1->id, contig1->id, scaff2->id, contig2->id, numCommonLocales);

	  for (i = 0; i < numCommonLocales; i++)
	  {
		fprintf( stderr, "fragIidLeftContig: %d (contig %d), fragIidRightContig: %d (contig %d)\n",
				 GetCIFragT(ScaffoldGraph->CIFrags, (int32) scaffoldGapInfoArray[i].fragIDLeftContig)->iid,
				 scaffoldGapInfoArray[i].fragLeftContigID, 
				 GetCIFragT(ScaffoldGraph->CIFrags, (int32) scaffoldGapInfoArray[i].fragIDRightContig)->iid,
				 scaffoldGapInfoArray[i].fragRightContigID);	  
	  }

	  if (scaff1->flags.bits.isDead || scaff2->flags.bits.isDead)
		continue;

	  // get the contig at the BEnd
	  contig2 = GetGraphNode( ScaffoldGraph->ContigGraph, scaff2->info.Scaffold.BEndCI);
	  assert(contig2 != NULL);

	  // find the locales that are at the apropos end of contig2
	  contig2GapEnd = ( contig2->offsetAEnd.mean < contig2->offsetBEnd.mean ? B_END : A_END);
	  
	  numCommonLocales = FindCommonLocales( contig1, contig1GapEnd, contig2, contig2GapEnd, &scaffoldGapInfoArray);

	  fprintf( stderr, "The A_END of scaff %d (contig %d) and the B_END end of scaff %d (contig %d) have %d common locales\n",
			   scaff1->id, contig1->id, scaff2->id, contig2->id, numCommonLocales);

	  for (i = 0; i < numCommonLocales; i++)
		fprintf( stderr, "fragIidLeftContig: %d (contig %d), fragIidRightContig: %d (contig %d)\n",
				 GetCIFragT(ScaffoldGraph->CIFrags, (int32) scaffoldGapInfoArray[i].fragIDLeftContig)->iid,
				 scaffoldGapInfoArray[i].fragLeftContigID, 
				 GetCIFragT(ScaffoldGraph->CIFrags, (int32) scaffoldGapInfoArray[i].fragIDRightContig)->iid,
				 scaffoldGapInfoArray[i].fragRightContigID);
	}

	// scaff1, BEnd
	contig1 = GetGraphNode( ScaffoldGraph->ContigGraph, scaff1->info.Scaffold.BEndCI);
	assert(contig1 != NULL);

	// find the locales that are at the apropos end of contig1
	contig1GapEnd = ( contig1->offsetAEnd.mean < contig1->offsetBEnd.mean ? B_END : A_END);
	
	for (sid2 = sid1 + 1; sid2 < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid2++) 
	{
	  CIScaffoldT *scaff2;

	  scaff2 = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid2);
	  // make sure the scaffold is there
	  assert(scaff2 != NULL);   
	  
	  //
	  // not interested in dead scaffold, not real scaffolds, or singleton
	  // scaffolds
	  
	  if ((isDeadCIScaffoldT(scaff2)) ||
		  (scaff2->type != REAL_SCAFFOLD) ||
		  (scaff2->info.Scaffold.numElements < 2))
	  {
		continue;
	  }
	  
	  // get the contig at the AEnd
	  contig2 = GetGraphNode( ScaffoldGraph->ContigGraph, scaff2->info.Scaffold.AEndCI);
	  assert(contig2 != NULL);

	  // find the locales that are at the apropos end of contig2
	  contig2GapEnd = ( contig2->offsetAEnd.mean < contig2->offsetBEnd.mean ? A_END : B_END);
	  
	  numCommonLocales = FindCommonLocales( contig1, contig1GapEnd, contig2, contig2GapEnd, &scaffoldGapInfoArray);

	  fprintf( stderr, "The B_END of scaff %d (contig %d) and the A_END end of scaff %d (contig %d) have %d common locales\n",
			   scaff1->id, contig1->id, scaff2->id, contig2->id, numCommonLocales);

	  for (i = 0; i < numCommonLocales; i++)
		fprintf( stderr, "fragIidLeftContig: %d (contig %d), fragIidRightContig: %d (contig %d)\n",
				 GetCIFragT(ScaffoldGraph->CIFrags, (int32) scaffoldGapInfoArray[i].fragIDLeftContig)->iid,
				 scaffoldGapInfoArray[i].fragLeftContigID, 
				 GetCIFragT(ScaffoldGraph->CIFrags, (int32) scaffoldGapInfoArray[i].fragIDRightContig)->iid,
				 scaffoldGapInfoArray[i].fragRightContigID);



	  // get the contig at the BEnd
	  contig2 = GetGraphNode( ScaffoldGraph->ContigGraph, scaff2->info.Scaffold.BEndCI);
	  assert(contig2 != NULL);

	  // find the locales that are at the apropos end of contig2
	  contig2GapEnd = ( contig2->offsetAEnd.mean < contig2->offsetBEnd.mean ? B_END : A_END);
	  
	  numCommonLocales = FindCommonLocales( contig1, contig1GapEnd, contig2, contig2GapEnd, &scaffoldGapInfoArray);

	  fprintf( stderr, "The B_END of scaff %d (contig %d) and the B_END end of scaff %d (contig %d) have %d common locales\n",
			   scaff1->id, contig1->id, scaff2->id, contig2->id, numCommonLocales);

	  for (i = 0; i < numCommonLocales; i++)
		fprintf( stderr, "fragIidLeftContig: %d (contig %d), fragIidRightContig: %d (contig %d)\n",
				 GetCIFragT(ScaffoldGraph->CIFrags, (int32) scaffoldGapInfoArray[i].fragIDLeftContig)->iid,
				 scaffoldGapInfoArray[i].fragLeftContigID, 
				 GetCIFragT(ScaffoldGraph->CIFrags, (int32) scaffoldGapInfoArray[i].fragIDRightContig)->iid,
				 scaffoldGapInfoArray[i].fragRightContigID);
	}
  }

  // temp hack
  if (0)
  {
	// really need the contigs relative scaffold positions based on frags locale positions
	double gapSize = 3 * 275.0;  
    CIScaffoldT *scaff1, *scaff2;
	
    scaff1 = GetGraphNode(ScaffoldGraph->ScaffoldGraph, 638);
    assert(scaff1 != NULL);   

    scaff2 = GetGraphNode(ScaffoldGraph->ScaffoldGraph, 698);
    assert(scaff2 != NULL);   

	fprintf( stderr, "merging scaffolds %d and %d\n", scaff1->id, scaff2->id);
	
	walkScaffolds( scaff1, 1, scaff2, 0, gapSize);

    scaff1 = GetGraphNode(ScaffoldGraph->ScaffoldGraph, 638);
    assert(scaff1 != NULL);   

    scaff2 = GetGraphNode(ScaffoldGraph->ScaffoldGraph, 698);
    assert(scaff2 != NULL);   

	scaff1->flags.bits.isDead = TRUE;
	scaff1->info.Scaffold.numElements = 0;
	scaff2->flags.bits.isDead = TRUE;	
	scaff2->info.Scaffold.numElements = 0;
  }
}


// Outline for Inter_Scaffold_Walking (2/27/01 MJF)
// 1. Step through all scaffolds finding locales within 2000bp of end
// 2. Store this information by BAC, sorted by iid
// 3. Prioritize walking order (by frag separation in locale?)
// 4. Attempt walks by using walkScaffolds

int Inter_Scaffold_Walking(void)
{
  CDS_CID_t sid;
  int i;
  int *AEnds, *BEnds;
  ScaffoldEndT scaffoldEnds[100000];
  int scaffoldEndCount = 0;
  int interScaffoldWalkAttempts = 0, interScaffoldWalkSuccesses = 0;
  
  fprintf(stderr,"------------------------------: interscaffold gap walking running\n");

  // step through all scaffolds finding ends that have BAC frags within ISGW_BASES_FROM_END_CUTOFF
  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++) 
  {
    CIScaffoldT *scaff;
	ContigT *contig;
	int contigEnd, numLocales;
	localeInfoT	*currLocaleInfo, *localeInfoContig = NULL;

    scaff = GetGraphNode( ScaffoldGraph->ScaffoldGraph, sid);
    // make sure the scaffold is there
    assert(scaff != NULL);   

    //
    // not interested in dead scaffold,  or unreal scaffolds
    // what about singletons? 		  (scaff->info.Scaffold.numElements < 2)
    if ( isDeadCIScaffoldT( scaff) || (scaff->type != REAL_SCAFFOLD))
	  continue;

	// get rid of small contigs on end that may prevent walking
	// removeSmallContigs( scaff, TRUE);
	trimScaffoldEnds( scaff );
	
    scaff = GetGraphNode( ScaffoldGraph->ScaffoldGraph, sid);

	// temp hack

	//if (scaff->bpLength.mean < 1000000)
	//continue;



	
	// scaff, AEnd
	contig = GetGraphNode( ScaffoldGraph->ContigGraph, scaff->info.Scaffold.AEndCI);
	assert (contig != NULL);
	
	// check to make sure this function is more general than its name, ie, it should be contains_bac
	if ( contains_fbac( contig ))
	{
	  // find the locales that are at the apropos end of contig
	  contigEnd = ( contig->offsetAEnd.mean < contig->offsetBEnd.mean ? A_END : B_END);
#define ISGW_BASES_FROM_END_CUTOFF 2000
	  numLocales = getLocalesInNode( contig, &localeInfoContig, contigEnd, ISGW_BASES_FROM_END_CUTOFF);	
	  currLocaleInfo = localeInfoContig;
	  
	  if (contig->id == 74288)
		fprintf( stderr, "for contig %d, currLocaleInfo->inSurrogate: %d\n", contig->id, currLocaleInfo->inSurrogate);

	  if (numLocales > 0)
	  {
		int numInserted = 0;
		
		// now insert info re the scaffold in scaffoldEnds
		for ( i = 0; i < numLocales; i++)
		{
		  if (currLocaleInfo->inSurrogate)  // don't use frags in surrogates already to base walks on
			continue;                       // it also messes things up since then frags can occur more than once
		  
		  if (1) // (currLocaleInfo->numLocaleFragments > 5 )
		  {
			scaffoldEnds[scaffoldEndCount].localeNumber = currLocaleInfo->localeNumber;
			scaffoldEnds[scaffoldEndCount].leftScaffID = scaff->id;
			scaffoldEnds[scaffoldEndCount].leftScaffEnd = A_END; // contigEnd;
			scaffoldEnds[scaffoldEndCount].leftExtremalFragID = currLocaleInfo->extremalFragID;
			scaffoldEnds[scaffoldEndCount].leftExtremalFragIid = currLocaleInfo->extremalFragIid;
			scaffoldEnds[scaffoldEndCount].leftContigID = currLocaleInfo->contigID;
			scaffoldEnds[scaffoldEndCount].leftBasesFromEnd = currLocaleInfo->basesFromEnd;
			scaffoldEnds[scaffoldEndCount].leftBacOffset = currLocaleInfo->BacOffset;
			scaffoldEnds[scaffoldEndCount].leftIidTrend = currLocaleInfo->iidTrend;
			scaffoldEndCount++;
			numInserted++;
		  }
		  currLocaleInfo = currLocaleInfo->next;		  
		}
		
		// clean up localeInfoContig
		currLocaleInfo = localeInfoContig->next;
		while (currLocaleInfo != NULL)
		{
		  safe_free (localeInfoContig);
		  localeInfoContig = currLocaleInfo;
		  currLocaleInfo = localeInfoContig->next;
		}
		safe_free (localeInfoContig);
		
		for ( i = 0; i < numInserted; i++)
		{
		  fprintf( stderr, "scaffoldEnds[%d].localeNumber = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].localeNumber);
		  fprintf( stderr, "scaffoldEnds[%d].leftScaffID = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftScaffID);
		  fprintf( stderr, "scaffoldEnds[%d].leftContigID = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftContigID);
		  fprintf( stderr, "scaffoldEnds[%d].leftExtremalFragID = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftExtremalFragID);
		  fprintf( stderr, "scaffoldEnds[%d].leftExtremalFragIid = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftExtremalFragIid);
		  fprintf( stderr, "scaffoldEnds[%d].leftInSurrogate = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftInSurrogate);
		  fprintf( stderr, "scaffoldEnds[%d].leftScaffEnd = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftScaffEnd);
		  fprintf( stderr, "\n");
		}
	  }
	}
	
	// scaff, BEnd
	contig = GetGraphNode( ScaffoldGraph->ContigGraph, scaff->info.Scaffold.BEndCI);
	assert (contig != NULL);
	
	// check to make sure this function is more general than its name, ie, it should be contains_bac
	if ( contains_fbac( contig ))
	{
	  // find the locales that are at the apropos end of contig
	  contigEnd = ( contig->offsetAEnd.mean < contig->offsetBEnd.mean ? B_END : A_END);

	  numLocales = getLocalesInNode( contig, &localeInfoContig, contigEnd, ISGW_BASES_FROM_END_CUTOFF);	
	  currLocaleInfo = localeInfoContig;

	  if (numLocales > 0)
	  {
		int numInserted = 0;
		
		// now insert info re the scaffold in scaffoldEnds
		for ( i = 0; i < numLocales; i++)
		{
		  if (currLocaleInfo->inSurrogate)  // don't use frags in surrogates already to base walks on
			continue;                       // it also messes things up since then frags can occur more than once
		  
		  if (1) // (currLocaleInfo->numLocaleFragments > 5 )
		  {
			scaffoldEnds[scaffoldEndCount].localeNumber = currLocaleInfo->localeNumber;
			scaffoldEnds[scaffoldEndCount].leftScaffID = scaff->id;
			scaffoldEnds[scaffoldEndCount].leftScaffEnd = B_END; // contigEnd;
			scaffoldEnds[scaffoldEndCount].leftExtremalFragID = currLocaleInfo->extremalFragID;
			scaffoldEnds[scaffoldEndCount].leftExtremalFragIid = currLocaleInfo->extremalFragIid;
			scaffoldEnds[scaffoldEndCount].leftContigID = currLocaleInfo->contigID;
			scaffoldEnds[scaffoldEndCount].leftBasesFromEnd = currLocaleInfo->basesFromEnd;
			scaffoldEnds[scaffoldEndCount].leftBacOffset = currLocaleInfo->BacOffset;
			scaffoldEnds[scaffoldEndCount].leftIidTrend = currLocaleInfo->iidTrend;
			scaffoldEndCount++;
			numInserted++;
		  } 
		  currLocaleInfo = currLocaleInfo->next;
		}
		
		// clean up localeInfoContig
		currLocaleInfo = localeInfoContig->next;
		while (currLocaleInfo != NULL)
		{
		  safe_free (localeInfoContig);
		  localeInfoContig = currLocaleInfo;
		  currLocaleInfo = localeInfoContig->next;
		}
		safe_free (localeInfoContig);
		
		for ( i = 0; i < numInserted; i++)
		{
		  fprintf( stderr, "scaffoldEnds[%d].localeNumber = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].localeNumber);
		  fprintf( stderr, "scaffoldEnds[%d].leftScaffID = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftScaffID);
		  fprintf( stderr, "scaffoldEnds[%d].leftContigID = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftContigID);
		  fprintf( stderr, "scaffoldEnds[%d].leftExtremalFragID = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftExtremalFragID);
		  fprintf( stderr, "scaffoldEnds[%d].leftExtremalFragIid = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftExtremalFragIid);
		  fprintf( stderr, "scaffoldEnds[%d].leftInSurrogate = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftInSurrogate);
		  fprintf( stderr, "scaffoldEnds[%d].leftScaffEnd = %d\n",
				   i, scaffoldEnds[scaffoldEndCount - numInserted + i].leftScaffEnd);
		  fprintf( stderr, "\n");
		}
	  }
	}
  }

  // now we need to sort find the gaps a particular BAC is capable of walking
  // so sort scaffoldEnds by localeNumber and BacOffset
  qsort( scaffoldEnds, scaffoldEndCount, sizeof(ScaffoldEndT), Sort_ScaffoldEnds_By_Locale); 

  fprintf( stderr, "after locale sorting\n");
  fprintf( stderr, "locale\t fragID\t fragIid\t contigID\t BacOffset\n");
  for ( i = 0; i < scaffoldEndCount; i++)
  {
	fprintf( stderr, "%d\t %d\t %d\t\t %d\t\t %d\n",
			 scaffoldEnds[i].localeNumber,
			 scaffoldEnds[i].leftExtremalFragID,
			 scaffoldEnds[i].leftExtremalFragIid,
			 scaffoldEnds[i].leftContigID,
			 scaffoldEnds[i].leftBacOffset);
  }

  // step through scaffoldEnds collecting gaps (leftContig, rightContig) that are walking candidates
  for ( i = 0; i < scaffoldEndCount - 1; i++)
  {
	if (scaffoldEnds[i].localeNumber == scaffoldEnds[i+1].localeNumber &&  // stay in locale
		scaffoldEnds[i].leftScaffID != scaffoldEnds[i+1].leftScaffID)      // scaffolds must differ
	{
	  scaffoldEnds[i].rightScaffID = scaffoldEnds[i+1].leftScaffID;
	  scaffoldEnds[i].rightScaffEnd = scaffoldEnds[i+1].leftScaffEnd;
	  scaffoldEnds[i].rightContigID = scaffoldEnds[i+1].leftContigID;
	  scaffoldEnds[i].rightExtremalFragID = scaffoldEnds[i+1].leftExtremalFragID;
	  scaffoldEnds[i].rightExtremalFragIid = scaffoldEnds[i+1].leftExtremalFragIid;
	  scaffoldEnds[i].rightBasesFromEnd = scaffoldEnds[i+1].leftBasesFromEnd;
	  scaffoldEnds[i].rightInSurrogate = scaffoldEnds[i+1].leftInSurrogate;
	  scaffoldEnds[i].rightBacOffset = scaffoldEnds[i+1].leftBacOffset;
	  scaffoldEnds[i].rightIidTrend = scaffoldEnds[i+1].leftIidTrend;
	}
	else
	{
	  scaffoldEnds[i].rightScaffID = NULLINDEX;
	  scaffoldEnds[i].rightScaffEnd = NULLINDEX;
	  scaffoldEnds[i].rightContigID = NULLINDEX;
	  scaffoldEnds[i].rightExtremalFragID = NULLINDEX;
	  scaffoldEnds[i].rightExtremalFragIid = NULLINDEX;
	  scaffoldEnds[i].rightBasesFromEnd = NULLINDEX;
	  scaffoldEnds[i].rightInSurrogate = NULLINDEX;
	  scaffoldEnds[i].rightBacOffset = NULLINDEX;
	}
  }

  // take care of last contig
  {
	scaffoldEnds[i].rightScaffID = NULLINDEX;
	scaffoldEnds[i].rightScaffEnd = NULLINDEX;
	scaffoldEnds[i].rightContigID = NULLINDEX;
	scaffoldEnds[i].rightExtremalFragID = NULLINDEX;
	scaffoldEnds[i].rightExtremalFragIid = NULLINDEX;
	scaffoldEnds[i].rightBasesFromEnd = NULLINDEX;
	scaffoldEnds[i].rightInSurrogate = NULLINDEX;
	scaffoldEnds[i].rightBacOffset = NULLINDEX;
  }

  fprintf( stderr, "before gap sorting\n");
  fprintf( stderr, "locale\t fragID\t fragIid  lcontigID (end)\t rcontigID (end)\n");
  for ( i = 0; i < scaffoldEndCount; i++)
  {
	fprintf( stderr, "%d\t %d\t %d\t  %d (%c)\t\t\t %d (%c)\n",
			 scaffoldEnds[i].localeNumber,
			 scaffoldEnds[i].leftExtremalFragID,
			 scaffoldEnds[i].leftExtremalFragIid,
			 scaffoldEnds[i].leftContigID,
			 scaffoldEnds[i].leftScaffEnd == A_END ? 'A' : 'B',
			 scaffoldEnds[i].rightContigID,
			 scaffoldEnds[i].rightScaffEnd == A_END ? 'A' : 'B');
  }

  // now we need to sort by gaps (to determine how many BACs can walk each gap)
  qsort( scaffoldEnds, scaffoldEndCount, sizeof(ScaffoldEndT), Sort_ScaffoldEnds_By_Gap); 

  fprintf( stderr, "after gap sorting\n");
  fprintf( stderr, "locale\t fragID\t fragIid\t lcontigID (end)\t rcontigID (end)\n");
  for ( i = 0; i < scaffoldEndCount; i++)
  {
	fprintf( stderr, "%d\t %d\t %d\t %d (%c)\t\t\t %d (%c)\n",
			 scaffoldEnds[i].localeNumber,
			 scaffoldEnds[i].leftExtremalFragID,
			 scaffoldEnds[i].leftExtremalFragIid,
			 scaffoldEnds[i].leftContigID,
			 scaffoldEnds[i].leftScaffEnd == A_END ? 'A' : 'B',
			 scaffoldEnds[i].rightContigID,
			 scaffoldEnds[i].rightScaffEnd == A_END ? 'A' : 'B');
  }

  // now determine if there are any pairs of contigs spanned by different BACs
  for ( i = 0; i < scaffoldEndCount; i++)
  {
	int j, duplicateCount = 1;
	
	if (scaffoldEnds[i].rightContigID == NULLINDEX)
	  continue;
	
	while (scaffoldEnds[i].leftContigID == scaffoldEnds[i + duplicateCount].leftContigID &&
		   scaffoldEnds[i].rightContigID == scaffoldEnds[i + duplicateCount].rightContigID)
	{
	  duplicateCount++;
	}
	
	for ( j = 0; j < duplicateCount; j++)
	  scaffoldEnds[i + j].spanningBACs = duplicateCount;
	
	i += duplicateCount - 1;
  }  

  // now we need to sort by spanning BACs across gaps
  qsort( scaffoldEnds, scaffoldEndCount, sizeof(ScaffoldEndT), Sort_Gaps_By_BAC_Count); 

  fprintf( stderr, "after spanning BAC sorting\n");
  fprintf( stderr, "locale\t lfragID  lfragIid\t lcontigID (end)\t  rfrag  rfragIid\t rcontigID (end)\t spanningBACs\n");
  for ( i = 0; i < scaffoldEndCount; i++)
  {
	fprintf( stderr, "%d\t %d\t %d\t %d (%c)\t\t %d\t %d\t %d (%c)\t\t %d\n",
			 scaffoldEnds[i].localeNumber,
			 scaffoldEnds[i].leftExtremalFragID,
			 scaffoldEnds[i].leftExtremalFragIid,
			 scaffoldEnds[i].leftContigID,
			 scaffoldEnds[i].leftScaffEnd == A_END ? 'A' : 'B',
			 scaffoldEnds[i].rightExtremalFragID,
			 scaffoldEnds[i].rightExtremalFragIid,
			 scaffoldEnds[i].rightContigID,
			 scaffoldEnds[i].rightScaffEnd == A_END ? 'A' : 'B',
			 scaffoldEnds[i].spanningBACs);
  }  

  // now make an array of all scaffold ends we're interested in
  // mark the ones that are to be involved in a walk as indicated by scaffoldEnds
  AEnds = (int *) safe_calloc( GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph), sizeof(int));
  BEnds = (int *) safe_calloc( GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph), sizeof(int));
  
  for ( i = 0; i < scaffoldEndCount; i++)
  {
	int leftEndUsed, rightEndUsed;

	fprintf( stderr, "scaffoldEnds[%d].leftScaffID: %d, scaffoldEnds[%d].rightScaffID: %d\n",
			 i, scaffoldEnds[i].leftScaffID, i, scaffoldEnds[i].rightScaffID);
	
	leftEndUsed = TRUE;
	rightEndUsed = TRUE;
	
	if ( scaffoldEnds[i].leftContigID == NULLINDEX || scaffoldEnds[i].rightContigID == NULLINDEX )
	{
	  scaffoldEnds[i].onWalkList = FALSE;
	  continue;
	}	

	// see if the left end is marked
	if (scaffoldEnds[i].leftScaffEnd == A_END)
	{
	  if ( AEnds[ scaffoldEnds[i].leftScaffID ] == FALSE) leftEndUsed = FALSE;
	}
	else
	{
	  if ( BEnds[ scaffoldEnds[i].leftScaffID ] == FALSE) leftEndUsed = FALSE;
	}

	// see if the right end is marked
	if (scaffoldEnds[i].rightScaffEnd == A_END)
	{
	  if ( AEnds[ scaffoldEnds[i].rightScaffID ] == FALSE) rightEndUsed = FALSE;
	}
	else
	{
	  if ( BEnds[ scaffoldEnds[i].rightScaffID ] == FALSE) rightEndUsed = FALSE;
	}

	// if either is already used continue on
	if ( leftEndUsed || rightEndUsed)
	{
	  scaffoldEnds[i].onWalkList = FALSE;
	  continue;
	}
	
	fprintf( stderr, "scaffoldEnds[%d].leftIidTrend: %d, scaffoldEnds[%d].rightIidTrend: %d\n",
			 i, scaffoldEnds[i].leftIidTrend, i, scaffoldEnds[i].rightIidTrend);

	// see if the iid trend on the ends is compatible, because of our sorting
	// the lower numbered frag iid is always the leftExtremalFragIid
	if (scaffoldEnds[i].leftIidTrend != 1 || scaffoldEnds[i].rightIidTrend != -1)  
	{
	  scaffoldEnds[i].onWalkList = FALSE;
	  continue;
	}
	
	fprintf( stderr, "marking scaffoldEnds[%d].leftScaffID: %d, scaffoldEnds[%d].rightScaffID: %d as used\n",
			 i, scaffoldEnds[i].leftScaffID, i, scaffoldEnds[i].rightScaffID);

	// mark both ends as used
	if (scaffoldEnds[i].leftScaffEnd == A_END)
	  AEnds[ scaffoldEnds[i].leftScaffID ] = TRUE;
	else
	  BEnds[ scaffoldEnds[i].leftScaffID ] = TRUE;

	if (scaffoldEnds[i].rightScaffEnd == A_END)
	  AEnds[ scaffoldEnds[i].rightScaffID ] = TRUE;
	else
	  BEnds[ scaffoldEnds[i].rightScaffID ] = TRUE;

	// mark the scaffoldEndT as one we should walk
	scaffoldEnds[i].onWalkList = TRUE;
  }  

  fprintf( stderr, "after marking scaffolds for merging\n");
  fprintf( stderr, "locale\t lfragID lfgIid\t lcontigID (end)\t rfragID rfgIid\t rcontigID (end)\t spanningBACs\t onWalkList\n");
  for ( i = 0; i < scaffoldEndCount; i++)
  {
	fprintf( stderr, "%d\t %d\t %d\t %d (%c)\t\t %d\t %d\t %d (%c)\t\t %d\t\t %d\n",
			 scaffoldEnds[i].localeNumber,
			 scaffoldEnds[i].leftExtremalFragID,
			 scaffoldEnds[i].leftExtremalFragIid,
			 scaffoldEnds[i].leftContigID,
			 scaffoldEnds[i].leftScaffEnd == A_END ? 'A' : 'B',
			 scaffoldEnds[i].rightExtremalFragID,
			 scaffoldEnds[i].rightExtremalFragIid,
			 scaffoldEnds[i].rightContigID,
			 scaffoldEnds[i].rightScaffEnd == A_END ? 'A' : 'B',
			 scaffoldEnds[i].spanningBACs,
			 scaffoldEnds[i].onWalkList);
  }  

  // now step though scaffoldEnds, trying to walk between scaffolds as per onWalkList
  for ( i = 0; i < scaffoldEndCount; i++)
  {
	CIFragT *leftExtremalFrag, *rightExtremalFrag;
	ContigT *leftContig, *rightContig;
	CIScaffoldT *leftScaffold, *rightScaffold;
	double gapSize;
	int invertLeftScaff, invertRightScaff, separationInLocale;
	 
	fprintf( stderr, "\n\n********* ********* ***********\n looking at walk from contig %d to %d\n",
			 scaffoldEnds[i].leftContigID, scaffoldEnds[i].rightContigID);

	if (!scaffoldEnds[i].onWalkList)
	{
	  fprintf( stderr, "skipping walk from contig %d to %d - not on walkList\n",
			   scaffoldEnds[i].leftContigID, scaffoldEnds[i].rightContigID);
	  continue;
	}	

	// grab the current contig & scaffold of the left extremal frag
	leftContig = GetGraphNode( ScaffoldGraph->ContigGraph, scaffoldEnds[i].leftContigID);
	if (leftContig->flags.bits.isDead)  // contig was used in a previous walk
	{
	  fprintf( stderr, "skipping walk from contig %d to %d - left contig dead\n",
			   scaffoldEnds[i].leftContigID, scaffoldEnds[i].rightContigID);
	  continue;	
	}
	
	leftScaffold = GetGraphNode( ScaffoldGraph->ScaffoldGraph, leftContig->scaffoldID);
	if (leftScaffold == NULL)
	{
	  fprintf( stderr, "skipping walk from contig %d to %d - left scaffold null\n",
			   scaffoldEnds[i].leftContigID, scaffoldEnds[i].rightContigID);
	  continue;
	}

	// grab the current contig & scaffold of the right extremal frag
	rightContig = GetGraphNode( ScaffoldGraph->ContigGraph, scaffoldEnds[i].rightContigID);
	if (rightContig->flags.bits.isDead)  // contig was used in a previous walk
	{
	  fprintf( stderr, "skipping walk from contig %d to %d - right contig dead\n",
			   scaffoldEnds[i].leftContigID, scaffoldEnds[i].rightContigID);
	  continue;
	}
	
	rightScaffold = GetGraphNode( ScaffoldGraph->ScaffoldGraph, rightContig->scaffoldID);
	if (rightScaffold == NULL)
	{
	  fprintf( stderr, "skipping walk from contig %d to %d - left scaffold null\n",
			   scaffoldEnds[i].leftContigID, scaffoldEnds[i].rightContigID);
	  continue;
	}

	// grab the extremal frags so we know how far apart they are in locale
	leftExtremalFrag = GetCIFragT( ScaffoldGraph->CIFrags, (int32) scaffoldEnds[i].leftExtremalFragID );
	AssertPtr( leftExtremalFrag );
	
	rightExtremalFrag = GetCIFragT( ScaffoldGraph->CIFrags, (int32) scaffoldEnds[i].rightExtremalFragID );
	AssertPtr( rightExtremalFrag );

	separationInLocale = rightExtremalFrag->localePos.bgn - leftExtremalFrag->localePos.end;

	gapSize = separationInLocale - scaffoldEnds[i].leftBasesFromEnd - scaffoldEnds[i].rightBasesFromEnd;
	fprintf( stderr, "figured gapSize of %f using leftExtremalFrag->iid %d and rightExtremalFrag->iid %d (sep: %d)\n",
			 gapSize, leftExtremalFrag->iid, rightExtremalFrag->iid, separationInLocale);

	// check to make sure the left contig is still on the end of its scaffold
	if (leftContig->id == leftScaffold->info.Scaffold.AEndCI ||
		leftContig->id == leftScaffold->info.Scaffold.BEndCI)
	{
	  int leftFragLeftEnd, leftFragRightEnd, leftFragScaffoldOrientation;
	  
	  GetFragmentPositionInScaffold(leftExtremalFrag, &leftFragLeftEnd, &leftFragRightEnd, &leftFragScaffoldOrientation);
	  if (leftFragScaffoldOrientation == 0)  // 0 means it is (5p->3p)
		// fprintf( stderr, "1 leftExtremalFrag %d will not be 5p->3p in scaffold\n", leftExtremalFrag->iid);
		invertLeftScaff = FALSE;
	  else
		invertLeftScaff = TRUE;
	}	
	else
	  continue;

	// check to make sure the right contig is still on the end of its scaffold
	if (rightContig->id == rightScaffold->info.Scaffold.AEndCI ||
		rightContig->id == rightScaffold->info.Scaffold.BEndCI)
	{
	  int rightFragLeftEnd, rightFragRightEnd, rightFragScaffoldOrientation;
	  
	  GetFragmentPositionInScaffold(rightExtremalFrag, &rightFragLeftEnd, &rightFragRightEnd, &rightFragScaffoldOrientation);
	  if (rightFragScaffoldOrientation == 0)  // 0 means it is (5p->3p)
		// fprintf( stderr, "1 rightExtremalFrag %d will not be 5p->3p in scaffold\n", rightExtremalFrag->iid);
		invertRightScaff = FALSE;
	  else
		invertRightScaff = TRUE;
	}
	else
	  continue;	

	// check to make sure we're not in the same scaffold (in case the other ends were used in a merge)
	if (leftContig->scaffoldID == rightContig->scaffoldID)
	  continue;

	if (gapSize < 100000)
	{
	  ScaffoldInfoT oldLeftScaffoldContigs, oldRightScaffoldContigs;
	  int fragDiff;
	  
	  SaveContigInformation( leftScaffold, &oldLeftScaffoldContigs);
	  SaveContigInformation( rightScaffold, &oldRightScaffoldContigs);	  

	  fprintf( stderr, 
			   "try inter-scaff walk from scaffold %d (contig %d, invert %d) to scaffold %d (contig %d, invert %d)\n",
			   leftScaffold->id, leftContig->id, invertLeftScaff,
			   rightScaffold->id, rightContig->id, invertRightScaff);
	  interScaffoldWalkAttempts++;
	  walkScaffolds( leftScaffold, invertLeftScaff, rightScaffold, invertRightScaff, gapSize);	

	  // check to see if there are contigs between the left and right contigs, or they overlap by 20
	  // if not, assume walk has failed, restore contigs to their original scaffold
	  // if success, continue with steps below

	  // now grab the new scaffold via the leftContig
	  leftContig = GetGraphNode( ScaffoldGraph->ContigGraph, scaffoldEnds[i].leftContigID);
	  rightContig = GetGraphNode( ScaffoldGraph->ContigGraph, scaffoldEnds[i].rightContigID);

	  fragDiff = scaffoldEnds[i].rightExtremalFragIid - scaffoldEnds[i].leftExtremalFragIid;
	  
	  // walk failure indication is when contigs are adjacent in scaffold as inserted
	  // but no intervening contigs are present
	  // exception is when fragIids are consecutive - walk succeeds by fiat with 20bp overlap
	  if (leftContig->BEndNext == rightContig->id && fragDiff != 1)
	  {
		CIScaffoldT *newScaffold;
		
		fprintf( stderr, "contigs %d and %d are adjacent in scaffold %d\n",
				 leftContig->id, rightContig->id, leftContig->scaffoldID);

		newScaffold = GetGraphNode( ScaffoldGraph->ScaffoldGraph, leftContig->scaffoldID);
		newScaffold->flags.bits.isDead = TRUE;
		newScaffold->info.Scaffold.numElements = 0;
		

		newScaffold = CreateNewScaffold();
		RestoreContigInformation( newScaffold, &oldLeftScaffoldContigs);

		newScaffold = CreateNewScaffold();
		RestoreContigInformation( newScaffold, &oldRightScaffoldContigs);		
	  }
	  else
		interScaffoldWalkSuccesses++;
	  
	  FreeContigInformation( &oldLeftScaffoldContigs );
	  FreeContigInformation( &oldRightScaffoldContigs );	
	}
	else
	  fprintf( stderr,
                   "skipping walk from contig %d to %d - gap too big\n",
                   leftContig->id, rightContig->id);
  }
  fprintf(stderr, " interScaffoldWalkAttempts: %d\n", interScaffoldWalkAttempts);
  if ( interScaffoldWalkAttempts > 0 )
	fprintf(stderr, "interScaffoldWalkSuccesses: %d (success rate: %f)\n", interScaffoldWalkSuccesses,
			(100.0 * interScaffoldWalkSuccesses) / interScaffoldWalkAttempts);
  fprintf(stderr, "------------------------------: interscaffold gap walking ending\n");
  
  return interScaffoldWalkSuccesses;
}

CIScaffoldT* CreateNewScaffold(void)
{
  CIScaffoldT* newScaffold;
  LengthT NullLength = {0.0, 0.0};

  // create & init a new scaffold
  newScaffold = (CIScaffoldT *) safe_malloc( sizeof (CIScaffoldT));
  InitializeScaffold( newScaffold, REAL_SCAFFOLD);
  newScaffold->info.Scaffold.AEndCI = NULLINDEX;
  newScaffold->info.Scaffold.BEndCI = NULLINDEX;
  newScaffold->info.Scaffold.numElements = 0;
  newScaffold->edgeHead = NULLINDEX;
  newScaffold->bpLength = NullLength; 
  newScaffold->id = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
  newScaffold->flags.bits.isDead = FALSE;
  newScaffold->aEndCoord = newScaffold->bEndCoord = -1;
  newScaffold->numEssentialA = newScaffold->numEssentialB = 0;
  newScaffold->essentialEdgeB = newScaffold->essentialEdgeA = NULLINDEX;
  AppendGraphNode( ScaffoldGraph->ScaffoldGraph, newScaffold);

  return newScaffold;
}




#define BASE_CUTOFF 50

void ComputeGapStatisitics(void)
{
  int sc, numGaps = 0, numGapsInScaffoldsWithBACFrags = 0;
  double basesInThisGap, basesInGaps = 0.0, basesInGapsInScaffoldsWithBACFrags = 0.0;
  int scaffoldContainsBACFrags;

  //
  // scan all the scaffolds
  //  
  for (sc = 0; sc < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sc++) 
  {
	CIScaffoldTIterator
	  CIsTemp;
	CIScaffoldT
	  * scaff;

    scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sc);
    // make sure the scaffold is there
    assert(scaff != NULL);
    
    //
    // not interested in dead scaffold, not real scaffolds, or singleton
    // scaffolds
    
    if ((isDeadCIScaffoldT(scaff)) ||
       	(scaff->type != REAL_SCAFFOLD) ||
		(scaff->info.Scaffold.numElements < 2))
	{
	  continue;
	}

    // fprintf(stderr,"=== Walking scaffold %d, size %d\n",sId,sortedScaffolds[sc].scaff_length);

	scaffoldContainsBACFrags = 0;
	InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE,
							FALSE, &CIsTemp);
	while (NextCIScaffoldTIterator(&CIsTemp))
	{
	  ChunkInstanceT
		* lchunkTemp,
		* rchunkTemp;
	  
	  // not walking off of scaffolds currently
	  if (CIsTemp.next == -1)
		break;
	  
	  //
	  // find the chunks in the gap by walking between the chunk <CIs.curr>
	  // and the chunk <CIs.next>
	  //
	  
	  lchunkTemp = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.curr);
	  rchunkTemp = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.next);
	  
	  assert(lchunkTemp != NULL);
	  assert(rchunkTemp != NULL);

	  if ( contains_fbac( lchunkTemp ) || contains_fbac( rchunkTemp ))
		scaffoldContainsBACFrags = 1;
	}

	InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE,
							FALSE, &CIsTemp);
	while (NextCIScaffoldTIterator(&CIsTemp))
	{
	  ChunkInstanceT
		* lchunkTemp,
		* rchunkTemp;
	  
	  // not walking off of scaffolds currently
	  if (CIsTemp.next == -1)
		break;
	  
	  //
	  // find the chunks in the gap by walking between the chunk <CIs.curr>
	  // and the chunk <CIs.next>
	  //
	  
	  lchunkTemp = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.curr);
	  rchunkTemp = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.next);
	  
	  assert(lchunkTemp != NULL);
	  assert(rchunkTemp != NULL);
	  
	  basesInThisGap = MIN (rchunkTemp->offsetAEnd.mean, rchunkTemp->offsetBEnd.mean) - 
		MAX (lchunkTemp->offsetAEnd.mean, lchunkTemp->offsetBEnd.mean);
	  // don't include "gaps" where contigs overlap
	  if (basesInThisGap > BASE_CUTOFF)
	  {
		numGaps++;
		basesInGaps += basesInThisGap;
	  }
	  
	  if (scaffoldContainsBACFrags)
	  {
		basesInThisGap = MIN (rchunkTemp->offsetAEnd.mean, rchunkTemp->offsetBEnd.mean) - 
		  MAX (lchunkTemp->offsetAEnd.mean, lchunkTemp->offsetBEnd.mean);
		// don't include "gaps" where contigs overlap
		if (basesInThisGap > BASE_CUTOFF)
		{
		  numGapsInScaffoldsWithBACFrags++;
		  basesInGapsInScaffoldsWithBACFrags += basesInThisGap;
		}
	  }	
	}	
  }
  
  fprintf( stderr, "Gap size cutoff: %d\n", BASE_CUTOFF);
  fprintf( stderr, "Number of gaps: %d\n", numGaps);
  fprintf( stderr, "Bases in gaps: %f\n", basesInGaps);
  fprintf( stderr, "Number of gaps in scaffolds with BAC frags: %d\n", numGapsInScaffoldsWithBACFrags);
  fprintf( stderr, "Bases in gaps in scaffolds with BAC frags: %f\n", basesInGapsInScaffoldsWithBACFrags);
}


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

        Module:  FbacREZ.c

   Description:  FBacREZ is a collection of functions that facilitate
                 the use of available FBAC information
                

    Programmer:  M. Flanigan

       Written:  7 Sep 00

  Last Revised:  

 **********************************************************************/

static char fileID[] = "$Id: FbacREZ.c,v 1.3 2005-03-22 19:07:34 jason_miller Exp $";

#define FBACDEBUG 2

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"

//
// AS_CNS
//
#include "MultiAlignment_CNS.h"
//
// AS_CGW
//
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ChiSquareTest_CGW.h"
#include "GraphCGW_T.h"

//
// AS_REZ
//
#include "DataTypesREZ.h"
#include "UtilsREZ.h"
#include "CommonREZ.h"
#include "SubgraphREZ.h"
#include "GapWalkerREZ.h"
#include "FbacREZ.h"
#include "RepeatRez.h"
#include "dpc_CNS.h"
#include "GWDriversREZ.h"

//
// vars
//

static VA_TYPE(char) *consensus1 = NULL;
static VA_TYPE(char) *consensus2 = NULL;
static VA_TYPE(char) *quality1 = NULL;
static VA_TYPE(char) *quality2 = NULL;

/*** mjf ***/ // cleanup
#define FBACMAXDIFF 5
#define BASES_FROM_END_CUTOFF 2000
#define BACWALK_THRESHOLD 0.5
// extern float Bayesian_Quality(CIEdgeT *, int32);

void ResetFrag( CIFragT *currFrag, int currFragOriginalContig, 
	       double currFragOriginalOffset5p, double currFragOriginalOffset3p);
int CheckWalkOverlaps( ChunkInsertInfoT *walkedContigs, int numChunks, int fixOverlapFailures);
int CheckRchunkContainment( ChunkInsertInfoT *walkedChunks, ChunkInsertInfoT *lastWalkedChunk);

// used to allow filter function to detect which locale we're filtering on
// int buildLocale;

/* returns TRUE if the node (Chunk instance or Contig) 
   has the hasFinishedBacFragments bit set */
int contains_fbac(NodeCGW_T* node)
{
  return (node->flags.bits.includesFinishedBacFragments);
}

int getLocalesInNode(NodeCGW_T* node, localeInfoT** localeInfo, int end, 
                     unsigned int basesFromEndCutoff)
{
  int i, ii;
  GraphCGW_T *graph = ScaffoldGraph->ContigGraph;
  MultiAlignT *ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, node->id, graph->type == CI_GRAPH); 
  // IntMultiPos *mp;
  int seen, 
	numLocales = 0,
        numUnitigs;
  localeInfoT * currentLocale, 
	* firstLocale = NULL,
	* lastLocale = NULL;
  int 
	basesFromEnd = 0;
  
  if( !contains_fbac(node) )
    return 0;

  numUnitigs = GetNumIntUnitigPoss( ma->u_list );

#if FBACDEBUG > 2
  fprintf( stderr, "\n\n\n********************************************************\n");
  fprintf(stderr,"Node with ID %d contains fbacs\n", node->id);
#endif

  numLocales = 0;

  for( i = 0; i < numUnitigs; i++)
  {
  	IntUnitigPos *upos = GetIntUnitigPos( ma->u_list, i);
	ChunkInstanceT *originalUnitig,
          *unitig = GetGraphNode( ScaffoldGraph->CIGraph, upos->ident);
	MultiAlignT *uma;
	int icntfrag, surrogateOffset, surrogateOrientation;
	int tempFragContigOffset5p = 0, tempFragContigOffset3p = 0;
	CIFragT *frag;

	// we don't need to examine unitigs more than basesFromEndCutoff away from apropos end
	if (end == A_END)
	{
	  if ((unitig->offsetAEnd.mean >= basesFromEndCutoff) && 
		  (unitig->offsetBEnd.mean >= basesFromEndCutoff))
	  {
		// fprintf( stderr, "(1)  contig %d, skipping unitig: %d (%f, %f)\n", node->id, unitig->id, 
		//	 unitig->offsetAEnd.mean, unitig->offsetBEnd.mean);
		continue;
	  }
	}
	else
	{
	  if ((node->bpLength.mean - unitig->offsetAEnd.mean >= basesFromEndCutoff) && 
		  (node->bpLength.mean - unitig->offsetBEnd.mean >= basesFromEndCutoff))
	  {
		// fprintf( stderr, "(2)  contig %d, skipping unitig: %d (%f, %f)\n", node->id, unitig->id, 
		//	 unitig->offsetAEnd.mean, unitig->offsetBEnd.mean);
		continue;
	  }
	}
	  
	  
	//if (unitig->info.CI.contigID != node->id)
	//{
	//fprintf( stderr, "suurogate: unitig->info.CI.contigID (%d) != node->id (%d)\n", unitig->info.CI.contigID, node->id);
	//}
	
	originalUnitig = unitig;
	if (unitig->flags.bits.isStoneSurrogate || unitig->flags.bits.isWalkSurrogate)
	{
	  if (((unitig->offsetAEnd.mean < unitig->offsetBEnd.mean) && end == A_END) ||
		  ((unitig->offsetAEnd.mean > unitig->offsetBEnd.mean) && end == B_END))
	  {
		surrogateOrientation = 1;
		surrogateOffset = (int) unitig->offsetAEnd.mean;
	  }
	  else
	  {
		surrogateOrientation = -1;
		surrogateOffset = (int) unitig->offsetBEnd.mean;
	  }	  
	  unitig = GetGraphNode( ScaffoldGraph->CIGraph, unitig->info.CI.baseID);
	}
	else
	{
	  surrogateOrientation = 0;
	  surrogateOffset = 0;
	}

	uma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, ScaffoldGraph->CIGraph->type == CI_GRAPH); 

	// now get info on the frags in the unitig
	for (icntfrag = 0; icntfrag < GetNumIntMultiPoss(uma->f_list); icntfrag++)
	{
	  IntMultiPos *imp = GetIntMultiPos(uma->f_list, icntfrag);

	  frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32) imp->source);

	  // keep track of whether the contig that this frag is in the assembly
	  // or if it was reached via a surrogate
	  if (frag->contigID != node->id && frag->locale != NULLINDEX)
	  {
		int originalUnitig_left_end, originalUnitig_right_end, originalUnitigContigOrientation;
		int frag_left_end, frag_right_end, fragmentContigOrientation;
		
		fprintf( stderr, "frag->contigID (%d) != node->id (%d)\n", frag->contigID, node->id);
		// continue;

		GetChunkPositionInContig(originalUnitig, &originalUnitig_left_end, &originalUnitig_right_end, 
								 &originalUnitigContigOrientation);
		fprintf( stderr, 
				 "surrogate: %d position in contig %d: surrogate_left_end: %d, surrogate_right_end: %d, orientation: %d\n", 
				 originalUnitig->id, node->id, originalUnitig_left_end, originalUnitig_right_end, 
				 originalUnitigContigOrientation);

		GetFragmentPositionInContigFromChunk(frag, &frag_left_end, &frag_right_end, 
											 &fragmentContigOrientation, 
											 originalUnitigContigOrientation,
											 originalUnitig_left_end, originalUnitig_right_end);
		fprintf(stderr, "frag iid: %d, frag_left_end: %d, frag_right_end: %d, fragmentContigOrientation: %d\n",
				frag->iid, frag_left_end, frag_right_end, fragmentContigOrientation);

		if (fragmentContigOrientation == 0)
		{
		  tempFragContigOffset5p = frag_left_end;
		  tempFragContigOffset3p = frag_right_end;
		}
		else
		{
		  tempFragContigOffset5p = frag_right_end;
		  tempFragContigOffset3p = frag_left_end;
		}
	  }
	  
	  if( frag->locale != NULLINDEX )
	  {
		// we know which end of the contig we are interested in
		if (end == A_END)
		{
		  if (surrogateOrientation == 0) // non-surrogate case
			basesFromEnd = min((int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean);
		  else
		  {
			basesFromEnd = min( tempFragContigOffset5p, tempFragContigOffset3p);
			// if (surrogateOrientation == 1) // surrogate is AB in contig
			  // basesFromEnd = surrogateOffset + min((int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean);
			// else // surrogate is BA in contig
			  // basesFromEnd = surrogateOffset + 
			  //			(originalUnitig->bpLength.mean - max((int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean));
			//basesFromEnd = max(tempFragContigOffset5p, tempFragContigOffset3p);
		  }
		}
		else if (end == B_END)
		{
		  if (surrogateOrientation == 0) // non-surrogate case
		  {
			basesFromEnd = min(node->bpLength.mean - (int) frag->contigOffset5p.mean,
							   node->bpLength.mean - (int) frag->contigOffset3p.mean);
		  }
		  else
		  {
			basesFromEnd = min(node->bpLength.mean - tempFragContigOffset5p,
							   node->bpLength.mean - tempFragContigOffset3p);
			//if (surrogateOrientation == 1) // surrogate is AB in contig
			//basesFromEnd = node->bpLength.mean - 
			//(surrogateOffset + max((int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean));
			//else // surrogate is BA in contig
			//basesFromEnd = node->bpLength.mean - 
			//(surrogateOffset + 
			// (originalUnitig->bpLength.mean - min((int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean)));
		  }
		}
		
#if FBACDEBUG > 10
		fprintf(stderr, "frag->iid: %d, locale: %d, location (5p, 3p): (%d, %d), bases from %s: %d\n", 
				frag->iid, frag->locale,
				(int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean,
				(end == A_END) ? "A_END" : "B_END", basesFromEnd); 
#endif
		
		// fprintf( stderr, "frag->iid: %d, basesFromEnd: %d\n", frag->iid, basesFromEnd);
		
		// don't bother with frags that are not within a certain distance from end
		if ( basesFromEnd > basesFromEndCutoff)
		{
		  // fprintf( stderr, "not going to look at frag: %d, locale: %d\n", frag->iid, frag->locale);
		  continue;
		}
		
		// determine if we've seen this locale before
		seen = FALSE;
		currentLocale = firstLocale;
		for (ii = 0; ii < numLocales; ii++)
		{
		  //fprintf(stderr, "currentLocale->localeNumber = %d, frag->locale = %d\n",
		  //currentLocale->localeNumber, frag->locale);
		  
		  if (currentLocale->localeNumber == frag->locale)
		  {
			seen = TRUE;
			break;
		  }
		  currentLocale = currentLocale->next;
		}
		
		if (!seen)
		{
		  localeInfoT *newLocale = malloc(sizeof(localeInfoT));
		  if (newLocale == NULL)
		  {
			fprintf(stderr, "could not malloc new localeInfo\n");
			exit(1);
		  }

		  numLocales++;
		  
#if FBACDEBUG > 2
		  fprintf(stderr, "adding locale %d to localeInfo( frag (5p, 3p): %d (%f, %f))\n", 
				  frag->locale, frag->iid, frag->contigOffset5p.mean, frag->contigOffset3p.mean);
#endif
		  // add it to the linked list
		  if (firstLocale == NULL)
			firstLocale = newLocale;
		  if (lastLocale != NULL)
			lastLocale->next = newLocale;
		  lastLocale = newLocale;
		  
		  // use the first frag we run across from this locale as the initial reference point
		  newLocale->localeNumber = frag->locale;
          newLocale->contigID = node->id;
		  newLocale->extremalFragID = (int32) imp->source;
		  newLocale->extremalFragIid = frag->iid;
		  newLocale->extremalFragUnitigID = originalUnitig->id;  // unitig may be a surrogate
		  newLocale->basesFromEnd = basesFromEnd; // frag->contigOffset5p.mean;
		  newLocale->next = NULL;
		  newLocale->inSurrogate = abs( surrogateOrientation );
		  newLocale->fragContigOffset5p = tempFragContigOffset5p;
		  newLocale->fragContigOffset3p = tempFragContigOffset3p;		  
		  newLocale->BacOffset = frag->localePos.bgn;

		  //fprintf(stderr, "frag->iid: %d, locale: %d, utg: %d orig_utg: %d, ctg (5p, 3p): (%d, %d), bases from %s: %d\n", 
		  //	  frag->iid, frag->locale, unitig->id, originalUnitig->id, 
		  //	  (int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean,
		  //	  (end == A_END) ? "A_END" : "B_END", basesFromEnd);

		  //newLocale->contigOffset = (int) frag->contigOffset5p.mean;
		  //newLocale->numFbacFrags = 0;
		  
		  currentLocale = newLocale;
		  currentLocale->numLocaleFragments = 0;
		  currentLocale->iidTrend = 0;
		}

		currentLocale->numLocaleFragments++;

		// if the new frag for this locale is closer to the end of the contig than the current one
		if (basesFromEnd < currentLocale->basesFromEnd)  
		{
		  if (currentLocale->extremalFragIid < frag->iid) // ... and its iid is greater
			currentLocale->iidTrend = 1;                  // the trend in iids is positive
		  else
			currentLocale->iidTrend = -1;                 // or negative
		}
		else   // else the new frag for this locale is further from the end of the contig than the current one
		{
		  if (currentLocale->extremalFragIid < frag->iid) // ... and its iid is greater
			currentLocale->iidTrend = -1;                  // the trend in iids is negative
		  else
			currentLocale->iidTrend = 1;                 // or positive
		}
		

		if (basesFromEnd < currentLocale->basesFromEnd)
		{
          currentLocale->contigID = node->id;
		  currentLocale->extremalFragID = (int32) imp->source;
		  currentLocale->extremalFragIid = frag->iid;
		  currentLocale->extremalFragUnitigID = originalUnitig->id;  // unitig may be a surrogate
		  currentLocale->basesFromEnd = basesFromEnd; // frag->contigOffset5p.mean;
		  currentLocale->inSurrogate = abs( surrogateOrientation );
		  currentLocale->fragContigOffset5p = tempFragContigOffset5p;
		  currentLocale->fragContigOffset3p = tempFragContigOffset3p;		  
		  currentLocale->BacOffset = frag->localePos.bgn;
		  
		  // fprintf(stderr, "frag->iid: %d, imp->source: %d, locale: %d, utg: %d orig_utg: %d, ctg (5p, 3p): (%d, %d), bases from %s: %d\n", 
		  //  frag->iid, (int32) imp->source, frag->locale, unitig->id, originalUnitig->id, 
		  //	  (int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean,
		  //	  (end == A_END) ? "A_END" : "B_END", basesFromEnd);
		}
	  }
	}
  }
#if FBACDEBUG > 2
  fprintf( stderr, "Node %d contains %d locales within %d bases of end %s\n", 
		   node->id, numLocales, basesFromEndCutoff, (end == A_END) ? "A_END" : "B_END");
#endif
  *localeInfo = firstLocale;
  return numLocales;
}

#if 1
/* returns   
   FBAC_CONSISTENT,
   if the node (Chunk instance or Contig) 
   is correctly covered by fbac fragments, 
   FBAC_INCONSISTENT if there is a discrapency
   NO_FBACS if the chunk has no NO_FBACS at all */

BacStatusREZ isFbacConsistent(NodeCGW_T* node)
{
  int i;
  MultiAlignT *ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, node->id, 0); 
  // MultiAlignT *ma = GetMultiAlignInStore(ScaffoldGraph->RezGraph->maStore, node->id);
  int32 numFrags;
  BacStatusREZ ret = BAC_CONSISTENT;

  float nodelength = node->bpLength.mean;
  IntMultiPos *mp;
  CIFragT *frag;
  localeInfoT *locales;
  int numLocales;
  fragmentInfoT *fragments;
  
  if( !contains_fbac(node) )
    return NO_BACS;
  
  // it doen't matter what end we choose, since going to look over whole node
  numLocales = getLocalesInNode( node, &locales, A_END, UINT_MAX);

  numFrags = GetNumIntMultiPoss(ma->f_list);

  if (numFrags == 1)
	return BAC_CONSISTENT;

  fragments = (fragmentInfoT *) malloc (sizeof( fragmentInfoT ) * numFrags);
  
  for (i = 0; i < numLocales; i++)
  {
	int currentLocale;
	int fragCnt;
	int contigDiff, localeDiff, iidDiff;
	int prevFragIid, prevFrag5pContigPos, prevFrag5pLocalePos;
	int fragsThisLocale;
	float localeDist, contigDist;  
	
	// should we check if all fbacs from a locale have the same orientation?
	
	currentLocale = locales->localeNumber;
	
#if FBACDEBUG > 0
	fprintf( stderr, "\n--------------------------------------------------------\n");
	fprintf( stderr, "Node with ID %d and length %f and numFrags %d and %d locales\n", 
			 node->id, nodelength, numFrags, numLocales);
	fprintf( stderr, "Looking for locale %d\n", currentLocale);
#endif

	/* in this loop we check whether the scaffold positions are consistent with
	   the locale positions of the fragments */

	fragsThisLocale = 0;
	for( fragCnt = 0; fragCnt < numFrags; fragCnt++)
	{
	  mp   = GetIntMultiPos(ma->f_list, fragCnt);
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32) mp->source);
	  
	  // fprintf( stderr, "frag->iid: %d, frag->locale: %d\n", frag->iid, frag->locale);

	  if( frag->locale != currentLocale )
		continue;

	  fragments[ fragsThisLocale ].fragIid = frag->iid;
	  fragments[ fragsThisLocale ].frag5pContigPos = (int) frag->contigOffset5p.mean;
	  fragments[ fragsThisLocale ].frag5pLocalePos = frag->localePos.bgn;
	  fragsThisLocale++;
	}
	
	qsort( fragments, fragsThisLocale, sizeof(fragmentInfoT), &compFragments);

	localeDist = fabs( fragments[fragsThisLocale - 1].frag5pLocalePos - 
					   fragments[0].frag5pLocalePos) + 550.0;
	contigDist = fabs( fragments[fragsThisLocale - 1].frag5pContigPos - 
					   fragments[0].frag5pContigPos) + 550.0;
	fprintf( stderr, "max frag iid: %d\n", fragments[fragsThisLocale - 1].fragIid);
	fprintf( stderr, "min frag iid: %d\n", fragments[0].fragIid);
	fprintf( stderr, "range: %d, numLocaleFrags: %d, difference: %d, locale dist: %f, contig dist: %f, diff: %f\n", 
			 fragments[fragsThisLocale - 1].fragIid - fragments[0].fragIid + 1, fragsThisLocale, 
			 fragments[fragsThisLocale - 1].fragIid - fragments[0].fragIid + 1 - fragsThisLocale,
			 localeDist, contigDist, fabs( localeDist - contigDist));
	fprintf( stderr, "100 * diff/locale dist: %.2f\n", 100 * fabs( localeDist - contigDist) / localeDist);
	
	for( fragCnt = 1; fragCnt < fragsThisLocale; fragCnt++)
	{
	  // find out where this frag is w/respect to the previous frag
	  iidDiff = abs( fragments[ fragCnt - 1 ].fragIid - fragments[ fragCnt ].fragIid );
	  contigDiff = abs( fragments[ fragCnt - 1 ].frag5pContigPos - fragments[ fragCnt ].frag5pContigPos );
	  localeDiff = abs( fragments[ fragCnt - 1 ].frag5pLocalePos - fragments[ fragCnt ].frag5pLocalePos );
	  
	  //fprintf( stderr, "for frag iid %d, prevFrag5pContigPos %.2f, prevFrag5pLocalePos: %d\n",
	  //	   frag->iid, prevFrag5pContigPos, prevFrag5pLocalePos);
	  
	  if ( (abs( contigDiff - localeDiff) > 30) ||
		   iidDiff > 1)
	  {
		fprintf( stderr, "contig: %d frag iid: %8d 5p contig: %8d 5p locale: %5d\n", node->id,
				 fragments[ fragCnt - 1 ].fragIid, fragments[ fragCnt - 1 ].frag5pContigPos,
				 fragments[ fragCnt - 1 ].frag5pLocalePos);
		fprintf( stderr, "contig: %d frag iid: %8d 5p contig: %8d 5p locale: %5d diffs: %4d %8d %8d %8d\n",
				 node->id, fragments[ fragCnt ].fragIid, fragments[ fragCnt ].frag5pContigPos,
				 fragments[ fragCnt ].frag5pLocalePos, iidDiff, contigDiff, localeDiff, abs( contigDiff - localeDiff));
		//fprintf( stderr, "diffs \t\t%8d \t\t\t%8d \t\t\t%8d \t%5d\n\n",
		//	 iidDiff, contigDiff, localeDiff, abs( contigDiff - localeDiff));
	  }
	  
	  // now set up for next frag
	  prevFragIid = frag->iid;
	  prevFrag5pContigPos = frag->contigOffset5p.mean;
	  prevFrag5pLocalePos = frag->localePos.bgn;
	}	
	locales = locales->next;
  }
  return ret;
}

#endif



/* returns TRUE if the edge connects two chunks that fullfill
   is_suitable_fbac_chunk (e.g. is_covered_by_fbac or contains_fbac
   and if the indicated overlap of this edge is consistent
   with the locpos fields of the chunks */
int edge_is_fbac_consistent(EdgeCGW_T* e)
{
  return TRUE;
}

int FindCommonLocales( ContigT * lcontig, int lcontigGapEnd, 
                       ContigT * rcontig, int rcontigGapEnd,
                       GapInfoT **gapInfoArray)
{
  int numLeftLocales, numRightLocales;
  int numCommonLocales = 0;
  int i, j;
  int  tmp_size;
  GapInfoT *gapInfoArrayTemp;

  localeInfoT
	*localeInfoLcontig = NULL,
	*localeInfoRcontig = NULL,
	*localeInfoLcontigTemp;

  if( !contains_fbac( lcontig ) || !contains_fbac( rcontig ))
  {
	numLeftLocales = 0;
	numRightLocales = 0;
  }  
  else
  {
	numLeftLocales = getLocalesInNode(lcontig, &localeInfoLcontig, lcontigGapEnd, BASES_FROM_END_CUTOFF);	  
	numRightLocales = getLocalesInNode(rcontig, &localeInfoRcontig, rcontigGapEnd, BASES_FROM_END_CUTOFF);  
  }
  
  tmp_size = min( numLeftLocales, numRightLocales);
  if  (tmp_size < 1)
      tmp_size = 1;
  gapInfoArrayTemp = (GapInfoT *) malloc (tmp_size * sizeof( GapInfoT));

  // find locales that are in both lists
  // locales are numbered starting at 1
  localeInfoLcontigTemp = localeInfoLcontig;
  for (i = 0; i < numLeftLocales; i++)
  {
	localeInfoT* localeInfoRcontigTemp = localeInfoRcontig;
	for (j = 0; j < numRightLocales; j++)
	{
	  // fprintf( stderr, "comparing left locale: %d to right locale: %d\n", 
	  //   localeInfoLcontigTemp->localeNumber, localeInfoRcontigTemp->localeNumber);
	  if (localeInfoLcontigTemp->localeNumber == localeInfoRcontigTemp->localeNumber) // &&   // same locale
		  // localeInfoLcontigTemp->numLocaleFragments > 5 &&  // with enough support
		  // localeInfoRcontigTemp->numLocaleFragments > 5 )
	  {
		gapInfoArrayTemp[numCommonLocales].localeNumber = localeInfoLcontigTemp->localeNumber;

		gapInfoArrayTemp[numCommonLocales].fragIDLeftContig = localeInfoLcontigTemp->extremalFragID;
		gapInfoArrayTemp[numCommonLocales].fragLeftContigID = localeInfoLcontigTemp->contigID;
		gapInfoArrayTemp[numCommonLocales].fragLeftContigOffset5p = localeInfoLcontigTemp->fragContigOffset5p;
		gapInfoArrayTemp[numCommonLocales].fragLeftContigOffset3p = localeInfoLcontigTemp->fragContigOffset3p;
		gapInfoArrayTemp[numCommonLocales].unitigIDLeftContig = localeInfoLcontigTemp->extremalFragUnitigID;
		gapInfoArrayTemp[numCommonLocales].fragInSurrogateLeftContig = localeInfoLcontigTemp->inSurrogate;
		gapInfoArrayTemp[numCommonLocales].basesFromEndLeftContig = localeInfoLcontigTemp->basesFromEnd;

		gapInfoArrayTemp[numCommonLocales].fragIDRightContig = localeInfoRcontigTemp->extremalFragID;
		gapInfoArrayTemp[numCommonLocales].fragRightContigID = localeInfoRcontigTemp->contigID;
		gapInfoArrayTemp[numCommonLocales].fragRightContigOffset5p = localeInfoRcontigTemp->fragContigOffset5p;
		gapInfoArrayTemp[numCommonLocales].fragRightContigOffset3p = localeInfoRcontigTemp->fragContigOffset3p;
		gapInfoArrayTemp[numCommonLocales].unitigIDRightContig = localeInfoRcontigTemp->extremalFragUnitigID;
		gapInfoArrayTemp[numCommonLocales].fragInSurrogateRightContig = localeInfoRcontigTemp->inSurrogate;
		gapInfoArrayTemp[numCommonLocales].basesFromEndRightContig = localeInfoRcontigTemp->basesFromEnd;		

		numCommonLocales++;
	  }
	  localeInfoRcontigTemp = localeInfoRcontigTemp->next;
	}
	localeInfoLcontigTemp = localeInfoLcontigTemp->next;
  }

  // free the linked lists (localeInfoLcontig and localeInfoRcontig)

  if (localeInfoLcontig)
  {
	localeInfoT *localeInfoTemp;
	for (i = 0; i < numLeftLocales - 1; i++)
	{
	  localeInfoTemp = localeInfoLcontig;
	  localeInfoLcontig = localeInfoLcontig->next;
	  free (localeInfoTemp);
	}
	free (localeInfoLcontig);
  }
  
  if (localeInfoRcontig)
  {
	localeInfoT *localeInfoTemp;
	for (i = 0; i < numRightLocales - 1; i++)
	{
	  localeInfoTemp = localeInfoRcontig;
	  localeInfoRcontig = localeInfoRcontig->next;
	  free (localeInfoTemp);
	}
	free (localeInfoRcontig);
  }  
  
#if DEBUG_GAP_WALKER > -1
  {
	fprintf(stderr, "numCommonLocales: %d\n", numCommonLocales);
  }
#endif

  *gapInfoArray = gapInfoArrayTemp;
  
  return numCommonLocales;
}

int compDistances( const void *s1, const void *s2)
{
  const GapInfoT * t1 = s1;
  const GapInfoT * t2 = s2;
  assert( t1 == s1 );
  assert( t2 == s2 );
  
  if (t1->distanceDifference < t2->distanceDifference)
	return -1;
  else if (t1->distanceDifference > t2->distanceDifference)
	return 1;
  else 
	return 0;
}

int compFragments( const void *s1, const void *s2)
{
  const fragmentInfoT * t1 = s1;
  const fragmentInfoT * t2 = s2;
  assert( t1 == s1 );
  assert( t2 == s2 );
  
  if (t1->fragIid < t2->fragIid)
	return -1;
  else if (t1->fragIid > t2->fragIid)
	return 1;
  else 
	return 0;
}

int compChunks( const void *s1, const void *s2)
{
  const ChunkInsertInfoT * t1 = s1;
  const ChunkInsertInfoT * t2 = s2;
  assert( t1 == s1 );
  assert( t2 == s2 );
  
  if (min (t1->aEndOffset.mean, t1->bEndOffset.mean) < min( t2->aEndOffset.mean, t2->bEndOffset.mean))
	return -1;
  else if (min (t1->aEndOffset.mean, t1->bEndOffset.mean) > min( t2->aEndOffset.mean, t2->bEndOffset.mean))
	return 1;
  else 
	return 0;
}

void SortCommonLocales( int numCommonLocales, GapInfoT **gapInfoArray, int intraScaffoldGap, LengthT gapEstimate)
{
  int leftFragLeftEnd, leftFragRightEnd;
  int rightFragLeftEnd, rightFragRightEnd;
  int leftFragScaffoldOrientation, rightFragScaffoldOrientation;
  int ilocale;
  GapInfoT *gapInfoArrayLocal = *gapInfoArray;
  CIFragT *leftFrag, *rightFrag;
  // float fragSeparationScaffold, fragSeparationLocale;

  // compute for each locale how well it agrees with our gap estimate
  for (ilocale = 0; ilocale < numCommonLocales; ilocale++)
  {
	leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, gapInfoArrayLocal[ilocale].fragIDLeftContig);
	GetFragmentPositionInScaffold(leftFrag, &leftFragLeftEnd, &leftFragRightEnd, &leftFragScaffoldOrientation);
	
	rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, gapInfoArrayLocal[ilocale].fragIDRightContig);
	GetFragmentPositionInScaffold(rightFrag, &rightFragLeftEnd, &rightFragRightEnd, &rightFragScaffoldOrientation);
	
	// find out separation in frag positions implied by positions in locale
	if (leftFrag->iid < rightFrag->iid)
	  gapInfoArrayLocal[ilocale].fragSeparationLocale = rightFrag->localePos.bgn - leftFrag->localePos.end;
	else
	  gapInfoArrayLocal[ilocale].fragSeparationLocale = leftFrag->localePos.bgn - rightFrag->localePos.end;

	if (intraScaffoldGap == 1)
	{
	  // find out separation in frag positions implied by positions in scaffold
	  gapInfoArrayLocal[ilocale].fragSeparationScaffold = rightFragLeftEnd - leftFragRightEnd;
	  gapInfoArrayLocal[ilocale].distanceDifference = 
		fabs( (double) (gapInfoArrayLocal[ilocale].fragSeparationScaffold - 
						gapInfoArrayLocal[ilocale].fragSeparationLocale));
	}
	else
	{
	  gapInfoArrayLocal[ilocale].distanceDifference = 
		fabs( (double) (gapEstimate.mean - gapInfoArrayLocal[ilocale].fragSeparationLocale));	
	}
	if (1)
	{
	  fprintf( stderr, "intraScaffoldGap: %d\n", intraScaffoldGap);
	  fprintf( stderr, "locale %d, fragSeparationScaffold: %f, fragSeparationLocale: %f\n",
		   gapInfoArrayLocal[ilocale].localeNumber, gapInfoArrayLocal[ilocale].fragSeparationScaffold, 
			   gapInfoArrayLocal[ilocale].fragSeparationLocale);
	  fprintf( stderr, " leftFrag->iid: %d, scaffold pos: %d, %d\n",
		   leftFrag->iid, leftFragLeftEnd, leftFragRightEnd);
	  fprintf( stderr, "rightFrag->iid: %d, scaffold pos: %d, %d\n",
		   rightFrag->iid, rightFragLeftEnd, rightFragRightEnd);
	  fprintf( stderr, "\n");
	}
  }
  
  // sort and arrange nodes of list
  if (numCommonLocales > 1)
	qsort(gapInfoArrayLocal, numCommonLocales, sizeof( GapInfoT ), &compDistances);
}

void SetFragPositions( int numCommonLocales, GapInfoT **gapInfoArray)
{
  GapInfoT *gapInfoArrayLocal = *gapInfoArray;
  int i;
  
  for ( i = 0; i < numCommonLocales; i++)
  {
	if (gapInfoArrayLocal[i].fragInSurrogateLeftContig == TRUE)
	{
	  // get the frag
	  CIFragT *leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, gapInfoArrayLocal[i].fragIDLeftContig);

	  // save its offsets in its current contig
	  gapInfoArrayLocal[i].fragOriginalLeftContigID = leftFrag->contigID;
	  gapInfoArrayLocal[i].fragOriginalLeftContigOffset5p = leftFrag->contigOffset5p.mean;
	  gapInfoArrayLocal[i].fragOriginalLeftContigOffset3p = leftFrag->contigOffset3p.mean;
	  
	  // set them in the frag array to what we want them to be
	  leftFrag->contigID = gapInfoArrayLocal[i].fragLeftContigID;
	  leftFrag->contigOffset5p.mean = gapInfoArrayLocal[i].fragLeftContigOffset5p;
	  leftFrag->contigOffset3p.mean = gapInfoArrayLocal[i].fragLeftContigOffset3p;
	  
	  SetVA_CIFragT(ScaffoldGraph->CIFrags, gapInfoArrayLocal[i].fragIDLeftContig, leftFrag);
	}

	if (gapInfoArrayLocal[i].fragInSurrogateRightContig == TRUE)
	{
	  // get the frag
	  CIFragT *rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, gapInfoArrayLocal[i].fragIDRightContig);

	  // save its offsets in its current contig
	  gapInfoArrayLocal[i].fragOriginalRightContigID = rightFrag->contigID;
	  gapInfoArrayLocal[i].fragOriginalRightContigOffset5p = rightFrag->contigOffset5p.mean;
	  gapInfoArrayLocal[i].fragOriginalRightContigOffset3p = rightFrag->contigOffset3p.mean;
	  
	  // set them in the frag array to what we want them to be
	  rightFrag->contigID = gapInfoArrayLocal[i].fragRightContigID;
	  rightFrag->contigOffset5p.mean = gapInfoArrayLocal[i].fragRightContigOffset5p;
	  rightFrag->contigOffset3p.mean = gapInfoArrayLocal[i].fragRightContigOffset3p;
	  
	  SetVA_CIFragT(ScaffoldGraph->CIFrags, gapInfoArrayLocal[i].fragIDRightContig, rightFrag);
	}
  }
}

void	ResetFragPositions( int numCommonLocales, GapInfoT **gapInfoArray)
{
  GapInfoT *gapInfoArrayLocal = *gapInfoArray;
  int i;
  
  for ( i = 0; i < numCommonLocales; i++)
  {
	if (gapInfoArrayLocal[i].fragInSurrogateLeftContig == TRUE)
	{
	  // get the frag
	  CIFragT *leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, gapInfoArrayLocal[i].fragIDLeftContig);

	  // set them in the frag array to what they used to be
	  leftFrag->contigID = gapInfoArrayLocal[i].fragOriginalLeftContigID;
	  leftFrag->contigOffset5p.mean = gapInfoArrayLocal[i].fragOriginalLeftContigOffset5p;
	  leftFrag->contigOffset3p.mean = gapInfoArrayLocal[i].fragOriginalLeftContigOffset3p;
	  
	  SetVA_CIFragT(ScaffoldGraph->CIFrags, gapInfoArrayLocal[i].fragIDLeftContig, leftFrag);
	}

	if (gapInfoArrayLocal[i].fragInSurrogateRightContig == TRUE)
	{
	  // get the frag
	  CIFragT *rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, gapInfoArrayLocal[i].fragIDRightContig);

	  // set them in the frag array to what they used to be
	  rightFrag->contigID = gapInfoArrayLocal[i].fragOriginalRightContigID;
	  rightFrag->contigOffset5p.mean = gapInfoArrayLocal[i].fragOriginalRightContigOffset5p;
	  rightFrag->contigOffset3p.mean = gapInfoArrayLocal[i].fragOriginalRightContigOffset3p;
	  
	  SetVA_CIFragT(ScaffoldGraph->CIFrags, gapInfoArrayLocal[i].fragIDRightContig, rightFrag);
	}
  }
}

int Set_Bac_Walking_Quality(int contigA, int contigB, int orient, float quality)
{
  ChunkOverlapSpecT spec;
  ChunkOverlapperT *chunkOverlapper = ScaffoldGraph->RezGraph->overlapper;
  ChunkOverlapCheckT *lookup;
  int isCanonical;

  isCanonical = InitCanonicalOverlapSpec(contigA,contigB, orient, &spec);
  lookup = LookupCanonicalOverlap(chunkOverlapper, &spec);
  if(!lookup)  // We didn't find anything
    return FALSE;

  // There was something in the table but the overlap length is 0
  if( lookup->overlap == 0 )
    return FALSE;

  // There is an overlap, set its quality
  lookup->quality = quality;
  return TRUE;
}

ChunkInsertInfoT* BuildLocaleOverlaps(GapInfoT *gapInfo, 
				      ChunkInstanceT *lchunk, ChunkInstanceT *rchunk, int *completeOverlapPath,
				      double *rchunkNewAEnd, double *rchunkNewBEnd,
				      ChunkOrientationType *tempOlapOrientation,
				      Overlap* rchunkOverlap, int fixOverlapFailures)
{
  int ifrag;
  int startFragIid, endFragIid;
  // ChunkOrientationType tempOlapOrientation;
  ChunkInsertInfoT *currentWalkedChunk, *walkedChunks = NULL;
  int numFrags;
  int succNodeBeg, succNodeEnd;
  ChunkInstanceT *succNode;
  int increment;
  int walkDone;
  int lchunkFragAssigned, lchunkFragOriginalUnitig;
  int rchunkFragAssigned, rchunkFragOriginalUnitig;
  int insertCounter = 0;
  static VA_TYPE(int32) *assignedLChunkFragments = NULL;
  static VA_TYPE(int32) *assignedRChunkFragments = NULL;


  // localeGraph should be made out of only Contigs that contain the locale we're working with
  // the current hack is to store that locale number in the global buildLocale
  // so it's accessable to the node filter (e.g., Contains_Locale_Frags)
  buildLocale = gapInfo->localeNumber;
  
  // find out the frag number in each that we want to include in overlap calcs
  // GetFragInfo(lchunk, locale, &minLchunkFragIid, &maxLchunkFragIid);
  // GetFragInfo(rchunk, locale, &minRchunkFragIid, &maxRchunkFragIid);	
  
  // now compute the overlaps we need
  // fprintf( stderr, "computing overlap for localeGraph->size: %d\n", localeGraph->size);

  // at the beginning the lchunk is the successor 
  if (lchunk->offsetAEnd.mean < lchunk->offsetBEnd.mean) 
  {
	succNodeBeg = (int) lchunk->offsetAEnd.mean;
	succNodeEnd = (int) lchunk->offsetBEnd.mean;
  }
  else
  {
	succNodeBeg = (int) lchunk->offsetBEnd.mean;
	succNodeEnd = (int) lchunk->offsetAEnd.mean;
  }

  startFragIid = GetCIFragT(ScaffoldGraph->CIFrags, (int32) gapInfo->fragIDLeftContig)->iid;
  endFragIid = GetCIFragT(ScaffoldGraph->CIFrags, (int32) gapInfo->fragIDRightContig)->iid;
  
#if DEBUG_GAP_WALKER > -1		  
  fprintf( stderr, "walking from startFragIid: %d (contig: %d) to endFragIid: %d (contig: %d)\n",
		   startFragIid, lchunk->id, endFragIid, rchunk->id);
  
  fprintf( stderr, "succNodeBeg = %d, succNodeEnd = %d (current pos of lchunk)\n", succNodeBeg, succNodeEnd);
#endif

  // if the left chunk's extreme fragment is in a surrogate, temporarily assign this frag to the surrogate
  if (0)
  {
	CIFragT *leftFrag;
	InfoByIID *info;

    // get the CIFragT for the frag in the lchunk, check if it is in a surrogate
	info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, startFragIid);
	assert(info->set);
	leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

	if (leftFrag->contigID != lchunk->id)
	{
	  fprintf( stderr, "leftFrag->contigID (%d) != lchunk->id (%d)\n", leftFrag->contigID, lchunk->id);
  
	  lchunkFragAssigned = TRUE;
	  lchunkFragOriginalUnitig = leftFrag->cid;
	  
	  if ( assignedLChunkFragments )
		ResetVA_int32( assignedLChunkFragments );
	  else
		assignedLChunkFragments = CreateVA_int32(1);
	  
	  SetVA_int32( assignedLChunkFragments, 0, &startFragIid);
	  
	  AssignFragsToResolvedCI(ScaffoldGraph->CIGraph, leftFrag->CIid, gapInfo->unitigIDLeftContig, assignedLChunkFragments);
	  
	  // get the CIFragT for the frag in the lchunk, now it should not be in a surrogate
	  info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, startFragIid);
	  assert(info->set);
	  leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
	  
	  fprintf( stderr, "after assignment, leftFrag->contigID (%d), lchunk->id (%d)\n", leftFrag->contigID, lchunk->id);
	}
  }
  
  // if the right chunk's extreme fragment is in a surrogate, temporarily assign this frag to the surrogate
  if (0)
  {
	CIFragT *rightFrag;
	InfoByIID *info;

    // get the CIFragT for the frag in the lchunk, check if it is in a surrogate
	info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, endFragIid);
	assert(info->set);
	rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

	if (rightFrag->contigID != rchunk->id)
	{
	  fprintf( stderr, "rightFrag->contigID (%d) != rchunk->id (%d)\n", rightFrag->contigID, rchunk->id);
	  
	  rchunkFragAssigned = TRUE;
	  rchunkFragOriginalUnitig = rightFrag->cid;
	  if ( assignedRChunkFragments )
		ResetVA_int32( assignedRChunkFragments );
	  else
		assignedRChunkFragments = CreateVA_int32(1);
	  
	  SetVA_int32( assignedRChunkFragments, 1, &endFragIid);
	  
	  AssignFragsToResolvedCI(ScaffoldGraph->CIGraph, rightFrag->cid, gapInfo->unitigIDRightContig, assignedRChunkFragments);
	  
	  // get the CIFragT for the frag in the rchunk, now it should not be in a surrogate
	  info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, startFragIid);
	  assert(info->set);
	  rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
	  
	  fprintf( stderr, "after assignment, rightFrag->contigID (%d), rchunk->id (%d)\n", rightFrag->contigID, rchunk->id);
	}
  }
  

  // check to see if they are adjacent (in iid) frags that do not overlap in the locale
  if (abs (startFragIid - endFragIid) == 1)
  {
	CIFragT *startFrag, *endFrag;
	InfoByIID *info;
	int overlapInLocale;
	
	info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, startFragIid);
	assert(info->set);
	startFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

	info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, endFragIid);
	assert(info->set);
	endFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

	overlapInLocale = TRUE;
	
	if ( startFrag->locale == endFrag->locale)
        {
	  if ( startFragIid < endFragIid)
	  {
            if (startFrag->localePos.end < endFrag->localePos.bgn) 
              overlapInLocale = FALSE;
	  }
	  else
            if (endFrag->localePos.end < startFrag->localePos.bgn) 
              overlapInLocale = FALSE;
        }
          
        if (!overlapInLocale)
	{
	  fprintf( stderr, "iid adjacent frags %d and %d do not overlap in their locale\n",
			   startFragIid, endFragIid);
#if 0
	  fprintf( stderr, "stopping walk because of condition 9.0 (consecutive iids do not overlap in locale)\n");
	  fprintf( stderr, "----------------\n");
	  *completeOverlapPath = -9;
	  return walkedChunks;
#endif
	}
  }

  if ( startFragIid < endFragIid)
  {
	numFrags = endFragIid - startFragIid;	
	increment = 1;
  }
  else
  {
	numFrags = startFragIid - endFragIid;	
	increment = -1;
  }
  
  ifrag = startFragIid;

  // first node in walked chunks is lchunk
  currentWalkedChunk = walkedChunks = (ChunkInsertInfoT *) malloc ( sizeof(ChunkInsertInfoT));
  currentWalkedChunk->next = NULL;
  currentWalkedChunk->contigID = lchunk->id;
  currentWalkedChunk->scaffoldID = lchunk->scaffoldID;
  currentWalkedChunk->aEndOffset.mean = lchunk->offsetAEnd.mean;
  currentWalkedChunk->aEndOffset.variance = lchunk->offsetAEnd.variance;
  currentWalkedChunk->bEndOffset.mean = lchunk->offsetBEnd.mean;
  currentWalkedChunk->bEndOffset.variance = lchunk->offsetBEnd.variance;  
  currentWalkedChunk->insertOrder = insertCounter++;
  currentWalkedChunk->hasMoved = FALSE;
  
  walkDone = FALSE;
  while ( !walkDone )
  {
	CIFragT *currFrag, *fragSucc1;
	int succ1, currUnitigID, currContigID;
	InfoByIID *info;
	int setFrag, currFragOriginalContig = NULLINDEX;
	double currFragOriginalOffset5p = 0., currFragOriginalOffset3p = 0.;
	// int fixPairwiseOverlapFailure = TRUE;
	Overlap dummyOlap;
	
	// get the CIFragT for ifrag, which is the first frag of interest in the current contig
	info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, ifrag);
	assert(info->set);
	currFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

	// save the id of the current contig in case findLastLocaleFragInContig returns a frag in a surrogate
	currContigID = currFrag->contigID;

	// find the last frag of interest in the current contig
	// we need the unitig id so we know where we are in the contig (eg, the frag is "in" a surrogate)
	ifrag = findLastLocaleFragInContig( currFrag->contigID, ifrag, endFragIid, increment, &currUnitigID);

	// get the CIFragT for (the possibly new) ifrag
	info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, ifrag);

	// we used to assert, but now we just end the walk gracefully
	// assert(info->set);
	if (info->set == 0)
	{
	  fprintf( stderr, "ifrag iid %d is not set!\n", ifrag);
	  // walk failed
	  fprintf( stderr, "stopping walk because of condition 7.0 (frag not set)\n");
	  fprintf( stderr, "----------------\n");
	  *completeOverlapPath = -7;
	  return walkedChunks;
	}

	currFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

	setFrag = 0;
	if (currFrag->contigID != currContigID)  // means frag is in a surrogate
	{
	  int surrogateLeftEnd, surrogateRightEnd, surrogateContigOrientation;
	  int fragLeftEnd, fragRightEnd, fragmentContigOrientation;
	  ChunkInstanceT *surrogate = GetGraphNode( ScaffoldGraph->CIGraph, currUnitigID );

	  setFrag = 1;

	  // save the frag's current contig and frag's position in such
	  currFragOriginalContig = currFrag->contigID;
	  currFragOriginalOffset5p = currFrag->contigOffset5p.mean;
	  currFragOriginalOffset3p = currFrag->contigOffset3p.mean;

	  GetChunkPositionInContig(surrogate, &surrogateLeftEnd, &surrogateRightEnd, 
				   &surrogateContigOrientation);
	  fprintf( stderr, 
		   "surrogate: %d position in contig %d: surrogateLeftEnd: %d, surrogateRightEnd: %d, orientation: %d\n", 
		   surrogate->id, currContigID, surrogateLeftEnd, surrogateRightEnd, 
		   surrogateContigOrientation);
	  
	  GetFragmentPositionInContigFromChunk(currFrag, &fragLeftEnd, &fragRightEnd, 
					       &fragmentContigOrientation,
					       surrogateContigOrientation, 
					       surrogateLeftEnd, surrogateRightEnd);
	  fprintf(stderr, "frag iid: %d, frag_left_end: %d, frag_right_end: %d, fragmentContigOrientation: %d\n",
		  currFrag->iid, fragLeftEnd, fragRightEnd, fragmentContigOrientation);

	  // set the frag to where it would be in currContig if not in surrogate
	  currFrag->contigID = currContigID;
	  if (fragmentContigOrientation == 0)
	  {
	    currFrag->contigOffset5p.mean = fragLeftEnd;
	    currFrag->contigOffset3p.mean = fragRightEnd;
	  }
	  else
	  {
	    currFrag->contigOffset5p.mean = fragRightEnd;
	    currFrag->contigOffset3p.mean = fragLeftEnd;
	  }	    

	  info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, currFrag->iid);
	  SetVA_CIFragT( ScaffoldGraph->CIFrags, info->fragIndex, currFrag);
	}

	// check to see that frag is within bounds
	if (FragOutOfBounds(currFrag, currUnitigID))
	{
	  fprintf( stderr, "frag %d is too far from an end!\n", currFrag->iid);
	  // walk failed
	  fprintf( stderr, "stopping walk because of condition 2.0 (fragment out of bounds)\n");
	  fprintf( stderr, "----------------\n");
	  *completeOverlapPath = -2;

	  if (setFrag)
	    ResetFrag( currFrag, currFragOriginalContig, currFragOriginalOffset5p, currFragOriginalOffset3p);
	  return walkedChunks;
	}

	// get the CIFragT for succ1
	succ1 = ifrag + increment;
	info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, succ1);

	// we used to assert, but now we just end the walk gracefully
	// assert(info->set);
	if (info->set == 0)
	{
	  fprintf( stderr, "succ1 iid %d is not set!\n", succ1);
	  // walk failed
	  fprintf( stderr, "stopping walk because of condition 7.0 (frag not set)\n");
	  fprintf( stderr, "----------------\n");
	  *completeOverlapPath = -7;
	  if (setFrag)
	    ResetFrag( currFrag, currFragOriginalContig, currFragOriginalOffset5p, currFragOriginalOffset3p);
	  return walkedChunks;
	}

	fragSucc1 = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

	// check to see that frag is within bounds
	if (FragOutOfBounds( fragSucc1, fragSucc1->cid))
	{
	  fprintf( stderr, "frag %d is too far from an end!\n", fragSucc1->iid);
	  // walk failed
	  fprintf( stderr, "stopping walk because of condition 2.1 (fragment out of bounds)\n");
	  fprintf( stderr, "----------------\n");
	  *completeOverlapPath = -2;
	  if (setFrag)
	    ResetFrag( currFrag, currFragOriginalContig, currFragOriginalOffset5p, currFragOriginalOffset3p);
	  return walkedChunks;
	}

	// even if this is true we still check below to see if rchunk overlaps properly with the current contig
	if (fragSucc1->contigID == rchunk->id)
	  walkDone = TRUE;
	
	// try to overlap with the immediate successor
	if (currFrag->contigID != fragSucc1->contigID)
	{
	  Overlap *tempOlap = OverlapContainingContigs( currFrag, fragSucc1, tempOlapOrientation,
													lchunk, rchunk);
	  fprintf( stderr, "just called OverlapContainingContigs with currFrag: %d, fragSucc1: %d\n",
			   currFrag->iid, fragSucc1->iid);

	  // potentially set an overlap of 20 and move on
	  if ( tempOlap == NULL && fixOverlapFailures == TRUE)
	  {
		ChunkInstanceT *current, *successor;

		fprintf( stderr, "adjusting walk because of condition 1 (failed overlap - pairwise)\n");
		fprintf( stderr, "----------------\n");

		current = GetGraphNode( ScaffoldGraph->ContigGraph, currFrag->contigID);
		successor = GetGraphNode( ScaffoldGraph->ContigGraph, fragSucc1->contigID);

		dummyOlap.length = 20;
		dummyOlap.begpos = current->bpLength.mean - 20;   // ahang
		dummyOlap.endpos = successor->bpLength.mean - 20; // bhang
		tempOlap = &dummyOlap;
	  }

	  if ( tempOlap == NULL )
	  {
		// walk failed
		fprintf( stderr, "stopping walk because of condition 1 (failed overlap - pairwise)\n");
		fprintf( stderr, "----------------\n");
		if (1)
		{
		  ChunkInstanceT *contig;
		  contig = GetGraphNode( ScaffoldGraph->ContigGraph, currFrag->contigID);
		  dumpContigInfo( contig );
		  contig = GetGraphNode( ScaffoldGraph->ContigGraph, fragSucc1->contigID);
		  dumpContigInfo( contig );
		}
		*completeOverlapPath = -1;
		if (setFrag)
		  ResetFrag( currFrag, currFragOriginalContig, currFragOriginalOffset5p, currFragOriginalOffset3p);
		return walkedChunks;
	  }
	  else if (0) // ((tempOlap->begpos > 0 && tempOlap->endpos < 0) ||   // the first chunk contains the second
		// (tempOlap->begpos < 0 && tempOlap->endpos > 0) )   // the second chunk contains the first
	  {
		// walk 
		fprintf( stderr, "stopping walk because of condition 4\n");
		fprintf( stderr, "----------------\n");
		*completeOverlapPath = -4;
		if (setFrag)
		  ResetFrag( currFrag, currFragOriginalContig, currFragOriginalOffset5p, currFragOriginalOffset3p);
		return walkedChunks;
	  }		
	  else
	  {
		// we used to have some code here for checking, see #ifdef CHECKCODE at bottom of file

		// get the goods on the node we just overlapped with
		succNode = GetGraphNode( ScaffoldGraph->RezGraph, fragSucc1->contigID);

		// make sure we don't rearrange contigs in the same scaffold
		if (succNode->scaffoldID != NULLINDEX &&
			succNode->scaffoldID == lchunk->scaffoldID &&
			succNode->id != rchunk->id &&
			succNode->bpLength.mean > 2000) // should check for links
		{
		  // don't walk unique chunks in some scaffold with lots of mates
		  fprintf( stderr, "stopping walk because of condition 3 (contig in this scaffold)\n");
		  fprintf( stderr, "succNode->id: %d, succNode->scaffoldID: %d, rchunk->id: %d\n", 
				   succNode->id, succNode->scaffoldID, rchunk->id);
		  fprintf( stderr, "----------------\n");
		  *completeOverlapPath = -3;
		  if (setFrag)
		    ResetFrag( currFrag, currFragOriginalContig, currFragOriginalOffset5p, currFragOriginalOffset3p);
		  return walkedChunks;
		}

		if (1)
		{
		  // don't take contigs bigger than a certain size from other scaffolds
		  if (succNode->scaffoldID != NULLINDEX &&
			  succNode->scaffoldID != lchunk->scaffoldID && 
			  succNode->bpLength.mean > 2000 &&            // should check for links, not length
			                                               // and steal big singletons, but not big others
			  GetGraphNode(ScaffoldGraph->ScaffoldGraph, succNode->scaffoldID)->info.Scaffold.numElements > 1)
		  {
			// don't walk unique chunks in some scaffold with lots of mates
			fprintf( stderr, "stopping walk because of condition 8 (contig in other scaffold)\n");
			fprintf( stderr, "succNode->id: %d, succNode->scaffoldID: %d, succNode->bpLength.mean: %f, rchunk->id: %d\n", 
					 succNode->id, succNode->scaffoldID, succNode->bpLength.mean, rchunk->id);
			fprintf( stderr, "----------------\n");
			*completeOverlapPath = -8;
			if (setFrag)
			  ResetFrag( currFrag, currFragOriginalContig, currFragOriginalOffset5p, currFragOriginalOffset3p);
			return walkedChunks;
		  }
		}
		
		// update where we are in the walk, tempOlap->begpos is the ahang
		// resetting at the same time who the successor node is
		succNodeBeg = succNodeBeg + tempOlap->begpos;
		succNodeEnd = succNodeBeg + (int) succNode->bpLength.mean;
		
		// save positions of the CIs
		if (1) // (fragSucc1->contigID != rchunk->id)
		{
		  //if (walkedChunks == NULL)
		  //{
		  //currentWalkedChunk = walkedChunks = (ChunkInsertInfoT *) malloc ( sizeof(ChunkInsertInfoT));
		  //}
		  //else
		  {
			ChunkInsertInfoT *newWalkedChunk = (ChunkInsertInfoT *) malloc ( sizeof(ChunkInsertInfoT));
			if (newWalkedChunk == NULL)
			{
			  fprintf( stderr, "malloc of newWalkedChunk failed...\n");
			  exit(1);
			}
			currentWalkedChunk->next = newWalkedChunk;
			currentWalkedChunk = currentWalkedChunk->next;
		  }
		  currentWalkedChunk->next = NULL;
		  currentWalkedChunk->contigID = fragSucc1->contigID;
		  currentWalkedChunk->scaffoldID = lchunk->scaffoldID;
		  currentWalkedChunk->insertOrder = insertCounter++;
		  currentWalkedChunk->hasMoved = FALSE;

#if DEBUG_GAP_WALKER > -1		  
		  fprintf( stderr, "just set %d as currentWalkedChunk->contigID\n", currentWalkedChunk->contigID);
#endif
		}

		if (*tempOlapOrientation == AB_AB || *tempOlapOrientation == BA_AB)
		{
		  // fprintf( stderr, "1 succNodeBeg: %f, succNodeEnd: %f\n", succNodeBeg, succNodeEnd);		  
		  currentWalkedChunk->aEndOffset.mean = succNodeBeg;
		  currentWalkedChunk->aEndOffset.variance = min( lchunk->offsetAEnd.variance, lchunk->offsetBEnd.variance) +
			ComputeFudgeVariance( succNodeBeg - min( lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean));
		  currentWalkedChunk->bEndOffset.mean = succNodeEnd;
		  currentWalkedChunk->bEndOffset.variance = min( lchunk->offsetAEnd.variance, lchunk->offsetBEnd.variance) +
			ComputeFudgeVariance( succNodeEnd - min( lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean));
		}
		else
		{
		  // fprintf( stderr, "2 succNodeBeg: %f, succNodeEnd: %f\n", succNodeBeg, succNodeEnd);
		  currentWalkedChunk->aEndOffset.mean = succNodeEnd;
		  currentWalkedChunk->aEndOffset.variance = min( lchunk->offsetAEnd.variance, lchunk->offsetBEnd.variance) +
			ComputeFudgeVariance( succNodeEnd - min( lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean));
		  currentWalkedChunk->bEndOffset.mean = succNodeBeg;
		  currentWalkedChunk->bEndOffset.variance = min( lchunk->offsetAEnd.variance, lchunk->offsetBEnd.variance) +
			ComputeFudgeVariance( succNodeBeg - min( lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean));
		}

		// temp checking code
		if (1)
		{
		  fprintf( stderr, "adding contig %d to walk\n", currentWalkedChunk->contigID);
		  fprintf( stderr, "(aEndOffset.mean, bEndOffset.mean): (%f, %f)\n", 
				   currentWalkedChunk->aEndOffset.mean, currentWalkedChunk->bEndOffset.mean);
		  fprintf( stderr, "(aEndOffset.variance, bEndOffset.variance): (%f, %f)\n", 
				   currentWalkedChunk->aEndOffset.variance, currentWalkedChunk->bEndOffset.variance);
		  fprintf( stderr, "\n");
		}
		
		// check overlaps other than pairwise
		if (walkDone)  // don't check if we are at rchunk, since its position is not correct
		{
		  ChunkInsertInfoT *firstWalkedChunk, *lastWalkedChunk;
		  ChunkInsertInfoT* sortedWalkedChunks = SortWalk(walkedChunks, &firstWalkedChunk, &lastWalkedChunk);
		  int numChunks = 1;
		  
		  while (walkedChunks->next != NULL)
		  {
			numChunks++;
			walkedChunks = walkedChunks->next;
		  }
		  
		  // walked Chunks is now an array sorted by min contig position
		  // but the links preserve the original walk order, look at insertOrder to find first
		  walkedChunks = sortedWalkedChunks;
		  currentWalkedChunk = lastWalkedChunk;

		  if ( CheckWalkOverlaps( walkedChunks, numChunks, fixOverlapFailures)) // , currentWalkedChunk))
		  {
			// walk failed
			fprintf( stderr, "stopping walk because of condition 5 (failed overlap - stack)\n");
			fprintf( stderr, "----------------\n");
			*completeOverlapPath = -5;
			{
			  ChunkInsertInfoT* chunks = walkedChunks;
			  while (chunks->next != NULL)
			  {
				dumpContigInfo ( GetGraphNode( ScaffoldGraph->ContigGraph, chunks->contigID));
				chunks = chunks->next;
			  }
			  dumpContigInfo ( GetGraphNode( ScaffoldGraph->ContigGraph, chunks->contigID));
			}
			if (setFrag)
			  ResetFrag( currFrag, currFragOriginalContig, currFragOriginalOffset5p, currFragOriginalOffset3p);
			return walkedChunks;
		  }
 /* 
		}
		
		if (walkDone)  // i.e., we are at rchunk
		{ 
 */
		  // check to make sure none of the walked chunks contains the rchunk
		  if ( CheckRchunkContainment( walkedChunks, lastWalkedChunk))
		  {
			// walk failed
			fprintf( stderr, "stopping walk because of condition 6 (rchunk contained)\n");
			fprintf( stderr, "----------------\n");
			*completeOverlapPath = -6;
			if (setFrag)
			  ResetFrag( currFrag, currFragOriginalContig, currFragOriginalOffset5p, currFragOriginalOffset3p);
			return walkedChunks;
		  }

		  *rchunkNewAEnd = lastWalkedChunk->aEndOffset.mean;
		  *rchunkNewBEnd = lastWalkedChunk->bEndOffset.mean;
		  
		  // we need this info for fixing the case where the rchunk has slid past the lchunk
		  rchunkOverlap->length = tempOlap->length;
		  rchunkOverlap->begpos = tempOlap->begpos;   // ahang
		  rchunkOverlap->endpos = tempOlap->endpos;   // bhang
		}  
	  }	  
//	  free( tempOlap );     is a local in dpcompare
	}
	if (setFrag)
	  ResetFrag( currFrag, currFragOriginalContig, currFragOriginalOffset5p, currFragOriginalOffset3p);
	ifrag += increment;
  }

#if DEBUG_GAP_WALKER > -1		  
  fprintf( stderr, "succNodeBeg = %d, succNodeEnd = %d (current pos of rchunk)\n", 
		   (int) rchunk->offsetAEnd.mean, (int) rchunk->offsetBEnd.mean);
  // walk succeeded
  fprintf( stderr, "walk succeeded\n");
#endif

  *completeOverlapPath = 1;
  return walkedChunks;
}

ChunkInsertInfoT* SortWalk( ChunkInsertInfoT* walkedChunks, ChunkInsertInfoT** firstWalkedChunk,
							ChunkInsertInfoT** lastWalkedChunk)
{
  ChunkInsertInfoT* oldChunks;
  ChunkInsertInfoT* newChunks = NULL;
  int i, j, numChunks = 1;
  
  oldChunks = walkedChunks;
  while (oldChunks->next != NULL)
  {
	numChunks++;
	oldChunks = oldChunks->next;
  }

  newChunks = (ChunkInsertInfoT *) malloc( numChunks * sizeof( ChunkInsertInfoT ));
  
  // make a copy of walkedChunks
  oldChunks = walkedChunks;
  for (i = 0; i < numChunks; i++)
  {
	newChunks[i] = oldChunks[0];
	newChunks[i].next = NULL;
	oldChunks = oldChunks->next;
  }
  
  // sort by min {a,b}EndOffset
  qsort( newChunks, numChunks, sizeof(ChunkInsertInfoT), &compChunks);
  
  // now preserve the walk order by setting pointer to the next contig, ie, the one inserted immediately after
  for (i = 0; i < numChunks; i++)
	for (j = 0; j < numChunks; j++)
	{
	  //fprintf( stderr, "newChunks[i].insertOrder + 1: %d, newChunks[j].insertOrder: %d\n", 
	  //	   newChunks[i].insertOrder + 1, newChunks[j].insertOrder);
	  if (newChunks[i].insertOrder + 1 == newChunks[j].insertOrder)
	  {
		// fprintf( stderr, "assigning %p to newChunks[%d].next\n", &newChunks[j], i);
		newChunks[i].next = &newChunks[j];	
	  }
	}

  for (i = 0; i < numChunks; i++)
  {
	if (newChunks[i].insertOrder == 0)
	  *firstWalkedChunk = &newChunks[i];
	if (newChunks[i].insertOrder == numChunks - 1)
	  *lastWalkedChunk = &newChunks[i];
  }
  
  return newChunks;
}


void ResetFrag( CIFragT *currFrag, int currFragOriginalContig, 
	       double currFragOriginalOffset5p, double currFragOriginalOffset3p)
{
  InfoByIID *info;

  currFrag->contigID = currFragOriginalContig;
  currFrag->contigOffset5p.mean = currFragOriginalOffset5p;
  currFrag->contigOffset3p.mean = currFragOriginalOffset3p;
  info = GetInfoByIID( ScaffoldGraph->iidToFragIndex, currFrag->iid);
  SetVA_CIFragT( ScaffoldGraph->CIFrags, info->fragIndex, currFrag);
}

// walkedContigs is an array as well as a linked list here
int CheckWalkOverlaps( ChunkInsertInfoT *walkedContigs,
                       int numChunks, int fixOverlapFailures) 
{
  Overlap *tempOlap1;
  ChunkInstanceT *contig1, *contig2;
  ChunkOrientationType overlapOrientation;
  int computedAhang;
  char *seq1, *seq2;
  int min_ahang, max_ahang;
  double erate, thresh;
  int minlen;
  double computedOverlap;
  // int fixStackOverlapFailure = TRUE;
  // int numChunks = 1;
  int i, j, k;
  int positiveBhang;
  
  // while (firstContig->next != NULL)
  // {
  // numChunks++;
  // firstContig = firstContig->next;
  // }

  // check all overlaps
  for (i = 0; i < numChunks; i++)
  {
	contig1 = GetGraphNode(ScaffoldGraph->ContigGraph, walkedContigs[i].contigID);
	positiveBhang = FALSE;

	// fprintf( stderr, "\n*****\n start CheckWalkOverlaps: walkedContigs[%d]: %d (%f, %f)\n", 
	//	 i, walkedContigs[i].contigID, walkedContigs[i].aEndOffset.mean, walkedContigs[i].bEndOffset.mean);
	  
	for (j = 0; j < numChunks && !positiveBhang; j++)
	{
	  if (i == j)
		continue;

	  // fprintf( stderr, "CheckWalkOverlaps j loop: walkedContigs[%d]: %d\n", j, walkedContigs[j].contigID);

	  // is walkedContigs[i] to the right of walkedContigs[j]?  if so, skip walkedContigs[j]
	  if (min (walkedContigs[i].aEndOffset.mean, walkedContigs[i].bEndOffset.mean) >
		  min (walkedContigs[j].aEndOffset.mean, walkedContigs[j].bEndOffset.mean))
	  {
		// fprintf( stderr, "skipping walkedContigs[%d]: %d\n", j, walkedContigs[j].contigID);
		continue;
	  }
	  
	  // see if we have positive bhang
	  if (max (walkedContigs[i].aEndOffset.mean, walkedContigs[i].bEndOffset.mean) < 
		  max (walkedContigs[j].aEndOffset.mean, walkedContigs[j].bEndOffset.mean))
		positiveBhang = TRUE;
	  
	  // have already done pairwise so skip ahead one if contigs are adjacent as indicated by insertOrder
	  if ( walkedContigs[i].insertOrder + 1 == walkedContigs[j].insertOrder &&
		   !walkedContigs[i].hasMoved && !walkedContigs[j].hasMoved)
	  {
		// fprintf( stderr, "checked already, skipping walkedContigs[%d]: %d\n", j, walkedContigs[j].contigID);
		continue;
	  }
	  
	  contig2 = GetGraphNode(ScaffoldGraph->ContigGraph, walkedContigs[j].contigID);
	
	  //fprintf( stderr, "CheckWalkOverlaps: checking contig %d (%f, %f) vs %d (%f, %f)\n", 
	  //	   walkedContigs[i].contigID, walkedContigs[i].aEndOffset.mean, walkedContigs[i].bEndOffset.mean,
	  //	   walkedContigs[j].contigID, walkedContigs[j].aEndOffset.mean, walkedContigs[j].bEndOffset.mean);

      // determine if contigs overlap
      computedOverlap = max ( walkedContigs[i].aEndOffset.mean, walkedContigs[i].bEndOffset.mean) -	
		min ( walkedContigs[j].aEndOffset.mean, walkedContigs[j].bEndOffset.mean);
	  
	  // fprintf( stderr, "computedOverlap: %f\n", computedOverlap);

      if (computedOverlap < (int) contig2->bpLength.mean)
		positiveBhang = TRUE;
	  
      if ( computedOverlap > 60 )
      {      
		// figure out orientation
		if ( walkedContigs[i].aEndOffset.mean < walkedContigs[i].bEndOffset.mean)
		{
		  if ( walkedContigs[j].aEndOffset.mean < walkedContigs[j].bEndOffset.mean)
			overlapOrientation = AB_AB;
		  else
			overlapOrientation = AB_BA;
		}
		else
		{
		  if ( walkedContigs[j].aEndOffset.mean < walkedContigs[j].bEndOffset.mean)
			overlapOrientation = BA_AB;
		  else
			overlapOrientation = BA_BA;
		}
		
		// compute expected ahang
		computedAhang = min ( walkedContigs[j].aEndOffset.mean, walkedContigs[j].bEndOffset.mean) -
		  min (walkedContigs[i].aEndOffset.mean, walkedContigs[i].bEndOffset.mean);
		min_ahang = computedAhang - 60;
		max_ahang = computedAhang + 60;
		
		if(consensus1 == NULL)
		{
		  consensus1 = CreateVA_char(2048);
		  consensus2 = CreateVA_char(2048);
		  quality1 = CreateVA_char(2048);
		  quality2 = CreateVA_char(2048);
		}
		
		// Get the consensus sequences for both contigs from the Store
		GetConsensus(ScaffoldGraph->RezGraph, contig1->id, consensus1, quality1);
		GetConsensus(ScaffoldGraph->RezGraph, contig2->id, consensus2, quality2);
		
		seq1 = Getchar(consensus1, 0);
		seq2 = Getchar(consensus2, 0);
		
		erate = 2 * CGW_DP_ERATE;
		thresh = CGW_DP_THRESH;
		minlen = CGW_DP_MINLEN;
		
		//fprintf( stderr, "CheckWalkOverlaps: contig1: %d, contig2: %d, computedAhang: %d, orientation: %c\n",
		//	 contig1->id, contig2->id, computedAhang, (char) overlapOrientation);
		
		// tempOlap1 points to a static down inside of DP_Compare
		tempOlap1 = OverlapSequences( seq1, seq2, overlapOrientation, min_ahang, max_ahang, 
									  erate, thresh, minlen, AS_FIND_ALIGN);
		
		if (tempOlap1 == NULL && fixOverlapFailures == TRUE)
		{
		  int newPos, delta;
		  ChunkInsertInfoT *tempContig = &walkedContigs[j];

		  fprintf( stderr, "fixing stack overlap failure\n");

		  // set the position of walkedContig2 so there's a 20 base pair overlap with walkedContig1
		  newPos = max( walkedContigs[i].aEndOffset.mean,  walkedContigs[i].bEndOffset.mean) - 20;
		  delta = newPos - min( tempContig->aEndOffset.mean, tempContig->bEndOffset.mean);
		  fprintf( stderr, "moving contig %d from (%f, %f) to (%f, %f)\n",
				   tempContig->contigID,
				   tempContig->aEndOffset.mean, tempContig->bEndOffset.mean,
				   tempContig->aEndOffset.mean + delta, tempContig->bEndOffset.mean + delta);			
		  tempContig->aEndOffset.mean += delta;
		  tempContig->bEndOffset.mean += delta;
		  tempContig->hasMoved = TRUE;
		  // also have to adjust the postion of everybody else in the walk

		  for (k = j + 1; k < numChunks; k++)
		  {
            // only move the chunks inserted after walkedContigs[j], not all those to the right of it
			if ( walkedContigs[k].insertOrder > walkedContigs[j].insertOrder)  
			{
			  fprintf( stderr, "moving contig %d from (%f, %f) to (%f, %f)\n",
					   walkedContigs[k].contigID,
					   walkedContigs[k].aEndOffset.mean, walkedContigs[k].bEndOffset.mean,
					   walkedContigs[k].aEndOffset.mean + delta, walkedContigs[k].bEndOffset.mean + delta);			
			  walkedContigs[k].aEndOffset.mean += delta;
			  walkedContigs[k].bEndOffset.mean += delta;
			  walkedContigs[k].hasMoved = TRUE;
			}
		  }
		}
		else if (tempOlap1 == NULL)
		{
		  // dumpContigInfo(contig1);
		  // dumpContigInfo(contig2);
		  return 1;
		}
      }
    }
  }
  return 0;
}

// // walkedContigs is an array as well as a linked list here
// int CheckWalkOverlaps_old( ChunkInsertInfoT *walkedContigs) // , ChunkInsertInfoT *lastWalkedContig)
// {
//   Overlap *tempOlap1;
//   ChunkInsertInfoT *firstContig = walkedContigs;
//   ChunkInsertInfoT *walkedContig1, *walkedContig2;;
//   ChunkInstanceT *contig1, *contig2;
//   ChunkOrientationType overlapOrientation;
//   int computedAhang, negativeBhang;
//   char *seq1, *seq2;
//   int min_ahang, max_ahang;
//   double erate, thresh;
//   int minlen;
//   double computedOverlap;
//   int fixStackOverlapFailure = TRUE;
//   int numChunks = 1;
//   
//   while (walkedContigs->next != NULL)
//   {
// 	numChunks++;
// 	walkedContigs = walkedContigs->next;
//   }
// 
//   walkedContig1 = firstContig;
// 
//   // check all overlaps
//   // while ( walkedContig1->next != lastWalkedContig )
//   while ( walkedContig1 != NULL )
//   {
//     contig1 = GetGraphNode(ScaffoldGraph->ContigGraph, walkedContig1->contigID);
// 	if (walkedContig1->next == NULL)
// 	  break;
// 	
//     fprintf( stderr, "\n*****\nCheckWalkOverlaps: contig1: %d (%f, %f)\n", 
// 			 contig1->id, walkedContig1->aEndOffset.mean, walkedContig1->bEndOffset.mean);
// 	
//     // have already done pairwise so skip ahead one if contigs are adjacent as indicated by insertOrder
//     walkedContig2 = walkedContig1->next;
// 	if (walkedContig1->insertOrder == walkedContig2->insertOrder - 1)
// 	{
// 	  if (walkedContig2->next == NULL)
// 		break;
// 	  walkedContig2 = walkedContig2->next;
// 	}
// 
// 	// did we move walkedContig1 to the right of walkedContig2?  if so, skip walkedContig2
// 	while (min (walkedContig1->aEndOffset.mean, walkedContig1->bEndOffset.mean) > 
// 		   min (walkedContig2->aEndOffset.mean, walkedContig2->bEndOffset.mean))
// 	{
// 	  if (walkedContig2->next == NULL)
// 		break;
// 	  walkedContig2 = walkedContig2->next;
// 	}
// 	if (min (walkedContig1->aEndOffset.mean, walkedContig1->bEndOffset.mean) > 
// 		min (walkedContig2->aEndOffset.mean, walkedContig2->bEndOffset.mean))
// 	  break;
// 
//     contig2 = GetGraphNode(ScaffoldGraph->ContigGraph, walkedContig2->contigID);
// 	
//     // determine if bhang is positive
//     computedOverlap = max ( walkedContig1->aEndOffset.mean, walkedContig1->bEndOffset.mean) -
//       min ( walkedContig2->aEndOffset.mean, walkedContig2->bEndOffset.mean);
//     if (computedOverlap < (int) contig2->bpLength.mean)
//       negativeBhang = FALSE;
//     else
//       negativeBhang = TRUE;
// 	
//     fprintf( stderr, "CheckWalkOverlaps: outer loop, contig2: %d (%f, %f)\n", 
// 			 contig2->id, walkedContig2->aEndOffset.mean, walkedContig2->bEndOffset.mean);
//     fprintf( stderr, "computedOverlap: %f, contig2->bpLength.mean: %d, negativeBhang: %d\n", 
// 			 computedOverlap, (int) contig2->bpLength.mean, negativeBhang);
//     do
//     {
//       // determine if contigs overlap
//       computedOverlap = max ( walkedContig1->aEndOffset.mean, walkedContig1->bEndOffset.mean) -
// 		min ( walkedContig2->aEndOffset.mean, walkedContig2->bEndOffset.mean);
// 	  
//       fprintf( stderr, "CheckWalkOverlaps: inner loop, contig2: %d, computedOverlap: %d\n", 
// 			   contig2->id, computedOverlap);
// 
//       if (computedOverlap < (int) contig2->bpLength.mean)
// 		negativeBhang = FALSE;
//       else
// 		negativeBhang = TRUE;
// 	  
//       if ( computedOverlap > 60 )
//       {      
// 		// figure out orientation
// 		if ( walkedContig1->aEndOffset.mean < walkedContig1->bEndOffset.mean)
// 		{
// 		  if ( walkedContig2->aEndOffset.mean < walkedContig2->bEndOffset.mean)
// 			overlapOrientation = AB_AB;
// 		  else
// 			overlapOrientation = AB_BA;
// 		}
// 		else
// 		{
// 		  if ( walkedContig2->aEndOffset.mean < walkedContig2->bEndOffset.mean)
// 			overlapOrientation = BA_AB;
// 		  else
// 			overlapOrientation = BA_BA;
// 		}
// 		
// 		// compute expected ahang
// 		computedAhang = min ( walkedContig2->aEndOffset.mean, walkedContig2->bEndOffset.mean) -
// 		  min (walkedContig1->aEndOffset.mean, walkedContig1->bEndOffset.mean);
// 		min_ahang = computedAhang - 60;
// 		max_ahang = computedAhang + 60;
// 		
// 		if(consensus1 == NULL)
// 		{
// 		  consensus1 = CreateVA_char(2048);
// 		  consensus2 = CreateVA_char(2048);
// 		  quality1 = CreateVA_char(2048);
// 		  quality2 = CreateVA_char(2048);
// 		}
// 		
// 		// Get the consensus sequences for both contigs from the Store
// 		GetConsensus(ScaffoldGraph->RezGraph, contig1->id, consensus1, quality1);
// 		GetConsensus(ScaffoldGraph->RezGraph, contig2->id, consensus2, quality2);
// 		
// 		seq1 = Getchar(consensus1, 0);
// 		seq2 = Getchar(consensus2, 0);
// 		
// 		erate = 2 * CGW_DP_ERATE;
// 		thresh = CGW_DP_THRESH;
// 		minlen = CGW_DP_MINLEN;
// 		
// 		fprintf( stderr, "CheckWalkOverlaps: contig1: %d, contig2: %d, computedAhang: %d, orientation: %c\n",
// 				 contig1->id, contig2->id, computedAhang, (char) overlapOrientation);
// 		
// 		// tempOlap1 points to a static down inside of DP_Compare
// 		tempOlap1 = OverlapSequences( seq1, seq2, overlapOrientation, min_ahang, max_ahang, 
// 									  erate, thresh, minlen, AS_FIND_ALIGN);
// 		
// 		if (tempOlap1 == NULL && fixStackOverlapFailure == TRUE)
// 		{
// 		  int newPos, delta;
// 		  ChunkInsertInfoT *tempContig = walkedContig2;
// 
// 		  fprintf( stderr, "fixing stack overlap failure\n");
// 
// 		  // set the position of walkedContig2 so there's a 20 base pair overlap with walkedContig1
// 		  newPos = max( walkedContig1->aEndOffset.mean,  walkedContig1->bEndOffset.mean) - 20;
// 		  delta = newPos - min( tempContig->aEndOffset.mean, tempContig->bEndOffset.mean);
// 		  fprintf( stderr, "moving contig %d from (%f, %f) to (%f, %f)\n",
// 				   tempContig->contigID,
// 				   tempContig->aEndOffset.mean, tempContig->bEndOffset.mean,
// 				   tempContig->aEndOffset.mean + delta, tempContig->bEndOffset.mean + delta);			
// 		  tempContig->aEndOffset.mean += delta;
// 		  tempContig->bEndOffset.mean += delta;
// 		  // also have to adjust the postion of everybody else in the walk
// 		  do
// 		  {
// 			tempContig = tempContig->next;
//             // only move the chunks inserted after walkedContig2, not all those to the right of it
// 			if (tempContig->insertOrder > walkedContig2->insertOrder)  
// 			{
// 			  fprintf( stderr, "moving contig %d from (%f, %f) to (%f, %f)\n",
// 					   tempContig->contigID,
// 					   tempContig->aEndOffset.mean, tempContig->bEndOffset.mean,
// 					   tempContig->aEndOffset.mean + delta, tempContig->bEndOffset.mean + delta);			
// 			  tempContig->aEndOffset.mean += delta;
// 			  tempContig->bEndOffset.mean += delta;
// 			}
// 		  } while (tempContig->next != NULL);
// 		}
// 		else if (tempOlap1 == NULL)
// 		{
// 		  // dumpContigInfo(contig1);
// 		  // dumpContigInfo(contig2);
// 		  return 1;
// 		}
//       }
// 	  if (walkedContig2->next != NULL)
// 	  {
// 		walkedContig2 = walkedContig2->next;
// 		contig2 = GetGraphNode(ScaffoldGraph->ContigGraph, walkedContig2->contigID);
// 	  }
// 	  else
// 		negativeBhang = FALSE;
//     } while (negativeBhang);
//     walkedContig1 = walkedContig1->next;
//   }
//   return 0;
// }



int CheckRchunkContainment( ChunkInsertInfoT *walkedChunks,
                            ChunkInsertInfoT *lastWalkedChunk)
{
  ChunkInsertInfoT *previousChunk = walkedChunks;

  // check all overlaps
  while ( previousChunk != lastWalkedChunk )
  {
	// check to make sure the rchunk is not contained by any chunk
	// on entry previousChunk points to lchunk
	if ( min ( lastWalkedChunk->aEndOffset.mean, lastWalkedChunk->bEndOffset.mean) >
		 min ( previousChunk->aEndOffset.mean, previousChunk->bEndOffset.mean)
		 &&
		 max ( lastWalkedChunk->aEndOffset.mean, lastWalkedChunk->bEndOffset.mean) <
		 max ( previousChunk->aEndOffset.mean, previousChunk->bEndOffset.mean))
	  return 1;
	previousChunk = previousChunk->next;
  }
  return 0;
}

// check whether a frag is within some cutoff of it's containing contig
// to do: check whether it's the appropriate end
int FragOutOfBounds(CIFragT* frag, int unitigID)
{
  double distToAEnd, distToBEnd;
  ContigT *contig;
  NodeCGW_T *unitig;
  int surrogateOffset;
  
  contig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);

  unitig = GetGraphNode( ScaffoldGraph->CIGraph, unitigID);

  if (unitig->flags.bits.isWalkSurrogate || unitig->flags.bits.isStoneSurrogate)
  {
	surrogateOffset = (int) max( unitig->offsetAEnd.mean, unitig->offsetAEnd.mean);
	// unitig = GetGraphNode( ScaffoldGraph->CIGraph, unitig->info.CI.baseID);
  }
  else
	surrogateOffset = 0;

  distToAEnd = surrogateOffset + min( frag->contigOffset5p.mean, frag->contigOffset3p.mean);
  distToBEnd = min( contig->bpLength.mean - (frag->contigOffset5p.mean + surrogateOffset), 
					   contig->bpLength.mean - (frag->contigOffset3p.mean + surrogateOffset));
  
  if ( min( distToAEnd, distToBEnd) > (double) BASES_FROM_END_CUTOFF)
  {
	fprintf( stderr, "frag->iid: %d, distance to end: %f\n", frag->iid, min( distToAEnd, distToBEnd));
	return 1;
  }
  else
	return 0;
}

// find the last fragment in this contig that is consecutive with ifrag 
int findLastLocaleFragInContig( int contigID, unsigned int currFragIid,
                                unsigned int endFragIid, int increment, 
                                int *currExtremeUnitigID)
{
  int numUnitigs, i;
  MultiAlignT *uma;
  CDS_IID_t currExtremeFragIid;

  // loop over all the unitigs (grabbing parent if necessary) looking for the extreme
  // fragment between ifrag and endFragID
  // first, step through the f_list of the contig
  // get multi-align for the contig
  uma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contigID, ScaffoldGraph->RezGraph->type == CI_GRAPH); 
  numUnitigs = GetNumIntUnitigPoss( uma->u_list );

  if (increment == 1)
    currExtremeFragIid = 0;
  else
    currExtremeFragIid = CDS_IID_MAX;

  for ( i = 0; i < numUnitigs; i++)
  {
  	IntUnitigPos *upos = GetIntUnitigPos( uma->u_list, i);
	ChunkInstanceT *unitig = GetGraphNode( ScaffoldGraph->CIGraph, upos->ident);
	ChunkInstanceT *originalUnitig;
	MultiAlignT *uma;
	int icntfrag;

	originalUnitig = unitig;
	if (unitig->flags.bits.isWalkSurrogate || unitig->flags.bits.isStoneSurrogate)
	{
	  unitig = GetGraphNode( ScaffoldGraph->CIGraph, unitig->info.CI.baseID);
	}
	uma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, ScaffoldGraph->CIGraph->type == CI_GRAPH); 
	// fprintf( stderr, "  unitig: %d\t num frags: %d\n", unitig->id, GetNumIntMultiPoss(uma->f_list));

	// now get info on the frags in the unitig
	for (icntfrag = 0; icntfrag < GetNumIntMultiPoss(uma->f_list); icntfrag++)
	{
	  IntMultiPos *imp = GetIntMultiPos(uma->f_list, icntfrag);

	  // fprintf( stderr, "  unitig: %d\t examining frag: %d\n", unitig->id, imp->ident);
	  if (increment == 1)
	  {
		if (imp->ident >= currFragIid && imp->ident < endFragIid)
		  if (imp->ident > currExtremeFragIid)
		  {
			currExtremeFragIid = imp->ident;
			*currExtremeUnitigID = originalUnitig->id;
		  }
	  }
	  else
	  {
		if (imp->ident <= currFragIid && imp->ident > endFragIid)
		  if (imp->ident < currExtremeFragIid)
		  {
			currExtremeFragIid = imp->ident;
			*currExtremeUnitigID = originalUnitig->id;
		  }
	  }
	}
  }
  return currExtremeFragIid;
}

void CheckOrientation( ContigT* lchunk, ContigT* rchunk,
                       LengthT* rchunk_delta, 
                       ChunkOrientationType olapOrientation,
                       Overlap* rchunkOverlap)
{
  ChunkOrientationType scaffoldOrientation = XX_XX;
  int lchunkOrientation, rchunkOrientation;
  
  if (lchunk->offsetAEnd.mean < lchunk->offsetBEnd.mean)
	lchunkOrientation = 0;
  else
	lchunkOrientation = 1;

  if (rchunk->offsetAEnd.mean < rchunk->offsetBEnd.mean)
	rchunkOrientation = 0;
  else
	rchunkOrientation = 1;

  if (lchunkOrientation == 0 && rchunkOrientation == 0)
	scaffoldOrientation = AB_AB;
  else if (lchunkOrientation == 0 && rchunkOrientation == 1)
	scaffoldOrientation = AB_BA;
  else if (lchunkOrientation == 1 && rchunkOrientation == 0)
	scaffoldOrientation = BA_AB;
  else if (lchunkOrientation == 1 && rchunkOrientation == 1)
	scaffoldOrientation = BA_BA;

  if (scaffoldOrientation != olapOrientation)
	fprintf( stderr, "Warning!!!! scaffoldOrientation != olapOrientation for lchunk: %d and rchunk: %d\n",
			 lchunk->id, rchunk->id);
  else  // they have the right orientation
	return;
  
  // these situations indicate sliding and can be fixed
  if ( ! ((scaffoldOrientation == AB_BA && olapOrientation == BA_AB) ||
		  (scaffoldOrientation == BA_AB && olapOrientation == AB_BA) ||
		  (scaffoldOrientation == AB_AB && olapOrientation == BA_BA) ||
		  (scaffoldOrientation == BA_BA && olapOrientation == AB_AB)))
  {
	fprintf( stderr, "Incompatible overlap orientations detected!\n");
	assert(0);
  }

  // olapOrientation is the real one
  if (scaffoldOrientation == AB_BA && olapOrientation == BA_AB)
  {
	double rchunkFinalBEndPosMean = lchunk->offsetAEnd.mean - rchunkOverlap->endpos;
	rchunk_delta->mean = rchunkFinalBEndPosMean - rchunk->offsetBEnd.mean;
  }
  else if (scaffoldOrientation == BA_AB && olapOrientation == AB_BA)
  {
	double rchunkFinalAEndPosMean = lchunk->offsetBEnd.mean - rchunkOverlap->endpos;
	rchunk_delta->mean = rchunkFinalAEndPosMean - rchunk->offsetAEnd.mean;
  }
  else if (scaffoldOrientation == AB_AB && olapOrientation == BA_BA)
  {
	double rchunkFinalAEndPosMean = lchunk->offsetAEnd.mean - rchunkOverlap->endpos;
	rchunk_delta->mean = rchunkFinalAEndPosMean - rchunk->offsetAEnd.mean;
  }
  else if (scaffoldOrientation == BA_BA && olapOrientation == AB_AB)
  {
	double rchunkFinalBEndPosMean = lchunk->offsetBEnd.mean - rchunkOverlap->endpos;
	rchunk_delta->mean = rchunkFinalBEndPosMean - rchunk->offsetBEnd.mean;
  }
}

// adjust the position of the right contig to where it should be after walk
// this prevents contigs inserted during the walk from being contained in the destination contig
// we overstimate the distance here just to be safe, the positions are adjusted again after the walk
void AdjustRightContigPos( ContigT* lchunk, ContigT* rchunk, LengthT rchunk_delta, CIScaffoldT* scaff)
{
  LengthT rchunkNewAEndOffset, rchunkNewBEndOffset;
  float64 lchunkMaxVariance = max( lchunk->offsetAEnd.variance, lchunk->offsetBEnd.variance);
  float64 lchunkMaxMean     = max( lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean);
  float64 rchunkMinVariance = min( rchunk->offsetAEnd.variance, rchunk->offsetBEnd.variance);
  float64 rchunkMinMean     = min( rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean);
  int startingContigIndex, scaffoldIndex;
  
  // rchunk_delta.mean has already been computed
  rchunk_delta.variance = lchunkMaxVariance                                    // the variance at the left of the gap
	+ ComputeFudgeVariance( rchunkMinMean + rchunk_delta.mean - lchunkMaxMean) // plus the variance across the gap
	- rchunkMinVariance;                                                       // minus the current variance for the rchunk

  // we used to just add a delta to rchunk and all contigs towards the B end
  // but if the lchunk and rchunk abut then the rchunk min might slide past the lchunk and violate the invariant 
  // now remove rchunk from the scaffold and reinsert him in his new position
  rchunkNewAEndOffset.mean = rchunk->offsetAEnd.mean + rchunk_delta.mean;
  rchunkNewAEndOffset.variance = rchunk->offsetAEnd.variance + rchunk_delta.variance;
  rchunkNewBEndOffset.mean = rchunk->offsetBEnd.mean + rchunk_delta.mean;
  rchunkNewBEndOffset.variance = rchunk->offsetBEnd.variance + rchunk_delta.variance;

  startingContigIndex = rchunk->BEndNext;
  scaffoldIndex = rchunk->scaffoldID;
  
  // take rchunk out
  RemoveCIFromScaffold( ScaffoldGraph, scaff, rchunk, FALSE);

  // now move everybody past rchunk into position
  if (startingContigIndex != NULLINDEX)
	AddDeltaToScaffoldOffsets(ScaffoldGraph,
							  scaffoldIndex,
							  startingContigIndex,  // was rchunk->id,
							  1,
							  0,
							  rchunk_delta);

  // now put rchunk back in
  InsertCIInScaffold(ScaffoldGraph, rchunk->id, scaff->id,
					 rchunkNewAEndOffset, rchunkNewBEndOffset, TRUE, NO_CONTIGGING);

  Force_Increasing_Variances_One_Scaffold( scaff->id );
}

void InsertWalkedChunks( ChunkInsertInfoT* chunksWalked,
                         ChunkInstanceT* lchunk, ChunkInstanceT* rchunk, 
                         int completeOverlapPath,
                         int *checkScaffold, int *checkScaffoldCount)
{
  // skip over lchunk, which is first entry in chunksWalked
  ChunkInsertInfoT *currentWalkedChunk = chunksWalked->next;
  CDS_CID_t rchunkId = rchunk->id;
  CDS_CID_t lchunkId = lchunk->id;
  int foundSufficientOverlap = FALSE;
  static int walkCIsInScaffold = 0;
  static int workingScaffold = -1;
  int containsDiscriminatorUniqueChunk, terminate;
  NodeCGW_T *unitig;
  
  if (completeOverlapPath != 1)
	assert(0);

  // this mean that we have only the rchunk in the walk and it's already been positioned
  if (chunksWalked->next == NULL)
	return;
  
  fprintf( stderr, "at scaffold %d, , walkCISInScaffold = %d\n", 
		   lchunk->scaffoldID, walkCIsInScaffold);

  if (lchunk->scaffoldID != workingScaffold)
  {
	walkCIsInScaffold = 0;	
	workingScaffold = lchunk->scaffoldID;
	fprintf( stderr, "at workingScaffold %d, walkCISInScaffold = %d\n", 
			 lchunk->scaffoldID, walkCIsInScaffold);
  }
  
  // if completeOverlapPath is 1, then the last chunk is the rchunk
  // if completeOverlapPath is 0, then the last chunk is just where the walk ended

  // translate coordinates for insertion from absolute to relative to lchunk->offsetAEnd
  while (currentWalkedChunk->next != NULL)
  {
	// fprintf( stderr, "moving walked chunk %d from (%f, %f) to ", currentWalkedChunk->contigID,
	//	 currentWalkedChunk->aEndOffset.mean, currentWalkedChunk->bEndOffset.mean);

	currentWalkedChunk->aEndOffset.mean -= lchunk->offsetAEnd.mean;
	currentWalkedChunk->aEndOffset.variance -= lchunk->offsetAEnd.variance;
	currentWalkedChunk->bEndOffset.mean -= lchunk->offsetAEnd.mean;
	currentWalkedChunk->bEndOffset.variance -= lchunk->offsetAEnd.variance;

	// fprintf( stderr, "(%f, %f)\n", currentWalkedChunk->aEndOffset.mean, currentWalkedChunk->bEndOffset.mean);

	currentWalkedChunk = currentWalkedChunk->next;
  }

  // reset, but skip the lchunk, which is the first member of chunksWalked
  currentWalkedChunk = chunksWalked->next;

  // do all except the last chunk
  while (currentWalkedChunk->next != NULL)
  {
	ContigT *walkedChunk = GetGraphNode(ScaffoldGraph->RezGraph, currentWalkedChunk->contigID);

	// we need to update lchunk, rchunk since they point into the contig array, which can be reallocated
	// eg, during SplitUnresolvedContig
	lchunk = GetGraphNode( ScaffoldGraph->ContigGraph, lchunkId);
	rchunk = GetGraphNode( ScaffoldGraph->ContigGraph, rchunkId);

	// don't insert chunks contained by the lchunk
	if (0)
	  if ( max( currentWalkedChunk->aEndOffset.mean, currentWalkedChunk->bEndOffset.mean) < 
		   max( lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean))
	  {
		currentWalkedChunk = currentWalkedChunk->next;
		continue;
	  }
	
	// the old way - just avoid contained contigs
	if (1)
	{
	  // we don't want to insert any chunks contained by the rchunk
	  // they make walking from the rchunk to the next gap difficult
	  // have to adjust to gap coordiante system

	  //fprintf( stderr, "lchunk->offsetAEnd.mean + min( currentWalkedChunk->{a,b}EndOffset.mean): %f\n",
	  //	   lchunk->offsetAEnd.mean + min( currentWalkedChunk->aEndOffset.mean, currentWalkedChunk->bEndOffset.mean));
	  //fprintf( stderr, "min( rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean): %f\n",
	  //	   min( rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean));
	  
	  if ( lchunk->offsetAEnd.mean + min( currentWalkedChunk->aEndOffset.mean, currentWalkedChunk->bEndOffset.mean) > 
		   min( rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean))
	  {
		currentWalkedChunk = currentWalkedChunk->next;
		continue;
	  }
	}
	
	// the new way - stop as soon as we have a chunk that overlaps enough
	if (0 && foundSufficientOverlap)
	{
		currentWalkedChunk = currentWalkedChunk->next;
		continue;
	}
	
	if ( max( currentWalkedChunk->aEndOffset.mean, currentWalkedChunk->bEndOffset.mean) 
		 - min( rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean) > 60)
	  foundSufficientOverlap = TRUE;
	  
	// fprintf( stderr, "during InsertWalkedChunks, foundSufficientOverlap:%d\n", foundSufficientOverlap);
	
	// determine whether contig should be moved or surrogated
	containsDiscriminatorUniqueChunk = FALSE;
	terminate = FALSE;
	unitig = GetGraphNode( ScaffoldGraph->CIGraph, walkedChunk->info.Contig.AEndCI);
	do 
	{
	  if (unitig->type == DISCRIMINATORUNIQUECHUNK_CGW)
	  {
		containsDiscriminatorUniqueChunk = TRUE;
		terminate = TRUE;
	  }
	  if (unitig->BEndNext != NULLINDEX)
		unitig = GetGraphNode( ScaffoldGraph->CIGraph, unitig->BEndNext);
	  else
		terminate = TRUE;
	} while (!terminate);
	

	// check on status of contig
	// first do contigs that are already scaffolded
	if (walkedChunk->scaffoldID != -1 || containsDiscriminatorUniqueChunk)
	// if (walkedChunk->type != UNRESOLVEDCHUNK_CGW)
	{
	  CIScaffoldT *walkedChunkCurrScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, walkedChunk->scaffoldID);

	  if (walkedChunk->scaffoldID != -1)
	  {
		fprintf( stderr, "1 contig %d is in scaffold %d, which has %d contigs\n", 
				 currentWalkedChunk->contigID, walkedChunk->scaffoldID, 
				 walkedChunkCurrScaffold->info.Scaffold.numElements);
		
		RemoveCIFromScaffold( ScaffoldGraph, walkedChunkCurrScaffold, walkedChunk, TRUE);
		CheckCIScaffoldTLength(ScaffoldGraph, walkedChunkCurrScaffold);

		// if we know that the scaffold only had 1 contig we should kill it
		if (walkedChunkCurrScaffold->info.Scaffold.numElements == 0)
		{
		  DeleteGraphNode( ScaffoldGraph->ScaffoldGraph, walkedChunkCurrScaffold);
		  fprintf( stderr, "Killed Scaffold %d and inserted its contig into Scaffold %d\n",
				   walkedChunkCurrScaffold->id, currentWalkedChunk->scaffoldID);
		}
	  }
	  else
	  {
		fprintf( stderr, "1 contig %d is in scaffold %d\n", 
				 currentWalkedChunk->contigID, currentWalkedChunk->scaffoldID);
	  }
	  
	  
	  
	  // the next two lines have bad effects if we try to move contigs within the same scaffold 
	  // and that scaffold gets broken up - since we're currently walking it
#if 0
	  else
	  {
		// mark the scaffold we stole contig from so it can be checked later for connectivity
		// fprintf( stderr, "*checkScaffoldCount before: %d\n", *checkScaffoldCount);
		checkScaffold[(*checkScaffoldCount)++] = walkedChunkCurrScaffold->id;
		// fprintf( stderr, "*checkScaffoldCount after: %d\n", *checkScaffoldCount);
		if (*checkScaffoldCount > MAX_SCAFFOLDS_CHECKED)
		{
		  fprintf( stderr, "Whoops! *checkScaffoldCount is now greater than MAX_SCAFFOLDS_CHECKED (%d)\n", 
				   *checkScaffoldCount, MAX_SCAFFOLDS_CHECKED);
		  assert(0);
	  }
#endif	  
	  
	  // set the chunks to their insert spot relative to lchunk
	  currentWalkedChunk->aEndOffset.mean += lchunk->offsetAEnd.mean;
	  currentWalkedChunk->aEndOffset.variance += lchunk->offsetAEnd.variance;
	  currentWalkedChunk->bEndOffset.mean += lchunk->offsetAEnd.mean;
	  currentWalkedChunk->bEndOffset.variance += lchunk->offsetAEnd.variance;
	  
	  fprintf( stderr, "inserting contig %d in scaffold %d at %f, %f\n", 
			   currentWalkedChunk->contigID, currentWalkedChunk->scaffoldID,
			   currentWalkedChunk->aEndOffset.mean, currentWalkedChunk->bEndOffset.mean);
	  fprintf( stderr, "during InsertWalkedChunks, rchunk %d: AEndOffset.mean: %f, BEndOffset.mean: %f\n",
			   rchunk->id, rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean);
	  
	  InsertCIInScaffold(ScaffoldGraph, currentWalkedChunk->contigID, currentWalkedChunk->scaffoldID,
						 currentWalkedChunk->aEndOffset, currentWalkedChunk->bEndOffset, TRUE, NO_CONTIGGING);
	}
	else
	{
	  MultiAlignT *uma;
	  int numUnitigs, i, splitCid;
	  
	  splitCid = SplitUnresolvedContig(ScaffoldGraph->RezGraph, currentWalkedChunk->contigID, 
					   NULL, FALSE);
	      
	  fprintf( stderr, "2 contig %d (surrogate: %d) is in scaffold %d, which has %d contigs\n", 
			   currentWalkedChunk->contigID, splitCid, walkedChunk->scaffoldID, 1);
	  
	  // set the chunks to their insert spot relative to lchunk
	  currentWalkedChunk->aEndOffset.mean += lchunk->offsetAEnd.mean;
	  currentWalkedChunk->aEndOffset.variance += lchunk->offsetAEnd.variance;
	  currentWalkedChunk->bEndOffset.mean += lchunk->offsetAEnd.mean;
	  currentWalkedChunk->bEndOffset.variance += lchunk->offsetAEnd.variance;
	  
	  fprintf( stderr, "inserting contig %d in scaffold %d at %f, %f\n", 
			   currentWalkedChunk->contigID, currentWalkedChunk->scaffoldID,
			   currentWalkedChunk->aEndOffset.mean, currentWalkedChunk->bEndOffset.mean);
	  fprintf( stderr, "during InsertWalkedChunks, rchunk %d: AEndOffset.mean: %f, BEndOffset.mean: %f\n",
			   rchunk->id, rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean);
	  
      // now we insert the surrogate into the scaffold
	  // remember that it has a different number
	  InsertCIInScaffold(ScaffoldGraph, splitCid, currentWalkedChunk->scaffoldID,
						 currentWalkedChunk->aEndOffset, currentWalkedChunk->bEndOffset, TRUE, NO_CONTIGGING);

	  // cycle through chunks in contig setting their walkSurrogate flag
	  uma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, splitCid, FALSE); 
	  // uma = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, splitCid);
	  numUnitigs = GetNumIntUnitigPoss( uma->u_list );

	  for ( i = 0; i < numUnitigs; i++)
	  {
		IntUnitigPos *upos = GetIntUnitigPos( uma->u_list, i);
		ChunkInstanceT *unitig = GetGraphNode( ScaffoldGraph->CIGraph, upos->ident);
		unitig->flags.bits.isSurrogate = TRUE;
		unitig->flags.bits.isWalkSurrogate = TRUE;

		fprintf( stderr, "just marked unitig %d as a walkSurrogate\n", unitig->id);
	  }
	}	
	currentWalkedChunk = currentWalkedChunk->next;
	walkCIsInScaffold++;
  }  

  if (rchunkId != rchunk->id)
  {
    fprintf( stderr, "rchunkID: " F_CID ", rchunk->id: " F_CID "\n",
             rchunkId, rchunk->id);
  }
  

  if (0)
  {
	if (completeOverlapPath == 0)
	{
	  //InsertCIInScaffold(ScaffoldGraph, currentWalkedChunk->contigID, currentWalkedChunk->scaffoldID,
	  //			   currentWalkedChunk->aEndOffset, currentWalkedChunk->bEndOffset, TRUE, NO_CONTIGGING);
	}  
	else
	  // check where the rchunk is now versus where we want to put it
	  // that is the amount we want to adjust the gap size by
	{
	  LengthT delta;	
	  float64 lchunkMaxVariance = max( lchunk->offsetAEnd.variance, lchunk->offsetBEnd.variance);
	  float64 lchunkMaxMean     = max( lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean);
	  float64 rchunkMinVariance = min( rchunk->offsetAEnd.variance, rchunk->offsetBEnd.variance);
	  float64 rchunkMinMean     = min( rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean);
	  int currentChunkBegPos    = min( currentWalkedChunk->aEndOffset.mean, currentWalkedChunk->bEndOffset.mean);
	  
	  delta.mean = currentChunkBegPos - rchunkMinMean;
	  delta.variance = lchunkMaxVariance                                   // the variance at the left of the gap
		+ ComputeFudgeVariance( currentChunkBegPos - lchunkMaxMean)        // plus the variance across the gap
		- rchunkMinVariance;                                               // the current variance for the rchunk
	  
	  fprintf( stderr, "in InsertWalkedChunks, lchunkMaxVariance = %f\n", lchunkMaxVariance);
	  fprintf( stderr, "in InsertWalkedChunks, currentChunkBegPos = %d\n", currentChunkBegPos);
	  fprintf( stderr, "in InsertWalkedChunks, lchunkMaxMean = %f\n", lchunkMaxMean);
	  fprintf( stderr, "in InsertWalkedChunks, ComputeFudgeVariance( currentChunkBegPos - lchunkMaxMean) = %f\n", 
			   ComputeFudgeVariance( currentChunkBegPos - lchunkMaxMean));
	  fprintf( stderr, "in InsertWalkedChunks, rchunkMinVariance = %f\n", rchunkMinVariance);
	  fprintf( stderr, "in InsertWalkedChunks, delta: mean = %f, variance = %f\n", 
			   delta.mean, delta.variance);
	  
	  Force_Increasing_Variances_One_Scaffold( rchunk->scaffoldID );
	  AddDeltaToScaffoldOffsets(ScaffoldGraph,
								rchunk->scaffoldID,
								rchunk->id,
								1,
								0,
								delta);
	  Force_Increasing_Variances_One_Scaffold( rchunk->scaffoldID );
	}
  }
}

void MergeWalkedScaffolds( CIScaffoldT* scaff, CIScaffoldT* nextScaff, 
						   double nextScaffExtremeContigAEnd, double nextScaffExtremeContigBEnd,
						   ContigT *scaffExtremeContig, ContigT *nextScaffExtremeContig)
{
  FragOrient orient;
  LengthT offset;



  // these orient values need to be checked!
  // as does walking off the A end of a scaffold in general!
  // basically, redo this routine




  if (scaff->info.Contig.AEndCI == scaffExtremeContig->id)
  {
	if (nextScaff->info.Contig.AEndCI == nextScaffExtremeContig->id)
	  orient = B_A;
	else
	  orient = A_B;
  }
  else
  {
	if (nextScaff->info.Contig.AEndCI == nextScaffExtremeContig->id)
	  orient = A_B;
	else
	  orient = B_A;
  }
  
  offset.mean = min (nextScaffExtremeContigAEnd, nextScaffExtremeContigBEnd);
  offset.variance = max( scaffExtremeContig->offsetAEnd.variance, scaffExtremeContig->offsetBEnd.variance) +
	ComputeFudgeVariance( min( nextScaffExtremeContigAEnd, nextScaffExtremeContigBEnd) -
						  max( scaffExtremeContig->offsetAEnd.mean, scaffExtremeContig->offsetBEnd.mean));

  // the orient below is either A_B or B_A
  // the offset is where the min end of the scaffold goes
  // the orient says whether the B_END of the A_END is the min end
  InsertScaffoldContentsIntoScaffold(ScaffoldGraph,
                                     scaff->id, nextScaff->id,
                                     orient, &offset, TRUE); // contig now
  
  DeleteGraphNode( ScaffoldGraph->ScaffoldGraph, nextScaff);
  fprintf( stderr, "Killed Scaffold %d and merged it with Scaffold %d\n",
		   nextScaff->id, scaff->id);
}

void printScaffoldContigs(CIScaffoldT *scaffold)
{
    CIScaffoldTIterator CIsTemp;
	ChunkInstanceT
	  * lchunkTemp = NULL,
	  * rchunkTemp = NULL;
    //
    // print all the gaps of this scaffold
    //
    InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE,
							FALSE, &CIsTemp);
    while (NextCIScaffoldTIterator(&CIsTemp))
	{
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

	  fprintf( stderr, "gap in scaffold: lchunk: %d (%f, %f), rchunk: %d (%f, %f)\n", 
			   lchunkTemp->id, lchunkTemp->offsetAEnd.mean, lchunkTemp->offsetBEnd.mean, 
			   rchunkTemp->id, rchunkTemp->offsetAEnd.mean, rchunkTemp->offsetBEnd.mean);
	}	
	fprintf( stderr, "gap in scaffold: lchunk: %d (%f, %f), rchunk: %d (%f, %f)\n", 
			 lchunkTemp->id, lchunkTemp->offsetAEnd.mean, lchunkTemp->offsetBEnd.mean, 
			 rchunkTemp->id, rchunkTemp->offsetAEnd.mean, rchunkTemp->offsetBEnd.mean);
}

Overlap* OverlapContainingContigs(CIFragT *frag1, CIFragT *frag2, ChunkOrientationType *overlapOrientation,
								  ChunkInstanceT *leftContig, ChunkInstanceT *rightContig)
{
  // ChunkOverlapCheckT tempOlap1;
  Overlap * tempOlap1;
  ChunkInstanceT *contig1, *contig2;
  int orientationFrag1, orientationFrag2;
  int computedAhang;
  char *seq1, *seq2;
  int min_ahang, max_ahang;
  double erate, thresh;
  int minlen;

  // grab the corresponding contigs
  contig1 = GetGraphNode(ScaffoldGraph->RezGraph, frag1->contigID);
  contig2 = GetGraphNode(ScaffoldGraph->RezGraph, frag2->contigID);

  // walking frags ascending by iid 
  if (frag1->iid < frag2->iid)
  {
	if (frag1->contigOffset5p.mean < frag1->contigOffset3p.mean)
	  orientationFrag1 = 0;
	else
	  orientationFrag1 = 1;
		
	if (frag2->contigOffset5p.mean < frag2->contigOffset3p.mean)
	  orientationFrag2 = 0;
	else
	  orientationFrag2 = 1;
  }
  // walking frags descending by iid 
  else
  {
	if (frag1->contigOffset5p.mean < frag1->contigOffset3p.mean)
	  orientationFrag1 = 1;
	else
	  orientationFrag1 = 0;
		
	if (frag2->contigOffset5p.mean < frag2->contigOffset3p.mean)
	  orientationFrag2 = 1;
	else
	  orientationFrag2 = 0;
  }
  
  if (orientationFrag1 == 0 && orientationFrag2 == 0)
	*overlapOrientation = AB_AB;
  else if (orientationFrag1 == 0 && orientationFrag2 == 1)
	*overlapOrientation = AB_BA;
  else if (orientationFrag1 == 1 && orientationFrag2 == 0)
	*overlapOrientation = BA_AB;
  else if (orientationFrag1 == 1 && orientationFrag2 == 1)
	*overlapOrientation = BA_BA;

  // check to make sure left contig has the same orientation we have just assigned it
  if (contig1->id == leftContig->id) 
  {
	if (contig1->offsetAEnd.mean < contig1->offsetBEnd.mean)  // left contig is AB
	{
	  if ( *overlapOrientation == BA_AB || *overlapOrientation == BA_BA) // bad
		return NULL;
	}
	else  // left contig is BA
	{
	  if ( *overlapOrientation == AB_AB || *overlapOrientation == AB_BA) // bad
		return NULL;
	}	  
  }

  // check to make sure right contig has the same orientation we have just assigned it
  if (contig2->id == rightContig->id) 
  {
	if (contig2->offsetAEnd.mean < contig2->offsetBEnd.mean)  // right contig is AB
	{
	  if ( *overlapOrientation == AB_BA || *overlapOrientation == BA_BA) // bad
		return NULL;
	}
	else  // right contig is BA
	{
	  if ( *overlapOrientation == AB_AB || *overlapOrientation == BA_AB) // bad
		return NULL;
	}
  }

#if DEBUG_GAP_WALKER > -1		  
  fprintf( stderr, "\ncomputing overlap for contig1: %d (frag %d) [%d, %d] and contig2: %d (frag: %d)  [%d, %d]\n", 
		   contig1->id, frag1->iid, frag1->aEndCoord, frag1->bEndCoord,
		   contig2->id, frag2->iid, frag2->aEndCoord, frag2->bEndCoord);
  fprintf( stderr, "orientation is %c\n", (char) *overlapOrientation);
#endif
  
  computedAhang = ComputeAhangSize( frag1, contig1, frag2, contig2, *overlapOrientation);
  fprintf( stderr, "computedAhang: %d\n", computedAhang);

  min_ahang = computedAhang - 30;
  max_ahang = computedAhang + 30;

  if(consensus1 == NULL)
  {
	consensus1 = CreateVA_char(2048);
	consensus2 = CreateVA_char(2048);
	quality1 = CreateVA_char(2048);
	quality2 = CreateVA_char(2048);
  }

  // Get the consensus sequences for both chunks from the Store
  GetConsensus(ScaffoldGraph->RezGraph, contig1->id, consensus1, quality1);
  GetConsensus(ScaffoldGraph->RezGraph, contig2->id, consensus2, quality2);

  seq1 = Getchar(consensus1, 0);
  seq2 = Getchar(consensus2, 0);

  erate = CGW_DP_ERATE;
  thresh = CGW_DP_THRESH;
  minlen = CGW_DP_MINLEN;

  // tempOlap1 is a static down inside of DP_Compare
  tempOlap1 = OverlapSequences( seq1, seq2, *overlapOrientation, min_ahang, max_ahang, erate, thresh, minlen, AS_FIND_ALIGN);

  // the negative insures that we use it before any of the other edges lying around
  // fprintf(stderr, "set a negative edge weight between %d and %d\n", contig1->id, contig2->id);
  // Set_Bac_Walking_Quality( contig1->id, contig2->id, *overlapOrientation, -0.02);
  
  if (tempOlap1 != NULL)
  {
	// fprintf(stderr, "%d, %d ahang: %d, bhang:%d\n", contig1->id, contig2->id, tempOlap1->begpos, tempOlap1->endpos);

	if (fabs (tempOlap1->begpos - computedAhang) > 100.0)
	{
	  // dumpContigInfo(contig1);
	  // dumpContigInfo(contig2);
	  fprintf( stderr, "**** Warning: computed ahangs differ significantly in OverlapContainingContigs\n");
	}
	return tempOlap1;
  }
  else
  {
	fprintf(stderr, "%d, %d do not overlap\n", contig1->id, contig2->id);
	// dumpContigInfo(contig1);
	// dumpContigInfo(contig2);

	return NULL;	
  }
}

int ComputeOverlapSize( CIFragT *currFrag, ChunkInstanceT *contigCurrFrag, 
					CIFragT *fragSucc, ChunkInstanceT *contigSucc, int overlapOrientation)
{
  float currAEndDist, currBEndDist;
  float succAEndDist, succBEndDist;
  float computedOverlap = 0.f;
  int shredOverlap;
  int currLength, succLength;
  int fragDirection;
  
  currLength = currFrag->localePos.end - currFrag->localePos.bgn;
  succLength = fragSucc->localePos.end - fragSucc->localePos.bgn;
  
  // find out how much the two frags overlap by in shredding
  if (currFrag->iid < fragSucc->iid)
  {
	fragDirection = POSDIR;
	shredOverlap = currFrag->localePos.end - fragSucc->localePos.bgn;
  }
  else
  {
	fragDirection = NEGDIR;
	shredOverlap = fragSucc->localePos.end - currFrag->localePos.bgn;
  }
  
  // now figure out how much the contigs overlap by
  // draw the contigs in their respective orientation
  // draw the frags in the contigs (note that the frags are always in the same orientation to one another)
  // here's the picture for the normal (AB_AB) situation, fragment numbers ascending (POSDIR)
  /* 
     |   currAEndDist   |5  currFrag  3|   currBEndDist  |
    A------------------------------------------------>B  contigCurrFrag
                        ---------------> currFrag (i)

           |   succAEndDist   |5  fragSucc  3|   succBEndDist  |
          A------------------------------------------------>B  contigSucc
                              ---------------> fragSucc (i+1)
  */

  // here's the picture for the normal (AB_AB) situation, fragment numbers descending (NEGDIR)
  /* 
     |   currAEndDist   |3  currFrag  5|   currBEndDist  |
    A------------------------------------------------>B  contigCurrFrag
                        <--------------- currFrag (i+1)

           |   succAEndDist   |3  fragSucc  5|   succBEndDist  |
          A------------------------------------------------>B  contigSucc
                              <--------------- fragSucc (i)
  */

  if (overlapOrientation == AB_AB)
  {
	if (fragDirection == POSDIR)
	{
	  succAEndDist = fragSucc->contigOffset5p.mean;
	  currBEndDist = contigCurrFrag->bpLength.mean - currFrag->contigOffset3p.mean;
	}
	else
	{
	  succAEndDist = fragSucc->contigOffset3p.mean;
	  currBEndDist = contigCurrFrag->bpLength.mean - currFrag->contigOffset5p.mean;
	}
	computedOverlap = succAEndDist + shredOverlap + currBEndDist;
  }
  else if (overlapOrientation == AB_BA)
  {
	if (fragDirection == POSDIR)
	{
	  succBEndDist = contigSucc->bpLength.mean - fragSucc->contigOffset5p.mean;
	  currBEndDist = contigCurrFrag->bpLength.mean - currFrag->contigOffset3p.mean;
	}
	else
	{
	  succBEndDist = contigSucc->bpLength.mean - fragSucc->contigOffset3p.mean;
	  currBEndDist = contigCurrFrag->bpLength.mean - currFrag->contigOffset5p.mean;
	}
	computedOverlap = succBEndDist + shredOverlap + currBEndDist;
  }
  else if (overlapOrientation == BA_AB)
  {
	currAEndDist = currFrag->contigOffset3p.mean;
	currBEndDist = contigCurrFrag->bpLength.mean - currFrag->contigOffset5p.mean;
	succAEndDist = fragSucc->contigOffset5p.mean;
	succBEndDist = contigSucc->bpLength.mean - fragSucc->contigOffset3p.mean;

	computedOverlap = min( currBEndDist + (currLength - shredOverlap), succAEndDist) + shredOverlap +
	  min( currAEndDist, succBEndDist + (succLength - shredOverlap));
  }
  else if (overlapOrientation == BA_BA)
  {
	currAEndDist = currFrag->contigOffset3p.mean;
	currBEndDist = contigCurrFrag->bpLength.mean - currFrag->contigOffset5p.mean;
	succAEndDist = fragSucc->contigOffset3p.mean;
	succBEndDist = contigSucc->bpLength.mean - fragSucc->contigOffset5p.mean;

	computedOverlap = min( currBEndDist + (currLength - shredOverlap), succBEndDist) + shredOverlap + 
	  min( currAEndDist, succAEndDist + (succLength - shredOverlap));
  }
  return ((int) computedOverlap);
}
  
int ComputeAhangSize( CIFragT *currFrag, ChunkInstanceT *contigCurrFrag, 
					  CIFragT *fragSucc, ChunkInstanceT *contigSucc, int overlapOrientation)
{
  float currAEndDist, currBEndDist;
  float succAEndDist, succBEndDist;
  float computedAhang = 0.f;
  int shredOverlap;
  int currLength, succLength;
  int fragDirection;
  
  currLength = currFrag->localePos.end - currFrag->localePos.bgn;
  succLength = fragSucc->localePos.end - fragSucc->localePos.bgn;
  
  // find out how much the two frags overlap by in shredding
  if (currFrag->iid < fragSucc->iid)
  {
	fragDirection = POSDIR;
	shredOverlap = currFrag->localePos.end - fragSucc->localePos.bgn;
  }
  else
  {
	fragDirection = NEGDIR;
	shredOverlap = fragSucc->localePos.end - currFrag->localePos.bgn;
  }

  if (shredOverlap <= 200)
    fprintf( stderr, "warning: adjacent shredded fragments %d and %d have an overlap of %d\n",
	     currFrag->iid, fragSucc->iid, shredOverlap);

  // now figure out how much the contigs overlap by
  // draw the contigs in their respective orientation
  // draw the frags in the contigs (note that the frags are always in the same orientation to one another)
  // here's the picture for the normal (AB_AB) situation, fragments POSDIR
  /* 
     |   currAEndDist   |5  currFrag  3|   currBEndDist  |
    A------------------------------------------------>B  contigCurrFrag
                        ---------------> currFrag

           |   succAEndDist   |5  fragSucc  3|   succBEndDist  |
          A------------------------------------------------>B  contigSucc
                              ---------------> fragSucc
  */
  // here's the picture for the normal (AB_AB) situation, fragments NEGDIR
  /* 
     |   currAEndDist   |3  currFrag  5|   currBEndDist  |
    A------------------------------------------------>B  contigCurrFrag
                        <--------------- currFrag

           |   succAEndDist   |3  fragSucc  5|   succBEndDist  |
          A------------------------------------------------>B  contigSucc
                              <--------------- fragSucc
  */

  if (overlapOrientation == AB_AB)
  {
	if (fragDirection == POSDIR)
	{
	  currAEndDist = currFrag->contigOffset5p.mean;
	  succAEndDist = fragSucc->contigOffset5p.mean;
	}
	else
	{
	  currAEndDist = currFrag->contigOffset3p.mean;
	  succAEndDist = fragSucc->contigOffset3p.mean;
	}
	computedAhang = currAEndDist + (currLength - shredOverlap) - succAEndDist;
  }
  else if (overlapOrientation == AB_BA)
  {
	if (fragDirection == POSDIR)
	{
	  currAEndDist = currFrag->contigOffset5p.mean;
	  succBEndDist = contigSucc->bpLength.mean - fragSucc->contigOffset5p.mean;
	}
	else
	{
	  currAEndDist = currFrag->contigOffset3p.mean;
	  succBEndDist = contigSucc->bpLength.mean - fragSucc->contigOffset3p.mean;
	}
	computedAhang = currAEndDist + (currLength - shredOverlap) - succBEndDist;
  }
  else if (overlapOrientation == BA_AB)
  {
	if (fragDirection == POSDIR)
	{
	  currBEndDist = contigCurrFrag->bpLength.mean - currFrag->contigOffset5p.mean;
	  succAEndDist = fragSucc->contigOffset5p.mean;
	}
	else
	{
	  currBEndDist = contigCurrFrag->bpLength.mean - currFrag->contigOffset3p.mean;
	  succAEndDist = fragSucc->contigOffset3p.mean;
	}	
	computedAhang = currBEndDist + (currLength - shredOverlap) - succAEndDist;
  }
  else if (overlapOrientation == BA_BA)
  {
	if (fragDirection == POSDIR)
	{
	  currBEndDist = contigCurrFrag->bpLength.mean - currFrag->contigOffset5p.mean;
	  succBEndDist = contigSucc->bpLength.mean - fragSucc->contigOffset5p.mean;
	}
	else
	{
	  currBEndDist = contigCurrFrag->bpLength.mean - currFrag->contigOffset3p.mean;
	  succBEndDist = contigSucc->bpLength.mean - fragSucc->contigOffset3p.mean;
	}
	computedAhang = currBEndDist + (currLength - shredOverlap) - succBEndDist;
  }
  return ((int) computedAhang);
}
  
void UnlinkLocaleOverlaps(int locale)
{
  int inode;
  chunk_subgraph *localeGraph = NULL;

  // localeGraph should be made out of only Contigs that contain the locale we're working with
  // the current hack is to store that locale number in the global currentLocale
  // currentLocale = locale;
  
  // now build the graph that contains only contigs that contain frags from the locale of interest
  // this should be replaced with calls to the "Bactig Store" when we get it
  // the store will tell us which frags are in which contigs
  localeGraph = Build_Subgraph(NULL, -1, -1, 0, 
							   Contains_Locale_Frags,
							   All_Edges);
  
  edgesUnlinkedListHead = NULL;
  edgesUnlinkedList = NULL;  

  for (inode = 0; inode < localeGraph->size; inode++)
  {
	CIEdgeT * e;
	GraphEdgeIterator edges;
	int32 cid, next_chunk_cid;
	// ChunkInstanceT *contig;

	cid = localeGraph->node[inode].cid;
	// contig = GetGraphNode(ScaffoldGraph->RezGraph, cid);
	InitGraphEdgeIterator(ScaffoldGraph->RezGraph, cid,
                              ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
	while((e = NextGraphEdgeIterator(&edges)) != NULL)
	{
	  assert(e != NULL);

	  //
	  // get the other end
	  //
	  if (cid == e->idA)
		next_chunk_cid = e->idB;
	  else
		next_chunk_cid = e->idA;

	  //
	  // check if the next chunk is in the subgraph
	  //
	  if (Belong_To(localeGraph, next_chunk_cid)) 
	  {
		unlinkedEdgeT *uedge;
		assert(localeGraph->table[next_chunk_cid] != NULL);
		
		// add this edge to the list we are gong to unlink
		uedge = (unlinkedEdgeT *) malloc( sizeof( unlinkedEdgeT ));
		uedge->next = NULL;
		uedge->edge = e;
		
		// if we have none in the list for this graph initialize
		if (edgesUnlinkedList == NULL)
		{
		  edgesUnlinkedList = uedge;
		  edgesUnlinkedListHead = uedge;
		}
		else  // add edge to list
		{
		  edgesUnlinkedList->next = uedge;
		}
	  }
	}
  }
	
  // now we know what edges to unlink for this contig
  edgesUnlinkedList = edgesUnlinkedListHead;
  if (edgesUnlinkedList != NULL)
  {
	while (edgesUnlinkedList->next != NULL)
	{
	  edgesUnlinkedList->eid = GetVAIndex_EdgeCGW_T(ScaffoldGraph->RezGraph->edges, edgesUnlinkedList->edge);
	  fprintf( stderr, "unlinking edge between: %d and %d\n", edgesUnlinkedList->edge->idA, edgesUnlinkedList->edge->idB);
	  UnlinkGraphEdge(ScaffoldGraph->RezGraph, edgesUnlinkedList->edge);
	  edgesUnlinkedList = edgesUnlinkedList->next;
	}
  }
}

void LinkLocaleOverlaps(void)
{
  unlinkedEdgeT *uedge = NULL;

  // if we did no unlinking return
  if (edgesUnlinkedListHead == NULL)
	return;
  
  // we know what edges we unlinked for all the contigs
  edgesUnlinkedList = edgesUnlinkedListHead;
  while (edgesUnlinkedList->next != NULL)
  {
	InsertGraphEdge( ScaffoldGraph->RezGraph, edgesUnlinkedList->eid, 0);
	uedge = edgesUnlinkedList;
	edgesUnlinkedList = edgesUnlinkedList->next;
	free(uedge);
  }
  free(edgesUnlinkedList);
  edgesUnlinkedList = NULL;
  edgesUnlinkedListHead = NULL;  
}

void DumpLocaleOverlaps(void)
{
  
  // we know what edges we unlinked for all the contigs
  edgesUnlinkedList = edgesUnlinkedListHead;
  while (edgesUnlinkedList->next != NULL)
  {
	fprintf( stderr, "edgesUnlinkedList, edge: %p, eid: %d, next.eid: %d\n", 
			 edgesUnlinkedList->edge, edgesUnlinkedList->eid, edgesUnlinkedList->next->eid);
	edgesUnlinkedList = edgesUnlinkedList->next;
  }
  fprintf( stderr, "edgesUnlinkedList, edge: %p, eid: %d\n", 
		   edgesUnlinkedList->edge, edgesUnlinkedList->eid);
}

int GetFragInfo (ChunkInstanceT* chunk, int locale, unsigned int* minFragIid, unsigned int *maxFragIid)
{
  int numFrags, ii;
  MultiAlignT *uma;
  IntMultiPos *ump;
  
  *minFragIid = CDS_INT32_MAX;
  *maxFragIid = 0;
  
  // then get multi-align for the node
  uma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, chunk->id, ScaffoldGraph->RezGraph->type == CI_GRAPH); 
  // uma = GetMultiAlignInStore(ScaffoldGraph->RezGraph->maStore, chunk->id);
  numFrags = GetNumIntMultiPoss(uma->f_list);
  
  for (ii = 0; ii < numFrags; ii++)
  {
	CIFragT *frag;  //
	ump = GetIntMultiPos(uma->f_list, ii);
	
	frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32) ump->source);
	
	// deal with frags from locales first
	//GetFragmentPositionInContigFromChunk( frag, &min, &max, &fragContigOrientation,
	//								  chunkContigOrientation, 
	//								  chunkLeftEnd, chunkRightEnd);
	
	if (frag->locale == locale)
	{
	  if (frag->iid < *minFragIid)
		*minFragIid = frag->iid;
	  if (frag->iid > *maxFragIid)
		*maxFragIid = frag->iid;
	}
	
	if (0)
	{
	  fprintf(GlobalData->gwlogfp,
			  "frag->iid %d: containing contig: %d, celsim coords: %d, %d, locale: %d\n",
			  frag->iid,
			  frag->contigID,
			  frag->aEndCoord,
			  frag->bEndCoord,
			  frag->locale);
	}
  }
  return *minFragIid;
}

static VA_TYPE(char) *consensus = NULL;
static VA_TYPE(char) *quality = NULL;

void dumpContigInfo(ChunkInstanceT *contig)
{
  int 
	offset, contigLeftEnd, contigRightEnd, contigOrientation;
  MultiAlignT 
	*ma;
  char *seq1;

  // return;
  
  fprintf( stderr, "*********************** contig analysis **************************\n");
  fprintf( stderr, "analyzing contig: %d\n", contig->id);
  fflush(stderr);

  if (contig->offsetAEnd.mean < contig->offsetBEnd.mean)
	contigOrientation = 0;
  else
	contigOrientation = 1;

  fprintf( stderr, "contig orientation: %d\t length: %d  ", contigOrientation, (int) contig->bpLength.mean);
  fprintf( stderr, "contig offsetAEnd: %d\t offsetBEnd: %d\n", (int) contig->offsetAEnd.mean, (int) contig->offsetBEnd.mean);

  offset = contigLeftEnd = min( (int) contig->offsetAEnd.mean, (int) contig->offsetBEnd.mean);
  contigRightEnd = max( (int) contig->offsetAEnd.mean, (int) contig->offsetBEnd.mean);

  // now some contig info
  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, ScaffoldGraph->RezGraph->type == CI_GRAPH); 
  // ma = GetMultiAlignInStore(ScaffoldGraph->RezGraph->maStore, contig->id);

  if(consensus == NULL)
  {
    consensus = CreateVA_char(2048);
    quality = CreateVA_char(2048);
  }
		
  // Get the consensus sequences for the contig from the Store
  GetConsensus(ScaffoldGraph->ContigGraph, contig->id, consensus, quality);
  seq1 = Getchar(consensus, 0);

  if (contigOrientation == 1)
	SequenceComplement(seq1, NULL);
  
  fprintf( stderr, "> contig %d consensus seq (flipped to reflect scaff orientation)\n", contig->id);
  fprintf( stderr, "%s\n", seq1);

  {
	char tempChar;
	
	tempChar = seq1[2500];
	seq1[2500] = '\0';
	
	fprintf( stderr, "> contig %d left end\n", contig->id);
	fprintf( stderr, "%s\n", seq1);

	seq1[2500] = tempChar;
  }

  contigRightEnd = strlen(seq1) - 2501;
  fprintf( stderr, "> contig %d right end\n", contig->id);
  fprintf( stderr, "%s\n", &seq1[contigRightEnd]);

#if 0
  numUnitigs = GetNumIntUnitigPoss(ma->u_list);
  fprintf( stderr, "number unitigs: %d\n", numUnitigs);

  for (i = 0; i < numUnitigs; i++)
  {
	IntUnitigPos *upos = GetIntUnitigPos( ma->u_list, i);
	ChunkInstanceT *unitig = GetGraphNode( ScaffoldGraph->CIGraph, upos->ident);
	MultiAlignT *uma;
	IntMultiPos *ump;
	int icntfrag;
	ChunkInstanceT *surrogateUnitig;

	uma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, ScaffoldGraph->CIGraph->type == CI_GRAPH); 
	// uma = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, unitig->id);
	fprintf( stderr, "  unitig: %d\t num frags: %ld surrogate: %d\n", unitig->id, GetNumIntMultiPoss(uma->f_list),
		 (unitig->flags.bits.isStoneSurrogate || unitig->flags.bits.isWalkSurrogate));

	surrogateUnitig = unitig;
	if (unitig->flags.bits.isStoneSurrogate || unitig->flags.bits.isWalkSurrogate)
	{
	  fprintf ( stderr, "  surrogate unitig offsetAEnd: %f, offsetBEnd: %f\n", 
		    unitig->offsetAEnd.mean, unitig->offsetBEnd.mean);
	  unitig = GetGraphNode( ScaffoldGraph->CIGraph, unitig->info.CI.baseID);
	  fprintf ( stderr, "  using original unitig: %d\n", unitig->id);
	  uma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, 
					      ScaffoldGraph->CIGraph->type == CI_GRAPH); 
	}

	// now print out info on the frags in the unitig
	for (icntfrag = 0; icntfrag < GetNumIntMultiPoss(uma->f_list); icntfrag++)
	{
	  CIFragT *frag;  //
	  IntMultiPos *imp = GetIntMultiPos(uma->f_list, icntfrag);
	  InfoByIID *info;

	  // get frags position in contig	
	  info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, imp->ident);
	  assert(info->set);
	  frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
	  
	  fprintf( stderr, "    frag: %6d\t contig pos (5p, 3p): %6d, %6d", 
			   imp->ident, (int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean);	  
	  fprintf( stderr, "   real pos (a, b): %6d, %6d, type: %c\n",
			   frag->aEndCoord, frag->bEndCoord, frag->type);
	}
  }
#endif
  // now get some edge info
  if (0)
  {
	CIEdgeT
	  * e;
	GraphEdgeIterator
	  edges;

	InitGraphEdgeIterator(ScaffoldGraph->RezGraph, contig->id,
                              ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
	while((e = NextGraphEdgeIterator(&edges)) != NULL)
	{
	  assert(e != NULL);
	  
	  // fprintf( stderr, "edge: %d to %d, orient: %c, length: %f\n", e->idA, e->idB, (char) e->orient, e->distance.mean);
	  PrintGraphEdge( stderr, ScaffoldGraph->RezGraph, "Analyzing edge", e, 0);
	}
  }
}

void dumpGapInfo(ChunkInstanceT *leftContig, ChunkInstanceT *rightContig)
{

  // return;
  
  fprintf( stderr, "*********************** gap analysis **************************\n");
  fprintf( stderr, "analyzing gap (%d, %d)\n", leftContig->id, rightContig->id);
  fflush(stderr);

  // DumpContig( GlobalData->stderrc, ScaffoldGraph, leftContig, TRUE);
  // DumpContig( GlobalData->stderrc, ScaffoldGraph, rightContig, TRUE);

  if (1)
  {
	CIEdgeT* e;
	GraphEdgeIterator edges;
	
	InitGraphEdgeIterator(ScaffoldGraph->RezGraph, leftContig->id,
						  ALL_END, ALL_EDGES, ITERATOR_VERBOSE | GRAPH_EDGE_RAW_ONLY, &edges);
	while((e = NextGraphEdgeIterator(&edges)) != NULL)
	{
	  CIFragT *leftFrag, *rightFrag, *tempFrag;	  
	  int leftFragLeftEnd, leftFragRightEnd, leftFragScaffoldOrientation;
	  int rightFragLeftEnd, rightFragRightEnd, rightFragScaffoldOrientation;
	  int leftContigLeftEnd, leftContigRightEnd, leftContigOrientation;
	  int rightContigLeftEnd, rightContigRightEnd, rightContigOrientation;
	  
	  assert(e != NULL);

 	  //fprintf( stderr, "e->idA: %d (frag: %d), e->idB: %d (frag: %d), e->distance.mean: %f\n", 
	  //	   e->idA, e->fragA, e->idB, e->fragB, e->distance.mean);

	  if (e->fragA == NULLINDEX || e->fragB == NULLINDEX)
		continue;
	  
	  if (e->fragA < 0 || e->fragB < 0)
		continue;
	  
	  leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, e->fragA);
	  if (leftFrag->contigID != leftContig->id)
		continue;
	  
	  rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, e->fragB);
	  if (rightFrag->contigID != rightContig->id)
		continue;

	  // if idA, idB != leftContig, rightContig
	  if (e->idA != leftContig->id)
	  {
		tempFrag = leftFrag;
		leftFrag = rightFrag;
		rightFrag = tempFrag;
	  }

	  GetFragmentPositionInScaffold( leftFrag, &leftFragLeftEnd, &leftFragRightEnd, &leftFragScaffoldOrientation);
	  GetContigPositionInScaffold( leftContig, &leftContigLeftEnd, &leftContigRightEnd, &leftContigOrientation);

	  GetFragmentPositionInScaffold( rightFrag, &rightFragLeftEnd, &rightFragRightEnd, &rightFragScaffoldOrientation);
	  GetContigPositionInScaffold( rightContig, &rightContigLeftEnd, &rightContigRightEnd, &rightContigOrientation);	  

	  // we are looking for 10k mates
	  if ( rightFragLeftEnd - leftFragRightEnd > 35000 ||
		   rightFragLeftEnd - leftFragRightEnd < 4000)
		continue;

	  fprintf( stderr, "%d to %d, o: %c, len: %f ", e->idA, e->idB, (char) e->orient, e->distance.mean);
	  fprintf( stderr, "lfrag: %d (ctg: %d, %d from end), rfrag: %d (ctg: %d, %d from end)\n",
			   leftFrag->iid, leftFrag->contigID, leftContigRightEnd - leftFragRightEnd,
			   rightFrag->iid, rightFrag->contigID, rightFragLeftEnd - rightContigLeftEnd);

	  PrintGraphEdge( stderr, ScaffoldGraph->ContigGraph, "\t", e, leftContig->id);

	  // fprintf( stderr, "edge: %d to %d, orient: %c, length: %f\n", e->idA, e->idB, (char) e->orient, e->distance.mean);
	  // PrintGraphEdge( stderr, ScaffoldGraph->RezGraph, "Analyzing edge", e, 0);

	}
  }
}

void GetContigPositionInScaffold(ChunkInstanceT *contig, int *left_end, int *right_end, 
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
}

// this routine should work for surrogates as well as normal unitigs
void GetCIPositionInScaffold(ChunkInstanceT *CI, int *left_end, int *right_end, 
							 int *CIScaffoldOrientation)
{
  int contigLeftEnd, contigRightEnd, contigScaffoldOrientation;
  NodeCGW_T *contig = GetGraphNode( ScaffoldGraph->ContigGraph, CI->info.CI.contigID);

  GetContigPositionInScaffold( contig, &contigLeftEnd, &contigRightEnd, &contigScaffoldOrientation);  

  if ( contigLeftEnd < contigRightEnd )  // this is relative to scaffold
  {
	if (CI->offsetAEnd.mean <= CI->offsetBEnd.mean) // this is relative to contig
	{
	  *left_end = contigLeftEnd + CI->offsetAEnd.mean;
	  *right_end = contigLeftEnd + CI->offsetBEnd.mean;
	  *CIScaffoldOrientation = 0;
	} 
	else 
	{
	  // *left_end = contigRightEnd - CI->offsetBEnd.mean;
	  // *right_end = contigRightEnd - CI->offsetAEnd.mean;
	  *left_end = contigLeftEnd + CI->offsetBEnd.mean;
	  *right_end = contigLeftEnd + CI->offsetAEnd.mean;
	  *CIScaffoldOrientation = 1;
	}
  }
  else  // contig is reversed in scaffold
  {
	if (CI->offsetAEnd.mean <= CI->offsetBEnd.mean) // this is relative to contig
	{
	  *left_end = contigRightEnd - CI->offsetBEnd.mean;
	  *right_end = contigRightEnd - CI->offsetAEnd.mean;
	  *CIScaffoldOrientation = 1;
	} 
	else 
	{
	  *left_end = contigRightEnd - CI->offsetAEnd.mean;
	  *right_end = contigRightEnd - CI->offsetBEnd.mean;
	  *CIScaffoldOrientation = 0;
	}
  }
}


void GetFragmentPositionInScaffoldFromContig(CIFragT *frag, int *left_end, int *right_end, 
											 int *fragmentScaffoldOrientation, 
											 int contigLeftEnd, int contigRightEnd, int contigScaffoldOrientation)
{
  if (contigScaffoldOrientation == 0)  // contig is direct in scaffold
  {
	if (frag->contigOffset5p.mean < frag->contigOffset3p.mean)  // frag is direct in contig
	{
	  *left_end = contigLeftEnd + frag->contigOffset5p.mean;
	  *right_end = contigLeftEnd + frag->contigOffset3p.mean;
	  *fragmentScaffoldOrientation = 0;
	}
	else  // frag is reversed in contig
	{
	  *left_end = contigLeftEnd + frag->contigOffset3p.mean;
	  *right_end = contigLeftEnd + frag->contigOffset5p.mean;
	  *fragmentScaffoldOrientation = 1;
	}
  }
  else   // contig is reversed in scaffold
  {
	if (frag->contigOffset5p.mean < frag->contigOffset3p.mean)  // frag is direct in contig
	{
	  *left_end = contigRightEnd - frag->contigOffset3p.mean;
	  *right_end = contigRightEnd - frag->contigOffset5p.mean;
	  *fragmentScaffoldOrientation = 1;
	}
	else  // frag is reversed in contig
	{
	  *left_end = contigRightEnd - frag->contigOffset5p.mean;
	  *right_end = contigRightEnd - frag->contigOffset3p.mean;
	  *fragmentScaffoldOrientation = 0;
	}
  }
}

void GetFragmentPositionInScaffoldFromCI(CIFragT *frag, int *left_end, int *right_end, 
										 int *fragmentScaffoldOrientation, 
										 int CILeftEnd, int CIRightEnd, int CIScaffoldOrientation)
{
  if (CIScaffoldOrientation == 0)  // CI is direct in scaffold
  {
	if (frag->offset5p.mean < frag->offset3p.mean)  // frag is direct in CI
	{
	  *left_end = CILeftEnd + frag->offset5p.mean;
	  *right_end = CILeftEnd + frag->offset3p.mean;
	  *fragmentScaffoldOrientation = 0;
	}
	else  // frag is reversed in CI
	{
	  *left_end = CILeftEnd + frag->offset3p.mean;
	  *right_end = CILeftEnd + frag->offset5p.mean;
	  *fragmentScaffoldOrientation = 1;
	}
  }
  else   // CI is reversed in scaffold
  {
	if (frag->offset5p.mean < frag->offset3p.mean)  // frag is direct in CI
	{
	  *left_end = CIRightEnd - frag->offset3p.mean;
	  *right_end = CIRightEnd - frag->offset5p.mean;
	  *fragmentScaffoldOrientation = 1;
	}
	else  // frag is reversed in CI
	{
	  *left_end = CIRightEnd - frag->offset5p.mean;
	  *right_end = CIRightEnd - frag->offset3p.mean;
	  *fragmentScaffoldOrientation = 0;
	}
  }
}

void GetFragmentPositionInScaffold(CIFragT *frag, int *left_end, int *right_end, 
								   int *fragmentScaffoldOrientation)
{
  ContigT *containingContig = GetGraphNode(ScaffoldGraph->ContigGraph, frag->contigID);
  int contigLeftEnd, contigRightEnd, contigScaffoldOrientation;
  
  GetContigPositionInScaffold( containingContig, &contigLeftEnd, &contigRightEnd, &contigScaffoldOrientation);

  GetFragmentPositionInScaffoldFromContig( frag, left_end, right_end, fragmentScaffoldOrientation,
									   contigLeftEnd, contigRightEnd, contigScaffoldOrientation);
}

void GetFragmentPositionInScaffold_new(CIFragT *frag, int *left_end, int *right_end, 
									   int *fragmentScaffoldOrientation, int operateOnContigs)
{
  if (operateOnContigs == TRUE)
  {
	ContigT *containingContig = GetGraphNode(ScaffoldGraph->ContigGraph, frag->contigID);
	int contigLeftEnd, contigRightEnd, contigScaffoldOrientation;
	
	GetContigPositionInScaffold( containingContig, &contigLeftEnd, &contigRightEnd, &contigScaffoldOrientation);
	
	GetFragmentPositionInScaffoldFromContig( frag, left_end, right_end, fragmentScaffoldOrientation,
											 contigLeftEnd, contigRightEnd, contigScaffoldOrientation);
  }
  else
  {
	NodeCGW_T *containingCI = GetGraphNode(ScaffoldGraph->CIGraph, frag->CIid);
	int CILeftEnd, CIRightEnd, CIScaffoldOrientation;
	
	GetCIPositionInScaffold( containingCI, &CILeftEnd, &CIRightEnd, &CIScaffoldOrientation);
	
	GetFragmentPositionInScaffoldFromCI( frag, left_end, right_end, fragmentScaffoldOrientation,
										 CILeftEnd, CIRightEnd, CIScaffoldOrientation);
  }
}



/*
----------------------------------------------------------------------------------
----------------------------------------------------------------------------------
Cam stuff
----------------------------------------------------------------------------------
----------------------------------------------------------------------------------
*/

#define NUM_COLORS 10

// an attempt to walk through the unitigs of a contig and grab their fragments,
// thus allowing us to detect a surrogate and grab the parent's frags
// the original localeCam operates off of a contig's f_list
void localeCam(char *middleName) 
{
  int
    i,
    numFrags;
  int32
    min,
    max;
  char
    filename[STR_LEN],
	*ScaffoldColor = "0ScaffoldColor: CFF0040 T2 S # Scaffolds",
	*BactigGapColor= "0BactigGapColor: CFF0040 T2 S # BactigGap",
	*FragShiftColor= "0FragShiftColor: CFF8040 T2 S # FragShift",
	*UniqueColor   = "0UniqueColor: CFFFF00 T2 S  # Unique",
	*UniqueColorStr= "UniqueColor",
	*RockColor     = "0RockColor: CFF0000 T2 S  # RockCI",
	*RockColorStr  = "RockColor",
	*StoneColor    = "0StoneColor: C77EF77 T2 S  # StoneCI",
	*StoneColorStr = "StoneColor",
	*WalkColor     = "0WalkColor: C8080FF T2 S  # WalkCI",
	*WalkColorStr  = "WalkColor",
    * Colour[NUM_COLORS] = {
      "C0040FF T2 S # Contig",
      "C0AAAA0 T2 S # non-BAC frag",
      "C0000F0 T2 S # locale_1",
      "CAAAA00 T2 S # locale_2",
      "C00FFF0 T2 S # locale_3",
      "CFF8000 T2 S # locale_4",
      "CA0A0FF T2 S # locale_5",
      "CA00A0F T2 S # locale_6",
      "CFFA0F0 T2 S # locale_7",
      "C00AAAA T1 S # locale_8"};
  FILE
    * cam_file;
  MultiAlignT 
	*ma;
  CIFragT 
	*frag;
  NodeCGW_T 
	*contig = NULL;
  int sc, scaffoldOffset, currentScaffoldID;
  int contigCnt, contigIDs[1000000];
#define MAX_FRAG_IID 1000000
  int fragsPresent[MAX_FRAG_IID], fragsMinPos[MAX_FRAG_IID], fragsMaxPos[MAX_FRAG_IID], 
	fragsScaffold[MAX_FRAG_IID], fragsLocale[MAX_FRAG_IID], fragsLocaleMin[MAX_FRAG_IID], fragsLocaleMax[MAX_FRAG_IID];
  int icnt;
  int lastFragIid, lastFragIidPerm, lastFragMinPos, lastFragMaxPos, lastFragLocale, lastFragScaffoldId;

  for (icnt = 0; icnt < MAX_FRAG_IID; icnt++)
	fragsPresent[icnt] = 0;
  
  //
  // open the cam file
  //
  sprintf(filename, "./cam/locales_%s.cam", middleName);
  system("mkdir cam");
  // sprintf(filename, "./cam/locale.%d.cam", locale);
  cam_file = file_open (filename, "w");
  assert(cam_file != NULL);
  
  //
  // output the colors
  //
  fprintf(cam_file, "%s\n", ScaffoldColor);
  fprintf(cam_file, "%s\n", BactigGapColor);
  fprintf(cam_file, "%s\n", FragShiftColor);
  fprintf(cam_file, "%s\n", UniqueColor);
  fprintf(cam_file, "%s\n", RockColor);
  fprintf(cam_file, "%s\n", StoneColor);
  fprintf(cam_file, "%s\n", WalkColor);
  for (i = 0; i < NUM_COLORS; i++)
    fprintf(cam_file, "%d: %s\n", i, Colour[i]);

  scaffoldOffset = 0;
  currentScaffoldID = -2;  

  for (sc = 0; sc < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sc++) 
  {
	CIScaffoldTIterator
	  CIs;
	CIScaffoldT
	  * scaff;
	int contigOffset, numUnitigs;
	int contigOrientation, contigLeftEnd, contigRightEnd;  
	char* unitigColor;
	int readCnt = 0;
	
	scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sc);
	
	// make sure the scaffold is there
    assert(scaff != NULL);

    if ((isDeadCIScaffoldT(scaff)) ||
       	(scaff->type != REAL_SCAFFOLD)) // || (scaff->info.Scaffold.numElements < 2))
    {
	  continue;
    }

	contigCnt = 0;
    InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE,	FALSE, &CIs);
    while (NextCIScaffoldTIterator(&CIs))
	{
	  contig = GetGraphNode(ScaffoldGraph->ContigGraph, CIs.curr);
	  contigIDs[contigCnt++] = contig->id;
	  
	  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
	  // ma = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id);

	  numUnitigs = GetNumIntUnitigPoss(ma->u_list);
	  
	  if (contig->offsetAEnd.mean < contig->offsetBEnd.mean)
		contigOrientation = 0;
	  else
		contigOrientation = 1;
	  
	  contigOffset = contigLeftEnd = min( (int) contig->offsetAEnd.mean, (int) contig->offsetBEnd.mean);
	  contigRightEnd = max( (int) contig->offsetAEnd.mean, (int) contig->offsetBEnd.mean);
	  
	  if (contig->scaffoldID >= 0)
		fprintf(cam_file,
				"%d: %d A%d %d R1 # contig: %d, containing scaffold: %d\n",
				contig->id,
				scaffoldOffset + contigLeftEnd,
				0,
				scaffoldOffset + contigRightEnd,
				contig->id,
				contig->scaffoldID);
	  
	  //if( !contains_fbac(contig))
	  //continue;
	  
	  // get the unitigs in this contig
	  numUnitigs = GetNumIntUnitigPoss(ma->u_list);
	  for (i = 0; i < numUnitigs; i++)
	  {
		IntUnitigPos *upos = GetIntUnitigPos( ma->u_list, i);
		ChunkInstanceT *originalUnitig, *unitig = GetGraphNode( ScaffoldGraph->CIGraph, upos->ident);
		MultiAlignT *uma;
		IntMultiPos *ump;
		int chunkContigOrientation, chunkLeftEnd, chunkRightEnd;
		int isSurrogate, ii;
		
		// find the orientation of this chunk in the contig
		GetChunkPositionInContig(unitig, &chunkLeftEnd, &chunkRightEnd, &chunkContigOrientation);

		if (0)
		  fprintf(stderr, 
				  "contig: %d, unitig: %d, unitig->info.CI.baseID: %d, contigOrientation: %d, chunkContigOrientation: %d\n",
				  contig->id, unitig->id, unitig->info.CI.baseID, contigOrientation, chunkContigOrientation);
		
		unitigColor = UniqueColorStr;
		if (unitig->flags.bits.isRock)
		  unitigColor = RockColorStr;
		else if (unitig->flags.bits.isStone || unitig->flags.bits.isStoneSurrogate)
		  unitigColor = StoneColorStr;
		else if (unitig->flags.bits.isWalk || unitig->flags.bits.isWalkSurrogate)
		  unitigColor = WalkColorStr;
		
		if (contigOrientation == 0)
		  fprintf(cam_file,
				  "%dUnitig: %d A0%s %d R2 # containing contig: %d, unitig->id: %d, baseID: %d\n",
				  unitig->id,   // the 1000000 prevents collisons in the Celamy namespace
				  scaffoldOffset + contigLeftEnd + chunkLeftEnd,
				  unitigColor,
				  scaffoldOffset + contigLeftEnd + chunkRightEnd,
				  unitig->info.CI.contigID,
				  unitig->id,
				  unitig->info.CI.baseID);
		else
		  fprintf(cam_file,
				  "%dUnitig: %d A0%s %d R2 # containing contig: %d, unitig->id: %d, baseID: %d contig reversed\n",
				  unitig->id,   // the 1000000 prevents collisons in the Celamy namespace
				  scaffoldOffset + contigLeftEnd + (int) contig->bpLength.mean - chunkRightEnd,
				  unitigColor,
				  scaffoldOffset + contigLeftEnd + (int) contig->bpLength.mean - chunkLeftEnd,
				  unitig->info.CI.contigID,
				  unitig->id,
				  unitig->info.CI.baseID);

		// if this is a surrogate switch pointer to parent, put frags in surrogate temporarily
		isSurrogate = FALSE;
		originalUnitig = unitig;
		if (unitig->info.CI.baseID > 0)
		{
		  isSurrogate = TRUE;
		  unitig = GetGraphNode( ScaffoldGraph->CIGraph, unitig->info.CI.baseID);
		}

		// then get multi-align for that unitig
		uma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, ScaffoldGraph->CIGraph->type == CI_GRAPH); 
		// uma = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, unitig->id);
		numFrags = GetNumIntMultiPoss(uma->f_list);
		
		if (0)
		{
		  fprintf(stderr, "contig: %d, parent unitig: %d, unitig->info.CI.baseID: %d, numFrags: %d\n",
				contig->id, unitig->id, unitig->info.CI.baseID, numFrags);
		}

		if (1)
		{
		  for (ii = 0; ii < numFrags; ii++)
		  {
			ump = GetIntMultiPos(uma->f_list, ii);
			frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32)ump->source);
		
			// deal with frags from locales first
			if (frag->locale != -1) 
			{
			  int fragContigOrientation;
			  int fragLocaleMod;
			  
			  GetFragmentPositionInContigFromChunk( frag, &min, &max, &fragContigOrientation,
													chunkContigOrientation, 
													chunkLeftEnd, chunkRightEnd);
			  fragLocaleMod = 1 + frag->locale % 8;
			  if (contigOrientation == 0)
			  {
				fprintf(cam_file,
						"%dFrag%d: %d A%d %d R6 # containing contig: %d, frag->iid: %d, baseID: %d, locale: %d\n",
						frag->iid,   // the 1000000 prevents collisons in the Celamy namespace
						originalUnitig->id,
						scaffoldOffset + contigLeftEnd + min,
						fragLocaleMod, // 1 + frag->locale,  // color 1 are the non-BAC frags
						scaffoldOffset + contigLeftEnd + max,
						frag->contigID,
						frag->iid,
						unitig->info.CI.baseID,
						frag->locale);
				fragsMinPos[frag->iid] = scaffoldOffset + contigLeftEnd + min;
				fragsMaxPos[frag->iid] = scaffoldOffset + contigLeftEnd + max;
			  }
			  else
			  {
				fprintf(cam_file,
						"%dFrag%d: %d A%d %d R6 # containing contig: %d, frag->iid: %d, baseID: %d, locale: %d\n",
						frag->iid,   // the 1000000 prevents collisons in the Celamy namespace
						originalUnitig->id,
						scaffoldOffset + contigLeftEnd + (int) contig->bpLength.mean - max,
						fragLocaleMod, // 1 + frag->locale,  // color 1 are the non-BAC frags
						scaffoldOffset + contigLeftEnd + (int) contig->bpLength.mean - min,
						frag->contigID,
						frag->iid,
						unitig->info.CI.baseID,
						frag->locale);
				fragsMinPos[frag->iid] = scaffoldOffset + contigLeftEnd + (int) contig->bpLength.mean - max;
				fragsMaxPos[frag->iid] = scaffoldOffset + contigLeftEnd + (int) contig->bpLength.mean - min;
			  }

			  fragsPresent[frag->iid] = 1;
			  fragsScaffold[frag->iid] = scaff->id;			  
			  fragsLocale[frag->iid] = frag->locale;
			  fragsLocaleMin[frag->iid] = frag->localePos.bgn;
			  fragsLocaleMax[frag->iid] = frag->localePos.end;
			}
			else  // look at non-locale frags if they are within a range from the end of the contig
			{
			  int 
                            contigOrientation,
                            fragOrientationInContig,
                            basesFromEnd;

			  int fragDisplayed;
			  
			  if (contig->offsetAEnd.mean < contig->offsetBEnd.mean)
				contigOrientation = 0;
			  else
				contigOrientation = 1;
			  
			  if (frag->contigOffset5p.mean < frag->contigOffset3p.mean)
				fragOrientationInContig = 0;
			  else
				fragOrientationInContig = 1;
			  
			  fragDisplayed = FALSE;

			  // do A end first
			  basesFromEnd = min((int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean);
			  if ( basesFromEnd < BASES_FROM_END_CUTOFF)
			  {
				if (contigOrientation == 0)
				{
				  if (fragOrientationInContig == 0)
				  {
					min = contigOffset + frag->contigOffset5p.mean;
					max = contigOffset + frag->contigOffset3p.mean;
				  }
				  else
				  {
					min = contigOffset + frag->contigOffset3p.mean;
					max = contigOffset + frag->contigOffset5p.mean;
				  }
				}
				else if (contigOrientation == 1)
				{
				  if (fragOrientationInContig == 0)
				  {
					min = contigOffset + contig->bpLength.mean - frag->contigOffset3p.mean;
					max = contigOffset + contig->bpLength.mean - frag->contigOffset5p.mean;
				  }
				  else
				  {
					min = contigOffset + contig->bpLength.mean - frag->contigOffset5p.mean;
					max = contigOffset + contig->bpLength.mean - frag->contigOffset3p.mean;
				  }			
				}
				fprintf(cam_file,
						"%dRead%d: %d A%d %d R4 # containing contig: %d frag->iid: %d\n",
						readCnt++,
						frag->iid,   // the 2000000 prevents collisons in the Celamy namespace
						scaffoldOffset + min,
						1,
						scaffoldOffset + max,
						frag->contigID,
						frag->iid);
				// note that we already displayed this frag
				fragDisplayed = TRUE;
			  }

			  // now do B end
			  basesFromEnd = min(contig->bpLength.mean - (int) frag->contigOffset5p.mean,
								 contig->bpLength.mean - (int) frag->contigOffset3p.mean);
			  if ( basesFromEnd < BASES_FROM_END_CUTOFF && !fragDisplayed)
			  {
				if (contigOrientation == 0)
				{
				  if (fragOrientationInContig == 0)
				  {
					min = contigOffset + frag->contigOffset5p.mean;
					max = contigOffset + frag->contigOffset3p.mean;
				  }
				  else
				  {
					min = contigOffset + frag->contigOffset3p.mean;
					max = contigOffset + frag->contigOffset5p.mean;
				  }
				}
				else if (contigOrientation == 1)
				{
				  if (fragOrientationInContig == 0)
				  {
					min = contigOffset + contig->bpLength.mean - frag->contigOffset3p.mean;
					max = contigOffset + contig->bpLength.mean - frag->contigOffset5p.mean;
				  }
				  else
				  {
					min = contigOffset + contig->bpLength.mean - frag->contigOffset5p.mean;
					max = contigOffset + contig->bpLength.mean - frag->contigOffset3p.mean;
				  }			
				}
				fprintf(cam_file,
						"%dRead%d: %d A%d %d R4 # containing contig: %d frag->iid: %d\n",
						readCnt++,
						frag->iid,   // the 3000000 prevents collisons in the Celamy namespace
						scaffoldOffset + min,
						1,
						scaffoldOffset + max,
						frag->contigID,
						frag->iid);
			  }		
			}		
		  }
		}
	  }
	}
	fprintf( cam_file, "LNK: ");
	for (i = 0; i < contigCnt; i++)
	  fprintf( cam_file, "%d ", contigIDs[i]);
	fprintf( cam_file, "A0ScaffoldColor\n");
	
	scaffoldOffset += max( (int) contig->offsetAEnd.mean, (int) contig->offsetBEnd.mean) + 20;
	if (0) fprintf( stderr, "setting scaffoldOffset to %d\n", scaffoldOffset);
  }

  icnt = 0;
  while (fragsPresent[icnt++] == 0);
  lastFragIid = lastFragIidPerm = icnt - 1;
  lastFragMinPos = fragsMinPos[lastFragIid];
  lastFragMaxPos = fragsMaxPos[lastFragIid];
  lastFragScaffoldId = fragsScaffold[lastFragIid];
  lastFragLocale = fragsLocale[lastFragIid];
	
  for (icnt = lastFragIid + 1; icnt < MAX_FRAG_IID; icnt++)
  {
	if (!fragsPresent[icnt])
	  continue;
	
	if ( (icnt - lastFragIid > 5) && lastFragScaffoldId == fragsScaffold[icnt] &&
		 lastFragLocale == fragsLocale[icnt])
	{
	  if (lastFragMaxPos < fragsMinPos[icnt])
		fprintf( cam_file, "%dBactigGap%d: %d A0BactigGapColor %d R8 # jump from frag %d to frag %d\n", 
				 lastFragIid, icnt, 
				 lastFragMaxPos, fragsMinPos[icnt],
				 lastFragIid, icnt);
	  else if (lastFragMinPos > fragsMaxPos[icnt])
		fprintf( cam_file, "%dBactigGap%d: %d A0BactigGapColor %d R8 # jump from frag %d to frag %d\n", 
				 lastFragIid, icnt, 
				 fragsMaxPos[icnt], lastFragMinPos, 
				 lastFragIid, icnt);
	}
	lastFragIid = icnt;
	lastFragMinPos = fragsMinPos[icnt];
	lastFragMaxPos = fragsMaxPos[icnt];
	lastFragScaffoldId = fragsScaffold[icnt];
	lastFragLocale = fragsLocale[icnt];
  }
  
  for (icnt = lastFragIidPerm + 1; icnt < MAX_FRAG_IID - 1; icnt++)
  {
	int localeDiff, scaffDiff;
	
	if ( !fragsPresent[icnt] || !fragsPresent[icnt + 1] )
	  continue;
	
	localeDiff = abs( fragsLocaleMin[icnt] - fragsLocaleMin[icnt+1] );
	scaffDiff = abs( fragsMinPos[icnt] - fragsMinPos[icnt+1] );

	if (abs(scaffDiff - localeDiff) > 10)
	{
	  if ( fragsMinPos[icnt] < fragsMinPos[icnt+1])  // icnt to left of icnt + 1
		fprintf( cam_file, "%dFragShift%d: %d A0FragShiftColor %d R8 # scaff length:%d locale length: %d\n",
				 icnt, icnt + 1, 
				 fragsMinPos[icnt], fragsMaxPos[icnt+1],
				 fragsMaxPos[icnt+1] - fragsMinPos[icnt], fragsLocaleMax[icnt+1] - fragsLocaleMin[icnt]);
	  else 
		fprintf( cam_file, "%dFragShift%d: %d A0FragShiftColor %d R8 # scaff length:%d locale length: %d\n",
				 icnt, icnt + 1, 
				 fragsMinPos[icnt+1], fragsMaxPos[icnt],
				 fragsMaxPos[icnt] - fragsMinPos[icnt+1], fragsLocaleMax[icnt+1] - fragsLocaleMin[icnt]);
	}
  }
  fclose (cam_file);
}

void localeSimCam() 
{
  int
    i,
    numFrags;
  int32
    min,
    max;
  char
    filename[STR_LEN],
    * Colour[NUM_COLORS] = {
      "CFFFF00 T2 S # chunk",
      "C0AAAA0 T2 S # non-BAC frag",
      "C0000F0 T2 S # locale_1",
      "CAAAA00 T2 S # locale_2",
      "C00FF00 T2 S # locale_3",
      "CFF8000 T2 S # locale_4",
      "CA0A0FF T2 S # locale_5",
      "CA00A0F T2 S # locale_6",
      "CFFA0F0 T2 S # locale_7",
      "C00AAAA T1 S # locale_8"};
  
  FILE
    * cam_file;
  MultiAlignT 
	*ma;
  IntMultiPos 
	*mp;
  CIFragT 
	*frag;
  // CIScaffoldTIterator CIs;
  NodeCGW_T 
	*contig;
  GraphNodeIterator
	ContigGraphIterator;

  //
  // open the cam file
  //
  sprintf(filename, "./cam/locales.sim.cam");
  // sprintf(filename, "./cam/locale.%d.cam", locale);
  cam_file = file_open (filename, "w");
  assert(cam_file != NULL);
  
  //
  // output the colors
  //
  for (i = 0; i < NUM_COLORS; i++)
    fprintf(cam_file, "%d: %s\n", i, Colour[i]);

  
  
  // InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE, FALSE, &CIs);
  // while(NextCIScaffoldTIterator(&CIs));
  
  InitGraphNodeIterator(&ContigGraphIterator, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((contig = NextGraphNodeIterator(&ContigGraphIterator)) != NULL)
  {
	int offset;

	ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, ScaffoldGraph->RezGraph->type == CI_GRAPH); 
	// ma = GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, contig->id);
	numFrags = GetNumIntMultiPoss(ma->f_list);
	
	min = min( contig->aEndCoord, contig->bEndCoord);
	max = max( contig->aEndCoord, contig->bEndCoord);  
	
	// if (contains_fbac(contig))
	{
	  fprintf(cam_file,
			  "%d: %d A%d %d # contig: %d, containing scaffold: %d\n",
			  contig->id,
			  min,
			  0,
			  max,
			  contig->id,
			  contig->scaffoldID);
	
	  offset = 0;
	  
	  for(i = 0; i < 1 /* numFrags */; i++)
	  {
		int fragLocaleMod;
		
		mp   = GetIntMultiPos(ma->f_list, i);
		frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32)mp->source);
		
		if( frag->locale != -1 )
		{
		  if (frag->aEndCoord < frag->bEndCoord)
		  {
			min = frag->aEndCoord;
			max = frag->bEndCoord;
		  }
		  else
		  {
			min = frag->bEndCoord;
			max = frag->aEndCoord;
		  }

		  fragLocaleMod = 1 + frag->locale % 8;
		  
		  fprintf(cam_file,
				  "%d: %d A%d %d # containing contig: %d, locale: %d\n",
				  frag->iid,
				  min,
				  fragLocaleMod,
				  max,
				  frag->contigID,
				  frag->locale);
		}
	  }
	}
  }
  fclose (cam_file);
}

void GetFragmentPositionInContigFromChunk(CIFragT *frag, int *left_end, int *right_end, 
										  int *fragmentContigOrientation, 
										  int chunkContigOrientation, int chunkLeftEnd, int chunkRightEnd)
{
  //fprintf(stderr, "for frag %d, chunkLeftEnd: %d, chunkRightEnd: %d\n", frag->iid, chunkLeftEnd, chunkRightEnd);
  //fprintf(stderr, "for frag %d, frag->offset5p.mean: %f, frag->offset3p.mean : %f\n", 
  //	  frag->iid, frag->offset5p.mean, frag->offset3p.mean);
  
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
}

void DumpContigGraph( int scaffID )
{
  CIEdgeT* e;
  GraphEdgeIterator edges;
  CIScaffoldT *scaffold;
  CIScaffoldTIterator Contigs;
  ChunkInstanceT *contig;
  FILE * out;
  char buffer[256];
  
  sprintf( buffer, "scaffold_graph_%d", scaffID);
  out = fopen( buffer, "w");
  if (out == NULL)
  {
	fprintf( stderr, "could not open file for writing: %s\n", buffer);
	exit(1);
  }

  // grab the scaffold of interest
  scaffold = GetGraphNode( ScaffoldGraph->ScaffoldGraph, scaffID);  

  // iterate through the contigs
  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &Contigs);
  while (NextCIScaffoldTIterator(&Contigs))
  {
	contig = GetGraphNode(ScaffoldGraph->RezGraph, Contigs.curr);

	// now get some edge info
	InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, contig->id,
						  ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &edges);
	while((e = NextGraphEdgeIterator(&edges)) != NULL)
	{
	  assert(e != NULL);
	  // if (e->flags.bits.isEssential && !e->flags.bits.isDeleted) note: isEssential is never set for some reason?
	  if (e->flags.bits.isUniquetoUnique && !e->flags.bits.isDeleted)
	  {
            // NodeCGW_T * left, *right;
            // left  = GetGraphNode(ScaffoldGraph->ContigGraph, e->idA);
            // right = GetGraphNode(ScaffoldGraph->ContigGraph, e->idB);
            if (e->idA == contig->id)
              fprintf( out, "%d %d %f\n",
                       e->idA, e->idB, (double) e->edgesContributing);
	  }
	}
  }
  fclose(out);
}

void CheckScaffoldOverlaps( void ) 
{
  ChunkInstanceT
    * lchunk,
    * rchunk;
  Overlap *tempOlap1;
  ChunkOrientationType overlapOrientation;
  int computedAhang;
  char *seq1, *seq2;
  int min_ahang, max_ahang;
  double erate, thresh;
  int minlen;
  int sc;
  
  fprintf( stderr, "=== Checking scaffold overlaps\n");
  
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
	
    // print some stats of the scaffold
#if 0
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
    fprintf( stderr, "=== Checking scaffold %d\n", sc);
	
    //
    //
    InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE,	FALSE, &CIsTemp);
    while (NextCIScaffoldTIterator(&CIsTemp))
	{
	  // not walking off of scaffolds currently
	  if (CIsTemp.next == NULLINDEX)
		break;
	  
      lchunk = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.curr);
	  rchunk = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.next);
	  assert(lchunk != NULL);
	  assert(rchunk != NULL);

	  while (min ( rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean) < 
			 max ( lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean))
	  {
		// determine chunks overlap size
		double computedOverlap = max ( lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean) -
		  min ( rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean);

		if (computedOverlap <= 40.0)
		  break;
		
		// figure out orientation
		if ( lchunk->offsetAEnd.mean < lchunk->offsetBEnd.mean)
		{
		  if ( rchunk->offsetAEnd.mean < rchunk->offsetBEnd.mean)
			overlapOrientation = AB_AB;
		  else
			overlapOrientation = AB_BA;
		}
		else
		{
		  if ( rchunk->offsetAEnd.mean < rchunk->offsetBEnd.mean)
			overlapOrientation = BA_AB;
		  else
			overlapOrientation = BA_BA;
		}
		
		// compute expected ahang
		computedAhang = min ( rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean) -
		  min (lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean);
		min_ahang = computedAhang - 30;
		max_ahang = computedAhang + 30;
		
		if(consensus1 == NULL)
		{
		  consensus1 = CreateVA_char(2048);
		  consensus2 = CreateVA_char(2048);
		  quality1 = CreateVA_char(2048);
		  quality2 = CreateVA_char(2048);
		}
		
		// Get the consensus sequences for both chunks from the Store
		GetConsensus(ScaffoldGraph->RezGraph, lchunk->id, consensus1, quality1);
		GetConsensus(ScaffoldGraph->RezGraph, rchunk->id, consensus2, quality2);
		
		seq1 = Getchar(consensus1, 0);
		seq2 = Getchar(consensus2, 0);
		
		erate = CGW_DP_ERATE;
		thresh = CGW_DP_THRESH;
		minlen = CGW_DP_MINLEN;
		
		// fprintf( stderr, "CheckScaffoldOverlaps: lchunk: %d, rchunk: %d, computedAhang: %d, orientation: %c\n",
		//	 lchunk->id, rchunk->id, computedAhang, (char) overlapOrientation);
		
		// tempOlap1 points to a static down inside of DP_Compare
		tempOlap1 = OverlapSequences( seq1, seq2, overlapOrientation, min_ahang, max_ahang, 
									  erate, thresh, minlen, AS_FIND_ALIGN);
		
		if (tempOlap1 == NULL)
		{
		  fprintf( stderr, "Warning: contig %d and %d do not have a valid overlap!!!\n",
			lchunk->id, rchunk->id);
		  fprintf( stderr, "contig %d position: ( %f, %f)\n", 
				   lchunk->id, lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean);
		  fprintf( stderr, "contig %d position: ( %f, %f)\n",		  
				   rchunk->id, rchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean);
		}
		rchunk = GetGraphNode(ScaffoldGraph->RezGraph, rchunk->BEndNext);
		assert(rchunk != NULL);
	  }
	}
  }

  fprintf( stderr, "=== Done checking scaffold overlaps\n");
  
  return;
}

int GetChunkPositionInContig(ChunkInstanceT *chunk, int *left_end, int *right_end, 
							 int *chunkContigOrientation)
{  
  if (chunk->offsetAEnd.mean <= chunk->offsetBEnd.mean) 
  {
	*left_end = chunk->offsetAEnd.mean;
	*right_end = chunk->offsetBEnd.mean;
	*chunkContigOrientation = 0;
  } 
  else 
  {
	*left_end = chunk->offsetBEnd.mean;
	*right_end = chunk->offsetAEnd.mean;
	*chunkContigOrientation = 1;
  }
  return 0;
}

void CheckAlteredScaffolds( int *checkScaffold, int checkScaffoldCount)
{
  int icnt;

  fprintf( stderr, "checking %d altered scaffolds for connectivity\n", checkScaffoldCount);
  fflush(stderr);
  for (icnt = 0; icnt < checkScaffoldCount; icnt++)
  {
	fprintf( stderr, "going to check scaffold %d for connectivity\n", checkScaffold[icnt]);
  }
  fflush(stderr);
  
  for (icnt = 0; icnt < checkScaffoldCount; icnt++)
  {
	CIScaffoldT *scaff= GetCIScaffoldT(ScaffoldGraph->CIScaffolds, checkScaffold[icnt]);
	
	fprintf( stderr, "checking scaffold %d for connectivity\n", scaff->id);
	
	CheckScaffoldConnectivityAndSplit( ScaffoldGraph, scaff, ALL_EDGES, FALSE);
  }
}

void removeSmallContigs( CIScaffoldT *scaff )
{
  CIScaffoldTIterator CIsTemp;
  int numContigs = scaff->info.Scaffold.numElements;
  int *allContigs = (int *) malloc( scaff->info.Scaffold.numElements * sizeof(int) );
  int *removedContigs = (int *) malloc( scaff->info.Scaffold.numElements * sizeof(int) );
  int icnt = 0;
  ChunkInstanceT *lcontig, *mcontig, *rcontig;
  int numToRemove = 0;
  
  fprintf( stderr, "Doing small contig removal on scaffold: %d\n", scaff->id);

  if (scaff->info.Scaffold.numElements < 3)
	return;
  
  // first find all the contigs in the gap
  InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE, FALSE, &CIsTemp);
  while (NextCIScaffoldTIterator(&CIsTemp)) 
  {
	rcontig = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.curr);
	assert(rcontig != NULL);

	fprintf( stderr, "in scaffold %d looking at contig %d, length: %f, contains bac frags: %d\n", 
		 scaff->id, rcontig->id, rcontig->bpLength.mean, contains_fbac(rcontig));
	
	if (rcontig->bpLength.mean < 2700) // && !contains_fbac(rcontig))
	{
	  numToRemove++;
	  removedContigs[icnt] = 1;
	}
	else
	  removedContigs[icnt] = 0;

	allContigs[icnt++] = CIsTemp.curr;
  } 

  if (scaff->info.Scaffold.numElements - numToRemove < 2)
	return;
  
  // now step through contigs examining triplets where we might want to remove the middle
  for (icnt = 1; icnt < numContigs - 1; icnt++)
  {
    int lcontigGapEnd, rcontigGapEnd;
    GapInfoT *gapInfoArray;
    int numCommonLocales;

    if (!removedContigs[icnt])
      continue;
    mcontig = GetGraphNode(ScaffoldGraph->RezGraph, allContigs[icnt]);

#if 1  // we used to worry about whether the flanking contigs had a common locale	
    fprintf( stderr, "in scaffold %d, triplet: (%d, %d, %d)\n", 
	     scaff->id, allContigs[icnt - 1], allContigs[icnt], allContigs[icnt + 1]);
    
    lcontig = GetGraphNode(ScaffoldGraph->RezGraph, allContigs[icnt - 1]);
    rcontig = GetGraphNode(ScaffoldGraph->RezGraph, allContigs[icnt + 1]);
    
    if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
      lcontigGapEnd = B_END;
    else
      lcontigGapEnd = A_END;
    
    if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
      rcontigGapEnd = A_END;
    else
      rcontigGapEnd = B_END;
    
    // find the locales that are at the apropos end of the contigs
    numCommonLocales = FindCommonLocales( lcontig, lcontigGapEnd, rcontig, rcontigGapEnd, &gapInfoArray);
    
    if ( numCommonLocales > 0)
#endif
	{
	  // fprintf( stderr, "in scaffold %d, contig %d is flanked by contigs (%d, %d) that have %d common locales\n",
	  //   scaff->id, mcontig->id, lcontig->id, rcontig->id, numCommonLocales);
	  fprintf( stderr, "in scaffold %d, removing contig %d\n", scaff->id, mcontig->id);
	  RemoveCIFromScaffold( ScaffoldGraph, scaff, mcontig, TRUE);
	  mcontig->flags.bits.isPotentialStone = 0;  // this prevents Art from throwing this guy back in

	  if (mcontig->info.Contig.numCI > 1)  // has to become a singleton scaffold
	  {
		CIScaffoldT* newScaffold = (CIScaffoldT *) malloc( sizeof (CIScaffoldT));
		int32 newScaffoldID;
		LengthT NullLength = {0.0, 0.0};

		assert (newScaffold != NULL);
		
		if (1)
		  newScaffold = CreateNewScaffold();
		else
		{
		  InitializeScaffold( newScaffold, REAL_SCAFFOLD);
		  newScaffold->info.Scaffold.AEndCI = NULLINDEX;
		  newScaffold->info.Scaffold.BEndCI = NULLINDEX;
		  newScaffold->info.Scaffold.numElements = 0;
		  newScaffold->edgeHead = NULLINDEX;
		  newScaffold->bpLength = NullLength;
		  newScaffoldID = newScaffold->id = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
		  newScaffold->flags.bits.isDead = FALSE;
		  newScaffold->aEndCoord = newScaffold->bEndCoord = -1;
		  newScaffold->numEssentialA = newScaffold->numEssentialB = 0;
		  newScaffold->essentialEdgeB = newScaffold->essentialEdgeA = NULLINDEX;
		  AppendGraphNode( ScaffoldGraph->ScaffoldGraph, newScaffold);
		}
		
		InsertCIInScaffold(ScaffoldGraph, mcontig->id, newScaffold->id,
						   NullLength, mcontig->bpLength, TRUE, FALSE);
	  }
	  else
	  {
		NodeCGW_T *unitig;
		
		assert( mcontig->info.Contig.AEndCI == mcontig->info.Contig.BEndCI);
		unitig = GetGraphNode( ScaffoldGraph->CIGraph, mcontig->info.Contig.AEndCI);
		SetNodeType(unitig, UNRESOLVEDCHUNK_CGW);
	  }
	}
  }
}

void trimScaffoldEnds( CIScaffoldT *scaff )
{
  ChunkInstanceT *contig, *mcontig;
  int removeContig;
  
  fprintf( stderr, "Doing end trimming on scaffold: %d\n", scaff->id);

  if (scaff->info.Scaffold.numElements < 2)
	return;
  
  // examine the contig on the A End
  removeContig = FALSE;
  contig = GetGraphNode( ScaffoldGraph->ContigGraph, scaff->info.Scaffold.AEndCI );

  // is it small and shredded BAC fragment free?
  if ( contig->bpLength.mean < 2700 && !contains_fbac(contig) )
  {
	// see if the next contig over contains BAC frags
	mcontig = GetGraphNode( ScaffoldGraph->ContigGraph, contig->BEndNext);
	
	if ( mcontig->bpLength.mean > 2700 && contains_fbac( mcontig ))
	  removeContig = TRUE;
  }
  
  if (removeContig == TRUE)  // get rid of contig
  {
	fprintf( stderr, "in scaffold %d, removing contig %d\n", scaff->id, contig->id);
	RemoveCIFromScaffold( ScaffoldGraph, scaff, contig, TRUE);  

	if (contig->info.Contig.numCI > 1)  // has to become a singleton scaffold
	{
	  CIScaffoldT* newScaffold = (CIScaffoldT *) malloc( sizeof (CIScaffoldT));
	  int32 newScaffoldID;
	  LengthT NullLength = {0.0, 0.0};
	  
	  assert (newScaffold != NULL);
	  
	  if (1)
		newScaffold = CreateNewScaffold();
	  else
	  {
		InitializeScaffold( newScaffold, REAL_SCAFFOLD);
		newScaffold->info.Scaffold.AEndCI = NULLINDEX;
		newScaffold->info.Scaffold.BEndCI = NULLINDEX;
		newScaffold->info.Scaffold.numElements = 0;
		newScaffold->edgeHead = NULLINDEX;
		newScaffold->bpLength = NullLength;
		newScaffoldID = newScaffold->id = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
		newScaffold->flags.bits.isDead = FALSE;
		newScaffold->aEndCoord = newScaffold->bEndCoord = -1;
		newScaffold->numEssentialA = newScaffold->numEssentialB = 0;
		newScaffold->essentialEdgeB = newScaffold->essentialEdgeA = NULLINDEX;
		AppendGraphNode( ScaffoldGraph->ScaffoldGraph, newScaffold);
	  }
	  InsertCIInScaffold(ScaffoldGraph, contig->id, newScaffold->id,
						 NullLength, contig->bpLength, TRUE, FALSE);
	}
	else
	{
	  NodeCGW_T *unitig;
		
	  assert( contig->info.Contig.AEndCI == contig->info.Contig.BEndCI);
	  unitig = GetGraphNode( ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);
	  SetNodeType(unitig, UNRESOLVEDCHUNK_CGW);
	}
  }

  // examine the contig on the B End
  removeContig = FALSE;
  contig = GetGraphNode( ScaffoldGraph->ContigGraph, scaff->info.Scaffold.BEndCI );

  // is it small and shredded BAC fragment free?
  if ( contig->bpLength.mean < 2700 && !contains_fbac(contig) )
  {
	// see if the next contig over contains BAC frags
	mcontig = GetGraphNode( ScaffoldGraph->ContigGraph, contig->AEndNext);
	
	if ( mcontig->bpLength.mean > 2700 && contains_fbac( mcontig ))
	  removeContig = TRUE;
  }
  
  if (removeContig == TRUE)  // get rid of contig
  {
	fprintf( stderr, "in scaffold %d, removing contig %d\n", scaff->id, contig->id);
	RemoveCIFromScaffold( ScaffoldGraph, scaff, contig, TRUE);  

	if (contig->info.Contig.numCI > 1)  // has to become a singleton scaffold
	{
	  CIScaffoldT* newScaffold = (CIScaffoldT *) malloc( sizeof (CIScaffoldT));
	  int32 newScaffoldID;
	  LengthT NullLength = {0.0, 0.0};
	  
	  assert (newScaffold != NULL);
	  
	  if (1)
		newScaffold = CreateNewScaffold();
	  else
	  {
		InitializeScaffold( newScaffold, REAL_SCAFFOLD);
		newScaffold->info.Scaffold.AEndCI = NULLINDEX;
		newScaffold->info.Scaffold.BEndCI = NULLINDEX;
		newScaffold->info.Scaffold.numElements = 0;
		newScaffold->edgeHead = NULLINDEX;
		newScaffold->bpLength = NullLength;
		newScaffoldID = newScaffold->id = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
		newScaffold->flags.bits.isDead = FALSE;
		newScaffold->aEndCoord = newScaffold->bEndCoord = -1;
		newScaffold->numEssentialA = newScaffold->numEssentialB = 0;
		newScaffold->essentialEdgeB = newScaffold->essentialEdgeA = NULLINDEX;
		AppendGraphNode( ScaffoldGraph->ScaffoldGraph, newScaffold);
	  }
	  InsertCIInScaffold(ScaffoldGraph, contig->id, newScaffold->id,
						 NullLength, contig->bpLength, TRUE, FALSE);
	}
	else
	{
	  NodeCGW_T *unitig;
		
	  assert( contig->info.Contig.AEndCI == contig->info.Contig.BEndCI);
	  unitig = GetGraphNode( ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);
	  SetNodeType(unitig, UNRESOLVEDCHUNK_CGW);
	}
  }
}

void analyzeGaps( CIScaffoldT *scaff)
{
  CIScaffoldTIterator CIsTemp;
  ChunkInstanceT *contig;
  
  fprintf( stderr, "analyzing gaps in scaffold: %d\n", scaff->id);

  if (scaff->info.Scaffold.numElements < 2)
	return;
  
  // first find all the contigs in the scaffold
  InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE, FALSE, &CIsTemp);
  while (NextCIScaffoldTIterator(&CIsTemp)) 
  {
	contig = GetGraphNode(ScaffoldGraph->RezGraph, CIsTemp.curr);
	assert(contig != NULL);
  }
}


void dumpWalkContigs ( GapInfoT gapInfo)
{
  int i, increment;
  CIFragT *leftFrag, *rightFrag;

  leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, gapInfo.fragIDLeftContig);  
  rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, gapInfo.fragIDRightContig);

  fprintf( stderr, "for locale %d, leftFrag iid: %d, rightFrag iid: %d\n",
	   gapInfo.localeNumber, gapInfo.fragIDLeftContig, gapInfo.fragIDRightContig);

  if ( leftFrag->iid < rightFrag->iid)
    increment = 1;
  else
    increment = -1;

  for (i = leftFrag->iid; i != rightFrag->iid + increment; i += increment)
  {
    CIFragT *frag;
    InfoByIID *info;    
    NodeCGW_T *contig;
    int lastContigID = -1;

    // get frags position in contig	
    info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, i);
    assert(info->set);
    frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

    if (frag->contigID != lastContigID)
    {
      fprintf( stderr, "    frag: %6d\t contig id: %d\t contig pos (5p, 3p): %6d, %6d\n", 
	       i, frag->contigID, (int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean);	  
      contig = GetGraphNode( ScaffoldGraph->RezGraph, frag->contigID);
      fprintf( stderr, "lastContigID: %d\n", lastContigID);
      lastContigID = frag->contigID;
      fprintf( stderr, "lastContigID: %d\n", lastContigID);
      dumpContigInfo( contig );
    }
    fflush(stderr);
  }
}

void walkScaffolds( CIScaffoldT *scaff1,
                    int invertScaff1, CIScaffoldT *scaff2,
                    int invertScaff2, double gapSize)
{
  CIScaffoldT* newScaffold;
  int32 newScaffoldID;
  int numContigsScaff1, numContigsScaff2;
  ChunkInstanceT *contig;
  NodeOrient orientContig;
  LengthT NullLength = {0.0, 0.0};
  LengthT scaff2Pos;
  int scaff1ID = scaff1->id, scaff2ID = scaff2->id;

  if (1)
	newScaffold = CreateNewScaffold();
  else
  {
	// create & init a new scaffold
	newScaffold = (CIScaffoldT *) malloc( sizeof (CIScaffoldT));
	assert (newScaffold != NULL);  
	InitializeScaffold( newScaffold, REAL_SCAFFOLD);
	newScaffold->info.Scaffold.AEndCI = NULLINDEX;
	newScaffold->info.Scaffold.BEndCI = NULLINDEX;
	newScaffold->info.Scaffold.numElements = 0;
	newScaffold->edgeHead = NULLINDEX;
	newScaffold->bpLength = NullLength;
	newScaffoldID = newScaffold->id = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
	newScaffold->flags.bits.isDead = FALSE;
	newScaffold->aEndCoord = newScaffold->bEndCoord = -1;
	newScaffold->numEssentialA = newScaffold->numEssentialB = 0;
	newScaffold->essentialEdgeB = newScaffold->essentialEdgeA = NULLINDEX;
	AppendGraphNode( ScaffoldGraph->ScaffoldGraph, newScaffold);
  }
  
  // we need to save info on the scaffolds in case the walk doesn't work out
  numContigsScaff1 = scaff1->info.Scaffold.numElements;
  numContigsScaff2 = scaff2->info.Scaffold.numElements;

  fprintf( stderr, "inserting scaff %d into scaff %d at %f\n", scaff1->id, newScaffold->id, 0.0);
  
  orientContig = (invertScaff1 == TRUE) ? B_A : A_B;
  {
    LengthT trashLength = NullLength;
    InsertScaffoldContentsIntoScaffold( ScaffoldGraph,
                                        newScaffold->id,
                                        scaff1->id,
                                        orientContig,
                                        &trashLength,
                                        TRUE ); // contig now
  }
  scaff1 = GetGraphNode( ScaffoldGraph->ScaffoldGraph, scaff1ID);
  scaff1->flags.bits.isDead = TRUE;
  scaff1->info.Scaffold.numElements = 0;
  
  // reget pointer since InsertScaffoldContentsIntoScaffold may have realloced nodes
  newScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, newScaffold->id);
  
  // add gapSize onto extreme contig position, don't worry about variance, it gets redone after contigging
  contig = GetGraphNode( ScaffoldGraph->ContigGraph, newScaffold->info.Scaffold.BEndCI);
  scaff2Pos = (contig->offsetAEnd.mean < contig->offsetBEnd.mean) ? contig->offsetBEnd : contig->offsetAEnd;
  scaff2Pos.mean += gapSize; 

  fprintf( stderr, "inserting scaff %d into scaff %d at %f\n", scaff2->id, newScaffold->id, scaff2Pos.mean);

  orientContig = (invertScaff2 == TRUE) ? B_A : A_B;
  InsertScaffoldContentsIntoScaffold( ScaffoldGraph,
                                      newScaffold->id,
                                      scaff2->id,
                                      orientContig,
                                      &scaff2Pos,
                                      TRUE); // contig now!
  scaff2 = GetGraphNode( ScaffoldGraph->ScaffoldGraph, scaff2ID);
  scaff2->flags.bits.isDead = TRUE;
  scaff2->info.Scaffold.numElements = 0;

  Intra_Scaffold_Path_Finding( newScaffold->id, 
                               AGGRESSIVE_WALKING_STD_DEVS,
                               // CONSERVATIVE_WALKING_STD_DEVS,
                               FALSE,
                               FALSE);
}

void SaveContigInformation( CIScaffoldT *scaffold, ScaffoldInfoT *scaffoldContigs)
{
  int contigCount, numContigs = scaffold->info.Scaffold.numElements;
  CIScaffoldTIterator Contigs;
  
  scaffoldContigs->numContigs = numContigs;
  scaffoldContigs->ContigIDs = (int *) malloc( numContigs * sizeof( int ));
  scaffoldContigs->AEndPositions = (LengthT *) malloc( numContigs * sizeof( LengthT ));
  scaffoldContigs->BEndPositions = (LengthT *) malloc( numContigs * sizeof( LengthT ));
  
  if ( scaffoldContigs->ContigIDs == NULL || 
	   scaffoldContigs->AEndPositions == NULL ||
	   scaffoldContigs->BEndPositions == NULL)
  {
	fprintf( stderr, "memory allocation failure in SaveContigInformation!\n");
	assert(0);
  }
  
  contigCount = 0;
  InitCIScaffoldTIterator( ScaffoldGraph, scaffold, TRUE, FALSE, &Contigs);
  while ( NextCIScaffoldTIterator(&Contigs))
  {
	ContigT *contig;
	
	contig = GetGraphNode( ScaffoldGraph->ContigGraph, Contigs.curr);
	scaffoldContigs->ContigIDs[contigCount] = contig->id;
	scaffoldContigs->AEndPositions[contigCount] = contig->offsetAEnd;
	scaffoldContigs->BEndPositions[contigCount++] = contig->offsetBEnd;
  }
}

void RestoreContigInformation( CIScaffoldT *scaffold, ScaffoldInfoT *scaffoldContigs)
{
  int contigCount;
  ContigT *contig;
  
  for (contigCount = 0; contigCount < scaffoldContigs->numContigs; contigCount++)
  {
	contig = GetGraphNode( ScaffoldGraph->ContigGraph, scaffoldContigs->ContigIDs[contigCount]);

	InsertCIInScaffold(ScaffoldGraph, contig->id, scaffold->id,
					   scaffoldContigs->AEndPositions[contigCount], scaffoldContigs->BEndPositions[contigCount],
					   TRUE, NO_CONTIGGING);
  }
}

void FreeContigInformation( ScaffoldInfoT *scaffoldContigs)
{
  free( scaffoldContigs->ContigIDs);
  free( scaffoldContigs->AEndPositions);
  free( scaffoldContigs->BEndPositions);  
}



#ifdef DEAD_CODE
*******************************************************************************
*******************************************************************************
*******************************************************************************
*******************************************************************************
*******************************************************************************


**** int localeCamOrig() 
**** {
****   int
****     color,
****     k,
****     o,
****     id,
****     i,
**** 	numFrags,
**** 	contigNumber;
****   CIEdgeT
****     * edge;
****   int32
****     low,
****     high,
****     min,
****     max,
****     a_end,
****     b_end,
****     cid,
****     other_cid;
****   ChunkInstanceT
****     * chunk;
****   char
****     unique[STR_LEN],
****     orientation[STR_LEN],
****     filename[STR_LEN],
****     * Colour[NUM_COLORS] = {
****       "CFFFF00 T2 S # contig",
****       "C0AAAA0 T2 S # non-BAC frag",
****       "C0000F0 T2 S # locale_1",
****       "CAAAA00 T2 S # locale_2",
****       "C00FF00 T2 S # locale_3",
****       "CFF8000 T2 S # locale_4",
****       "CA0A0FF T2 S # locale_5",
****       "CA00A0F T2 S # locale_6",
****       "CFFA0F0 T2 S # locale_7",
****       "C00AAAA T1 S # locale_8"};
****   FILE
****     * cam_file;
****   GraphCGW_T 
**** 	*graph = ScaffoldGraph->RezGraph;
****   MultiAlignT 
**** 	*ma;
****   IntMultiPos 
**** 	*mp;
****   CIFragT 
**** 	*frag;
****   CIScaffoldTIterator
**** 	CIs;
****   CIScaffoldT 
**** 	*scaff1;
****   NodeCGW_T 
**** 	*contig;
****   GraphNodeIterator
**** 	ContigGraphIterator;
**** 
****   //
****   // open the cam file
****   //
****   sprintf(filename, "./cam/locales.cam");
****   // sprintf(filename, "./cam/locale.%d.cam", locale);
****   cam_file = file_open (filename, "w");
****   assert(cam_file != NULL);
****   
****   //
****   // output the colors
****   //
****   for (i = 0; i < NUM_COLORS; i++)
****     fprintf(cam_file, "%d: %s\n", i, Colour[i]);
**** 
****   
****   
****   // InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE, FALSE, &CIs);
****   // while(NextCIScaffoldTIterator(&CIs));
****   
****   InitGraphNodeIterator(&ContigGraphIterator, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
****   while (contig = NextGraphNodeIterator(&ContigGraphIterator))
****   {
**** 	int offset;
**** 
**** 	if (contig->scaffoldID == -1)
**** 	  continue;
**** 	
**** 	ma = GetMultiAlignInStore(graph->maStore, contig->id);
**** 	numFrags = GetNumIntMultiPoss(ma->f_list);
**** 	
**** 	min = min( (int) contig->offsetAEnd.mean, (int) contig->offsetBEnd.mean);
**** 	max = max( (int) contig->offsetAEnd.mean, (int) contig->offsetBEnd.mean);
**** 	
**** 	fprintf(cam_file,
**** 			"%d: %d A%d %d # contig: %d, containing scaffold: %d\n",
**** 			contig->id,
**** 			min,
**** 			0,
**** 			max,
**** 			contig->id,
**** 			contig->scaffoldID);
**** 	
**** 	if( !contains_fbac(contig))
**** 	  continue;
**** 	
**** 	offset = min( (int) contig->offsetAEnd.mean, (int) contig->offsetBEnd.mean);
**** 	for (i = 0; i < numFrags; i++)
**** 	{
**** 	  mp   = GetIntMultiPos(ma->f_list, i);
**** 	  frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32)mp->source);
**** 
**** 
**** 	  /*
		**** 	  if (contig->id == 98) 
		**** 		fprintf( stderr, "%d, %d, contig: %d, frag: %d, frag->locale: %d, celsim coords: %d, %d\n",
		**** 				 min ((int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean),
		**** 				 max ((int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean),
		**** 				 contig->id,
		**** 				 frag->iid,
		**** 				 frag->locale,
		**** 				 frag->aEndCoord,
		**** 				 frag->bEndCoord);
		**** 	  */
**** 
**** 
**** 	  
**** 	  if( frag->locale != -1 )
**** 	  {
**** 		if (frag->contigOffset5p.mean < frag->contigOffset3p.mean)
**** 		{
**** 		  min = offset + frag->contigOffset5p.mean;
**** 		  max = offset + frag->contigOffset3p.mean;
**** 		}
**** 		else
**** 		{
**** 		  min = offset + frag->contigOffset3p.mean;
**** 		  max = offset + frag->contigOffset5p.mean;
**** 		}
**** 		
**** 		fprintf(cam_file,
**** 				"%d: %d A%d %d # containing contig: %d, celsim coords: %d, %d\n",
**** 				frag->iid,
**** 				min,
**** 				1 + frag->locale,  // color 1 are the non-BAC frags
**** 				max,
**** 				frag->contigID,
**** 				frag->aEndCoord,
**** 				frag->bEndCoord);
**** 	  }
**** 	  else
**** 	  {
**** 		int 
**** 		  contigOrientation,
**** 		  fragOrientationInContig,
**** 		  basesFromEnd,
**** 		  end;
**** 		
**** 		if (contig->offsetAEnd.mean < contig->offsetBEnd.mean)
**** 		  contigOrientation = 0;
**** 		else
**** 		  contigOrientation = 1;
**** 		
**** 		if (frag->contigOffset5p.mean < frag->contigOffset3p.mean)
**** 		  fragOrientationInContig = 0;
**** 		else
**** 		  fragOrientationInContig = 1;
**** 
**** 		// do A end first
**** 		basesFromEnd = min((int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean);
**** 		if ( basesFromEnd < )
**** 		{
**** 		  if (contigOrientation == 0)
**** 		  {
**** 			if (fragOrientationInContig == 0)
**** 			{
**** 			  min = offset + frag->contigOffset5p.mean;
**** 			  max = offset + frag->contigOffset3p.mean;
**** 			}
**** 			else
**** 			{
**** 			  min = offset + frag->contigOffset3p.mean;
**** 			  max = offset + frag->contigOffset5p.mean;
**** 			}
**** 		  }
**** 		  else if (contigOrientation == 1)
**** 		  {
**** 			if (fragOrientationInContig == 0)
**** 			{
**** 			  min = offset + contig->bpLength.mean - frag->contigOffset3p.mean;
**** 			  max = offset + contig->bpLength.mean - frag->contigOffset5p.mean;
**** 			}
**** 			else
**** 			{
**** 			  min = offset + contig->bpLength.mean - frag->contigOffset5p.mean;
**** 			  max = offset + contig->bpLength.mean - frag->contigOffset3p.mean;
**** 			}			
**** 		  }
**** 		  fprintf(cam_file,
**** 				  "%d: %d A%d %d # containing contig: %d celsim coords: %d, %d\n",
**** 				  1000000 + frag->iid,   // the 1000000 prevents collisons in the Celamy namespace
**** 				  min,
**** 				  1,
**** 				  max,
**** 				  frag->contigID,
**** 				  frag->aEndCoord,
**** 				  frag->bEndCoord);
**** 		}
**** 
**** 		// now do B end
**** 		basesFromEnd = min(contig->bpLength.mean - (int) frag->contigOffset5p.mean,
**** 						   contig->bpLength.mean - (int) frag->contigOffset3p.mean);
**** 		if ( basesFromEnd < 1000)
**** 		{
**** 		  if (contigOrientation == 0)
**** 		  {
**** 			if (fragOrientationInContig == 0)
**** 			{
**** 			  min = offset + frag->contigOffset5p.mean;
**** 			  max = offset + frag->contigOffset3p.mean;
**** 			}
**** 			else
**** 			{
**** 			  min = offset + frag->contigOffset3p.mean;
**** 			  max = offset + frag->contigOffset5p.mean;
**** 			}
**** 		  }
**** 		  else if (contigOrientation == 1)
**** 		  {
**** 			if (fragOrientationInContig == 0)
**** 			{
**** 			  min = offset + contig->bpLength.mean - frag->contigOffset3p.mean;
**** 			  max = offset + contig->bpLength.mean - frag->contigOffset5p.mean;
**** 			}
**** 			else
**** 			{
**** 			  min = offset + contig->bpLength.mean - frag->contigOffset5p.mean;
**** 			  max = offset + contig->bpLength.mean - frag->contigOffset3p.mean;
**** 			}			
**** 		  }
**** 		  fprintf(cam_file,
**** 				  "%d: %d A%d %d # containing contig: %d celsim coords: %d, %d\n",
**** 				  2000000 + frag->iid,   // the 2000000 prevents collisons in the Celamy namespace
**** 				  min,
**** 				  1,
**** 				  max,
**** 				  frag->contigID,
**** 				  frag->aEndCoord,
**** 				  frag->bEndCoord);
**** 		}		
**** 	  }		
**** 	}
****   }
****   fclose (cam_file);
**** }
#endif

#ifdef CHECKCODE
		if (0)
		{
		  int end;
		  
		  if (*tempOlapOrientation == AB_AB || *tempOlapOrientation == AB_BA)
			end = B_END;
		  else
			end = A_END;

		  {
			GraphEdgeIterator edges;
			CIEdgeT *edge;
			
			fprintf( stderr, "Edges before insert----------------------------------------------------\n");
			InitGraphEdgeIterator(ScaffoldGraph->RezGraph, currFrag->contigID, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
			while( edge = NextGraphEdgeIterator(&edges) )
			  PrintGraphEdge( stderr, ScaffoldGraph->RezGraph, "relativ from idA ", edge, edge->idA);
		  }
		  
		  fprintf( stderr, "Inserting an overlap between %d and %d, orientation: %c\n", 
				   currFrag->contigID, fragSucc1->contigID, (char) *tempOlapOrientation);
		  InsertOverlapInHashTable( tempOlap, currFrag->contigID, fragSucc1->contigID, *tempOlapOrientation);
		  
		  {
			GraphEdgeIterator edges;
			CIEdgeT *edge;
			
			fprintf( stderr, "Edges after insert----------------------------------------------------\n");
			InitGraphEdgeIterator(ScaffoldGraph->RezGraph, currFrag->contigID, end, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
			while( edge = NextGraphEdgeIterator(&edges) )
			PrintGraphEdge( stderr, ScaffoldGraph->RezGraph, "relativ from idA ", edge, edge->idA);
		  }

#endif


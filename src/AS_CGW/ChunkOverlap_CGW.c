
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
static char CM_ID[] = "$Id: ChunkOverlap_CGW.c,v 1.14 2007-05-14 09:27:11 brianwalenz Exp $";

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_HashCommon.h"
#include "AS_UTL_Var.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"    // For DeleteCIOverlapEdge
#include "AS_ALN_aligners.h"
#include "ChunkOverlap_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "UtilsREZ.h"
#include "dpc_CNS.h"

//#define DEBUG_CHUNKOVERLAP 

#define GREEDYDEBUG 0
#define USE_NEW_DP_COMPARE


// this is the initial range we use to compute 
// overlaps with a different min/max range
#define BEGENDSLOP 10


// Forward declarations

int InsertChunkOverlap(ChunkOverlapperT *chunkOverlapper,
                       ChunkOverlapCheckT *olap);

int InsertOverlapInHashTable(Overlap *tempOlap, CDS_CID_t cidA, CDS_CID_t cidB,
                             ChunkOrientationType orientation);

static OverlapMesg *ComputeCanonicalOverlapWithTrace(GraphCGW_T *graph, 
                                                     ChunkOverlapCheckT *canOlap,
                                                     InternalFragMesg *AFR,
                                                     InternalFragMesg *BFR,
                                                     FILE* fp, int iterate);


/* compute the overlap length given the result of DP_Compare */
static int adapt_overlap(OverlapMesg* O, ChunkOverlapCheckT* canOlap,
                         CDS_COORD_t lengthA, CDS_COORD_t lengthB,
                         CDS_COORD_t intended);

/* compute the beg, end and opposite parameters for passing to DP_Compare */
static void prepare_overlap(ChunkOverlapCheckT *canOlap,
                            CDS_COORD_t *intended,
                            CDS_COORD_t * beg, CDS_COORD_t * end,
                            int* opposite,
                            CDS_COORD_t lengthA, CDS_COORD_t lengthB);

/* v*** Functions supporting hash table ***v */

/* Type of function to compare hash table entries */
int CanOlapCmp(const void *cO1, const void *cO2){
  ChunkOverlapSpecT *c1 = (ChunkOverlapSpecT *)cO1;
  ChunkOverlapSpecT *c2 = (ChunkOverlapSpecT *)cO2;
  int diff;

  diff = c1->cidA - c2->cidA;
  if(diff)
    return diff;

  diff = c1->cidB - c2->cidB;
  if(diff)
    return diff;

  if(c1->orientation == c2->orientation)
    return 0;
  else return (c1 - c2);
  

}
int CanOlapHash(const void *cO, int length){
  return  Hash_AS((uint8 *)cO, (uint32) length, 37);
}

/* ^*** Functions supporting hash table ***^ */

/************************************************************************/

ChunkOverlapperT *CreateChunkOverlapper(void){
  ChunkOverlapperT *chunkOverlapper;
  
  // Allocate a chunkOverlapper
  chunkOverlapper = (ChunkOverlapperT *)safe_malloc(sizeof(ChunkOverlapperT));
  AssertPtr(chunkOverlapper);


  // Create a symbol table for the sparse chunk-chunk overlap relationships
  chunkOverlapper->hashTable = CreateHashTable_AS(1000, CanOlapHash, CanOlapCmp);

  // Create the VA for the overlap results
  chunkOverlapper->ChunkOverlaps = AllocateChunkOverlapCheckTHeap(1000);

  return chunkOverlapper;
}

/************************************************************************/

void DestroyChunkOverlapper(ChunkOverlapperT *chunkOverlapper){
  // Blow away the symbol table
  DeleteHashTable_AS(chunkOverlapper->hashTable);
  // Free the chunkOverlapper
  FreeChunkOverlapCheckTHeap(chunkOverlapper->ChunkOverlaps);

}

/************************************************************************/

void  SaveChunkOverlapperToStream(ChunkOverlapperT *chunkOverlapper, FILE *stream){

  /* We will save the ChunkOverlaps to a file, and rebuild the heap after
     we read the file.  In other words, we don't worry about hashTable persistence. */

  HashTable_Iterator_AS iterator;
  void *key, *value;
  int64 numOverlaps = 0;

  // Iterate over all hashtable elements, just to count them

  InitializeHashTable_Iterator_AS(chunkOverlapper->hashTable, &iterator);

  while(NextHashTable_Iterator_AS(&iterator, &key, &value)){
    numOverlaps++;
  }

  AS_UTL_safeWrite(stream, &numOverlaps, "SaveChunkOverlapperToStream", sizeof(int64), 1);

  // Iterate over all hashtable elements, writing

  InitializeHashTable_Iterator_AS(chunkOverlapper->hashTable, &iterator);

  while(NextHashTable_Iterator_AS(&iterator, &key, &value)){
    ChunkOverlapCheckT *olap = (ChunkOverlapCheckT*) value;    
    
    AS_UTL_safeWrite(stream, olap, "SaveChunkOverlapperToStream", sizeof(ChunkOverlapCheckT), 1);
  }
}
/************************************************************************/

ChunkOverlapperT *  LoadChunkOverlapperFromStream(FILE *stream){
  int64 numOverlaps;
  int status;
  int64 overlap;
  ChunkOverlapCheckT olap;
  ChunkOverlapperT *chunkOverlapper;

  // Allocate a chunkOverlapper
  chunkOverlapper = (ChunkOverlapperT *)safe_malloc(sizeof(ChunkOverlapperT));
  AssertPtr(chunkOverlapper);


  // open the chunkStore at chunkStorepath
  

  status = AS_UTL_safeRead(stream, &numOverlaps, "LoadChunkOverlapperFromStream", sizeof(int64), 1);
  assert(status == 1);


  // Create a symbol table for the sparse chunk-chunk overlap relationships
  chunkOverlapper->hashTable = CreateHashTable_AS(1000, CanOlapHash, CanOlapCmp);

  // Create the VA for the overlap results
  chunkOverlapper->ChunkOverlaps = AllocateChunkOverlapCheckTHeap(1000);

  for(overlap = 0; overlap < numOverlaps; overlap++){
    status = AS_UTL_safeRead(stream, &olap, "LoadChunkOverlapperFromStream", sizeof(ChunkOverlapCheckT), 1);
    assert(status == 1);
    assert(olap.errorRate > 0.0);
    InsertChunkOverlap(chunkOverlapper, &olap);
  }

  return chunkOverlapper;

}



/************************************************************************/

void PrintChunkOverlapSpec(FILE *fp, ChunkOverlapSpecT *spec){
  fprintf(fp," Overlap " F_CID " " F_CID " %c\n",
	  spec->cidA, spec->cidB, spec->orientation);
}
/************************************************************************/

void InitOverlapSpec(CDS_CID_t cidA, CDS_CID_t cidB,
                     ChunkOrientationType orientation,
                     ChunkOverlapSpecT *spec){
  spec->cidA = cidA;
  spec->cidB = cidB;
  spec->orientation = orientation;
}


/************************************************************************
 * A canonical overlap hs the following properties 
 *       if(orientation is symmetric -- AB_BA or BA_AB )
 *       	then the overlap is stored with the lower numbered chunk first
 *	if(orientation is antinormal -- BA_BA)
 *	        then the overlap is stored as an normal overlap with the chunks reversed
 *************************************************************************/
int InitCanonicalOverlapSpec(CDS_CID_t cidA, CDS_CID_t cidB,
                             ChunkOrientationType orientation,
                             ChunkOverlapSpecT *spec){
  int canonical = TRUE;
  
  switch(orientation){
      
    case BA_BA:
      spec->orientation = AB_AB;
      spec->cidA = cidB;
      spec->cidB = cidA;
      canonical = FALSE;
      break;
    case AB_AB:
      spec->orientation = orientation;
      spec->cidA = cidA;
      spec->cidB = cidB;
      break;
    case BA_AB:
    case AB_BA:
      spec->orientation = orientation;
      if(cidA > cidB){
	spec->cidA = cidB;
	spec->cidB = cidA;
	canonical = FALSE;
      }else{
	spec->cidA = cidA;
	spec->cidB = cidB;
      }
      break;
    default:
      assert(0);
      break;
  }

#if 0
  fprintf(GlobalData->stderrc,"* ChunkOverlap canonical:%d\n\t",
	  canonical);
  PrintChunkOverlapSpec(GlobalData->stderrc, spec);
#endif

  return canonical;
}


/************************************************************************/
/* Given a graph edge, create an overlap in the hashtable */

/* Given a graph edge, create an overlap in the hashtable and mark it as computed */
void CreateChunkOverlapFromEdge(GraphCGW_T *graph,
                                EdgeCGW_T *edge,
                                int bayesian){
  ChunkOverlapCheckT olap;
  double delta = sqrt(edge->distance.variance) * 3.0;
  InitCanonicalOverlapSpec(edge->idA, edge->idB, edge->orient, &olap.spec);
  olap.computed = TRUE;
  olap.overlap = -edge->distance.mean;
  olap.minOverlap = -edge->distance.mean - delta;
  olap.maxOverlap = -edge->distance.mean + delta;;
  olap.fromCGB = FALSE;
  olap.cgbMinOverlap = 0;
  olap.cgbMaxOverlap = 0;
  olap.hasBayesianQuality = bayesian;
  olap.errorRate = CGW_DP_ERATE;
  olap.quality = edge->quality;
  olap.ahg = olap.bhg = 0;
  olap.min_offset = olap.max_offset = 0;
#if 0
  fprintf(GlobalData->stderrc,"* CreateChunkOverlapFromEdge bayes:%d\n",
          olap.hasBayesianQuality);
  PrintChunkOverlapSpec(GlobalData->stderrc, &olap.spec);
#endif
  InsertChunkOverlap(graph->overlapper, &olap);
}

#define DELTA 20

/************************************************************************/
/* Given a graph edge, create an overlap in the hashtable */
void FillChunkOverlapWithEdge(EdgeCGW_T *edge, ChunkOverlapCheckT *olap){
  double delta = sqrt(edge->distance.variance) * 3.0;
  InitCanonicalOverlapSpec(edge->idA, edge->idB, edge->orient, &olap->spec);
  olap->computed = FALSE;
  olap->overlap = -edge->distance.mean;

  // might be unsafe for big variances after tandem mark propagation
  olap->minOverlap = (CDS_COORD_t) -edge->distance.mean - delta;
  olap->maxOverlap = (CDS_COORD_t) -edge->distance.mean + delta;
 
  
  olap->minOverlap = (CDS_COORD_t) -edge->distance.mean - DELTA;
  olap->maxOverlap = (CDS_COORD_t) -edge->distance.mean + DELTA;
  olap->fromCGB = FALSE;
  olap->cgbMinOverlap = 0;
  olap->cgbMaxOverlap = 0;
  olap->hasBayesianQuality = FALSE;
  olap->errorRate = CGW_DP_ERATE;
  olap->quality = edge->quality;
  olap->ahg = olap->bhg = 0;
  olap->min_offset = olap->max_offset = 0;
}


/************************************************************************/
/* Given a graph edge, create an overlap in the hashtable */
void FillChunkOverlapWithUOM(ChunkOverlapCheckT *olap, UnitigOverlapMesg *uom_mesg){
  InitCanonicalOverlapSpec(uom_mesg->chunk1, uom_mesg->chunk2, uom_mesg->orient, &olap->spec);
  olap->computed = FALSE;
  olap->overlap = uom_mesg->best_overlap_length;
  olap->minOverlap = uom_mesg->min_overlap_length;
  olap->maxOverlap = uom_mesg->max_overlap_length;
  olap->hasBayesianQuality = !ScaffoldGraph->alignOverlaps;
  olap->errorRate = CGW_DP_ERATE;
  olap->quality = uom_mesg->quality;
  olap->ahg = olap->bhg = 0;
  olap->min_offset = olap->max_offset = 0;
}


/************************************************************************/
// Returns FALSE if none found
// Returns TRUE and overwrites *olap if found
//
int LookupOverlap(GraphCGW_T *graph,
		  CDS_CID_t cidA, CDS_CID_t cidB,
		  ChunkOrientationType orientation,
		  ChunkOverlapCheckT *olap){
  ChunkOverlapperT *chunkOverlapper = graph->overlapper;
  int isCanonical ;
  ChunkOverlapSpecT spec;
  ChunkOverlapCheckT *lookup;
  isCanonical = InitCanonicalOverlapSpec(cidA, cidB, orientation, &spec);
  lookup = LookupCanonicalOverlap(chunkOverlapper, &spec);
  if(!lookup)  // We didn't find anything
    return FALSE;

  *olap = *lookup;

  if(isCanonical){  // we're done
    return TRUE;
  }
  // Otherwise we have to fix up the retrieved canonical overlap for the non-canonical query
  //
  olap->spec.orientation = orientation;
  olap->spec.cidA = cidA;
  olap->spec.cidB = cidB;
  if(olap->BContainsA || olap->AContainsB)
    {
      int swap;
      NodeCGW_T *a = GetGraphNode(graph, cidA);
      NodeCGW_T *b = GetGraphNode(graph, cidB);

      swap = olap->BContainsA;
      olap->BContainsA = olap->AContainsB;
      olap->AContainsB = swap;
      /// NEW!!!!
      if(olap->AContainsB){

	switch(orientation){
          case AB_AB:
          case AB_BA:
            olap->overlap = b->bpLength.mean + olap->bhg;
            break;
          case BA_AB:
          case BA_BA:
            olap->overlap = b->bpLength.mean - olap->ahg;
            break;
          default:
            assert(0);
            break;
	}
      }else if(olap->BContainsA){
	switch(orientation){
          case AB_AB:
          case AB_BA:
            olap->overlap = a->bpLength.mean - olap->bhg;
            break;
          case BA_AB:
          case BA_BA:
            olap->overlap = a->bpLength.mean + olap->ahg;
            break;
          default:
            assert(0);
            break;
	}
      }
      /// END NEW!
    }

  return TRUE;

}



/************************************************************************/

ChunkOverlapCheckT *LookupCanonicalOverlap(ChunkOverlapperT *chunkOverlapper,
                                           ChunkOverlapSpecT *spec){
  ChunkOverlapCheckT *olap;
  
  olap = (ChunkOverlapCheckT *)LookupInHashTable_AS(chunkOverlapper->hashTable,
                                                    spec, sizeof(*spec));
  
  return olap;
}

/************************************************************************/
int InsertChunkOverlap(ChunkOverlapperT *chunkOverlapper,
                       ChunkOverlapCheckT *olap){
  ChunkOverlapCheckT *nolap =
    GetChunkOverlapCheckTHeapItem(chunkOverlapper->ChunkOverlaps);
  *nolap = *olap;
  assert(nolap->overlap==0||nolap->errorRate >= 0.0);  
  return InsertInHashTable_AS(chunkOverlapper->hashTable, &nolap->spec,
                              sizeof(olap->spec), nolap);

}

/************************************************************************/
int DeleteChunkOverlap(ChunkOverlapperT *chunkOverlapper,
                       ChunkOverlapCheckT *olap){
  return DeleteFromHashTable_AS(chunkOverlapper->hashTable, &olap->spec,
                                sizeof(olap->spec));
}

// use this routine to put an overlap into hash table
int InsertOverlapInHashTable(Overlap *tempOlap,
                             CDS_CID_t cidA, CDS_CID_t cidB,
                             ChunkOrientationType orientation)
{
  ChunkOverlapCheckT olap;
  CDS_CID_t edgeIndex;
  
  olap.spec.cidA = cidA;
  olap.spec.cidB = cidB;
  olap.spec.orientation = orientation;
  
  olap.errorRate = CGW_DP_ERATE;
  olap.computed = 1;
  olap.fromCGB = 0;
  olap.hasBayesianQuality = 0;
  if (tempOlap->begpos < 0 && tempOlap->endpos > 0)  // begpos is ahang, endpos is bhang
    {
      olap.AContainsB = 0;
      olap.BContainsA = 1;
    }
  else if (tempOlap->begpos > 0 && tempOlap->endpos < 0)
    {
      olap.AContainsB = 1;	  
      olap.BContainsA = 0;
    }
  else
    {
      olap.AContainsB = 0;	  
      olap.BContainsA = 0;
    }
  olap.suspicious = 0;

  // olap.overlap is interpreted as the edge length in InsertComputedOverlapEdge
  {
    ChunkInstanceT *contigB = GetGraphNode(ScaffoldGraph->RezGraph, cidB);
    olap.overlap = contigB->bpLength.mean - tempOlap->endpos;	
  }
  
  olap.ahg = tempOlap->begpos;
  olap.bhg = tempOlap->endpos;
  olap.quality = -0.02;
  // olap.min_offset, olap.max_offset not relevant
  
  // if some overlap has already been computed for these chunks, get rid of it
  DeleteChunkOverlap(ScaffoldGraph->RezGraph->overlapper, &olap);

  // now add to hash table
  if(InsertChunkOverlap(ScaffoldGraph->RezGraph->overlapper,
                        &olap) != HASH_SUCCESS)
    assert(0);

  // now add to graph
  edgeIndex = InsertComputedOverlapEdge( ScaffoldGraph->RezGraph, &olap);
  
  fprintf(GlobalData->stderrc, "edgeIndex = " F_CID "\n", edgeIndex);
  
  return 1;
}

/************************************************************************/

void CollectChunkOverlap(GraphCGW_T *graph,
                         CDS_CID_t cidA, CDS_CID_t cidB,
                         ChunkOrientationType orientation, 
                         float   meanOverlap, float   deltaOverlap,
                         float   quality, int bayesian,
                         int fromCGB, 
			 int verbose){
  ChunkOverlapperT *chunkOverlapper = graph->overlapper;
  ChunkOverlapCheckT canOlap, *olap;
  CDS_COORD_t delta;
  CDS_COORD_t minOverlap,maxOverlap;

  delta = (CDS_COORD_t)(3.0 * deltaOverlap);
  delta = MAX(delta, 10);
  minOverlap = MAX(meanOverlap - delta, 0);
  maxOverlap = meanOverlap + delta;
  
  //  assert(maxOverlap - minOverlap >= 20);

  if(maxOverlap < 0){
    // There is no chance that these overlap!
    return;
  }

  // Create a canonical representation of the overlap
  InitCanonicalOverlapSpec(cidA, cidB, orientation, &canOlap.spec);

  if(verbose){
    fprintf(GlobalData->stderrc,"* CollectChunkOverlap! %c (" F_CID "," F_CID ") [" F_COORD "," F_COORD "] %s\n", 
	    canOlap.spec.orientation,
            canOlap.spec.cidA, canOlap.spec.cidB,
            minOverlap, maxOverlap, 
	    (fromCGB?"(Overlap)":""));
    PrintChunkOverlapSpec(GlobalData->stderrc, &canOlap.spec);
  }
  // Lookup to see if we've already stored such an overlap
  olap = LookupCanonicalOverlap(chunkOverlapper, &canOlap.spec);
  if(!olap){
    // If we don't find it, add it.
    if(verbose){
      fprintf(GlobalData->stderrc,"* No Overlap found ...inserting!\n");
      PrintChunkOverlapSpec(GlobalData->stderrc, &canOlap.spec);
    }
    canOlap.computed = FALSE;
    canOlap.overlap = FALSE;
    canOlap.quality = 1.0;
    canOlap.hasBayesianQuality = FALSE;
    canOlap.minOverlap = minOverlap;
    canOlap.maxOverlap = maxOverlap;
    canOlap.fromCGB = fromCGB;
    if(fromCGB && bayesian){
      canOlap.computed = TRUE;
      canOlap.quality = quality;
      canOlap.hasBayesianQuality = TRUE;
      canOlap.overlap = (canOlap.minOverlap + canOlap.maxOverlap)/2;
    }
    canOlap.cgbMinOverlap = minOverlap;
    canOlap.cgbMaxOverlap = maxOverlap;
    canOlap.errorRate = CGW_DP_ERATE;

    // Add it to the symbol table
    if(InsertChunkOverlap(chunkOverlapper, &canOlap) != HASH_SUCCESS)
      assert(0);

  }else{ // we found an overlap
    if(verbose){
      fprintf(GlobalData->stderrc,"* Overlap (" F_CID "," F_CID ") found ...updating! [" F_COORD "," F_COORD "] %s computed:%d\n", 
              olap->spec.cidA, olap->spec.cidB,
              olap->minOverlap, olap->maxOverlap,
              (olap->fromCGB?"(Overlap)":""), olap->computed);
    }
    // We found one.  So now we need to update the existing record so that
    // it reflects the maximal interval we're interested in overlapping.
    if(olap->computed){
      // If we've already computed this one, only recompute if the range is expanded
      if(olap->overlap == 0 &&
	 (minOverlap < olap->minOverlap ||
          maxOverlap > olap->maxOverlap)){
	olap->computed = FALSE; // Recompute
	olap->hasBayesianQuality = FALSE;
	olap->minOverlap = MIN(minOverlap, olap->minOverlap);
	olap->maxOverlap = MAX(maxOverlap, olap->maxOverlap);
      }


    }else{
      if(fromCGB){  // a CGB overlap clobbers whatever is there
        if(!olap->fromCGB){ // Not from the chunk graph builder
          olap->cgbMinOverlap = minOverlap;
          olap->cgbMaxOverlap = maxOverlap;
          olap->overlap = (minOverlap + maxOverlap)/2;
          olap->fromCGB = fromCGB;
        }else{              // From the chunk graph builder
          if(quality < olap->quality){
            olap->cgbMinOverlap = minOverlap;
            olap->cgbMaxOverlap = maxOverlap;
          }
        }	
	olap->quality = quality;
	olap->hasBayesianQuality = bayesian;
	if(bayesian){
	  olap->overlap = (olap->cgbMinOverlap + olap->cgbMaxOverlap)/2;
	  olap->computed = TRUE;
	}

      }
      olap->minOverlap = MIN(minOverlap, olap->minOverlap);
      olap->maxOverlap = MAX(maxOverlap, olap->maxOverlap);
    }
#ifdef DEBUG_CHUNKOVERLAP
    fprintf(GlobalData->stderrc,"* Finished updating! [" F_COORD "," F_COORD "] %s\n", 
	    olap->minOverlap, olap->maxOverlap,
            (olap->fromCGB?"(Overlap)":""));
#endif

  }
}

/************************************************************************/
static VA_TYPE(char) *consensusA = NULL;
static VA_TYPE(char) *consensusB = NULL;
static VA_TYPE(char) *qualityA = NULL;
static VA_TYPE(char) *qualityB = NULL;


void ComputeCanonicalOverlap(GraphCGW_T *graph, 
			     ChunkOverlapCheckT *canOlap){
  CDS_COORD_t lengthA, lengthB;
  CDS_COORD_t beg, end;
  int opposite;
  CDS_COORD_t intended;

  if(consensusA == NULL){
    fprintf(GlobalData->stderrc,"* Allocating static arrays in ComputeCanonicalOverlap \n");
    consensusA = CreateVA_char(2048);
    consensusB = CreateVA_char(2048);
    qualityA = CreateVA_char(2048);
    qualityB = CreateVA_char(2048);
  }


  // Get the consensus sequences for both chunks from the ChunkStore
  lengthA = GetConsensus(graph, canOlap->spec.cidA, consensusA, qualityA);
  lengthB = GetConsensus(graph, canOlap->spec.cidB, consensusB, qualityB);


  prepare_overlap(canOlap,&intended, &beg,&end,&opposite,lengthA, lengthB);

#ifdef DEBUG_CHUNKOVERLAP
  fprintf(GlobalData->stderrc,"* ComputeOverlap " F_CID ", " F_CID " orient = %c min:" F_COORD " max:" F_COORD " beg:" F_COORD " end:" F_COORD " opposite:%d\n",
          canOlap->spec.cidA, canOlap->spec.cidB, canOlap->spec.orientation,
          canOlap->minOverlap, canOlap->maxOverlap,
          beg, end, opposite);
#endif
  // Return value is length of chunk sequence/quality
  // overlap 'em

  {
    InternalFragMesg AFR, BFR;
    OverlapMesg     *O;
    CDS_COORD_t where;
    AFR.sequence   = Getchar(consensusA, 0);
    AFR.quality    = Getchar(qualityA, 0);
    AFR.eaccession = 0;
    AFR.iaccession = canOlap->spec.cidA;
    BFR.sequence   = Getchar(consensusB, 0);
    BFR.quality    = Getchar(qualityB, 0);
    BFR.eaccession = 0;
    BFR.iaccession = canOlap->spec.cidB;

#ifdef DEBUG_CHUNKOVERLAP
    fprintf(GlobalData->stderrc,"* Calling DP_Compare with %s sequences lengths " F_SIZE_T "," F_SIZE_T " " F_COORD "," F_COORD " beg,end =[" F_COORD "," F_COORD "] %c %c\n",
            (!strcmp(AFR.sequence, BFR.sequence)?"EQUAL":"DIFFERENT"),
            strlen(AFR.sequence), strlen(BFR.sequence),
            lengthA, lengthB,
            beg, end, *AFR.sequence, *BFR.sequence);
#endif
    //CGW_DP_ERATE=0.1 changed parameter to canOlap->errRate
    // a new field in ChunkOverlapCheckT
    O = DP_Compare_AS(&AFR,&BFR, 
                      beg, end, opposite,
                      canOlap->errorRate,CGW_DP_THRESH, CGW_DP_MINLEN,
#if 0
                      AS_FIND_OVERLAP ,  // Faster, doesn't compute an alignment, slightly less accurate
#else
                      AS_FIND_ALIGN_NO_TRACE,    // Slower, computes an alignment, throw away delta encoding
#endif
                      &where);



    if(O == NULL){  // Didn't find an overlap
      if(canOlap->fromCGB){

#ifdef DEBUG_CHUNKOVERLAP
        fprintf(GlobalData->stderrc,
                "**** MISSED **** No Overlap %c " F_CID " and " F_CID " looked on (" F_COORD "," F_COORD ") overlapper found (" F_COORD "," F_COORD ")\n",
                canOlap->spec.orientation,
                canOlap->spec.cidA, canOlap->spec.cidB,
                canOlap->minOverlap, canOlap->maxOverlap,
                canOlap->cgbMinOverlap, canOlap->cgbMaxOverlap );
        PrintChunkOverlapSpec(GlobalData->stderrc, &canOlap->spec);
#endif

        DeleteGraphOverlapEdge(graph, canOlap->spec.cidA, canOlap->spec.cidB, canOlap->spec.orientation);

      }
      canOlap->computed = TRUE;
      canOlap->overlap = FALSE;
      canOlap->ahg = canOlap->bhg = 0;
      canOlap->suspicious = FALSE;
    }else{         // Found one....
	    
      adapt_overlap(O,canOlap,lengthA,lengthB,intended);
#ifdef DEBUG_CHUNKOVERLAP
      fprintf(GlobalData->stderrc,
              "**** FOUND ****\n");
#endif
      // did we get the overlap we expected?
      if( ! canOlap->suspicious ){	      
        if(!canOlap->fromCGB){
#if 0
          //#ifdef DEBUG_CHUNKOVERLAP
          fprintf(GlobalData->stderrc,"! %s OVERLAP %c between " F_CID " and " F_CID " overlap:" F_COORD " ahg:" F_COORD " bhg:" F_COORD " quality:%f (beg:" F_COORD " end:" F_COORD " opposite:%d)\n",
                  (!canOlap->fromCGB?"****NEW****":"****CGB****"),
                  canOlap->spec.orientation,
                  canOlap->spec.cidA, canOlap->spec.cidB,
                  canOlap->overlap, 
                  ahg,bhg,
                  O->quality,
                  beg,end,opposite);
#endif
        }
#ifdef DEBUG_CHUNKOVERLAP
        // If we found an overlap that differs from what the overlapper found,
        // flag it
        if(canOlap->fromCGB && 
           (canOlap->overlap < canOlap->cgbMinOverlap - 5  || canOlap->overlap > canOlap->cgbMaxOverlap + 5)){
          ChunkInstanceT *a, *b;
          a = GetGraphNode(graph, canOlap->spec.cidA);
          b = GetGraphNode(graph, canOlap->spec.cidB);
		
          if(a->flags.bits.isUnique && b->flags.bits.isUnique){
            fprintf(GlobalData->stderrc,"*** DIFFERS ** overlap " F_COORD " is outside [" F_COORD "," F_COORD "]\n",
                    canOlap->overlap, canOlap->cgbMinOverlap,
                    canOlap->cgbMaxOverlap);
          }
        }
#endif
        // If direct computation yields a smaller overlap, replace
        // the cgb overlap.  dp_align is more aggressive than the overlapper, otherwise
        // we should always update the mean.  We assume the variance won't change much.
        if(canOlap->fromCGB &&
           canOlap->overlap < (canOlap->cgbMinOverlap + canOlap->cgbMaxOverlap)/2){
          EdgeCGW_T *edge = FindGraphOverlapEdge(graph, canOlap->spec.cidA, canOlap->spec.cidB, canOlap->spec.orientation);
          AssertPtr(edge);
          /*
            PrintGraphEdge(GlobalData->stderrc,graph,
            "\nUpdating Overlap * ", edge, edge->idA);
            fprintf(GlobalData->stderrc,
            "* New value is " F_COORD " (cgbmin:" F_COORD " cgbmax:" F_COORD "\n\n",
            -canOlap->overlap,
            canOlap->cgbMinOverlap,
            canOlap->cgbMaxOverlap);
          */
          edge->distance.mean = -canOlap->overlap;
        }
	      
        //	    Print_Overlap_AS(GlobalData->stderrc,&AFR,&BFR,O);
        canOlap->computed = TRUE;
        canOlap->ahg = O->ahg;
        canOlap->bhg = O->bhg;
        canOlap->quality = O->quality;
        canOlap->min_offset = O->min_offset;
        canOlap->max_offset = O->max_offset;
      }
    }
  }
  
}

void ComputeCanonicalOverlap_new(GraphCGW_T *graph,
                                 ChunkOverlapCheckT *canOlap)
{
  CDS_COORD_t lengthA, lengthB;
  // CDS_COORD_t beg, end;
  // int opposite;
  // CDS_COORD_t intended;
  ChunkOverlapperT *chunkOverlapper = graph->overlapper;
  Overlap * tempOlap1;
  ChunkOverlapSpecT inSpec;


  if(consensusA == NULL)
    {
      fprintf(GlobalData->stderrc,"* Allocating static arrays in ComputeCanonicalOverlap \n");
      consensusA = CreateVA_char(2048);
      consensusB = CreateVA_char(2048);
      qualityA = CreateVA_char(2048);
      qualityB = CreateVA_char(2048);
    }
  
  // Save the input spec 
  inSpec = canOlap->spec;

  canOlap->BContainsA = FALSE;
  canOlap->AContainsB = FALSE;
  canOlap->computed = TRUE;
  canOlap->overlap = FALSE;
  canOlap->ahg = canOlap->bhg = 0;

  if(canOlap->maxOverlap < 0) // no point doing the expensive part if there can be no overlap
    return;


  // Get the consensus sequences for both chunks from the ChunkStore
  lengthA = GetConsensus(graph, canOlap->spec.cidA, consensusA, qualityA);


  lengthB = GetConsensus(graph, canOlap->spec.cidB, consensusB, qualityB);

  if(canOlap->minOverlap > (lengthA+lengthB-CGW_DP_MINLEN)) // no point doing the expensive part if there can be no overlap
    return;

  // prepare_overlap(canOlap, &intended, &beg, &end, &opposite, lengthA, lengthB);

#ifdef DEBUG_CHUNKOVERLAP
  fprintf(GlobalData->stderrc,"* ComputeOverlap " F_CID ", " F_CID " orient = %c min:" F_COORD " max:" F_COORD " beg:" F_COORD " end:" F_COORD " opposite:%d\n",
          canOlap->spec.cidA, canOlap->spec.cidB,
          canOlap->spec.orientation,
          canOlap->minOverlap, canOlap->maxOverlap,
          beg, end, opposite);
#endif
  // Return value is length of chunk sequence/quality
  // overlap 'em  
  {
    char *seq1, *seq2;
    CDS_COORD_t min_ahang, max_ahang;
    seq1   = Getchar(consensusA, 0);
    seq2   = Getchar(consensusB, 0);
	
#ifdef DEBUG_CHUNKOVERLAP
    fprintf(GlobalData->stderrc,"* Calling DP_Compare with %s sequences lengths " F_SIZE_T "," F_SIZE_T " " F_COORD "," F_COORD " beg,end =[" F_COORD "," F_COORD "] %c %c\n",
            (!strcmp(seq1, seq2) ? "EQUAL" : "DIFFERENT"),
            strlen(seq1), strlen(seq2),
            lengthA, lengthB,
            beg, end, *seq1, *seq2);
#endif
	
    min_ahang = lengthA - canOlap->maxOverlap;
    max_ahang = lengthA - canOlap->minOverlap;

    // tempOlap1 is a static down inside of DP_Compare, don't free it
    tempOlap1 =
      OverlapSequences(seq1, seq2, canOlap->spec.orientation,
                       min_ahang, max_ahang,
                       canOlap->errorRate,
                       CGW_DP_THRESH, CGW_DP_MINLEN,
                       (GlobalData->aligner == DP_Compare) ?
                       AS_FIND_ALIGN : AS_FIND_LOCAL_OVERLAP);

    if (tempOlap1 ) {     // Found one....
	  
      // adapt_overlap(O,canOlap,lengthA,lengthB,intended);
	  
      if( tempOlap1->begpos < 0 && tempOlap1->endpos > 0) // ahang is neg and bhang is pos
        canOlap->BContainsA = TRUE;
      else if( tempOlap1->begpos > 0 && tempOlap1->endpos < 0) // ahang is pos and bhang is neg
        canOlap->AContainsB = TRUE;
	  
#ifdef DEBUG_CHUNKOVERLAP
      fprintf(GlobalData->stderrc,
              "**** FOUND ****\n");
#endif

      //	    Print_Overlap_AS(GlobalData->stderrc,&AFR,&BFR,O);
      canOlap->computed = TRUE;
      canOlap->ahg = tempOlap1->begpos;
      canOlap->bhg = tempOlap1->endpos;

      //  Make the overlap field be the number of bases from the tail of
      //  the A sequence to the beginning of the B sequence
      canOlap -> overlap = tempOlap1 -> length;
      if  (canOlap -> ahg < 0)
        canOlap -> overlap -= canOlap -> ahg;
      if  (canOlap -> bhg < 0)
        canOlap -> overlap -= canOlap -> bhg;
	  
      // fields are no longer used in DP_Compare (as opposed to DP_Compare_AS)
      canOlap->quality = 0.0;
      // canOlap->min_offset = O->min_offset;
      // canOlap->max_offset = O->max_offset;

      // here we set the suspicious flag based on what the overlap is
      // if the sequences have slid (e.g., an AB_AB has become a BA_BA)
      // then we change the orientation and set the suspicious flag

      // not dealing with containments here - they can go in either way?

      if (canOlap->ahg < 0 && canOlap->bhg < 0){
#ifndef USE_NEW_DP_COMPARE   // use new computecanonical in pseudo-emulation mode
        // temp hack
        canOlap->computed = TRUE;
        canOlap->overlap = FALSE;
        canOlap->ahg = canOlap->bhg = 0;
#else
        // Try to delete the overlap from the hashtable.  It may or may not
        // be there.  We will (re)insert it later.  If we didn't do this,
        // the hashtable would be corrupted, since the chains in the buckets
        // are ordered by (ida,idb,orientation), so we can't screw around with
        // these, without removing the entry and reinserting it.
        if(HASH_SUCCESS != DeleteChunkOverlap(chunkOverlapper, canOlap)){
#ifdef DEBUG_CHUNKOVERLAP
          fprintf(stderr,"* Couldn't delete (" F_CID "," F_CID ",%c)\n",
                  canOlap->spec.cidA, canOlap->spec.cidB,
                  canOlap->spec.orientation);
#endif
        }
        canOlap->suspicious = TRUE;
        canOlap -> overlap = tempOlap1 -> length;
        switch  (canOlap -> spec . orientation)
          {
            case  AB_AB :
              // we want to switch to a non-canonical orientation
              // ...canOlap->spec.orientation = BA_BA; but since
              // we can't, we switch order of operands instead
              canOlap -> spec. cidA = inSpec . cidB;
              canOlap -> spec. cidB = inSpec . cidA;
              canOlap -> ahg = - canOlap -> ahg;
              canOlap -> bhg = - canOlap -> bhg;
              break;

            case  AB_BA :
              canOlap -> spec. orientation = BA_AB; 
              canOlap -> ahg = - tempOlap1 -> endpos;
              canOlap -> bhg = - tempOlap1 -> begpos;
              break;

            case  BA_AB :
              canOlap -> spec. orientation = AB_BA; 
              canOlap -> ahg = - tempOlap1 -> endpos;
              canOlap -> bhg = - tempOlap1 -> begpos;
              break;

            default :
              fprintf (GlobalData->stderrc, "Non_canonical orientation = %c\n",
                       canOlap -> spec . orientation);
              assert (FALSE);
          }

        fprintf(GlobalData->stderrc,">>> Fixing up suspicious overlap (" F_CID "," F_CID ",%c) (ahg:" F_COORD " bhg:" F_COORD ") to (" F_CID "," F_CID ",%c) (ahg:" F_COORD " bhg:" F_COORD ") len: " F_COORD "\n",
                inSpec.cidA, inSpec.cidB,
                inSpec.orientation,
                tempOlap1->begpos, tempOlap1->endpos,
                canOlap->spec.cidA, canOlap->spec.cidB,
                canOlap->spec.orientation,
                canOlap->ahg, canOlap->bhg,
                canOlap->overlap);

        // Add it to the symbol table

#ifndef DEBUG_CHUNKOVERLAP
        InsertChunkOverlap(chunkOverlapper, canOlap);
#else
        // This insertion may not succeed, if the modified
        // overlap is already in the hashtable
        //
        if(InsertChunkOverlap(chunkOverlapper, canOlap) != HASH_SUCCESS){
          fprintf(GlobalData->stderrc,"* (" F_CID "," F_CID ",%c) already exists...not inserting in hash\n",
                  canOlap->spec.cidB, canOlap->spec.cidA,
                  canOlap->spec.orientation);
        }else{
          fprintf(GlobalData->stderrc,"* (" F_CID "," F_CID ",%c) doesn't exist...inserting in hash\n",
                  canOlap->spec.cidB, canOlap->spec.cidA,
                  canOlap->spec.orientation);
        }
#endif

#endif
      }
    }
  }
}

/* Insert a computed overlap as a CIEdgeT into the Scaffold Graph. */

CDS_CID_t InsertComputedOverlapEdge(GraphCGW_T *graph,
                                    ChunkOverlapCheckT *olap){
  CDS_CID_t eid;
  int fudge;
  int isRepeatOverlap = FALSE; 
  int isDoveTail = FALSE;
  LengthT overlap;
  ChunkOrientationType orient = olap->spec.orientation;
  EdgeCGW_T *existing = FindGraphOverlapEdge(graph, olap->spec.cidA, olap->spec.cidB, orient);
  int verbose = FALSE;

  overlap.mean   = -olap->overlap;
  overlap.variance = MAX(1.0, ComputeFudgeVariance((double)olap->overlap));
  fudge = sqrt(overlap.variance);

  isDoveTail = !(olap->AContainsB || olap->BContainsA);
  if(isDoveTail){
    isRepeatOverlap = IsRepeatOverlap(graph,
                                      olap->spec.cidA, 
                                      olap->spec.cidB,
                                      orient, 
                                      overlap);
    isDoveTail = !isRepeatOverlap;
  }


  // If there is an existing overlap edge, don't insert this one!
  if(existing){
    double diff = abs(existing->distance.mean - overlap.mean);
    if(diff < 5.0){ // this is the same edge
      CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, existing);
      return eid;
    }
  }
  if(verbose)
    fprintf(GlobalData->stderrc,"* Calling AddGraphEdge# (" F_CID "," F_CID ",%c) %g +/- %g (" F_COORD ")\n",
            olap->spec.cidA, olap->spec.cidB,
            orient, overlap.mean,
            overlap.variance, olap->overlap);

  eid = AddGraphEdge(graph, 
                     olap->spec.cidA, 
                     olap->spec.cidB, 
                     NULLINDEX, NULLINDEX, // frags
                     NULLINDEX,  // dist
                     overlap, 
                     olap->quality,
                     fudge,
                     orient, 
                     FALSE, // inducedByUnknownOrientation
                     FALSE, // isGuide
                     FALSE, // isMayJoin
                     FALSE, // isMustJoin
                     isDoveTail,  // isOverlap
                     isRepeatOverlap, // isRepeatOverlap
                     FALSE,           // is tandem
                     olap->AContainsB,                 // isAContainsB
                     olap->BContainsA,                 // isBContainsA
                     FALSE,                        // isTransChunk
                     FALSE, FALSE,  // extremalA and extremalB
                     UNKNOWN_EDGE_STATUS,
                     FALSE,
                     TRUE /* insert*/ );

  {
    EdgeCGW_T *e = GetGraphEdge(graph, eid);

    if(verbose)PrintGraphEdge(GlobalData->stderrc,graph, " Added  ", e, e->idA);

  }
  return eid;
}



static void prepare_overlap(ChunkOverlapCheckT *canOlap,
                            CDS_COORD_t *intended,
                            CDS_COORD_t * beg, CDS_COORD_t * end,
                            int* opposite,
                            CDS_COORD_t lengthA, CDS_COORD_t lengthB)
{

  /* We have to deal with the orientation of the chunks */
  switch(canOlap->spec.orientation){
    case BA_BA:// ANTINORMAL
      // there can be no ANTINORMAL canonical overlap
      assert(0);
      break;
    case AB_AB: // NORMAL
      /* We use the chunks as retrieved from the chunk store */
      // The normal canonical overlaps that can occur are (first frag is A, second frag is B)
      // 1)   =============>
      //      +++++=============>
      // 
      // DP_Compare should return the orientation NORMAL, and NOT flip the fragments

      // 2)   ==================>
      //      ++++========>
      //
      // DP_Compare should return the orientation NORMAL, and NOT flip the fragments
    

      // 3)   -----=====>
      //      ================>
      //
      // DP_Compare should return the orientation NORMAL, and DO flip the fragments
     
      // we do not complement the second fragment
      // and set beg,end to the interval around the marked  region 
      // (+++ means positiv, --- negativ)
      // we need the intended to check whether DP_Compare computed
      // the overlap we intended to find.
      *beg = lengthA - canOlap->maxOverlap;
      *end = lengthA - canOlap->minOverlap;
      *intended = lengthA - canOlap->overlap;
      *opposite = FALSE;
      break;
    case AB_BA: // INNIE
      /* We need to reverse chunk B */
      // The innie overlaps that can occur are (first frag is A, second frag is B)
      // 1)   =============>
      //      +++++<===========
      // 
      // DP_Compare should return the orientation INNIE, and NOT flip the fragments
   
      // 2)   ================>
      //      +++++<======
      //
      // DP_Compare should return the orientation INNIE, and NOT flip the fragments

      // 3)   ------=======>
      //      <===============
      //
      // DP_Compare should return the orientation OUTTIE , and DO flip the fragments
    
      // we DO complement the second fragment
      // and set beg,end to the interval around the marked  region 
      // (+++ means positiv, --- negativ)
      *beg = lengthA - canOlap->maxOverlap;
      *end = lengthA - canOlap->minOverlap;
      *intended = lengthA - canOlap->overlap;
      *opposite = TRUE;
      break;
    case BA_AB: // OUTTIE
      // The outtie overlaps that can occur are (first frag is A, second frag is B)
      // since DP_Compare cannot handle ANTINORMALS we compute a different
      // beg, end which means that means we really ask DP_Compare for an INNIE
      //
      // 1)   <=============---
      //           ===========>
      // 
      // DP_Compare should return the orientation OUTTIE, and flip the fragments
   
      // 2)   <================
      //           ======>+++++
      //
      // DP_Compare should return the orientation INNIE, and NOT flip the fragments

      // 3)       <=======----
      //      ===============>
      //
      // DP_Compare should return the orientation OUTTIE , and flip the fragments
    
      // we DO complement the second fragment
      // and set beg,end to the interval around the marked  region 
      // (+++ means positiv, --- negativ)
      //        *beg = lengthA - canOlap->maxOverlap;
      //        *end = lengthA - canOlap->minOverlap;
      //        *intended = lengthA - canOlap->overlap;
      *beg = -(lengthB - canOlap->minOverlap);
      *end = -(lengthB - canOlap->maxOverlap);
      *intended = -(lengthB - canOlap->overlap);
      *opposite = TRUE;
      break;
    default:
      assert(0);
  }

#ifdef DEBUG    
  fprintf(GlobalData->stderrc,"In prepare overlap. beg = " F_COORD ", end = " F_COORD "\n",*beg,*end);
#endif
  
  // Add some slop
  *beg -= BEGENDSLOP;
  *end += BEGENDSLOP;  

}
/* 
   make sense of the results of DP_Compare.
   - Compute the overlap length
   - set containment flags
   - mark an overlap as suspicious if we got not back what we exspected
*/

#define SLOP 5

static int adapt_overlap(OverlapMesg* O, ChunkOverlapCheckT* canOlap,
                         CDS_COORD_t lengthA, CDS_COORD_t lengthB,
                         CDS_COORD_t intended)
{
  // did DP_Compare flip the frgment IDs ?
  int flipped  = (O->aifrag != canOlap->spec.cidA);
  // did DP_Compare flip the orientation ?
  int oflipped = (O->orientation != canOlap->spec.orientation);

  CDS_COORD_t overlap;
  CDS_COORD_t ahg = O->ahg;
  CDS_COORD_t bhg = O->bhg;
  
  assert( ahg >= 0 );

  canOlap->BContainsA = FALSE;
  canOlap->AContainsB = FALSE;

  // for OUTTIES it should be ((lengthA - bhg) + ( lengthB - ahg))/2;
  // which is the same as below
  overlap= ((lengthA - ahg) + ( lengthB - bhg))/2;

  if( bhg < 0) 
    {
      if( flipped )
	canOlap->BContainsA = TRUE;
      else
	canOlap->AContainsB = TRUE;
    }  
  
  canOlap->suspicious = FALSE;
#if GREEDYDEBUG > 1
  fprintf(GlobalData->stderrc,"In adapt overlap. Computed overlap = " F_COORD "\n",overlap);
#endif

  /* In the cas of an OUTTIE, if DP_Compare flips the two fragment IDs it 
     should not have a different orientation
     In the case of an INNIE if DP_Compare does flip the fragment IDs
     then it also should flip the orientation */
  if( canOlap->spec.orientation == BA_AB ) // OUTTIE
    canOlap->suspicious = (flipped == oflipped) ;
  else
    if( canOlap->spec.orientation == AB_BA ) // INNIE
      canOlap->suspicious = (flipped != oflipped) ;
    else
      if( oflipped ) // in the NORMAL case the orientation should not be flipped
        canOlap->suspicious = TRUE;
  
  /* if our intended beg was negativ, the computed should also be 
     and vice versa */
  if( intended < -SLOP && ! flipped )
    canOlap->suspicious = TRUE; 
  
  if( intended > SLOP && flipped)
    canOlap->suspicious = TRUE; 
  
#if 0
  // if the overlap is suspicious we adapt its length to 0
  // since a lot  of code assumes that as an indication of
  // a failed overlap
  if( canOlap->suspicious )
    canOlap->overlap = 0;
  else
    canOlap->overlap = overlap;
#else
  canOlap->overlap = overlap;
#endif

  return flipped;
  
}




/************************************************************************/
OverlapMesg *ComputeCanonicalOverlapWithTrace(GraphCGW_T *graph, 
					      ChunkOverlapCheckT *canOlap,
					      InternalFragMesg *AFR,
                                              InternalFragMesg *BFR,
					      FILE* fp, int iterate){
  CDS_COORD_t lengthA, lengthB;
  CDS_COORD_t beg, end;
  int opposite;
  CDS_COORD_t intended;
  OverlapMesg     *O;

  if(consensusA == NULL){
    consensusA = CreateVA_char(2048);
    consensusB = CreateVA_char(2048);
    qualityA = CreateVA_char(2048);
    qualityB = CreateVA_char(2048);
  }

  if( fp == NULL )
    fp = GlobalData->stderrc;

  // Get the consensus sequences for both chunks from the ChunkStore
  lengthA = GetConsensus(graph, canOlap->spec.cidA, consensusA, qualityA);
  lengthB = GetConsensus(graph, canOlap->spec.cidB, consensusB, qualityB);

  prepare_overlap(canOlap,&intended,&beg,&end,&opposite,lengthA, lengthB);
  
  // Return value is length of chunk sequence/quality
  // **consensus will point to consensus data
  // **quality   will point to quality data
  // overlap 'em
  
  {
    CDS_COORD_t where;

    AFR->sequence   =  Getchar(consensusA, 0);
    AFR->quality    =  Getchar(qualityA, 0);
    AFR->eaccession = 0;
    AFR->iaccession = canOlap->spec.cidA;
    BFR->sequence   =  Getchar(consensusB, 0);
    BFR->quality    =  Getchar(qualityB, 0);
    BFR->eaccession = 0;
    BFR->iaccession = canOlap->spec.cidB;
    

#if GREEDYDEBUG > 1
    fprintf(fp,"* Calling DP_Compare with trace with %s sequences lengths " F_SIZE_T "," F_SIZE_T " " F_COORD "," F_COORD " beg,end =[" F_COORD "," F_COORD "]\n",
	    (!strcmp(AFR->sequence, BFR->sequence)?"EQUAL":"DIFFERENT"),
	    strlen(AFR->sequence), strlen(BFR->sequence),
            lengthA, lengthB,
	    beg, end);
#endif
    //CGW_DP_ERATE=0.1 changed parameter to canOlap->errRate
    // a new field in ChunkOverlapCheckT

    if( ! iterate )
      O = DP_Compare_AS(AFR,BFR, 
			beg, end, opposite,
			canOlap->errorRate,CGW_DP_THRESH, CGW_DP_MINLEN,
			AS_FIND_ALIGN,
			&where);
    else
      {
	CDS_COORD_t minBeg = beg;
	CDS_COORD_t maxEnd = end;
	int goOn = TRUE;

	beg = intended-BEGENDSLOP;
	end = intended+BEGENDSLOP;
	
	O = NULL;
	while( O == NULL && goOn )
	  {
	    O = DP_Compare_AS(AFR,BFR, 
			      beg, end, opposite,
			      canOlap->errorRate,CGW_DP_THRESH, CGW_DP_MINLEN,
			      AS_FIND_ALIGN,
			      &where);
#if GREEDYDEBUG > 0
	    if (O != NULL)
	      fprintf(GlobalData->stderrc,"* TANDEM REPEAT succeded with beg=" F_COORD " , end=" F_COORD "\n",beg,end);
#endif
	    
	    if( (beg == minBeg) && (end == maxEnd) )
	      goOn = FALSE;
	    beg = MAX(minBeg,beg-BEGENDSLOP);
	    end = MIN(maxEnd,end+BEGENDSLOP);
	  }
	
      }
    
#if GREEDYDEBUG > 0
    fprintf(GlobalData->stderrc,"* queried overlap " F_CID ", " F_CID " orient = %c min:" F_COORD " max:" F_COORD " beg:" F_COORD " end:" F_COORD " intended:" F_COORD " opposite:%d\n",
	    canOlap->spec.cidA, canOlap->spec.cidB, canOlap->spec.orientation,
	    canOlap->minOverlap, canOlap->maxOverlap,
	    beg, end, intended, opposite);
    if( canOlap->minOverlap != canOlap->maxOverlap )
      fprintf(GlobalData->stderrc,"* MIN/MAX different min=" F_COORD ", max = " F_COORD "\n",
              canOlap->minOverlap,canOlap->maxOverlap);
    if( O != NULL )
      fprintf(GlobalData->stderrc,"$ computed overlap (DP_Compare) idA = " F_CID " , idB = " F_CID ", orientation = %c, ahg = " F_COORD ", bhg = " F_COORD "\n",
              O->aifrag,O->bifrag,O->orientation,O->ahg,O->bhg);
#endif


    if( O != NULL)
      { 
#if GREEDYDEBUG > 0
	int flipped = adapt_overlap(O,canOlap,lengthA,lengthB, intended);
        
	fprintf(GlobalData->stderrc,"* queried overlap after adaption " F_CID ", " F_CID " orient = %c min:" F_COORD " max:" F_COORD " overlap:" F_COORD " beg:" F_COORD " end:" F_COORD " opposite:%d\n",
		canOlap->spec.cidA, canOlap->spec.cidB,
                canOlap->spec.orientation,
		canOlap->minOverlap, canOlap->maxOverlap, canOlap->overlap,
		beg, end, opposite);

	if ( flipped )
	  fprintf(fp,"FLIPPED\n");
#endif
	// handle containment overlaps correctly

#if GREEDYDEBUG > 4
	Print_Overlap_AS(fp,AFR,BFR,O);
#endif     
      }

    return O;
    
  }
}





/************************************************************************/
// Returns FALSE if none found
// Returns TRUE  if found and sets quality bit if necessary
 
int LookupQualityOverlap(GraphCGW_T *graph, 
                         EdgeCGW_T *edge,
                         ChunkOrientationType orientation,
                         ChunkOverlapCheckT *olap, QualityFuncT qfunc,
                         float* quality, FILE* fp){
  int isCanonical;
  ChunkOverlapSpecT spec;
  ChunkOverlapCheckT *lookup;
  ChunkOverlapperT *chunkOverlapper = graph->overlapper;
  OverlapMesg* omesg = NULL;

  CDS_CID_t cidA = edge->idA;
  CDS_CID_t cidB = edge->idB;

  if(fp == NULL)
    fp = GlobalData->stderrc;

  isCanonical = InitCanonicalOverlapSpec(cidA, cidB, orientation, &spec);
  lookup = LookupCanonicalOverlap(chunkOverlapper, &spec);
  if(!lookup)  // We didn't find anything
    return FALSE;

  // There was sth in the table but the overlap length is 0
  if( lookup->overlap == 0 )
    return FALSE;

  // Otherwise we found something
  // First we check whether the this overlap was already computed for the
  // given quality function
  // Do this if-statement for each quality-bit
  // in ChunkOverlapCheckT

  *olap = *lookup;

  if( lookup->hasBayesianQuality != FALSE)
    {
#if GREEDYDEBUG > 1
      fprintf(fp,"BAYESIAN flag set\n");
#endif
      *quality = lookup->quality;
      return TRUE;
    }

  // if we are here we have to recompute the alignment 
  // and invoke the appropriate
  // quality function
  switch( qfunc ){
    case BAYESIAN :
      {
        InternalFragMesg IFG1,IFG2;
        static OverlapMesg omesgBuffer;
        CDS_COORD_t length;
        float normalQuality;
        // float QVQuality;
        // Compute the DP_Compare alignment and fill the IFG1 and IFG2 with the overlapping
        // pieces of the two chunks (resp. the quality values
        // if the edge has the tandem bit set we try different error rates
        if( edge->flags.bits.hasTandemOverlap)
          {
            omesg = ComputeCanonicalOverlapWithTrace(graph, lookup, &IFG1, &IFG2, fp, TRUE);
          }
        else
          {
            omesg = ComputeCanonicalOverlapWithTrace(graph, lookup, &IFG1, &IFG2, fp, FALSE);
          }

        if( omesg == NULL )
          return FALSE;

        if( omesg->delta == NULL )
          {
            fprintf(GlobalData->stderrc, "\n\n\n***** tripped the hack for the case of omesg->delta == NULL\n\n\n");
            return FALSE;
          }

        omesgBuffer.aifrag       = omesg->aifrag;
        omesgBuffer.bifrag       = omesg->bifrag;
        omesgBuffer.ahg          = omesg->ahg;
        omesgBuffer.bhg          = omesg->bhg;
        omesgBuffer.orientation  = omesg->orientation;
        omesgBuffer.overlap_type = omesg->overlap_type;
        omesgBuffer.quality      = omesg->quality;
        omesgBuffer.min_offset   = omesg->min_offset;
        omesgBuffer.max_offset   = omesg->max_offset;
        omesgBuffer.polymorph_ct = omesg->polymorph_ct;

        omesgBuffer.delta = (signed char*) safe_calloc(sizeof(signed char),strlen((char*)omesg->delta)+1);
        strcpy((char*)omesgBuffer.delta,(char*)omesg->delta);
        // compute the quality value and
        // set the quality and the bit in the ChunkOverlapCheckT
        // we do this with and without quality realigning
        compute_bayesian_quality(&IFG1,&IFG2,&omesgBuffer,0,&length,fp);

        // Realign it using quality values
        normalQuality = omesgBuffer.quality;
#if GREEDYDEBUG > 1
        fprintf(fp,"Quality between " F_CID " and " F_CID " WITHOUT QV Realigning = %f\n",
                cidA,cidB,normalQuality);
#endif   

        // Postponed. First integrate greedy walking
        // Then see whetehr QV Realigning makes it better.
        /*
          omesg = QV_ReAligner_AS(&IFG1,&IFG2,&omesgBuffer);
          compute_bayesian_quality(&IFG1,&IFG2,omesg,0,&length,fp);
          QVQuality = omesg->quality;
          fprintf(fp,"Quality WITH QV Realigning = %f\n",QVQuality);
          if( QVQuality != normalQuality)
	  fprintf(fp,"CHANGED quality values\n");
        */

        *quality = omesgBuffer.quality;
        lookup->quality = *quality;
        lookup->hasBayesianQuality = TRUE;
        *olap = *lookup;
      }
  }
  

  // We check whether DP_compare has changed cidA and cidB
  // with respect to the original cidA cidB. Hence we have to
  // pay attention, whether Saul has switched the ids to make
  // them canonical
  // If so, we allow only a sloppy positiv ahang
  // otherwise DP_Compare has computed a different
  // overlap. Hence we return FALSE
  
  if(  lookup->suspicious  )
    {
#if GREEDYDEBUG > -1
      fprintf(fp,"SUSPICIOUS OVERLAP omesg cidA = " F_CID " , canOlap " F_CID " cidA = " F_CID ", omesg Orientation = %c\n",
              omesg->aifrag,olap->spec.cidA,cidA,omesg->orientation);
#endif
      return FALSE;
    }
  else
    {
#if GREEDYDEBUG > -1
      fprintf(fp,"OVERLAP omesg cidA = " F_CID " , canOlap " F_CID " cidA = " F_CID ", omesg Orientation = %c\n",
              omesg->aifrag,olap->spec.cidA,cidA,omesg->orientation);
#endif
    }

  if(isCanonical){  // we're done
    return TRUE;
  }
  // Otherwise we have to fix up the retrieved canonical overlap for the non-canonical query
  //
  olap->spec.orientation = orientation;
  olap->spec.cidA = cidA;
  olap->spec.cidB = cidB;
  {
    int swap;
    swap = olap->BContainsA;
    olap->BContainsA = olap->AContainsB;
    olap->AContainsB = swap;
  }
    
  return TRUE;

}




/************************************************************************/
// Returns FALSE if none found
// Returns TRUE  if found and sets quality bit if necessary

int ComputeQualityOverlap(GraphCGW_T *graph, 
			  EdgeCGW_T *edge,
			  ChunkOrientationType orientation,
			  ChunkOverlapCheckT *olap, QualityFuncT qfunc,
			  float* quality, FILE* fp){
  int isCanonical;
  ChunkOverlapSpecT spec;
  ChunkOverlapCheckT lookup;
  OverlapMesg* omesg = NULL;

  CDS_CID_t cidA = edge->idA;
  CDS_CID_t cidB = edge->idB;


  if(fp == NULL)
    fp = GlobalData->stderrc;

  isCanonical = InitCanonicalOverlapSpec(cidA, cidB, orientation, &spec);

  FillChunkOverlapWithEdge(edge,&lookup);
  *olap = lookup;

  
  // if we are here we have to recompute the alignment 
  // and invoke the appropriate
  // quality function
  switch( qfunc ){
    case BAYESIAN :
      {
        InternalFragMesg IFG1,IFG2;
        static OverlapMesg omesgBuffer;
        CDS_COORD_t length;
        float normalQuality;
        // float QVQuality;
        // Compute the DP_Compare alignment and fill the IFG1 and IFG2 with the overlapping
        // pieces of the two chunks (resp. the quality values
        // if the edge has the tandem bit set we try different error rates
        if( edge->flags.bits.hasTandemOverlap)
          {
            omesg = ComputeCanonicalOverlapWithTrace(graph, &lookup, &IFG1, &IFG2, fp, TRUE);
          }
        else
          {
            omesg = ComputeCanonicalOverlapWithTrace(graph, &lookup, &IFG1, &IFG2, fp, FALSE);
          }

        if( omesg == NULL )
          return FALSE;

        omesgBuffer.aifrag       = omesg->aifrag;
        omesgBuffer.bifrag       = omesg->bifrag;
        omesgBuffer.ahg          = omesg->ahg;
        omesgBuffer.bhg          = omesg->bhg;
        omesgBuffer.orientation  = omesg->orientation;
        omesgBuffer.overlap_type = omesg->overlap_type;
        omesgBuffer.quality      = omesg->quality;
        omesgBuffer.min_offset   = omesg->min_offset;
        omesgBuffer.max_offset   = omesg->max_offset;
        omesgBuffer.polymorph_ct = omesg->polymorph_ct;

        omesgBuffer.delta = (signed char*) safe_calloc(sizeof(signed char),strlen((char*)omesg->delta)+1);
        strcpy((char*)omesgBuffer.delta,(char*)omesg->delta);
        // compute the quality value and
        // set the quality and the bit in the ChunkOverlapCheckT
        // we do this with and without quality realigning
        compute_bayesian_quality(&IFG1,&IFG2,&omesgBuffer,0,&length,fp);
        // Realign it using quality values
        normalQuality = omesgBuffer.quality;
#if GREEDYDEBUG > 1
        fprintf(fp,"Quality between " F_CID " and " F_CID " WITHOUT QV Realigning = %f\n",
                cidA,cidB,normalQuality);
#endif   

        // Postponed. First integrate greedy walking
        // Then see whether QV Realigning makes it better.
        /*
          omesg = QV_ReAligner_AS(&IFG1,&IFG2,&omesgBuffer);
          compute_bayesian_quality(&IFG1,&IFG2,omesg,0,&length,fp);
          QVQuality = omesg->quality;
          fprintf(fp,"Quality WITH QV Realigning = %f\n",QVQuality);
          if( QVQuality != normalQuality)
	  fprintf(fp,"CHANGED quality values\n");
        */

        *quality = omesgBuffer.quality;
        lookup.quality = *quality;
        lookup.hasBayesianQuality = TRUE;
        *olap = lookup;
      }
  }
  
  // We check whether DP_compare has changed cidA and cidB
  // with respect to the original cidA cidB. Hence we have to
  // pay attention, whether Saul has switched the ids to make
  // them canonical
  // If so, we allow only a sloppy positiv ahang
  // otherwise DP_Compare has computed a different
  // overlap. Hence we return FALSE


  
  if(  lookup.suspicious )
    {
      fprintf(fp,"SUSPICIOUS OVERLAP omesg cidA = " F_CID " , canOlap " F_CID " cidA = " F_CID ", omesg Orientation = %c\n",
              omesg->aifrag,olap->spec.cidA,cidA,omesg->orientation);
      return FALSE;
    }
  else
    fprintf(fp,"OVERLAP omesg cidA = " F_CID " , canOlap " F_CID " cidA = " F_CID ", omesg Orientation = %c\n",
            omesg->aifrag,olap->spec.cidA,cidA,omesg->orientation);


  if(isCanonical){  // we're done
    return TRUE;
  }
  // Otherwise we have to fix up the retrieved canonical overlap for the non-canonical query
  //

  olap->spec.orientation = orientation;
  olap->spec.cidA = cidA;
  olap->spec.cidB = cidB;
  {
    int swap;
    swap = olap->BContainsA;
    olap->BContainsA = olap->AContainsB;
    olap->AContainsB = swap;
  }

  return TRUE;

}





/************************************************************************/
// Returns FALSE if none found
// Returns TRUE  if found and sets quality bit if necessary

int ComputeUOMQualityOverlap(GraphCGW_T *graph, 
			     UnitigOverlapMesg *uom_mesg,
			     ChunkOverlapCheckT *olap,
			     float* quality){
  int isCanonical;
  ChunkOverlapSpecT spec;
  ChunkOverlapCheckT lookup;

  CDS_CID_t cidA = uom_mesg->chunk1;
  CDS_CID_t cidB = uom_mesg->chunk2;
  OverlapMesg* omesg;

#if GREEDYDEBUG > 1
  fprintf(GlobalData->stderrc,"\nENTERING ComputeUOMQualityOverlap to compute overlap between cidA = " F_CID " and cidB = " F_CID " orient = %c\n",
          cidA,cidB,uom_mesg->orient);
#endif


  isCanonical = InitCanonicalOverlapSpec(cidA, cidB, uom_mesg->orient, &spec);
 
  FillChunkOverlapWithUOM(&lookup,uom_mesg);
  *olap = lookup;

  // if we are here we have to recompute the alignment 
  // and invoke the appropriate
  // quality function
  {
    InternalFragMesg IFG1,IFG2;
    static OverlapMesg omesgBuffer;
    CDS_COORD_t length;
    float normalQuality;
    // Compute the DP_Compare alignment and fill the IFG1 and IFG2 with the overlapping
    // pieces of the two chunks (resp. the quality values
    // if the edge has the tandem bit set we try different error rates

    if( uom_mesg->min_overlap_length != uom_mesg->max_overlap_length)
      {
	omesg = ComputeCanonicalOverlapWithTrace(graph, &lookup, &IFG1, &IFG2, NULL, TRUE);
      }
    else
      {
	omesg = ComputeCanonicalOverlapWithTrace(graph, &lookup, &IFG1, &IFG2, NULL, FALSE);
      }
    
    if( omesg == NULL )
      {
#if GREEDYDEBUG > 1
	fprintf(GlobalData->stderrc,"FAILED to find overlap for type %c\n",
                uom_mesg->overlap_type);
#endif
	return FALSE;
      }    
    
    omesgBuffer.aifrag       = omesg->aifrag;
    omesgBuffer.bifrag       = omesg->bifrag;
    omesgBuffer.ahg          = omesg->ahg;
    omesgBuffer.bhg          = omesg->bhg;
    omesgBuffer.orientation  = omesg->orientation;
    omesgBuffer.overlap_type = omesg->overlap_type;
    omesgBuffer.quality      = omesg->quality;
    omesgBuffer.min_offset   = omesg->min_offset;
    omesgBuffer.max_offset   = omesg->max_offset;
    omesgBuffer.polymorph_ct = omesg->polymorph_ct;
    
    omesgBuffer.delta = (signed char*) safe_calloc(sizeof(signed char),strlen((char*)omesg->delta)+1);
    strcpy((char*)omesgBuffer.delta,(char*)omesg->delta);
    // compute the quality value and
    // set the quality and the bit in the ChunkOverlapCheckT
    // we do this with and without quality realigning

    compute_bayesian_quality(&IFG1,&IFG2,&omesgBuffer,0,&length,NULL);
    normalQuality = omesgBuffer.quality;

    *quality = omesgBuffer.quality;
    lookup.quality = *quality;
    lookup.hasBayesianQuality = TRUE;
    *olap = lookup;
  }

  if(  lookup.suspicious  )
    {
#if GREEDYDEBUG > 1
      fprintf(GlobalData->stderrc,"SUSPICIOUS OVERLAP omesg cidA = " F_CID " , canOlap " F_CID " cidA = " F_CID ", omesg Orientation = %c\n",
              omesg->aifrag,olap->spec.cidA,cidA,omesg->orientation);
#endif
      return FALSE;
    }
  else
    {
#if GREEDYDEBUG > 1
      fprintf(GlobalData->stderrc,"OVERLAP omesg cidA = " F_CID " , canOlap " F_CID " cidA = " F_CID ", omesg Orientation = %c\n",
              omesg->aifrag,olap->spec.cidA,cidA,omesg->orientation);
#endif
    }
  

  if(isCanonical){  // we're done
    return TRUE;
  }
  // Otherwise we have to fix up the retrieved canonical overlap for the non-canonical query
  //
  olap->spec.orientation = uom_mesg->orient;
  olap->spec.cidA = cidA;
  olap->spec.cidB = cidB;
  {
    int swap;
    swap = olap->BContainsA;
    olap->BContainsA = olap->AContainsB;
    olap->AContainsB = swap;
  }
  return TRUE;
  
}







/************************************************************************/
// This is the top level routine that computes all new potential overlaps.
//
void DumpOverlaps(GraphCGW_T *graph){
  HashTable_Iterator_AS iterator;
  void *key, *value;

  StartTimerT(&GlobalData->OverlapTimer);

  fprintf(GlobalData->stderrc,"* DumpOverlaps ************\n");

  // Iterate over all hashtable elements, computing overlaps

  InitializeHashTable_Iterator_AS(graph->overlapper->hashTable, &iterator);

  while(NextHashTable_Iterator_AS(&iterator, &key, &value)){
    ChunkOverlapCheckT *olap = (ChunkOverlapCheckT*) value;    
    
    fprintf(GlobalData->stderrc,"* (" F_CID "," F_CID ",%c) olap:" F_COORD " ahg:" F_COORD " bhg:" F_COORD "  fromcgb:%d\n",
            olap->spec.cidA, olap->spec.cidB,
            olap->spec.orientation, olap->overlap,
	    olap->ahg, olap->bhg, olap->fromCGB);
  }
}




/************************************************************************/
// This is the top level routine that computes all new potential overlaps.
//
#define NUM_SECTIONS (5)

void ComputeOverlaps(GraphCGW_T *graph, int addEdgeMates,
                     int recomputeCGBOverlaps)
{
  int i = 0;
  HashTable_Iterator_AS iterator;
  void *key, *value;
  int sectionOuter, sectionOuterMin, sectionOuterMax;
  int sectionInner, sectionInnerMin, sectionInnerMax;
  int numOverlaps = 0;

  StartTimerT(&GlobalData->OverlapTimer);

  fprintf(GlobalData->stderrc,"* ComputeOverlaps ************\n");

  // Iterate over all hashtable elements, computing overlaps
  for ( sectionOuter = 0; sectionOuter < NUM_SECTIONS; sectionOuter++)
    {
      sectionOuterMin = sectionOuter * (GetNumGraphNodes(graph)) / NUM_SECTIONS;
      sectionOuterMax = (sectionOuter + 1) * (GetNumGraphNodes(graph)) / NUM_SECTIONS;

      for ( sectionInner = 0; sectionInner < NUM_SECTIONS; sectionInner++)
	{
	  sectionInnerMin = sectionInner * (GetNumGraphNodes(graph)) / NUM_SECTIONS;
	  sectionInnerMax = (sectionInner + 1) * (GetNumGraphNodes(graph)) / NUM_SECTIONS;
	  
	  fprintf(GlobalData->stderrc,"ComputeOverlaps section (o %d,i %d) outer:[%d,%d) inner:[%d,%d)\n",
                  sectionOuter,  sectionInner,
                  sectionOuterMin, sectionOuterMax,
                  sectionInnerMin, sectionInnerMax);
	  fflush(GlobalData->stderrc);
	  
	  InitializeHashTable_Iterator_AS(graph->overlapper->hashTable, &iterator);

	  while(NextHashTable_Iterator_AS(&iterator, &key, &value))
            {
              ChunkOverlapCheckT *olap = (ChunkOverlapCheckT*) value;    
		
#if 0
              fprintf(GlobalData->stderrc,"* olap is (" F_CID "," F_CID ",%c) bayes:%d computed:%d\n",
                      olap->spec.cidA, 
                      olap->spec.cidB, 
                      olap->spec.orientation,
                      olap->hasBayesianQuality,
                      olap->computed);
#endif
		
              assert(key == value);
		
              {
                int inSection = FALSE;
                CDS_CID_t smaller, bigger;
		  
                smaller = MIN( olap->spec.cidA, olap->spec.cidB);
                bigger = MAX( olap->spec.cidA, olap->spec.cidB);
		  
                inSection = (smaller < sectionOuterMax && smaller >= sectionOuterMin) &&
                  (bigger < sectionInnerMax && bigger >= sectionInnerMin);
		  
                // Only do the overlaps where the larger of the ids is within range.
                // The overlaps order of appearance is sorted by the smaller of the indices.
                if(olap->computed || !inSection)
                  continue;
		  
#if 0
                fprintf(stderr,"* (" F_CID "," F_CID ") \n",
                        smaller, bigger);
#endif
              }
		
		
              if(!olap->computed &&  /* We haven't computed this one, and it isn't from cgb, or recomputeCGBOverlaps is set */
                 (!olap->fromCGB || recomputeCGBOverlaps)){
                ChunkOrientationType orientation = olap->spec.orientation;
		  
                // set errRate to old value
                olap->errorRate = CGW_DP_ERATE;
		  
                // first we trust that overlap
                olap->suspicious = FALSE;
		  
                if(olap->maxOverlap < 0){ // Dummy!  Who put this overlap in the table?  No overlap is possible.....SAK
                  fprintf(stderr,"* maxOverlap < 0\n");
                  olap->overlap = 0;
                  continue;
                }
                if((++i % 100000) == 0){
                  fprintf(GlobalData->stderrc,
                          "* ComputeOverlaps %d  (" F_CID "," F_CID ",%c)\n",
                          i, olap->spec.cidA, olap->spec.cidB,
                          olap->spec.orientation);
                }		

                numOverlaps++;

#ifdef USE_NEW_DP_COMPARE
                {
		  ChunkOverlapSpecT inSpec;
		  inSpec = olap->spec;
		  ComputeCanonicalOverlap_new(graph, olap);
		  if(olap->suspicious)
                    {
                      int lengthA, lengthB;
                      lengthA = GetConsensus(graph, olap->spec.cidA, consensusA, qualityA);
                      lengthB = GetConsensus(graph, olap->spec.cidB, consensusB, qualityB);

                      fprintf(GlobalData->stderrc,"* CO: SUSPICIOUS Overlap found! Looked for (" F_CID "," F_CID ",%c)[" F_COORD "," F_COORD "]"
                              "found (" F_CID "," F_CID ",%c) " F_COORD "; contig lengths as found (%d,%d)\n",
                              inSpec.cidA, inSpec.cidB, orientation, olap->minOverlap, olap->maxOverlap,
                              olap->spec.cidA, olap->spec.cidB, olap->spec.orientation, olap->overlap,
                              lengthA,lengthB);
                    }
                }
#else
                ComputeCanonicalOverlap(graph, olap);
#endif
		  
#ifdef DEBUG_CHUNKOVERLAP
                fprintf(GlobalData->stderrc,"* addEdgeMates = %d fromCGB = %d (" F_CID "," F_CID ",%c) olap:" F_COORD "\n",
                        addEdgeMates, olap->fromCGB,
                        olap->spec.cidA, olap->spec.cidB,
                        olap->spec.orientation, olap->overlap);
#endif
                if(addEdgeMates && !olap->fromCGB && olap->overlap)
		  {
                    //fprintf(GlobalData->stderrc,"* Calling AddGraphEdge$\n");
                    InsertComputedOverlapEdge(graph, olap);
		  }	  
              }
            }
	}
    }
  
  StopTimerT(&GlobalData->OverlapTimer);

  fprintf(GlobalData->stderrc,"* CGW Overlapper took %g seconds to compute %d overlaps\n",
	  TotalTimerT(&GlobalData->OverlapTimer, NULL), numOverlaps);


}


/* Added by Knut Reinert 05/26/99 */
/* Bug fixed by Knut Reinert 07/28/99 */

static int checkChunkOverlapCheckTOld(ChunkOverlapCheckT *co1,
                                      CDS_COORD_t minOverlap,
				      CDS_COORD_t maxOverlap,
                                      float errorRate){
  if( co1->errorRate != errorRate )
    return FALSE;
  if( co1->minOverlap != minOverlap )
    return FALSE;
  if( co1->maxOverlap != maxOverlap )
    return FALSE;

  return TRUE;
}

static int checkChunkOverlapCheckT(ChunkOverlapCheckT *co1,
                                   CDS_COORD_t minOverlap,
				   CDS_COORD_t maxOverlap,
                                   float errorRate)
{
  if( co1->errorRate != errorRate )
    return FALSE;
  if( minOverlap >= co1->minOverlap && maxOverlap <= co1->maxOverlap && 
      ( minOverlap <= co1->overlap && maxOverlap >= co1->overlap ))
    return TRUE;
  
  return FALSE;
}

static int checkChunkOverlapCheckTLoose(ChunkOverlapCheckT *co1,
                                        CDS_COORD_t minOverlap,
                                        CDS_COORD_t maxOverlap,
                                        float errorRate)
{
  if( co1->errorRate != errorRate )
    return FALSE;
  if( minOverlap < co1->minOverlap || maxOverlap > co1->maxOverlap)
    return FALSE;
  
  return TRUE;
}



CDS_COORD_t SmallOverlapExists(GraphCGW_T *graph,
                               CDS_CID_t cidA, CDS_CID_t cidB,
                               ChunkOrientationType orientation,
                               CDS_COORD_t minOverlap){
  ChunkOverlapCheckT olap;

  olap = OverlapChunks(graph, cidA, cidB, orientation,
                       minOverlap, CGW_DP_MINLEN + 5, CGW_DP_ERATE, FALSE);

  return olap.overlap;

}

CDS_COORD_t  LargeOverlapExists(GraphCGW_T *graph,
                                CDS_CID_t cidA, CDS_CID_t cidB,
                                ChunkOrientationType orientation,
                                CDS_COORD_t minOverlap,
                                CDS_COORD_t maxOverlap){
  ChunkOverlapCheckT olap;

  olap = OverlapChunks(graph, cidA, cidB, orientation,
                       minOverlap, maxOverlap, CGW_DP_ERATE, FALSE);

  return olap.overlap;

}

BranchPointResult OverlapChunksWithBPDetection(GraphCGW_T *graph,
                                               CDS_CID_t cidA, CDS_CID_t cidB,
                                               ChunkOrientationType orientation, 
                                               CDS_COORD_t minOverlap,
                                               CDS_COORD_t maxOverlap){

  BranchPointResult bpResult;
  ChunkOverlapCheckT olap;
  SeqInterval singleCoverageA, singleCoverageB;
  MultiAlignT *maA = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, cidA, graph->type == CI_GRAPH);
  MultiAlignT *maB = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, cidB, graph->type == CI_GRAPH);
  //  MultiAlignT *maA = GetMultiAlignInStore(graph->maStore, cidA);
  //  MultiAlignT *maB = GetMultiAlignInStore(graph->maStore, cidB);
  int endA, endB;
  CDS_COORD_t lengthA, lengthB;
  CDS_COORD_t singleLengthA, singleLengthB;
  CDS_COORD_t beg, end;
  int opposite;
#if 0
  int hasMultiCoverageA, hasMultiCoverageB;
  SeqInterval multiCoverageA, multiCoverageB;
  CDS_COORD_t multiLengthA, multiLengthB;
#endif

  InitCanonicalOverlapSpec(cidA,cidB,orientation,&olap.spec);
  
  /* Create consensusA,B and qualityA,B, if necessary */
  if(consensusA == NULL){
    consensusA = CreateVA_char(2048);
    consensusB = CreateVA_char(2048);
    qualityA = CreateVA_char(2048);
    qualityB = CreateVA_char(2048);
  }else{
    ResetVA_char(consensusA);
    ResetVA_char(consensusB);
    ResetVA_char(qualityA);
    ResetVA_char(qualityB);
  }

  // This is a waste, but get the length of the two contigs
  // Get the consensus sequences for both chunks from the ChunkStore
  lengthA = strlen(Getchar(maA->consensus,0));
  lengthB = strlen(Getchar(maB->consensus,0));

  bpResult.apnt = NULLINDEX;
  bpResult.bpnt = NULLINDEX;

  /* We have to deal with the orientation of the chunks */
  switch(olap.spec.orientation){
    case BA_BA:
      assert(0);
      break;
    case AB_AB:
      /* We use the chunks as retrieved from the chunk store */
      opposite = FALSE;
      endA = B_END;
      endB = A_END;
      break;
    case AB_BA:
      /* We need to reverse chunk B */
      opposite = TRUE;
      endA = B_END;
      endB = B_END;
      break;
    case BA_AB:
      /* We need to reverse both chunks, and then reverse chunk B ==> reverse chunkB */
      opposite = TRUE;
      endA = A_END;
      endB = A_END;
      break;
    default:
      assert(0);
  }

#if 0
  /* Identify the multiple coverage portions of the two contigs */
  hasMultiCoverageA = GetGappedMultipleCoverageInterval(graph, cidA, &multiCoverageA, endA);
  if(!hasMultiCoverageA){
    fprintf(GlobalData->stderrc,"# Contig " F_CID " has NO multi-coverage region...returning\n", cidA);
    return bpResult;
  }
  hasMultiCoverageB = GetGappedMultipleCoverageInterval(graph, cidB, &multiCoverageB, endB);
  if(!hasMultiCoverageB){
    fprintf(GlobalData->stderrc,"# Contig " F_CID " has NO multi-coverage region...returning\n", cidB);
    return bpResult;
  }

  multiLengthA = multiCoverageA.end - multiCoverageA.bgn;  
  multiLengthB = multiCoverageB.end - multiCoverageB.bgn;
  fprintf(GlobalData->stderrc,"* MultiCoverageA [" F_COORD "," F_COORD "]\n",
	  multiCoverageA.bgn, multiCoverageA.end);
  fprintf(GlobalData->stderrc,"* MultiCoverageB [" F_COORD "," F_COORD "]\n",
	  multiCoverageB.bgn, multiCoverageB.end);
#endif
  
  /* Output fasta records for single coverage portion of contig */
  if(endA == A_END){
    singleCoverageA.bgn = 0;
    singleCoverageA.end = MIN(500, lengthA);
  }else{
    singleCoverageA.bgn = MAX(0,lengthA - 500);
    singleCoverageA.end = lengthA;
  }
  if(endB == A_END){
    singleCoverageB.bgn = 0;
    singleCoverageB.end = MIN(500, lengthB);
  }else{
    singleCoverageB.bgn = MAX(0,lengthB - 500);
    singleCoverageB.end = lengthB;
  }
  
#if 0
  // Extract relative portions of consensus sequence 
  GetMultiAlignUngappedConsensusFromInterval(maA,multiCoverageA,consensusA, qualityA);
  GetMultiAlignUngappedConsensusFromInterval(maB,multiCoverageB,consensusB, qualityB);

  fprintf(GlobalData->stderrc,"* OverlapChunksWithMultipleCoverage A: " F_CID " [" F_COORD "," F_COORD "] of [" F_COORD "," F_COORD "]   B: " F_CID "  [" F_COORD "," F_COORD "] of [" F_COORD "," F_COORD "]\n",
	  cidA, multiCoverageA.bgn, multiCoverageA.end,0, lengthA,
	  cidB, multiCoverageB.bgn, multiCoverageB.end,0, lengthB);

  fflush(NULL);
  olap.overlap = 0;
  olap.minOverlap = minOverlap;
  olap.maxOverlap = MAX(multiLengthA, multiLengthB);

  /* Set up the beg,end, opposite arguments for the call to dp_compare */
  switch(olap.spec.orientation){
    case BA_BA:
      assert(0);
      break;
    case AB_AB:
      /* We use the chunks as retrieved from the chunk store */
      beg = multiLengthA - olap.maxOverlap;
      end = multiLengthA - olap.minOverlap;
      opposite = FALSE;
      break;
    case AB_BA:
      /* We need to reverse chunk B */
      beg = multiLengthA - olap.maxOverlap;
      end = multiLengthA - olap.minOverlap;
      opposite = TRUE;
      break;
    case BA_AB:
      /* We need to reverse both chunks, and then reverse chunk B ==> reverse chunkB */
      beg = -(multiLengthB - olap.minOverlap );
      end = -(multiLengthB - olap.maxOverlap );
      opposite = TRUE;
      break;
    default:
      assert(0);
  }
  // Add some slop
  beg -= 20;
  end += 20;
#else

  // Extract relative portions of consensus sequence 
  GetMultiAlignUngappedConsensusFromInterval(maA,singleCoverageA,consensusA, qualityA);
  GetMultiAlignUngappedConsensusFromInterval(maB,singleCoverageB,consensusB, qualityB);

  fprintf(GlobalData->stderrc,"* OverlapChunksWithBPDetection A: " F_CID " [" F_COORD "," F_COORD "] of [" F_COORD "," F_COORD "]   B: " F_CID "  [" F_COORD "," F_COORD "] of [" F_COORD "," F_COORD "]\n",
	  cidA, singleCoverageA.bgn, singleCoverageA.end,0, lengthA,
	  cidB, singleCoverageB.bgn, singleCoverageB.end,0, lengthB);

  fflush(NULL);

  singleLengthA = singleCoverageA.end - singleCoverageA.bgn;
  singleLengthB = singleCoverageB.end - singleCoverageB.bgn;
  olap.overlap = 0;
  olap.minOverlap = minOverlap;
  olap.maxOverlap = MAX(singleCoverageA.end - singleCoverageA.bgn, 
			singleCoverageB.end - singleCoverageB.bgn);

  /* Set up the beg,end, opposite arguments for the call to dp_compare */
  switch(olap.spec.orientation){
    case BA_BA:
      assert(0);
      break;
    case AB_AB:
      /* We use the chunks as retrieved from the chunk store */
      beg = singleLengthA - olap.maxOverlap;
      end = singleLengthA - olap.minOverlap;
      opposite = FALSE;
      break;
    case AB_BA:
      /* We need to reverse chunk B */
      beg = singleLengthA - olap.maxOverlap;
      end = singleLengthA - olap.minOverlap;
      opposite = TRUE;
      break;
    case BA_AB:
      /* We need to reverse both chunks, and then reverse chunk B ==> reverse chunkB */
      beg = -(singleLengthB - olap.minOverlap );
      end = -(singleLengthB - olap.maxOverlap );
      opposite = TRUE;
      break;
    default:
      assert(0);
  }
  // Add some slop
  beg -= 20;
  end += 20;

#endif  
  
  // Return value is length of chunk sequence/quality
  // **consensus will point to consensus data
  // **quality   will point to quality data
  // overlap 'em
  
  {
    BranchPointResult *BPResult = NULL;
    InternalFragMesg AFR,BFR;
    AFR.sequence   =  Getchar(consensusA, 0);
    AFR.quality    =  Getchar(qualityA, 0);
    AFR.eaccession = 0;
    AFR.iaccession = olap.spec.cidA;
    BFR.sequence   =  Getchar(consensusB, 0);
    BFR.quality    =  Getchar(qualityB, 0);
    BFR.eaccession = 0;
    BFR.iaccession = olap.spec.cidB;
    
#if 1
    fprintf(GlobalData->stderrc,"* Calling BPnt_Compare  with %s sequences lengths " F_SIZE_T "," F_SIZE_T " beg,end =[" F_COORD "," F_COORD "]\n",
	    (!strcmp(AFR.sequence, BFR.sequence)?"EQUAL":"DIFFERENT"),
	    strlen(AFR.sequence), strlen(BFR.sequence),
	    beg, end);
#endif
    //CGW_DP_ERATE=0.1 changed parameter to canOlap->errRate
    // a new field in ChunkOverlapCheckT

    BPResult = BPnt_Compare_AS(&AFR,&BFR, 
			       beg, end, opposite,
			       0.10,CGW_DP_THRESH, 
			       25,25);
    if(BPResult){
      bpResult = *BPResult;
    }


    return bpResult;

  }
}




ChunkOverlapCheckT OverlapChunks( GraphCGW_T *graph,
                                  CDS_CID_t cidA, CDS_CID_t cidB,
                                  ChunkOrientationType orientation, 
                                  CDS_COORD_t minOverlap,
                                  CDS_COORD_t maxOverlap,
                                  float errorRate,
                                  int insertGraphEdges)
{
  /* this function takes two chunks cidA and cidB, their orientation
     and an assumed minimum and maximum overlap for which it checks.
     It then tries to lookup whether the overlap was computed earlier
     or, if not, it computes the overlap and inserts the result in the
     hash table only if there was not such symbol there before.
  */
 
  int recompute = FALSE;
  int insert    = FALSE;
  int isCanonical;
  // If we init olap the return value is stored
  // here indicating whether the orientation of the two chunks
  // was canonical in Saul's definition or not.

  ChunkOverlapCheckT *lookup;
  // This pointer holds the return value of LookupCanonicalOverlap

  ChunkOverlapCheckT olap;
  // This is a temporary variable. The return value will be in lookup
  // or the pointer returned by the lookup following the insert


  isCanonical = InitCanonicalOverlapSpec(cidA,cidB,orientation,&olap.spec);
  // initalize olap with the chunk IDs and their orientation and record 
  // whether the orientation was already canonical or not
  
  lookup = LookupCanonicalOverlap(graph->overlapper,&olap.spec);
  // lookup whether the overlap had already been computed

  
  if( lookup != NULL ){
    olap = *lookup;
    if( checkChunkOverlapCheckT(lookup,minOverlap,maxOverlap,errorRate) == FALSE )
      recompute = TRUE;
    insert = FALSE;
  }
  else
    {
      recompute = TRUE;
      insert = insertGraphEdges; 
    }

  if( recompute == TRUE ){
    // compute new overlap and store it into the table
    // If it has an overlap add an edge mate to the CG
    // and return TRUE
    olap.computed      = FALSE;
    // olap.overlap       = 0;
    olap.overlap       = (minOverlap + maxOverlap) / 2;
    olap.minOverlap    = minOverlap;
    olap.maxOverlap    = maxOverlap;
    olap.fromCGB       = FALSE;
    olap.cgbMinOverlap = minOverlap;
    olap.cgbMaxOverlap = maxOverlap; 
    olap.errorRate     = errorRate;
    olap.suspicious = FALSE;

    {
      ChunkOverlapCheckT olap_copy;

      olap_copy = olap;

#ifdef USE_NEW_DP_COMPARE
      ComputeCanonicalOverlap_new(graph, &olap);
#else
      ComputeCanonicalOverlap_new(graph, &olap_copy);
	  
      // We compute the canonical overlap
      ComputeCanonicalOverlap(graph, &olap);
	  
      if (abs(olap.overlap - olap_copy.overlap) > 10 ||
          abs(olap.ahg - olap_copy.ahg) > 10 ||
          abs(olap.bhg - olap_copy.bhg) > 10)
        {
          char *seq1, *seq2;
          seq1   = Getchar(consensusA, 0);
          seq2   = Getchar(consensusB, 0);

          fprintf(GlobalData->stderrc, "nodes: " F_CID ", " F_CID ", orientation: %c\n",
                  olap.spec.cidA, olap.spec.cidB, orientation);
          fprintf(GlobalData->stderrc, "minOverlap, maxOverlap: " F_COORD ", " F_COORD "\n",
                  minOverlap, maxOverlap);
          fprintf( GlobalData->stderrc, "olap: overlap, ahg, bhg: " F_COORD ", " F_COORD ", " F_COORD "\n",
                   olap.overlap, olap.ahg, olap.bhg);
          fprintf( GlobalData->stderrc, "olap_copy: overlap, ahg, bhg: " F_COORD ", " F_COORD ", " F_COORD "\n",
                   olap_copy.overlap, olap_copy.ahg, olap_copy.bhg);
        }

      olap = olap_copy;
#endif
    }
	
    if(insert)
      { // Insert new entries in hashtable, if requested
        //
        int suspicious = olap.suspicious;
        if(olap.suspicious){
          olap.suspicious = FALSE;
          lookup = LookupCanonicalOverlap(graph->overlapper,&olap.spec); // see if it is already in the table
          if(!lookup){
            fprintf(GlobalData->stderrc,"* Inserting into hash (" F_CID "," F_CID ",%c) " F_COORD "\n",
                    olap.spec.cidA, olap.spec.cidB,
                    olap.spec.orientation, olap.overlap);
            if(InsertChunkOverlap(graph->overlapper, &olap) !=
               HASH_SUCCESS)
              assert(0);
          }
        }else{
          if(InsertChunkOverlap(graph->overlapper, &olap) != HASH_SUCCESS)
            assert(0);
        }
        lookup = LookupCanonicalOverlap(graph->overlapper,&olap.spec);
        assert(lookup != NULL);
        // ComputeCanonicalOverlap does not return the olap, so we look it up again.

        olap = *lookup;
        olap.suspicious = suspicious;
      }


    if(olap.overlap && insert){ // Insert all non-zero overlaps, if requested
      /*
        fprintf(GlobalData->stderrc,"* Inserting into graph (" F_CID "," F_CID ",%c) " F_COORD " ahang: " F_COORD ", bhang: " F_COORD "\n",
        olap.spec.cidA, olap.spec.cidB,
        olap.spec.orientation, olap.overlap,
        olap.ahg, olap.bhg);
      */
      InsertComputedOverlapEdge(graph, &olap);
    }
    /* Make sure the orientation of the edge we return is IDENTICAL with the spec returned */
    // if the input was not canonical we set the cid's and orientation
    // back to the input value (see als LookupOverlap)
    if( !olap.suspicious  && ! isCanonical ){
#if 0
      NodeCGW_T *a = GetGraphNode(graph, cidA);
      NodeCGW_T *b = GetGraphNode(graph, cidB);
#endif
      int swap;

      olap.spec.orientation = orientation;

      // If this is non-canonical,s wap things back
      olap.spec.cidA = cidA;
      olap.spec.cidB = cidB;
      swap = olap.BContainsA;
      olap.BContainsA = olap.AContainsB;
      olap.AContainsB = swap;
      swap = olap.ahg;
      olap.ahg = olap.bhg;
      olap.bhg = swap;

#if 0
      // NEW!!!!
      if(olap.AContainsB){
        assert(olap.ahg > 0 && olap.bhg < 0);
        switch(orientation){
          case AB_AB:
          case AB_BA:
            olap.overlap = b->bpLength.mean - olap.bhg;
            break;
          case BA_AB:
          case BA_BA:
            olap.overlap = b->bpLength.mean + olap.ahg;
            break;
        }
      }else if(olap.BContainsA){
        assert(olap.bhg > 0 && olap.ahg < 0);
		
        switch(orientation){
          case AB_AB:
          case AB_BA:
            olap.overlap = a->bpLength.mean - olap.ahg;
            break;
          case BA_AB:
          case BA_BA:
            olap.overlap = a->bpLength.mean + olap.bhg;
            break;
        }
      }
      // END NEW!
#else
      /*
        The following adjustments are unnecessary as they change the
        overlap value from the average of the two ways it can be
        calculated to one or the other. The average is more stable &
        reproducible.

        // RENEW!!!!
        if(olap.AContainsB){
        assert(olap.ahg > 0 && olap.bhg < 0);
        olap.overlap = b->bpLength.mean - olap.bhg;
        }else if(olap.BContainsA){
        assert(olap.bhg > 0 && olap.ahg < 0);
        olap.overlap = a->bpLength.mean - olap.ahg;
        }
        // END RENEW!
        */
#endif
    }
  }

  CheckScaffoldGraphCache(ScaffoldGraph); // flush the cache if it has gotten too big  

  if(olap.overlap==0){olap.quality=0;}
  return olap;
}

#if 0
Overlap* OverlapChunksNew( GraphCGW_T *graph,
                           CDS_CID_t cidA, CDS_CID_t cidB,
                           ChunkOrientationType orientation, 
                           CDS_COORD_t minOverlap, CDS_COORD_t maxOverlap,
                           float errorRate,
                           int insertGraphEdges)
{
  /* this function takes two chunks cidA and cidB, their orientation
     and an assumed minimum and maximum overlap for which it checks.
     It then tries to lookup whether the overlap was computed earlier
     or, if not, it computes the overlap and inserts the result in the
     hash table only if there was not such symbol there before.
  */

  int recompute = FALSE;
  int insert    = FALSE;
  int isCanonical;
  // If we init olap the return value is stored
  // here indicating whether the orientation of the two chunks
  // was canonical in Saul's definition or not.

  static int numOverlaps = 0;
  static int numHits = 0;

  ChunkOverlapCheckT *lookup;
  // This pointer holds the return value of LookupCanonicalOverlap

  ChunkOverlapCheckT olap;
  // This is a temporary variable. The return value will be in lookup
  // or the pointer returned by the lookup following the insert

  // make a static to mimic the behavoir of the one returned by DP_Compare
  static Overlap olapReturn;

  isCanonical = InitCanonicalOverlapSpec(cidA,cidB,orientation,&olap.spec);
  // initalize olap with the chunk IDs and their orientation and record 
  // whether the orientation was already canonical or not
  
  lookup = LookupCanonicalOverlap(graph->overlapper,&olap.spec);
  // lookup whether the overlap had already been computed

  
  if( lookup != NULL ){
    olap = *lookup;

    // this indicates that CreateChunkOverlapFromEdge inserted w/o ahg & bhg info
    if (lookup->ahg == 0 && lookup->bhg == 0)  
      {
        recompute = TRUE;
        if(DeleteChunkOverlap(graph->overlapper, &olap) != HASH_SUCCESS)
          assert(0);
        // RemoveComputedOverlapEdge(graph, &olap);
        insert = TRUE;
      }
    else
      {
        if( checkChunkOverlapCheckTLoose(lookup,minOverlap,maxOverlap,errorRate) == FALSE )
          recompute = TRUE;
        insert = FALSE;
      }
  }
  else
    {
      recompute = TRUE;
      insert = insertGraphEdges; 
    }

  // hit rate measuring
  numOverlaps++;
  if (recompute == FALSE)
    numHits++;
  
  if (!(numOverlaps % 10))
    fprintf( stderr, "numOverlaps: %d, numHits: %d\n",
             numOverlaps, numHits);

  if( recompute == TRUE )
    {
      // compute new overlap and store it into the table
      // If it has an overlap add an edge mate to the CG
      // and return TRUE
      olap.computed      = FALSE;
      // olap.overlap       = 0;
      olap.overlap       = (minOverlap + maxOverlap) / 2;
      olap.minOverlap    = minOverlap;
      olap.maxOverlap    = maxOverlap;
      olap.fromCGB       = FALSE;
      olap.cgbMinOverlap = minOverlap;
      olap.cgbMaxOverlap = maxOverlap; 
      olap.errorRate     = errorRate;
      olap.suspicious = FALSE;

      ComputeCanonicalOverlap_new( graph, &olap);
	
      if (insert)
        { // Insert new entries in hashtable, if requested
          //
          int suspicious = olap.suspicious;
          if (olap.suspicious)
            {
              olap.suspicious = FALSE;
              lookup = LookupCanonicalOverlap(graph->overlapper,&olap.spec); // see if it is already in the table
              if(!lookup)
		{
		  fprintf(GlobalData->stderrc,"* Inserting into hash (" F_CID "," F_CID ",%c) " F_COORD "\n",
                          olap.spec.cidA, olap.spec.cidB,
                          olap.spec.orientation, olap.overlap);
                  if(InsertChunkOverlap(graph->overlapper, &olap) !=
                     HASH_SUCCESS)
                    assert(0);
		}
            }
	  else
            {
              if(InsertChunkOverlap(graph->overlapper, &olap) != HASH_SUCCESS)
                assert(0);
            }
          lookup = LookupCanonicalOverlap(graph->overlapper,&olap.spec);
          assert(lookup != NULL);
          // ComputeCanonicalOverlap does not return the olap, so we look it up again.
	  
          olap = *lookup;
          olap.suspicious = suspicious;
        }

      if (olap.overlap && insert)
	{ // Insert all non-zero overlaps, if requested
	  /*
            fprintf(GlobalData->stderrc,"* Inserting into graph (" F_CID "," F_CID ",%c) " F_COORD " ahang: " F_COORD ", bhang: " F_COORD "\n",
            olap.spec.cidA, olap.spec.cidB,
            olap.spec.orientation, olap.overlap,
            olap.ahg, olap.bhg);
	  */
          InsertComputedOverlapEdge(graph, &olap);
        }
    }

  /* Make sure the orientation of the edge we return is IDENTICAL with the spec returned */
  // if the input was not canonical we set the cid's and orientation
  // back to the input value (see als LookupOverlap)
  if( !olap.suspicious  && ! isCanonical )
    {
      NodeCGW_T *a = GetGraphNode(graph, cidA);
      NodeCGW_T *b = GetGraphNode(graph, cidB);
      int swap;
	
      olap.spec.orientation = orientation;
	
      // If this is non-canonical,s wap things back
      olap.spec.cidA = cidA;
      olap.spec.cidB = cidB;
      swap = olap.BContainsA;
      olap.BContainsA = olap.AContainsB;
      olap.AContainsB = swap;
      swap = olap.ahg;
      olap.ahg = olap.bhg;
      olap.bhg = swap;

      // NEW!!!!
	
      if(olap.AContainsB)
	{
	  assert(olap.ahg > 0 && olap.bhg < 0);
	  switch(orientation)
            {
              case AB_AB:
              case AB_BA:
                olap.overlap = b->bpLength.mean + olap.bhg;
                break;
              case BA_AB:
              case BA_BA:
                olap.overlap = b->bpLength.mean - olap.ahg;
                break;
              default:
                assert(0);
                break;
            }
	}
      else if (olap.BContainsA)
	{
	  assert(olap.bhg > 0 && olap.ahg < 0);
	  switch(orientation)
            {
              case AB_AB:
              case AB_BA:
                olap.overlap = a->bpLength.mean - olap.bhg;
                break;
              case BA_AB:
              case BA_BA:
                olap.overlap = a->bpLength.mean + olap.ahg;
                break;
              default:
                assert(0);
                break;
            }
	}
      // END NEW!
    }
  
  CheckScaffoldGraphCache(ScaffoldGraph); // flush the cache if it has gotten too big

  // return olap;
  // convert olap (of type ChunkOverlapCheckT) to olapReturn (of type Overlap)
  if (olap.overlap != 0)
    {
      olapReturn.begpos = olap.ahg;
      olapReturn.endpos = olap.bhg;
      olapReturn.length = olap.overlap;  
      return (&olapReturn);
    }
  else
    return(NULL);
}
#endif


size_t ReportMemorySize_CO(ChunkOverlapperT *chunkOverlapper, FILE *stream){
  size_t totalMemorySize  = 0;
  totalMemorySize += ReportMemorySize_HT(chunkOverlapper->hashTable, "overlap",stream);
  totalMemorySize += ReportMemorySize_HP(chunkOverlapper->ChunkOverlaps, "heap",stream);
  return totalMemorySize;
}


/* This function takes two sequences, their orientation
   and an assumed minimum and maximum ahang for which it checks.
   
   This version uses the new DPCompare.
*/

Overlap* OverlapSequences( char *seq1, char *seq2,
                           ChunkOrientationType orientation, 
                           CDS_COORD_t min_ahang, CDS_COORD_t max_ahang,
                           double erate, double thresh, CDS_COORD_t minlen,
                           CompareOptions what)
{
  Overlap *omesg;
  int flip = 0;
  
  // if the orientation is BA_AB or BA_BA, we need to reverse complement the first contig
  if (orientation == BA_AB || orientation == BA_BA)
    Complement_Seq( seq1 );
  
  // if the orientation is AB_BA or BA_BA, we need to set the flip variable for the second contig
  if (orientation == AB_BA || orientation == BA_BA)
    flip = 1;

  // min_ahang and end are essentially bounds on the a-hang
  omesg = GlobalData->aligner(seq1, seq2,
                              min_ahang, max_ahang, (int) flip,
                              erate, thresh, minlen,
                              what);

  // return seq1 to its original state
  if (orientation == BA_AB || orientation == BA_BA)
    Complement_Seq( seq1 );

  // to protect against overlaps shorter than minlen that could be extended
  // by end-gaps to minlen long ...
  if(omesg!=NULL&&omesg->length<=minlen)return NULL;

  // omesg->begpos is the a-hang, omesg->endpos is the b-hang
  return omesg;
}




static VA_TYPE(char) *consensus1 = NULL;
static VA_TYPE(char) *consensus2 = NULL;
static VA_TYPE(char) *quality1 = NULL;
static VA_TYPE(char) *quality2 = NULL;

Overlap* OverlapContigs(NodeCGW_T *contig1, NodeCGW_T *contig2, 
                        ChunkOrientationType *overlapOrientation,
                        CDS_COORD_t minAhang, CDS_COORD_t maxAhang,
                        int computeAhang)
{
  Overlap *tempOlap1;
  char *seq1, *seq2;
  double erate, thresh;
  CDS_COORD_t minlen;

  /* 
     fprintf( GlobalData->stderrc, "\ncomputing overlap for contig1: " F_CID " and contig2: " F_CID "\n", 
     contig1->id, contig2->id);
     fprintf( GlobalData->stderrc, "orientation is %c\n", (char) *overlapOrientation);
  */

  erate = CGW_DP_ERATE;
  thresh = CGW_DP_THRESH;
  minlen = CGW_DP_MINLEN;

  // if computeAhang is TRUE, allow a lot of slide in ahang
  if (computeAhang == TRUE)
    {
      minAhang = - (CDS_COORD_t) contig2->bpLength.mean;
      maxAhang = (CDS_COORD_t) contig1->bpLength.mean;
      minlen -= 3;  // we subtract 3 because of an asymmetry in DP_COMPARE re AB_BA and BA_AB
    }
  /* 
     fprintf( stderr, "computeAhang is %s\n", ((computeAhang == TRUE) ? "TRUE" : "FALSE"));
     fprintf( stderr, "minAhang is " F_COORD "\n", minAhang);
     fprintf( stderr, "maxAhang is " F_COORD "\n", maxAhang);
  */
  if(consensus1 == NULL)
    {
      consensus1 = CreateVA_char(2048);
      consensus2 = CreateVA_char(2048);
      quality1 = CreateVA_char(2048);
      quality2 = CreateVA_char(2048);
    }else{
      ResetVA_char(consensus1);
      ResetVA_char(consensus2);
      ResetVA_char(quality1);
      ResetVA_char(quality2);
    }
  // Get the consensus sequences for both chunks from the Store
  GetConsensus(ScaffoldGraph->RezGraph, contig1->id, consensus1, quality1);
  GetConsensus(ScaffoldGraph->RezGraph, contig2->id, consensus2, quality2);

  seq1 = Getchar(consensus1, 0);
  seq2 = Getchar(consensus2, 0);

  // tempOlap1 is a static down inside of DP_Compare
  tempOlap1 =
    OverlapSequences( seq1, seq2,
                      *overlapOrientation, minAhang, maxAhang, 
                      erate, thresh, minlen,
                      (GlobalData->aligner == Local_Overlap_AS_forCNS) ?
                      AS_FIND_LOCAL_OVERLAP : AS_FIND_ALIGN);

  if (tempOlap1 != NULL)
    {
      /*
        fprintf(GlobalData->stderrc, F_CID ", " F_CID " ahang: " F_COORD ", bhang:" F_COORD "\n",
        contig1->id, contig2->id, tempOlap1->begpos, tempOlap1->endpos);
      */
      return tempOlap1;
    }
  else
    {
      fprintf(GlobalData->stderrc, F_CID ", " F_CID " do not overlap\n",
              contig1->id, contig2->id);
      // dumpContigInfo(contig1);
      // dumpContigInfo(contig2);

      return NULL;	
    }
}


Overlap* OverlapContigsLocal(NodeCGW_T *contig1, NodeCGW_T *contig2, 
                             ChunkOrientationType overlapOrientation,
                             CDS_COORD_t minAhang, CDS_COORD_t maxAhang,
                             int computeAhang)
{
  Overlap * tempOlap1;
  char *seq1, *seq2;
  double erate, thresh;
  CDS_COORD_t minlen;

  erate = CGW_DP_ERATE;
  thresh = CGW_DP_THRESH;
  minlen = CGW_DP_MINLEN;

  // if computeAhang is TRUE, allow a lot of slide in ahang
  if (computeAhang == TRUE)
    {
      minAhang = - (CDS_COORD_t) contig2->bpLength.mean;
      maxAhang = (CDS_COORD_t) contig1->bpLength.mean;
      minlen -= 3;  // we subtract 3 because of an asymmetry in DP_COMPARE re AB_BA and BA_AB
    }

  if(consensus1 == NULL)
    {
      consensus1 = CreateVA_char(2048);
      consensus2 = CreateVA_char(2048);
      quality1 = CreateVA_char(2048);
      quality2 = CreateVA_char(2048);
    }else{
      ResetVA_char(consensus1);
      ResetVA_char(consensus2);
      ResetVA_char(quality1);
      ResetVA_char(quality2);
    }
  // Get the consensus sequences for both chunks from the Store
  GetConsensus(ScaffoldGraph->RezGraph, contig1->id, consensus1, quality1);
  GetConsensus(ScaffoldGraph->RezGraph, contig2->id, consensus2, quality2);

  seq1 = Getchar(consensus1, 0);
  seq2 = Getchar(consensus2, 0);

  // tempOlap1 is a static down inside of DP_Compare
  tempOlap1 = Local_Overlap_AS_forCNS(seq1, seq2,
                                      minAhang,
                                      maxAhang,
                                      (overlapOrientation == AB_BA ||
                                       overlapOrientation == BA_AB),
                                      erate, thresh, minlen,
                                      AS_FIND_LOCAL_OVERLAP);
  return tempOlap1;
}

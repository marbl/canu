
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
static char *rcsid = "$Id: ChunkOverlap_CGW.c,v 1.58 2012-06-10 05:52:34 brianwalenz Exp $";

#include "ChunkOverlap_CGW.h"
#include "AS_UTL_reverseComplement.h"
#include "ScaffoldGraph_CGW.h"    // For DeleteCIOverlapEdge

#undef DEBUG_OVERLAP_SEQUENCES

// this is the initial range we use to compute
// overlaps with a different min/max range
#define BEGENDSLOP 10

// This will enable the use of Local_Overlap_AS_forCNS as a fallback
// when DP_Compare fails.  It will slow down cgw by over an order of
// magnitude (cgw on some microbes was taking more than 12 hours).
//
// In p.cnpt3, enabling this found 3 more overlaps.
//
//      3 only found with Local_Overlap
//     20 only found with DP_Compare
//    618 found by both
//   2711 found by none
//   3352 tries
//
#undef USE_LOCAL_OVERLAP_AS_FALLBACK


/* ChunkOverlap_CGW provides tools for invoking Gene's dpalign tool to compute
   overlaps between chunks.  Such overlaps are first 'collected' and then 'computed'.
   Then, they may be 'looked up'.

   This sparse 'database' of relationships between chunks (both overlaps and lack of overlaps) is
   stored in a symbol table based on AS_UTL_Hash to facilitate quick lookups and insertions.
   Key to this database is a canonical overlap representation.  The canonical storage for a
   potential overlap is as follows:
   if(orientation is symettric -- AB_BA or BA_AB )
   then the overlap is stored with the lower numbered chunk first
   if(orientation is antinormal -- BA_BA)
   then the overlap is stored as an normal overlap with the chunks reversed

   The idea of the collection phase, is to determine the maximal overlap range between
   two chunks that is worth computing.  So, during the construction of the raw
   mate edges for the extended chunk graph, any mate link that implies a possible overlap is
   collected for later evaluation.  If multiple overlaps are collected for the same chunk pair
   and orientation, the maximal overlap interval is collected.

   Once a set of inter-chunk relationships have been collected, Gene's dpalign tool is invoked
   on the consensus sequence for each chunk, and the results are stored in the database.  As
   an option, these overlap edges can also be added to the extended chunk graph.  The implied
   overlaps that are detected are thus added to the set of raw edge mates prior to merging.

   Once the extended chunk graph has been merged, we look for potential overlaps between
   unique chunks that are implied by the confirmed edges are checked.  These overlaps are NOT
   added tot he extended chunk graph, but are simply stored in the database for later
   retrieval in the scaffold construction code.

   Saul Kravitz
   April 1999


   Additions made by Knut Reinert
   09/10/99

   We make use of the above hash table to compute and store quality values with the overlaps
   Initially no meaningful quality value is stored in the ChunkOverlapCheckT::overlap
   We add bits to the ChunkOverlapCheckT struct to indicate whether there was a quality
   computation, and which function computed the quality value.
   Whenever a new method for computing the quality of an overlap should be tested,
   one should augment ChunkOverlapCheckT by the appropriate bit and add a value
   to the enum QualityFuncsT.

*/



static VA_TYPE(char) *consensusA = NULL;
static VA_TYPE(char) *consensusB = NULL;
static VA_TYPE(char) *qualityA = NULL;
static VA_TYPE(char) *qualityB = NULL;




//  Hash table support functions
//
static
int CanOlapCmp(uint64 cO1, uint64 cO2){
  ChunkOverlapSpecT *c1 = (ChunkOverlapSpecT *)(INTPTR)cO1;
  ChunkOverlapSpecT *c2 = (ChunkOverlapSpecT *)(INTPTR)cO2;

  if (c1->cidA < c2->cidA)
    return(-1);
  if (c1->cidA > c2->cidA)
    return(1);

  if (c1->cidB < c2->cidB)
    return(-1);
  if (c1->cidB > c2->cidB)
    return(1);

  if (c1->orientation.toLetter() < c2->orientation.toLetter())
    return(-1);
  if (c1->orientation.toLetter() > c2->orientation.toLetter())
    return(1);

  return(0);
}

static
int CanOlapHash(uint64 cO, uint32 length){
  ChunkOverlapSpecT *spec = (ChunkOverlapSpecT *)(INTPTR)cO;
  uint64             arr[3];

  assert(length == sizeof(ChunkOverlapSpecT));

  //  Why?!  Because of the chance that the unused bytes in ChunkOverlapSpecT are
  //  corrupting the hash values.

  arr[0] = spec->cidA;
  arr[1] = spec->cidB;
  arr[2] = spec->orientation.toLetter();

  return  Hash_AS((uint8 *)arr, sizeof(uint64) * 3, 37);
}







//external
ChunkOverlapperT *CreateChunkOverlapper(void){
  ChunkOverlapperT *chunkOverlapper = (ChunkOverlapperT *)safe_malloc(sizeof(ChunkOverlapperT));
  chunkOverlapper->hashTable = CreateGenericHashTable_AS(CanOlapHash, CanOlapCmp);
  chunkOverlapper->ChunkOverlaps = AllocateHeap_AS(sizeof(ChunkOverlapCheckT));
  return chunkOverlapper;
}

//external
void DestroyChunkOverlapper(ChunkOverlapperT *chunkOverlapper){
  DeleteHashTable_AS(chunkOverlapper->hashTable);
  FreeHeap_AS(chunkOverlapper->ChunkOverlaps);
  safe_free(chunkOverlapper);
}



static
void
printChunkOverlapCheckT(char *label, ChunkOverlapCheckT *olap) {
  fprintf(stderr, "%s %d,%d,%c - min/max %d/%d %d/%d erate %f flags %d%d%d%d%d overlap %d hang %d,%d qual %f offset %d,%d\n",
          label,
          olap->spec.cidA, olap->spec.cidB, olap->spec.orientation.toLetter(),
          olap->minOverlap, olap->maxOverlap,
          olap->cgbMinOverlap, olap->cgbMaxOverlap,
          olap->errorRate,
          olap->computed, olap->fromCGB, olap->AContainsB, olap->BContainsA, olap->suspicious,
          olap->overlap, olap->ahg, olap->bhg, olap->quality, olap->min_offset, olap->max_offset);
}




static
int
ExistsChunkOverlap(ChunkOverlapperT *chunkOverlapper,
                   ChunkOverlapCheckT *olap) {
  return(ExistsInHashTable_AS(chunkOverlapper->hashTable,
                              (uint64)(INTPTR)&olap->spec,
                              sizeof(ChunkOverlapSpecT)));
}


static
void
DeleteChunkOverlap(ChunkOverlapperT *chunkOverlapper,
                   ChunkOverlapCheckT *olap){

  ChunkOverlapCheckT  *del = LookupCanonicalOverlap(chunkOverlapper, &olap->spec);

  if (del == NULL)
    return;

  if (DeleteFromHashTable_AS(chunkOverlapper->hashTable,
                             (uint64)(INTPTR)&olap->spec,
                             sizeof(ChunkOverlapSpecT)) != HASH_SUCCESS) {
    //  This is triggered if the overlap DOESN'T exist in the table.  Which, since
    //  we just checked that it does exist, should never occur.
    printChunkOverlapCheckT("WARNING:  Failed to delete overlap", olap);
    assert(0);
  }

  if (ExistsChunkOverlap(chunkOverlapper, olap) == TRUE) {
    //  This is triggered if the overlap STILL exists in the table.
    printChunkOverlapCheckT("WARNING:  Deleted overlap still exists", olap);
    assert(0);
  }

  //  And finally, nuke the original data.
  memset(del, 0xff, sizeof(ChunkOverlapCheckT));
}

//external
int InsertChunkOverlap(ChunkOverlapperT *chunkOverlapper,
                       ChunkOverlapCheckT *olap){
  ChunkOverlapCheckT *nolap = (ChunkOverlapCheckT *)GetHeapItem_AS(chunkOverlapper->ChunkOverlaps);

  *nolap = *olap;

  assert((nolap->overlap == 0) || (nolap->errorRate >= 0.0));

  int exists = ExistsInHashTable_AS(chunkOverlapper->hashTable,
                                    (uint64)(INTPTR)&nolap->spec,
                                    sizeof(ChunkOverlapSpecT));

  if (exists == HASH_SUCCESS) {
    ChunkOverlapCheckT  *old = LookupCanonicalOverlap(chunkOverlapper, &nolap->spec);

    fprintf(stderr, "WARNING:  InsertChunkOverlap()-- Chunk overlap already exists.\n");
    printChunkOverlapCheckT("NEW", nolap);
    printChunkOverlapCheckT("OLD", old);
  }

  int inserted = InsertInHashTable_AS(chunkOverlapper->hashTable,
                                    (uint64)(INTPTR)&nolap->spec,
                                    sizeof(ChunkOverlapSpecT),
                                    (uint64)(INTPTR)nolap,
                                    0);

  //  Either the overlap didn't exist (in which case we should have inserted it), or
  //  the overlap existed already (in which case we should not have inserted it).
  //
  assert((exists == HASH_FAILURE) || (inserted == HASH_FAILURE));

  return(inserted);
}


//external
void  SaveChunkOverlapperToStream(ChunkOverlapperT *chunkOverlapper, FILE *stream){
  HashTable_Iterator_AS iterator;
  uint64 key, value;
  uint32 valuetype;
  int64 numOverlaps = 0;

  //  Iterate over all hashtable elements, just to count them

  InitializeHashTable_Iterator_AS(chunkOverlapper->hashTable, &iterator);

  while (NextHashTable_Iterator_AS(&iterator, &key, &value, &valuetype))
    numOverlaps++;

  AS_UTL_safeWrite(stream, &numOverlaps, "SaveChunkOverlapperToStream", sizeof(int64), 1);

  // Iterate over all hashtable elements, writing

  InitializeHashTable_Iterator_AS(chunkOverlapper->hashTable, &iterator);

  while(NextHashTable_Iterator_AS(&iterator, &key, &value, &valuetype)){
    ChunkOverlapCheckT *olap = (ChunkOverlapCheckT*)(INTPTR)value;

    AS_UTL_safeWrite(stream, olap, "SaveChunkOverlapperToStream", sizeof(ChunkOverlapCheckT), 1);
  }
}


//external
ChunkOverlapperT *  LoadChunkOverlapperFromStream(FILE *stream){
  ChunkOverlapperT  *chunkOverlapper = CreateChunkOverlapper();
  int64              numOverlaps;

  int status = AS_UTL_safeRead(stream, &numOverlaps, "LoadChunkOverlapperFromStream", sizeof(int64), 1);

  assert(status == 1);

  for (int64 overlap = 0; overlap < numOverlaps; overlap++) {
    ChunkOverlapCheckT olap;

    status = AS_UTL_safeRead(stream, &olap, "LoadChunkOverlapperFromStream", sizeof(ChunkOverlapCheckT), 1);

    assert(status == 1);
    assert(olap.errorRate > 0.0);

    if (InsertChunkOverlap(chunkOverlapper, &olap) != HASH_SUCCESS) {
      fprintf(stderr, "LoadChunkOverlapperFromStream()-- WARNING:  Duplicate chunk overlap!\n");
      assert(0);
    }
  }

  return chunkOverlapper;
}




/************************************************************************
 * A canonical overlap hs the following properties
 *       if(orientation is symmetric -- AB_BA or BA_AB )
 *       	then the overlap is stored with the lower numbered chunk first
 *	if(orientation is antinormal -- BA_BA)
 *	        then the overlap is stored as an normal overlap with the chunks reversed
 *************************************************************************/

//external
int
InitCanonicalOverlapSpec(CDS_CID_t           cidA,
                         CDS_CID_t           cidB,
                         PairOrient          orientation,
                         ChunkOverlapSpecT  *spec){
  int canonical = TRUE;

  assert(orientation.isUnknown() == false);

  memset(spec, 0, sizeof(ChunkOverlapSpecT));

  if (orientation.isBA_BA()) {
    spec->orientation.setIsAB_AB();
    spec->cidA = cidB;
    spec->cidB = cidA;
    canonical = FALSE;
  }

  if (orientation.isAB_AB()) {
    spec->orientation = orientation;
    spec->cidA = cidA;
    spec->cidB = cidB;
  }

  if (orientation.isBA_AB() || orientation.isAB_BA()) {
    spec->orientation = orientation;
    if(cidA > cidB){
      spec->cidA = cidB;
      spec->cidB = cidA;
      canonical = FALSE;
    }else{
      spec->cidA = cidA;
      spec->cidB = cidB;
    }
  }

  return canonical;
}



/* Given a graph edge, create an overlap in the hashtable and mark it as computed */
//external
void
CreateChunkOverlapFromEdge(GraphCGW_T *graph, EdgeCGW_T *edge){
  ChunkOverlapCheckT olap;
  double             delta = sqrt(edge->distance.variance) * 3.0;

  assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));

  memset(&olap, 0, sizeof(ChunkOverlapCheckT));

  InitCanonicalOverlapSpec(edge->idA, edge->idB, edge->orient, &olap.spec);

  olap.computed      = TRUE;
  olap.overlap       = -edge->distance.mean;
  olap.minOverlap    = -edge->distance.mean - delta;
  olap.maxOverlap    = -edge->distance.mean + delta;;
  olap.fromCGB       = FALSE;
  olap.cgbMinOverlap = 0;
  olap.cgbMaxOverlap = 0;
  olap.errorRate     = AS_CGW_ERROR_RATE;
  olap.quality       = edge->quality;
  olap.ahg           = 0;
  olap.bhg           = 0;
  olap.min_offset    = 0;
  olap.max_offset    = 0;

  //  Delete any overlap that is there already.  BPW is't sure if we want to REPLACE
  //  any overlap that exists with the edge-based overlap, or leave it alone.  I'm
  //  thinking that we DO want to replace, just because, for whatever reason, the
  //  client didn't use the existing overlap, or didn't find it, or just didn't look.
  //  Plus, it wasn't clear if the original intent was to replace or ignore.
  //
  //DeleteChunkOverlap(ScaffoldGraph->ChunkOverlaps, &olap);

  //  And add the new one.
  if (InsertChunkOverlap(ScaffoldGraph->ChunkOverlaps, &olap) != HASH_SUCCESS) {
    ChunkOverlapCheckT  *old = LookupCanonicalOverlap(ScaffoldGraph->ChunkOverlaps, &olap.spec);

    fprintf(stderr, "WARNING:  CreateChunkOverlapFromEdge()-- Chunk overlap already exists.  Keeping old overlap.\n");
    printChunkOverlapCheckT("NEW", &olap);
    printChunkOverlapCheckT("OLD", old);
    //assert(0);
  }
}



/* Given a graph edge, create an overlap in the hashtable */
//external
void FillChunkOverlapWithEdge(EdgeCGW_T *edge, ChunkOverlapCheckT *olap){
  double delta = sqrt(edge->distance.variance) * 3.0;
  assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
  InitCanonicalOverlapSpec(edge->idA, edge->idB, edge->orient, &olap->spec);
  olap->computed = FALSE;
  olap->overlap = -edge->distance.mean;

  // might be unsafe for big variances after tandem mark propagation
  olap->minOverlap = (int32) -edge->distance.mean - delta;
  olap->maxOverlap = (int32) -edge->distance.mean + delta;

  olap->minOverlap = (int32) -edge->distance.mean - 20;
  olap->maxOverlap = (int32) -edge->distance.mean + 20;
  olap->fromCGB = FALSE;
  olap->cgbMinOverlap = 0;
  olap->cgbMaxOverlap = 0;
  olap->errorRate = AS_CGW_ERROR_RATE;
  olap->quality = edge->quality;
  olap->ahg = olap->bhg = 0;
  olap->min_offset = olap->max_offset = 0;
}



/* Given a graph edge, create an overlap in the hashtable */
static
void FillChunkOverlapWithOVL(GraphCGW_T   *graph,
                             OverlapMesg  *ovl) {

  ChunkOverlapCheckT  olap;

  CIFragT            *cifraga, *cifragb;

  CDS_CID_t           cia;
  CDS_CID_t           cib;

  int                 bega, enda, lena;
  int                 begb, endb, lenb;

  assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));

  //  Find the chunks with the two fragments

  cifraga = GetCIFragT(ScaffoldGraph->CIFrags, ovl->aifrag);
  cia     = cifraga->cid; // cifrag.CIid;
  bega    = MIN(cifraga->offset5p.mean, cifraga->offset3p.mean);
  enda    = MAX(cifraga->offset5p.mean, cifraga->offset3p.mean);

  cifragb = GetCIFragT(ScaffoldGraph->CIFrags, ovl->bifrag);
  cib     = cifragb->cid; // cifrag.CIid;
  begb    = MIN(cifragb->offset5p.mean, cifragb->offset3p.mean);
  endb    = MAX(cifragb->offset5p.mean, cifragb->offset3p.mean);

  if (consensusA == NULL) {
    //  These are globals used by ComputeCanonicalOverlap_new().
    consensusA = CreateVA_char(2048);
    consensusB = CreateVA_char(2048);
    qualityA   = CreateVA_char(2048);
    qualityB   = CreateVA_char(2048);
  }

  //  Length of the ungapped consensus of each unitig.
  lena = GetConsensus(graph, cia, consensusA, qualityA);
  lenb = GetConsensus(graph, cib, consensusB, qualityB);

  //  Based on the fragment placement in the two unitigs, determine
  //  the size of the overlap we are expecting.
  //
  //  first orientation is one of
  //     AS_READ_ORIENT_INNIE
  //     AS_READ_ORIENT_NORMAL
  //     AS_READ_ORIENT_ANTINORMAL
  //     AS_READ_ORIENT_OUTTIE
  //  which are different from the overlap orientations.

  int  oo = AS_READ_ORIENT_UNKNOWN;

  if (ovl->orientation.isInnie())   oo = AS_READ_ORIENT_INNIE;
  if (ovl->orientation.isNormal())  oo = AS_READ_ORIENT_NORMAL;
  if (ovl->orientation.isAnti())    oo = AS_READ_ORIENT_ANTINORMAL;
  if (ovl->orientation.isOuttie())  oo = AS_READ_ORIENT_OUTTIE;

  PairOrient orient = ciEdgeOrientFromFragment(oo, getCIFragOrient(cifraga), getCIFragOrient(cifragb));

  //  Adjust positions, orientation.  We place unitig A starting at 0, then
  //  figure out where to place unitig B based on the overlap orientation.

  int utgabeg=0, utgaend=lena;
  int utgbbeg=0, utgbend=lenb;
  int olapsize = 0;
  int olapmin  = 0;
  int olapmax  = 0;

  assert(orient.isUnknown() == false);

  if (orient.isAB_AB()) {
    utgbbeg = (bega)        + ovl->ahg - (begb);
    utgbend = utgbbeg + lenb;
  }

  if (orient.isAB_BA()) {
    utgbbeg = (bega)        + ovl->ahg - (lenb - endb);
    utgbend = utgbbeg + lenb;
  }

  if (orient.isBA_AB()) {
    utgbbeg = (lena - enda) + ovl->ahg - (begb);
    utgbend = utgbbeg + lenb;
  }

  if (orient.isBA_BA()) {
    utgbbeg = (lena - enda) + ovl->ahg - (lenb - endb);
    utgbend = utgbbeg + lenb;
  }

  if (utgabeg < utgbbeg)
    if (utgaend < utgbend)
      olapsize = utgaend - utgbbeg;
    else
      olapsize = lenb;
  else
    if (utgaend < utgbend)
      olapsize = lena;
    else
      olapsize = utgbend - utgabeg;

  olapmin = olapsize - 0.1 * olapsize;
  olapmax = olapsize + 0.1 * olapsize;

  olap = OverlapChunks(graph,
                       cia, cib,
                       orient,
                       olapmin,
                       olapmax,
                       AS_CGW_ERROR_RATE,
                       TRUE);

#if 0
  fprintf(stderr, "TRY: utga %d-%d utgb %d-%d ovl %d (%d,%d) frgori %c utgori %c\n",
          utgabeg, utgaend,
          utgbbeg, utgbend,
          olapsize, olapmin, olapmax,
          ovl->orientation,
          orient);

  if (olap.overlap > 0)
    fprintf(stderr, "FillChunkOverlapWithOVL()-- frg "F_IID" in utg "F_IID" <-> frg "F_IID" in utg "F_IID" ovl %d min %d,%d hang %d,%d\n",
            ovl->aifrag, cia, ovl->bifrag, cib,
            olap.overlap, olapmin, olapmax, olap.ahg, olap.bhg);
  else
    fprintf(stderr, "FillChunkOverlapWithOVL()-- frg "F_IID" in utg "F_IID" <-> frg "F_IID" in utg "F_IID" ovl %d min %d,%d FAILED\n",
            ovl->aifrag, cia, ovl->bifrag, cib,
            olapsize, olapmin, olapmax);
#endif
}


//external
ChunkOverlapCheckT *LookupCanonicalOverlap(ChunkOverlapperT *chunkOverlapper,
                                           ChunkOverlapSpecT *spec){
  assert(spec->orientation.isValid() == true);
  return (ChunkOverlapCheckT *)(INTPTR)LookupValueInHashTable_AS(chunkOverlapper->hashTable,
                                                                 (uint64)(INTPTR)spec,
                                                                 sizeof(ChunkOverlapSpecT));
}


// Returns FALSE if none found
// Returns TRUE and overwrites *olap if found
//
//external
int
LookupOverlap(GraphCGW_T *graph,
              CDS_CID_t cidA, CDS_CID_t cidB,
              PairOrient orientation,
              ChunkOverlapCheckT *olap) {
  ChunkOverlapSpecT   spec;
  int                 isCanonical = InitCanonicalOverlapSpec(cidA, cidB, orientation, &spec);
  ChunkOverlapCheckT *lookup      = LookupCanonicalOverlap(ScaffoldGraph->ChunkOverlaps, &spec);

  if (lookup == NULL)
    // We didn't find anything
    return(FALSE);

  *olap = *lookup;

  if (isCanonical)
    return(TRUE);

  // Otherwise we have to fix up the retrieved canonical overlap for the non-canonical query
  //
  olap->spec.orientation = orientation;
  olap->spec.cidA = cidA;
  olap->spec.cidB = cidB;

  if(olap->BContainsA || olap->AContainsB) {
    NodeCGW_T *a = GetGraphNode(graph, cidA);
    NodeCGW_T *b = GetGraphNode(graph, cidB);

    int swap = olap->BContainsA;
    olap->BContainsA = olap->AContainsB;
    olap->AContainsB = swap;

    assert(orientation.isUnknown() == false);

    if (olap->AContainsB) {
      if (orientation.isAB_AB() || orientation.isAB_BA())
        olap->overlap = b->bpLength.mean + olap->bhg;
      if (orientation.isBA_AB() || orientation.isBA_BA())
        olap->overlap = b->bpLength.mean - olap->ahg;

    } else if (olap->BContainsA) {
      if (orientation.isAB_AB() || orientation.isAB_BA())
        olap->overlap = a->bpLength.mean - olap->bhg;
      if (orientation.isBA_AB() || orientation.isBA_BA())
        olap->overlap = a->bpLength.mean + olap->ahg;
    }
  }

  return(TRUE);
}



/* Insert a computed overlap as a CIEdgeT into the Scaffold Graph. */

//external
CDS_CID_t InsertComputedOverlapEdge(GraphCGW_T *graph,
                                    ChunkOverlapCheckT *olap){
  CDS_CID_t eid;
  int fudge;
  int isDoveTail = FALSE;
  LengthT overlap;
  PairOrient orient = olap->spec.orientation;
  EdgeCGW_T *existing = FindGraphOverlapEdge(graph, olap->spec.cidA, olap->spec.cidB, orient);

  overlap.mean   = -olap->overlap;
  overlap.variance = MAX(1.0, ComputeFudgeVariance((double)olap->overlap));
  fudge = sqrt(overlap.variance);

  isDoveTail = !(olap->AContainsB || olap->BContainsA);


  // If there is an existing overlap edge, don't insert this one if it is the same!
  if (existing) {
    double diff = existing->distance.mean - overlap.mean;

    if ((-5.0 < diff) &&
        (diff < 5.0))
      return(GetVAIndex_EdgeCGW_T(graph->edges, existing));
  }

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
                     isDoveTail,  // isOverlap
                     olap->AContainsB,                 // isAContainsB
                     olap->BContainsA,                 // isBContainsA
                     FALSE,                        // isTransChunk
                     FALSE, FALSE,  // extremalA and extremalB
                     UNKNOWN_EDGE_STATUS,
                     FALSE,
                     TRUE /* insert*/ );

  return eid;
}






//external
void CollectChunkOverlap(GraphCGW_T *graph,
                         CDS_CID_t cidA, CDS_CID_t cidB,
                         PairOrient orientation,
                         double   meanOverlap, double   deltaOverlap,
                         double   quality, int bayesian,
                         int fromCGB,
			 int verbose){
  ChunkOverlapCheckT canOlap, *olap;
  int32 delta;
  int32 minOverlap,maxOverlap;

  delta = (int32)(3.0 * deltaOverlap);
  delta = MAX(delta, 10);
  minOverlap = MAX(meanOverlap - delta, 0);
  maxOverlap = meanOverlap + delta;

  if(maxOverlap < 0){
    // There is no chance that these overlap!
    return;
  }

  memset(&canOlap, 0, sizeof(ChunkOverlapCheckT));

  // Create a canonical representation of the overlap
  InitCanonicalOverlapSpec(cidA, cidB, orientation, &canOlap.spec);

  // Lookup to see if we've already stored such an overlap
  olap = LookupCanonicalOverlap(ScaffoldGraph->ChunkOverlaps, &canOlap.spec);
  if(!olap){
    assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
    canOlap.computed = FALSE;
    canOlap.overlap = FALSE;
    canOlap.quality = 1.0;
    canOlap.minOverlap = minOverlap;
    canOlap.maxOverlap = maxOverlap;
    canOlap.fromCGB = fromCGB;
    if(fromCGB && bayesian){
      canOlap.computed = TRUE;
      canOlap.quality = quality;
      canOlap.overlap = (canOlap.minOverlap + canOlap.maxOverlap)/2;
    }
    canOlap.cgbMinOverlap = minOverlap;
    canOlap.cgbMaxOverlap = maxOverlap;
    canOlap.errorRate = AS_CGW_ERROR_RATE;

    // Add it to the symbol table
    if(InsertChunkOverlap(ScaffoldGraph->ChunkOverlaps, &canOlap) != HASH_SUCCESS)
      assert(0);

  }else{ // we found an overlap
    // We found one.  So now we need to update the existing record so that
    // it reflects the maximal interval we're interested in overlapping.
    if(olap->computed){
      // If we've already computed this one, only recompute if the range is expanded
      if(olap->overlap == 0 &&
	 (minOverlap < olap->minOverlap ||
          maxOverlap > olap->maxOverlap)){
	olap->computed = FALSE; // Recompute
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
	if(bayesian){
	  olap->overlap = (olap->cgbMinOverlap + olap->cgbMaxOverlap)/2;
	  olap->computed = TRUE;
	}

      }
      olap->minOverlap = MIN(minOverlap, olap->minOverlap);
      olap->maxOverlap = MAX(maxOverlap, olap->maxOverlap);
    }
  }
}






//external
ALNoverlap* OverlapSequences(char *seq1, char *seq2,
                             PairOrient orientation,
                             int32 min_ahang, int32 max_ahang,
                             double erate, double thresh, int32 minlen,
                             uint32 tryLocal)
{
  ALNoverlap *dp_omesg = NULL;
  ALNoverlap *lo_omesg = NULL;
  int      flip = 0;
  int      len1 = strlen(seq1);
  int      len2 = strlen(seq2);

  //  This function takes two sequences, their orientation and an
  //  assumed minimum and maximum ahang for which it checks.

  // if the orientation is BA_AB or BA_BA, we need to reverse
  // complement the first contig
  if (orientation.isBA_AB() || orientation.isBA_BA())
    reverseComplementSequence( seq1, len1 );

  // if the orientation is AB_BA or BA_BA, we need to set the flip
  // variable for the second contig
  if (orientation.isAB_BA() || orientation.isBA_BA())
    flip = 1;

  // min_ahang and end are essentially bounds on the a-hang

  dp_omesg = DP_Compare(seq1, seq2,
                        min_ahang, max_ahang,
                        len1, len2,
                        flip,
                        erate, thresh, minlen,
                        AS_FIND_ALIGN);

  if ((dp_omesg != NULL) && (dp_omesg->length <= minlen)) {
#ifdef DEBUG_OVERLAP_SEQUENCES
    fprintf(stderr, "OverlapSequences()-- Found overlap with DP_Compare   begpos=%d endpos=%d length=%d diffs=%d comp=%d/%d ISSHORT\n",
            dp_omesg->begpos, dp_omesg->endpos, dp_omesg->length, dp_omesg->diffs, dp_omesg->comp, flip);
#endif
    dp_omesg = NULL;
  }

  if ((dp_omesg != NULL) && ((double)dp_omesg->diffs / dp_omesg->length > erate)) {
#ifdef DEBUG_OVERLAP_SEQUENCES
    fprintf(stderr, "OverlapSequences()-- Found overlap with DP_Compare   begpos=%d endpos=%d length=%d diffs=%d comp=%d/%d ISCRAP\n",
            dp_omesg->begpos, dp_omesg->endpos, dp_omesg->length, dp_omesg->diffs, dp_omesg->comp, flip);
#endif
    dp_omesg = NULL;
  }

#ifdef DEBUG_OVERLAP_SEQUENCES
  if (dp_omesg != NULL)
    fprintf(stderr, "OverlapSequences()-- Found overlap with DP_Compare   begpos=%d endpos=%d length=%d diffs=%d comp=%d/%d\n",
            dp_omesg->begpos, dp_omesg->endpos, dp_omesg->length, dp_omesg->diffs, dp_omesg->comp, flip);
#endif


#ifdef USE_LOCAL_OVERLAP_AS_FALLBACK
  if (!dp_omesg)
    lo_omesg = Local_Overlap_AS_forCNS(seq1, seq2,
                                       min_ahang, max_ahang,
                                       len1, len2,
                                       flip,
                                       erate, thresh, minlen,
                                       AS_FIND_LOCAL_OVERLAP);
#else
if (tryLocal == TRUE && !dp_omesg) {
    lo_omesg = Local_Overlap_AS_forCNS(seq1, seq2,
                                       min_ahang, max_ahang,
                                       len1, len2,
                                       flip,
                                       erate, thresh, minlen,
                                       AS_FIND_LOCAL_OVERLAP);
}
#endif

  if (lo_omesg != NULL) {
    if ((lo_omesg != NULL) && (lo_omesg->length <= minlen)) {
      lo_omesg = NULL;
    }
    if ((lo_omesg != NULL) && ((double)lo_omesg->diffs / lo_omesg->length > erate)) {
      lo_omesg = NULL;
    }

#ifdef DEBUG_OVERLAP_SEQUENCES
    if (lo_omesg != NULL)
      fprintf(stderr, "OverlapSequences()-- Found overlap with Local_Overlap   begpos=%d endpos=%d length=%d diffs=%d comp=%d/%d\n",
              lo_omesg->begpos, lo_omesg->endpos, lo_omesg->length, lo_omesg->diffs, lo_omesg->comp, flip);
#endif
  }
  
  if (orientation.isBA_AB() || orientation.isBA_BA())
    reverseComplementSequence( seq1, len1 );

  return((dp_omesg) ? dp_omesg : lo_omesg);
}



static
void
ComputeCanonicalOverlap_new(GraphCGW_T *graph, ChunkOverlapCheckT *canOlap) {

  if (consensusA == NULL) {
    consensusA = CreateVA_char(2048);
    consensusB = CreateVA_char(2048);
    qualityA   = CreateVA_char(2048);
    qualityB   = CreateVA_char(2048);
  }

  //  Reset, make it look like there is no overlap.

  canOlap->BContainsA = FALSE;
  canOlap->AContainsB = FALSE;
  canOlap->computed   = TRUE;
  canOlap->overlap    = FALSE;
  canOlap->ahg        = 0;
  canOlap->bhg        = 0;

  //  Save a copy of the spec supplied, then reset.  The copies are made because 'canOlap' is
  //  probably a reference to an overlap in the store.  If we were to modify the spec, we screw up
  //  the hash function, and also cannot even remove the original overlap from the store.

  ChunkOverlapCheckT inOlap = *canOlap;  //  Copy of the original input, will be removed from the store
  ChunkOverlapCheckT nnOlap = *canOlap;  //  Working copy, will be added to the store

  if (nnOlap.maxOverlap < 0)
    //  No point doing the expensive part if there can be no overlap
    return;

  // Get the consensus sequences for both chunks from the ChunkStore
  int32 lengthA = GetConsensus(graph, nnOlap.spec.cidA, consensusA, qualityA);
  int32 lengthB = GetConsensus(graph, nnOlap.spec.cidB, consensusB, qualityB);

  if (nnOlap.minOverlap > (lengthA + lengthB - CGW_DP_MINLEN))
    //  No point doing the expensive part if there can be no overlap
    return;

  // Return value is length of chunk sequence/quality
  // overlap 'em

  char *seq1   = Getchar(consensusA, 0);
  char *seq2   = Getchar(consensusB, 0);

  int32 min_ahang = lengthA - nnOlap.maxOverlap;
  int32 max_ahang = lengthA - nnOlap.minOverlap;

  // tempOlap1 is a static down inside of DP_Compare, don't free it
  ALNoverlap *tempOlap1 = OverlapSequences(seq1, seq2, nnOlap.spec.orientation,
                                           min_ahang, max_ahang,
                                        nnOlap.errorRate,
                                        CGW_DP_THRESH, CGW_DP_MINLEN);

  if (tempOlap1 == NULL)
    //  Didn't find an overlap.  Bail.
    return;

  if (tempOlap1->begpos < 0 && tempOlap1->endpos > 0)
    // ahang is neg and bhang is pos
    nnOlap.BContainsA = TRUE;

  else if (tempOlap1->begpos > 0 && tempOlap1->endpos < 0)
    // ahang is pos and bhang is neg
    nnOlap.AContainsB = TRUE;

  //Print_Overlap_AS(stderr,&AFR,&BFR,O);  Last appears in AS_ALN_qvaligner.c r1.24
  nnOlap.ahg = tempOlap1->begpos;
  nnOlap.bhg = tempOlap1->endpos;

  //  Make the overlap field be the number of bases from the tail of
  //  the A sequence to the beginning of the B sequence
  nnOlap.overlap = tempOlap1->length;

  if  (nnOlap.ahg < 0)
    nnOlap.overlap -= nnOlap.ahg;
  if  (nnOlap.bhg < 0)
    nnOlap.overlap -= nnOlap.bhg;

  // fields are no longer used in DP_Compare (as opposed to DP_Compare_AS)
  nnOlap.quality = 0.0;

  // not dealing with containments here - they can go in either way?  (BPW: ??)

  //  If the overlap found has both negative hangs, we swap the orders.  In the AB_AB case, this
  //  would result in a non-canonical overlap (BA_BA), and we're forced to swap fragments instead.
  //
  //  The suspicious flag is set, indicating the sequences have slid.
  //
  //  Finally, we replace the overlap in the overlap table.
  //
  if (nnOlap.ahg < 0 && nnOlap.bhg < 0) {
    nnOlap.suspicious = TRUE;
    nnOlap.overlap    = tempOlap1->length;

    assert(nnOlap.spec.orientation.isUnknown() == false);
    assert(nnOlap.spec.orientation.isBA_BA()   == false);  //  Not canonical!?

    if        (nnOlap.spec.orientation.isAB_AB()) {
      // we want to switch to a non-canonical orientation ...nnOlap.spec.orientation = BA_BA; but
      // since we can't, we switch order of operands instead
      nnOlap.spec.cidA = inOlap.spec.cidB;
      nnOlap.spec.cidB = inOlap.spec.cidA;
      nnOlap.ahg = -nnOlap.ahg;
      nnOlap.bhg = -nnOlap.bhg;

    } else if (nnOlap.spec.orientation.isAB_BA()) {
      nnOlap.spec.orientation.setIsBA_AB();
      nnOlap.ahg = -tempOlap1->endpos;
      nnOlap.bhg = -tempOlap1->begpos;

    } else if (nnOlap.spec.orientation.isBA_AB()) {
      nnOlap.spec.orientation.setIsAB_BA();
      nnOlap.ahg = -tempOlap1->endpos;
      nnOlap.bhg = -tempOlap1->begpos;
    }

    fprintf(stderr,">>> Fixing up suspicious overlap ("F_CID","F_CID",%c) (ahg:"F_S32" bhg:"F_S32") to ("F_CID","F_CID",%c) (ahg:"F_S32" bhg:"F_S32") len: "F_S32"\n",
            inOlap.spec.cidA, inOlap.spec.cidB, inOlap.spec.orientation.toLetter(), tempOlap1->begpos, tempOlap1->endpos,
            nnOlap.spec.cidA, nnOlap.spec.cidB, nnOlap.spec.orientation.toLetter(), nnOlap.ahg,        nnOlap.bhg,
            nnOlap.overlap);

    DeleteChunkOverlap(ScaffoldGraph->ChunkOverlaps, &inOlap);  //  Delete the old overlap
    DeleteChunkOverlap(ScaffoldGraph->ChunkOverlaps, &nnOlap);  //  New one shouldn't exist, but we'll delete it anyway

    if (InsertChunkOverlap(ScaffoldGraph->ChunkOverlaps, &nnOlap) != HASH_SUCCESS) {
      fprintf(stderr, "ComputeCanonicalOverlap_new()-- Failed to insert.\n");
      assert(0);
    }
  }

  *canOlap = nnOlap;
}



static
int checkChunkOverlapCheckT(ChunkOverlapCheckT *co1,
                            int32 minOverlap,
                            int32 maxOverlap,
                            double errorRate)
{
  if( co1->errorRate != errorRate )
    return FALSE;
  if( minOverlap >= co1->minOverlap && maxOverlap <= co1->maxOverlap &&
      ( minOverlap <= co1->overlap && maxOverlap >= co1->overlap ))
    return TRUE;

  return FALSE;
}



//external
ChunkOverlapCheckT OverlapChunks( GraphCGW_T *graph,
                                  CDS_CID_t cidA, CDS_CID_t cidB,
                                  PairOrient orientation,
                                  int32 minOverlap,
                                  int32 maxOverlap,
                                  double errorRate,
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

  memset(&olap, 0, sizeof(ChunkOverlapCheckT));

  isCanonical = InitCanonicalOverlapSpec(cidA,cidB,orientation,&olap.spec);
  // initalize olap with the chunk IDs and their orientation and record
  // whether the orientation was already canonical or not

  lookup = LookupCanonicalOverlap(ScaffoldGraph->ChunkOverlaps,&olap.spec);
  // lookup whether the overlap had already been computed


  if( lookup != NULL ){
    olap = *lookup;
    if( checkChunkOverlapCheckT(&olap,minOverlap,maxOverlap,errorRate) == FALSE )
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

    ComputeCanonicalOverlap_new(graph, &olap);

    if(insert)
      { // Insert new entries in hashtable, if requested
        //
        int suspicious = olap.suspicious;
        if(olap.suspicious){
          olap.suspicious = FALSE;
          lookup = LookupCanonicalOverlap(ScaffoldGraph->ChunkOverlaps,&olap.spec); // see if it is already in the table
          if(!lookup){
            if(InsertChunkOverlap(ScaffoldGraph->ChunkOverlaps, &olap) != HASH_SUCCESS) {
              fprintf(stderr, "Failed to insert overlap into hash table.\n");
              assert(0);
            }
          }
        }else{
          if(InsertChunkOverlap(ScaffoldGraph->ChunkOverlaps, &olap) != HASH_SUCCESS) {
            fprintf(stderr, "Failed to insert overlap into hash table.\n");
            assert(0);
          }
        }
        lookup = LookupCanonicalOverlap(ScaffoldGraph->ChunkOverlaps,&olap.spec);
        assert(lookup != NULL);
        // ComputeCanonicalOverlap does not return the olap, so we look it up again.

        olap = *lookup;
        olap.suspicious = suspicious;
      }


    if(olap.overlap && insert){ // Insert all non-zero overlaps, if requested
      InsertComputedOverlapEdge(graph, &olap);
    }
    /* Make sure the orientation of the edge we return is IDENTICAL with the spec returned */
    // if the input was not canonical we set the cid's and orientation
    // back to the input value (see als LookupOverlap)
    if( !olap.suspicious  && ! isCanonical ){
      int swap;

      olap.spec.orientation = orientation;

      // If this is non-canonical, swap things back
      olap.spec.cidA = cidA;
      olap.spec.cidB = cidB;
      swap = olap.BContainsA;
      olap.BContainsA = olap.AContainsB;
      olap.AContainsB = swap;
      swap = olap.ahg;
      olap.ahg = olap.bhg;
      olap.bhg = swap;
    }
  }

  if(olap.overlap==0){olap.quality=0;}
  return olap;
}




//external
ALNoverlap* OverlapContigs(NodeCGW_T *contig1, NodeCGW_T *contig2,
                           PairOrient *overlapOrientation,
                           int32 minAhang, int32 maxAhang,
                           int computeAhang,
                           uint32 tryLocal, uint32 tryRev)
{
  ALNoverlap *tempOlap1;
  char *seq1, *seq2;
  double erate, thresh;
  int32 minlen;

static VA_TYPE(char) *consensus1 = NULL;
static VA_TYPE(char) *consensus2 = NULL;
static VA_TYPE(char) *quality1 = NULL;
static VA_TYPE(char) *quality2 = NULL;

  assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
  erate = AS_CGW_ERROR_RATE;
  thresh = CGW_DP_THRESH;
  minlen = CGW_DP_MINLEN;

  // if computeAhang is TRUE, allow a lot of slide in ahang
  if (computeAhang == TRUE)
    {
      minAhang = - (int32) contig2->bpLength.mean;
      maxAhang = (int32) contig1->bpLength.mean;
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
  int len1 = GetConsensus(ScaffoldGraph->ContigGraph, contig1->id, consensus1, quality1) - 1;
  int len2 = GetConsensus(ScaffoldGraph->ContigGraph, contig2->id, consensus2, quality2) - 1;

  seq1 = Getchar(consensus1, 0);
  seq2 = Getchar(consensus2, 0);
  
  // tempOlap1 is a static down inside of DP_Compare
  tempOlap1 =
    OverlapSequences( seq1, seq2,
                      *overlapOrientation, minAhang, maxAhang,
                      erate, thresh, minlen, tryLocal);

  if (tryRev && tempOlap1 == NULL) {
     reverseComplementSequence(seq1, len1);
     reverseComplementSequence(seq2, len2);
     tempOlap1 =
       OverlapSequences( seq1, seq2,
                         *overlapOrientation, 
                         MIN(len1-minAhang-len2, len1-maxAhang-len2), MAX(len1-minAhang-len2, len1-maxAhang-len2),
                         erate, thresh, minlen, tryLocal);
     reverseComplementSequence(seq1, len1);
     reverseComplementSequence(seq2, len2);
  }

  return(tempOlap1);
}












// This is the top level routine that computes all new potential overlaps.
//
#define NUM_SECTIONS (5)

//external
void
ComputeOverlaps(GraphCGW_T *graph, int addEdgeMates, int recomputeCGBOverlaps) {
  HashTable_Iterator_AS iterator;
  uint64 key, value;
  uint32 valuetype;

  InitializeHashTable_Iterator_AS(ScaffoldGraph->ChunkOverlaps->hashTable, &iterator);

  while(NextHashTable_Iterator_AS(&iterator, &key, &value, &valuetype)) {

    //  VERY IMPORTANT.  Do NOT directly use the overlap stored in the hash table.  If we recompute
    //  it (ComputeCanonicalOverlap_new) we can and do screw up the hash table.  This function
    //  occasionally changes the hash key on us, once the key changes, we cannot delete the original
    //  overlap (because we fail to find it in the hash table now) and end up with duplicate entries
    //  in the table.
    //
    ChunkOverlapCheckT olap = *(ChunkOverlapCheckT*)(INTPTR)value;

    assert(key == value);

    if (olap.computed)
      continue;

    if (olap.fromCGB && !recomputeCGBOverlaps)
      continue;

    // set errRate to old value
    assert((0.0 <= AS_CGW_ERROR_RATE) && (AS_CGW_ERROR_RATE <= AS_MAX_ERROR_RATE));
    olap.errorRate = AS_CGW_ERROR_RATE;

    // first we trust that overlap
    olap.suspicious = FALSE;

    if(olap.maxOverlap < 0){ // Dummy!  Who put this overlap in the table?  No overlap is possible.....SAK
      olap.overlap = 0;
      continue;
    }

    ChunkOverlapSpecT inSpec = olap.spec;

    ComputeCanonicalOverlap_new(graph, &olap);

    if (olap.suspicious) {
      int lengthA = GetConsensus(graph, olap.spec.cidA, consensusA, qualityA);
      int lengthB = GetConsensus(graph, olap.spec.cidB, consensusB, qualityB);

      fprintf(stderr,"* CO: SUSPICIOUS Overlap found! Looked for ("F_CID","F_CID",%c)["F_S32","F_S32"] found ("F_CID","F_CID",%c) "F_S32"; contig lengths as found (%d,%d)\n",
              inSpec.cidA,    inSpec.cidB,    inSpec.orientation.toLetter(),    olap.minOverlap, olap.maxOverlap,
              olap.spec.cidA, olap.spec.cidB, olap.spec.orientation.toLetter(), olap.overlap,
              lengthA,lengthB);
    }

    if (addEdgeMates && !olap.fromCGB && olap.overlap)
      InsertComputedOverlapEdge(graph, &olap);
  }
}



//external
void AddUnitigOverlaps(GraphCGW_T *graph,
                       char       *ovlFileName) {

  if (ovlFileName == NULL)
    return;

  fprintf(stderr, "AddUnitigOverlaps()-- reading overlaps from '%s'\n", ovlFileName);

  errno = 0;
  FILE *F = fopen(ovlFileName, "r");
  if (errno) {
    fprintf(stderr, "AddUnitigOverlaps()--   Couldn't open '%s' to read unitig overlaps.\n", ovlFileName);
    exit(1);
    return;
  }

  GenericMesg  *pmesg = NULL;

  while (ReadProtoMesg_AS(F, &pmesg) != EOF) {
    if (pmesg->t == MESG_OVL) {
      FillChunkOverlapWithOVL(graph, (OverlapMesg *)pmesg->m);
    } else {
      fprintf(stderr, "AddUnitigOverlaps()-- unexpected message type '%s'\n", MessageTypeName[pmesg->t]);
    }
  }

  fclose(F);
}

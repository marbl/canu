
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
/* 	$Id: GraphCGW_T.h,v 1.7 2006-03-28 02:51:31 ahalpern Exp $	 */

/**************************************************************************
 *  GraphCGW
 *  
 *  Saul A. Kravitz 8/99
 *
 *  This issgc the proposed data structures and operations for the
 *  scaffold graph that will be used for the final phases of the
 *  Chunk Graph walker.
 *************************************************************************/
#ifndef GRAPH_CGW_H
#define GRAPH_CGW_H

#include "AS_UTL_Var.h"
#include "InputDataTypes_CGW.h"
#include "ChunkOverlap_CGW.h"
//#include "MultiAlignStore_CNS.h"
#include "AS_ALN_aligners.h"
#include "dpc_CNS.h"

typedef enum {
        INVALID_EDGE_STATUS = 0,
        UNKNOWN_EDGE_STATUS = 1,
        UNTRUSTED_EDGE_STATUS = 2,
	TENTATIVE_UNTRUSTED_EDGE_STATUS = 4,
	TENTATIVE_TRUSTED_EDGE_STATUS = 8,
        TRUSTED_EDGE_STATUS = 16,
	LARGE_VARIANCE_EDGE_STATUS = 32,
	INTER_SCAFFOLD_EDGE_STATUS = 64
}EdgeStatus;

EdgeStatus AS_CGW_SafeConvert_uintToEdgeStatus(unsigned int input);


typedef struct {
	CDS_CID_t idA;
	CDS_CID_t idB;
	ChunkOrientationType orient; // Orientation of idA <-> idB, BE CAREFUL IF YOU WANT idB<->idA
	//
	int32 edgesContributing;
        float quality;   // Used to order edges by decreasing quality (quality = 0.0 is the best, 1.0 is the worst)

        union{
	  int64 all;  // Used for overlap, repeatoverlap, tandemOverlap, etc
	  struct {
	    unsigned int isInferred:1; // Is this an inferred edge.
	    unsigned int isTentative:1; // Is this inferred edge finalized.
	    unsigned int isLeastSquares:1; // High confidence edge based on
	    // least squares calculation.
	    unsigned int isEssential:1; // Is this edge essential.
	    unsigned int wasEssential:1; // wass this edge essential before.
	    unsigned int isActive:1; // Is this edge currently of interest.
	    unsigned int isConfirmed:1; // Is this edge confirmed.
	    unsigned int isContigConfirming:1; /* This edge indicates to the contigger that the two
						elements should be merged into a contig */
	    unsigned int isUniquetoUnique:1; /* Is this edge from one unique CI
						to another unique CI. */
	    unsigned int isTransitivelyRemoved:1; /* Has this edge been transitively removed */
	    unsigned int isInferredRemoved:1;  /*Has this edge been
						 transitively removed by an inferred edge. */
	    unsigned int isRedundantRemoved:1; /*Has this edge been
						 removed by another edge between the sam pair of unique CIs. */
	    unsigned int isDeleted:1; // Is edge deleted.
	    unsigned int isPossibleChimera:1;         /* An edgemate consisting of a single raw edgemate
							 and an overlap, where the same read participates 
							 both the mate and overlap relationship */

	    /* Flags relating to the type of relationship that induced this edge */

	    unsigned int inducedByUnknownOrientation:1; /* One of 4 raw edges induced by a LKG with AS_UNKNOWN orientation */
	    unsigned int hasContributingOverlap:1;  /* EdgeMate includes contribution from an overlap, not repeat only */
	    unsigned int hasRepeatOverlap:1;        /* Has overlaps beyond branch points ==> repeat only */
	    unsigned int hasTandemOverlap:1;  /* Overlapper reported a min-max range of overlaps */
	    unsigned int aContainsB:1;        /* From CGB for multiply contained fragments */
	    unsigned int bContainsA:1;        /* From CGB for multiply contained fragments */
	    unsigned int mustOverlap:1;       /* Marked on merged edges when mate-link variance is signficantly less than overlap length */
	    unsigned int hasGuide:1;                /* Contains one or more guide edges */
	    unsigned int hasSTSGuide:1;             /* Contains an STS Guide */
	    unsigned int hasMayJoin:1;             /* Contains a may join constraint */
	    unsigned int hasMustJoin:1;             /* Contains a must join constraint */
	    unsigned int hasTransChunk:1;           /* was a transitively removed edge in cgb */
	    unsigned int hasContainmentOverlap:1;  /* Implies a containment */
	    unsigned int isRaw:1;                   /* True for raw edges, false for merged edges */

	    unsigned int hasExtremalAFrag:1; /* Used in merging raw edge mates, meaningless elsewhere */
	    unsigned int hasExtremalBFrag:1; /* Used in merging raw edge mates, meaningless elsewhere */
            unsigned int rangeTruncated:1;    /* TRUE if we looked for an overlap and didn't find any,
						 and thus truncated the range of distances */
	    unsigned int inAssembly:1;
	    unsigned int isBogus:1;  // determined from simulator annotations, overloaded by scaffold merging code (see CIScaffold_Merge_CGW.c)
	    unsigned int isProbablyBogus:1; // determined from distance and distIndex empirically, overloaded for one-sided edges in CIScaffold_Merge
	    unsigned int hasConfirmingPath:1; // Has this edge another path that confirms its length & var (used in GapWalkerREZ.c)
	    

	    unsigned int edgeStatus:7; 
	    
	    unsigned int isMarkedForDeletion:1;   // We plan to delete this guy
	    unsigned int MeanChangedByWalking:1;
            unsigned int highQualityA:1;           // One of the top-ranked edges incident on node idA
            unsigned int highQualityB:1;           // One of the top-ranked edges incident on node idB
	    unsigned int isSloppy:1;
	    unsigned int isBridge:1;               // Bridge edge in a scaffold

	  }bits;
	}flags;
	//	
	LengthT distance; // gap/overlap length

  /* vvvvv Iterator does not fill in fields below this line vvvvvv */
	//
	CDS_CID_t nextALink; // next edge involving cidA, -1 if none
	CDS_CID_t nextBLink; // next edge involving cidB, -1 if none
	CDS_CID_t prevALink; // prev edge involving cidA, -1 if none
	CDS_CID_t prevBLink; // prev edge involving cidB, -1 if none

        float32 minDistance;  /* negative implies potential overlap 
			       * This Field is overloaded to store the distance.mean when
			       * the flag MeanChangedByWalk is true.  The function RestoreEdgeMeans
			       * restores the value and unsets the flag.
			       */

  /*** We need these to reference back to the fragment pair that induced the edge */
      CDS_CID_t fragA;    // The fragment in chunk/contigA, or NULLINDEX
      CDS_CID_t fragB;    // The fragment in chunk/contigB or NULLINDEX
      CDS_CID_t distIndex; // Index of underlying distance record or NULLINDEX
      CDS_CID_t nextRawEdge; // index to next raw edge in the chain.  These are maintained in a singly linked list.
      CDS_CID_t topLevelEdge; /* If this is a raw edge, references the 'owner' or top-level edge in which the raw
			     edge is linked */
      CDS_CID_t referenceEdge;  /*** Reference to inducing edge */
      
}EdgeCGW_T;

typedef struct
{
  CDS_COORD_t samples;
  CDS_CID_t frags;
  CDS_CID_t mates;
} MateInfoT;

typedef enum {
	DISCRIMINATORUNIQUECHUNK_CGW,  // Initially, all scaffolded chunks are these
	UNRESOLVEDCHUNK_CGW,     // Not discriminator unique
	UNIQUECHUNK_CGW,               // These are unique chunks that were not discriminator unique
	RESOLVEDREPEATCHUNK_CGW,        // A subset of a repeat chunk that has been instantiated
	//
	CONTIG_CGW,                     // A contig that has subsumed >=1 chunk
	UNIQUECONTIG_CGW,
	RESOLVEDCONTIG_CGW,
	UNRESOLVEDCONTIG_CGW,
	//
	REAL_SCAFFOLD,     // the genuine article
	OUTPUT_SCAFFOLD,    // an artefact generated for output
	SCRATCH_SCAFFOLD    // a temporary scaffold 
} ChunkInstanceType;


/* This enum is used to encode the classification of unitigs based on the
   simulator coordinates of the fragments */
typedef enum {
  UU_CGBTYPE = 0,
  UR_CGBTYPE = 1,
  RU_CGBTYPE = 2,
  RR_CGBTYPE = 3,
  XX_CGBTYPE = 4
}CGB_Type;


/* This enum is used to encode the presence/absence of tandem repeat overlaps
   on the end of the ChunkInstance (only meaningful if isCI = TRUE)
*/
typedef enum {
  NO_TANDEM_OVERLAP =       (unsigned)0x0,
  AEND_TANDEM_OVERLAP =     (unsigned)0x1,
  BEND_TANDEM_OVERLAP =     (unsigned)0x2,
  BOTH_END_TANDEM_OVERLAP = (unsigned)0x3
}TandemOverlapType;


/* ChunkInstanceT:
   This structure is overloaded to hold ChunkInstances, Contigs and Scaffolds.
   The enum ChunkInstanceType and the flag bits isCI, isContig, and isScaffold
   should be maintained consistently as follows:
   isCI = TRUE (isContig = FALSE, isScaffold = FALSE)
	DISCRIMINATORUNIQUECHUNK_CGW,  // Initially, all scaffolded chunks are these
	UNRESOLVEDCHUNK_CGW,     // Not discriminator unique
	UNIQUECHUNK_CGW,               // These are unique chunks that were not discriminator unique
	RESOLVEDREPEATCHUNK_CGW,        // A subset of a repeat chunk that has been instantiated
   isContig = TRUE (isCI = FALSE, isScaffold = FALSE)
	CONTIG_CGW,                     // A contig that has subsumed >=1 chunk
	UNIQUECONTIG_CGW,
	RESOLVEDCONTIG_CGW,
	UNRESOLVEDCONTIG_CGW,
   isScaffold = TRUE (isCI = FALSE, isContig = FALSE)
	REAL_SCAFFOLD,     // the genuine article
	OUTPUT_SCAFFOLD,    // an artefact generated for output
	SCRATCH_SCAFFOLD    // a temporary scaffold 
*/	

        
typedef struct{
  ChunkInstanceType type; //
  CDS_CID_t id;        // Instance ID
  CDS_CID_t outputID;  // InstanceID (dense encoding over LIVE instances)
  CDS_CID_t scaffoldID; // scaffold ID
  CDS_CID_t prevScaffoldID; // previous scaffold ID from last iteration
  int32     indexInScaffold; // Relative position from A end of Scaffold (not kept current)
  CDS_CID_t smoothExpectedCID; // Used to avoid cycles in smoothing the trasitively
  // reduced unique unitig graph.

  int32 numEssentialA; // Number of essential edges off the A end.
  int32 numEssentialB; // Number of essential edges off the B end.
  CDS_CID_t essentialEdgeA; // Essential edge off the A end.  (also used for free list maintenance)
  CDS_CID_t essentialEdgeB; // Essential edge off the B end.
  //
  // Chunk of Scaffolds / Contigs
  // If this chunk instance is in a contig, the links and positions are to contig neighbors
  // If this chunk is not in a contig, the links and positions are within a scaffold
  CDS_CID_t AEndNext;  // index to predecessor on A end of Chunk/Contig in Scaffold
  CDS_CID_t BEndNext;  // index to predecessor on B end of Chunk/Contig in Scaffold
  //
  LengthT  bpLength;
  LengthT  offsetAEnd;     // Offset of A end of CI relative to A end of Contig/CIScaffold
  LengthT  offsetBEnd;     // Offset of B end of CI relative to A end of Contig/CIScaffold

  /***** Simulator Coordinates ****/
  CDS_COORD_t aEndCoord, bEndCoord;

  union{  // ChunkInstanceType discriminates
    struct CIINFO_TAG {
      CDS_CID_t contigID;   // contigID -- if -1, this chunkInstance not merged into a contig
      //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/
      // Both of these are redundant, since we can get
      // to the same info by looking in the multiAlignT for this CI
      // at some point we should nuke them
      CDS_CID_t headOfFragments; /* Index of first FragInfoT record belonging to this  chunkInstance
				These will be linked together.  */
      int32 numFragments;
      /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
      CDS_COORD_t  branchPointA;
      CDS_COORD_t  branchPointB;    
      int32        coverageStat;
      CDS_CID_t    baseID;    /* If this is a  RESOLVEDREPEAT, the id of the original 
			     CI from which it was spawned */
      int32 numInstances; /* Number of actual or surrogate instances in scaffolds
			     If this is not a RESOLVEDREPEAT, numInstances should be = 0 */
      union{
	struct {            // if numInstances is <=2
	  CDS_CID_t instance1;
	  CDS_CID_t instance2;
	}in_line;
	VA_TYPE(CDS_CID_t) *va; // if numInstances is > 2
      }instances;
#ifdef DEBUG_DATA
      CDS_CID_t source;
#endif
    }CI;
    struct CONTIGINFO_TAG{
      // All of this is redundant to what is incldued in the multiAlignment
      CDS_CID_t AEndCI;   //Index of Chunk Instance at A end of Contig
      CDS_CID_t BEndCI;   //Index of Chunk Instance at B end of Contig
      cds_int32 numCI;    //Number of CI in contig
      CDS_COORD_t branchPointA;
      CDS_COORD_t branchPointB;    
    }Contig;
    struct CISCAFFOLD_TAG{
        CDS_CID_t AEndCI; // Index of Chunk Instance at A end of Scaffold
	CDS_CID_t BEndCI; // Index of Chunk Instance at B End of Scaffold
        cds_int32 numElements; // If containsCIs, these are CIs, else they are Contigs
      /*** Info computed by RecomputeScaffoldPositions ***/
      float32 leastSquareError;   // Measure of chi-squared of computed positions
      cds_int32 numLeastSquareClones; // Relates to degrees of freedom for chi-square calculation
      /*** Info computed by MarkInternaledges ***/
      cds_int32 internalEdges;  // Number of merged edges (not including UNTRUSTED) that are internal to scaffold
      cds_int32 confirmedInternalEdges; // Number of merged edges confirmed by current CI positions
    }Scaffold;
  }info;

  //
  union{
    struct {
      unsigned int isUnique:1;
      unsigned int smoothSeenAlready:1; // Used for cycle detection
      /* This field has one of the following values (first refers to essentials, 2nd to contains):
	 XX_CGBTYPE   -- No annotation on unitig in cgb output
	 RU_CGBTYPE   -- CGB reported "ru"
	 UR_CGBTYPE   -- CGB reported "ur"
	 RR_CGBTYPE   -- CGB reported "rr"
	 UU_CGBTYPE   -- CGB reported "uu"
      */
      unsigned int cgbType:3; 
      unsigned int isDead:1;      
      unsigned int isFree:1;      
      unsigned int containsCIs:1;        // Scaffold contains either CIs or Contigs

      /* This field has one of the following values:
	 NO_TANDEM_OVERLAP
	 AEND_TANDEM_OVERLAP
	 BEND_TANDEM_OVERLAP
	 BOTH_END_TANDEM_OVERLAP = AEND_TANDEM_OVERLAP | BEND_TANDEM_OVERLAP
      */
      unsigned int tandemOverlaps:2;
      unsigned int isCI:1;
      unsigned int isContig:1;
      unsigned int isScaffold:1;
      /* The following is TRUE for CIs that are surrogate CIs and
	 for Contigs that contain a single surrogate CI.  Set at creation
	 time in the Split functions */
      unsigned int isSurrogate:1; 
      unsigned int beingContigged:1;
      unsigned int includesFinishedBacFragments:1;
      /* The following is used by gap walking to see which scaffolds have been walked.
	 Initialized to zero when a scaffold is created, managed by gw afterward */
      unsigned int walkedAlready:1;
      unsigned int walkedTooShort:1;
      unsigned int walkedTooLong:1;
      unsigned int walkMaxedOut:1;
      unsigned int walkedTrivial:1;
      unsigned int isStoneSurrogate:1;
      unsigned int isWalkSurrogate:1;
      unsigned int failedToContig:1;
      unsigned int isChaff:1;  // a singleton unitig/contig that is not placed or used as a surrogate
      unsigned int isMisplaced:1; // Contig placed out of order wrt simulator coordinates
      unsigned int isStone:1;
      unsigned int isWalk:1;
      unsigned int isRock:1;
      unsigned int isPotentialRock:1; // has at least 2 celera reads with external mates
      unsigned int isPotentialStone:1; // has at least 1 celera reads with external mates
    }bits;
    int32 all;
  }flags;

  CDS_CID_t edgeHead;  // Pointer to linked list of edges  in edges;

  float32 microhetScore; /* Score from Knut&Aaron's microhet detecter, valid for CIs only! Could be in union above, but CI
			  variant is biggest one.  So I reused an unused field at the top level SAK */


  CDS_CID_t setID;

}NodeCGW_T;

VA_DEF(NodeCGW_T)
VA_DEF(EdgeCGW_T)


/* GraphCGW_T
   This is the basis structure for holding CGW graphs, of which there are 3:
      Graph of CIs
      Graph of Contigs
      Graph of Scaffolds
*/

typedef enum{
  CI_GRAPH =       'c',
  CONTIG_GRAPH =   'C',
  SCAFFOLD_GRAPH = 'S'
}GraphType;


typedef struct{
  GraphType type;
  VA_TYPE(NodeCGW_T) *nodes;
  VA_TYPE(EdgeCGW_T) *edges;
  cds_int32 numActiveNodes;
  cds_int32 numActiveEdges;
  CDS_CID_t freeEdgeHead;
  CDS_CID_t tobeFreeEdgeHead; // staging area for edges waiting to be moved to the free list
  CDS_CID_t freeNodeHead;
  CDS_CID_t tobeFreeNodeHead; // staging area for nodes waiting to be moved to the free list
  CDS_CID_t deadNodeHead;
  ChunkOverlapperT *overlapper;
  // Storage for consensus
  //  MultiAlignStoreT *maStore;
}GraphCGW_T;


typedef struct{
  CDS_CID_t id;
  CDS_COORD_t minGap;
  CDS_COORD_t maxGap;
  ChunkOrientationType orient;
}RevivedEdgeT;

static void InitializeGraph(GraphCGW_T *graph){
  graph->nodes = NULL;
  graph->edges = NULL;
  graph->numActiveNodes = 0;
  graph->numActiveEdges = 0;
  graph->freeEdgeHead = NULLINDEX;
  graph->tobeFreeEdgeHead = NULLINDEX;
  graph->freeNodeHead = NULLINDEX;
  graph->tobeFreeNodeHead = NULLINDEX;
  graph->deadNodeHead = NULLINDEX;
  graph->overlapper = NULL;
}



static int32 GetNumGraphNodes(GraphCGW_T *graph){
  return (int32) GetNumNodeCGW_Ts(graph->nodes);
}
static int32 GetNumGraphEdges(GraphCGW_T *graph){
  return (int32) GetNumEdgeCGW_Ts(graph->edges);
}
/* Accessors */
static EdgeCGW_T *GetGraphEdge(GraphCGW_T *graph, CDS_CID_t edgeID){
  return GetEdgeCGW_T(graph->edges, edgeID);
}
static NodeCGW_T *GetGraphNode(GraphCGW_T *graph, CDS_CID_t nodeID){
  return GetNodeCGW_T(graph->nodes, nodeID);
}

/* Append */
static void AppendGraphNode(GraphCGW_T *graph, NodeCGW_T *node){
  AppendNodeCGW_T(graph->nodes, node);
}
static void AppendGraphEdge(GraphCGW_T *graph, EdgeCGW_T *edge){
  AppendEdgeCGW_T(graph->edges, edge);
}


/* Create New */
static NodeCGW_T *CreateNewGraphNode(GraphCGW_T *graph){
  NodeCGW_T node;
  node.id = GetNumGraphNodes(graph);
  node.flags.all = 0;
  node.flags.bits.cgbType = XX_CGBTYPE;
  node.edgeHead = NULLINDEX;
  node.setID = 0;

  switch(graph->type){
  case CI_GRAPH:
    node.type = UNRESOLVEDCHUNK_CGW;
    node.flags.bits.isCI = TRUE;
    node.info.CI.contigID = NULLINDEX;
    node.info.CI.headOfFragments = NULLINDEX;
    node.info.CI.numFragments = 0;
    node.info.CI.branchPointA = node.info.CI.branchPointB = 0;
    node.microhetScore = -1.0;
    node.info.CI.coverageStat = 0;
    node.info.CI.numInstances = 0;
    node.info.CI.instances.va = NULL;
    break;
  case CONTIG_GRAPH:
    node.type = CONTIG_CGW;
    node.flags.bits.isContig = TRUE;
    node.info.Contig.AEndCI = node.info.Contig.BEndCI = NULLINDEX;
    node.info.Contig.numCI = 0;
    break;
  case SCAFFOLD_GRAPH:
    node.type = REAL_SCAFFOLD;
    node.flags.bits.isScaffold = TRUE;
    node.info.Scaffold.AEndCI =     node.info.Scaffold.BEndCI = NULLINDEX;
    node.info.Scaffold.numElements = 0;
    break;
  default:
    assert(0);
  }
  node.scaffoldID = NULLINDEX;
  node.smoothExpectedCID = NULLINDEX;
  node.numEssentialA = node.numEssentialB = 0;
  node.essentialEdgeA = node.essentialEdgeB = NULLINDEX;
  node.indexInScaffold = NULLINDEX;
  node.microhetScore = -1.0;
  node.edgeHead = NULLINDEX;
  node.smoothExpectedCID = 0;
  node.AEndNext = node.BEndNext = NULLINDEX;
  node.bpLength.mean = node.bpLength.variance = 0.0;
  node.offsetAEnd.mean = node.offsetAEnd.variance = 0.0;
  node.offsetBEnd.mean = node.offsetBEnd.variance = 0.0;
  node.aEndCoord = node.bEndCoord = NULLINDEX;

  AppendNodeCGW_T(graph->nodes, &node);
  return(GetGraphNode(graph, node.id));
}
static EdgeCGW_T *CreateNewGraphEdge(GraphCGW_T *graph ){
  EdgeCGW_T edge;
  AppendEdgeCGW_T(graph->edges, &edge);
  return(GetGraphEdge(graph, GetNumGraphEdges(graph) - 1));
}

/* Constructor and Destructor */
GraphCGW_T *CreateGraphCGW(GraphType type, cds_int32 numNodes,
                           cds_int32 numEdges);
void DeleteGraphCGW(GraphCGW_T *graph);

/* Persistence */
void SaveGraphCGWToStream(GraphCGW_T *graph, FILE *stream);
GraphCGW_T *LoadGraphCGWFromStream(FILE *stream);


static EdgeStatus GetEdgeStatus(EdgeCGW_T *edge){
  return (EdgeStatus) edge->flags.bits.edgeStatus;
}


#if 0
/*
Create a new node by extracting fragments from  an existing node.
The node must be either a contig or a CI.
 
The fragsToExtract argument is a list of fragments, all of which must belong
to the node, that are to be extracted and become the basis for the new node.
The multi-alignment of the new node will include these fragments,as well as
all unitigs that are present in the base node's multi-alignment.
 
Following the split, the new node will have copies of all raw overlap edges
incident on the base node, and the raw link edges will be distributed between
the two as a function of the split of the fragments.  Re-merging of edges
incident on both nodes completes the operation.
*/
CDS_CID_t CreateNewNodeByExtractingEdges(GraphCGW_T *graph,
                                         CDS_CID_t id,
                                         VA_TYPE(CDS_CID_t) *fragsToExtract);
#endif

static void SetEdgeStatus(GraphCGW_T *graph,
                          EdgeCGW_T *edge,
                          EdgeStatus status){
  edge->flags.bits.edgeStatus = (unsigned int) status;
  if(edge->flags.bits.isRaw)
    return;
  // Propagate to raw edges that are attached
  while(NULL != (edge = GetGraphEdge(graph, edge->nextRawEdge))){
    edge->flags.bits.edgeStatus = (unsigned int) status;
  }
}


/*  PropagateEdgeStatusToFrag
    Iterate through all raw edges constituting this edge, and mark
    their fragments mate status
*/
void PropagateEdgeStatusToFrag(GraphCGW_T *graph, EdgeCGW_T *edge);


static ChunkOrientationType InvertEdgeOrient
(const ChunkOrientationType orient){
  ChunkOrientationType new_orient = orient;
  switch(orient){
  case AB_BA:
    new_orient = BA_AB; break;
  case BA_AB:
    new_orient = AB_BA; break;
  case XX_XX:
    new_orient = orient; break;
  case AB_AB:
    new_orient = BA_BA; break;
  case BA_BA:
    new_orient = AB_AB; break;
  default:
    assert(0);
  }
  return new_orient;
}


static ChunkOrientationType FlipEdgeOrient(ChunkOrientationType orient){
  ChunkOrientationType new_orient = orient;
  switch(orient){
  case AB_BA:
  case BA_AB:
  case XX_XX:
    new_orient = orient; break;
  case AB_AB:
    new_orient = BA_BA; break;
  case BA_BA:
    new_orient = AB_AB; break;
  default:
    assert(0);
  }
  return new_orient;
}


// Given orientation relating X and Y as ori(X)_ori(Y), return ori(Y)_ori(X)
static ChunkOrientationType EdgeOrientSwap
(const ChunkOrientationType orient){
  ChunkOrientationType new_orient = orient;
  switch(orient){
  case AB_BA:
    new_orient = BA_AB; break;
  case BA_AB:
    new_orient = AB_BA; break;
  case XX_XX:
    new_orient = orient; break;
  case AB_AB:
    new_orient = AB_AB; break;
  case BA_BA:
    new_orient = BA_BA; break;
  default:
    assert(0);
  }
  return new_orient;
}


static NodeOrient FlipNodeOrient(NodeOrient orient){
  NodeOrient new_orient = orient;
  switch(orient){
  case A_B:
    new_orient = B_A; break;
  case B_A:
    new_orient = A_B; break;
  case X_X:
    new_orient = orient; break;
  default:
    assert(0);
  }
  return new_orient;
}


static int32 IsSurrogateNode(NodeCGW_T *node){
  return node->flags.bits.isSurrogate;
}


static NodeOrient GetNodeOrient(NodeCGW_T *CI){
  if(CI->offsetBEnd.mean > CI->offsetAEnd.mean){
    return A_B;
  }
  //else
  return B_A;
}


static ChunkOrientationType GetEdgeOrientationWRT(EdgeCGW_T* edge,
                                                  CDS_CID_t wrtCI){
  AssertPtr(edge);
  if(edge->idA == wrtCI)
    return edge->orient;
  //else
  return FlipEdgeOrient(edge->orient);
}

// which side of the seed is a chunk on in the scaffold
typedef enum {LEFT_SIDE, RIGHT_SIDE, THE_SEED} SeedSideType;

static SeedSideType GetChunkSeedSide( ChunkOrientationType edgeOrient){
  SeedSideType ret = RIGHT_SIDE;
  switch(edgeOrient){
  case AB_AB:
  case AB_BA:
    ret = RIGHT_SIDE; break;
  case BA_BA: 
  case BA_AB: 
    ret = LEFT_SIDE; break;
  default:
    assert(0);
    break;
  }
  return ret;
}

/* Always consider the seed in the AB orientation, so we flip
   the other chunk accordingly to make it so */
static ChunkOrient GetRelativeChunkOrientation( ChunkOrientationType edgeOrient){
  ChunkOrient ret = A_B;
  switch(edgeOrient){
  case AB_AB:
  case BA_BA:  // flip so that it is really AB_AB with the seed on the right
    ret = A_B;
    break;
  case AB_BA:
  case BA_AB: // flip so that it is really BA_AB with the seed on the left
    ret = B_A;
    break;
  default:
    assert(0);
    break;
  }
  return ret;
}

static ChunkOrientationType GetChunkPairOrientation(ChunkOrient orientA,
                                                    ChunkOrient orientB){
  ChunkOrientationType ret = BA_AB;
  int code = 0;

  if(orientA == A_B)
    code += 2;

  if(orientB == A_B)
    code += 1;

  switch(code){
    case 0: // BA_BA
      ret = BA_BA; break;
    case 1: // BA_AB
      ret = BA_AB; break;
    case 2: // AB_BA
      ret = AB_BA; break;
    case 3: // BA_AB
      ret = AB_AB; break;
    default:
      assert(0);
  }
  return ret;
}





static unsigned int GetNodeTandemOverlaps(NodeCGW_T *ci){
  return ((unsigned int)ci->flags.bits.tandemOverlaps);
}

static ChunkInstanceType GetNodeType(NodeCGW_T *ci){
  return ci->type;
}


// Reset to initial, unscaffolded value
static void ResetNodeType(NodeCGW_T *ci){
  if(ci->flags.bits.isCI){
    switch(ci->type){
    case DISCRIMINATORUNIQUECHUNK_CGW:
      ci->flags.bits.isUnique = 1;
      break;
    case UNIQUECHUNK_CGW:
      ci->type = UNRESOLVEDCHUNK_CGW;
      ci->flags.bits.isUnique = 0;
      break;
    case UNRESOLVEDCHUNK_CGW:
    case RESOLVEDREPEATCHUNK_CGW:
      assert(0); // Shouldn't happen
      break;
    default:
      assert(0);
      break;
    }
  }else if(ci->flags.bits.isContig){
    assert(ci->type == CONTIG_CGW);
    assert(0);     // don't know how to handle this yet
  }else if(ci->flags.bits.isScaffold){
    switch(ci->type){
    case REAL_SCAFFOLD:     // the genuine article
    case OUTPUT_SCAFFOLD:    // an artefact generated for output
    case SCRATCH_SCAFFOLD:    // a temporary scaffold 
    default:
      assert(0); // shouldn't be here
      break;
    }
  }else{
    assert(0);
  }
}


static void SetNodeType(NodeCGW_T *ci, ChunkInstanceType type){
  if(ci->flags.bits.isCI){
    switch(type){
    case DISCRIMINATORUNIQUECHUNK_CGW:
    case UNIQUECHUNK_CGW:
      ci->flags.bits.isUnique = 1;
      break;
    case UNRESOLVEDCHUNK_CGW:
    case RESOLVEDREPEATCHUNK_CGW:
      ci->flags.bits.isUnique = 0;
      break;
    default:
      assert(0);
      break;
    }
  }else if(ci->flags.bits.isContig){
    assert(type == CONTIG_CGW);
  }else if(ci->flags.bits.isScaffold){
    switch(type){
    case REAL_SCAFFOLD:     // the genuine article
    case OUTPUT_SCAFFOLD:    // an artefact generated for output
    case SCRATCH_SCAFFOLD:    // a temporary scaffold 
      break;
    default:
      assert(0);
    }
  }else{
    assert(0);
  }
  ci->type = type;
}

/* isBacOnlyEdge:  true if all of the non-overlap edges are from BACs */
int isBacOnlyEdge(GraphCGW_T *graph, EdgeCGW_T *edge);

/* EdgeDegree: */
void EdgeDegree(GraphCGW_T *graph, EdgeCGW_T *edge,
                int32 *totalDegree, int32 *noBacDegree);

static int isContainmentEdge(EdgeCGW_T *edge){
  return edge->flags.bits.hasContainmentOverlap;
}

static int isMustOverlapEdge(EdgeCGW_T *edge){
  return edge->flags.bits.mustOverlap;
}

static int isInferredEdge(EdgeCGW_T *edge){
  return edge->flags.bits.isInferred;
}

static int isOverlapEdge(EdgeCGW_T *edge){
  return (edge->flags.bits.hasContributingOverlap || 
    edge->flags.bits.hasRepeatOverlap || 
    edge->flags.bits.hasTandemOverlap ||
    edge->flags.bits.aContainsB ||
    edge->flags.bits.bContainsA);
}

static int edgeContainsCI(EdgeCGW_T *edge, CDS_CID_t id){
  assert((id == edge->idA) || (id == edge->idB));
  return((id == edge->idA) ? (int)edge->flags.bits.bContainsA :
	 (int)edge->flags.bits.aContainsB);
}

static int isTransChunkEdge(EdgeCGW_T *edge){
  return (edge->flags.bits.hasTransChunk);
}


// a 10k mate has std of 1.0-1.5k
// a 50k mate has std of 5.0k
// a BE has std of 10+k
//
// for test run on mosquito to treat 10ks libraries as guides
// #define SLOPPY_EDGE_VARIANCE_THRESHHOLD (38025)    // sigma = 195
//
// for TVG run to use fosmids in initial scaffolding  ALD
// #define SLOPPY_EDGE_VARIANCE_THRESHHOLD (14.0e+6)    // sigma = 3700
//
#define SLOPPY_EDGE_VARIANCE_THRESHHOLD (4.0e+6)    // sigma = 2000

static int isSloppyEdge(EdgeCGW_T *edge){
  if(edge->flags.bits.isInferred){
    return 0;
  }
  if(edge->flags.bits.hasGuide || 
    edge->flags.bits.hasSTSGuide || 
      edge->flags.bits.isSloppy){
     return 1;
  }
     
  // for older files
  return (edge->distance.variance > SLOPPY_EDGE_VARIANCE_THRESHHOLD);  
}

static int isSingletonOverlapEdge(EdgeCGW_T *edge){
  return ( isOverlapEdge(edge) &&   edge->edgesContributing <= 1);
}

static int isProbablyBogusEdge(EdgeCGW_T *edge){
  return edge->flags.bits.isProbablyBogus;
}

static int isConfirmedEdge(EdgeCGW_T *edge){
#if 0  // If we want to ignore repeat overlaps, use this
  if(edge->edgesContributing > 2)
    return TRUE;
  if(edge->edgesContributing == 2){
    if( !isOverlapEdgeMate(edge) )
      return TRUE;
    else if(edge->hasContributingOverlap) // this means it does NOT have a repeat or tandem overlap
      return TRUE;
  }
#endif
  return(edge->edgesContributing >=2);
}


static void setEssentialEdgeStatus(EdgeCGW_T *edge, int status){
  if(status){
    assert(!isSloppyEdge(edge));
  }
  edge->flags.bits.isEssential = status;
}

static int getEssentialEdgeStatus(EdgeCGW_T *edge){
  return edge->flags.bits.isEssential;
}

/***************************************************************************/
/* Node Iterator */

typedef struct{
  // indices of nodes
  CDS_CID_t next;
  CDS_CID_t curr;
  CDS_CID_t prev;
  //
  int uniqueOnly;
  int verbose;
  GraphCGW_T *graph;
}GraphNodeIterator;

#define GRAPH_NODE_DEFAULT 0
#define GRAPH_NODE_UNIQUE_ONLY 1
#define GRAPH_NODE_VERBOSE 2

static void InitGraphNodeIterator(GraphNodeIterator *iterator,
				  GraphCGW_T *graph,
				  int flags){
  iterator->curr = iterator->prev = NULLINDEX;
  iterator->next = 0;
  iterator->graph = graph;
  iterator->uniqueOnly = flags &  GRAPH_NODE_UNIQUE_ONLY;
  iterator->verbose = flags & GRAPH_NODE_VERBOSE;
}
				  
static NodeCGW_T *NextGraphNodeIterator(GraphNodeIterator *e){
  NodeCGW_T *retNode = NULL;
  
  if(e->verbose)
    fprintf(stderr,
            "* NextGraphNodeIterator prev:" F_CID " curr:" F_CID " next:" F_CID "\n",
	    e->prev, e->curr, e->next);

  if(e->next == NULLINDEX){
    if(e->curr != NULLINDEX){ // do this once
      e->prev = e->curr;
      e->curr = e->next;
    }
    if(e->verbose)
      fprintf(stderr,"* Fell off the end ---next = NULLINDEX\n");
    return retNode;
  }


  while(retNode == (NodeCGW_T *)NULL &&
	(e->next != NULLINDEX)){
    NodeCGW_T *node = GetGraphNode(e->graph, e->next);
    int isUniqueEnough = TRUE;
    int isInitialized = FALSE;
    if(!node){
      if(e->verbose)
	fprintf(stderr,"* Fell off the end at index " F_CID "\n", e->next);
      break;
    }
    {
      CDS_CID_t idx = GetVAIndex_NodeCGW_T(e->graph->nodes, node);
      isInitialized = node->flags.bits.isDead || node->flags.bits.isFree ||
        (node->id == idx);
    }
    if(isInitialized){
      if(e->uniqueOnly){
	if(node->flags.bits.isContig){
	  isUniqueEnough = (node->scaffoldID > NULLINDEX ||
                            node->flags.bits.isUnique);
	  if(e->verbose & !isUniqueEnough)
	    fprintf(stderr,
                    "* Skipping contig " F_CID " since not scaffolded\n",
		    node->id);
	}else if(node->flags.bits.isScaffold){
	  isUniqueEnough = TRUE;
	}else{
	  assert(node->flags.bits.isCI);
	  isUniqueEnough = (node->type == DISCRIMINATORUNIQUECHUNK_CGW ||
                            node->type == UNIQUECHUNK_CGW ||
                            node->flags.bits.isUnique);
	  if(e->verbose & !isUniqueEnough)
	    fprintf(stderr,"* Skipping CI " F_CID " since not scaffolded\n",
		    node->id);
	}
      }
    }
    if(isInitialized &&
       (!node->flags.bits.isDead) &&
       (!node->flags.bits.isFree) && 
        isUniqueEnough ){
      retNode = node;
      if(e->verbose)
	fprintf(stderr,"* Returning node " F_CID "\n", node->id);
    }
    e->prev = e->curr;
    e->curr = e->next;
    e->next++;
  }
  return retNode;
}

/***************************************************************************/
/* EDGE ITERATOR */
typedef struct {
  // Indices of MERGED Edges
  CDS_CID_t next;  /* Index into edges */
  CDS_CID_t curr;  /* Index into edges */
  CDS_CID_t prev;  /* Index into edges */
  // Indices of Raw Edges
  CDS_CID_t nextRaw;  /* Index into edges */
  CDS_CID_t currRaw;  /* Index into edges */
  CDS_CID_t prevRaw;  /* Index into edges */
  int includeContainment; /* Include the containment edges?  Default is no */
  int confirmedOnly;
  int rawOnly;
  int end; /* A_END, B_END or ALL_END */
  unsigned int edgeStatusSet; /* See EdgeStatus and TRUSTED_EDGES and ALL_EDGES */
  int verbose;
  CDS_CID_t cid;
  GraphCGW_T *graph;
}GraphEdgeIterator;




#define GRAPH_EDGE_DEFAULT 0
#define GRAPH_EDGE_CONFIRMED_ONLY 1
#define GRAPH_EDGE_RAW_ONLY 2
#define GRAPH_EDGE_VERBOSE 4
//#define GRAPH_EDGE_OVERLAP_ONLY 8
#define GRAPH_EDGE_EXCLUDE_CONTAINMENT 16

static void InitGraphEdgeIterator(GraphCGW_T *graph,
				  CDS_CID_t cid,
				  int end, 
				  int edgeStatusSet,
				  int flags,
				  GraphEdgeIterator *e){
  NodeCGW_T *CI;

  assert(graph && e);

  CI = GetGraphNode(graph, cid);
  e->curr = e->prev = NULLINDEX ;
  e->next = CI->edgeHead;
  e->currRaw = e->nextRaw = e->prevRaw = NULLINDEX;
  e->graph = graph;
  assert(end & ALL_END);
  e->confirmedOnly = flags & GRAPH_EDGE_CONFIRMED_ONLY;
  e->includeContainment = !(flags & GRAPH_EDGE_EXCLUDE_CONTAINMENT);
  e->rawOnly = flags & GRAPH_EDGE_RAW_ONLY;;
  e->verbose = flags & GRAPH_EDGE_VERBOSE;
  e->end = end;
  e->edgeStatusSet = edgeStatusSet;
  e->cid = cid;
  if(e->verbose)
  fprintf(stderr,
          "* Iterator for CI " F_CID " end %d  head = " F_CID " confirmed:%d raw:%d \n",
	  cid, e->end, e->next,e->confirmedOnly, e->rawOnly);
}


static  EdgeCGW_T *NextGraphEdgeIterator(GraphEdgeIterator *e){
  int isA;
  EdgeCGW_T *r = (EdgeCGW_T *)NULL;
  EdgeCGW_T *retEdge = (EdgeCGW_T *)NULL;
  assert(e->graph && (e->end & ALL_END));

  if(e->verbose)
    fprintf(stderr,
            "* NextGraphEdgeIterator nextRaw:" F_CID " prev:" F_CID " curr:" F_CID " next:" F_CID "\n",
	  e->nextRaw, e->prev, e->curr, e->next);

  if(e->nextRaw == NULLINDEX){
    if(e->currRaw != NULLINDEX){ // do this once
      e->prevRaw = e->currRaw;
      e->currRaw = e->nextRaw;
    }
    if(e->next == NULLINDEX){
      if(e->curr != NULLINDEX){ // do this once
	e->prev = e->curr;
	e->curr = e->next;
      }
      if(e->verbose)
	fprintf(stderr,"* Fell off the end \n");
      return retEdge;
    }
  }

  while(retEdge == (EdgeCGW_T *)NULL &&
	(e->nextRaw != NULLINDEX || e->next != NULLINDEX)){
    ChunkOrientationType orient;
  
    if(e->verbose)
      fprintf(stderr,"* In loop (" F_CID "," F_CID "," F_CID ")\n",
	      e->prev, e->curr, e->next);

    /* We are iterating through the raw edges of a merged edge */
    if(e->nextRaw != NULLINDEX){
      // These edges are guaranteed to be the right orientation and end points
      if(e->verbose)
	fprintf(stderr,"* Getting a raw edge\n");

      r = GetGraphEdge(e->graph, e->nextRaw);
      AssertPtr(r);

      assert(!r->flags.bits.isDeleted &&
             (r->idA == e->cid ||
	      r->idB == e->cid));

      e->prevRaw = e->currRaw;
      e->currRaw = e->nextRaw;
      e->nextRaw = r->nextRawEdge;
      retEdge = r;
    }else{
      /* We are iterating through the merged and unmerged raw edges */
      r = GetGraphEdge(e->graph, e->next);
      AssertPtr(r);

      assert(r->idA == e->cid ||
	     r->idB == e->cid);

      assert(!r->flags.bits.isDeleted);

      // Flip orientation to accomodate canonical graph form
      orient = GetEdgeOrientationWRT(r, e->cid);

      /* Check for correct end and confirmed status */
      switch(orient){
        /* EdgeMate from the A-End */
	case BA_BA:
	case BA_AB:
	  if((e->end & A_END) &&
	     (e->confirmedOnly == 0 || isConfirmedEdge(r)) &&
	     (e->includeContainment || !isContainmentEdge(r))){
	    retEdge = r;
	  }else{
	    if(e->verbose)
	      fprintf(stderr,
                      "* Skipping edge (" F_CID "," F_CID ") with orient:%c orientWRT %c (e->end&A_END) = %d confirmedOnly = %d failedContainment:%d\n",
		      r->idA, r->idB, r->orient, orient,e->end&A_END, e->confirmedOnly,(isContainmentEdge(r) && e->includeContainment == FALSE));
	  }
	  break;
	case AB_BA:
	case AB_AB:
	  if((e->end & B_END) &&
	     (e->confirmedOnly == 0 || isConfirmedEdge(r)) &&
	     (e->includeContainment || !isContainmentEdge(r))){
	    retEdge = r;
	  }else{
	    if(e->verbose)
	      fprintf(stderr,"* Skipping edge (" F_CID "," F_CID ") with orient:%c orientWRT %c (e->end&B_END) = %d confirmedOnly = %d failedContainment:%d\n",
		      r->idA, r->idB, r->orient, orient,e->end&B_END, e->confirmedOnly,(isContainmentEdge(r) && e->includeContainment == FALSE));
	  }
	  break;
	default:
	  assert(0);
      }
      
      e->prev = e->curr;
      e->curr = e->next;
      isA = (r->idA == e->cid);
      if(isA){
	e->next = r->nextALink;
      }else{
	e->next = r->nextBLink;
      }

      // Found a top level (not a raw under a merged) that is inappropriate
      // Only top level edges are marked with status values
      if(retEdge &&
         !((uint32)GetEdgeStatus(retEdge) & (uint32)e->edgeStatusSet)){
	  if(e->verbose)
	    fprintf(stderr,"* Looking for status 0x%x only, found an edge with status 0x%x\n",
		    e->edgeStatusSet, GetEdgeStatus(retEdge));
	  retEdge = NULL;    // This stops us from exiting the loop */
      }
      /* If we are iterating over raw edges, dive into teh raw edge list */
      /* assert((!retEdge->flags.bits.isRaw &&
	      (retEdge->nextRawEdge != NULLINDEX)) ||
	     (retEdge->flags.bits.isRaw &&
	     (retEdge->nextRawEdge == NULLINDEX))); */
      if(retEdge && !retEdge->flags.bits.isRaw && e->rawOnly){
	  if(e->verbose)
	    fprintf(stderr,"* Looking for raw only, found a merged edge\n");
	  e->nextRaw = r->nextRawEdge;
	  retEdge = NULL;    // This stops us from exiting the loop */
      }
    }
  }
  if(retEdge && e->verbose){
    fprintf(stderr,"* Found CIEdge (" F_CID "," F_CID ") with orient %c and e->end = %d status:%d (" F_CID "," F_CID "," F_CID ")\n",
	    retEdge->idA, retEdge->idB,
	    retEdge->orient, e->end,
	    GetEdgeStatus(retEdge),
	    e->prev, e->curr, e->next);
  }
  assert(e->curr == NULLINDEX || e->curr != e->next);

  return retEdge;
}

/****************************************************************************/

/* Diagnostic */
size_t ReportMemorySizeGraphCGW(GraphCGW_T *graph, FILE *stream);




/* Operations on Edges */

/* MergeNodeGraphEdges
 *   Merge the edges incident on a particular node
 *
 */
void MergeNodeGraphEdges(GraphCGW_T *graph, NodeCGW_T *node,
                         int includeGuides, int mergeAll, int debug);

// Merge all of the edges
void MergeAllGraphEdges(GraphCGW_T *graph, int includeGuides);



/* MergeGraphEdges
   Input:  GraphCGW_T *graph    the edges manipulated are all in graph->edges
           VA_TYPE(int) *inputEdges  of references to edges that link the
           same two IDs, to be merged. The merged CIEdgeTs reference the
           raw edges they incorporate by a singly linked list via the
           nextRawEdge field.  The merged edges are marked as not raw.
           The merged edges are APPENDED to the edges array.

   The return value is the number of edges that this routine appended to
   the edges array.

   The client of this routine must INSERT the resulting edges into the
   appropriate graph.
*/
int MergeGraphEdges(GraphCGW_T *graph, VA_TYPE(CDS_CID_t) *inputEdges);


/* FindGraphEdgeChain
   Starting from the CIEdge with index eid, identify a chain of edges
   of length > 1 that connect the same pair of CIs.  The indices of these
   edges are collected in the chain VA, and the entire chain is UNLINKED
   from the graph.
   
   The return value is the length of the identified chain.  If 1 is returned,
   the graph has not been perturbed, unless extractSingletons is TRUE.
*/
CDS_CID_t FindGraphEdgeChain(GraphCGW_T *graph, CDS_CID_t eid,
                             VA_TYPE(CDS_CID_t) *chain,
                             int extractSingletons,
                             int includeGuides);

int GraphEdgeSanity(GraphCGW_T *graph, CDS_CID_t eid);

static void SetGraphEdgeStatus(GraphCGW_T *graph, EdgeCGW_T *edge,
                               EdgeStatus status){
  edge->flags.bits.edgeStatus = (unsigned int) status;
  if(edge->flags.bits.isRaw)
    return;
  // Propagate to raw edges that are attached
  while(NULL != (edge = GetGraphEdge(graph, edge->nextRawEdge))){
    edge->flags.bits.edgeStatus = (unsigned int) status;
  }
}






// Initialize the status flags for the given edge.
void InitGraphEdgeFlags(GraphCGW_T *graph, EdgeCGW_T *edge);

int32 InsertGraphEdge(GraphCGW_T *graph,  CDS_CID_t cedgeID, int verbose);
int32 InsertComputedOverlapEdge(GraphCGW_T *graph,
				ChunkOverlapCheckT *olap);
void InsertGraphEdgeInList(GraphCGW_T *graph, CDS_CID_t CIedgeID,
                           CDS_CID_t sid, int verbose);
void PrintGraphEdge(FILE *logfp, GraphCGW_T *graph, char *label,
                    EdgeCGW_T *edge, CDS_CID_t cid);

void PrintContigEdgeInScfContext(FILE *logfp, GraphCGW_T *graph, char *label, EdgeCGW_T *edge, int cid);



CDS_CID_t AddGraphEdge( GraphCGW_T *graph,
                        CDS_CID_t cidA, CDS_CID_t cidB, 
                        CDS_CID_t fidA, CDS_CID_t fidB,
                        CDS_CID_t dist,
                        LengthT distance,
                        float32 quality,
                        CDS_COORD_t fudgeDistance,
                        OrientType orientation,
                        int isInducedByUnknownOrientation,
                        int isGuide,    // Add it to GuideMates and flag it
                        int isSTSGuide,
                        int isMayJoin,
                        int isMustJoin,
                        int isOverlap,
                        int isRepeatOverlap,
                        int isTandemOverlap,
                        int isTransChunk,
                        int isAContainsB,
                        int isBContainsA,
                        int isExtremalA,
                        int isExtremalB,
                        EdgeStatus status,
                        int collectOverlap,
                        int insert); // insert the edge in the graph, if true

// Delete a node and all incident edges
void DeleteGraphNode(GraphCGW_T *graph, NodeCGW_T *node);

// Unlink a graph edge from the graph, without marking it as deleted
void UnlinkGraphEdge(GraphCGW_T *graph, EdgeCGW_T *edge);

// Delete a top level CIEdgeT by unlinking it from the graph - remember
// this will unlink any dangling raw edges as well - and mark it as deleted.
void  DeleteGraphEdge(GraphCGW_T *graph,  EdgeCGW_T *edge);

// Find an overlap edge..returns NULL if none found
EdgeCGW_T *FindGraphOverlapEdge(GraphCGW_T *graph,
                                CDS_CID_t idA, CDS_CID_t idB,
                                ChunkOrientationType orient);

// Find an edge..returns NULL if none found
EdgeCGW_T *FindGraphEdge(GraphCGW_T *graph,
                         CDS_CID_t idA, CDS_CID_t idB,
                         ChunkOrientationType orient);

// Find an overlap edge (assumed to exist) between the two CIs.
// Unlink it from the graph
void  DeleteGraphOverlapEdge(GraphCGW_T *graph,
                             CDS_CID_t idB, CDS_CID_t idA,
                             ChunkOrientationType orient);

// Move the edge to the free list
void FreeGraphEdgeByEID(GraphCGW_T *graph,  CDS_CID_t eid);

static void FreeGraphEdge(GraphCGW_T *graph,  EdgeCGW_T *edge){
  CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, edge);
  
  FreeGraphEdgeByEID(graph, eid);
}

// Get the edge from the free list
EdgeCGW_T *GetFreeGraphEdge(GraphCGW_T *graph);


// Clean hash table
// Move deleted nodes and edges to their respective free lists
void RecycleDeletedGraphElements(GraphCGW_T *graph);


void DumpGraph(GraphCGW_T *graph, FILE *stream);

// Compute all of the potential overlaps that have been collected
void ComputeOverlaps(GraphCGW_T *graph, int addEdgeMates,
                     int recomputeCGBOverlaps);

// Dump Overlaps
void DumpOverlaps(GraphCGW_T *graph);

int IsRepeatOverlap(GraphCGW_T *graph,
                    CDS_CID_t cid1, CDS_CID_t cid2,
                    ChunkOrientationType orient, LengthT overlap);

/* Check that edge with index eid is properly wired in the graph:
   - we find it when looking for it in both lists which it is supposed
     to be a member
   - we can walk backwards from the edge to the heads of the two lists
*/

// If the chunks DON'T Overlap, returns NULL

ChunkOverlapCheckT OverlapChunks(GraphCGW_T *graph,
                                 CDS_CID_t cidA, CDS_CID_t cidB,
                                 ChunkOrientationType orientation, 
                                 CDS_COORD_t minOverlap,
                                 CDS_COORD_t maxOverlap,
                                 float errRate,
                                 int insertGraphEdges);

Overlap*  OverlapSequences(char * seq1, char * seq2,
                           ChunkOrientationType orientation, 
                           CDS_COORD_t min_ahang, CDS_COORD_t max_ahang,
                           double erate, double thresh, CDS_COORD_t minlen,
                           CompareOptions what);

Overlap* OverlapContigs(NodeCGW_T *contig1, NodeCGW_T *contig2, 
                        ChunkOrientationType *overlapOrientation,
                        CDS_COORD_t minAhang, CDS_COORD_t maxAhang,
                        int computeAhang);
/*
Overlap* OverlapContigsLocal(NodeCGW_T *contig1, NodeCGW_T *contig2, 
                             ChunkOrientationType overlapOrientation,
                             int minAhang, int maxAhang, int computeAhang);
*/

/* OverlapChunksWithMultipleCoverage:
   This is a special version of OverlapChunks that lops off the
   single coverage ends of the chunks prior to overlapping */
BranchPointResult OverlapChunksWithBPDetection(GraphCGW_T *graph,
                                               CDS_CID_t cidA, CDS_CID_t cidB,
                                               ChunkOrientationType orientation, 
                                               CDS_COORD_t minOverlap,
                                               CDS_COORD_t maxOverlap);

int32  SmallOverlapExists(GraphCGW_T *graph,
                          CDS_CID_t cidA, CDS_CID_t cidB,
                          ChunkOrientationType orientation,
                          CDS_COORD_t minOverlap);

int32  LargeOverlapExists(GraphCGW_T *graph,
                          CDS_CID_t cidA, CDS_CID_t cidB,
                          ChunkOrientationType orientation,
                          CDS_COORD_t minOverlap, CDS_COORD_t maxOverlap);

void CollectChunkOverlap(GraphCGW_T *graph,
                         CDS_CID_t cidA, CDS_CID_t cidB,
                         ChunkOrientationType orientation, 
                         float32 meanOverlap, float32 deltaOverlap,
                         float32 quality, int bayesian,
                         int fromCGB, int verbose);


int32 GetGappedMultipleCoverageInterval(GraphCGW_T *graph, CDS_CID_t cid,
                                        SeqInterval *interval, int end);


void CheckEdgesAgainstOverlapper(GraphCGW_T *graph);
int AddImplicitOverlaps_(GraphCGW_T *graph);
// If the chunks DON'T Overlap, returns NULL


/* Given a graph edge, create an overlap in the hashtable */
void CreateChunkOverlapFromEdge(GraphCGW_T *graph,
                                EdgeCGW_T *edge,
                                int bayesian);


void FillChunkOverlapWithEdge(EdgeCGW_T *edge, ChunkOverlapCheckT *olap);

int LookupOverlap(GraphCGW_T *graph,
		  CDS_CID_t cidA, CDS_CID_t cidB,
		  ChunkOrientationType orientation,
		  ChunkOverlapCheckT *olap);

/*
  The function LookupQualityOverlap takes the same arguments as LookupOverlap
  with an additional argument for the quality function and the quality value. 
  It first looks up the overlap
  in the hash table. If it has the appropriate quality bit set the quality
  is returned. If not, the quality is computed using the appropriate quality
  function and stored in the hash table. In addition, the appropriate bit
  is set indicating that computation.
*/
int LookupQualityOverlap(GraphCGW_T *graph,
			 EdgeCGW_T *edge,
			 ChunkOrientationType orientation,
			 ChunkOverlapCheckT *olap, QualityFuncT qfunc,
			 float* quality, FILE* log);


int ComputeQualityOverlap(GraphCGW_T *graph, 
			  EdgeCGW_T *edge,
			  ChunkOrientationType orientation,
			  ChunkOverlapCheckT *olap, QualityFuncT qfunc,
			  float* quality, FILE* log);


int ComputeUOMQualityOverlap(GraphCGW_T *graph, 
			     UnitigOverlapMesg *uom_mesg,
			     ChunkOverlapCheckT *olap,
			     float* quality);
/* Update all the fragments belonging to the multiAlignment for Chunk cid
   so that their membership and offsets in their CIFragT record are recorded
   properly. Call this on either a CI or a Contig after it is created */
void UpdateNodeFragments(GraphCGW_T *graph, CDS_CID_t cid,
                         int markFragmentsPlaced, int markUnitigAndContig);

void UpdateNodeUnitigs(MultiAlignT *ma, NodeCGW_T *contig);


/* Compute the offset and orientation of a fragment in its chunk
   orientIsOpposite == TRUE
     Offset is from 5p end of fragment to the end of the chunk in the
     direction of the 3p end of the fragment.
   orientIsOpposite == FALSE
     Offset is from 5p end of fragment to the end of the chunk in the
     direction of the 5p end of the fragment.
*/
/* orientIsOpposite:
     TRUE if offset should be calculated from 5' towards end of
       chunk closest to 3' end of fragment.  
     FALSE if offset should be calculated from 5' towards end of
       chunk closest to 5' end.
*/
int FragOffsetAndOrientation(CIFragT     *frag,
			     NodeCGW_T *chunk,
			     LengthT    *chunkOffset, // output
			     FragOrient *chunkOrient, // output
			     int32 *extremal,         // output
			     int32 orientIsOpposite);

// GraphEdgeStat is used to collect stats on graph edge building
typedef struct {
  int32 totalFragments;
  int32 totalMatePairs;
  int32 totalExternalMatePairs;
  int32 totalBacPairs;
  int32 totalUUBacPairs;
}GraphEdgeStatT;

static void InitGraphEdgeStatT(GraphEdgeStatT *stat){
  stat->totalFragments = 0;
  stat->totalMatePairs = 0;
  stat->totalExternalMatePairs = 0;
  stat->totalBacPairs = 0;
  stat->totalUUBacPairs = 0;
}

/*
  Create All raw link-based graph edges directly from fragment links
  and multi-alignments
*/
void  BuildGraphEdgesDirectly(GraphCGW_T *graph);

// Create the raw link-based edges incident on a particular graph node
void  BuildGraphEdgesFromMultiAlign(GraphCGW_T *graph, NodeCGW_T *node,
                                    MultiAlignT *ma, GraphEdgeStatT *stats,
                                    int buildAll);

void    UpdateScaffoldSimCoordinates(NodeCGW_T *scaffold);

void    UpdateContigSimCoordinates(NodeCGW_T *contig);


// We want tandem overlaps to have large variances, roughly the
// twoce variance of a 10k mate link.
// Sigma of a 10k link is 1/3 of 20% of 10000
#define TANDEM_OVERLAP_STDDEV (2.0 * 0.20 * 10000.0/3.0)
#define TANDEM_OVERLAP_VARIANCE (TANDEM_OVERLAP_STDDEV * TANDEM_OVERLAP_STDDEV)

void PropagateTandemMarks(GraphCGW_T *graph);

int  MarkTandemEdge(GraphCGW_T *graph, EdgeCGW_T *edge);



/**********************

   Split an unresolved CI, moving a subset of its fragments to the new node.
   Returns the index of the new node, or NULLINDEX if failure.
   The new CI gets copies of all of the overlap edges of its parent, and
   inherits some of the mate edges, as a function of fragments selected.
   The fragments listed are marked for membership (via their CIid field) in
   the new element
   An empty fragments array is handled the same as fragments == NULL.
 */
int32 SplitUnresolvedCI(GraphCGW_T *graph, CDS_CID_t nodeID,
                        VA_TYPE(CDS_CID_t) *fragments);

/**********************
   Split an unresolved Contig.
   This involves splitting its underling CI, and creating a new Contig
   containing the split CI.
   The new CI gets copies of all of the overlap edges of its parent, and
   inherits some of the mate edges, as a function of fragments selected.
   Returns the index of the new node, or NULLINDEX if failure.
   The fragments listed are marked for membership (via their contigID field)
   int he new element
   An empty fragments array is handled the same as fragments == NULL.
   If copyAllOverlaps == TRUE, all overlap edges that are NOT indicative of
   the split Contig being contained, are duplicated and assigned to the new
   contig.  If FALSE, only overlaps indicative of the contig CONTAINING
   another contig are duplicated and assigned.
 */
int32 SplitUnresolvedContig(GraphCGW_T *graph, CDS_CID_t nodeID,
                            VA_TYPE(CDS_CID_t) *fragments,
                            int32 copyAllOverlaps);

/* existsContainmentRelationship */
UnitigOverlapType existsContainmentRelationship(NodeCGW_T *ci,
                                                NodeCGW_T *otherCI);

/*
  ComputeMatePairStatistics:
   Compute the mate pair distance distributions on either contigs
   (operateOnContigs == TRUE) or Unitigs, optionally updating the
   nominal values stored in the DistT records.

   Should also bucketize the data.
*/

#define UNITIG_OPERATIONS   0
#define CONTIG_OPERATIONS   1
#define SCAFFOLD_OPERATIONS 2

void ComputeMatePairStatistics( int operateOnNodes,
                                int minSamplesForOverride,
                                char *instance_label);
void ComputeMatePairStatisticsRestricted( int operateOnNodes,
                                          int minSamplesForOverride,
                                          char *instance_label);



/* Find the original contig ID of a surrogate contig */
CDS_CID_t GetOriginalContigID(CDS_CID_t contigID);

/* IsDefinitielyUniqueContig
 Returns TRUE if this is either a multi-unitig contig, or a single unitig
 contig where the unitig has a coverageState above threshhold
*/
int32 IsDefinitelyUniqueContig(NodeCGW_T *contig);

int32 IsShakyContigAtScaffoldEnd(NodeCGW_T *contig);

/* Restore the Edge Means we squirreled away during walking */
void RestoreEdgeMeans(GraphCGW_T *graph);

void CheckGraph(GraphCGW_T *graph);

void ReallocGraphEdges(GraphCGW_T *graph, int32 numEdges);

void InitGraphEdge(EdgeCGW_T *edge);

void AssignFragsToResolvedCI(GraphCGW_T *graph,
                             CDS_CID_t fromID, CDS_CID_t toID,
                             VA_TYPE(CDS_CID_t) *fragments);

#endif


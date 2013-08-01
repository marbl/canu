
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

#ifndef GRAPH_CGW_H
#define GRAPH_CGW_H

static const char *rcsid_GRAPH_CGW_H = "$Id$";

#include "AS_UTL_Var.H"
#include "AS_CGW_dataTypes.H"
#include "InputDataTypes_CGW.H"
#include "MultiAlign.H"

#include <vector>
#include <set>

using namespace std;

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
  PairOrient orient; // Orientation of idA <-> idB, BE CAREFUL IF YOU WANT idB<->idA
  //
  int32 edgesContributing;
  double quality;   // Used to order edges by decreasing quality (quality = 0.0 is the best, 1.0 is the worst)

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
      // 8 bits used

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

      // 16 bits used

      unsigned int aContainsB:1;        /* From CGB for multiply contained fragments */
      unsigned int bContainsA:1;        /* From CGB for multiply contained fragments */
      unsigned int mustOverlap:1;       /* Marked on merged edges when mate-link variance is signficantly less than overlap length */

      // 24 bits used

      unsigned int hasTransChunk:1;           /* was a transitively removed edge in cgb */
      unsigned int hasContainmentOverlap:1;  /* Implies a containment */
      unsigned int isRaw:1;                   /* True for raw edges, false for merged edges */

      unsigned int hasExtremalAFrag:1; /* Used in merging raw edge mates, meaningless elsewhere */
      unsigned int hasExtremalBFrag:1; /* Used in merging raw edge mates, meaningless elsewhere */
      unsigned int rangeTruncated:1;    /* TRUE if we looked for an overlap and didn't find any,
                                           and thus truncated the range of distances */
      unsigned int inAssembly:1;

      // 32 bits used

      unsigned int isBogus:1;  // determined from simulator annotations, overloaded by scaffold merging code (see CIScaffold_Merge_CGW.c)
      unsigned int isProbablyBogus:1; // determined from distance and distIndex empirically, overloaded for one-sided edges in CIScaffold_Merge
      unsigned int hasConfirmingPath:1; // Has this edge another path that confirms its length & var (used in GapWalkerREZ.c)


      unsigned int edgeStatus:7;

      // 42 bits used

      unsigned int isMarkedForDeletion:1;   // We plan to delete this guy
      unsigned int MeanChangedByWalking:1;
      unsigned int highQualityA:1;           // One of the top-ranked edges incident on node idA
      unsigned int highQualityB:1;           // One of the top-ranked edges incident on node idB
      unsigned int isSloppy:1;
      unsigned int isBridge:1;               // Bridge edge in a scaffold

      // 48 bits used total
    }bits;
  }flags;
  //
  LengthT distance; // gap/overlap length

  /* vvvvv Iterator does not fill in fields below this line vvvvvv */
  //
  CDS_CID_t nextALinkUNUSED; // next edge involving cidA, -1 if none
  CDS_CID_t nextBLinkUNUSED; // next edge involving cidB, -1 if none
  CDS_CID_t prevALinkUNUSED; // prev edge involving cidA, -1 if none
  CDS_CID_t prevBLinkUNUSED; // prev edge involving cidB, -1 if none

  double  minDistance;  /* negative implies potential overlap
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

typedef EdgeCGW_T CIEdgeT;
typedef EdgeCGW_T ContigEdgeT;
typedef EdgeCGW_T SEdgeT;



class EdgeCGWLabel_T {
public:
  CDS_CID_t   idA;
  CDS_CID_t   idB;
  PairOrient  orient;
  LengthT     distance;

  bool  operator<(EdgeCGWLabel_T const &that) const {
    if (idA < that.idA)   return(true);
    if (idA > that.idA)   return(false);

    if (idB < that.idB)   return(true);
    if (idB > that.idB)   return(false);

    if (orient.toLetter() < that.orient.toLetter())   return(true);
    if (orient.toLetter() > that.orient.toLetter())   return(false);

    if (distance.mean < that.distance.mean)   return(true);
    if (distance.mean > that.distance.mean)   return(false);

    return(false);
  }
};


#ifdef TRACK_MATE_PAIR_TEST
//  For linking matepair tests with merge results
//
struct mergeTestResult {
public:
  mergeTestResult();
  ~mergeTestResult();

  void     report(void) {
    fprintf(stderr, "%s edgeID %u %u +- %u weight %u -- scaffold %u (%ubp %u contigs h/g/b %.2f/%.2f/%.2f) -- scaffold %u (%ubp %u contigs h/g/b %.2f/%.2f/%.2f) -- MERGE (h/g/b %.2f/%.2f/%.2f) -- BEFORE %.3f %.2f/%.2f -- AFTER %.3f %.2f/%.2f\n",
            mergeAccepted ? "PASS" : "FAIL",
            edgeID, edgeLength, edgeVariance, edgeWeight,
            scaffoldAid, scaffoldAlength, scaffoldAcontigs, Ahappy, Agap, Abad,
            scaffoldBid, scaffoldBlength, scaffoldBcontigs, Bhappy, Bgap, Bbad,
            Mhappy, Mgap, Mbad,
            beforeSatisfied, beforeGood, beforeBad,
            afterSatisfied, afterGood, afterBad);
  };

  uint32   edgeID;
  uint32   edgeLength;
  uint32   edgeVariance;
  uint32   edgeWeight;

  uint32   scaffoldAid;
  uint32   scaffoldAlength;
  uint32   scaffoldAcontigs;

  uint32   scaffoldBid;
  uint32   scaffoldBlength;
  uint32   scaffoldBcontigs;

  double   Ahappy;
  double   Agap;
  double   Abad;

  double   Bhappy;
  double   Bgap;
  double   Bbad;

  double   Mhappy;
  double   Mgap;
  double   Mbad;

  double   beforeSatisfied;
  double   beforeGood;
  double   beforeBad;

  double   afterSatisfied;
  double   afterGood;
  double   afterBad;

  bool     mergeAccepted;
};
#endif




typedef struct
{
  int32 samples;
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

typedef UnitigFUR ChunkFUR;


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


VA_DEF(CDS_CID_t)


typedef struct{
  ChunkInstanceType type; //

  CDS_CID_t id;        // Instance ID
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
  LengthT  offsetDelta;

  union{  // ChunkInstanceType discriminates
    struct CIINFO_TAG {
      CDS_CID_t    contigID;   // contigID -- if -1, this chunkInstance not merged into a contig
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
      CDS_CID_t source;
    }CI;

    struct CONTIGINFO_TAG{
      CDS_CID_t AEndCI;     //  Index of Chunk Instance at A end of Contig
      CDS_CID_t BEndCI;     //  Index of Chunk Instance at B end of Contig
      int32     numCI;      //  Number of CI in contig
      CDS_CID_t contigNum;  //  Used in CIScaffoldT_Biconnected_CGW.c
    }Contig;

    struct CISCAFFOLD_TAG{
      CDS_CID_t AEndCI; // Index of Chunk Instance at A end of Scaffold
      CDS_CID_t BEndCI; // Index of Chunk Instance at B End of Scaffold
      int32 numElements; // If containsCIs, these are CIs, else they are Contigs
      /*** Info computed by RecomputeScaffoldPositions ***/
      double  leastSquareError;   // Measure of chi-squared of computed positions
      int32 numLeastSquareClones; // Relates to degrees of freedom for chi-square calculation
      /*** Info computed by MarkInternaledges ***/
      int32 internalEdges;  // Number of merged edges (not including UNTRUSTED) that are internal to scaffold
      int32 confirmedInternalEdges; // Number of merged edges confirmed by current CI positions
    }Scaffold;

  }info;

  //
  union{
    struct {
      unsigned int isUnique:1;
      unsigned int smoothSeenAlready:1; // Used for cycle detection

      unsigned int isDead:1;
      unsigned int isFree:1;
      unsigned int containsCIs:1;        // Scaffold contains either CIs or Contigs

      unsigned int isCI:1;
      unsigned int isContig:1;
      unsigned int isScaffold:1;

      /* The following is TRUE for CIs that are surrogate CIs and
         for Contigs that contain a single surrogate CI.  Set at creation
         time in the Split functions */
      unsigned int isSurrogate:1;
      unsigned int beingContigged:1;

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
      unsigned int isStone:1;
      unsigned int isWalk:1;
      unsigned int isRock:1;
      unsigned int isPotentialRock:1; // has at least 2 celera reads with external mates
      unsigned int isPotentialStone:1; // has at least 1 celera reads with external mates
      
      /* THe following is used to closure read placement to identify candidate closure reads
       * Initialized to match the closure input status of the read or contig
       */
      unsigned int isClosure:1;
      unsigned int isJiggled:1;

      //  When rebuilding scaffold edges for only newly created scaffolds, we need to know if the
      //  edge is between two new scaffolds, or between a new and an old.  This flags marks the new
      //  scaffolds.
      unsigned int isHavingEdgesBuilt:1;

      unsigned int unused:5;
    }bits;
    int32 all;
  }flags;

  CDS_CID_t edgeHeadUNUSERD;  // Pointer to linked list of edges  in edges;
  CDS_CID_t setID;

}NodeCGW_T;

typedef NodeCGW_T ChunkInstanceT;
typedef NodeCGW_T ContigT;
typedef NodeCGW_T CIScaffoldT;

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

  VA_TYPE(NodeCGW_T)      *nodes;
  VA_TYPE(EdgeCGW_T)      *edges;

  vector<set<CDS_CID_t,bool(*)(CDS_CID_t,CDS_CID_t)> >  edgeLists;

  int32 numActiveNodes;
  int32 numActiveEdges;

  CDS_CID_t freeEdgeHead;
  CDS_CID_t tobeFreeEdgeHead; // staging area for edges waiting to be moved to the free list
  CDS_CID_t freeNodeHead;
  CDS_CID_t tobeFreeNodeHead; // staging area for nodes waiting to be moved to the free list
  CDS_CID_t deadNodeHead;     //  UNUSED
}GraphCGW_T;


typedef struct{
  CDS_CID_t id;
  int32 minGap;
  int32 maxGap;
  PairOrient orient;
}RevivedEdgeT;


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


void ResizeEdgeList(GraphCGW_T *graph);

static void SetGraphNode(GraphCGW_T *graph, NodeCGW_T *node) {
  SetNodeCGW_T(graph->nodes, node->id, node);
  ResizeEdgeList(graph);
}

/* Append */
static void AppendGraphNode(GraphCGW_T *graph, NodeCGW_T *node){
  AppendNodeCGW_T(graph->nodes, node);
  ResizeEdgeList(graph);
}
static void AppendGraphEdge(GraphCGW_T *graph, EdgeCGW_T *edge){
  AppendEdgeCGW_T(graph->edges, edge);
}


/* Create New */
static NodeCGW_T *CreateNewGraphNode(GraphCGW_T *graph){
  NodeCGW_T node;
  memset(&node, 0, sizeof(NodeCGW_T));

  node.id         = GetNumGraphNodes(graph);
  node.flags.all  = 0;
  node.setID      = 0;

  switch(graph->type){
    case CI_GRAPH:
      //fprintf(stderr, "CreateNewGraphNode()-- CI %d\n", node.id);
      node.type                      = UNRESOLVEDCHUNK_CGW;
      node.flags.bits.isCI           = TRUE;
      node.info.CI.contigID          = NULLINDEX;
      node.info.CI.numInstances      = 0;
      node.info.CI.instances.va      = NULL;
      break;
    case CONTIG_GRAPH:
      //fprintf(stderr, "CreateNewGraphNode()-- Contig %d\n", node.id);
      node.type                      = CONTIG_CGW;
      node.flags.bits.isContig       = TRUE;
      node.info.Contig.AEndCI        = NULLINDEX;
      node.info.Contig.BEndCI        = NULLINDEX;
      node.info.Contig.numCI         = 0;
      break;
    case SCAFFOLD_GRAPH:
      //fprintf(stderr, "CreateNewGraphNode()-- Scaffold %d\n", node.id);
      node.type                      = REAL_SCAFFOLD;
      node.flags.bits.isScaffold     = TRUE;
      node.info.Scaffold.AEndCI      = NULLINDEX;
      node.info.Scaffold.BEndCI      = NULLINDEX;
      node.info.Scaffold.numElements = 0;
      break;
    default:
      assert(0);
  }

  node.scaffoldID          = NULLINDEX;
  node.smoothExpectedCID   = NULLINDEX;
  node.numEssentialA       = 0;
  node.numEssentialB       = 0;
  node.essentialEdgeA      = NULLINDEX;
  node.essentialEdgeB      = NULLINDEX;
  node.indexInScaffold     = NULLINDEX;
  node.smoothExpectedCID   = 0;
  node.AEndNext            = NULLINDEX;
  node.BEndNext            = NULLINDEX;

  node.bpLength.mean       = 0.0;
  node.bpLength.variance   = 0.0;

  node.offsetAEnd.mean     = 0.0;
  node.offsetAEnd.variance = 0.0;

  node.offsetBEnd.mean     = 0.0;
  node.offsetBEnd.variance = 0.0;

  //node.aEndCoord = node.bEndCoord = NULLINDEX;

  AppendGraphNode(graph, &node);

  return(GetGraphNode(graph, node.id));
}

/* Constructor and Destructor */
GraphCGW_T *CreateGraphCGW(GraphType type, int32 numNodes, int32 numEdges);
void DeleteGraphCGW(GraphCGW_T *graph);

/* Persistence */
void        SaveGraphCGWToStream(GraphCGW_T *graph, FILE *stream);
GraphCGW_T *LoadGraphCGWFromStream(FILE *stream);
void        RebuildGraphEdges(GraphCGW_T *graph);


static EdgeStatus GetEdgeStatus(EdgeCGW_T *edge){
  return (EdgeStatus) edge->flags.bits.edgeStatus;
}


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


static int32 IsSurrogateNode(NodeCGW_T *node){
  return node->flags.bits.isSurrogate;
}


static SequenceOrient GetNodeOrient(NodeCGW_T *CI){
  SequenceOrient  orient;

  if(CI->offsetBEnd.mean > CI->offsetAEnd.mean)
    orient.setIsForward();
  else
    orient.setIsReverse();

  return(orient);
}


static PairOrient GetEdgeOrientationWRT(EdgeCGW_T* edge, CDS_CID_t wrtCI){
  PairOrient ret = edge->orient;
  if(edge->idA != wrtCI)
    ret.flip();
  return(ret);
}

// which side of the seed is a chunk on in the scaffold
typedef enum {LEFT_SIDE, RIGHT_SIDE, THE_SEED} SeedSideType;

static SeedSideType GetChunkSeedSide( PairOrient edgeOrient){
  assert(edgeOrient.isUnknown() == false);
  if (edgeOrient.isAB_AB() || edgeOrient.isAB_BA())
    return(RIGHT_SIDE);
  return(LEFT_SIDE);
}

/* Always consider the seed in the AB orientation, so we flip
   the other chunk accordingly to make it so */
static SequenceOrient GetRelativeChunkOrientation( PairOrient edgeOrient){
  SequenceOrient  ret;

  assert(edgeOrient.isUnknown() == false);

  if (edgeOrient.isAB_AB() || edgeOrient.isBA_BA())
    ret.setIsForward();
  else
    ret.setIsReverse();
  return(ret);
}

static PairOrient GetChunkPairOrientation(SequenceOrient orientA,
                                          SequenceOrient orientB){
  PairOrient ret;
  int code = 0;

  if(orientA.isForward())
    code += 2;

  if(orientB.isForward())
    code += 1;

  switch(code){
    case 0: // BA_BA
      ret.setIsBA_BA(); break;
    case 1: // BA_AB
      ret.setIsBA_AB(); break;
    case 2: // AB_BA
      ret.setIsAB_BA(); break;
    case 3: // BA_AB
      ret.setIsAB_AB(); break;
    default:
      assert(0);
  }

  return ret;
}


PairOrient
ciEdgeOrientFromFragment(int             orient,
                         SequenceOrient  ciOrient,
                         SequenceOrient  mciOrient);


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

/* EdgeDegree: */
int32 EdgeDegree(GraphCGW_T *graph, EdgeCGW_T *edge);

static int isContainmentEdge(const EdgeCGW_T *edge){
  return edge->flags.bits.hasContainmentOverlap;
}

static int isMustOverlapEdge(const EdgeCGW_T *edge){
  return edge->flags.bits.mustOverlap;
}

static int isInferredEdge(const EdgeCGW_T *edge){
  return edge->flags.bits.isInferred;
}

static int isOverlapEdge(const EdgeCGW_T *edge){
  return (edge->flags.bits.hasContributingOverlap ||
          edge->flags.bits.aContainsB ||
          edge->flags.bits.bContainsA);
}

static int edgeContainsCI(const EdgeCGW_T *edge, CDS_CID_t id){
  assert((id == edge->idA) || (id == edge->idB));
  return((id == edge->idA) ? (int)edge->flags.bits.bContainsA :
	 (int)edge->flags.bits.aContainsB);
}

static int isTransChunkEdge(const EdgeCGW_T *edge){
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
#define SLOPPY_EDGE_VARIANCE_THRESHHOLD (2.5e+7)    // sigma = 5000
// #define SLOPPY_EDGE_VARIANCE_THRESHHOLD (14.0e+6)    // sigma = 3700
// #define SLOPPY_EDGE_VARIANCE_THRESHHOLD (4.0e+6)    // sigma = 2000

static int isSloppyEdge(EdgeCGW_T const *edge){
  if(edge->flags.bits.isInferred)
    return 0;
  if(edge->flags.bits.isSloppy)
    return 1;
  // for older files
  return (edge->distance.variance > SLOPPY_EDGE_VARIANCE_THRESHHOLD);
}

static int isSingletonOverlapEdge(EdgeCGW_T *edge){
  return ( isOverlapEdge(edge) &&   edge->edgesContributing <= 1);
}

static int isProbablyBogusEdge(EdgeCGW_T *edge){
  return edge->flags.bits.isProbablyBogus;
}

#define MIN_EDGES 2

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
  return(edge->edgesContributing >= MIN_EDGES);
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
    fprintf(stderr, "* NextGraphNodeIterator prev:" F_CID " curr:" F_CID " next:" F_CID "\n",
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
      isInitialized = node->flags.bits.isDead || node->flags.bits.isFree || (node->id == idx);
    }
    if(isInitialized){
      if(e->uniqueOnly){
        if(node->flags.bits.isContig){
          isUniqueEnough = (node->scaffoldID > NULLINDEX || node->flags.bits.isUnique);
          if(e->verbose & !isUniqueEnough)
            fprintf(stderr, "* Skipping contig " F_CID " since not scaffolded\n", node->id);
        }else if(node->flags.bits.isScaffold){
          isUniqueEnough = TRUE;
        }else{
          assert(node->flags.bits.isCI);
          isUniqueEnough = (node->type == DISCRIMINATORUNIQUECHUNK_CGW ||
                            node->type == UNIQUECHUNK_CGW ||
                            node->flags.bits.isUnique);
          if(e->verbose & !isUniqueEnough)
            fprintf(stderr,"* Skipping CI " F_CID " since not scaffolded\n", node->id);
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

#include "GraphEdgeIterator.H"


/****************************************************************************/

/* Diagnostic */
size_t ReportMemorySizeGraphCGW(GraphCGW_T *graph, FILE *stream);




/* Operations on Edges */

void
MergeAllGraphEdges(GraphCGW_T        *graph,
                   vector<CDS_CID_t> &rawEdges,
                   bool               includeGuides,
                   bool               mergeAll);



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
int MergeGraphEdges(GraphCGW_T *graph, vector<CDS_CID_t> &inputEdges);



void DumpGraphEdges(GraphCGW_T *graph, char *outname);

void GraphEdgeSanity(GraphCGW_T *graph, CDS_CID_t eid);

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

void InsertGraphEdgeInList(GraphCGW_T *graph, CDS_CID_t CIedgeID, CDS_CID_t sid);
void InsertGraphEdge(GraphCGW_T *graph,  CDS_CID_t cedgeID);

void PrintGraphEdge(FILE *fp, GraphCGW_T *graph, const char *label,
                    EdgeCGW_T *edge, CDS_CID_t cid);

void PrintContigEdgeInScfContext(FILE *fp, GraphCGW_T *graph, char *label, EdgeCGW_T *edge, int cid);



CDS_CID_t AddGraphEdge( GraphCGW_T *graph,
                        CDS_CID_t cidA, CDS_CID_t cidB,
                        CDS_CID_t fidA, CDS_CID_t fidB,
                        CDS_CID_t dist,
                        LengthT distance,
                        double   quality,
                        int32 fudgeDistance,
                        PairOrient orientation,
                        int isInducedByUnknownOrientation,
                        int isOverlap,
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
                                PairOrient orient);

// Find an edge..returns NULL if none found
EdgeCGW_T *FindGraphEdge(GraphCGW_T *graph,
                         CDS_CID_t idA, CDS_CID_t idB,
                         PairOrient orient);

// Find an overlap edge (assumed to exist) between the two CIs.
// Unlink it from the graph
void  DeleteGraphOverlapEdge(GraphCGW_T *graph,
                             CDS_CID_t idB, CDS_CID_t idA,
                             PairOrient orient);

// Move the edge to the free list
void FreeGraphEdge(GraphCGW_T *graph, EdgeCGW_T *edge);

// Get the edge from the free list
EdgeCGW_T *GetFreeGraphEdge(GraphCGW_T *graph);


// Clean hash table
// Move deleted nodes and edges to their respective free lists
void RecycleDeletedGraphElements(GraphCGW_T *graph);



void CheckEdgesAgainstOverlapper(GraphCGW_T *graph);



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
void
FragOffsetAndOrientation(CIFragT          *frag,
                         ChunkInstanceT   *chunk,
                         LengthT          *chunkOffset,    // output
                         SequenceOrient   *chunkOrient,    // output
                         int32            *extremal,       // output
                         int32             orientIsOpposite);

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
void
BuildGraphEdgesDirectly(GraphCGW_T *graph, vector<CDS_CID_t> &rawEdges);

// Create the raw link-based edges incident on a particular graph node
void
BuildGraphEdgesFromMultiAlign(GraphCGW_T         *graph,
                              NodeCGW_T          *node,
                              MultiAlignT        *ma,
                              GraphEdgeStatT     *stats,
                              int                 buildAll,
                              vector<CDS_CID_t>  *rawEdges);

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

void ComputeMatePairDetailedStatus(void);
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

void ReallocGraphEdges(GraphCGW_T *graph, int32 numEdges);

void InitGraphEdge(EdgeCGW_T *edge);

void AssignFragsToResolvedCI(GraphCGW_T *graph,
                             CDS_CID_t fromID, CDS_CID_t toID,
                             VA_TYPE(CDS_CID_t) *fragments);

// QC check on surrogate sanity
void CheckSurrogateUnitigs();


#endif


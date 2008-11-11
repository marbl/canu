
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

#ifndef _AS_CGB_BUBBLE_POPPER_H_
#define _AS_CGB_BUBBLE_POPPER_H_

static const char *rcsid__AS_CGB_BUBBLE_POPPER_H_ = "$Id: AS_CGB_Bubble_Popper.h,v 1.8 2008-11-11 16:22:26 brianwalenz Exp $";

/* The bubble popper will REJECT bubbles with more than this many fragments! */
#define POPPER_MAX_BUBBLE_SIZE 200

/* The parameters for the AS_ALN_affine_overlap() function (should
   eventually become a CGB parameters). */
#define POPPER_ALN_ERATE    AS_READ_ERROR_RATE
#define POPPER_ALN_THRESH   1e-6
#define POPPER_ALN_MIN_LEN  AS_OVERLAP_MIN_LEN
#define POPPER_ALN_FN       Local_Overlap_AS

#if AS_CGB_BUBBLE_DIST_OUTPUT
#define POPPER_ALN_TYPE     AS_FIND_LOCAL_ALIGN
#else
#define POPPER_ALN_TYPE     AS_FIND_LOCAL_OVERLAP
#endif


/* This provides a minimum discriminator score, below which bubbles
   are rejected. */
#define POPPER_MIN_DISCRIMINATOR  0.0

typedef struct BubblePopper {

  /* Permanent information fields. */

  BubGraph_t bg;
  TChunkMesg *chunks;
  TChunkFrag *chunkFrags;
  float globalArrivalRate;
  GateKeeperStore *gkpStore;

  /* Permanent statistics fields. */

  int numBubblesProcessed;
  int numBubblesCollapsed;
  int numOlapsComputed;
  int numOlapsSuccessful;
  int numOlapsRetained;
  int numFragsInBubbles;
  int numFragsInCollapsedBubbles;
  int totalDistSpannedByBubbles;
  int totalDistSpannedByCollapsedBubbles;
  int numRejectedByDiscriminator;

  FILE *nfragsFile;		/* Distribution output files */
  FILE *nfragspopFile;
  FILE *sdiffFile;

  /* Recycled scratch space for each bubble. */

  fragRecord  rsp;
  IntFragment_ID curBubSize;
  int *topDistArray;		/* Scratch space for longest path algorithm */
  BG_E_Iter *dfsStack;	        /* Scratch space for traversals (dfs and
				   longest path) */
  int *vidToBid;	        /* Maps VIDs to IDs in the bubFrags array */
  IntFragment_ID *bubFrags;	/* Fragments in bubble */
  InternalFragMesg *bubMesgs;	/* Frags in message form for aligner */
  int *adj;			/* Adjacency array */
  int numOlaps;
  OverlapMesg *bubOlaps;        /* Temporary array of generated overlaps */
  int allocatedByCreate;	/* True if create method was used for
				   allocation */
} BubblePopper;


/*
 * METHODS
 */

/* Performs the bulk of the work of initializing a BubblePopper (all
   the create method in "AS_CGB_Bubble.h" really does is generate bg). */
void
BP_init(BubblePopper_t bp, BubGraph_t bg, TChunkMesg chunks[],
	TChunkFrag cfrgs[], float global_arrival_rate,
        GateKeeperStore *gkpStore,
        const char * fileprefix
        );

/* Resets the scratch space, for use before processing a new bubble. */
void
BP_reset(BubblePopper_t bp);

/* Returns the Chunk Array */
TChunkMesg *
BP_chunks(BubblePopper_t bp);

/* Returns the ChunkFrag Array */
TChunkMesg *
BP_chunkFrags(BubblePopper_t bp);

/* Returns the Chunk with the specific ID. */
AChunkMesg *
BP_getChunk(BubblePopper_t bp, IntChunk_ID cid);

/* Returns the number of fragments in the bubble. */
IntFragment_ID
BP_numFrags(BubblePopper_t bp);

/* Returns the VID of fragment with BID bid in the bubble.
   Fragments are sorted in increasing order by DFS coordinate. */
IntFragment_ID
BP_getFrag(BubblePopper_t bp, IntFragment_ID bid);

/* Given a VID, returns the BID of the fragment.  Undefined if the fragment
   is not in the bubble. */
int
BP_getBID(BubblePopper_t bp, IntFragment_ID vid);

/* Tries to find an overlap between two fragments.  If an overlap is found,
   then TRUE is returned, the overlap is placed in the bubOlaps array, and
   the adjacency array is updated.  Otherwise, FALSE is returned. */
int
BP_findOverlap(BubblePopper_t bp, IntFragment_ID bid1, IntFragment_ID bid2);

/* Operations on the directed adjacency array.  An Entry (V1, V2)
   indicates an from V1 to V2.  Also, by default, (V1, V1) is 1 even
   though the graph has no loops.

   An entry of 0 indicates no edge.
   1 indicates a direct edge.
   2 indicates transitive connectivity. */
int
BP_getAdj(BubblePopper_t bp, IntFragment_ID bid1, IntFragment_ID bid2);

int
BP_getAdj_VID(BubblePopper_t bp, IntFragment_ID vid1, IntFragment_ID vid2);

void
BP_setAdj(BubblePopper_t bp, IntFragment_ID bid1, IntFragment_ID bid2, int v);

void
BP_setAdj_VID(BubblePopper_t bp, IntFragment_ID vid1, IntFragment_ID vid2,
	      int v);

#endif

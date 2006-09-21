
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

#ifndef FRAGMENTPLACEMENT_H
#define FRAGMENTPLACEMENT_H

#define MAX_SIGMA_SLOP 3.0

#define EXTERNAL_MATES_ONLY 1
#define ALL_MATES 0

/*  my_getChunkInstanceID() is a clone of  pushdgetChunkInstanceID() from AS_REZ/PlaceFragments.c */
int my_getChunkInstanceID(ChunkInstanceT *chunk, int index);


void GetRelationBetweenFragsOnChunks(CIFragT *frag,ChunkInstanceT*fragChunk,
				     CIFragT *mate,ChunkInstanceT*mateChunk,
				     LengthT*separation,OrientType* edgeOri);


/* FragAndMateAreCompatible() is based partly on PlaceFragmentBasedOnMate()
   from AS_REZ/PlaceFragments.c */
int FragAndMateAreCompatible(CIFragT *frag, ChunkInstanceT *fragChunk, 
			     CIFragT *mate, ChunkInstanceT *mateChunk,
			     OrientType expectedOri);




// stuff that was in localUnitigging initially

// new iterator for looping over the fragments in a chunk -- note need to call cleanup routine when done

typedef struct FragIterator{
  /* we need to have:
     an ma for the  top chunk
     if (type == contig) 
     an iterator over chunks in top ma
     an ma for sub chunks
     an iterator over fragments in the (sub) ma

     the fragment iterator should be set to go after init

     each call to iterator should 
     compute next in fragment iterator
     if non null, return
     else
     iterate to next sub chunk
     if null, return failure
     else 
     restart fragment iterator     
     compute next


     return the current
  */

  uint32 isUtg:1;
  uint32 isCtg:1;
  uint32 includeSurrogateFrgs:1;
  CDS_CID_t id;
  MultiAlignT *ma;
  cds_int32 fOrder; /* which fragment within multialign is next */
  ContigTIterator subchunks;
  struct FragIterator *subchunkIterator;

} CGWFragIterator;

void InitCIFragTInChunkIterator(CGWFragIterator* frags,NodeCGW_T *chunk, int includeSurrogates);
int NextCIFragTInChunkIterator(CGWFragIterator* frags, CIFragT**nextfrg);
void CleanupCIFragTInChunkIterator(CGWFragIterator* frags);




// New iterator for finding the mates of a fragment

typedef struct  {
  CDS_CID_t thisFragIID;
  CDS_CID_t nextLink;
  NodeCGW_T *node;
  GateKeeperLinkRecordIterator GKPLinks;
  uint32 external_only:1;
  uint32 getLinksFromStore:1;
} CGWMateIterator;

// set up iterator
void InitCGWMateIterator(CGWMateIterator* mates,CDS_CID_t fragIID, int external_only,ChunkInstanceT *node);
// determine whether a fragment is internal to a node
static int check_internal(NodeCGW_T *node,CDS_CID_t frgIID);
// return the next link (external only as appropriate)
int NextCGWMateIterator(CGWMateIterator* mates,CDS_CID_t * linkIID);


// returns true if frg has a mate in scaffold sid
int matePlacedIn(CIFragT *frg, CDS_CID_t sid);

// returns true if all of frg's placed mates are in scaffold sid; mate and mateChunk point to 
// the last found mate and its chunk
int matePlacedOnlyIn(CIFragT *frg, CDS_CID_t sid, CIFragT **mate, ChunkInstanceT **mateChunk);

// returns the scaffold containing a given fragment
int scaffoldOf(CDS_CID_t fiid);

// append f_list to a  multialignt's f_list and update seqdb entry
void PlaceFragmentsInMultiAlignT(CDS_CID_t toID, int isUnitig,
				 VA_TYPE(IntMultiPos) *f_list);

// based on AssignFragsToResolvedCI() from GraphCGW_T; make adjustments 
// to everything  to reflect assignment of fragments to a surrogate CI
void ReallyAssignFragsToResolvedCI(GraphCGW_T *graph,
				   CDS_CID_t fromID, CDS_CID_t toID,
				   VA_TYPE(CDS_CID_t) *fragments);




#endif


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

   The function LookupQualityOverlap takes the same arguments as LookupOverlap
   with an additional argument for the quality function and the quality value.
   It first looks up the overlap
   in the hash table. If it has the appropriate quality bit set the quality is returned.
   If not, the quality is computed using the appropriate quality function and stored
   in the hash table. In addition, the appropriate bit is set indicating that computation.


*/


#ifndef CHUNKOVERLAP_H
#define CHUNKOVERLAP_H

static const char *rcsid_CHUNKOVERLAP_H = "$Id: ChunkOverlap_CGW.h,v 1.10 2008-10-08 22:02:55 brianwalenz Exp $";

#include "AS_CGW_dataTypes.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_Hash.h"
#include "MultiAlignStore_CNS.h"

/* this enum indicates with which quality function the quality of an overlap
   should be computed */
typedef enum { BAYESIAN } QualityFuncT;


/* This structure comprises the 'symbol' that is the key for database lookups */
typedef struct {
  CDS_CID_t cidA;
  CDS_CID_t cidB;
  ChunkOrientationType orientation;
}ChunkOverlapSpecT;

/* This is the value stored in the symbol table.  Note that it's first field is the key
   for the symbol table.  Thus we kill two birds with one stone.  These structures are allocated
   from a heap */

typedef struct {
  ChunkOverlapSpecT spec;

  CDS_COORD_t  minOverlap;
  CDS_COORD_t  maxOverlap;
  CDS_COORD_t  cgbMinOverlap;
  CDS_COORD_t  cgbMaxOverlap;
  float        errorRate;

  // This is what we found
  uint32       computed:1;
  uint32       fromCGB:1;
  uint32       hasBayesianQuality:1;
  uint32       AContainsB:1;
  uint32       BContainsA:1;
  uint32       suspicious:1;
  uint32       unused:26;

  CDS_COORD_t  overlap;  // The overlaplength if there is an overlap, 0 otherwise
  CDS_COORD_t  ahg;
  CDS_COORD_t  bhg;
  float        quality;
  CDS_COORD_t  min_offset;
  CDS_COORD_t  max_offset;
}ChunkOverlapCheckT;

typedef struct {
  HashTable_AS *hashTable;
  Heap_AS      *ChunkOverlaps;  //  Heap of ChunkOverlapCheckT
}ChunkOverlapperT;

// Constructor
ChunkOverlapperT *CreateChunkOverlapper(void);

// Destructor
void DestroyChunkOverlapper(ChunkOverlapperT *chunkOverlapper);


// Save to File
void  SaveChunkOverlapperToStream(ChunkOverlapperT *chunkOverlapper,
                                  FILE *stream);

// Load from File
ChunkOverlapperT *LoadChunkOverlapperFromStream(FILE *stream);


int InitCanonicalOverlapSpec(CDS_CID_t cidA, CDS_CID_t cidB,
                             ChunkOrientationType orientation,
                             ChunkOverlapSpecT *spec);

ChunkOverlapCheckT *LookupCanonicalOverlap(ChunkOverlapperT *chunkOverlapper,
                                           ChunkOverlapSpecT *spec);

int InsertChunkOverlap(ChunkOverlapperT *chunkOverlapper,
                       ChunkOverlapCheckT *olap);

#endif


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

        Module:  FbacREZ.h

   Description:  Declaration of common data types for FbacREZ.c

    Programmer:  M. Flanigan

       Written:  May 4, 2000

  Last Revised:  

 **********************************************************************/

/*********************************************************************
   CVS_ID: $Id: FbacREZ.h,v 1.3 2005-03-22 19:07:35 jason_miller Exp $
 *********************************************************************/

#ifndef FBAC_REZ_H
#define FBAC_REZ_H

#define POSDIR 1
#define NEGDIR -1

#define LEFTEND 0
#define RIGHTEND 1

#define MAX_SCAFFOLDS_CHECKED 1024*1024

#include "InputDataTypes_CGW.h"  /* for LengthT */

typedef enum {
  BAC_CONSISTENT,
  BAC_INCONSISTENT,
  NO_BACS
} BacStatusREZ;

int contains_fbac(NodeCGW_T* n);
static int buildLocale;

/* returns   
   BAC_CONSISTENT,
   if the node (Chunk instance or Contig) 
   is correctly covered by bac fragments, 
   BAC_INCONSISTENT if there is a discrapency
   NO_BACS if the chunk has no NO_BACS at all */
BacStatusREZ isFbacConsistent(NodeCGW_T* node);

// returns TRUE if the edge connects two chunks that fullfill
// is_suitable_fbac_chunk (e.g. is_covered_by_fbac_ or contains_fbac
// and if the indicated overlap of this edge is consistent
// with the locpos fields of the chunks
int edge_is_fbac_consistent(EdgeCGW_T* n);


typedef struct fragmentInfo 
{
	  int fragIid;
	  int frag5pContigPos;
	  int frag5pLocalePos;
} fragmentInfoT;  

typedef struct localeInfo 
{
	  int localeNumber;
	  int32 extremalFragID;
	  int32 extremalFragIid;
	  int contigID;
	  int extremalFragUnitigID;    // this is the unitig (possibly a surrogate) used to get at extremalFragID
	  int basesFromEnd;  // LengthT fragContigOffset5p;
	  struct localeInfo * next;        // we link structs together in getLocalesInNode
	  int inSurrogate;                 // set if extremal fragment in in surrogate
	  int fragContigOffset5p, fragContigOffset3p;
	  int BacOffset;
	  int numLocaleFragments;
	  int iidTrend;    // tracks whether the iids are increasing or decreasing towards end of contig/scaffold
} localeInfoT;  

typedef struct ChunkInsertInfo 
{
	  int32 scaffoldID;
	  int32 contigID;
	  LengthT aEndOffset, bEndOffset;
	  int insertOrder;
	  struct ChunkInsertInfo *next;
	  int hasMoved;
} ChunkInsertInfoT;  

typedef struct GapInfo
{
	  int localeNumber;
	  int fragIDLeftContig;           // the id of the frag on the left side of the gap
	  int fragIDRightContig;          // the id of the frag on the right side of the gap
	  int fragLeftContigID;           // the contigID for the frag on the left side of the gap
	  int fragRightContigID;          // the contigID for the frag on the right side of the gap
	  int fragOriginalLeftContigID;   // the original contigID for the frag on the left side of the gap
	  int fragOriginalRightContigID;  // the original contigID for the frag on the right side of the gap
	  int unitigIDLeftContig;
	  int unitigIDRightContig;	  
	  float distanceDifference;       // the difference of the next two fields
	  float fragSeparationScaffold;   // the distance the extremal frags are apart in the scaffold
	  float fragSeparationLocale;     // the distance the extremal frags are apart in the locale
	  int fragInSurrogateLeftContig;
	  int fragInSurrogateRightContig;	  
	  int basesFromEndLeftContig;
	  int basesFromEndRightContig;
	  int fragLeftContigOffset5p;
	  int fragLeftContigOffset3p;	  
	  int fragRightContigOffset5p;
	  int fragRightContigOffset3p;
	  double fragOriginalLeftContigOffset5p;
	  double fragOriginalLeftContigOffset3p;	  
	  double fragOriginalRightContigOffset5p;
	  double fragOriginalRightContigOffset3p;
} GapInfoT;

typedef struct ScaffoldInfo
{
	  int numContigs;
	  int *ContigIDs;
	  LengthT *AEndPositions;
	  LengthT *BEndPositions;
} ScaffoldInfoT;

typedef struct unlinkedEdge 
{
	  int eid;
	  CIEdgeT *edge;
	  struct unlinkedEdge *next;
} unlinkedEdgeT;

static unlinkedEdgeT *edgesUnlinkedList = NULL, *edgesUnlinkedListHead = NULL;

ChunkInsertInfoT* BuildLocaleOverlaps(GapInfoT *gapInfo,
                                      ChunkInstanceT *lchunk,
                                      ChunkInstanceT *rchunk,
                                      int *completeOverlapPath,
                                      double *rchunkNewAEnd,
                                      double *rchunkNewBEnd,
                                      ChunkOrientationType *tempOlapOrientation,
                                      Overlap* rchunkOverlap,
                                      int fixLocaleOverlaps);

void SortCommonLocales( int numCommonLocales,
                        GapInfoT **gapInfoArray,
                        int intraScaffoldGap,
                        LengthT gapEstimate);

int getLocalesInNode(NodeCGW_T* node, localeInfoT** localeInfo,
                     int end, unsigned int basesFromEndCutoff);

int findLastLocaleFragInContig( int contigID, unsigned int currFragIid,
                                unsigned int endFragIid, int increment, 
                                int *currExtremeUnitigID);

int FragOutOfBounds(CIFragT* frag, int unitigID);

void CheckOrientation( ContigT* lchunk, ContigT* rchunk,
                       LengthT* rchunk_delta, 
                       ChunkOrientationType olapOrientation,
                       Overlap* rchunkOverlap);

void InsertWalkedChunks( ChunkInsertInfoT* chunksWalked,
                         ChunkInstanceT* lchunk, ChunkInstanceT* rchunk, 
                         int completeOverlapPath,
                         int *checkScaffold, int *checkScaffoldCount );

void MergeWalkedScaffolds( CIScaffoldT* scaff, CIScaffoldT* nextScaff, 
                           double nextScaffExtremeContigAEnd,
                           double nextScaffExtremeContigBEnd,
                           ContigT *scaffExtremeContig,
                           ContigT *nextScaffExtremeContig);

void printScaffoldContigs(CIScaffoldT *scaffold);

int ComputeAhangSize( CIFragT *currFrag, ChunkInstanceT *contigCurrFrag, 
                      CIFragT *fragSucc, ChunkInstanceT *contigSucc,
                      int overlapOrientation);

void dumpContigInfo(ChunkInstanceT *contig);

void LinkLocaleOverlaps(void);

void UnlinkLocaleOverlaps(int locale);

void DumpLocaleOverlaps(void);

void GetFragmentPositionInScaffoldFromContig(CIFragT *frag,
                                             int *left_end, int *right_end, 
                                             int *fragmentScaffoldOrientation, 
                                             int contigLeftEnd,
                                             int contigRightEnd,
                                             int contigScaffoldOrientation);

void GetFragmentPositionInScaffold(CIFragT *frag,
                                   int *left_end, int *right_end, 
                                   int *fragmentScaffoldOrientation);

void GetFragmentPositionInContigFromChunk(CIFragT *frag,
                                          int *left_end, int *right_end, 
                                          int *fragmentContigOrientation, 
                                          int chunkContigOrientation,
                                          int chunkLeftEnd, int chunkRightEnd);

void GetCIPositionInScaffold(ChunkInstanceT *CI,
                             int *left_end, int *right_end, 
                             int *CIScaffoldOrientation);

void GetFragmentPositionInScaffoldFromCI(CIFragT *frag,
                                         int *left_end, int *right_end, 
                                         int *fragmentScaffoldOrientation, 
                                         int CILeftEnd, int CIRightEnd,
                                         int CIScaffoldOrientation);

void localeCam(char *);

void localeSimCam();

void GetContigPositionInScaffold(ChunkInstanceT *contig,
                                 int *left_end, int *right_end, 
                                 int *contigScaffoldOrientation);

void AdjustRightContigPos( ContigT* lchunk, ContigT* rchunk,
                           LengthT rchunk_delta, CIScaffoldT* scaff);

int FindCommonLocales( ContigT * lcontig, int lcontigGapEnd, 
                       ContigT * rcontig, int rcontigGapEnd,
                       GapInfoT **gapInfoArray);

Overlap* OverlapContainingContigs(CIFragT *frag1, CIFragT *frag2,
                                  ChunkOrientationType *overlapOrientation,
                                  ChunkInstanceT *leftContig,
                                  ChunkInstanceT *rightContig);

int compFragments( const void *s1, const void *s2);

ChunkInsertInfoT* SortWalk( ChunkInsertInfoT* walkedChunks,
                            ChunkInsertInfoT** firstWalkedChunk,
                            ChunkInsertInfoT** lastWalkedChunk);

void SaveContigInformation( NodeCGW_T *scaffold,
                            ScaffoldInfoT *scaffoldContigs);

void RestoreContigInformation( NodeCGW_T *scaffold,
                               ScaffoldInfoT *scaffoldContigs);

void FreeContigInformation( ScaffoldInfoT *scaffoldContigs);

int GetChunkPositionInContig(ChunkInstanceT *chunk,
                             int *left_end, int *right_end, 
                             int *chunkContigOrientation);

void trimScaffoldEnds( CIScaffoldT *scaff );

void walkScaffolds( CIScaffoldT *scaff1,
                    int invertScaff1, CIScaffoldT *scaff2,
                    int invertScaff2, double gapSize);

void	ResetFragPositions( int numCommonLocales, GapInfoT **gapInfoArray);

void dumpGapInfo(ChunkInstanceT *leftContig, ChunkInstanceT *rightContig);

void SetFragPositions( int numCommonLocales, GapInfoT **gapInfoArray);

void removeSmallContigs( CIScaffoldT *scaff );

int CheckRchunkContainment( ChunkInsertInfoT *walkedChunks,
                            ChunkInsertInfoT *lastWalkedChunk);

int CheckWalkOverlaps( ChunkInsertInfoT *walkedContigs,
                       int numChunks, int fixOverlapFailures);

#endif

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

#undef DEBUG_ECR
#undef DIAG_PRINTS

//#error I'm probably not opening the frag store for read/write, see the XXX

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "FbacREZ.h"
#include "PublicAPI_CNS.h"
#include "AS_PER_SafeIO.h"
#include "ChiSquareTest_CGW.h"

#define NUM_STDDEV_CUTOFF 5.0
#define CONTIG_BASES 1000
#define MAX_EXTENDABLE_FRAGS 100

typedef struct extendableFragT {
  int fragIid;
  int extension;
  int addedBases;
  int basesToNextFrag;
  int fragOnEnd;
  int unitigID;
} extendableFrag;

typedef struct fragPositionsT {
  int bgn;
  int end;
} fragPositions;

int findFirstFrag(ContigT *contig, int *fragIid,
                  int *extension, int *basesToNextFrag);
int findLastFrag(ContigT *contig, int *fragIid,
                 int *extension, int *basesToNextFrag);
int dumpFragInfo(int fragIid);
int examineGap(ContigT *lcontig, int lFragIid,
               ContigT *rcontig, int rFragIid, 
               int gapNumber, int *ahang,
               int* olapLengthOut, int *bhang, int *currDiffs,
               int *lcontigBasesIntact, int *rcontigBasesIntact,
               int *closedGapDelta,
               int lBasesToNextFrag, int rBasesToNextFrag,
               int *leftFragFlapLength, int *rightFragFlapLength);
void initVarArrays(void);
void SequenceComplement(char *sequence, char *quality);
LengthT FindGapLength(ChunkInstanceT * lchunk,
                      ChunkInstanceT * rchunk, int verbose);
int findFirstUnitig(ContigT *contig, int *unitigID);
int findLastUnitig(ContigT *contig, int *unitigID);
void collectContigStats(int *allContigLengths, int *contigValid);
void produceContigStatsConventional(void);
void produceGapStats(int numGaps, int *closedGap,
                     int *closedGapDelta, int *originalGaps);
void produceScaffoldStats(int *alteredScaffoldLengths);
void produceContigStats(int numGaps, int *closedGap, int *closedGapDelta,
                        int *lcontigIdGap, int *rcontigIdGap,
                        int *lcontigLength, int *rcontigLength,
                        int *contigValid, int *allContigLengths,
                        int numDeletedContigs);
int findRightExtendableFrags(ContigT *contig, int *basesToNextFrag,
                             extendableFrag *extFragsArray);
void extendUnitigs(NodeCGW_T *unitig, int fragIid,
                   extendableFrag extFrag, int scaffoldLeftEnd);
void extendContig(ContigT *contig, int extendAEnd);
int findFirstExtendableFrags(ContigT *contig, extendableFrag *extFragsArray);
int findLastExtendableFrags(ContigT *contig, extendableFrag *extFragsArray);
int GetNewUnitigMultiAlign(NodeCGW_T *unitig, fragPositions *fragPoss,
                           int extendedFragIid);
int setCgwClearRange(int fragIid, int frag3pDelta);
// int getAlteredFragPositions(NodeCGW_T *unitig, fragPositions *fragPoss,
//                             int alteredFragIid, int extension);
void getAlteredFragPositions(NodeCGW_T *unitig, fragPositions **fragPoss,
                             int alteredFragIid, int extension);
void bubbleSortIUMs(IntMultiPos *f_list, int numIMPs);
void adjustUnitigCoords(NodeCGW_T *contig);

void DumpContigMultiAlignInfo (char *label, MultiAlignT *cma, int contigID);
void DumpUnitigInfo(char *label, NodeCGW_T *unitig);
void DumpContigUngappedOffsets(char *label, int contigID);

void leftShiftIUM(IntMultiPos *f_list, int numFrags, int extendedFragIid);
void rightShiftIUM(IntMultiPos *f_list, int numFrags, int extendedFragIid);
void saveFragAndUnitigData(int lFragIid, int rFragIid);
void restoreFragAndUnitigData(int lFragIid, int rFragIid);

void printGapSizes();
int writeEcrCheckpoint(int *numGapsInScaffold,
                       int *numGapsClosedInScaffold,
                       int *numSmallGapsInScaffold,
                       int *numSmallGapsClosedInScaffold,
                       int *numLargeGapsInScaffold,
                       int *numLargeGapsClosedInScaffold);
int loadEcrCheckpoint(int ckptNum, int *numGapsInScaffold,
                      int *numGapsClosedInScaffold,
                      int *numSmallGapsInScaffold,
                      int *numSmallGapsClosedInScaffold,
                      int *numLargeGapsInScaffold,
                      int *numLargeGapsClosedInScaffold);
int extendCgwClearRange(int fragIid, int frag3pDelta);
void SynchUnitigTWithMultiAlignT(NodeCGW_T *unitig);
int revertToCnsClearRange(int fragIid);


extern int                     totalContigsBaseChange;
extern ReadStructp             fsread;
extern VA_TYPE(char)          *lContigConsensus;
extern VA_TYPE(char)          *rContigConsensus;
extern VA_TYPE(char)          *lContigQuality;
extern VA_TYPE(char)          *rContigQuality;
extern VA_TYPE(char)          *reformed_consensus;
extern VA_TYPE(char)          *reformed_quality;
extern VA_TYPE(int32)         *reformed_deltas;
extern VA_TYPE(IntElementPos) *ContigPositions;
extern VA_TYPE(IntElementPos) *UnitigPositions;

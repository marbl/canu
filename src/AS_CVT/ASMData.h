
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
#ifndef ASMDATA_H
#define ASMDATA_H

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_PER_asmStore.h"

typedef enum
{
  ASM_GnuplotOutput,
  ASM_CliqueFinderOutput
} ASM_OutputType;


typedef enum
{
  ASM_Stretched = 0,
  ASM_Compressed,
  ASM_Normal,
  ASM_Antinormal,
  ASM_OuttieLR,        // breakpoint pair to left & right of outties
  ASM_OuttieM,         // region between outtie matepairs
  NumBreakpointTypes
} BreakpointType;


typedef struct
{
  Fragment_ID leftUID;
  SeqInterval pts[4];
} ASM_Quad;
VA_DEF(ASM_Quad)


typedef struct
{
  Fragment_ID leftUID;
  Fragment_ID rightUID;
  Dist_ID     distUID;
  Scaffold_ID containerUID;
  SeqInterval fivePrimes;
  CDS_COORD_t basePairs;
  float       stddevs;
  OrientType  orient;
} ASM_MatePair;
VA_DEF(ASM_MatePair)

typedef struct
{
  VA_TYPE(ASM_MatePair) * innie;
  VA_TYPE(ASM_MatePair) * normal;
  VA_TYPE(ASM_MatePair) * antinormal;
  VA_TYPE(ASM_MatePair) * outtie;
  VA_TYPE(ASM_MatePair) * elsewhere;
  VA_TYPE(ASM_MatePair) * inDegenerate;
  VA_TYPE(ASM_MatePair) * inSurrogate;
  VA_TYPE(ASM_MatePair) * inChaff;
  VA_TYPE(ASM_MatePair) * missing;
} CloneData;
typedef CloneData * CloneDatap;


AssemblyStore * CreateAssemblyStoreFromASMFile(FILE * fi,
                                               char * storePath,
                                               char * gkpStorePath,
                                               char * frgStorePath);
void PrintFragmentScaffoldCoordinates(AssemblyStore * asmStore,
                                      int doInstances,
                                      int doSingleSurrogates,
                                      int doDegenerates,
                                      int doChaff,
                                      int doLinks,
                                      int doUnreferenced,
                                      FILE * fo);
void PrintReadsPlaced(AssemblyStore * asmStore,
                      int doInstances,
                      int doSingleSurrogates,
                      int doDegenerates,
                      int doChaff,
                      int doUnreferenced,
                      FILE * fo);

void DeleteCloneData(CloneData * cd);
CloneData * CreateCloneData(void);
CloneData * GetScaffoldCloneData(AssemblyStore * asmStore,
                                 int32 scaffoldIndex,
                                 int canonicalOnly, int isDegenerate);
void PrintCloneData(Scaffold_ID containerUID,
                    CloneData * cd, char * which, FILE * fo);
void PrintAsmStoreIntraTampaData(Scaffold_ID containerUID,
                                 CloneData * cd,
                                 int includeSPS,
                                 FILE * fo);
void PrintChromosomeElsewheres(AssemblyStore * asmStore,
                               MapStore * mapStore,
                               Scaffold_ID containerUID,
                               FILE * fo);
void PrintScaffoldElsewheres(AssemblyStore * asmStore,
                             int32 scaffoldIndex,
                             int canonicalOnly,
                             int isDegenerate,
                             FILE * fo);
void PrintBadCloneData(Scaffold_ID containerUID,
                       CloneData * cd, FILE * fo);
void PrintMissingCloneData(Scaffold_ID containerUID,
                           CloneData * cd, FILE * fo);

void CopyGateKeeperLNKStore(AssemblyStore * asmStore);
int SetZeroInstanceFragments(AssemblyStore * asmStore);
void SetMultiInstanceFragments(AssemblyStore * asmStore);
void PrintSurrogateCoordinates(AssemblyStore * asmStore, FILE * fo);
void PrintSurrogateSequenceCoordinates(AssemblyStore * asmStore, FILE * fo);
void PrintUnreferencedFrags(AssemblyStore * asmStore, FILE * fo);
void PrintMissingFrags(AssemblyStore * asmStore, FILE * fo);
void PrintDeletedFrags(AssemblyStore * asmStore, FILE * fo);
void PrintScaffoldLength(AssemblyStore * asmStore, int32 index, FILE * fo);
void PrintDegenerateLength(AssemblyStore * asmStore, int32 index, FILE * fo);
void PrintQuadrilaterals(VA_TYPE(ASM_Quad) * quads,
                         ASM_OutputType ot,
                         FILE * fo);
void PrintLibraries(AssemblyStore * asmStore, FILE * fo);
VA_TYPE(ASM_Quad) * IdentifyBadMateQuads(AssemblyStore * asmStore,
                                         CloneData * cd,
                                         char * fragTypes,
                                         BreakpointType bpType,
                                         float32 numStddevs);
void CreateBMPFilename(char * fname,
                       CDS_UID_t uid,
                       char * fragTypes,
                       BreakpointType bpType,
                       float32 numStddevs);
void PrintScaffoldContigCoordinates(AssemblyStore * asmStore,
                                    uint32 index,
                                    int doGnuplotOutput,
                                    int doAssemblyCoords,
                                    FILE * fo);

MapStore * CreateMapStoreFromFiles(AssemblyStore * asmStore,
                                   FILE * chromFp,
                                   FILE * fragFp,
                                   char * mapStorePath);
void PrintMapStore(AssemblyStore * asmStore,
                   MapStore * mapStore,
                   FILE * fo);
CloneData * GetChromosomeCloneData(AssemblyStore * asmStore,
                                   MapStore * mapStore,
                                   Scaffold_ID containerUID);
void WriteBinaryCloneData(Scaffold_ID containerUID, CloneData * cd);

void PrintIntraCloneData(Scaffold_ID containerUID, CloneData * cd, FILE * fo);

void PrintATACScaffoldGenomicAxis(AssemblyStore * asmStore,
                                  int32 index,
                                  char * parent,
                                  CDS_COORD_t * offset,
                                  FILE * fo);
void PrintDeflineATACAxes(AssemblyStore * asmStore,
                          FILE * dfp,
                          char * parent,
                          FILE * fo);
void PrintATACSurrogates(AssemblyStore * asmStore,
                         char * parent,
                         FILE * fo);
void InitializeATACFile(AssemblyStore * asmStore, FILE * fo);
void PrintBACFragCloneIntervals(AssemblyStore * asmStore,
                                MapStore * mapStore,
                                Scaffold_ID containerUID);
void PrintFastaFragmentCoordinates(AssemblyStore * asmStore, FILE * fo);

#endif // ASMDATA_H

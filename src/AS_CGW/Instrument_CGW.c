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
static char CM_ID[] = "$Id: Instrument_CGW.c,v 1.26 2007-08-18 11:42:07 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "AS_global.h"
#include "Instrument_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "FbacREZ.h"
#include "UtilsREZ.h"
#include "GraphCGW_T.h"

//#define DEBUG
//#define DEBUG_BRK
//#define DEBUG2
//#define LIST_TERMINAL_TYPES

#ifdef LIST_TERMINAL_TYPES
CDS_CID_t   ContigFirstFragIID;
FragType    ContigFirstFragType;
CDS_CID_t   ContigLastFragIID;
FragType    ContigLastFragType;
CDS_COORD_t ContigLastBP;
CDS_COORD_t UnitigOffset;
#endif

#define READ_TRIM_BASES                 30
#define INITIAL_CONTIGS_PER_SCAFFOLD  1000

//#define INSTRUMENT_CUTOFF         3.0
#define INSTRUMENT_CUTOFF         CGW_CUTOFF

#define INTERVAL(a)  ((a)->mu + INSTRUMENT_CUTOFF * (a)->sigma)


// should we print mate info for all clones or only long ones?
#define USE_LONG_MATES 0
#define USE_ALL_MATES 1

int do_surrogate_tracking=1;

// if printMateUIDs==1, labels of celamy lines for mate pairs will contain UIDs, not IIDs
int printMateUIDs=0;

typedef struct
{
  CDS_CID_t   iid;
  CDS_COORD_t length;
} IID_Size;
VA_DEF(IID_Size);

#undef DUMP_MATE_PAIRS 
//#define DUMP_MATE_PAIRS
  
#define TOLERANCE  (0.1f)
  
int DoSimpleScaffoldChecks(FILE * fp,
                           ScaffoldGraphT * graph,
                           CIScaffoldT * scaff)
{
  CIScaffoldTIterator ciIterator;
  ChunkInstanceT * ci;
  float   lastEnd = 0;
  int32 ciCount = 0;
  CDS_CID_t lastCIID = NULLINDEX;
  
  InitCIScaffoldTIterator(graph, scaff, TRUE, FALSE, &ciIterator);
  for(ciCount = 0, ci = NextCIScaffoldTIterator(&ciIterator);
      ci != NULL;
      ciCount++, ci = NextCIScaffoldTIterator(&ciIterator))
    {
      float   minCoord = MIN(ci->offsetAEnd.mean, ci->offsetBEnd.mean);
      float   maxCoord = MAX(ci->offsetAEnd.mean, ci->offsetBEnd.mean);
      float   minVariance =
        MIN(ci->offsetAEnd.variance,
            ci->offsetBEnd.variance);
      float   maxVariance =
        MAX(ci->offsetAEnd.variance,
            ci->offsetBEnd.variance);
    
      if(ciCount == 0)
        {
          if(minCoord > TOLERANCE || minCoord < -TOLERANCE)
            {
              fprintf(fp, "scf " F_CID ": 1st contig (id " F_CID ") min offset is %f\n",
                      scaff->id, ci->id, minCoord);
            }
        }
      else if(minCoord < TOLERANCE)
        {
          fprintf(fp, "scf " F_CID ": non-1st CI (" F_CID ") min offset is %f\n",
                  scaff->id, ci->id, minCoord);
        }
      else if(minVariance < TOLERANCE)
        {
          fprintf(fp, "scf " F_CID ": CI (" F_CID ") variance (%f) is zero or negative\n",
                  scaff->id, ci->id, minVariance);
        }

      if(ciCount == scaff->info.Scaffold.numElements)
        {
          if(maxCoord > scaff->bpLength.mean + TOLERANCE ||
             maxCoord < scaff->bpLength.mean - TOLERANCE)
            {
              fprintf(fp, "scf " F_CID " length (%f) isn't last CI's (" F_CID ") end (%f)\n",
                      scaff->id, scaff->bpLength.mean,
                      ci->id, maxCoord);
            }
          if(maxVariance > scaff->bpLength.variance + TOLERANCE ||
             maxVariance < scaff->bpLength.variance - TOLERANCE)
            {
              fprintf(fp,
                      "scf " F_CID " variance (%f) isn't last CI's (" F_CID ") variance (%f)\n",
                      scaff->id, scaff->bpLength.variance,
                      ci->id, maxVariance);
            }
        }
      else if(minVariance < 1.f + TOLERANCE &&
              ((ci->offsetAEnd.variance < 1.f + TOLERANCE &&
                ci->offsetAEnd.mean > TOLERANCE) ||
               (ci->offsetBEnd.variance < 1.f + TOLERANCE &&
                ci->offsetBEnd.mean > TOLERANCE)))
        {
          fprintf(fp, "scf " F_CID ": CI (" F_CID ") variance (%f) is less than 1\n",
                  scaff->id, ci->id, minVariance);
        }

      if(minCoord < lastEnd - 20.f - TOLERANCE)
        {
          fprintf(fp,
                  "scf " F_CID ": large negative gap %f between CIs " F_CID " and " F_CID "\n",
                  scaff->id, (minCoord - lastEnd), lastCIID, ci->id);
        }
    
      if(maxCoord > scaff->bpLength.mean + TOLERANCE)
        {
          fprintf(fp,
                  "scf " F_CID ": CI (" F_CID ") end (%f) is past scaffold end (%f)\n",
                  scaff->id, ci->id, maxCoord, scaff->bpLength.mean);
        }

      if(maxVariance > scaff->bpLength.variance + TOLERANCE)
        {
          fprintf(fp,
                  "scf " F_CID ": CI (" F_CID ") position variance (%f) exceeds scaffold length variance (%f)\n",
                  scaff->id, ci->id, maxVariance, scaff->bpLength.variance);
        }
    
      lastCIID = ci->id;
      lastEnd = maxCoord;
    }
  return 0;
}

int DoSimpleScaffoldGraphChecks(ScaffoldGraphT * graph,
                                char * filename)
{
  GraphNodeIterator scaffolds;
  CIScaffoldT * scaff;
  FILE * fp;
  char myFilename[1024];

  sprintf(myFilename, "%s.preCkp%d", filename, graph->checkPointIteration);
  fp = fopen(myFilename, "w");
  assert(fp != NULL);
  /*
    Iterate over all real & live scaffolds
    Iterate over all contigs
    Catch all gaps < -20
    Catch 1st contig with bgn != 0
    Catch last contig with end != scaffold end
  */
  // loop over all scaffolds in the graph
  InitGraphNodeIterator(&scaffolds,
                        graph->ScaffoldGraph,
                        GRAPH_NODE_DEFAULT);
  while(NULL != (scaff = NextGraphNodeIterator(&scaffolds)))
    {
      if(scaff->flags.bits.isDead == FALSE && scaff->type == REAL_SCAFFOLD)
        {
          DoSimpleScaffoldChecks(fp, graph, scaff);
        }
    }
  fclose(fp);
  return 0;
}
  
/********************************************************************
          Functions for hashing, qsorting, et cetera
********************************************************************/


/*
  qsort comparison function
  for mate detail order (sort mate pairs left 5p to right 5p)
*/
static int md1Compare(const void * A, const void * B)
{
  const MateDetail *a = (const MateDetail *)A;
  const MateDetail *b = (const MateDetail *)B;
  return(a->fragOffset5p - b->fragOffset5p);
}


/*
  qsort comparison function
  for IIDs & sizes - sort by size (large to small)
*/
static int sizeCompare(const IID_Size * a, const IID_Size * b)
{
  return(b->length - a->length);
}


/*
  qsort comparison function
  for mate detail order (sort (fragChunkIID, mateChunkIID), fragOffset5p
*/
static int md2Compare(const MateDetail * a, const MateDetail * b)
{
  if(a->fragChunkIID < b->fragChunkIID)
    return -1;
  else if(a->fragChunkIID > b->fragChunkIID)
    return 1;
  else
    {
      if(a->mateChunkIID < b->mateChunkIID)
        return -1;
      else if(a->mateChunkIID > b->mateChunkIID)
        return 1;
      else
        {
          return(a->fragOffset5p - b->fragOffset5p);
        }
    }
}


/* qsort comparison function
   for use in first step of building cpIndex array
*/
static int cp1Compare(const CP_Index * a, const CP_Index * b)
{
  return(a->contig1 - b->contig1);
}


/* qsort comparison function
   for use in second step of building cpIndex array
*/
static int cp2Compare(const CP_Index * a, const CP_Index *b)
{
  if(a->contig1 < b->contig1)
    return -1;
  else if(a->contig1 > b->contig1)
    return 1;
  else
    return(a->contig2 - b->contig2);
}


/********************************************************************
  Functions for freeing dynamically allocated members of structures
********************************************************************/
void FreeMateStatusPositions(MateStatusPositions * msp)
{
  if(msp)
    {
      int ori1;

      for(ori1 = 0; ori1 < NUM_ORIENTATIONS_INSTR; ori1++)
        {
          int ori2;
      
          if(msp->happy[ori1])
            DeleteVA_MateDetail(msp->happy[ori1]);
          for(ori2 = 0; ori2 < NUM_ORIENTATIONS_INSTR; ori2++)
            {
              if(msp->misoriented[ori1][ori2])
                DeleteVA_MateDetail(msp->misoriented[ori1][ori2]);
            }
          if(msp->misseparatedClose[ori1])
            DeleteVA_MateDetail(msp->misseparatedClose[ori1]);
          if(msp->misseparatedFar[ori1])
            DeleteVA_MateDetail(msp->misseparatedFar[ori1]);
        }
    
      if(msp->inter)
        DeleteVA_FragDetail(msp->inter);

      safe_free(msp);
    }
}


void FreeMateStatusPositionsSet(MateStatusPositionsSet * msps)
{
  if(msps)
    {
      FreeMateStatusPositions(msps->intra);
      FreeMateStatusPositions(msps->inter);
    }
}

void DestroyMateStatusPositionsSet(MateStatusPositionsSet * msps)
{
  if(msps)
    {
      FreeMateStatusPositionsSet(msps);
      safe_free(msps);
    }
}


void FreeMateInstrumenter(MateInstrumenter * mi)
{
  if(mi)
    {
      DestroyMateStatusPositionsSet(mi->mateStatus);
      if(mi->noMate)
        DeleteVA_FragDetail(mi->noMate);
    }
}


void FreeSurrogateTracker(SurrogateTracker * st)
{
  if(st)
    {
      if(st->surrogateFragHT)
        DeleteHashTable_AS(st->surrogateFragHT);
      if(st->surrogateFragLocs)
        safe_free(st->surrogateFragLocs);
    }
}


void FreeInstrumenterBookkeeping(InstrumenterBookkeeping * bk)
{
  if(bk)
    {
      if(bk->fragHT)
        DeleteHashTable_AS(bk->fragHT);
      if(bk->fragArray)
        DeleteVA_CDS_CID_t(bk->fragArray);
      if(bk->wExtMates)
        DeleteVA_MateDetail(bk->wExtMates);
    }
}


void FreeUnitigInstrumenter(UnitigInstrumenter * ui)
{
  if(ui)
    {
      FreeMateInstrumenter(&(ui->mates));
      FreeInstrumenterBookkeeping(&(ui->bookkeeping));
      if(ui->breakpoints)
        DeleteVA_InstrumenterBreakpoint(ui->breakpoints);
    }
}


void FreeContigInstrumenter(ContigInstrumenter * ci)
{
  if(ci)
    {
      if(ci->unitigSizes)
        DeleteVA_float(ci->unitigSizes);
      if(ci->surrogateSizes)
        DeleteVA_float(ci->surrogateSizes);
      FreeUnitigInstrumenter(&(ci->reusableUI));
      FreeUnitigInstrumenter(&(ci->unitig));
      FreeMateInstrumenter(&(ci->mates));
      FreeInstrumenterBookkeeping(&(ci->bookkeeping));
      if(ci->breakpoints)
        DeleteVA_InstrumenterBreakpoint(ci->breakpoints);
    }
}


void FreeScaffoldInstrumenter(ScaffoldInstrumenter * si)
{
  if(si)
    {
      if(si->scaffoldGapSizes)
        DeleteVA_float(si->scaffoldGapSizes);
      if(si->inferredEdgeStddevs)
        DeleteVA_float(si->inferredEdgeStddevs);
      if(si->contigSizes)
        DeleteVA_float(si->contigSizes);
      FreeContigInstrumenter(&(si->reusableCI));
      FreeContigInstrumenter(&(si->contig));
      FreeMateInstrumenter(&(si->mates));
      FreeInstrumenterBookkeeping(&(si->bookkeeping));
      if(si->breakpoints)
        DeleteVA_InstrumenterBreakpoint(si->breakpoints);
      FreeSurrogateTracker(&(si->surrogateTracker));

      if(si->cpHT)
        DeleteHashTable_AS(si->cpHT);
      if(si->cpArray)
        DeleteVA_ContigPlacement(si->cpArray);

      if(si->anchoredHT)
        DeleteHashTable_AS(si->anchoredHT);
      /*
        if(si->icps)
        DeleteVA_IntContigPairs(si->icps);
        if(si->options & INST_OPT_CONTIG_PAIRS)
        {
        if(si->contigPairs)
        DeleteVA_InstrumenterContigPair(si->contigPairs);
        if(si->cpIndex)
        DeleteVA_CP_Index(si->cpIndex);
        }
      */
    }
}


void FreeScaffoldGraphInstrumenter(ScaffoldGraphInstrumenter * sgi)
{
  if(sgi)
    {
      if(sgi->singletonScaffoldSizes)
        DeleteVA_float(sgi->singletonScaffoldSizes);
      if(sgi->unitigsPerSingletonScaffold)
        DeleteVA_int32(sgi->unitigsPerSingletonScaffold);
      if(sgi->degenerateScaffoldSizes)
        DeleteVA_float(sgi->degenerateScaffoldSizes);
      if(sgi->scaffoldSizes)
        DeleteVA_float(sgi->scaffoldSizes);
      FreeScaffoldInstrumenter(&(sgi->scaffold));
      FreeInstrumenterBookkeeping(&(sgi->bookkeeping));
    }
}


/********************************************************************
  Functions for freeing dynamically allocated members of structures
        and freeing the dynamically allocated structures
********************************************************************/
void DestroyMateInstrumenter(MateInstrumenter * mi)
{
  if(mi)
    {
      FreeMateInstrumenter(mi);
      safe_free(mi);
    }
}


void DestroyUnitigInstrumenter(UnitigInstrumenter * ui)
{
  if(ui)
    {
      FreeUnitigInstrumenter(ui);
      safe_free(ui);
    }
}


void DestroyContigInstrumenter(ContigInstrumenter * ci)
{
  if(ci)
    {
      FreeContigInstrumenter(ci);
      safe_free(ci);
    }
}


void DestroyScaffoldInstrumenter(ScaffoldInstrumenter * si)
{
  if(si)
    {
      FreeScaffoldInstrumenter(si);
      safe_free(si);
    }
}


void DestroyScaffoldGraphInstrumenter(ScaffoldGraphInstrumenter * sgi)
{
  if(sgi)
    {
      FreeScaffoldGraphInstrumenter(sgi);
      safe_free(sgi);
    }
}


/********************************************************************
         Functions for allocating or reseting arrays and
                    other structure members
********************************************************************/
void InitializeFragDetailArray(VA_TYPE(FragDetail) ** fda)
{
  if(*fda == NULL)
    *fda = CreateVA_FragDetail(1000);
  else
    ResetVA_FragDetail(*fda);
}  


void InitializeMateDetailArray(VA_TYPE(MateDetail) ** mda)
{
  if(*mda == NULL)
    *mda = CreateVA_MateDetail(1000);
  else
    ResetVA_MateDetail(*mda);
}  


void InitializeMateStatusPositions(MateStatusPositions * msp)
{
  int ori1;
  
  for(ori1 = 0; ori1 < NUM_ORIENTATIONS_INSTR; ori1++)
    {
      int ori2;
    
      InitializeMateDetailArray(&(msp->happy[ori1]));
      for(ori2 = 0; ori2 < NUM_ORIENTATIONS_INSTR; ori2++)
        InitializeMateDetailArray(&(msp->misoriented[ori1][ori2]));
      InitializeMateDetailArray(&(msp->misseparatedClose[ori1]));
      InitializeMateDetailArray(&(msp->misseparatedFar[ori1]));
    }
  
  InitializeFragDetailArray(&(msp->inter));
}


void InitializeMateStatusPositionsSet(MateStatusPositionsSet * msps)
{
  InitializeMateStatusPositions(msps->intra);
  InitializeMateStatusPositions(msps->inter);
}


void ResetMateStats(MateStats * ms)
{
  ms->reads = 0;
  ms->externalReads = 0;
  ms->externalFrags = 0;
}


void ResetMateStatsSet(MateStatsSet * mss)
{
  ResetMateStats(&(mss->happy));
  ResetMateStats(&(mss->misoriented));
  ResetMateStats(&(mss->misseparatedClose));
  ResetMateStats(&(mss->misseparatedFar));
  ResetMateStats(&(mss->inter));
}


void ResetMateInstrumenterCounts(MateInstrumenter * mi)
{
  ResetMateStatsSet(&(mi->intra));
  ResetMateStatsSet(&(mi->inter));
  ResetMateStats(&(mi->mateless));
}


MateStatusPositions * CreateMateStatusPositions(void)
{
  MateStatusPositions * msp;
  
  msp = safe_calloc(1, sizeof(MateStatusPositions));

  InitializeMateStatusPositions(msp);
  return msp;
}


MateStatusPositionsSet * CreateMateStatusPositionsSet(void)
{
  MateStatusPositionsSet * msps;

  msps = safe_calloc(1, sizeof(MateStatusPositionsSet));

  if((msps->intra = CreateMateStatusPositions()) == NULL)
    {
      fprintf(stderr, "Failed to allocate mate status positions arrays.\n");
      return NULL;
    }
  if((msps->inter = CreateMateStatusPositions()) == NULL)
    {
      fprintf(stderr, "Failed to allocate mate status positions arrays.\n");
      return NULL;
    }
  return msps;
}


int InitializeMateInstrumenter(ScaffoldGraphT * graph,
                               MateInstrumenter * mi)
{
  if(mi->mateStatus == NULL)
    {
      if((mi->mateStatus = CreateMateStatusPositionsSet()) == NULL)
        {
          fprintf(stderr, "Failed to allocate mate status positions\n");
          return 1;
        }
    }
  else
    {
      InitializeMateStatusPositionsSet(mi->mateStatus);
    }

  ResetMateInstrumenterCounts(mi);

  if(mi->noMate == NULL)
    {
      mi->noMate = CreateVA_FragDetail(1000);
      if(mi->noMate == NULL)
        {
          fprintf(stderr, "Failed to create no mate variable array\n");
          return 1;
        }
    }
  else
    {
      ResetVA_FragDetail(mi->noMate);
    }
  return 0;
}


MateInstrumenter * CreateMateInstrumenter(ScaffoldGraphT * graph,
                                          uint32 options)
{
  MateInstrumenter * mi = safe_calloc(1, sizeof(MateInstrumenter));
  
  mi->options = options;
  
  if(InitializeMateInstrumenter(graph, mi))
    {
      fprintf(stderr, "Failed to initialize MateInstrumenter!\n");
      safe_free(mi);
      return NULL;
    }
  return mi;
}


int InitializeSurrogateTracker(ScaffoldGraphT * graph,
                               SurrogateTracker * st)
{
  if(st->surrogateFragHT == NULL)
    {
      // assume 1/10 # of frags?
      st->numAllocatedLocs = MAX(50, GetNumCIFragTs(graph->CIFrags) / 100);
      st->surrogateFragHT = CreateScalarHashTable_AS(st->numAllocatedLocs);
      if(st->surrogateFragHT == NULL)
        {
          fprintf(stderr, "Failed to allocate surrogate fragment hashtable\n");
          return 1;
        }
    }
  else
    {
      ResetHashTable_AS(st->surrogateFragHT);
    }
  
  if(st->surrogateFragLocs == NULL)
    {
      st->numAllocatedLocs = MAX(50, GetNumCIFragTs(graph->CIFrags) / 100);
      st->numUsedLocs = 0;
      st->surrogateFragLocs = safe_calloc(st->numAllocatedLocs,
                                          sizeof(SurrogateFragLocation));
    }
  else
    {
      memset(st->surrogateFragLocs,
             0,
             st->numAllocatedLocs * sizeof(SurrogateFragLocation));
      st->numUsedLocs = 0;
    }
  return 0;
}


int InitializeInstrumenterBookkeeping(ScaffoldGraphT * graph,
                                      InstrumenterBookkeeping * bk,
                                      InstrumenterLevel level)
{
  int32 numFrags = 0;
  int32 numWithExternalMates = 0;

  if(bk->fragHT == NULL ||
     bk->fragArray == NULL ||
     bk->wExtMates == NULL)
    {

#define LARGEST_TO_MEAN_RATIO   5.f
    
      // guestimates for max frags per contig & scaffold
      // largest contg/scaffold is x times larger than mean?
      switch(level)
        {
          case InstrumenterUnitigLevel:
          case InstrumenterContigLevel:
            if(graph->numContigs == 0)
              {
                fprintf(stderr,
                        "*** Inititializing contig instrumenter bookkeeping, "
                        "but graph has no contigs! ***\n");
                // this shouldn't be the case...
                numFrags = 128;
              }
            else
              {
                numFrags = LARGEST_TO_MEAN_RATIO *
                  (GetNumCIFragTs(graph->CIFrags) / MAX(1, graph->numContigs));
              }
            break;
          case InstrumenterScaffoldGraphLevel:
          case InstrumenterScaffoldLevel:
            if(graph->numLiveScaffolds == 0)
              {
                fprintf(stderr,
                        "*** Inititializing scaffold instrumenter bookkeeping, "
                        "but graph has no live scaffolds! ***\n");
                numFrags = LARGEST_TO_MEAN_RATIO *
                  (GetNumCIFragTs(graph->CIFrags) / MAX(1, graph->numContigs));
              }
            else
              {
                numFrags = LARGEST_TO_MEAN_RATIO *
                  (GetNumCIFragTs(graph->CIFrags) / graph->numLiveScaffolds);
              }
            break;
            /*
              numFrags = GetNumCIFragTs(graph->CIFrags);
              break;
            */
        }
      numFrags = MAX(3, MIN(numFrags, GetNumCIFragTs(graph->CIFrags)));
      numWithExternalMates = numFrags / 3;
    }

  if(bk->fragHT == NULL)
    {
      bk->fragHT = CreateScalarHashTable_AS(numFrags);
      if(bk->fragHT == NULL)
        {
          fprintf(stderr, "Failed to allocate fragment hashtable\n");
          return 1;
        }
    }
  else
    {
      ResetHashTable_AS(bk->fragHT);
    }

  if(bk->fragArray == NULL)
    {
      bk->fragArray = CreateVA_CDS_CID_t(numFrags);
      if(bk->fragArray == NULL)
        {
          fprintf(stderr, "Failed to allocate fragment array\n");
          return 1;
        }
    }
  else
    {
      ResetVA_CDS_CID_t(bk->fragArray);
    }
  
  if(bk->wExtMates == NULL)
    {
      bk->wExtMates = CreateVA_MateDetail(numWithExternalMates);
    }
  else
    {
      ResetVA_MateDetail(bk->wExtMates);
    }

  return 0;
}


int InitializeUnitigInstrumenter(ScaffoldGraphT * graph,
                                 UnitigInstrumenter * ui)
{
  ui->id = NULLINDEX;
  ui->isSurrogate = 0;
  ui->leftEnd = ui->rightEnd = ui->orientation = 0;

  ui->numReads = ui->numExtReads = ui->numExtFrags = 0;

  ui->mates.options = ui->options;
  if(InitializeMateInstrumenter(graph, &(ui->mates)))
    {
      fprintf(stderr, "Failed to initialize mate instrumenter\n");
      return 1;
    }

  if(InitializeInstrumenterBookkeeping(graph,
                                       &(ui->bookkeeping),
                                       InstrumenterUnitigLevel))
    {
      fprintf(stderr, "Failed to initialize unitig bookkeeping data\n");
      return 1;
    }

  // allocate or reset breakpoints
  if(ui->options & INST_OPT_BREAKPOINTS)
    {
      if(ui->breakpoints == NULL)
        {
          ui->breakpoints = CreateVA_InstrumenterBreakpoint(100);
        }
      else
        {
          ResetVA_InstrumenterBreakpoint(ui->breakpoints);
        }
    }
  else
    ui->breakpoints = NULL;

  return 0;
}


int InitializeContigInstrumenter(ScaffoldGraphT * graph,
                                 ContigInstrumenter * ci)
{
  ci->id = NULLINDEX;
  ci->leftEnd = ci->rightEnd = ci->orientation = 0;

  // allocate or reset unitig lengths
  if(ci->unitigSizes == NULL)
    {
      ci->unitigSizes = CreateVA_float(1000);
    }
  else
    {
      ResetVA_float(ci->unitigSizes);
    }

  // allocate of reset surrogate lengths
  if(ci->surrogateSizes == NULL)
    {
      ci->surrogateSizes = CreateVA_float(100);
    }
  else
    {
      ResetVA_float(ci->surrogateSizes);
    }
  
  ci->numReads = ci->numExtReads = ci->numExtFrags = 0;

  // initialize unitig instrumenters
  ci->reusableUI.options = ci->options;
  if(InitializeUnitigInstrumenter(graph, &(ci->reusableUI)))
    {
      fprintf(stderr, "Failed to initialize unitig instrumenter\n");
      return 1;
    }
  ci->unitig.options = ci->options;
  if(InitializeUnitigInstrumenter(graph, &(ci->unitig)))
    {
      fprintf(stderr, "Failed to initialize unitig instrumenter\n");
      return 1;
    }
  
  ci->mates.options = ci->options;
  if(InitializeMateInstrumenter(graph, &(ci->mates)))
    {
      fprintf(stderr, "Failed to initialize mate instrumenter\n");
      return 1;
    }

  if(InitializeInstrumenterBookkeeping(graph,
                                       &(ci->bookkeeping),
                                       InstrumenterContigLevel))
    {
      fprintf(stderr, "Failed to initialize contig bookkeeping data\n");
      return 1;
    }

  // allocate or reset breakpoints
  if(ci->options & INST_OPT_BREAKPOINTS)
    {
      if(ci->breakpoints == NULL)
        {
          ci->breakpoints = CreateVA_InstrumenterBreakpoint(100);
        }
      else
        {
          ResetVA_InstrumenterBreakpoint(ci->breakpoints);
        }
    }
  else
    ci->breakpoints = NULL;

  return 0;
}


int InitializeScaffoldInstrumenter(ScaffoldGraphT * graph,
                                   ScaffoldInstrumenter * si)
{
  si->id = NULLINDEX;

  si->size = 0.0f;
  
  if(si->scaffoldGapSizes == NULL)
    {
      si->scaffoldGapSizes = CreateVA_float(100);
    }
  else
    {
      ResetVA_float(si->scaffoldGapSizes);
    }

  if(si->inferredEdgeStddevs == NULL)
    {
      si->inferredEdgeStddevs = CreateVA_float(100);
    }
  else
    {
      ResetVA_float(si->inferredEdgeStddevs);
    }

  if(si->contigSizes == NULL)
    {
      si->contigSizes = CreateVA_float(100);
    }
  else
    {
      ResetVA_float(si->contigSizes);
    }

  /*
  // Initialize the aggregate unitig instrumenter
  si->unitig.options = si->options;
  if(si->unitig.options & INST_OPT_INTER_MATES)
  si->unitig.options ^= INST_OPT_INTER_MATES;
  if(InitializeUnitigInstrumenter(graph, &(si->unitig)))
  {
  fprintf(stderr, "Failed to initialize aggregate unitig instrumenter\n");
  return 1;
  }
  */

  // Initialize the reusable contig instrumenter
  si->reusableCI.options = si->options;
  if(si->reusableCI.options & INST_OPT_INTER_MATES)
    si->reusableCI.options ^= INST_OPT_INTER_MATES;
  if(InitializeContigInstrumenter(graph, &(si->reusableCI)))
    {
      fprintf(stderr, "Failed to initialize reusable contig instrumenter\n");
      return 1;
    }

  // Initialize the aggregate contig instrumenter
  si->contig.options = si->options;
  if(si->contig.options & INST_OPT_INTER_MATES)
    si->contig.options ^= INST_OPT_INTER_MATES;
  if(InitializeContigInstrumenter(graph, &(si->contig)))
    {
      fprintf(stderr, "Failed to initialize aggregate contig instrumenter\n");
      return 1;
    }

  si->mates.options = si->options;
  if(InitializeMateInstrumenter(graph, &(si->mates)))
    {
      fprintf(stderr, "Failed to initialize mate status positions arrays\n");
      return 1;
    }

  if(InitializeInstrumenterBookkeeping(graph,
                                       &(si->bookkeeping),
                                       InstrumenterScaffoldLevel))
    {
      fprintf(stderr, "Failed to initialize scaffold bookkeeping data\n");
      return 1;
    }

  // allocate or reset breakpoints
  if(si->options & INST_OPT_BREAKPOINTS)
    {
      if(si->breakpoints == NULL)
        {
          si->breakpoints = CreateVA_InstrumenterBreakpoint(100);
        }
      else
        {
          ResetVA_InstrumenterBreakpoint(si->breakpoints);
        }
    }
  else
    si->breakpoints = NULL;

  if(InitializeSurrogateTracker(graph, &(si->surrogateTracker)))
    {
      fprintf(stderr, "Failed to initialize surrogate tracker\n");
      return 1;
    }

  if(si->cpHT == NULL)
    {
      if((si->cpHT = CreateScalarHashTable_AS(1000)) == NULL)
        {
          fprintf(stderr, "Failed to allocate contig pair hashtable\n");
          return 1;
        }
    }
  else
    {
      ResetHashTable_AS(si->cpHT);
    }

  // NOTE: If this gets resized, cpHT needs to be repopulated

  if (si->cpArray == NULL) {
    si->cpArray = CreateVA_ContigPlacement(1000);
  }
  ResetVA_ContigPlacement(si->cpArray);

  if(si->anchoredHT == NULL)
    {
      if((si->anchoredHT = CreateScalarHashTable_AS(1000)) == NULL)
        {
          fprintf(stderr, "Failed to allocate contig anchoring hashtable\n");
          return 1;
        }
    }
  else
    {
      ResetHashTable_AS(si->anchoredHT);
    }
  
  return 0;
}


/********************************************************************
              Functions for creating instrumenters
********************************************************************/
UnitigInstrumenter * CreateUnitigInstrumenter(ScaffoldGraphT * graph,
                                              uint32 options)
{
  UnitigInstrumenter * ui;
  ui = (UnitigInstrumenter *) safe_calloc(1, sizeof(UnitigInstrumenter));

  ui->options = options;
  if(ui->options & INST_OPT_INTER_MATES)
    ui->options ^= INST_OPT_INTER_MATES;
  
  if(InitializeUnitigInstrumenter(graph, ui))
    {
      fprintf(stderr, "Failed to initialize UnitigInstrumenter!\n");
      safe_free(ui);
      return NULL;
    }
  return ui;
}


ContigInstrumenter * CreateContigInstrumenter(ScaffoldGraphT * graph,
                                              uint32 options)
{
  ContigInstrumenter * ci;
  ci = (ContigInstrumenter *) safe_calloc(1, sizeof(ContigInstrumenter));

  ci->options = options;
  if(ci->options & INST_OPT_INTER_MATES)
    ci->options ^= INST_OPT_INTER_MATES;
  
  if(InitializeContigInstrumenter(graph, ci))
    {
      fprintf(stderr, "Failed to initialize ContigInstrumenter!\n");
      safe_free(ci);
      return NULL;
    }
  return ci;
}


ScaffoldInstrumenter * CreateScaffoldInstrumenter(ScaffoldGraphT * graph,
                                                  uint32 options)
{
  ScaffoldInstrumenter * si;
  si = (ScaffoldInstrumenter *) safe_calloc(1, sizeof(ScaffoldInstrumenter));

  si->options = options;
  if(InitializeScaffoldInstrumenter(graph, si))
    {
      fprintf(stderr, "Failed to initialize ScaffoldInstrumenter!\n");
      safe_free(si);
      return NULL;
    }
  return si;
}


int InitializeScaffoldGraphInstrumenter(ScaffoldGraphT * graph,
                                        ScaffoldGraphInstrumenter * sgi)
{
  if(sgi->singletonScaffoldSizes == NULL)
    {
      sgi->singletonScaffoldSizes = CreateVA_float(10000);
    }
  else
    {
      ResetVA_float(sgi->singletonScaffoldSizes);
    }

  if(sgi->unitigsPerSingletonScaffold == NULL)
    {
      sgi->unitigsPerSingletonScaffold = CreateVA_int32(10000);
    }
  else
    {
      ResetVA_int32(sgi->unitigsPerSingletonScaffold);
    }

  if(sgi->degenerateScaffoldSizes == NULL)
    {
      sgi->degenerateScaffoldSizes = CreateVA_float(10000);
    }
  else
    {
      ResetVA_float(sgi->degenerateScaffoldSizes);
    }

  sgi->numDegenerateScaffoldsWithoutReads = 0;

  if(sgi->scaffoldSizes == NULL)
    {
      sgi->scaffoldSizes = CreateVA_float(10000);
    }
  else
    {
      ResetVA_float(sgi->scaffoldSizes);
    }

  sgi->scaffold.options = sgi->options;
  if(InitializeScaffoldInstrumenter(graph, &(sgi->scaffold)))
    {
      fprintf(stderr, "Failed to initialize scaffold instrumenter\n");
      return 1;
    }

  if(InitializeInstrumenterBookkeeping(graph,
                                       &(sgi->bookkeeping),
                                       InstrumenterScaffoldGraphLevel))
    {
      fprintf(stderr, "Failed to initialize scaffold graph bookkeeping data\n");
      return 1;
    }
  return 0;
}


ScaffoldGraphInstrumenter *
CreateScaffoldGraphInstrumenter(ScaffoldGraphT * graph, uint32 options)
{
  ScaffoldGraphInstrumenter * sgi;
  sgi =
    (ScaffoldGraphInstrumenter *)safe_calloc(1, sizeof(ScaffoldGraphInstrumenter));
  
  sgi->options = options;
  if(InitializeScaffoldGraphInstrumenter(graph, sgi))
    {
      fprintf(stderr, "Failed to initialize ScaffoldGraphInstrumenter!\n");
      safe_free(sgi);
      return NULL;
    }
  return sgi;
}


void FindRockStoneUnitigs(ScaffoldGraphT * graph)
{
  GraphNodeIterator unitigIterator;
  ChunkInstanceT * unitig;
  int32 numRockStones = 0;

  InitGraphNodeIterator(&unitigIterator, graph->CIGraph, GRAPH_NODE_DEFAULT);
  while((unitig = NextGraphNodeIterator(&unitigIterator)) != NULL)
    {
      if(unitig->flags.bits.isStone && unitig->flags.bits.isRock)
        fprintf(stderr, "%d. unitig " F_CID " is both a rock and a stone.\n",
                ++numRockStones, unitig->id);
    }
}

/********************************************************************
                 Functions for detecting breakpoints
********************************************************************/
int AppendBreakpointSet(VA_TYPE(InstrumenterBreakpoint) * bps,
                        InstrumenterBreakpoint * bp1,
                        InstrumenterBreakpoint * bp2,
                        InstrumenterBreakpoint * bp3)
{
  static CDS_CID_t persistentIID = 0;

  bp1->iid = persistentIID++;
  
  // bp1 is left, bp2 is right, bp3 is in between or NULL
  // see if two breakpoint intervals overlap
  if(bp1->end2 >= bp2->end1)
    {
      bp1->end2 = bp2->end2;
      bp1->section = BP_ALL;
      AppendVA_InstrumenterBreakpoint(bps, bp1);
    }
  else
    {
      bp1->section = BP_LEFT;
      bp2->iid = bp1->iid;
      bp2->section = BP_RIGHT;
      bp2->type = bp1->type;
      bp2->pairs = bp1->pairs;
      bp2->contig1 = bp1->contig1;
      bp2->contig2 = bp1->contig2;
      AppendVA_InstrumenterBreakpoint(bps, bp1);
      AppendVA_InstrumenterBreakpoint(bps, bp2);
      if(bp3)
        {
          bp3->iid = bp1->iid;
          bp3->section = BP_MIDDLE;
          bp3->type = bp1->type;
          bp3->pairs = bp1->pairs;
          bp3->contig1 = bp1->contig1;
          bp3->contig2 = bp1->contig2;
          AppendVA_InstrumenterBreakpoint(bps, bp3);
        }
    }
  return 0;
}


int CreateBreakpointIntervalsFromMateDetail(MateDetail * md,
                                            DistT * dptr,
                                            InstBreakpointType problem,
                                            InstrumenterBreakpoint * bp1,
                                            InstrumenterBreakpoint * bp2,
                                            InstrumenterBreakpoint * bp3)
{
  switch( problem )
    {
      case INST_BP_TOO_CLOSE:
        /*
          --->     <---
          (   bp    )
        */
        bp1->end1 = bp2->end1 = md->fragOffset5p + READ_TRIM_BASES;
        bp1->end2 = bp2->end2 = md->mateOffset5p - READ_TRIM_BASES;
        bp1->contigEnd = bp1->end2;
        break;
      case INST_BP_TOO_FAR:
        /*
          --->                    <---
          (    bp1    ) (    bp2   )
          bp1 & bp2 may overlap
        */
        bp1->end1 = md->fragOffset5p + READ_TRIM_BASES;
        bp1->end2 = md->fragOffset5p + READ_TRIM_BASES + INTERVAL(dptr);
        bp2->end1 = md->mateOffset5p - READ_TRIM_BASES - INTERVAL(dptr);
        bp2->end2 = md->mateOffset5p - READ_TRIM_BASES;
        bp1->contigEnd = bp2->end2;
        break;
      case INST_BP_NORMAL:
        /*
          --->            --->
          (   bp1   )     (   bp2   )
          bp1 & bp2 may overlap
        */
        bp1->end1 = md->fragOffset5p + READ_TRIM_BASES;
        bp1->end2 = md->fragOffset5p + READ_TRIM_BASES + INTERVAL(dptr);
        bp2->end1 = md->mateOffset5p + READ_TRIM_BASES;
        bp2->end2 = md->mateOffset5p + READ_TRIM_BASES + INTERVAL(dptr);
        bp1->contigEnd = bp2->end1;
        break;
      case INST_BP_ANTINORMAL:
        /*
          <---            <---
          (   bp1   )     (   bp2   )
          bp1 & bp2 may overlap
        */
        bp1->end1 = md->fragOffset5p - READ_TRIM_BASES - INTERVAL(dptr);
        bp1->end2 = md->fragOffset5p - READ_TRIM_BASES;
        bp2->end1 = md->mateOffset5p - READ_TRIM_BASES - INTERVAL(dptr);
        bp2->end2 = md->mateOffset5p - READ_TRIM_BASES;
        bp1->contigEnd = bp2->end2;
        break;
      case INST_BP_OUTTIE:
        /*
          <---     --->
          (   bp1   )       (   bp2   )
          (  bp3  )
          bp1 & bp2 may overlap over the set of pairs
        */
        bp1->end1 = md->fragOffset5p - READ_TRIM_BASES - INTERVAL(dptr);
        bp1->end2 = md->fragOffset5p - READ_TRIM_BASES;
        bp2->end1 = md->mateOffset5p + READ_TRIM_BASES;
        bp2->end2 = md->mateOffset5p + READ_TRIM_BASES + INTERVAL(dptr);
        bp3->end1 = bp1->end2;
        bp3->end2 = bp2->end1;
        bp1->contigEnd = bp2->end1;
        break;
      case INST_BP_INNIE:
        /*
          --->   <---
          (   bp1   )       (   bp2   )
          (  bp3  )
          bp1 & bp2 may overlap over the set of pairs
        */
        bp1->end1 = md->fragOffset5p + READ_TRIM_BASES - INTERVAL(dptr);
        bp1->end2 = md->fragOffset5p + READ_TRIM_BASES;
        bp2->end1 = md->mateOffset5p - READ_TRIM_BASES;
        bp2->end2 = md->mateOffset5p - READ_TRIM_BASES + INTERVAL(dptr);
        bp3->end1 = bp1->end2;
        bp3->end2 = bp2->end1;
        bp1->contigEnd = bp2->end1;
        break;
      default:
        return 1;
        break;
    }
  // track some details using only breakpoint interval 1
  bp1->contig1 = md->fragChunkIID;
  bp1->contig2 = md->mateChunkIID;
  bp1->type = problem;
  bp1->mates = NULL;
  bp1->pairs = 1;
  return 0;
}


int BreakpointsOverlap(InstBreakpointType problem,
                       InstrumenterBreakpoint * currBP1,
                       InstrumenterBreakpoint * currBP2,
                       InstrumenterBreakpoint * newBP1,
                       InstrumenterBreakpoint * newBP2)
{
  switch(problem)
    {
      case INST_BP_TOO_CLOSE:
        /*
          --->     <---
          (   bp    )
                   
          Not in same breakpoint if
          newBP1->end1 is to the right of bp1->end2
        */
        return((newBP1->end1 > currBP1->end2) ? 0 : 1);
        break;
      case INST_BP_TOO_FAR:
        /*
          --->                    <---
          (    bp1    ) (    bp2   )
          bp1 & bp2 may overlap
          When proceeding left to right, mate pair is not in breakpoint when:
          A: left fragment's interval doesn't intersect current bp1
          i.e., md->fragOffset5p >= bp1->end2
          B: right fragment's interval doesn't intersect current bp2
          i.e., md->mateOffset5p - (u + ns) >= bp2->end2
          or, md->mateOffset5p < bp2->end1
          NOTE: This suggests that bp1 & bp2 should be pursued as separate
          breakpoints
        */
      case INST_BP_NORMAL:
        /*
          --->            --->
          (   bp1   )     (   bp2   )
          bp1 & bp2 may overlap
        */
      case INST_BP_ANTINORMAL:
        /*
          <---            <---
          (   bp1   )     (   bp2   )
          bp1 & bp2 may overlap
        */
      case INST_BP_OUTTIE:
        /*
          <---     --->
          (   bp1   )       (   bp2   )
          (  bp3  )
          bp1 & bp2 may overlap over the set of pairs
        */
      case INST_BP_INNIE:
        /*
          --->   <---
          (   bp1   )       (   bp2   )
          (  bp3  )
          bp1 & bp2 may overlap over the set of pairs
        */
        /*
          For all of the cases other than too close,
          Not in same breakpoint if
          newBP1->end1 is to the right of bp1->end2, or
          newBP2->end1 is to the right of bp2->end2, or
          newBP2->end2 is to the left of bp2->end1
        */
        return(((newBP1->end1 > currBP1->end2) ? 0 :
                ((newBP2->end1 > currBP2->end2) ? 0 :
                 ((newBP2->end2 < currBP2->end1) ? 0 : 1))));
        break;
      default:
        break;
    }
  return 1;
}


void NarrowBreakpointInterval(InstrumenterBreakpoint * destBP,
                              InstrumenterBreakpoint * sourceBP)
{
  destBP->end1 = MAX(destBP->end1, sourceBP->end1);
  destBP->end2 = MIN(destBP->end2, sourceBP->end2);
}


void UpdateBreakpointSet(InstBreakpointType problem,
                         InstrumenterBreakpoint * inOutBP1,
                         InstrumenterBreakpoint * inOutBP2,
                         InstrumenterBreakpoint * inOutBP3,
                         InstrumenterBreakpoint * inBP1,
                         InstrumenterBreakpoint * inBP2,
                         InstrumenterBreakpoint * inBP3)
{
  switch(problem)
    {
      case INST_BP_TOO_CLOSE:
        /*
          --->     <---
          (   bp    )
          Want to narrow bp
        */
        NarrowBreakpointInterval(inOutBP1, inBP1);
        break;
      case INST_BP_TOO_FAR:
        /*
          --->                    <---
          (    bp1    ) (    bp2   )
          bp1 & bp2 may overlap
        */
      case INST_BP_NORMAL:
        /*
          --->            --->
          (   bp1   )     (   bp2   )
          bp1 & bp2 may overlap
        */
      case INST_BP_ANTINORMAL:
        /*
          <---            <---
          (   bp1   )     (   bp2   )
          bp1 & bp2 may overlap
        */
        NarrowBreakpointInterval(inOutBP1, inBP1);
        NarrowBreakpointInterval(inOutBP2, inBP2);
        break;
      case INST_BP_OUTTIE:
        /*
          <---     --->
          (   bp1   )       (   bp2   )
          (  bp3  )
          bp1 & bp2 may overlap over the set of pairs
        */
      case INST_BP_INNIE:
        /*
          --->   <---
          (   bp1   )       (   bp2   )
          (  bp3  )
          bp1 & bp2 may overlap over the set of pairs
        */
        NarrowBreakpointInterval(inOutBP1, inBP1);
        NarrowBreakpointInterval(inOutBP2, inBP2);
        NarrowBreakpointInterval(inOutBP3, inBP3);
        break;
      default:
        assert(0);
        break;
    }

  // keep track of the right-most position & contig involved
  if(inBP1->contigEnd > inOutBP1->contigEnd)
    {
      inOutBP1->contig2 = inBP1->contig2;
      inOutBP1->contigEnd = inBP1->contigEnd;
    }
  inOutBP1->pairs++;
  return;
}


/*
  Function to detect a given type of breakpoint given a set of
  that type of bad mate intervals. Linear in the number of input
  intervals

  MateDetails in mda keep only 5p end of each fragment
  frag is left-most, mate is right-most, and problem specifies
  orientation of frag & mate

  currently assuming mates should be innie unless problem is INST_BP_INNIE,
  in which case mates should be outtie
*/
int DetectBreakpointType(ScaffoldGraphT * graph,
                         VA_TYPE(InstrumenterBreakpoint) * bps,
                         CDS_COORD_t bgn,
                         CDS_COORD_t end,
                         VA_TYPE(MateDetail) * mda,
                         InstBreakpointType problem)
{
  // don't bother if there aren't enough intervals
  if(GetNumVA_MateDetail(mda) >= INST_MIN_BREAK_MATES)
    {
      InstrumenterBreakpoint bp1;
      InstrumenterBreakpoint bp2;
      InstrumenterBreakpoint bp3;
      int32 i;
      MateDetail * md;
      DistT * dptr;

      // sort intervals by frag's (left-most) offset5p
      qsort(GetVA_MateDetail(mda, 0),
            GetNumVA_MateDetail(mda),
            sizeof(MateDetail),
            md1Compare );

      // loop over intervals
      for(i = 0, bp1.pairs = 0; i < GetNumVA_MateDetail(mda); i++)
        {
          InstrumenterBreakpoint newBP1;
          InstrumenterBreakpoint newBP2;
          InstrumenterBreakpoint newBP3;
          md = GetVA_MateDetail(mda, i);
          dptr = GetDistT(graph->Dists, md->libIID);

          CreateBreakpointIntervalsFromMateDetail(md, dptr, problem,
                                                  &newBP1, &newBP2, &newBP3);
      
          // if a new possible breakpoint should be initiated
          if(i == 0 || ! BreakpointsOverlap(problem, &bp1, &bp2, &newBP1, &newBP2))
            {
#ifdef DEBUG_BRK
              if(problem == INST_BP_OUTTIE)
                {
                  PrintMateDetailAndDist(md, dptr, "", stderr);
                  fprintf(stderr, "bp1. ");
                  PrintBreakpoint(&newBP1, "", stderr);
                  fprintf(stderr, "bp2. ");
                  PrintBreakpoint(&newBP2, "", stderr);
                  fprintf(stderr, "New/End-of-Previous breakpoint\n");
                }
#endif
              // if there was a good one, save it
              if(bp1.pairs >= INST_MIN_BREAK_MATES)
                {
#ifdef DEBUG_BRK
                  fprintf(stderr, "Got a keeper!\n");
#endif
                  /*
                    bp1.end1 = MAX(bp1.end1, bgn);
                    bp1.end2 = MIN(bp1.end2, end);
                    bp2.end1 = MAX(bp2.end1, bgn);
                    bp2.end2 = MIN(bp2.end2, end);
                  */
                  AppendBreakpointSet(bps, &bp1, &bp2,
                                      ((problem == INST_BP_INNIE ||
                                        problem == INST_BP_OUTTIE) ?
                                       &bp3 : NULL));
                }
              // start considering new breakpoing interval
              bp1 = newBP1;
              bp2 = newBP2;
              bp3 = newBP3;
            }
          else
            {
              // if here, continuing in same breakpoint
              UpdateBreakpointSet(problem,
                                  &bp1, &bp2, &bp3,
                                  &newBP1, &newBP2, &newBP3);
            }
        }
      // if there's an unfinished breakpoint interval, save it
      if(bp1.pairs >= INST_MIN_BREAK_MATES)
        {
          /*
            bp1.end1 = MAX(bp1.end1, bgn);
            bp1.end2 = MIN(bp1.end2, end);
            bp2.end1 = MAX(bp2.end1, bgn);
            bp2.end2 = MIN(bp2.end2, end);
          */
          AppendBreakpointSet(bps, &bp1, &bp2,
                              ((problem == INST_BP_INNIE ||
                                problem == INST_BP_OUTTIE) ?
                               &bp3 : NULL));
        }
    }
  return 0;
}


int DetectRoughIntraContigBreakpoints(ScaffoldGraphT * graph,
                                      ContigInstrumenter * ci)
{
  // should be innie, are too close
  DetectBreakpointType(graph, ci->breakpoints,
                       0, abs(ci->rightEnd - ci->leftEnd),
                       ci->mates.mateStatus->intra->misseparatedClose[INNIE_INSTR],
                       INST_BP_TOO_CLOSE);

  // should be innie, are too far apart
  DetectBreakpointType(graph, ci->breakpoints,
                       0, abs(ci->rightEnd - ci->leftEnd),
                       ci->mates.mateStatus->intra->misseparatedFar[INNIE_INSTR],
                       INST_BP_TOO_FAR);

  // should be innie, are normal
  DetectBreakpointType(graph, ci->breakpoints,
                       0, abs(ci->rightEnd - ci->leftEnd),
                       ci->mates.mateStatus->intra->misoriented[INNIE_INSTR][NORMAL_INSTR],
                       INST_BP_NORMAL);

  // should be innie, are antinormal
  DetectBreakpointType(graph, ci->breakpoints,
                       0, abs(ci->rightEnd - ci->leftEnd),
                       ci->mates.mateStatus->intra->misoriented[INNIE_INSTR][ANTINORMAL_INSTR],
                       INST_BP_ANTINORMAL);

  // should be innie, are outtie
  DetectBreakpointType(graph, ci->breakpoints,
                       0, abs(ci->rightEnd - ci->leftEnd),
                       ci->mates.mateStatus->intra->misoriented[INNIE_INSTR][OUTTIE_INSTR],
                       INST_BP_OUTTIE);
  return 0;
}


int DetectRoughInterContigBreakpoints(ScaffoldGraphT * graph,
                                      ScaffoldInstrumenter * si)
{
  // should be innie, are too close
  DetectBreakpointType(graph, si->breakpoints,
                       0, (CDS_COORD_t) si->size,
                       si->mates.mateStatus->inter->misseparatedClose[INNIE_INSTR],
                       INST_BP_TOO_CLOSE);

  // should be innie, are too far apart
  DetectBreakpointType(graph, si->breakpoints,
                       0, (CDS_COORD_t) si->size,
                       si->mates.mateStatus->inter->misseparatedFar[INNIE_INSTR],
                       INST_BP_TOO_FAR);

  // should be innie, are normal
  DetectBreakpointType(graph, si->breakpoints,
                       0, (CDS_COORD_t) si->size,
                       si->mates.mateStatus->inter->misoriented[INNIE_INSTR][NORMAL_INSTR],
                       INST_BP_NORMAL);

  // should be innie, are antinormal
  DetectBreakpointType(graph, si->breakpoints,
                       0, (CDS_COORD_t) si->size,
                       si->mates.mateStatus->inter->misoriented[INNIE_INSTR][ANTINORMAL_INSTR],
                       INST_BP_ANTINORMAL);

  // should be innie, are outtie
  DetectBreakpointType(graph, si->breakpoints,
                       0, (CDS_COORD_t) si->size,
                       si->mates.mateStatus->inter->misoriented[INNIE_INSTR][OUTTIE_INSTR],
                       INST_BP_OUTTIE);
  return 0;
}


/********************************************************************
                  Functions for computing stats
********************************************************************/
void AddFloatToInstrumenterStatistics(InstrumenterStatistics * is,
                                      float * var)
{
  is->mean += (*var);
  is->sumOfSquares += (*var) * (*var);
  is->min = MIN(is->min, *var);
  is->max = MAX(is->max, *var);
}


void AddCoordToInstrumenterStatistics(InstrumenterStatistics * is,
                                      CDS_COORD_t * var)
{
  is->mean += (*var);
  is->sumOfSquares += (*var) * (*var);
  is->min = MIN(is->min, *var);
  is->max = MAX(is->max, *var);
}


void InstrumentFloatStatistics(VA_TYPE(float) * va,
                               InstrumenterStatistics * is,
                               int separateNegatives)
{
  int32 i;
  float * var;

  memset(is, 0, sizeof(InstrumenterStatistics));
  is->min = FLT_MAX;
  is->minNegative = FLT_MAX;

  // use mean to hold sum
  for(i = 0; i < GetNumVA_float(va); i++)
    {
      var = GetVA_float(va, i);

      if((*var) < 0.0f)
        {
          is->numNegatives++;
          is->minNegative = MIN(is->minNegative, *var);
          if(!separateNegatives)
            AddFloatToInstrumenterStatistics(is, var);
        }
      else
        {
          is->numPositives++;
          AddFloatToInstrumenterStatistics(is, var);
        }
    }

  is->min = (is->min == FLT_MAX) ? 0.0f : is->min;
  if(separateNegatives)
    is->num = is->numPositives;
  else
    is->num = is->numPositives + is->numNegatives;

  if(is->num <= 1)
    {
      is->stddev = 0;
    }
  else
    {
      is->mean = is->mean / is->num;
      is->stddev = sqrt((is->sumOfSquares - is->num * is->mean * is->mean) /
                        (is->num - 1.0));
    }
}


void InstrumentCoordStatistics(VA_TYPE(CDS_COORD_t) * va,
                               InstrumenterStatistics * is,
                               int separateNegatives)
{
  int32 i;
  CDS_COORD_t * var;

  memset(is, 0, sizeof(InstrumenterStatistics));
  is->min = FLT_MAX;
  is->minNegative = FLT_MAX;

  // use mean to hold sum
  for(i = 0; i < GetNumVA_CDS_COORD_t(va); i++)
    {
      var = GetVA_CDS_COORD_t(va, i);

      if((*var) < 0.0f)
        {
          is->numNegatives++;
          is->minNegative = MIN(is->minNegative, *var);
          if(!separateNegatives)
            AddCoordToInstrumenterStatistics(is, var);
        }
      else
        {
          is->numPositives++;
          AddCoordToInstrumenterStatistics(is, var);
        }
    }

  is->min = (is->min == FLT_MAX) ? 0.0f : is->min;
  if(separateNegatives)
    is->num = is->numPositives;
  else
    is->num = is->numPositives + is->numNegatives;

  if(is->num <= 1)
    {
      is->stddev = 0;
    }
  else
    {
      is->mean = is->mean / is->num;

      is->stddev = sqrt((is->sumOfSquares - is->num * is->mean * is->mean) /
                        (is->num - 1.0));
    }
}


void AddFragDetailsToMateStats(MateStats * ms,
                               VA_TYPE(FragDetail) * fda)
{
  int32 i;
  for(i = 0; i < GetNumVA_FragDetail(fda); i++)
    {
      FragDetail * fd = GetVA_FragDetail(fda, i);
      switch(fd->type)
        {
          case AS_READ:
          case AS_TRNR:
            ms->reads++;
            break;
          case AS_EXTR:
            ms->externalReads++;
            break;
          default:
            ms->externalFrags++;
            break;
        }
    }
}


void AddMateDetailsToMateStatsSet(MateStats * ms, VA_TYPE(MateDetail) * mda)
{
  int32 i;
  for(i = 0; i < GetNumVA_MateDetail(mda); i++)
    {
      MateDetail * md = GetVA_MateDetail(mda, i);
      switch(md->type)
        {
          case AS_READ:
          case AS_TRNR:
            ms->reads++;
            break;
          case AS_EXTR:
            ms->externalReads++;
            break;
          default:
            ms->externalFrags++;
            break;
        }
    }
}


void ComputeMateStats(MateStatsSet * mss, MateStatusPositions * msp)
{
  int ori1;

  AddFragDetailsToMateStats(&(mss->inter), msp->inter);
  for(ori1 = 0; ori1 < NUM_ORIENTATIONS_INSTR; ori1++)
    {
      int ori2;
    
      AddMateDetailsToMateStatsSet(&(mss->happy), msp->happy[ori1]);
      AddMateDetailsToMateStatsSet(&(mss->misseparatedClose),
                                   msp->misseparatedClose[ori1]);
      AddMateDetailsToMateStatsSet(&(mss->misseparatedFar),
                                   msp->misseparatedFar[ori1]);
    
      for(ori2 = 0; ori2 < NUM_ORIENTATIONS_INSTR; ori2++)
        AddMateDetailsToMateStatsSet(&(mss->misoriented),
                                     msp->misoriented[ori1][ori2]);
    }
}


void ComputeMateInstrumenterStats(MateInstrumenter * mi)
{
  ResetMateInstrumenterCounts(mi);

  if(mi->options & INST_OPT_INTRA_MATES)
    ComputeMateStats(&(mi->intra), mi->mateStatus->intra);

  if(mi->options & INST_OPT_INTER_MATES)
    ComputeMateStats(&(mi->inter), mi->mateStatus->inter);
  
  AddFragDetailsToMateStats(&(mi->mateless), mi->noMate);
}


void ComputeUnitigInstrumenterStats(ScaffoldGraphT * graph,
                                    UnitigInstrumenter * ui)
{
  ComputeMateInstrumenterStats(&(ui->mates));
}


void ComputeContigInstrumenterStats(ScaffoldGraphT * graph,
                                    ContigInstrumenter * ci)
{
  InstrumentFloatStatistics(ci->unitigSizes,
                            &(ci->unitigSizeStats),
                            0);
  InstrumentFloatStatistics(ci->surrogateSizes,
                            &(ci->surrogateSizeStats),
                            0);
  ComputeUnitigInstrumenterStats(graph, &(ci->unitig));
  ComputeMateInstrumenterStats(&(ci->mates));
  if(ci->options & INST_OPT_BREAKPOINTS)
    DetectRoughIntraContigBreakpoints(graph, ci);
}


void ComputeScaffoldInstrumenterStats(ScaffoldGraphT * graph,
                                      ScaffoldInstrumenter * si)
{
  InstrumentFloatStatistics(si->scaffoldGapSizes,
                            &(si->scaffoldGapSizeStats),
                            1);
  InstrumentFloatStatistics(si->contigSizes,
                            &(si->contigSizeStats),
                            0);
  ComputeContigInstrumenterStats(graph, &(si->contig));
  ComputeMateInstrumenterStats(&(si->mates));
  if(si->options & INST_OPT_BREAKPOINTS)
    DetectRoughInterContigBreakpoints(graph, si);
}


void ComputeScaffoldGraphInstrumenterStats(ScaffoldGraphT * graph,
                                           ScaffoldGraphInstrumenter * sgi)
{
  if(sgi->options & INST_OPT_FRAGMENTS)
    {
      CDS_CID_t i;
      CDS_CID_t j;
    
      sgi->numNotInUnitigs = 0;
      sgi->numNotInContigs = 0;
      sgi->numNotPlaced = 0;
      sgi->numChaff = 0;
      sgi->numInUnresolvedChunks = 0;
    
      sgi->numFragments = GetNumCIFragTs(graph->CIFrags);
      for( j = 0, i = 0; i < sgi->numFragments; i++)
        {
          CIFragT * frag = GetCIFragT(graph->CIFrags, i);
      
          sgi->numNotInUnitigs += (frag->cid == NULLINDEX) ? 1 : 0;
          sgi->numNotInContigs += (frag->CIid == NULLINDEX) ? 1 : 0;
          sgi->numNotPlaced += (frag->flags.bits.isPlaced) ? 0 : 1;
          sgi->numChaff += (frag->flags.bits.isChaff) ? 1 : 0;
      
        }
    
      // fragments in unresolved unitigs
      for(i = 0; i < GetNumGraphNodes(graph->CIGraph); i++)
        {
          ChunkInstanceT * chunk = GetGraphNode(graph->CIGraph, i);
          if(chunk->scaffoldID == NULLINDEX)
            sgi->numInUnresolvedChunks += chunk->info.CI.numFragments;
        }
    }
  
  InstrumentFloatStatistics(sgi->singletonScaffoldSizes,
                            &(sgi->singletonScaffoldSizeStats),
                            0);
  InstrumentCoordStatistics(sgi->unitigsPerSingletonScaffold,
                            &(sgi->unitigsPerSingletonStats),
                            0);
  InstrumentFloatStatistics(sgi->degenerateScaffoldSizes,
                            &(sgi->degenerateScaffoldSizeStats),
                            0);
  InstrumentFloatStatistics(sgi->scaffoldSizes,
                            &(sgi->scaffoldSizeStats),
                            0);
  ComputeScaffoldInstrumenterStats(graph, &(sgi->scaffold));
}


/********************************************************************
                     Functions for printing
********************************************************************/
void PrintConsensus(VA_TYPE(char) * consensus, FILE * printTo)
{
  int32 i;
  char * seq = GetVA_char(consensus, 0);
  int32 length = GetNumVA_char(consensus);

  for(i = 0; i < length - 1; i++)
    {
      if(i > 0 && i % 70 == 0)
        fprintf(printTo, "\n");
      fprintf(printTo, "%c", seq[i]);
    }
  fprintf(printTo, "\n");
}


void PrintFragment(CIFragT * frag, CDS_CID_t index, FILE * printTo)
{
  fprintf(printTo, "Fragment iid " F_CID ", index " F_CID "\n",
          frag->iid, index);

  if(frag->flags.bits.hasMate){
    fprintf(printTo, "  index of mate "F_CID"\n", frag->mateOf);
    fprintf(printTo, "  iid of mate " F_CID "\n",
	    GetCIFragT(ScaffoldGraph->CIFrags, frag->mateOf)->iid);
  } else {
    fprintf(printTo, "  unmated");
  }
  switch(frag->type)
    {
      case AS_READ:
        fprintf(printTo, "  type: AS_READ; ");
        break;
      case AS_EXTR:
        fprintf(printTo, "  type: AS_EXTR; ");
        break;
      case AS_TRNR:
        fprintf(printTo, "  type: AS_TRNR; ");
        break;
      default:
        fprintf(printTo, "  type: -other-; ");
        break;
    }
  fprintf(printTo, "%s, %s, %s\n",
          (frag->flags.bits.isPlaced) ? "placed" : "not placed",
          (frag->flags.bits.isSingleton) ? "singleton" : "not singleton",
          (frag->flags.bits.isChaff) ? "chaff" : "not chaff");
  fprintf(printTo, "  cid " F_CID "; CIid " F_CID "; 5p,3p: %d,%d\n",
          frag->cid, frag->CIid,
          (int) frag->offset5p.mean, (int) frag->offset3p.mean);
  fprintf(printTo, "  contigID " F_CID "; 5p,3p: %d,%d\n",
          frag->contigID,
          (int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean);
}


void PrintContigPlacement(ContigPlacement * cp,
                          CDS_COORD_t inFromLeft,
                          char * prefix,
                          FILE * printTo)
{
  if(inFromLeft != NULLINDEX)
    fprintf(printTo, "%s" F_COORD ". ", prefix, inFromLeft);
  else
    fprintf(printTo, "%s ", prefix);
  fprintf(printTo,
          "IID:" F_CID ", Offset:%.f, Length:%.f, Orientation:%s\n",
          cp->id, cp->offset, cp->length,
          (cp->orient == A_B) ? "A_B" : "B_A");
}


void PrintContigPair(IntContigPairs * cp, char * prefix, FILE * printTo)
{
  fprintf(printTo,
          "%scontig1 = " F_CID ", contig2 = " F_CID ", mean = %f, stddev = %f, orientation = %c\n",
          prefix, cp->contig1, cp->contig2,
          cp->mean, cp->stddev, cp->orient);
}


void PrintInstrumenterStats(InstrumenterStatistics * is,
                            char * prefix,
                            FILE * printTo)
{
  fprintf(printTo, "%sNumber of negatives: %d\n",
          prefix, is->numNegatives);
  fprintf(printTo, "%sNumber of positives: %d\n",
          prefix, is->numPositives);
  
  fprintf(printTo, "%sMin: %.2f, Mean: %.2f, Max: %.2f, Stddev: %.2f",
          prefix, is->min, is->mean, is->max, is->stddev);
  
  if(is->numNegatives != 0 && is->num == is->numPositives)
    {
      fprintf(printTo, ", Min negative: %.2f\n", is->minNegative);
      fprintf(printTo, "%s(Negative values not used in statistics)\n",
              prefix);
    }
  else
    {
      fprintf(printTo, "\n");
    }
}


void PrintMateDetailAndDist(MateDetail * md,
                            DistT * dptr,
                            char * prefix,
                            FILE * printTo)
{
  fprintf(printTo,
          "%sIn (" F_CID "," F_CID ") fid: " F_CID ", 5p: %.f\tmid: " F_CID ", 5p: %.f, type: %c, dist: %.2f, stddev: %.2f\n",
          prefix,
          md->fragChunkIID, md->mateChunkIID,
          md->fragIID, md->fragOffset5p,
          md->mateIID, md->mateOffset5p,
          md->type, dptr->mu, dptr->sigma);
}



#ifdef TRACK_3P

void safelyAppendInstInfo(char **locs,int32 utgIID, int *lenloc, int *lenUsed){
  ChunkInstanceT *mateUtg;
  int32 mateCtg;
  int32 mateScf;
  char teststring[100];
  int testsize;
   
  mateUtg = GetGraphNode(ScaffoldGraph->CIGraph,utgIID);
  mateCtg = mateUtg->info.CI.contigID;
  mateScf = mateUtg->scaffoldID;
  testsize = snprintf(teststring,99," surroTig " F_CID " ctg " F_CID " scf " F_CID,
		      mateUtg->id,mateCtg,mateScf);
  assert(testsize <= 100); /* test against truncation */
  assert(testsize >0); /* test against other error */
  if(*lenUsed+testsize>*lenloc){
    *lenloc+=1000;
    *locs = (char *) safe_realloc(*locs, *lenloc * sizeof(char));
  }
  strcat(*locs,teststring);
  *lenUsed+=testsize-1; /* -1 because snprintf includes the '\0' in its return,
			   but strcat effectively adds one less since the '\0'
			   at the previous string end is overwritten */
}


void PrintExternalMateDetailAndDist(MateDetail * md,
                                    DistT * dptr,
                                    char * prefix,
                                    FILE * printTo,
                                    int printtype)
{
  if(printtype==PRINTTABLE){
    fprintf(printTo,
	    "%sIn (" F_CID "," F_CID ") fid: " F_CID ", 5p: %.f\tmid: " F_CID ", 5p: %.f, type: %c, dist: %.2f, stddev: %.2f\n",
	    prefix,
	    md->fragChunkIID, md->mateChunkIID,
	    md->fragIID, md->fragOffset5p,
	    md->mateIID, md->mateOffset5p,
	    md->type, dptr->mu, dptr->sigma);
  } else {
    int fragLeftEnd,fragRightEnd,fragOri;
    int32 mateChunk,mateScf,mateCtg;
    assert(PRINTCELAMY==printtype);

    if(! USE_ALL_MATES && ! (USE_LONG_MATES && dptr->mu > 15000))
      return;

    if(md->fragOffset5p<md->fragOffset3p){
      fragLeftEnd = md->fragOffset5p;
      fragRightEnd = md->fragOffset3p;
      fragOri = A_B;
    } else {
      fragRightEnd = md->fragOffset5p;
      fragLeftEnd = md->fragOffset3p;
      fragOri = B_A;
    }
    mateChunk = md->mateChunkIID;
    {
      ChunkInstanceT *unitig = GetGraphNode(ScaffoldGraph->CIGraph,mateChunk);
      int numInst = unitig->info.CI.numInstances;
      if(unitig->info.CI.numInstances>0){
	static char *locs=NULL;
	static int lenloc=0;
	int lenUsed = 0;
	if(locs==NULL){
	  lenloc = 1000;
	  locs = (char *) safe_malloc(lenloc*sizeof(char));
	}
	locs[0]='\0';
	if(numInst<=2){
	  safelyAppendInstInfo(&locs,unitig->info.CI.instances.in_line.instance1,&lenloc,&lenUsed);
	  if(numInst==2){
	    safelyAppendInstInfo(&locs,unitig->info.CI.instances.in_line.instance2,&lenloc,&lenUsed);
	  }
	} else {
	  int i,n;
          int32 *inst_list;
	  n=unitig->info.CI.numInstances;
	  assert(n == GetNumint32s(unitig->info.CI.instances.va));
	  inst_list = Getint32(unitig->info.CI.instances.va,0);
	  for(i=0;i<n;i++){
	    safelyAppendInstInfo(&locs,inst_list[i],&lenloc,&lenUsed);
	  }
	}

	fprintf(printTo,F_CID "Fragment: " F_COORD " %s " F_COORD " R50 # Externally-mated fragment " F_CID " ori:%s lib %f +/- %f Mate info: BaseCI " F_CID " Instances(%d): %s\n",
		md->fragIID,fragLeftEnd,
		fragOri==A_B ? "A7CMColor" : "A8CMColor" ,
		fragRightEnd,md->fragIID,fragOri==A_B?"A_B":"B_A",
		dptr->mu,dptr->sigma,
		unitig->id,
		unitig->info.CI.numInstances,
		locs);
      } else {
	mateCtg = unitig->info.CI.contigID;
	mateScf = unitig->scaffoldID;
	fprintf(printTo,F_CID "Fragment: " F_COORD " %s " F_COORD " R50 # Externally-mated fragment " F_CID " ori:%s lib %f +/- %f mateChunk " F_CID " mateCtg %d mateScf %d\n",
		md->fragIID,fragLeftEnd,
		fragOri==A_B ? "A7CMColor" : "A8CMColor" ,
		fragRightEnd,md->fragIID,fragOri==A_B?"A_B":"B_A",
		dptr->mu,dptr->sigma,
		mateChunk,mateCtg,mateScf);
      }
    }
  }
}
#endif


/*
  For debugging - call from gdb, for instance
*/
void PrintMateDetailsAndDists(ScaffoldGraphT * graph,
                              VA_TYPE(MateDetail) * mda,
                              char * prefix,
                              FILE * printTo)
{
  int32 i;
  
  for(i = 0; i < GetNumVA_MateDetail(mda); i++)
    {
      MateDetail * md = GetVA_MateDetail(mda, i);
      DistT * dptr = GetDistT(graph->Dists, md->libIID);
      PrintMateDetailAndDist(md, dptr, prefix, printTo);
    }
}

#ifdef TRACK_3P
/*
  For debugging - call from gdb, for instance
*/

void PrintExternalMateDetailsAndDists(ScaffoldGraphT * graph,
                                      VA_TYPE(MateDetail) * mda,
                                      char * prefix,
                                      FILE * printTo,
                                      int printtype)
{
  int i;
  
  for(i = 0; i < GetNumVA_MateDetail(mda); i++)
    {
      MateDetail * md = GetVA_MateDetail(mda, i);
      DistT * dptr = GetDistT(graph->Dists, md->libIID);
      PrintExternalMateDetailAndDist(md, dptr, prefix, printTo,printtype);
    }
}
#endif


void PrintMateComparison(int32 numBefore,
                         int32 totalBefore,
                         int32 numAfter,
                         int32 totalAfter,
                         char * string,
                         int moreIsBetter,
                         int printRawNumbers,
                         FILE * printTo)
{
  float ratioBefore = (totalBefore > 0) ?
    ((float) numBefore) / totalBefore : 0.0f;
  float ratioAfter = (totalAfter > 0) ?
    ((float) numAfter) / totalAfter : 0.0f;
  
  fprintf(printTo, string);
  if(totalBefore < 1.f || totalAfter < 1.f)
    fprintf(printTo, "%s", "-      ");
  else if(ratioAfter > ratioBefore + 0.000001)
    fprintf(printTo, "%s", (moreIsBetter) ? "better " : "worse  ");
  else if( ratioAfter < ratioBefore - 0.00001)
    fprintf(printTo, "%s", (moreIsBetter) ? "worse  " : "better ");
  else
    fprintf(printTo, "same   ");

  if(printRawNumbers)
    fprintf(printTo, "(%8d = %7.3f%% before, %8d = %7.3f%% after)\n",
            (int) numBefore, 100.f * ratioBefore,
            (int) numAfter, 100.f * ratioAfter);
  else
    fprintf(printTo, "(%7.3f%% before, %7.3f%% after)\n",
            100.f * ratioBefore, 100.f * ratioAfter);
}


int32 GetMateStatsSum(MateStats * ms)
{
  return(ms->reads + ms->externalReads + ms->externalFrags);
}


void PrintMateStatsSet(MateStatsSet * mss,
                       InstrumenterLevel level,
                       int inter,
                       char * prefix,
                       FILE * printTo)
{
  char * unit;
  char * nextUnit;
  char * relation;
  char * otherRelation;

  if(level == InstrumenterUnitigLevel)
    {
      unit = "unitig";
      nextUnit = "contig";
    }
  else
    {
      unit = "contig";
      nextUnit = "scaffold";
    }

  if(inter)
    {
      relation = "inter";
      otherRelation = "intra";
    }
  else
    {
      relation = "intra";
      otherRelation = "inter";
    }

  fprintf(printTo, "%sCelera reads:\n", prefix);
  fprintf(printTo, "%s\t%10d happy %s-%s mate pairs\n",
          prefix, mss->happy.reads, relation, unit);
  fprintf(printTo, "%s\t%10d mis-oriented %s-%s mate pairs\n",
          prefix, mss->misoriented.reads, relation, unit);
  fprintf(printTo,
          "%s\t%10d mis-separated too close %s-%s mate pairs\n",
          prefix, mss->misseparatedClose.reads, relation, unit);
  fprintf(printTo,
          "%s\t%10d mis-separated too far %s-%s mate pairs\n",
          prefix, mss->misseparatedFar.reads, relation, unit);
  /*
    fprintf(printTo, "%s\t%10d %s-%s read mates\n",
    prefix, mss->inter.reads, otherRelation, unit);
  */

  fprintf(printTo, "%sExternal reads:\n", prefix);
  fprintf(printTo, "%s\t%10d happy %s-%s mate pairs\n",
          prefix, mss->happy.externalReads, relation, unit);
  fprintf(printTo, "%s\t%10d mis-oriented %s-%s mate pairs\n",
          prefix, mss->misoriented.externalReads, relation, unit);
  fprintf(printTo,
          "%s\t%10d mis-separated too close %s-%s mate pairs\n",
          prefix, mss->misseparatedClose.externalReads, relation, unit);
  fprintf(printTo,
          "%s\t%10d mis-separated too far %s-%s mate pairs\n",
          prefix, mss->misseparatedFar.externalReads, relation, unit);
  /*
    fprintf(printTo, "%s\t%10d %s-%s mates\n",
    prefix, mss->inter.externalReads, otherRelation, unit);
  */

  fprintf(printTo, "%sExternal frags:\n", prefix);
  fprintf(printTo, "%s\t%10d happy %s-%s mate pairs\n",
          prefix, mss->happy.externalFrags, relation, unit);
  fprintf(printTo, "%s\t%10d mis-oriented %s-%s mate pairs\n",
          prefix, mss->misoriented.externalFrags, relation, unit);
  fprintf(printTo,
          "%s\t%10d mis-separated too close %s-%s mate pairs\n",
          prefix, mss->misseparatedClose.externalFrags, relation, unit);
  fprintf(printTo,
          "%s\t%10d mis-separated too far %s-%s mate pairs\n",
          prefix, mss->misseparatedFar.externalFrags, relation, unit);
  /*
    fprintf(printTo, "%s\t%10d %s-%s mates\n",
    prefix, mss->inter.externalFrags, otherRelation, unit);
  */

  fprintf(printTo, "%sTotal:\n", prefix);
  fprintf(printTo, "%s\t%10d happy %s-%s mate pairs\n",
          prefix, GetMateStatsSum(&(mss->happy)),
          relation, unit);
  fprintf(printTo, "%s\t%10d mis-oriented %s-%s mate pairs\n",
          prefix, GetMateStatsSum(&(mss->misoriented)),
          relation, unit);
  fprintf(printTo,
          "%s\t%10d mis-separated too close %s-%s mate pairs\n",
          prefix, GetMateStatsSum(&(mss->misseparatedClose)),
          relation, unit);
  fprintf(printTo,
          "%s\t%10d mis-separated too far %s-%s mate pairs\n",
          prefix, GetMateStatsSum(&(mss->misseparatedFar)),
          relation, unit);
  /*
    fprintf(printTo, "%s\t%10d %s-%s mates\n",
    prefix, GetMateStatsSum(&(mss->inter)),
    otherRelation, unit);
  */
}


void PrintMateInstrumenter(MateInstrumenter * mi,
                           InstrumenterLevel level,
                           char * prefix,
                           FILE * printTo)
{
  if(mi->options & INST_OPT_INTRA_MATES)
    PrintMateStatsSet(&(mi->intra), level, 0, prefix, printTo);
  
  if(mi->options & INST_OPT_INTER_MATES)
    PrintMateStatsSet(&(mi->inter), level, 1, prefix, printTo);

  // print mateless
  fprintf(printTo, "%sMateless:\n", prefix);
  fprintf(printTo, "%s\t%10d reads\n", prefix, mi->mateless.reads);
  fprintf(printTo, "%s\t%10d external reads\n",
          prefix, mi->mateless.externalReads);
  fprintf(printTo, "%s\t%10d external frags\n",
          prefix, mi->mateless.externalFrags);
  fprintf(printTo, "%s\t%10d total\n", prefix, GetMateStatsSum(&(mi->mateless)));
}


void PrintInstrumenterBookkeeping(InstrumenterBookkeeping * bk,
                                  int printMissingMates,
                                  char * prefix,
                                  FILE * printTo)
{
}


void PrintBreakpoint(InstrumenterBreakpoint * bp,
                     char * prefix,
                     FILE * printTo)
{
  fprintf(printTo, "%s" F_CID " (%c)\t" F_CID "\t" F_CID "\t(" F_COORD "," F_COORD ")\t%6d\t\t",
          prefix,
          bp->iid,
          ((bp->section == BP_ALL) ? 'a' :
           ((bp->section == BP_LEFT) ? 'l' :
            (bp->section == BP_RIGHT) ? 'r' : 'm')),
          bp->contig1, bp->contig2,
          bp->end1, bp->end2, bp->pairs);
  switch(bp->type)
    {
      case INST_BP_TOO_CLOSE:
        fprintf(printTo, "too close\n");
        break;
      case INST_BP_TOO_FAR:
        fprintf(printTo, "too far\n");
        break;
      case INST_BP_NORMAL:
        fprintf(printTo, "mis-oriented normal\n");
        break;
      case INST_BP_ANTINORMAL:
        fprintf(printTo, "mis-oriented antinormal\n");
        break;
      case INST_BP_OUTTIE:
        fprintf(printTo, "mis-oriented outtie\n");
        break;
      case INST_BP_INNIE:
        fprintf(printTo, "mis-oriented innie\n");
        break;
      case INST_BP_UNKNOWN:
        fprintf(printTo, "unknown!\n");
        break;
    }
}


void PrintBreakpoints(VA_TYPE(InstrumenterBreakpoint) * bps,
                      char * prefix,
                      FILE * printTo)
{
  int32 i;
  if(GetNumVA_InstrumenterBreakpoint(bps) > 0)
    {
      fprintf(printTo,
              "%sIID\tContig1\tContig2\tInterval\tMate Pairs\tType\n", prefix);
      for(i = 0; i < GetNumVA_InstrumenterBreakpoint(bps); i++)
        {
          PrintBreakpoint(GetVA_InstrumenterBreakpoint(bps, i),
                          prefix,
                          printTo);
        }
    }
  else
    {
      fprintf(printTo,"%sNone\n", prefix);
    }
}


void PrintUnitigInstrumenter(ScaffoldGraphT * graph,
                             UnitigInstrumenter * ui,
                             InstrumenterVerbosity verbose,
                             char * prefix,
                             FILE * printTo)
{
  if(ui->id != NULLINDEX)
    {
      fprintf(printTo, "%sStatistics for unitig " F_CID "\n", prefix, ui->id);
      fprintf(printTo, "%sSize: " F_COORD "\n", prefix, ui->rightEnd - ui->leftEnd);
    }
  
  fprintf(printTo,
          "\n%s%d reads, %d external reads, %d external frags\n",
          prefix, ui->numReads, ui->numExtReads, ui->numExtFrags);

  fprintf(printTo, "\n%sMate summary:\n", prefix);
  PrintMateInstrumenter(&(ui->mates), InstrumenterUnitigLevel,
                        prefix, printTo);

  fprintf(printTo, "\n%sExternal mates:\n", prefix);
  PrintInstrumenterBookkeeping(&(ui->bookkeeping), FALSE, prefix, printTo);

  if(ui->options & INST_OPT_BREAKPOINTS)
    {
      fprintf(printTo, "\n%sBreakpoints:\n", prefix);
      PrintBreakpoints(ui->breakpoints, prefix, printTo);
    }

  if(verbose >= InstrumenterVerbose5)
    {
      int32 i;
      ChunkInstanceT * unitig;
    
      if((unitig = GetGraphNode(graph->RezGraph, ui->id)) == NULL)
        {
          fprintf(stderr, "Unitig " F_CID " does not exist in the graph!\n", ui->id);
          return;
        }
      else
        {
          MultiAlignT * cma = LoadMultiAlignTFromSequenceDB(graph->sequenceDB,
                                                            unitig->id, FALSE);
          if(cma == NULL)
            {
              fprintf(stderr,
                      "Failed to load MultiAlignT of unitig " F_CID "\n", unitig->id);
              return;
            }

          fprintf(printTo, "\n%sGapped sequence:\n", prefix);
          fprintf(printTo, "\n>Unitig " F_CID "\n", unitig->id);
          PrintConsensus(cma->consensus, printTo);
          UnloadMultiAlignTFromSequenceDB(graph->sequenceDB, unitig->id, FALSE);
        }

      fprintf(printTo, "\n%sConstituent fragments:\n", prefix);
      for(i = 0; i < GetNumVA_CDS_CID_t(ui->bookkeeping.fragArray); i++)
        {
          // Don't use the convenience function, since we want to print the index
          CDS_CID_t * iid = GetVA_CDS_CID_t(ui->bookkeeping.fragArray, i);
          InfoByIID * info = GetInfoByIID(graph->iidToFragIndex, *iid);
          CIFragT * frag = GetCIFragT(graph->CIFrags, info->fragIndex);
          PrintFragment(frag, info->fragIndex, printTo);
        }
    }
  
  fprintf(printTo, "\n");
  fflush(printTo);
}


void PrintContigInstrumenter(ScaffoldGraphT * graph,
                             ContigInstrumenter * ci,
                             InstrumenterVerbosity verbose,
                             char * prefix,
                             FILE * printTo)
{
  if(ci->id != NULLINDEX)
    {
      fprintf(printTo, "%sStatistics for contig " F_CID "\n", prefix, ci->id);
      fprintf(printTo, "%sSize: " F_COORD "\n", prefix, ci->rightEnd - ci->leftEnd);
      fprintf(printTo, "%sLeft end: " F_COORD ", Right end: " F_COORD ", Orientation: %s\n",
              prefix, ci->leftEnd, ci->rightEnd,
              (ci->orientation == A_B) ? "A_B" : "B_A");
    }
  
  fprintf(printTo, "\n%sUnitig sizes:\n", prefix);
  PrintInstrumenterStats(&(ci->unitigSizeStats), prefix, printTo);
  fprintf(printTo, "\n%sSurrogate sizes:\n", prefix);
  PrintInstrumenterStats(&(ci->surrogateSizeStats), prefix, printTo);

  fprintf(printTo,
          "\n%s%d reads, %d external reads, %d external frags\n",
          prefix, ci->numReads, ci->numExtReads, ci->numExtFrags);

  fprintf(printTo, "\n%sMate summary:\n", prefix);
  PrintMateInstrumenter(&(ci->mates), InstrumenterContigLevel,
                        prefix, printTo);

  fprintf(printTo, "\n%sExternal mates:\n", prefix);
  PrintInstrumenterBookkeeping(&(ci->bookkeeping), FALSE, prefix, printTo);

  if(ci->options & INST_OPT_BREAKPOINTS)
    {
      fprintf(printTo, "\n%sBreakpoints:\n", prefix);
      PrintBreakpoints(ci->breakpoints, prefix, printTo);
    }

  if(verbose >= InstrumenterVerbose5)
    {
      int32 i;
      ContigT * contig;
    
      if((contig = GetGraphNode(graph->RezGraph, ci->id)) == NULL)
        {
          fprintf(stderr, "Contig " F_CID " does not exist in the graph!\n", ci->id);
          return;
        }
      else
        {
          MultiAlignT * cma = LoadMultiAlignTFromSequenceDB(graph->sequenceDB,
                                                            contig->id, FALSE);
          if(cma == NULL)
            {
              fprintf(stderr,
                      "Failed to load MultiAlignT of contig " F_CID "\n", contig->id);
              return;
            }

          fprintf(printTo, "\n%sGapped sequence:\n", prefix);
          fprintf(printTo, "\n>Contig " F_CID "\n", contig->id);
          PrintConsensus(cma->consensus, printTo);
          UnloadMultiAlignTFromSequenceDB(graph->sequenceDB, contig->id, FALSE);
        }

      fprintf(printTo, "\n%sConstituent fragments:\n", prefix);
      for(i = 0; i < GetNumVA_CDS_CID_t(ci->bookkeeping.fragArray); i++)
        {
          // Don't use the convenience function, since we want to print the index
          CDS_CID_t * iid = GetVA_CDS_CID_t(ci->bookkeeping.fragArray, i);
          InfoByIID * info = GetInfoByIID(graph->iidToFragIndex, *iid);
          CIFragT * frag = GetCIFragT(graph->CIFrags, info->fragIndex);
          PrintFragment(frag, info->fragIndex, printTo);
        }
    }
  
  fprintf(printTo, "\n");
  fflush(printTo);
}


void PrintInferredStddevs(VA_TYPE(float) * stddevs,
                          char * prefix,
                          FILE * printTo)
{
  int32 incRun = 0;
  int32 decRun = 0;
  int32 maxIncRun = 0;
  int32 maxDecRun = 0;
  int increasing = 0;
  float * prevStddev = NULL;
  float * currStddev = NULL;
  
  if(GetNumVA_float(stddevs) > 0)
    {
      int32 i;
      fprintf(printTo, "%s\t", prefix);
      for( i = 0; i < GetNumVA_float(stddevs); i++)
        {
          currStddev = GetVA_float(stddevs, i);
          fprintf(printTo, "%.0lf ", *currStddev);
          if(prevStddev != NULL)
            {
              if(*currStddev > *prevStddev)
                {
                  incRun++;
                  maxIncRun = MAX(maxIncRun, incRun);
                  if(!increasing)
                    {
                      increasing = 1;
                      decRun = 0;
                    }
                }
              else
                {
                  decRun++;
                  maxDecRun = MAX(maxDecRun, decRun);
                  if(increasing)
                    {
                      increasing = 0;
                      incRun = 0;
                    }
                }
            }
          prevStddev = currStddev;
        }
      fprintf(printTo, "(%d,%d)\n", maxDecRun, maxIncRun);
    }
  else
    {
      fprintf(printTo, "%sNone.\n", prefix);
    }
}


void PrintScaffoldGaps(ScaffoldInstrumenter * si,
                       FILE * printTo)
{
  int32 numGapSizes = GetNumVA_float(si->scaffoldGapSizes);

  if(numGapSizes > 0)
    {
      int32 i;
      float * gapSizes = GetVA_float(si->scaffoldGapSizes, 0);
      for(i = 0; i < numGapSizes; i++)
        {
          fprintf(printTo, F_CID "\t%d\n", si->id, (int) gapSizes[i]);
        }
    }
}


void PrintUnanchoredContigIDs(ScaffoldInstrumenter * si,
                              char * prefix,
                              FILE * printTo)
{
  CDS_CID_t i;
  // loop over contig pairs & look up in anchoredHT

  for(i = 0; i < GetNumVA_ContigPlacement(si->cpArray); i++)
    {
      ContigPlacement * cp = GetVA_ContigPlacement(si->cpArray, i);
      if(!ExistsInHashTable_AS(si->anchoredHT, (uint64)cp->id, 0))
        {
          PrintContigPlacement(cp, i, prefix, printTo);
        }
    }
  fprintf(printTo, "\n");
}


void PrintScaffoldInstrumenter(ScaffoldGraphT * graph,
                               ScaffoldInstrumenter * si,
                               InstrumenterVerbosity verbose,
                               char * prefix,
                               FILE * printTo)
{
  char nextPrefix[1024];

  sprintf(nextPrefix, "%s\t", prefix);
  if(si->id != NULLINDEX)
    {
      fprintf(printTo, "%sStatistics for scaffold " F_CID "\n", prefix, si->id);
      fprintf(printTo, "\n%sSize: %.0f\n", prefix, si->size);
    }

  fprintf(printTo, "\n%sScaffold gap sizes:\n", prefix);
  PrintInstrumenterStats(&(si->scaffoldGapSizeStats), prefix, printTo);
  
  fprintf(printTo, "\n%sContig sizes:\n", prefix);
  PrintInstrumenterStats(&(si->contigSizeStats), prefix, printTo);

  fprintf(printTo, "\n%sUnanchored contigs:\n", prefix);
  PrintUnanchoredContigIDs(si, prefix, printTo);
  
  if(verbose > InstrumenterVerbose3)
    {
      fprintf(printTo, "\n%sUnitig summary:\n", prefix);
      PrintUnitigInstrumenter(graph, &(si->contig.unitig),
                              MIN(verbose, InstrumenterVerbose4),
                              nextPrefix, printTo);
    }
  
  if(verbose >= InstrumenterVerbose2)
    {
      fprintf(printTo, "\n%sContig summary:\n", prefix);
      PrintContigInstrumenter(graph, &(si->contig),
                              MIN(verbose, InstrumenterVerbose4),
                              nextPrefix, printTo);
    }
  
  fprintf(printTo, "\n%sMate summary:\n", prefix);
  PrintMateInstrumenter(&(si->mates), InstrumenterScaffoldLevel,
                        prefix, printTo);

  fprintf(printTo, "\n%sExternal mates:\n", prefix);
  PrintInstrumenterBookkeeping(&(si->bookkeeping), FALSE, prefix, printTo);

  fprintf(printTo, "\n%sInferred edge stddevs: ", prefix);
  PrintInferredStddevs(si->inferredEdgeStddevs, prefix, printTo);

  if(si->options & INST_OPT_BREAKPOINTS)
    {
      fprintf(printTo, "\n%sBreakpoints:\n", prefix);
      PrintBreakpoints(si->breakpoints, prefix, printTo);
    }
  fprintf(printTo, "\n");
  fflush(printTo);
}


void PrintScaffoldGraphInstrumenter(ScaffoldGraphT * graph,
                                    ScaffoldGraphInstrumenter * sgi,
                                    InstrumenterVerbosity verbose,
                                    FILE * printTo)

{
  fprintf(printTo, "Scaffold graph summary:\n");

  if(sgi->options & INST_OPT_FRAGMENTS)
    {
      fprintf(printTo, "\n%10d fragments\n", sgi->numFragments);

      fprintf(printTo, "%10d not in unitigs\n", sgi->numNotInUnitigs);
      fprintf(printTo, "%10d not in contigs\n", sgi->numNotInContigs);
      fprintf(printTo, "%10d not placed\n", sgi->numNotPlaced);
      fprintf(printTo, "%10d chaff\n", sgi->numChaff);
      fprintf(printTo, "%10d in unresolved chunks & not in scaffolds\n",
              sgi->numInUnresolvedChunks);
      fprintf(printTo, "%10d fragments accessible via scaffolds\n",
              sgi->numFragments - sgi->numInUnresolvedChunks);
    }

  fprintf(printTo, "\nSingleton scaffold sizes:\n");
  PrintInstrumenterStats(&(sgi->singletonScaffoldSizeStats), "", printTo);

  fprintf(printTo, "\nUnitigs per singleton scaffold:\n");
  PrintInstrumenterStats(&(sgi->unitigsPerSingletonStats), "", printTo);
  
  fprintf(printTo, "\nDegenerate scaffold sizes:\n");
  PrintInstrumenterStats(&(sgi->degenerateScaffoldSizeStats), "", printTo);

  fprintf(printTo, "\nScaffold sizes:\n");
  PrintInstrumenterStats(&(sgi->scaffoldSizeStats), "", printTo);

  fprintf(printTo, "\nScaffold summary:\n");
  PrintScaffoldInstrumenter(graph, &(sgi->scaffold), verbose, "\t", printTo);

  fprintf(printTo, "\nExternal mates:\n");
  PrintInstrumenterBookkeeping(&(sgi->bookkeeping), TRUE, "", printTo);
  
  fflush(printTo);
}


/********************************************************************
                    Functions for accessing
********************************************************************/
void CopyMateStats(MateStats * dest,
                   MateStats * src)
{
  dest->reads = src->reads;
  dest->externalReads = src->externalReads;
  dest->externalFrags = src->externalFrags;
}


void CopyMateStatsSet(MateStatsSet * dest,
                      MateStatsSet * src)
{
  CopyMateStats(&(dest->happy), &(src->happy));
  CopyMateStats(&(dest->misoriented), &(src->misoriented));
  CopyMateStats(&(dest->misseparatedClose), &(src->misseparatedClose));
  CopyMateStats(&(dest->misseparatedFar), &(src->misseparatedFar));
  CopyMateStats(&(dest->inter), &(src->inter));
}


void GetMateInstrumenterFromScaffoldInstrumenter(MateInstrumenter * mi,
                                                 ScaffoldInstrumenter * si)
{
  mi->options = si->options;

  CopyMateStatsSet(&(mi->intra), &(si->mates.intra));
  CopyMateStatsSet(&(mi->inter), &(si->mates.inter));
  CopyMateStats(&(mi->mateless), &(si->mates.mateless));
}


void GetMateInstrumenterFromScaffoldGraphInstrumenter(
                                                      MateInstrumenter * mi,
                                                      ScaffoldGraphInstrumenter * sgi)
{
  GetMateInstrumenterFromScaffoldInstrumenter(mi, &(sgi->scaffold));
}


CIFragT * getFragByIID(ScaffoldGraphT * graph,
                       CDS_CID_t iid)
{
  InfoByIID * info = GetInfoByIID(graph->iidToFragIndex, iid);
  return(GetCIFragT(graph->CIFrags, info->fragIndex));
}


CIFragT * getBookkeepingFrag(ScaffoldGraphT * graph,
                             VA_TYPE(CDS_CID_t) * fragArray,
                             CDS_CID_t index)
{
  return(getFragByIID(graph, *(GetVA_CDS_CID_t(fragArray, index))));
}


/********************************************************************
            Functions for adding/comparing instrumenters
********************************************************************/
void AddInstrumenterBreakpoints(VA_TYPE(InstrumenterBreakpoint) * dest,
                                VA_TYPE(InstrumenterBreakpoint) * src)
{
  int32 i;
  
  for(i = 0; i < GetNumVA_InstrumenterBreakpoint(src); i++)
    AppendVA_InstrumenterBreakpoint(dest,
                                    GetVA_InstrumenterBreakpoint(src, i));
}


void AddFragDetails(VA_TYPE(FragDetail) * dest,
                    VA_TYPE(FragDetail) * src)
{
  int32 i;
  
  for(i = 0; i < GetNumVA_FragDetail(src); i++)
    AppendVA_FragDetail(dest, GetVA_FragDetail(src, i));
}


void AddMateDetails(VA_TYPE(MateDetail) * dest,
                    VA_TYPE(MateDetail) * src)
{
  int32 i;
  
  for(i = 0; i < GetNumVA_MateDetail(src); i++)
    AppendVA_MateDetail(dest, GetVA_MateDetail(src, i));
}


int AddMateStatusPositions(MateStatusPositions * dest,
                           MateStatusPositions * src)
{
  int ori1;
  int ori2;

#ifdef DEBUG
  fprintf(stderr, "Adding mate status positions\n");
#endif
  
  AddFragDetails(dest->inter, src->inter);
  
  for(ori1 = 0; ori1 < NUM_ORIENTATIONS_INSTR; ori1++)
    {
      AddMateDetails(dest->happy[ori1], src->happy[ori1]);
    
      for(ori2 = 0; ori2 < NUM_ORIENTATIONS_INSTR; ori2++)
        {
          AddMateDetails(dest->misoriented[ori1][ori2],
                         src->misoriented[ori1][ori2]);
        }
      AddMateDetails(dest->misseparatedClose[ori1],
                     src->misseparatedClose[ori1]);
    
      AddMateDetails(dest->misseparatedFar[ori1],
                     src->misseparatedFar[ori1]);
    }
  return 0;
}


void AddMateStatusPositionsSets(MateStatusPositionsSet * dest,
                                MateStatusPositionsSet * src)
{
  AddMateStatusPositions(dest->intra, src->intra);
  AddMateStatusPositions(dest->inter, src->inter);
}


int AddMateInstrumenters(MateInstrumenter * dest,
                         MateInstrumenter * src)
{
  int32 i;

  AddMateStatusPositionsSets(dest->mateStatus, src->mateStatus);

  for(i = 0; i < GetNumVA_FragDetail(src->noMate); i++)
    {
      AppendVA_FragDetail(dest->noMate, GetVA_FragDetail(src->noMate, i));
    }
  return 0;
}


/*
  if as_is flag is non-zero it means adding like-level bookkeeping
  if zero, ignore src->fragHT
*/
int AddInstrumenterBookkeeping(ScaffoldGraphT * graph,
                               InstrumenterBookkeeping * dest,
                               InstrumenterBookkeeping * src,
                               int as_is)
{
#ifdef DEBUG
  fprintf(stderr, "Adding instrumenter bookkeeping, %s\n",
          (as_is) ? "as is" : "not as is");
#endif

  if(as_is)
    {
      int32 i;
      // do a straight copy, but from array
      for(i = 0; i < GetNumVA_CDS_CID_t(src->fragArray); i++)
        {
          CDS_CID_t * fragIID = GetVA_CDS_CID_t(src->fragArray, i);
      
          if(!ExistsInHashTable_AS(dest->fragHT, (uint64)fragIID, 0))
            {
              CIFragT * frag = getFragByIID(graph, *fragIID);
              InsertInHashTable_AS(dest->fragHT, (uint64)frag->iid, 0, (uint64)frag, 0);
              AppendVA_CDS_CID_t(dest->fragArray, fragIID);
            }
        }
      for(i = 0; i < GetNumVA_MateDetail(src->wExtMates); i++)
        {
          AppendVA_MateDetail(dest->wExtMates,
                              GetVA_MateDetail(src->wExtMates, i));
        }
    }
  else
    {
      int32 i;
      //add fragments with external mates to destination's fragHT
      for(i = 0; i < GetNumVA_MateDetail(src->wExtMates); i++)
        {
          MateDetail * md = GetVA_MateDetail(src->wExtMates, i);
      
          if(!ExistsInHashTable_AS(dest->fragHT, (uint64)md->fragIID, 0))
            {
              // need to find the CIFragT - stable pointer
              CIFragT * frag = getFragByIID(graph, md->fragIID);
        
              InsertInHashTable_AS(dest->fragHT, (uint64)frag->iid, 0, (uint64)frag, 0);
              AppendVA_CDS_CID_t(dest->fragArray, &(md->fragIID));
            }
        }
    }
  return 0;
}


int AddUnitigInstrumenters(ScaffoldGraphT * graph,
                           UnitigInstrumenter * dest,
                           UnitigInstrumenter * src)
{
#ifdef DEBUG
  fprintf(stderr, "Adding unitig instrumenters\n");
#endif
  
  dest->numReads += src->numReads;
  dest->numExtReads += src->numExtReads;
  dest->numExtFrags += src->numExtFrags;

  AddMateInstrumenters(&(dest->mates), &(src->mates));

  AddInstrumenterBookkeeping(graph,
                             &(dest->bookkeeping), &(src->bookkeeping),
                             TRUE);

  if((dest->options & INST_OPT_BREAKPOINTS) &&
     (src->options & INST_OPT_BREAKPOINTS))
    AddInstrumenterBreakpoints(dest->breakpoints, src->breakpoints);
  return 0;
}


int AddUnitigToContigInstrumenter(ScaffoldGraphT * graph,
                                  ContigInstrumenter * ci,
                                  UnitigInstrumenter * ui)
{
#ifdef DEBUG
  fprintf(stderr, "Adding unitig instrumenter to contig instrumenter\n");
#endif

  // add the size to the unitig sizes
  {
    float unitigSize = fabs(ui->leftEnd - ui->rightEnd);
    if(ui->isSurrogate)
      AppendVA_float(ci->surrogateSizes, &unitigSize);
    else
      AppendVA_float(ci->unitigSizes, &unitigSize);
  }

  ci->numReads += ui->numReads;
  ci->numExtReads += ui->numExtReads;
  ci->numExtFrags += ui->numExtFrags;
  
  // add unitig data
  AddUnitigInstrumenters(graph, &(ci->unitig), ui);

  // mate data
  AddMateInstrumenters(&(ci->mates), &(ui->mates));

  // bookkeeping
  AddInstrumenterBookkeeping(graph,
                             &(ci->bookkeeping), &(ui->bookkeeping),
                             FALSE);
  return 0;
}


int AddContigInstrumenters(ScaffoldGraphT * graph,
                           ContigInstrumenter * dest,
                           ContigInstrumenter * src)
{
  int32 i;

#ifdef DEBUG
  fprintf(stderr, "Adding contig instrumenters\n");
#endif

  for(i = 0; i < GetNumVA_float(src->unitigSizes); i++)
    {
      AppendVA_float(dest->unitigSizes,
                           GetVA_float(src->unitigSizes, i));
    }
  for(i = 0; i < GetNumVA_float(src->surrogateSizes); i++)
    {
      AppendVA_float(dest->surrogateSizes,
                           GetVA_float(src->surrogateSizes, i));
    }

  dest->numReads += src->numReads;
  dest->numExtReads += src->numExtReads;
  dest->numExtFrags += src->numExtFrags;

  // add unitig data
  AddUnitigInstrumenters(graph, &(dest->unitig), &(src->unitig));

  AddMateInstrumenters(&(dest->mates), &(src->mates));

  AddInstrumenterBookkeeping(graph,
                             &(dest->bookkeeping), &(src->bookkeeping),
                             TRUE);

  if((dest->options & INST_OPT_BREAKPOINTS) &&
     (src->options & INST_OPT_BREAKPOINTS))
    AddInstrumenterBreakpoints(dest->breakpoints, src->breakpoints);
  return 0;
}


int AddContigToScaffoldInstrumenter(ScaffoldGraphT * graph,
                                    ScaffoldInstrumenter * si,
                                    ContigInstrumenter * ci)
{
  float size;

#ifdef DEBUG
  fprintf(stderr, "Adding contig instrumenter to scaffold instrumenter\n");
#endif

  // accumulate contig sizes
  size = fabs((float) ci->leftEnd - (float) ci->rightEnd);
  AppendVA_float(si->contigSizes, &size);

  // figure out scaffold size
  si->size = MAX(si->size, MAX(ci->leftEnd, ci->rightEnd));

  // add contig data
  AddContigInstrumenters(graph, &(si->contig), ci);

  // mate data
  AddMateInstrumenters(&(si->mates), &(ci->mates));

  // bookkeeping
  AddInstrumenterBookkeeping(graph,
                             &(si->bookkeeping), &(ci->bookkeeping),
                             FALSE);
  return 0;
}


void AddMateStats(MateStats * dest, MateStats * src)
{
  dest->reads += src->reads;
  dest->externalReads += src->externalReads;
  dest->externalFrags += src->externalFrags;
}

void AddMateStatsSet(MateStatsSet * dest, MateStatsSet * src)
{
  AddMateStats(&(dest->happy), &(src->happy));
  AddMateStats(&(dest->misoriented), &(src->misoriented));
  AddMateStats(&(dest->misseparatedClose), &(src->misseparatedClose));
  AddMateStats(&(dest->misseparatedFar), &(src->misseparatedFar));
}


void AddMateInstrumenterCounts(MateInstrumenter * dest,
                               MateInstrumenter * src)
{
  AddMateStatsSet(&(dest->intra), &(src->intra));
  AddMateStatsSet(&(dest->inter), &(src->inter));
  AddMateStats(&(dest->mateless), &(src->mateless));
}


int32 GetMateStatsBad(MateStatsSet * mss)
{
  return(GetMateStatsSum(&(mss->misoriented)) +
         GetMateStatsSum(&(mss->misseparatedClose)) +
         GetMateStatsSum(&(mss->misseparatedFar)));
}


int32 GetMateStatsHappy(MateStatsSet * mss)
{
  return GetMateStatsSum(&(mss->happy));
}


/*
  returns: InstrumenterResult
  intra_inter is 1 if intra only, 2 if inter only, anything else for both
*/
InstrumenterResult CompareMateInstrumenters(MateInstrumenter * miBefore,
                                            MateInstrumenter * miAfter,
                                            InstrumenterVerbosity verbose,
                                            FILE * printTo)
{
  int32 totalHappyBeforeIntra;
  int32 totalUnhappyBeforeIntra;
  int32 totalBeforeIntra;
  int32 totalHappyBeforeInter;
  int32 totalUnhappyBeforeInter;
  int32 totalBeforeInter;
  int32 totalHappyBefore;
  int32 totalUnhappyBefore;
  int32 totalBefore;
  
  int32 totalHappyAfterIntra;
  int32 totalUnhappyAfterIntra;
  int32 totalAfterIntra;
  int32 totalHappyAfterInter;
  int32 totalUnhappyAfterInter;
  int32 totalAfterInter;
  int32 totalHappyAfter;
  int32 totalUnhappyAfter;
  int32 totalAfter;
  float delta;
  InstrumenterResult retVal;

  // before counts
  totalHappyBeforeIntra = GetMateStatsHappy(&(miBefore->intra));
  totalUnhappyBeforeIntra = GetMateStatsBad(&(miBefore->intra));
  totalBeforeIntra = totalHappyBeforeIntra + totalUnhappyBeforeIntra;

  totalHappyBeforeInter = GetMateStatsHappy(&(miBefore->inter));
  totalUnhappyBeforeInter = GetMateStatsBad(&(miBefore->inter));
  totalBeforeInter = totalHappyBeforeInter + totalUnhappyBeforeInter;

  totalHappyBefore = totalHappyBeforeIntra +
    totalHappyBeforeInter;
  totalUnhappyBefore = totalUnhappyBeforeIntra +
    totalUnhappyBeforeInter;
  
  totalBefore = totalBeforeIntra + totalBeforeInter;
  
  // after counts
  totalHappyAfterIntra = GetMateStatsHappy(&(miAfter->intra));
  totalUnhappyAfterIntra = GetMateStatsBad(&(miAfter->intra));
  totalAfterIntra = totalHappyAfterIntra + totalUnhappyAfterIntra;

  totalHappyAfterInter = GetMateStatsHappy(&(miAfter->inter));
  totalUnhappyAfterInter = GetMateStatsBad(&(miAfter->inter));
  totalAfterInter = totalHappyAfterInter + totalUnhappyAfterInter;

  totalHappyAfter = totalHappyAfterIntra +
    totalHappyAfterInter;
  totalUnhappyAfter = totalUnhappyAfterIntra +
    totalUnhappyAfterInter;
  
  totalAfter = totalAfterIntra + totalAfterInter;
  
  // overall happy & unhappy
  if(verbose >= InstrumenterVerbose1)
    {
      // compare the ratios of happy to mis-separated & mis-oriented
      fprintf(printTo, "* MateInstrumenter comparison:\n");

      if(miBefore->options & INST_OPT_ALL_MATES)
        {
          PrintMateComparison(totalHappyBefore, totalBefore,
                              totalHappyAfter, totalAfter,
                              "Happy (all)       :                  ",
                              1, 1, printTo);
    
          // overall unhappy
          PrintMateComparison(totalUnhappyBefore, totalBefore,
                              totalUnhappyAfter, totalAfter,
                              "Unhappy (all)     :                  ",
                              0, 1, printTo);
        }
    }

  // happy/unhappy intra-inter contig
  if(verbose >= InstrumenterVerbose2)
    {
      // intra
      if(miBefore->options & INST_OPT_INTRA_MATES)
        {
          // happy intracontig
          PrintMateComparison(totalHappyBeforeIntra,
                              totalBeforeIntra,
                              totalHappyAfterIntra,
                              totalAfterIntra,
                              "Happy intra-contig:                  ",
                              1, 1, printTo);
      
          // unhappy intracontig
          PrintMateComparison(totalUnhappyBeforeIntra,
                              totalBeforeIntra,
                              totalUnhappyAfterIntra,
                              totalAfterIntra,
                              "Unhappy intra-contig:                ",
                              0, 1, printTo);
        }

      // inter
      if(miBefore->options & INST_OPT_INTER_MATES)
        {
          // happy intercontig
          PrintMateComparison(totalHappyBeforeInter,
                              totalBeforeInter,
                              totalHappyAfterInter,
                              totalAfterInter,
                              "Happy inter-contig:                  ",
                              1, 1, printTo);
      
          // unhappy intercontig
          PrintMateComparison(totalUnhappyBeforeInter,
                              totalBeforeInter,
                              totalUnhappyAfterInter,
                              totalAfterInter,
                              "Unhappy inter-contig:                ",
                              0, 1, printTo);
        }
    }

  // details of unhappy
  if(verbose >= InstrumenterVerbose3)
    {
      // intra
      if(miBefore->options & INST_OPT_INTRA_MATES)
        {
          // misseparated close:
          PrintMateComparison(
                              GetMateStatsSum(&(miBefore->intra.misseparatedClose)),
                              totalBeforeIntra,
                              GetMateStatsSum(&(miAfter->intra.misseparatedClose)),
                              totalAfterIntra,
                              "Misseparated too close intra-contig: ",
                              0, 1, printTo);
          // misseparated far:
          PrintMateComparison(
                              GetMateStatsSum(&(miBefore->intra.misseparatedFar)),
                              totalBeforeIntra,
                              GetMateStatsSum(&(miAfter->intra.misseparatedFar)),
                              totalAfterIntra,
                              "Misseparated too far intra-contig:   ",
                              0, 1, printTo);
          // misoriented:
          PrintMateComparison(
                              GetMateStatsSum(&(miBefore->intra.misoriented)),
                              totalBeforeIntra,
                              GetMateStatsSum(&(miAfter->intra.misoriented)),
                              totalAfterIntra,
                              "Misoriented intra-contig:            ",
                              0, 1, printTo);
        }

      // inter
      if(miBefore->options & INST_OPT_INTER_MATES)
        {
          // misseparated close:
          PrintMateComparison(
                              GetMateStatsSum(&(miBefore->inter.misseparatedClose)),
                              totalBeforeInter,
                              GetMateStatsSum(&(miAfter->inter.misseparatedClose)),
                              totalAfterInter,
                              "Misseparated too close inter-contig: ",
                              0, 1, printTo);
          // misseparated far:
          PrintMateComparison(
                              GetMateStatsSum(&(miBefore->inter.misseparatedFar)),
                              totalBeforeInter,
                              GetMateStatsSum(&(miAfter->inter.misseparatedFar)),
                              totalAfterInter,
                              "Misseparated too far inter-contig:   ",
                              0, 1, printTo);
          // misoriented:
          PrintMateComparison(
                              GetMateStatsSum(&(miBefore->inter.misoriented)),
                              totalBeforeInter,
                              GetMateStatsSum(&(miAfter->inter.misoriented)),        
                              totalAfterInter,
                              "Misoriented inter-contig:            ",
                              0, 1, printTo);
        }
    }

  // good to bad ratio
  if(verbose >= InstrumenterVerbose1)
    {
      // intra
      if(miBefore->options & INST_OPT_INTRA_MATES)
        {
          PrintMateComparison(totalUnhappyBeforeIntra,
                              totalHappyBeforeIntra,
                              totalUnhappyAfterIntra,
                              totalHappyAfterIntra,
                              "Unhappy/Happy ratio intra-contig:    ",
                              0, 0, printTo);
        }

      // inter
      if(miBefore->options & INST_OPT_INTER_MATES)
        {
          PrintMateComparison(totalUnhappyBeforeInter,
                              totalHappyBeforeInter,
                              totalUnhappyAfterInter,
                              totalHappyAfterInter,
                              "Unhappy/Happy ratio inter-contig:    ",
                              0, 0, printTo);
        }

      // all
      if(miBefore->options & INST_OPT_ALL_MATES)
        {
          PrintMateComparison(totalUnhappyBefore,
                              totalHappyBefore,
                              totalUnhappyAfter,
                              totalHappyAfter,
                              "Unhappy/Happy ratio:                 ",
                              0, 0, printTo);
        }
    }

  // set the return value
  retVal = InstrumenterSame;
  if(!(miBefore->options & INST_OPT_INTER_MATES))
    {
      delta = ((totalBeforeIntra == 0) ?
               1.f : ((float) totalHappyBeforeIntra) / totalBeforeIntra) -
        ((totalAfterIntra == 0) ?
         1.f : ((float) totalHappyAfterIntra) / totalAfterIntra);
      if(totalBeforeIntra < .9f || totalAfterIntra < .9f)
        retVal = InstrumenterIndeterminate;
    }
  else if(!(miBefore->options & INST_OPT_INTRA_MATES))
    {
      delta = ((totalBeforeInter == 0) ?
               1.f : ((float) totalHappyBeforeInter) / totalBeforeInter) -
        ((totalAfterInter == 0) ?
         1.f : ((float) totalHappyAfterInter) / totalAfterInter);
      if(totalBeforeInter < .9f || totalAfterInter < .9f)
        retVal = InstrumenterIndeterminate;
    }
  else
    {
      delta = ((totalBefore == 0) ?
               1.f : ((float) totalHappyBefore) / totalBefore) -
        ((totalAfter == 0) ?
         1.f : ((float) totalHappyAfter) / totalAfter);
      if(totalBefore < .9f || totalAfter < .9f)
        retVal = InstrumenterIndeterminate;
    }
  
  // set return value
  if(retVal != InstrumenterIndeterminate)
    {
      if(delta > 0.00001)
        retVal = InstrumenterWorse;
      else if(delta < -0.00001)
        retVal = InstrumenterBetter;
    }            
  
  return retVal;
}


int AddScaffoldInstrumenters(ScaffoldGraphT * graph,
                             ScaffoldInstrumenter * dest,
                             ScaffoldInstrumenter * src)
{
  int32 i;

#ifdef DEBUG
  fprintf(stderr, "Adding scaffold instrumenters\n");
#endif

  for(i = 0; i < GetNumVA_float(src->scaffoldGapSizes); i++)
    {
      AppendVA_float(dest->scaffoldGapSizes,
                           GetVA_float(src->scaffoldGapSizes, i));
    }

  for(i = 0; i < GetNumVA_float(src->contigSizes); i++)
    {
      AppendVA_float(dest->contigSizes,
                           GetVA_float(src->contigSizes, i));
    }

  AddContigInstrumenters(graph, &(dest->contig), &(src->contig));

  AddMateInstrumenters(&(dest->mates), &(src->mates));

  AddInstrumenterBookkeeping(graph,
                             &(dest->bookkeeping), &(src->bookkeeping),
                             TRUE);
  return 0;
}


int AddScaffoldToScaffoldGraphInstrumenter(ScaffoldGraphT * graph,
                                           ScaffoldGraphInstrumenter * sgi,
                                           ScaffoldInstrumenter * si)
{
#ifdef DEBUG
  fprintf(stderr,
          "Adding scaffold instrumenter to scaffold graph instrumenter.\n");
#endif

  // distinguish between degenerate, singleton, & other scaffolds
  if(GetNumVA_float(si->contigSizes) == 1)
    {
      if(GetNumVA_float(si->contig.unitigSizes) == 1)
        {
          // degenerate
          AppendVA_float(sgi->degenerateScaffoldSizes, &(si->size));
        }
      else
        {
          // singleton
          int32 numUnitigs = GetNumVA_float(si->contig.unitigSizes);
          AppendVA_float(sgi->singletonScaffoldSizes, &(si->size));
          AppendVA_int32(sgi->unitigsPerSingletonScaffold, &numUnitigs);
        }
    }
  else
    {
      // not degenerate or singleton
      AppendVA_float(sgi->scaffoldSizes, &(si->size));
    }

  AddScaffoldInstrumenters(graph, &(sgi->scaffold), si);

  AddInstrumenterBookkeeping(graph,
                             &(sgi->bookkeeping), &(si->bookkeeping),
                             TRUE);
  
  // doing nothing with surrogateTracker yet
  return 0;
}


void SortMateDetails(VA_TYPE(MateDetail) * mda)
{
  if(GetNumVA_MateDetail(mda) > 1)
    qsort(GetVA_MateDetail(mda, 0),
          GetNumVA_MateDetail(mda),
          sizeof(MateDetail),
          (int (*) (const void *, const void *)) md2Compare );
}

/*
  This is a placeholder - simple function
  Should evolve into more of a clustering function
*/
int AppendOrientedCP(ScaffoldGraphT * graph,
                     VA_TYPE(MateDetail) * mdaa[NUM_ORIENTATIONS_INSTR],
                     int * indices,
                     CP_Index * cpip,
                     VA_TYPE(InstrumenterContigPair) * cps)
{
  MateDetail * md;
  InstrumenterContigPair cp;
  int i;

  for(i = 0; i < NUM_ORIENTATIONS_INSTR; i++)
    {
      VA_TYPE(MateDetail) * mda = mdaa[i];
    
      cp.numPairs = 0;
      while((md = GetVA_MateDetail(mda, indices[i])) != NULL &&
            md->fragChunkIID == cpip->contig1 &&
            md->mateChunkIID == cpip->contig2)
        {
          if(cp.numPairs == 0)
            {
              // start a new set
        
              cp.numPairs++;
              indices[i]++;
            }
          else
            {
        
            }
        }
    }
  return 0;
}
  
int AppendMisorientedCP(ScaffoldGraphT * graph,
                        VA_TYPE(MateDetail) * mdaa[NUM_ORIENTATIONS_INSTR],
                        int * indices,
                        CP_Index * cpip,
                        VA_TYPE(InstrumenterContigPair) * cps)
{
  return 0;
}

int AppendContigPairSets(ScaffoldGraphT * graph,
                         VA_TYPE(CP_Index) * cpi,
                         VA_TYPE(InstrumenterContigPair) * cps,
                         MateStatusPositions * msp)
{
  int i, j;
  int hI[NUM_ORIENTATIONS_INSTR];
  int oI[NUM_ORIENTATIONS_INSTR][NUM_ORIENTATIONS_INSTR];
  int cI[NUM_ORIENTATIONS_INSTR];
  int fI[NUM_ORIENTATIONS_INSTR];

  // set indexes into mate detail arrays...
  memset(hI, 0, NUM_ORIENTATIONS_INSTR * sizeof(int));
  for(i = 0; i < NUM_ORIENTATIONS_INSTR; i++)
    memset(oI[i], 0, NUM_ORIENTATIONS_INSTR* sizeof(int));
  memset(cI, 0, NUM_ORIENTATIONS_INSTR * sizeof(int));
  memset(fI, 0, NUM_ORIENTATIONS_INSTR* sizeof(int));

  for(i = 0; i < GetNumVA_CP_Index(cpi); i++)
    {
      CP_Index * cpip = GetVA_CP_Index(cpi, i);

      AppendOrientedCP(graph, msp->happy, hI, cpip, cps);
      AppendOrientedCP(graph, msp->misseparatedClose, cI, cpip, cps);
      AppendOrientedCP(graph, msp->misseparatedFar, fI, cpip, cps);
      for(j = 0; j < NUM_ORIENTATIONS_INSTR; j++)
        AppendMisorientedCP(graph, msp->misoriented[j], oI[j], cpip, cps);
    }
  return 0;
}


/********************************************************************
                     Functions for instrumenting
********************************************************************/
/*
  Initial version ignores fragments in surrogates;
  as of Feb, 2004, ALH trying to take surrogates into account ...
*/
int AddFragmentToUnitigInstrumenter(ScaffoldGraphT * graph,
                                    MultiAlignT * uma,
                                    CDS_CID_t fi,
                                    UnitigInstrumenter * ui)
{
  IntMultiPos * imp = GetIntMultiPos(uma->f_list, fi);
  InfoByIID * info = GetInfoByIID(graph->iidToFragIndex, imp->ident);
  CIFragT * frag = GetCIFragT(graph->CIFrags, info->fragIndex);

#ifdef DEBUG2
  fprintf(stderr, "Adding fragment " F_CID " (index = " F_CID ") to unitig instrumenter\n",
          frag->iid, info->fragIndex);
#endif

  switch(frag->type)
    {
      case AS_READ:
      case AS_EXTR:
      case AS_TRNR:
        {
          ui->numReads +=
            (frag->type == AS_READ ||
             frag->type == AS_TRNR) ? 1 : 0;
          ui->numExtReads += (frag->type == AS_EXTR) ? 1 : 0;

          if(frag->mateOf != NULLINDEX)
            {
              if(!ExistsInHashTable_AS(ui->bookkeeping.fragHT, (uint64)&(frag->iid), 0))
                {
                  if(InsertInHashTable_AS(ui->bookkeeping.fragHT,
                                          (uint64)frag->iid, 0,
                                          (uint64)frag, 0) == HASH_FAILURE)
                    {
                      fprintf(stderr, "Failed to insert frag into hashtable.\n");
                      return 1;
                    }
                  AppendVA_CDS_CID_t(ui->bookkeeping.fragArray, &(frag->iid));
                }
            }
          else
            {
              FragDetail fragDetail;
              // no mate, log it
              fragDetail.iid = frag->iid;
              fragDetail.type = frag->type;
              fragDetail.offset5p = frag->contigOffset5p.mean;
              AppendVA_FragDetail(ui->mates.noMate, &fragDetail);
            }
        }
        break;
      default:
        fprintf(stderr, "Unknown fragment type %c encountered\n",
                (char) frag->type);
        break;
    }
  return 0;
}


int AddFragmentToSurrogateTracker(ScaffoldGraphT * graph,
				  HashTable_AS *cpHT,
                                  CDS_CID_t contigID,
                                  IntMultiPos * imp,
                                  float aEnd,
                                  float bEnd,
                                  SurrogateTracker * st)
{
  // Don't bother, if we've run out of room in the array
  if(st->numAllocatedLocs == st->numUsedLocs)
    {
      fprintf(stderr, "Ran out of space for tracking surrogate fragments\n");
      return 1;
    }

  // the goal is to compute the fragment's coordinates relative to the
  // contig, and store them appropriately in the SurrogateTracker;
  // so, we want aEnd,bEnd to be the ends of the unitig relative to 
  // the contig ...
  {
    SurrogateFragLocation * sflp;
    int addToHashTable = 1;
    InfoByIID * info = GetInfoByIID(graph->iidToFragIndex, imp->ident);
    CIFragT * frag = GetCIFragT(graph->CIFrags, info->fragIndex);

    /*
      Look up the IID in the surrogate ht
      if present, follow linked list until next is NULL
      add an entry to the array & change the NULL to point to it
      if not present, add to hashtable
    */
    if((sflp = (SurrogateFragLocation *)LookupValueInHashTable_AS(st->surrogateFragHT,
                                                                  (uint64)frag->iid, 0)))
      {
        // found entry for fragment. follow linked list to the last one
        while(sflp->nextIndex != 0)
          {
            sflp = &(st->surrogateFragLocs[sflp->nextIndex]);
          }
        sflp->nextIndex = st->numUsedLocs;
        addToHashTable = 0;
      }

    /* aEnd is A end coordinate of unitig in contig
       bEnd is B end coordinate of unitig in contig
       imp has coordinates of fragment in unitig
    */
    // populate the new surrogate fragment entry in the array
    st->surrogateFragLocs[st->numUsedLocs].contig = contigID;
#if 0
    st->surrogateFragLocs[st->numUsedLocs].offset5p =
      imp->position.bgn + (aEnd < bEnd) ? aEnd : bEnd;
    st->surrogateFragLocs[st->numUsedLocs].offset3p =
      imp->position.end + (aEnd < bEnd) ? aEnd : bEnd;
    st->surrogateFragLocs[st->numUsedLocs].nextIndex = 0;
#else
    if(aEnd<bEnd){
      st->surrogateFragLocs[st->numUsedLocs].offset5p =
	aEnd + imp->position.bgn;
      st->surrogateFragLocs[st->numUsedLocs].offset3p =
	aEnd + imp->position.end;
    }else{
      st->surrogateFragLocs[st->numUsedLocs].offset5p =
	aEnd - imp->position.bgn;
      st->surrogateFragLocs[st->numUsedLocs].offset3p =
	aEnd - imp->position.end;
    }
#endif

    // if the fragment wasn't in the hashtable, add it
    if(addToHashTable)
      {
        if(InsertInHashTable_AS(st->surrogateFragHT,
                                (uint64)frag->iid, 0,
                                (uint64)&st->surrogateFragLocs[st->numUsedLocs], 0)
           != HASH_SUCCESS)
          {
            fprintf(stderr,
                    "Failed to insert surrogate fragment into hashtable.\n");
            return 1;
          }
      }
    st->numUsedLocs++;
  }
  return 0;
}


void CheckMateLinkStatus(unsigned int innieMates,
                         DistT * dptr,
                         float frag5p,
                         FragOrient fragOrient,
                         float mate5p,
                         FragOrient mateOrient,
                         InstrumentOrientations * orientShouldBe,
                         InstrumentOrientations * orientIs,
                         InstrumentDistStatus * distStatus)
{
  // determine what the orientation should be
  *orientShouldBe = (innieMates) ? INNIE_INSTR : OUTTIE_INSTR;

  // determine what the orientation is
  if(fragOrient == mateOrient)
    *orientIs = (fragOrient == A_B) ? NORMAL_INSTR : ANTINORMAL_INSTR;
  else
    {
      if((fragOrient == A_B && frag5p < mate5p) ||
         (mateOrient == A_B && mate5p < frag5p))
        *orientIs = INNIE_INSTR;
      else
        *orientIs = OUTTIE_INSTR;
    }

  // if orientations agree, determine distance
  if(*orientShouldBe == *orientIs)
    {
      float dist = fabs(mate5p - frag5p);

      if(dist > dptr->mu - INSTRUMENT_CUTOFF * dptr->sigma)
        {
          if(dist < INTERVAL(dptr))
            *distStatus = DISTANCE_OKAY_INSTR;
          else
            *distStatus = DISTANCE_TOO_FAR_INSTR;
        }
      else
        *distStatus = DISTANCE_TOO_CLOSE_INSTR;
    }
  // else, don't bother with *distStatus
}


int GetFragmentPositionInFauxScaffold(HashTable_AS * cpHT,
                                      CIFragT * frag,
                                      CDS_COORD_t * fragLeftEnd,
                                      CDS_COORD_t * fragRightEnd,
                                      int * fragOrientInScaffold)
{
  ContigPlacement * cp;

  cp = (ContigPlacement *)LookupValueInHashTable_AS(cpHT, (uint64)frag->contigID, 0);
  if(cp == NULL)
    {
      fprintf(stderr, "Fragment " F_CID "'s contig " F_CID " is not in hashtable!\n",
              frag->iid, frag->contigID);
      return 1;
    }

  // 0 means A_B, non-0 means B_A
  if(cp->orient == A_B)
    {
      *fragLeftEnd = cp->offset +
        MIN(frag->contigOffset5p.mean, frag->contigOffset3p.mean);
      *fragRightEnd = cp->offset +
        MAX(frag->contigOffset5p.mean, frag->contigOffset3p.mean);
      *fragOrientInScaffold =
        (frag->contigOffset5p.mean < frag->contigOffset3p.mean) ? 0: 1;
    }
  else
    {
      *fragLeftEnd = cp->offset + cp->length -
        MAX(frag->contigOffset5p.mean, frag->contigOffset3p.mean);
      *fragRightEnd = cp->offset + cp->length -
        MIN(frag->contigOffset5p.mean, frag->contigOffset3p.mean);
      *fragOrientInScaffold =
        (frag->contigOffset5p.mean > frag->contigOffset3p.mean) ? 0: 1;
    }
  return 0;
}


void GetFragmentPosition(HashTable_AS * cpHT,
                         CIFragT * frag,
                         CDS_COORD_t * frag5p,
                         CDS_COORD_t * frag3p,
                         FragOrient * fragOrient)
{
  if(cpHT == NULL)
    {
      // get fragment position & orientation in contig
      *frag5p = (CDS_COORD_t) frag->contigOffset5p.mean;
      *frag3p = (CDS_COORD_t) frag->contigOffset3p.mean;
      *fragOrient = ((frag->contigOffset5p.mean < frag->contigOffset3p.mean) ?
                     A_B : B_A);
    }
  else
    {
      CDS_COORD_t fragLeftEnd;
      CDS_COORD_t fragRightEnd;
      int fragOrientInScaffold;
      // get fragment position & orientation in scaffold
      GetFragmentPositionInFauxScaffold(cpHT,
                                        frag,
                                        &fragLeftEnd,
                                        &fragRightEnd,
                                        &fragOrientInScaffold);
      *fragOrient = (fragOrientInScaffold == 0) ? A_B : B_A;
      *frag5p = (*fragOrient == A_B) ? fragLeftEnd : fragRightEnd;
      *frag3p = (*fragOrient == B_A) ? fragLeftEnd : fragRightEnd;
    }
}




int GetFragment5pPositionInFauxScaffoldGivenCtgPsn(HashTable_AS * cpHT,
						   int32 contigIID,
						   CDS_COORD_t ctg5p,
						   CDS_COORD_t *scf5p)
{
  ContigPlacement * cp;
  static int firstTime=1;

  cp = (ContigPlacement *)LookupValueInHashTable_AS(cpHT, (uint64)contigIID, 0);
  if(cp == NULL)
    {
      fprintf(stderr, "Contig %u is not in hashtable!\n",
              contigIID);
      return 1;
    }

  if(cp->orient == A_B)
    {
      *scf5p = (CDS_COORD_t) (cp->offset + ctg5p);
    }
  else
    {
      *scf5p = (CDS_COORD_t) (cp->offset + cp->length - ctg5p);
    }


  return 0;
}


void GetFragment5pPositionGivenCtgPsn(HashTable_AS * cpHT,
				      int32 contigIID,
				      CDS_COORD_t ctg5p,
				      CDS_COORD_t *newPsn)
{
  if(cpHT == NULL)
    {
      // get fragment position & orientation in contig
      *newPsn = ctg5p;
    }
  else
    {
      // get fragment position & orientation in scaffold
      GetFragment5pPositionInFauxScaffoldGivenCtgPsn(cpHT,contigIID,
                                                     ctg5p,
                                                     newPsn);
    }
}




int GetSurrogatePositionInFauxScaffoldFromSFL(HashTable_AS * cpHT,
                                              int32 contigID,
                                              CDS_COORD_t *frag5p,
                                              CDS_COORD_t *frag3p)
{
  ContigPlacement * cp;

  // on entry, *frag5p should be set to the position on the contig
  // of the 5p end of the fragment; likewise for *frag3p
  CDS_COORD_t AEndOnCtg = *frag5p;
  CDS_COORD_t BEndOnCtg = *frag3p;

  cp = (ContigPlacement *)LookupValueInHashTable_AS(cpHT, (uint64)contigID, 0);
  if(cp == NULL)
    {
      fprintf(stderr, "A surrogate fragment's contig %u is not in hashtable!\n",
              contigID);
      return 1;
    }


  // 0 means A_B, non-0 means B_A
  if(cp->orient == A_B)
    {
      *frag5p = (CDS_COORD_t) cp->offset + AEndOnCtg; 
      *frag3p = (CDS_COORD_t) cp->offset + BEndOnCtg; 
    }
  else
    {
      *frag5p = (CDS_COORD_t) (cp->offset + cp->length) - AEndOnCtg;
      *frag3p = (CDS_COORD_t) (cp->offset + cp->length) - BEndOnCtg;
    }
#if 0
  fprintf(stderr,"Looking for surrogate position: from [" F_COORD "," F_COORD "] to [" F_COORD "," F_COORD "] using offset %g, orient %s , length %g\n",
	  AEndOnCtg,BEndOnCtg,*frag5p,*frag3p,
	  cp->offset, cp->orient == A_B ? "A_B" : "B_A",cp->length);
#endif
  return 0;
}


void GetSurrogatePositionFromSFL(HashTable_AS * cpHT,
				 SurrogateFragLocation *sflp,
				 CDS_COORD_t *frag5p,
				 CDS_COORD_t *frag3p){
  *frag5p = (CDS_COORD_t) sflp->offset5p;
  *frag3p = (CDS_COORD_t) sflp->offset3p;
  if( cpHT != NULL ){
    // get fragment position in scaffold
    GetSurrogatePositionInFauxScaffoldFromSFL(cpHT,
                                              sflp->contig,
                                              frag5p,
                                              frag3p);
  }
#if 0
  fprintf(stderr,"Looking for surrogate position: from [" F_COORD "," F_COORD "] to [" F_COORD "," F_COORD "] using cpHT %x\n",
	  sflp->offset5p,sflp->offset3p,*frag5p,*frag3p,cpHT);
#endif
}


void PrintScaffoldMateDetail(HashTable_AS * cpHT,
                             MateDetail * md,
                             ChunkOrientationType oShouldBe,
                             ChunkOrientationType oIs,
                             char * category,
                             CDS_CID_t id,
                             FILE * printTo,
			     int printType)
{
  CDS_CID_t fragIID;
  CDS_CID_t mateIID;
  CDS_CID_t fragChunkIID;
  CDS_CID_t mateChunkIID;
  CDS_COORD_t frag5p,frag3p;
  CDS_COORD_t mate5p,mate3p;
  FragOrient fragOrient;
  FragOrient mateOrient;
  char oString[1024];
  CIFragT *tmpfrg;

  tmpfrg = GetCIFragT(ScaffoldGraph->CIFrags,
		      GetInfoByIID(ScaffoldGraph->iidToFragIndex,
				   md->fragIID)->fragIndex);

  
  {
    DistT *dst=GetDistT(ScaffoldGraph->Dists,tmpfrg->dist);
    assert(dst!=NULL);
    if(! USE_ALL_MATES && ! (USE_LONG_MATES && dst->mu > 15000))
      return;
  }

  if(!ExistsInHashTable_AS(cpHT,(uint64)tmpfrg->contigID,0)){
    // this should apply only to surrogate fragments

    /* intra-ctg mate pair; coordinates are relative to the contig */
    /* inter-ctg mate pair; coordinates are relative to the scaffold */
    
    //    GetFragment5pPositionGivenCtgPsn(cpHT,
    GetFragment5pPositionGivenCtgPsn( (md->fragChunkIID==md->mateChunkIID) ? cpHT : NULL,
				      md->fragChunkIID,
				      md->fragOffset5p,
				      &frag5p);
  } else {
    GetFragmentPosition(cpHT,
			tmpfrg,
			&frag5p,
			&frag3p,
			&fragOrient);
  }

  tmpfrg = GetCIFragT(ScaffoldGraph->CIFrags,
		      GetInfoByIID(ScaffoldGraph->iidToFragIndex,
				   md->mateIID)->fragIndex);
  if(!ExistsInHashTable_AS(cpHT,(uint64)tmpfrg->contigID,0)){
    // this should apply only to surrogate fragments

    /* intra-ctg mate pair; coordinates are relative to the contig */
    /* inter-ctg mate pair; coordinates are relative to the scaffold */

    //    GetFragment5pPositionGivenCtgPsn(cpHT,
    GetFragment5pPositionGivenCtgPsn( (md->fragChunkIID==md->mateChunkIID) ? cpHT : NULL,
				      md->mateChunkIID,
				      md->mateOffset5p,
				      &mate5p);
  } else {
    GetFragmentPosition(cpHT,
			tmpfrg,
			&mate5p,
			&mate3p,
			&mateOrient);
  }

  if(mate5p < frag5p)
    {
      CDS_COORD_t temp5p;
      fragIID = md->mateIID;
      mateIID = md->fragIID;
      fragChunkIID = md->mateChunkIID;
      mateChunkIID = md->fragChunkIID;
      temp5p = frag5p;
      frag5p = mate5p;
      mate5p = temp5p;
      oShouldBe = FlipEdgeOrient(oShouldBe);
      oIs = FlipEdgeOrient(oIs);
    }
  else
    {
      fragIID = md->fragIID;
      mateIID = md->mateIID;
      fragChunkIID = md->fragChunkIID;
      mateChunkIID = md->mateChunkIID;
    }

#if 0
  fprintf(stderr,"pair (" F_CID "," F_CID ") assigned 5p positions " F_COORD " " F_COORD " based on offsets %g %g on ctgs " F_CID " " F_CID "\n",
	  fragIID,mateIID,frag5p,mate5p,
	  ( md->fragChunkIID == fragChunkIID) ? md->fragOffset5p : md->mateOffset5p,	  
	  ( md->fragChunkIID == fragChunkIID) ? md->mateOffset5p : md->fragOffset5p,	  
	  fragChunkIID,mateChunkIID);
#endif

  if(oShouldBe == oIs)
    {
      sprintf(oString, "%s",
              oIs == AB_AB ? "AB_AB" :
              (oIs == AB_BA ? "AB_BA" :
               (oIs == BA_AB ? "BA_AB" : "BA_BA")));
    }
  else
    {
      sprintf(oString, "%s:%s",
              oShouldBe == AB_AB ? "AB_AB" :
              (oShouldBe == AB_BA ? "AB_BA" :
               (oShouldBe == BA_AB ? "BA_AB" : "BA_BA")),
              oIs == AB_AB ? "AB_AB" :
              (oIs == AB_BA ? "AB_BA" :
               (oIs == BA_AB ? "BA_AB" : "BA_BA")));
    }
  if(printType == PRINTTABLE ){
    fprintf(printTo,
	    "%s\t%s\t" F_CID "\t" F_CID "\t%c\t"
	    "" F_CID "\t" F_COORD "\t" F_CID "\t"
	    "" F_CID "\t" F_COORD "\t" F_CID "\n",
	    category, oString, id, md->libIID, md->type,
	    fragIID, frag5p, fragChunkIID,
	    mateIID, mate5p, mateChunkIID);
  } else {

    char markString[300] = "";
    char catString[50];
    int row;

#define MARK_INTER_SCAFFOLD_MATES
#ifdef MARK_INTER_SCAFFOLD_MATES
    int frgScf,mateScf;
    InfoByIID * info;
    CIFragT * frag;
    int32 scfId;
#endif
    
    // new version here
    if(strcmp(category,"SATISFIED")==0){
      strcpy(catString,"A1CMColor");
      row=110+(int)( abs(frag5p-mate5p)/1000.);
    } else if(strcmp(category,"MISORIENTED")==0){
      if(oIs == BA_BA){
	strcpy(catString,"A2CMColor");
	row=70;
      } else if(oIs ==AB_AB) {
	strcpy(catString,"A5CMColor");
	row=80;
      } else {
	strcpy(catString,"A6CMColor");
	row=90;
      }
    } else if(strcmp(category,"TOO_CLOSE")==0){ 
      strcpy(catString,"A3CMColor");
      row=100;
    } else if(strcmp(category,"TOO_FAR")==0){
      strcpy(catString,"A4CMColor");
      row=100;
    }

#ifdef MARK_INTER_SCAFFOLD_MATES
    frgScf = GetGraphNode(ScaffoldGraph->ContigGraph,
			  fragChunkIID)->scaffoldID;
    mateScf = GetGraphNode(ScaffoldGraph->ContigGraph,
			   mateChunkIID)->scaffoldID;
    if(frgScf==mateScf){
      markString[0]='\0';
    }else{
      sprintf(markString,"%s " F_COORD " A0InterScfColor " F_COORD " ",catString,frag5p+(mate5p-frag5p)/3,frag5p+2*(mate5p-frag5p)/3);
    }
#endif

    if(printMateUIDs){
      uint64 fragUID;
      uint64 mateUID;
      GateKeeperFragmentRecord gkpFrag;

      if(getGateKeeperFragment(ScaffoldGraph->gkpStore,fragIID,&gkpFrag)!=0){
	assert(0);
      }
      fragUID=gkpFrag.readUID;
      if(getGateKeeperFragment(ScaffoldGraph->gkpStore,mateIID,&gkpFrag)!=0){
	assert(0);
      }
      mateUID=gkpFrag.readUID;
      fprintf(printTo,F_UID "Mate" F_UID ": " F_COORD " %s%s " F_COORD " R%d # %s " F_CID " " F_CID "\n",
	      fragUID,mateUID,frag5p,markString,catString,mate5p,row,category,fragIID,mateIID);
    } else {
      fprintf(printTo,F_CID "Mate" F_CID ": " F_COORD " %s%s " F_COORD " R%d # %s " F_CID " " F_CID "\n",
	      fragIID,mateIID,frag5p,markString,catString,mate5p,row,category,fragIID,mateIID);
    }


  }
}


void PrintScaffoldMateDetailArray(HashTable_AS * cpHT,
                                  VA_TYPE(MateDetail) * mda,
                                  ChunkOrientationType oShouldBe,
                                  ChunkOrientationType oIs,
                                  char * category,
                                  CDS_CID_t id,
                                  FILE * printTo,
				  int printType)
{
  int32 i;
  for(i = 0; i < GetNumVA_MateDetail(mda); i++)
    PrintScaffoldMateDetail(cpHT,
                            GetVA_MateDetail(mda, i),
                            oShouldBe, oIs, category, id, printTo,printType);
}


void PrintScaffoldInstrumenterMateDetails(ScaffoldInstrumenter * si,
                                          FILE * printTo,
					  int printType)
{
  int ori1, ori2;
  MateStatusPositions * msp = si->mates.mateStatus->intra;
  // VA_TYPE(MateDetail) * wExtMates = si->bookkeeping.wExtMates;

  while(msp != NULL)
    {
      for(ori1 = 0; ori1 < NUM_ORIENTATIONS_INSTR; ori1++)
        {
          ChunkOrientationType oShouldBe, oIs;
          switch(ori1)
            {
              case INNIE_INSTR:
                oShouldBe = AB_BA;
                break;
              case OUTTIE_INSTR:
                oShouldBe = BA_AB;
                break;
              case NORMAL_INSTR:
                oShouldBe = AB_AB;
                break;
              case ANTINORMAL_INSTR:
                oShouldBe = BA_BA;
                break;
              default:
                oShouldBe = XX_XX;
                break;
            }
      
          oIs = oShouldBe;
          /*
            PrintScaffoldMateDetailArray(si->cpHT, wExtMates, oShouldBe, oIs,
            (doingScaffold) ? "INTER-SCAFFOLD" : "INTER-CONTIG",
            si->id, printTo,printType);
          */
      
          PrintScaffoldMateDetailArray(si->cpHT, msp->happy[ori1], oShouldBe, oIs,
                                       "SATISFIED", si->id, printTo,printType);
      
          for(ori2 = 0; ori2 < NUM_ORIENTATIONS_INSTR; ori2++)
            {
              switch(ori2)
                {
                  case INNIE_INSTR:
                    oIs = AB_BA;
                    break;
                  case OUTTIE_INSTR:
                    oIs = BA_AB;
                    break;
                  case NORMAL_INSTR:
                    oIs = AB_AB;
                    break;
                  case ANTINORMAL_INSTR:
                    oIs = BA_BA;
                    break;
                  default:
                    oIs = XX_XX;
                    break;
                }
              PrintScaffoldMateDetailArray(si->cpHT,
                                           msp->misoriented[ori1][ori2], oShouldBe, oIs,
                                           "MISORIENTED", si->id, printTo,printType);
            }
          oIs = oShouldBe;
          PrintScaffoldMateDetailArray(si->cpHT,
                                       msp->misseparatedClose[ori1],
                                       oShouldBe, oIs,
                                       "TOO_CLOSE", si->id, printTo,printType);
          PrintScaffoldMateDetailArray(si->cpHT,
                                       msp->misseparatedFar[ori1],
                                       oShouldBe, oIs,
                                       "TOO_FAR", si->id, printTo,printType);
        }
      if(msp == si->mates.mateStatus->intra)
        msp = si->mates.mateStatus->inter;
      else
        msp = NULL;
    }
}


void PrintUnmatedDetails(ScaffoldInstrumenter * si,
			 FILE * printTo,
			 int printType)
{
  assert(printType==PRINTCELAMY); // nothing else coded yet ..
  assert(si->cpHT != NULL);
  {//  stuff based on CheckFragmentMatePairs() ...

    int i;
    for(i = 0; i < GetNumVA_FragDetail(si->mates.noMate); i++){
      CIFragT * frag;
      FragDetail *fd;
      int32 fragIID;
      CDS_COORD_t fragLeftEnd, fragRightEnd,frag5p, frag3p;
      FragOrient fragOrient;
      
      fd=GetVA_FragDetail(si->mates.noMate,i);
      fragIID = fd->iid;
      //      fprintf(stderr,"Got frag detail %x iid %d\n",fd,fragIID);
      frag = getFragByIID(ScaffoldGraph,fragIID);
      //      fprintf(stderr,"Got frag %x\n",frag);


      GetFragmentPosition(si->cpHT,frag,&frag5p,&frag3p,&fragOrient);
      //      fprintf(stderr,"Got frag positions\n",frag);
      
      if(frag5p<frag3p){
	fragLeftEnd=frag5p;
	fragRightEnd=frag3p;
      }else{
	fragLeftEnd=frag3p;
	fragRightEnd=frag5p;
      }

      fprintf(printTo,F_CID "Fragment: " F_COORD " A9CMColor " F_COORD " R45 # Unmated Fragment " F_CID " ori:%s\n",
	      frag->iid,fragLeftEnd,fragRightEnd,frag->iid,
	      fragOrient==0 ? "A_B" : "B_A");
    }
  }
  return;
}


int CheckFragmentMatePairs(ScaffoldGraphT * graph,
                           InstrumenterBookkeeping * bookkeeping,
                           SurrogateTracker * st,
                           HashTable_AS * cpHT,
                           HashTable_AS * anchoredHT,
                           MateStatusPositionsSet * msps,
                           uint32 options,
                           int32 * numMiso,
                           int32 * numFar,
                           int32 * numClose,
                           CDS_CID_t chunkIID)
{
  int32 i;
  int doingContig = (cpHT == NULL) ? 1 : 0;
  MateStatusPositions * msp = (doingContig) ? msps->intra : msps->inter;

  //  fprintf(stderr,"CheckFragmentMatePairs applied to %d cases\n",
  //	  GetNumVA_int32(bookkeeping->fragArray));

  for(i = 0; i < GetNumVA_CDS_CID_t(bookkeeping->fragArray); i++)
    {
      CIFragT * frag;
      CIFragT * graphMate;
      CIFragT * lookupMate;
      CIFragT   mockMate;
      InstrumentOrientations orientShouldBe;
      InstrumentOrientations orientIs;
      InstrumentDistStatus distStatus;
      CDS_COORD_t frag5p,frag3p;
      FragOrient fragOrient;
      CDS_COORD_t mate5p,mate3p;
      FragOrient mateOrient;

      // get current fragment & its mate
      frag = getBookkeepingFrag(graph, bookkeeping->fragArray, i);
      GetFragmentPosition(cpHT, frag, &frag5p, &frag3p, &fragOrient);
    
      /*    fprintf(stderr,"checking mates of %d CIid %d doingContig %d\n",
            frag->iid,frag->CIid,doingContig);*/

      graphMate = GetCIFragT(graph->CIFrags, frag->mateOf);

      // see if the mate is in a unitig
      if((lookupMate = (CIFragT *)LookupValueInHashTable_AS(bookkeeping->fragHT, (uint64)graphMate->iid, 0)) == NULL)
        {
          // if here, mate is not in a non-surrogate unitig in this node
          // see if it's in the set of surrogate fragments
          SurrogateFragLocation * sflp;
          if((sflp = (SurrogateFragLocation *)LookupValueInHashTable_AS(st->surrogateFragHT,
                                                                   (uint64)graphMate->iid, 0)) != NULL
             && (!doingContig /* i.e. working on scf */ || 
                 sflp->contig == chunkIID /* surrogate is in same contig */)
             )
            {
              // found mate in a surrogate, make a mock fragment with usable coords
              if((doingContig && (options & INST_OPT_INTRA_MATES)) ||
                 (!doingContig && (options & INST_OPT_INTER_MATES))){
                lookupMate = &mockMate;
                //	  fprintf(stderr,"found mate in a surrogate, make a mock fragment with usable coords\n");
        
                mockMate.iid = graphMate->iid;
                do
                  {
                    // see if it's a good pair: coords should be in appropriate reference
#if 0
                    mate5p = sflp->offset5p; // this is relative to a contig!
                    mateOrient = (sflp->offset5p < sflp->offset3p) ? A_B : B_A;
                    mockMate.CIid = sflp->contig;
#else
                    // get position; if analyzing contig, this is relative to
                    // the contig; if analyzing scaffold, this is relative to
                    // the scaffold.
                    GetSurrogatePositionFromSFL(cpHT,sflp,&mate5p,&mate3p);
                    mateOrient = ( mate5p < mate3p ) ? A_B : B_A;
                    mockMate.contigID = sflp->contig;
#if 0
                    fprintf(stderr,"Found a surrogate mate (%d, frag = %d), was ctg %d [%d,%d], mapped to [%d,%d], orientation %c\n",
                            mockMate.iid,frag->iid,sflp->contig,sflp->offset5p,sflp->offset3p,mate5p,mate3p,mateOrient);
                    fprintf(stderr,"  paired to %d at 5p %d, orient %c contig %d doingContig %d\n",
                            frag->iid,frag5p,fragOrient,frag->CIid,doingContig);
#endif
#endif
                    CheckMateLinkStatus(frag->flags.bits.innieMate,
                                        GetDistT(graph->Dists,frag->dist),
                                        frag5p,
                                        fragOrient,
                                        mate5p,
                                        mateOrient,
                                        &orientShouldBe,
                                        &orientIs,
                                        &distStatus);
#if 0
                    fprintf(stderr,"   CheckMateLinkStatus gives orientIs %d\n",
                            orientIs);
#endif


                    if(orientShouldBe == orientIs &&
                       distStatus == DISTANCE_OKAY_INSTR)
                      {
                        break; // preference is given to good ones
                      }
                    sflp = &(st->surrogateFragLocs[sflp->nextIndex]);
                  } while(sflp->nextIndex != 0);
                /* here, we've found the best or the last mate position/orientation
                   NOTE: The problem is, there may be a 'good' link to an
                   instance of the mate fragment in a surrogate in another
                   scaffold
                */
              }else{
                lookupMate = NULL;
                //	  fprintf(stderr,"found mate in a surrogate, but don't make a mock fragment with usable coords\n");
              }

            }
          else
            {
              // if here, the mate isn't present, even in an surrogate
              MateDetail md;
              md.fragIID = frag->iid;
              md.fragOffset5p = frag5p;
              md.fragChunkIID = chunkIID;
              md.libIID = frag->dist;
              md.type = frag->type;
              md.mateIID = graphMate->iid;
              md.mateOffset5p = -1.f;
#ifdef TRACK_3P
              md.fragOffset3p = frag3p;
              md.mateOffset3p = -1.f;
#endif
              md.mateChunkIID = getFragByIID(ScaffoldGraph,md.mateIID)->cid;
              AppendVA_MateDetail(bookkeeping->wExtMates, &(md));
              lookupMate = NULL;
            }
        }
      else if((doingContig && (options & INST_OPT_INTRA_MATES)) ||
              (!doingContig && (options & INST_OPT_INTER_MATES)))
        {
          if(graphMate->iid > frag->iid)
            continue;

          GetFragmentPosition(cpHT, graphMate, &mate5p, &mate3p, &mateOrient);
      
          // check orientation, separation...
          CheckMateLinkStatus(frag->flags.bits.innieMate,
                              GetDistT(graph->Dists,frag->dist),
                              frag5p,
                              fragOrient,
                              mate5p,
                              mateOrient,
                              &orientShouldBe,
                              &orientIs,
                              &distStatus);
        }
      else
        {
          // mate is there, but we don't care
          lookupMate = NULL;
        }

      // if lookupMate == NULL, then mate not found even in surrogate
      if(lookupMate != NULL)
        {
          MateDetail matePair;
          // populate with fragOffset5p <= mateOffset5p
          // for subsequent intra-contig breakpoint detection
          if(frag5p < mate5p)
            {
              matePair.fragOffset5p = frag5p;
              matePair.mateOffset5p = mate5p;
#ifdef TRACK_3P
              matePair.fragOffset3p = frag3p;
              matePair.mateOffset3p = mate3p;
#endif
              matePair.fragIID = frag->iid;
              matePair.mateIID = lookupMate->iid;
              matePair.fragChunkIID = frag->contigID;
              matePair.mateChunkIID = lookupMate->contigID;
            }
          else
            {
              matePair.fragOffset5p = mate5p;
              matePair.mateOffset5p = frag5p;
#ifdef TRACK_3P
              matePair.fragOffset3p = mate3p;
              matePair.mateOffset3p = frag3p;
#endif
              matePair.fragIID = lookupMate->iid;
              matePair.mateIID = frag->iid;
              matePair.fragChunkIID = lookupMate->contigID;
              matePair.mateChunkIID = frag->contigID;
            }
          matePair.libIID = frag->dist;
          matePair.type = frag->type;
      
          if(anchoredHT && frag->contigID != lookupMate->contigID)
            {
              if(!ExistsInHashTable_AS(anchoredHT, (uint64)frag->contigID, 0))
                InsertInHashTable_AS(anchoredHT,
                                     (uint64)frag->contigID, 0,
                                     (uint64)frag->contigID, 0);
              if(!ExistsInHashTable_AS(anchoredHT, (uint64)lookupMate->contigID, 0))
                InsertInHashTable_AS(anchoredHT,
                                     (uint64)lookupMate->contigID, 0,
                                     (uint64)lookupMate->contigID, 0);
            }
         
          if(orientShouldBe == orientIs)
            {
              switch(distStatus)
                {
                  case DISTANCE_OKAY_INSTR:
                    AppendVA_MateDetail(msp->happy[orientIs], &matePair);
                    break;
                  case DISTANCE_TOO_CLOSE_INSTR:
                    *numClose++;
                    AppendVA_MateDetail(msp->misseparatedClose[orientIs], &matePair);
                    break;
                  case DISTANCE_TOO_FAR_INSTR:
                    *numFar++;
                    AppendVA_MateDetail(msp->misseparatedFar[orientIs], &matePair);
                    break;
                }
            }
          else
            {
              *numMiso++;
              AppendVA_MateDetail(msp->misoriented[orientShouldBe][orientIs],
                                  &matePair);
            }
        } // if(lookupMate != NULL) - mate was found in unitig or surrogate
    } // loop over all fragments in bookkeeping

  return 0;
}


int InstrumentUnitig(ScaffoldGraphT * graph,
		     HashTable_AS *cpHT,
                     SurrogateTracker * st,
                     ChunkInstanceT * unitig,
                     float contigAEnd,
                     float contigBEnd,
                     UnitigInstrumenter * ui)
{
  MultiAlignT * uma;
  CDS_CID_t fi;
  int32 numMiso;
  int32 numFar;
  int32 numClose;
  uint32 ctgID;

#ifdef DEBUG
  fprintf(stderr, 
	  "\t\tInstrumenting unitig " F_CID
	  " [" F_COORD "," F_COORD "] relative to contig [" F_COORD "," F_COORD "]\n",
	  unitig->id,
	  unitig->offsetAEnd.mean,
	  unitig->offsetBEnd.mean,
	  contigAEnd,contigBEnd);
#endif

  InitializeUnitigInstrumenter(graph, ui);
  
  // get the position, whether or not it's a surrogate
  ui->leftEnd = MIN(unitig->offsetAEnd.mean, unitig->offsetBEnd.mean) + 0.5f;
  ui->rightEnd = MAX(unitig->offsetAEnd.mean, unitig->offsetBEnd.mean) + 0.5f;
  ui->orientation =
    (unitig->offsetAEnd.mean < unitig->offsetBEnd.mean) ? A_B : B_A;
  ctgID = unitig->info.CI.contigID;
  
  // surrogate? make sure we have a real unitig
  if(unitig->flags.bits.isStoneSurrogate || unitig->flags.bits.isWalkSurrogate)
    {
      ui->isSurrogate = 1;
      unitig = GetGraphNode(graph->CIGraph, unitig->info.CI.baseID);
      if(unitig == NULL)
        {
          fprintf(stderr, "Surrogate's unitig " F_CID " does not exist in the graph!\n",
                  unitig->info.CI.baseID);
          return 1;
        }
    }

  // get the multialignment - lists fragments
  uma = LoadMultiAlignTFromSequenceDB(graph->sequenceDB, unitig->id, TRUE);
  if(uma == NULL)
    {
      fprintf(stderr, "Failed to load MultiAlignT of unitig " F_CID "\n", unitig->id);
      return 1;
    }
  
  if(ui->isSurrogate)
    {
      if(do_surrogate_tracking){
        float utgAEndOnCtg,utgBEndOnCtg;
        utgAEndOnCtg = ((ui->orientation == A_B) ? ui->leftEnd : ui->rightEnd );
        utgBEndOnCtg = ((ui->orientation == B_A) ? ui->leftEnd : ui->rightEnd );

        // iterate over fragments in surrogate
        for(fi = 0; fi < GetNumIntMultiPoss(uma->f_list); fi++)
          {
            // add to surrogate set - position & orientation in contig
            AddFragmentToSurrogateTracker(graph, cpHT, ctgID,
                                          GetIntMultiPos(uma->f_list, fi),
                                          utgAEndOnCtg,utgBEndOnCtg,
                                          st);
          }
      }
    }
  else
    {
#ifdef LIST_TERMINAL_TYPES
      CDS_CID_t firstFragIID;
      FragType  firstFragType;
      CDS_CID_t lastFragIID;
      FragType  lastFragType;
      CDS_COORD_t lastBP = 0;
      IntMultiPos * imp;
#endif
      // iterate over fragments
      for(fi = 0; fi < GetNumIntMultiPoss(uma->f_list); fi++)
        {
#ifdef LIST_TERMINAL_TYPES
          imp = GetIntMultiPos(uma->f_list, fi);
          if(lastBP == 0)
            {
              firstFragIID = imp->ident;
              firstFragType = imp->type;
            }
          if(lastBP < MAX(imp->position.bgn, imp->position.end))
            {
              lastFragIID = imp->ident;
              lastFragType = imp->type;
            }
          lastBP = MAX(lastBP, MAX(imp->position.bgn, imp->position.end));
      
          if(ContigLastBP == 0)
            {
              ContigFirstFragIID = imp->ident;
              ContigFirstFragType = imp->type;
            }
          if(ContigLastBP < UnitigOffset + MAX(imp->position.bgn,
                                               imp->position.end))
            {
              ContigLastFragIID = imp->ident;
              ContigLastFragType = imp->type;
            }
          ContigLastBP = UnitigOffset + MAX(imp->position.bgn,
                                            imp->position.end);
#endif
          AddFragmentToUnitigInstrumenter(graph, uma, fi, ui);
        }
#ifdef LIST_TERMINAL_TYPES
      fprintf(stdout, "Terminal fragments for unitig " F_CID ":\t" F_CID ", %c\t" F_CID ", %c\n",
              uma->id, firstFragIID, firstFragType, lastFragIID, lastFragType);
#endif
    }

  // unload the multialignment
  UnloadMultiAlignTFromSequenceDB(graph->sequenceDB, unitig->id, TRUE);

  CheckFragmentMatePairs(graph,
                         &(ui->bookkeeping),
                         st,
                         NULL,
                         NULL,
                         ui->mates.mateStatus,
                         ui->options,
                         &numMiso,
                         &numFar,
                         &numClose,
                         ui->id);
  
  return 0;
}


/* For a contig
   loop over unitigs in contig
   - add unitig size to variable array
   if surrogate, count
   for each fragment in unitig/surrogate
   check that frag is in correct unitig & contig
   count type
   if read or bac end & has mate, add to fragHT.
   if in surrogate, add to surrogateFragHT & surrogateFragLocs
       
   after all unitigs have been processed,
   
   loop over elements 'a' in fragHT
   look for 'a' in surrogateFragHT
   get cifrag 'a' & look for mate 'b' in fragHT
   if 'a' is in surrogateFragHT:
   if 'b' is missing, no biggy
   else check distance & orientation in contig & record
   else
   if 'b' is missing,
   if 'a' is too far from contig end, record
   else add 'a' to list of external_mates
   else check distance & orientation in contig & record

   loop over unitigs & generate stats

   report #reads, #bac ends, #locales
*/
int InstrumentContig(ScaffoldGraphT * graph,
		     HashTable_AS *cpHT,
                     SurrogateTracker * st,
                     ChunkInstanceT * contig,
                     ContigInstrumenter * ci,
                     float aEnd,
                     float bEnd)
{
  ChunkInstanceT * unitig;
  ContigTIterator unitigIterator;
  int32 numMiso;
  int32 numFar;
  int32 numClose;

#ifdef DEBUG
  fprintf(stderr, "\tInstrumenting contig " F_CID "\n", contig->id);
#endif
  if(graph == NULL || ci == NULL)
    {
      fprintf(stderr, "graph or contig instrumenter is NULL!\n");
      return 1;
    }
  InitializeContigInstrumenter(graph, ci);

  ci->id = contig->id;
  ci->leftEnd = MIN(aEnd, bEnd) + 0.5f;
  ci->rightEnd = MAX(aEnd, bEnd) + 0.5f;
  ci->orientation = (aEnd < bEnd) ? A_B : B_A;

#ifdef LIST_TERMINAL_TYPES
  ContigLastBP = 0;
#endif
  
  // Iterate over unitigs in contig & add data to contig instrumenter
  InitContigTIterator(graph, contig->id, TRUE, FALSE, &unitigIterator);
  while((unitig = NextContigTIterator(&unitigIterator)) != NULL)
    {
#ifdef LIST_TERMINAL_TYPES
      UnitigOffset = MIN(unitig->offsetAEnd.mean, unitig->offsetBEnd.mean);
#endif
      InstrumentUnitig(graph, cpHT, st, unitig, aEnd, bEnd, &(ci->reusableUI));
      AddUnitigToContigInstrumenter(graph, ci, &(ci->reusableUI));
    }
#ifdef LIST_TERMINAL_TYPES
  fprintf(stdout, "Terminal fragments for contig " F_CID ":\t" F_CID ", %c\t" F_CID ", %c\n",
          contig->id,
          ContigFirstFragIID, ContigFirstFragType,
          ContigLastFragIID, ContigLastFragType);
#endif
  
  CheckFragmentMatePairs(graph,
                         &(ci->bookkeeping),
                         st,
                         NULL,
                         NULL,
                         ci->mates.mateStatus,
                         ci->options,
                         &numMiso,
                         &numFar,
                         &numClose,
                         ci->id);

  // detect intra-contig breakpoints
  /*
    if(ci->options & INST_OPT_BREAKPOINTS &&
    (numClose >= INST_MIN_BREAK_MATES ||
    numFar   >= INST_MIN_BREAK_MATES ||
    numMiso  >= INST_MIN_BREAK_MATES))
    DetectRoughIntraContigBreakpoints(graph, ci);
  */

  return 0;
}


void PopulateICP(IntContigPairs * icp, CDS_CID_t id, CIEdgeT * edge)
{
  icp->contig1 = id;
  icp->contig2 = (edge->idA == id) ? edge->idB : edge->idA;
  icp->mean = edge->distance.mean;
  if (edge->distance.variance < 0) {
      fprintf(stderr,"Negative variance in sqrt: ctg1:%ld ctg2:%ld mean:%lf var:%lf\n",
              id, icp->contig2, icp->mean, edge->distance.variance);
      //assert(0);
  }
  icp->stddev = sqrt(edge->distance.variance);

  switch(edge->orient)
    {
      case AB_AB:
        icp->orient = (edge->idA == id) ? AB_AB : BA_BA;
        break;
      case BA_BA:
        icp->orient = (edge->idA == id) ? BA_BA : AB_AB;
        break;
      default:
        icp->orient = edge->orient;
        break;
    }
}


int GetOppositeEndOfOtherCI(CIEdgeT * edge, CDS_CID_t thisID)
{
  switch(edge->orient)
    {
      case AS_NORMAL:
        return((edge->idA != thisID) ? A_END : B_END);
        break;
      case AS_ANTI:
        return((edge->idA != thisID) ? B_END : A_END);
        break;
      case AS_INNIE:
        return A_END;
        break;
      case AS_OUTTIE:
        return B_END;
        break;
      default:
        return NO_END;
        break;
    }
  return NO_END;
}


int AddICP(VA_TYPE(IntContigPairs) * icps,
           CDS_CID_t * thisID,
           int * thisEnd,
           CIEdgeT * edge)
{
  IntContigPairs icp;
  // append the edge
  PopulateICP(&icp, *thisID, edge);
  AppendVA_IntContigPairs(icps, &icp);
  *thisEnd = GetOppositeEndOfOtherCI(edge, *thisID);
  if(*thisEnd == NO_END)
    {
      // fprintf(stderr, "\n");
      DeleteVA_IntContigPairs(icps);
    }
  *thisID = (edge->idA == *thisID) ? edge->idB : edge->idA;
  // fprintf(stderr, F_CID "(%c) ", *thisID, (*thisEnd == A_END) ? 'A' : 'B');
  return 0;
}


/*
  NOTE: this needs to be implemented so that the scaffold instrumenter
  maintains a pointer to icps, so it can be reused
*/
void FinishIntScaffoldMesg(IntScaffoldMesg * isf,
                           VA_TYPE(IntContigPairs) * icps)
{
  isf->num_contig_pairs = GetNumVA_IntContigPairs(icps);
  if(isf->num_contig_pairs == 0)
    isf->contig_pairs = NULL;
  else
    {
      isf->contig_pairs = safe_malloc(isf->num_contig_pairs * sizeof(IntContigPairs));
      memcpy(isf->contig_pairs, GetVA_IntContigPairs(icps, 0),
             isf->num_contig_pairs * sizeof(IntContigPairs));
    }
}


void FreeIntScaffoldMesg(IntScaffoldMesg * isf)
{
  if(isf->contig_pairs)
    safe_free(isf->contig_pairs);
  isf->num_contig_pairs = 0;
}


int BuildFauxIntScaffoldMesgFromScaffold(ScaffoldGraphT * graph,
                                         CIScaffoldT * scaffold,
                                         ScaffoldInstrumenter * si,
                                         IntScaffoldMesg * isf)
{
  CIScaffoldTIterator CIsTemp;
  VA_TYPE(IntContigPairs) * icps = CreateVA_IntContigPairs(1000);
  assert(icps != NULL);
  
  isf->iaccession = scaffold->id;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIsTemp);
  while( NextCIScaffoldTIterator(&CIsTemp) && CIsTemp.next != NULLINDEX)
    {
      float gapSize;
      float currVariance;
      ChunkOrientationType pairwiseOrient;
      ContigT * lContig;
      ContigT * rContig;
      EdgeCGW_T * edge;
      EdgeCGW_T myEdge;
      CDS_CID_t thisID;
      int thisEnd;

      // get the left contig
      if((lContig = GetGraphNode(graph->RezGraph, CIsTemp.curr)) == NULL)
        {
          fprintf(stderr, "Left contig " F_CID " does not exist in the graph!\n",
                  CIsTemp.curr);
          return 1;
        }
      // get the right contig
      if((rContig = GetGraphNode(graph->RezGraph, CIsTemp.next)) == NULL)
        {
          fprintf(stderr, "Right contig " F_CID " does not exist in the graph!\n",
                  CIsTemp.next);
          return 1;
        }

      // capture gap information between curr & next in si
      if(lContig->offsetAEnd.mean < lContig->offsetBEnd.mean)
        {
          if(rContig->offsetAEnd.mean < rContig->offsetBEnd.mean)
            {
              pairwiseOrient = AB_AB;
              gapSize = rContig->offsetAEnd.mean - lContig->offsetBEnd.mean;
              currVariance =
                rContig->offsetAEnd.variance - lContig->offsetBEnd.variance;
            }
          else
            {
              pairwiseOrient = AB_BA;
              gapSize = rContig->offsetBEnd.mean - lContig->offsetBEnd.mean;
              currVariance =
                rContig->offsetBEnd.variance - lContig->offsetBEnd.variance;
            }
        }
      else
        {
          if(rContig->offsetAEnd.mean < rContig->offsetBEnd.mean)
            {
              pairwiseOrient = BA_AB;
              gapSize = rContig->offsetAEnd.mean - lContig->offsetAEnd.mean;
              currVariance =
                rContig->offsetAEnd.variance - lContig->offsetAEnd.variance;
            }
          else
            {
              pairwiseOrient = BA_BA;
              gapSize = rContig->offsetBEnd.mean - lContig->offsetAEnd.mean;
              currVariance =
                rContig->offsetBEnd.variance - lContig->offsetAEnd.variance;
            }
        }

      /*
        if(gapSize < 0.0
        && ((-gapSize < lContig->bpLength.mean + 0.5 && -gapSize + 0.5 > lContig->bpLength.mean)
        || (-gapSize < rContig->bpLength.mean + 0.5 && -gapSize + 0.5 > rContig->bpLength.mean))
        )
        {
        fprintf(stderr, "\n*****Found contigLength = -gapSize:\n");
        fprintf(stderr, "***** ID: " F_CID ", length %.2f, gap = %.2f, ID: " F_CID ", length %.2f\n\n",
        CIsTemp.curr, lContig->bpLength.mean, gapSize, CIsTemp.next, rContig->bpLength.mean);
        }
      */
    
      // set some temporary variables that get changed in AddICP
      thisID = CIsTemp.curr;
      thisEnd = (pairwiseOrient == AB_AB || pairwiseOrient == AB_BA) ? 2 : 1;
    
      {
        myEdge.idA = CIsTemp.curr;
        myEdge.idB = CIsTemp.next;
        myEdge.distance.mean = gapSize;
        myEdge.distance.variance = currVariance;
        myEdge.orient = pairwiseOrient;
        edge = &myEdge;
      }
      // add to set of ICPs
      AddICP(icps, &thisID, &thisEnd, edge);
    }

  // At this point, if there are no contig pairs, need to create one
  if(GetNumVA_IntContigPairs(icps) == 0)
    {
      EdgeCGW_T myEdge;
      IntContigPairs icp;
    
      myEdge.idA = myEdge.idB = CIsTemp.curr;
      myEdge.distance.mean = myEdge.distance.variance = 0.0f;
      myEdge.orient = AB_AB;
      PopulateICP(&icp, CIsTemp.curr, &myEdge);
      AppendVA_IntContigPairs(icps, &icp);
    }
  FinishIntScaffoldMesg(isf, icps);
  DeleteVA_IntContigPairs(icps);
  return 0;
}


int AddCPToHashTable(HashTable_AS * ht,
                     ContigPlacement * cp)
{
  if(InsertInHashTable_AS(ht,
                          (uint64)cp->id, 0,
                          (uint64)cp, 0) == HASH_FAILURE)
    {
      fprintf(stderr, "Failed to insert contig position into hashtable.\n");
      return 1;
    }
  return 0;
}


int AddContigPlacementToScaffoldInstrumenter(ScaffoldInstrumenter * si,
                                             ContigPlacement * cp)
{
  if(GetNumVA_ContigPlacement(si->cpArray) > 0)
    {
      ContigPlacement * originCP = GetVA_ContigPlacement(si->cpArray, 0);
      AppendVA_ContigPlacement(si->cpArray, cp);
      if(originCP != GetVA_ContigPlacement(si->cpArray, 0))
        {
          // if here, need to repopulate hashtable
          int32 i;
          ResetHashTable_AS(si->cpHT);
          for(i = 0; i < GetNumVA_ContigPlacement(si->cpArray); i++)
            {
              AddCPToHashTable(si->cpHT, GetVA_ContigPlacement(si->cpArray, i));
            }
        }
      else
        {
          AddCPToHashTable(si->cpHT,
                           GetVA_ContigPlacement(si->cpArray,
                                                 GetNumVA_ContigPlacement(si->cpArray) - 1));
        }
    }
  else
    {
      AppendVA_ContigPlacement(si->cpArray, cp);
      AddCPToHashTable(si->cpHT,
                       GetVA_ContigPlacement(si->cpArray,
                                             GetNumVA_ContigPlacement(si->cpArray) - 1));
    }
  return 0;
}


int InstrumentScaffoldNextContig(ScaffoldGraphT * graph,
                                 ScaffoldInstrumenter * si,
                                 IntScaffoldMesg * isf,
                                 InstrumenterVerbosity verbose,
                                 FILE * printTo)
{
  ContigT * contig;
  ContigPlacement cp;
  CDS_CID_t contigIndex = GetNumVA_ContigPlacement(si->cpArray);
  CDS_CID_t contigID;

  contigID = (contigIndex == 0) ? 
    isf->contig_pairs[0].contig1 :
    isf->contig_pairs[contigIndex - 1].contig2;
  
  // get the contig
  if((contig = GetGraphNode(graph->RezGraph, contigID)) == NULL)
    {
      fprintf(stderr, "Contig " F_CID " does not exist in the graph!\n", contigID);
      return 1;
    }
      
  // set the contig's position in the contigplacement array
  cp.id = contigID;
  cp.length = contig->bpLength.mean;
  if(contigIndex == 0)
    {
      cp.offset = 0.f;
      cp.orient = (isf->contig_pairs[0].orient == AB_AB ||
                   isf->contig_pairs[0].orient == AB_BA) ? A_B : B_A;
    }
  else
    {
      ContigPlacement * prevCP = GetVA_ContigPlacement(si->cpArray,
                                                       contigIndex - 1);
      cp.offset = prevCP->offset + prevCP->length +
        isf->contig_pairs[contigIndex - 1].mean;
      cp.orient = (isf->contig_pairs[contigIndex - 1].orient == AB_AB ||
                   isf->contig_pairs[contigIndex - 1].orient == BA_AB) ?
        A_B : B_A;
    }

  // append the contig placement to the array & hashtable
  // array may be realloc'd, so potentially repopulate hashtable
  AddContigPlacementToScaffoldInstrumenter(si, &cp);
  
  // instrument the contig
  if(InstrumentContig(graph, si->cpHT, &(si->surrogateTracker),
                      contig, &(si->reusableCI),
                      (cp.orient == A_B) ? cp.offset : cp.offset + cp.length,
                      (cp.orient == A_B) ? cp.offset + cp.length: cp.offset))
    {
      fprintf(stderr, "Failed to instrument contig " F_CID "\n", contig->id);
      return 1;
    }

  ComputeContigInstrumenterStats(graph, &(si->reusableCI));
  
  if(printTo != NULL && verbose >= InstrumenterVerbose4)
    {
      PrintContigInstrumenter(graph,
                              &(si->reusableCI),
                              verbose,
                              "\t",
                              printTo);
    }

  AddContigToScaffoldInstrumenter(graph, si, &(si->reusableCI));

  return 0;
}


int InstrumentScaffoldNextGapAndContig(ScaffoldGraphT * graph,
                                       ScaffoldInstrumenter * si,
                                       IntScaffoldMesg * isf,
                                       InstrumenterVerbosity verbose,
                                       FILE * printTo)
{
  CDS_CID_t pairIndex = GetNumVA_float(si->scaffoldGapSizes);
  EdgeCGW_T * edge;

  // capture the gap
  AppendVA_float(si->scaffoldGapSizes,
                       &(isf->contig_pairs[pairIndex].mean));

  // capture inferred edge stddevs
  if((edge = FindGraphEdge(graph->RezGraph,
                           isf->contig_pairs[pairIndex].contig1,
                           isf->contig_pairs[pairIndex].contig2,
                           isf->contig_pairs[pairIndex].orient)) == NULL)
    {
      AppendVA_float(si->inferredEdgeStddevs,
                           &(isf->contig_pairs[pairIndex].stddev));
    }

  InstrumentScaffoldNextContig(graph, si, isf, verbose, printTo);
  return 0;
}


/*
  for each contig,
  get the contig from the graph
  instrument it
  populate contigplacement item with id, offset, length, orient
  partly based on contigplacement of previous contig
  add to array & check pointer to first item
  if(pointer changed)
  repopulate entire hashtable
  else
  add latest contig to hashtable
  track the gap size
  populate the contig pairs array
*/
int InstrumentScaffoldInitialContigPair(ScaffoldGraphT * graph,
                                        ScaffoldInstrumenter * si,
                                        IntScaffoldMesg * isf,
                                        InstrumenterVerbosity verbose,
                                        FILE * printTo)
{
  InstrumentScaffoldNextContig(graph, si, isf, verbose, printTo);
  InstrumentScaffoldNextGapAndContig(graph, si, isf, verbose, printTo);
  return 0;
}


int InstrumentIntScaffoldMesg(ScaffoldGraphT * graph,
                              ScaffoldInstrumenter * si,
                              IntScaffoldMesg * isf,
                              InstrumenterVerbosity verbose,
                              FILE * printTo)
{
  CDS_CID_t i;
  int32 numMiso;
  int32 numFar;
  int32 numClose;
  
  InitializeScaffoldInstrumenter(graph, si);

  si->id = isf->iaccession;

#ifdef DEBUG2
  fprintf(stderr, "Contig pairs in scaffold " F_CID "\n", isf->iaccession);
  for(i = 0; i < isf->num_contig_pairs; i++)
    {
      PrintContigPair(&(isf->contig_pairs[i]), "", stderr);
    }
#endif

  if(isf->num_contig_pairs == 1 &&
     isf->contig_pairs[0].contig1 == isf->contig_pairs[0].contig2)
    {
      // singleton scaffold
      InstrumentScaffoldNextContig(graph, si, isf, verbose, printTo);
    }
  else
    {
      // multi-contig scaffold
      InstrumentScaffoldInitialContigPair(graph, si, isf, verbose, printTo);
      for(i = 1; i < isf->num_contig_pairs; i++)
        {
          InstrumentScaffoldNextGapAndContig(graph, si, isf, verbose, printTo);
        }
    }

#ifdef DEBUG
  fprintf(stderr, "Post-processing contig data from Scaffold " F_CID "\n",
          scaffold->id);
#endif
  // convert surrogate fragment positions from contig to scaffold coords

  //  fprintf(stderr,"Checking mate pairs within scaffold; cpHT = %x\n",si->cpHT);
  CheckFragmentMatePairs(graph,
                         &(si->bookkeeping),
                         &(si->surrogateTracker),
                         si->cpHT,
                         si->anchoredHT,
                         si->mates.mateStatus,
                         si->options,
                         &numMiso,
                         &numFar,
                         &numClose,
                         si->id);
  //  fprintf(stderr,"End checking mate pairs within scaffold\n");

#ifdef DUMP_MATE_PAIRS
  {
    FILE * fp;

    /*
      char siFile[1024];
      sprintf(siFile, "s" F_CID "Mates.txt", si->id);
      fprintf(GlobalData->stderrc, "Writing mates in scaffold " F_CID " to %s\n",
      si->id, siFile);
      fp = fopen(siFile, "w");
    */

    fp = fopen("sMates.txt", "a");
    assert(fp != NULL);
    PrintScaffoldInstrumenterMateDetails(si, fp, PRINTTABLE);
    
    fclose(fp);
  }
#endif
  return 0;
}


/*
  Scaffold:
  # gaps
  # negative
  # non-negative
  MIN (non-negative)
  max
  mean (excluding negatives)
  stddev
  # contigs, stats on contig sizes, contigs / scaffold
  min, max, mean, stddev
  # contigs with no reads/bac ends
  # unitigs, stats on unitig sizes, unitigs / contig
  min, max, mean, stddev
  # unitigs with no reads/bac ends
  # surrogates
  mates:
  happy intra-contig (in unitigs vs in surrogates)
  happy inter-contig, intra-scaffold (in unitigs vs in surrogates)
  missing - intra-scaffold gaps
  mis-oriented:
  innie vs. normal/anti-normal or outtie
  outtie vs. normal/anti-normal or innie
  normal/anti-normal vs. innie or outtie
  mis-separated
  fragments:
  reads in unitigs, reads in surrogates & stats
  bac ends...
  locales

  Contig:
  # unitigs, stats on unitig sizes, unitigs / contig
  min, max, mean, stddev
  # unitigs with no reads/bac ends
  # unitigs with no external data
  # surrogates
  mates:
  happy intra-contig (in unitigs vs in surrogates)
  happy inter-contig, intra-scaffold (in unitigs vs in surrogates)
  missing - intra-scaffold gaps
  mis-oriented:
  innie vs. normal/anti-normal or outtie
  outtie vs. normal/anti-normal or innie
  normal/anti-normal vs. innie or outtie
  mis-separated
  fragments:
  reads in unitigs, reads in surrogates & stats
  bac ends...
  locales
*/
int InstrumentScaffold(ScaffoldGraphT * graph,
                       CIScaffoldT * scaffold,
                       ScaffoldInstrumenter * si,
                       InstrumenterVerbosity verbose,
                       FILE * printTo)
{
  IntScaffoldMesg isf;

  if(graph == NULL ||
     scaffold == NULL ||
     si == NULL)
    {
      fprintf(stderr, "graph or scaffold or instrumenter is NULL!\n");
      return 1;
    }
  InitializeScaffoldInstrumenter(graph, si);

  if(printTo && verbose >= InstrumenterVerbose3)
    {
      fprintf(printTo, "Instrumenting Scaffold " F_CID "\n", scaffold->id);
    }

  // build a faux scaffold message - facilitates code reuse
  BuildFauxIntScaffoldMesgFromScaffold(graph, scaffold, si, &isf);
  
  InstrumentIntScaffoldMesg(graph, si, &isf, verbose, printTo);
  FreeIntScaffoldMesg(&isf);
  ComputeScaffoldInstrumenterStats(graph, si);
  
  if(printTo)
    {
      if(verbose >= InstrumenterVerbose3)
        {
          PrintScaffoldInstrumenter(graph, si, verbose, "\t", printTo);
        }
      else if(verbose == InstrumenterSilent)
        {
          PrintScaffoldGaps(si, printTo);
        }
    }
  return 0;
}


int FinishMissingMateList(ScaffoldGraphInstrumenter * sgi)
{
  MateDetail * wExtMates = GetVA_MateDetail(sgi->bookkeeping.wExtMates, 0);
  int32 numMatePairs = GetNumVA_MateDetail(sgi->bookkeeping.wExtMates);
  HashTable_AS * mateDetailHT;
  CDS_CID_t i;

  if(numMatePairs == 0)
    return 0;

  mateDetailHT = CreateScalarHashTable_AS(numMatePairs);
  if(mateDetailHT == NULL)
    {
      fprintf(stderr, "Failed to allocate hashtable of %d elements\n",
              numMatePairs);
      return 1;
    }

  // loop through all missing mates, add to hashtable, look up mate, ....
  for(i = 0; i < numMatePairs; i++)
    {
      MateDetail * mate;
      if((mate = (MateDetail *)LookupValueInHashTable_AS(mateDetailHT, (uint64)wExtMates[i].mateIID, 0)) == NULL)
        {
          if(wExtMates[i].mateChunkIID != NULLINDEX)
            {
              // fragment pair listed more than twice
            }
          InsertInHashTable_AS(mateDetailHT,
                               (uint64)wExtMates[i].fragIID, 0,
                               (uint64)&wExtMates[i], 0);
        }
      else
        {
          // populate both entries (mate and wExtMates[i])
          wExtMates[i].mateOffset5p = mate->fragOffset5p;
          wExtMates[i].mateChunkIID = mate->fragChunkIID;
          mate->mateOffset5p = wExtMates[i].fragOffset5p;
#ifdef TRACK_3P
          wExtMates[i].mateOffset3p = mate->fragOffset3p;
          mate->mateOffset3p = wExtMates[i].fragOffset3p;
#endif
          mate->mateChunkIID = wExtMates[i].fragChunkIID;
        }
    }
  
  if(numMatePairs > 1)
    qsort(wExtMates, numMatePairs, sizeof(MateDetail),
          (int (*) (const void *, const void *)) md2Compare);

  DeleteHashTable_AS(mateDetailHT);
  return 0;
}

/*
  Scaffold Graph:
  # singleton scaffolds
  min, max, mean, stddev
  # degenerate scaffolds, # without reads
  # non-degenerate singleton scaffolds, # unitigs/per
  # scaffolds, stats on sizes (and/or .cgm of scaffold sizes)
  min, max, mean, stddev
  # gaps, #gaps/scaffold, stats on gap sizes
  # negative
  # non-negative
  MIN (non-negative)
  max
  mean (excluding negatives)
  stddev
  # contigs, stats on contig sizes, contigs / scaffold
  min, max, mean, stddev
  # contigs with no reads/bac ends
  # unitigs, stats on unitig sizes, unitigs / contig
  min, max, mean, stddev
  # unitigs with no reads/bac ends
  # surrogates
  mates - by library:
  happy intra-contig (in unitigs vs in surrogates)
  happy inter-contig, intra-scaffold (in unitigs vs in surrogates)
  missing - intra-scaffold gaps
  mis-oriented:
  innie vs. normal/anti-normal or outtie
  outtie vs. normal/anti-normal or innie
  normal/anti-normal vs. innie or outtie
  mis-separated
  fragments:
  reads in unitigs, reads in surrogates
  bac ends...
  locales
*/
int InstrumentScaffoldGraph(ScaffoldGraphT * graph,
                            ScaffoldGraphInstrumenter * sgi,
                            CDS_COORD_t lowerLimit,
                            CDS_COORD_t upperLimit,
                            InstrumenterVerbosity verbose,
                            FILE * printTo)
{
  GraphNodeIterator scaffolds;
  CIScaffoldT * scaff;
  ScaffoldInstrumenter si;  // reused for each scaffold
  VA_TYPE(IID_Size) * iidSizes;
  int32 i;
  IID_Size * iidSize;

  if(graph == NULL || sgi == NULL)
    {
      fprintf(stderr, "graph or scaffold graph instrumenter is NULL!\n");
      return 1;
    }

  if(printTo && verbose >= InstrumenterVerbose1)
    {
      fprintf(printTo, "** Instrumenting Scaffold Graph\n");
    }
  // setup scaffold graph instrumenter
  InitializeScaffoldGraphInstrumenter(graph, sgi);

  // setup the reusable scaffold instrumenter
  memset(&si, 0, sizeof(ScaffoldInstrumenter));
  si.options = sgi->options;
  InitializeScaffoldInstrumenter(graph, &si);

  // allocate an array to hold iids & sizes
  iidSizes =
    CreateVA_IID_Size(GetNumVA_NodeCGW_T(graph->ScaffoldGraph->nodes));
  
  // loop over all scaffolds in the graph
  InitGraphNodeIterator(&scaffolds,
                        graph->ScaffoldGraph,
                        GRAPH_NODE_DEFAULT);
  while(NULL != (scaff = NextGraphNodeIterator(&scaffolds)))
    {
      if(scaff->flags.bits.isDead == FALSE && scaff->type == REAL_SCAFFOLD)
        {
          IID_Size iidSize;
          iidSize.iid = scaff->id;
          iidSize.length = (CDS_COORD_t) scaff->bpLength.mean;
          AppendVA_IID_Size(iidSizes, &iidSize);
        }
    }

  // sort the iidSizes by size - largest to smallest
  if(GetNumVA_IID_Size(iidSizes) > 1)
    {
      qsort(GetVA_IID_Size(iidSizes, 0),
            GetNumVA_IID_Size(iidSizes),
            sizeof(IID_Size),
            (int (*) (const void *, const void *)) sizeCompare );
    }

  {
    CDS_COORD_t largest, smallest;
    fprintf(stderr, "Instrumenting %d real, live scaffolds\n",
            (int) GetNumVA_IID_Size(iidSizes));
    iidSize = GetVA_IID_Size(iidSizes, 0);
    smallest = iidSize->length;
    iidSize = GetVA_IID_Size(iidSizes, GetNumVA_IID_Size(iidSizes) - 1);
    largest = MAX(smallest, iidSize->length);
    smallest = MIN(smallest, iidSize->length);
    fprintf(stderr, "Smallest: " F_COORD "bp, Largest: " F_COORD "bp\n", smallest, largest);
  }
  
  if(printTo && verbose == InstrumenterSilent)
    fprintf(printTo, "Scaffold Gap Sizes\n");
  
  for(i = 0; i < GetNumVA_IID_Size(iidSizes); i++)
    {
      iidSize = GetVA_IID_Size(iidSizes, i);

      if(iidSize->length >= lowerLimit && iidSize->length <= upperLimit)
        {
          scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, iidSize->iid);
          /*
            fprintf(stderr, "\r%d\t" F_CID "\t%15.0fbp",
            i + 1, scaff->id, scaff->bpLength.mean);
          */
          if(InstrumentScaffold(graph, scaff, &si, verbose, printTo))
            {
              fprintf(stderr,
                      "Failed to instrument scaffold " F_CID "\n",scaff->scaffoldID);
              return 1;
            }
      
          // consolidate scaffold data with scaffoldgraph data
          AddScaffoldToScaffoldGraphInstrumenter(graph, sgi, &si);
        }
    }
  fprintf(stderr, "\n");

  DeleteVA_IID_Size(iidSizes);
  
  // free scaffold instrumenter
  FreeScaffoldInstrumenter(&si);
  
  ComputeScaffoldGraphInstrumenterStats(graph, sgi);

  // bookkeeping wExtMates should list each missing mate fragment twice
  FinishMissingMateList(sgi);
    
  if(printTo && verbose >= InstrumenterVerbose1)
    {
      PrintScaffoldGraphInstrumenter(graph, sgi, verbose, printTo);
    }

  return 0;
}


/* Function to build IntScaffoldMesg from one of:
   if(terminalID == NULLINDEX)
   start at CI id's end & follow the edge until a branch is encountered

   BuildFauxIntScaffoldMesg()
*/   
int BuildFauxIntScaffoldMesg(ScaffoldGraphT * graph,
                             CDS_CID_t id, int end,
                             int terminalID, CIEdgeT * edge,
                             IntScaffoldMesg * isf)
{
  IntContigPairs icp;
  int nextEnd;
  int32 numEssential;
  ChunkInstanceT * nextCI;
  VA_TYPE(IntContigPairs) * icps = CreateVA_IntContigPairs(1000);
  assert(icps != NULL);

  // need a basis for constructing a faux or real scaffold
  if(terminalID == NULLINDEX && edge == NULL)
    return 1;
  
  // add the first pair of CIs to the list
  if(!edge)
    {
      GraphEdgeIterator edges;
    
      nextCI = GetGraphNode(graph->RezGraph, id);
      InitGraphEdgeIterator(graph->RezGraph,
                            nextCI->id,
                            end,
                            ALL_EDGES,
                            GRAPH_EDGE_DEFAULT,
                            &edges);
      while((edge = NextGraphEdgeIterator(&edges))!= NULL &&
            !getEssentialEdgeStatus(edge));
    }
  if(edge == NULL)
    return 1;
  PopulateICP(&icp, id, edge);
  AppendVA_IntContigPairs(icps, &icp);

  // get the next CI & the end of it that extends the path
  nextCI = GetGraphNode(graph->RezGraph, icp.contig2);
  nextEnd = (icp.orient == AB_AB || icp.orient == BA_AB) ? B_END : A_END;
  numEssential = (nextEnd == A_END) ?
    nextCI->numEssentialA : nextCI->numEssentialB;
  while(numEssential == 1)
    {
      CIEdgeT * nextEdge;
      GraphEdgeIterator edges;
    
      InitGraphEdgeIterator(graph->RezGraph,
                            nextCI->id,
                            nextEnd,
                            ALL_EDGES,
                            GRAPH_EDGE_DEFAULT,
                            &edges);
      while((nextEdge = NextGraphEdgeIterator(&edges))!= NULL &&
            !getEssentialEdgeStatus(nextEdge));
    
      PopulateICP(&icp, nextCI->id, nextEdge);
      AppendVA_IntContigPairs(icps, &icp);

      if(nextCI->id == terminalID)
        break;
    
      // get the next CI & the end of it that extends the path
      nextCI = GetGraphNode(graph->RezGraph, icp.contig2);
      nextEnd = (icp.orient == AB_AB || icp.orient == BA_AB) ? B_END : A_END;
      numEssential = (nextEnd == A_END) ?
        nextCI->numEssentialA : nextCI->numEssentialB;
    }
  FinishIntScaffoldMesg(isf, icps);
  DeleteVA_IntContigPairs(icps);
  return 0;
}


/*
  returns number of contigs instrumented (excluding thisCI)
*/
int InstrumentContigEnd(ScaffoldGraphT * graph,
                        ScaffoldInstrumenter * si,
                        ChunkInstanceT * thisCI,
                        int end)
{
  MateInstrumenter * mi;
  CIEdgeT * edge;
  GraphEdgeIterator edges;
  int32 numContigs = 1;

  // allocate a mate instrumenter to populate & return
  mi = CreateMateInstrumenter(graph, si->options);

  // loop over one end's essential edges & instrument each 'scaffold'
  InitGraphEdgeIterator(graph->RezGraph,
                        thisCI->id,
                        end,
                        ALL_EDGES,
                        GRAPH_EDGE_DEFAULT,
                        &edges);
  while((edge = NextGraphEdgeIterator(&edges))!= NULL)
    {
      IntScaffoldMesg isf;
      if(!getEssentialEdgeStatus(edge))
        continue;

  
      // instrument down this essential edge off this end of thisCI
      BuildFauxIntScaffoldMesg(graph, thisCI->id, end,
                               NULLINDEX, edge,
                               &isf);
      numContigs += isf.num_contig_pairs;
    
      InstrumentIntScaffoldMesg(graph, si, &isf, 0, NULL);
      FreeIntScaffoldMesg(&isf);
    
      // accumulate mate statuses
      AddMateInstrumenters(mi, &(si->mates));
    }
  // wipe the si's mate status data
  InitializeMateInstrumenter(graph, &(si->mates));
  AddMateInstrumenters(&(si->mates), mi);
  
  ComputeScaffoldInstrumenterStats(graph, si);
  
  DestroyMateInstrumenter(mi);
  return numContigs;
}


/*
  returns number of contigs instrumented (excluding thisCI)
*/
int InstrumentContigEndPartial(ScaffoldGraphT * graph,
                               ScaffoldInstrumenter * si,
                               ChunkInstanceT * thisCI,
                               int end,
                               int32 numContigs)
{
  MateInstrumenter * mi;
  CIEdgeT * edge;
  GraphEdgeIterator edges;

  // allocate a mate instrumenter to use
  if((mi = CreateMateInstrumenter(graph, si->options)) == NULL)
    {
      fprintf(stderr, "Failed to allocate mate instrumenter\n");
      return 1;
    }

  // loop over one end's essential edges & instrument each 'scaffold'
  InitGraphEdgeIterator(graph->RezGraph,
                        thisCI->id,
                        end,
                        ALL_EDGES,
                        GRAPH_EDGE_DEFAULT,
                        &edges);
  while((edge = NextGraphEdgeIterator(&edges))!= NULL)
    {
      IntScaffoldMesg isf;
      if(!getEssentialEdgeStatus(edge))
        continue;
  
      // instrument down this essential edge off this end of thisCI
      BuildFauxIntScaffoldMesg(graph, thisCI->id, end,
                               NULLINDEX, edge,
                               &isf);

      // if numContigs can't be instrumented, return failure
      if(isf.num_contig_pairs <= numContigs - 1)
        return 1;

      isf.num_contig_pairs = numContigs - 1;
      InstrumentIntScaffoldMesg(graph, si, &isf, 0, NULL);
      FreeIntScaffoldMesg(&isf);
    
      // accumulate mate statuses
      AddMateInstrumenters(mi, &(si->mates));
    }
  // wipe the si's mate status data
  InitializeMateInstrumenter(graph, &(si->mates));
  AddMateInstrumenters(&(si->mates), mi);
  
  ComputeScaffoldInstrumenterStats(graph, si);
  
  DestroyMateInstrumenter(mi);
  return 0;
}


int InstrumentContigPath(ScaffoldGraphT * graph,
                         ScaffoldInstrumenter * si,
                         CDS_CID_t firstID,
                         int firstEnd,
                         CDS_CID_t lastID)
{
  // build a faux scaffold from first to last & instrument
  // reuse BuildFauxIntScaffoldMesg code, but do so inefficiently
  // count number of contigs from firstCI to lastCI
  int thisEnd = firstEnd;
  IntScaffoldMesg isf;
  CDS_CID_t thisID = firstID;
  int done = 0;
  static IntScaffold_ID persistentScaffoldID = 0;
  VA_TYPE(IntContigPairs) * icps = CreateVA_IntContigPairs(1000);
  assert(icps != NULL);

  isf.iaccession = persistentScaffoldID++;
  
  // fprintf(stderr,  F_CID "(%c): ", firstID, (firstEnd == A_END) ? 'A' : 'B');
  while(!done)
    {
      ChunkInstanceT * thisCI = GetGraphNode(graph->RezGraph, thisID);
      CIEdgeT * edge;
      GraphEdgeIterator edges;

      /*
        if((thisEnd == 1 && thisCI->numEssentialA > 1) ||
        (thisEnd == 2 && thisCI->numEssentialB > 1))
        {
        fprintf(stderr, "Branch encountered. Instrumenting aborted.\n");
        return 1;
        }
      */
    
      // loop over one end's essential edges & instrument each 'scaffold'
      InitGraphEdgeIterator(graph->RezGraph,
                            thisCI->id,
                            thisEnd,
                            ALL_EDGES,
                            GRAPH_EDGE_DEFAULT,
                            &edges);
      if((thisEnd == 1 && thisCI->numEssentialA > 1) ||
         (thisEnd == 2 && thisCI->numEssentialB > 1))
        {
          int foundIt = 0;
          while((edge = NextGraphEdgeIterator(&edges))!= NULL)
            {
              if(getEssentialEdgeStatus(edge) &&
                 ((edge->idA == thisID && edge->idB == lastID) ||
                  (edge->idB == thisID && edge->idA == lastID)))
                {
                  foundIt = 1;
                  break;
                }
            }
          if(foundIt == 0)
            {
              fprintf(stderr,
                      "Unresolvable branch encountered. Instrumenting aborted.\n");
              return 1;
            }
          else
            {
              if(AddICP(icps, &thisID, &thisEnd, edge))
                {
                  fprintf(stderr, "Failed to add contig pairs!\n");
                  return 1;
                }
            }
        }
      else
        {
          while((edge = NextGraphEdgeIterator(&edges))!= NULL)
            {
              if(getEssentialEdgeStatus(edge))
                {
                  if(AddICP(icps, &thisID, &thisEnd, edge))
                    {
                      fprintf(stderr, "Failed to add contig pairs!\n");
                      return 1;
                    }
                }
            }
        }
      done = (thisID == lastID) ? 1 : 0;
    }
  FinishIntScaffoldMesg(&isf, icps);
  DeleteVA_IntContigPairs(icps);
  
  // instrument the scaffoldmesg
  // fprintf(stderr, "\n");
  InstrumentIntScaffoldMesg(graph, si, &isf, 0, NULL);
  FreeIntScaffoldMesg(&isf);
  ComputeScaffoldInstrumenterStats(graph, si);
  return 0;
}


void PrintEssentialEdges(ScaffoldGraphT * graph,
                         CDS_CID_t chunkID, int end)
{
  ChunkInstanceT * thisCI = GetGraphNode(graph->RezGraph, chunkID);
  CIEdgeT * edge;
  GraphEdgeIterator edges;
  int thisEnd = end;
  
  fprintf(stderr, F_CID "(%c): ", chunkID, (end == A_END) ? 'A' : 'B');

  // loop over one end's essential edges & instrument each 'scaffold'
  InitGraphEdgeIterator(graph->RezGraph,
                        thisCI->id,
                        thisEnd,
                        ALL_EDGES,
                        GRAPH_EDGE_DEFAULT,
                        &edges);
  while((edge = NextGraphEdgeIterator(&edges))!= NULL)
    {
      if(getEssentialEdgeStatus(edge))
        {
          thisEnd = GetOppositeEndOfOtherCI(edge, thisCI->id);
          if(thisEnd == NO_END)
            {
              fprintf(stderr, "\n");
              return;
            }
          fprintf(stderr, F_CID "(%.0f, %d) ",
                  (thisCI->id == edge->idA) ? edge->idB : edge->idA,
                  edge->distance.mean, thisEnd);
        }
    }
  fprintf(stderr, "\n");
}


/*
  Allocate char array to track scaffold IDs seen
  Loop over all CIs
  if ci's scaffoldID < numScaffoldIDs & not seen,
  mark it as seen
  go to left end of scaffold & build message & instrument
  if bad scaffold
  renumber each CI's scaffold to isolate it
*/
int AdjustCIScaffoldLabels(ScaffoldGraphT * graph,
                           int32 * numScaffoldIDs)
{
  char * scaffoldSeen;
  ScaffoldInstrumenter * si;
  int32 myNumScaffoldIDs = *numScaffoldIDs;
  ChunkInstanceT * firstCI;
  int firstEnd;
  GraphNodeIterator ciIterator;

  si = CreateScaffoldInstrumenter(graph, INST_OPT_INTER_MATES);
  assert(si != NULL);
  
  scaffoldSeen = (char *) safe_calloc( myNumScaffoldIDs, sizeof(char));

  // loop over all CIs
  InitGraphNodeIterator(&ciIterator,
                        graph->RezGraph,
                        GRAPH_NODE_DEFAULT);
  while(NULL != (firstCI = NextGraphNodeIterator(&ciIterator)))
    {
      int32 numCIs = 1;
    
      // only examine old CIs not yet seen
      if(firstCI->scaffoldID < *numScaffoldIDs &&
         !scaffoldSeen[firstCI->scaffoldID])
        {
          // go to left-most CI (off Aend of this CI)
          firstEnd = A_END;

          // if there's one essential edge to the left, we're not at the left end
          if(firstCI->numEssentialA == 1)
            {
              int atLeft = 0;
              CIEdgeT * edge;
              GraphEdgeIterator edges;

              // get to left-most (relative to the first CI encountered)
              while(!atLeft)
                {
                  // iterate over edges off relevant end to get next CI
                  InitGraphEdgeIterator(graph->RezGraph,
                                        firstCI->id,
                                        firstEnd,
                                        ALL_EDGES,
                                        GRAPH_EDGE_DEFAULT,
                                        &edges);
                  while((edge = NextGraphEdgeIterator(&edges))!= NULL)
                    {
                      if(getEssentialEdgeStatus(edge))
                        {
                          ChunkInstanceT * nextCI;
                          int nextEnd;

                          // get the CI to check it's scaffold ID
                          nextCI = GetGraphNode(graph->RezGraph,
                                                (edge->idA == firstCI->id) ?
                                                edge->idB : edge->idA);
                          if(nextCI->scaffoldID != firstCI->scaffoldID)
                            {
                              atLeft = 1;
                              break;
                            }
                          else
                            {
                              nextEnd = GetOppositeEndOfOtherCI(edge, firstCI->id);
                              if(nextEnd == NO_END)
                                return 1;
                              firstCI = nextCI;
                              firstEnd = nextEnd;
                              numCIs++;
                              break;
                            }
                        }
                    }
                }
            }

          /* at this point firstCI is the left-most CI of a scaffold
             and thisEnd is end away from the scaffold
             So, create scaffold message starting with firstCI &
             moving in direction opposite of thisEnd
          */
          if(numCIs > 1)
            {
              IntScaffoldMesg isf;
              ChunkInstanceT * thisCI = firstCI;
              int thisEnd = (firstEnd == A_END) ? B_END : A_END;
              ChunkInstanceT * nextCI;
              int done = 0;
              int foundNext;
              VA_TYPE(IntContigPairs) * icps = CreateVA_IntContigPairs(numCIs - 1);

              isf.iaccession = thisCI->scaffoldID;
              while(!done)
                {
                  CIEdgeT * edge;
                  GraphEdgeIterator edges;

                  foundNext = 0;
                  // loop over one end's essential edges & instrument each 'scaffold'
                  InitGraphEdgeIterator(graph->RezGraph,
                                        thisCI->id,
                                        thisEnd,
                                        ALL_EDGES,
                                        GRAPH_EDGE_DEFAULT,
                                        &edges);
                  while((edge = NextGraphEdgeIterator(&edges))!= NULL)
                    {
                      if(getEssentialEdgeStatus(edge))
                        {
                          // make sure scaffoldID of other CI is same
                          nextCI = GetGraphNode(graph->RezGraph,
                                                (edge->idA == thisCI->id) ?
                                                edge->idB : edge->idA);
                          if(nextCI->scaffoldID == thisCI->scaffoldID)
                            {
                              CDS_CID_t thisID = thisCI->id;
                              if(AddICP(icps, &thisID, &thisEnd, edge))
                                {
                                  fprintf(stderr, "Failed to add contig pairs!\n");
                                  return 1;
                                }
                              else
                                {
                                  thisCI = nextCI;
                                  foundNext = 1;
                                  break;
                                }
                            }
                        }
                    }
                  done = (foundNext == 1) ? 0 : 1;
                }

              // create the internal scaffold message
              FinishIntScaffoldMesg(&isf, icps);
  
              // instrument the scaffoldmesg
              InstrumentIntScaffoldMesg(graph, si, &isf, 0, NULL);
              FreeIntScaffoldMesg(&isf);
              ComputeScaffoldInstrumenterStats(graph, si);

              /* if the stats aren't okay, iterate over contigs &
                 assign them new scaffold IDs
              */
              {
                int32 badInterMates = GetMateStatsBad(&(si->mates.inter));
                int32 allInterMates =
                  badInterMates + GetMateStatsHappy(&(si->mates.inter));

                if(allInterMates > 0.5 &&
                   ((float) badInterMates) / allInterMates > 0.05 &&
                   isf.contig_pairs[0].contig1 != isf.contig_pairs[0].contig2)
                  {
                    int32 q;

                    fprintf(stderr, "**** Splitting scaffold " F_CID " into %d contigs:\n",
                            isf.iaccession, isf.num_contig_pairs + 1);
                    PrintScaffoldInstrumenter(graph, si, InstrumenterVerbose2, "\t", stderr);
                    for(q = 0; q < isf.num_contig_pairs; q++)
                      {
                        ContigT * contig = GetGraphNode(graph->RezGraph,
                                                        isf.contig_pairs[q].contig2);
                        contig->scaffoldID = myNumScaffoldIDs++;
                      }
                  }
              }
              DeleteVA_IntContigPairs(icps);
            }
        }
    }
  
  safe_free(scaffoldSeen);
  DestroyScaffoldInstrumenter(si);
  *numScaffoldIDs = myNumScaffoldIDs;
  return 0;
}


void PopulateICPContigs(ScaffoldGraphT * graph,
                        IntScaffoldMesg * ism,
                        CIScaffoldT * scaffold,
                        SEdgeT * sEdge,
                        LengthT lengthToAdd,
                        CDS_CID_t index)
{
  int isA = (sEdge->idA == scaffold->id);
  int sIsA2B = (sEdge->orient == AB_AB ||
                (isA && sEdge->orient == AB_BA) ||
                (!isA && sEdge->orient == BA_AB));
  double varFrom = sIsA2B ? 0.0 : scaffold->bpLength.variance;
  double meanFrom = sIsA2B ? 0.0 : scaffold->bpLength.mean;
  CDS_CID_t cIndex;
  CIScaffoldTIterator ciIterator;
  ChunkInstanceT * ci;
                                    
  InitCIScaffoldTIterator(graph, scaffold, sIsA2B, FALSE, &ciIterator);
  cIndex = 0;
  while((ci = NextCIScaffoldTIterator(&ciIterator)) != NULL)
    {
      ism->contig_pairs[cIndex + index].contig1 = ci->id;
    
      if((sIsA2B && ci->offsetAEnd.mean < ci->offsetBEnd.mean) ||
         (!sIsA2B && ci->offsetAEnd.mean > ci->offsetBEnd.mean))
        {
          // contig is A2B
          ism->contig_pairs[cIndex + index].orient = AB_AB;
          ism->contig_pairs[cIndex + index].mean =
            lengthToAdd.mean + (isA ? 0.0 : sEdge->distance.mean) +
            fabs(meanFrom - ci->offsetAEnd.mean);
          ism->contig_pairs[cIndex + index].stddev =
            sqrt(lengthToAdd.variance +
                 (isA ? 0.0 : sEdge->distance.variance) +
                 fabs(varFrom - ci->offsetAEnd.variance));
        }
      else
        {
          // contig is B2A
          ism->contig_pairs[cIndex + index].orient = BA_BA;
          ism->contig_pairs[cIndex + index].mean =
            lengthToAdd.mean + (isA ? 0.0 : sEdge->distance.mean) +
            fabs(meanFrom - ci->offsetBEnd.mean);
          ism->contig_pairs[cIndex + index].stddev =
            sqrt(lengthToAdd.variance +
                 (isA ? 0.0 : sEdge->distance.variance) +
                 fabs(varFrom - ci->offsetBEnd.variance));
        }
      cIndex++;
    }
}


static int ICPCompare(const IntContigPairs * a, const IntContigPairs * b)
{
  return (int) (a->mean - b->mean);
}


int InstrumentScaffoldPair(ScaffoldGraphT * graph,
                           SEdgeT * sEdge,
                           ScaffoldInstrumenter * si,
                           InstrumenterVerbosity verbose,
                           FILE * printTo)
{
  CIScaffoldT * scaffoldA = GetCIScaffoldT(graph->CIScaffolds, sEdge->idA);
  CIScaffoldT * scaffoldB = GetCIScaffoldT(graph->CIScaffolds, sEdge->idB);
  CDS_CID_t i;
  IntScaffoldMesg ism;
  LengthT offset = {0, 0};

  assert(scaffoldA != NULL && scaffoldB != NULL);

  ism.iaccession = 0;
  ism.num_contig_pairs = (scaffoldA->info.Scaffold.numElements +
                          scaffoldB->info.Scaffold.numElements) - 1;
  ism.contig_pairs =
    (IntContigPairs *) safe_malloc((ism.num_contig_pairs + 1) *
                                   sizeof(IntContigPairs));

  // populate contig_pairs with contigs, not contig pairs:
  {
    // first, determine starting position of A scaffold;
    // normally, this is 0;
    LengthT offset = {0, 0};
    // but if the overlap is longer than A (does this happen
    // only when A is contained by B?), then we need to put
    // the beginning of B at 0, and the beginning of A at
    // - length of overlap - length of A
    if(scaffoldA->bpLength.mean<-sEdge->distance.mean){
      offset.mean = -sEdge->distance.mean - scaffoldA->bpLength.mean;
      offset.variance = sEdge->distance.variance;
    }

    PopulateICPContigs(graph, &ism, scaffoldA, sEdge, offset, 0);

    offset.mean+=scaffoldA->bpLength.mean;
    offset.variance+=scaffoldA->bpLength.variance;

    PopulateICPContigs(graph, &ism, scaffoldB, sEdge, scaffoldA->bpLength,
		       scaffoldA->info.Scaffold.numElements);
  }
  // sort by start
  qsort(ism.contig_pairs,
        ism.num_contig_pairs + 1,
        sizeof(IntContigPairs),
        (int (*) (const void *, const void *)) ICPCompare );

  // convert to contig pairs
  for(i = 0; i < ism.num_contig_pairs; i++)
    {
      ChunkInstanceT * ci = GetGraphNode(graph->RezGraph,
                                         ism.contig_pairs[i].contig1);
      ism.contig_pairs[i].contig2 = ism.contig_pairs[i+1].contig1;
      if(ism.contig_pairs[i].orient == AB_AB)
        ism.contig_pairs[i].orient =
          (ism.contig_pairs[i+1].orient == AB_AB ? AB_AB : AB_BA);
      else
        ism.contig_pairs[i].orient =
          (ism.contig_pairs[i+1].orient == AB_AB ? BA_AB : BA_BA);
    
      ism.contig_pairs[i].mean = ism.contig_pairs[i+1].mean -
        ism.contig_pairs[i].mean - ci->bpLength.mean;
      ism.contig_pairs[i].stddev =
        sqrt(MAX(400., 
                 ism.contig_pairs[i+1].stddev * ism.contig_pairs[i+1].stddev -
                 ism.contig_pairs[i].stddev * ism.contig_pairs[i].stddev -
                 ci->bpLength.variance));
    }
  
#if 0
  fprintf(stderr,"Instrumenting a scaffold pair ... %d to %d, orient %s to %s, mean %g\n",
	  scaffoldA->id,scaffoldB->id,
	  ( sEdge->orient == AB_AB || sEdge->orient == AB_BA ) ? "AB" : "BA",
	  ( sEdge->orient == AB_AB || sEdge->orient == BA_AB ) ? "AB" : "BA",
	  sEdge->distance.mean);
#endif

  InstrumentIntScaffoldMesg(graph, si, &ism, verbose, printTo);
  ComputeScaffoldInstrumenterStats(graph, si);
  safe_free(ism.contig_pairs);
  return 0;
}

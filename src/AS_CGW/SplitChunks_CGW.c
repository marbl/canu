
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
static char CM_ID[] = "$Id: SplitChunks_CGW.c,v 1.7 2006-05-18 18:30:31 vrainish Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "InputDataTypes_CGW.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Input_CGW.h"
#include "Globals_CNS.h"
#include "MultiAlignment_CNS.h"
#include "SplitChunks_CGW.h"

//#define DEBUG
extern void CheckUnitigs(CDS_CID_t start, CDS_CID_t end);


/*
  Initial splitting of unitigs is performed as follows:

  For all, shortDiscriminatorUnique(+), or just DiscriminatorUnique unitigs?
  Only for unitigs >= min distance mean + CGW_CUTOFF * stddev in length
  
  1. create sequence coverage map of reads (AS_READ | AS_EXTR | AS_TRNR?)
    VR instead of explicit fragment type macro AS_FA_READ is used 
     - trim 30bp from each end of each read in creating map
     if, beyond initial 1x coverage & until last 1x coverage, there is
       at least one location with <=1x coverage then continue
  2. create good & bad clone coverage maps
     if at some point in seq coverage map where coverage <= 1 (not on ends)
     gcc <= 1 && bcc >= 2 split at that point/interval.
     - remember to scan entire unitig for possibility of multiple break points

  Splitting - make sure it will work for contigs & unitigs
              and for contig calling fn for unitig:
    Inputs: chunk index, interval where to break, graph, chunk type
    Outputs: if unitig, need to report info for contig breaking
  FOR UNITIGS
  1. original unitig becomes 'right' unitig
     add 2 new unitigs : chimeric middle & left
  2. from left to right, populate new unitigs/depopulate original
     until interval reached: frags go to left unitig
     until interval left: frags go to chimeric middle unitig
  3. finish unitig creation - A stat, consensus, unitig edges

  take right unitig & repeat checks for good/bad/sequence coverage
  & possibly split again
 */


// in creating coverage maps, trim bases off each end of fragment
#define READ_TRIM_BASES         30
#define MAX_SEQUENCE_COVERAGE    1
#define MIN_BAD_CLONE_COVERAGE   3
#define MAX_GOOD_CLONE_COVERAGE  0
#define LN_2                      .693147f

VA_DEF(SeqInterval)


#define SOURCE_LENGTH 100

typedef struct
{
  IntUnitigMesg    ium;
  /*
    ium fields for deltas, sequence, & quality are variable arrays
    for compatibility with the consensus interface
  */
  VA_TYPE(int32) * deltas;
  VA_TYPE(char)  * sequence;
  VA_TYPE(char)  * quality;
  CDS_COORD_t      minPos;
  int32            numFragsAllocated;
  int32            numRandomFragments;
} IUMStruct;


void PrintPositions(ScaffoldGraphT * graph, MultiAlignT * ma, int isUnitig)
{
  int i;
  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++)
  {
    IntMultiPos * imp = GetIntMultiPos(ma->f_list, i);
    InfoByIID * info = GetInfoByIID(graph->iidToFragIndex,
                                    imp->ident);
    CIFragT   * frag = GetCIFragT(graph->CIFrags,
                                  info->fragIndex);

    fprintf(stdout, F_IID ": (" F_COORD "," F_COORD ")\t",
            imp->ident, min(imp->position.bgn, imp->position.end), max(imp->position.bgn, imp->position.end));
    if(isUnitig)
      fprintf(stdout, "(%f,%f)\t",
              min(frag->offset5p.mean, frag->offset3p.mean),
              max(frag->offset5p.mean, frag->offset3p.mean));
    else
      fprintf(stdout, "(%f,%f)\t",
              min(frag->contigOffset5p.mean, frag->contigOffset3p.mean),
              max(frag->contigOffset5p.mean, frag->contigOffset3p.mean));
    if(imp->contained)
      fprintf(stdout, F_IID "\n", imp->contained);
    else
      fprintf(stdout, "\n");
  }
}


static void IncrementMapInterval(VA_TYPE(cds_uint16) * map,
                                 CDS_COORD_t minPos,
                                 CDS_COORD_t maxPos)
{
  CDS_COORD_t i;
  
  for(i = minPos; i < maxPos; i++)
  {
    cds_uint16 * temp = GetVA_cds_uint16(map, i);
    (*temp)++;
  }
}


static int CreateReadCoverageMap(ScaffoldGraphT * graph,
                                 VA_TYPE(cds_uint16) * rc,
                                 MultiAlignT * ma,
                                 int isUnitig)
{
  cds_uint32 i;
  
  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++)
  {
    IntMultiPos * imp = GetIntMultiPos(ma->f_list, i);

    // only add 'reads' to sequence coverage map
    
    if (AS_FA_READ(imp->type))  
    {
      InfoByIID * info = GetInfoByIID(graph->iidToFragIndex,
                                      imp->ident);
      CIFragT   * frag = GetCIFragT(graph->CIFrags,
                                    info->fragIndex);
      CDS_COORD_t minPos;
      CDS_COORD_t maxPos;

      // if unitig is forward vs. reverse, add to correct fragment location
      minPos = (isUnitig) ?
        (min(frag->offset5p.mean,
             frag->offset3p.mean) + READ_TRIM_BASES) :
        (min(frag->contigOffset5p.mean,
             frag->contigOffset3p.mean) + READ_TRIM_BASES);
      maxPos = (isUnitig) ?
        (max(frag->offset5p.mean,
             frag->offset3p.mean) - READ_TRIM_BASES) :
        (max(frag->contigOffset5p.mean,
             frag->contigOffset3p.mean) - READ_TRIM_BASES);
      IncrementMapInterval(rc, minPos, maxPos);
    }
  }
  return 0;
}


static int AddLinkToMaps(ScaffoldGraphT * graph,
                         VA_TYPE(cds_uint16) * gcc,
                         VA_TYPE(cds_uint16) * bcc,
                         CIFragT * frag,
                         CIFragT * mfrag,
                         OrientType orient,
                         CDS_CID_t distID,
                         CDS_COORD_t length,
                         int isUnitig)
{
  DistT * dist = GetDistT(graph->Dists, distID);
  CDS_COORD_t minPos = (isUnitig) ?
    (min(frag->offset5p.mean, mfrag->offset5p.mean)) :
    (min(frag->contigOffset5p.mean, mfrag->contigOffset5p.mean));
  CDS_COORD_t maxPos = (isUnitig) ?
    (max(frag->offset5p.mean, mfrag->offset5p.mean)) :
    (max(frag->contigOffset5p.mean, mfrag->contigOffset5p.mean));
  
  // if they're in the same unitig
  if((isUnitig && frag->cid == mfrag->cid) ||
     (!isUnitig && frag->contigID == mfrag->contigID))
  {
    // for pairs in the same unitig, just process the lesser
    if(frag->iid < mfrag->iid)
    {
      // if orientation is the same, the pair is bad
      if((isUnitig &&
          getCIFragOrient(frag) == getCIFragOrient(mfrag)) ||
         (!isUnitig &&
          GetContigFragOrient(frag) == GetContigFragOrient(mfrag)))
      {
        // bad pair - increment intervals
        if(orient == AS_INNIE)
        {
          // link is innie
          if((isUnitig && getCIFragOrient(frag) == A_B) ||
             (!isUnitig && GetContigFragOrient(frag) == A_B))
          {
            // both are A_B oriented, so intervals > 5p are suspect
            /* Increment from minPos to somewhere & maxPos to somewhere
               minPos                   maxPos
                  |                         |
                  -------->                 -------->
                  --- increment ---?        --- increment ---?
            */
            // note that maxPos < length
            IncrementMapInterval(bcc,
                                 minPos,
                                 min(minPos + dist->mean +
                                     CGW_CUTOFF * dist->stddev,
                                     maxPos));
            IncrementMapInterval(bcc,
                                 maxPos,
                                 min(maxPos + dist->mean +
                                     CGW_CUTOFF * dist->stddev,
                                     length - 1));
          }
          else
          {
            // both are B_A oriented, so intervals < 5p are suspect
            /* Increment from minPos to somewhere & maxPos to somewhere
                             minPos                   maxPos
                                |                         |
                        <--------                 <--------
               ?--- increment ---        ?--- increment ---
            */
            IncrementMapInterval(bcc,
                                 max(0,
                                     minPos - dist->mean -
                                     CGW_CUTOFF * dist->stddev),
                                 minPos);
            IncrementMapInterval(bcc,
                                 max(minPos,
                                     maxPos - dist->mean -
                                     CGW_CUTOFF * dist->stddev),
                                 maxPos);
          }
        }
        else
        {
          // link is outtie
          if((isUnitig && getCIFragOrient(frag) == A_B) ||
             (!isUnitig && GetContigFragOrient(frag) == A_B))
          {
            // both are A_B oriented, so
            /* Increment from minPos to somewhere & maxPos to somewhere
                           minPos                   maxPos
                              |                         |
                              -------->                 -------->
            ?--- increment ---        ?--- increment ---
            */
            IncrementMapInterval(bcc,
                                 max(0,
                                     minPos -
                                     dist->mean - CGW_CUTOFF * dist->stddev),
                                 minPos);
            IncrementMapInterval(bcc,
                                 max(minPos,
                                     maxPos -
                                     dist->mean - CGW_CUTOFF * dist->stddev),
                                 maxPos);
          }
          else
          {
            // both are B_A oriented, so
            /* Increment from minPos to somewhere & maxPos to somewhere
                 minPos                      maxPos
                    |                            |
            <--------                    <--------
                     --- increment ---?           --- increment ---?
            */
            IncrementMapInterval(bcc,
                                 minPos,
                                 min(minPos +
                                     dist->mean + CGW_CUTOFF * dist->stddev,
                                     maxPos));
            IncrementMapInterval(bcc,
                                 maxPos,
                                 min(length - 1,
                                     maxPos +
                                     dist->mean + CGW_CUTOFF * dist->stddev));
          }
        }
      }
      else
      {
        // fragments are oriented differently - may be okay
        CDS_COORD_t distance = maxPos - minPos;
        
        // if the distance is wrong, the pair is bad
        if(distance < dist->mean - CGW_CUTOFF * dist->stddev ||
           distance > dist->mean + CGW_CUTOFF * dist->stddev)
        {
          // bad pair
          // same intervals for either innie or outtie
          /*
              innie:
               minPos                                maxPos
                  |                                      |
                  ------->                       <--------
                  --- increment ---?    ?--- increment ---

              outtie:
                  minPos                                maxPos
                     |                                      |
              <-------                                      -------->
                      --- increment ---?  ?--- increment ---
          */
          IncrementMapInterval(bcc,
                               minPos,
                               min(maxPos,
                                   minPos +
                                   dist->mean + CGW_CUTOFF * dist->stddev));
          IncrementMapInterval(bcc,
                               max(minPos +
                                   dist->mean + CGW_CUTOFF * dist->stddev,
                                   maxPos -
                                   dist->mean - CGW_CUTOFF * dist->stddev),
                               maxPos);
        }
        else
        {
          // good pair
          IncrementMapInterval(gcc, minPos, maxPos);
        }
      }
    }
  }
  else
  {
    // in different unitigs, should be close to end of unitig
    if(orient == AS_INNIE)
    {
      if((isUnitig && getCIFragOrient(frag) == A_B &&
          frag->offset5p.mean <
          length - dist->mean - CGW_CUTOFF * dist->stddev) ||
          (!isUnitig && GetContigFragOrient(frag) == A_B &&
          frag->contigOffset5p.mean <
           length - dist->mean - CGW_CUTOFF * dist->stddev))
      {
        /*
                     --------->
                     --- increment ---
         */
        IncrementMapInterval(bcc,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean),
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean) +
                             dist->mean + CGW_CUTOFF * dist->stddev);
      }
      else if((isUnitig && getCIFragOrient(frag) == B_A &&
               frag->offset5p.mean >
               dist->mean + CGW_CUTOFF * dist->stddev) ||
              (!isUnitig && GetContigFragOrient(frag) == B_A &&
               frag->contigOffset5p.mean >
               dist->mean + CGW_CUTOFF * dist->stddev))
      {
        /*
                             <----------
                       --- increment ---
         */
        IncrementMapInterval(bcc,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean) -
                             dist->mean - CGW_CUTOFF * dist->stddev,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean));
      }
    }
    else
    {
      // outtie
      if((isUnitig && getCIFragOrient(frag) == B_A &&
          frag->offset5p.mean <
          length - dist->mean - CGW_CUTOFF * dist->stddev) ||
         (!isUnitig && GetContigFragOrient(frag) == B_A &&
          frag->contigOffset5p.mean <
          length - dist->mean - CGW_CUTOFF * dist->stddev))
      {
        /*
                   <----------
                              --- increment ---
         */
        IncrementMapInterval(bcc,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean),
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean) +
                             dist->mean + CGW_CUTOFF * dist->stddev);
                             
      }
      else if((isUnitig && getCIFragOrient(frag) == A_B &&
               frag->offset5p.mean >
               dist->mean + CGW_CUTOFF * dist->stddev) ||
              (!isUnitig && GetContigFragOrient(frag) == A_B &&
               frag->contigOffset5p.mean >
               dist->mean + CGW_CUTOFF * dist->stddev))
      {
      /*
                          --------->
         --- increment ---
       */
        IncrementMapInterval(bcc,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean) -
                             dist->mean - CGW_CUTOFF * dist->stddev,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean));
      }
    }
  }
  return 0;
}


static int CreateCloneCoverageMaps(ScaffoldGraphT * graph,
                                   VA_TYPE(cds_uint16) * gcc,
                                   VA_TYPE(cds_uint16) * bcc,
                                   MultiAlignT * ma,
                                   int isUnitig)
{
  cds_uint32 i;
  IntMultiPos * imp;
  InfoByIID * info;
  CIFragT * frag;
  
  // loop over fragments in unitig
  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++)
  {
    imp = GetIntMultiPos(ma->f_list, i);
    info = GetInfoByIID(graph->iidToFragIndex, imp->ident);
    frag = GetCIFragT(graph->CIFrags, info->fragIndex);

    // this is unlike ComputeMatePairStatistics, since we're interested
    // even in pairs that are in different unitigs/contigs
    if(frag->numLinks > 0)
    {
      if(frag->flags.bits.getLinksFromStore)
      {
        // going to the gatekeeper store should be the exception
        GateKeeperLinkRecordIterator gkpli;
        GateKeeperLinkRecord gkpl;
        CDS_CID_t linkHead = frag->linkHead;

        CreateGateKeeperLinkRecordIterator(graph->gkpStore.lnkStore,
                                           linkHead,
                                           frag->iid,
                                           &gkpli);
        while(NextGateKeeperLinkRecordIterator(&gkpli, &gkpl))
        {
          /* don't double count if clone is good
             do process bad clones from each fragment's perspective
             but don't cover the same bases twice?
             This is handled in AddLinkToMaps
          */
          if(gkpl.type != AS_REREAD)
          {
            InfoByIID * info = GetInfoByIID(graph->iidToFragIndex,
                                            ((frag->iid == gkpl.frag1) ?
                                             gkpl.frag2 : gkpl.frag1));
            AddLinkToMaps(graph, gcc, bcc, frag,
                          GetCIFragT(graph->CIFrags, info->fragIndex),
                          ((gkpl.orientation==AS_GKP_INNIE) ?
                           AS_INNIE : AS_OUTTIE),
                          (CDS_CID_t) gkpl.distance,
                          GetMultiAlignLength(ma),
                          isUnitig);
          }
        }
      }
      else
      {
        // get lone mate fragment
        if(frag->mateOf != NULLINDEX && frag->linkType != AS_REREAD)
        {
          AddLinkToMaps(graph, gcc, bcc, frag,
                        GetCIFragT(graph->CIFrags, frag->mateOf),
                        ((frag->flags.bits.innieMate) ? AS_INNIE : AS_OUTTIE),
                        frag->dist,
                        GetMultiAlignLength(ma),
                        isUnitig);
        }
      }
    }
  }
  return 0;
}


static int WriteMapToCelagramFile(VA_TYPE(cds_uint16) * map,
                                  CDS_COORD_t length,
                                  char * filestub,
                                  char * comment,
                                  CDS_CID_t index)
{
  char filename[1024];
  FILE * fp;
  CDS_COORD_t b;
  cds_uint16 n;
  
  sprintf(filename, "%s.%010" F_CIDP ".cgm", filestub, index);
  if((fp = fopen(filename, "w")) == NULL)
  {
    fprintf(GlobalData->stderrc, "Failed to open file %s for writing.\n",
            filename);
    return 1;
  }
  fprintf(fp, "%s " F_CID "\n", comment, index);
  for(b = 0; b < length; b++)
  {
    for(n = 0; n < *(GetVA_cds_uint16(map,b)); n++)
      fprintf(fp, F_COORD "\n", b);
  }
  fclose(fp);
  return 0;
}


static IUMStruct * CreateIUMStruct(void)
{
  IUMStruct * is;
  
  is = (IUMStruct *) safe_calloc(1, sizeof(IUMStruct));
  assert(is != NULL);
  
#ifdef AS_ENABLE_SOURCE
  is->ium.source = safe_malloc(SOURCE_LENGTH);
  memset(is->ium.source, (int) ' ', SOURCE_LENGTH);
  is->ium.source[SOURCE_LENGTH - 1] = '\0';
#endif
  is->deltas = CreateVA_int32(1);
  is->sequence = CreateVA_char(200000);
  is->quality = CreateVA_char(200000);
  is->minPos = CDS_COORD_MAX;
  return is;
}


static void ResetIUMStruct(IUMStruct * is)
{
  int i;
  for(i = 0; i < is->ium.num_frags; i++)
  {
#ifdef AS_ENABLE_SOURCE
    /*
    if(is->ium.f_list[i].source)
      free(is->ium.f_list[i].source);
    */
#endif
    is->ium.f_list[i].delta_length = 0;
    is->ium.f_list[i].delta = NULL;
  }
  is->ium.consensus = NULL;
  is->ium.quality = NULL;
  is->ium.num_frags = 0;

  is->minPos = CDS_COORD_MAX;
  is->numRandomFragments = 0;
}


static void FreeIUMStruct(IUMStruct * is)
{
  if(is)
  {
    if(is->deltas)
      DeleteVA_int32(is->deltas);
    if(is->sequence)
      DeleteVA_char(is->sequence);
    if(is->quality)
      DeleteVA_char(is->quality);
#ifdef AS_ENABLE_SOURCE
    if(is->ium.source)
      free(is->ium.source);
#endif
    if(is->ium.f_list)
    {
#ifdef AS_ENABLE_SOURCE
      /*
      int i;
      for(i = 0; i < is->ium.num_frags; i++)
      {
        if(is->ium.f_list[i].source)
          free(is->ium.f_list[i].source);
      }
      */
#endif
      free(is->ium.f_list);
    }
    free(is);
  }
}


/*
  Don't worry about containing-contained IMPs.
  That's what AdjustContainedOrder() is for. If it were handled here,
  there could be pathological cases where several IMPs have identical
  bgn & end positions with only some contained by others. qsort won't
  do all the comparisons and some contained IMPs might end up listed
  before their containing IMPs.
 */
static int positionCompare(const IntMultiPos * a, const IntMultiPos * b)
{
  return min(a->position.bgn, a->position.end) -
    min(b->position.bgn, b->position.end);
}


static int AddIMPToIUMStruct(IUMStruct * is, IntMultiPos * imp)
{
#define USE_UNGAPPED_POSITION
#ifdef USE_UNGAPPED_POSITION  
  InfoByIID * info = GetInfoByIID(ScaffoldGraph->iidToFragIndex,
				  imp->ident);
  CIFragT   * frag = GetCIFragT(ScaffoldGraph->CIFrags,
				info->fragIndex);
#endif
  if(is->ium.num_frags == is->numFragsAllocated)
  {
    is->ium.f_list = (IntMultiPos *) safe_realloc(is->ium.f_list,
                                             (is->numFragsAllocated + 10) *
                                             sizeof(IntMultiPos));
    assert(is->ium.f_list != NULL);
    is->numFragsAllocated += 10;
  }

  // copy the fields - not source, or delta, so no memcpy...
  is->ium.f_list[is->ium.num_frags].type = imp->type;
  is->ium.f_list[is->ium.num_frags].ident = imp->ident;
  is->ium.f_list[is->ium.num_frags].contained = imp->contained;
#ifdef USE_UNGAPPED_POSITION
  is->ium.f_list[is->ium.num_frags].position.bgn = frag->offset5p.mean;
  is->ium.f_list[is->ium.num_frags].position.end = frag->offset3p.mean;
#else
  is->ium.f_list[is->ium.num_frags].position = imp->position;
#endif
  is->ium.f_list[is->ium.num_frags].delta_length = 0;
  is->ium.f_list[is->ium.num_frags].delta = NULL;
  
  // update the min & max
#ifdef USE_UNGAPPED_POSITION
  //  is->minPos = min(frag->offset5p.mean,frag->offset3p.mean); /* 3/22/05 ALH: this looks wrong! */
  is->minPos = min(is->minPos,min(frag->offset5p.mean,frag->offset3p.mean));
#else
  is->minPos = min(is->minPos,min(imp->position.bgn, imp->position.end));
#endif

  is->ium.num_frags++;
  is->numRandomFragments += (AS_FA_RANDOM(imp->type)) ? 1 : 0; 
  return 0;
}


float EstimateGlobalFragmentArrivalRate(ChunkInstanceT * ci,
                                        MultiAlignT * ma)
{
  int i;
  int numRF = 0;
  cds_int32 rho = 0;

  // estimate random fragments/bp
  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++)
  {
    IntMultiPos * imp = GetIntMultiPos(ma->f_list, i);
    if (AS_FA_RANDOM(imp->type)) 
    {
      rho = max(rho, min(imp->position.bgn, imp->position.end));
      numRF++;
    }
  }
#ifdef DEBUG
  fprintf(GlobalData->stderrc,
          "EGFAR: rho = %d, numRandomFrags = %d, CI coverageStat = %d\n",
          rho, numRF, ci->info.CI.coverageStat);
#endif

  return((ci->info.CI.coverageStat + LN_2 * (numRF - 1)) / rho);
}


/* Keep contained fragments just after containing fragment
   1. Remove contained imps from impList
   2. Append them to the end of impList so that
      even contained imp A that contains another imp B will come before B
   3. Move contained imps up through impList to position them just
      after their containing imp

   How?
   1. copy contained imps into its own array
      while compacting non-contained imps in impList &
      adding non-contained imps to a hashtable
   2. looping until no more contained imps get moved to the end of impList,
      make sure contained imps that somehow don't have containing imps
      in the impList get moved to the end, too (orphans).
   3. Move contained imps up to the position after their containing imp
      Move orphaned contained imps up to a position where they should
      overlap the preceding imp
*/   
void AdjustContainedOrder(IntMultiPos * impList, int numIMPs)
{
  HashTable_AS * containing;
  IntMultiPos * contained;
  int numContained;
  int numContaining;
  int numMoved;
  int totalMoved;
  int numContainedLeft;
  int index;
  
  containing = CreateHashTable_int32_AS(numIMPs);
  contained = (IntMultiPos *) safe_malloc(numIMPs * sizeof(IntMultiPos));
  assert(containing != NULL && contained != NULL);
  
  for(index = 0, numContained = 0, numContaining = 0; index < numIMPs; index++)
  {
    if(impList[index].contained != 0)
    {
      // copy to contained array
      contained[numContained++] = impList[index];
    }
    else
    {
      // move up and add to hashtable
      impList[numContaining] = impList[index];
      InsertInHashTable_AS(containing,
                           (void *) &(impList[numContaining].ident),
                           sizeof(impList[numContaining].ident),
                           (void *) &(impList[numContaining].ident));
      numContaining++;
    }
  }

  /*
    Until all IMPs have been moved
    Iterate over contained IMP list
      If imp status is Waiting, look up containing IMP in hashtable
        if it's in the hashtable, change to MarkedToMove
    Iterate over contained IMP list
      If imp status is MarkedToMove, insert in impList just after
        containing imp and change status to AlreadyMoved & add to hashtable
   */
  totalMoved = 0;
  numMoved = 100;
  numContainedLeft = numContained;
  while(numMoved > 0)
  {
    numMoved = 0;
    // identify contained imps whose containing imp is in impList
    for(index = 0; index < numContainedLeft; index++)
    {
      // if this imp hasn't moved & it's containing imp is in the hashtable
      // then move it to the end of the impList & put it in the hashtable
      if(LookupInHashTable_AS(containing,
                              (void *) &(contained[index].contained),
                              sizeof(CDS_IID_t)))
      {
        impList[numContaining + totalMoved] = contained[index];
        InsertInHashTable_AS(containing,
                             (void *) &(impList[numContaining +
                                               totalMoved].ident),
                             sizeof(impList[numContaining +
                                               totalMoved].ident),
                             (void *) &(impList[numContaining +
                                               totalMoved].ident));
        numMoved++;
        totalMoved++;
      }
      else
      {
        // move toward the head of the array
        contained[index - numMoved] = contained[index];
      }
    }
    numContainedLeft -= numMoved;
  }
  DeleteHashTable_AS(containing);

  // if numContainedLeft != 0, there is a contained missing its containing
  /*
  for(index = 0; index < numContainedLeft; index++)
  {
    fprintf(GlobalData->stderrc,
            "Contained imp " F_IID " will be in separate unitig from containing imp " F_IID "\n",
            contained[index].ident,
            contained[index].contained);
    contained[index].contained = 0;
    impList[numContaining + totalMoved + index] = contained[index];
  }
  */
  free(contained);

  assert(numContainedLeft == 0);
  
  // Now move containeds & orphaned containeds up to behind their containing
  for(index = numContaining; index < numIMPs; index++)
  {
    int index2;
    IntMultiPos imp = impList[index];
      
    // move up to behind its containing
    if(impList[index].contained != 0)
    {
      for(index2 = index - 1;
          impList[index2].ident != imp.contained &&
            impList[index2].contained != imp.contained;
          index2--)
      {
        impList[index2 + 1] = impList[index2];
      }
    }
    else
    {
      // move up to behind something it should overlap
      for(index2 = index - 1;
          index2 >= 0 &&
            (min(impList[index2].position.bgn, impList[index2].position.end) >
             min(imp.position.bgn, imp.position.end) ||
             max(impList[index2].position.bgn, impList[index2].position.end) <
             min(imp.position.bgn, imp.position.end));
          index2--)
      {
        impList[index2 + 1] = impList[index2];
      }
    }
    impList[index2 + 1] = imp;
  }
}


static int StoreIUMStruct(ScaffoldGraphT * graph,
                          IUMStruct * is,
                          int isUnitig,
                          UnitigStatus status,
                          float egfar,
                          int cgbType)
{
  int i;
  cds_int32 rho = 0;
  
  is->ium.iaccession = ((isUnitig) ? GetNumGraphNodes(graph->CIGraph) :
                        GetNumGraphNodes(graph->ContigGraph));

  is->ium.status = status;

  // re-position fragments & compute (estimated) coverage statistic
  is->ium.length = 0;
  for(i = 0; i < is->ium.num_frags; i++)
  {
    is->ium.f_list[i].position.bgn -= is->minPos;
    is->ium.f_list[i].position.end -= is->minPos;
    if(is->ium.f_list[i].position.end > is->ium.f_list[i].position.bgn)
    {
      rho = max(rho,is->ium.f_list[i].position.bgn);
      is->ium.length = max(is->ium.length, is->ium.f_list[i].position.end);
      is->ium.f_list[i].position.end--;
    }
    else
    {
      rho = max(rho,is->ium.f_list[i].position.end);
      is->ium.length = max(is->ium.length, is->ium.f_list[i].position.bgn);
      is->ium.f_list[i].position.bgn--;
    }
    // check that containing fragment is present, if specified
    // NOTE: ugly search
    if(is->ium.f_list[i].contained)
    {
      int j;
      for(j = 0; j < is->ium.num_frags; j++)
      {
        if(is->ium.f_list[j].ident == is->ium.f_list[i].contained)
          break;
      }
      if(j == is->ium.num_frags)
        is->ium.f_list[i].contained = 0;
    }
  }
  is->ium.a_branch_point = is->ium.b_branch_point = 0;


  /* egfar is the estimated global fragment arrival rate based on the original
     chunk.
     rho is the sum of the a-hangs (= largest starting position)
     Based on compute_coverage_statistic() in AS_CGB_cgb.h
  */
  is->ium.coverage_stat = ((egfar > 0.f && is->numRandomFragments > 0) ?
                           (rho * egfar - LN_2 *
                            (is->numRandomFragments - 1)) :
                           0.f);

#ifdef DEBUG
  fprintf(GlobalData->stderrc, "accession = %d\n", is -> ium . iaccession);
  fprintf(GlobalData->stderrc, "length = " F_COORD ", egfar = %f, rho = %d, num random frags = %d, coverage stat = %f\n",
          is->ium.length,
          egfar, rho, is->numRandomFragments,
          is->ium.coverage_stat);
#endif
  
  is->ium.forced = 0;

  // get the multi-alignment - this fills in some unitig fields
  qsort(is->ium.f_list,
        is->ium.num_frags,
        sizeof(IntMultiPos),
        (int (*) (const void *, const void *)) positionCompare );
  AdjustContainedOrder(is->ium.f_list, is->ium.num_frags);
  
  if(MultiAlignUnitig(&(is->ium),
                      ScaffoldGraph->fragStore,
                      is->sequence,
                      is->quality,
                      is->deltas,
                      CNS_STATS_ONLY,
                      0,
                      Local_Overlap_AS_forCNS,      //  DP_Compare
                      NULL)) {
    fprintf(GlobalData->stderrc,
            "FATAL ERROR: MultiAlignUnitig call failed in unitig splitting\n");
    assert(FALSE);
  }
  
  // add ium to the system
  {
    MultiAlignT *ma = CreateMultiAlignTFromIUM(&(is->ium),
                                               GetNumCIFragTs(graph->CIFrags),
                                               FALSE);
    CDS_COORD_t length = GetMultiAlignUngappedLength(ma);
    ChunkInstanceT * ci;
    int hasFbacFrags = 0;
    
    // need to point fragments to their new unitig/contig
    for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++)
    {
      IntMultiPos * imp = GetIntMultiPos(ma->f_list, i);
      InfoByIID * info = GetInfoByIID(graph->iidToFragIndex,
                                      imp->ident);
      CIFragT   * frag = GetCIFragT(graph->CIFrags,
                                    info->fragIndex);

      hasFbacFrags = (frag->type == AS_FBAC) ? 1 : hasFbacFrags;
#ifdef AS_ENABLE_SOURCE
      imp->sourceInt = info->fragIndex;
#endif
    }
    
    // This adds a reference only if keepInCache is true...
    InsertMultiAlignTInSequenceDB(graph->sequenceDB, is->ium.iaccession,
                                  TRUE, ma, TRUE);
    DuplicateEntryInSequenceDB(graph->sequenceDB,
                               is->ium.iaccession, TRUE,
                               is->ium.iaccession, FALSE, FALSE);
    ProcessIUM_ScaffoldGraph(&(is->ium), length, TRUE);
    /*
    ma = LoadMultiAlignTFromSequenceDB(graph->sequenceDB,
                                       is->ium.iaccession,
                                       TRUE);
    */
    
    ci = GetGraphNode((isUnitig ? graph->CIGraph : graph->ContigGraph),
                      is->ium.iaccession);
    ci->flags.bits.cgbType = cgbType;
    ci->aEndCoord = ci->bEndCoord = 0;
    ci->flags.bits.includesFinishedBacFragments = hasFbacFrags;
/*
    if(is->ium.num_frags < 2)
      ci->flags.bits.isChaff = TRUE;
*/

    UpdateNodeFragments((isUnitig ? graph->CIGraph : graph->ContigGraph),
                        is->ium.iaccession,
                        ci->type == DISCRIMINATORUNIQUECHUNK_CGW,
                        TRUE);
    
    UnloadMultiAlignTFromSequenceDB(graph->sequenceDB,
                                    is->ium.iaccession, TRUE);
  }

  return 0;
}


// NOTE: csis intervals are in ungapped coordinates
static int SplitChunkByIntervals(ScaffoldGraphT * graph,
                                 CDS_CID_t ciID,
                                 MultiAlignT * ma,
                                 VA_TYPE(SeqInterval) * csis,
                                 int isUnitig)
{
  int currI = 0;
  SeqInterval * currInterval = GetVA_SeqInterval(csis,currI);
  int i;
  CDS_CID_t  lastIID = 0;
  IUMStruct * lastPlacement = NULL;
  IUMStruct * good = CreateIUMStruct();
  IUMStruct * bad = CreateIUMStruct();
  float egfar;
  ChunkInstanceT * ci = GetGraphNode(graph->CIGraph, ciID);

  // not implemented for contigs yet
  assert(isUnitig);

  // fprintf(stdout, "%10" F_IIDP ": Pre-sort\n", ma->id);
  // PrintPositions(graph, ma, 1);
  // NOTE: is this safe?
  // sort imps by starting position
  qsort(ma->f_list->Elements,
        GetNumIntMultiPoss(ma->f_list),
        sizeof(IntMultiPos),
        (int (*) (const void *, const void *)) positionCompare );
  egfar = EstimateGlobalFragmentArrivalRate(ci, ma);

  // need to have contained fragments listed right after containing fragment
  // fprintf(stdout, "%10" F_IIDP ": Pre-adjust\n", ma->id);
  // PrintPositions(graph, ma, 1);
  AdjustContainedOrder((IntMultiPos *) ma->f_list->Elements,
                       GetNumIntMultiPoss(ma->f_list));
  // fprintf(stdout, "%10" F_IIDP ": Final\n", ma->id);
  // PrintPositions(graph, ma, 1);
  
#ifdef DEBUG
  fprintf(GlobalData->stderrc, "egfar = %f\n", egfar);
#endif
  
  // feedback for log file
  fprintf(GlobalData->stderrc,
          "Splitting %s " F_CID " into as many as %d %s at intervals:\n",
          (isUnitig ? "unitig" : "contig"), ma->id,
          (int) (2 * GetNumVA_SeqInterval(csis) + 1),
          (isUnitig ? "unitigs" : "contigs"));

  for(currI = 0; currI < GetNumVA_SeqInterval(csis); currI++)
  {
    currInterval = GetVA_SeqInterval(csis, currI);
    fprintf(GlobalData->stderrc, "\t" F_COORD "," F_COORD "\n",
            currInterval->bgn, currInterval->end);
  }

  /*
    loop over fragments, placing them in different chunks
    Two chunks are populated at the same stage - a good & a bad - until the
    last chimeric interval is passed, in which case only a good is populated.
    So, when the current bad interval is exceeded, the next good/bad pair
    must be used. However, intervals may be close together to the point that
    the very next chimeric interval may not be the furthest chimeric interval
    into which the fragment falls...
  */
  currI = 0;
  currInterval = GetVA_SeqInterval(csis, currI);
  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++)
  {
    IntMultiPos * imp = GetIntMultiPos(ma->f_list, i);
    InfoByIID * info = GetInfoByIID(graph->iidToFragIndex,
                                    imp->ident);
    CIFragT   * frag = GetCIFragT(graph->CIFrags,
                                  info->fragIndex);
    CDS_COORD_t minPos =
      ((isUnitig) ?
       min(frag->offset5p.mean, frag->offset3p.mean) :
       min(frag->contigOffset5p.mean, frag->contigOffset3p.mean));
    CDS_COORD_t maxPos =
      ((isUnitig) ?
       max(frag->offset5p.mean, frag->offset3p.mean) :
       max(frag->contigOffset5p.mean, frag->contigOffset3p.mean));

#ifdef DEBUG
    fprintf(GlobalData->stderrc, "Interval: " F_COORD ", " F_COORD "\tfragment(" F_IID "): " F_COORD ", " F_COORD " (" F_COORD ", " F_COORD ")\n",
            currInterval->bgn, currInterval->end, imp->ident,
            minPos, maxPos, imp->position.bgn, imp->position.end);
#endif
    fprintf(stderr, "Interval: " F_COORD ", " F_COORD "\tfragment(" F_IID "): " F_COORD ", " F_COORD " (" F_COORD ", " F_COORD ")\n",
            currInterval->bgn, currInterval->end, imp->ident,
            minPos, maxPos, imp->position.bgn, imp->position.end);

    /* keep contained with containing - avoids problem of consensus failures
       if overlap of contained & next fragment is very short
    */
    if(imp->contained)
    {
      AddIMPToIUMStruct(lastPlacement, imp);
      continue;
    }

    lastIID = imp->ident;
    // determine where the fragment goes wrt the chimeric interval
    if(minPos >= currInterval->bgn)
    {
      if(minPos > currInterval->end)
      {
        /* Three possibilities here:
           1. we've run out of intervals & the rest of the fragments
              go in the last good interval (handled by the following if())
           2. the fragment is not beyond the next bad interval
           3. the fragment is beyond the next bad interval
        */
        if(currI < GetNumVA_SeqInterval(csis))
        {
          SeqInterval * lastInterval = NULL;
          /* chimeric intervals may be close together - close enough that
             the current fragment may be in the next 'bad' interval rather
             than the next 'good' interval - or even further down...
             we ARE done with the current good/bad IUMs.
          */
          // old bad & old good are done
          if(good->ium.num_frags > 0)
          {
            StoreIUMStruct(graph, good, isUnitig, AS_UNASSIGNED, egfar,
                           ci->flags.bits.cgbType);
            // refresh ci, since the VarArray may have been resized
            ci = GetGraphNode(graph->CIGraph, ciID);
          }
            
          if(bad->ium.num_frags > 0)
          {
            StoreIUMStruct(graph, bad, isUnitig, AS_UNASSIGNED, egfar,
                           ci->flags.bits.cgbType);
            // refresh ci, since the VarArray may have been resized
            ci = GetGraphNode(graph->CIGraph, ciID);
          }

          // reset the old
          ResetIUMStruct(good);
          ResetIUMStruct(bad);
          
          // advance to the next relevant interval
          while(currI < GetNumVA_SeqInterval(csis) &&
                minPos > currInterval->end)
          {
            lastInterval = currInterval;
            currInterval = GetVA_SeqInterval(csis, ++currI);
          }

          if(currInterval)
          {
            // currInterval is relevant bad interval
            // place frag in good or bad
            if(minPos >= currInterval->bgn)
            {
              if(minPos > currInterval->end)
              {
                AddIMPToIUMStruct((lastPlacement = good), imp);
              }
              else
              {
                AddIMPToIUMStruct((lastPlacement = bad), imp);
              }
            }
            else
            {
              if(maxPos >= currInterval->bgn)
              {
                AddIMPToIUMStruct((lastPlacement = bad), imp);
              }
              else
              {
                AddIMPToIUMStruct((lastPlacement = good), imp);
              }
            }
          }
          else
          {
            // ran out of chimeric intervals - in last good interval
            AddIMPToIUMStruct((lastPlacement = good), imp);
            currInterval = lastInterval;
            currI = GetNumVA_SeqInterval(csis);
          }
        }
        else
        {
          AddIMPToIUMStruct((lastPlacement = good), imp);
        }
      }
      else
      {
        // goes in bad
        AddIMPToIUMStruct((lastPlacement = bad), imp);
      }
    }
    else
    {
      if(maxPos >= currInterval->bgn)
      {
        // goes in bad
        AddIMPToIUMStruct((lastPlacement = bad), imp);
      }
      else
      {
        // goes in good
        AddIMPToIUMStruct((lastPlacement = good), imp);
      }
    }
  }
  
  // now make the last good chunk a chunk
  if(good->ium.num_frags > 0)
  {
    StoreIUMStruct(graph, good, isUnitig, AS_UNASSIGNED, egfar,
                   ci->flags.bits.cgbType);
    // refresh ci, since the VarArray may have been resized
    ci = GetGraphNode(graph->CIGraph, ciID);
  }
  
  // if the last bad section went to the end, it's still here
  if(bad->ium.num_frags > 0)
  {
    StoreIUMStruct(graph, bad, isUnitig, AS_UNASSIGNED, egfar,
                   ci->flags.bits.cgbType);
    // refresh ci, since the VarArray may have been resized
    ci = GetGraphNode(graph->CIGraph, ciID);
  }
  
  // delete the original unitig & contig
  if(ci->info.CI.contigID != NULLINDEX)
  {
    ChunkInstanceT * contig = GetGraphNode(graph->ContigGraph,
                                           ci->info.CI.contigID);
    if(contig)
      DeleteGraphNode(graph->ContigGraph, contig);
  }
  DeleteGraphNode(graph->CIGraph, ci);
  
  return 0;
}


int SplitInputUnitigs(ScaffoldGraphT * graph)
{
  VA_TYPE(cds_uint16) * rc = NULL;   // read coverage map
  VA_TYPE(cds_uint16) * gcc = NULL;  // good clone coverage
  VA_TYPE(cds_uint16) * bcc = NULL;  // bad clone coverage
  VA_TYPE(SeqInterval) * csis = NULL; // chimeric sequence intervals
  CDS_COORD_t minLength = CDS_COORD_MAX;
  int i;
  int numCIs = GetNumGraphNodes(graph->CIGraph);

  rc = CreateVA_cds_uint16(10000);
  gcc = CreateVA_cds_uint16(10000);
  bcc = CreateVA_cds_uint16(10000);
  csis = CreateVA_SeqInterval(100);
  assert(rc != NULL && gcc != NULL && bcc != NULL && csis != NULL);

  // set terminate_cond to 1 - this is needed so that if MultiAlignUnitig()
  // fails, the program will exit with a non-zero
  terminate_cond = 1;
  
  // determine minimum unitig length of interest
  // NOTE: consider adding a multipler to dptr->upper
  for(i = 0; i < GetNumDistTs(graph->Dists); i++)
  {
    DistT * dptr = GetDistT(graph->Dists,i);
    minLength = min(minLength,dptr->mean + CGW_CUTOFF * dptr->stddev);
  }

#ifdef DEBUG
  fprintf(GlobalData->stderrc, "%10d unitigs.\n", numCIs);
  fprintf(GlobalData->stderrc, "minimum length: " F_COORD "\n", minLength);
  minLength = max(minLength, 2000 + CGW_CUTOFF * 200);
#endif
  
  // loop over number of original unitigs - not new ones being generated
  // otherwise would use an iterator
  for(i = 0; i < numCIs; i++)
  {
    ChunkInstanceT * ci = GetGraphNode(graph->CIGraph, i);
    MultiAlignT * ma = LoadMultiAlignTFromSequenceDB(graph->sequenceDB,
                                                     ci->id,
                                                     TRUE);
#ifdef DEBUG2
    if(i > 0 && i % 1000 == 0)
      fprintf(GlobalData->stderrc, "%d\r", i);
#endif
    
    // NOTE: add discriminator statistic checks?
    if(GetMultiAlignLength(ma) >= minLength)
    {
      CDS_COORD_t minBase, maxBase;
      CDS_COORD_t curBase;

      // create read coverage map for unitig
      // function intended to work for either contig or unitig
      ResetVA_cds_uint16(rc);
      EnableRangeVA_cds_uint16(rc, GetMultiAlignLength(ma));
#ifdef DEBUG2
      fprintf(GlobalData->stderrc, "Unitig %10d, length %d\n",
              i, (int) GetMultiAlignLength(ma));
      fprintf(GlobalData->stderrc, "\tVA_num = %d\n",
              (int) GetNumVA_cds_uint16(rc));
#endif
      CreateReadCoverageMap(graph, rc, ma, TRUE);

#ifdef DEBUG3
      // abuse the celagram file concept
      WriteMapToCelagramFile(rc, (CDS_COORD_t) GetMultiAlignLength(ma),
                             "rc", "Read coverage for unitig", i);
#endif
      
      // locate region within which to look for possible chimeric points
      // i.e., ignore initial & trailing 0/1 values
      for(minBase = READ_TRIM_BASES;
          minBase < GetMultiAlignLength(ma) - READ_TRIM_BASES &&
          *(GetVA_cds_uint16(rc,minBase)) <= 1;
          minBase++);
      for(maxBase = GetMultiAlignLength(ma) - READ_TRIM_BASES;
          maxBase > READ_TRIM_BASES &&
          *(GetVA_cds_uint16(rc,maxBase)) <= 1;
          maxBase--);

      // see if there is a candidate interval
      for( curBase = minBase; curBase < maxBase; curBase++)
        if(*(GetVA_cds_uint16(rc,curBase)) <= MAX_SEQUENCE_COVERAGE)
          break;

#ifdef DEBUG2
      fprintf(GlobalData->stderrc,
              "minBase = " F_COORD ", maxBase = " F_COORD ", curBase = " F_COORD "\n",
              minBase, maxBase, curBase);
#endif
      
      // see if above loop ended in a candidate interval
      if(curBase < maxBase)
      {
        CDS_COORD_t checkBase;
        SeqInterval interval;
        int inInterval = 0;
        
#ifdef DEBUG2
        fprintf(GlobalData->stderrc, "Candidate interval found!\n");
#endif
        
        // create gcc & bcc maps for unitig
        ResetVA_cds_uint16(gcc);
        EnableRangeVA_cds_uint16(gcc, GetMultiAlignLength(ma));
        ResetVA_cds_uint16(bcc);
        EnableRangeVA_cds_uint16(bcc, GetMultiAlignLength(ma));

        // reset the number of splitting intervals to 0
        ResetVA_SeqInterval(csis);
        EnableRangeVA_SeqInterval(csis, 0);
        
        // function intended to work for either contig or unitig
        CreateCloneCoverageMaps(graph, gcc, bcc, ma, TRUE);
        
        // identify & count chimeric sequence intervals
        for(checkBase = curBase; checkBase < maxBase; checkBase++)
        {
          if(*(GetVA_cds_uint16(rc,checkBase)) <= MAX_SEQUENCE_COVERAGE &&
             *(GetVA_cds_uint16(bcc,checkBase)) >= MIN_BAD_CLONE_COVERAGE &&
             *(GetVA_cds_uint16(gcc,checkBase)) <= MAX_GOOD_CLONE_COVERAGE)
          {
            if(inInterval)
            {
              // continuing in interval
              interval.end = checkBase;
            }
            else
            {
              // starting interval
              interval.bgn = interval.end = checkBase;
              inInterval = 1;
            }
          }
          else
          {
            if(inInterval)
            {
              // ended interval
              /* if it is more than a minimum distance from the last one,
                 add it. otherwise, combine the two - since no fragment can
                 possibly be entirely within the intervening 'good' interval
              */
              if(GetNumVA_SeqInterval(csis) > 0)
              {
                SeqInterval * tempSI =
                  GetVA_SeqInterval(csis, GetNumVA_SeqInterval(csis) - 1);
                if(tempSI->end < interval.bgn - AS_FRAG_MIN_LEN)
                  AppendVA_SeqInterval(csis, &interval);
                else
                  tempSI->end = interval.end;
              }
              else
                AppendVA_SeqInterval(csis, &interval);
              inInterval = 0;
            }
            // otherwise continuing not in interval - do nothing
          }
        }

        if(GetNumVA_SeqInterval(csis) > 0)
        {
#ifdef DEBUG2
          WriteMapToCelagramFile(rc,
                                 (cds_uint32) GetMultiAlignLength(ma),
                                 "rc",
                                 "Read coverage for unitig", i);
          WriteMapToCelagramFile(gcc,
                                 (cds_uint32) GetMultiAlignLength(ma),
                                 "gcc",
                                 "Good clone coverage for unitig", i);
          WriteMapToCelagramFile(bcc,
                                 (cds_uint32) GetMultiAlignLength(ma),
                                 "bcc",
                                 "Bad clone coverage for unitig", i);
#endif
          // compute global fragment arrival rate
          SplitChunkByIntervals(graph, ci->id, ma, csis, TRUE);
          // refresh ci, since the VarArray may have been resized
          ci = GetGraphNode(graph->CIGraph, i);
        }
      }
    }    
    UnloadMultiAlignTFromSequenceDB(graph->sequenceDB, ci->id, TRUE);
  }

  DeleteVA_SeqInterval(csis);
  DeleteVA_cds_uint16(rc);
  DeleteVA_cds_uint16(gcc);
  DeleteVA_cds_uint16(bcc);

  //CheckUnitigs(0, GetNumGraphNodes(graph->CIGraph));
  
  return 0;
}


/*
  Operate on contigs.
  Create variable array with 2*numContig entries
  Create hashtable with 2*numContig entries
    key is contigID, contigEnd pair
    value is pointer into array
  Create variable array for satisfied mate interval map
  Iterate over contigs
    reset array & hashtable
    iterate over fragments in contig
      access mate of this fragment
      if(in same contig)
        if fragIID < mateIID, add link to clone coverage map
      else
        determine if edge extends off AEnd or BEnd of contig
        lookup other contigID, end pair in hashtable
        create/update entry with 3' end of fragment
    sort elements in array by left offset
    iterate over entries in sorted array
      track A) 

        
  Should return variable array of contig/interval objects
 */

typedef struct
{
  CDS_CID_t id;
  int end;
} CELKey;

typedef struct
{
  CELKey      key;
  int32       numEdges;
  CDS_COORD_t leftCoord;
  CDS_COORD_t rightCoord;
} ChunkEndLink;
VA_DEF(ChunkEndLink)
  
int CELHashFn(const void * item, int length)
{
  return Hash_AS((uint8 *) item, length, 37);
}

int CELCompareFn(const void * item1, const void * item2)
{
  const ChunkEndLink * cel1 = (const ChunkEndLink *) item1;
  const ChunkEndLink * cel2 = (const ChunkEndLink *) item2;

 if(cel1->key.id == cel2->key.id)
 {
   if(cel1->key.id == cel2->key.id)
     return 0;
   else
     return (cel1->key.end == A_END) ? -1 : 1;
 }
 else
   return cel1->key.id - cel2->key.id;
}

                                             
void ProcessLinkForChimeraDetection(ScaffoldGraphT * graph,
                                    CIFragT * frag,
                                    CIFragT * mfrag,
                                    OrientType orient,
                                    CDS_CID_t distID,
                                    ChunkInstanceT * ci,
                                    MultiAlignT * ma,
                                    HashTable_AS * ht,
                                    VA_TYPE(ChunkEndLink) * cels,
                                    VA_TYPE(cds_uint16) * gcc,
                                    VA_TYPE(cds_uint16) * bcc)
{
  if(mfrag->contigID == frag->contigID)
  {
    AddLinkToMaps(graph, gcc, bcc, frag, mfrag,
                  orient, distID, GetMultiAlignLength(ma), FALSE);
  }
  else
  {
    // ChunkEndLink
    ChunkEndLink cel;
    ChunkEndLink * celp;
    cel.key.id = mfrag->contigID;
    cel.key.end = (frag->contigOffset3p.mean <
                   frag->contigOffset5p.mean) ? A_END : B_END;
    celp = LookupInHashTable_AS(ht, (void *) &(cel.key), sizeof(CELKey));
    if(celp != NULL)
    {
      // update
      celp->numEdges++;
      celp->leftCoord = min(celp->leftCoord,
                            frag->contigOffset3p.mean);
      celp->rightCoord = max(celp->rightCoord,
                             frag->contigOffset3p.mean);
    }
    else
    {
      cel.numEdges = 1;
      cel.leftCoord = frag->contigOffset3p.mean;
      cel.rightCoord = frag->contigOffset3p.mean;
      AppendVA_ChunkEndLink(cels, &cel);
      celp = GetVA_ChunkEndLink(cels, GetNumVA_ChunkEndLink(cels) - 1);
      InsertInHashTable_AS(ht,
                           (void *) &(celp->key),
                           sizeof(CELKey),
                           (void *) celp);
    }
  }
}


static int CELCompare(const ChunkEndLink * a, const ChunkEndLink * b)
{
  if(a->numEdges >= 2 && b->numEdges >= 2)
  {
    if(a->leftCoord == b->leftCoord)
      return a->rightCoord - b->rightCoord;
    else
      return a->leftCoord - b->leftCoord;
  }
  else
  {
    if(a->numEdges < 2 && b->numEdges < 2)
      return 0;
    else
      return (b->numEdges < 2) ? -1 : 1;
  }
}


void PrintChunkEndLinks(CDS_CID_t id,
                        VA_TYPE(ChunkEndLink) * cels,
                        FILE * fp)
{
  int i;
  fprintf(fp, F_CID " chunk end links:\n", id);
  for(i = 0; i < GetNumVA_ChunkEndLink(cels); i++)
  {
    ChunkEndLink * celp = GetVA_ChunkEndLink(cels, i);
    if(celp->numEdges >= 2)
      fprintf(fp, "\t" F_CID "-%c: (" F_COORD "," F_COORD "), %d edges\n",
              celp->key.id,
              (celp->key.end == A_END) ? 'A' : 'B',
              celp->leftCoord, celp->rightCoord, celp->numEdges);
  }
}

VA_TYPE(SplitInterval) * DetectChimericChunksInGraph(ScaffoldGraphT * graph)
{
  VA_TYPE(SplitInterval) * sis;
  int numCIs;
  VA_TYPE(cds_uint16) * rc;   // read coverage
  VA_TYPE(cds_uint16) * gcc;  // good clone coverage
  VA_TYPE(cds_uint16) * bcc;  // bad clone coverage
  HashTable_AS * ht;
  VA_TYPE(ChunkEndLink) * cels;
  int i;
  ChunkEndLink * celp;
  CDS_COORD_t maxRightCoord;
  int num1A = 0;
  int num1B = 0;
  
  if(graph->RezGraph == NULL ||
     (numCIs = GetNumGraphNodes(graph->RezGraph)) == 0)
  {
    fprintf(GlobalData->stderrc,
            "No contigs in graph for chimera detection.\n");
    return NULL;
  }
  rc = CreateVA_cds_uint16(10000);
  gcc = CreateVA_cds_uint16(10000);
  bcc = CreateVA_cds_uint16(10000);
  sis = CreateVA_SplitInterval(100);
  ht = CreateHashTable_AS(numCIs * 2, CELHashFn, CELCompareFn);
  cels = CreateVA_ChunkEndLink(numCIs * 2);
  
  assert(cels != NULL &&
         ht != NULL &&
         sis != NULL &&
         gcc != NULL &&
         bcc != NULL &&
         rc != NULL);

  // loop over contigs to find chimeric ones
  for(i = 0; i < numCIs; i++)
  {
    ChunkInstanceT * ci = GetGraphNode(graph->RezGraph, i);
    MultiAlignT * ma;
    int fIndex;
    int cIndex;

    if(ci->flags.bits.isDead)
      continue;

/*
    if(i != 0 && i % 10 == 0)
      fprintf(stderr, F_CID "\r", i);
*/
 
    ma = LoadMultiAlignTFromSequenceDB(graph->sequenceDB, ci->id, FALSE);
    
    if(GetMultiAlignLength(ma) < 2000)
      continue;

    // create gcc & bcc maps for unitig
    ResetVA_cds_uint16(rc);
    EnableRangeVA_cds_uint16(rc, GetMultiAlignLength(ma));
    ResetVA_cds_uint16(gcc);
    EnableRangeVA_cds_uint16(gcc, GetMultiAlignLength(ma));
    ResetVA_cds_uint16(bcc);
    EnableRangeVA_cds_uint16(bcc, GetMultiAlignLength(ma));
    ResetVA_ChunkEndLink(cels);
    ResetHashTable_AS(ht);
    
    // iterate over fragments to find inter-contig mates
    for(fIndex = 0; fIndex < GetNumIntMultiPoss(ma->f_list); fIndex++)
    {
      IntMultiPos * imp = GetIntMultiPos(ma->f_list, fIndex);
      InfoByIID * info = GetInfoByIID(graph->iidToFragIndex,
                                      imp->ident);
      CIFragT   * frag = GetCIFragT(graph->CIFrags,
                                    info->fragIndex);

      assert(frag->contigID == ci->id);
      if(frag->numLinks > 0)
      {
        if(frag->flags.bits.getLinksFromStore)
        {
          // going to the gatekeeper store should be the exception
          GateKeeperLinkRecordIterator gkpli;
          GateKeeperLinkRecord gkpl;
          CDS_CID_t linkHead = frag->linkHead;
          
          CreateGateKeeperLinkRecordIterator(graph->gkpStore.lnkStore,
                                             linkHead,
                                             frag->iid,
                                             &gkpli);
          while(NextGateKeeperLinkRecordIterator(&gkpli, &gkpl))
          {
            if(gkpl.type != AS_REREAD)
            {
              InfoByIID * info = GetInfoByIID(graph->iidToFragIndex,
                                              ((frag->iid == gkpl.frag1) ?
                                               gkpl.frag2 : gkpl.frag1));
              ProcessLinkForChimeraDetection(graph, frag, 
                                             GetCIFragT(graph->CIFrags,
                                                        info->fragIndex),
                                             (gkpl.orientation==AS_GKP_INNIE) ?
                                              AS_INNIE : AS_OUTTIE,
                                             (CDS_CID_t) gkpl.distance,
                                             ci, ma, ht, cels, gcc, bcc);
            }
          }
        }
        else
        {
          // get lone mate fragment
          if(frag->mateOf != NULLINDEX && frag->linkType != AS_REREAD)
          {
            ProcessLinkForChimeraDetection(graph, frag,
                                           GetCIFragT(graph->CIFrags,
                                                      frag->mateOf),
                                           (frag->flags.bits.innieMate) ?
                                           AS_INNIE : AS_OUTTIE,
                                           frag->dist,
                                           ci, ma, ht, cels, gcc, bcc);
          }
        }
      }
    }

    // process the cels array - sort, step through
    if(GetNumVA_ChunkEndLink(cels) < 2)
      continue;

    qsort(GetVA_ChunkEndLink(cels, 0),
          GetNumVA_ChunkEndLink(cels),
          sizeof(ChunkEndLink),
          (int (*) (const void *, const void *)) CELCompare);

    /*
    if(GetMultiAlignLength(ma) > 10000)
      PrintChunkEndLinks(ci->id, cels, GlobalData->stderrc);
    */
    
    celp = GetVA_ChunkEndLink(cels, 0);
    maxRightCoord = celp->rightCoord;
    if(celp->key.end == A_END)
    {
      num1A = 1;
      num1B = 0;
    }
    else
    {
      num1A = 0;
      num1B = 1;
    }

    // iterate left to right through chunk end links
    for(cIndex = 1; cIndex < GetNumVA_ChunkEndLink(cels); cIndex++)
    {
      celp = GetVA_ChunkEndLink(cels, cIndex);
      if(celp->numEdges < 2)
        break;
      if(num1A >= CHIMERA_SET_THRESHOLD &&
         num1B >= CHIMERA_SET_THRESHOLD &&
         celp->leftCoord > maxRightCoord)
      {
        cds_int32 c2Index;
        cds_int32 num2A = 0;
        cds_int32 num2B = 0;

        for(c2Index = cIndex;
            c2Index < GetNumVA_ChunkEndLink(cels);
            c2Index++)
        {
          ChunkEndLink * cel2p = GetVA_ChunkEndLink(cels, c2Index);
          if(cel2p->numEdges < 2)
            break;
          if(cel2p->key.end == A_END)
            num2A++;
          else
            num2B++;
          if(num2A >= CHIMERA_SET_THRESHOLD &&
             num2B >= CHIMERA_SET_THRESHOLD)
          {
            // maybe append an interval...
            // CDS_COORD_t k;
            SeqInterval interval;
            int inInterval = FALSE;
            
            CreateReadCoverageMap(graph, rc, ma, TRUE);

            interval.bgn = maxRightCoord;
            interval.end = celp->leftCoord;
            {
              SplitInterval si;
              
              PrintChunkEndLinks(ci->id, cels, GlobalData->stderrc);
              WriteMapToCelagramFile(rc,
                                     (CDS_COORD_t) GetMultiAlignLength(ma),
                                     "rc",
                                     "Read coverage for unitig", i);
              WriteMapToCelagramFile(gcc,
                                     (CDS_COORD_t) GetMultiAlignLength(ma),
                                     "gcc",
                                     "Good clone coverage for unitig", i);
              WriteMapToCelagramFile(bcc,
                                     (CDS_COORD_t) GetMultiAlignLength(ma),
                                     "bcc",
                                     "Bad clone coverage for unitig", i);
              
              si.id = ci->id;
              si.interval = interval;
              AppendVA_SplitInterval(sis, &si);
              fprintf(GlobalData->stderrc,
                      "Chimeric interval in contig " F_CID ": (" F_COORD "," F_COORD ")\n",
                      si.id, si.interval.bgn, si.interval.end);
              inInterval = FALSE;
              break;
            }

            /*
            // check against good clone coverage
            for(k = maxRightCoord; k < celp->leftCoord; k++)
            {
              if(*(GetVA_cds_uint16(rc,k)) <= MAX_SEQUENCE_COVERAGE &&
                 *(GetVA_cds_uint16(gcc, k)) <= MAX_GOOD_CLONE_COVERAGE)
              {
                if(inInterval)
                {
                  interval.end = k;
                }
                else
                {
                  interval.bgn = interval.end = k;
                  inInterval = TRUE;
                }
              }
              else
              {
                if(inInterval)
                {
                  SplitInterval si;
                  
                  WriteMapToCelagramFile(rc,
                                         (CDS_COORD_t) GetMultiAlignLength(ma),
                                         "rc",
                                         "Read coverage for unitig", i);
                  WriteMapToCelagramFile(gcc,
                                         (CDS_COORD_t) GetMultiAlignLength(ma),
                                         "gcc",
                                         "Good clone coverage for unitig", i);
                  WriteMapToCelagramFile(bcc,
                                         (CDS_COORD_t) GetMultiAlignLength(ma),
                                         "bcc",
                                         "Bad clone coverage for unitig", i);
          
                  si.id = ci->id;
                  si.interval = interval;
                  AppendVA_SplitInterval(sis, &si);
                  fprintf(GlobalData->stderrc,
                          "Chimeric interval in contig " F_CID ": (" F_COORD "," F_COORD ")\n",
                          si.id, si.interval.bgn, si.interval.end);
                  inInterval = FALSE;
                }
              }
            }
            */
          }
        }
      }
      
      if(celp->key.end == A_END)
        num1A++;
      else
        num1B++;
      
      // always adjust this, so that we don't get overlapping intervals
      maxRightCoord = max(maxRightCoord, celp->rightCoord);
    }
  }

  DeleteVA_cds_uint16(rc);
  DeleteVA_cds_uint16(gcc);
  DeleteVA_cds_uint16(bcc);
  
  return sis;
}

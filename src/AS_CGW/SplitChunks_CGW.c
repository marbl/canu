
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

static char *rcsid = "$Id: SplitChunks_CGW.c,v 1.34 2008-11-07 06:13:55 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "MultiAlignment_CNS.h"
#include "ScaffoldGraph_CGW.h"
#include "Input_CGW.h"

//  Initial splitting of unitigs is performed as follows:
//
//  For all, shortDiscriminatorUnique(+), or just DiscriminatorUnique
//  unitigs?  Only for unitigs >= min distance mean + CGW_CUTOFF *
//  stddev in length
//
//  1. create sequence coverage map of reads (AS_READ | AS_EXTR |
//  AS_TRNR?) - trim 30bp from each end of each read in creating map
//  if, beyond initial 1x coverage & until last 1x coverage, there is
//  at least one location with <=1x coverage then continue
//
//  2. create good & bad clone coverage maps.  if at some point in seq
//  coverage map where coverage <= 1 (not on ends) gcc <= 1 && bcc >=
//  2 split at that point/interval.  - remember to scan entire unitig
//  for possibility of multiple break points
//
//  Splitting - make sure it will work for contigs & unitigs and for
//  contig calling fn for unitig:
//
//  Inputs: chunk index, interval where to break, graph, chunk type
//
//  Outputs: if unitig, need to report info for contig breaking
//
//  FOR UNITIGS
//
//  1. original unitig becomes 'right' unitig.  add 2 new unitigs :
//  chimeric middle & left
//
//  2. from left to right, populate new unitigs/depopulate original
//  until interval reached: frags go to left unitig
//  until interval left: frags go to chimeric middle unitig
//
//  3. finish unitig creation - A stat, consensus, unitig edges
//
//  take right unitig & repeat checks for good/bad/sequence coverage &
//  possibly split again


#define READ_TRIM_BASES         30
#define MAX_SEQUENCE_COVERAGE    1
#define MIN_BAD_CLONE_COVERAGE   3
#define MAX_GOOD_CLONE_COVERAGE  0
#define LN_2                     0.693147

#define SOURCE_LENGTH            100

typedef struct {
  CDS_CID_t   id;
  SeqInterval interval;
} SplitInterval;

VA_DEF(SplitInterval);
VA_DEF(uint16);
VA_DEF(SeqInterval);




typedef struct {
  IntUnitigMesg    ium;

  VA_TYPE(int32)  *deltas;
  VA_TYPE(char)   *sequence;
  VA_TYPE(char)   *quality;

  CDS_COORD_t      minPos;
  int32            numFragsAllocated;
  int32            numRandomFragments;
} IUMStruct;


static
IUMStruct *
CreateIUMStruct(void) {
  IUMStruct *is = (IUMStruct *) safe_calloc(1, sizeof(IUMStruct));

  is->deltas    = CreateVA_int32(1);
  is->sequence  = CreateVA_char(200000);
  is->quality   = CreateVA_char(200000);

  is->minPos    = CDS_COORD_MAX;

  return is;
}


static
void
ResetIUMStruct(IUMStruct *is) {
  int i;

  for(i = 0; i < is->ium.num_frags; i++) {
    is->ium.f_list[i].delta_length = 0;
    is->ium.f_list[i].delta = NULL;
  }

  is->ium.consensus      = NULL;
  is->ium.quality        = NULL;
  is->ium.num_frags      = 0;

  is->minPos             = CDS_COORD_MAX;
  is->numRandomFragments = 0;
}


static
void
FreeIUMStruct(IUMStruct *is) {

  if (is == NULL)
    return;

  DeleteVA_int32(is->deltas);
  DeleteVA_char(is->sequence);
  DeleteVA_char(is->quality);

  safe_free(is->ium.f_list);

  safe_free(is);
}


static
void
AddIMPToIUMStruct(IUMStruct *is, IntMultiPos *imp) {
  InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, imp->ident);
  CIFragT   *frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

  if (is->ium.num_frags >= is->numFragsAllocated) {
    is->numFragsAllocated += 10;
    is->ium.f_list         = (IntMultiPos *)safe_realloc(is->ium.f_list, is->numFragsAllocated * sizeof(IntMultiPos));
  }

  is->ium.f_list[is->ium.num_frags].type = imp->type;
  is->ium.f_list[is->ium.num_frags].ident = imp->ident;
  is->ium.f_list[is->ium.num_frags].contained = imp->contained;
  is->ium.f_list[is->ium.num_frags].position.bgn = frag->offset5p.mean;
  is->ium.f_list[is->ium.num_frags].position.end = frag->offset3p.mean;
  is->ium.f_list[is->ium.num_frags].delta_length = 0;
  is->ium.f_list[is->ium.num_frags].delta = NULL;

  is->ium.num_frags++;

  is->minPos = MIN(is->minPos,MIN(frag->offset5p.mean,frag->offset3p.mean));
  is->numRandomFragments += (AS_FA_RANDOM(imp->type)) ? 1 : 0;
}


static
void
PrintPositions(ScaffoldGraphT *graph, MultiAlignT *ma, int isUnitig, char *label) {
  int i;

  fprintf(stderr, "%10" F_IIDP ": %s\n", ma->maID, label);

  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp  = GetIntMultiPos(ma->f_list, i);
    InfoByIID   *info = GetInfoByIID(graph->iidToFragIndex, imp->ident);
    CIFragT     *frag = GetCIFragT(graph->CIFrags, info->fragIndex);

    fprintf(stderr, F_IID ": (" F_COORD "," F_COORD ")\t(%f,%f)\t"F_IID"\n",
            imp->ident,
            MIN(imp->position.bgn, imp->position.end),
            MAX(imp->position.bgn, imp->position.end),
            (isUnitig) ? MIN(frag->offset5p.mean, frag->offset3p.mean) : MIN(frag->contigOffset5p.mean, frag->contigOffset3p.mean),
            (isUnitig) ? MAX(frag->offset5p.mean, frag->offset3p.mean) : MAX(frag->contigOffset5p.mean, frag->contigOffset3p.mean),
            imp->contained);
  }
}


static
void
IncrementMapInterval(VA_TYPE(uint16) *map, CDS_COORD_t minPos, CDS_COORD_t maxPos) {
  CDS_COORD_t i;
  uint16     *t = GetVA_uint16(map, 0);

  for (i=minPos; i<maxPos; i++)
    t[i]++;
}


static
void
AddLinkToMaps(ScaffoldGraphT *graph,
              VA_TYPE(uint16) *gcc,
              VA_TYPE(uint16) *bcc,
              CIFragT *frag,
              CIFragT *mfrag,
              OrientType orient,
              CDS_CID_t distID,
              CDS_COORD_t length,
              int isUnitig) {

  DistT *dist = GetDistT(graph->Dists, distID);

  CDS_COORD_t minPos = ((isUnitig) ?
                        (MIN(frag->offset5p.mean, mfrag->offset5p.mean)) :
                        (MIN(frag->contigOffset5p.mean, mfrag->contigOffset5p.mean)));
  CDS_COORD_t maxPos = ((isUnitig) ?
                        (MAX(frag->offset5p.mean, mfrag->offset5p.mean)) :
                        (MAX(frag->contigOffset5p.mean, mfrag->contigOffset5p.mean)));

  // if they're in the same unitig
  if((isUnitig && frag->cid == mfrag->cid) ||
     (!isUnitig && frag->contigID == mfrag->contigID)) {
    // for pairs in the same unitig, just process the lesser
    if(frag->iid < mfrag->iid) {
      // if orientation is the same, the pair is bad
      if((isUnitig &&
          getCIFragOrient(frag) == getCIFragOrient(mfrag)) ||
         (!isUnitig &&
          GetContigFragOrient(frag) == GetContigFragOrient(mfrag))) {
        // bad pair - increment intervals
        if(orient == AS_INNIE) {
          // link is innie
          if((isUnitig && getCIFragOrient(frag) == A_B) ||
             (!isUnitig && GetContigFragOrient(frag) == A_B)) {
            // both are A_B oriented, so intervals > 5p are suspect
            // Increment from minPos to somewhere & maxPos to somewhere
            // minPos                   maxPos
            //     |                         |
            //     -------->                 -------->
            //     --- increment ---?        --- increment ---?
            //
            // note that maxPos < length
            IncrementMapInterval(bcc,
                                 minPos,
                                 MIN(minPos + dist->mu +
                                     CGW_CUTOFF * dist->sigma,
                                     maxPos));
            IncrementMapInterval(bcc,
                                 maxPos,
                                 MIN(maxPos + dist->mu +
                                     CGW_CUTOFF * dist->sigma,
                                     length - 1));
          } else {
            // both are B_A oriented, so intervals < 5p are suspect
            // Increment from minPos to somewhere & maxPos to somewhere
            //              minPos                   maxPos
            //                  |                         |
            //          <--------                 <--------
            // ?--- increment ---        ?--- increment ---
            //
            IncrementMapInterval(bcc,
                                 MAX(0,
                                     minPos - dist->mu -
                                     CGW_CUTOFF * dist->sigma),
                                 minPos);
            IncrementMapInterval(bcc,
                                 MAX(minPos,
                                     maxPos - dist->mu -
                                     CGW_CUTOFF * dist->sigma),
                                 maxPos);
          }
        } else {
          // link is outtie
          if((isUnitig && getCIFragOrient(frag) == A_B) ||
             (!isUnitig && GetContigFragOrient(frag) == A_B)) {
            // both are A_B oriented, so
            // Increment from minPos to somewhere & maxPos to somewhere
            //               minPos                    maxPos
            //                   |                         |
            //                   -------->                 -------->
            // ?--- increment ---        ?--- increment ---
            //
            IncrementMapInterval(bcc,
                                 MAX(0,
                                     minPos -
                                     dist->mu - CGW_CUTOFF * dist->sigma),
                                 minPos);
            IncrementMapInterval(bcc,
                                 MAX(minPos,
                                     maxPos -
                                     dist->mu - CGW_CUTOFF * dist->sigma),
                                 maxPos);
          } else {
            // both are B_A oriented, so
            // Increment from minPos to somewhere & maxPos to somewhere
            //     minPos                       maxPos
            //         |                            |
            // <--------                    <--------
            //          --- increment ---?           --- increment ---?
            //
            IncrementMapInterval(bcc,
                                 minPos,
                                 MIN(minPos +
                                     dist->mu + CGW_CUTOFF * dist->sigma,
                                     maxPos));
            IncrementMapInterval(bcc,
                                 maxPos,
                                 MIN(length - 1,
                                     maxPos +
                                     dist->mu + CGW_CUTOFF * dist->sigma));
          }
        }
      } else {
        // fragments are oriented differently - may be okay
        CDS_COORD_t distance = maxPos - minPos;

        // if the distance is wrong, the pair is bad
        if(distance < dist->mu - CGW_CUTOFF * dist->sigma ||
           distance > dist->mu + CGW_CUTOFF * dist->sigma) {
          // bad pair
          // same intervals for either innie or outtie
          //
          //  innie:
          //  minPos                                 maxPos
          //      |                                      |
          //      ------->                       <--------
          //      --- increment ---?    ?--- increment ---
          //
          //  outtie:
          //     minPos                                 maxPos
          //         |                                      |
          //  <-------                                      -------->
          //          --- increment ---?  ?--- increment ---
          //
          IncrementMapInterval(bcc,
                               minPos,
                               MIN(maxPos,
                                   minPos +
                                   dist->mu + CGW_CUTOFF * dist->sigma));
          IncrementMapInterval(bcc,
                               MAX(minPos +
                                   dist->mu + CGW_CUTOFF * dist->sigma,
                                   maxPos -
                                   dist->mu - CGW_CUTOFF * dist->sigma),
                               maxPos);
        } else {
          // good pair
          IncrementMapInterval(gcc, minPos, maxPos);
        }
      }
    }
  } else {
    // in different unitigs, should be close to end of unitig
    if(orient == AS_INNIE) {
      if((isUnitig && getCIFragOrient(frag) == A_B &&
          frag->offset5p.mean <
          length - dist->mu - CGW_CUTOFF * dist->sigma) ||
         (!isUnitig && GetContigFragOrient(frag) == A_B &&
          frag->contigOffset5p.mean <
          length - dist->mu - CGW_CUTOFF * dist->sigma)) {
        //
        //  --------->
        //  --- increment ---
        //
        IncrementMapInterval(bcc,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean),
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean) +
                             dist->mu + CGW_CUTOFF * dist->sigma);
      } else if((isUnitig && getCIFragOrient(frag) == B_A &&
                 frag->offset5p.mean >
                 dist->mu + CGW_CUTOFF * dist->sigma) ||
                (!isUnitig && GetContigFragOrient(frag) == B_A &&
                 frag->contigOffset5p.mean >
                 dist->mu + CGW_CUTOFF * dist->sigma)) {
        //
        //        <----------
        //  --- increment ---
        //
        IncrementMapInterval(bcc,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean) -
                             dist->mu - CGW_CUTOFF * dist->sigma,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean));
      }
    } else {
      // outtie
      if((isUnitig && getCIFragOrient(frag) == B_A &&
          frag->offset5p.mean <
          length - dist->mu - CGW_CUTOFF * dist->sigma) ||
         (!isUnitig && GetContigFragOrient(frag) == B_A &&
          frag->contigOffset5p.mean <
          length - dist->mu - CGW_CUTOFF * dist->sigma)) {
        //
        //  <----------
        //             --- increment ---
        //
        IncrementMapInterval(bcc,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean),
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean) +
                             dist->mu + CGW_CUTOFF * dist->sigma);

      } else if((isUnitig && getCIFragOrient(frag) == A_B &&
                 frag->offset5p.mean >
                 dist->mu + CGW_CUTOFF * dist->sigma) ||
                (!isUnitig && GetContigFragOrient(frag) == A_B &&
                 frag->contigOffset5p.mean >
                 dist->mu + CGW_CUTOFF * dist->sigma)) {
        //
        //                   --------->
        //  --- increment ---
        //
        IncrementMapInterval(bcc,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean) -
                             dist->mu - CGW_CUTOFF * dist->sigma,
                             ((isUnitig) ?
                              frag->offset5p.mean :
                              frag->contigOffset5p.mean));
      }
    }
  }
}


static
void
CreateReadCoverageMap(ScaffoldGraphT *graph,
                      VA_TYPE(uint16) *rc,
                      MultiAlignT *ma,
                      int isUnitig) {
  uint32 i;

  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    // only add 'reads' to sequence coverage map

    if (AS_FA_READ(imp->type)) {
      InfoByIID *info = GetInfoByIID(graph->iidToFragIndex, imp->ident);
      CIFragT   *frag = GetCIFragT(graph->CIFrags, info->fragIndex);

      // if unitig is forward vs. reverse, add to correct fragment location

      CDS_COORD_t minPos = ((isUnitig) ?
                            (MIN(frag->offset5p.mean, frag->offset3p.mean) + READ_TRIM_BASES) :
                            (MIN(frag->contigOffset5p.mean, frag->contigOffset3p.mean) + READ_TRIM_BASES));

      CDS_COORD_t maxPos = ((isUnitig) ?
                            (MAX(frag->offset5p.mean, frag->offset3p.mean) - READ_TRIM_BASES) :
                            (MAX(frag->contigOffset5p.mean, frag->contigOffset3p.mean) - READ_TRIM_BASES));

      IncrementMapInterval(rc, minPos, maxPos);
    }
  }
}


static
void
CreateCloneCoverageMaps(ScaffoldGraphT *graph,
                        VA_TYPE(uint16) *gcc,
                        VA_TYPE(uint16) *bcc,
                        MultiAlignT *ma,
                        int isUnitig) {
  uint32 i;

  // loop over fragments in unitig
  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp  = GetIntMultiPos(ma->f_list, i);
    InfoByIID   *info = GetInfoByIID(graph->iidToFragIndex, imp->ident);
    CIFragT     *frag = GetCIFragT(graph->CIFrags, info->fragIndex);

    // this is unlike ComputeMatePairStatistics, since we're interested
    // even in pairs that are in different unitigs/contigs

    // get lone mate fragment
    if ((frag->flags.bits.hasMate > 0) &&
        (frag->mateOf != NULLINDEX) &&
        (frag->flags.bits.linkType != AS_REREAD))
      AddLinkToMaps(graph, gcc, bcc, frag,
                    GetCIFragT(graph->CIFrags, frag->mateOf),
                    ((frag->flags.bits.innieMate) ? AS_INNIE : AS_OUTTIE),
                    frag->dist,
                    GetMultiAlignLength(ma),
                    isUnitig);
  }
}







static
int
positionCompare(const void *A, const void *B) {
  const IntMultiPos *a = (const IntMultiPos *)A;
  const IntMultiPos *b = (const IntMultiPos *)B;
  return(MIN(a->position.bgn, a->position.end) - MIN(b->position.bgn, b->position.end));
}

//  There are two methods for fixing containment ordering.
//
//  AdjustContainedOrder() is very aggressive, moving all contained
//  fragments to be immediately after their container.  This fact is
//  used in SplitChunkByIntervals()
//
//  FixContainedOrder() just makes sure that all containers are after
//  their containee, and does so without changing the position-based
//  sorting.

static
void
AdjustContainedOrder(IntMultiPos *impList, int numIMPs) {
  HashTable_AS *containing;
  IntMultiPos *contained;
  int numContained;
  int numContaining;
  int numMoved;
  int totalMoved;
  int numContainedLeft;
  int index;

  // Keep contained fragments just after containing fragment
  //
  // 1. Remove contained imps from impList
  //
  // 2. Append them to the end of impList so that even contained imp A
  // that contains another imp B will come before B
  //
  // 3. Move contained imps up through impList to position them just
  // after their containing imp
  //
  // How?
  //
  // 1. copy contained imps into its own array.  while compacting
  // non-contained imps in impList & adding non-contained imps to a
  // hashtable
  //
  // 2. looping until no more contained imps get moved to the end of
  // impList, make sure contained imps that somehow don't have
  // containing imps in the impList get moved to the end, too (orphans).
  //
  // 3. Move contained imps up to the position after their containing
  // imp
  //
  // 4. Move orphaned contained imps up to a position where they should
  // overlap the preceding imp
  //

  containing = CreateScalarHashTable_AS(numIMPs);
  contained = (IntMultiPos *) safe_malloc(numIMPs * sizeof(IntMultiPos));
  assert(containing != NULL && contained != NULL);

  for(index = 0, numContained = 0, numContaining = 0; index < numIMPs; index++) {
    if(impList[index].contained != 0) {
      // copy to contained array
      contained[numContained++] = impList[index];
    } else {
      // move up and add to hashtable
      impList[numContaining] = impList[index];
      InsertInHashTable_AS(containing, impList[numContaining].ident, 0, 0, 0);
      numContaining++;
    }
  }

  //
  // Until all IMPs have been moved
  //
  // Iterate over contained IMP list
  //
  // If imp status is Waiting, look up containing IMP in hashtable if
  // it's in the hashtable, change to MarkedToMove
  //
  // Iterate over contained IMP list
  //
  // If imp status is MarkedToMove, insert in impList just after
  // containing imp and change status to AlreadyMoved & add to
  // hashtable
  //
  totalMoved = 0;
  numMoved = 100;
  numContainedLeft = numContained;
  while(numMoved > 0) {
    numMoved = 0;
    // identify contained imps whose containing imp is in impList
    for(index = 0; index < numContainedLeft; index++) {
      // if this imp hasn't moved & it's containing imp is in the hashtable
      // then move it to the end of the impList & put it in the hashtable
      if(ExistsInHashTable_AS(containing, contained[index].contained, 0)) {
        impList[numContaining + totalMoved] = contained[index];
        InsertInHashTable_AS(containing, impList[numContaining + totalMoved].ident, 0, 0, 0);
        numMoved++;
        totalMoved++;
      } else {
        // move toward the head of the array
        contained[index - numMoved] = contained[index];
      }
    }
    numContainedLeft -= numMoved;
  }
  DeleteHashTable_AS(containing);

  safe_free(contained);

  assert(numContainedLeft == 0);

  // Now move containeds & orphaned containeds up to behind their containing
  for(index = numContaining; index < numIMPs; index++) {
    int index2;
    IntMultiPos imp = impList[index];

    // move up to behind its containing
    if(impList[index].contained != 0) {
      for(index2 = index - 1;
          impList[index2].ident != imp.contained &&
            impList[index2].contained != imp.contained;
          index2--) {
        impList[index2 + 1] = impList[index2];
      }
    } else {
      // move up to behind something it should overlap
      for(index2 = index - 1;
          index2 >= 0 &&
            (MIN(impList[index2].position.bgn, impList[index2].position.end) >
             MIN(imp.position.bgn, imp.position.end) ||
             MAX(impList[index2].position.bgn, impList[index2].position.end) <
             MIN(imp.position.bgn, imp.position.end));
          index2--) {
        impList[index2 + 1] = impList[index2];
      }
    }
    impList[index2 + 1] = imp;
  }
}



static
void
FixContainedOrder(IntMultiPos *impList, int numIMPs) {
  int i=0, j=0;

#if 0
  VERBOSE_MULTIALIGN_OUTPUT = 1;
#endif

#if 0
  fprintf(stderr, "BEFORE\n");
  for (i=0; i<numIMPs; i++)
    fprintf(stderr, "%d\t"F_IID"\t%d\t%d\t"F_IID"\n",
            i,
            impList[i].ident,
            impList[i].position.bgn,
            impList[i].position.end,
            impList[i].contained);
#endif

  while (i < numIMPs - 1) {
    j = i + 1;

#if 0
    fprintf(stderr, "TEST %d\t"F_IID"\t%d\t%d\t"F_IID"\n",
            i,
            impList[i].ident,
            impList[i].position.bgn,
            impList[i].position.end,
            impList[i].contained);
#endif

    int  ipos = MIN(impList[i].position.bgn, impList[i].position.end);
    int  jpos = MIN(impList[j].position.bgn, impList[j].position.end);

    while ((ipos == jpos) && (j < numIMPs)) {
      j++;
      jpos = MIN(impList[j].position.bgn, impList[j].position.end);
    }

    if (i+1 != j) {

      //  Got a set of reads that all start at the same position.  We
      //  don't expect this to be a very large set, and we do lots of
      //  linear searches.

#if 0
      {
        int xx;
        fprintf(stderr, "Found potential containments to fix.\n");
        for (xx=i; xx<j; xx++)
          fprintf(stderr, "%d\t"F_IID"\t%d\t%d\t"F_IID"\n",
                  xx,
                  impList[xx].ident,
                  impList[xx].position.bgn,
                  impList[xx].position.end,
                  impList[xx].contained);
      }
#endif

      int  beg = i;
      int  elt = i;

      //  Move all the non-contained reads to the start.  Scan the
      //  list (variable 'elt') move any uncontained reads to position
      //  'beg'.  Advance 'beg' to the next uncontained read.
      //
      while ((beg < j) && (impList[beg].contained == 0))
        beg++;
      for (elt = beg+1; elt < j; elt++) {
        if (impList[elt].contained == 0) {
#if 0
          fprintf(stderr, "swap 1 container %d to %d\n", elt, beg);
#endif
          IntMultiPos t = impList[beg];
          impList[beg]  = impList[elt];
          impList[elt]  = t;
          beg++;
        }
      }

      //  beg == first read that is contained
      //  elt == read we are examining

      elt = j - 1;

      //  Starting with the last read in the i-j range, check that
      //  it's container is present in the set of properly ordered
      //  containers [i...beg].  If it is, move that read to the set
      //  of properly ordered containers.  If its container is in the
      //  range [beg...elt] we need to first move the container into
      //  the set of properly ordered containers.  And if its
      //  container isn't found in the range [i...elt] then the
      //  container must be before our range.
      //
      while (beg < elt) {
        int scan = 0;
        int swap = 0;

        assert(impList[elt].contained != 0);

        for (scan=i; scan<elt; scan++) {

          //  If we found the container, move elements and stop scanning.
          //
          if (impList[scan].ident == impList[elt].contained) {
            if (scan < beg) {
              //  Container is properly ordered.  Move 'elt' into the
              //  set of containers.
#if 0
              fprintf(stderr, "swap 2 containee %d to %d\n", elt, beg);
#endif
              IntMultiPos t = impList[beg];
              impList[beg]  = impList[elt];
              impList[elt]  = t;
              swap = 1;
              beg++;
            } else {
              //  Our container is here, but it's not yet in the
              //  correct spot.  Swap the two elements, and restart
              //  the loop to process the container.
#if 0
              fprintf(stderr, "swap 3 container %d to %d\n", scan, elt);
#endif
              IntMultiPos t = impList[scan];
              impList[scan] = impList[elt];
              impList[elt]  = t;
              swap = 1;
            }
            scan = j;
          }  //  end of if container found
        }  //  end of scan loop

        //  Move to the next element.  If we swapped elements, this
        //  will move us back to the element we swapped in, and we'll
        //  next process it.  If we didn't swap elements, this element
        //  is contained in something not in the [i...j] range, and is thus
        //  eligible to be moved into the set of valid containers.
        //
        if (swap == 0) {
          IntMultiPos t = impList[beg];
          impList[beg]  = impList[elt];
          impList[elt]  = t;
          beg++;
        }
      }

#if 0
      {
        int xx;
        fprintf(stderr, "Final containments fixed?\n");
        for (xx=i; xx<j; xx++)
          fprintf(stderr, "%d\t"F_IID"\t%d\t%d\t"F_IID"\n",
                  xx,
                  impList[xx].ident,
                  impList[xx].position.bgn,
                  impList[xx].position.end,
                  impList[xx].contained);
      }
#endif

    }  //  Done fixing the containments in range [i...j]

    i = j;
  }  //  Move to the next range

#if 0
  fprintf(stderr, "AFTER\n");
  for (i=0; i<numIMPs; i++)
    fprintf(stderr, F_IID"\t%d\t%d\t"F_IID"\n",
            impList[i].ident,
            impList[i].position.bgn,
            impList[i].position.end,
            impList[i].contained);
#endif
}




static
void
StoreIUMStruct(ScaffoldGraphT *graph,
               IUMStruct *is,
               int isUnitig,
               UnitigStatus status,
               float egfar) {
  int i;
  int32 rho = 0;

  is->ium.iaccession = ((isUnitig) ? GetNumGraphNodes(graph->CIGraph) :
                        GetNumGraphNodes(graph->ContigGraph));

  is->ium.status = status;

  // re-position fragments & compute (estimated) coverage statistic
  is->ium.length = 0;
  for(i = 0; i < is->ium.num_frags; i++) {
    is->ium.f_list[i].position.bgn -= is->minPos;
    is->ium.f_list[i].position.end -= is->minPos;
    if(is->ium.f_list[i].position.end > is->ium.f_list[i].position.bgn) {
      rho = MAX(rho,is->ium.f_list[i].position.bgn);
      is->ium.length = MAX(is->ium.length, is->ium.f_list[i].position.end);
      is->ium.f_list[i].position.end--;
    } else {
      rho = MAX(rho,is->ium.f_list[i].position.end);
      is->ium.length = MAX(is->ium.length, is->ium.f_list[i].position.bgn);
      is->ium.f_list[i].position.bgn--;
    }
    // check that containing fragment is present, if specified
    // NOTE: ugly search
    if(is->ium.f_list[i].contained) {
      int j;
      for(j = 0; j < is->ium.num_frags; j++) {
        if(is->ium.f_list[j].ident == is->ium.f_list[i].contained)
          break;
      }
      if(j == is->ium.num_frags)
        is->ium.f_list[i].contained = 0;
    }
  }


  // egfar is the estimated global fragment arrival rate based on the original
  // chunk.
  //
  // rho is the sum of the a-hangs (= largest starting position)
  //
  // Based on compute_coverage_statistic() in AS_CGB_cgb.h
  //
  is->ium.coverage_stat = 0.0;
  if ((egfar > 0.f && is->numRandomFragments > 0))
    is->ium.coverage_stat = rho * egfar - LN_2 * (is->numRandomFragments - 1);

  is->ium.forced = 0;

  qsort(is->ium.f_list, is->ium.num_frags, sizeof(IntMultiPos), positionCompare );

  FixContainedOrder(is->ium.f_list, is->ium.num_frags);

  // get the multi-alignment - this fills in some unitig fields

  int unitigFail = MultiAlignUnitig(&(is->ium),
                                    ScaffoldGraph->gkpStore,
                                    is->sequence,
                                    is->quality,
                                    is->deltas,
                                    CNS_STATS_ONLY,
                                    0,
                                    Local_Overlap_AS_forCNS, // DP_Compare
                                    NULL);

  //  Whoops!  Failed!  Like consensus does in MultiAlignUnitig() we
  //  now (2007-12-07) will try again, allowing negative hangs.
  //
  if (unitigFail) {
    int  ov = VERBOSE_MULTIALIGN_OUTPUT, oh=allow_neg_hang;

    VERBOSE_MULTIALIGN_OUTPUT = 1;
    allow_neg_hang            = 1;
    unitigFail = MultiAlignUnitig(&(is->ium),
                                  ScaffoldGraph->gkpStore,
                                  is->sequence,
                                  is->quality,
                                  is->deltas,
                                  CNS_STATS_ONLY,
                                  0,
                                  Local_Overlap_AS_forCNS, // DP_Compare
                                  NULL);
    VERBOSE_MULTIALIGN_OUTPUT = ov;
    allow_neg_hang            = oh;
  }

  //  Whoops!  Failed!  Unlike consensus does, we'll also try to
  //  force the fragment (2007-12-07).
  //
#if 0
  //  Not sure if the forcing was doing anything appropriate.  It
  //  would try to get an overlap (and a trace), fail, but then
  //  proceed to apply the (bogus) trace to align the fragment.
  //  (2008-11-05)
  //
  if (unitigFail) {
    int  ov = VERBOSE_MULTIALIGN_OUTPUT, of=allow_forced_frags;

    VERBOSE_MULTIALIGN_OUTPUT = 1;
    allow_forced_frags        = 1;
    unitigFail = MultiAlignUnitig(&(is->ium),
                                  ScaffoldGraph->gkpStore,
                                  is->sequence,
                                  is->quality,
                                  is->deltas,
                                  CNS_STATS_ONLY,
                                  0,
                                  Local_Overlap_AS_forCNS, // DP_Compare
                                  NULL);
    VERBOSE_MULTIALIGN_OUTPUT = ov;
    allow_forced_frags        = of;
  }
#endif

  if (unitigFail) {
    fprintf(GlobalData->stderrc, "FATAL ERROR: MultiAlignUnitig call failed in unitig splitting.\n");
    assert(FALSE);
  }


  // add ium to the system
  {
    MultiAlignT *ma = CreateMultiAlignTFromIUM(&(is->ium), GetNumCIFragTs(graph->CIFrags), FALSE);
    CDS_COORD_t length = GetMultiAlignUngappedLength(ma);
    ChunkInstanceT *ci;

    // need to point fragments to their new unitig/contig

    for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++) {
      IntMultiPos *imp  = GetIntMultiPos(ma->f_list, i);
      InfoByIID   *info = GetInfoByIID(graph->iidToFragIndex, imp->ident);

      imp->sourceInt = info->fragIndex;
    }


    //  Insert both a unitig and a contig

    insertMultiAlignTInSequenceDB(graph->sequenceDB, is->ium.iaccession, TRUE,
                                  ma,
                                  TRUE);
    insertMultiAlignTInSequenceDB(graph->sequenceDB, is->ium.iaccession, FALSE,
                                  CopyMultiAlignT(NULL, ma),
                                  FALSE);

    ProcessIUM_ScaffoldGraph(&(is->ium), length, TRUE);

    ci = GetGraphNode((isUnitig ? graph->CIGraph : graph->ContigGraph), is->ium.iaccession);

    UpdateNodeFragments((isUnitig ? graph->CIGraph : graph->ContigGraph),
                        is->ium.iaccession,
                        ci->type == DISCRIMINATORUNIQUECHUNK_CGW,
                        TRUE);
  }
}






static
float
EstimateGlobalFragmentArrivalRate(ChunkInstanceT *ci, MultiAlignT *ma) {
  int    i, numRF=0, rho=0;

  // estimate random fragments/bp

  for (i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    if (AS_FA_RANDOM(imp->type)) {
      rho = MAX(rho, MIN(imp->position.bgn, imp->position.end));
      numRF++;
    }
  }

  return((ci->info.CI.coverageStat + LN_2 * (numRF - 1)) / rho);
}


// NOTE: csis intervals are in ungapped coordinates
static
void
SplitChunkByIntervals(ScaffoldGraphT *graph,
                      CDS_CID_t ciID,
                      MultiAlignT *ma,
                      VA_TYPE(SeqInterval) *csis,
                      int isUnitig) {
  int             currI = 0;
  SeqInterval    *currInterval = GetVA_SeqInterval(csis,currI);
  int             i;
  CDS_CID_t       lastIID = 0;
  IUMStruct      *lastPlacement = NULL;
  IUMStruct      *good = CreateIUMStruct();
  IUMStruct      *bad  = CreateIUMStruct();
  float           egfar;
  ChunkInstanceT *ci = GetGraphNode(graph->CIGraph, ciID);

  assert(isUnitig);  // not implemented for contigs yet

  qsort(ma->f_list->Elements, GetNumIntMultiPoss(ma->f_list), sizeof(IntMultiPos), positionCompare);

  egfar = EstimateGlobalFragmentArrivalRate(ci, ma);

  AdjustContainedOrder((IntMultiPos *) ma->f_list->Elements, GetNumIntMultiPoss(ma->f_list));

#if 1
  // feedback for log file
  fprintf(GlobalData->stderrc, "Splitting %s " F_CID " into as many as %d %s at intervals:",
          (isUnitig ? "unitig" : "contig"), ma->maID,
          (int) (2 * GetNumVA_SeqInterval(csis) + 1),
          (isUnitig ? "unitigs" : "contigs"));

  for(currI = 0; currI < GetNumVA_SeqInterval(csis); currI++) {
    currInterval = GetVA_SeqInterval(csis, currI);
    fprintf(GlobalData->stderrc, "\t" F_COORD "," F_COORD, currInterval->bgn, currInterval->end);
  }
  fprintf(GlobalData->stderrc, "\n");
#endif

  //
  // loop over fragments, placing them in different chunks
  //
  // Two chunks are populated at the same stage - a good & a bad -
  // until the last chimeric interval is passed, in which case only a
  // good is populated.
  //
  // So, when the current bad interval is exceeded, the next good/bad
  // pair must be used. However, intervals may be close together to
  // the point that the very next chimeric interval may not be the
  // furthest chimeric interval into which the fragment falls...
  //
  currI = 0;
  currInterval = GetVA_SeqInterval(csis, currI);
  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp  = GetIntMultiPos(ma->f_list, i);
    InfoByIID   *info = GetInfoByIID(graph->iidToFragIndex, imp->ident);
    CIFragT     *frag = GetCIFragT(graph->CIFrags, info->fragIndex);

    CDS_COORD_t  minPos = ((isUnitig) ?
                           MIN(frag->offset5p.mean, frag->offset3p.mean) :
                           MIN(frag->contigOffset5p.mean, frag->contigOffset3p.mean));
    CDS_COORD_t  maxPos = ((isUnitig) ?
                           MAX(frag->offset5p.mean, frag->offset3p.mean) :
                           MAX(frag->contigOffset5p.mean, frag->contigOffset3p.mean));

    //  keep contained with containing - avoids problem of consensus
    //  failures if overlap of contained & next fragment is very short
    //
    if(imp->contained) {
      AddIMPToIUMStruct(lastPlacement, imp);
      continue;
    }

    lastIID = imp->ident;
    // determine where the fragment goes wrt the chimeric interval
    if(minPos >= currInterval->bgn) {
      if(minPos > currInterval->end) {

        // Three possibilities here:
        //
        // 1. we've run out of intervals & the rest of the fragments
        // go in the last good interval (handled by the following if())
        //
        // 2. the fragment is not beyond the next bad interval
        //
        // 3. the fragment is beyond the next bad interval
        //
        if(currI < GetNumVA_SeqInterval(csis)) {
          SeqInterval *lastInterval = NULL;
          // chimeric intervals may be close together - close enough that
          // the current fragment may be in the next 'bad' interval rather
          // than the next 'good' interval - or even further down...
          // we ARE done with the current good/bad IUMs.
          //
          // old bad & old good are done
          if(good->ium.num_frags > 0) {
            good->ium.unique_rept = AS_FORCED_NONE;
            StoreIUMStruct(graph, good, isUnitig, AS_UNASSIGNED, egfar);
            ci = GetGraphNode(graph->CIGraph, ciID);
          }

          if(bad->ium.num_frags > 0) {
            bad->ium.unique_rept = AS_FORCED_NONE;
            StoreIUMStruct(graph, bad, isUnitig, AS_UNASSIGNED, egfar);
            ci = GetGraphNode(graph->CIGraph, ciID);
          }

          // reset the old
          ResetIUMStruct(good);
          ResetIUMStruct(bad);

          // advance to the next relevant interval
          while(currI < GetNumVA_SeqInterval(csis) && minPos > currInterval->end) {
            lastInterval = currInterval;
            currInterval = GetVA_SeqInterval(csis, ++currI);
          }

          if(currInterval) {
            // currInterval is relevant bad interval
            // place frag in good or bad
            if(minPos >= currInterval->bgn) {
              if(minPos > currInterval->end) {
                AddIMPToIUMStruct((lastPlacement = good), imp);
              } else {
                AddIMPToIUMStruct((lastPlacement = bad), imp);
              }
            } else {
              if(maxPos >= currInterval->bgn) {
                AddIMPToIUMStruct((lastPlacement = bad), imp);
              } else {
                AddIMPToIUMStruct((lastPlacement = good), imp);
              }
            }
          } else {
            // ran out of chimeric intervals - in last good interval
            AddIMPToIUMStruct((lastPlacement = good), imp);
            currInterval = lastInterval;
            currI = GetNumVA_SeqInterval(csis);
          }
        } else {
          AddIMPToIUMStruct((lastPlacement = good), imp);
        }
      } else {
        // goes in bad
        AddIMPToIUMStruct((lastPlacement = bad), imp);
      }
    } else {
      if(maxPos >= currInterval->bgn) {
        // goes in bad
        AddIMPToIUMStruct((lastPlacement = bad), imp);
      } else {
        // goes in good
        AddIMPToIUMStruct((lastPlacement = good), imp);
      }
    }
  }

  // now make the last good chunk a chunk
  if(good->ium.num_frags > 0) {
    good->ium.unique_rept = AS_FORCED_NONE;
    StoreIUMStruct(graph, good, isUnitig, AS_UNASSIGNED, egfar);
    ci = GetGraphNode(graph->CIGraph, ciID);
  }

  // if the last bad section went to the end, it's still here
  if(bad->ium.num_frags > 0) {
    bad->ium.unique_rept = AS_FORCED_NONE;
    StoreIUMStruct(graph, bad, isUnitig, AS_UNASSIGNED, egfar);
    ci = GetGraphNode(graph->CIGraph, ciID);
  }

  // delete the original unitig & contig
  if(ci->info.CI.contigID != NULLINDEX) {
    ChunkInstanceT *contig = GetGraphNode(graph->ContigGraph, ci->info.CI.contigID);
    if(contig)
      DeleteGraphNode(graph->ContigGraph, contig);
  }

  DeleteGraphNode(graph->CIGraph, ci);
}


void
SplitInputUnitigs(ScaffoldGraphT *graph) {
  VA_TYPE(uint16) *rc = NULL;   // read coverage map
  VA_TYPE(uint16) *gcc = NULL;  // good clone coverage
  VA_TYPE(uint16) *bcc = NULL;  // bad clone coverage
  VA_TYPE(SeqInterval) *csis = NULL; // chimeric sequence intervals
  CDS_COORD_t minLength = CDS_COORD_MAX;
  int i;
  int numCIs = GetNumGraphNodes(graph->CIGraph);

  rc = CreateVA_uint16(10000);
  gcc = CreateVA_uint16(10000);
  bcc = CreateVA_uint16(10000);
  csis = CreateVA_SeqInterval(100);
  assert(rc != NULL && gcc != NULL && bcc != NULL && csis != NULL);

  // determine minimum unitig length of interest
  //
  for (i=1; i<GetNumDistTs(graph->Dists); i++) {
    DistT *dptr = GetDistT(graph->Dists,i);
    minLength = MIN(minLength,dptr->mu + CGW_CUTOFF * dptr->sigma);
  }

  // loop over number of original unitigs - not new ones being
  // generated
  //
  for (i=0; i<numCIs; i++) {
    ChunkInstanceT *ci = GetGraphNode(graph->CIGraph, i);
    MultiAlignT *ma = loadMultiAlignTFromSequenceDB(graph->sequenceDB, ci->id, TRUE);

    // NOTE: add discriminator statistic checks?

    if(GetMultiAlignLength(ma) >= minLength) {
      CDS_COORD_t minBase, maxBase;
      CDS_COORD_t curBase;

      // create read coverage map for unitig

      ResetVA_uint16(rc);
      EnableRangeVA_uint16(rc, GetMultiAlignLength(ma));
      CreateReadCoverageMap(graph, rc, ma, TRUE);

      // locate region within which to look for possible chimeric
      // points i.e., ignore initial & trailing 0/1 values

      for(minBase = READ_TRIM_BASES;
          minBase < GetMultiAlignLength(ma) - READ_TRIM_BASES && *(GetVA_uint16(rc,minBase)) <= 1;
          minBase++)
        ;

      for(maxBase = GetMultiAlignLength(ma) - READ_TRIM_BASES;
          maxBase > READ_TRIM_BASES && *(GetVA_uint16(rc,maxBase)) <= 1;
          maxBase--)
        ;

      // see if there is a candidate interval

      for( curBase = minBase; curBase < maxBase; curBase++)
        if(*(GetVA_uint16(rc,curBase)) <= MAX_SEQUENCE_COVERAGE)
          break;

      // see if above loop ended in a candidate interval

      if(curBase < maxBase) {
        CDS_COORD_t checkBase;
        SeqInterval interval;
        int inInterval = 0;

        // create gcc & bcc maps for unitig
        ResetVA_uint16(gcc);
        EnableRangeVA_uint16(gcc, GetMultiAlignLength(ma));

        ResetVA_uint16(bcc);
        EnableRangeVA_uint16(bcc, GetMultiAlignLength(ma));

        // reset the number of splitting intervals to 0
        ResetVA_SeqInterval(csis);
        EnableRangeVA_SeqInterval(csis, 0);

        CreateCloneCoverageMaps(graph, gcc, bcc, ma, TRUE);

        // identify & count chimeric sequence intervals
        for(checkBase = curBase; checkBase < maxBase; checkBase++) {
          if(*(GetVA_uint16(rc,checkBase)) <= MAX_SEQUENCE_COVERAGE &&
             *(GetVA_uint16(bcc,checkBase)) >= MIN_BAD_CLONE_COVERAGE &&
             *(GetVA_uint16(gcc,checkBase)) <= MAX_GOOD_CLONE_COVERAGE) {
            if(inInterval) {
              // continuing in interval
              interval.end = checkBase;
            } else {
              // starting interval
              interval.bgn = interval.end = checkBase;
              inInterval = 1;
            }
          } else {
            if(inInterval) {
              // ended interval
              // if it is more than a minimum distance from the last one,
              // add it. otherwise, combine the two - since no fragment can
              // possibly be entirely within the intervening 'good' interval
              //
              if(GetNumVA_SeqInterval(csis) > 0) {
                SeqInterval *tempSI = GetVA_SeqInterval(csis, GetNumVA_SeqInterval(csis) - 1);

                if(tempSI->end < interval.bgn - AS_FRAG_MIN_LEN)
                  AppendVA_SeqInterval(csis, &interval);
                else
                  tempSI->end = interval.end;
              } else
                AppendVA_SeqInterval(csis, &interval);
              inInterval = 0;
            }
            // otherwise continuing not in interval - do nothing
          }
        }

        if(GetNumVA_SeqInterval(csis) > 0) {
          SplitChunkByIntervals(graph, ci->id, ma, csis, TRUE);
          ci = GetGraphNode(graph->CIGraph, i);
        }
      }
    }
  }

  DeleteVA_SeqInterval(csis);
  DeleteVA_uint16(rc);
  DeleteVA_uint16(gcc);
  DeleteVA_uint16(bcc);

  //CheckUnitigs(0, GetNumGraphNodes(graph->CIGraph));
}

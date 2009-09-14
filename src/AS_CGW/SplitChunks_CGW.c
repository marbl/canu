
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

static char *rcsid = "$Id: SplitChunks_CGW.c,v 1.51 2009-09-14 16:09:04 brianwalenz Exp $";

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


#define READ_TRIM_BASES          (AS_OVERLAP_MIN_LEN / 2 - 1)   //  was 30
#define MAX_SEQUENCE_COVERAGE    1
#define MIN_BAD_CLONE_COVERAGE   3
#define MAX_GOOD_CLONE_COVERAGE  0
#define LN_2                     0.693147

VA_DEF(uint16);
VA_DEF(SeqInterval);

typedef struct {
  IntUnitigMesg    ium;
  int32            isGood;
  int32            intbgn;
  int32            intend;
  int32            utgbgn;
  int32            utgend;
} IUMStruct;




static
void
IncrementMapInterval(char             *name,
                     VA_TYPE(uint16)  *map,
                     int32             minPos,
                     int32             maxPos) {
  uint16     *t = GetVA_uint16(map, 0);
  uint32      l = GetNumuint16s(map);

  //fprintf(stderr, "ADD %s %d-%d\n", name, minPos, maxPos);

  for (int32 i=minPos; (i<maxPos) && (i<l); i++)
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
              int32 length,
              int isUnitig) {

  DistT *dist    = GetDistT(graph->Dists, distID);
  int32  distMin = dist->mu - CGW_CUTOFF * dist->sigma;
  int32  distMax = dist->mu + CGW_CUTOFF * dist->sigma;

  int32 minPos = ((isUnitig) ?
                  (MIN(frag->offset5p.mean, mfrag->offset5p.mean)) :
                  (MIN(frag->contigOffset5p.mean, mfrag->contigOffset5p.mean)));
  int32 maxPos = ((isUnitig) ?
                  (MAX(frag->offset5p.mean, mfrag->offset5p.mean)) :
                  (MAX(frag->contigOffset5p.mean, mfrag->contigOffset5p.mean)));

  // if they're in the same unitig
  if((isUnitig && frag->cid == mfrag->cid) ||
     (!isUnitig && frag->contigID == mfrag->contigID)) {
    // for pairs in the same unitig, just process the lesser
    if(frag->read_iid < mfrag->read_iid) {
      // if orientation is the same, the pair is bad
      if((isUnitig && getCIFragOrient(frag) == getCIFragOrient(mfrag)) ||
         (!isUnitig && GetContigFragOrient(frag) == GetContigFragOrient(mfrag))) {
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
            IncrementMapInterval("BCC", bcc, minPos, MIN(minPos + distMax, maxPos));
            IncrementMapInterval("BCC", bcc, maxPos, MIN(maxPos + distMax, length - 1));
          } else {
            // both are B_A oriented, so intervals < 5p are suspect
            // Increment from minPos to somewhere & maxPos to somewhere
            //              minPos                   maxPos
            //                  |                         |
            //          <--------                 <--------
            // ?--- increment ---        ?--- increment ---
            //
            IncrementMapInterval("BCC", bcc, MAX(0, minPos - distMax), minPos);
            IncrementMapInterval("BCC", bcc, MAX(minPos, maxPos - distMax), maxPos);
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
            IncrementMapInterval("BCC", bcc, MAX(0, minPos - distMax), minPos);
            IncrementMapInterval("BCC", bcc, MAX(minPos, maxPos - distMax), maxPos);
          } else {
            // both are B_A oriented, so
            // Increment from minPos to somewhere & maxPos to somewhere
            //     minPos                       maxPos
            //         |                            |
            // <--------                    <--------
            //          --- increment ---?           --- increment ---?
            //
            IncrementMapInterval("BCC", bcc, minPos, MIN(minPos + distMax, maxPos));
            IncrementMapInterval("BCC", bcc, maxPos, MIN(length - 1, maxPos + distMax));
          }
        }
      } else {
        // fragments are oriented differently - may be okay
        int32 distance = maxPos - minPos;

        // if the distance is wrong, the pair is bad
        if(distance < distMin ||
           distance > distMax) {
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
          IncrementMapInterval("BCC", bcc, minPos, MIN(maxPos, minPos + distMax));
          IncrementMapInterval("BCC", bcc, MAX(minPos + distMax, maxPos - distMax), maxPos);
        } else {
          // good pair
          IncrementMapInterval("GCC", gcc, minPos, maxPos);
        }
      }
    }
  } else {
    // in different unitigs, should be close to end of unitig
    if(orient == AS_INNIE) {
      if((isUnitig && getCIFragOrient(frag) == A_B && frag->offset5p.mean < length - distMax) ||
         (!isUnitig && GetContigFragOrient(frag) == A_B && frag->contigOffset5p.mean < length - distMax)) {
        //
        //  --------->
        //  --- increment ---
        //
        IncrementMapInterval("BCC", bcc,
                             ((isUnitig) ? frag->offset5p.mean : frag->contigOffset5p.mean),
                             ((isUnitig) ? frag->offset5p.mean : frag->contigOffset5p.mean) + distMax);
      } else if((isUnitig && getCIFragOrient(frag) == B_A && frag->offset5p.mean > distMax) ||
                (!isUnitig && GetContigFragOrient(frag) == B_A && frag->contigOffset5p.mean > distMax)) {
        //
        //        <----------
        //  --- increment ---
        //
        IncrementMapInterval("BCC", bcc,
                             ((isUnitig) ? frag->offset5p.mean : frag->contigOffset5p.mean) - distMax,
                             ((isUnitig) ? frag->offset5p.mean : frag->contigOffset5p.mean));
      }
    } else {
      // outtie
      if((isUnitig && getCIFragOrient(frag) == B_A && frag->offset5p.mean < length - distMax) ||
         (!isUnitig && GetContigFragOrient(frag) == B_A && frag->contigOffset5p.mean < length - distMax)) {
        //
        //  <----------
        //             --- increment ---
        //
        IncrementMapInterval("BCC", bcc,
                             ((isUnitig) ? frag->offset5p.mean : frag->contigOffset5p.mean),
                             ((isUnitig) ? frag->offset5p.mean : frag->contigOffset5p.mean) + distMax);

      } else if((isUnitig && getCIFragOrient(frag) == A_B && frag->offset5p.mean > distMax) ||
                (!isUnitig && GetContigFragOrient(frag) == A_B && frag->contigOffset5p.mean > distMax)) {
        //
        //                   --------->
        //  --- increment ---
        //
        IncrementMapInterval("BCC", bcc,
                             ((isUnitig) ? frag->offset5p.mean : frag->contigOffset5p.mean) - distMax,
                             ((isUnitig) ? frag->offset5p.mean : frag->contigOffset5p.mean));
      }
    }
  }
}


static
void
CreateReadCoverageMap(ScaffoldGraphT *graph,
                      VA_TYPE(uint16) *rc,
                      MultiAlignT *ma,
                      int32 minSplit,
                      int isUnitig) {

  ResetVA_uint16(rc);
  EnableRangeVA_uint16(rc, GetMultiAlignLength(ma));

  uint32 n      = GetNumIntMultiPoss(ma->f_list);

  int32  length = GetMultiAlignLength(ma);

  for (uint32 i=0; i<n; i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    if (!AS_FA_READ(imp->type))
      continue;

    CIFragT   *frag = GetCIFragT(graph->CIFrags, imp->ident);

    int32 minPos = ((isUnitig) ?
                    (MIN(frag->offset5p.mean, frag->offset3p.mean) + READ_TRIM_BASES) :
                    (MIN(frag->contigOffset5p.mean, frag->contigOffset3p.mean) + READ_TRIM_BASES));

    int32 maxPos = ((isUnitig) ?
                    (MAX(frag->offset5p.mean, frag->offset3p.mean) - READ_TRIM_BASES) :
                    (MAX(frag->contigOffset5p.mean, frag->contigOffset3p.mean) - READ_TRIM_BASES));

    //  Undo the offsets if we're near the end of a unitig where we expect coverage to be low.

    if ((i < 5) || (maxPos < minSplit)) {
      minPos -= READ_TRIM_BASES;
      maxPos += READ_TRIM_BASES;
    }

    if ((i >= n-5) || (minPos + minSplit > length)) {
      minPos -= READ_TRIM_BASES;
      maxPos += READ_TRIM_BASES;
    }

    IncrementMapInterval("RC ", rc, minPos, maxPos);
  }
}


static
void
CreateCloneCoverageMaps(ScaffoldGraphT *graph,
                        VA_TYPE(uint16) *gcc,
                        VA_TYPE(uint16) *bcc,
                        MultiAlignT *ma,
                        int isUnitig) {

  ResetVA_uint16(gcc);
  EnableRangeVA_uint16(gcc, GetMultiAlignLength(ma));

  ResetVA_uint16(bcc);
  EnableRangeVA_uint16(bcc, GetMultiAlignLength(ma));

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp  = GetIntMultiPos(ma->f_list, i);
    CIFragT     *frag = GetCIFragT(graph->CIFrags, imp->ident);

    // this is unlike ComputeMatePairStatistics, since we're interested
    // even in pairs that are in different unitigs/contigs

    if ((frag->flags.bits.hasMate > 0) &&
        (frag->mate_iid           > 0))
      AddLinkToMaps(graph, gcc, bcc, frag,
                    GetCIFragT(graph->CIFrags, frag->mate_iid),
                    ((frag->flags.bits.innieMate) ? AS_INNIE : AS_OUTTIE),
                    frag->dist,
                    GetMultiAlignLength(ma),
                    isUnitig);
  }
}





static
void
StoreIUMStruct(ScaffoldGraphT *graph,
               IUMStruct *is,
               int isUnitig,
               float egfar) {
  int32 rho = 0;

  is->ium.iaccession = ((isUnitig) ? GetNumGraphNodes(graph->CIGraph) : GetNumGraphNodes(graph->ContigGraph));

  is->ium.status      = AS_UNASSIGNED;
  is->ium.forced      = 0;
  is->ium.length      = 0;
  is->ium.unique_rept = AS_FORCED_NONE;

  int32 minPos  = INT32_MAX;
  int32 numRand = 0;

  for (int i=0; i<is->ium.num_frags; i++) {
    minPos = MIN(minPos, is->ium.f_list[i].position.bgn);
    minPos = MIN(minPos, is->ium.f_list[i].position.end);

    numRand += (AS_FA_RANDOM(is->ium.f_list[i].type)) ? 1 : 0;
  }

  for (int i=0; i<is->ium.num_frags; i++) {
    is->ium.f_list[i].position.bgn -= minPos;
    is->ium.f_list[i].position.end -= minPos;

    if (is->ium.f_list[i].position.end > is->ium.f_list[i].position.bgn) {
      rho            = MAX(rho,            is->ium.f_list[i].position.bgn);
      is->ium.length = MAX(is->ium.length, is->ium.f_list[i].position.end);
    } else {
      rho            = MAX(rho,            is->ium.f_list[i].position.end);
      is->ium.length = MAX(is->ium.length, is->ium.f_list[i].position.bgn);
    }
  }


  // egfar is the estimated global fragment arrival rate based on the original chunk.
  //
  // rho is the sum of the a-hangs (= largest starting position)
  //
  // Based on compute_coverage_statistic() in AS_CGB_cgb.h
  //
  is->ium.coverage_stat = 0.0;
  if ((egfar > 0.0) &&
      (numRand > 0))
    is->ium.coverage_stat = rho * egfar - LN_2 * (numRand - 1);

  // get the multi-alignment - this fills in some unitig fields

  VA_TYPE(int32)   *deltas     = CreateVA_int32(1);
  VA_TYPE(char)    *sequence   = CreateVA_char(200000);
  VA_TYPE(char)    *quality    = CreateVA_char(200000);

  int unitigSuccess = MultiAlignUnitig(&is->ium,
                                       ScaffoldGraph->gkpStore,
                                       sequence,
                                       quality,
                                       deltas,
                                       CNS_STATS_ONLY,
                                       NULL);

  if (unitigSuccess == 0) {
    GenericMesg pmesg;

    pmesg.t = MESG_IUM;
    pmesg.m = &is->ium;

    fprintf(stderr, "================================================================================\n");
    WriteProtoMesg_AS(stderr, &pmesg);
    fprintf(stderr, "================================================================================\n");
    fprintf(stderr, "FATAL ERROR: MultiAlignUnitig call failed in unitig splitting.\n");
    assert(FALSE);
  }

  //  Add ium to the system

  MultiAlignT    *ma = CreateMultiAlignTFromIUM(&is->ium, GetNumCIFragTs(graph->CIFrags), FALSE);

  //  Insert both a unitig and a contig

  insertMultiAlignTInSequenceDB(graph->sequenceDB, is->ium.iaccession, TRUE,
                                ma,
                                TRUE);
  insertMultiAlignTInSequenceDB(graph->sequenceDB, is->ium.iaccession, FALSE,
                                CopyMultiAlignT(NULL, ma),
                                FALSE);

  ProcessIUM_ScaffoldGraph(&is->ium, GetMultiAlignUngappedLength(ma), TRUE);

  ChunkInstanceT *ci = GetGraphNode((isUnitig ? graph->CIGraph : graph->ContigGraph), is->ium.iaccession);

  UpdateNodeFragments((isUnitig ? graph->CIGraph : graph->ContigGraph),
                      is->ium.iaccession,
                      ci->type == DISCRIMINATORUNIQUECHUNK_CGW,
                      TRUE);

  DeleteVA_int32(deltas);
  DeleteVA_char(sequence);
  DeleteVA_char(quality);
}






static
float
EstimateGlobalFragmentArrivalRate(ChunkInstanceT *ci, MultiAlignT *ma) {
  int32  numRF = 0;
  int32  rho   = 0;

  // estimate random fragments/bp

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    if (AS_FA_RANDOM(imp->type)) {
      rho = MAX(rho, MIN(imp->position.bgn, imp->position.end));
      numRF++;
    }
  }

  return((ci->info.CI.coverageStat + LN_2 * (numRF - 1)) / rho);
}



static
int
positionCompare(const void *A, const void *B) {
  const IntMultiPos *a = (const IntMultiPos *)A;
  const IntMultiPos *b = (const IntMultiPos *)B;

  //  A before B -> negative
  int  diff = MIN(a->position.bgn, a->position.end) - MIN(b->position.bgn, b->position.end);

  if (diff == 0) {
    //  A longer than B -> negative
    int alen = MAX(a->position.bgn, a->position.end) - MIN(a->position.bgn, a->position.end);
    int blen = MAX(b->position.bgn, b->position.end) - MIN(b->position.bgn, b->position.end);
    diff = blen - alen;
  }

  return(diff);
}



// NOTE: csis intervals are in ungapped coordinates
static
void
SplitChunkByIntervals(ScaffoldGraphT *graph,
                      CDS_CID_t ciID,
                      MultiAlignT *ma,
                      VA_TYPE(SeqInterval) *csis,
                      int isUnitig) {
  ChunkInstanceT *ci = GetGraphNode(graph->CIGraph, ciID);

  assert(isUnitig);  // not implemented for contigs yet

  float egfar = EstimateGlobalFragmentArrivalRate(ci, ma);

  fprintf(stderr, "Splitting %s "F_CID " into as many as %d %s at intervals:",
          (isUnitig ? "unitig" : "contig"), ma->maID,
          (int) (2 * GetNumVA_SeqInterval(csis) + 1),
          (isUnitig ? "unitigs" : "contigs"));

  for (uint32 i=0; i<GetNumVA_SeqInterval(csis); i++) {
    SeqInterval *I = GetVA_SeqInterval(csis, i);
    fprintf(stderr, "\t"F_S32","F_S32, I->bgn, I->end);
  }

  fprintf(stderr, "\n");

  
  //  If a fragment even touches a bad interval, it gets placed in that bad interval.
  //
  //  This causes problems in cases such as:
  //
  //        |     bad      |    good     |
  //      --------------------
  //                         --
  //                      ----------------
  //                               -------
  //
  //  If the little fragment is not marked as contained, it will be placed in the good interval,
  //  resulting in a disconnected unitig.  This can (and does) happen due to consensus needing to
  //  align the second to last fragment with a negative ahang.
  //
  //  We get around this by explicitly checking if the next fragment is completely contained in the
  //  unitig associated with the previous interval.

  int32            utgNum = 2 * GetNumVA_SeqInterval(csis) + 1;
  IUMStruct       *utg    = (IUMStruct *)safe_calloc(utgNum, sizeof(IUMStruct));

  for (int i=0; i<utgNum; i++) {
    utg[i].ium.f_list = (IntMultiPos *)safe_calloc(GetNumIntMultiPoss(ma->f_list), sizeof(IntMultiPos));
    utg[i].isGood = 1;

    utg[i].intbgn = INT32_MIN;
    utg[i].intend = INT32_MAX;

    utg[i].utgbgn = INT32_MAX;
    utg[i].utgend = INT32_MIN;
  }


  for (uint32 currI=0; currI < GetNumVA_SeqInterval(csis); currI++) {
    int cp = 2 * currI;
    int ci = 2 * currI + 1;
    int cn = 2 * currI + 2;

    utg[ci].isGood = 0;
    utg[cp].intend = utg[ci].intbgn = GetVA_SeqInterval(csis, currI)->bgn;
    utg[ci].intend = utg[cn].intbgn = GetVA_SeqInterval(csis, currI)->end;
  }

  //  Sort fragments based on their start position.  This is needed so that contained fragments get
  //  placed in the correct unitig.

  qsort(ma->f_list->Elements, GetNumIntMultiPoss(ma->f_list), sizeof(IntMultiPos), positionCompare);

  //  Then place fragments in new unitigs.

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp  = GetIntMultiPos(ma->f_list, i);
    CIFragT     *frag = GetCIFragT(graph->CIFrags, imp->ident);

    //  Determine the position of the fragment in the unitig or contig

    int32  minPos = ((isUnitig) ?
                     MIN(frag->offset5p.mean, frag->offset3p.mean) :
                     MIN(frag->contigOffset5p.mean, frag->contigOffset3p.mean));
    int32  maxPos = ((isUnitig) ?
                     MAX(frag->offset5p.mean, frag->offset3p.mean) :
                     MAX(frag->contigOffset5p.mean, frag->contigOffset3p.mean));

    //  Find the interval it goes in.  We don't expect to have very many intervals, so this isn't
    //  as terrible as it looks.

    int32  interval = -1;

    //  Do we belong to a toxic interval?  Anything that touches it belongs in it.
    //
    for (int i=0; (interval == -1) && (i < utgNum); i++)
      if ((utg[i].isGood == 0) && (maxPos >= utg[i].intbgn) && (minPos <= utg[i].intend))
        interval = i;

    //  If not in an interval, it must intersect exactly one good interval (otherwise, it would be
    //  spanning a bad interval).  This lets us use the same test as for 'toxic'; simply touching a
    //  good interval puts the fragment in it.
    //
    for (int i=0; (interval == -1) && (i < utgNum); i++)
      if ((utg[i].isGood == 1) && (maxPos >= utg[i].intbgn) && (minPos <= utg[i].intend))
        interval = i;

    assert(interval >= 0);
    assert(interval <  utgNum);

    //  Adjust the placement in a few cases

    //  1) If this fragment is not overlapping with the unitig that exists in this interval, check
    //  if the next interval contains it.  Place it where the thicker overlap exists.
    //
    //  -------------BADINTERVAL---
    //  ---------
    //    ---------------
    //        --
    //            -- <- is in 'good' but there is no sequence to support it
    //
    if ((interval < utgNum - 1) &&
        (utg[interval].isGood == 1) &&
        (utg[interval].utgend - AS_OVERLAP_MIN_LEN < minPos)) {
      int32  thisOvl = utg[interval].utgend - minPos;
      int32  thatOvl = maxPos - utg[interval+1].utgbgn;

      assert((thisOvl > 0) || (thatOvl > 0));

      if (thatOvl > thisOvl)
        interval++;
    }

    //  2) Similar to the above, but looking backwards.  Move a fragment if it is completely
    //  contained in a previous bad interval, and if the current interval is empty.
    //
    //  ----------BAD1----BAD2----------
    //   --------------            (in bad1)
    //       --------------------  (in bad1)
    //                   -------   (in bad2)
    //                         --- (**)
    //
    //  The ** fragment does NOT need to be moved.  It is not contained.  The only way that this
    //  unitig could have been originally connected is by an overlap from ** to a later fragment
    //  (OK, except if that overlap is high error noise).  Because this fragment isn't contained, it
    //  will have the thickest overlap to any later fragment.
    //
    if ((interval > 0) &&
        (utg[interval].utgend < utg[interval].utgbgn)) {
      for (int32 newint=interval-1; newint>=0; newint--) {
        if (utg[newint].isGood == 1)
          continue;
        if ((maxPos <= utg[newint].utgend) &&
            (utg[newint].utgbgn <= minPos))
          interval = newint;
      }
    }

    //  3) If we are a in a bad interval, but the previous good interval is empty then that good
    //  interval will never get fragments; we're already past it (fragments are sorted).  Instead of
    //  making two adjacent bad intervals, we merge them together.
    //
    //  The first bad interval we need to worry about is interval 3.
    //
    if ((interval > 2) &&
        (utg[interval].isGood == false) &&
        (utg[interval-1].utgend < utg[interval-1].utgbgn)) {
      int32 newint = interval - 2;

      //  Search for the first non-empty unitig before us.
      while ((newint >= 0) && (utg[newint].utgend < utg[newint].utgbgn))
        newint--;

      //  If that non-empty unitig is a bad one, reassign the fragment.
      if ((newint > 0) &&
          (utg[newint].isGood == 0)) {
        fprintf(stderr, "REASSIGN bad fragment from interval %d to interval %d\n", interval, newint);
        interval = newint;
      }
    }

    //  Place it.

    {
      CIFragT        *frag = GetCIFragT(ScaffoldGraph->CIFrags, imp->ident);
      IntUnitigMesg  *ium  = &utg[interval].ium;

      ium->f_list[ium->num_frags] = *imp;

      ium->f_list[ium->num_frags].position.bgn = ((isUnitig) ? frag->offset5p.mean : frag->contigOffset5p.mean);
      ium->f_list[ium->num_frags].position.end = ((isUnitig) ? frag->offset3p.mean : frag->contigOffset3p.mean);

      ium->num_frags++;
    }

    utg[interval].utgbgn = MIN(utg[interval].utgbgn, minPos);
    utg[interval].utgend = MAX(utg[interval].utgend, maxPos);

#ifdef DEBUG
    fprintf(stderr, "FRAG %d (%d,%d) placed in interval %d (%s) (now %d,%d)\n",
            imp->ident, minPos, maxPos,
            interval,
            utg[interval].isGood ? "GOOD" : "BAD",
            utg[interval].utgbgn, utg[interval].utgend);
#endif
  }

  //  Now, add all the new unitigs.

  for (int i=0; i<utgNum; i++) {
#ifdef DEBUG
    if (utg[i].ium.num_frags == 0)
      fprintf(stderr, "EMPTY UNITIG\n");
    if (utg[i].ium.num_frags > 0)
      fprintf(stderr, "NEW UNITIG with %d fragments\n", utg[i].ium.num_frags);
#endif
    if (utg[i].ium.num_frags > 0)
      StoreIUMStruct(graph, utg + i, isUnitig, egfar);

    safe_free(utg[i].ium.f_list);
  }

  safe_free(utg);

  //  Delete the original unitig and contig.

  ci = GetGraphNode(graph->CIGraph, ciID);

  if (ci->info.CI.contigID != NULLINDEX) {
    ChunkInstanceT *contig = GetGraphNode(graph->ContigGraph, ci->info.CI.contigID);
    if (contig)
      DeleteGraphNode(graph->ContigGraph, contig);
  }

  DeleteGraphNode(graph->CIGraph, ci);
}




//static
void
SplitInputUnitigs(ScaffoldGraphT *graph) {
  VA_TYPE(uint16)       *rc  = CreateVA_uint16(10000);      // read coverage map
  VA_TYPE(uint16)       *gcc = CreateVA_uint16(10000);      // good clone coverage
  VA_TYPE(uint16)       *bcc = CreateVA_uint16(10000);      // bad clone coverage
  VA_TYPE(SeqInterval)  *csis = CreateVA_SeqInterval(100);  // chimeric sequence intervals

  int32 minLength = INT32_MAX;
  int32 minSplit  = INT32_MAX;  //  Not exactly the "minimum split"...

  for (uint32 i=1; i<GetNumDistTs(graph->Dists); i++) {
    DistT *dptr = GetDistT(graph->Dists,i);

    fprintf(stderr, "LIB %d mu = %f sigma = %f\n", i, dptr->mu, dptr->sigma);

    if (dptr->mu > 0) {
      minLength = MIN(minLength, dptr->mu + CGW_CUTOFF * dptr->sigma);
      minSplit  = MIN(minSplit,  dptr->mu - CGW_CUTOFF * dptr->sigma);
    }
  }

  fprintf(stderr, "minLength = %d\n", minLength);
  fprintf(stderr, "minSplit  = %d\n", minSplit);

  //  Over all unitigs present now (we add unitigs as this loop progresses)

  int32 numUnitigsBefore = GetNumGraphNodes(graph->CIGraph);

  for (int32 i=0; i<numUnitigsBefore; i++) {
    ChunkInstanceT *ci = GetGraphNode(graph->CIGraph, i);
    MultiAlignT *ma = loadMultiAlignTFromSequenceDB(graph->sequenceDB, ci->id, TRUE);

    // NOTE: add discriminator statistic checks?

    if(GetMultiAlignLength(ma) < minLength)
      continue;

    int32 minBase, maxBase;
    int32 curBase;

    CreateReadCoverageMap(graph, rc, ma, minSplit, TRUE);

    //  Locate region within which to look for possible chimeric points i.e., ignore initial &
    //  trailing 0/1 values

    for(minBase = READ_TRIM_BASES;
        minBase < GetMultiAlignLength(ma) - READ_TRIM_BASES && *(GetVA_uint16(rc,minBase)) <= 1;
        minBase++)
      ;

    for(maxBase = GetMultiAlignLength(ma) - READ_TRIM_BASES;
        maxBase > READ_TRIM_BASES && *(GetVA_uint16(rc,maxBase)) <= 1;
        maxBase--)
      ;

    //  Find a candidate interval

    for( curBase = minBase; curBase < maxBase; curBase++)
      if(*(GetVA_uint16(rc,curBase)) <= MAX_SEQUENCE_COVERAGE)
        break;

    //  A good candidate interval?

    if(curBase >= maxBase)
      continue;

    SeqInterval interval;
    int         inInterval = 0;

    CreateCloneCoverageMaps(graph, gcc, bcc, ma, TRUE);

    ResetVA_SeqInterval(csis);
    EnableRangeVA_SeqInterval(csis, 0);

    //  Identify & count chimeric sequence intervals

    for (int32 checkBase = curBase; checkBase < maxBase; checkBase++) {
      if (*(GetVA_uint16(rc, checkBase)) <= MAX_SEQUENCE_COVERAGE &&
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

            if(tempSI->end < interval.bgn - AS_READ_MIN_LEN)
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

  DeleteVA_SeqInterval(csis);
  DeleteVA_uint16(rc);
  DeleteVA_uint16(gcc);
  DeleteVA_uint16(bcc);

  //CheckUnitigs(0, GetNumGraphNodes(graph->CIGraph));
}

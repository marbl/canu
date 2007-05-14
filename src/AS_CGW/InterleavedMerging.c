
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

#include "AS_global.h"
#include "AS_UTL_Var.h"

#include "InterleavedMerging.h"

#include "AS_CGW_dataTypes.h"
#include "InputDataTypes_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "GraphCGW_T.h"
#include "UtilsREZ.h"
#include "ChiSquareTest_CGW.h"

#include "CA_ALN_local.h"
#include "CA_ALN_scafcomp.h"

#define MIN_GAP_LENGTH 50

#define INTERLEAVE_CUTOFF   3.5

#undef DEBUG1
#undef PRINT_OVERLAPS

#define CONNECTEDNESS_CHECKS  

//  Define this to enable checking of the stretching / compression
//  of interleaved scaffold merges.  See the detailed comment at the
//  code.
#undef CHECK_INTERLEAVE_DISTANCE



/*
  What this code is supposed to do:
  Take two scaffolds & a negative edge between them & determine if the
  two scaffolds can be merged - either by overlapping one or more contigs
  or by interleaving with no contig overlaps

  How to use this code:
  1. Creat a reusable ScaffoldAlignmentInterface object:
  
  ScaffoldAlignmentInterface * sai = CreateScaffoldAlignmentInterface();

    
  2. Populate the sai with a couple scaffolds & edge:
  
  PopulateScaffoldAlignmentInterface(scaffoldA, scaffoldB, edge, sai);

  This populates the sai object with scaffold & edge information as
  well as contig overlaps that may be involved
      
  This function returns 0 if successful, non-0 if failure

    
  3. Call Align_Scaffold to see if scaffolds can be aligned
  
  Align_Scaffold(sai->segmentList,
  sai->numSegs,
  sai->varWin,
  sai->scaffoldA->scaffold,
  sai->scaffoldB->scaffold,
  &(sai->best),
  sai->scaffoldA->bandBeg,
  sai->scaffoldA->bandEnd);

  This function returns:
  a) non-NULL if the scaffolds can be merged with one or more
  contig overlaps
  b) NULL and sai->best >= 0 if the scaffolds can be merged
  via interleaving with no contig overlapping
  c) NULL and sai->best < 0 if the scaffolds cannot be merged

  The non-NULL returned pointer can be ignored.

      
  4. If the scaffolds are mergeable, call MakeScaffoldAlignmentAdjustments
  to adjust contig positions in the scaffolds and adjust edge distance
  so that subsequent call to InsertScaffoldContentsIntoScaffold will
  merge the scaffolds correctly

  MakeScaffoldAlignmentAdjustments(scaffoldA, scaffoldB, edge, sai);

  This function returns:
  a) non-NULL pointer to an edge (statically allocated in the function)
  with adjusted distance, if everything is well
  b) NULL if adjustments failed




  Notes on call to Gene's/Aaron's Align_Scaffold function:
  
  Segment *Align_Scaffold(Segment *seglist, int numsegs, int varwin,
  Scaffold *AF, Scaffold *BF, int *best)
                          
  returns a pointer to a (linked) list of Segment,
  namely the set of contig/contig overlaps that were used,
  or NULL if unsuccessful
    
  takes as arguments:
  a pointer to a (linked) list of Segment,
  the set of contig/contig overlaps that are allowed
  (but not required) to be part of the solution
  the number of segments in this (input) linked list
  the number of standard deviations away from the mean you
  consider acceptable (ie. probably 3)
  a pointer to the (appropriately specified) first scaffold
  a pointer to the (appropriately specified) second scaffold
  a pointer to an integer that will receive the alignment score
  of a successful alignment

  Pre-processing must generate all arguments to this function from
  a pointer to the graph
  a pointer to scaffold A
  a pointer to scaffold B
  a pointer to the edge between A and B

  Post-processing must examine the contents of the Segment returned and
  possibly adjust the edge and contig positions in one or both scaffolds
*/

typedef struct
{
  int minIndex;
  int maxIndex;
  CDS_COORD_t minCoord;
  CDS_COORD_t maxCoord;
} ContigSetInterval;

typedef struct
{
  int index;
  int firstOverlap;
  int lastOverlap;
  ContigSetInterval a;
  ContigSetInterval b;
} COSData;

VA_DEF(COSData);




int GetNumSegmentsInList(Segment * segmentList)
{
  int numSegments = 0;
  Segment * segment = segmentList;
  while(segment != NULL)
    {
      segment = segment->next;
      numSegments++;
    }
  return numSegments;
}


void PrintSegment(FILE * fp, Segment * segment)
{
  fprintf(fp, "  A: %d, B: %d, beg,end,len: (" F_COORD ", " F_COORD ", " F_COORD ")\n",
          segment->a_contig, segment->b_contig,
          segment->overlap->begpos, segment->overlap->endpos,
          segment->overlap->length);
}

void PrintSegmentList(FILE * fp, Segment * segmentList)
{
  Segment * segment;
  for(segment = segmentList; segment != NULL; segment = segment->next)
    PrintSegment(fp, segment);
}

void PrintScafCompScaffold(FILE * fp, Scaffold * scaffold)
{
  int i;
  fprintf(fp, "  length: " F_COORD "\n", scaffold->length);
  fprintf(fp, "  %d gaps:\n", scaffold->num_gaps);
  fprintf(fp, "    tig   left: %10" F_COORDP ",   length: %10" F_COORDP "\n",
          scaffold->ctgs[0].lft_end, scaffold->ctgs[0].length);
  for(i = 0; i < scaffold->num_gaps; i++)
    {
      fprintf(fp, "    gap length: %10" F_COORDP ",   stddev: %.2f\n",
              scaffold->gaps[i].gap_length,
              scaffold->gaps[i].gap_var);
      fprintf(fp, "    tig   left: %10" F_COORDP ",   length: %10" F_COORDP "\n",
              scaffold->ctgs[i+1].lft_end, scaffold->ctgs[i+1].length);
    }
}

void PrintContigElement(FILE * fp, ContigElement * ce)
{
  fprintf(fp, "    (%d, " F_CID ") len: %.f, min: %.f, max: %.f, orient: %c\n",
          ce->index, ce->id, ce->length, ce->minCoord, ce->maxCoord,
          ce->orient);
}


void PrintContigElements(FILE * fp, VA_TYPE(ContigElement) * contigs)
{
  int i;

  for(i = 0; i < GetNumVA_ContigElement(contigs); i++)
    {
      PrintContigElement(fp, GetVA_ContigElement(contigs, i));
    }
}

void PrintScaffoldStuff(FILE * fp, ScaffoldStuff * ss)
{
  PrintScafCompScaffold(fp, ss->scaffold);
  fprintf(fp, "  Contig elements:\n");
  PrintContigElements(fp, ss->contigs);
}

void PrintLocal_Overlaps(FILE * fp, VA_TYPE(Local_Overlap) * los)
{
  int i;
  Local_Overlap * lo = GetVA_Local_Overlap(los, 0);
  for(i = 0; i < GetNumVA_Local_Overlap(los); i++)
    {
      fprintf(fp, "  beg: " F_COORD ", end: " F_COORD ", len: " F_COORD "\n",
              lo[i].begpos, lo[i].endpos, lo[i].length);
    }
}

void PrintScaffoldAlignmentInterface(FILE * fp,
                                     ScaffoldAlignmentInterface * sai)
{
  fprintf(fp, "--------------- ScaffoldAlignmentInterface ----------------\n");
  fprintf(fp, "band begin,end, best: " F_COORD ", " F_COORD ", %d\n",
          sai->scaffoldA->bandBeg, sai->scaffoldA->bandEnd, sai->best);
  fprintf(fp, "\n");
  
  fprintf(fp, "Segment list (%d segments):\n", sai->numSegs);
  PrintSegmentList(fp, sai->segmentList);
  fprintf(fp, "\n");

  fprintf(fp, "Scaffold A:\n");
  PrintScaffoldStuff(fp, sai->scaffoldA);
  fprintf(fp, "\n");
  
  fprintf(fp, "Scaffold B:\n");
  PrintScaffoldStuff(fp, sai->scaffoldB);
  fprintf(fp, "\n");
}

  
void DeleteScaffoldPools(ScaffoldPools * sp)
{
  if(sp)
    {
      if(sp->gapPool)
        DeleteVA_Scaffold_Gap(sp->gapPool);
      if(sp->tigPool)
        DeleteVA_Scaffold_Tig(sp->tigPool);
      safe_free(sp);
    }
}


void DeleteScaffoldStuff(ScaffoldStuff * ss)
{
  if(ss)
    {
      if(ss->scaffold)
        safe_free(ss->scaffold);
      if(ss->pools)
        DeleteScaffoldPools(ss->pools);
      if(ss->contigs)
        DeleteVA_ContigElement(ss->contigs);
      if(ss->edgeContigs)
        DeleteVA_ContigElement(ss->edgeContigs);
    }
}

void DeleteSegmentList(Segment * segmentList)
{
  Segment * curr;
  Segment * next;
  for(curr = segmentList; curr != NULL; curr = next)
    {
      next = curr->next;
      safe_free(curr->overlap);
      safe_free(curr);
    }
}

Segment* DuplicateSegmentList(Segment * segmentList)
{
  Segment * curr;
  Segment * head=NULL;
  Segment * tail=NULL;
  Local_Overlap * ovl=NULL;
  for(curr = segmentList; curr != NULL; curr = curr->next)
    {
      Segment *this = (Segment *)safe_malloc(sizeof(Segment));
      ovl=(Local_Overlap*)safe_malloc(sizeof(Local_Overlap));
      *ovl=*(curr->overlap);
      *this = *curr;
      this->overlap=ovl;
      if(head==NULL){
        tail=head=this;
      } else {
        tail->next=this;
        tail=this;
      }
    }
  return head;
}

void DeleteScaffoldAlignmentInterface(ScaffoldAlignmentInterface * sai)
{
  if(sai)
    {
      if(sai->scaffoldA)
        DeleteScaffoldStuff(sai->scaffoldA);
    
      if(sai->scaffoldB)
        DeleteScaffoldStuff(sai->scaffoldB);

      DeleteSegmentList(sai->segmentList);

      if(sai->scaffInst)
        DestroyScaffoldInstrumenter(sai->scaffInst);

      safe_free(sai);
    }
}

ScaffoldPools * CreateScaffoldPools(void)
{
  ScaffoldPools * sp = safe_calloc(1, sizeof(ScaffoldPools));

  if(sp == NULL)
    return NULL;

  sp->gapPool = CreateVA_Scaffold_Gap(100);
  sp->tigPool = CreateVA_Scaffold_Tig(100);

  if(sp->gapPool == NULL ||
     sp->tigPool == NULL)
    {
      DeleteScaffoldPools(sp);
      return NULL;
    }
  return sp;
}


ScaffoldStuff * CreateScaffoldStuff(void)
{
  ScaffoldStuff * ss = safe_calloc(1, sizeof(ScaffoldStuff));

  if(ss == NULL)
    return NULL;

  ss->scaffold = safe_calloc(1, sizeof(Scaffold));
  ss->pools = CreateScaffoldPools();
  ss->contigs = CreateVA_ContigElement(100);
  ss->edgeContigs = CreateVA_ContigElement(100);

  if(ss->scaffold == NULL ||
     ss->pools == NULL ||
     ss->contigs == NULL ||
     ss->edgeContigs == NULL)
    {
      DeleteScaffoldStuff(ss);
      return NULL;
    }

  return ss;
}

ScaffoldAlignmentInterface * CreateScaffoldAlignmentInterface(void)
{
  ScaffoldAlignmentInterface * sai =
    safe_calloc(1, sizeof(ScaffoldAlignmentInterface));

  if(sai == NULL)
    return NULL;

  sai->scaffoldA = CreateScaffoldStuff();
  sai->scaffoldB = CreateScaffoldStuff();
  sai->scaffInst = CreateScaffoldInstrumenter(ScaffoldGraph, INST_OPT_ALL);
  if(sai->scaffoldA == NULL ||
     sai->scaffoldB == NULL ||
     sai->scaffInst == NULL)
    {
      fprintf(stderr, "Failed to allocate scaffold alignment interface structures!  (Out of memory?)\n");
      assert(0);
    }
  return sai;
}

void ResetScaffoldAlignmentScaffold(Scaffold * scf)
{
  scf->packed_seq = NULL;
  scf->num_gaps = 0;
  scf->gaps = NULL;
  scf->ctgs = NULL;
}


void ResetScaffoldPools(ScaffoldPools * sp)
{
  ResetVA_Scaffold_Gap(sp->gapPool);
  ResetVA_Scaffold_Tig(sp->tigPool);
}

void ResetScaffoldStuff(ScaffoldStuff * ss)
{
  ResetScaffoldAlignmentScaffold(ss->scaffold);
  ResetScaffoldPools(ss->pools);
  ResetVA_ContigElement(ss->contigs);
  ResetVA_ContigElement(ss->edgeContigs);
}

void ResetScaffoldAlignmentInterface(ScaffoldAlignmentInterface * sai)
{
  assert(sai != NULL);

  sai->numSegs = 0;
  sai->varWin = 3.0;
  sai->best = -1;
  
  ResetScaffoldStuff(sai->scaffoldA);
  ResetScaffoldStuff(sai->scaffoldB);

  DeleteSegmentList(sai->segmentList);
  sai->segmentList = NULL;
}


/*
  For Align_Scaffolds call, A & B have to appear to have NORMAL orientations.
  Two coordinate systems are used:
  
  1. treat scaffoldA is being to the left of scaffoldB
  left end of B = 0 with positive to the right
  this is used to determine which contigs may be 'in the edge'
  contigs in A that have a min coord >= 0 || max coord <= m may overlap
  such contigs in B
  -n         0          +m
  A     ---------------------->
  B               ----------------------->

  2. positions of contigs within scaffold are assigned with
  leftmost (start of) scaffold set to 0, increasing left to right
  This is used by Align_Scaffolds & to adjust contig positions later
  0                            j
  A    ------  ------  -----  ------(>)

  0                    k
  B    -----  ------- ------(>)
  
*/
void PopulateScaffoldStuff(ScaffoldStuff * ss,
                           CIScaffoldT * scaffold,
                           SEdgeT * sEdge,
                           LengthT osLength)
{
  CIScaffoldTIterator contigIterator;
  ChunkInstanceT * contig;
  LengthT thisLeft;
  LengthT thisRight;
  LengthT lastRight;
  int contigCount = 0;
  int isA = (sEdge->idA == scaffold->id) ? TRUE : FALSE;
  float thinEdge = -(sEdge->distance.mean +
                     INTERLEAVE_CUTOFF *
                     sqrt((double) sEdge->distance.variance));
  float thickEdge = -(sEdge->distance.mean -
                      INTERLEAVE_CUTOFF *
                      sqrt((double) sEdge->distance.variance));
  float osMin;
  float osMax;
  double varDelta;

#ifdef DEBUG1
  fprintf(GlobalData->stderrc, "****PopulateScaffoldStuff:\n");
  fprintf(GlobalData->stderrc, "  scaffold id:" F_CID ", length: %d, isA:%c\n",
          scaffold->id, (int) scaffold->bpLength.mean, isA ? 'T' : 'F');
  fprintf(GlobalData->stderrc, "  edge: %s, %d - %d\n",
          sEdge->orient == AB_AB ? "AB_AB" :
          (sEdge->orient == AB_BA ? "AB_BA" :
           (sEdge->orient == BA_AB ? "BA_AB" : "BA_BA")),
          (int) thinEdge, (int) thickEdge);
#endif
  /*
    loop through contigs in scaffold, left to right
    compute coordinates left to right
    also measure distances from left end of scaffoldB
  */
  if(sEdge->orient == AB_AB ||
     (isA && sEdge->orient == AB_BA) ||
     (!isA && sEdge->orient == BA_AB))
    ss->orient = A_B;
  else
    ss->orient = B_A;

  /*
    compute range in which to identify contigs that may be involved
    in overlaps. limit by the length of the other scaffold

    NOTE: increase variance going away from the end of scaffold
    'in the edge'
  */
  if(isA)
    {
      osMin = CGW_DP_MINLEN;
      osMax = osLength.mean +
        INTERLEAVE_CUTOFF * sqrt((double) osLength.variance) - CGW_DP_MINLEN;
      varDelta = (ss->orient == A_B) ? scaffold->bpLength.variance : 0.0;
    }
  else
    {
      osMin = thinEdge - (osLength.mean + INTERLEAVE_CUTOFF *
                          sqrt((double) osLength.variance)) + CGW_DP_MINLEN;
      osMax = thickEdge - CGW_DP_MINLEN;
      varDelta = (ss->orient == B_A) ? scaffold->bpLength.variance : 0.0;
    }
  
#ifdef DEBUG1
  fprintf(GlobalData->stderrc, "    scaffold is oriented %s\n",
          ss->orient == A_B ? "A_B" : "B_A");
  fprintf(GlobalData->stderrc, "    scaffold overlap limits: %d,%d\n",
          (int) osMin, (int) osMax);
#endif
  
  /*
    bandBeg must be less than or equal to bandEnd
    initialize maxAHang & flag as to whether or not minAHang has been set
    Check these values after looping over all contigs

    bandBeg:
    minCoord           maxCoord
    (-------            )
    (--------          )
    scfA: ------------------------------------
    |
    |
    scfB:                  ------------------------------------------------
    0

    thisLeft.mean - distance from left end of scfA to left end of contig
    bandBeg = thisLeft.mean - minCoord

      
    bandEnd:
    minCoord           maxCoord
    (          -------)
    (          -------)
    scfA: ------------------------------------
    |
    |
    scfB:                  ------------------------------------------------
    0

    thisLeft.mean - distance from left end of scfA to left end of contig
    bandEnd = thisLeft.mean - (maxCoord - contigLength)
  */
  ss->bandBeg = scaffold->bpLength.mean - CGW_DP_MINLEN;
  ss->bandEnd = CGW_DP_MINLEN;
  
  InitCIScaffoldTIterator(ScaffoldGraph, scaffold,
                          (ss->orient == A_B), FALSE, &contigIterator);
  while((contig = NextCIScaffoldTIterator(&contigIterator)) != NULL)
    {
      Scaffold_Tig tig;
      ContigElement ce;

      ce.index = contigCount;
      ce.id = contig->id;
      ce.length = contig->bpLength.mean;

      if(ss->orient == A_B)
        {
          // scaffold is -------------->
          if(contig->offsetAEnd.mean < contig->offsetBEnd.mean)  // contig: -->
            {
              thisLeft.mean = contig->offsetAEnd.mean;
              thisLeft.variance = fabs(varDelta - contig->offsetAEnd.variance);
              thisRight.mean = contig->offsetBEnd.mean;
              thisRight.variance = fabs(varDelta - contig->offsetBEnd.variance);
              ce.orient = A_B;
            }
          else // contig: <--
            {
              thisLeft.mean = contig->offsetBEnd.mean;
              thisLeft.variance = fabs(varDelta - contig->offsetBEnd.variance);
              thisRight.mean = contig->offsetAEnd.mean;
              thisRight.variance = fabs(varDelta - contig->offsetAEnd.variance);
              ce.orient = B_A;
            }
        }
      else
        {
          // scaffold is <--------------
          if(contig->offsetAEnd.mean < contig->offsetBEnd.mean) // contig: <--
            {
              thisLeft.mean = scaffold->bpLength.mean - contig->offsetBEnd.mean;
              thisLeft.variance = fabs(varDelta - contig->offsetBEnd.variance);
              thisRight.mean = scaffold->bpLength.mean - contig->offsetAEnd.mean;
              thisRight.variance = fabs(varDelta - contig->offsetAEnd.variance);
              ce.orient = B_A;
            }
          else // contig: -->
            {
              thisLeft.mean = scaffold->bpLength.mean - contig->offsetAEnd.mean;
              thisLeft.variance = fabs(varDelta - contig->offsetAEnd.variance);
              thisRight.mean = scaffold->bpLength.mean - contig->offsetBEnd.mean;
              thisRight.variance = fabs(varDelta - contig->offsetBEnd.variance);
              ce.orient = A_B;
            }
        }

      // store gap for Align_Scaffold interface
      if(contigCount > 0)
        {
          Scaffold_Gap gap;
          gap.gap_length = (int) ((thisLeft.mean - lastRight.mean) + 0.5);
          // yes, gap_var is actually used as stddev in Align_Scaffold()
          gap.gap_var = MAX(1.,sqrt(fabs((double) thisLeft.variance -
                                         lastRight.variance)));
#ifdef DEBUG1
          fprintf(GlobalData->stderrc, "\tgap %d: (%d,%d)\n",
                  contigCount, (int) gap.gap_length, (int) gap.gap_var);
#endif
          AppendVA_Scaffold_Gap(ss->pools->gapPool, &gap);
        }
    
      // compute distance measured from left end of scaffoldB
      // with edge slop factored into scaffoldA's distances
      if(isA)
        {
          ce.minCoord = thinEdge + thisLeft.mean -
            INTERLEAVE_CUTOFF * sqrt((double) thisLeft.variance) -
            scaffold->bpLength.mean;
          ce.maxCoord = thickEdge + thisRight.mean +
            INTERLEAVE_CUTOFF * sqrt((double) thisRight.variance) -
            scaffold->bpLength.mean;
        
          /*
            compute min & max ahang - contigs are iterated left to right
            ahangs are in or off the left end of scaffoldA

            contig or gap may specify an ahang
            NOTE: this is a simple initial implementation - should be pretty close
          */
          if(ce.minCoord <= 0 && ce.maxCoord >= CGW_DP_MINLEN)
            {
              // the interval this contig is in specifies an ahang range
              ss->bandBeg = MIN(ss->bandBeg, .5 +
                                thisLeft.mean + MAX(0,
                                                    contig->bpLength.mean -
                                                    ce.maxCoord));
              ss->bandEnd = MAX(ss->bandEnd, .5 +
                                thisLeft.mean + MIN(-ce.minCoord,
                                                    contig->bpLength.mean -
                                                    CGW_DP_MINLEN));
            }
          if(contigCount > 0)
            {
              Scaffold_Gap * gap =
                GetVA_Scaffold_Gap(ss->pools->gapPool, contigCount - 1);
              float minGapCoord = thinEdge + lastRight.mean -
                INTERLEAVE_CUTOFF * sqrt((double) lastRight.variance) -
                scaffold->bpLength.mean + 0.5;
              float maxGapCoord = thickEdge + thisLeft.mean +
                INTERLEAVE_CUTOFF * sqrt((double) thisLeft.variance) -
                scaffold->bpLength.mean + 0.5;
              if(minGapCoord < 0 && maxGapCoord > 0)
                {
                  // the interval the last gap is in specifies an ahang range
                  ss->bandBeg = MIN(ss->bandBeg, .5 +
                                    lastRight.mean + MAX(0,
                                                         gap->gap_length - maxGapCoord));
                  ss->bandEnd = MAX(ss->bandEnd, .5 +
                                    lastRight.mean + MIN(-minGapCoord, gap->gap_length));
                }
            }
        }
      else
        {
          ce.minCoord = thisLeft.mean - INTERLEAVE_CUTOFF *
            sqrt((double) thisLeft.variance);
          ce.maxCoord = thisRight.mean + INTERLEAVE_CUTOFF *
            sqrt((double) thisRight.variance);
        }

      /*
        store contig element for detecting contig overlaps and post-processing
        range must intersect overlap interval - defined above
      */
      if(ce.maxCoord >= osMin && ce.minCoord <= osMax)
        AppendVA_ContigElement(ss->edgeContigs, &ce);
      AppendVA_ContigElement(ss->contigs, &ce);

#ifdef DEBUG1
      PrintContigElement(GlobalData->stderrc, &ce);
      if(isA)
        fprintf(GlobalData->stderrc, "  minAHang: " F_COORD ", maxAHang:" F_COORD "\n",
                ss->bandBeg, ss->bandEnd);
#endif
    
      // store all contigs for Align_Scaffold interface
      tig.length = (int) (contig->bpLength.mean + .5);
      tig.lft_end = (int) (thisLeft.mean + .5);
      tig.insert_pnt = 0;
      AppendVA_Scaffold_Tig(ss->pools->tigPool, &tig);
    
      contigCount++;
      lastRight = thisRight;
    }

  // make sure minAHang & maxAHang are acceptable
  if(isA)
    {
      if(ss->bandBeg > ss->bandEnd)
        {
          ss->bandBeg = scaffold->bpLength.mean - INTERLEAVE_CUTOFF *
            sqrt((double) scaffold->bpLength.variance) - thickEdge;
          ss->bandEnd = scaffold->bpLength.mean + INTERLEAVE_CUTOFF *
            sqrt((double) scaffold->bpLength.variance) - thinEdge;
        }
      ss->bandBeg = MAX(ss->bandBeg, -(osLength.mean - CGW_DP_MINLEN));
      ss->bandEnd = MIN(ss->bandEnd, scaffold->bpLength.mean - CGW_DP_MINLEN);
    }
  
  // populate Scaffold structure members
  ss->scaffold->length = scaffold->bpLength.mean;
  ss->scaffold->packed_seq = NULL;
  ss->scaffold->num_gaps = GetNumVA_Scaffold_Gap(ss->pools->gapPool);
  ss->scaffold->gaps = GetVA_Scaffold_Gap(ss->pools->gapPool, 0);
  ss->scaffold->ctgs = GetVA_Scaffold_Tig(ss->pools->tigPool, 0);

  // above, the determination of ss->bandBeg and ss->bandEnd is pretty
  // elaborate, and seemed at one point to be the cause of some trouble with
  // some improperly (un)connected scaffolds; a much simpler computation 
  // that often suffices for Align_Scaffolds(), which itself does some
  // adjustments similar to the above, is the following simpler version.
  // All that said, we have made some progress towards cleaning up the 
  // problem with disconnected scaffolds, so perhaps it is worth trying
  // to use the version above once again.  -- ALH, 9/28/04

#define SIMPLE_BAND_COMPUTE
#ifdef SIMPLE_BAND_COMPUTE
  if(isA){
    ss->bandBeg =
      scaffold->bpLength.mean 
      + sEdge->distance.mean 
      - INTERLEAVE_CUTOFF * sqrt((double) sEdge->distance.variance);
    if(ss->bandBeg > scaffold->bpLength.mean)
      ss->bandBeg = scaffold->bpLength.mean;
    if(ss->bandBeg < -osLength.mean)
      ss->bandBeg = -osLength.mean;

    ss->bandEnd =
      scaffold->bpLength.mean 
      + sEdge->distance.mean 
      + INTERLEAVE_CUTOFF * sqrt((double) sEdge->distance.variance);
    if(ss->bandEnd > scaffold->bpLength.mean)
      ss->bandEnd = scaffold->bpLength.mean;
    if(ss->bandEnd < -osLength.mean)
      ss->bandEnd = -osLength.mean;
    assert(ss->bandBeg <= ss->bandEnd);

    fprintf(stderr,
            "Setting band [%d,%d] in A scaffold (PopulateScaffoldStuff())\n",
            ss->bandBeg,ss->bandEnd);
  }
#endif
}

Overlap * LookForChunkOverlapFromContigElements(ContigElement * ceA,
                                                ContigElement * ceB,
                                                SEdgeT * sEdge)
{
  NodeOrient orientA;
  NodeOrient orientB;
  static Overlap myOverlap;
  Overlap * retOverlap = NULL;
  CDS_COORD_t minOverlap;
  CDS_COORD_t maxOverlap;
  ChunkOverlapCheckT chunkOverlap;
  ChunkOrientationType overlapOrient;
  ChunkInstanceT * contigA;
  ChunkInstanceT * contigB;
  CDS_COORD_t minLengthA;
  CDS_COORD_t maxLengthA;
  CDS_COORD_t minLengthB;
  CDS_COORD_t maxLengthB;

  contigA = GetGraphNode(ScaffoldGraph->RezGraph, ceA->id);
  contigB = GetGraphNode(ScaffoldGraph->RezGraph, ceB->id);
  
  minLengthA = contigA->bpLength.mean - 5. * sqrt(contigA->bpLength.variance);
  maxLengthA = contigA->bpLength.mean + 5. * sqrt(contigA->bpLength.variance);
  minLengthB = contigB->bpLength.mean - 5. * sqrt(contigB->bpLength.variance);
  maxLengthB = contigB->bpLength.mean + 5. * sqrt(contigB->bpLength.variance);

  // contig orientation already factored into the ContigElement
  orientA = ceA->orient; 
  orientB = ceB->orient; 

#ifdef DEBUG1
  fprintf(stderr,"scaffold orientation: %c ; contigA %c contigB %c\n",
          sEdge->orient,orientA,orientB); 
 
  // now take care of orientation of contigs within scaffold 
  if(contigA->offsetAEnd.mean > contigA->offsetBEnd.mean){ 
    fprintf(stderr,"Reverse orientation for " F_CID " in scf\n",contigA->id); 
  } 
  if(contigB->offsetAEnd.mean > contigB->offsetBEnd.mean){ 
    fprintf(stderr,"Reverse orientation for " F_CID " in scf\n",contigB->id); 
  } 
#endif
  
  /*
    consider two orientations of overlap
    contigA:   --------
    contigB:       --------

    contigA:       --------
    contigB:   --------

    and B containing A
    contigA:       ----
    contigB:   ------------

    can also have A contain B when the ends are close
    contigA:   ------------
    contigB:          ----
  */


  // first, consider min coord of A to the left of min coord of ceB, 

  // WAS:  if(ceA->minCoord < ceB->maxCoord - minLengthB)
  if(ceA->minCoord+minLengthA < ceB->maxCoord )
    {
      if(orientA == A_B)
        overlapOrient = (orientB == A_B) ? AB_AB : AB_BA;
      else
        overlapOrient = (orientB == A_B) ? BA_AB : BA_BA;
    
      // min overlap: push ceA far to left wrt ceB
      minOverlap = MAX(CGW_MISSED_OVERLAP,
                       ceA->minCoord + minLengthA -
                       (ceB->maxCoord - minLengthB) + .5);
      minOverlap = MIN(minLengthA, MIN(minLengthB, minOverlap));
    
      // max overlap: push ceA far to right up to leftmost right end of ceB
      // this allows for overlaps up to maxLengthA ... IF ceA can be pushed
      // far enough ... but not further given constraint that ceA stick
      // out to the left.  Oh, ... but actually, we want to handle ceB 
      // contains ceA here too, which means that we just want to disallow
      // the case of ceA sticking out to the right.  So, slippage allowing,
      // we allow overlap up to the length of the longer of A or B.
      maxOverlap = MIN(MAX(maxLengthA,maxLengthB),ceA->maxCoord-ceB->minCoord+.5);
      maxOverlap = MAX(CGW_MISSED_OVERLAP, maxOverlap);

      chunkOverlap = OverlapChunks(ScaffoldGraph->RezGraph,
                                   ceA->id, ceB->id, overlapOrient,
                                   minOverlap, maxOverlap,
                                   CGW_DP_ERATE, FALSE);
    
#ifdef PRINT_OVERLAPS    
      fprintf(stderr,"A_B: Trying to find overlap between " F_CID " and "
              F_CID " overlap range " F_COORD "," F_COORD
              " orientation %c --> " F_COORD "\n",
              ceA->id,ceB->id,minOverlap,maxOverlap,
              overlapOrient,chunkOverlap.overlap);
#endif
    
      if(chunkOverlap.overlap != 0)
        {
          if(chunkOverlap.AContainsB || chunkOverlap.BContainsA)
            {
              // if both contigs are reverse(anti-normal) then
              // OverlapChunks switches them and makes b a so switch back
              if (chunkOverlap.spec.cidA == ceB->id)
                { 
                  myOverlap.begpos = chunkOverlap.bhg;
                  myOverlap.endpos = chunkOverlap.ahg;
                } else {
                  myOverlap.begpos = chunkOverlap.ahg;
                  myOverlap.endpos = chunkOverlap.bhg;
                }
            }
          else
            {
              myOverlap.begpos = ceA->length - chunkOverlap.overlap;
              myOverlap.endpos = ceB->length - chunkOverlap.overlap;
            }
          myOverlap.length = chunkOverlap.overlap;
          myOverlap.diffs = 0;
          retOverlap = &myOverlap;

#ifdef PRINT_OVERLAPS    
          fprintf(stderr,"B_A: Trying to find overlap between " F_CID " and "
                  F_CID " overlap range " F_COORD "," F_COORD
                  " orientation %c --> " F_COORD "\n",
                  ceA->id,ceB->id,minOverlap,maxOverlap,
                  overlapOrient,chunkOverlap.overlap);
#endif

        }
    }

  // consider ceA to the right of ceB
  if(retOverlap == NULL && ceB->minCoord + minLengthB < ceA->maxCoord)
    {
      if(orientA == A_B)
        overlapOrient = (orientB == A_B) ? BA_BA : BA_AB;
      else
        overlapOrient = (orientB == A_B) ? AB_BA : AB_AB;
    
      // min overlap: push ceA far to the right wrt ceB
      minOverlap = MAX(CGW_MISSED_OVERLAP,
                       ceB->minCoord + minLengthB -
                       (ceA->maxCoord - minLengthA) + .5);
      minOverlap = MIN(minLengthA, MIN(minLengthB, minOverlap));
      //max overlap: push ceA far to left up to rightmost left end of ceB
      // this allows for overlaps up to maxLengthA ... IF ceA can
      // be pushed far enough ... but not further because we restrict
      // ourselves to ceA sticking out further to the right than ceB;
      // but again, we actually can allow for containment, so allow
      // max length of A or B
      maxOverlap = MIN(MAX(maxLengthA,maxLengthB),ceA->maxCoord-ceB->minCoord+.5);
      maxOverlap = MAX(CGW_MISSED_OVERLAP, maxOverlap);
    
      chunkOverlap = OverlapChunks(ScaffoldGraph->RezGraph,
                                   ceA->id, ceB->id, overlapOrient,
                                   minOverlap, maxOverlap,
                                   CGW_DP_ERATE, FALSE);

#ifdef PRINT_OVERLAPS    
      fprintf(stderr,"B_A: Trying to find overlap between " F_CID " and "
              F_CID " overlap range " F_COORD "," F_COORD
              " orientation %c --> " F_COORD "\n",
              ceA->id,ceB->id,minOverlap,maxOverlap,
              overlapOrient,chunkOverlap.overlap);
#endif

      if(chunkOverlap.overlap != 0)
        {
          if(chunkOverlap.AContainsB || chunkOverlap.BContainsA)
            {
              if (chunkOverlap.spec.cidA == ceB->id)
                { 
                  myOverlap.begpos = chunkOverlap.bhg;
                  myOverlap.endpos = chunkOverlap.ahg;
                } else {
                  myOverlap.begpos = chunkOverlap.ahg;
                  myOverlap.endpos = chunkOverlap.bhg;
                }
              myOverlap.length = chunkOverlap.overlap;
            }
          else
            {

              //     non-canonical return:
              //
              //     -------        B
              //        -------     A
              //

	
              myOverlap.begpos = -(ceB->length - chunkOverlap.overlap);
              myOverlap.endpos = -(ceA->length - chunkOverlap.overlap);
              myOverlap.length = chunkOverlap.overlap - myOverlap.begpos - myOverlap.endpos;
            }
          myOverlap.diffs = 0;
          retOverlap = &myOverlap;
        }
    }

  // consider ceB containing ceA? ... No, we did that above ...

#ifdef PRINT_OVERLAPS
  if(retOverlap != NULL)
    fprintf(GlobalData->stderrc, "Contig overlap (scaff:" F_CID ", contig:" F_CID ") & (scaff:" F_CID ", contig " F_CID "): ahg:" F_COORD ", bhg:" F_COORD ", length:" F_COORD ", quality:%.2f orient:%s\n",
            sEdge->idA, ceA->id, sEdge->idB, ceB->id,
            retOverlap->begpos, retOverlap->endpos, retOverlap->length,
            (float) retOverlap->diffs / (float) retOverlap->length,
            overlapOrient == AB_AB ? "AB_AB" :
            (overlapOrient == AB_BA ? "AB_BA" :
             (overlapOrient == BA_AB ? "BA_AB" : "BA_BA")));
#endif
  return retOverlap;
}


int PopulateScaffoldAlignmentInterface(CIScaffoldT * scaffoldA,
                                       CIScaffoldT * scaffoldB,
                                       SEdgeT * sEdge,
                                       ScaffoldAlignmentInterface * sai)
{
  int indexA;
  CDS_CID_t idA;
  CDS_CID_t idB;
  ChunkOrientationType orient;

  if(sEdge->distance.mean -
     INTERLEAVE_CUTOFF * sqrt((double) sEdge->distance.variance) >= 0.0)
    {
      fprintf(GlobalData->stderrc, "PopulateScaffoldAlignmentInterface called with non-negative edge!\n");
      return 1;
    }

  idA = sEdge->idA;
  idB = sEdge->idB;
  orient = sEdge->orient;

  // adjust sEdge for ease of use - change back later
  // assume scaffoldA is to 'left' of scaffoldB (i.e., non-negative ahang)
  if(sEdge->idA != scaffoldA->id)
    {
      sEdge->idA = idB;
      sEdge->idB = idA;
      sEdge->orient = FlipEdgeOrient(sEdge->orient);
    }

  if(sEdge->distance.mean -
     INTERLEAVE_CUTOFF * sqrt((double) sEdge->distance.variance) >= 0.0)
    {
      fprintf(GlobalData->stderrc, "PopulateScaffoldAlignmentInterface called with non-negative edge!\n");
      return 1;
    }

  // need to do this before adding contig overlaps
  ResetScaffoldAlignmentInterface(sai);
  PopulateScaffoldStuff(sai->scaffoldA, scaffoldA, sEdge, scaffoldB->bpLength);
  PopulateScaffoldStuff(sai->scaffoldB, scaffoldB, sEdge, scaffoldA->bpLength);

  // now look for overlaps between edgeContigs
  // loop over all edgeContigs in scaffoldA
  // loop from edge out, to catch low variance overlaps first
  // this may catch missing, required overlaps to quit quickly
  for(indexA = GetNumVA_ContigElement(sai->scaffoldA->edgeContigs) - 1;
      indexA >= 0;
      indexA--)
    {
      ContigElement * ceA =
        GetVA_ContigElement(sai->scaffoldA->edgeContigs, indexA);
      int indexB;

      // loop over all edgeContigs in scaffloldB
      for(indexB = 0;
          indexB < GetNumVA_ContigElement(sai->scaffoldB->edgeContigs);
          indexB++)
        {
          ContigElement * ceB =
            GetVA_ContigElement(sai->scaffoldB->edgeContigs, indexB);

          // look for overlap if they intersect
          if(ceA->maxCoord >= ceB->minCoord && ceB->maxCoord >= ceA->minCoord)
            {
              Overlap * overlap = NULL;
              overlap = LookForChunkOverlapFromContigElements(ceA, ceB, sEdge);
              if(overlap)
                {
                  // MODIFICATIONS Nov 17 2003 by ALH: 
                  // ARRGH!!!: Align_Scaffold() assumes
                  // that the seglist is actually in reverse order ... that is,
                  // overlaps between later contigs precede overlaps between
                  // earlier contigs ... so, we need to push new overlaps onto 
                  // the head rather than tail
                  if(sai->segmentList == NULL)
                    {
                      sai->segmentList = safe_calloc(1, sizeof(Segment));
                      assert(sai->segmentList != NULL);
                      sai->segmentList->next = NULL;
                    }
                  else
                    {
                      Segment *s =  safe_calloc(1, sizeof(Segment));
                      assert(s!=NULL);
                      s->next=sai->segmentList;
                      sai->segmentList=s;
                    }
                  sai->segmentList->overlap = safe_calloc(1, sizeof(Local_Overlap));
                  assert(sai->segmentList->overlap != NULL);
                  sai->segmentList->overlap->begpos = overlap->begpos;
                  sai->segmentList->overlap->endpos = overlap->endpos;
                  sai->segmentList->overlap->length = overlap->length;
                  sai->segmentList->a_contig = ceA->index;
                  sai->segmentList->b_contig = ceB->index;

                  sai->numSegs++;


                  // The following are used only for the XFIG diagrams produced by CA_ALN_scafcomp routines;
                  // there are four cases to be concerned with:

                  //            -------    A     alow = begpos; ahgh = begpos + length 
                  //                -----  B     blow = 0 ; bhgh = length
                  //  
                  //            -------    A     alow = 0; ahgh = Alen + endpos = length + begpos + endpos 
                  //         -----         B     blow = -begpos ; bhgh = length + endpos
                  //
                  //         ------------  A     alow = begpos; ahgh = length + begpos + endpos
                  //            -----      B     blow = 0 = ?; bhgh = length + endpos
                  // 
                  //            -----      A     alow = 0; ahgh = length + begpos
                  //         ------------  B     blow = -begpos; bhgh = length
                  //

                  sai->segmentList->alow = MAX(0,overlap->begpos);  
                  sai->segmentList->blow = MAX(0,-(overlap->begpos));  
                  sai->segmentList->ahgh = overlap->length + overlap->begpos + MIN(0,overlap->endpos);
                  sai->segmentList->bhgh = overlap->length + MIN(0,overlap->endpos);

                }
              else
                {
#ifdef DEBUG1
                  fprintf(GlobalData->stderrc,
                          "There is no overlap between (" F_CID ":"
                          F_COORD "," F_COORD ") and (" F_CID ":"
                          F_COORD "," F_COORD ")\n",
                          ceA->id, ceA->minCoord, ceA->maxCoord,
                          ceB->id, ceB->minCoord, ceB->maxCoord);
#endif
                }
            }
          else
            {
              assert(ceA->maxCoord < ceB->minCoord + CGW_DP_MINLEN + 1||
                     ceB->maxCoord < ceA->minCoord + CGW_DP_MINLEN + 1);
#ifdef DEBUG1
              fprintf(GlobalData->stderrc,
                      "Not looking for overlap between (" F_CID ":"
                      F_COORD "," F_COORD ") and (" F_CID ":"
                      F_COORD "," F_COORD ")\n",
                      ceA->id, ceA->minCoord, ceA->maxCoord,
                      ceB->id, ceB->minCoord, ceB->maxCoord);
#endif
            }
        }
    }

  // change sEdge back to what it was
  sEdge->idA = idA;
  sEdge->idB = idB;
  sEdge->orient = orient;

#ifdef DEBUG1
  PrintScaffoldAlignmentInterface(GlobalData->stderrc, sai);
  fprintf(GlobalData->stderrc, "Scaffold A CGW data structure:\n");
  DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffoldA, FALSE);
  DumpACIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffoldA, FALSE);

  fprintf(GlobalData->stderrc, "Scaffold B CGW data structure:\n");
  DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffoldB, FALSE);
  DumpACIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffoldB, FALSE);

  fprintf(GlobalData->stderrc, "\nsEdge cgw data structure:\n");
  PrintSEdgeT(GlobalData->stderrc, ScaffoldGraph, "sEdge", sEdge, scaffoldA->id);
#endif

  return 0;
}


static int segmentCompare(const Segment * a, const Segment * b)
{
  if(a->a_contig == b->a_contig)
    return(a->b_contig - b->b_contig);
  else
    return(a->a_contig - b->a_contig);
}


/*
  compute expansion of gaps[index] given current placement of
  contigs[index] and contigs[index+1]
*/
CDS_COORD_t ComputeCurrentGapExpansion(Scaffold_Tig * contigs,
                                       Scaffold_Gap * gaps,
                                       int index)
{
  return MAX(0, contigs[index+1].lft_end - gaps[index].gap_length -
             contigs[index].lft_end - contigs[index].length);
}


/*
  compute expansion of gap between contigs1[index1-1] and contigs1[index1]
  needed to accomodate contigs2[index2]
  assumes we're proceeding left to right
*/
CDS_COORD_t ComputeAdditionalLeftGapExpansion(Scaffold_Tig * contigs1,
                                              int index1,
                                              Scaffold_Tig * contigs2,
                                              int index2)
{
  return MAX(0, contigs2[index2].length + 2 * MIN_GAP_LENGTH -
             (contigs1[index1].lft_end - MAX(contigs1[index1-1].lft_end +
                                             contigs1[index1-1].length,
                                             contigs2[index2-1].lft_end +
                                             contigs2[index2-1].length)));
}

/*
  compute expansion of gap between contigs1[index1] and contigs1[index1+1]
  needed to accomodate contigs2[index2]
  assumes we're proceeding right to left
*/
CDS_COORD_t ComputeAdditionalRightGapExpansion(Scaffold_Tig * contigs1,
                                               int index1,
                                               Scaffold_Tig * contigs2,
                                               int index2)
{
  return MAX(0, contigs2[index2].length + 2 * MIN_GAP_LENGTH -
             (MIN(contigs1[index1+1].lft_end, contigs2[index2+1].lft_end) -
              (contigs1[index1].lft_end + contigs1[index1].length)));

}


/*
  compression of gaps1[index1-1] to fit contigs1[index1] to the left of
  contigs2[index2]
*/
CDS_COORD_t ComputeLeftGapCompression(Scaffold_Tig * contigs1,
                                      Scaffold_Gap * gaps1,
                                      int index1,
                                      Scaffold_Tig * contigs2,
                                      int index2,
                                      CDS_COORD_t gap2Expansion)
{
  return MAX(0, (contigs1[index1-1].lft_end + contigs1[index1-1].length +
                 gaps1[index1-1].gap_length) -
             (contigs2[index2].lft_end + gap2Expansion - MIN_GAP_LENGTH -
              contigs1[index1].length));
}

/*
  compression of gaps1[index1] to fit contigs1[index1] to the right of
  contigs2[index2]
*/
CDS_COORD_t ComputeRightGapCompression(Scaffold_Tig * contigs1,
                                       Scaffold_Gap * gaps1,
                                       int index1,
                                       Scaffold_Tig * contigs2,
                                       int index2,
                                       CDS_COORD_t gap2Expansion)
{
  return MAX(0, (contigs2[index2].lft_end + contigs2[index2].length +
                 MIN_GAP_LENGTH - gap2Expansion) -
             (contigs1[index1+1].lft_end - gaps1[index1].gap_length -
              contigs1[index1].length));
}


int AdjustNonOverlappingContigsLeftToRight(Scaffold_Tig * contigsA,
                                           Scaffold_Gap * gapsA,
                                           int minIndexA,
                                           int maxIndexA,
                                           Scaffold_Tig * contigsB,
                                           Scaffold_Gap * gapsB,
                                           int minIndexB,
                                           int maxIndexB)
{
  int ia;
  int ib;
  // Detect unwanted overlaps & adjust locally as we go
  // never move ia-1 or ib-1
  ia = minIndexA;
  ib = minIndexB;

  assert(ia != 0 && ib != 0);
  
  while(ia <= maxIndexA || ib <= maxIndexB)
    {
      if(ib <= maxIndexB)
        {
          // tentatively place contigsB[ib]
          contigsB[ib].lft_end =
            MAX(contigsB[ib-1].lft_end + contigsB[ib-1].length +
                gapsB[ib-1].gap_length,
                contigsA[ia-1].lft_end + contigsA[ia-1].length + MIN_GAP_LENGTH);
          if(ia > maxIndexA)
            {
              // done with all contigsA
              ib++;
              continue;
            }
        }
    
      if(ia <= maxIndexA)
        {
          // tentatively place contigsA[ia]
          contigsA[ia].lft_end =
            MAX(contigsA[ia-1].lft_end + contigsA[ia-1].length +
                gapsA[ia-1].gap_length,
                contigsB[ib-1].lft_end + contigsB[ib-1].length + MIN_GAP_LENGTH);
          if(ib > maxIndexB)
            {
              // done with all contigsB
              ia++;
              continue;
            }
        }
    
      /*
        if we're here, ia & ib are still in play
        and they're completely to the right of all contigs ia-1 & ib-1
        choose one to go to the right of the other
        for simplicity, consider just the gap to the left of ia & ib
        Case 1:
        Contigia-1  gapia-1    contigia
        ScaffoldA          --------           -----------------
      
        ScaffoldB       ---------            ----------------
        Contigib-1   gapib-1     contigib
      */
      if(contigsB[ib].lft_end + contigsB[ib].length < contigsA[ia].lft_end)
        ib++;
      else if(contigsA[ia].lft_end + contigsA[ia].length < contigsB[ib].lft_end)
        ia++;
      else
        {
          float preExpansionA;  // expansion of gap A already
          float addExpansionA;  // additional expansion of gap A to fit contig B
          float compressionA;   // compression of gap A to fit contig A in gap B
          float preExpansionB;
          float addExpansionB;
          float compressionB;
          // yes, gap_var is actually used as stddev in Align_Scaffold()
          float stddevGapA = gapsA[ia-1].gap_var; 
          float stddevGapB = gapsB[ib-1].gap_var;
          float deltaAB;
          float deltaBA;
      
          // expansion of gaps from tentative placement
          preExpansionA = ComputeCurrentGapExpansion(contigsA, gapsA, ia - 1);
          preExpansionB = ComputeCurrentGapExpansion(contigsB, gapsB, ib - 1);
      
          // expansion of gaps to accomodate other contig fitting
          // contigsB[ib] in gapsA[ia-1]
          //
          addExpansionA = ComputeAdditionalLeftGapExpansion(contigsA, ia, contigsB, ib);
          addExpansionB = ComputeAdditionalLeftGapExpansion(contigsB, ib, contigsA, ia);


          //  compression of gapsX[ix-1] to put contigsX[ix] to the left
          //  of contigsY[iy]
          //
          compressionA = ComputeLeftGapCompression(contigsA, gapsA, ia, contigsB, ib, addExpansionB);
          compressionB = ComputeLeftGapCompression(contigsB, gapsB, ib, contigsA, ia, addExpansionA);

          // deltaAB is stddevs involved in placing A to the left of B
          deltaAB = (preExpansionB + addExpansionB) / stddevGapB +
            compressionA / stddevGapA;
          deltaBA = (preExpansionA + addExpansionA) / stddevGapA +
            compressionB / stddevGapB;
          if(deltaAB > deltaBA)
            {
              // place contigsB[ib] in gapsA[ia-1]
              contigsB[ib].lft_end =
                MAX(contigsA[ia-1].lft_end + contigsA[ia-1].length + MIN_GAP_LENGTH,
                    contigsB[ib-1].lft_end + contigsB[ib-1].length + gapsB[ib-1].gap_length);
              ib++;
            }
          else
            {
              contigsA[ia].lft_end =
                MAX(contigsB[ib-1].lft_end + contigsB[ib-1].length + MIN_GAP_LENGTH,
                    contigsA[ia-1].lft_end + contigsA[ia-1].length + gapsA[ia-1].gap_length);
              ia++;
            }
        }
    }
  return 0;
}


void PlaceContigsBetweenOverlapSets(COSData * cosLeft,
                                    COSData * cosRight,
                                    Scaffold_Tig * contigsA,
                                    Scaffold_Gap * gapsA,
                                    Scaffold_Tig * contigsB,
                                    Scaffold_Gap * gapsB)
{
  /*
    Identify rightmost contigA, contigB in cosLeft
    these are already placed
    Identify leftmost contigA, contigB in cosRight

    Iterate over contigs between contig overlap sets & place them
  */
  int ia, ib;
  float sumStddevsA = 0;
  float sumStddevsB = 0;

  /*
    Make initial placement of contigs to determine difference between
    scafflolds in distance between overlap sets
  */
  for(ia = cosLeft->a.maxIndex + 1; ia < cosRight->a.minIndex; ia++)
    {
      contigsA[ia].lft_end = contigsA[ia-1].lft_end +
        contigsA[ia-1].length + gapsA[ia-1].gap_length;
      // yes, gap_var is actually used as stddev in Align_Scaffold()
      sumStddevsA += gapsA[ia-1].gap_var;
    }
  for(ib = cosLeft->b.maxIndex + 1; ib < cosRight->b.minIndex; ib++)
    {
      contigsB[ib].lft_end = contigsB[ib-1].lft_end +
        contigsB[ib-1].length + gapsB[ib-1].gap_length;
      // yes, gap_var is actually used as stddev in Align_Scaffold()
      sumStddevsB += gapsB[ib-1].gap_var;
    }

  //Compute expansion of gaps in one of the scaffolds
  {
    CDS_COORD_t initialPositionA;
    CDS_COORD_t initialPositionB;
    CDS_COORD_t spreadInA; // = negative spread in scaffold B

    initialPositionA = contigsA[cosRight->a.minIndex-1].lft_end +
      contigsA[cosRight->a.minIndex-1].length +
      gapsA[cosRight->a.minIndex-1].gap_length;
    // yes, gap_var is actually used as stddev in Align_Scaffold()
    sumStddevsA += gapsA[cosRight->a.minIndex-1].gap_var;

    initialPositionB = contigsB[cosRight->b.minIndex-1].lft_end +
      contigsB[cosRight->b.minIndex-1].length +
      gapsB[cosRight->b.minIndex-1].gap_length;
    // yes, gap_var is actually used as stddev in Align_Scaffold()
    sumStddevsB += gapsB[cosRight->b.minIndex-1].gap_var;
    
    /*
      initialPositionA - initialPositionB is distance A is to the right of B
      by tentative placement from gap spacing
      a.lft_end - b.lft_end is distance A is to right of B by overlap
      if(former > latter)
      then spread is in scaffoldA
      else
      spread is in scaffoldB
    */

    spreadInA = (initialPositionA - initialPositionB) -
      (contigsA[cosRight->a.minIndex].lft_end -
       contigsB[cosRight->b.minIndex].lft_end);
    if(spreadInA > 0)
      {
        for(ia = cosLeft->a.maxIndex + 1; ia < cosRight->a.minIndex; ia++)
          {
            // warning - changing gap_length field in the gap structure
            // yes, gap_var is actually used as stddev in Align_Scaffold()
            gapsA[ia-1].gap_length += spreadInA * gapsA[ia-1].gap_var / sumStddevsA;
            contigsA[ia].lft_end = contigsA[ia-1].lft_end + contigsA[ia-1].length + gapsA[ia-1].gap_length;
          }
      }
    else
      {
        for(ib = cosLeft->b.maxIndex + 1; ib < cosRight->b.minIndex; ib++)
          {
            // warning - changing gap_length field in the gap structure
            // yes, gap_var is actually used as stddev in Align_Scaffold()
            gapsB[ib-1].gap_length -= spreadInA * gapsB[ib-1].gap_var / sumStddevsB;
            contigsB[ib].lft_end = contigsB[ib-1].lft_end + contigsB[ib-1].length + gapsB[ib-1].gap_length;
          }
      }
  }

  AdjustNonOverlappingContigsLeftToRight(contigsA, gapsA,
                                         cosLeft->a.maxIndex + 1,
                                         cosRight->a.minIndex - 1,
                                         contigsB, gapsB,
                                         cosLeft->b.maxIndex + 1,
                                         cosRight->b.minIndex - 1);
}


int PlaceContigsInOverlapSet(COSData * cos,
                             Scaffold_Tig * contigsA,
                             Scaffold_Gap * gapsA,
                             Scaffold_Tig * contigsB,
                             Scaffold_Gap * gapsB)
{
  int ia = cos->a.minIndex;
  int ib = cos->b.minIndex;
  CDS_COORD_t offset;

  if(contigsA[ia].lft_end == 0)
    {
      // contigA is leftmost, place it to get offset
      if(ia == 0)
        offset = 0;
      else if(ib > 0)
        offset = MAX(contigsA[ia-1].lft_end + contigsA[ia-1].length +
                     gapsA[ia-1].gap_length,
                     contigsB[ib-1].lft_end + contigsB[ib-1].length +
                     MIN_GAP_LENGTH);
      else
        offset = contigsA[ia-1].lft_end + contigsA[ia-1].length +
          gapsA[ia-1].gap_length;
    }
  else
    {
      //  You may be able to get around this assert by undef
      //  ALLOW_NEG_GAP_BACKUP in CA_ALN_scafcomp.c.  You may also want
      //  to restart from scratch instead of a checkpoint after doing
      //  that.
      //
      assert(contigsB[ib].lft_end == 0);
      if(ib == 0)
        offset = 0;
      else if(ia > 0)
        offset = MAX(contigsB[ib-1].lft_end + contigsB[ib-1].length +
                     gapsB[ib-1].gap_length,
                     contigsA[ia-1].lft_end + contigsA[ia-1].length +
                     MIN_GAP_LENGTH);
      else
        offset = contigsB[ib-1].lft_end + contigsB[ib-1].length +
          gapsB[ib-1].gap_length;
    }

  for(ia = cos->a.minIndex; ia <= cos->a.maxIndex; ia++)
    contigsA[ia].lft_end += offset;
  for(ib = cos->b.minIndex; ib <= cos->b.maxIndex; ib++)
    contigsB[ib].lft_end += offset;
  
  return 0;
}

#define MIN_OVERLAP_SET             0
#define NO_OVERLAP_SET             -1
#define SKIPPED_CONTIG             -2
#define SWITCHED_CONTIG            -3


int MarkSkippedContigsInOverlapSet(COSData * cos,
                                   Scaffold_Tig * contigs,
                                   int isA)
{
  int i;
  int numSkipped = 0;
  ContigSetInterval * csi = isA ? &(cos->a) : &(cos->b);
  
  for(i = csi->minIndex + 1; i < csi->maxIndex; i++)
    {
      if(contigs[i].insert_pnt != cos->index)
        {
          assert(contigs[i].insert_pnt == NO_OVERLAP_SET);
          //fprintf(GlobalData->stderrc, "** Contig out of order in overlap set!\n");
          contigs[i].insert_pnt = SKIPPED_CONTIG;
          numSkipped++;
        }
    }
  return numSkipped;
}


int * CreateContigMap(Scaffold_Tig * contigs, int numContigs)
{
  int * contigMap = safe_malloc(numContigs * sizeof(int));

  if(contigMap)
    {
      int i;
      for(i = 0; i < numContigs; i++)
        contigMap[i] = i;
    }
  return contigMap;
}


void InitializeContigSetInterval(ContigSetInterval * csi,
                                 Scaffold_Tig * contigs,
                                 int contigIndex)
{
  csi->minIndex = csi->maxIndex = contigIndex;
  // initialize to max possible, Update function will adjust
  csi->minCoord = contigs[contigIndex].length;
  // initialize to min possible, Update function will adjust
  csi->maxCoord = -contigs[contigIndex].length ;
}


void InitializeCOSDataItem(COSData * cos,
                           Segment * overlaps,
                           int overlapIndex,
                           Scaffold_Tig * contigsA,
                           Scaffold_Tig * contigsB)
{
  cos->firstOverlap = cos->lastOverlap = overlapIndex;
  InitializeContigSetInterval(&(cos->a),
                              contigsA,
                              overlaps[overlapIndex].a_contig);
  InitializeContigSetInterval(&(cos->b),
                              contigsB,
                              overlaps[overlapIndex].b_contig);
}


void UpdateContigSetInterval(ContigSetInterval * csi,
                             Scaffold_Tig * contigs,
                             int contigIndex)
{
  csi->minIndex = MIN(csi->minIndex, contigIndex);
  csi->maxIndex = MAX(csi->maxIndex, contigIndex);
  csi->minCoord = MIN(csi->minCoord, contigs[contigIndex].lft_end);
  csi->maxCoord = MAX(csi->maxCoord, contigs[contigIndex].lft_end +
                      contigs[contigIndex].length);
}


void UpdateCOSDataItem(COSData * cos,
                       Segment * overlaps,
                       int overlapIndex,
                       Scaffold_Tig * contigsA,
                       Scaffold_Tig * contigsB)
{
  cos->lastOverlap = overlapIndex;
  UpdateContigSetInterval(&(cos->a), contigsA,
                          overlaps[overlapIndex].a_contig);
  UpdateContigSetInterval(&(cos->b), contigsB,
                          overlaps[overlapIndex].b_contig);
}


void AdjustContigOverlapSetOffsets(COSData * cos,
                                   Scaffold_Tig * contigsA,
                                   Scaffold_Tig * contigsB,
                                   CDS_COORD_t offset)
{
  int i;
  cos->a.minCoord += offset;
  cos->a.maxCoord += offset;
  cos->b.minCoord += offset;
  cos->b.maxCoord += offset;

  for(i = cos->a.minIndex; i <= cos->a.maxIndex; i++)
    contigsA[i].lft_end += offset;
  for(i = cos->b.minIndex; i <= cos->b.maxIndex; i++)
    contigsB[i].lft_end += offset;
}


// returns number of contigs that need to be reordered
int ExamineContigOverlapSets(ScaffoldAlignmentInterface * sai,
                             Segment * overlaps,
                             int numOverlaps,
                             VA_TYPE(COSData) * overlapSet)
{
  Scaffold_Tig * contigsA =
    GetVA_Scaffold_Tig(sai->scaffoldA->pools->tigPool, 0);
  int numContigsA = GetNumVA_Scaffold_Tig(sai->scaffoldA->pools->tigPool);

  Scaffold_Tig * contigsB =
    GetVA_Scaffold_Tig(sai->scaffoldB->pools->tigPool, 0);
  int numContigsB = GetNumVA_Scaffold_Tig(sai->scaffoldB->pools->tigPool);

  int i;
  int numSkippedA;
  int numSkippedB;
  COSData cos;

  // overlaps are already sorted by a_contig, b_contig
  
  /*
    An overlap set implies a new contig will be created upon scaffold merging
    Iterate through overlaps & traverse overlap sets
    Figure out how to detect implicit reorderings

    Overlap set will have following pattern of a_contig, b_contig indices

    2,1
    ----
    3,2 |
    3,3 |
    4,2 |- Overlap set
    4,3 |
    4,4 |
    ---
    6,5

    To traverse an overlap set, keep min & max contig index in A & B
    while a_contig is within minA:maxA or b_contig is within minB:maxB
    you're in the same overlap set

    Two potential problems:
    1) a contig is skipped in an overlap set
    such a contig must be moved to before or after the set,
    whichever minimizes the sum of stddevs affected
    2) overlap intervals imply a scaffold reordering within the overlap set
  */

  // overload the insert_pnt field of Scaffold_Tig for identifying
  // which overlap set each contig is in. Initialize all to -1
  //
  for(i=0; i<numContigsA; i++)
    contigsA[i].insert_pnt = NO_OVERLAP_SET;
  for(i=0; i<numContigsB; i++)
    contigsB[i].insert_pnt = NO_OVERLAP_SET;

  // loop over overlaps to identify contig overlap sets
  numSkippedA = numSkippedB = 0;
  for(cos.index = MIN_OVERLAP_SET, i = 0; i < numOverlaps; cos.index++)
    {
      InitializeCOSDataItem(&cos, overlaps, i, contigsA, contigsB);
      // while this contig is in the set, label & increment i
      for( ;
           i < numOverlaps &&
             ((overlaps[i].a_contig >= cos.a.minIndex && overlaps[i].a_contig <= cos.a.maxIndex) ||
              (overlaps[i].b_contig >= cos.b.minIndex && overlaps[i].b_contig <= cos.b.maxIndex));
           i++)
        {
          /*
            begpos = ahang
            we want the leftmost relative coordinate to be 0, so
            0      begpos      0     -begpos
            |     |            |     |
            a_contig:    ----------               ---------
            b_contig:          ---------    ---------
          */
          if(contigsA[overlaps[i].a_contig].insert_pnt == NO_OVERLAP_SET)
            {
              if(contigsB[overlaps[i].b_contig].insert_pnt == NO_OVERLAP_SET)
                {
                  // new overlap set
                  contigsA[overlaps[i].a_contig].lft_end = MAX(0, -overlaps[i].overlap->begpos);
                  contigsB[overlaps[i].b_contig].lft_end = MAX(0, overlaps[i].overlap->begpos);
                }
              else
                {
                  // b left end has been set
                  contigsA[overlaps[i].a_contig].lft_end = contigsB[overlaps[i].b_contig].lft_end - overlaps[i].overlap->begpos;
                }
            }
          else if(contigsB[overlaps[i].b_contig].insert_pnt == NO_OVERLAP_SET)
            {
              // a left end has been set
              contigsB[overlaps[i].b_contig].lft_end = contigsA[overlaps[i].a_contig].lft_end + overlaps[i].overlap->begpos;
            }
          else
            {
              /*
                both already overlap in overlap set -
                don't check consistency of coordinates
                NOTE: may want to add a check
              */
            }

          UpdateCOSDataItem(&cos, overlaps, i, contigsA, contigsB);
          overlaps[i].alow = cos.index;
          contigsA[overlaps[i].a_contig].insert_pnt = cos.index;
          contigsB[overlaps[i].b_contig].insert_pnt = cos.index;
        }

      // this may not be necessary, but we'll do it anyway
      if(cos.a.minCoord != 0 && cos.b.minCoord != 0)
        {
          // adjust cos & contig coordinates so min = 0
          AdjustContigOverlapSetOffsets(&cos, contigsA, contigsB, -MIN(cos.a.minCoord, cos.b.minCoord));
        }
      // find skipped contigs in this overlap set
      numSkippedA += MarkSkippedContigsInOverlapSet(&cos, contigsA, TRUE);
      numSkippedB += MarkSkippedContigsInOverlapSet(&cos, contigsB, FALSE);

      AppendVA_COSData(overlapSet, &cos);
    }

  return numSkippedA + numSkippedB;
}


CDS_COORD_t ComputeLeftMostLeftEnd(Scaffold_Tig * contigs1,
                                   Scaffold_Gap * gaps1,
                                   int index1,
                                   Scaffold_Tig * contigs2,
                                   int index2)
{
  return MIN(contigs1[index1+1].lft_end - gaps1[index1].gap_length,
             contigs2[index2+1].lft_end - MIN_GAP_LENGTH) -
    contigs1[index1].length;
}


void PlaceContigsLeftOfFirstOverlapSet(COSData * cos,
                                       Scaffold_Tig * contigsA,
                                       Scaffold_Gap * gapsA,
                                       Scaffold_Tig * contigsB,
                                       Scaffold_Gap * gapsB)
{
  int ia = cos->a.minIndex - 1;
  int ib = cos->b.minIndex - 1;

  // Detect unwanted overlaps & adjust locally as we go
  // don't move ia+1 or ib+1
  while(ia >= 0 || ib >= 0)
    {
      // set tentative position for contig ia and/or ib;
      if(ib >= 0)
        {
          contigsB[ib].lft_end =
            ComputeLeftMostLeftEnd(contigsB, gapsB, ib, contigsA, ia);
          if(ia < 0)
            {
              // we're done with A contigs, keep going with B contigs
              ib--;
              continue;
            }
        }

      if(ia >= 0)
        {
          contigsA[ia].lft_end =
            ComputeLeftMostLeftEnd(contigsA, gapsA, ia, contigsB, ib);
          if(ib < 0)
            {
              // we're done with B contigs, keep going with A contigs
              ia--;
              continue;
            }
        }

      /*
        if we're here, ia >= 0 && ib >= 0 and
        they're tentatively placed entirely to the left of contigs ia+1 & ib+1
        if their coordinates don't intersect, place the rightmost one.
      */
      if(contigsA[ia].lft_end >
         contigsB[ib].lft_end + contigsB[ib].length + MIN_GAP_LENGTH)
        ia--;
      else if(contigsB[ib].lft_end >
              contigsA[ia].lft_end + contigsA[ia].length + MIN_GAP_LENGTH)
        ib--;
      else
        {
          float preExpansionA;  // expansion of gap A already
          float addExpansionA;  // additional expansion of gap A to fit contig B
          float compressionA;   // compression of gap A to fit contig A in gap B
          float preExpansionB;
          float addExpansionB;
          float compressionB;
          // yes, gap_var is actually used as stddev in Align_Scaffold()
          float stddevGapA = gapsA[ia].gap_var;
          float stddevGapB = gapsB[ib].gap_var;
          float deltaAB;
          float deltaBA;

          /*
            the two contigs  overlap each other.
            Choose one to go to the right of the other
            for simplicity, consider just the gap to the right of ia & ib
            Case 1:
            contigia  gapia    contigia+1
            ScaffoldA      --------           -----------------
        
            ScaffoldB       ---------            ----------------
            contigib   gapib     contigib+1

            Make the choice based on minimal stretching of gap sizes in terms of
            stddevs
          */
      
          // calculate expansion of gaps needed to accomodate other contig
          // = (gap stretch already) + (additional stretch to fit other contig)
          preExpansionA = ComputeCurrentGapExpansion(contigsA, gapsA, ia);
          preExpansionB = ComputeCurrentGapExpansion(contigsB, gapsB, ib);

          addExpansionA = ComputeAdditionalRightGapExpansion(contigsA, ia, contigsB, ib);
          addExpansionB = ComputeAdditionalRightGapExpansion(contigsB, ib, contigsA, ia);

          // calculate compression of gaps needed to place this contig
          compressionA = ComputeRightGapCompression(contigsA, gapsA, ia, contigsB, ib, addExpansionB);
          compressionB = ComputeRightGapCompression(contigsB, gapsB, ib, contigsA, ia, addExpansionA);

          // deltaAB is stddevs involved in placing B to right of A
          deltaAB = (preExpansionA + addExpansionA) / stddevGapA + compressionB / stddevGapB;
          deltaBA = (preExpansionB + addExpansionB) / stddevGapB + compressionA / stddevGapA;
          if(deltaBA > deltaAB)
            {
              // place contigB to the right of contigA
              contigsB[ib].lft_end =
                MIN(contigsB[ib+1].lft_end - gapsB[ib].gap_length,
                    contigsA[ia+1].lft_end - MIN_GAP_LENGTH) - contigsB[ib].length;
              ib--;
            }
          else
            {
              contigsA[ia].lft_end =
                MIN(contigsA[ia+1].lft_end - gapsA[ia].gap_length,
                    contigsB[ib+1].lft_end - MIN_GAP_LENGTH) - contigsA[ia].length;
              ia--;
            }
        }
    }
}

int AdjustScaffoldContigPositions(CIScaffoldT * scaffold,
                                  Scaffold_Tig * contigs,
                                  int numContigs,
                                  int isAB)
{
  float offset;
  float scaffoldLength;
  CIScaffoldTIterator contigIterator;

  offset = contigs[0].lft_end;
  scaffoldLength = contigs[numContigs - 1].lft_end + contigs[numContigs - 1].length - offset;
  scaffold->bpLength.mean = scaffoldLength;
  
  InitCIScaffoldTIterator(ScaffoldGraph, scaffold,
                          TRUE, FALSE, &contigIterator);

  if(isAB)
    {
      // contigs in contigs array and scaffold are in same order
      int i;
      ChunkInstanceT * contig;
    
      i = 0;
      while((contig = NextCIScaffoldTIterator(&contigIterator)) != NULL)
        {
#ifndef RAT_RUN
          assert(contigs[i].length < contig->bpLength.mean + 5 &&
                 contigs[i].length > contig->bpLength.mean - 5);
#endif
          if(contig->offsetAEnd.mean < contig->offsetBEnd.mean)
            {
              // AB scaffold and AB contig
              contig->offsetAEnd.mean = contigs[i].lft_end - offset;
              contig->offsetBEnd.mean =
                contigs[i].lft_end + contigs[i].length - offset;
            }
          else
            {
              contig->offsetBEnd.mean = contigs[i].lft_end - offset;
              contig->offsetAEnd.mean =
                contigs[i].lft_end + contigs[i].length - offset;
            }
          i++;
        }
    }
  else
    {
      // contigs in contigs array and scaffold are in reverse order
      int i;
      ChunkInstanceT * contig;
    
      i = numContigs - 1;
      while((contig = NextCIScaffoldTIterator(&contigIterator)) != NULL)
        {
#ifndef RAT_RUN
          assert(contigs[i].length < contig->bpLength.mean + 5 &&
                 contigs[i].length > contig->bpLength.mean - 5);
#endif
          if(contig->offsetAEnd.mean < contig->offsetBEnd.mean)
            {
              // BA scaffold and AB contig
              contig->offsetAEnd.mean = scaffoldLength -
                (contigs[i].lft_end + contigs[i].length - offset);
              contig->offsetBEnd.mean =
                scaffoldLength - (contigs[i].lft_end - offset);
            }
          else
            {
              contig->offsetBEnd.mean = scaffoldLength -
                (contigs[i].lft_end + contigs[i].length - offset);
              contig->offsetAEnd.mean =
                scaffoldLength - (contigs[i].lft_end - offset);
            }
          i--;
        }
    }
  return 0;
}


void AdjustForPureInterleaving(Scaffold_Tig * contigsA,
                               Scaffold_Gap * gapsA,
                               int numContigsA,
                               Scaffold_Tig * contigsB,
                               Scaffold_Gap * gapsB,
                               int numContigsB,
                               SEdgeT * sEdge)
{
  // get at least one contig in each scaffold set then call other function
  int ia;
  int ib;
  double edgeStddev = sqrt((double) sEdge->distance.variance);
  double deltaAB;
  double deltaBA;
  double adjust =  -sEdge->distance.mean - (contigsA[numContigsA-1].lft_end +
                                            contigsA[numContigsA-1].length);

#ifdef DEBUG1
  fprintf(stderr, "Adjusting positions of Scaffold A contigs by %f, edge length %f len(A) %d len(B) %d\n",
	  adjust,sEdge->distance.mean,
	  (contigsA[numContigsA-1].lft_end +contigsA[numContigsA-1].length),
	  (contigsB[numContigsB-1].lft_end +contigsB[numContigsB-1].length));
#endif

  // adjust contigs in scaffold A
  for(ia = 0; ia < numContigsA; ia++)
    contigsA[ia].lft_end += adjust;
  
  for(ia = 0, ib = 0; ia == 0 || ib == 0; )
    {
      CDS_COORD_t aLeft, bLeft;
      // contigs ia & ib occupy same space. place one or the other
      if(ia == numContigsA)
        {
          // no more A contigs, ib == 0
          contigsB[0].lft_end =
            MAX(contigsB[0].lft_end,
                contigsA[ia-1].lft_end + contigsA[ia-1].length + MIN_GAP_LENGTH);
          ib++;
          break;
        }
      if(ib == numContigsB)
        {
          // no more B contigs, ia == 0
          contigsA[0].lft_end =
            MAX(contigsA[0].lft_end,
                contigsB[ib-1].lft_end + contigsB[ib-1].length + MIN_GAP_LENGTH);
          ia++;
          break;
        }

      aLeft = contigsA[ia].lft_end;
      /*
        if(ia > 0)
        aLeft = contigsA[ia-1].lft_end + contigsA[ia-1].length + gapsA[ia-1].gap_length;
        else */
      if(ib > 0)
        aLeft = MAX(aLeft, contigsB[ib-1].lft_end + contigsB[ib-1].length + MIN_GAP_LENGTH);

      bLeft = contigsB[ib].lft_end;
      /*
        if(ib > 0)
        bLeft = contigsB[ib-1].lft_end + contigsB[ib-1].length + gapsB[ib-1].gap_length;
        else
      */
      if(ia > 0)
        bLeft = MAX(bLeft, contigsA[ia-1].lft_end + contigsA[ia-1].length + MIN_GAP_LENGTH);
                
      // edge compression for ia - ib
      if(ia == 0)
        deltaAB = MAX(0, (aLeft - bLeft - MIN_GAP_LENGTH - 
                          contigsA[ia].length) / edgeStddev);
      else
        deltaAB = MAX(0, (aLeft -
                          MAX(bLeft - MIN_GAP_LENGTH -
                              contigsA[ia].length,
                              contigsA[ia-1].lft_end + contigsA[ia-1].length +
                              MIN_GAP_LENGTH)) / edgeStddev);

      // add possible gaps[ia] compression
      if(ia < numContigsA - 1)
        deltaAB += MAX(0, (aLeft - contigsA[ia+1].lft_end -
                           gapsA[ia].gap_length -
                           contigsA[ia].length) / gapsA[ia].gap_var);
       
      /*
      // add possible gap stretch for ia - ib
      if(ia < numContigsA - 1 &&
      contigsB[ib].lft_end < contigsA[ia+1].lft_end + contigsA[ia+1].length)
      deltaAB += MAX(0, (contigsB[ib].length + 2 * MIN_GAP_LENGTH -
      gapsA[ia].gap_length) / gapsA[ia].gap_var);
      */

      // edge stretch for ib - ia
      if(ib == 0)
        deltaBA = MAX(0, (bLeft - aLeft - MIN_GAP_LENGTH - 
                          contigsB[ib].length) / edgeStddev);
      else
        deltaBA = MAX(0, (bLeft -
                          MAX(aLeft - MIN_GAP_LENGTH -
                              contigsB[ib].length,
                              contigsB[ib-1].lft_end + contigsB[ib-1].length +
                              MIN_GAP_LENGTH)) / edgeStddev);

      // add possible gapsB[ib] compression
      if(ib < numContigsB - 1)
        deltaBA += MAX(0, (bLeft - contigsB[ib+1].lft_end -
                           gapsB[ib].gap_length -
                           contigsB[ib].length) / gapsB[ib].gap_var);

      /*
      // add possible gap stretch for ib - ia
      if(ib < numContigsB - 1 &&
      contigsA[ia].lft_end < contigsB[ib+1].lft_end + contigsB[ib+1].length)
      deltaBA += MAX(0, (contigsA[ia].length + 2 * MIN_GAP_LENGTH -
      gapsB[ib].gap_length) / gapsB[ib].gap_var);
      */

      if(deltaAB < deltaBA)
        {
          if(ib > 0) // ia == 0
            contigsA[ia].lft_end = aLeft;
          /*
            MAX(contigsA[ia].lft_end,
            contigsB[ib-1].lft_end + contigsB[ib-1].length + MIN_GAP_LENGTH);
          */
          // else ia is already okay relative to contigsA[ia-1]
          ia++;
        }
      else
        {
          if(ia > 0) // ib == 0
            contigsB[ib].lft_end = bLeft;
          /*
            MAX(contigsB[ib].lft_end,
            contigsA[ia-1].lft_end + contigsA[ia-1].length + MIN_GAP_LENGTH);
          */
          // else ib is already okay relative to contigsB[ib-1]
          ib++;
        }
    }
  
  AdjustNonOverlappingContigsLeftToRight(contigsA,
                                         gapsA,
                                         ia,
                                         numContigsA - 1,
                                         contigsB,
                                         gapsB,
                                         ib,
                                         numContigsB - 1);
}


void PrintScaffold_Tigs(FILE * fp, Scaffold_Tig * contigs, int numContigs)
{
  int i;
  for(i = 0; i < numContigs; i++)
    {
      fprintf(fp, "%d: length: " F_COORD ", left end: " F_COORD "\n",
              i, contigs[i].length, contigs[i].lft_end);
    }
}


void PrintSegments(FILE * fp, Segment * segments)
{
  if(segments == NULL)
    {
      fprintf(fp, "None.\n");
    }
  else
    {
      Segment * thisSegment;
      for(thisSegment = segments;
          thisSegment != NULL;
          thisSegment = thisSegment->next)
        {
          fprintf(fp, "contigA: %d, contigB: %d, begpos: " F_COORD ", endpos: " F_COORD ", length: " F_COORD "\n",
                  thisSegment->a_contig, thisSegment->b_contig,
                  thisSegment->overlap->begpos, thisSegment->overlap->endpos,
                  thisSegment->overlap->length);
        }
    }
}


/*
  push contigs left or right so contig overlaps in segments will be
  found back later when scaffolds are actually merged
  also push contigs left or right so non-overlapping contigs don't overlap
  scaffolds in sai are populated so overlap betwen them is AB_AB
  ScaffoldStuff orient field indicates whether scaffold needs to be flipped
  best is best aHang
    
  if best ahang indicates different edge distance, adjust sEdge
  
  work through segment list (overlaps) & anchor overlapping contigs
  wrt each other & adjust positions of contigs between overlapping sets


  Returns non-NULL pointer to static edge with new distance.mean if
  all went well, otherwise returns NULL 
*/
SEdgeT * MakeScaffoldAlignmentAdjustments(CIScaffoldT * scaffoldA,
                                          CIScaffoldT * scaffoldB,
                                          SEdgeT * sEdge,
                                          ScaffoldAlignmentInterface * sai)
{
  Scaffold_Gap * gapsA = GetVA_Scaffold_Gap(sai->scaffoldA->pools->gapPool, 0);
  Scaffold_Gap * gapsB = GetVA_Scaffold_Gap(sai->scaffoldB->pools->gapPool, 0);
  Scaffold_Tig * contigsA = GetVA_Scaffold_Tig(sai->scaffoldA->pools->tigPool, 0);
  int numContigsA = GetNumVA_Scaffold_Tig(sai->scaffoldA->pools->tigPool);
  int lastPlacedA = 0;
  Scaffold_Tig * contigsB = GetVA_Scaffold_Tig(sai->scaffoldB->pools->tigPool, 0);
  int numContigsB = GetNumVA_Scaffold_Tig(sai->scaffoldB->pools->tigPool);
  int lastPlacedB = 0;
  Segment * overlaps;
  int numOverlaps;
  int i;
  VA_TYPE(COSData) * cosData = NULL;
  CDS_CID_t idA;
  CDS_CID_t idB;
  ChunkOrientationType orient;
  static SEdgeT mySEdge;
  double newEdgeMean;

  idA = sEdge->idA;
  idB = sEdge->idB;
  orient = sEdge->orient;

  if(sEdge->idA != scaffoldA->id)
    {
      sEdge->idA = idB;
      sEdge->idB = idA;
      sEdge->orient = FlipEdgeOrient(sEdge->orient);
    }

  // organize overlaps so they're easier to use
  if(sai->segmentList != NULL)
    {
      int numReordered;
    
      // copy linked list of overlaps to sorted array

      // count the number of overlaps in the linked list
      for(numOverlaps = 0, overlaps = sai->segmentList;
          overlaps != NULL;
          overlaps = overlaps->next, numOverlaps++);
    
      overlaps = (Segment *) safe_malloc(numOverlaps * sizeof(Segment));
      for(numOverlaps = 0;
          sai->segmentList != NULL;
          sai->segmentList = sai->segmentList->next, numOverlaps++)
        overlaps[numOverlaps] = *(sai->segmentList);
    
      // sort them left to right
      if(numOverlaps > 1)
        qsort(overlaps, numOverlaps, sizeof(Segment),
              (int (*) (const void *, const void *)) segmentCompare);

      cosData = CreateVA_COSData(numOverlaps);
      assert(cosData != NULL);

      // identify sets of overlapping contigs
      if((numReordered = ExamineContigOverlapSets(sai,
                                                  overlaps, numOverlaps,
                                                  cosData)) != 0)
        {
          fprintf(GlobalData->stderrc, "*** WARNING ***\n"
                  "\tScaffolds " F_CID " and " F_CID " can't be merged with edge (%.2f,%.2f)\n"
                  "\tbecause %d contigs would have to be re-ordered\n"
                  "\t(consider modifying the code to handle this)\n",
                  scaffoldA->id, scaffoldB->id,
                  sEdge->distance.mean, sEdge->distance.variance,
                  numReordered);

          // change sEdge back to what it was
          sEdge->idA = idA;
          sEdge->idB = idB;
          sEdge->orient = orient;

          return NULL;
        }

      //  Dump debug on cosData
#ifdef DEBUG1
      for(i = 1; i < GetNumVA_COSData(cosData); i++) {
        COSData *cos = GetVA_COSData(cosData, i);

        fprintf(stderr, "COSdata[%2d] - index=%d firstOverlap=%d lastOverlap=%d\n",
                i,
                cos->index,
                cos->firstOverlap,
                cos->lastOverlap);
        fprintf(stderr, "              Ainterval min %d(%d) max %d(%d)\n",
                cos->a.minIndex, cos->a.minCoord,
                cos->b.minIndex, cos->b.minCoord);
        fprintf(stderr, "              Binterval min %d(%d) max %d(%d)\n",
                cos->a.minIndex, cos->a.minCoord,
                cos->b.minIndex, cos->b.minCoord);
      }
#endif

      // Set positions of contigs to the left of the first overlap set
      PlaceContigsLeftOfFirstOverlapSet(GetVA_COSData(cosData,0),
                                        contigsA, gapsA,
                                        contigsB, gapsB);

      // loop until we run out of overlap sets
      for(i = 1; i < GetNumVA_COSData(cosData); i++)
        {
          PlaceContigsBetweenOverlapSets(GetVA_COSData(cosData, i - 1),
                                         GetVA_COSData(cosData, i),
                                         contigsA, gapsA,
                                         contigsB, gapsB);
          PlaceContigsInOverlapSet(GetVA_COSData(cosData, i),
                                   contigsA, gapsA,
                                   contigsB, gapsB);
        }

      lastPlacedA = (GetVA_COSData(cosData,
                                   GetNumVA_COSData(cosData) - 1))->a.maxIndex;
      lastPlacedB = (GetVA_COSData(cosData,
                                   GetNumVA_COSData(cosData) - 1))->b.maxIndex;

      DeleteVA_COSData(cosData);
      safe_free(overlaps);
    
      // deal with contigs to the right of the last placed contigs
      AdjustNonOverlappingContigsLeftToRight(contigsA,
                                             gapsA,
                                             lastPlacedA + 1,
                                             numContigsA - 1,
                                             contigsB,
                                             gapsB,
                                             lastPlacedB + 1,
                                             numContigsB - 1);
    }
  else
    {
      // pure interleaving
      AdjustForPureInterleaving(contigsA, gapsA, numContigsA,
                                contigsB, gapsB, numContigsB,
                                sEdge);
    }

  // check if edge distance is acceptable
  newEdgeMean = contigsB[0].lft_end -
    (contigsA[numContigsA-1].lft_end + contigsA[numContigsA-1].length);


  //  What this does is test whether the proposed layout is consistent
  //  with the original edge.  The problem with the test as written is
  //  as follows: assume two scaffolds with several contigs each,
  //  interleaved in a fashion that leaves each individual gap ok but
  //  slightly stretches both scaffolds in the region of their overlap.
  //  In this case, the thickness of the overlap will increase by the
  //  amount that the A scaffold is stretched; when this stretched
  //  value is compared back to the original estimate of the edge
  //  depth, they may turn out to be incompatible (leading to a failure
  //  of the above test) without there really being anything wrong.
  //
  //  Now it may be that we simply want to be VERY cautious about any
  //  interleaving, especially pure interleaving with no contig
  //  overlaps to support it, but this particular piece of logic doesnt
  //  seem well-suited to distinguishing the good from the bad.  If we
  //  dont really want to do aggressive interleaving, lets turn it off
  //  rather than wasting a lot of time doing computations we are going
  //  to discard.  --Aaron
  //
#ifdef CHECK_INTERLEAVE_DISTANCE
  if(sEdge->distance.mean + INTERLEAVE_CUTOFF * sqrt(sEdge->distance.variance) < newEdgeMean ||
     sEdge->distance.mean - INTERLEAVE_CUTOFF * sqrt(sEdge->distance.variance) > newEdgeMean)
    {
      // change sEdge back to what it was
      sEdge->idA = idA;
      sEdge->idB = idB;
      sEdge->orient = orient;

      fprintf(GlobalData->stderrc,
              "WARNING - Interleaved scaffold adjustments stretch or compress edge too much!\n");
    
      fprintf(GlobalData->stderrc, "Original Edge:\n");
      PrintSEdgeT(GlobalData->stderrc, ScaffoldGraph, "", sEdge, sEdge->idA);
      fprintf(GlobalData->stderrc, "\nNew mean would be %.f\n", newEdgeMean);

      fprintf(GlobalData->stderrc, "Contig overlaps:\n");
      PrintSegments(GlobalData->stderrc, sai->segmentList);

      /*
        fprintf(GlobalData->stderrc, "\nOriginal scaffolds:\n");
        DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffoldA, FALSE);
        DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffoldB, FALSE);

        fprintf(GlobalData->stderrc, "\nAdjusted contig positions:");
        fprintf(GlobalData->stderrc, "Scaffold 'A'\n");
        PrintScaffold_Tigs(GlobalData->stderrc, contigsA, numContigsA);
        fprintf(GlobalData->stderrc, "\nScaffold 'B'\n");
        PrintScaffold_Tigs(GlobalData->stderrc, contigsB, numContigsB);
      */
      return NULL;
    }
#endif

    
  // now map the adjustments over to the cgw scaffold/contig data structures
  // adjust scaffoldA


#ifdef CONNECTEDNESS_CHECKS
  // not sure whether the MarkInternalEdgeStatus calls below are helpful or not?
  MarkInternalEdgeStatus(ScaffoldGraph,scaffoldA, 
			 PAIRWISECHI2THRESHOLD_CGW,
			 1000.0 * SLOPPY_EDGE_VARIANCE_THRESHHOLD,
			 TRUE, TRUE, 0, TRUE);
  MarkInternalEdgeStatus(ScaffoldGraph,scaffoldB, 
			 PAIRWISECHI2THRESHOLD_CGW,
			 1000.0 * SLOPPY_EDGE_VARIANCE_THRESHHOLD,
			 TRUE, TRUE, 0, TRUE);

  if(!IsScaffoldInternallyConnected(ScaffoldGraph,scaffoldA,ALL_TRUSTED_EDGES))
    fprintf(stderr,"WARNING: Interleaved merging: scaffoldA %d pre-adjustment is not internally connected\n",scaffoldA->id);

  if(!IsScaffoldInternallyConnected(ScaffoldGraph,scaffoldB,ALL_TRUSTED_EDGES))
    fprintf(stderr,"WARNING: Interleaved merging: scaffoldB %d pre-adjustment is not internally connected\n",scaffoldB->id);

  if(!IsScaffold2EdgeConnected(ScaffoldGraph,scaffoldA))
    fprintf(stderr,"WARNING: Interleaved merging: scaffoldA %d pre-adjustment is not 2-edge connected\n",scaffoldA->id);

  if(!IsScaffold2EdgeConnected(ScaffoldGraph,scaffoldB))
    fprintf(stderr,"WARNING: Interleaved merging: scaffoldB %d pre-adjustment is not 2-edge connected\n",scaffoldB->id);
#endif

  AdjustScaffoldContigPositions(scaffoldA, contigsA, numContigsA,
                                (sEdge->orient == AB_AB ||
                                 sEdge->orient == AB_BA));
  AdjustScaffoldContigPositions(scaffoldB, contigsB, numContigsB,
                                (sEdge->orient == AB_AB ||
                                 sEdge->orient == BA_AB));

#ifdef CONNECTEDNESS_CHECKS
  MarkInternalEdgeStatus(ScaffoldGraph,scaffoldA, 
			 PAIRWISECHI2THRESHOLD_CGW,
			 1000.0 * SLOPPY_EDGE_VARIANCE_THRESHHOLD,
			 TRUE, TRUE, 0, TRUE);
  MarkInternalEdgeStatus(ScaffoldGraph,scaffoldB, 
			 PAIRWISECHI2THRESHOLD_CGW,
			 1000.0 * SLOPPY_EDGE_VARIANCE_THRESHHOLD,
			 TRUE, TRUE, 0, TRUE);

  if(!IsScaffoldInternallyConnected(ScaffoldGraph,scaffoldA,ALL_TRUSTED_EDGES))
    fprintf(stderr,"WARNING: Interleaved merging: scaffoldA %d post-adjustment is not internally connected\n",scaffoldA->id);

  if(!IsScaffoldInternallyConnected(ScaffoldGraph,scaffoldB,ALL_TRUSTED_EDGES))
    fprintf(stderr,"WARNING: Interleaved merging: scaffoldB %d post-adjustment is not internally connected\n",scaffoldB->id);

  if(!IsScaffold2EdgeConnected(ScaffoldGraph,scaffoldA))
    fprintf(stderr,"WARNING: Interleaved merging: scaffoldA %d post-adjustment is not 2-edge connected\n",scaffoldA->id);

  if(!IsScaffold2EdgeConnected(ScaffoldGraph,scaffoldB))
    fprintf(stderr,"WARNING: Interleaved merging: scaffoldB %d post-adjustment is not 2-edge connected\n",scaffoldB->id);
#endif

  // change sEdge back to what it was
  sEdge->idA = idA;
  sEdge->idB = idB;
  sEdge->orient = orient;

  mySEdge = *sEdge;
  // adjust edge so everything fits together nicely
  mySEdge.distance.mean = newEdgeMean;

#ifdef DEBUG1
  fprintf(GlobalData->stderrc, "Post-adjustment: Scaffold A CGW data structure:\n");
  DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffoldA, FALSE);
  fprintf(GlobalData->stderrc, "Post-adjustment: Scaffold B CGW data structure:\n");
  DumpCIScaffold(GlobalData->stderrc, ScaffoldGraph, scaffoldB, FALSE);
  fprintf(GlobalData->stderrc, "\nsEdge cgw data structure:\n");
  PrintSEdgeT(GlobalData->stderrc, ScaffoldGraph, "sEdge", sEdge,
              scaffoldA->id);
#endif


  return &mySEdge;
}


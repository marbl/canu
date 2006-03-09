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

#include "eCR.h"

void saveDefaultLocalAlignerVariables(void);
void restoreDefaultLocalAlignerVariables(void);


// externable variables for controlling use of Local_Overlap_AS_forCNS

// [initialized value is 12 -- no more than this many segments in the chain]
extern int MaxGaps;

// [ init value is 200; this could be set to the amount you extend the clear 
// range of seq b, plus 10 for good measure]
extern int MaxBegGap;

// [ init value is 200; this could be set to the amount you extend the 
// clear range of seq a, plus 10 for good measure]
extern int MaxEndGap;

// [ initial value is 1000 (should have almost no effect)
// and defines the largest gap between segments in the chain]
// Also: allowed size of gap within the alignment
// -- forcing relatively good alignments, compared to those
// allowed in bubble-smoothing where indel polymorphisms are expected
extern int MaxInteriorGap;

// boolean to cause the size of an "end gap" to
// be evaluated with regard to the clear range extension
extern int asymmetricEnds;


static int DefaultMaxBegGap;
static int DefaultMaxEndGap;
static int DefaultMaxGaps;
static int DefaultMaxInteriorGap;
static int DefaultAsymmetricEnds;

void
saveDefaultLocalAlignerVariables(void) {
  DefaultMaxBegGap      = MaxBegGap;
  DefaultMaxEndGap      = MaxEndGap;
  DefaultMaxGaps        = MaxGaps;
  DefaultMaxInteriorGap = MaxInteriorGap;
  DefaultAsymmetricEnds = asymmetricEnds;
}

void
restoreDefaultLocalAlignerVariables(void) {
  MaxBegGap      = DefaultMaxBegGap;
  MaxEndGap      = DefaultMaxEndGap;
  MaxGaps        = DefaultMaxGaps;
  MaxInteriorGap = DefaultMaxInteriorGap;
  asymmetricEnds = DefaultAsymmetricEnds;
}


#if 0
//  Used to be before CreateAContigInScaffold
MaxGaps = 5;
MaxBegGap = 200;
MaxEndGap = 200;
MaxInteriorGap = 30;
asymmetricEnds = TRUE;
#endif



int
examineGap(ContigT *lcontig, int lFragIid, ContigT *rcontig, int rFragIid, 
           int gapNumber, int *ahang, int *olapLengthOut, int *bhang, int *currDiffs,
           int *lcontigBasesIntact, int *rcontigBasesIntact,
           int *closedGapDelta, int lBasesToNextFrag, int rBasesToNextFrag,
           int *leftFragFlapLength, int *rightFragFlapLength) {

  CIFragT *lFrag = NULL, *rFrag = NULL;
  static VA_TYPE(char) *ungappedSequence=NULL, *ungappedQuality=NULL;
  char lFragSeqBuffer[AS_BACTIG_MAX_LEN+1], lqltbuffer[AS_BACTIG_MAX_LEN+1];
  char rFragSeqBuffer[AS_BACTIG_MAX_LEN+1], rqltbuffer[AS_BACTIG_MAX_LEN+1];
  char lcompBuffer[AS_BACTIG_MAX_LEN+1], rcompBuffer[AS_BACTIG_MAX_LEN+1];
  uint lclr_bgn, lclr_end;
  uint rclr_bgn, rclr_end;
  char tmp_char, *lSequence, *rSequence;
  NodeOrient lContigOrientation, rContigOrientation;
  InfoByIID *info;
  int temp, len;
  int i;
  int lcontigBaseStart, lcontigBasesUsed, rcontigBasesUsed;
  int lFragContigOverlapLength, rFragContigOverlapLength;

  //  bpw new
  // set some variables to control Local_Overlap_AS_forCNS
  MaxGaps        = 5;
  MaxInteriorGap = 30;
  asymmetricEnds = TRUE;

  if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
    lContigOrientation = A_B;
  else
    lContigOrientation = B_A;
  
  if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
    rContigOrientation = A_B;
  else
    rContigOrientation = B_A;

  if (lFragIid != -1) {
    info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, lFragIid);
    assert(info->set);
    lFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
  }
  
  if (rFragIid != -1) {
    info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);
    assert(info->set);
    rFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
  }
  
  if (ungappedSequence== NULL) {
    ungappedSequence = CreateVA_char(0);
    ungappedQuality = CreateVA_char(0);
  } else {
    ResetVA_char(ungappedSequence);
    ResetVA_char(ungappedQuality);
  }

  if (lFragIid != -1) {
    getFragStore(ScaffoldGraph->fragStore, lFragIid, FRAG_S_ALL, fsread);
    getClearRegion_ReadStruct(fsread, &lclr_bgn, &lclr_end, READSTRUCT_CNS);
    getSequence_ReadStruct(fsread, lFragSeqBuffer, lqltbuffer, AS_BACTIG_MAX_LEN);
  }
  
  if (rFragIid != -1) {
    getFragStore(ScaffoldGraph->fragStore, rFragIid, FRAG_S_ALL, fsread);
    getClearRegion_ReadStruct(fsread, &rclr_bgn, &rclr_end, READSTRUCT_CNS);
    getSequence_ReadStruct(fsread, rFragSeqBuffer, rqltbuffer, AS_BACTIG_MAX_LEN);
  }
  
  // Get the consensus sequences for both chunks from the Store
  GetConsensus(ScaffoldGraph->ContigGraph, lcontig->id, lContigConsensus, lContigQuality);
  GetConsensus(ScaffoldGraph->ContigGraph, rcontig->id, rContigConsensus, rContigQuality);
  
  lSequence = Getchar(lContigConsensus, 0);
  rSequence = Getchar(rContigConsensus, 0);
  
  // ----------------------> lContigOrientation == A_B
  //                  -----> frag is 5p->3p into gap, aligned with contig
  
  // <---------------------- lContigOrientation == B_A
  //                  -----> frag is 5p->3p into gap, aligned opposite to contig
  
  // the frag is oriented opposite to the contig in this case
  if (lContigOrientation == B_A)
    SequenceComplement(lSequence, NULL);
  
  // print out info on the left contig and fragment

  if (lFragIid != -1) {
#if 0
    char tmp_char = lFragSeqBuffer[lclr_end];
    lFragSeqBuffer[lclr_end] = '\0';
    fprintf(stderr, " last 50 bases of lfragIid clr range: %s\n", &lFragSeqBuffer[lclr_end - 50]);
    lFragSeqBuffer[lclr_end] = tmp_char;
#endif  
    fprintf(stderr, "for frag %d, lclr_bgn: %d, lclr_end: %d, strlen(lFragSeqBuffer): " F_SIZE_T "\n", 
            lFragIid, lclr_bgn, lclr_end, strlen(lFragSeqBuffer));
    // fprintf(stderr, " lfrag: %s\n", lFragSeqBuffer);
  }
	
#if 0
  fprintf(stderr, "  last 50 bases of lContig consensus: %s\n\n", &lSequence[ strlen(lSequence) - 50]);
#endif
  
  // ----------------------> rContigOrientation == A_B
  // <-----                  frag is 5p->3p into gap, aligned opposite to contig
  
  // <---------------------- rContigOrientation == B_A
  // <-----                  frag is 5p->3p into gap, aligned with contig
  
  // now do right contig
  if (rContigOrientation == B_A)  // the frag is oriented opposite to the contig in this case
    {
      SequenceComplement(rSequence, NULL);  // flip contig sequence to its orientation in scaffold
    }
  
  if (rFragIid != -1) {
    // we want to flip the frag in either case
    SequenceComplement(rFragSeqBuffer, NULL);
    len = strlen(rFragSeqBuffer);
    temp = len - rclr_bgn;  // new rclr_end
    rclr_bgn = len - rclr_end;
    rclr_end = temp;
  }
  
  // print out info on the right contig and fragment
  if (rFragIid != -1) {
#if 0
    char tmp_char = rFragSeqBuffer[rclr_bgn + 50];
    rFragSeqBuffer[rclr_bgn + 50] = '\0';
    fprintf(stderr, "first 50 bases of rFragIid clr range: %s\n", &rFragSeqBuffer[rclr_bgn]);
    rFragSeqBuffer[rclr_bgn + 50] = tmp_char;
#endif
	  
    fprintf(stderr, "for frag %d, rclr_bgn: %d, rclr_end: %d, strlen(rFragSeqBuffer): " F_SIZE_T "\n", 
            rFragIid, rclr_bgn, rclr_end, strlen(rFragSeqBuffer));
    // fprintf(stderr, " rfrag: %s\n", rFragSeqBuffer);
  }
	
#if 0
  tmp_char = rSequence[50];
  rSequence[50] = '\0';
  fprintf(stderr, " first 50 bases of rContig consensus: %50s\n\n", rSequence);
  rSequence[50] = tmp_char;
#endif
  
  // we use frag sequence from where the clear range ends to the end of the frag
  // ----------------------> lContigOrientation == A_B
  //               ----->    frag is 5p->3p into gap, aligned with contig
  // <---------------------- lContigOrientation == B_A
  //               ----->    frag is 5p->3p into gap, aligned opposite to contig
  
  if (lFragIid != -1) {
    if (lContigOrientation == A_B)
      lFragContigOverlapLength = (int) (lcontig->bpLength.mean - lFrag->contigOffset3p.mean);
    else
      lFragContigOverlapLength = (int) (lFrag->contigOffset3p.mean);
  } else {
    lFragContigOverlapLength = 0;
  }

  // grab the last CONTIG_BASES bases of the lcontig consensus sequence
  lcontigBasesUsed = min(CONTIG_BASES - lFragContigOverlapLength, 
                         (int) lcontig->bpLength.mean - lFragContigOverlapLength);
  // lcontigBasesUsed = 100.0;  // temp hack, but it is sometimes better to do 100 than 1000.  why???
  
  // lcontigBaseStart is the base where we start using the consensus sequence in lcompBuffer
  // and thus also the number of bases from the contig that are intact
  lcontigBaseStart = strlen(lSequence) - lcontigBasesUsed - lFragContigOverlapLength;
  
  // grab the bases from the contig, ie, those not from the non-clear range of the frag
  for (i = 0; i < lcontigBasesUsed; i++) {
    lcompBuffer[ i ] = lSequence[ lcontigBaseStart + i];  // a bit ugly
  }
  lcompBuffer[ i ] = '\0';

  // now tack on the 3p clr range extension to the bases of the contig consensus sequence
  if (lFragIid != -1) {
    // if (lcontigBasesUsed < lFragContigOverlapLength)  // this means that the frag and contig overlap by
    // lFragContigOverlapLength = lcontigBasesUsed;    // more bases than we are using from the contig


    // basesToNextFrag is the number of bases back to the first frag that gets us to 2x
    // used to be, but Aaron thought it looked fishy
    //MaxEndGap = strlen(lFragSeqBuffer) - lclr_end - lFragContigOverlapLength + lBasesToNextFrag + 20;  // 20 is slop
    MaxEndGap = strlen(lFragSeqBuffer) - lclr_end + lBasesToNextFrag + 20;  // 20 is slop
    fprintf(stderr,"## MaxEndGap %d\n",MaxEndGap);
    for (i = lclr_end; i < strlen(lFragSeqBuffer); i++)
      lcompBuffer[ lcontigBasesUsed + i - lclr_end ] = lFragSeqBuffer[ i ];
    lcompBuffer[ lcontigBasesUsed + i - lclr_end ] = '\0';
  } else {
    MaxEndGap = 100;
  }
  
  // we use frag sequence from where the clear range ends to the end of the frag
  // ----------------------> rContigOrientation == A_B
  //    <-----               frag is 5p->3p into gap, aligned opposite to contig

  // <---------------------- rContigOrientation == B_A
  //      <-----             frag is 5p->3p into gap, aligned with contig

  if (rFragIid != -1) {
    if (rContigOrientation == A_B)
      rFragContigOverlapLength = (int) (rFrag->contigOffset3p.mean);
    else
      rFragContigOverlapLength = (int) (rcontig->bpLength.mean - rFrag->contigOffset3p.mean);
  } else {
    rFragContigOverlapLength = 0;
  }
 
  if (rFragIid != -1) {
    // basesToNextFrag is the number of bases back to the first frag that gets us to 2x
    // used to be, but Aaron thought it looked fishy
    //MaxBegGap = rclr_bgn - rFragContigOverlapLength + rBasesToNextFrag + 20;  // 20 is slop	
    MaxBegGap = rclr_bgn + rBasesToNextFrag + 20;  // 20 is slop	
    fprintf(stderr,"## MaxBegGap %d\n",MaxBegGap);
  } else {
    MaxBegGap = 100;
  }
  
  // now if we have a right frag, grab the "5p" clr range extension - remember the frag has been flipped
  if (rFragIid != -1) {
    for (i = 0; i < rclr_bgn; i++)
      rcompBuffer[ i ] = rFragSeqBuffer[ i ];
  } else {
    rclr_bgn = 0;  // need this for the next loop
  }
  
  // grab the first CONTIG_BASES bases of the rcontig consensus sequence
  // the rcontig consensus has been flipped if necessary
  rcontigBasesUsed = min(CONTIG_BASES - rFragContigOverlapLength, 
                         (int) rcontig->bpLength.mean - rFragContigOverlapLength);
  for (i = 0; i < rcontigBasesUsed; i++)
    rcompBuffer[ rclr_bgn + i ] = rSequence[ i + rFragContigOverlapLength ]; //Aaron ad rFragContigOverlapLength
  rcompBuffer[ rclr_bgn + i ] = '\0';

  fprintf(stderr, "> lcompBuffer gap %d (len: " F_SIZE_T "): \n%s\n", gapNumber, strlen(lcompBuffer), lcompBuffer);
  fprintf(stderr, "> rcompBuffer gap %d (len: " F_SIZE_T "): \n%s\n", gapNumber, strlen(rcompBuffer), rcompBuffer);
  
  // now lcompBuffer and rcompBuffer hold the sequence of the fragments in the correct strand
  // now prepare for call to Local_Overlap_AS_forCNS
  {
    int beg, end, opposite = FALSE;
    double erate, thresh, minlen;
    CompareOptions what;
    Overlap *overlap;
    int basesAdded;
    LengthT gapSize;
    char *rcompBufferTrimmed = NULL;
	
    // MaxGaps = 5;  // we don't want a lot of segments
	
    // stole from MultiAlignment_CNS.c
#define CNS_DP_ERATE .12 // was .20 during testing
#define CNS_DP_THRESH 1e-6
#define CNS_DP_MINLEN 30

    beg    = -strlen (rcompBuffer);
    end    = strlen (lcompBuffer);
    erate  = CNS_DP_ERATE;
    thresh = CNS_DP_THRESH;
    minlen = CNS_DP_MINLEN;
    what   = AS_FIND_LOCAL_ALIGN;        //  was: what = AS_FIND_LOCAL_ALIGN_NO_TRACE;
	
    fprintf(stderr, "MaxBegGap: %d\n", MaxBegGap);
    fprintf(stderr, "MaxEndGap: %d\n", MaxEndGap);
	
    assert(strlen(rcompBuffer)>0);

    overlap = Local_Overlap_AS_forCNS(lcompBuffer,
                                      rcompBuffer,
                                      -strlen(rcompBuffer),
                                      strlen(lcompBuffer),
                                      opposite,
                                      erate,
                                      thresh,
                                      minlen,
                                      what);

    if (0 && overlap != NULL) {
      fprintf(stderr, "initial ahang: %d, bhang:%d, length: %d, diffs: %d, diffs / length %%: %f\n",
              overlap->begpos, overlap->endpos, overlap->length, overlap->diffs, 
              100.0 * overlap->diffs / overlap->length);
      Print_Overlap(stderr, lcompBuffer, rcompBuffer, overlap);
    }

    // not interested in overlaps with negative ahangs or bhangs
    if (overlap != NULL && overlap->begpos > 0 && overlap->endpos > 0) {
      /*
        \
        \
        \  left flap, right frag
        \  
        -------------------
        -------------------
        \   right flap, left frag
        \
        \
      */
          
      // do flap trimming
      // left flap, right frag
      *rightFragFlapLength = 0;
      if (overlap->begpos == abs(overlap->trace[0]) - 1 && overlap->trace[0] < 0) {
        while (overlap->begpos == abs(overlap->trace[ *rightFragFlapLength ]) - 1 && 
               overlap->trace[ *rightFragFlapLength ] < 0)
          (*rightFragFlapLength)++;
      }
	  
      // right flap, left frag
      {
        // first find the number of elements in trace
        int numTraceElements = 0, icnt = 0;
        while (overlap->trace[ icnt++ ] != 0) numTraceElements++;
		
        *leftFragFlapLength = 0;
        if (overlap->trace[ numTraceElements - 1 ] > 0 && overlap->endpos > 0) {
          icnt = numTraceElements - 1;
          while ((icnt >= 0) &&
                 (overlap->trace[ icnt ] == strlen(rcompBuffer) - overlap->endpos + 1)) {
            (*leftFragFlapLength)++; icnt--;
          }
        }
      }
	
      lcompBuffer[ strlen(lcompBuffer) - *leftFragFlapLength ] = '\0';
      // MaxEndGap -= *leftFragFlapLength;
      rcompBufferTrimmed = &rcompBuffer[ *rightFragFlapLength ];
      // MaxBegGap = *rightFragFlapLength;

      // now do overlap again after trimming to make sure it is still
      // there, sometimes trimming makes them go away
#if 1
      overlap = Local_Overlap_AS_forCNS(lcompBuffer,
                                        rcompBufferTrimmed,
                                        -strlen(rcompBufferTrimmed),
                                        strlen(lcompBuffer),
                                        opposite,
                                        erate,
                                        thresh,
                                        minlen,
                                        what);
      if (!overlap)
        fprintf(stderr, "lost overlap to flap trimming!\n");
#endif
    }

    if (1 && overlap != NULL && overlap->begpos > 0 && overlap->endpos > 0) {
      fprintf(stderr, "post-flap trimming ahang: %d, bhang:%d, length: %d, diffs: %d, diffs / length %%: %f\n",
              overlap->begpos, overlap->endpos, overlap->length, overlap->diffs, 
              100.0 * overlap->diffs / overlap->length);
      Print_Overlap(stderr, lcompBuffer, rcompBufferTrimmed, overlap);
    }
	
    // not interested in overlaps with negative ahangs or bhangs
    if (overlap != NULL && overlap->begpos > 0 && overlap->endpos > 0) {
      int baseChangeLeftContig, baseChangeRightContig;

      fprintf(stderr, "found overlap between frags %d and %d, length = %d\n", 
              lFragIid, rFragIid, overlap->length);
      fprintf(stderr, "ahang + overlap->length: %d, strlen(lcompBuffer): " F_SIZE_T ", diff: " F_SIZE_T "\n",
              overlap->begpos + overlap->length, strlen(lcompBuffer),
              overlap->begpos + overlap->length - strlen(lcompBuffer));
      fprintf(stderr, "overlap->length + bhang: %d, strlen(rcompBufferTrimmed): " F_SIZE_T ", diff: " F_SIZE_T "\n",
              overlap->length + overlap->endpos, strlen(rcompBufferTrimmed), 
              overlap->length + overlap->endpos - strlen(rcompBufferTrimmed));
	  
      *lcontigBasesIntact = lcontigBaseStart;
      *ahang = max (0, overlap->begpos);
      *olapLengthOut = overlap->length;
      *bhang = max (0, overlap->endpos);
      *currDiffs = overlap->diffs;
      *rcontigBasesIntact = max((int) rcontig->bpLength.mean - CONTIG_BASES, 0);

      // calculate how many bases would be changed if gap was closed
      // take the whole lcompBuffer, subtract lcontigBasesUsed and lFragContigOverlapLength from ahang
      // and add in the length of the overlap
      // take the whole rcompBuffer, subtract rcontigBasesUsed and rFragContigOverlapLength from bhang

      baseChangeLeftContig  = overlap->begpos - lcontigBasesUsed - lFragContigOverlapLength;
      baseChangeRightContig = overlap->endpos - rcontigBasesUsed - rFragContigOverlapLength;

      fprintf(stderr, "lcontigBasesIntact: %d\n", *lcontigBasesIntact);
      fprintf(stderr, "overlap->begpos: %d\n", overlap->begpos);
      fprintf(stderr, "overlap->length: %d\n", overlap->length);
      fprintf(stderr, "overlap->endpos: %d\n", overlap->endpos);
      fprintf(stderr, "rcontigBasesIntact: %d\n", *rcontigBasesIntact);

      fprintf(stderr, "base change  left contig: %d\n", baseChangeLeftContig);
      fprintf(stderr, "base change right contig: %d\n", baseChangeRightContig);

      basesAdded = baseChangeLeftContig + overlap->length + baseChangeRightContig;

      gapSize = FindGapLength(lcontig, rcontig, FALSE);

      totalContigsBaseChange += basesAdded;

      *closedGapDelta = basesAdded - (int) gapSize.mean;
	  
      fprintf(stderr, "would fill gap %d of size %d with %d bases, net change: %d\n",
              gapNumber, (int) gapSize.mean, basesAdded, basesAdded - (int) gapSize.mean);
      fprintf(stderr, "lcontig->bpLength.mean: %f, baseChangeLeftContig: %d\n",
              lcontig->bpLength.mean, baseChangeLeftContig);
      fprintf(stderr, "rcontig->bpLength.mean: %f, baseChangeRightContig: %d\n",
              rcontig->bpLength.mean, baseChangeRightContig);

      fprintf(stderr, "new contig size: %d\n", 
              (int) lcontig->bpLength.mean + baseChangeLeftContig + 
              (int) rcontig->bpLength.mean + baseChangeRightContig + overlap->length);
      fprintf(stderr, "totalContigsBaseChange: %d\n", totalContigsBaseChange);
	  
      // *closedGapSizeDiff = basesAdded - (int) gapSize.mean;

      restoreDefaultLocalAlignerVariables();
      return 1;
    } else {
      fprintf(stderr, "no overlap found between frags %d and %d\n", lFragIid, rFragIid);
      restoreDefaultLocalAlignerVariables();
      return 0;
    }
  }
}

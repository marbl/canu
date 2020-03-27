
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

/*
 * https://github.com/PacificBiosciences/FALCON/blob/master/src/c/falcon.c
 *
 * Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted (subject to the limitations in the
 * disclaimer below) provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above
 *  copyright notice, this list of conditions and the following
 *  disclaimer in the documentation and/or other materials provided
 *  with the distribution.
 *
 *  * Neither the name of Pacific Biosciences nor the names of its
 *  contributors may be used to endorse or promote products derived
 *  from this software without specific prior written permission.
 *
 * NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
 * GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
 * BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#include "falconConsensus.H"
#include "edlib.H"

#undef  DEBUG_ALIGN
#undef  DEBUG_ALIGN_VERBOSE

static
alignTagList *
getAlignTags(char       *Qalign,   int32 Qbgn,  int32 Qlen, int32 UNUSED(Qid),    //  read
             char       *Talign,   int32 Tbgn,  int32 Tlen,                       //  template
             int32       alignLen) {
  int32   i        = Qbgn - 1;   //  Position in query, not really used.
  int32   j        = Tbgn - 1;   //  Position in template
  int32   p_j      = -1;

  uint32  jj       = 0;          //  Number of non-gap bases in Q aligned to a gap in T
  uint32  p_jj     = 0;

  char    p_q_base = '.';

  alignTagList  *tags = new alignTagList(alignLen);

  for (int32 k=0; k < alignLen; k++) {
    if (Qalign[k] != '-') {
      i++;
      jj++;
    }

    if (Talign[k] != '-') {
      j++;
      jj = 0;
    }

    assert(i >= 0);
    assert(i <  Qlen);
    assert(j >= 0);
    assert(j <  Tlen);

    if ((jj   >= uint16MAX) ||
        (p_jj >= uint16MAX))
      continue;

    tags->setTag(j, p_j, jj, p_jj, Qalign[k], p_q_base);

#ifdef DEBUG_ALIGN_VERBOSE
    fprintf(stderr, "set tag j %5d p_j %5d jj %5d p_jj %5d base %c p_q_base %c\n",
            j, p_j, jj, p_jj, Qalign[k], p_q_base);
#endif

    p_j       = j;
    p_jj      = jj;
    p_q_base  = Qalign[k];
  }

  return(tags);
}



alignTagList **
alignReadsToTemplate(falconInput    *evidence,
                     uint32          evidenceLen,
                     double          minOlapIdentity,
                     uint32          minOlapLength,
                     bool            restrictToOverlap) {

  double         maxDifference = 1.0 - minOlapIdentity;
  alignTagList **tagList = new alignTagList * [evidenceLen];

  //  I don't remember where this was causing problems, but reads longer than the template were.  So truncate them.

  for (uint32 j=0; j<evidenceLen; j++)
    if (evidence[j].readLength > evidence[0].readLength) {
      evidence[j].readLength = evidence[0].readLength;
      evidence[j].read[evidence[j].readLength]  = 0;
    }

  //  Set everything to an empty list.  Makes aborting the algnment loop much easier.

  for (uint32 j=0; j<evidenceLen; j++)
    tagList[j] = NULL;


#pragma omp parallel for schedule(dynamic)
  for (uint32 j=0; j<evidenceLen; j++) {
    if (evidence[j].readLength < minOlapLength)
      continue;

    int32 tolerance =  (int32)ceil(min(evidence[j].readLength, evidence[0].readLength) * maxDifference * 1.1);

    int32  alignBgn = (restrictToOverlap == true) ? evidence[j].placedBgn : 0;
    int32  alignEnd = (restrictToOverlap == true) ? evidence[j].placedEnd : evidence[0].readLength;

    assert(alignEnd > alignBgn);

    //  Extend the region we align to by ... some amount.
    //  For simplicity, we'll use 10% of the read length.

    int32  expansion = 0.1 * evidence[j].readLength;

  again:
    alignBgn -= expansion;
    alignEnd += expansion;

    if (alignBgn < 0)                         alignBgn = 0;
    if (alignEnd > evidence[0].readLength)    alignEnd = evidence[0].readLength;

#ifdef DEBUG_ALIGN
    fprintf(stderr, "ALIGN to %d-%d length %d\n",
            alignBgn, alignEnd, evidence[0].readLength);
#endif

    EdlibAlignResult align = edlibAlign(evidence[j].read,            evidence[j].readLength,
                                        evidence[0].read + alignBgn, alignEnd - alignBgn,
                                        edlibNewAlignConfig(tolerance, EDLIB_MODE_HW, EDLIB_TASK_PATH));

#ifdef DEBUG_ALIGN
    for (int32 l=0; l<align.numLocations; l++)
      fprintf(stderr, "read%u #%u location %d to template %d-%d length %d diff %f\n",
              evidence[j].ident,
              j,
              l,
              align.startLocations[l],
              align.endLocations[l],
              align.endLocations[l] - align.startLocations[l],
              (float)align.editDistance / (align.endLocations[l] - align.startLocations[l]));
#endif

    if (align.numLocations == 0) {
      edlibFreeAlignResult(align);
#ifdef DEBUG_ALIGN
      fprintf(stderr, "read %7u failed to map\n", j);
#endif
      continue;
    }

    int32  alignLen  = align.endLocations[0] - align.startLocations[0];
    double alignDiff = align.editDistance / (double)alignLen;

    if (alignLen < minOlapLength) {
      edlibFreeAlignResult(align);
#ifdef DEBUG_ALIGN
      fprintf(stderr, "read %7u failed to map - short\n", j);
#endif
      continue;
    }

    if (alignDiff >= maxDifference) {
      edlibFreeAlignResult(align);
#ifdef DEBUG_ALIGN
      fprintf(stderr, "read %7u failed to map - different\n", j);
#endif
      continue;
    }

    int32  rBgn = 0;
    int32  rEnd = evidence[j].readLength;

    int32  tBgn = alignBgn + align.startLocations[0];
    int32  tEnd = alignBgn + align.endLocations[0] + 1;    //  Edlib returns position of last base aligned

    if ((alignBgn > 0) &&
        (tBgn <= alignBgn)) {
      edlibFreeAlignResult(align);
#ifdef DEBUG_ALIGN
      fprintf(stderr, "bumped into start align %d-%d mapped %d-%d\n", alignBgn, alignEnd, tBgn, tEnd);
#endif
      goto again;
    }

    if ((alignEnd < evidence[0].readLength) &&
        (tEnd >= alignEnd)) {
      edlibFreeAlignResult(align);
#ifdef DEBUG_ALIGN
      fprintf(stderr, "bumped into end align %d-%d mapped %d-%d\n", alignBgn, alignEnd, tBgn, tEnd);
#endif
      goto again;
    }

    char *tAln = new char [align.alignmentLength + 1];
    char *rAln = new char [align.alignmentLength + 1];

    edlibAlignmentToStrings(align.alignment,
                            align.alignmentLength,
                            tBgn, tEnd,
                            rBgn, rEnd,
                            evidence[0].read, evidence[j].read,
                            tAln, rAln);

    //  Strip leading/trailing gaps on template sequence.

    uint32 fBase = 0;                        //  First non-gap in the alignment
    uint32 lBase = align.alignmentLength;    //  Last base in the alignment (actually, first gap in the gaps at the end, but that was too long for a variable name)

    while ((fBase < align.alignmentLength) && (tAln[fBase] == '-'))
      fBase++;

    while ((lBase > fBase) && (tAln[lBase-1] == '-'))
      lBase--;

    rBgn += fBase;
    rEnd -= align.alignmentLength - lBase;

    assert(rBgn >= 0);      assert(rEnd <= evidence[j].readLength);
    assert(tBgn >= 0);      assert(tEnd <= evidence[0].readLength);

    rAln[lBase] = 0;   //  Truncate the alignments before the gaps.
    tAln[lBase] = 0;

#ifdef DEBUG_ALIGN
    fprintf(stderr, "mapped %5u %5u-%5u to template %6u-%6u trimmed by %6u-%6u %s %s\n",
            evidence[j].ident,
            rBgn - fBase, rEnd + align.alignmentLength - lBase,
            tBgn, tEnd,
            fBase, align.alignmentLength - lBase,
            rAln + lBase - 10,
            tAln + lBase - 10);
#endif

    tagList[j] = getAlignTags(rAln + fBase, rBgn, evidence[j].readLength, j,
                              tAln + fBase, tBgn, evidence[0].readLength,
                              lBase - fBase);

    delete [] tAln;
    delete [] rAln;

    edlibFreeAlignResult(align);
  }

  return(tagList);
}

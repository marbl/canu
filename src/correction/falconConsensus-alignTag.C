
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

#include <stdarg.h>

#undef  DEBUG_ALIGN_VERBOSE


FILE **thrlog = nullptr;

void
openAlignLogFiles(char const *outputPrefix) {
  char  N[FILENAME_MAX+1] = {0};

  if (thrlog != nullptr)
    return;

  thrlog = new FILE * [getNumThreads()];

  for (uint32 ii=0; ii<getNumThreads(); ii++) {
    snprintf(N, FILENAME_MAX, "%s.align_%02u.log", outputPrefix, ii);

    thrlog[ii] = AS_UTL_openOutputFile(N);
  }
}


void
closeAlignLogFiles(char const *prefix) {

  if (thrlog == nullptr)
    return;

  for (uint32 ii=0; ii<getNumThreads(); ii++)
    AS_UTL_closeFile(thrlog[ii]);

  delete [] thrlog;
}


void
printAllLog(char const *format, ...) {
  va_list ap;

  if (thrlog == nullptr)
    return;

  va_start(ap, format);

  for (uint32 ii=0; ii<getNumThreads(); ii++)
    vfprintf(thrlog[ii], format, ap);

  va_end(ap);
}


void
printLog(falconInput *evidence, uint32 j, char const *format, ...) {
  va_list ap;

  if (thrlog == nullptr)
    return;

  uint32  tid = getThreadNum();

  va_start(ap, format);
  fprintf(thrlog[tid], "read%08u-evidence%08u len %-6u ", evidence[0].ident, evidence[j].ident, evidence[j].readLength);
  vfprintf(thrlog[tid], format, ap);
  va_end(ap);
}


void
printLog(char const *format, ...) {
  va_list ap;

  if (thrlog == nullptr)
    return;

  uint32  tid = getThreadNum();

  va_start(ap, format);
  vfprintf(thrlog[tid], format, ap);
  va_end(ap);
}



uint32
stripLeadingGaps(char *tAln,
                 char *rAln, int32 alignLen) {
  int32  fBase    = 0;
  bool   skipping = true;

  //  Skip initial gaps.
  while ((fBase < alignLen) && (tAln[fBase] == '-'))
    fBase++;

  //  Count the number of ACGT in a row, then the number of gaps we find.
  //  If there are more ACGT than gaps, stop skipping.
  while (skipping == true) {
    int32  nBase    = 0;
    int32  nGap     = 0;

    while ((fBase+nBase < alignLen) && (tAln[fBase+nBase] != '-'))
      nBase++;
    while ((fBase+nBase+nGap < alignLen) && (tAln[fBase+nBase+nGap] == '-'))
      nGap++;

    if ((nBase < 16) || (nBase < nGap))
      fBase += nBase + nGap;
    else
      skipping = false;

    if (fBase >= alignLen)
      skipping = false;
  }

  return(fBase);
}



uint32
stripTrailingGaps(char *tAln,
                  char *rAln, int32 alignLen) {
  int32  lBase    = alignLen - 1;
  bool   skipping = true;

  //  Skip initial gaps.
  while ((lBase >= 0) && (tAln[lBase] == '-'))
    lBase--;

  //  Count the number of ACGT in a row, then the number of gaps we find.
  //  If there are more ACGT than gaps, stop skipping.
  while (skipping == true) {
    int32  nBase    = 0;
    int32  nGap     = 0;

    while ((lBase-nBase >= 0) && (tAln[lBase-nBase] != '-'))
      nBase++;
    while ((lBase-nBase-nGap >= 0) && (tAln[lBase-nBase-nGap] == '-'))
      nGap++;

    if ((nBase < 16) || (nBase < nGap))
      lBase -= nBase + nGap;
    else
      skipping = false;

    if (lBase < 0)
      skipping = false;
  }

  return(alignLen - 1 - lBase);
}



static
alignTagList *
getAlignTags(char       *Qalign,   int32 Qbgn,  int32 Qlen,    //  read
             char       *Talign,   int32 Tbgn,  int32 Tlen,    //  template
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

    if ((jj   >= uint16max) ||
        (p_jj >= uint16max))
      continue;

    tags->setTag(j, p_j, jj, p_jj, Qalign[k], p_q_base);

#ifdef DEBUG_ALIGN_VERBOSE
    printLog("set tag j %5d p_j %5d jj %5d p_jj %5d base %c p_q_base %c\n",
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

  printAllLog("STARTING read %u length %u with %u evidence reads.\n",
              evidence[0].ident, evidence[0].readLength,
              evidenceLen);

  //  I don't remember where this was causing problems, but reads longer than the template were.  So truncate them.

  for (uint32 j=0; j<evidenceLen; j++)
    if (evidence[j].readLength > evidence[0].readLength) {
      evidence[j].readLength = evidence[0].readLength;
      evidence[j].read[evidence[j].readLength]  = 0;
    }

  //  Set everything to an empty list.  Makes aborting the algnment loop much easier.

  for (uint32 j=0; j<evidenceLen; j++)
    tagList[j] = NULL;

  //  The first evidence read is ourself, and we know how to add that alignment immediately.

  tagList[0] = getAlignTags(evidence[0].read, 0, evidence[0].readLength,
                            evidence[0].read, 0, evidence[0].readLength,
                            evidence[0].readLength);

  //  For the rest of the evidence: align, parse, add.

#pragma omp parallel for schedule(dynamic)
  for (uint32 j=1; j<evidenceLen; j++) {
    if (evidence[j].readLength < minOlapLength)
      continue;

    int32 tolerance =  (int32)ceil(std::min(evidence[j].readLength, evidence[0].readLength) * maxDifference * 1.1);

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

    EdlibAlignResult align = edlibAlign(evidence[j].read,            evidence[j].readLength,
                                        evidence[0].read + alignBgn, alignEnd - alignBgn,
                                        edlibNewAlignConfig(tolerance, EDLIB_MODE_HW, EDLIB_TASK_PATH));

    if (align.numLocations == 0) {
      printLog(evidence, j, "ALIGN to template %7d-%-7d   tolerance %5d -- FAILED TO ALIGN\n",
               alignBgn, alignEnd, tolerance);
      edlibFreeAlignResult(align);
      continue;
    }

    int32  tBgn = alignBgn + align.startLocations[0];
    int32  tEnd = alignBgn + align.endLocations[0] + 1;    //  Edlib returns position of last base aligned
    int32  tLen = tEnd - tBgn;
    double tDif = 100.0 * align.editDistance / tLen;

    int32  rBgn = 0;
    int32  rEnd = evidence[j].readLength;

    for (int32 ll=0; ll<align.numLocations; ll++)   //  Report ALL aligns, not just the first
      printLog(evidence, j, "    #%-2u  template %7d-%-7d  evidence %7d-%-7d  %6.3f%%\n",
               ll,
               alignBgn + align.startLocations[ll], alignBgn + align.endLocations[ll] + 1,
               rBgn, rEnd,
               100.0 * align.editDistance / (align.endLocations[ll] - align.startLocations[ll]));

    if (tLen < minOlapLength) {
      printLog(evidence, j, "ALIGNED  template %7d-%-7d   tolerance %5d -- TOO SHORT %u < %u\n",
               tBgn, tEnd, tolerance, tLen, minOlapLength);
      edlibFreeAlignResult(align);
      continue;
    }

    if (tDif >= 100.0 * maxDifference) {
      printLog(evidence, j, "ALIGNED  template %7d-%-7d   tolerance %5d -- TOO DIFFERENT %.3f%% >= %.3f%%\n",
               alignBgn, alignEnd, tolerance, tDif, 100.0 * maxDifference);
      edlibFreeAlignResult(align);
      continue;
    }

    if ((alignBgn > 0) && (tBgn <= alignBgn)) {
      printLog(evidence, j, "ALIGNED  template %7d-%-7d   tolerance %5d -- HIT BEGIN\n",
               tBgn, tEnd, tolerance);
      edlibFreeAlignResult(align);
      goto again;
    }

    if ((alignEnd < evidence[0].readLength) && (tEnd >= alignEnd)) {
      printLog(evidence, j, "ALIGNED  template %7d-%-7d   tolerance %5d -- HIT END\n",
               tBgn, tEnd, tolerance);
      edlibFreeAlignResult(align);
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

    //  Strip leading/trailing gaps on template sequence.  These are the
    //  number of alignment letters to ignore on either end.

    int32 bSkip = stripLeadingGaps (tAln, rAln, align.alignmentLength);   //  Alignment letters to skip
    int32 eSkip = stripTrailingGaps(tAln, rAln, align.alignmentLength);   //  because they're junk gaps.

    //  We need to know the position in the template that the (trimmed)
    //  alignment starts at, so we count how many non-gap letters are skipped
    //  in either sequence at either end.

    int32  tbSkip = 0, teSkip = 0;   //  Sequence bases to skip in (t)emplate or evidence (r)ead,
    int32  rbSkip = 0, reSkip = 0;   //  at the (b)egin or (e)nd of the sequence.

    for (int32 ii=0; ii<bSkip; ii++) {
      if (tAln[ii] != '-')   tbSkip++;
      if (rAln[ii] != '-')   rbSkip++;
    }

    for (int32 ii=0; ii<eSkip; ii++) {
      if (tAln[align.alignmentLength - ii - 1] != '-')   teSkip++;
      if (rAln[align.alignmentLength - ii - 1] != '-')   reSkip++;
    }

#if 1
    printLog(evidence, j, "ALIGNED  template %7d-%-7d  evidence %7d-%-7d  %6.3f%%  skip %d %d endgaps %d %d\n",
             tBgn + tbSkip, tEnd - teSkip,
             rBgn + rbSkip, rEnd - reSkip,
             bSkip, eSkip,
             rbSkip, reSkip);
#endif

#if 0
    printLog(evidence, j, "    full template %100.100s...\n", tAln);
    printLog(evidence, j, "    full evidence %100.100s...\n", rAln);

    printLog(evidence, j, "    full template ...%100.100s\n", tAln + align.alignmentLength - 100);
    printLog(evidence, j, "    full evidence ...%100.100s\n", rAln + align.alignmentLength - 100);
#endif

    tAln[align.alignmentLength - eSkip] = 0;   //  Truncate the alignment strings for printing.
    rAln[align.alignmentLength - eSkip] = 0;

#if 0
    printLog(evidence, j, "         template %100.100s...\n", tAln + bSkip);
    printLog(evidence, j, "         evidence %100.100s...\n", rAln + bSkip);

    printLog(evidence, j, "         template ...%100.100s\n", tAln + align.alignmentLength - eSkip - 100);
    printLog(evidence, j, "         evidence ...%100.100s\n", rAln + align.alignmentLength - eSkip - 100);
#endif

    //  Convert the alignment into tags.

    if (align.alignmentLength > bSkip + eSkip + minOlapLength)
      tagList[j] = getAlignTags(rAln + bSkip, rBgn + rbSkip, evidence[j].readLength,
                                tAln + bSkip, tBgn + tbSkip, evidence[0].readLength,
                                align.alignmentLength - bSkip - eSkip);

    delete [] tAln;
    delete [] rAln;

    edlibFreeAlignResult(align);
  }

  return(tagList);
}

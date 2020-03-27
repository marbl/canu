
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
#include "falconConsensus-alignTag.H"
#include "falconConsensus-msa.H"

#include "intervalList.H"


#undef DEBUG
#undef DEBUG_VERBOSE


falconData *
falconConsensus::getConsensus(uint32         tagsLen,                //  Number of evidence reads
                              alignTagList **tags,                   //  Alignment tags
                              uint32         templateLen) {          //  Length of template read

  //  If no tags, return an empty result.

  if (tagsLen == 0)
    return(new falconData);

  msa.resize(templateLen);

  //  For each alignment position, insert the alignment tag to msa

  int32  t_pos   = 0;

  for (uint32 i=0; i<tagsLen; i++) {
    if (tags[i] == NULL)
      continue;

    for (uint32 j=0; j<tags[i]->numberOfTags(); j++) {
      alignTag *tag = (*tags[i])[j];

      if (tag->delta == 0) {
        t_pos = tag->t_pos;
        msa[t_pos]->coverage++;
      }

#ifdef DEBUG
      fprintf(stderr, "Processing position %d in sequence %d (in msa it is column %d with cov %d) with delta %d and current size is %d\n", j, i, t_pos, msa[t_pos]->coverage, tag->delta, msa[t_pos]->size);
#endif

      // Assume t_pos was set on earlier iteration.
      // (Otherwise, use its initial value, which might be an error. ~cd)

      assert(tag->delta < uint16MAX);

      msa[t_pos]->increaseDeltaGroup(tag->delta);

      uint32 base = 4;

      switch (tag->q_base) {
        case 'A':  base = 0;  break;
        case 'C':  base = 1;  break;
        case 'G':  base = 2;  break;
        case 'T':  base = 3;  break;
        case '-':  base = 4;  break;
        default :  base = 4;  break;
      }

      if (j > 0)    assert(tag->p_t_pos >= 0);


      //  Update the column

      assert(tag->delta < msa[t_pos]->deltaLen);
      align_tag_col_t  &col = msa[t_pos]->delta[tag->delta]->base[base];

      bool updated = false;

      col.count += 1;

      //  Search for a matching column.  If found, add one.  If not found, make a new entry.

      for (int32 kk=0; kk<col.n_link; kk++) {
        if ((tag->p_t_pos   == col.p_t_pos[kk]) &&
            (tag->p_delta   == col.p_delta[kk]) &&
            (tag->p_q_base  == col.p_q_base[kk])) {
          col.link_count[kk]++;
          updated = true;
          break;
        }
      }

      if (updated == false)
        col.addEntry(tag);

#ifdef DEBUG
      fprintf(stderr, "Updating column from seq %d at position %d in column %d base pos %d base %d to be %c and length is %d\n", i, j, t_pos, base, tag->p_t_pos, tag->p_q_base, msa[t_pos]->deltaLen);
#endif
    }

    updateRSS();

    delete tags[i];
    tags[i] = NULL;
  }

  //  Done with the tags.

  delete [] tags;
  tags = NULL;

  // propogate score throught the alignment links, setup backtracking information

  align_tag_col_t *g_best_aln_col = NULL;
  int32            g_best_t_pos   = -1;
  double           g_best_score   = -1;  //  Might be a magic value.

  //  Over every template base,
  //  And every delta position,
  //  And every base at that position
  //  Search links to previous columns, remember the highest scoring one,
  //  Then remember the highest scoring link for each

  for (uint32 i=0; i<templateLen; i++) {
    for (uint32 j=0; j<msa[i]->deltaLen; j++) {
      for (uint32 kk=0; kk<5; kk++) {
        align_tag_col_t *aln_col = msa[i]->delta[j]->base + kk;

        aln_col->score    = -1;  //  Probably needs to be the same magic value as above.

        double best_score = -1;  //  Magic too?

        //fprintf(stderr, "Processing consensus template %d which as %d delta and on base %d i pulled up col %d with %d links and best %d %d %d\n",
        //        i, j, kk, aln_col, aln_col->n_link, aln_col->best_p_t_pos, aln_col->best_p_delta, aln_col->best_p_q_base);

        //  Search links to previous columns, remember the highest scoring one.

        for (uint32 ck=0; ck<aln_col->n_link; ck++) {
          int32 pi  = aln_col->p_t_pos[ck];
          int32 pj  = aln_col->p_delta[ck];
          int32 pkk = 4;

          switch (aln_col->p_q_base[ck]) {
            case 'A': pkk = 0; break;
            case 'C': pkk = 1; break;
            case 'G': pkk = 2; break;
            case 'T': pkk = 3; break;
            case '-': pkk = 4; break;
            default : pkk = 4; break;
          }

          //  Score is just our link weight, possibly with the previous column's score, and
          //  penalizing for coverage.

          double score = aln_col->link_count[ck] - msa[i]->coverage * 0.5;

          if ((aln_col->p_t_pos[ck] != -1) &&
              (pj <= msa[pi]->deltaLen))
            score += msa[pi]->delta[pj]->base[pkk].score;

          //  Save best score.

#ifdef DEBUG_VERBOSE
          fprintf(stderr, "best_score %f at pi %d pj %d pkk %d -- score %f\n", score, pi, pj, pkk, score);
#endif

          if (best_score < score) {
            aln_col->best_p_t_pos    = pi;
            aln_col->best_p_delta    = pj;
            aln_col->best_p_q_base   = pkk;
            best_score               = score;

#ifdef DEBUG
            fprintf(stderr, "best_score %f at pi %d pj %d pkk %d\n", score, pi, pj, pkk);
#endif
          }
        }  //  Over all links

        aln_col->score = best_score;

        if (g_best_score < best_score) {
          g_best_aln_col = aln_col;
          g_best_t_pos   = i;
          g_best_score   = best_score;
        }
      }
    }
  }

  updateRSS();

  //  If there is no best found, then there was no evidence aligned to the template.  Prior versions
  //  used to fail here:
  //    assert(g_best_t_pos   != -1);
  //    assert(g_best_aln_col != NULL);
  //  but BPW thinks we can just fall through and return an empty sequence.  It should actually never
  //  happen - the original falcon_sense had the asserts and worked fine; this version triggered the
  //  asserts when short overlaps were falsely removed.  But really, I think it's a coordination
  //  between this code and the thing that generates layouts.

  //  Reconstruct the sequences.

  falconData *fd = new falconData(templateLen * 2 + 1);

  int32      i  = g_best_t_pos;
  int32      j  = 0;
  uint32     kk = (g_best_aln_col == NULL) ? 0 : g_best_aln_col->best_p_q_base;

  while ((i != -1) && (fd->len < templateLen * 2)) {
    char  bb = '-';

    switch (kk) {
      case 0: bb = (msa[i]->coverage <= minOutputCoverage) ? 'a' : 'A'; break;
      case 1: bb = (msa[i]->coverage <= minOutputCoverage) ? 'c' : 'C'; break;
      case 2: bb = (msa[i]->coverage <= minOutputCoverage) ? 'g' : 'G'; break;
      case 3: bb = (msa[i]->coverage <= minOutputCoverage) ? 't' : 'T'; break;
      case 4: bb =                                                 '-'; break;
    }

    if (bb != '-') {
      fd->seq[fd->len] = bb;
      fd->eqv[fd->len] = (msa[i]->coverage == g_best_aln_col->count) ? (40) : (-10 * log((msa[i]->coverage - g_best_aln_col->count + 1) / (double)msa[i]->coverage));
      fd->pos[fd->len] = i;

#ifdef DEBUG_VERBOSE
      //fprintf(stderr, "seq %5u pos %5u '%c' cov %3u eqv %4d\n",
      //        fd->len, i, bb, msa[i]->coverage, fd->eqv[fd->len]);
      fprintf(stderr, "seq %5u pos %5u '%c' cov %3u\n",
              fd->len, i, bb, msa[i]->coverage);
#endif

      if (fd->eqv[fd->len] > 40)
        fd->eqv[fd->len] = 40;

      assert(0 <= fd->eqv[fd->len]);

      fd->len++;
    }

    i   = g_best_aln_col->best_p_t_pos;
    j   = g_best_aln_col->best_p_delta;
    kk  = g_best_aln_col->best_p_q_base;

    if (i != -1)
      g_best_aln_col = msa[i]->delta[j]->base + kk;
  }

  fd->seq[fd->len] = 0;

  // reverse the sequence

  reverse(fd->seq, fd->seq + fd->len);
  reverse(fd->eqv, fd->eqv + fd->len);
  reverse(fd->pos, fd->pos + fd->len);

  updateRSS();

  return(fd);
}



falconData *
falconConsensus::generateConsensus(falconInput   *evidence,
                                   uint32         evidenceLen) {

  setRSS();

  alignTagList **tags = alignReadsToTemplate(evidence, evidenceLen, minOlapIdentity, minOlapLength, restrictToOverlap);

  updateRSS();

  return(getConsensus(evidenceLen, tags, evidence[0].readLength));
}



uint64
falconConsensus::estimateMemoryUsage(uint32 UNUSED(evidenceLen),
                                     uint64        nBasesInOlaps,
                                     uint32        templateLen) {

  //  For evidence, each aligned base makes an alignTag, then 2 bytes for the read itself.
  //  This _should_ be a vast over-estimate, but it is just barely the actual size.
  //
  //  Then during consensus, each base in the template allocates:
  //     an msa_delta_group_t           each of which allocates:
  //     at least 8 msa_base_group_t    each of which allocates:    (assume 16 max)
  //     at least 8 64-bit words.                                   (assume 16 max)
  //
  //  Based on a single long nanopore read, using 16 instead of 8 is an overestimate.  I don't
  //  understand what makes these grow.

  uint64  perEvidence = sizeof(alignTag) + 2;
  uint64  perTemplate = (sizeof(msa_delta_group_t) +
                         16 * (sizeof(msa_base_group_t) +
                               24 * (sizeof(int32) + sizeof(uint16) + sizeof(char) + sizeof(uint16))));
  uint64  slush       = 500 * 1024 * 1024;

  //fprintf(stderr, "evidence  %4lu x %9lu bases = %9lu %9lu MB\n",
  //        perEvidence, nBasesInOlaps, nBasesInOlaps * perEvidence, nBasesInOlaps * perEvidence >> 20);
  //fprintf(stderr, "template  %4lu x %9u bases = %9lu %9lu MB\n",
  //        perTemplate, templateLen,   templateLen * perTemplate,   templateLen * perTemplate >> 20);

  return(nBasesInOlaps * perEvidence + templateLen * perTemplate + slush);
}



uint64
falconConsensus::analyzeLength(tgTig            *layout,
                               uint32           &correctedLength,
                               uint64           &memoryRequired) {

  //  Estimate the length of the corrected read, using overlap position and depth.

  intervalList<int32>   coverage;
  uint64                basesInOlaps = 0;

  for (uint32 ii=0; ii<layout->numberOfChildren(); ii++) {
    tgPosition *pos = layout->getChild(ii);

    coverage.add(pos->_min, pos->_max - pos->_min);

    basesInOlaps += pos->_max - pos->_min;
  }

  intervalDepth<int32>  depth(coverage);

  int32    bgn       = INT32_MAX;
  int32    corLen    = 0;

  for (uint32 dd=0; dd<depth.numberOfIntervals(); dd++) {
    if (depth.depth(dd) < minOutputCoverage) {
      bgn = INT32_MAX;
      continue;
    }

    if (bgn == INT32_MAX)
      bgn = depth.lo(dd);

    if (corLen < depth.hi(dd) - bgn)
      corLen = depth.hi(dd) - bgn;
  }

  //  Save our results.

  correctedLength = corLen;
  memoryRequired  = estimateMemoryUsage(layout->numberOfChildren(),
                                        basesInOlaps,
                                        layout->length());

  return(memoryRequired);
}


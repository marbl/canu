// This file is part of A2Amapper.
// Copyright (c) 2005 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "bio++.H"
#include "atac.H"


//  While loading matches, we compute the mapped length and covered
//  length.

matchList::matchList(char *filename, char matchOrRun, bool saveLine) {

  if ((matchOrRun != 'u') && (matchOrRun != 'x') && (matchOrRun != 'r') && (matchOrRun != 'm')) {
    fprintf(stderr, "matchList::matchList()-- Invalid value '%c' for matchOrRun, should be:\n");
    fprintf(stderr, "                         'u' -- ungapped matches, with mismatches\n");
    fprintf(stderr, "                         'x' -- ungapped matches, exact\n");
    fprintf(stderr, "                         'm' -- ungapped matches, both 'u' and 'x'\n");
    fprintf(stderr, "                         'r' -- runs (gapped)\n");
    exit(1);
  }

  errno = 0;
  FILE   *inFile = fopen(filename, "r");
  if (errno)
    fprintf(stderr, "matchList::matchList()-- failed to load %s: %s\n", filename, strerror(errno)), exit(1);

  //  Read the preamble, look for our data sources.  This leaves us with
  //  the first match in the inLine, and fills in file1 and file2.
  //
  char    inLine[1024];
  readHeader(inLine, inFile, _file1, _file2, 0L);

  //  Open some FastAWrappers for each of the files -- we use these
  //  only to get the length of the sequence.
  //
  _seq1 = new FastAWrapper(_file1);
  _seq2 = new FastAWrapper(_file2);

  _seq1->openIndex();
  _seq2->openIndex();

  _matchesLen = 0;
  _matchesMax = 32 * 1048576;
  _matches    = new match_t [_matchesMax];

  //  For the coverage to work correctly, we need to either have one
  //  intervalList per input sequence, or build a table of the chained
  //  sequence positions.
  //
  u64bit  *offset1 = new u64bit [_seq1->getNumberOfSequences()];
  u64bit  *offset2 = new u64bit [_seq2->getNumberOfSequences()];

  u32bit   length1 = 0;
  u32bit   length2 = 0;

  offset1[0] = 1000000;
  for (u32bit i=1; i<_seq1->getNumberOfSequences(); i++)
    offset1[i] = offset1[i-1] + _seq1->sequenceLength(i-1) + 1;

  offset2[0] = 1000000;
  for (u32bit i=1; i<_seq2->getNumberOfSequences(); i++)
    offset2[i] = offset2[i-1] + _seq2->sequenceLength(i-1) + 1;

  for (u32bit i=0; i<_seq1->getNumberOfSequences(); i++)
    length1 += _seq1->sequenceLength(i);

  for (u32bit i=0; i<_seq2->getNumberOfSequences(); i++)
    length2 += _seq2->sequenceLength(i);

  intervalList  intervalA;
  intervalList  intervalB;

  while (!feof(inFile)) {
    if (inLine[0] == 'M') {
      splitToWords  S(inLine);

      if ((S[1][0] == matchOrRun) ||
          ((matchOrRun == 'm') && ((S[1][0] == 'u') || (S[1][0] == 'x')))) {

        //  Save the name/label
        if (_name1[0] == 0)
          decodeMatchNames(S, _name1, _name2);

        u32bit  iid1=0, pos1=0, len1=0, fwd1=0;
        u32bit  iid2=0, pos2=0, len2=0, fwd2=0;
        decodeMatch(S, iid1, pos1, len1, fwd1, iid2, pos2, len2, fwd2);

        if ((pos1) > _seq1->sequenceLength(iid1) || (pos1 + len1) > _seq1->sequenceLength(iid1)) {
          chomp(inLine);
          fprintf(stderr, "Match longer than sequence in 1: "u32bitFMT" %s\n", _seq1->sequenceLength(iid1), inLine);
        }

        if ((pos2) > _seq2->sequenceLength(iid2) || (pos2 + len2) > _seq2->sequenceLength(iid2)) {
          chomp(inLine);
          fprintf(stderr, "Match longer than sequence in 2: "u32bitFMT" %s\n", _seq2->sequenceLength(iid2), inLine);
        }

        if ((iid1 >= _seq1->getNumberOfSequences()) || (iid2 >= _seq2->getNumberOfSequences())) {
          chomp(inLine);
          fprintf(stderr, "Match references invalid sequence iid: %s\n", inLine);
        } else {
          intervalA.add(offset1[iid1] + (u64bit)pos1, (u64bit)len1);
          intervalB.add(offset2[iid2] + (u64bit)pos2, (u64bit)len2);

          //  Add it to our list of matches
          //
          if (_matchesLen > _matchesMax) {
            fprintf(stderr, "SORRY!  I don't feel like reallocating matches.  Increase\n");
            fprintf(stderr, "the preallocated size in %s\n", __FILE__);
            exit(1);
          }

          _matches[_matchesLen].matchiid = _matchesLen;
          _matches[_matchesLen].matchuid = strtou32bit(S[2], 0L);

          _matches[_matchesLen].iid1 = iid1;
          _matches[_matchesLen].pos1 = pos1;
          _matches[_matchesLen].len1 = len1;
          _matches[_matchesLen].fwd1 = fwd1;

          _matches[_matchesLen].iid2 = iid2;
          _matches[_matchesLen].pos2 = pos2;
          _matches[_matchesLen].len2 = len2;
          _matches[_matchesLen].fwd2 = fwd2;

          _matchesLen++;
        }
      }
    }

    fgets(inLine, 1024, inFile);
  }

  fprintf(stderr, "numberOfItems  %c "u64bitFMT"\n",
          matchOrRun, (u64bit)_matchesLen);

  fprintf(stderr, "totalLength    A "u64bitFMT" B "u64bitFMT"\n",
          (u64bit)length1,
          (u64bit)length2);

  fprintf(stderr, "intervalLength A "u64bitFMT" B "u64bitFMT"\n",
          (u64bit)intervalA.sumOfLengths(),
          (u64bit)intervalB.sumOfLengths());

  intervalA.merge();
  intervalB.merge();

  fprintf(stderr, "coveredLength  A "u64bitFMT" B "u64bitFMT"\n",
          (u64bit)intervalA.sumOfLengths(),
          (u64bit)intervalB.sumOfLengths());
}




static
int
sort1_(const void *a, const void *b) {
  const match_t *A = (const match_t *)a;
  const match_t *B = (const match_t *)b;

  if (A->iid1 < B->iid1)  return(-1);
  if (A->iid1 > B->iid1)  return(1);
  if (A->pos1 < B->pos1)  return(-1);
  if (A->pos1 > B->pos1)  return(1);
  if (A->len1 > B->len1)  return(-1);
  if (A->len1 < B->len1)  return(1);
#if 0
  //  disabled so that clumpMaker can use these sorts
  if (A->fwd1 > B->fwd1)  return(-1);
  if (A->fwd1 < B->fwd1)  return(1);
#endif
  if (A->iid2 < B->iid2)  return(-1);
  if (A->iid2 > B->iid2)  return(1);
  if (A->pos2 < B->pos2)  return(-1);
  if (A->pos2 > B->pos2)  return(1);
  if (A->len2 > B->len2)  return(-1);
  if (A->len2 < B->len2)  return(1);
  return(0);
}

static
int
sort2_(const void *a, const void *b) {
  const match_t *A = (const match_t *)a;
  const match_t *B = (const match_t *)b;

  if (A->iid2 < B->iid2)  return(-1);
  if (A->iid2 > B->iid2)  return(1);
  if (A->pos2 < B->pos2)  return(-1);
  if (A->pos2 > B->pos2)  return(1);
  if (A->len2 > B->len2)  return(-1);
  if (A->len2 < B->len2)  return(1);
#if 0
  //  disabled so that clumpMaker can use these sorts
  if (A->fwd2 > B->fwd2)  return(-1);
  if (A->fwd2 < B->fwd2)  return(1);
#endif
  if (A->iid1 < B->iid1)  return(-1);
  if (A->iid1 > B->iid1)  return(1);
  if (A->pos1 < B->pos1)  return(-1);
  if (A->pos1 > B->pos1)  return(1);
  if (A->len1 > B->len1)  return(-1);
  if (A->len1 < B->len1)  return(1);

  return(0);
}


void
matchList::sort1(u32bit first, u32bit len) {
  if (len == 0) len = _matchesLen;
  qsort(_matches + first, len, sizeof(match_t), sort1_);
}


void
matchList::sort2(u32bit first, u32bit len) {
  if (len == 0) len = _matchesLen;
  qsort(_matches + first, len, sizeof(match_t), sort2_);
}

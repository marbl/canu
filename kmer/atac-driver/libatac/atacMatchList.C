// This file is part of A2Amapper.
// Copyright (c) 2005, 2006 J. Craig Venter Institute
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


atacMatchList::atacMatchList(char const *filename, char matchOrRun, FILE *headerOut) {

  _labelA[0]  = 0;
  _labelB[0]  = 0;
  _fileA[0]   = 0;
  _fileB[0]   = 0;
  _seqA       = 0L;
  _seqB       = 0L;
  _matchesLen = 0;
  _matchesMax = 0;
  _matches    = 0L;

  if (filename == 0L)
    return;

  if ((matchOrRun != 'u') && (matchOrRun != 'x') && (matchOrRun != 'r') && (matchOrRun != 'm')) {
    fprintf(stderr, "atacMatchList::atacMatchList()-- Invalid value '%c' for matchOrRun, should be:\n", matchOrRun);
    fprintf(stderr, "                                 'u' -- ungapped matches, with mismatches\n");
    fprintf(stderr, "                                 'x' -- ungapped matches, exact\n");
    fprintf(stderr, "                                 'm' -- ungapped matches, both 'u' and 'x'\n");
    fprintf(stderr, "                                 'r' -- runs (gapped)\n");
    exit(1);
  }

  FILE *inFile = stdin;
  if ((filename != 0L) && (strcmp(filename, "-") != 0)) {
    errno = 0;
    inFile = fopen(filename, "r");
    if (errno)
      fprintf(stderr, "atacMatchList::atacMatchList()-- failed to load %s: %s\n", filename, strerror(errno)), exit(1);
  }

  //  Read the preamble, look for our data sources.  This leaves us with
  //  the first match in the inLine, and fills in fileA and fileB.
  //
  char    inLine[1024];
  readHeader(inLine, inFile, _fileA, _fileB, headerOut);

  //  Open some FastAWrappers for each of the files -- we use these
  //  only to get the length of the sequence.
  //
  _seqA = 0L;
  _seqB = 0L;
  if (_fileA && _fileA[0]) {
    _seqA = new FastAWrapper(_fileA);
    _seqA->openIndex();
  }
  if (_fileB && _fileB[0]) {
    _seqB = new FastAWrapper(_fileB);
    _seqB->openIndex();
  }

  _labelA[0] = 0;
  _labelB[0] = 0;

  _matchesLen = 0;
  _matchesMax = 1 * 1048576;
  _matches    = new atacMatch [_matchesMax];

  while (!feof(inFile)) {
    if (inLine[0] == 'M') {
      splitToWords  S(inLine);

      if ((S[1][0] == matchOrRun) ||
          ((matchOrRun == 'm') && ((S[1][0] == 'u') || (S[1][0] == 'x')))) {

        //  Save the name/label
        if (_labelA[0] == 0) {
          decodeAtacName(S[4], _labelA);
          decodeAtacName(S[8], _labelB);
        }

        u32bit  iid1=0, pos1=0, len1=0, fwd1=0;
        u32bit  iid2=0, pos2=0, len2=0, fwd2=0;
        decodeMatch(S, iid1, pos1, len1, fwd1, iid2, pos2, len2, fwd2);

        bool matchOK = true;
        if (_seqA && _seqB) {
          if ((pos1) > _seqA->sequenceLength(iid1) || (pos1 + len1) > _seqA->sequenceLength(iid1)) {
            chomp(inLine);
            fprintf(stderr, "Match longer than sequence (by "u32bitFMT"bp) in 1: seqLen="u32bitFMTW(8)" %s\n",
                    pos1 + len1 - _seqA->sequenceLength(iid1),
                    _seqA->sequenceLength(iid1), inLine);
            matchOK = false;
          }

          if ((pos2) > _seqB->sequenceLength(iid2) || (pos2 + len2) > _seqB->sequenceLength(iid2)) {
            chomp(inLine);
            fprintf(stderr, "Match longer than sequence (by "u32bitFMT"bp) in 2: seqLen="u32bitFMTW(8)" %s\n",
                    pos2 + len2 - _seqB->sequenceLength(iid2),
                    _seqB->sequenceLength(iid2), inLine);
            matchOK = false;
          }

          if ((iid1 >= _seqA->getNumberOfSequences()) || (iid2 >= _seqB->getNumberOfSequences())) {
            chomp(inLine);
            fprintf(stderr, "Match references invalid sequence iid: %s\n", inLine);
            matchOK = false;
          }
        }

        if (matchOK) {

          //  Add it to our list of matches
          //
          if (_matchesLen > _matchesMax) {
            _matchesMax *= 2;
            atacMatch *n = new atacMatch [_matchesMax];
            memcpy(n, _matches, sizeof(atacMatch));
            delete [] _matches;
            _matches = n;
          }

          strncpy(_matches[_matchesLen].matchuid,  S[2], 16);
          strncpy(_matches[_matchesLen].parentuid, S[3], 16);

          _matches[_matchesLen].matchuid[15]  = 0;
          _matches[_matchesLen].parentuid[15] = 0;

          _matches[_matchesLen].matchiid = _matchesLen;

          _matches[_matchesLen].type[0] = 0;
          _matches[_matchesLen].type[1] = 0;
          _matches[_matchesLen].type[2] = 0;
          _matches[_matchesLen].type[3] = 0;

          _matches[_matchesLen].type[0] = S[1][0];
          _matches[_matchesLen].type[1] = S[1][1];
          if (S[1][1])
            _matches[_matchesLen].type[2] = S[1][2];

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
}


//  atacMatchOrder



void
atacMatchOrder::mergeMatches(atacMatch *l, atacMatch *r, u32bit mergeuid) {
  atacMatch   n;

  //  Create a new match record for the merged match.  We could
  //  probably do this inplace in l.

  //  Copy all the defaults from L first.  This copies most of the stuff.
  //
  memcpy(&n, l, sizeof(atacMatch));

  sprintf(n.matchuid, "merge"u32bitFMT, mergeuid);

  n.len1 = (r->pos1 + r->len1) - (l->pos1);
  n.len2 = n.len1;

  if (r->fwd2 == false)
    n.pos2 = r->pos2;

  n.fwd2 = r->fwd2;

  //  Update l with the new contents.

  memcpy(l, &n, sizeof(atacMatch));

  //  Remove the r match from our set.  The hardest part is figuring
  //  out what index the r match is at.  The easiest way to do that is
  //  the most inefficient (start at zero, when we find the r match,
  //  start updating).  The quickest way (given we want an array)
  //  makes us trust our index.
  //
  _matchesLen--;
  for (u32bit idx = index(r->matchiid); idx < _matchesLen; idx++) {
    _matches[idx]                           = _matches[idx+1];
    _matchIIDtoIdx[_matches[idx]->matchiid] = idx;
  }
}



static
int
sortA_(const void *a, const void *b) {
  const atacMatch *A = *(const atacMatch **)a;
  const atacMatch *B = *(const atacMatch **)b;

  if (A->iid1 < B->iid1)  return(-1);
  if (A->iid1 > B->iid1)  return(1);
  if (A->pos1 < B->pos1)  return(-1);
  if (A->pos1 > B->pos1)  return(1);
  if (A->len1 > B->len1)  return(-1);
  if (A->len1 < B->len1)  return(1);
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
sortB_(const void *a, const void *b) {
  const atacMatch *A = *(const atacMatch **)a;
  const atacMatch *B = *(const atacMatch **)b;

  if (A->iid2 < B->iid2)  return(-1);
  if (A->iid2 > B->iid2)  return(1);
  if (A->pos2 < B->pos2)  return(-1);
  if (A->pos2 > B->pos2)  return(1);
  if (A->len2 > B->len2)  return(-1);
  if (A->len2 < B->len2)  return(1);
  if (A->iid1 < B->iid1)  return(-1);
  if (A->iid1 > B->iid1)  return(1);
  if (A->pos1 < B->pos1)  return(-1);
  if (A->pos1 > B->pos1)  return(1);
  if (A->len1 > B->len1)  return(-1);
  if (A->len1 < B->len1)  return(1);

  return(0);
}

static
int
sortdiagonal_(const void *a, const void *b) {
  const atacMatch *A = *(const atacMatch **)a;
  const atacMatch *B = *(const atacMatch **)b;

  if (A->iid2 < B->iid2)  return(-1);
  if (A->iid2 > B->iid2)  return(1);
  if (A->iid1 < B->iid1)  return(-1);
  if (A->iid1 > B->iid1)  return(1);
  if (A->fwd2 < B->fwd2)  return(-1);
  if (A->fwd2 > B->fwd2)  return(1);

  //  We're now in the same sequence pair with the same orientation.

  //  So much easier if we use signed math.

  //  This works for forward matches
  s32bit dA = (s32bit)A->pos2 - (s32bit)A->pos1;
  s32bit dB = (s32bit)B->pos2 - (s32bit)B->pos1;

  if (A->fwd2 == 0) {
    //  OK, so not the greatest diagonal computation ever.  We end up
    //  with a gigantic discontinuity at the origin, but we don't
    //  care, just as long as the diagonals are distinct.
    //
    dA = (s32bit)A->pos2 - (1000000000 - (s32bit)(A->pos2 + A->len2));
    dB = (s32bit)B->pos2 - (1000000000 - (s32bit)(B->pos2 + B->len2));
  }

  if (dA < dB)  return(-1);
  if (dA > dB)  return(1);

  //  This is just candy; might make things easier later
  if (A->pos1 < B->pos1)  return(-1);
  if (A->pos1 > B->pos1)  return(1);
  if (A->len1 > B->len1)  return(-1);
  if (A->len1 < B->len1)  return(1);

  return(0);
}

static
int
sortmatchuid_(const void *a, const void *b) {
  const atacMatch *A = *(const atacMatch **)a;
  const atacMatch *B = *(const atacMatch **)b;

  int  r = strcmp(A->matchuid, B->matchuid);
  if (r < 0)  return(-1);
  if (r > 0)  return(1);
  r = strcmp(A->parentuid, B->parentuid);
  if (r < 0)  return(-1);
  if (r > 0)  return(1);

  return(0);
}

static
int
sortparentuid_(const void *a, const void *b) {
  const atacMatch *A = *(const atacMatch **)a;
  const atacMatch *B = *(const atacMatch **)b;

  int  r = strcmp(A->parentuid, B->parentuid);
  if (r < 0)  return(-1);
  if (r > 0)  return(1);
  r = strcmp(A->matchuid, B->matchuid);
  if (r < 0)  return(-1);
  if (r > 0)  return(1);
  
  return(0);
}


void
atacMatchOrder::sortA(u32bit first, u32bit len) {
  if (len == 0) len = _matchesLen;
  qsort(_matches + first, len, sizeof(atacMatch*), sortA_);
  updateIndex();
}

void
atacMatchOrder::sortB(u32bit first, u32bit len) {
  if (len == 0) len = _matchesLen;
  qsort(_matches + first, len, sizeof(atacMatch*), sortB_);
  updateIndex();
}

void
atacMatchOrder::sortDiagonal(u32bit first, u32bit len) {
  if (len == 0) len = _matchesLen;
  qsort(_matches + first, len, sizeof(atacMatch*), sortdiagonal_);
  updateIndex();
}

void
atacMatchOrder::sortMatchUID(u32bit first, u32bit len) {
  if (len == 0) len = _matchesLen;
  qsort(_matches + first, len, sizeof(atacMatch*), sortmatchuid_);
  updateIndex();
}

void
atacMatchOrder::sortParentUID(u32bit first, u32bit len) {
  if (len == 0) len = _matchesLen;
  qsort(_matches + first, len, sizeof(atacMatch*), sortparentuid_);
  updateIndex();
}

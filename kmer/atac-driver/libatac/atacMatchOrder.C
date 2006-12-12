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
  const atacMatch *A = *(const atacMatch * const *)a;
  const atacMatch *B = *(const atacMatch * const *)b;

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
  const atacMatch *A = *(const atacMatch * const *)a;
  const atacMatch *B = *(const atacMatch * const *)b;

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
  const atacMatch *A = *(const atacMatch * const *)a;
  const atacMatch *B = *(const atacMatch * const *)b;

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
  const atacMatch *A = *(const atacMatch * const *)a;
  const atacMatch *B = *(const atacMatch * const *)b;

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
  const atacMatch *A = *(const atacMatch * const *)a;
  const atacMatch *B = *(const atacMatch * const *)b;

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

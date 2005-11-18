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

#include "overlap.H"

int
sortMatches1(const void *a, const void *b) {
  const match_t *A = *((const match_t * const *)a);
  const match_t *B = *((const match_t * const *)b);

  if (A->iid1 < B->iid1)  return(-1);
  if (A->iid1 > B->iid1)  return(1);
  if (A->pos1 < B->pos1)  return(-1);
  if (A->pos1 > B->pos1)  return(1);
  if (A->len1 > B->len1)  return(-1);
  if (A->len1 < B->len1)  return(1);
  if (A->fwd1 > B->fwd1)  return(-1);
  if (A->fwd1 < B->fwd1)  return(1);
  return(0);
}

int
sortMatches2(const void *a, const void *b) {
  const match_t *A = *((const match_t * const *)a);
  const match_t *B = *((const match_t * const *)b);

  if (A->iid2 < B->iid2)  return(-1);
  if (A->iid2 > B->iid2)  return(1);
  if (A->pos2 < B->pos2)  return(-1);
  if (A->pos2 > B->pos2)  return(1);
  if (A->len2 > B->len2)  return(-1);
  if (A->len2 < B->len2)  return(1);
  if (A->fwd2 > B->fwd2)  return(-1);
  if (A->fwd2 < B->fwd2)  return(1);
  return(0);
}



int
spanCompare(const void *a, const void *b) {
  const span_t *A = *((const span_t * const *)a);
  const span_t *B = *((const span_t * const *)b);

  if (A->_iid < B->_iid)  return(-1);
  if (A->_iid > B->_iid)  return(1);
  if (A->_beg < B->_beg)  return(-1);
  if (A->_beg > B->_beg)  return(1);
  return(0);
}

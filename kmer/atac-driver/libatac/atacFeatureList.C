// This file is part of A2Amapper.
// Copyright (c) 2006 J. Craig Venter Institute
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


atacFeatureList::atacFeatureList() {
  _featuresLen = 0;
  _featuresMax = 256;
  _features    = new atacFeature [_featuresMax];
}

atacFeatureList::~atacFeatureList() {
  delete [] _features;
}

void
atacFeatureList::add(atacFeature &m) {

  if (_featuresLen >= _featuresMax) {
    _featuresMax <<= 2;
    atacFeature  *A = new atacFeature [_featuresMax];
    memcpy(A, _features, sizeof(atacFeature) * _featuresLen);
    delete [] _features;
    _features = A;
  }

  memcpy(&_features[_featuresLen], &m, sizeof(atacFeature));

  _features[_featuresLen].featureiid = _featuresLen++;
}



static
int
sort_(const void *a, const void *b) {
  const atacFeature *A = (const atacFeature *)a;
  const atacFeature *B = (const atacFeature *)b;

  if (A->iid < B->iid)  return(-1);
  if (A->iid > B->iid)  return(1);
  if (A->pos < B->pos)  return(-1);
  if (A->pos > B->pos)  return(1);
  if (A->len > B->len)  return(-1);
  if (A->len < B->len)  return(1);
  return(0);
}

static
int
sortfeatureuid_(const void *a, const void *b) {
  const atacFeature *A = (const atacFeature *)a;
  const atacFeature *B = (const atacFeature *)b;

  int  r = strcmp(A->featureuid, B->featureuid);
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
  const atacFeature *A = (const atacFeature *)a;
  const atacFeature *B = (const atacFeature *)b;

  int  r = strcmp(A->parentuid, B->parentuid);
  if (r < 0)  return(-1);
  if (r > 0)  return(1);
  r = strcmp(A->featureuid, B->featureuid);
  if (r < 0)  return(-1);
  if (r > 0)  return(1);
  
  return(0);
}


void
atacFeatureList::sort(u32bit first, u32bit len) {
  if (len == 0) len = _featuresLen;
  qsort(_features + first, len, sizeof(atacFeature), sort_);
}

void
atacFeatureList::sortFeatureUID(u32bit first, u32bit len) {
  if (len == 0) len = _featuresLen;
  qsort(_features + first, len, sizeof(atacFeature), sortfeatureuid_);
}

void
atacFeatureList::sortParentUID(u32bit first, u32bit len) {
  if (len == 0) len = _featuresLen;
  qsort(_features + first, len, sizeof(atacFeature), sortparentuid_);
}

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


static
u32bit
decodeAtacName(char *atac,
               char *label) {
  if (label) {
    while (*atac && (*atac != ':'))
      *label++ = *atac++;
    *label = 0;
  } else {
    while (*atac && (*atac != ':'))
      atac++;
  }
  if (*atac)
    return(strtou32bit(atac+1, 0L));
  return(~u32bitZERO);
}


atacFeature::atacFeature(char *line) {
  decode(line);
}


atacFeature::atacFeature(char *fuid,
                         char *puid,
                         u32bit fiid,
                         char *t,
                         u32bit i, u32bit p, u32bit l) {

  strcpy(featureuid, fuid);
  strcpy(parentuid, puid);

  featureuid[15] = 0;
  parentuid[15]  = 0;

  featureiid = fiid;

  type[0] = 0;
  type[1] = 0;
  type[2] = 0;
  type[3] = 0;
  strcpy(type, t);

  iid = i;
  pos = p;
  len = l;
}


void
atacFeature::decode(char *line) {

  splitToWords  W(line);

  strcpy(featureuid, W[2]);
  strcpy(parentuid, W[3]);

  featureuid[15] = 0;
  parentuid[15]  = 0;

  featureiid = 0;

  type[0] = 0;
  type[1] = 0;
  type[2] = 0;
  type[3] = 0;
  strcpy(type, W[1]);

  iid = decodeAtacName(W[4], 0L);
  pos = strtou32bit(W[5], 0L);
  len = strtou32bit(W[6], 0L);
}


bool
atacFeature::sanity(seqFile *A, char *inLine) {

  bool  featureOK = true;

  if (A) {
    if ((pos) > A->sequenceLength(iid) || (pos + len) > A->sequenceLength(iid)) {
      chomp(inLine);
      fprintf(stderr, "Feature longer than sequence (by "u32bitFMT"bp): seqLen="u32bitFMTW(8)" %s\n",
              pos + len - A->sequenceLength(iid),
              A->sequenceLength(iid), inLine);
      featureOK = false;
    }

    if (iid >= A->getNumberOfSequences()) {
      chomp(inLine);
      fprintf(stderr, "Feature references invalid sequence iid: %s\n", inLine);
      featureOK = false;
    }
  }

  return(featureOK);
}

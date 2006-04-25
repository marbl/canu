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


void
decodeMatch(splitToWords &W,
            u32bit &iid1, u32bit &pos1, u32bit &len1, u32bit &fwd1,
            u32bit &iid2, u32bit &pos2, u32bit &len2, u32bit &fwd2) {

  iid1 = decodeAtacName(W[4], 0L);
  pos1 = strtou32bit(W[5], 0L);
  len1 = strtou32bit(W[6], 0L);
  fwd1 = (W[7][0] == '-') ? 0 : 1;
  iid2 = decodeAtacName(W[8], 0L);
  pos2 = strtou32bit(W[9], 0L);
  len2 = strtou32bit(W[10], 0L);
  fwd2 = (W[11][0] == '-') ? 0 : 1;
}


void
decodeFeature(splitToWords &W,
              u32bit &iid, u32bit &pos, u32bit &len) {

  iid = decodeAtacName(W[4], 0L);
  pos = strtou32bit(W[5], 0L);
  len = strtou32bit(W[6], 0L);
}


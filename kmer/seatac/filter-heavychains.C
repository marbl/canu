// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Copyright (c) 2005 The J. Craig Venter Institute
// Author: Clark Mobarry
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

#include "util++.H"
#include "heavychains.H"


extern "C" {
  void    *construct(char *options);
  void     destruct(void *handle);
  void     addHit(void *handle,
                  char    orientation,
                  u32bit  id1,
                  u32bit  pos1,
                  u32bit  len1,
                  u32bit  id2,
                  u32bit  pos2,
                  u32bit  len2,
                  u32bit  filled);
  void     filter(void *handle);
  u64bit   output(void *handle, FILE *file, u64bit matchid);

  void    *constructStats(char *options);
  void     destructStats(void *handle);
  void     addStats(void  *handle, void *sp);
  void     showStats(void *handle, FILE *file);
}


void*
construct(char *options) {
  int    beVerbose    = 0;
  char  *assemblyId1  = "UNK";
  char  *assemblyId2  = "UNK";
  double minScore     = 100.0;   // Default minimum of bp filled in a good run.
  int    maxJump      = 100000;  // Default maximum intra-run jump allowed in a good run.

  //  Parse the options to find the parameters
  //
  splitToWords  W(options);

  u32bit arg = 0;
  while (arg < W.numWords()) {
    if        (strcmp(W.getWord(arg), "-v") == 0) {
      beVerbose++;
    } else if (strcmp(W.getWord(arg), "-s") == 0) {
      minScore = atof(W.getWord(++arg));
    } else if (strcmp(W.getWord(arg), "-j") == 0) {
      maxJump = atoi(W.getWord(++arg));
    } else if (strcmp(W.getWord(arg), "-1") == 0) {
      assemblyId1 = W.getWord(++arg);
    } else if (strcmp(W.getWord(arg), "-2") == 0) {
      assemblyId2 = W.getWord(++arg);
    }

    arg++;
  }

  return((void *)(new StrandPair(beVerbose, assemblyId1, assemblyId2, maxJump, minScore)));
}

void
destruct(void *handle) {
  delete (StrandPair *)handle;
}

void
addHit(void   *handle,
       char    orientation,
       u32bit  id1,
       u32bit  pos1,
       u32bit  len1,
       u32bit  id2,
       u32bit  pos2,
       u32bit  len2,
       u32bit  filled) {
  ((StrandPair *)handle)->addHit(orientation, id1, pos1, len1, id2, pos2, len2, filled);
}

void
filter(void *handle) {
  ((StrandPair *)handle)->process();
}


u64bit
output(void *handle, FILE *file, u64bit matchid) {
  return(((StrandPair *)handle)->print(file, matchid));
}






void*
constructStats(char *options) {
  int    beVerbose    = 0;
  char  *assemblyId1  = 0L;
  char  *assemblyId2  = 0L;
  double minScore     = 100.0;   // Default minimum of bp filled in a good run.
  int    maxJump      = 100000;  // Default maximum intra-run jump allowed in a good run.

  //  Parse the options to find the parameters
  //
  splitToWords  W(options);

  u32bit arg = 0;
  while (arg < W.numWords()) {
    if        (strcmp(W.getWord(arg), "-v") == 0) {
      beVerbose++;
    } else if (strcmp(W.getWord(arg), "-s") == 0) {
      minScore = atof(W.getWord(++arg));
    } else if (strcmp(W.getWord(arg), "-j") == 0) {
      maxJump = atoi(W.getWord(++arg));
    } else if (strcmp(W.getWord(arg), "-1") == 0) {
      assemblyId1 = W.getWord(++arg);
    } else if (strcmp(W.getWord(arg), "-2") == 0) {
      assemblyId2 = W.getWord(++arg);
    }

    arg++;
  }

  return((void *)(new TheStats(beVerbose, assemblyId1, assemblyId2, maxJump, minScore)));
}

void
destructStats(void *handle) {
  delete (TheStats *)handle;
}

void
addStats(void   *handle, void *sp) {
  ((TheStats *)handle)->add((StrandPair *)sp);
}

void
showStats(void *handle, FILE *file) {
  ((TheStats *)handle)->show(file);
}

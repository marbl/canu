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

#ifndef ATAC_FEATURE_H
#define ATAC_FEATURE_H

#include <string>
#include <map>
using namespace std;

#include "bio++.H"

//  A barebones feature.
//
//  F type featureuid parentuid LABEL:IID POS
//
//  The idea is to mark some region as containing some feature.  You
//  can use the featureiid to map to a chunk of non-atac data, e.g., a
//  strucure containing information about the feature.
//
//  An EST feature could be represented as:
//    F est  est4     . B35:3 423551 10421
//    F exon exon7 est4 B35:3 423551 346
//    F exon exon8 est4 B35:3 425931 146
//    F exon exon9 est4 B35:3 433426 546
//
//  There is nothing there that immediately links these atac features
//  to their respective est/exon data structures.  This is a Good
//  Thing (tm) because usually we don't have any form of ID with those
//  data structures, so we'd be using an offset or something, which
//  would change if the features are filtered.
//
class atacFeature {
public:
  atacFeature() {};
  atacFeature(char *line);
  atacFeature(char *fuid,
              char *puid,
              u32bit fiid,
              char *t,
              u32bit i, u32bit p, u32bit l);

  void   decode(char *line);

  bool   sanity(seqFile *A, char *inLine);

  char    featureuid[16];
  char    parentuid[16];
  u32bit  featureiid;
  char    type[4];
  u32bit  iid, pos, len;

  void  print(FILE *f, char const *label) const {
    fprintf(f, "F %s %s %s %s:"u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
            type,
            featureuid, parentuid,
            label, iid, pos, len);
  };
};

#endif  //  ATAC_FEATURE_H

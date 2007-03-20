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

#ifndef ATAC_MATCH_H
#define ATAC_MATCH_H

#include <string>
#include <map>
using namespace std;

#include "bio++.H"

class atacMatch {
public:
  atacMatch() {
    matchuid[0] = 0;
    parentuid[0] = 0;
    matchiid = 0;
    type[0] = 0;
    iid1 = pos1 = len1 = fwd1 = 0;
    iid2 = pos2 = len2 = fwd2 = 0;
  };
  atacMatch(char *line);
  atacMatch(char *muid,
            char *puid,
            u32bit miid,
            char *t,
            u32bit i1, u32bit p1, u32bit l1, u32bit f1,
            u32bit i2, u32bit p2, u32bit l2, u32bit f2);

  void  decode(char *line);

  //  Sanity check the match record -- make sure it's within the
  //  sequence itself.
  //
  bool   sanity(FastABase *A, FastABase *B, char *inLine);


  char    matchuid[16];     //  external id
  char    parentuid[16];    //  external parent id
  u32bit  matchiid;         //  internal id, usually pointing to an entry in atacMatchList
  char    type[4];          //  right now, only need one byte, but we keep things aligned
  u32bit  iid1, pos1, len1, fwd1;
  u32bit  iid2, pos2, len2, fwd2;

  void  print(FILE *f, char const *label1, char const *label2) const {
    fprintf(f, "M %s %s %s %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
            type,
            matchuid, parentuid,
            label1, iid1, pos1, len1, fwd1 ? 1 : -1,
            label2, iid2, pos2, len2, fwd2 ? 1 : -1);
  };
};

#endif  //  ATAC_MATCH_H

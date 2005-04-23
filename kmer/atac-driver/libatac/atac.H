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

#ifndef ATAC_COMMON_H
#define ATAC_COMMON_H

#include "bio++.H"

void   readHeader(char *inLine, FILE *in, char *file1, char *file2, FILE *out);
void   decodeMatch(splitToWords &W,
                   u32bit &iid1, u32bit &pos1, u32bit &len1, u32bit &ori1,
                   u32bit &iid2, u32bit &pos2, u32bit &len2, u32bit &ori2);


class match_t {
public:
  u32bit  matchiid;  //  internal id, usually pointing to an entry in matchList
  u32bit  matchuid;  //  external id
  u32bit  iid1, pos1, len1, ori1;
  u32bit  iid2, pos2, len2, ori2;
};


class matchList {
public:
  matchList(char *filename);
  ~matchList() {
    delete _seq1;
    delete _seq2;
    delete [] _matches;
  };

  u32bit         numMatches(void) {
    return(_matchesLen);
  };
  match_t       *getMatch(u32bit i) {
    return(_matches + i);
  };

  match_t       *operator[](u32bit i) {
    return(_matches + i);
  };

  //  Sort by the A assembly or B assembly location
  void           sort1(u32bit first=0, u32bit len=0);
  void           sort2(u32bit first=0, u32bit len=0);

  char          _file1[1024];
  char          _file2[1024];

  FastAWrapper *_seq1;
  FastAWrapper *_seq2;

  u32bit        _matchesLen;
  u32bit        _matchesMax;
  match_t      *_matches;
};


#endif  //  ATAC_COMMON_H

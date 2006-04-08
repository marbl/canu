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

#include <string>
#include <map>
using namespace std;

#include "bio++.H"

void   readHeader(char *inLine,
                  FILE *in,
                  char *file1=0L,
                  char *file2=0L,
                  FILE *out=0L,
                  map<string,string> *params=0L);

void   decodeMatchNames(splitToWords &W,
                        char *name1,
                        char *nam2);
void   decodeMatch(splitToWords &W,
                   u32bit &iid1, u32bit &pos1, u32bit &len1, u32bit &fwd1,
                   u32bit &iid2, u32bit &pos2, u32bit &len2, u32bit &fwd2);


class atacMatch {
public:
  char    matchuid[16];     //  external id
  char    parentuid[16];    //  external parent id
  u32bit  matchiid;         //  internal id, usually pointing to an entry in atacMatchList
  u32bit  pad;              //  nothing, just keep things aligned
  u32bit  iid1, pos1, len1, fwd1;
  u32bit  iid2, pos2, len2, fwd2;

  void  print(FILE *f, char *label1, char *label2) {
    fprintf(f, "M u %s %s %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
            matchuid, parentuid,
            label1, iid1, pos1, len1, fwd1 ? 1 : -1,
            label2, iid2, pos2, len2, fwd2 ? 1 : -1);
  }
};


class atacMatchList {
public:
  atacMatchList(char *filename, char matchOrRun = 'm', bool saveLine=false);
  ~atacMatchList() {
    delete _seq1;
    delete _seq2;
    delete [] _matches;
  };

  u32bit         numMatches(void) {
    return(_matchesLen);
  };
  atacMatch     *getMatch(u32bit i) {
    return(_matches + i);
  };

  atacMatch     *operator[](u32bit i) {
    return(_matches + i);
  };

  //  Sort by the A assembly or B assembly location
  void           sort1(u32bit first=0, u32bit len=0);
  void           sort2(u32bit first=0, u32bit len=0);
  void           sortDiagonal(u32bit first=0, u32bit len=0);
  void           sortMatchUID(u32bit first=0, u32bit len=0);
  void           sortParentUID(u32bit first=0, u32bit len=0);

  char          _file1[1024];
  char          _file2[1024];

  //  The name / label of each of the sequences
  char          _name1[256];
  char          _name2[256];

  FastAWrapper *_seq1;
  FastAWrapper *_seq2;

  u32bit        _matchesLen;
  u32bit        _matchesMax;
  atacMatch    *_matches;
};


#endif  //  ATAC_COMMON_H

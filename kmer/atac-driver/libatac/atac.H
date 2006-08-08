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

u32bit decodeAtacName(char *name, char *label);
void   decodeMatch(splitToWords &W,
                   u32bit &iid1, u32bit &pos1, u32bit &len1, u32bit &fwd1,
                   u32bit &iid2, u32bit &pos2, u32bit &len2, u32bit &fwd2);
void   decodeFeature(splitToWords &W,
                     u32bit &iid, u32bit &pos, u32bit &len);


class atacMatch {
public:
  char    matchuid[16];     //  external id
  char    parentuid[16];    //  external parent id
  u32bit  matchiid;         //  internal id, usually pointing to an entry in atacMatchList
  char    type[4];          //  right now, only need one byte, but we keep things aligned
  u32bit  iid1, pos1, len1, fwd1;
  u32bit  iid2, pos2, len2, fwd2;

  void  print(FILE *f, char *label1, char *label2) {
    fprintf(f, "M %s %s %s %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
            type,
            matchuid, parentuid,
            label1, iid1, pos1, len1, fwd1 ? 1 : -1,
            label2, iid2, pos2, len2, fwd2 ? 1 : -1);
  };
};


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
  char    featureuid[16];
  char    parentuid[16];
  u32bit  featureiid;
  char    type[4];
  u32bit  iid, pos, len;

  void  print(FILE *f, char *label) {
    fprintf(f, "F %s %s %s %s:"u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
            type,
            featureuid, parentuid,
            label, iid, pos, len);
  };
};


class atacMatchList {
public:
  atacMatchList(char *filename, char matchOrRun = 'm', FILE *headerOut=0L);
  ~atacMatchList() {
    delete _seq1;
    delete _seq2;
    delete [] _matches;
  };

  u32bit         numMatches(void)      { return(_matchesLen);  };
  u32bit         numberOfMatches(void) { return(_matchesLen);  };
  atacMatch     *getMatch(u32bit i)    { return(_matches + i); };
  atacMatch     *operator[](u32bit i)  { return(_matches + i); };

  char          *labelA(void) { return(_label1); };
  char          *labelB(void) { return(_label2); };

  char          *assemblyFileA(void) { return(_file1); };
  char          *assemblyFileB(void) { return(_file2); };

  FastAWrapper  *fastaA(void) { return(_seq1); };
  FastAWrapper  *fastaB(void) { return(_seq2); };

  //  Sort by the A assembly or B assembly location
  void           sort1(u32bit first=0, u32bit len=0);
  void           sort2(u32bit first=0, u32bit len=0);
  void           sortDiagonal(u32bit first=0, u32bit len=0);
  void           sortMatchUID(u32bit first=0, u32bit len=0);
  void           sortParentUID(u32bit first=0, u32bit len=0);

  //  Merge the r match into the l match.  The l match gets the result.
  void           mergeMatches(atacMatch *l, atacMatch *r, u32bit mergeuid);

private:
  char          _label1[256];  //  The label of each of the sequences
  char          _label2[256];

  char          _file1[1024];  //  The name of our genome files
  char          _file2[1024];

  FastAWrapper *_seq1;
  FastAWrapper *_seq2;

  u32bit        _matchesLen;
  u32bit        _matchesMax;
  atacMatch    *_matches;
};




class atacFeatureList {
public:
  atacFeatureList(char *filename, char *featureType=0L, bool saveLine=false);
  ~atacFeatureList() {
    delete _seq;
    delete [] _features;
  };

  //u32bit         numFeatures(void)      { return(_featuresLen);  };
  u32bit         numberOfFeatures(void) { return(_featuresLen);  };
  atacFeature   *getFeature(u32bit i)   { return(_features + i); };
  atacFeature   *operator[](u32bit i)   { return(_features + i); };

  char          *label(void)        { return(_label); };
  char          *assemblyFile(void) { return(_file); };

  //  Sort by the A assembly or B assembly location
  void           sort(u32bit first=0, u32bit len=0);
  void           sortFeatureUID(u32bit first=0, u32bit len=0);
  void           sortParentUID(u32bit first=0, u32bit len=0);

private:
  char          _label[256];   //  The label of the sequence
  char          _file[1024];  //  The name of our genome files
  FastAWrapper *_seq;

  u32bit        _featuresLen;
  u32bit        _featuresMax;
  atacFeature  *_features;
};


#endif  //  ATAC_COMMON_H

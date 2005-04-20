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

//  Looks for 1's surrounded by U's
void
findIsolatedUnique(annoList *AL, u32bit ALlen) {
  bool   only1 = true;
  u32bit sumA = 0, tA=0;
  u32bit sumB = 0, tB=0;

  for (u32bit i=1; i<ALlen; i++) {
    if (AL[i].type == 'U') {
      if (only1 && (!tA || !tB)) {
        sumA += tA;
        sumB += tB;
      }
      tA = 0;
      tB = 0;
      only1 = true;
    } else if (AL[i].type != '1') {
      only1 = false;
    } else {
      tA += AL[i].len2a;
      tB += AL[i].len2b;
    }
  }

  fprintf(stderr, "isolated Unique: map1: "u32bitFMT" map2: "u32bitFMT"\n", sumA, sumB);
}



u32bit encodeType(u32bit type) {
  u32bit  t = 9;
  if      (type == 'U') t = 0;
  else if (type == '1') t = 1;
  else if (type == 'Y') t = 2;
  else if (type == 'N') t = 3;
  else if (type == '?') t = 4;
  else if (type == '!') t = 5;
  if (t == 9)
    fprintf(stderr, "got invalid type; "u32bitFMT" -- %c\n", type, (char)type), exit(1);
  return(t);
}



//  Looks for 1's at the end of a Y, that have the same match id
void
findExtended(annoList *AL, u32bit ALlen) {

  class stat_s {
  public:
    u32bit  count;
    u32bit  len;

    stat_s() {
      count = 0;
      len   = 0;
    };

    stat_s &operator+=(u32bit x) {
      count++;
      len += x;
      return(*this);
    };
    stat_s &operator+=(stat_s &x) {
      count += x.count;
      len   += x.len;
      return(*this);
    };

    void print(char *msg) {
      fprintf(stderr, "%s:  "u32bitFMT" len:"u32bitFMT"\n", msg, count, len);
    };
  };

  stat_s  count1[6][6][2][2][2];
  stat_s  count2[6][6][2][2][2];

  //  Look forward for the next event:
  //    non-1
  //    the iid changed
  //
  //  End is the match after the 1, beg is the match before the 1.
  //
  for (u32bit i=0; i<ALlen; i++) {
    if (AL[i].type == '1') {
      u32bit  beg = i-1;
      u32bit  gap = i;
      u32bit  end = i;
      u32bit  len = 0;

      if (AL[gap].iid2a) {
        while ((AL[gap].iid2a == AL[end].iid2a) && (AL[end].type == '1')) {
          len += AL[end].len1;
          end++;
        }
      } else {
        while ((AL[gap].iid2b == AL[end].iid2b) && (AL[end].type == '1')) {
          len += AL[end].len1;
          end++;
        }
      }

      if (beg == gap)
        fprintf(stderr, "beg == gap?\n"), exit(1);
      if (gap == end)
        fprintf(stderr, "end == gap?\n"), exit(1);

      u32bit tbeg = encodeType(AL[beg].type);
      u32bit tend = encodeType(AL[end].type);

      bool  mbeg, mend, moth;
      if (AL[gap].iid2a) {
        mbeg = (AL[beg].iid2a == AL[gap].iid2a);
        mend = (AL[end].iid2a == AL[gap].iid2a);
        moth = (AL[beg].iid2b == AL[end].iid2b);
        count1[tbeg][tend][mbeg][mend][moth] += len;
      } else {
        mbeg = (AL[beg].iid2b == AL[gap].iid2b);
        mend = (AL[end].iid2b == AL[gap].iid2b);
        moth = (AL[beg].iid2a == AL[end].iid2a);
        count2[tbeg][tend][mbeg][mend][moth] += len;
      }
    }
  }


  char  label[6] = {'U', '1', 'Y', 'N', '?', '!'};

#define TYPE_U 0
#define TYPE_1 1
#define TYPE_Y 2
#define TYPE_N 3
#define TYPE_Q 4
#define TYPE_E 5

  //  If the other iid is the same, then these are all interesting cases
  //  Maybe not.
#if 0
  stat_s   oth_cnst_1;
  stat_s   oth_cnst_2;
  stat_s   oth_diff_1;
  stat_s   oth_diff_2;

  for (u32bit i=0; i<6; i++)
    for (u32bit j=0; j<6; j++)
      for (u32bit k=0; k<2; k++)
        for (u32bit l=0; l<2; l++) {
          oth_cnst_1 += count1[i][j][k][l][1];
          oth_cnst_2 += count2[i][j][k][l][1];
          oth_diff_1 += count1[i][j][k][l][0];
          oth_diff_2 += count2[i][j][k][l][0];
        }
  fprintf(stderr, "count1 other iid constant:  "u32bitFMT" len:"u32bitFMT"\n", oth_cnst_1.count, oth_cnst_1.len);
  fprintf(stderr, "count2 other iid constant:  "u32bitFMT" len:"u32bitFMT"\n", oth_cnst_2.count, oth_cnst_2.len);
  fprintf(stderr, "count1 other iid different: "u32bitFMT" len:"u32bitFMT"\n", oth_diff_1.count, oth_diff_1.len);
  fprintf(stderr, "count2 other iid different: "u32bitFMT" len:"u32bitFMT"\n", oth_diff_2.count, oth_diff_2.len);
#endif

#if 0
  for (u32bit tbeg=0; tbeg<6; tbeg++)
    for (u32bit tend=0; tend<6; tend++)
      for (u32bit mbeg=0; mbeg<2; mbeg++)
        for (u32bit mend=0; mend<2; mend++)
          for (u32bit moth=0; moth<2; moth++) {
          }
#endif


  //  Look for things with an extra piece in the middle -- mbeg and mend
  //
  stat_s   extraMid1[7];
  stat_s   extraMid2[7];

  for (u32bit tbeg=0; tbeg<6; tbeg++)
    for (u32bit tend=0; tend<6; tend++)
      for (u32bit mbeg=0; mbeg<2; mbeg++)
        for (u32bit mend=0; mend<2; mend++)
          for (u32bit moth=0; moth<2; moth++) {
            if (mbeg && mend) {
              extraMid1[6] += count1[tbeg][tend][mbeg][mend][moth];
              extraMid2[6] += count2[tbeg][tend][mbeg][mend][moth];

              if (tbeg == tend) {
                extraMid1[tbeg] += count1[tbeg][tend][mbeg][mend][moth];
                extraMid2[tbeg] += count2[tbeg][tend][mbeg][mend][moth];
              }
            }
          }

  extraMid1[6].print("extra middle 1");
  extraMid2[6].print("extra middle 2");
  for (u32bit i=0; i<6; i++) {
    char l[64];
    sprintf(l, "extra middle 1 %c", label[i]);
    extraMid1[i].print(l);
    sprintf(l, "extra middle 2 %c", label[i]);
    extraMid2[i].print(l);
  }

  //  Look for true extensions
  //
  stat_s  extY1;
  stat_s  extN1;
  stat_s  extY2;
  stat_s  extN2;

  for (u32bit tbeg=0; tbeg<6; tbeg++)
    for (u32bit tend=0; tend<6; tend++)
      for (u32bit mbeg=0; mbeg<2; mbeg++)
        for (u32bit mend=0; mend<2; mend++)
          for (u32bit moth=0; moth<2; moth++) {
            if (mbeg && !mend) {
              if (tbeg == TYPE_Y) {
                extY1 += count1[tbeg][tend][mbeg][mend][moth];
              } else {
                extN1 += count1[tbeg][tend][mbeg][mend][moth];
              }
            }
            if (!mbeg && mend) {
              if (tend == TYPE_Y) {
                extY2 += count2[tbeg][tend][mbeg][mend][moth];
              } else {
                extN2 += count2[tbeg][tend][mbeg][mend][moth];
              }
            }
          }
  extY1.print("extension Y 1");
  extN1.print("extension N 1");
  extY2.print("extension Y 2");
  extN2.print("extension N 2");


  fprintf(stderr, "----------------------------------------\n");
  fprintf(stderr, "BEGIN IID SAME, END IID SAME, OTHER IID SAME\n");
  fprintf(stderr, "----------------------------------------\n");
  fprintf(stderr, "END IID\n");
  fprintf(stderr, "----------------------------------------\n");
  fprintf(stderr, "OTHER IID\n");
}

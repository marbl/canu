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
#include <time.h>

#include "bio++.H"
#include "util++.H"
#include "atac-common.H"

//  Compute some simple statistics on a set of matches

int
main(int argc, char **argv) {
  char    inLine[1024];

  char    file1[1024];
  char    file2[1024];

  intervalList  intervalA;
  intervalList  intervalB;

  u32bit        skippedCount = 0;
  u64bit        skippedLengthA = 0;
  u64bit        skippedLengthB = 0;

  //  Read the preamble, look for our data sources.  This leaves us with
  //  the first match in the inLine, and fills in file1 and file2.
  //
  readHeader(inLine, stdin, file1, file2, 0L);

  fprintf(stderr, "Opening '%s' for sequence one.\n", file1);
  fprintf(stderr, "Opening '%s' for sequence two.\n", file2);

  //  Open some FastAWrappers for each of the files -- we use these
  //  only to get the length of the sequence.
  //
  FastAWrapper  *C1 = new FastAWrapper(file1);
  FastAWrapper  *C2 = new FastAWrapper(file2);

  C1->openIndex();
  C2->openIndex();

  //  For the coverage to work correctly, we need to either have one
  //  intervalList per input sequence, or build a table of the chained
  //  sequence positions.
  //
  u64bit  *offset1 = new u64bit [C1->getNumberOfSequences()];
  u64bit  *offset2 = new u64bit [C2->getNumberOfSequences()];

  offset1[0] = 1000000;
  for (u32bit i=1; i<C1->getNumberOfSequences(); i++)
    offset1[i] = offset1[i-1] + C1->sequenceLength(i-1) + 1;

  offset2[0] = 1000000;
  for (u32bit i=1; i<C2->getNumberOfSequences(); i++)
    offset2[i] = offset2[i-1] + C2->sequenceLength(i-1) + 1;

  while (!feof(stdin)) {
    if (inLine[0] == 'M') {
      splitToWords  S(inLine);

      if ((S[1][0] == 'u') || (S[1][0] == 'x')) {
        //if ((S[1][0] == 'r')) {
        u32bit  iid1=0, pos1=0, len1=0, ori1=0;
        u32bit  iid2=0, pos2=0, len2=0, ori2=0;
        decodeMatch(S, iid1, pos1, len1, ori1, iid2, pos2, len2, ori2);

        if ((pos1 + len1) > C1->sequenceLength(iid1)) {
          chomp(inLine);
          fprintf(stderr, "Too Long in 1: "u32bitFMT" %s\n", C1->sequenceLength(iid1), inLine);
        }

        if ((pos2 + len2) > C2->sequenceLength(iid2)) {
          chomp(inLine);
          fprintf(stderr, "Too Long in 2: "u32bitFMT" %s\n", C2->sequenceLength(iid2), inLine);
        }


        if ((iid1 >= C1->getNumberOfSequences()) || (iid2 >= C2->getNumberOfSequences())) {
          //  Hmmm.  Skip it.
          skippedCount++;
          skippedLengthA += len1;
          skippedLengthB += len2;
        } else {
          intervalA.add(offset1[iid1] + (u64bit)pos1, (u64bit)len1);
          intervalB.add(offset2[iid2] + (u64bit)pos2, (u64bit)len2);
        }
      }
    }

    fgets(inLine, 1024, stdin);
  }

  fprintf(stderr, "skipped "u32bitFMT" matches with length "u64bitFMT" and "u64bitFMT"\n",
          skippedCount, skippedLengthA, skippedLengthB);

  fprintf(stderr, "intervalLength A "u64bitFMT" B "u64bitFMT"\n",
          (u64bit)intervalA.sumOfLengths(),
          (u64bit)intervalB.sumOfLengths());

  intervalA.merge();
  intervalB.merge();

  fprintf(stderr, "coveredLength  A "u64bitFMT" B "u64bitFMT"\n",
          (u64bit)intervalA.sumOfLengths(),
          (u64bit)intervalB.sumOfLengths());

  return(0);
}

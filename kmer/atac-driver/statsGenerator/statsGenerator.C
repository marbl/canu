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

  char   *file1 = 0L;
  char   *file2 = 0L;

  u64bit  coveredLengthA = 0;
  u64bit  coveredLengthB = 0;


  //  Read the preamble, look for our data sources.  This leaves us with
  //  the first match in the inLine, and fills in file1 and file2.
  //
  readHeader(inLine, stdin, file1, file2, 0L);

  //  Open some FastAWrappers for each of the files -- we use these
  //  only to get the length of the sequence.
  //
  FastAWrapper  *C1 = new FastAWrapper(file1);
  FastAWrapper  *C2 = new FastAWrapper(file2);

  C1->openIndex();
  C2->openIndex();


  while (!feof(stdin)) {
    if (inLine[0] == 'M') {
      splitToWords  S(inLine);

      if ((S[1][0] == 'u') || (S[1][0] == 'x')) {
        u32bit  iid1=0, pos1=0, len1=0, ori1=0;
        u32bit  iid2=0, pos2=0, len2=0, ori2=0;
        decodeMatch(S, iid1, pos1, len1, ori1, iid2, pos2, len2, ori2);

        coveredLengthA += len1;
        coveredLengthB += len2;
      }
    }

    fgets(inLine, 1024, stdin);
  }

  fprintf(stderr, "coveredLengthA:    "u64bitFMT"\n", coveredLengthA);
  fprintf(stderr, "coveredLengthB:    "u64bitFMT"\n", coveredLengthB);
}


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
#include "atac-common.H"

//#define ANNOTATE

//  Generates a histogram of the exact match block sizes
//  Counts to global number of mismatches
//  Annotates each match with the number of mismatches
//  Checks for identities outside matches




void
updateExactBlockHistogram(u32bit *blockHistogram, u32bit blockMatches) {

  if (blockMatches > 8 * 1024 * 1024)
    blockHistogram[0]++;
  else
    blockHistogram[blockMatches]++;
}



int
main(int argc, char *argv[]) {

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-h") == 0) {
      //  Generate a histogram of exact-match lengths
    } else if (strcmp(argv[arg], "-a") == 0) {
      //  Annotate each match with the percent error, compute
      //  the global percent error.
    } else if (strcmp(argv[arg], "-e") == 0) {
      //  Generate a histogram of the percent error in each match
    } else if (strcmp(argv[arg], "-c") == 0) {
      //  Check the edges of each match to ensure there isn't a match
    } else {
      fprintf(stderr, "usage: %s [-h exact-match-histogram] [-a] [-e error-histogram] [-c]\n", argv[0]);
      fprintf(stderr, "  -h:     histogram of the length of the exact match blocks\n");
      fprintf(stderr, "  -a:     annotate each match with the percent error, write to stdout\n");
      fprintf(stderr, "  -e:     histogram of the error rate of each match\n");
      fprintf(stderr, "  -c:     check that the next base on each side is a mismatch\n");
      exit(1);
    }
  }

  char  inLine[1024] = {0};
  char  file1[1024]  = {0};
  char  file2[1024]  = {0};

  u32bit   globalSequence = 0;
  u32bit   globalMismatches = 0;
  u32bit   blockMatches = 0;
  u32bit  *blockHistogram = new u32bit [8 * 1024 * 1024];

  for (u32bit x=0; x<8*1024*1024; x++)
    blockHistogram[x] = 0;

  //  Read the preamble, look for our data sources.  This leaves us with
  //  the first match in the inLine, and fills in file1 and file2.
  //
#ifdef ANNOTATE
  readHeader(inLine, stdin, file1, file2, stdout);
#else
  readHeader(inLine, stdin, file1, file2, 0L);
#endif

  //  Open some FastACache's for each of the files
  //
  FastACache  *C1 = new FastACache(file1, 1, false, false);
  FastACache  *C2 = new FastACache(file2, 1, false, false);

  //  While not end-of-file, read the matches.  If we get something
  //  that isn't a match, just silently emit it.
  //
  //  Need to exclude RUNS
  //
  while (!feof(stdin)) {
    if ((inLine[0] == 'M') && (inLine[2] != 'r')) {
      splitToWords  W(inLine);

      //  Parse out the sequence iid from the atac iid
      //
      u32bit  iid1=0, pos1=0, len1=0, ori1=0;
      u32bit  iid2=0, pos2=0, len2=0, ori2=0;

      decodeMatch(W, iid1, pos1, len1, ori1, iid2, pos2, len2, ori2);

      //  Grab those sequences from the cache
      //
      FastASequenceInCore  *S1 = C1->getSequence(iid1);
      FastASequenceInCore  *S2 = C2->getSequence(iid2);

      FastAAccessor A1(S1, false);
      FastAAccessor A2(S2, (ori1 != ori2));

      A1.setReverseComplementRange(pos1, len1);
      A2.setReverseComplementRange(pos2, len2);

      u32bit extraMatches = 0;
      u32bit localMismatches = 0;

      //  Check for matches on either side of the region.
      //  (but only if there is stuff on the other side)

      if (A1.setPosition(pos1 - 1) && A2.setPosition(pos2 - 1))
        if (validSymbol[*A1] &&
            validSymbol[*A2] &&
            IUPACidentity[*A1][*A2])
          extraMatches++;

      if (A1.setPosition(pos1+len1) && A2.setPosition(pos2+len2))
        if (validSymbol[*A1] &&
            validSymbol[*A2] &&
            IUPACidentity[*A1][*A2])
          extraMatches++;


      //  WARN if we found extra identities
      //
      if (extraMatches > 0) {
        A1.setPosition(pos1);
        A2.setPosition(pos2);

        fprintf(stderr, "WARNING: found extra matches in %s\n", inLine);

        for (u32bit ii=0; ii<len1; ii++, ++A1)
          fprintf(stdout, "%c", *A1);
        fprintf(stdout, "\n");

        for (u32bit ii=0; ii<len1; ii++, ++A2)
          fprintf(stdout, "%c", *A2);
        fprintf(stdout, "\n");
      }


      A1.setPosition(pos1);
      A2.setPosition(pos2);

      for (u32bit ii=0; ii<len1; ii++, ++A1, ++A2) {

        //
        //  do stuff here
        //

        ////////////////////////////////////////
        //
        //  Count global matches / mismatches
        //
        globalSequence++;
        if (validSymbol[*A1] &&
            validSymbol[*A2] &&
            !IUPACidentity[*A1][*A2]) {
          globalMismatches++;
          localMismatches++;
        }

        ////////////////////////////////////////
        //
        //  Histogram of exact match block lengths
        //
        if (validSymbol[*A1] &&
            validSymbol[*A2] &&
            IUPACidentity[*A1][*A2]) {
          blockMatches++;
        } else {
          updateExactBlockHistogram(blockHistogram, blockMatches);
          blockMatches = 0;
        }
      }

      ////////////////////////////////////////
      //
      //  Finish off stuff
      //
      updateExactBlockHistogram(blockHistogram, blockMatches);
      blockMatches = 0;

#ifdef ANNOTATE
      chomp(inLine);
      fputs(inLine, stdout);
      fprintf(stdout, " > /extramatches="u32bitFMT" /mismatches="u32bitFMT"\n",
              extraMatches, localMismatches);
    } else {
      fputs(inLine, stdout);
#endif
    }

    fgets(inLine, 1024, stdin);
  }


  ////////////////////////////////////////
  //
  //  Report stuff
  //
  fprintf(stderr, "globalSequence   = "u32bitFMT"\n", globalSequence);
  fprintf(stderr, "globalMismatches = "u32bitFMT"\n", globalMismatches);

  FILE *O = fopen("MismatchCounter.block.histogram.out", "w");
  for (u32bit i=0; i<8 * 1024 * 1024; i++)
    fprintf(O, u32bitFMT" "u32bitFMT"\n", i, blockHistogram[i]);
  fclose(O);

  return(0);
}

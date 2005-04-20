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
#include "atac.H"

//  Compute some simple statistics on a set of matches

void mappedMultiply1(matchList &matches, char *prefix);
void mappedMultiply2(matchList &matches, char *prefix);

int
u32bitcompare(const void *a, const void *b) {
  const u32bit A = *((const u32bit *)a);
  const u32bit B = *((const u32bit *)b);
  if (A < B) return(-1);
  if (A > B) return(1);
  return(0);
}

int
main(int argc, char **argv) {

  if (argc != 3) {
    fprintf(stderr, "usage: %s <file.atac> <outprefix>\n", argv[0]);
    exit(1);
  }

  //  matchList also computes the length and coverage of the matches.
  //
  matchList     matches(argv[1]);
  FILE         *out;
  char         *prefix = argv[2];
  char          filename[1024];

  //  Generate an Nx plot, and a match length histogram
  //  We block the histogram at 1kb, 0 is for things too big.
  //
  u32bit    histogramMax    = 1000000;
  u32bit    histogramBlock  = 100;
  u32bit   *lengthHistogram = new u32bit [histogramMax / histogramBlock];
  u32bit   *n50             = new u32bit [matches.numMatches()];

  for (u32bit i=0; i<histogramMax / histogramBlock; i++)
    lengthHistogram[i] = 0;

  for (u32bit i=0; i<matches.numMatches(); i++) {
    match_t  *m = matches.getMatch(i);

    //  Update the length histogram
    //
    if (m->len1 > histogramMax) {
      lengthHistogram[0]++;
    } else {
      lengthHistogram[m->len1 / histogramBlock]++;
    }

    //  Save the length for the n50 plot
    n50[i] = m->len1;
  }

  //  Dump the histogram
  sprintf(filename, "%s.matchlengthhistogram", prefix);
  out = fopen(filename, "w");
  for (u32bit i=0; i<histogramMax / histogramBlock; i++)
    fprintf(out, u32bitFMT" "u32bitFMT"\n", i, lengthHistogram[i]);
  fclose(out);

  //  Compute the total length of the sequence
  u32bit totalLength = 0;
  for (u32bit i=0; i<matches._seq1->getNumberOfSequences(); i++)
    totalLength += matches._seq1->sequenceLength(i);

  //  Sort the n50 list of lengths
  qsort(n50, matches.numMatches(), sizeof(u32bit), u32bitcompare);

  //  It's slow and obvious and, yes, there is a better way.  Dump the
  //  Nx plot as it's being generated.

  sprintf(filename, "%s.Nx", prefix);
  out = fopen(filename, "w");

  for (u32bit n=1; n<100; n++) {
    u32bit  limit = totalLength / 100 * n;
    u32bit  iter  = 0;
    u32bit  sum   = 0;

    while ((sum < limit) && (iter < matches.numMatches())) {
      sum += n50[iter++];
    }

    fprintf(out, u32bitFMT" "u32bitFMT" "u32bitFMT"\n", n, n50[iter-1], iter);
  }

  fclose(out);

  //  plot [0:200][0:50000] "gunk.matchlengthhistogram" using 2 with lines, "gunk.Nx" with lines, "funk.matchlengthhistogram" using 2 with lines, "funk.Nx" with lines


  ////////////////////////////////////////
  //
  //  Decide how many sequences in A are mapped to more than one thing
  //  in B.  For each sequence in A, we compute the percent of it
  //  mapped, the largest and smallest percent mapped to a single
  //  sequence, and the number of sequences it mapped to.
  //
  //  A:4 num-mapped-to percent-mapped largest-mapped smallest-mapped
  //
  mappedMultiply1(matches, prefix);
  mappedMultiply2(matches, prefix);

  return(0);
}



//  XXX:  Sadly, this function is needed twice, one for each direction.

void
mappedMultiply1(matchList &matches, char *prefix) {
  u32bit   beg = 0;
  u32bit   end = 0;
  char     filename[1024];
  FILE    *out;

  sprintf(filename, "%s.multiple1", prefix);
  out = fopen(filename, "w");

  matches.sort1();

  while (beg < matches.numMatches()) {

    //  Figure out how many matches we have for this sequence
    end = beg;
    while ((end < matches.numMatches()) && (matches[beg]->iid1 == matches[end]->iid1))
      end++;

    //  Sort that range by the other index
    matches.sort2(beg, end-beg);

    //  Count the number of sequences this sequence maps to
    u32bit numTargets = 1;

    for (u32bit i=beg+1; i<end; i++) {
      if (matches[i-1]->iid2 != matches[i]->iid2)
        numTargets++;
    }

    //  Build an interval list (of the regions in us) for each sequence we map to
    intervalList  AL;
    intervalList  IL[numTargets];
    u32bit target = 0;

    for (u32bit i=beg; i<end; i++) {
      AL.add(matches[i]->pos1, matches[i]->len1);
      IL[target].add(matches[i]->pos1, matches[i]->len1);

      if ((i+1 < end) &&
          (matches[i]->iid2 != matches[i+1]->iid2))
        target++;
    }

    //  Squash the lists, compute the min and max coverage

    double minc = 100.0;
    double maxc =   0.0;

    for (u32bit i=0; i<numTargets; i++) {
      IL[i].merge();

      double c = 100.0 * IL[i].sumOfLengths() / matches._seq1->sequenceLength(matches[beg]->iid1);

      if (minc > c)  minc = c;
      if (maxc < c)  maxc = c;
    }

    //  report

    AL.merge();

    fprintf(out, u32bitFMT" "u32bitFMT" "u32bitFMT" %8.4f %8.4f %8.4f\n",
            matches[beg]->iid1, numTargets,
            end - beg,
            minc,
            maxc,
            100.0 * AL.sumOfLengths() / matches._seq1->sequenceLength(matches[beg]->iid1) );

    beg = end+1;
  }

  fclose(out);
}




void
mappedMultiply2(matchList &matches, char *prefix) {
  u32bit   beg = 0;
  u32bit   end = 0;
  char     filename[1024];
  FILE    *out;

  sprintf(filename, "%s.multiple2", prefix);
  out = fopen(filename, "w");

  matches.sort2();

  while (beg < matches.numMatches()) {

    //  Figure out how many matches we have for this sequence
    end = beg;
    while ((end < matches.numMatches()) && (matches[beg]->iid2 == matches[end]->iid2))
      end++;

    //  Sort that range by the other index
    matches.sort1(beg, end-beg);

    //  Count the number of sequences this sequence maps to
    u32bit numTargets = 1;

    for (u32bit i=beg+1; i<end; i++) {
      if (matches[i-1]->iid1 != matches[i]->iid1)
        numTargets++;
    }

    //  Build an interval list (of the regions in us) for each sequence we map to
    intervalList  AL;
    intervalList  IL[numTargets];
    u32bit target = 0;

    for (u32bit i=beg; i<end; i++) {
      AL.add(matches[i]->pos2, matches[i]->len2);
      IL[target].add(matches[i]->pos2, matches[i]->len2);

      if ((i+1 < end) &&
          (matches[i]->iid1 != matches[i+1]->iid1))
        target++;
    }

    //  Squash the lists, compute the min and max coverage

    double minc = 100.0;
    double maxc =   0.0;

    for (u32bit i=0; i<numTargets; i++) {
      IL[i].merge();

      double c = 100.0 * IL[i].sumOfLengths() / matches._seq2->sequenceLength(matches[beg]->iid2);

      if (minc > c)  minc = c;
      if (maxc < c)  maxc = c;
    }

    //  report

    AL.merge();

    fprintf(out, u32bitFMT" "u32bitFMT" "u32bitFMT" %8.4f %8.4f %8.4f\n",
            matches[beg]->iid2, numTargets,
            end - beg,
            minc,
            maxc,
            100.0 * AL.sumOfLengths() / matches._seq2->sequenceLength(matches[beg]->iid2) );

    beg = end+1;
  }

  fclose(out);
}

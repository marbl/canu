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


int
main(int argc, char **argv) {

  if (argc != 4) {
    fprintf(stderr, "usage: %s <matches-1> <matches-2> <out-prefix>\n", argv[0]);
    exit(1);
  }
  atacFile       *AF1 = new atacFile(argv[1]);
  atacFile       *AF2 = new atacFile(argv[2]);

  atacMatchList  *M1 = AF1->matches();
  atacMatchList  *M2 = AF2->matches();

  char           *OP = argv[3];


  //  We want to annotate the two assembies with:
  //    a) mapped by both, the same
  //    b) mapped by both, differently
  //    c) mapped by the first, unmapped by the second
  //    d) mapped by the second, unmapped by the first
  //    e) unmapped by both
  //
  //  If unmapped, we could further annotate with the reason it was
  //  unmapped -- not found, or found multiple times.
  //
  //  Our annotation datastructure is a tree of spans.  Each span is a
  //  sequence, and an interval on that sequence.  We assume that the
  //  tree contains the spans for the whole sequence, that is, that we
  //  never need to increase a span, just split.
  //
  spanTree    *S1 = new spanTree();
  spanTree    *S2 = new spanTree();

  //  Initialize the tree of spans by inserting a single span for each
  //  sequence in the file.
  //
  for (u32bit i=0; i<AF1->fastaA()->getNumberOfSequences(); i++)
    S1->addNewSpan(i, AF1->fastaA()->getSequenceLength(i));
  for (u32bit i=0; i<AF1->fastaB()->getNumberOfSequences(); i++)
    S2->addNewSpan(i, AF1->fastaB()->getSequenceLength(i));

  //  Add every match to the spanTrees.

  for (u32bit i=0; i<M1->numberOfMatches(); i++) {
    S1->addMatch(M1->getMatch(i), 0, 0);
    S2->addMatch(M1->getMatch(i), 1, 0);
  }
  for (u32bit i=0; i<M2->numberOfMatches(); i++) {
    S1->addMatch(M2->getMatch(i), 0, 1);
    S2->addMatch(M2->getMatch(i), 1, 1);
  }

  //  Dump each spanTree: For each span, we need to check that
  //    it has matches?
  //    only one match, or only matches from one mapping?
  //    matches from both mappings?  need to check that
  //     the span in the other tree also has the same matches
  //
  //  Doesn't handle weird stuff like this span (on sequence 1)
  //  mapping onto seq2 correctly, but the span in seq2 having an
  //  extra match to somewhere else in seq1.
  //
  //  we want to find the single span in the other spanTree that
  //  corresponds to this span.  once we do that, we can verify that
  //  all the matches are the same.
  //
  //  because we are gapless matches, we can, for each match,
  //  compute the exact location this span should occur on the other
  //  sequence.  then, do a lookup() to get that span, or just
  //  verify that everybody is the same location.

  char  outname[1024];
  FILE *outfile;

  overlapStats  statsA;
  u32bit     ALmax = (u32bit)dict_count(S1->_tree);
  u32bit     ALlen = 0;
  annoList  *AL    = new annoList [ ALmax ];

  sprintf(outname, "%s.map1annotation", OP);
  errno = 0;
  outfile = fopen(outname, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outname, strerror(errno));
  process1(outfile, S1, M1, M2, statsA, AL, ALlen, ALmax);
  fclose(outfile);

  overlapStats  statsB;
  u32bit     BLmax = (u32bit)dict_count(S1->_tree);
  u32bit     BLlen = 0;
  annoList  *BL    = new annoList [ ALmax ];

  sprintf(outname, "%s.map2annotation", OP);
  errno = 0;
  outfile = fopen(outname, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outname, strerror(errno));
  process2(outfile, S2, M1, M2, statsB, BL, BLlen, BLmax);
  fclose(outfile);

  fprintf(stderr, "unmapped:           A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", statsA.unmapped.getSum(),     statsB.unmapped.getSum());
  fprintf(stderr, "unique mapping 1:   A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", statsA.map1unique.getSum(),   statsB.map1unique.getSum());
  fprintf(stderr, "unique mapping 2:   A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", statsA.map2unique.getSum(),   statsB.map2unique.getSum());
  fprintf(stderr, "different:          A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", statsA.different.getSum(),    statsB.different.getSum());
  fprintf(stderr, "wild diff:          A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", statsA.wilddiff.getSum(),     statsB.wilddiff.getSum());
  fprintf(stderr, "same:               A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", statsA.same.getSum(),         statsB.same.getSum());
  fprintf(stderr, "inconsistent:       A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", statsA.inconsistent.getSum(), statsB.inconsistent.getSum());

  //  Dump the histograms for each of the labelings
  //
  sprintf(outname, "%s.asm1histogram", OP);
  statsA.writeHistogram(outname);
  sprintf(outname, "%s.asm2histogram", OP);
  statsB.writeHistogram(outname);

  //  Draw some pretty pictures
  //
  sprintf(outname, "%s.histogram.gnuplot", OP);
  errno = 0;
  outfile = fopen(outname, "w");
  if (errno)
    fprintf(stderr, "failed to open '%s': %s\n", outname, strerror(errno)), exit(1);
  fprintf(outfile, "set terminal postscript color\n");
  fprintf(outfile, "set output \"%s.unmapped.histogram.ps\"\n", OP);
  fprintf(outfile, "set ylabel \"number of regions\"\n");
  fprintf(outfile, "set xlabel \"length of region\"\n");
  fprintf(outfile, "plot [0:10000][0:400] \\\n");
  fprintf(outfile, "          \"%s.asm1histogram.unmapped\" using 2 title \"assembly 1 unmapped\" with lines, \\\n", OP);
  fprintf(outfile, "          \"%s.asm2histogram.unmapped\" using 2 title \"assembly 2 unmapped\" with lines\n", OP);
  fprintf(outfile, "set output \"%s.same.histogram.ps\"\n", OP);
  fprintf(outfile, "plot [0:20000][0:2000] \\\n");
  fprintf(outfile, "          \"%s.asm1histogram.same\" using 2 title \"assembly 1 same\" with lines, \\\n", OP);
  fprintf(outfile, "          \"%s.asm2histogram.same\" using 2 title \"assembly 2 same\" with lines\n", OP);
  fprintf(outfile, "set output \"%s.histogram.ps\"\n", OP);
  fprintf(outfile, "plot [0:2000][0:100] \\\n");
  fprintf(outfile, "          \"%s.asm1histogram.different\" using 2 title \"assembly 1 different\" with lines, \\\n", OP);
  fprintf(outfile, "          \"%s.asm2histogram.different\" using 2 title \"assembly 2 different\" with lines, \\\n", OP);
  fprintf(outfile, "          \"%s.asm1histogram.wilddiff\" using 2 title \"assembly 1 wildly diff\" with lines, \\\n", OP);
  fprintf(outfile, "          \"%s.asm2histogram.wilddiff\" using 2 title \"assembly 2 wildly diff\" with lines\n", OP);
  fprintf(outfile, "set output \"%s.unique.histogram.ps\"\n", OP);
  fprintf(outfile, "plot [0:2000][0:100] \\\n");
  fprintf(outfile, "          \"%s.asm1histogram.map1unique\" using 2 title \"map 1, assembly 1 unique\" with lines, \\\n", OP);
  fprintf(outfile, "          \"%s.asm1histogram.map2unique\" using 2 title \"map 2, assembly 1 unique\" with lines, \\\n", OP);
  fprintf(outfile, "          \"%s.asm2histogram.map1unique\" using 2 title \"map 1, assembly 2 unique\" with lines, \\\n", OP);
  fprintf(outfile, "          \"%s.asm2histogram.map2unique\" using 2 title \"map 2, assembly 2 unique\" with lines\n", OP);
  fclose(outfile);

  sprintf(outname, "gnuplot < %s.histogram.gnuplot", OP);
  if (system(outname))
    fprintf(stderr, "Failed to '%s'\n", outname);

#if 0
  findIsolatedUnique(AL, ALlen);
  findExtended(AL, ALlen);
#endif

  //  Deleting the spanTrees takes a long time, so we don't bother with any cleanup.
  return(0);
}

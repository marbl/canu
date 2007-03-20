// This file is part of A2Amapper.
// Copyright (c) 2006 J. Craig Venter Institute
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
#include "atac.H"


//  Attempts to align unmapped regions.
//
//  For each unmapped region, we extract the corresponding sequences,
//  sim4db them together, parse the output to make atac-format
//  matches, but of a lower confidence.
//
//  IDX1   ------------------------------------------
//         ||||||||                        ||||||||||
//  IDX2   -A------------\          /--------------B-
//                        \        /
//
//  The nasty case is that IDX1 could be doubly mapped, once by A and
//  once by B.  So we also need to label those regions that are mapped
//  multiple times as an even lower confidence.
//
//  We probably should bias the alignment towards the anchored edge,
//  implying I should use something other than sim4db here.
//
//  We end up with three confidence classes:
//  1)  mapped by atac itself, 1-to-1 matches
//  2)  mapped by sim4db above, with no conflict, between anchors
//  3)  same as 2, but conflicting
//
//  Why sim4db?  It's splicing model might introduce some noise on the
//  ends (which we'll clean up), but more importantly, the splicing
//  allows us to skip over large blocks of whatever (rearrangement,
//  tandem repeat, etc).  And it's also in my source tree and I know
//  how to use it.


this is unfinished crap



//  The below is the main from writing unmatched regions




int
main(int argc, char *argv[]) {
  FILE         *Aoutput = 0L;
  FILE         *Boutput = 0L;
  char         *matchesFile = 0L;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-a") == 0) {
      errno = 0;
      Aoutput = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s': %s\n", argv[arg], strerror(errno)), exit(1);
    } else if (strcmp(argv[arg], "-b") == 0) {
      errno = 0;
      Boutput = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s': %s\n", argv[arg], strerror(errno)), exit(1);
    } else if (strcmp(argv[arg], "-m") == 0) {
      matchesFile = argv[++arg];
    } else {
      fprintf(stderr, "usage: %s -a Aunmatched.fasta -b B.unmatched.fasta < matches\n", argv[0]);
      exit(1);
    }
    arg++;
  }

  if ((Aoutput == 0L) || (Boutput == 0L) || (matchesFile == 0L)) {
    fprintf(stderr, "usage: %s -a Aunmatched.fasta -b B.unmatched.fasta < matches\n", argv[0]);
    exit(1);
  }

  atacMatchList  ML1(matchesFile, 'm', false);
  atacMatchList  ML2(matchesFile, 'm', false);

  ML1.sort1();  //  Sorted by first index
  ML2.sort2();  //  Sorted by second index

  FastABase            *W1 = ML1._seq1;
  FastABase            *W2 = ML1._seq2;


  //  For every match,
  //    find the match before and the match after, on both axes
  //    




  //  Extract unmapped in sequence 1

  ML.sort1();
  W = ML._seq1;
  W->find(ML[0]->iid1);
  S = W->getSequence();
  for (u32bit i=1; i<ML.numMatches(); i++) {
    atacMatch *l = ML[i-1];
    atacMatch *r = ML[i];

    if (l->iid1 != r->iid1)
      continue;

    if (l->iid1 != S->getIID()) {
      delete S;
      W->find(l->iid1);
      S = W->getSequence();
    }

    //  Extract from (l->pos1 + l->len1) to (r->pos1), if it's longer than 20bp
    //
    if (l->pos1 + l->len1 + 20 < r->pos1)
      writeGaplessSequence(Aoutput,
                           S,
                           l->pos1 + l->len1,
                           r->pos1);
  }

  //  Extract unmapped in sequence 2

  ML.sort2();
  W = ML._seq2;
  W->find(ML[0]->iid2);
  S = W->getSequence();
  for (u32bit i=1; i<ML.numMatches(); i++) {
    atacMatch *l = ML[i-1];
    atacMatch *r = ML[i];

    if (l->iid2 != r->iid2)
      continue;

    if (l->iid2 != S->getIID()) {
      delete S;
      W->find(l->iid2);
      S = W->getSequence();
    }

    //  Extract from (l->pos2 + l->len2) to (r->pos2), if it's longer than 20bp
    //
    if (l->pos2 + l->len2 + 20 < r->pos2)
      writeGaplessSequence(Boutput,
                           S,
                           l->pos2 + l->len2,
                           r->pos2);
  }


  fclose(Aoutput);
  fclose(Boutput);

  return(0);
}

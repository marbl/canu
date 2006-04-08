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

//  Reads a set of atac matches, computes the percent identity of the
//  regions, and warns if any identites are low.

int
main(int argc, char *argv[]) {
  char         *matchesFile   = 0L;
  double        identityLimit = 0.9;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      matchesFile = argv[++arg];
    } else if (strcmp(argv[arg], "-i") == 0) {
      identityLimit = atof(argv[++arg]);
      if (identityLimit > 1.0)
        identityLimit /= 100;
    } else {
      fprintf(stderr, "usage: %s -m matches\n", argv[0]);
      exit(1);
    }
    arg++;
  }

  if (matchesFile == 0L)
      fprintf(stderr, "usage: %s -m matches\n", argv[0]), exit(1);

  atacMatchList  ML(matchesFile, 'm', false);
  FastACache     Acache(ML._file1, 0, true, false);
  FastACache     Bcache(ML._file2, 0, true, false);

  for (u32bit i=0; i<ML.numMatches(); i++) {
    atacMatch            *m = ML.getMatch(i);

    u32bit identities = 0;

    char  *a = Acache.getSequence(m->iid1)->sequence() + m->pos1;
    char  *b = Bcache.getSequence(m->iid2)->sequence() + m->pos2;

    if (m->fwd2) {
      for (u32bit p=0; p<m->len1; p++) {
        if (toUpper[a[p]] == toUpper[b[p]])
          identities++;
      }
    } else {
      for (u32bit p=0, q=m->len2-1; p<m->len1; p++, q--) {
        if (toUpper[a[p]] == toUpper[complementSymbol[b[q]]])
          identities++;
      }
    }

    if ((double)identities / m->len1 < identityLimit) {
      fprintf(stderr, "match "u32bitFMT" is only %6.2f%% identity:  ",
              i, 100.0 * identities / m->len1);
      m->print(stderr, ML._name1, ML._name2);
      if (m->len1 < 200) {
        char   tmp[1000];

        strncpy(tmp, a, m->len1);
        tmp[m->len1] = 0;
        fprintf(stderr, "  %s\n", tmp);

        strncpy(tmp, b, m->len1);
        tmp[m->len1] = 0;
        fprintf(stderr, "  %s\n", tmp);
      }
    }
  }

  return(0);
}

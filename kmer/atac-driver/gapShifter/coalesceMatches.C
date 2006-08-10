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

//  Reads a set of matches, coalesces those on the same diagonal.
//  Does not preserve runs.
//
//  No args, reads stdin, writes stdout.

int
main(int argc, char *argv[]) {
  atacMatchList  ML("-", 'm', stdout);
  atacMatch     *l = 0L;
  atacMatch     *r = 0L;

  ML.sortDiagonal();

  for (u32bit i=1; i<ML.numMatches(); i++) {
    l = ML[i-1];
    r = ML[i];

    if ((l->iid1 == r->iid1) &&
        (l->iid2 == r->iid2) &&
        (l->fwd1 == r->fwd1) &&
        (l->fwd2 == r->fwd2) &&
        (l->pos1 + l->len1 == r->pos1) &&
        (l->pos2 + l->len2 == r->pos2) &&
        (strcmp(l->type, r->type) == 0) &&
        (strcmp(l->parentuid, r->parentuid) == 0)) {

      fprintf(stderr, "MERGE:\n");
      l->print(stderr, ML.labelA(), ML.labelB());
      r->print(stderr, ML.labelA(), ML.labelB());

      l->len1 += r->len1;
      l->len2 += r->len2;

      l->print(stderr, ML.labelA(), ML.labelB());
    } else {
      l->print(stdout, ML.labelA(), ML.labelB());
      l = 0L;
    }
  }

  if (l)
    l->print(stdout, ML.labelA(), ML.labelB());

  return(0);
}

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
//  Does not preserve runs or comments or much else in the input
//  file(s).

int
main(int argc, char *argv[]) {

  //  No args, reads stdin, writes stdout.

  atacMatchList  ML("-", 'm');

  //  Sort by the ID's; we'll coalesce 
  ML.sortDiagonal();

  for (u32bit i=0; i<ML.numMatches(); i++) {
    atacMatch *m1 = ML[i];

    fprintf(stdout, "M u %s %s B35LC:"u32bitFMT" "u32bitFMT" "u32bitFMT" 1  HUREF2:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
            m1->matchuid, m1->parentuid,
            m1->iid1, m1->pos1, m1->len1,
            m1->iid2, m1->pos2, m1->len2, m1->fwd2 ? 1 : -1);
  }

  return(0);
}

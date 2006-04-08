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

#ifndef OVERLAP_H
#define OVERLAP_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "bio++.H"
#include "util++.H"
#include "atac.H"

int sortMatches1(const void *a, const void *b);
int sortMatches2(const void *a, const void *b);
int spanCompare(const void *a, const void *b);

#include "overlap-span.H"
#include "overlap-matchTree.H"
#include "overlap-spanTree.H"
#include "overlap-annoList.H"
#include "overlap-stats.H"

void
process1(FILE             *outfile,
         spanTree         *S,
         atacMatchList    *M1,
         atacMatchList    *M2,
         overlapStats     &stats,
         annoList         *AL,
         u32bit           &ALlen,
         u32bit           &ALmax);

void
process2(FILE             *outfile,
         spanTree         *S,
         atacMatchList    *M1,
         atacMatchList    *M2,
         overlapStats     &stats,
         annoList         *AL,
         u32bit           &ALlen,
         u32bit           &ALmax);

void
printAnno(FILE *F, annoList *AL, u32bit &ALlen,
          char label,
          u32bit axis,
          span_t *span,
          u32bit match1=u32bitZERO, atacMatch *m1=0L,
          u32bit match2=u32bitZERO, atacMatch *m2=0L);

void
findIsolatedUnique(annoList *AL, u32bit ALlen);

void
findExtended(annoList *AL, u32bit ALlen);


#endif  //  OVERLAP_H

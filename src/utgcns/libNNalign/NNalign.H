
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

#ifndef ALIGNERS_H
#define ALIGNERS_H

static const char *rcsid_ALIGNERS_H = "$Id$";

#include "AS_global.H"

//
//  BPW's brute-force dynamic programming
//

typedef struct {
  int32   matches;
  int32   alignLen;
  int32   begI, begJ;
  int32   endI, endJ;
  int32   lenA, lenB;

  double  pIdentity;
  double  pCoverageA;
  double  pCoverageB;
} alignLinker_s;


void
alignLinker(char           *alignA,
            char           *alignB,
            char           *stringA,
            char           *stringB,
            alignLinker_s  *a,
            int             endToEnd,
            int             allowNs,
            int             ahang, int bhang);



#endif

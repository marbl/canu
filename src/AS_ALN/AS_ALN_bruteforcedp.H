
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007-2008, J. Craig Venter Institute
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

#ifndef AS_ALN_BRUTEFORCEDP
#define AS_ALN_BRUTEFORCEDP

static const char *rcsid_AS_ALN_BRUTEFORCEDP = "$Id$";


typedef struct {
  unsigned int  score  : 30;
  unsigned int  action : 2;
} dpCell;


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
            dpCell        (*M)[AS_READ_MAX_NORMAL_LEN + 1],
            alignLinker_s  *a,
            int             endToEnd,
            int             allowNs,
            int             ahang, int bhang);

#endif  //  AS_ALN_BRUTEFORCEDP

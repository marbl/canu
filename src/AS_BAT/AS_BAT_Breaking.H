
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

#ifndef INCLUDE_AS_BAT_BREAKING
#define INCLUDE_AS_BAT_BREAKING

static const char *rcsid_INCLUDE_AS_BAT_BREAKING = "$Id$";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"

//  What to do with contained fragments at the split point?
//
//  ejectContains -- eject all contains from the unitig, client is expected to
//  place them again.
//
//  keepContains -- contains that are associated with the fragEnd are moved to
//  the new unitig.  There are five cases:
//
//  A #-----------------*
//  1   ------                           - contained in A
//  2        -------                     - contained in A, dovetail to B
//  B          -----------------------
//  3             ------                 - contained in A, contained in B
//  4                 -------            - dovetail to A, contained in B
//  5                       -------      - contained in B
//
//  If the break point is at #:
//    keepContains == true  -- 1, 2, 3 will move with A.
//    keepContains == false -- 1, 2    will move with A.
//
//  If the break point is at *:
//    keepContains == true  -- 1, 2, 3 will move with A.
//    keepContains == false -- 1, 2    will move with A.
//

class breakPoint {
public:
  breakPoint() {
    fragEnd             = FragmentEnd();
    ejectContains       = false;
    keepContains        = false;
    break5p             = false;
    break3p             = false;
    rememberLastTig     = false;
    searchDiscontinuous = false;
  };

  breakPoint(uint32 isectFrg, bool isect3p, bool eject, bool keep) {
    fragEnd             = FragmentEnd(isectFrg, isect3p);
    ejectContains       = eject;
    keepContains        = keep;
    break5p             = false;
    break3p             = false;
    rememberLastTig     = false;
    searchDiscontinuous = false;
  };

  bool operator<(breakPoint const that) const {
    return(fragEnd < that.fragEnd);
  };

public:
  FragmentEnd fragEnd;

  bool        ejectContains;
  bool        keepContains;

  bool        break5p;
  bool        break3p;
  bool        rememberLastTig;
  bool        searchDiscontinuous;
};



bool
breakUnitigAt(UnitigVector &unitigs, Unitig *tig, vector<breakPoint> &breaks, bool doDelete);

#endif //  INCLUDE_AS_BAT_BREAKING

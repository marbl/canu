
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

static const char *rcsid_INCLUDE_AS_BAT_BREAKING = "$Id: AS_BAT_Breaking.H,v 1.1 2010-12-06 08:03:48 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"


struct UnitigBreakPoint {
  FragmentEnd fragEnd;          // frag id and which end to break on
  SeqInterval fragPos;          // coordinates in unitig (used to get fwd/rev)

  int         fragsBefore;      // number of frags before the break fragment
  int         fragsAfter;       // number of frags after the break fragment

  FragmentEnd inEnd;            // frag id of the incoming unitig

  int         inSize;           // the size of the incoming unitig
  int         inFrags;          // the number of fragments in incoming unitig

  UnitigBreakPoint() {
    fragEnd      = FragmentEnd();
    fragPos.bgn  = 0;
    fragPos.end  = 0;

    fragsBefore  = 0;
    fragsAfter   = 0;

    inEnd        = FragmentEnd();
    inSize       = 0;
    inFrags      = 0;
  };

  UnitigBreakPoint(uint32 isectFrg, bool isect3p, SeqInterval isectPos,
                   uint32 frgsBefore, uint32 frgsAfter,
                   uint32 invaderFrg, bool invader3p,
                   uint32 invaderTigLen,
                   uint32 invaderTigFrags) {
    fragEnd      = FragmentEnd(isectFrg, isect3p);
    fragPos      = isectPos;

    fragsBefore  = frgsBefore;
    fragsAfter   = frgsAfter;

    inEnd        = FragmentEnd(invaderFrg, invader3p);
    inSize       = invaderTigLen;
    inFrags      = invaderTigFrags;
  };
};


typedef std::vector<UnitigBreakPoint> UnitigBreakPoints;

void
filterBreakPoints(Unitig *, UnitigBreakPoints &);

UnitigVector *
breakUnitigAt(Unitig *, UnitigBreakPoints &);

#endif //  INCLUDE_AS_BAT_BREAKING

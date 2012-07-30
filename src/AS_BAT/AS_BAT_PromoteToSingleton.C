
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2011, The Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id: AS_BAT_PromoteToSingleton.C,v 1.3 2012-07-30 01:21:01 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"

//  If we are not reconstructing repeats, promote all the unplaced fragments to new unitigs.
//  Oodles of possibilities here; promote everything to a singleton unitig, promote only
//  the non-contained, then place contains, then promote what is left over, etc.

void
promoteToSingleton(UnitigVector &unitigs, bool enablePromoteToSingleton) {

  for (uint32 fi=1; fi<=FI->numFragments(); fi++) {
    if (Unitig::fragIn(fi) != 0)
      //  Placed already
      continue;

    if (FI->fragmentLength(fi) == 0)
      //  Deleted.
      continue;

    if (enablePromoteToSingleton == false) {
      writeLog("promoteToSingleton()--  Repeat fragment "F_U32" removed from assembly.\n", fi);
      FI->markAsIgnore(fi);
      continue;
    }

    Unitig *utg = unitigs.newUnitig(false);
    ufNode  frag;
  
    frag.ident             = fi;
    frag.contained         = 0;
    frag.parent            = 0;
    frag.ahang             = 0;
    frag.bhang             = 0;
    frag.position.bgn      = 0;
    frag.position.end      = FI->fragmentLength(fi);
    frag.containment_depth = 0;

    utg->addFrag(frag, 0, false);
  }
}

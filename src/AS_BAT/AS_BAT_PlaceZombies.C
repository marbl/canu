

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

static const char *rcsid = "$Id: AS_BAT_PlaceZombies.C,v 1.2 2010-12-06 08:03:48 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_PlaceZombies.H"


//  Zombies are either caused by a bug (hope not) or by a contained fragment being contained in something that
//  is eventually contained in the original fragment.  A contains B, B contains C....and C contains A.

void
placeZombies(UnitigVector &unitigs) {

  fprintf(logFile, "==> SEARCHING FOR ZOMBIES\n");

  uint32 *inUnitig   = new uint32 [FI->numFragments()+1];
  int     numZombies = 0;

  //  Mark fragments as dead.
  //
  for (uint32 i=0; i<FI->numFragments()+1; i++)
    inUnitig[i] = noUnitig;

  //  ZZZzzzzaapaapppp!  IT'S ALIVE!
  //
  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    for (uint32 fi=0; fi<utg->ufpath.size(); fi++) {
      ufNode  *frag = &utg->ufpath[fi];

      inUnitig[frag->ident] = utg->id();
    }
  }

  //  Anything still dead?
  //
  for (uint32 i=0; i<FI->numFragments()+1; i++) {
    if (FI->fragmentLength(i) == 0)
      //  Deleted fragment
      continue;

    if (inUnitig[i] != noUnitig)
      //  Valid fragment in a unitig, other errors caught in checkUnitigMembership()
      continue;

    Unitig      *utg = new Unitig(false);
    ufNode       frg;

    frg.ident             = i;
    frg.contained         = 0;
    frg.containment_depth = 0;

    frg.position.bgn      = 0;
    frg.position.end      = FI->fragmentLength(i);

    utg->addFrag(frg, 0, false);
    unitigs.push_back(utg);

    fprintf(logFile, "placeZombies()-- unitig %d created from zombie fragment %d\n",
            utg->id(), i);
    numZombies++;
  }

  fprintf(logFile, "RESURRECTED %d ZOMBIE FRAGMENT%s.\n", numZombies, (numZombies != 1) ? "s" : "");

  delete [] inUnitig;
}

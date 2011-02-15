

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

static const char *rcsid = "$Id: AS_BAT_PlaceZombies.C,v 1.4 2011-02-15 08:10:11 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_PlaceZombies.H"


//  Zombies are caused by a fragment being contained in a fragment that is eventually contained in
//  the original fragment -- circular containmnents.
//
//  Here we detect Zombies, and reset their best container to something that is already placed.

void
placeZombies(UnitigVector &unitigs, double erate, double elimit) {

  fprintf(logFile, "==> SEARCHING FOR ZOMBIES\n");

  uint32 *inUnitig   = new uint32 [FI->numFragments()+1];
  int     numZombies = 0;

  //  Mark fragments as dead, then unmark them if they are in a real living unitig.

  for (uint32 i=0; i<FI->numFragments()+1; i++)
    inUnitig[i] = noUnitig;

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    for (uint32 fi=0; fi<utg->ufpath.size(); fi++)
      inUnitig[utg->ufpath[fi].ident] = utg->id();
  }

  //  For anything not in a living unitig, reload the overlaps and find a new container.
  //  (NOT IMPLEMENTED - for now we just move these to new singleton unitigs).

  for (uint32 i=0; i<FI->numFragments()+1; i++) {
    if (FI->fragmentLength(i) == 0)
      //  Deleted fragment
      continue;

    if (inUnitig[i] != noUnitig)
      //  Valid fragment in a unitig
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



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

static const char *rcsid = "$Id: AS_BOG_PlaceZombies.cc,v 1.2 2010-09-28 09:17:54 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"





//  This is a huge hack to get around a bug somewhere before us.  We
//  seem to not be placing all fragments.  So, run through all
//  fragments, and for anything not placed, toss it into a new
//  unitig.  Let scaffolder figure out where to put it.
//
//  Notice the similarity between this code and checkUnitigMembership().
//
void
UnitigGraph::placeZombies(void) {

  fprintf(logFile, "==> SEARCHING FOR ZOMBIES\n");

  uint32 *inUnitig   = new uint32 [_fi->numFragments()+1];
  int     numZombies = 0;

  //  Mark fragments as dead.
  //
  for (uint32 i=0; i<_fi->numFragments()+1; i++)
    inUnitig[i] = noUnitig;

  //  ZZZzzzzaapaapppp!  IT'S ALIVE!
  //
  for (uint32 ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];

    if (utg == NULL)
      continue;

    for (uint32 fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode  *frag = &(*utg->dovetail_path_ptr)[fi];

      inUnitig[frag->ident] = utg->id();
    }
  }

  //  Anything still dead?
  //
  for (uint32 i=0; i<_fi->numFragments()+1; i++) {
    if (_fi->fragmentLength(i) == 0)
      //  Deleted fragment
      continue;

    if (inUnitig[i] != noUnitig)
      //  Valid fragment in a unitig, other errors caught in checkUnitigMembership()
      continue;

    //  Ha!  Gotcha!  You're now a resurrected brain eating
    //  zomibie?  Some day we'll figure out how to put you in
    //  properly.  For now, enjoy the ride.

    Unitig      *utg = new Unitig(false);
    DoveTailNode frg;

    frg.ident             = i;
    frg.contained         = 0;
    frg.containment_depth = 0;

    frg.position.bgn      = 0;
    frg.position.end      = _fi->fragmentLength(i);

    utg->addFrag(frg, 0, false);
    unitigs->push_back(utg);

    fprintf(logFile, "placeZombies()-- unitig %d created from zombie fragment %d\n",
            utg->id(), i);
    numZombies++;
  }

  fprintf(logFile, "RESURRECTED %d ZOMBIE FRAGMENT%s.\n", numZombies, (numZombies != 1) ? "s" : "");

  delete [] inUnitig;
}



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

static const char *rcsid = "$Id: AS_BOG_PlaceZombies.cc,v 1.1 2010-09-23 09:34:50 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"

#undef max




//  This is a huge hack to get around a bug somewhere before us.  We
//  seem to not be placing all fragments.  So, run through all
//  fragments, and for anything not placed, toss it into a new
//  unitig.  Let scaffolder figure out where to put it.
//
//  Notice the similarity between this code and checkUnitigMembership().
//
void
UnitigGraph::placeZombies(void) {

  fprintf(stderr, "==> SEARCHING FOR ZOMBIES\n");

  uint32 *inUnitig   = new uint32 [_fi->numFragments()+1];
  int     numZombies = 0;

  //  Sorry, poor fella that has more than 987,654,321 unitigs.  I
  //  can't imagine the rest of the pipeline would run though.
  //
  //  Mark fragments as dead.
  //
  for (uint32 i=0; i<_fi->numFragments()+1; i++)
    inUnitig[i] = noUnitig;

  //  ZZZzzzzaapaapppp!  IT'S ALIVE!
  //
  for (int  ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];

    if (utg)
      for (DoveTailIter it=utg->dovetail_path_ptr->begin(); it != utg->dovetail_path_ptr->end(); it++)
        inUnitig[it->ident] = utg->id();
  }

  //  Anything still dead?
  //
  for (uint32 i=0; i<_fi->numFragments()+1; i++) {
    if (_fi->fragmentLength(i) > 0) {
      if (inUnitig[i] == 0) {
        //  We'll catch this error inna second in checkUnitigMembership().
      } else if (inUnitig[i] != noUnitig) {
        //  We'll count this inna second there too.
      } else {
        //  Ha!  Gotcha!  You're now a resurrected brain eating
        //  zomibie?  Some day we'll figure out how to put you in
        //  properly.  For now, enjoy the ride.

        Unitig *utg = new Unitig(false);

        DoveTailNode frag;

        frag.ident             = i;
        frag.contained         = 0;
        frag.containment_depth = 0;

        frag.position.bgn      = 0;
        frag.position.end      = _fi->fragmentLength(i);

        utg->addFrag(frag, 0, false);
        unitigs->push_back(utg);

        numZombies++;
      }
    }
  }

  if (numZombies > 0)
    fprintf(stderr, "RESURRECTED %d ZOMBIE FRAGMENTS.\n", numZombies);

  delete [] inUnitig;
}

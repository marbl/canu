
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

static const char *rcsid = "$Id: AS_BAT_ReconstructRepeats.C,v 1.2 2012-02-15 03:41:08 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_ChunkGraph.H"

#include "AS_BAT_PopulateUnitig.H"
#include "AS_BAT_PlaceContains.H"

#include "AS_BAT_Unitig.H"


//  estimate read error rate from best overlaps (per library?)
//  use that error rate below when rebuilding repeats

void
reconstructRepeats(UnitigVector &unitigs,
                   double        erateGraph,
                   double        elimitGraph) {

  //  Similar to mate extension, we build a set<> of all the unplaced fragments, then
  //  construct a new BOG and CG from which we construct unitigs.

  BestOverlapGraph  *OGsave = OG;
  ChunkGraph        *CGsave = CG;

  set<AS_IID>        unplaced;

  for (uint32 fi=1; fi<=FI->numFragments(); fi++)
    if (Unitig::fragIn(fi) == 0)
      unplaced.insert(fi);

  OG = new BestOverlapGraph(erateGraph / 2.0, elimitGraph / 2.0, &unplaced);
  CG = new ChunkGraph(&unplaced);

  fprintf(logFile, "==> BUILDING REPEAT UNITIGS from %d fragments.\n", unplaced.size());

  for (uint32 fi=CG->nextFragByChunkLength(); fi>0; fi=CG->nextFragByChunkLength())
    populateUnitig(unitigs, fi);

  fprintf(logFile, "==> BUILDING REPEAT UNITIGS catching missed fragments.\n");

  for (uint32 fi=1; fi <= FI->numFragments(); fi++)
    populateUnitig(unitigs, fi);

  fprintf(logFile, "==> BUILDING REPEAT UNITIGS placing contained fragments.\n");

  placeContainsUsingBestOverlaps(unitigs);

  delete OG;
  delete CG;

  OG = OGsave;
  CG = CGsave;
}

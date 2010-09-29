
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

static const char *rcsid = "$Id: AS_BOG_Instrumentation.cc,v 1.3 2010-09-29 22:04:29 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"



void
UnitigGraph::checkUnitigMembership(void) {
  int nutg = 0;
  int nfrg = 0;

  fprintf(logFile, "checkUnitigMembership()--  numfrags=%d\n", _fi->numFragments());

  uint32 *inUnitig = new uint32 [_fi->numFragments()+1];
  uint32  logSize[20] = {0};

  for (uint32 i=0; i<_fi->numFragments()+1; i++)
    inUnitig[i] = noUnitig;

  for (uint32 ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];
    uint32   len = 0;

    if (utg) {
      nutg++;

      for (DoveTailIter it=utg->dovetail_path_ptr->begin(); it != utg->dovetail_path_ptr->end(); it++) {
        nfrg++;

        if (it->ident > _fi->numFragments())
          fprintf(logFile, "HUH?  ident=%d numfrags=%d\n", it->ident, _fi->numFragments());

        inUnitig[it->ident] = ti;

        len = MAX(len, it->position.bgn);
        len = MAX(len, it->position.end);
      }

      logSize[ (uint32)(log10(len) / log10(2)) ]++;
    }
  }

  int lost = 0;
  int found = 0;

  for (uint32 i=0; i<_fi->numFragments()+1; i++) {
    if (_fi->fragmentLength(i) > 0) {
      if (inUnitig[i] == 0) {
        fprintf(logFile, "ERROR frag %d is in unitig 0!\n", i);
      } else if (inUnitig[i] != noUnitig) {
        found++;
      } else {
        fprintf(logFile, "ERROR frag %d disappeared!\n", i);
        lost++;
      }
    }
  }

  fprintf(logFile, "checkUnitigMembership()-- nutg=%d nfrg=%d lost=%d found=%d\n", nutg, nfrg, lost, found);

  fprintf(logFile, "checkUnitigMembership()-- log2 length histogram:\n");
  for (uint32 i=5; i<20; i++)
    fprintf(logFile, "checkUnitigMembership()-- %2u (%9u-%9u) %u\n", i, (uint32)1 << i, (uint32)1 << (i+1), logSize[i]);

  assert(lost == 0);
}


//  For every unitig, report the best overlaps contained in the
//  unitig, and all overlaps contained in the unitig.
void
UnitigGraph::reportOverlapsUsed(const char *filename) {

  if (logFileFlagSet(LOG_OVERLAPS_USED) == 0)
    return;

  FILE *F = fopen(filename, "w");

  if (F == NULL)
    return;

  for (uint32  ti=0; ti<unitigs->size(); ti++) {
    Unitig  *utg = (*unitigs)[ti];

    if (utg == NULL)
      continue;

    for (uint32 fi=0; fi<utg->dovetail_path_ptr->size(); fi++) {
      DoveTailNode  *frg = &(*utg->dovetail_path_ptr)[fi];

      //  Where is our best overlap?  Contained or dovetail?

      BestEdgeOverlap *bestedge5 = bog_ptr->getBestEdgeOverlap(frg->ident, FIVE_PRIME);
      BestEdgeOverlap *bestedge3 = bog_ptr->getBestEdgeOverlap(frg->ident, THREE_PRIME);

      uint32           bestident5 = 0;
      uint32           bestident3 = 0;

      if (bestedge5)
        bestident5 = bestedge5->frag_b_id;

      if (bestedge3)
        bestident3 = bestedge3->frag_b_id;

      //  Now search ahead, reporting any overlap to any fragment.
      //
      for (uint32 oi=fi+1; oi<utg->dovetail_path_ptr->size(); oi++) {
        DoveTailNode  *ooo = &(*utg->dovetail_path_ptr)[oi];

        int frgbgn = MIN(frg->position.bgn, frg->position.end);
        int frgend = MAX(frg->position.bgn, frg->position.end);

        int ooobgn = MIN(ooo->position.bgn, ooo->position.end);
        int oooend = MAX(ooo->position.bgn, ooo->position.end);

        if ((frgbgn <= ooobgn) && (ooobgn + 40 < frgend)) {
          BestContainment *bestcont  = bog_ptr->getBestContainer(ooo->ident);

          uint32           bestident  = 0;
          if (bestcont)
            bestident = bestcont->container;

          bool isBest = ((frg->ident == bestident) ||
                         (ooo->ident == bestident5) ||
                         (ooo->ident == bestident3));

          fprintf(F, "%d\t%d%s\n", frg->ident, ooo->ident, (isBest) ? ((bestident) ? "\tbc" : "\tbe") : "");
        }

        if (frgend < ooobgn)
          break;
      }
    }
  }

  fclose(F);
}


void
UnitigGraph::reportUnitigs(const char *filename) {

  if (logFileFlagSet(LOG_INTERMEDIATE_UNITIGS) == 0)
    return;

  char tigStorePath[FILENAME_MAX];
  sprintf(tigStorePath, "debug.%s.tigStore", filename);

  writeIUMtoFile(tigStorePath, tigStorePath, 90000, false);
}


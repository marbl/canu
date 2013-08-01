
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

static const char *rcsid = "$Id$";

#include "AS_BOG_Datatypes.H"
#include "AS_BOG_UnitigGraph.H"
#include "AS_BOG_BestOverlapGraph.H"

#include "MultiAlignStore.H"



void
UnitigGraph::checkUnitigMembership(void) {
  int nutg = 0;
  int nfrg = 0;

  fprintf(logFile, "checkUnitigMembership()--  numfrags=%d\n", FI->numFragments());

  uint32 *inUnitig = new uint32 [FI->numFragments()+1];
  uint32  logSizeMax  = 0;
  uint32  logSize[64] = {0};

  for (uint32 i=0; i<FI->numFragments()+1; i++)
    inUnitig[i] = noUnitig;

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];
    uint32   len = 0;

    if (tig) {
      nutg++;

      for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
        ufNode  *frg = &tig->ufpath[fi];
        nfrg++;

        if (frg->ident > FI->numFragments())
          fprintf(logFile, "HUH?  ident=%d numfrags=%d\n", frg->ident, FI->numFragments());

        inUnitig[frg->ident] = ti;

        len = MAX(len, frg->position.bgn);
        len = MAX(len, frg->position.end);
      }

      uint32  ls = (uint32)(log10(len) / log10(2));
      logSizeMax = (logSizeMax < ls) ? ls : logSizeMax;
      logSize[ls]++;
    }
  }

  int lost = 0;
  int found = 0;

  for (uint32 i=0; i<FI->numFragments()+1; i++) {
    if (FI->fragmentLength(i) > 0) {
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
  for (uint32 i=5; i<=logSizeMax; i++)
    fprintf(logFile, "checkUnitigMembership()-- %2u (%9u-%9u) %u\n", i, (uint32)1 << i, (uint32)1 << (i+1), logSize[i]);

  assert(lost == 0);

  delete [] inUnitig;
}


//  For every unitig, report the best overlaps contained in the
//  unitig, and all overlaps contained in the unitig.
void
UnitigGraph::reportOverlapsUsed(const char *prefix, const char *name) {

  if (logFileFlagSet(LOG_OVERLAPS_USED) == 0)
    return;

  char ovlPath[FILENAME_MAX];
  sprintf(ovlPath, "%s.%03u.%s.overlaps", prefix, logFileOrder, name);

  FILE *F = fopen(ovlPath, "w");

  if (F == NULL)
    return;

  for (uint32  ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    for (uint32 fi=0; fi<utg->ufpath.size(); fi++) {
      ufNode  *frg = &utg->ufpath[fi];

      //  Where is our best overlap?  Contained or dovetail?

      BestEdgeOverlap *bestedge5 = OG->getBestEdgeOverlap(frg->ident, false);
      BestEdgeOverlap *bestedge3 = OG->getBestEdgeOverlap(frg->ident, true);

      uint32           bestident5 = 0;
      uint32           bestident3 = 0;

      if (bestedge5)
        bestident5 = bestedge5->fragId();

      if (bestedge3)
        bestident3 = bestedge3->fragId();

      //  Now search ahead, reporting any overlap to any fragment.
      //
      for (uint32 oi=fi+1; oi<utg->ufpath.size(); oi++) {
        ufNode  *ooo = &utg->ufpath[oi];

        int frgbgn = MIN(frg->position.bgn, frg->position.end);
        int frgend = MAX(frg->position.bgn, frg->position.end);

        int ooobgn = MIN(ooo->position.bgn, ooo->position.end);
        int oooend = MAX(ooo->position.bgn, ooo->position.end);

        if ((frgbgn <= ooobgn) && (ooobgn + 40 < frgend)) {
          BestContainment *bestcont  = OG->getBestContainer(ooo->ident);

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
UnitigGraph::reportUnitigs(const char *prefix, const char *name) {

  if (logFileFlagSet(LOG_INTERMEDIATE_UNITIGS) == 0)
    return;

  uint32  numFragsT  = 0;
  uint32  numFragsP  = 0;

  //  Compute average frags per partition.
  for (uint32  ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    numFragsT += utg->ufpath.size();
  }

  numFragsP = numFragsT / 7;

  char tigStorePath[FILENAME_MAX];
  sprintf(tigStorePath, "%s.%03u.%s.tigStore", prefix, logFileOrder, name);

  writeIUMtoFile(tigStorePath, tigStorePath, numFragsP, false);
}


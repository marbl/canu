
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_Instrumentation.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-27
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2014-DEC-23
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static const char *rcsid = "$Id$";

#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_SetParentAndHang.H"
#include "AS_BAT_Outputs.H"

void
checkUnitigMembership(UnitigVector &unitigs) {
  int nutg = 0;
  int nfrg = 0;

  writeLog("checkUnitigMembership()--  numfrags=%d\n", FI->numFragments());

  uint32 *inUnitig = new uint32 [FI->numFragments()+1];
  uint32  logSizeMax  = 0;
  uint32  logSize[64] = {0};

  for (uint32 i=0; i<FI->numFragments()+1; i++)
    inUnitig[i] = noUnitig;

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];
    int32    len = 0;

    if (tig) {
      nutg++;

      for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
        ufNode  *frg = &tig->ufpath[fi];
        nfrg++;

        if (frg->ident > FI->numFragments())
          writeLog("HUH?  ident=%d numfrags=%d\n", frg->ident, FI->numFragments());

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
        writeLog("ERROR frag %d is in unitig 0!\n", i);
      } else if (inUnitig[i] != noUnitig) {
        found++;
      } else {
        writeLog("ERROR frag %d disappeared!\n", i);
        lost++;
      }
    }
  }

  writeLog("checkUnitigMembership()-- nutg=%d nfrg=%d lost=%d found=%d\n", nutg, nfrg, lost, found);

  writeLog("checkUnitigMembership()-- log2 length histogram:\n");
  for (uint32 i=5; i<=logSizeMax; i++)
    writeLog("checkUnitigMembership()-- %2u (%9u-%9u) %u\n", i, (uint32)1 << i, (uint32)1 << (i+1), logSize[i]);

  assert(lost == 0);

  delete [] inUnitig;
}


//  For every unitig, report the best overlaps contained in the
//  unitig, and all overlaps contained in the unitig.
void
reportOverlapsUsed(UnitigVector &unitigs, const char *prefix, const char *name) {

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
          if (bestcont->isContained)
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
reportUnitigs(UnitigVector &unitigs, const char *prefix, const char *name) {

  if (logFileFlagSet(LOG_INTERMEDIATE_UNITIGS) == 0)
    return;

  uint32  numFragsT  = 0;
  uint32  numFragsP  = 0;
  uint64  utgLen     = 0;

  //  Compute average frags per partition.
  for (uint32  ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    numFragsT += utg->ufpath.size();

    if (utg->ufpath.size() > 2)
      utgLen    += utg->getLength();
  }

  if      (utgLen < 16 * 1024 * 1024)
    numFragsP = numFragsT / 7;
  else if (utgLen < 64 * 1024 * 1024)
    numFragsP = numFragsT / 63;
  else
    numFragsP = numFragsT / 127;

  char tigStorePath[FILENAME_MAX];
  sprintf(tigStorePath, "%s.%03u.%s.tigStore", prefix, logFileOrder, name);

  //  Failing to do this results in consensus running about 40 times slower.  Three hours instead of
  //  five minutes.
  setParentAndHang(unitigs);

  writeUnitigsToStore(unitigs, tigStorePath, tigStorePath, numFragsP, false);
}


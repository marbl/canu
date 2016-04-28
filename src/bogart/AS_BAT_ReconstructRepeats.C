
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
 *    src/AS_BAT/AS_BAT_ReconstructRepeats.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-JAN-05 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2015-APR-24
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_ChunkGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_PlaceFragUsingOverlaps.H"

#include "AS_BAT_PopulateUnitig.H"
#include "AS_BAT_PlaceContains.H"



//  estimate read error rate from best overlaps (per library?)
//  use that error rate below when rebuilding repeats

void
reconstructRepeats(UnitigVector &unitigs,
                   double        erateGraph,
                   double        deviationGraph) {

  //  Build a set<> of all the unplaced fragments, then construct a new BOG and CG from which we
  //  construct unitigs.

  BestOverlapGraph  *OGsave = OG;
  ChunkGraph        *CGsave = CG;

  set<uint32>        unplaced;

  for (uint32 fi=1; fi<=FI->numFragments(); fi++)
    if (Unitig::fragIn(fi) == 0)
      unplaced.insert(fi);

  OG = new BestOverlapGraph(erateGraph / 2.0, deviationGraph, &unplaced);
  CG = new ChunkGraph(&unplaced);

  writeLog("==> BUILDING REPEAT UNITIGS from %d fragments.\n", unplaced.size());

  for (uint32 fi=CG->nextFragByChunkLength(); fi>0; fi=CG->nextFragByChunkLength())
    populateUnitig(unitigs, fi);

  writeLog("==> BUILDING REPEAT UNITIGS catching missed fragments.\n");

  for (uint32 fi=1; fi <= FI->numFragments(); fi++)
    populateUnitig(unitigs, fi);

  writeLog("==> BUILDING REPEAT UNITIGS placing contained fragments.\n");

  placeUnplacedUsingAllOverlaps(unitigs, "PREFIX");

  delete OG;
  delete CG;

  OG = OGsave;
  CG = CGsave;
}

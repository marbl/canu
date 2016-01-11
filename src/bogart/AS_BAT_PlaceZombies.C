
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
 *    src/AS_BAT/AS_BAT_PlaceZombies.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2015-APR-24
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_PlaceZombies.H"


//  Zombies are caused by a fragment being contained in a fragment that is eventually contained in
//  the original fragment -- circular containmnents.
//
//  Here we detect Zombies, and reset their best container to something that is already placed.

void
placeZombies(UnitigVector &unitigs, double erate) {

  writeLog("==> SEARCHING FOR ZOMBIES\n");

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

    Unitig      *utg = unitigs.newUnitig(false);
    ufNode       frg;

    frg.ident             = i;
    frg.contained         = 0;
    frg.parent            = 0;

    frg.ahang             = 0;
    frg.bhang             = 0;

    frg.position.bgn      = 0;
    frg.position.end      = FI->fragmentLength(i);

    frg.containment_depth = 0;

    utg->addFrag(frg, 0, false);

    writeLog("placeZombies()-- unitig %d created from zombie fragment %d\n",
            utg->id(), i);
    numZombies++;
  }

  writeLog("RESURRECTED %d ZOMBIE FRAGMENT%s.\n", numZombies, (numZombies != 1) ? "s" : "");

  delete [] inUnitig;
}

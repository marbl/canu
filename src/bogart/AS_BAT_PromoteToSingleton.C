
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
 *    src/AS_BAT/AS_BAT_PromoteToSingleton.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-JAN-05 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_Unitig.H"

//  If we are not reconstructing repeats, promote all the unplaced fragments to new unitigs.
//  Oodles of possibilities here; promote everything to a singleton unitig, promote only
//  the non-contained, then place contains, then promote what is left over, etc.

void
promoteToSingleton(UnitigVector &unitigs) {

  for (uint32 fi=1; fi<=FI->numFragments(); fi++) {
    if (Unitig::fragIn(fi) != 0)
      //  Placed already
      continue;

    if (FI->fragmentLength(fi) == 0)
      //  Deleted.
      continue;

    Unitig *utg = unitigs.newUnitig(false);
    ufNode  frag;

    frag.ident             = fi;
    frag.contained         = 0;
    frag.parent            = 0;
    frag.ahang             = 0;
    frag.bhang             = 0;
    frag.position.bgn      = 0;
    frag.position.end      = FI->fragmentLength(fi);

    utg->addFrag(frag, 0, false);
  }
}

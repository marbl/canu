
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
 *    src/AS_BAT/AS_BAT_PopulateUnitig.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
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

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_PopulateUnitig.H"



void
populateUnitig(Unitig           *unitig,
               BestEdgeOverlap  *bestnext) {

  assert(unitig->getLength() > 0);

  if ((bestnext == NULL) || (bestnext->fragId() == 0))
    //  Nothing to add!
    return;

  ufNode  frag    = unitig->ufpath.back();

  //  The ID of the last fragment in the unitig, and the end we should walk off of it.
  int32   lastID  = frag.ident;
  bool    last3p  = (frag.position.bgn < frag.position.end);

  uint32  nAdded  = 0;

  //  While there are fragments to add AND those fragments to add are not already in a unitig,
  //  construct a reverse-edge, and add the fragment.

  while ((bestnext->fragId() != 0) &&
         (Unitig::fragIn(bestnext->fragId()) == 0)) {
    BestEdgeOverlap  bestprev;

    //  Reverse nextedge (points from the unitig to the next fragment to add) so that it points from
    //  the next fragment to add back to something in the unitig.  If the fragments are
    //  innie/outtie, we need to reverse the overlap to maintain that the A fragment is forward.

    if (last3p == bestnext->frag3p())
      bestprev.set(lastID, last3p, bestnext->bhang(), bestnext->ahang(), bestnext->evalue());
    else
      bestprev.set(lastID, last3p, -bestnext->ahang(), -bestnext->bhang(), bestnext->evalue());

    //  We just made 'bestprev' pointing from read 'bestnext->fragId()' end 'bestnext->frag3p()'
    //  back to read 'lastID' end 'last3p'.  Compute the placement.

    if (unitig->placeFrag(frag, bestnext->fragId(), bestnext->frag3p(), &bestprev)) {
      unitig->addFrag(frag, 0, false);
      nAdded++;

    } else {
      writeLog("ERROR:  Failed to place frag %d into BOG path.\n", frag.ident);
      assert(0);
    }

    //  Set up for the next fragmnet

    lastID  = frag.ident;
    last3p  = (frag.position.bgn < frag.position.end);

    bestnext = OG->getBestEdgeOverlap(lastID, last3p);
  }

  if (logFileFlagSet(LOG_BUILD_UNITIG))
    if (bestnext->fragId() == 0)
      writeLog("Stopped adding at frag %u/%c' because no next best edge.  Added %u reads.\n",
               lastID, (last3p) ? '3' : '5',
               nAdded);
    else
      writeLog("Stopped adding at frag %u/%c' beacuse next best frag %u/%c' is in unitig %u.  Added %u reads.\n",
               lastID, (last3p) ? '3' : '5',
               bestnext->fragId(), bestnext->frag3p() ? '3' : '5',
               Unitig::fragIn(bestnext->fragId()),
               nAdded);
}




void
populateUnitig(UnitigVector &unitigs,
               int32 fi) {

  if ((FI->fragmentLength(fi) == 0) ||  //  Skip deleted
      (Unitig::fragIn(fi) != 0) ||      //  Skip placed
      (OG->isContained(fi) == true))    //  Skip contained
    return;

  Unitig *utg = unitigs.newUnitig(logFileFlagSet(LOG_BUILD_UNITIG));

  //  Add a first fragment -- to be 'compatable' with the old code, the first fragment is added
  //  reversed, we walk off of its 5' end, flip it, and add the 3' walk.

  ufNode  frag;

  frag.ident             = fi;
  frag.contained         = 0;
  frag.parent            = 0;
  frag.ahang             = 0;
  frag.bhang             = 0;
  frag.position.bgn      = FI->fragmentLength(fi);
  frag.position.end      = 0;

  utg->addFrag(frag, 0, logFileFlagSet(LOG_BUILD_UNITIG));

  //  Add fragments as long as there is a path to follow...from the 3' end of the first fragment.

  BestEdgeOverlap  *bestedge5 = OG->getBestEdgeOverlap(fi, false);
  BestEdgeOverlap  *bestedge3 = OG->getBestEdgeOverlap(fi, true);

  assert(bestedge5->ahang() <= 0);  //  Best Edges must be dovetail, which makes this test
  assert(bestedge5->bhang() <= 0);  //  much simpler.
  assert(bestedge3->ahang() >= 0);
  assert(bestedge3->bhang() >= 0);

  //  If this fragment is not covered by the two best overlaps we are finished.  We will not follow
  //  the paths out.  This indicates either low coverage, or a chimeric fragment.  If it is low
  //  coverage, then the best overlaps will be mutual and we'll recover the same path.  If it is a
  //  chimeric fragment the overlaps will not be mutual and we will skip this fragment.
  //
  //  The amount of our fragment that is covered by the two best overlaps is
  //
  //    (fragLen + bestedge5->bhang()) + (fragLen - bestedge3->ahang())
  //
  //  If that is not significantly longer than the fragment length, then we will not use this
  //  fragment as a seed for unitig construction.
  //

  if (OG->isSuspicious(fi))
    return;

#if 0
  uint32  covered = FI->fragmentLength(fi) + bestedge5->bhang() + FI->fragmentLength(fi) - bestedge3->ahang();

  //  This breaks unitigs at 0x best-coverage regions.  There might be a contain that spans (joins)
  //  the two best overlaps to verify the fragment, but we can't easily tell right now.
  if (covered < FI->fragmentLength(fi) + AS_OVERLAP_MIN_LEN / 2) {
    writeLog("Stopping unitig construction of suspicious frag %d in unitig %d\n",
            utg->ufpath.back().ident, utg->id());
    return;
  }
#endif

  if (logFileFlagSet(LOG_BUILD_UNITIG))
    writeLog("Adding 5' edges off of frag %d in unitig %d\n",
            utg->ufpath.back().ident, utg->id());

  if (bestedge5->fragId())
    populateUnitig(utg, bestedge5);

  utg->reverseComplement(false);

  if (logFileFlagSet(LOG_BUILD_UNITIG))
    writeLog("Adding 3' edges off of frag %d in unitig %d\n",
            utg->ufpath.back().ident, utg->id());

  if (bestedge3->fragId())
    populateUnitig(utg, bestedge3);

  //  Enabling this reverse complement is known to degrade the assembly.  It is not known WHY it
  //  degrades the assembly.
  //
  //utg->reverseComplement(false);
}

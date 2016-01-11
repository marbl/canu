
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
 *    src/AS_BAT/AS_BAT_Breaking.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-DEC-06 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_BREAKING
#define INCLUDE_AS_BAT_BREAKING

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"

//  What to do with contained fragments at the split point?
//
//  ejectContains -- eject all contains from the unitig, client is expected to
//  place them again.
//
//  keepContains -- contains that are associated with the fragEnd are moved to
//  the new unitig.  There are five cases:
//
//  A #-----------------*
//  1   ------                           - contained in A
//  2        -------                     - contained in A, dovetail to B
//  B          -----------------------
//  3             ------                 - contained in A, contained in B
//  4                 -------            - dovetail to A, contained in B
//  5                       -------      - contained in B
//
//  If the break point is at #:
//    keepContains == true  -- 1, 2, 3 will move with A.
//    keepContains == false -- 1, 2    will move with A.
//
//  If the break point is at *:
//    keepContains == true  -- 1, 2, 3 will move with A.
//    keepContains == false -- 1, 2    will move with A.
//

class breakPoint {
public:
  breakPoint() {
    fragEnd             = FragmentEnd();
    ejectContains       = false;
    keepContains        = false;
    break5p             = false;
    break3p             = false;
    rememberLastTig     = false;
    searchDiscontinuous = false;
  };

  breakPoint(uint32 isectFrg, bool isect3p, bool eject, bool keep) {
    fragEnd             = FragmentEnd(isectFrg, isect3p);
    ejectContains       = eject;
    keepContains        = keep;
    break5p             = false;
    break3p             = false;
    rememberLastTig     = false;
    searchDiscontinuous = false;
  };

  bool operator<(breakPoint const that) const {
    return(fragEnd < that.fragEnd);
  };

public:
  FragmentEnd fragEnd;

  bool        ejectContains;
  bool        keepContains;

  bool        break5p;
  bool        break3p;
  bool        rememberLastTig;
  bool        searchDiscontinuous;
};



bool
breakUnitigAt(UnitigVector &unitigs, Unitig *tig, vector<breakPoint> &breaks, bool doDelete);

#endif //  INCLUDE_AS_BAT_BREAKING

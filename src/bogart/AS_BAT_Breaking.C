
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
 *    src/AS_BAT/AS_BAT_Breaking.C
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Breaking.H"

#include "AS_BAT_BestOverlapGraph.H"

#define LOG_ADDUNITIG_BREAKING  1
#define LOG_ADDFRAG_BREAKING    0

//  The four cases we must handle:
//
//  -------------------------A
//        --
//              ---------------------------
//                  ------
//                    B----------------------------
//                        ---
//                                  ------
//
//  When at A:
//    keepContains == true  -- Remember lastUnitig, and the coordinate of A.
//                             Any fragment that ends before A is placed in lastUnitig.
//                             When a fragment begins after A, we can forget lastUnitig.
//
//    keepContains == false -- Remember lastUnitig. and the coordinate of A.
//                             Search forward until the first fragment that starts after A,
//                             search for any gaps in the layout.  Move all fragments before
//                             the last gap to lastUnitig.  Gap is defined as being less
//                             than an overlap of overlap.
//
//  When at B:
//    keepContains == true  -- Do nothing special.
//
//    keepContains == false -- Same as "A: keepContains == true".
//
//  UGLIES:
//
//    1) When A and B are on the same fragment, both breakPoints must have set keepContains
//       consistently.  This could arise if the fragment is a chimera.  In this case, we definitely
//       do not want to keep contains.  So, instead of asserting (we'll probably do that, just to
//       see if this occurs) we'll set keepContains only if ALL breakpoints request it.
//
//    2) Rapid fire breakpoints.  When breakpoints closely follow each other -- specificially when
//       the fragments they break on are overlapping -- the meaning of keepContains becomes horribly
//       confused.  The algorithm will keep lastTig updated to the last break point seen.  The
//       picture above (without trying) shows this case.  The tiny fragment second from the end is
//       contained in both A and B.  It will end up with B if that is keepContains=true, regardless
//       of what A said.


static
bool
processBreakpoints(vector<breakPoint>       &breaks,
                   map<uint32, breakPoint>  &breakMap) {
  bool ejectContains = false;

  for (uint32 bp=0; bp<breaks.size(); bp++) {
    uint32  fid = breaks[bp].fragEnd.fragId();

    //writeLog("found break %u/%c' eject=%d keep=%d\n",
    //        breaks[bp].fragEnd.fragId(), (breaks[bp].fragEnd.frag3p() ? '3' : '5'),
    //        breaks[bp].ejectContains, breaks[bp].keepContains);

    if (breakMap.find(fid) == breakMap.end()) {
      breakMap[fid] = breaks[bp];

      ejectContains              |= breaks[bp].ejectContains;
      breakMap[fid].break3p       = (breaks[bp].fragEnd.frag3p() == true);
      breakMap[fid].break5p       = (breaks[bp].fragEnd.frag3p() == false);

    } else {
      ejectContains              |= breaks[bp].ejectContains;
      breakMap[fid].break3p      |= (breaks[bp].fragEnd.frag3p() == true);
      breakMap[fid].break5p      |= (breaks[bp].fragEnd.frag3p() == false);
      breakMap[fid].keepContains &= breaks[bp].keepContains;
    }
  }

  return(ejectContains);
}


static
bool
isBreakpoint(ufNode &frg,
             map<uint32, breakPoint> &breakMap,
             bool &break5p, bool &break3p,
             bool &rememberLastTig,
             bool &searchDiscontinuous) {

  if (breakMap.find(frg.ident) == breakMap.end())
    return(false);

  breakPoint &bp = breakMap[frg.ident];

  assert(bp.fragEnd.fragId() == frg.ident);

  //writeLog("BREAK %u %d/%d\n",
  //        bp.fragEnd.fragId(), bp.break3p, bp.break5p);

  break3p             = bp.break3p;
  break5p             = bp.break5p;
  rememberLastTig     = false;
  searchDiscontinuous = false;

  bool    frgReversed = (frg.position.end < frg.position.bgn);
  bool    isFarEnd    = false;

  if (((break3p == true) && (frgReversed == false)) ||
      ((break5p == true) && (frgReversed == true)))
    isFarEnd = true;

  //  Remember the lastTig if we are case A (isFarEnd == true) and keepContains is true,
  //  of if we are case B (isFarEnd == false) and keepContains is false.
  if (isFarEnd == bp.keepContains)
    rememberLastTig = true;

  //  Do the painful search for disconnects if we are case A and keepContains is false.
  if ((isFarEnd == true) && (bp.keepContains == false))
    searchDiscontinuous = true;

  return(true);
}







bool
breakUnitigAt(UnitigVector       &unitigs,
              Unitig             *tig,
              vector<breakPoint> &breaks,
              bool                doDelete) {
  uint32  newTigs = 0;

  if (breaks.empty())
    return(false);

  // we cannot predict the number of new unitigs created from the number of break points.  to
  // prevent infinite splitting loops, we need to count the number of new unitigs we create here,
  // and return 'no work done' if only one new unitig is created (aka, if we just copied all
  // fragments from the old unitig to the new unitig)

  Unitig       *lastTig    = NULL;  //  For saving contained frags in the last unitig constructed.
  int32         lastOffset = 0;
  int32         lastLength = 0;

  Unitig       *saveTig    = NULL;  //  The last unitig constructed.
  int32         saveOffset = 0;

  Unitig       *currTig    = NULL;  //  The current unitig being constructed.  Not guaranteed to exist.
  int32         currOffset = 0;

  uint32        tigsCreated = 0;

  //  Examine the break points.  Merge multiple break points on a single fragment together.  Put the
  //  break points into a map -- without this, we'd need to make sure the break points are sorted in
  //  the same order as the fragments in the unitig.

  map<uint32, breakPoint>  breakMap;

  bool ejectContains = processBreakpoints(breaks, breakMap);

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode  frg         = tig->ufpath[fi];
    bool    frgReversed = (frg.position.end < frg.position.bgn);

    bool    break5p             = false;
    bool    break3p             = false;
    bool    rememberLastTig     = false;
    bool    searchDiscontinuous = false;

    if (ejectContains && OG->isContained(frg.ident)) {
      Unitig::removeFrag(frg.ident);
      continue;
    }

    //  Current fragment is the one we want to break on.  Figure out which end(s) to break on -- we
    //  might (stupidly) have requested to break on both ends -- and if we should be remembering
    //  lastTig.

    if (isBreakpoint(frg,
                     breakMap,
                     break5p, break3p,
                     rememberLastTig,
                     searchDiscontinuous)) {
      //writeLog("NEW BREAK at frag %u 5p=%d 3p=%d remember=%d search=%d\n",
      //        frg.ident, break5p, break3p, rememberLastTig, searchDiscontinuous);

      lastTig    = NULL;
      lastOffset = 0;
      lastLength = 0;

      if ((rememberLastTig) && (saveTig != NULL)) {
        lastTig    = saveTig;
        lastOffset = saveOffset;
        lastLength = saveOffset + saveTig->getLength();
      }

      //  Should we clear saveTig?  This is always set to the last unitig created.  It is
      //  possible to get two breakpoints with no unitig creation (break at the high end of one
      //  fragment, and the low end of the next fragment) in which case we still want to have
      //  the last unitig...and not NULL...saved.  So, no, don't clear.
      //
      //saveTig    = NULL;
      //saveOffset = 0;

      if (searchDiscontinuous) {
      }
    }


    //  If both break5p and break3p, this fragment is ejected to a singleton.  This is an easy case,
    //  and we'll get it out of the way.

    if ((break5p == true) &&
        (break3p == true)) {
      newTigs++;
      saveTig    = currTig    = unitigs.newUnitig(false);  //  LOG_ADDUNITIG_BREAKING);
      saveOffset = currOffset = (frgReversed) ? -frg.position.end : -frg.position.bgn;

      currTig->addFrag(frg, currOffset, LOG_ADDFRAG_BREAKING);

      currTig = NULL;

      continue;
    }

    //  If neither break and we're saving contains, add to the last unitig.

    if ((break5p == false) &&
        (break3p == false) &&
        (lastTig) &&
        (frg.position.bgn < lastLength) &&
        (frg.position.end < lastLength)) {
      lastTig->addFrag(frg, lastOffset, LOG_ADDFRAG_BREAKING);
      continue;
    }

    //  If neither break, just add the fragment to the existing unitig.

    if ((break5p == false) &&
        (break3p == false)) {
      if (currTig == NULL) {
        newTigs++;
        saveTig    = currTig    = unitigs.newUnitig(false);  //  LOG_ADDUNITIG_BREAKING);
        saveOffset = currOffset = (frgReversed) ? -frg.position.end : -frg.position.bgn;
      }

      currTig->addFrag(frg, currOffset, LOG_ADDFRAG_BREAKING);
      continue;
    }

    //  Breaking at the left end of the fragment.  This fragment starts a new unitig.

    if ((break5p && (frgReversed == false)) ||
        (break3p && (frgReversed == true))) {
      newTigs++;
      saveTig    = currTig    = unitigs.newUnitig(false);  //  LOG_ADDUNITIG_BREAKING);
      saveOffset = currOffset = (frgReversed) ? -frg.position.end : -frg.position.bgn;

      currTig->addFrag(frg, currOffset, LOG_ADDFRAG_BREAKING);

      continue;
    }

    //  Breaking at the right end of the fragment.  This fragment ends the existing unitig (which
    //  might not even exist).

    if ((break5p && (frgReversed == true)) ||
        (break3p && (frgReversed == false))) {
      if (currTig == NULL) {
        newTigs++;
        saveTig    = currTig    = unitigs.newUnitig(false);  //  LOG_ADDUNITIG_BREAKING);
        saveOffset = currOffset = (frgReversed) ? -frg.position.end : -frg.position.bgn;
      }

      currTig->addFrag(frg, currOffset, LOG_ADDFRAG_BREAKING);

      currTig = NULL;

      continue;
    }
  }

  if (newTigs == 0)
    return(false);

  if (doDelete) {
    unitigs[tig->id()] = NULL;
    delete tig;
  }

  return(true);
}


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
 *    src/AS_BAT/AS_BAT_Unitig.C
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

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

static std::map<uint32,int>* containPartialOrder;

uint32* Unitig::_inUnitig     = NULL;
uint32* Unitig::_pathPosition = NULL;


#warning WHAT REALLLY HAPPENS IF NO BACKBONE NODE, OR NO PREVIOUS BACKBONE NODE

ufNode Unitig::getLastBackboneNode(void) {
  for (int32 fi=ufpath.size()-1; fi >= 0; fi--) {
    ufNode  &node = ufpath[fi];

    if (node.contained)
      continue;

    return(node);
  }

  writeLog("Unitig::getLastBackboneNode()--  WARNING:  unitig %d has no backbone nodes, all contained!\n", id());
  ufNode last;
  memset(&last, 0, sizeof(ufNode));
  return(last);
}


ufNode Unitig::getLastBackboneNode(uint32 &prevID) {
  ufNode last;

  memset(&last, 0, sizeof(ufNode));

  prevID = 0;

  for (int32 fi=ufpath.size()-1; (fi >= 0) && (prevID == 0); fi--) {
    ufNode  *node = &ufpath[fi];

    if (node->contained)
      continue;

    if (last.ident == 0)
      //  Save the last dovetail node, but keep looking....
      last = *node;
    else
      //  ...for the next to last ID.
      prevID = node->ident;
  }

  return(last);
}




void Unitig::reverseComplement(bool doSort) {

  //  If there are contained fragments, we need to sort by position to place them correctly after
  //  their containers.  If there are no contained fragments, sorting can break the initial unitig
  //  building.  When two frags start at position zero, we'll exchange the order.  Initial unitig
  //  building depends on having the first fragment added become the last fragment in the unitig
  //  after reversing.

  for (uint32 fi=0; fi<ufpath.size(); fi++) {
    ufNode  *frg = &ufpath[fi];

    frg->position.bgn = getLength() - frg->position.bgn;
    frg->position.end = getLength() - frg->position.end;

    //if (frg->contained != 0)
    //  doSort = true;

    assert(frg->position.bgn >= 0);
    assert(frg->position.end >= 0);
  }

  //  We've updated the positions of everything.  Now, sort or reverse the list, and rebuild the
  //  pathPosition map.

  if (doSort) {
    sort();
  } else {
    std::reverse(ufpath.begin(), ufpath.end());

    for (uint32 fi=0; fi<ufpath.size(); fi++)
      _pathPosition[ufpath[fi].ident] = fi;
  }
}



int
ufNodeCmp(const void *a, const void *b){
  ufNode *impa = (ufNode *)a;
  ufNode *impb = (ufNode *)b;

  int32 abgn = (impa->position.bgn < impa->position.end) ? impa->position.bgn : impa->position.end;
  int32 aend = (impa->position.bgn < impa->position.end) ? impa->position.end : impa->position.bgn;

  int32 bbgn = (impb->position.bgn < impb->position.end) ? impb->position.bgn : impb->position.end;
  int32 frag3p = (impb->position.bgn < impb->position.end) ? impb->position.end : impb->position.bgn;

  //  NEWSORT does not work.  When bubbles are popped, we add non-contained fragments to
  //  a unitig, but just stick them at the end of the list.  NEWSORT would then maintain
  //  this ordering, which is an error.
  //
#undef NEWSORT

#ifdef NEWSORT
  bool  aIsCont = OG->isContained(impa->ident);
  bool  bIsCont = OG->isContained(impb->ident);

  if ((aIsCont == false) && (bIsCont == false))
    //  Both dovetail nodes, keep same order
    return((int)impa->containment_depth - (int)impb->containment_depth);
#endif

  if (abgn != bbgn)
    //  Return negative for the one that starts first.
    return(abgn - bbgn);

  if (aend != frag3p)
    //  Return negative for the one that ends last.
    return(frag3p - aend);

#ifdef NEWSORT
  if (bIsCont == true)
    //  b is contained in a, so it comes after a.
    return(-1);

 if (aIsCont == true)
    //  a is contained in b, so it comes after b.
    return(1);
#endif

  //  Both contained, fallback on depth added, negative for earliest added
  return((int)impa->containment_depth - (int)impb->containment_depth);
}


void
Unitig::sort(void) {

#ifdef NEWSORT
  for (int fi=0; fi<ufpath.size(); fi++) {
    ufNode *f = &(ufpath[fi]);

    if (OG->isContained(f->ident) == false)
      f->containment_depth = fi;
  }
#endif

  qsort( &(ufpath.front()), getNumFrags(), sizeof(ufNode), &ufNodeCmp );

  for (uint32 fi=0; fi<ufpath.size(); fi++)
    _pathPosition[ufpath[fi].ident] = fi;
}

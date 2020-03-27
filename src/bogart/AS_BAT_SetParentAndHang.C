
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_TigVector.H"

#include "AS_BAT_SetParentAndHang.H"

#include <map>

using namespace std;

void
setParentAndHang(TigVector &tigs) {

  return;

  map<uint32,bool>  forward;
  map<uint32,bool>  allreads;

  //  Just for stats, build a map fo the reads in the unitig.


  for (uint32 ti=0; ti<tigs.size(); ti++) {
    Unitig        *tig = tigs[ti];

    if (tig == NULL)
      continue;

    if (tig->ufpath.size() == 0)
      continue;

    //  Reset parent and hangs, build a map of the reads in the unitig.

    for (uint32 fi=1; fi<tig->ufpath.size(); fi++) {
      ufNode *frg = &tig->ufpath[fi];

      frg->parent          = 0;
      frg->ahang           = 0;
      frg->bhang           = 0;

      allreads[frg->ident] = true;
    }

    //  For each read, set parent/hangs using the edges.

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode *frg  = &tig->ufpath[fi];



    //  Remember that we've placed this read, and if it was forward or reverse.
    forward[frg->ident] = (frg->position.bgn < frg->position.end);

    //  If the first read, there is no parent possible.
    if (ti == 0)
      continue;

    //  Otherwise, find the thickest overlap to any read already placed in the unitig.

    uint32         olapsLen = 0;
    BAToverlap    *olaps = OC->getOverlaps(frg->ident, olapsLen);

    uint32         tt     = UINT32_MAX;
    uint32         ttLen  = 0;
    double         ttErr  = DBL_MAX;

    int32          ah     = 0;
    int32          bh     = 0;

    uint32         notPresent = 0;  //  Potential parent isn't in the unitig
    uint32         notPlaced  = 0;  //  Potential parent isn't placed yet
    uint32         negHang    = 0;  //  Potential parent has a negative hang to a placed read
    uint32         goodOlap   = 0;

    for (uint32 oo=0; oo<olapsLen; oo++) {

      if (allreads.count(olaps[oo].b_iid) == 0) {
        notPresent++;
        continue;
      }

      if (forward.count(olaps[oo].b_iid) == 0) {       //  Potential parent not placed yet
        notPlaced++;
        continue;
      }

      uint32  l = RI->overlapLength(olaps[oo].a_iid, olaps[oo].b_iid, olaps[oo].a_hang, olaps[oo].b_hang);

      //  Compute the hangs, so we can ignore those that would place this read before the parent.
      //  This is a flaw somewhere in bogart, and should be caught and fixed earlier.

      //  Consensus is expecting the have the hangs for the parent read, not this read, and some
      //  fiddling is needed to flip the overlap for this:
      //    First, swap the reads so it's b-vs-a.
      //    Then, flip the overlap if the b read is in the unitig flipped.

      int32 ah = (olaps[oo].flipped == false) ? (-olaps[oo].a_hang) : (olaps[oo].b_hang);
      int32 bh = (olaps[oo].flipped == false) ? (-olaps[oo].b_hang) : (olaps[oo].a_hang);

      if (forward[olaps[oo].b_iid] == false) {
        swap(ah, bh);
        ah = -ah;
        bh = -bh;
      }

      //  If the ahang is negative, we flubbed up somewhere, and want to place this read before
      //  the parent (even though positions say to place it after, because we sorted by position).

      if (ah < 0) {
        //fprintf(stderr, "ERROR: read %u in tig %u has negative ahang from parent read %u, ejected.\n",
        //        frg->ident, ti, olaps[oo].b_iid);
        negHang++;
        continue;
      }

      //  The overlap is good.  Count it as such.

      goodOlap++;

      //  If the overlap is worse than the one we already have, we don't care.

      if ((l < ttLen) ||                    //  Too short
          (ttErr < olaps[oo].erate())) {    //  Too noisy
        continue;
      }

      tt    = oo;
      ttLen = l;
      ttErr = olaps[oo].erate();
    }

    //  If no thickest overlap, we screwed up somewhere.  Complain and eject the read.

    if (tt == UINT32_MAX) {
      fprintf(stderr, "ERROR: read %u in tig %u has no overlap to any previous read, ejected.  %u overlaps total.  %u negative hang.  %u to read not in tig.  %u to read later in tig.  %u good overlaps.\n",
              frg->ident, tig->id(), olapsLen, negHang, notPresent, notPlaced, goodOlap);
      continue;
    }

    frg->parent = olaps[tt].b_iid;
    frg->ahang  = ah;
    frg->bhang  = bh;



    }  //  Over all reads
  }  //  Over all tigs
}

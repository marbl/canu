
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
 *    src/AS_BAT/AS_BAT_Joining.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010,2012-2013 J. Craig Venter Institute, and
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

#include "AS_BAT_Joining.H"

class joinEntry {
public:
  joinEntry() {
    frFragID = 0;
    toFragID = 0;

    frFirst  = false;

    joinLen  = 0;
  };

  joinEntry(uint32 frFragID_,
            bool   frFirst_,
            uint32 toFragID_,
            bool   toFlip_,
            uint32 joinLen_) {
    frFragID = frFragID_;
    toFragID = toFragID_;

    frFirst  = frFirst_;
    toFlip   = toFlip_;

    joinLen  = joinLen_;
  };

  ~joinEntry() {
  };

  uint32            frFragID;
  uint32            toFragID;

  bool              frFirst;
  bool              toFlip;

  uint32            joinLen;

  bool operator<(const joinEntry &that) const { return(joinLen < that.joinLen); };
  bool operator>(const joinEntry &that) const { return(joinLen > that.joinLen); };
};





//  A cheap test if this unitig looks like it is a bubble in some other unitig.  If both best edges
//  off the end are to the same other unitig, this is probably a bubble.
static
bool
joinUnitigs_looksLikeBubble(Unitig *fr) {
  ufNode          *bgn  = &fr->ufpath[0];
  ufNode          *end  = &fr->ufpath[fr->ufpath.size() - 1];

  BestEdgeOverlap *bestEdge5 = OG->getBestEdgeOverlap(bgn->ident, (bgn->position.end < bgn->position.bgn));
  BestEdgeOverlap *bestEdge3 = OG->getBestEdgeOverlap(end->ident, (end->position.bgn < end->position.end));

  uint32           unitig5 = fr->fragIn(bestEdge5->fragId());
  uint32           unitig3 = fr->fragIn(bestEdge3->fragId());

  if (unitig5 == unitig3) {
    //writeLog("joinUnitigs_looksLikeBubble()-- unitig %u looks like a bubble in unitig %u\n",
    //         fr->id(), unitig5);
    return(true);
  }

  return(false);
}



//  Examine the first (few?) fragments of a unitig, evaluate if they indicate a join should be made.
static
bool
joinUnitigs_examineEnd(UnitigVector      &unitigs,
                       Unitig            *fr,
                       uint32             idx,
                       bool               frFirstEnd,
                       vector<joinEntry> &joins) {
  uint32           frgIdx  = (frFirstEnd) ? (idx) : (fr->ufpath.size() - 1 - idx);
  ufNode          *frg     = &fr->ufpath[frgIdx];
  bool             frgRev  = (frg->position.end < frg->position.bgn);

  //  Grab the best edge for this end frag.  The last arg requests the 3' end if true.
  //
  //  If we're looking at the first read, we want to get:
  //    5' - if the frag is forward
  //    3' - if the frag is reverse (frgRev == true)
  //
  //  If we're looking at the lat read, we want to get:
  //    5' - if the frag is reverse
  //    3' - if the frag is forward  (frgRev == false)
  //
  BestEdgeOverlap *bestEdge    = OG->getBestEdgeOverlap(frg->ident, (frgRev == frFirstEnd));

  uint32      tgtId = bestEdge->fragId();
  bool        tgt3p = bestEdge->frag3p();

  if (tgtId == 0)
    //  No best edge?  Skip it.
    return(false);

  //  Grab the unitig for that best edge.

  uint32   toID  = fr->fragIn(tgtId);
  Unitig  *to    = unitigs[toID];

  if (to->ufpath.size() == 1)
    //  Joining to something teeny?  Don't bother checking further.
    return(false);

  if (to->id() == fr->id())
    //  Join to myself?  Nope.
    return(false);

  //  Grab the read we have an edge to, an compute the overlapping length and left over length.

  ufNode  *tgt    = &to->ufpath[to->pathPosition(tgtId)];
  bool     tgtRev = (tgt->position.end < tgt->position.bgn);

  //  If tgt3p (we overlap to the 3' end) is the same as tgtRev (read is reverse) then the unitig is oriented
  //  correctly.  Otherwise, positions need to be reverse-complemented.


  bool     toFlip = false;

  if ((frFirstEnd == true) && (tgt3p == false) && (tgtRev == false))
    //  source read is at the start, overlap to 5' and the read is forward, need to flip the target unitig
    toFlip = true;

  if ((frFirstEnd == true) && (tgt3p == true) && (tgtRev == true))
    //  source read is at the start, overlap to 3' and the read is reverse, need to flip the target unitig
    toFlip = true;


  if ((frFirstEnd == false) && (tgt3p == false) && (tgtRev == true))
    //  source read is at the end, overlap to 5' and the read is reverse, need to flip the target unitig
    toFlip = true;

  if ((frFirstEnd == false) && (tgt3p == true) && (tgtRev == false))
    //  source read is at the end, overlap to 3' and the read is forward, need to flip the target unitig
    toFlip = true;


  uint32   toMin = MIN(tgt->position.bgn, tgt->position.end);
  uint32   toMax = MAX(tgt->position.bgn, tgt->position.end);
  uint32   toLen = to->getLength();
  uint32   frLen = fr->getLength();

  if (toFlip) {
    toMin = toLen - MAX(tgt->position.bgn, tgt->position.end);
    toMax = toLen - MIN(tgt->position.bgn, tgt->position.end);
  }

  assert(toMin < toMax);

  //  Our two unitigs are of length frLen and toLen.  We are appending some portion of 'to' onto
  //  'fr', and 'discarding' the rest.  If the 'discarded' piece is larger than the 'fr' unitig, we
  //  don't want to do the join.
  //
  //  We err on the side of the discarded piece.

  uint32   joinLen = 0;
  uint32   discLen = 0;

  if (frFirstEnd == true) {
    joinLen = toMin + frLen;  //  Prepend the start of 'to' onto 'fr'.
    discLen = toLen - toMin;

  } else {
    joinLen = frLen + toLen - toMax;  //  Append the end of 'to' onto 'fr'.
    discLen = toMax;
  }

  //  If the discard is bigger than us, we do damage by joining.

  if (discLen > frLen)
    return(false);

  //  The joined should be much larger and the discarded much smaller.

  uint32    maxLen = MAX(frLen, toLen);
  uint32    minLen = MIN(frLen, toLen);

  double    joinChange = (double)joinLen / maxLen;
  double    discChange = (double)discLen / minLen;

  bool      isBad = false;

  if ((joinChange < 1.10) ||
      (0.75       < discChange))
    //  Bad if we didn't really change sizes.
    isBad = true;

  if ((1.0        < joinChange) &&
      (discChange < 0.5))
    //  But good if discard is tiny.  This occurs if we merge a small with a big.  The join change
    //  is somewhat small (1.05 say) yet most of the smaller unitig is used.
    isBad = false;

  if (isBad) {
    writeLog("joinUnitigs_examineEnd()-- join unitig %6u (%7ubp) frag %6u %s <-> unitig %6u (%7ubp) frag %6u %s <-> length %5.2f %7u and %5.2f %7u BAD\n",
             fr->id(), fr->getLength(), frg->ident, (frgRev) ? "rev" : "fwd",
             to->id(), to->getLength(), tgt->ident, (tgtRev) ? "rev" : "fwd",
             joinChange, joinLen,
             discChange, discLen);
    return(false);
  }

  //  OK, join.

  writeLog("joinUnitigs_examineEnd()-- join unitig %6u (%7ubp) frag %6u %s <-> unitig %6u (%7ubp) frag %6u %s <-> length %5.2f %7u and %5.2f %7u\n",
           fr->id(), fr->getLength(), frg->ident, (frgRev) ? "rev" : "fwd",
           to->id(), to->getLength(), tgt->ident, (tgtRev) ? "rev" : "fwd",
           joinChange, joinLen,
           discChange, discLen);

  joins.push_back(joinEntry(frg->ident, frFirstEnd, tgt->ident, toFlip, joinLen));

  return(true);
}









static
void
joinUnitigs_append(UnitigVector &unitigs, joinEntry *join) {
  uint32    frId = Unitig::fragIn(join->frFragID);
  uint32    toId = Unitig::fragIn(join->toFragID);

  Unitig   *fr   = unitigs[frId];
  Unitig   *to   = unitigs[toId];

  uint32    frIdx = Unitig::pathPosition(join->frFragID);
  uint32    toIdx = Unitig::pathPosition(join->toFragID);

  //  The 'fr' unitig is assumed to be forward, and assumed to be the one we join to.

  //  Compute the offset for our append.  We just need to compute where the join fragment would
  //  appear in the unitig.  The join fragment MUST be the first thing in the frUnitig.

  //int32 offset = MIN(frF.position.bgn, frF.position.end);

  //  Over all fragments in the frUnitig, add them to either the joinUnitig or the discUnitig.

  Unitig *joinUnitig = unitigs.newUnitig(false);
  Unitig *discUnitig = unitigs.newUnitig(false);

  //  Reverse the 'to' unitig if needed.

  if (join->toFlip)
    to->reverseComplement(true);

  //  If we're joining off the 5' end of the fr untiig, add the to reads first.

  if (join->frFirst == true) {
    uint32 ii=0;

    for (; ii < toIdx; ii++)
      joinUnitig->addFrag(to->ufpath[ii], 0, false);

    for (; ii < to->ufpath.size(); ii++)
      discUnitig->addFrag(to->ufpath[ii], 0, false);
  }

  //  Now add all the fr unitig reads.

  for (uint32 ii=0; ii < fr->ufpath.size(); ii++)
    joinUnitig->addFrag(to->ufpath[ii], 0, false);

  //  If we're not joining off the 5' end, add the to unitig reads last.

  if (join->frFirst == false) {
    uint32 ii = 0;

    for (; ii < toIdx; ii++)
      discUnitig->addFrag(to->ufpath[ii], 0, false);

    for (; ii < to->ufpath.size(); ii++)
      joinUnitig->addFrag(to->ufpath[ii], 0, false);
  }

  //  Delete the donor unitigs.

  delete fr;
  delete to;

  unitigs[frId] = NULL;
  unitigs[toId] = NULL;

  //  And make sure the new unitigs are consistent.

  joinUnitig->sort();
  discUnitig->sort();
}





void
joinUnitigs(UnitigVector &unitigs, bool enableJoining) {

  if (enableJoining == false)
    return;

  writeLog("==> JOINING SPLIT UNITIGS\n");

  //  Sort unitigs by joined size.  Sort.  Join the largest first.

  vector<joinEntry>  joins;

  //  Over all unitigs, evaluate if a unitig is a candidate for merging onto something.

  for (uint32 frID=0; frID<unitigs.size(); frID++) {
    Unitig        *fr   = unitigs[frID];

    if (fr == NULL)
      //  Ain't no unitig here, mister!
      continue;

    if (fr->ufpath.size() < 2)
      //  Ain't no real unitig here, mister!
      continue;

    //  Do we look like a bubble?

    if (joinUnitigs_looksLikeBubble(fr))
      continue;

    //  The for loop tries reads close to the end - but we don't support joining these.

    for (uint32 ii=0; (ii < 1) && (ii < fr->ufpath.size()); ii++)
      if (joinUnitigs_examineEnd(unitigs, fr, ii, true,  joins))
        break;


    for (uint32 ii=0; (ii < 1) && (ii < fr->ufpath.size()); ii++)
      if (joinUnitigs_examineEnd(unitigs, fr, ii, false, joins))
        break;
  }  //  Over all unitigs.


  writeLog("Found %d pairs of unitigs to join.\n", (int)joins.size());

  std::sort(joins.begin(), joins.end(), greater<joinEntry>());


  return;


  for (uint32 j=0; j<joins.size(); j++) {
    joinEntry  *join = &joins[j];

    //joinUnitigs_append(unitigs, join);
  }
}


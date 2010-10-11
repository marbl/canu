
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static const char *rcsid = "$Id: AS_BOG_IntersectBubble.cc,v 1.7 2010-10-11 03:43:44 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include "MultiAlignStore.h"

#define MAX_OVERLAPS_PER_FRAG   (16 * 1024 * 1024)



//  Sometimes, we don't include fragments into a unitig, even though the best overlap off of both
//  ends is to the same unitig.  This will include those.
//
void
UnitigGraph::popIntersectionBubbles(OverlapStore *ovlStoreUniq, OverlapStore *ovlStoreRept) {
  uint32      ovlMax = MAX_OVERLAPS_PER_FRAG;
  uint32      ovlLen = 0;
  OVSoverlap *ovl    = new OVSoverlap [ovlMax];
  uint32     *ovlCnt = new uint32     [AS_READ_MAX_NORMAL_LEN];

  uint32      nBubblePopped    = 0;
  uint32      nBubbleTooBig    = 0;
  uint32      nBubbleConflict  = 0;
  uint32      nBubbleNoEdge    = 0;

  fprintf(logFile, "==> SEARCHING FOR BUBBLES\n");

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig        *shortTig = unitigs[ti];
    Unitig        *mergeTig = NULL;

    if ((shortTig == NULL) ||
        (shortTig->ufpath.size() >= 30))
      continue;

    uint32         otherUtg     = noUnitig;
    uint32         conflicts    = 0;
    uint32         self         = 0;
    uint32         nonmated     = 0;
    uint32         matedcont    = 0;
    uint32         spurs        = 0;

    uint32         diffOrient   = 0;
    uint32         tooLong      = 0;
    uint32         tigLong      = 0;
    uint32         tigShort     = 0;

    uint32         tooDifferent = 0;

    int32          minNewPos    = INT32_MAX;
    int32          maxNewPos    = INT32_MIN;

    if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
      fprintf(logFile, "popBubbles()-- try unitig %d of length %d with %d fragments\n",
              shortTig->id(),
              shortTig->getLength(),
              shortTig->ufpath.size());

    for (uint32 fi=0; fi<shortTig->ufpath.size(); fi++) {
      ufNode *frg = &shortTig->ufpath[fi];

      int32  frgID = frg->ident;
      int32  utgID = shortTig->id();

      BestEdgeOverlap *bestedge5 = OG->getBestEdgeOverlap(frg->ident, false);
      BestEdgeOverlap *bestedge3 = OG->getBestEdgeOverlap(frg->ident, true);
      BestContainment *bestcont  = OG->getBestContainer(frg->ident);

      if (FI->mateIID(frgID) > 0) {
        if (bestcont)
          matedcont++;
      } else {
        nonmated++;
      }

      if (bestcont)
        continue;

      if (bestedge5->fragId() == 0) {
        spurs++;
      } else {
        uint32 ou5 = shortTig->fragIn(bestedge5->fragId());

        assert(ou5 > 0);

        if (ou5 == shortTig->id()) {
          self++;
        } else {
          if ((otherUtg == noUnitig) && (ou5 != 0))
            otherUtg = ou5;
          if (otherUtg != ou5)
            conflicts++;
        }
      }

      if (bestedge3->fragId() == 0) {
        spurs++;
      } else {
        uint32 ou3 = shortTig->fragIn(bestedge3->fragId());

        assert(ou3 > 0);

        if (ou3 == shortTig->id()) {
          self++;
        } else {
          if ((otherUtg == noUnitig) && (ou3 != 0))
            otherUtg = ou3;
          if (otherUtg != ou3)
            conflicts++;
        }
      }

      if (otherUtg == noUnitig)
        //  Didn't find a unitig to merge this fragment into.
        continue;

      if (conflicts > 0)
        //  Found multiple unitigs to merge this fragment into, or multiple unitigs to merge this unitig into.
        continue;

      //  Check sanity of the new placement.

      mergeTig = unitigs[otherUtg];

      ufNode place5;
      ufNode place3;

      place5.ident = frgID;
      place3.ident = frgID;

      int32 bidx5 = -1;
      int32 bidx3 = -1;

      mergeTig->placeFrag(place5, bidx5, bestedge5,
                          place3, bidx3, bestedge3);

      //  Either or both ends can fail to place in the new unitig -- edges can be to a fragment in
      //  the tig we're testing.  This doesn't matter for finding the min/max, but does matter
      //  when we decide if the placements are about the correct size.

      //  Update the min/max placement for the whole tig.  At the end of everyhing we'll
      //  make sure the min/max agree with the size of the tig.

      int32  min5 = INT32_MAX, max5 = INT32_MIN;
      int32  min3 = INT32_MAX, max3 = INT32_MIN;

      int32  minU = INT32_MAX;
      int32  maxU = INT32_MIN;

      if (bidx5 != -1) {
        if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
          fprintf(logFile, "popBubbles()-- place frag %d using 5' edge (%d,%c) in unitig %d at %d,%d\n",
                  frgID, bestedge5->fragId(), bestedge5->frag3p() ? '3' : '5', otherUtg, place5.position.bgn, place5.position.end);
        min5 = MIN(place5.position.bgn, place5.position.end);
        max5 = MAX(place5.position.bgn, place5.position.end);
      }

      if (bidx3 != -1) {
        if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
          fprintf(logFile, "popBubbles()-- place frag %d using 3' edge (%d,%c) in unitig %d at %d,%d\n",
                  frgID, bestedge3->fragId(), bestedge3->frag3p() ? '3' : '5', otherUtg, place3.position.bgn, place3.position.end);
        min3 = MIN(place3.position.bgn, place3.position.end);
        max3 = MAX(place3.position.bgn, place3.position.end);
      }

      minU = MIN(min5, min3);
      maxU = MAX(max5, max3);

      minNewPos = MIN(minU, minNewPos);
      maxNewPos = MAX(maxU, maxNewPos);

      //  Check that the two placements agree with each other.  Same orientation, more or less the
      //  same location, more or less the correct length.  For the location, we only need to test
      //  that the min and max positions are roughly the same as the fragment length.  If the two
      //  edges place the fragment in different locations, then the min/max values of the placement
      //  will be too large (never too small).

      if ((bidx5 != -1) && (bidx3 != -1)) {
        if ((min5 < max5) != (min3 < max3))
          //  Orientation bad.
          diffOrient++;

        if (maxU - minU > 1.25 * FI->fragmentLength(frgID)) {
          //  Location bad.
          if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
            fprintf(logFile, "popBubbles()--   too long1 %d - %d = %d > 1.06 * %d = %.2f\n",
                    maxU, minU, maxU - minU, FI->fragmentLength(frgID), 1.06 * FI->fragmentLength(frgID));
          tooLong++;
        }

        if (max5 - min5 > 1.25 * FI->fragmentLength(frgID)) {
          //  Length bad.
          if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
            fprintf(logFile, "popBubbles()--   too long2 %d - %d = %d > 1.06 * %d = %.2f\n",
                    max5, min5, max5 - min5, FI->fragmentLength(frgID), 1.06 * FI->fragmentLength(frgID));
          tooLong++;
        }

        if (max3 - min3 > 1.25 * FI->fragmentLength(frgID)) {
          //  Length bad.
          if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
            fprintf(logFile, "popBubbles()--   too long3 %d - %d = %d > 1.06 * %d = %.2f\n",
                    max3, min3, max3 - min3, FI->fragmentLength(frgID), 1.06 * FI->fragmentLength(frgID));
          tooLong++;
        }
      }
    }  //  Over all fragments in the source unitig/

    //
    //  If we are bad already, just stop.  If we pass these tests, continue on to checking overlaps.
    //

    if (otherUtg == noUnitig) {
      if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
        fprintf(logFile, "popBubbles()-- unitig %d has NO EDGES\n",
                shortTig->id());
      nBubbleNoEdge++;
      continue;
    }

    if (maxNewPos - minNewPos > 1.25 * shortTig->getLength()) {
      //  Bad placement; edges indicate we blew the unitig apart.
      if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
        fprintf(logFile, "popBubbles()--   tig long %d - %d = %d > 1.06 * %d = %.2f\n",
                maxNewPos, minNewPos, maxNewPos - minNewPos, shortTig->getLength(), 1.06 * shortTig->getLength());
      tigLong++;
    }

    if (maxNewPos - minNewPos < 0.75 * shortTig->getLength()) {
      //  Bad placement; edges indicate we compressed the unitig (usually by placing only one fragment)
      if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
        fprintf(logFile, "popBubbles()--   tig short %d - %d = %d < 1.06 * %d = %.2f\n",
                maxNewPos, minNewPos, maxNewPos - minNewPos, shortTig->getLength(), 1.06 * shortTig->getLength());
      tigShort++;
    }

    if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
      fprintf(logFile, "popBubbles()-- unitig %d CONFLICTS %d SPURS %d SELF %d len %d frags %u matedcont %d nonmated %d diffOrient %d tooLong %d tigLong %d tigShort %d\n",
              shortTig->id(), conflicts, spurs, self, shortTig->getLength(), (uint32)shortTig->ufpath.size(), matedcont, nonmated, diffOrient, tooLong, tigLong, tigShort);

    if ((spurs        > 0) ||
#if 0
        (self         > 6) ||  //  These are possible too aggressive.  They were originally
        (matedcont    > 6) ||  //  used before CHECK_OVERLAPS existed.
#endif
        (diffOrient   > 0) ||
        (tooLong      > 0) ||
        (tigLong      > 0) ||
        (tigShort     > 0) ||
        (tooDifferent > 0) ||
        (conflicts    > 0)) {
      nBubbleConflict++;
      continue;
    }

    if (mergeTig == NULL)
      //  Didn't find any place to put this short unitig.  It's not a bubble, just a short unitig
      //  with no overlaps anywhere.
      continue;

    //
    //  Grab the overlaps for this read, paint the number of times we overlap some fragment
    //  already in this unitig.  If we have any significant blocks with no coverage, the bubble
    //  is probably too big for consensus.
    //

#define CHECK_OVERLAPS
#ifdef CHECK_OVERLAPS
    for (uint32 fi=0; fi<shortTig->ufpath.size(); fi++) {
      ufNode *frg = &shortTig->ufpath[fi];

      int32  frgID = frg->ident;
      int32  utgID = shortTig->id();

      BestEdgeOverlap *bestedge5 = OG->getBestEdgeOverlap(frg->ident, false);
      BestEdgeOverlap *bestedge3 = OG->getBestEdgeOverlap(frg->ident, true);
      BestContainment *bestcont  = OG->getBestContainer(frg->ident);

      if (bestcont)
        continue;

      ovlLen  = 0;

      if (ovlStoreUniq) {
        AS_OVS_setRangeOverlapStore(ovlStoreUniq, frgID, frgID);
        ovlLen += AS_OVS_readOverlapsFromStore(ovlStoreUniq, ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
      }

      if (ovlStoreRept) {
        AS_OVS_setRangeOverlapStore(ovlStoreRept, frgID, frgID);
        ovlLen += AS_OVS_readOverlapsFromStore(ovlStoreRept, ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
      }

      memset(ovlCnt, 0, sizeof(uint32) * AS_READ_MAX_NORMAL_LEN);

      uint32 alen = FI->fragmentLength(frgID);

      for (uint32 o=0; o<ovlLen; o++) {
        int32 a_iid = ovl[o].a_iid;
        int32 b_iid = ovl[o].b_iid;

        assert(a_iid == frgID);

        if (shortTig->fragIn(b_iid) != mergeTig->id())
          //  Ignore overlaps to fragments not in this unitig.  We might want to ignore overlaps
          //  to fragments in this unitig but not at this position, though that assumes we
          //  actually got the placement correct....and it only matters for repeats.
          continue;

        uint32        mrgPos = mergeTig->pathPosition(b_iid);
        ufNode *mrg    = &mergeTig->ufpath[mrgPos];

        if ((mrg->position.bgn < minNewPos - AS_OVERLAP_MIN_LEN) &&
            (mrg->position.end < minNewPos - AS_OVERLAP_MIN_LEN)) {
          //  This overlapping fragment is before the position we are supposed to be merging to, skip it.
          if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
            fprintf(logFile, "frag %d ignores overlap to frag %d - outside range\n", frgID, b_iid);
          continue;
        }

        if ((mrg->position.bgn > maxNewPos - AS_OVERLAP_MIN_LEN) &&
            (mrg->position.end > maxNewPos - AS_OVERLAP_MIN_LEN)) {
          //  This overlapping fragment is after the position we are supposed to be merging to, skip it.
          if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
            fprintf(logFile, "frag %d ignores overlap to frag %d - outside range\n", frgID, b_iid);
          continue;
        }

        uint32 bgn = 0;
        uint32 end = 0;

        int32 a_hang = ovl[o].dat.ovl.a_hang;
        int32 b_hang = ovl[o].dat.ovl.b_hang;

        if (a_hang < 0) {
          //  b_hang < 0      ?     ----------  :     ----
          //                  ?  ----------     :  ----------
          //
          bgn = 0;
          end = (b_hang < 0) ? (alen + b_hang) : (alen);
        } else {
          //  b_hang < 0      ?  ----------              :  ----------
          //                  ?     ----                 :     ----------
          //
          bgn = a_hang;
          end = (b_hang < 0) ? (alen + b_hang) : alen;
        }

        for (uint32 x=bgn; x<end; x++)
          ovlCnt[x]++;
      }

      //  Score the overlap coverage.  Allow a small amount of total zero, but don't
      //  allow any significant blocks of zero.

      uint32 zTotal    = 0;
      uint32 zInternal = 0;
      uint32 zSum      = 0;

      for (uint32 x=0; x<alen; x++) {
        if (ovlCnt[x] == 0) {
          zTotal++;
          zSum++;
        } else {
          if (zSum > zInternal)
            zInternal = zSum;
          zSum = 0;
        }
      }

      if ((zTotal > 0.10 * alen) ||
          (zInternal > AS_OVERLAP_MIN_LEN)) {
        if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG)) {
          fprintf(logFile, "frag %d too different with zTotal=%d (limit=%d) and zInternal=%d (ovl=%d)\n",
                  frgID, zTotal, (int)(0.10 * alen), zInternal, AS_OVERLAP_MIN_LEN);
          for (uint32 x=0; x<alen; x++)
            fprintf(logFile, "%c", (ovlCnt[x] < 10) ? '0' + ovlCnt[x] : '*');
          fprintf(logFile, "\n");
        }
        tooDifferent++;
      }
    }  //  over all frags

    //  If there are fragments with missing overlaps, don't merge.

    if (logFileFlagSet(LOG_INTERSECTION_BUBBLES_DEBUG))
      fprintf(logFile, "popBubbles()-- unitig %d tooDifferent %d\n",
              shortTig->id(), tooDifferent);

    if (tooDifferent > 0) {
      nBubbleTooBig++;
      continue;
    }
#endif // CHECK_OVERLAPS


    //  Merge this unitig into otherUtg.

    if (logFileFlagSet(LOG_INTERSECTION_BUBBLES))
      fprintf(logFile, "popBubbles()-- merge unitig %d (len %d) into unitig %d (len %d)\n",
              shortTig->id(), shortTig->getLength(),
              mergeTig->id(), mergeTig->getLength());
    nBubblePopped++;

    //  Every once in a while, we get a unitig that is misordered, caused by conflicting best overlaps.
    //  For example:
    //
    //  A          <---------                    <----------
    //  B                ------------>                  ----------->
    //  C    <--------                                          <------------
    //
    //  Where the C fragment has edges off both A and B.  If the layout is the one on the right,
    //  when we go to move the A fragment to mergeTig, it has no edges to any fragment there (since
    //  the picture on the left shows it has edges to B and C), and so we cannot place A until
    //  either B or C is placed.

    bool  tryAgain  = true;
    bool  allPlaced = true;
    bool  isStuck   = false;

    while (tryAgain) {
      tryAgain  = false;
      allPlaced = true;
      isStuck   = true;

      for (uint32 fi=0; fi<shortTig->ufpath.size(); fi++) {
        ufNode  *frag = &shortTig->ufpath[fi];

        if (mergeTig->fragIn(frag->ident) == mergeTig->id())
          //  Already placed in mergeTig.
          continue;

        allPlaced = false;

        if (OG->getBestContainer(frag->ident))
          mergeTig->addContainedFrag(frag->ident,
                                     OG->getBestContainer(frag->ident),
                                     logFileFlagSet(LOG_INTERSECTION_BUBBLES));
        else
          mergeTig->addAndPlaceFrag(frag->ident,
                                    OG->getBestEdgeOverlap(frag->ident, false),
                                    OG->getBestEdgeOverlap(frag->ident, true),
                                    logFileFlagSet(LOG_INTERSECTION_BUBBLES));

        if (mergeTig->fragIn(frag->ident) == mergeTig->id()) {
          //  Placed something, making progress!
          if (logFileFlagSet(LOG_INTERSECTION_BUBBLES))
            fprintf(logFile, "popBubbles()-- Moved frag %d from unitig %d to unitig %d (isStuck <- false)\n",
                    frag->ident, shortTig->id(), mergeTig->id());
          isStuck = false;
        } else {
          //  Failed to place, gotta do the loop again.
          if (logFileFlagSet(LOG_INTERSECTION_BUBBLES))
            fprintf(logFile, "popBubbles()-- Failed to move frag %d from unitig %d to unitig %d (tryAgain <- true)\n",
                    frag->ident, shortTig->id(), mergeTig->id());
          tryAgain = true;
        }
      }

      if ((allPlaced == false) && (isStuck == true)) {
        if (logFileFlagSet(LOG_INTERSECTION_BUBBLES))
          fprintf(logFile, "popBubbles()--  Failed to completely merge unitig %d into unitig %d.\n",
                  shortTig->id(), mergeTig->id());
        tryAgain = false;
      }
    }

    mergeTig->sort();

    //  Mark this unitig as dead if we fully merged it.

    if (isStuck == false) {
      unitigs[ti] = NULL;
      delete shortTig;
    }
  }  //  over all unitigs

  delete [] ovl;
  delete [] ovlCnt;

  fprintf(logFile, "==> SEARCHING FOR BUBBLES done, %u popped, %u had conflicting placement, %u had no edges, %u were too dissimilar.\n",
          nBubblePopped, nBubbleConflict, nBubbleNoEdge, nBubbleTooBig);
}  //  bubble popping scope


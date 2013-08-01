
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

static const char *rcsid = "$Id: AS_BOG_MoveContains.cc,v 1.1 2010-10-01 13:43:49 brianwalenz Exp $";

#include "AS_BOG_BestOverlapGraph.H"
#include "AS_BOG_UnitigGraph.H"
#include "AS_BOG_MateLocation.H"


//  Make sure that contained fragments are in the same unitig
//  as their container.  Due to sorting, contained fragments
//  can come much later in the unitig:
//
//  ------------1
//    -------------2
//       --------------3
//         ----4 (contained in 1, too much error keeps it out of 2 and 3)
//
//  So, our first pass is to move contained fragments around.
//
void UnitigGraph::moveContains(void) {

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *thisUnitig = unitigs[ti];

    if ((thisUnitig == NULL) ||
        (thisUnitig->ufpath.size() < 2))
      continue;

    MateLocation positions(thisUnitig);

    ufNode               *frags         = new ufNode [thisUnitig->ufpath.size()];
    uint32                fragsLen      = 0;

    if (logFileFlagSet(LOG_MATE_SPLIT_UNHAPPY_CONTAINS))
      fprintf(logFile, "moveContain unitig %d\n", thisUnitig->id());

    for (uint32 fi=0; fi<thisUnitig->ufpath.size(); fi++) {
      ufNode  *frg = &thisUnitig->ufpath[fi];

      BestContainment   *bestcont   = OG->getBestContainer(frg->ident);
      MateLocationEntry  mloc       = positions.getById(frg->ident);

      uint32  thisFrgID = frg->ident;
      uint32  contFrgID = (bestcont) ? bestcont->container : 0;
      uint32  mateFrgID = FI->mateIID(frg->ident);

      uint32  thisUtgID = thisUnitig->fragIn(thisFrgID);
      uint32  contUtgID = thisUnitig->fragIn(contFrgID);
      uint32  mateUtgID = thisUnitig->fragIn(mateFrgID);

      //  id1 != 0 -> we found the fragment in the mate happiness table
      //  isBad -> and the mate is unhappy.
      //
      //  What's id1 vs id2 in MateLocationEntry?  Dunno.  All I
      //  know is that if there is no mate present, one of those
      //  will be 0.  (Similar test used above too.)
      //
      bool    isMated    = (mateFrgID > 0);
      bool    isGrumpy   = ((isMated) && (mloc.mleFrgID1 != 0) && (mloc.mleFrgID2 != 0) && (mloc.isGrumpy == true));

      //
      //  Figure out what to do.
      //

      bool    moveToContainer = false;
      bool    moveToSingleton = false;

      if        ((frg->contained == 0) && (bestcont == NULL)) {
        //  CASE 1:  Not contained.  Leave the fragment here.
        //fprintf(logFile, "case1 frag %d fragsLen %d\n", thisFrgID, fragsLen);

      } else if (isMated == false) {
        //  CASE 2: Contained but not mated.  Move to be with the
        //  container (if the container isn't here).
        //fprintf(logFile, "case2 frag %d contID %d fragsLen %d\n", thisFrgID, contUtgID, fragsLen);

        if (thisUtgID != contUtgID)
          moveToContainer = true;

      } else if ((isGrumpy == true) && (thisUtgID == mateUtgID)) {
        //  CASE 3: Not happy, and the frag and mate are together.
        //  Kick out to a singleton.

        //fprintf(logFile, "case3 frag %d utg %d mate %d utg %d cont %d utg %d fragsLen %d\n",
        //        thisFrgID, thisUtgID, mateFrgID, mateUtgID, contFrgID, contUtgID, fragsLen);

        if (thisUtgID == mateUtgID)
          moveToSingleton = true;

      } else {

        //  This makes for some ugly code (we break the nice if else
        //  if else structure we had going on) but the next two cases
        //  need to know if there is an overlap to the rest of the
        //  unitig.

        bool  hasOverlap   = (thisUtgID == contUtgID);
        bool  allContained = false;


        if (hasOverlap == false) {
          if (fragsLen == 0) {
            //  The first fragment.  Check fragments after to see if
            //  there is an overlap (note only frags with an overlap
            //  in the layout are tested).  In rare cases, we ejected
            //  the container, and left a containee with no overlap to
            //  fragments remaining.
            //
            //  Note that this checks if there is an overlap to the
            //  very first non-contained (aka dovetail) fragment ONLY.
            //  If there isn't an overlap to the first non-contained
            //  fragment, then that fragment will likely NOT align
            //  correctly.

            uint32 ft = fi + 1;

#warning 2x BUGS IN COMPARISON HERE

            //  Skip all the contains.
            while ((ft < thisUnitig->ufpath.size()) &&
                   (OG->isContained(thisUnitig->ufpath[ft].ident) == true) &&
                   (MAX(frg->position.bgn, frg->position.end) < MIN(thisUnitig->ufpath[ft].position.bgn, thisUnitig->ufpath[ft].position.end)))
              ft++;

            //  If the frag is not contained (we could be the
            //  container), and overlaps in the layout, see if there
            //  is a real overlap.
            if ((ft < thisUnitig->ufpath.size()) &&
                (OG->isContained(thisUnitig->ufpath[ft].ident) == false) &&
                (MAX(frg->position.bgn, frg->position.end) < MIN(thisUnitig->ufpath[ft].position.bgn, thisUnitig->ufpath[ft].position.end)))
              hasOverlap = OG->containHaveEdgeTo(thisFrgID, thisUnitig->ufpath[ft].ident);
          } else {
            //  Not the first fragment, search for an overlap to an
            //  already placed frag.

            uint32  ft = fi;

            do {
              ft--;

              //  OK to overlap to a contained frag; he could be our
              //  container.

              hasOverlap = OG->containHaveEdgeTo(thisFrgID, thisUnitig->ufpath[ft].ident);

              //  Stop if we found an overlap, or we just checked the
              //  first frag in the unitig, or we no longer overlap in
              //  the layout.
            } while ((hasOverlap == false) &&
                     (ft > 0) &&
                     (MIN(frg->position.bgn, frg->position.end) < MAX(thisUnitig->ufpath[ft].position.bgn, thisUnitig->ufpath[ft].position.end)));
          }
        }  //  end of hasOverlap


        //  An unbelievabe special case.  When the unitig is just a
        //  single container fragment (and any contained frags under
        //  it) rule 4 breaks.  The first fragment has no overlap (all
        //  later reads are contained) and so we want to eject it to a
        //  new unitig.  Since there are multiple fragments in this
        //  unitig, the ejection occurs.  Later, all the contains get
        //  moved to the new unitig.  And we repeat.  To prevent, we
        //  abort the ejection if the unitig is all contained in one
        //  fragment.
        //
        if (fragsLen == 0) {
          allContained = true;

          for (uint32 ft = fi + 1; ((allContained == true) && (ft < thisUnitig->ufpath.size())); ft++)
            allContained = OG->isContained(thisUnitig->ufpath[ft].ident);
        }



        if (isGrumpy == true) {
          //  CASE 4: Not happy and not with the mate.  This one is a
          //  bit of a decision.
          //
          //  If an overlap exists to the rest of the unitig, we'll
          //  leave it here.  We'll also leave it here if it is the
          //  rest of the unitig is all contained in this fragment.
          //
          //  If no overlap, and the mate and container are in the
          //  same unitig, we'll just eject.  That also implies the
          //  other unitig is somewhat large, at least as big as the
          //  insert size.
          //
          //  Otherwise, we'll move to the container and cross our
          //  fingers we place it correctly.  The alternative is to
          //  eject, and hope that we didn't also eject the mate to a
          //  singleton.

          //fprintf(logFile, "case4 frag %d utg %d mate %d utg %d cont %d utg %d fragsLen %d\n",
          //        thisFrgID, thisUtgID, mateFrgID, mateUtgID, contFrgID, contUtgID, fragsLen);

          if ((hasOverlap == false) && (allContained == false))
            if (mateUtgID == contUtgID)
              moveToSingleton = true;
            else
              moveToContainer = true;

        } else {
          //  CASE 5: Happy!  If with container, or an overlap exists to
          //  some earlier fragment, leave it here.  Otherwise, eject it
          //  to a singleton.  The fragment is ejected instead of moved
          //  to be with its container since we don't know which is
          //  correct - the mate or the overlap.
          //
          //  If not happy, we've already made sure that the mate is not
          //  here (that was case 3).

          //fprintf(logFile, "case5 frag %d utg %d mate %d utg %d cont %d utg %d fragsLen %d\n",
          //        thisFrgID, thisUtgID, mateFrgID, mateUtgID, contFrgID, contUtgID, fragsLen);

          //  If no overlap (so not with container or no overlap to
          //  other frags) eject.
          if ((hasOverlap == false) && (allContained == false))
            moveToSingleton = true;
        }
      }  //  End of cases

      //
      //  Do it.
      //

      if (moveToContainer == true) {
        //  Move the fragment to be with its container.

        Unitig         *thatUnitig = unitigs[contUtgID];
        ufNode          containee  = *frg;

        assert(thatUnitig->id() == contUtgID);

        //  Nuke the fragment in the current list
        frg->ident        = 999999999;
        frg->contained    = 999999999;
        frg->position.bgn = 0;
        frg->position.end = 0;

        assert(thatUnitig->id() == contUtgID);

        if (logFileFlagSet(LOG_MATE_SPLIT_UNHAPPY_CONTAINS))
          fprintf(logFile, "Moving contained fragment %d from unitig %d to be with its container %d in unitig %d\n",
                  thisFrgID, thisUtgID, contFrgID, contUtgID);

        assert(bestcont->container == contFrgID);

        thatUnitig->addContainedFrag(thisFrgID, bestcont, logFileFlagSet(LOG_MATE_SPLIT_UNHAPPY_CONTAINS));
        assert(thatUnitig->id() == Unitig::fragIn(thisFrgID));

      } else if ((moveToSingleton == true) && (thisUnitig->getNumFrags() != 1)) {
        //  Eject the fragment to a singleton (unless we ARE the singleton)
        Unitig        *singUnitig  = new Unitig(logFileFlagSet(LOG_MATE_SPLIT_UNHAPPY_CONTAINS));
        ufNode         containee  = *frg;

        //  Nuke the fragment in the current list
        frg->ident        = 999999999;
        frg->contained    = 999999999;
        frg->position.bgn = 0;
        frg->position.end = 0;

        if (logFileFlagSet(LOG_MATE_SPLIT_UNHAPPY_CONTAINS))
          fprintf(logFile, "Ejecting unhappy contained fragment %d from unitig %d into new unitig %d\n",
                  thisFrgID, thisUtgID, singUnitig->id());

        containee.contained = 0;

        singUnitig->addFrag(containee, -MIN(containee.position.bgn, containee.position.end), logFileFlagSet(LOG_MATE_SPLIT_UNHAPPY_CONTAINS));

        unitigs.push_back(singUnitig);
        thisUnitig = unitigs[ti];  //  Reset the pointer; unitigs might be reallocated

      } else {
        //  Leave fragment here.  Copy the fragment to the list -- if
        //  we need to rebuild the unitig (because fragments were
        //  removed), the list is used, otherwise, we have already
        //  made the changes needed.
        //
        //  Also, very important, update our containment mark.  If our
        //  container was moved, but we stayed put because of a happy
        //  mate, we're still marked as being contained.  Rather than
        //  put this check in all the places where we stay put in the
        //  above if-else-else-else, it's here.

        if ((frg->contained) && (thisUtgID != contUtgID))
          frg->contained = 0;

        frags[fragsLen] = *frg;
        fragsLen++;
      }

    }  //  over all frags

    //  Now, rebuild this unitig if we made changes.

    if (fragsLen != thisUnitig->ufpath.size()) {
      if (logFileFlagSet(LOG_MATE_SPLIT_UNHAPPY_CONTAINS))
        fprintf(logFile, "Rebuild unitig %d after removing contained fragments.\n", thisUnitig->id());

      thisUnitig->ufpath.clear();

      //  Occasionally, we move all fragments out of the original unitig.  Might be worth checking
      //  if that makes sense!!
      //
#warning EMPTIED OUT A UNITIG
      if (fragsLen > 0) {
        //  No need to resort.  Offsets only need adjustment if the first fragment is thrown out.
        //  If not, splitOffset will be zero.
        //
        int splitOffset = -MIN(frags[0].position.bgn, frags[0].position.end);

        //  This is where we clean up from the splitting not dealing with contained fragments -- we
        //  force the first frag to be uncontained.
        //
        frags[0].contained = 0;

        for (uint32 i=0; i<fragsLen; i++)
          thisUnitig->addFrag(frags[i], splitOffset, logFileFlagSet(LOG_MATE_SPLIT_UNHAPPY_CONTAINS));
      }
    }

    delete [] frags;
    frags = NULL;

  }  //  Over all unitigs
}


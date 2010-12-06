
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

static const char *rcsid = "$Id: AS_BAT_Breaking.C,v 1.2 2010-12-06 08:03:48 brianwalenz Exp $";

#include "AS_BAT_Breaking.H"

#include "AS_BAT_BestOverlapGraph.H"

//  The original version was filtering breakpoints.  It was accepting any break point with more than
//  MIN_BREAK_FRAGS fragments and longer than MIN_BREAK_LENGTH.  The shorter ones in between two
//  large break points were (I suspect) analyzed to see if many short break points were piling up in
//  one region.  If so, one was selected and accepted into the list of final break points.

static const int MIN_BREAK_FRAGS = 1;
static const int MIN_BREAK_LENGTH = 500;


void
filterBreakPoints(Unitig *tig,
                               UnitigBreakPoints &breaks) {

  UnitigBreakPoints  newBreaks;

  for (uint32 i=0; i<breaks.size(); i++)
    if ((breaks[i].inFrags > MIN_BREAK_FRAGS) &&
        (breaks[i].inSize  > MIN_BREAK_LENGTH))
      newBreaks.push_back(breaks[i]);

  breaks = newBreaks;
}




UnitigVector *
breakUnitigAt(Unitig *tig,
                           UnitigBreakPoints &breaks) {

  if (breaks.empty())
    return NULL;

  //  Contained fragments provide a huge headache when splitting a unitig.  We just don't know where
  //  to put them:
  //
  //  ------------------------------>
  //            <------------------------------------
  //                  --------->
  //
  //  If we split at either point, we cannot easily tell if we should keep the contain with the
  //  forward or with the reverse fragment.
  //
  //  To resolve this, we simply do not place contained fragments in a new unitig (and we must
  //  explicitly remove them from the original unitig).  After all splits are done, we'll go back
  //  and place the contains again.

  UnitigVector *newTigs = new UnitigVector();
  Unitig       *newTig  = NULL;
  int32         offset  = 0;
  uint32        bpCur   = 0;  //  Current break point in the breaks list

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode  frg         = tig->ufpath[fi];
    bool    frgReversed = (frg.position.end < frg.position.bgn);

    bool    break5p     = false;
    bool    break3p     = false;

    if (OG->isContained(frg.ident)) {
      Unitig::removeFrag(frg.ident);
      continue;
    }

    //  Current fragment is the one we want to break on.  Figure out which ends to break on.

    if (frg.ident == breaks[bpCur].fragEnd.fragId()) {
      break3p   = (breaks[bpCur].fragEnd.frag3p() == true);
      break5p   = (breaks[bpCur].fragEnd.frag3p() == false);

      bpCur++;

      while (((bpCur < breaks.size()) &&
              (frg.ident == breaks[bpCur].fragEnd.fragId()))) {
        break3p |= (breaks[bpCur].fragEnd.frag3p() == true);
        break5p |= (breaks[bpCur].fragEnd.frag3p() == false);
        bpCur++;
      }
    }

    //  Move this fragment to a new unitig, possibly creating a new one for it, and maybe even
    //  forgetting about the new unitig after the fragment is added.

    if        (break5p && break3p) {
      //  Break at both ends, this fragment ends up as a singleton.

      newTig = new Unitig(logFileFlagSet(LOG_INTERSECTION_BREAKING));  //  always make a new unitig, we put a frag in it now
      newTigs->push_back(newTig);
      offset = (frgReversed) ? -frg.position.end : -frg.position.bgn;

      newTig->addFrag(frg, offset, logFileFlagSet(LOG_INTERSECTION_BREAKING));

      newTig = NULL;


    } else if ((break5p && (frgReversed == false)) ||
               (break3p && (frgReversed == true))) {
      //  Break at left end, this fragment starts a new unitig.

      newTig = new Unitig(logFileFlagSet(LOG_INTERSECTION_BREAKING));  //  always make a new unitig, we put a frag in it now
      newTigs->push_back(newTig);
      offset = (frgReversed) ? -frg.position.end : -frg.position.bgn;

      newTig->addFrag(frg, offset, logFileFlagSet(LOG_INTERSECTION_BREAKING));


    } else if ((break5p && (frgReversed == true)) ||
               (break3p && (frgReversed == false))) {
      //  Break at right end, this fragment ends the existing unitig (which might not even exist).

      if (newTig == NULL) {
        newTig = new Unitig(logFileFlagSet(LOG_INTERSECTION_BREAKING));
        newTigs->push_back(newTig);
        offset = (frgReversed) ? -frg.position.end : -frg.position.bgn;
      }

      newTig->addFrag(frg, offset, logFileFlagSet(LOG_INTERSECTION_BREAKING));

      newTig = NULL;


    } else {
      //  Don't break, just add to the existing unitig (which might not even exist).

      if (newTig == NULL) {
        newTig = new Unitig(logFileFlagSet(LOG_INTERSECTION_BREAKING));
        newTigs->push_back(newTig);
        offset = (frgReversed) ? -frg.position.end : -frg.position.bgn;
      }

      newTig->addFrag(frg, offset, logFileFlagSet(LOG_INTERSECTION_BREAKING));
    }
  }


  //  Finally, add back in contains.


  return(newTigs);
}

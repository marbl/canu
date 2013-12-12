
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007-2008, J. Craig Venter Institute
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

static const char *rcsid = "$Id$";

#include "AS_global.H"
#include "AS_ALN_bruteforcedp.H"
#include "AS_UTL_reverseComplement.H"

#undef DEBUG

#define STOP             0  //  IMPORTANT; default is STOP
#define MATCH            1
#define GAPA             2
#define GAPB             3

//  3, -2, -2 original
//  5, -4, -3 suggested
//  3, -6, -4 got around one problem, but made many more
//
#define MATCHSCORE       3
#define GAPSCORE        -6
#define MISMATCHSCORE   -4

#define SLOP             10

//  Scores in the dynamic programming matrix are unsigned ints, currently 30 bits wide.
//
//  The smallest score is 0.
//  The zero score is 2^29 (512 million)
//  The largest score is 2^30-1 (1 billion).
//
//  Edges that we don't want to end up at are initialized to DP_NEGT (2^28, 256 million) which will
//  take an awful lot of alignment to make it decent, or to underflow.
//
//  Edges that we want to end up at are set to DP_ZERO (2^29, 512 million).

#define DP_NEGT          (1 << 28)  //  A very negative value
#define DP_ZERO          (1 << 29)  //  Zero


//
//  ahang, bhang represent any sequence to EXCLUDE from the alignmet.  There is a little bit of slop
//  in this exclusion.
//
//  Setting both to the length of the sequence will try to find an alignment over all bases.  The
//  highest scoring alignment is returned, which is likely not an alignment that uses all bases --
//  the only requirement (assuming endToEnd is set) is that the alignment reaches the end of one
//  sequence.
//

void
alignLinker(char           *alignA,
            char           *alignB,
            char           *stringA,
            char           *stringB,
            alignLinker_s  *a,
            int             endToEnd,
            int             allowNs,
            int             ahang,
            int             bhang) {

  int32 lenA = strlen(stringA);
  int32 lenB = strlen(stringB);

  memset(a, 0, sizeof(alignLinker_s));

  dpActions  ACT;

  uint32     SCO[4][AS_READ_MAX_NORMAL_LEN * 2];

  memset(SCO, 0, sizeof(uint32) * AS_READ_MAX_NORMAL_LEN * 2);

  uint32    *lastCol = SCO[0];
  uint32    *thisCol = SCO[1];
  uint32    *iFinal  = SCO[2];
  uint32    *jFinal  = SCO[3];

  if ((lenA > AS_READ_MAX_NORMAL_LEN) || (lenB > AS_READ_MAX_NORMAL_LEN)) {
    fprintf(stderr, "alignLinker()-- Reads too long.  %d or %d > %d\n", lenA, lenB, AS_READ_MAX_NORMAL_LEN);
    return;
  }

  //  Definition of the box we want to do dynamic programming in.
  int32 ibgn = 1;
  int32 iend = lenA;
  int32 jbgn = 1;
  int32 jend = lenB;

  //  Set the edges.

  for (int32 i=ibgn-1; i<=iend+1; i++) {
    iFinal[i]  = DP_ZERO;
  }

  for (int32 j=jbgn-1; j<=jend+1; j++) {
    lastCol[j] = DP_ZERO;
    thisCol[j] = DP_ZERO;
    jFinal[j]  = DP_ZERO;
  }

  //  Catch an invalid use case

  if ((endToEnd == true) &&
      ((ahang != 0) || (bhang != 0))) {
    fprintf(stderr, "Invalid code path.  I shouldn't be here.  Please report to atg@jcvi.org, thanks.\n");
    assert(0);
  }

  int32 scoreMax  = 0;

  int32 endI=0, curI=0;
  int32 endJ=0, curJ=0;

#ifdef DEBUG
  fprintf(stderr, "%d,%d - %d,%d -- ahang,bhang %d,%d  alen,blen %d,%d\n",
          ibgn, jbgn, iend, jend, ahang, bhang, lenA, lenB);
#endif

  for (int32 i=ibgn; i<=iend; i++) {
    for (int32 j=jbgn; j<=jend; j++) {

      //  Pick the max of these

      //  i-1 -> lastCol
      //  i   -> thisCol

      uint32 ul = lastCol[j-1] + ((stringA[i-1] == stringB[j-1]) ? MATCHSCORE : MISMATCHSCORE);
      uint32 lf = lastCol[j  ] + GAPSCORE;  //  Gap in B
      uint32 up = thisCol[j-1] + GAPSCORE;  //  Gap in A

      //  When computing unitig consensus, the N letter is used to indicate a gap.  This might be a
      //  match, or it might be a mismatch.  We just don't know.  Don't count it as anything.  This
      //  should be extended to allow ambiguitity codes -- but that then needs base calling support.
      //
      if (allowNs)
        if ((stringA[i-1] == 'N') ||
            (stringA[i-1] == 'n') ||
            (stringB[j-1] == 'N') ||
            (stringB[j-1] == 'n')) {
          ul = lastCol[j-1] + MATCHSCORE;
        }

      //  For unitig consensus, if the consensus sequence (stringA) is lowercase, count it as a
      //  match, otherwise ignore the gap it induces in stringB.
      //  
#warning stop complaing that this is ambiguous, please.
      if (stringA[i-1] == 'a')
        if (stringB[j-1] == 'A')
          ul = lastCol[j-1] + MATCHSCORE;
        else {
          ul = lastCol[j-1];
          lf = lastCol[j];
        }

      if (stringA[i-1] == 'c')
        if (stringB[j-1] == 'C')
          ul = lastCol[j-1] + MATCHSCORE;
        else {
          ul = lastCol[j-1];
          lf = lastCol[j];
        }

      if (stringA[i-1] == 'g')
        if (stringB[j-1] == 'G')
          ul = lastCol[j-1] + MATCHSCORE;
        else {
          ul = lastCol[j-1];
          lf = lastCol[j];
        }

      if (stringA[i-1] == 't')
        if (stringB[j-1] == 'T')
          ul = lastCol[j-1] + MATCHSCORE;
        else {
          ul = lastCol[j-1];
          lf = lastCol[j];
        }


      if (endToEnd) {
        //  Set score to the smallest value possible; we will then ALWAYS pick an action below.
        //
        thisCol[j]  = 0;
        ACT.set(i, j, MATCH);
        assert(ACT.get(i, j) == MATCH);
      } else {
        //  (i,j) is the beginning of a subsequence, our default behavior.  If we're looking for
        //  local alignments, make the alignment stop here.
        //
        thisCol[j]  = DP_ZERO;
        ACT.set(i, j, STOP);
        assert(ACT.get(i, j) == STOP);
      }

      if (thisCol[j] < ul) {
        thisCol[j]  = ul;
        ACT.set(i, j, MATCH);
        assert(ACT.get(i, j) == MATCH);
      }

      if (thisCol[j] < lf) {
        thisCol[j]  = lf;
        ACT.set(i, j, GAPB);
        assert(ACT.get(i, j) == GAPB);
      }

      if (thisCol[j] < up) {
        thisCol[j]  = up;
        ACT.set(i, j, GAPA);
        assert(ACT.get(i, j) == GAPA);
      }

      if (scoreMax < thisCol[j]) {
        scoreMax  = thisCol[j];
        endI = curI = i;
        endJ = curJ = j;
      }


      if (j == jend)
        jFinal[i] = thisCol[j];
      if (i == iend)
        iFinal[j] = thisCol[j];  //  Not really needed (just a copy of thisCol), but simplifies later code
    }  //  Over j


    //  Swap the columns
    uint32  *tc = thisCol;
    thisCol = lastCol;
    lastCol = tc;
  }  //  Over i

  //  If we're not looking for local alignments, scan the end points for the best value.  If we are
  //  looking for local alignments, we've already found and remembered the best end point.
  //
  //  lastCol - [iend][]
  //  jFinal  - [][jend]
  //

  if (endToEnd) {
    //fprintf(stderr, "ENDTOEND  curI %u curJ %u\n", curI, curJ);
    scoreMax    = 0;
    endI = curI = 0;
    endJ = curJ = 0;

    for (int32 i=ibgn, j=jend; i<=iend; i++) {
      if (scoreMax < jFinal[i]) {
        scoreMax = jFinal[i];
        endI = curI = i;
        endJ = curJ = j;
        //fprintf(stderr, "IscoreMax = %d at %d,%d\n", scoreMax - DP_ZERO, i, j);
      }
    }

    for (int32 j=jbgn, i=iend; j<=jend; j++) {
      if (scoreMax < iFinal[j]) {
        scoreMax = iFinal[j];
        endI = curI = i;
        endJ = curJ = j;
        //fprintf(stderr, "JscoreMax = %d at %d,%d\n", scoreMax - DP_ZERO, i, j);
      }
    }

    //  Not sure what this is for.
    //M[endI][endJ].score = 0;
  }

  //fprintf(stderr, "FINAL  curI %u curJ %u\n", curI, curJ);

  int32  alignLen  = 0;

  int32  nGapA     = 0;
  int32  nGapB     = 0;
  int32  nMatch    = 0;
  int32  nMismatch = 0;

  int32  terminate = 0;

  while (terminate == 0) {
    uint32 act = ACT.get(curI, curJ);

    //fprintf(stderr, "ACTION %2u curI %u curJ %u\n", act, curI, curJ);

    switch (act) {
      case STOP:
        terminate = 1;
        break;
      case MATCH:
        alignA[alignLen] = stringA[curI-1];
        alignB[alignLen] = stringB[curJ-1];

        if (alignA[alignLen] == alignB[alignLen]) {
          alignA[alignLen] = tolower(alignA[alignLen]);
          alignB[alignLen] = tolower(alignB[alignLen]);
          nMatch++;
        } else {
          alignA[alignLen] = toupper(alignA[alignLen]);
          alignB[alignLen] = toupper(alignB[alignLen]);
          nMismatch++;
        }

        curI--;
        curJ--;
        alignLen++;
        break;
      case GAPA:
        alignA[alignLen] = '-';
        alignB[alignLen] = stringB[curJ-1];
        curJ--;
        alignLen++;
        nGapA++;
        break;
      case GAPB:
        alignA[alignLen] = stringA[curI-1];
        alignB[alignLen] = '-';
        curI--;
        alignLen++;
        nGapB++;
        break;
    }
  }

  alignA[alignLen] = 0;
  alignB[alignLen] = 0;

  reverse(alignA, alignB, alignLen);

  a->matches    = nMatch;
  a->alignLen   = alignLen;
  a->begI       = curI;
  a->begJ       = curJ;
  a->endI       = endI;
  a->endJ       = endJ;
  a->lenA       = lenA;
  a->lenB       = lenB;
  a->pIdentity  = (double)(nMatch) / (double)(nGapA + nGapB + nMatch + nMismatch);
  a->pCoverageA = (double)(a->endI - a->begI) / (double)(lenA);
  a->pCoverageB = (double)(a->endJ - a->begJ) / (double)(lenB);
}

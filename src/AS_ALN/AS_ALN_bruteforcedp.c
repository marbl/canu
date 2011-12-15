
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

static const char *rcsid = "$Id: AS_ALN_bruteforcedp.c,v 1.19 2011-12-15 19:54:30 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_ALN_bruteforcedp.h"
#include "AS_UTL_reverseComplement.h"

#undef DEBUG

#define MATCH            0
#define GAPA             1
#define GAPB             2
#define STOP             3

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
            dpCell        (*M)[AS_READ_MAX_NORMAL_LEN + 1],
            alignLinker_s  *a,
            int             endToEnd,
            int             allowNs,
            int             ahang,
            int             bhang) {

  int32 lenA = strlen(stringA);
  int32 lenB = strlen(stringB);

  memset(a, 0, sizeof(alignLinker_s));

  if ((lenA > AS_READ_MAX_NORMAL_LEN) || (lenB > AS_READ_MAX_NORMAL_LEN))
    return;

  //  Definition of the box we want to do dynamic programming in.
  int32 ibgn = 1;
  int32 iend = lenA;
  int32 jbgn = 1;
  int32 jend = lenB;

  //  Set the edges.

  for (int32 i=ibgn-1, j=jbgn-1; i<=iend+1; i++) {
    M[i][j].score  = DP_ZERO;
    M[i][j].action = STOP;
  }

  for (int32 i=ibgn-1, j=jbgn-1; j<=jend+1; j++) {
    M[i][j].score  = DP_ZERO;
    M[i][j].action = STOP;
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

      int ul = M[i-1][j-1].score + ((stringA[i-1] == stringB[j-1]) ? MATCHSCORE : MISMATCHSCORE);
      int lf = M[i-1][j].score + GAPSCORE;  //  Gap in B
      int up = M[i][j-1].score + GAPSCORE;  //  Gap in A

      //  When computing unitig consensus, the N letter is used to indicate a gap.  This might be a
      //  match, or it might be a mismatch.  We just don't know.  Don't count it as anything.  This
      //  should be extended to allow ambiguitity codes -- but that then needs base calling support.
      //
      if (allowNs)
        if ((stringA[i-1] == 'N') ||
            (stringA[i-1] == 'n') ||
            (stringB[j-1] == 'N') ||
            (stringB[j-1] == 'n'))
          ul = M[i-1][j-1].score + MATCHSCORE;

      //  For unitig consensus, if the consensus sequence (stringA) is lowercase, count it as a
      //  match, otherwise ignore the gap it induces in stringB.
      //  
      if (stringA[i-1] == 'a')
        if (stringB[j-1] == 'A')
          ul = M[i-1][j-1].score + MATCHSCORE;
        else {
          ul = M[i-1][j-1].score;
          lf = M[i-1][j].score;
        }

      if (stringA[i-1] == 'c')
        if (stringB[j-1] == 'C')
          ul = M[i-1][j-1].score + MATCHSCORE;
        else {
          ul = M[i-1][j-1].score;
          lf = M[i-1][j].score;
        }

      if (stringA[i-1] == 'g')
        if (stringB[j-1] == 'G')
          ul = M[i-1][j-1].score + MATCHSCORE;
        else {
          ul = M[i-1][j-1].score;
          lf = M[i-1][j].score;
        }

      if (stringA[i-1] == 't')
        if (stringB[j-1] == 'T')
          ul = M[i-1][j-1].score + MATCHSCORE;
        else {
          ul = M[i-1][j-1].score;
          lf = M[i-1][j].score;
        }

      if (endToEnd) {
        //  Set score to the smallest value possible; we will then ALWAYS pick an action below.
        //
        M[i][j].score  = 0;
        M[i][j].action = MATCH;
      } else {
        //  (i,j) is the beginning of a subsequence, our default behavior.  If we're looking for
        //  local alignments, make the alignment stop here.
        //
        M[i][j].score  = DP_ZERO;
        M[i][j].action = STOP;
      }

      if (M[i][j].score < ul) {
        M[i][j].score  = ul;
        M[i][j].action = MATCH;
      }

      if (M[i][j].score < lf) {
        M[i][j].score  = lf;
        M[i][j].action = GAPB;
      }

      if (M[i][j].score < up) {
        M[i][j].score  = up;
        M[i][j].action = GAPA;
      }

      if (scoreMax < M[i][j].score) {
        scoreMax  = M[i][j].score;
        endI = curI = i;
        endJ = curJ = j;
      }
    }
  }

  //  If we're not looking for local alignments, scan the end points for the best value.  If we are
  //  looking for local alignments, we've already found and remembered the best end point.

  if (endToEnd) {
    scoreMax    = 0;
    endI = curI = 0;
    endJ = curJ = 0;

    for (int32 i=ibgn, j=jend; i<=iend; i++) {
      if (scoreMax < M[i][j].score) {
        scoreMax  = M[i][j].score;
        endI = curI = i;
        endJ = curJ = j;
        //fprintf(stderr, "IscoreMax = %d at %d,%d\n", scoreMax - DP_ZERO, i, j);
      }
    }

    for (int32 j=jbgn, i=iend; j<=jend; j++) {
      if (scoreMax < M[i][j].score) {
        scoreMax  = M[i][j].score;
        endI = curI = i;
        endJ = curJ = j;
        //fprintf(stderr, "JscoreMax = %d at %d,%d\n", scoreMax - DP_ZERO, i, j);
      }
    }

    M[endI][endJ].score = 0;
  }

  int32  alignLen  = 0;

  int32  nGapA     = 0;
  int32  nGapB     = 0;
  int32  nMatch    = 0;
  int32  nMismatch = 0;

  int32  terminate = 0;

  while (terminate == 0) {
    switch (M[curI][curJ].action) {
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

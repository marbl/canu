
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

static const char *rcsid = "$Id: AS_ALN_bruteforcedp.c,v 1.4 2008-11-12 12:44:46 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_ALN_bruteforcedp.h"
#include "AS_UTL_reverseComplement.h"


#define MATCH            0
#define GAPA             1
#define GAPB             2
#define STOP             3

#define MATCHSCORE       3
#define GAPSCORE        -2
#define MISMATCHSCORE   -2


void
alignLinker(char           *alignA,
            char           *alignB,
            char           *stringA,
            char           *stringB,
            dpCell        (*M)[AS_READ_MAX_LEN],
            alignLinker_s  *a,
            int             endToEnd, int abeg, int aend) {

  int   lenA = strlen(stringA);
  int   lenB = strlen(stringB);

  int   i, j;

  abeg = MAX(0, abeg);

  for (i=0; i<lenA+1; i++) {
    M[i][0].score  = 1 << 29;
    M[i][0].action = STOP;

    M[i][lenB].score  = 0;
    M[i][lenB].action = STOP;
  }

  for (j=0; j<lenB+1; j++) {
    M[abeg][j].score  = 1 << 29;
    M[abeg][j].action = STOP;

    M[lenA][j].score  = 0;
    M[lenA][j].action = STOP;
  }

  int   scoreMax  = 0;

  int   endI=0, curI=0;
  int   endJ=0, curJ=0;

  for (i=abeg+1; i<=lenA; i++){
    int  jb = 1;
    int  je = lenB;

#if 0
    //  Attempt at further banding; doesn't work.
    //
    if (endToEnd) {
      if (i <= aend) {
        jb = 1;
        je = 3 * (i-abeg) / 2;
      } else {
        jb = 2 * (i-aend) / 3;
        je = 3 * (i-abeg) / 2;
      }

      if (jb < 1)
        jb = 1;

      if (jb > lenB)
        jb = lenB;
      if (je > lenB)
        je = lenB;

      M[i-1][je].score  = 1 << 29;
      M[i-1][je].action = STOP;

      M[i-1][jb-1].score  = 1 << 29;
      M[i-1][jb-1].action = STOP;

      M[i][jb-1].score  = 1 << 29;
      M[i][jb-1].action = STOP;
    }
#endif

    for (j=jb; j<=je; j++){

      //  Pick the max of these

      int ul = M[i-1][j-1].score + ((stringA[i-1] == stringB[j-1]) ? MATCHSCORE : MISMATCHSCORE);
      int lf = M[i-1][j].score + GAPSCORE;
      int up = M[i][j-1].score + GAPSCORE;

      //  (i,j) is the beginning of a subsequence, our default
      //  behavior.  If we're looking for local alignments, make the
      //  alignment stop here.  Otherwise, it's a match.
      //
      M[i][j].score  = 1 << 29;
      M[i][j].action = (endToEnd) ? MATCH : STOP;

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

  //  If we're not looking for local alignments, scan the end points
  //  for the best value.
  //
 findAnother:
  if (endToEnd) {
    scoreMax = 0;
    endI     = 0;
    endJ     = 0;

    //  Instead of these, it is somewhat easy to look for an overlap
    //  with a given ahang,bhang.

    //  Search for the best dovetail overlap
    for (i=lenA, j=0; j<=lenB; j++) {
      if (scoreMax < M[i][j].score) {
        scoreMax  = M[i][j].score;
        endI = curI = i;
        endJ = curJ = j;
      }
    }

    //  Search for the best contained overlap
    for (i=MAX(0, abeg), j=lenB; i<=lenA; i++) {
      if (scoreMax < M[i][j].score) {
        scoreMax  = M[i][j].score;
        endI = curI = i;
        endJ = curJ = j;
      }
    }

    M[endI][endJ].score = 0;
  }

  assert(scoreMax > 0);

  int  alignLen  = 0;
  int  matches   = 0;
  int  terminate = 0;

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
          matches++;
        } else {
          alignA[alignLen] = toupper(alignA[alignLen]);
          alignB[alignLen] = toupper(alignB[alignLen]);
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
        break;
      case GAPB:
        alignA[alignLen] = stringA[curI-1];
        alignB[alignLen] = '-';
        curI--;
        alignLen++;
        break;
    }
  }

  if (endToEnd) {
    //  Life would have been so much easier if we got ahang and bhang,
    //  instead of getting a requested start position in the A
    //  sequence.
    //
    //  One could probably fix this by reversing the sequences then
    //  aligning.
    //
    if ((curI < abeg) || (aend < curI)) {
      fprintf(stderr, "ANOTHER A %d-%d B %d-%d beg,end %d,%d\n", curI, endI, curJ, endJ, abeg, aend);
      goto findAnother;
    }
  }

  alignA[alignLen] = 0;
  alignB[alignLen] = 0;

  reverse(alignA, alignB, alignLen);

  a->matches  = matches;
  a->alignLen = alignLen;
  a->begI     = curI;
  a->begJ     = curJ;
  a->endI     = endI;
  a->endJ     = endJ;
  a->lenA     = lenA;
  a->lenB     = lenB;
}

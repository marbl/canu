// This file is part of A2Amapper.
// Copyright (c) 2005 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "bio++.H"
#include "atac.H"


u32bit  gapLimit = 5;


//  Reads a set of matches and shifts the location of all the gaps to
//  one side or the other.
//
//  Example:
//
//    A middle gap:
//    GGGGGGGGGGATATATATATATATATATATATATGGGGGGGGG
//    GGGGGGGGGGATAT--ATATATATATATATATATGGGGGGGGG
//
//    A left-most gap:
//    GGGGGGGGGGATATATATATATATATATATATATGGGGGGGGG
//    GGGGGGGGGG--ATATATATATATATATATATATGGGGGGGGG
//
//    A right-most gap:
//    GGGGGGGGGGATATATATATATATATATATATATGGGGGGGGG
//    GGGGGGGGGGATATATATATATATATATATAT--GGGGGGGGG
//
//  Shifting is done for both assembly-axes.


////  SWITCHED to use atacMatchList, untested!


u32bit
shiftRight(atacMatch *m1, atacMatch *m2, FastACache *C1, FastACache *C2) {
  u32bit shifted = 0;

  //  Don't do anything if they're on different sequences
  if ((m1->iid1 != m2->iid1) || (m1->iid2 != m2->iid2))
    return(0);

  //  Don't do anything if the orientation of the two matches is
  //  different.  This is probably not a gap we want to muck with.
  //
  if (m1->fwd2 != m2->fwd2)
    return(0);

  //  Don't do anything if the length is zero
  if ((m1->len1 == 0) ||
      (m1->len2 == 0) ||
      (m2->len1 == 0) ||
      (m2->len2 == 0))
    return(0);

  //  Don't do anything if both sequences are not consecutive
  if ((m1->pos1 + m1->len1 != m2->pos1) &&
      (m1->pos2 + m1->len2 != m2->pos2))
    return(0);

  //  Don't do anything if they are overlapping in the first sequence
  if (m1->pos1 + m1->len1 > m2->pos1)
    return(0);

  //  Don't do anything if the gap is big -- if we are reverse, then m2 is before m1 on
  //  the other assembly!
  if (m1->fwd1 == 0) {
    if ((m1->pos1 + m1->len1 == m2->pos1) && (m2->pos2 + m2->len2 + gapLimit <= m1->pos2))
      return(0);
    if ((m1->pos2 + m1->len2 == m2->pos2) && (m2->pos1 + m2->len1 + gapLimit <= m1->pos1))
      return(0);
  } else {
    if ((m1->pos1 + m1->len1 == m2->pos1) && (m1->pos2 + m1->len2 + gapLimit <= m2->pos2))
      return(0);
    if ((m1->pos2 + m1->len2 == m2->pos2) && (m1->pos1 + m1->len1 + gapLimit <= m2->pos1))
      return(0);
  }

  atacMatch m1c;
  atacMatch m2c;

  memcpy(&m1c, m1, sizeof(atacMatch));
  memcpy(&m2c, m2, sizeof(atacMatch));

  //  Grab those sequences from the wrapper
  //
  FastASequenceInCore  *S1 = C1->getSequence(m1->iid1);
  FastASequenceInCore  *S2 = C2->getSequence(m1->iid2);

  FastAAccessor A1(S1, false);
  FastAAccessor A2(S2, (m1->fwd1 != m1->fwd2));
  A1.setRange(m1->pos1, m1->len1);
  A2.setRange(m1->pos2, m1->len2);

  FastAAccessor B1(S1, false);
  FastAAccessor B2(S2, (m2->fwd1 != m2->fwd2));
  B1.setRange(m2->pos1, m2->len1);
  B2.setRange(m2->pos2, m2->len2);

  A1.setPosition(m1->pos1 + m1->len1 - 1);
  A2.setPosition(m1->pos2 + m1->len2 - 1);
  ++A1;
  ++A2;

  B1.setPosition(m2->pos1);  //  Used only for output.
  B2.setPosition(m2->pos2);

  //  We want to extend m1 to the right, this will shift the gap
  //  to the right-most position (relative to the forward
  //  genomic).
  //
  //  While there is a match at the end of m1, extend m1 to the right,
  //  and chop the head off m2.
  //
  while (A1.isValid() && A2.isValid() && B1.isValid() && B2.isValid() &&
         (m2->len1 > 0) &&
         (m2->len2 > 0) &&
         validSymbol[(int)*A1] &&
         validSymbol[(int)*A2] &&
         IUPACidentity[(int)*A1][(int)*A2]) {
    A1.extendRight(1);
    A2.extendRight(1);

    B1.extendLeft(-1);
    B2.extendLeft(-1);

    //  Only used to make sure that this sequence isn't turned
    //  into a black hole.  Probably not needed.
    m2->len1--;
    m2->len2--;

    ++A1;
    ++A2;

    ++B1;
    ++B2;

    shifted++;
  }

  m1->pos1 = A1.getRangeBegin();
  m1->len1 = A1.getRangeLength();
  m1->pos2 = A2.getRangeBegin();
  m1->len2 = A2.getRangeLength();
  m2->pos1 = B1.getRangeBegin();
  m2->len1 = B1.getRangeLength();
  m2->pos2 = B2.getRangeBegin();
  m2->len2 = B2.getRangeLength();

#if 0
  if (shifted) {
    fprintf(stderr, "inp "u32bitFMT" "u32bitFMT" "u32bitFMT" 1  "u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
            m1c.iid1, m1c.pos1, m1c.len1,
            m1c.iid2, m1c.pos2, m1c.len2, m1c.fwd2 ? 1 : -1);
    fprintf(stderr, "inp "u32bitFMT" "u32bitFMT" "u32bitFMT" 1  "u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
            m2c.iid1, m2c.pos1, m2c.len1,
            m2c.iid2, m2c.pos2, m2c.len2, m2c.fwd2 ? 1 : -1);
    fprintf(stderr, "shifted "u32bitFMT"\n", shifted);
    fprintf(stderr, "out "u32bitFMT" "u32bitFMT" "u32bitFMT" 1  "u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
            m1->iid1, m1->pos1, m1->len1,
            m1->iid2, m1->pos2, m1->len2, m1->fwd2 ? 1 : -1);
    fprintf(stderr, "out "u32bitFMT" "u32bitFMT" "u32bitFMT" 1  "u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
            m2->iid1, m2->pos1, m2->len1,
            m2->iid2, m2->pos2, m2->len2, m2->fwd2 ? 1 : -1);
  }
#endif

  u32bit errors = 0;

  if (m1c.pos1 != m1->pos1) {
    errors++;
    fprintf(stderr, "WARNING:  begin of assembly 1 moved!\n");
  }
  if (m2c.pos1 + m2c.len1 != m2->pos1 + m2->len1) {
    errors++;
    fprintf(stderr, "WARNING:  end of assembly 1 moved!\n");
  }

  if (m1->fwd2 == 0) {
    if (m2c.pos2 != m2->pos2) {
      errors++;
      fprintf(stderr, "WARNING:  begin of assembly 2 moved (rc)!\n");
    }

    if (m1c.pos2 + m1c.len2 != m1->pos2 + m1->len2) {
      errors++;
      fprintf(stderr, "WARNING:  end of assembly 2 moved (rc)!\n");
    }
  } else {
    if (m1c.pos2 != m1->pos2) {
      errors++;
      fprintf(stderr, "WARNING:  begin of assembly 2 moved!\n");
    }

    if (m2c.pos2 + m2c.len2 != m2->pos2 + m2->len2) {
      errors++;
      fprintf(stderr, "WARNING:  end of assembly 1 moved!\n");
    }
  }

  if (errors)
    exit(1);

  return(shifted);
}





int
main(int argc, char *argv[]) {

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-l") == 0) {
    } else if (strcmp(argv[arg], "-r") == 0) {
    } else if (strcmp(argv[arg], "-g") == 0) {
      gapLimit = strtou32bit(argv[++arg], 0L);
    } else {
      fprintf(stderr, "usage: %s [-l | -r] [-g limit] < matches > matches\n", argv[0]);
      fprintf(stderr, "  -l:     Shift left (UNIMPLEMENTED!)\n");
      fprintf(stderr, "  -r:     Shift right\n");
      exit(1);
    }
    arg++;
  }

  fprintf(stderr, "gapLimit is "u32bitFMT"\n", gapLimit);

  atacMatchList  ML("-", 'm', false);
  ML.sort1();

  FastACache  *C1 = new FastACache(ML.assemblyFileA(),    2, true, false);
  FastACache  *C2 = new FastACache(ML.assemblyFileB(), 1024, true, false);

  u32bit gapsShifted = 0;

  for (u32bit i=1; i<ML.numMatches(); i++)
    if (shiftRight(ML[i-1], ML[i], C1, C2)) {
      gapsShifted++;
      //fprintf(stderr, "Sort 1: shifted "u32bitFMT" out of "u32bitFMT" (%6.2f%%)\r", gapsShifted, i, (double)gapsShifted / (double)i * 100.0);
      //fflush(stderr);
    }
  fprintf(stderr, "Sort 1: shifted "u32bitFMT" out of "u32bitFMT" (%6.2f%%)\n", gapsShifted, ML.numMatches(), (double)gapsShifted / (double)ML.numMatches() * 100.0);


  //  Dump all the matches.
  //
  for (u32bit i=0; i<ML.numMatches(); i++) {
    atacMatch *m1 = ML[i];
    if ((m1->len1 > 0) && (m1->len2 > 0))
      fprintf(stdout, "M u "u32bitFMT" . B35LC:"u32bitFMT" "u32bitFMT" "u32bitFMT" 1  HUREF2:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
              i,
              m1->iid1, m1->pos1, m1->len1,
              m1->iid2, m1->pos2, m1->len2, m1->fwd2 ? 1 : -1);
  }

  return(0);
}

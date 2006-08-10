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

#define DEBUG_SHIFTING

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


u32bit
shiftGap(atacMatch *ma, atacMatch *mb,
         FastACache *C1, FastACache *C2,
         u32bit gapLimit,
         bool shiftRight) {

  u32bit shifted = 0;

  //  Don't do anything if they're on different sequences
  //
  if ((ma->iid1 != mb->iid1) || (ma->iid2 != mb->iid2))
    return(0);

  //  Don't do anything if the orientation of the two matches is
  //  different.  This is probably not a gap we want to muck with.
  //
  if (ma->fwd2 != mb->fwd2)
    return(0);

  //  Don't do anything if the length is zero
  //
  if ((ma->len1 == 0) ||
      (ma->len2 == 0) ||
      (mb->len1 == 0) ||
      (mb->len2 == 0))
    return(0);

  //  Don't do anything if both sequences are not consecutive -- if
  //  both are consecutive they should have been merged already, so
  //  this leaves us with exactly one sequence with a gap.
  //
  if ((ma->pos1 + ma->len1 != mb->pos1) &&
      (ma->pos2 + ma->len2 != mb->pos2))
    return(0);

  //  And just to be sure, don't do anything if there really is no gap
  //  -- why does this occur??
  //
  if ((ma->pos1 + ma->len1 == mb->pos1) &&
      (ma->pos2 + ma->len2 == mb->pos2))
    return(0);

  //  Don't do anything if they are overlapping
  //
  if (ma->pos1 + ma->len1 > mb->pos1)
    return(0);
  if (ma->pos2 + ma->len2 > mb->pos2)
    return(0);

  //  Don't do anything if the gap is big -- if we are reverse, then mb is before ma on
  //  the other assembly!
  //
  if (ma->fwd2 == 0) {
    if ((ma->pos1 + ma->len1 == mb->pos1) && (mb->pos2 + mb->len2 + gapLimit < ma->pos2))
      return(0);
    if ((ma->pos2 + ma->len2 == mb->pos2) && (mb->pos1 + mb->len1 + gapLimit < ma->pos1))
      return(0);
  } else {
    if ((ma->pos1 + ma->len1 == mb->pos1) && (ma->pos2 + ma->len2 + gapLimit < mb->pos2))
      return(0);
    if ((ma->pos2 + ma->len2 == mb->pos2) && (ma->pos1 + ma->len1 + gapLimit < mb->pos1))
      return(0);
  }

  atacMatch mac;
  atacMatch mbc;

  memcpy(&mac, ma, sizeof(atacMatch));
  memcpy(&mbc, mb, sizeof(atacMatch));

  //  Grab those sequences from the wrapper
  //
  FastASequenceInCore  *S1 = C1->getSequence(ma->iid1);
  FastASequenceInCore  *S2 = C2->getSequence(ma->iid2);

  FastAAccessor A1(S1, false);
  FastAAccessor A2(S2, (ma->fwd1 != ma->fwd2));
  A1.setRange(ma->pos1, ma->len1);
  A2.setRange(ma->pos2, ma->len2);

  FastAAccessor B1(S1, false);
  FastAAccessor B2(S2, (mb->fwd1 != mb->fwd2));
  B1.setRange(mb->pos1, mb->len1);
  B2.setRange(mb->pos2, mb->len2);


  //  We want to extend ma to the right, this will shift the gap
  //  to the right-most position (relative to the forward
  //  genomic).
  //
  //  While there is a match at the end of ma, extend ma to the right,
  //  and chop the head off mb.
  //
  if (shiftRight) {

    //  A wants to be the first thing after ma -- the first base in the
    //  gap.  Set the position to the last thing in the range, then use
    //  the increment operator to extend past that.  We can't directly
    //  go here due to reverse complement issues.
    //
    A1.setPosition(ma->pos1 + ma->len1 - 1);
    A2.setPosition(ma->pos2 + ma->len2 - 1);
    ++A1;
    ++A2;

    //  But B can be set to the first thing in the match with no
    //  problem.
    //
    B1.setPosition(mb->pos1);  //  Used only for output.
    B2.setPosition(mb->pos2);

    //  While we're still in sequence (isValid()) and we haven't
    //  obliterated the match we're shifting the gap into, and we can
    //  extend the other match (being both validSymbols and identity),
    //  shift the gap to the right.
    //
    while (A1.isValid() && A2.isValid() && B1.isValid() && B2.isValid() &&
           (mb->len1 > 0) &&
           (mb->len2 > 0) &&
           validSymbol[(int)*A1] &&
           validSymbol[(int)*A2] &&
           IUPACidentity[(int)*A1][(int)*A2]) {
      A1.extendRight(1);
      A2.extendRight(1);

      B1.extendLeft(-1);
      B2.extendLeft(-1);

      //  Only used to make sure that this sequence isn't turned
      //  into a black hole.  Probably not needed.
      mb->len1--;
      mb->len2--;

      ++A1;
      ++A2;

      ++B1;
      ++B2;

      shifted++;
    }
  } else {

    //  Similar to above, A wants to be the last thing in the match,
    //  and B the first thing in the gap.
    //
    A1.setPosition(ma->pos1 + ma->len1 - 1);
    A2.setPosition(ma->pos2 + ma->len2 - 1);

    B1.setPosition(mb->pos1);
    B2.setPosition(mb->pos2);
    --B1;
    --B2;

    while (A1.isValid() && A2.isValid() && B1.isValid() && B2.isValid() &&
           (ma->len1 > 0) &&
           (ma->len2 > 0) &&
           validSymbol[(int)*B1] &&
           validSymbol[(int)*B2] &&
           IUPACidentity[(int)*B1][(int)*B2]) {
      A1.extendRight(-1);
      A2.extendRight(-1);

      B1.extendLeft(1);
      B2.extendLeft(1);

      //  Only used to make sure that this sequence isn't turned
      //  into a black hole.  Probably not needed.
      ma->len1--;
      ma->len2--;

      --A1;
      --A2;

      --B1;
      --B2;

      shifted++;
    }
  }

  ma->pos1 = A1.getRangeBegin();
  ma->len1 = A1.getRangeLength();
  ma->pos2 = A2.getRangeBegin();
  ma->len2 = A2.getRangeLength();
  mb->pos1 = B1.getRangeBegin();
  mb->len1 = B1.getRangeLength();
  mb->pos2 = B2.getRangeBegin();
  mb->len2 = B2.getRangeLength();


  u32bit errors = 0;

  if (mac.pos1 != ma->pos1) {
    errors++;
    fprintf(stderr, "WARNING:  begin of assembly 1 moved!\n");
  }
  if (mbc.pos1 + mbc.len1 != mb->pos1 + mb->len1) {
    errors++;
    fprintf(stderr, "WARNING:  end of assembly 1 moved!\n");
  }

  if (ma->fwd2 == 0) {
    if (mbc.pos2 != mb->pos2) {
      errors++;
      fprintf(stderr, "WARNING:  begin of assembly 2 moved (rc)!\n");
    }

    if (mac.pos2 + mac.len2 != ma->pos2 + ma->len2) {
      errors++;
      fprintf(stderr, "WARNING:  end of assembly 2 moved (rc)!\n");
    }
  } else {
    if (mac.pos2 != ma->pos2) {
      errors++;
      fprintf(stderr, "WARNING:  begin of assembly 2 moved!\n");
    }

    if (mbc.pos2 + mbc.len2 != mb->pos2 + mb->len2) {
      errors++;
      fprintf(stderr, "WARNING:  end of assembly 1 moved!\n");
    }
  }

  if ((errors > 0) || (shifted > 999999)) {
    atacMatch *l = 0L;
    atacMatch *r = 0L;
    char    str1[1024];
    char    str2[1024];
    char    str3[1024];
    u32bit  i;
    u32bit  len;
    u32bit  pos;

#define MAXPRINT 90

    fprintf(stderr, "----------------------------------------\n");

    mac.print(stderr, "A", "B");
    mbc.print(stderr, "A", "B");

    fprintf(stderr, u32bitFMT"-"u32bitFMT" -- "u32bitFMT"-"u32bitFMT"\n",
            mac.pos1, mac.pos1 + mac.len1,
            mbc.pos1, mbc.pos1 + mbc.len1);
    fprintf(stderr, u32bitFMT"-"u32bitFMT" -- "u32bitFMT"-"u32bitFMT"\n",
            mac.pos2, mac.pos2 + mac.len2,
            mbc.pos2, mbc.pos2 + mbc.len2);

    l = &mac;
    r = &mbc;

    str1[0] = 0;
    str2[0] = 0;
    str3[0] = 0;

    len = l->len1;
    if (len > 0) {
      if (len > MAXPRINT)
        len = MAXPRINT;
      pos = l->pos1 + l->len1 - len;
      A1.setRange(pos, len);
      A1.setPosition(pos);
      for (i=0; i<len; i++, ++A1)
        str1[i] = *A1;
      str1[i] = 0;
    }

    len = r->pos1 - l->pos1 - l->len1;
    if (len > 0) {
      pos = l->pos1 + l->len1;
      A1.setRange(pos, len);
      A1.setPosition(pos);
      for (i=0; i<len; i++, ++A1)
        str2[i] = *A1;
      str2[i] = 0;
    }

    len = r->len1;
    if (len > 0) {
      if (len > MAXPRINT)
        len = MAXPRINT;
      pos = r->pos1;
      A1.setRange(pos, len);
      A1.setPosition(pos);
      for (i=0; i<len; i++, ++A1)
        str3[i] = *A1;
      str3[i] = 0;
    }

    fprintf(stderr, "SEQA: %s -- %s -- %s\n", str1, str2, str3);


    if (l->fwd2 == 0) {
      l = &mbc;
      r = &mac;
    }

    str1[0] = 0;
    str2[0] = 0;
    str3[0] = 0;

    len = l->len2;
    if (len > 0) {
      if (len > MAXPRINT)
        len = MAXPRINT;
      pos = l->pos2 + l->len2 - len;
      A2.setRange(pos, len);
      A2.setPosition(pos);
      for (i=0; i<len; i++, ++A2)
        str1[i] = *A2;
      str1[i] = 0;
    }

    len = r->pos2 - l->pos2 - l->len2;
    if (len > 0) {
      pos = l->pos2 + l->len2;
      A2.setRange(pos, len);
      A2.setPosition(pos);
      for (i=0; i<len; i++, ++A2)
        str2[i] = *A2;
      str2[i] = 0;
    }

    len = r->len2;
    if (len > 0) {
      if (len > MAXPRINT)
        len = MAXPRINT;
      pos = r->pos2;
      A2.setRange(pos, len);
      A2.setPosition(pos);
      for (i=0; i<len; i++, ++A2)
        str3[i] = *A2;
      str3[i] = 0;
    }

    fprintf(stderr, "SEQB: %s -- %s -- %s\n", str1, str2, str3);


    fprintf(stderr, "shifted "u32bitFMT"\n", shifted);


    ma->print(stderr, "A", "B");
    mb->print(stderr, "A", "B");

    fprintf(stderr, u32bitFMT"-"u32bitFMT" -- "u32bitFMT"-"u32bitFMT"\n",
            ma->pos1, ma->pos1 + ma->len1,
            mb->pos1, mb->pos1 + mb->len1);
    fprintf(stderr, u32bitFMT"-"u32bitFMT" -- "u32bitFMT"-"u32bitFMT"\n",
            ma->pos2, ma->pos2 + ma->len2,
            mb->pos2, mb->pos2 + mb->len2);

    l = ma;
    r = mb;

    str1[0] = 0;
    str2[0] = 0;
    str3[0] = 0;

    len = l->len1;
    if (len > 0) {
      if (len > MAXPRINT)
        len = MAXPRINT;
      pos = l->pos1 + l->len1 - len;
      A1.setRange(pos, len);
      A1.setPosition(pos);
      for (i=0; i<len; i++, ++A1)
        str1[i] = *A1;
      str1[i] = 0;
    }

    len = r->pos1 - l->pos1 - l->len1;
    if (len > 0) {
      pos = l->pos1 + l->len1;
      A1.setRange(pos, len);
      A1.setPosition(pos);
      for (i=0; i<len; i++, ++A1)
        str2[i] = *A1;
      str2[i] = 0;
    }

    len = r->len1;
    if (len > 0) {
      if (len > MAXPRINT)
        len = MAXPRINT;
      pos = r->pos1;
      A1.setRange(pos, len);
      A1.setPosition(pos);
      for (i=0; i<len; i++, ++A1)
        str3[i] = *A1;
      str3[i] = 0;
    }

    fprintf(stderr, "SEQA: %s -- %s -- %s\n", str1, str2, str3);


    if (l->fwd2 == 0) {
      l = mb;
      r = ma;
    }

    str1[0] = 0;
    str2[0] = 0;
    str3[0] = 0;

    len = l->len2;
    if (len > 0) {
      if (len > MAXPRINT)
        len = MAXPRINT;
      pos = l->pos2 + l->len2 - len;
      A2.setRange(pos, len);
      A2.setPosition(pos);
      for (i=0; i<len; i++, ++A2)
        str1[i] = *A2;
      str1[i] = 0;
    }

    len = r->pos2 - l->pos2 - l->len2;
    if (len > 0) {
      pos = l->pos2 + l->len2;
      A2.setRange(pos, len);
      A2.setPosition(pos);
      for (i=0; i<len; i++, ++A2)
        str2[i] = *A2;
      str2[i] = 0;
    }

    len = r->len2;
    if (len > 0) {
      if (len > MAXPRINT)
        len = MAXPRINT;
      pos = r->pos2;
      A2.setRange(pos, len);
      A2.setPosition(pos);
      for (i=0; i<len; i++, ++A2)
        str3[i] = *A2;
      str3[i] = 0;
    }

    fprintf(stderr, "SEQB: %s -- %s -- %s\n", str1, str2, str3);
  }

  if (errors)
    exit(1);

  return(shifted);
}





int
main(int argc, char *argv[]) {
  bool    shiftRight = true;
  u32bit  gapLimit   = 5;
  bool    error      = false;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-l") == 0) {
      shiftRight = false;
    } else if (strcmp(argv[arg], "-r") == 0) {
      shiftRight = true;
    } else if (strcmp(argv[arg], "-g") == 0) {
      gapLimit = strtou32bit(argv[++arg], 0L);
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      error = true;
    }
    arg++;
  }

  if (error) {
    fprintf(stderr, "usage: %s [-l | -r] [-g limit] < matches > matches\n", argv[0]);
    fprintf(stderr, "  -l:     Shift left\n");
    fprintf(stderr, "  -r:     Shift right (default)\n");
    fprintf(stderr, "  -g:     Maxinum size of gap to shift\n");
    exit(1);
  }

  fprintf(stderr, "gapLimit is "u32bitFMT"\n", gapLimit);

  atacMatchList  ML("-", 'm', stdout);
  ML.sort1();

  FastACache  *C1 = new FastACache(ML.assemblyFileA(),    2, true, false);
  FastACache  *C2 = new FastACache(ML.assemblyFileB(), 1024, true, false);


  u32bit gapsShifted = 0;

  for (u32bit i=1; i<ML.numMatches(); i++)
    if (shiftGap(ML[i-1], ML[i], C1, C2, gapLimit, shiftRight)) {
      gapsShifted++;
      fprintf(stderr, "shiftR: shifted "u32bitFMT" out of "u32bitFMT" (%6.2f%%)\r", gapsShifted, i, (double)gapsShifted / (double)i * 100.0);
      fflush(stderr);
    }
  fprintf(stderr, "shiftR: "u32bitFMT" out of "u32bitFMT" (%6.2f%%)\n", gapsShifted, ML.numMatches(), (double)gapsShifted / (double)ML.numMatches() * 100.0);


  for (u32bit i=0; i<ML.numMatches(); i++) {
    atacMatch *ma = ML[i];
    if ((ma->len1 > 0) && (ma->len2 > 0))
      fprintf(stdout, "M u %s %s %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" 1  %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
              ma->matchuid, ma->parentuid,
              ML.labelA(), ma->iid1, ma->pos1, ma->len1,
              ML.labelB(), ma->iid2, ma->pos2, ma->len2, ma->fwd2 ? 1 : -1);
  }

  return(0);
}

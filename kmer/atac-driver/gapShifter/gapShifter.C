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

#define MAXPRINT 90

#if 0
#define REPORT_RESULTS
#define REPORT_UNSHIFTABLE
#define REPORT_SHIFTING
#endif

//  Global statistics, reset on each iteration
//
u32bit   numShifted;     //  valid to be shifted, and were shifted
u32bit   numNotShifted;  //  but valid to be shifted
u32bit   numDiffSeq;     //  all the rest are not valid to be shifted
u32bit   numDiffOri;
u32bit   numZeroLen;
u32bit   numOutOfOrder;
u32bit   numNotAdjacent;
u32bit   numNoGap;
u32bit   numGapTooBig;
u32bit   numOverlapping;
u32bit   amountShifted[1024];

FILE    *logFile = 0L;

//  XXXXXX  outer loop needs to skip empty matches!!

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



//  Returns true if these two matches have a potentially shiftable gap
//  between them.  Potentially shiftable means that the matches are
//  contiguous on one axis, and consecutive (no matches between) on
//  the other axis.
//
//  Assume the matches are sorted by the first sequence.
//
bool
isPotentiallyShiftable(atacMatch *ma,
                       atacMatch *mb,
                       atacMatchOrder &MOB,
                       u32bit gapLimit) {

#ifdef REPORT_UNSHIFTABLE
  fprintf(stderr, "isPotentiallyShiftable()\n");
  ma->print(stderr, "A", "B");
  mb->print(stderr, "A", "B");
#endif

  //  Not shiftable if on different sequences
  //
  if ((ma->iid1 != mb->iid1) ||
      (ma->iid2 != mb->iid2)) {
#ifdef REPORT_UNSHIFTABLE
    fprintf(stderr, "UNSHIFTABLE different sequences\n");
#endif
    numDiffSeq++;
    return(false);
  }

  //  Not shiftable if the orientation of the two matches is
  //  different.  This is probably not a gap we want to muck with.
  //
  if ((ma->fwd1 != mb->fwd1) ||
      (ma->fwd2 != mb->fwd2)) {
#ifdef REPORT_UNSHIFTABLE
    fprintf(stderr, "UNSHIFTABLE different orientation\n");
#endif
    numDiffOri++;
    return(false);
  }

  //  Not shiftable if any length is zero.  This isn't a gap, it's a
  //  dead match.
  //
  if ((ma->len1 == 0) ||
      (ma->len2 == 0) ||
      (mb->len1 == 0) ||
      (mb->len2 == 0)) {
#ifdef REPORT_UNSHIFTABLE
    fprintf(stderr, "UNSHIFTABLE zero length\n");
#endif
    numZeroLen++;
    return(false);
  }

  atacMatch  *bl = (ma->fwd2) ? ma : mb;  //  the left match on B, relative to forward orientation
  atacMatch  *br = (ma->fwd2) ? mb : ma;  //  the right

  //  Not shiftable if the B matches are out of order
  //
  if (bl->pos2 > br->pos2) {
#ifdef REPORT_UNSHIFTABLE
    fprintf(stderr, "UNSHIFTABLE misordered on B\n");
#endif
    numOutOfOrder++;
    return(false);
  }

  u32bit magap = mb->pos1 - (ma->pos1 + ma->len1);
  u32bit mbgap = br->pos2 - (bl->pos2 + bl->len2);

  //  Not shiftable if there is no zero size gap
  //
  if ((magap > 0) && (mbgap > 0)) {
#ifdef REPORT_UNSHIFTABLE
    fprintf(stderr, "UNSHIFTABLE no zero size gap ("u32bitFMT", "u32bitFMT")\n", magap, mbgap);
#endif
    numNotAdjacent++;
    return(false);
  }

  //  Not shiftabe if there is no gap
  //
  if ((magap == 0) && (mbgap == 0)) {
#ifdef REPORT_UNSHIFTABLE
    fprintf(stderr, "UNSHIFTABLE no gap on both sequences ("u32bitFMT", "u32bitFMT")\n", magap, mbgap);
#endif
    numNoGap++;
    return(false);
  }

  //  Not shiftable if the gap is big
  //
  if ((magap > gapLimit) || (mbgap > gapLimit)) {
#ifdef REPORT_UNSHIFTABLE
    fprintf(stderr, "UNSHIFTABLE gap too big ("u32bitFMT", "u32bitFMT")\n", magap, mbgap);
#endif
    numGapTooBig++;
    return(false);
  }

  //  Not shiftable if they overlap
  //
  if (ma->pos1 + ma->len1 > mb->pos1) {
#ifdef REPORT_UNSHIFTABLE
    fprintf(stderr, "UNSHIFTABLE overlap on sequence A\n");
#endif
    numOverlapping++;
    return(false);
  }
  if (bl->pos2 + bl->len2 > br->pos2) {
#ifdef REPORT_UNSHIFTABLE
    fprintf(stderr, "UNSHIFTABLE overlap on sequence B\n");
#endif
    numOverlapping++;
    return(false);
  }

  u32bit iid1 = ma->matchiid;
  u32bit iid2 = mb->matchiid;

  //  Check that there isn't another match stuck in the middle on
  //  the B axis.
  //
  if (ma->fwd2 == true) {
    if ((MOB.index(iid1)+1) != MOB.index(iid2)) {
      fprintf(stderr, "WARNING:  Match inbetween!  (forward)\n");

      fprintf(stderr, "iid1 "u32bitFMT", iid2 "u32bitFMT"\n", iid1, iid2);
      ma->print(stderr, "A", "B");
      mb->print(stderr, "A", "B");

      fprintf(stderr, "before, iid1, after\n");
      MOB[MOB.index(iid1)-1]->print(stderr, "A", "B");
      MOB[MOB.index(iid1)  ]->print(stderr, "A", "B");
      MOB[MOB.index(iid1)+1]->print(stderr, "A", "B");

      fprintf(stderr, "before, iid2, after\n");
      MOB[MOB.index(iid2)-1]->print(stderr, "A", "B");
      MOB[MOB.index(iid2)  ]->print(stderr, "A", "B");
      MOB[MOB.index(iid2)+1]->print(stderr, "A", "B");

      return(false);
    }
  } else {
    if ((MOB.index(iid1)-1) != MOB.index(iid2)) {
      fprintf(stderr, "WARNING:  Match inbetween!  (reverse-complement)\n");

      fprintf(stderr, "iid1 "u32bitFMT", iid2 "u32bitFMT"\n", iid1, iid2);
      ma->print(stderr, "A", "B");
      mb->print(stderr, "A", "B");

      fprintf(stderr, "before, iid1, after\n");
      MOB[MOB.index(iid1)-1]->print(stderr, "A", "B");
      MOB[MOB.index(iid1)  ]->print(stderr, "A", "B");
      MOB[MOB.index(iid1)+1]->print(stderr, "A", "B");

      fprintf(stderr, "before, iid2, after\n");
      MOB[MOB.index(iid2)-1]->print(stderr, "A", "B");
      MOB[MOB.index(iid2)  ]->print(stderr, "A", "B");
      MOB[MOB.index(iid2)+1]->print(stderr, "A", "B");

      return(false);
    }
  }

  return(true);
}



void
dumpAgap(atacMatch *ma, atacMatch *mb,
         atacMatchOrder &MOB,
         FastACache *C1, FastACache *C2,
         u32bit gapLimit,
         bool shiftRight) {
}


void
dumpBgap(atacMatch *ma, atacMatch *mb,
         atacMatchOrder &MOB,
         FastACache *C1, FastACache *C2,
         u32bit gapLimit,
         bool shiftRight) {
}


//  Returns the beginning of the sequence from pos to pos+len
char *
getSequenceBeg(char *str, u32bit pos, u32bit len, FastAAccessor &it) {
  u32bit i = 0;

  it.setRange(pos, len);

  if (len > MAXPRINT)
    len = MAXPRINT;

  if (len > 0) {
    it.setPosition(pos);

    for (i=0; i<len; i++, ++it)
      str[i] = *it;
  }
  str[i] = 0;

  return(str);
}


//  Returns all sequence from pos to pos+len
char *
getSequenceAll(char *str, u32bit pos, u32bit len, FastAAccessor &it) {
  u32bit i = 0;

  it.setRange(pos, len);

  if (len > 0) {
    it.setPosition(pos);

    for (i=0; i<len; i++, ++it)
      str[i] = *it;
  }
  str[i] = 0;

  return(str);
}


//  Returns the end of the sequence from pos to pos+len
char *
getSequenceEnd(char *str, u32bit pos, u32bit len, FastAAccessor &it) {
  u32bit i = 0;

  it.setRange(pos, len);

  pos += len;
  if (len > MAXPRINT)
    len = MAXPRINT;

  if (len > 0) {
    it.setPosition(pos - len);

    for (i=0; i<len; i++, ++it)
      str[i] = *it;
  }
  str[i] = 0;

  return(str);
}



u32bit
shiftGap(atacFile &AF,
         atacMatchList &ML,
         atacMatch *ma, atacMatch *mb,
         atacMatchOrder &MOB,
         FastACache *C1, FastACache *C2,
         u32bit gapLimit,
         bool shiftRight) {

#ifdef REPORT_RESULTS
  fprintf(stderr, "----------------------------------------\n");
#endif

  if (isPotentiallyShiftable(ma, mb, MOB, gapLimit) == false)
    return(0);

  //  Save a copy of the original matches
  atacMatch macopy = *ma;
  atacMatch mbcopy = *mb;

  //  Grab the sequences we use, and make the accesors
  //
  seqInCore  *s1 = C1->getSequenceInCore(ma->iid1);
  seqInCore  *s2 = C2->getSequenceInCore(ma->iid2);

  FastAAccessor mas1(s1, ma->fwd1 == false);
  FastAAccessor mas2(s2, ma->fwd2 == false);
  FastAAccessor mbs1(s1, mb->fwd1 == false);
  FastAAccessor mbs2(s2, mb->fwd2 == false);

  mas1.setRange(ma->pos1, ma->len1);
  mas2.setRange(ma->pos2, ma->len2);
  mbs1.setRange(mb->pos1, mb->len1);
  mbs2.setRange(mb->pos2, mb->len2);

  u32bit  shifted = 0;

  //  We want to extend ma to the right, this will shift the gap to
  //  the right-most position (relative to the forward genomic).
  //
  //  While there is a match after ma, extend ma to the right, and
  //  decrease mb from the left.
  //
  if (shiftRight == false) {

    //  Similar to above.  The accessor hides most of the pain caused
    //  by reverse complement.

    mas1.setPosition(ma->pos1 + ma->len1 - 1);
    mas2.setPosition(ma->pos2 + ma->len2 - 1);

    mbs1.setPosition(mb->pos1);  --mbs1;
    mbs2.setPosition(mb->pos2);  --mbs2;


#ifdef REPORT_DEBUG
    //  Dump out some sequence to see where we really are
    //
    fprintf(stderr, "A: ");
    for (u32bit i=0; i<50; i++) {
      fprintf(stderr, "%c", *mas1);
      --mas1;
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "B: ");
    for (u32bit i=0; i<50; i++) {
      fprintf(stderr, "%c", *mbs1);
      --mbs1;
    }
    fprintf(stderr, "\n");

    //  Reset the iterators
    //
    mas1.setPosition(ma->pos1 + ma->len1 - 1);
    mas2.setPosition(ma->pos2 + ma->len2 - 1);

    mbs1.setPosition(mb->pos1);  --mbs1;
    mbs2.setPosition(mb->pos2);  --mbs2;
#endif



    while (mas1.isValid() &&
           mas2.isValid() &&
           mbs1.isValid() &&
           mbs2.isValid() &&
           (ma->len1 > 0) &&
           (ma->len2 > 0) &&
           validSymbol[(int)*mbs1] &&
           validSymbol[(int)*mbs2] &&
           IUPACidentity[(int)*mbs1][(int)*mbs2]) {

#ifdef REPORT_SHIFTING
      fprintf(stderr, "EXTENDrev:  MA %c/%c ----- %c/%c MB\n",
              *mas1, *mas2, *mbs1, *mbs2);
#endif

      mas1.extendRight(-1);  ma->len1--;  --mas1;
      mas2.extendRight(-1);  ma->len2--;  --mas2;

      mbs1.extendLeft(1);  mb->len1++;  --mbs1;
      mbs2.extendLeft(1);  mb->len2++;  --mbs2;

      shifted++;
    }
  } else {

    //  A wants to be the first thing after ma -- the first base in
    //  the gap.  Set the position to the last thing in the range,
    //  then use the increment operator to extend past that.  The spec
    //  on FastAAccessor says we can't directly go somewhere outside
    //  the range.
    //
    mas1.setPosition(ma->pos1 + ma->len1 - 1);  ++mas1;
    mas2.setPosition(ma->pos2 + ma->len2 - 1);  ++mas2;

    //  B can be set to the first thing in the match with no problem.
    //
    mbs1.setPosition(mb->pos1);
    mbs2.setPosition(mb->pos2);

    //  While we're still in sequence (isValid()) and we haven't
    //  obliterated the match we're shifting the gap into, and we can
    //  extend the other match (being both validSymbols and identity),
    //  shift the gap to the right.
    //
    while (mas1.isValid() &&
           mas2.isValid() &&
           mbs1.isValid() &&
           mbs2.isValid() &&
           (mb->len1 > 0) &&
           (mb->len2 > 0) &&
           validSymbol[(int)*mas1] &&
           validSymbol[(int)*mas2] &&
           IUPACidentity[(int)*mas1][(int)*mas2]) {

#ifdef REPORT_SHIFTING
      fprintf(stderr, "EXTENDfwd:  MA %c/%c ----- %c/%c MB\n",
              *mas1, *mas2, *mbs1, *mbs2);
#endif

      mas1.extendRight(1);  ma->len1++;  ++mas1;
      mas2.extendRight(1);  ma->len2++;  ++mas2;

      mbs1.extendLeft(-1);  mb->len1--;  ++mbs1;
      mbs2.extendLeft(-1);  mb->len2--;  ++mbs2;

      shifted++;
    }
  }

  //  Finally, update the two matches with the shifted results.

  ma->pos1 = mas1.getRangeBegin();
  ma->len1 = mas1.getRangeLength();
  ma->pos2 = mas2.getRangeBegin();
  ma->len2 = mas2.getRangeLength();
  mb->pos1 = mbs1.getRangeBegin();
  mb->len1 = mbs1.getRangeLength();
  mb->pos2 = mbs2.getRangeBegin();
  mb->len2 = mbs2.getRangeLength();


  //
  //  The rest is just error checking.
  //

  if (shifted)
    numShifted++;
  else
    numNotShifted++;

  if (shifted < 1024)
    amountShifted[shifted]++;

  //  leftmatch origend newend rightmatch origbegin newbegin
  if (shifted && logFile) {
    if (ma->fwd2) {
      //  Forward matches are easy.
      //
      fprintf(logFile, "%s\t%s\t%s:"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t->\t"u32bitFMT"\t"u32bitFMT"\t",
              ma->matchuid, mb->matchuid,
              AF.labelA(), ma->iid1,
              macopy.pos1 + macopy.len1, ma->pos1 + ma->len1,
              mbcopy.pos1, mb->pos1);
      fprintf(logFile, "%s:"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t->\t"u32bitFMT"\t"u32bitFMT"\n",
              AF.labelB(), ma->iid2,
              macopy.pos2 + macopy.len2, ma->pos2 + ma->len2,
              mbcopy.pos2, mb->pos2);
    } else {
      //  Reverse matches are painful.  The gap on B is between the
      //  right edge of mb, and the the left edge of ma.
      //
      fprintf(logFile, "%s\t%s\t%s:"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t->\t"u32bitFMT"\t"u32bitFMT"\t",
              ma->matchuid, mb->matchuid,
              AF.labelA(), ma->iid1,
              macopy.pos1 + macopy.len1, ma->pos1 + ma->len1,
              mbcopy.pos1, mb->pos1);
      fprintf(logFile, "%s:"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t->\t"u32bitFMT"\t"u32bitFMT"\n",
              AF.labelB(), ma->iid2,
              mbcopy.pos2 + mbcopy.len2, mb->pos2 + mb->len2,
              macopy.pos2, ma->pos2);
    }
  }


#ifdef REPORT_RESULTS
  if (shifted)
    fprintf(stderr, "SHIFTED "u32bitFMT" bases.\n", shifted);
  else
    fprintf(stderr, "NOT SHIFTED.\n");

  fprintf(stderr, u32bitFMTW(9)"-"u32bitFMTW(9)" -- "u32bitFMTW(9)"-"u32bitFMTW(9)"\n",
          macopy.pos1, macopy.pos1 + macopy.len1,
          mbcopy.pos1, mbcopy.pos1 + mbcopy.len1);
  fprintf(stderr, u32bitFMTW(9)"-"u32bitFMTW(9)" -- "u32bitFMTW(9)"-"u32bitFMTW(9)"\n",
          macopy.pos2, macopy.pos2 + macopy.len2,
          mbcopy.pos2, mbcopy.pos2 + mbcopy.len2);
  fprintf(stderr, "shifted "u32bitFMT" bases (fwd1=%d fwd2=%d  fwd1=%d fwd2=%d)\n", shifted, ma->fwd1, ma->fwd2, mb->fwd1, mb->fwd2);
  fprintf(stderr, u32bitFMTW(9)"-"u32bitFMTW(9)" -- "u32bitFMTW(9)"-"u32bitFMTW(9)"\n",
          ma->pos1, ma->pos1 + ma->len1,
          mb->pos1, mb->pos1 + mb->len1);
  fprintf(stderr, u32bitFMTW(9)"-"u32bitFMTW(9)" -- "u32bitFMTW(9)"-"u32bitFMTW(9)"\n",
          ma->pos2, ma->pos2 + ma->len2,
          mb->pos2, mb->pos2 + mb->len2);
#endif

  u32bit errors = 0;

  if (macopy.pos1 != ma->pos1)
    fprintf(stderr, "WARNING:  begin of assembly 1 moved!\n"), errors++;
  if (mbcopy.pos1 + mbcopy.len1 != mb->pos1 + mb->len1)
    fprintf(stderr, "WARNING:  end of assembly 1 moved!\n"), errors++;

  if ((ma->fwd2 == true)  && (macopy.pos2 != ma->pos2))
    fprintf(stderr, "WARNING:  begin of assembly 2 moved!\n"), errors++;
  if ((ma->fwd2 == false) && (mbcopy.pos2 != mb->pos2))
    fprintf(stderr, "WARNING:  begin of assembly 2 moved (rc)!\n"), errors++;

  if ((ma->fwd2 == true)  && (mbcopy.pos2 + mbcopy.len2 != mb->pos2 + mb->len2))
    fprintf(stderr, "WARNING:  end of assembly 1 moved!\n"), errors++;
  if ((ma->fwd2 == false) && (macopy.pos2 + macopy.len2 != ma->pos2 + ma->len2))
    fprintf(stderr, "WARNING:  end of assembly 2 moved (rc)!\n"), errors++;

  //  For debugging, claim there were errors if we shifted something.
#ifdef REPORT_RESULTS
  errors++;
#endif

  if (errors > 0) {
    atacMatch *l = 0L;
    atacMatch *r = 0L;
    char    str1[1024];
    char    str2[1024];
    char    str3[1024];

    //  Print the sequence.  We could print each piece separately
    //  (and, indeed, we tried that initially) but that's difficult
    //  because we need to remember which match is first on B, 

    macopy.print(stderr, "A", "B");
    mbcopy.print(stderr, "A", "B");
    l = &macopy;
    r = &mbcopy;

    getSequenceEnd(str1, l->pos1, l->len1, mas1);
    getSequenceAll(str2, l->pos1 + l->len1, r->pos1 - l->pos1 - l->len1, mas1);
    getSequenceBeg(str3, r->pos1, r->len1, mas1);
    fprintf(stderr, "SEQA: %s -- %s -- %s\n", str1, str2, str3);

    if (macopy.fwd2) {
      //  We're forward, so l is really first on B.
      getSequenceEnd(str1, l->pos2, l->len2, mas2);
      getSequenceAll(str2, l->pos2 + l->len2, r->pos2 - l->pos2 - l->len2, mas2);
      getSequenceBeg(str3, r->pos2, r->len2, mas2);
    } else {
      //  Nope, reverse complement, so r is really first on B.  This
      //  only changes how we get the gap.
      //
      getSequenceEnd(str1, l->pos2, l->len2, mas2);
      getSequenceAll(str2, r->pos2 + r->len2, l->pos2 - r->pos2 - r->len2, mas2);
      getSequenceBeg(str3, r->pos2, r->len2, mas2);
    }
    fprintf(stderr, "SEQB: %s -- %s -- %s\n", str1, str2, str3);


    //  Do the same thing (same getSequence calls) for the after picture.

    ma->print(stderr, "A", "B");
    mb->print(stderr, "A", "B");
    l = ma;
    r = mb;

    getSequenceEnd(str1, l->pos1, l->len1, mas1);
    getSequenceAll(str2, l->pos1 + l->len1, r->pos1 - l->pos1 - l->len1, mas1);
    getSequenceBeg(str3, r->pos1, r->len1, mas1);
    fprintf(stderr, "SEQA: %s -- %s -- %s\n", str1, str2, str3);

    if (ma->fwd2) {
      getSequenceEnd(str1, l->pos2, l->len2, mas2);
      getSequenceAll(str2, l->pos2 + l->len2, r->pos2 - l->pos2 - l->len2, mas2);
      getSequenceBeg(str3, r->pos2, r->len2, mas2);
    } else {
      getSequenceEnd(str1, l->pos2, l->len2, mas2);
      getSequenceAll(str2, r->pos2 + r->len2, l->pos2 - r->pos2 - r->len2, mas2);
      getSequenceBeg(str3, r->pos2, r->len2, mas2);
    }
    fprintf(stderr, "SEQB: %s -- %s -- %s\n", str1, str2, str3);
  }

  //if (errors)
  //  exit(1);

  return(shifted);
}





int
main(int argc, char *argv[]) {

  if (argc == 1) {
    fprintf(stderr, "usage: %s [options] < matches > matches\n", argv[0]);
    fprintf(stderr, "  Instead of the usual switch based options to enable behavior\n");
    fprintf(stderr, "  gapShifter iterates of a list of shift directions and sizes.\n");
    fprintf(stderr, "    l      -- shift gaps to the left\n");
    fprintf(stderr, "    r      -- shift gaps to the right\n");
    fprintf(stderr, "    #      -- set the maximum size of a gap to shift\n");
    fprintf(stderr, "    log x  -- open a logfile 'x' for results of the next shift\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  for example\n");
    fprintf(stderr, "    gapShifter 1 l r l r 10 l r l log X r < some.atac > shifted.atac\n");
    fprintf(stderr, "  would shift 1bp gaps to the left, then to the right, then left,\n");
    fprintf(stderr, "  then set the gap size to 10bp and repeat.  The last shift is logged\n");
    fprintf(stderr, "  into fle 'X'.\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  This is useful since shifting gaps can obliterate matches, but possibly.\n");
    fprintf(stderr, "  when both left and right shifts are used.\n");
    fprintf(stderr, "      GCTAATTAGACG\n");
    fprintf(stderr, "      GCT-AT-AGACG\n");
    fprintf(stderr, "  The second gap can be shifted to the left, and the first gap can be\n");
    fprintf(stderr, "  shifted right, resulting in\n");
    fprintf(stderr, "      GCTAATTAGACG\n");
    fprintf(stderr, "      GCTA--TAGACG\n");
    fprintf(stderr, "  Thus, two one base gaps were merged into a two base gap, which might\n");
    fprintf(stderr, "  then be able to be shifted.  e.g.:\n");
    fprintf(stderr, "      atgatcatcttatc\n");
    fprintf(stderr, "      at---c-t--tatc\n");
    exit(1);
  }

  atacFile        AF("-");
  atacMatchList  &ML = *AF.matches();
  atacMatchOrder  MOA(ML);
  atacMatchOrder  MOB(ML);

  MOA.sortA();
  MOB.sortB();

  //  second to last == loadAll
  //  last           == report loading
  //
  FastACache  *C1 = new FastACache(AF.assemblyFileA(),    2, false, false);
  FastACache  *C2 = new FastACache(AF.assemblyFileB(), 1024, false, false);

  bool    shiftRight   = true;
  u32bit  gapLimit     = 5;

  char   *logFileName  = 0L;

  int arg=1;
  while (arg < argc) {
    bool doShift = false;

    if        (strcmp(argv[arg], "log") == 0) {
      logFileName = argv[++arg];
      errno = 0;
      logFile = fopen(logFileName, "w");
      if (errno)
        fprintf(stderr, "gapShifter: can't open log file '%s': %s\n", logFileName, strerror(errno)), exit(1);
    } else if (strcmp(argv[arg], "l") == 0) {
      shiftRight   = false;
      doShift      = true;
    } else if (strcmp(argv[arg], "r") == 0) {
      shiftRight   = true;
      doShift      = true;
    } else {
      gapLimit     = strtou32bit(argv[arg], 0L);
    }

    if (doShift) {

      for (u32bit x=0; x<1024; x++)
        amountShifted[x] = 0;

      numShifted = 0;
      numNotShifted = 0;
      numDiffSeq = 0;
      numDiffOri = 0;
      numZeroLen = 0;
      numOutOfOrder = 0;
      numNotAdjacent = 0;
      numNoGap = 0;
      numGapTooBig = 0;
      numOverlapping = 0;

      fprintf(stderr, "Shifting gaps of length at most "u32bitFMT" bases, to the %s.\n", gapLimit, (shiftRight) ? "right" : "left");

      u32bit gapsShifted = 0;
      for (u32bit i=1; i<ML.numMatches(); i++) {

        if (shiftGap(AF, ML, MOA[i-1], MOA[i], MOB, C1, C2, gapLimit, shiftRight)) {
          gapsShifted++;
          //fprintf(stderr, "shifted "u32bitFMT" out of "u32bitFMT" (%6.2f%%)\r", gapsShifted, i, (double)gapsShifted / (double)i * 100.0);
          //fflush(stderr);
        }
      }

      fprintf(stderr, "numShifted = "u32bitFMT"\n", numShifted);
      fprintf(stderr, "numNotShifted = "u32bitFMT"\n", numNotShifted);
      fprintf(stderr, "numDiffSeq = "u32bitFMT"\n", numDiffSeq);
      fprintf(stderr, "numDiffOri = "u32bitFMT"\n", numDiffOri);
      fprintf(stderr, "numZeroLen = "u32bitFMT"\n", numZeroLen);
      fprintf(stderr, "numOutOfOrder = "u32bitFMT"\n", numOutOfOrder);
      fprintf(stderr, "numNotAdjacent = "u32bitFMT"\n", numNotAdjacent);
      fprintf(stderr, "numNoGap = "u32bitFMT"\n", numNoGap);
      fprintf(stderr, "numGapTooBig = "u32bitFMT"\n", numGapTooBig);
      fprintf(stderr, "numOverlapping = "u32bitFMT"\n", numOverlapping);

      for (u32bit x=0; x<50; x++)
        fprintf(stderr, "amountShifted["u32bitFMT"] = "u32bitFMT" (number of gaps shifted by [number of bases])\n", x, amountShifted[x]);

      fprintf(stderr, "shifted "u32bitFMT" out of "u32bitFMT" (%6.2f%%)\n", gapsShifted, ML.numMatches(), (double)gapsShifted / (double)ML.numMatches() * 100.0);

      if (logFile) {
        fclose(logFile);
        logFileName = 0L;
        logFile     = 0L;
      }
    }

    arg++;
  }


  for (u32bit i=0; i<ML.numMatches(); i++) {
    atacMatch *ma = ML[i];
    if ((ma->len1 > 0) && (ma->len2 > 0))
      fprintf(stdout, "M u %s %s %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" 1  %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
              ma->matchuid, ma->parentuid,
              AF.labelA(), ma->iid1, ma->pos1, ma->len1,
              AF.labelB(), ma->iid2, ma->pos2, ma->len2, ma->fwd2 ? 1 : -1);
  }

  return(0);
}

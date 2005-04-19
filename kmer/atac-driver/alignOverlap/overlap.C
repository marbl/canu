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
#include <time.h>

#include "bio++.H"
#include "util++.H"
#include "atac-common.H"

int sortMatches1(const void *a, const void *b);
int sortMatches2(const void *a, const void *b);
int spanCompare(const void *a, const void *b);

#include "overlap-match.H"
#include "overlap-matchList.H"
#include "overlap-span.H"
#include "overlap-matchTree.H"
#include "overlap-spanTree.H"
#include "overlap-annoList.H"

#include "overlap-sort.C"

u32bit      ALmax = 0;
u32bit      ALlen = 0;
annoList   *AL    = 0L;


void
printAnno(FILE *F,
          char label,
          u32bit axis,
          span_t *span,
          u32bit match1=u32bitZERO, match_t *m1=0L,
          u32bit match2=u32bitZERO, match_t *m2=0L) {

  //  If we're just given match1, make it match2 if it is the second mapping
  //
  if ((match1 >> COLORSHIFT) && (match2 == u32bitZERO)) {
    match2 = match1; m2 = m1;
    match1 = 0;      m1 = 0;
  }

  u32bit len  = span->_end - span->_beg;

  //  axis is 1 or 2; if we're the first axis (B35 centric) make a
  //  list of the matches for later processing

  if (axis == 1)
    AL[ALlen++].add(label, span->_iid, span->_beg, len,
                    match1 & COLORMASK, m1,
                    match2 & COLORMASK, m2);

  fprintf(F, "%c "u32bitFMTW(4)":"u32bitFMTW(09)"-"u32bitFMTW(09)"["u32bitFMTW(6)"] ",
          label,
          span->_iid, span->_beg, span->_end, len);

  if (m1) {
    fprintf(F, u32bitFMTW(07)" ", m1->matchuid);
    u32bit off1 = span->_beg - m1->pos1;

    if (axis == 1) {
      u32bit  sta = m1->pos2 + off1;
      u32bit  end = m1->pos2 + off1 + len;

      if (m1->ori2 == 0) {
        sta = m1->pos2 + m1->len2 - off1;
        end = m1->pos2 + m1->len2 - off1 - len;
      }

      fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", m1->iid2, sta, end);
    } else {
      fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", m1->iid1, m1->pos1 + off1, m1->pos1 + off1 + len);
    }
  } else {
    fprintf(F, u32bitFMTW(07)" ", u32bitZERO);
    fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", u32bitZERO, u32bitZERO, u32bitZERO);
  }

  if (m2) {
    fprintf(F, u32bitFMTW(07)" ", m2->matchuid);
    u32bit off2 = span->_beg - m2->pos1;

    if (axis == 1) {
      u32bit  sta = m2->pos2 + off2;
      u32bit  end = m2->pos2 + off2 + len;
      
      if (m2->ori2 == 0) {
        sta = m2->pos2 + m2->len2 - off2;
        end = m2->pos2 + m2->len2 - off2 - len;
      }

      fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", m2->iid2, sta, end);
    } else {
      fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", m2->iid1, m2->pos1 + off2, m2->pos1 + off2 + len);
    }
  } else {
    fprintf(F, u32bitFMTW(07)" ", u32bitZERO);
    fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", u32bitZERO, u32bitZERO, u32bitZERO);
  }

  fprintf(F, "\n");
}


#include "overlap-find.C"



int
main(int argc, char **argv) {

  if (argc != 4) {
    fprintf(stderr, "usage: %s <matches-1> <matches-2> <out-prefix>\n", argv[0]);
    exit(1);
  }

  matchList  *M1 = new matchList(argv[1]);
  matchList  *M2 = new matchList(argv[2]);
  char       *OP = argv[3];

  //  Construct a pair of matchTrees for one of the maps.  We need
  //  to search by either side of the match, thus two trees.
  //
  //  Our trees have pointers into the matchList.  We make a temporary
  //  list of pointers so we can sort for loading into the tree.
  //
  //matchTree  *T1a = new matchTree(M1, 0);
  //matchTree  *T1b = new matchTree(M1, 1);


  //  We want to annotate the two assembies with:
  //    a) mapped by both, the same
  //    b) mapped by both, differently
  //    c) mapped by the first, unmapped by the second
  //    d) mapped by the second, unmapped by the first
  //    e) unmapped by both
  //
  //  If unmapped, we could further annotate with the reason it was
  //  unmapped -- not found, or found multiple times.
  //
  //  Our annotation datastructure is a tree of spans.  Each span is a
  //  sequence, and an interval on that sequence.  We assume that the
  //  tree contains the spans for the whole sequence, that is, that we
  //  never need to increase a span, just split.
  //
  spanTree    *S1 = new spanTree();
  spanTree    *S2 = new spanTree();

  //  Initialize the tree of spans by inserting a single span for each
  //  sequence in the file.
  //
  for (u32bit i=0; i<M1->_seq1->getNumberOfSequences(); i++)
    S1->addNewSpan(i, M1->_seq1->sequenceLength(i));
  for (u32bit i=0; i<M1->_seq2->getNumberOfSequences(); i++)
    S2->addNewSpan(i, M1->_seq2->sequenceLength(i));

  //  Add every match to the spanTrees.

  for (u32bit i=0; i<M1->_matchesLen; i++) {
    S1->addMatch(M1->_matches+i, 0, 0);
    S2->addMatch(M1->_matches+i, 1, 0);
  }
  for (u32bit i=0; i<M2->_matchesLen; i++) {
    S1->addMatch(M2->_matches+i, 0, 1);
    S2->addMatch(M2->_matches+i, 1, 1);
  }

  //  Dump each spanTree: For each span, we need to check that
  //    it has matches?
  //    only one match, or only matches from one mapping?
  //    matches from both mappings?  need to check that
  //     the span in the other tree also has the same matches

  //  Statistics and Histograms
  //
  //  Index 1 is the assembly, 2 is the mapping.  Stats count the
  //  number of bases covered, histograms are of the block sizes.
  //
  u32bit   unmapped[2]     = {0,0};
  u32bit   unique[2][2]    = {{0,0},{0,0}};
  u32bit   different[2]    = {0,0};
  u32bit   wilddiff[2]     = {0,0};
  u32bit   same[2]         = {0,0};
  u32bit   inconsistent[2] = {0,0};

  //  Histograms.  
  //
  u32bit   histogramMax = 1048576;
  u32bit  *unmappedH[2];
  u32bit  *uniqueH[2][2];
  u32bit  *differentH[2];
  u32bit  *wilddiffH[2];
  u32bit  *sameH[2];
  u32bit  *inconsistentH[2];

  for (u32bit i=0; i<2; i++) {
    unmappedH[i]     = new u32bit [histogramMax];
    uniqueH[i][0]    = new u32bit [histogramMax];
    uniqueH[i][1]    = new u32bit [histogramMax];
    differentH[i]    = new u32bit [histogramMax];
    wilddiffH[i]     = new u32bit [histogramMax];
    sameH[i]         = new u32bit [histogramMax];
    inconsistentH[i] = new u32bit [histogramMax];

    bzero(unmappedH[i],     sizeof(u32bit) * histogramMax);
    bzero(uniqueH[i][0],    sizeof(u32bit) * histogramMax);
    bzero(uniqueH[i][1],    sizeof(u32bit) * histogramMax);
    bzero(differentH[i],    sizeof(u32bit) * histogramMax);
    bzero(wilddiffH[i],     sizeof(u32bit) * histogramMax);
    bzero(sameH[i],         sizeof(u32bit) * histogramMax);
    bzero(inconsistentH[i], sizeof(u32bit) * histogramMax);
  }

#define updateHistogram(H,x)  H[ ((x) > histogramMax) ? 0 : x ]++

  //  Stats, debugging mostly
  u32bit   same1sf=0, same1ef=0, same1sn=0, same1en=0;
  u32bit   same2sf=0, same2ef=0, same2sn=0, same2en=0;
  u32bit   same1f=0, same1n=0;
  u32bit   same2f=0, same2n=0;

  //  We write the annotation to here
  char  outname[1024];
  sprintf(outname, "%s.Aannotation", OP);
  FILE *Aanno = fopen(outname, "w");
  sprintf(outname, "%s.Bannotation", OP);
  FILE *Banno = fopen(outname, "w");

  //  We also build a list of the annotation for classification later
  ALmax = (u32bit)dict_count(S1->_tree);
  ALlen = 0;
  AL    = new annoList [ ALmax ];


  //  Doesn't handle weird stuff like this span (on sequence 1)
  //  mapping onto seq2 correctly, but the span in seq2 having an
  //  extra match to somewhere else in seq1.

  dnode_t  *node1 = 0L;
  dnode_t  *node2 = 0L;

  //  we want to find the single span in the other spanTree that
  //  corresponds to this span.  once we do that, we can verify that
  //  all the matches are the same.
  //
  //  because we are gapless matches, we can, for each match,
  //  compute the exact location this span should occur on the other
  //  sequence.  then, do a lookup() to get that span, or just
  //  verify that everybody is the same location.

  node1 = dict_first(S1->_tree);
  while (node1) {
    span_t  *span = (span_t *)dnode_getkey(node1);
    u32bit   spanLen = span->_end - span->_beg;

    if (span->_matchesLen == 0) {
      unmapped[0] += spanLen;
      updateHistogram(unmappedH[0], spanLen);

      printAnno(Aanno, 'U', 1, span);
    } else if (span->_matchesLen == 1) {
      u32bit idx = span->_matches[0] >> COLORSHIFT;
      unique[0][idx] += spanLen;
      updateHistogram(uniqueH[0][idx], spanLen);

      u32bit    match = span->_matches[0];
      match_t  *m;

      if (match >> COLORSHIFT) {
        m = &M2->_matches[match & COLORMASK];
      } else {
        m = &M1->_matches[match & COLORMASK];
      }

      printAnno(Aanno, '1', 1, span, match, m);
    } else if ((span->_matchesLen == 2) &&
               ((span->_matches[0] >> COLORSHIFT) == (span->_matches[1] >> COLORSHIFT))) {
      inconsistent[0] += spanLen;

      printAnno(Aanno, '?', 1, span);
    } else if (span->_matchesLen == 2) {
      u32bit match1 = span->_matches[0];
      u32bit match2 = span->_matches[1];

      if (match1 >> COLORSHIFT) {
        match1 = span->_matches[1];
        match2 = span->_matches[0];
      }

      match_t  *m1 = &M1->_matches[match1 & COLORMASK];
      match_t  *m2 = &M2->_matches[match2 & COLORMASK];

#if 0
      fprintf(Aanno, "m1: "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
              m1->matchuid,
              m1->iid1, m1->pos1, m1->pos1+m1->len1,
              m1->iid2, m1->pos2, m1->pos2+m1->len2);
      fprintf(Aanno, "m2: "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
              m2->matchuid,
              m2->iid1, m2->pos1, m2->pos1+m2->len1,
              m2->iid2, m2->pos2, m2->pos2+m2->len2);
#endif

      if (m1->iid2 == m2->iid2) {
        u32bit off1 = 0, pos1l = 0, pos1r = 0;
        u32bit off2 = 0, pos2l = 0, pos2r = 0;

        //  If we're a forward match, compare that the begin position of the
        //  match on the other assembly is the same.
        //
        //  But, if we're a reverse match, we need to compare the
        //  'end' location for a match
        //
        off1  = span->_beg - m1->pos1;
        pos1l = m1->pos2 + off1;
        pos1r = m1->pos2 + m1->len2 - off1;

        off2  = span->_beg - m2->pos1;
        pos2l = m2->pos2 + off2;
        pos2r = m2->pos2 + m2->len2 - off2;

        //  We matched both
        //  We matched the start, and we're normal
        //  We matched the end, and we're normal
        if (m1->ori2 == 0) {
          if ((pos1l == pos2l) && (pos1r == pos2r)) same1n++;
          else if (pos1l == pos2l) same1sn++;
          else if (pos1r == pos2r) same1en++;
        } else {
          if ((pos1l == pos2l) && (pos1r == pos2r)) same1f++;
          else if (pos1l == pos2l) same1sf++;
          else if (pos1r == pos2r) same1ef++;
        }

        if ((pos1l == pos2l) || (pos1r == pos2r)) {
          same[0] += spanLen;
          updateHistogram(sameH[0], spanLen);

          printAnno(Aanno, 'Y', 1, span, match1, m1, match2, m2);
        } else {
          different[0] += spanLen;
          updateHistogram(differentH[0], spanLen);

          printAnno(Aanno, 'N', 1, span, match1, m1, match2, m2);
        }
      } else {
        //  Wildly different matches!  Mapped to different scaffolds!
        wilddiff[0] += spanLen;
        updateHistogram(wilddiffH[0], spanLen);
        printAnno(Aanno, '!', 1, span, match1, m1, match2, m2);
      }
    } else {
      inconsistent[0] += spanLen;
      updateHistogram(inconsistentH[0], spanLen);

      printAnno(Aanno, '?', 1, span);
    }

    node1 = dict_next(S1->_tree, node1);
  }

#if 0
  node2 = dict_first(S2->_tree);
  while (node2) {
    span_t  *span = (span_t *)dnode_getkey(node2);
    u32bit   spanLen = span->_end - span->_beg;

    if (span->_matchesLen == 0) {
      unmapped[1] += spanLen;
      updateHistogram(unmappedH[1], spanLen);

      printAnno(Banno, 'U', 2, span);
    } else if (span->_matchesLen == 1) {
      u32bit idx = span->_matches[0] >> COLORSHIFT;
      unique[1][idx] += spanLen;
      updateHistogram(uniqueH[1][idx], spanLen);

      u32bit    match = span->_matches[0];
      match_t  *m;

      if (match >> COLORSHIFT) {
        m = &M2->_matches[match & COLORMASK];
      } else {
        m = &M1->_matches[match & COLORMASK];
      }

      printAnno(Banno, '1', 2, span, match, m);
    } else if ((span->_matchesLen == 2) &&
               ((span->_matches[0] >> COLORSHIFT) == (span->_matches[1] >> COLORSHIFT))) {
      inconsistent[1] += spanLen;

      printAnno(Banno, '?', 2, span);
    } else if (span->_matchesLen == 2) {
      u32bit match1 = span->_matches[0];
      u32bit match2 = span->_matches[1];

      if (match1 >> COLORSHIFT) {
        match1 = span->_matches[1];
        match2 = span->_matches[0];
      }

      match1 &= COLORMASK;
      match2 &= COLORMASK;

      match_t  *m1 = &M1->_matches[match1];
      match_t  *m2 = &M2->_matches[match2];

      if (m1->iid1 == m2->iid1) {
        u32bit off1 = 0, pos1l = 0, pos1r = 0;
        u32bit off2 = 0, pos2l = 0, pos2r = 0;

        //  If we're a forward match, compare that the begin position of the
        //  match on the other assembly is the same.
        //
        //  But, if we're a reverse match, we need to compare the
        //  'end' location for a match
        //
        off1  = span->_beg - m1->pos2;
        pos1l = m1->pos1 + off1;
        pos1r = m1->pos1 + m1->len1 - off1;

        off2  = span->_beg - m2->pos2;
        pos2l = m2->pos1 + off2;
        pos2r = m2->pos1 + m2->len1 - off2;

        //  We matched both
        //  We matched the start, and we're normal
        //  We matched the end, and we're normal
        if (m1->ori2 == 0) {
          if ((pos1l == pos2l) && (pos1r == pos2r)) same2n++;
          else if (pos1l == pos2l) same2sn++;
          else if (pos1r == pos2r) same2en++;
        } else {
          if ((pos1l == pos2l) && (pos1r == pos2r)) same2f++;
          else if (pos1l == pos2l) same2sf++;
          else if (pos1r == pos2r) same2ef++;
        }


        if ((pos1l == pos2l) || (pos1r == pos2r)) {
          same[1] += spanLen;
          updateHistogram(sameH[1], spanLen);

          printAnno(Banno, 'Y', 2, span, match1, m1, match2, m2);
        } else {
          different[1] += spanLen;
          updateHistogram(differentH[1], spanLen);

          printAnno(Banno, 'N', 2, span, match1, m1, match2, m2);
        }
      } else {
        //  Wildly different matches!  Mapped to different scaffolds!
        wilddiff[1] += spanLen;
        updateHistogram(wilddiffH[1], spanLen);
        printAnno(Aanno, '!', 2, span, match1, m1, match2, m2);
      }
    } else {
      inconsistent[1] += spanLen;
      updateHistogram(inconsistentH[1], spanLen);

      printAnno(Banno, '?', 2, span);
    }

    node2 = dict_next(S2->_tree, node2);
  }
#endif


  fclose(Aanno);
  fclose(Banno);

  fprintf(stderr, "unmapped:           A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", unmapped[0],     unmapped[1]);
  fprintf(stderr, "unique mapping 1:   A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", unique[0][0],    unique[1][0]);
  fprintf(stderr, "unique mapping 2:   A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", unique[0][1],    unique[1][1]);
  fprintf(stderr, "different:          A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", different[0],    different[1]);
  fprintf(stderr, "wild diff:          A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", wilddiff[0],     wilddiff[1]);
  fprintf(stderr, "same:               A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", same[0],         same[1]);
  fprintf(stderr, "inconsistent:       A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", inconsistent[0], inconsistent[1]);

  fprintf(stderr, "\n");
  fprintf(stderr, "same1 both:         normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same1n,  same1f);
  fprintf(stderr, "same1 start:        normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same1sn, same1sf);
  fprintf(stderr, "same1 end:          normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same1en, same1ef);
  fprintf(stderr, "\n");
  fprintf(stderr, "same2 both:         normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same2n,  same2f);
  fprintf(stderr, "same2 start:        normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same2sn, same2sf);
  fprintf(stderr, "same2 end:          normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same2en, same2ef);



  //
  //  run through the annoList, count interesting bits
  //
#if 0
  findIsolatedUnique();
  findExtended();
#endif




#if 0
  FILE *hout = fopen("overlap.histogram", "w");
  for (u32bit i=0; i<histogramMax; i++) {
    fprintf(hout, u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)"\n",
            i,
            unmappedH[0][i], uniqueH[0][0][i], uniqueH[0][1][i], differentH[0][i], wilddiffH[0][i], sameH[0][i], inconsistentH[0][i],
            unmappedH[1][i], uniqueH[1][0][i], uniqueH[1][1][i], differentH[1][i], wilddiffH[1][i], sameH[1][i], inconsistentH[1][i]);
  }
  fclose(hout);
#endif


#if 0

plot [0:500] [0:1500] \
"overlap.histogram" using  2 title "unmapped A" with lines, \
"overlap.histogram" using  9 title "unmapped B" with lines, \
"overlap.histogram" using  5 title "different A" with lines, \
"overlap.histogram" using 12 title "different B" with lines, \
"overlap.histogram" using  6 title "wild diff A" with lines, \
"overlap.histogram" using 13 title "wild diff B" with lines, \
"overlap.histogram" using  8 title "inconsistent A" with lines, \
"overlap.histogram" using 15 title "inconsistent B" with lines

plot [0:500] [0:1500] \
"overlap.histogram" using  7 title "same A" with lines, \
"overlap.histogram" using 14 title "same B" with lines, \

plot [0:500] [0:1500] \
"overlap.histogram" using  3 title "unique A, mapping 1" with lines, \
"overlap.histogram" using  4 title "unique A, mapping 2" with lines, \
"overlap.histogram" using 10 title "unique B, mapping 1" with lines, \
"overlap.histogram" using 11 title "unique B, mapping 2" with lines

#endif


  //  Deleting the spanTrees takes a long time, so we don't bother with any cleanup.
  return(0);
}

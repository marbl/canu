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

//  Given two ATAC-format match files, compute the amount they agree / disagree.

class match_t {
public:
  u32bit  matchid;
  u32bit  iid1, pos1, len1, ori1;
  u32bit  iid2, pos2, len2, ori2;
};




int
sortMatches1(const void *a, const void *b) {
  const match_t *A = *((const match_t * const *)a);
  const match_t *B = *((const match_t * const *)b);

#if 0
  if (debugSort)
    fprintf(stderr, "sortMatches1 "u32bitFMT" "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
            A->iid1, A->pos1, A->len1,
            B->iid1, B->pos1, B->len1);
#endif

  if (A->iid1 < B->iid1)  return(-1);
  if (A->iid1 > B->iid1)  return(1);
  if (A->pos1 < B->pos1)  return(-1);
  if (A->pos1 > B->pos1)  return(1);
  if (A->len1 > B->len1)  return(-1);
  if (A->len1 < B->len1)  return(1);
  if (A->ori1 > B->ori1)  return(-1);
  if (A->ori1 < B->ori1)  return(1);
  return(0);
}

int
sortMatches2(const void *a, const void *b) {
  const match_t *A = *((const match_t * const *)a);
  const match_t *B = *((const match_t * const *)b);

#if 0
  if (debugSort)
    fprintf(stderr, "sortMatches1 "u32bitFMT" "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
            A->iid2, A->pos2, A->len2,
            B->iid2, B->pos2, B->len2);
#endif

  if (A->iid2 < B->iid2)  return(-1);
  if (A->iid2 > B->iid2)  return(1);
  if (A->pos2 < B->pos2)  return(-1);
  if (A->pos2 > B->pos2)  return(1);
  if (A->len2 > B->len2)  return(-1);
  if (A->len2 < B->len2)  return(1);
  if (A->ori2 > B->ori2)  return(-1);
  if (A->ori2 < B->ori2)  return(1);
  return(0);
}




#define COLORSHIFT 24
#define COLORMASK  0x00ffffff

class span_t {
public:
  u32bit  _iid;
  u32bit  _beg;
  u32bit  _end;
  u32bit  _matchesLen;
  u32bit  _matchesMax;
  u32bit *_matches;

  span_t::span_t(u32bit iid, u32bit beg, u32bit end) {
    _iid = iid;
    _beg = beg;
    _end = end;
    _matchesLen = 0;
    _matchesMax = 0;
    _matches    = 0L;
  };

  //  The top X bits of the _matches is for storing the color.  This
  //  does cut down the number of matches we can store.  Human-Human
  //  is ~1 million matches.

  void   addMatch(u32bit matchid, u32bit color) {
    if (_matchesLen >= _matchesMax) {
      if (_matchesMax == 0)
        _matchesMax = 2;
      _matchesMax *= 2;
      u32bit *X = new u32bit [_matchesMax];
      memcpy(X, _matches, sizeof(u32bit) * _matchesLen);
      delete [] _matches;
      _matches = X;
    }

    if (matchid >> COLORSHIFT)
      fprintf(stderr, "ERROR!  span_t::addMatch()-- match id too big, decrease the color space.\n"), exit(1);

    _matches[_matchesLen++] = (color << COLORSHIFT) | (matchid);
  };

  //  Split this span at position, return two new spans
  //
  void   split(u32bit position, span_t* &l, span_t* &r) {
    l = new span_t(_iid, _beg, position);
    r = new span_t(_iid, position, _end);

    l->_matchesLen = _matchesLen;
    l->_matchesMax = _matchesMax;
    l->_matches    = new u32bit [_matchesMax];
    memcpy(l->_matches, _matches, sizeof(u32bit) * _matchesLen);

    r->_matchesLen = _matchesLen;
    r->_matchesMax = _matchesMax;
    r->_matches    = new u32bit [_matchesMax];
    memcpy(r->_matches, _matches, sizeof(u32bit) * _matchesLen);
  };
};



int
spanCompare(const void *a, const void *b) {
  const span_t *A = *((const span_t * const *)a);
  const span_t *B = *((const span_t * const *)b);

  if (A->_iid < B->_iid)  return(-1);
  if (A->_iid > B->_iid)  return(1);
  if (A->_beg < B->_beg)  return(-1);
  if (A->_beg > B->_beg)  return(1);
  return(0);
}




//  Loads a set of matches from a file
class matchList {
public:
  matchList(char *filename);
  ~matchList() {
    delete _seq1;
    delete _seq2;
    delete [] _matches;
  };

  char          _file1[1024];
  char          _file2[1024];

  FastAWrapper *_seq1;
  FastAWrapper *_seq2;

  u32bit        _matchesLen;
  u32bit        _matchesMax;
  match_t      *_matches;
};



//  Contructs a search tree from a matchList
class matchTree {
public:
  matchTree(matchList *L, u32bit side);
  ~matchTree() {
    dict_free_nodes(_tree);
    dict_free(_tree);
  };

  dict_t        *_tree;
  dict_load_t    _load;
};




#include "spanTree.H"


matchTree::matchTree(matchList *L, u32bit side) {

  //  Construct a list of pointers to the matchList data
  //
  //  kazlib was modified to be qsort() compatible and so it passes a
  //  pointer to whatever it is sorting.  Since kazlib operates on
  //  pointers anyway, this means that it passes the compare function
  //  a pointer to a pointer to the object.
  //
  //  Which really fails in this case.  We have a list of pointers to
  //  objects that we sort, then want to load.
  //
  //  Uhhh, no, this is correct.  We give kazlib a pointer to the
  //  object, it gives the compare function a pointer to that pointer.
  //
  //  qsort() below sorts pointers to objects, and does the same.

  match_t  **matchPointers = new match_t * [L->_matchesLen];
  for (u32bit i=0; i<L->_matchesLen; i++)
    matchPointers[i] = L->_matches + i;

  //  Choose a comparison function based on the side we want

  int (*sortMatches)(const void *, const void *) = sortMatches1;
  if (side == 1)
    sortMatches = sortMatches2;
  
  //  Sort

  qsort(matchPointers, L->_matchesLen, sizeof(match_t *), sortMatches);

  //  Load the tree (use DICTCOUNT_T_MAX for max nodes)

  _tree = dict_create(L->_matchesLen, sortMatches);
  dict_allow_dupes(_tree);

  dict_load_begin(&_load, _tree);

  for (u32bit i=0; i<L->_matchesLen; i++) {
#if 0
    fprintf(stderr, "Load "u32bitFMT" "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
            matchPointers[i]->iid1, matchPointers[i]->pos1, matchPointers[i]->len1,
            matchPointers[i]->iid2, matchPointers[i]->pos2, matchPointers[i]->len2);
#endif

    dnode_t   *node = (dnode_t *)malloc(sizeof(dnode_t));
    dnode_init(node, 0L);
    dict_load_next(&_load, node, matchPointers[i]);
  }

  dict_load_end(&_load);

  //  Clean up
  delete [] matchPointers;
}



matchList::matchList(char *filename) {

  errno = 0;
  FILE   *inFile = fopen(filename, "r");
  if (errno)
    fprintf(stderr, "matchList::matchList()-- failed to load %s: %s\n", filename, strerror(errno)), exit(1);


  //  While loading matches, we compute the mapped length and covered
  //  length.

  fprintf(stderr, "Loading matches from %s\n", filename);
  
  //  Read the preamble, look for our data sources.  This leaves us with
  //  the first match in the inLine, and fills in file1 and file2.
  //
  char    inLine[1024];
  readHeader(inLine, inFile, _file1, _file2, 0L);

  fprintf(stderr, "Opening '%s' for sequence one.\n", _file1);
  fprintf(stderr, "Opening '%s' for sequence two.\n", _file2);

  //  Open some FastAWrappers for each of the files -- we use these
  //  only to get the length of the sequence.
  //
  _seq1 = new FastAWrapper(_file1);
  _seq2 = new FastAWrapper(_file2);

  _seq1->openIndex();
  _seq2->openIndex();

  _matchesLen = 0;
  _matchesMax = 2 * 1048576;
  _matches    = new match_t [_matchesMax];

  //  For the coverage to work correctly, we need to either have one
  //  intervalList per input sequence, or build a table of the chained
  //  sequence positions.
  //
  u64bit  *offset1 = new u64bit [_seq1->getNumberOfSequences()];
  u64bit  *offset2 = new u64bit [_seq2->getNumberOfSequences()];

  offset1[0] = 1000000;
  for (u32bit i=1; i<_seq1->getNumberOfSequences(); i++)
    offset1[i] = offset1[i-1] + _seq1->sequenceLength(i-1) + 1;

  offset2[0] = 1000000;
  for (u32bit i=1; i<_seq2->getNumberOfSequences(); i++)
    offset2[i] = offset2[i-1] + _seq2->sequenceLength(i-1) + 1;

  u32bit        skippedCount = 0;
  u64bit        skippedLengthA = 0;
  u64bit        skippedLengthB = 0;

  intervalList  intervalA;
  intervalList  intervalB;

  while (!feof(inFile)) {
    if (inLine[0] == 'M') {
      splitToWords  S(inLine);

      if ((S[1][0] == 'u') || (S[1][0] == 'x')) {
        //if ((S[1][0] == 'r')) {
        u32bit  iid1=0, pos1=0, len1=0, ori1=0;
        u32bit  iid2=0, pos2=0, len2=0, ori2=0;
        decodeMatch(S, iid1, pos1, len1, ori1, iid2, pos2, len2, ori2);


        if ((pos1 + len1) > _seq1->sequenceLength(iid1)) {
          chomp(inLine);
          fprintf(stderr, "Too Long in 1: "u32bitFMT" %s\n", _seq1->sequenceLength(iid1), inLine);
        }

        if ((pos2 + len2) > _seq2->sequenceLength(iid2)) {
          chomp(inLine);
          fprintf(stderr, "Too Long in 2: "u32bitFMT" %s\n", _seq2->sequenceLength(iid2), inLine);
        }


        if ((iid1 >= _seq1->getNumberOfSequences()) || (iid2 >= _seq2->getNumberOfSequences())) {
          //  Hmmm.  Skip it.
          skippedCount++;
          skippedLengthA += len1;
          skippedLengthB += len2;
        } else {
          intervalA.add(offset1[iid1] + (u64bit)pos1, (u64bit)len1);
          intervalB.add(offset2[iid2] + (u64bit)pos2, (u64bit)len2);

          //  Add it to our list of matches
          //
          if (_matchesLen > _matchesMax) {
            fprintf(stderr, "SORRY!  I don't feel like reallocating matches.  Increase\n");
            fprintf(stderr, "the preallocated size in %s\n", __FILE__);
            exit(1);
          }

          _matches[_matchesLen].matchid = _matchesLen;

          _matches[_matchesLen].iid1 = iid1;
          _matches[_matchesLen].pos1 = pos1;
          _matches[_matchesLen].len1 = len1;
          _matches[_matchesLen].ori1 = ori1;

          _matches[_matchesLen].iid2 = iid2;
          _matches[_matchesLen].pos2 = pos2;
          _matches[_matchesLen].len2 = len2;
          _matches[_matchesLen].ori2 = ori2;

          _matchesLen++;
        }
      }
    }

    fgets(inLine, 1024, inFile);
  }

  fprintf(stderr, "skipped "u32bitFMT" matches with length "u64bitFMT" and "u64bitFMT"\n",
          skippedCount, skippedLengthA, skippedLengthB);

  fprintf(stderr, "intervalLength A "u64bitFMT" B "u64bitFMT"\n",
          (u64bit)intervalA.sumOfLengths(),
          (u64bit)intervalB.sumOfLengths());

  intervalA.merge();
  intervalB.merge();

  fprintf(stderr, "coveredLength  A "u64bitFMT" B "u64bitFMT"\n",
          (u64bit)intervalA.sumOfLengths(),
          (u64bit)intervalB.sumOfLengths());
}





void
printAnno(FILE *F,
          char label,
          u32bit axis,
          span_t *span,
          u32bit match1=u32bitZERO, u32bit off1=u32bitZERO, match_t *m1=0L,
          u32bit match2=u32bitZERO, u32bit off2=u32bitZERO, match_t *m2=0L) {

  fprintf(F, "%c "u32bitFMT" "u32bitFMT"-"u32bitFMT"  ",
          label,
          span->_iid, span->_beg, span->_end);

  if (match1) {
    fprintf(F, u32bitFMT","u32bitFMT" ", match1 >> COLORSHIFT, match1 & COLORMASK);

    if (m1) {
      if (axis == 1) {
        u32bit  sta = m1->pos2 + off2;
        u32bit  end = m1->pos2 + off2 + (span->_end - span->_beg);

        if (m1->ori2) {
          sta = m1->pos2 + m1->len2 - off2 - (span->_end - span->_beg);
          end = m1->pos2 + m1->len2 - off2;
        }

        fprintf(F, "("u32bitFMT"-"u32bitFMT") ", sta, end);
      } else {
        fprintf(F, "("u32bitFMT"-"u32bitFMT") ", m1->pos1 + off1, m1->pos1 + off1 + span->_end - span->_beg);
      }
    }
  }

  if (match2) {
    fprintf(F, u32bitFMT","u32bitFMT" ", match2 >> COLORSHIFT, match2 & COLORMASK);

    if (m2)
      if (axis == 1) {
        u32bit  sta = m2->pos2 + off2;
        u32bit  end = m2->pos2 + off2 + (span->_end - span->_beg);

        if (m2->ori2) {
          sta = m2->pos2 + m2->len2 - off2 - (span->_end - span->_beg);
          end = m2->pos2 + m2->len2 - off2;
        }

        fprintf(F, "("u32bitFMT"-"u32bitFMT") ", sta, end);
      } else {
        fprintf(F, "("u32bitFMT"-"u32bitFMT") ", m2->pos1 + off1, m2->pos1 + off1 + span->_end - span->_beg);
      }
  }

  fprintf(F, "\n");
}


#if 0
          fprintf(stderr, "span: "u32bitFMT"-"u32bitFMT" ("u32bitFMT")  off1="u32bitFMT" off2="u32bitFMT"\n",
                  span->_beg, span->_end, spanLen, off1, off2);
          fprintf(stderr, "DIFF "u32bitFMT" "u32bitFMT" "u32bitFMT" ("u32bitFMT") -(%d)- " u32bitFMT" "u32bitFMT" "u32bitFMT" ("u32bitFMT"-"u32bitFMT")\n",
                  m1->iid1, m1->pos1, m1->len1,
                  m1->pos1+off1,
                  m1->ori2,
                  m1->iid2, m1->pos2, m1->len2,
                  pos1l, pos1r);
          fprintf(stderr, "     "u32bitFMT" "u32bitFMT" "u32bitFMT" ("u32bitFMT") -(%d)- "u32bitFMT" "u32bitFMT" "u32bitFMT" ("u32bitFMT"-"u32bitFMT")\n",
                  m2->iid1, m2->pos1, m2->len1,
                  m2->pos1+off2,
                  M1->_matches[match2].ori2,
                  m2->iid2, m2->pos2, m2->len2,
                  pos2l, pos2r);
#endif








int
main(int argc, char **argv) {

  if (argc != 3) {
    fprintf(stderr, "usage: %s <matches-1> <matches-2>\n", argv[0]);
    exit(1);
  }

  matchList  *M1 = new matchList(argv[1]);
  matchList  *M2 = new matchList(argv[2]);

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
  for (u32bit i=0; i<M1->_seq1->getNumberOfSequences(); i++) {
    S1->addNewSpan(i, M1->_seq1->sequenceLength(i));
  }
  for (u32bit i=0; i<M1->_seq2->getNumberOfSequences(); i++) {
    S2->addNewSpan(i, M1->_seq2->sequenceLength(i));
  }

  //  Add every match to the spanTrees.

  for (u32bit i=0; i<M1->_matchesLen; i++) {
    S1->addMatch(M1->_matches+i, 0, 0);
    S2->addMatch(M1->_matches+i, 1, 0);
  }
  for (u32bit i=0; i<M2->_matchesLen; i++) {
    S1->addMatch(M2->_matches+i, 0, 1);
    S2->addMatch(M2->_matches+i, 1, 1);
  }

  //fprintf(stderr, "S1:"u32bitFMT" S2:"u32bitFMT"\n", S1->size(), S2->size());

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
  u32bit   same[2]         = {0,0};
  u32bit   inconsistent[2] = {0,0};

  //  Histograms.  
  //
  u32bit   histogramMax = 1048576;
  u32bit  *unmappedH[2];
  u32bit  *uniqueH[2][2];
  u32bit  *differentH[2];
  u32bit  *sameH[2];
  u32bit  *inconsistentH[2];

  for (u32bit i=0; i<2; i++) {
    unmappedH[i]     = new u32bit [histogramMax];
    uniqueH[i][0]    = new u32bit [histogramMax];
    uniqueH[i][1]    = new u32bit [histogramMax];
    differentH[i]    = new u32bit [histogramMax];
    sameH[i]         = new u32bit [histogramMax];
    inconsistentH[i] = new u32bit [histogramMax];

    bzero(unmappedH[i], sizeof(u32bit) * histogramMax);
    bzero(uniqueH[i][0], sizeof(u32bit) * histogramMax);
    bzero(uniqueH[i][1], sizeof(u32bit) * histogramMax);
    bzero(differentH[i], sizeof(u32bit) * histogramMax);
    bzero(sameH[i], sizeof(u32bit) * histogramMax);
    bzero(inconsistentH[i], sizeof(u32bit) * histogramMax);
  }

#define updateHistogram(H,x)  H[ ((x) > histogramMax) ? 0 : x ]++

  //  Stats, debugging mostly
  u32bit   same1sf=0, same1ef=0, same1sn=0, same1en=0;
  u32bit   same2sf=0, same2ef=0, same2sn=0, same2en=0;
  u32bit   same1f=0, same1n=0;
  u32bit   same2f=0, same2n=0;

  //  We write the annotation to here
  FILE *Aanno = fopen("overlap.Aannotation", "w");
  FILE *Banno = fopen("overlap.Bannotation", "w");


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
  //

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

      printAnno(Aanno, '1', 1, span, span->_matches[0]);
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

          printAnno(Aanno, 'Y', 1, span, match1, off1, m1, match2, off2, m2);
        } else {
          different[0] += spanLen;
          updateHistogram(differentH[0], spanLen);

          printAnno(Aanno, 'N', 1, span, match1, off1, m1, match2, off2, m2);
        }
      }
    } else {
      inconsistent[0] += spanLen;
      updateHistogram(inconsistentH[0], spanLen);

      printAnno(Aanno, '?', 1, span);
    }

    node1 = dict_next(S1->_tree, node1);
  }

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

      printAnno(Banno, '1', 2, span, span->_matches[0]);
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

          printAnno(Banno, 'Y', 2, span, match1, off1, m1, match2, off2, m2);
        } else {
          different[1] += spanLen;
          updateHistogram(differentH[1], spanLen);

          printAnno(Banno, 'N', 2, span, match1, off1, m1, match2, off2, m2);
        }
      }
    } else {
      inconsistent[1] += spanLen;
      updateHistogram(inconsistentH[1], spanLen);

      printAnno(Banno, '?', 2, span);
    }

    node2 = dict_next(S2->_tree, node2);
  }

  fclose(Aanno);
  fclose(Banno);

  fprintf(stderr, "unmapped:           A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", unmapped[0], unmapped[1]);
  fprintf(stderr, "unique mapping 1:   A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", unique[0][0], unique[1][0]);
  fprintf(stderr, "unique mapping 2:   A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", unique[0][1], unique[1][1]);
  fprintf(stderr, "different:          A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", different[0], different[1]);
  fprintf(stderr, "same:               A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", same[0], same[1]);
  fprintf(stderr, "inconsistent:       A:"u32bitFMTW(10)" B:"u32bitFMTW(10)"\n", inconsistent[0], inconsistent[1]);

  fprintf(stderr, "\n");
  fprintf(stderr, "same1 both:         normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same1n, same1f);
  fprintf(stderr, "same1 start:        normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same1sn, same1sf);
  fprintf(stderr, "same1 end:          normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same1en, same1ef);
  fprintf(stderr, "\n");
  fprintf(stderr, "same2 both:         normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same2n, same2f);
  fprintf(stderr, "same2 start:        normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same2sn, same2sf);
  fprintf(stderr, "same2 end:          normal:"u32bitFMTW(10)" flipped:"u32bitFMTW(10)"\n", same2en, same2ef);


  FILE *hout = fopen("overlap.histogram", "w");
  for (u32bit i=0; i<histogramMax; i++) {
    fprintf(hout, u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)" "u32bitFMTW(5)"\n",
            i,
            unmappedH[0][i], uniqueH[0][0][i], uniqueH[0][1][i], differentH[0][i], sameH[0][i], inconsistentH[0][i],
            unmappedH[1][i], uniqueH[1][0][i], uniqueH[1][1][i], differentH[1][i], sameH[1][i], inconsistentH[1][i]);
  }
  fclose(hout);

#if 0

plot [0:500] [0:1500] \
"overlap.histogram" using  2 title "unmapped A" with lines, \
"overlap.histogram" using  3 title "unique A, mapping 1" with lines, \
"overlap.histogram" using  4 title "unique A, mapping 2" with lines, \
"overlap.histogram" using  5 title "different A" with lines, \
"overlap.histogram" using  6 title "same A" with lines, \
"overlap.histogram" using  7 title "inconsistent A" with lines, \
"overlap.histogram" using  8 title "unmapped B" with lines, \
"overlap.histogram" using  9 title "unique B, mapping 1" with lines, \
"overlap.histogram" using 10 title "unique B, mapping 2" with lines, \
"overlap.histogram" using 11 title "different B" with lines, \
"overlap.histogram" using 12 title "same B" with lines, \
"overlap.histogram" using 13 title "inconsistent B" with lines

plot [0:500] [0:1500] \
"overlap.histogram" using  3 title "unique A, mapping 1" with lines, \
"overlap.histogram" using  4 title "unique A, mapping 2" with lines, \
"overlap.histogram" using  9 title "unique B, mapping 1" with lines, \
"overlap.histogram" using 10 title "unique B, mapping 2" with lines

#endif


  //  Deleting the spanTrees takes a long time, so we don't bother with any cleanup.
  return(0);
}

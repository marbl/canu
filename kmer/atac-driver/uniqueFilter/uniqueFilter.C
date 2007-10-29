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

#include "util++.H"
#include "atac.H"

//  Kaz Kylheku <kaz@ashi.footprints.net> library.
#include "kazlib/dict.h"
#include "kazlib/except.h"
#include "kazlib/hash.h"
#include "kazlib/list.h"
#include "kazlib/sfx.h"

//  Filters out matches that have non-unique pieces.  Does not discard
//  the whole match, but just trims out the non-unique section.
//
//  Original implementation in Python by Clark Mobarry:
//
//    sort the matches in X
//    apply the mask to the X
//    sort the matchs in Y
//    apply the mask to the Y
//    output the matches
//
//  if we keep the coverage intervals in core, we get around sorting
//  the matches.  how big can they be -- especially if we only keep
//  things with > 1 coverage!  But, we also can't use an elegant
//  algorithm for trimming/splitting.


//  We can abuse this to subtract matches from other matches.  The
//  operation done for removing non-unique is to first find any
//  overlapping intervals in the set, then subract those from the
//  input.  If we instead just find all intervals in a set of matches,
//  and subtract those from the input, we get subtraction.


//  Reads the input, builds an interval list of the regions
//  that have coverage > 1.
//
struct coverage1_s {
  u32bit  axis;
  u32bit  position;
  int     increment;
};

struct coverage2_s {
  u32bit  axis;
  u32bit  beg;
  u32bit  end;
  u32bit  coverage;
};

struct match_s {
  u32bit  iid1, pos1, len1, ori1;
  u32bit  iid2, pos2, len2, ori2;
};

int
sortCoverage1(const void *a, const void *b) {
  const coverage1_s *A = *((const coverage1_s * const *)a);
  const coverage1_s *B = *((const coverage1_s * const *)b);

  if (A->axis      < B->axis)       return(-1);
  if (A->axis      > B->axis)       return(1);
  if (A->position  < B->position)   return(-1);
  if (A->position  > B->position)   return(1);
  if (A->increment > B->increment)  return(-1);
  if (A->increment < B->increment)  return(1);
  return(0);
}


//  Not a complete comparison, but we only use this for an interval
//  list.
//
int
sortCoverage2(const void *a, const void *b) {
  const coverage2_s *A = *((const coverage2_s * const *)a);
  const coverage2_s *B = *((const coverage2_s * const *)b);

  if (A->axis     < B->axis)  return(-1);
  if (A->axis     > B->axis)  return(1);
  if (A->beg      < B->beg)   return(-1);
  if (A->beg      > B->beg)   return(1);
  return(0);
}


//  Same as sortCoverage2, but the array being sorted is not
//  an array of pointers.
//
int
sortCoverage3(const void *a, const void *b) {
  const coverage2_s *A = (const coverage2_s *)a;
  const coverage2_s *B = (const coverage2_s *)b;

  if (A->axis     < B->axis)  return(-1);
  if (A->axis     > B->axis)  return(1);
  if (A->beg      < B->beg)   return(-1);
  if (A->beg      > B->beg)   return(1);
  return(0);
}




//  An interval list, searchable
//
class coverageIntervals {
public:
  dict_t       *_il;
  dict_load_t   _load;
public:
  coverageIntervals() {
    _il = dict_create(DICTCOUNT_T_MAX, sortCoverage2);
  };

  ~coverageIntervals() {
    dict_free(_il);
    pfree();
  };


  //  We want to return the first node that is before our thing.  Our comparison
  //  tests the start position.  If there is no first node before our thing,
  //  return the next node.
  //
  dnode_t  *lookup(void *thing) {
    dnode_t  *it = dict_upper_bound(_il, thing);
    if (it == 0L)
      it = dict_lower_bound(_il, thing);
    return(it);
  };


  void  addInterval(int axis, int beg, int end, int coverage) {
    dnode_t      *node = (dnode_t *)palloc(sizeof(dnode_t));
    coverage2_s  *cov  = (coverage2_s *)palloc(sizeof(coverage2_s));

    cov->axis     = axis;
    cov->beg      = beg;
    cov->end      = end;
    cov->coverage = coverage;

    //  initialize the node with the value
    dnode_init(node, 0L);

    //  insert the node into the tree using the key
    dict_insert(_il, node, (void *)cov);
  };

  void  beginLoad(void) {
    dict_load_begin(&_load, _il);
  };
  void  endLoad(void) {
    dict_load_end(&_load);
  };
  void  loadInterval(int axis, int beg, int end, int coverage) {
    dnode_t      *node = (dnode_t *)palloc(sizeof(dnode_t));
    coverage2_s  *cov  = (coverage2_s *)palloc(sizeof(coverage2_s));

    cov->axis     = axis;
    cov->beg      = beg;
    cov->end      = end;
    cov->coverage = coverage;

    //  initialize the node with the value
    dnode_init(node, 0L);

    //  insert the node into the tree using the key
    dict_load_next(&_load, node, (void *)cov);
  };
};




void
offsetsToCoverage(u32bit minCov, bigQueue *I, coverageIntervals *L) {
  u32bit  axis     = ~u32bitZERO;
  u32bit  position = ~u32bitZERO;
  u32bit  coverage = 0;
  u64bit  covered  = 0;


  L->beginLoad();
  speedCounter D(" %8.0f matches treed  (%8.2f matches/sec)\r", 1, 511, false);
  while (I->next()) {
    coverage1_s  *cov1 = (coverage1_s *)I->get();

    if ((cov1->axis != axis) && (coverage != 0))
      fprintf(stderr, "Sorting error -- have coverage at the end of an axis.\n"), exit(1);

    int length = cov1->position - position;

    if ((coverage >= minCov) && (length > 0)) {
      D.tick();
      L->loadInterval(axis, position, position+length, coverage);
      covered += length;
    }

    //  Occasionally, we get stung by insisting to use unsigned
    //  numbers.  This is one of them.
    //
    if ((coverage == 0) && (cov1->increment == -1))
      fprintf(stderr, "Sorting error -- have negative coverage (axis="u32bitFMT" position="u32bitFMT")!\n",
              axis, position), exit(1);

    coverage   += cov1->increment;
    axis        = cov1->axis;
    position    = cov1->position;
  }
  D.finish();
  L->endLoad();

  fprintf(stderr, "offsetsToCoverage()-- Found "u64bitFMT" bases at coverage "u32bitFMT" or greater.\n",
          covered, minCov);
}





void
findCoverageIntervals(char const *fileName,
                      u32bit      minCov,
                      coverageIntervals *Fint,
                      coverageIntervals *Rint) {
  bigQueue  F(sortCoverage1, 0L, 0L, 0L, sizeof(coverage1_s), 128, 0L);
  bigQueue  R(sortCoverage1, 0L, 0L, 0L, sizeof(coverage1_s), 128, 0L);

  //
  //  Read the input file, building a bigQueue of the interval offsets
  //


  atacFileStream  AF(fileName);
  atacMatch      *m = AF.nextMatch('u');

  while (m) {
    coverage1_s *fbeg = (coverage1_s *)malloc(sizeof(coverage1_s));
    coverage1_s *fend = (coverage1_s *)malloc(sizeof(coverage1_s));
    coverage1_s *rbeg = (coverage1_s *)malloc(sizeof(coverage1_s));
    coverage1_s *rend = (coverage1_s *)malloc(sizeof(coverage1_s));

    fbeg->axis      = m->iid1;
    fbeg->position  = m->pos1;
    fbeg->increment = 1;

    fend->axis      = m->iid1;
    fend->position  = m->pos1 + m->len1;
    fend->increment = -1;

    rbeg->axis      = m->iid2;
    rbeg->position  = m->pos2;
    rbeg->increment = 1;

    rend->axis      = m->iid2;
    rend->position  = m->pos2 + m->len2;
    rend->increment = -1;

    F.add(fbeg);
    F.add(fend);
    R.add(rbeg);
    R.add(rend);

    m = AF.nextMatch('u');
  }

  //  Sort each bigQueue
  //
  F.sort();
  R.sort();

  //  Convert the interval offsets into a coverage interval list
  //
  offsetsToCoverage(minCov, &F, Fint);
  offsetsToCoverage(minCov, &R, Rint);
}










void
intersectTest(match_s            *matches,
              u32bit              matchesLen,
              coverageIntervals  *Fint,
              coverageIntervals  *Rint,
              u32bit              matchNumber) {
  bool errors = false;

  for (u32bit i=0; i<matchesLen; i++) {
    coverage2_s   thing;

    //  Query the tree for the first interval intersecting iid1
    //
    thing.axis     = matches[i].iid1;
    thing.beg      = matches[i].pos1;
    thing.end      = matches[i].pos1 + matches[i].len1;
    thing.coverage = 0;
    dnode_t *node1 = Fint->lookup(&thing);

    thing.axis     = matches[i].iid2;
    thing.beg      = matches[i].pos2;
    thing.end      = matches[i].pos2 + matches[i].len2;
    thing.coverage = 0;
    dnode_t *node2 = Rint->lookup(&thing);


    //  Keep iterating until the node returned from the tree
    //  is empty, or it is after our region
    //
    while (node1 && node2) {
      const coverage2_s *key1 = 0L, *key2 = 0L;

      bool  isect1=false, before1=false;
      bool  isect2=false, before2=false;

      if (node1) {
        key1   = (const coverage2_s *)dnode_getkey(node1);

        isect1 = ((key1->axis == matches[i].iid1) &&
                  (matches[i].pos1 < key1->end) &&
                  (key1->beg       < matches[i].pos1 + matches[i].len1));

        before1 = ((key1->axis  < matches[i].iid1) ||
                   ((key1->axis == matches[i].iid1) && (key1->beg < matches[i].pos1 + matches[i].len1)));
      }


      if (node2) {
        key2   = (const coverage2_s *)dnode_getkey(node2);

        isect2 = ((key2->axis == matches[i].iid2) &&
                  (matches[i].pos2 < key2->end) &&
                  (key2->beg       <  matches[i].pos2 + matches[i].len2));

        before2 = ((key2->axis  < matches[i].iid2) ||
                   ((key2->axis == matches[i].iid2) && (key2->beg < matches[i].pos2 + matches[i].len2)));
      }


      if (isect1) {
        fprintf(stderr, "Got fwd intersection on i="u32bitFMT" matchNumber="u32bitFMT"\n", i, matchNumber);
        fprintf(stdout, "--"u32bitFMT" "u32bitFMT" 1    "u32bitFMT" "u32bitFMT" %d\n",
                matches[i].pos1, matches[i].pos1 + matches[i].len1,
                matches[i].pos2, matches[i].pos2 + matches[i].len2, matches[i].ori2 ? 1 : -1);
        fprintf(stdout, "--key1 beg="u32bitFMT" end="u32bitFMT"\n",
                key1->beg, key1->end);
        errors = true;
      }
      if (isect2) {
        fprintf(stderr, "Got rev intersection on i="u32bitFMT" matchNumber="u32bitFMT"\n", i, matchNumber);
        fprintf(stdout, "--"u32bitFMT" "u32bitFMT" 1    "u32bitFMT" "u32bitFMT" %d\n",
                matches[i].pos1, matches[i].pos1 + matches[i].len1,
                matches[i].pos2, matches[i].pos2 + matches[i].len2, matches[i].ori2 ? 1 : -1);
        fprintf(stdout, "--key2 beg="u32bitFMT" end="u32bitFMT"\n",
                key2->beg, key2->end);
        errors = true;
      }


      //  If we intersected or were before, move to the next, otherwise,
      //  stop
      //
      if (isect1 || before1)
        node1 = dict_next(Fint->_il, node1);
      else
        node1 = 0L;

      if (isect2 || before2)
        node2 = dict_next(Rint->_il, node2);
      else
        node2 = 0L;
    }
  }

  if (errors)
    abort();
}









//  This is used all over the place.
//
#define D08D2 u32bitFMTW(8)" "u32bitFMTW(8)
#define KEY1THING "key1 = "u32bitFMTW(8)" "u32bitFMTW(8)" "u32bitFMTW(8)"    thing = "D08D2" "D08D2"\n"
#define KEY2THING "key2 = "u32bitFMTW(8)" "u32bitFMTW(8)" "u32bitFMTW(8)"    thing = "D08D2" "D08D2"\n"

int
main(int argc, char **argv) {
  char   *inputName = 0L;
  char   *subtractName = 0L;

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-i") == 0) {
      inputName = argv[++arg];
    } else if (strcmp(argv[arg], "-s") == 0) {
      subtractName = argv[++arg];
    } else {
      fprintf(stderr, "usage: %s [-h] [-s subtractFile] [-i inputFile]\n", argv[0]);
      fprintf(stderr, "  -s     instead of finding regions to remove by looking\n");
      fprintf(stderr, "         for duplicatd regions in inputFile, load them\n");
      fprintf(stderr, "         from subtractFile.\n");
      exit(1);
    }
    arg++;
  }

  if (inputName == 0L)
    fprintf(stderr, "usage: %s [-i inputfile] [-o outputfile] [-h]\n", argv[0]), exit(1);

  coverageIntervals   *Fint = new coverageIntervals;
  coverageIntervals   *Rint = new coverageIntervals;

  if (subtractName)
    findCoverageIntervals(subtractName, 1, Fint, Rint);
  else
    findCoverageIntervals(inputName, 2, Fint, Rint);



  //  The original implementation would then sort the matches in X,
  //  merge the sorted intervals and the sorted matches together, resort
  //  the matches by Y, and merge with Rcov.  That is a lot of work to
  //  avoid keeping two interval lists in memory.
  //
  //  We build an in-core interval list for both assemblies, and
  //  then stream the matches by it.

  //  we need to ask the interval list:
  //    return the intervals that are covered by this interval


  u32bit        matchesLen = 0;
  u32bit        matchesMax = 1024;
  match_s      *matches    = new match_s [matchesMax];
  match_s       extent;

  u32bit        matchNumber = 0;

  atacFileStream  AF(inputName);
  atacMatch      *m = AF.nextMatch('u');

  while (m) {
    matches[0].iid1 = m->iid1;
    matches[0].pos1 = m->pos1;
    matches[0].len1 = m->len1;
    matches[0].ori1 = m->fwd1;
    matches[0].iid2 = m->iid2;
    matches[0].pos2 = m->pos2;
    matches[0].len2 = m->len2;
    matches[0].ori2 = m->fwd2;

    //  Save the original, we'll use this to test for intersections.
    //
    memcpy(&extent, matches, sizeof(match_s));

    bool fwd = (matches[0].ori1 == matches[0].ori2);


    //  Query the tree for the first interval intersecting iid1
    //    XXX: Should this be upper_bound instead?
    //
    //  A scratch interval used for querying the list
    //
    coverage2_s   thing;

    thing.axis          = extent.iid1;
    thing.beg           = extent.pos1;
    thing.end           = extent.pos1 + extent.len1;
    thing.coverage      = 0;
    dnode_t *node1start = Fint->lookup(&thing);
    dnode_t *node1      = node1start;


    thing.axis          = extent.iid2;
    thing.beg           = extent.pos2;
    thing.end           = extent.pos2 + extent.len2;
    thing.coverage      = 0;
    dnode_t *node2start = Rint->lookup(&thing);
    dnode_t *node2      = node2start;

    //  while the node intersects the match, trim or split it, then
    //  get the next node.
    //
    //  any way I tried this, it's ugly.
    //
    //  if there is one match, then it is in [0].  If there is more
    //  than one match, then the trimmed match is in [0], but the
    //  split matches are in [1] on.


    //  XXX  the problem is that we split off things, then move
    //  to the next match, without checking previously split
    //  things against this match
    //

    //  Keep iterating until the node returned from the tree
    //  is empty, or it is after our region
    //
    while (node1 || node2) {

      bool  before1=false;
      bool  before2=false;

      bool  modified=false;

      const coverage2_s *key1 = 0L;
      const coverage2_s *key2 = 0L;

      if (node1) {
        key1    = (const coverage2_s *)dnode_getkey(node1);
        before1 = ((key1->axis  < extent.iid1) ||
                   ((key1->axis == extent.iid1) && (key1->beg < extent.pos1 + extent.len1)));


        //  Three cases: (1) we trim off the front, (2) trim off the
        //  back or (3) split.  And, (4) delete the whole damn thing.
        //
        //  Further complicated by having multiple things to try.
        //
        //  If anything is modified, reset the node to the start
        //

        for (u32bit i=0; i<matchesLen; i++) {
          if (matches[i].len1 == 0)
            continue;

          if ((key1->beg <= matches[i].pos1) &&
              ((matches[i].pos1 + matches[i].len1) <= key1->end)) {
            modified = true;

            //  Trim the whole thing?
            //
            matches[i].pos1 = 0;
            matches[i].len1 = 0;
            matches[i].pos2 = 0;
            matches[i].len2 = 0;
          } else if ((matches[i].pos1 < key1->beg) &&
                     (key1->end < (matches[i].pos1 + matches[i].len1))) {
            modified = true;

            //  Contained.  Split it.
            //

            //  The left half
            //
            int newLen = key1->beg - matches[i].pos1;

            matches[matchesLen].iid1 = matches[i].iid1;
            matches[matchesLen].pos1 = matches[i].pos1;
            matches[matchesLen].len1 = newLen;
            matches[matchesLen].ori1 = matches[i].ori1;

            if (fwd) {
              matches[matchesLen].iid2 = matches[i].iid2;
              matches[matchesLen].pos2 = matches[i].pos2;
              matches[matchesLen].len2 = newLen;
              matches[matchesLen].ori2 = matches[i].ori2;
            } else {
              matches[matchesLen].iid2 = matches[i].iid2;
              matches[matchesLen].pos2 = matches[i].pos2 + matches[i].len2 - newLen;
              matches[matchesLen].len2 = newLen;
              matches[matchesLen].ori2 = matches[i].ori2;
            }

            matchesLen++;

            //  The right half
            //
            newLen = matches[i].pos1 + matches[i].len1 - key1->end;

            matches[matchesLen].iid1 = matches[i].iid1;
            matches[matchesLen].pos1 = key1->end;
            matches[matchesLen].len1 = newLen;
            matches[matchesLen].ori1 = matches[i].ori1;

            if (fwd) {
              matches[matchesLen].iid2 = matches[i].iid2;
              matches[matchesLen].pos2 = matches[i].pos2 + (key1->end - matches[i].pos1);
              matches[matchesLen].len2 = newLen;
              matches[matchesLen].ori2 = matches[i].ori2;
            } else {
              matches[matchesLen].iid2 = matches[i].iid2;
              matches[matchesLen].pos2 = matches[i].pos2;
              matches[matchesLen].len2 = newLen;
              matches[matchesLen].ori2 = matches[i].ori2;
            }

            matchesLen++;

            //  Invalidate this match
            //
            matches[i].pos1 = 0;
            matches[i].len1 = 0;
            matches[i].pos2 = 0;
            matches[i].len2 = 0;
          } else if ((key1->beg <= matches[i].pos1) &&
                     (matches[i].pos1 < key1->end)) {
            modified = true;

            //  Trim the begin?
            //

            int trimLen = key1->end - matches[i].pos1;
            matches[i].pos1 += trimLen;
            matches[i].len1 -= trimLen;

            if (fwd == true)
              matches[i].pos2 += trimLen;
            matches[i].len2 -= trimLen;

          } else if ((key1->beg < (matches[i].pos1 + matches[i].len1)) &&
                     ((matches[i].pos1 + matches[i].len1) <= key1->end)) {
            modified = true;

            //  Trim the end?
            //

            int trimLen = matches[i].pos1 + matches[i].len1 - key1->beg;
            matches[i].len1 -= trimLen;

            if (fwd == false)
              matches[i].pos2 += trimLen;
            matches[i].len2 -= trimLen;
          }
        } 
      }  // isect



      if (node2) {
        key2    = (const coverage2_s *)dnode_getkey(node2);
        before2 = ((key2->axis  < extent.iid2) ||
                   ((key2->axis == extent.iid2) && (key2->beg < extent.pos2 + extent.len2)));

        for (u32bit i=0; i<matchesLen; i++) {
          if (matches[i].len2 == 0)
            continue;

          if ((key2->beg <= matches[i].pos2) &&
              ((matches[i].pos2 + matches[i].len2) <= key2->end)) {
            modified = true;

            //  Trim the whole thing?
            //
            matches[i].pos1 = 0;
            matches[i].len1 = 0;
            matches[i].pos2 = 0;
            matches[i].len2 = 0;
          } else if ((matches[i].pos2 < key2->beg) &&
                     (key2->end < (matches[i].pos2 + matches[i].len2))) {
            modified = true;

            //  Contained.  Split it.
            //

            //  The left (forward strand) half
            //

            if (fwd) {
              int newLen = key2->beg - matches[i].pos2;

              matches[matchesLen].iid1 = matches[i].iid1;
              matches[matchesLen].pos1 = matches[i].pos1;
              matches[matchesLen].len1 = newLen;
              matches[matchesLen].ori1 = matches[i].ori1;

              matches[matchesLen].iid2 = matches[i].iid2;
              matches[matchesLen].pos2 = matches[i].pos2;
              matches[matchesLen].len2 = newLen;
              matches[matchesLen].ori2 = matches[i].ori2;
            } else {
              int newLen = matches[i].pos2 + matches[i].len2 - key2->end;

              matches[matchesLen].iid1 = matches[i].iid1;
              matches[matchesLen].pos1 = matches[i].pos1;
              matches[matchesLen].len1 = newLen;
              matches[matchesLen].ori1 = matches[i].ori1;

              matches[matchesLen].iid2 = matches[i].iid2;
              matches[matchesLen].pos2 = key2->end;
              matches[matchesLen].len2 = newLen;
              matches[matchesLen].ori2 = matches[i].ori2;
            }

            matchesLen++;

            //  The right (forward strand) half
            //
            if (fwd) {
              int newLen = matches[i].pos2 + matches[i].len2 - key2->end;

              matches[matchesLen].iid1 = matches[i].iid1;
              matches[matchesLen].pos1 = matches[i].pos1 + key2->end - matches[i].pos2;
              matches[matchesLen].len1 = newLen;
              matches[matchesLen].ori1 = matches[i].ori1;

              matches[matchesLen].iid2 = matches[i].iid2;
              matches[matchesLen].pos2 = key2->end;
              matches[matchesLen].len2 = newLen;
              matches[matchesLen].ori2 = matches[i].ori2;
            } else {
              int newLen = key2->beg - matches[i].pos2;

              matches[matchesLen].iid1 = matches[i].iid1;
              matches[matchesLen].pos1 = matches[i].pos1 + matches[i].pos2 + matches[i].len2 - key2->beg;
              matches[matchesLen].len1 = newLen;
              matches[matchesLen].ori1 = matches[i].ori1;

              matches[matchesLen].iid2 = matches[i].iid2;
              matches[matchesLen].pos2 = matches[i].pos2;
              matches[matchesLen].len2 = newLen;
              matches[matchesLen].ori2 = matches[i].ori2;
            }

            matchesLen++;

            //  Invalidate this match
            //
            matches[i].pos1 = 0;
            matches[i].len1 = 0;
            matches[i].pos2 = 0;
            matches[i].len2 = 0;
          } else if ((key2->beg <= matches[i].pos2) &&
                     (matches[i].pos2 < key2->end)) {
            modified = true;

            //  Trim the begin?  fwdOK, revOK
            //

            int trimLen = key2->end - matches[i].pos2;
            matches[i].pos2 += trimLen;
            matches[i].len2 -= trimLen;

            if (fwd == true)
              matches[i].pos1 += trimLen;
            matches[i].len1 -= trimLen;

          } else if ((key2->beg < (matches[i].pos2 + matches[i].len2)) &&
                     ((matches[i].pos2 + matches[i].len2) <= key2->end)) {
            modified = true;

            //  Trim the end?
            //

            int trimLen = matches[i].pos2 + matches[i].len2 - key2->beg;
            matches[i].len1 -= trimLen;

            if (fwd == false)
              matches[i].pos1 += trimLen;
            matches[i].len2 -= trimLen;
          }
        } 
      }

      //  If we intersected or were before, move to the next, otherwise,
      //  stop.
      //
      if (modified)
        node1 = node1start;
      else if (before1)
        node1 = dict_next(Fint->_il, node1);
      else
        node1 = 0L;

      if (modified)
        node2 = node2start;
      else if (before2)
        node2 = dict_next(Rint->_il, node2);
      else
        node2 = 0L;
    }  //  end of while (node1 || node2)


    //  Nobody should be outside the extent
    //
    for (u32bit i=0; i<matchesLen; i++) {
      if ((matches[i].len1 > 0) && (matches[i].len2 > 0)) {
        if ((matches[i].pos1 < extent.pos1) ||
            (matches[i].pos1 + matches[i].len1 > extent.pos1 + extent.len1) ||
            (matches[i].pos2 < extent.pos2) ||
            (matches[i].pos2 + matches[i].len2 > extent.pos2 + extent.len2)) {
          fprintf(stderr, "match "u32bitFMT" is outside the extent!\n", i);
          abort();
        }
      }
    }

    //  Print out all the modified matches
    //
    for (u32bit i=0; i<matchesLen; i++) {
      if ((matches[i].len1 > 0) && (matches[i].len2 > 0)) {
        fprintf(stdout, "M %s %s."u32bitFMT" . %s "u32bitFMT" "u32bitFMT" 1 %s "u32bitFMT" "u32bitFMT" %d\n",
                m->matchuid, m->parentuid, i,
                AF.labelA(), matches[i].pos1, matches[i].len1,
                AF.labelB(), matches[i].pos2, matches[i].len2, matches[i].ori2 ? 1 : -1);
      }
    }

    //  Check that the modified matches do not intersect anything in
    //  the tree.
    //
    intersectTest(matches, matchesLen, Fint, Rint, matchNumber);

    matchNumber++;
  }
}

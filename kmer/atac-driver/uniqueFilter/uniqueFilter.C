#include <stdio.h>
#include <stdlib.h>

#include "util++.H"
#include "atac-common.H"

//  Filters out matches that have non-unique pieces.  Does not discard
//  the whole match, but just trims out the non-unique section.
//
//  Original implementation in Python by Clark Mobarry.

char   *_inputName = 0L;
char   *_outputName = 0L;

//  Reads the input, builds an interval list of the regions
//  that have coverage > 1.
//
struct coverage1_s {
  int  axis;
  int  position;
  int  increment;
};

struct coverage2_s {
  int  axis;
  int  beg;
  int  end;
  int  coverage;
};

struct match_s {
  u32bit  iid1, pos1, len1, ori1;
  u32bit  iid2, pos2, len2, ori2;
};

int
sortCoverage1(const void *a, const void *b) {
  coverage1_s *A = *((coverage1_s **)a);
  coverage1_s *B = *((coverage1_s **)b);

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
  coverage2_s *A = *((coverage2_s **)a);
  coverage2_s *B = *((coverage2_s **)b);

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
  dict_t    *_il;

public:
  coverageIntervals() {
    _il = dict_create(DICTCOUNT_T_MAX, sortCoverage2);
  };

  ~coverageIntervals() {
    //  XXX: probably also need to release all the nodes
    dict_free(_il);
  };

  void        addInterval(int axis, int beg, int end, int coverage) {
    dnode_t      *node = (dnode_t *)malloc(sizeof(dnode_t));
    coverage2_s  *cov  = (coverage2_s *)malloc(sizeof(coverage2_s));

    cov->axis     = axis;
    cov->beg      = beg;
    cov->end      = end;
    cov->coverage = coverage;

    //  initialize the node with the value
    dnode_init(node, 0L);

    //  insert the node into the tree using the key
    dict_insert(_il, node, (void *)cov);
  };
};




void
offsetsToCoverage(bigQueue *I, bigQueue *O, coverageIntervals *L) {
  int  axis     = -1;
  int  position = -1;
  int  coverage = 0;

  fprintf(stderr, "offsetsToCoverage()-- begin\n");

  while (I->next()) {
    coverage1_s  *cov1 = (coverage1_s *)I->get();

    if ((cov1->axis != axis) && (coverage != 0))
      fprintf(stderr, "Sorting error -- have coverage at the end of an axis.\n"), exit(1);

    int length = cov1->position - position;

    //  Add this interval to the interval list.  We only really care
    //  about coverage > 1 though
    //
    if ((coverage > 0) && (length > 0)) {
      coverage2_s  *cov2 = (coverage2_s *)malloc(sizeof(coverage2_s));

      cov2->axis      = axis;
      cov2->beg       = position;
      cov2->end       = position + length;
      cov2->coverage  = coverage;

      O->add(cov2);
    }

    if ((coverage > 1) && (length > 0)) {
      fprintf(stderr, "axis=%8d %8d-%8d cov=%d\n", axis, position, position+length, coverage);

      L->addInterval(axis, position, position+length, coverage);
    }

    coverage   += cov1->increment;
    axis        = cov1->axis;
    position    = cov1->position;

    if (coverage < 0)
      fprintf(stderr, "Sorting error -- have negative coverage (axis=%d position=%d)!\n",
              axis, position), exit(1);
  }

  fprintf(stderr, "offsetsToCoverage()-- end\n");
}





void
findCoverageIntervals(bigQueue *Fcov, coverageIntervals *Fint,
                      bigQueue *Rcov, coverageIntervals *Rint) {
  bigQueue  F(sortCoverage1, 0L, 0L, 0L, sizeof(coverage1_s), 128, 0L);
  bigQueue  R(sortCoverage1, 0L, 0L, 0L, sizeof(coverage1_s), 128, 0L);
  char      inLine[1024];
  FILE     *inFile = fopen(_inputName, "r");

  fprintf(stderr, "findCoverageIntervals()--\n");

  fgets(inLine, 1024, inFile);

  //
  //  Read the input file, building a bigQueue of the interval offsets
  //

  while (!feof(inFile)) {
    if (inLine[0] == 'M') {
      splitToWords  S(inLine);

      if ((S[1][0] == 'u') || (S[1][0] == 'x')) {
        coverage1_s *fbeg = (coverage1_s *)malloc(sizeof(coverage1_s));
        coverage1_s *fend = (coverage1_s *)malloc(sizeof(coverage1_s));
        coverage1_s *rbeg = (coverage1_s *)malloc(sizeof(coverage1_s));
        coverage1_s *rend = (coverage1_s *)malloc(sizeof(coverage1_s));

        u32bit  iid1=0, pos1=0, len1=0, ori1=0;
        u32bit  iid2=0, pos2=0, len2=0, ori2=0;
        decodeMatch(S, iid1, pos1, len1, ori1, iid2, pos2, len2, ori2);

        fbeg->axis      = iid1;
        fbeg->position  = pos1;
        fbeg->increment = 1;

        fend->axis      = iid1;
        fend->position  = pos1 + len1;
        fend->increment = -1;

        rbeg->axis      = iid2;
        rbeg->position  = pos2;
        rbeg->increment = 1;

        rend->axis      = iid2;
        rend->position  = pos2 + len2;
        rend->increment = -1;

        F.add(fbeg);
        F.add(fend);
        R.add(rbeg);
        R.add(rend);
      }
    }

    fgets(inLine, 1024, inFile);
  }

  fclose(inFile);

  //  Sort each bigQueue
  //
  fprintf(stderr, "findCoverageIntervals()-- begin sorting\n");

  F.sort();
  R.sort();

  //  Convert the interval offsets into a coverage interval list
  //
  offsetsToCoverage(&F, Fcov, Fint);
  offsetsToCoverage(&R, Rcov, Rint);
}


//  sort the matches in X
//  apply the mask to the X
//  sort the matchs in Y
//  apply the mask to the Y
//  output the matches
//
//  if we keep the coverage intervals in core, we get around sorting
//  the matches.  how big can they be -- especially if we only keep
//  things with > 1 coverage!
//
//  need an interval list that answers "tell me any intervals from X
//  to Y", build one interval list for each sequence in the input
//  (yikes -- 100k of them?)


int
main(int argc, char **argv) {

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-i") == 0) {
      _inputName = argv[++arg];
    } else if (strcmp(argv[arg], "-o") == 0) {
      _outputName = argv[++arg];
    } else {
      fprintf(stderr, "usage: %s [-i inputfile] [-o outputfile] [-h]\n", argv[0]);
      exit(1);
    }
    arg++;
  }

  if (_inputName == 0L)
    fprintf(stderr, "usage: %s [-i inputfile] [-o outputfile] [-h]\n", argv[0]), exit(1);

  bigQueue            *Fcov = new bigQueue(0L, 0L, 0L, 0L, sizeof(coverage2_s), 128, 0L);
  bigQueue            *Rcov = new bigQueue(0L, 0L, 0L, 0L, sizeof(coverage2_s), 128, 0L);
  coverageIntervals   *Fint = new coverageIntervals;
  coverageIntervals   *Rint = new coverageIntervals;

  //  For silly reasons, we build both the coverage list and the
  //  interval list.  We do this so we can easily get back to the
  //  original out-of-core algorithm.

  findCoverageIntervals(Fcov, Fint, Rcov, Rint);



  //  The original implementation would then sort the matches in X,
  //  merge the sorted Fcov and the sorted matches together, resort
  //  the matches by Y, and merge with Rcov.  That is a lot of work to
  //  avoid keeping two interval lists in memory.
  //
  //  We build an in-core interval list for both assemblies, and
  //  then stream the matches by it.

  //  we need to ask the interval list:
  //    return the intervals that are covered by this interval

  char      inLine[1024];
  FILE     *inFile = fopen(_inputName, "r");

  readHeader(inLine, inFile, 0L, 0L, stdout);

  //  A scratch interval used for querying the list
  //
  coverage2_s   thing;

  u32bit        matchesLen = 0;
  u32bit        matchesMax = 1024;
  match_s      *matches    = new match_s [matchesMax];
  match_s       extent;

  while (!feof(inFile)) {
    if (inLine[0] == 'M') {
      splitToWords  S(inLine);

      matchesLen = 1;
      decodeMatch(S,
                  matches[0].iid1, matches[0].pos1, matches[0].len1, matches[0].ori1,
                  matches[0].iid2, matches[0].pos2, matches[0].len2, matches[0].ori2);

      memcpy(&extent, matches, sizeof(match_s));

      bool fwd = (matches[0].ori1 == matches[0].ori2);

      chomp(inLine);
      fprintf(stderr, "%s\n", inLine);

      //  Query the tree for the first interval intersecting iid1
      //
      thing.axis     = matches[0].iid1;
      thing.beg      = matches[0].pos1;
      thing.end      = matches[0].pos1 + matches[0].len1;
      thing.coverage = 0;
      dnode_t *node1 = dict_lower_bound(Fint->_il, &thing);

      thing.axis     = matches[0].iid2;
      thing.beg      = matches[0].pos2;
      thing.end      = matches[0].pos2 + matches[0].len2;
      thing.coverage = 0;
      dnode_t *node2 = dict_lower_bound(Rint->_il, &thing);

#if 0
      if (node1) {
        const coverage2_s *key1 = (const coverage2_s *)dnode_getkey(node1);
        fprintf(stderr, "Got node1 = %08d %08d %08d\n",
                key1->axis, key1->beg, key1->end);
      }

      if (node2) {
        const coverage2_s *key1 = (const coverage2_s *)dnode_getkey(node2);
        fprintf(stderr, "Got node2 = %08d %08d %08d\n",
                key1->axis, key1->beg, key1->end);
      }
#endif

      //  while the node intersects the match, trim or split it, then
      //  get the next node.
      //
      //  any way I tried this, it's ugly.
      //
      //  if there is one match, then it is in [0].  If there is more
      //  than one match, then the trimmed match is in [0], but the
      //  split matches are in [1] on.
      //
      while (node1 || node2) {
        const coverage2_s *key1 = 0L, *key2 = 0L;

        bool  isect1 = false;
        bool  isect2 = false;

        if (node1) {
          key1   = (const coverage2_s *)dnode_getkey(node1);
          isect1 = ((key1->axis == extent.iid1) &&
                    (key1->end  >= extent.pos1) &&
                    (key1->beg  <= extent.pos1 + extent.len1));
        }

        if (isect1) {
          fprintf(stderr, "Intersect: key1 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                  key1->axis, key1->beg, key1->end,
                  extent.pos1, extent.pos1 + extent.len1,
                  extent.pos2, extent.pos2 + extent.len2);

          //  Three cases: (1) we trim off the front, (2) trim off the
          //  back or (3) split.  And, (4) delete the whole damn thing.
          //
          //  Further complicated by having multiple things to try.

          for (u32bit i=0; i<matchesLen; i++) {
            if ((key1->beg <= matches[i].pos1) &&
                ((matches[i].pos1 + matches[i].len1) <= key1->end)) {

              //  Trim the whole thing?
              //
              fprintf(stderr, "remove1    key1 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                      key1->axis, key1->beg, key1->end,
                      matches[i].pos1, matches[i].pos1 + matches[i].len1,
                      matches[i].pos2, matches[i].pos2 + matches[i].len2);

              matches[i].pos1 = 0;
              matches[i].len1 = 0;
              matches[i].pos2 = 0;
              matches[i].len2 = 0;
            } else if ((matches[i].pos1 < key1->beg) &&
                       (key1->end < (matches[i].pos1 + matches[i].len1))) {
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

              fprintf(stderr, "cont1      key1 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                      key1->axis, key1->beg, key1->end,
                      matches[matchesLen].pos1, matches[matchesLen].pos1 + matches[matchesLen].len1,
                      matches[matchesLen].pos2, matches[matchesLen].pos2 + matches[matchesLen].len2);

              matchesLen++;

              //  The right half
              //
              newLen = matches[i].pos1 + matches[i].len1 - key1->end;

              matches[matchesLen].iid1 = matches[i].iid1;
              matches[matchesLen].pos1 = matches[i].pos1 + (key1->end - matches[i].pos1);
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

              fprintf(stderr, "cont1      key1 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                      key1->axis, key1->beg, key1->end,
                      matches[matchesLen].pos1, matches[matchesLen].pos1 + matches[matchesLen].len1,
                      matches[matchesLen].pos2, matches[matchesLen].pos2 + matches[matchesLen].len2);

              matchesLen++;

              //  Invalidate this match
              //
              matches[i].pos1 = 0;
              matches[i].len1 = 0;
              matches[i].pos2 = 0;
              matches[i].len2 = 0;
            } else if ((key1->beg <= matches[i].pos1) &&
                       (matches[i].pos1 < key1->end)) {

              //  Trim the begin?
              //

              int trimLen = key1->end - matches[i].pos1;
              matches[i].pos1 += trimLen;
              matches[i].len1 -= trimLen;

              if (fwd == true)
                matches[i].pos2 += trimLen;
              matches[i].len2 -= trimLen;

              fprintf(stderr, "begin1     key1 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                      key1->axis, key1->beg, key1->end,
                      matches[i].pos1, matches[i].pos1 + matches[i].len1,
                      matches[i].pos2, matches[i].pos2 + matches[i].len2);
            } else if ((key1->beg < (matches[i].pos1 + matches[i].len1)) &&
                       ((matches[i].pos1 + matches[i].len1) <= key1->end)) {

              //  Trim the end?
              //

              int trimLen = matches[i].pos1 + matches[i].len1 - key1->beg;
              matches[i].len1 -= trimLen;

              if (fwd == false)
                matches[i].pos2 += trimLen;
              matches[i].len2 -= trimLen;

              fprintf(stderr, "end1       key1 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                      key1->axis, key1->beg, key1->end,
                      matches[i].pos1, matches[i].pos1 + matches[i].len1,
                      matches[i].pos2, matches[i].pos2 + matches[i].len2);
            } else {
              //fprintf(stderr, "no intersection\n");
              //exit(1);
            }
          } 

          node1 = dict_next(Fint->_il, node1);
        } else {
          node1 = 0L;
        }

        if (node2) {
          key2   = (const coverage2_s *)dnode_getkey(node2);
          isect2 = ((key2->axis == extent.iid2) &&
                    (key2->end  >= extent.pos2) &&
                    (key2->beg  <= extent.pos2 + extent.len2));
        }

        if (isect2) {
          fprintf(stderr, "Intersect: key2 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                  key2->axis, key2->beg, key2->end,
                  extent.pos1, extent.pos1 + extent.len1,
                  extent.pos2, extent.pos2 + extent.len2);

          for (u32bit i=0; i<matchesLen; i++) {
            if ((key2->beg <= matches[i].pos2) &&
                ((matches[i].pos2 + matches[i].len2) <= key2->end)) {

              //  Trim the whole thing?
              //
              fprintf(stderr, "remove2    key2 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                      key2->axis, key2->beg, key2->end,
                      matches[i].pos1, matches[i].pos1 + matches[i].len1,
                      matches[i].pos2, matches[i].pos2 + matches[i].len2);

              matches[i].pos1 = 0;
              matches[i].len1 = 0;
              matches[i].pos2 = 0;
              matches[i].len2 = 0;
            } else if ((matches[i].pos2 < key2->beg) &&
                (key2->end < (matches[i].pos2 + matches[i].len2))) {

              //  Contained.  Split it.
              //

              //  The left half
              //
              int newLen = key2->beg - matches[i].pos2;

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
                matches[matchesLen].pos2 = matches[i].pos2 + matches[i].len2 - (key2->end - matches[i].pos2) - newLen;
                matches[matchesLen].len2 = newLen;
                matches[matchesLen].ori2 = matches[i].ori2;
              }

              fprintf(stderr, "cont2      key2 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                      key2->axis, key2->beg, key2->end,
                      matches[matchesLen].pos1, matches[matchesLen].pos1 + matches[matchesLen].len1,
                      matches[matchesLen].pos2, matches[matchesLen].pos2 + matches[matchesLen].len2);

              matchesLen++;

              //  The right half
              //
              newLen = matches[i].pos2 + matches[i].len2 - key2->end;

              matches[matchesLen].iid1 = matches[i].iid1;
              matches[matchesLen].pos1 = matches[i].pos1 + (key2->end - matches[i].pos2);
              matches[matchesLen].len1 = newLen;
              matches[matchesLen].ori1 = matches[i].ori1;

              if (fwd) {
                matches[matchesLen].iid2 = matches[i].iid2;
                matches[matchesLen].pos2 = matches[i].pos2 + (key2->end - matches[i].pos2);
                matches[matchesLen].len2 = newLen;
                matches[matchesLen].ori2 = matches[i].ori2;
              } else {
                matches[matchesLen].iid2 = matches[i].iid2;
                matches[matchesLen].pos2 = matches[i].pos2 + matches[i].len2 - newLen;
                matches[matchesLen].len2 = newLen;
                matches[matchesLen].ori2 = matches[i].ori2;
              }

              fprintf(stderr, "cont2      key2 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                      key2->axis, key2->beg, key2->end,
                      matches[matchesLen].pos1, matches[matchesLen].pos1 + matches[matchesLen].len1,
                      matches[matchesLen].pos2, matches[matchesLen].pos2 + matches[matchesLen].len2);

              matchesLen++;

              //  Invalidate this match
              //
              matches[i].pos1 = 0;
              matches[i].len1 = 0;
              matches[i].pos2 = 0;
              matches[i].len2 = 0;
            } else if ((key2->beg <= matches[i].pos2) &&
                       (matches[i].pos2 < key2->end)) {

              //  Trim the begin?  fwdOK, revOK
              //

              int trimLen = key2->end - matches[i].pos2;
              matches[i].pos2 += trimLen;
              matches[i].len2 -= trimLen;

              if (fwd == true)
                matches[i].pos1 += trimLen;
              matches[i].len1 -= trimLen;

              fprintf(stderr, "begin2     key2 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                      key2->axis, key2->beg, key2->end,
                      matches[i].pos1, matches[i].pos1 + matches[i].len1,
                      matches[i].pos2, matches[i].pos2 + matches[i].len2);
            } else if ((key2->beg < (matches[i].pos2 + matches[i].len2)) &&
                       ((matches[i].pos2 + matches[i].len2) <= key2->end)) {

              //  Trim the end?  fwdOK
              //

              int trimLen = matches[i].pos2 + matches[i].len2 - key2->beg;
              matches[i].len1 -= trimLen;

              if (fwd == false)
                matches[i].pos1 += trimLen;
              matches[i].len2 -= trimLen;

              fprintf(stderr, "end2       key2 = %08d %08d %08d    thing = %08d %08d  %08d %08d\n",
                      key2->axis, key2->beg, key2->end,
                      matches[i].pos1, matches[i].pos1 + matches[i].len1,
                      matches[i].pos2, matches[i].pos2 + matches[i].len2);
            } else {
              //fprintf(stderr, "no intersection\n");
              //exit(2);
            }
          } 

          node2 = dict_next(Rint->_il, node2);
        } else {
          node2 = 0L;
        }
      }  //  end of while (node1 || node2)

      //  Print out all the modified matches
      //
      for (u32bit i=0; i<matchesLen; i++) {
        if ((matches[i].len1 > 0) && (matches[i].len2 > 0)) {
          fprintf(stdout, "M %s %s.%d . %s %d %d %d %s %d %d %d\n",
                  S[1], S[2], i,
                  S[4], matches[i].pos1, matches[i].len1, matches[i].ori1 ? 1 : -1,
                  S[8], matches[i].pos2, matches[i].len2, matches[i].ori2 ? 1 : -1);
        }
      }
      fflush(stdout);
      fflush(stderr);
      fprintf(stderr, "----------------------------------------\n");
      fflush(stdout);
      fflush(stderr);
    } else {
      fflush(stdout);
      fflush(stderr);
      fputs(inLine, stdout);
      fflush(stdout);
      fflush(stderr);
    }

    fgets(inLine, 1024, inFile);
  }

  fclose(inFile);



  delete Fcov;
  delete Rcov;
}

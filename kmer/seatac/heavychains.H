#ifndef STRANDPAIR_H
#define STRANDPAIR_H

#include <stdio.h>
#include <stdlib.h>

#include <vector>
using namespace std;

#include "util.h"

//
//  The StrandPair does the heavy chains filtering, while the TheStats
//  collects statistics on all StrandPairs.  TheStats is also
//  responsible for reporting the options used by StrandPair.
//


struct Match {
  int xlo;
  int ylo;
  int xhi;
  int yhi;

  double selfS; // The intrinsic score of the Match.
  double S;     // The computed score of the Match?

  // We need two times the number of dimensions of scores.  That is
  // one score starting from each corner of the bounding box of the
  // space.
  //
  double neS;
  double swS;
  double nwS;
  double seS;

  int filled;  // Is this the same as selfS ?

  char ori;
};



class StrandPair {
public:
  StrandPair(bool   verbose,
             char  *assemblyid1,
             char  *assemblyid2,
             int    maxjump,
             double minscore) {
    beVerbose    = verbose;
    strncpy(assemblyId1, assemblyid1, 31);
    strncpy(assemblyId2, assemblyid2, 31);
    maxJump      = maxjump;
    minScore     = minscore;

    Plen = 0;
    Pmax = 1024;
    P    = new Match [Pmax];

    clear();
  };

  ~StrandPair(void) {
    if (beVerbose > 1)
      fprintf(stderr, "StrandPair::StrandPair()-- delete %s vs %s with %d hits\n", assemblyId1, assemblyId2, Plen);
    delete [] P;
  };

  void addHit(char   direction,
              u32bit id1,
              u32bit xlo,
              u32bit xln,
              u32bit id2,
              u32bit ylo,
              u32bit yln,
              u32bit filled);

  void process(void);

  u64bit print(FILE *outF, u64bit matchid);

  void clear(void) {
    iid1        = ~u32bitZERO;
    iid2        = ~u32bitZERO;

    _next       = 0L;

    sumlen1     = 0.0;
    sumlen2     = 0.0;
    maxlen1     = 0.0;
    maxlen2     = 0.0;
    maxScoreFwd = 0.0; // Maximum forward chain score for the strand pair.
    maxScoreRev = 0.0; // Maximum reverse chain score for the strand pair.

    Plen = 0;
  };

  //  The StrandPairManager (in filter-strandpair.H) is lazy and makes
  //  StrandPairs keep track of the next one.
  //
  StrandPair  *_next;
  StrandPair  *next(void)             { return(_next); };
  void         addNext(StrandPair *n) { _next = n;     };

  u32bit       sequenceIID1(void) { return(iid1); };
  //u32bit       sequenceIID2(void) { return(iid2); };

  double getsumlen1(void) const { return(sumlen1); };
  double getsumlen2(void) const { return(sumlen2); };
  double getmaxlen1(void) const { return(maxlen1); };
  double getmaxlen2(void) const { return(maxlen2); };
  double getmaxScoreFwd(void) const { return(maxScoreFwd); };
  double getmaxScoreRev(void) const { return(maxScoreRev); };

private:
  //  This used to use a vector<Match>, but DPTree wants a pointer to the
  //  array of matches
  //
  int    Plen;
  int    Pmax;
  Match *P;

  void   Padd(Match *m) {
    if (Plen >= Pmax) {
      Pmax *= 2;
      Match *n = new Match [Pmax];
      memcpy(n, P, sizeof(Match) * Plen);
      delete [] P;
      P = n;
    }
    memcpy(P+Plen, m, sizeof(Match));
    Plen++;
  };

private:
  u32bit     iid1;
  u32bit     iid2;
  u32bit     beVerbose;
  char       assemblyId1[32];
  char       assemblyId2[32];
  int        maxJump;          // Default maximum intra-run jump allowed in a good run.
  double     minScore;         // Default minimum of bp filled in a good run.

  // The following are only known after StrandPair::print().

  double     sumlen1;
  double     sumlen2;
  double     maxlen1;
  double     maxlen2;
  double     maxScoreFwd; // Maximum forward chain score for the strand pair.
  double     maxScoreRev; // Maximum reverse chain score for the strand pair.
};






class TheStats {
public:
  TheStats(bool   verbose,
           char  *assemblyid1,
           char  *assemblyid2,
           int    maxjump,
           double minscore) {
    beVerbose    = verbose;
    strncpy(assemblyId1, assemblyid1, 31);  //  Note the cap 'i'
    strncpy(assemblyId2, assemblyid2, 31);
    maxJump      = maxjump;
    minScore     = minscore;

    sumlen1        = 0.0;
    sumlen2        = 0.0;
    sumMaxLen1     = 0.0;
    sumMaxLen2     = 0.0;
    sumMaxScoreFwd = 0.0;
    sumMaxScoreRev = 0.0;
  };

  void add(StrandPair *sp) {
    sumlen1        += sp->getsumlen1();
    sumlen2        += sp->getsumlen2();
    sumMaxLen1     += sp->getmaxlen1();
    sumMaxLen2     += sp->getmaxlen2();
    sumMaxScoreFwd += sp->getmaxScoreFwd();
    sumMaxScoreRev += sp->getmaxScoreRev();
  };

  void show(FILE *outfile) {
    fprintf(outfile, "/assemblyId1=%s\n",         assemblyId1);
    fprintf(outfile, "/assemblyId2=%s\n",         assemblyId2);
    fprintf(outfile, "/heavyMaxJump=%d\n",        maxJump);
    fprintf(outfile, "/heavyMinFill=%f\n",        minScore);
    fprintf(outfile, "/heavySumLen1=%f\n",        sumlen1);
    fprintf(outfile, "/heavySumLen2=%f\n",        sumlen2);
    fprintf(outfile, "/heavySumMaxLen1=%f\n",     sumMaxLen1);
    fprintf(outfile, "/heavySumMaxLen2=%f\n",     sumMaxLen2);
    fprintf(outfile, "/heavySumMaxScoreFwd=%f\n", sumMaxScoreFwd);
    fprintf(outfile, "/heavySumMaxScoreRev=%f\n", sumMaxScoreRev);
  };

private:
  //  Parameters to the filter
  int    beVerbose;
  char   assemblyId1[32];
  char   assemblyId2[32];
  int    maxJump;
  double minScore;

  double sumlen1;
  double sumlen2;
  double sumMaxLen1;
  double sumMaxLen2;
  double sumMaxScoreFwd;
  double sumMaxScoreRev;
};





struct Interval {
  int     lo;
  int     hi;
  double  S;

  Interval() {}; 
  // This is an explicit redefinition of the default constructor.
};

int x_compar(const void *x,const void *y);
int y_compar(const void *x,const void *y);

class DPTree {
  Interval *node;
  Match *match;
  int node_size;
  int match_size; // The number of matches stored in the tree.

  // DP parameters
  int MaxJump;
  
  struct kd_node {
    bool Xy;
    int start,stop; // The indices to define a segment of the vector.
    int intv;       // some index
    kd_node() {};
    // This is an explicit redefinition of the default constructor.

    inline int nmatches() const {return stop-start;}
    inline int midpoint() const {return (start+stop+1)/2;}
    // This is the midpoint of the match list, not the midpoint of the region.

    inline bool lesser_is_node() const {return nmatches() > 2;}
    inline bool greater_is_node() const {return nmatches() > 3;}

    inline bool X() const {return Xy;}
    inline bool Y() const {return !Xy;}

    kd_node lesser_node() const {
      kd_node ret(*this); // This is an explicit invocation of the copy constructor filled with "this".
      // Equivalent to: kd_node ret = (*this);
      ret.intv = ret.intv*2+1;
      ret.stop = midpoint();
      ret.Xy = !ret.Xy;
      return ret;
    }
    kd_node greater_node() const {
      kd_node ret(*this); // This is an explicit invocation of the copy constructor filled with "this".
      // Equivalent to: kd_node ret = (*this);
      ret.intv = ret.intv*2+2;
      ret.start = midpoint();
      ret.Xy = !ret.Xy;
      return ret;
    }
    inline kd_node child_node(int x) const {
      if (x < midpoint())
	return lesser_node();
      else
	return greater_node();
    }

    // root is always real, and everyone else must have at least 2 points
    inline bool is_real() const {return intv==0 || nmatches() > 1;}
  };

  kd_node root() const {
    kd_node ret;
    ret.start = 0;
    ret.stop = match_size;
    ret.intv = 0;
    ret.Xy = true;
    return ret; // Return the structure by value.
  }

public:

  ~DPTree() {
    delete[] node;
  }
  
  DPTree(unsigned siz, Match *p): node(NULL),match(p),
				node_size(0),match_size(siz) {
    node_size = 2;
    for(unsigned sz =  match_size; sz>1 ; sz/=2, node_size *= 2);
    node = new Interval[node_size];
  }

  inline void setParams(int mj) {
    MaxJump = mj;
  }

  double treeScore() {
    init();
    if (match_size > 0)
      privScore(root(),root());

    return node[root().intv].S;
  }

private:

  inline double pairScore(const Match &pl,const Match &ph) const {
    const int dx = ph.xlo - pl.xlo;
    const int dy = ph.ylo - pl.ylo; // causality difference

    const int ix = ph.xlo - pl.xhi;
    const int iy = ph.ylo - pl.yhi;

    const int smaller_jump = (ix < iy)?ix:iy; // will be < 0 if they intersect
    const int larger_jump = (ix < iy)?iy:ix;  // must be < MaxJump for an interaction

    int intersection = smaller_jump * (smaller_jump < 0); 

    return (dx >= 0 && dy >= 0 && larger_jump < MaxJump) 
      * (pl.S + intersection );
  }

  inline bool pruneScore(kd_node f,
			 const Match &p) const {
    int d,jd;

    if (f.X()) {
      jd = p.xlo - node[f.intv].hi;
      d = p.xlo - node[f.intv].lo;
    }
    else {
      jd = p.ylo - node[f.intv].hi;
      d = p.ylo - node[f.intv].lo;
    }
    // returns true if we really need to check this score
    return (d >= 0 && jd < MaxJump) && (node[f.intv].S >= p.S);
  }

  double privScore(kd_node flo,kd_node fhi) {
    // no longer double recursive -- just iterate through fhi
    for(int x = fhi.start; x < fhi.stop; ++x) {
      match[x].S = 0;
      matchScore(flo,match[x]);
      match[x].S += match[x].selfS;

      for(kd_node tmp(fhi); tmp.is_real() ; tmp = tmp.child_node(x)) {
	if (node[tmp.intv].S < match[x].S) node[tmp.intv].S = match[x].S;
      }
    }
    return node[fhi.intv].S;
  }

  double matchScore(kd_node flo,Match &p) {
    double score = 0;

    if ( (flo.X() && node[flo.intv].lo <= p.xlo ||
	  flo.Y() && node[flo.intv].lo <= p.ylo)
	 && pruneScore(flo,p) ) {

      if (flo.greater_is_node())
	score = matchScore(flo.greater_node(),p);
      else
	score = pairScore(match[flo.stop-1],p);
      if (p.S < score) p.S = score;

      if (flo.lesser_is_node())
	score = matchScore(flo.lesser_node(),p);
      else
	score = pairScore(match[flo.start],p);

      if (p.S < score) p.S = score;

    }
    return p.S;
  }

  void init() {
    if (match_size > 0){ sort_nodes(root());}

    int minx=0,miny=0,maxx=0,maxy=0; // initial values will be overwritten
    get_bbox(root(),minx,miny,maxx,maxy);

    for (int i=0; i < node_size; ++i) node[i].S = -1;
    for (int i=0; i < match_size; ++i) match[i].S = -1;
  }

  void sort_nodes(kd_node fs) {
    if (fs.intv >= node_size) {
      fprintf(stderr,"overflow %d %d\n",fs.intv,node_size);
    }

    qsort( match+fs.start, fs.nmatches(), sizeof(Match),
	   (fs.X()?x_compar:y_compar) );

    if (fs.greater_is_node()) sort_nodes(fs.greater_node());
    if (fs.lesser_is_node()) sort_nodes(fs.lesser_node());
  }

  void get_bbox(kd_node fs,int &minx,int &miny,int &maxx,int &maxy) {
    int lminx,lminy,lmaxx,lmaxy;
    int gminx,gminy,gmaxx,gmaxy;

    if (fs.lesser_is_node()) {
      get_bbox(fs.lesser_node(),lminx,lminy,lmaxx,lmaxy);
    }
    else {
      lminx = match[fs.start].xlo;
      lmaxx = match[fs.start].xhi;
      lminy = match[fs.start].ylo;
      lmaxy = match[fs.start].yhi;
    }

    if (fs.greater_is_node()) {
      get_bbox(fs.greater_node(),gminx,gminy,gmaxx,gmaxy);
    }
    else {
      gminx = match[fs.stop-1].xlo;
      gmaxx = match[fs.stop-1].xhi;
      gminy = match[fs.stop-1].ylo;
      gmaxy = match[fs.stop-1].yhi;
    }

    miny = (lminy < gminy)?lminy:gminy;
    minx = (lminx < gminx)?lminx:gminx;
    maxy = (lmaxy > gmaxy)?lmaxy:gmaxy;
    maxx = (lmaxx > gmaxx)?lmaxx:gmaxx;

    if (fs.X()) {
      node[fs.intv].lo = minx;
      node[fs.intv].hi = maxx;
    }
    else {
      node[fs.intv].lo = miny;
      node[fs.intv].hi = maxy;
    }

  }
};





#endif  //  STRANDPAIR_H

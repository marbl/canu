// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Clark Mobarry
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


#ifndef __CYGWIN__
#include <libgen.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;

#include "heavyChainsMod.h"

// The input file must come in "sequence axis pair segments".  

// Assuming Brian's EST mapper file format:
// mkdir tmp
// nice -10 sort -y -T ./tmp -k 3n -k 4n -k 7n -k 8n inpfile.k20.u1.*.matches  > inpfile.k20.u1.matches.sorted
// heavyChains -v -1 id1 -2 id2 -S 100 -J 100000 inpfile.k20.u1.matches.sorted inpfile.k20.u1.matches.heavy
// heavyChains -v -g /assemblyId1=B34 -g /assemblyId2=B33A -g /heavyMinFill=100 -g /heavyMaxJump=100000 B34vsB33A-C16.atac.ckpBeforeHeavyChains B34vsB33A-C16.atac.ckpAfterHeavyChains.031210 > heavyChains.out 2>&1 &
//
// runHeavy.py inpfile.matches

// Ross Lippert stores the command parameters in static globals.

// Proposed interface with seatac.
// so = load shared library (*.so)
// filterObj = new so->heavyChains(options)
//
//    filterObj->addHit(int xlo, int ylo, int xlen, int ylen, int orienation)
// filterObj->filter()
// filterObj->write(FILE *fp)
// delete filterObj
// 


struct Interval {
  int lo,hi;
  double S;
  Interval() {}; 
  // This is an explicit redefinition of the default constructor.
};

// The following would need to parameterized for a general kD tree.
// Could we use one function with a static variable to remember the
// sorting direction?
int x_compar(const void *x,const void *y) {
  const Match &p1=*((const Match*)x);
  const Match &p2=*((const Match*)y);
  if (p1.xhi < p2.xhi) return -1;
  if (p1.xhi > p2.xhi) return  1;
  return 0;
}

int y_compar(const void *x,const void *y) {
  const Match &p1=*((const Match*)x);
  const Match &p2=*((const Match*)y);
  if (p1.yhi < p2.yhi) return -1;
  if (p1.yhi > p2.yhi) return  1;
  return 0;
}

class DPTree {
  Interval *node;
  Match *match;
  int node_size;
  int match_size; // The number of matches stored in the tree.

  // DP parameters
  int MaxJump;
  
  struct kd_node {
    //private:
    bool Xy;
    //public:
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
    if (match_size>0)
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
    assert(intersection <= 0);

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
}; // end class DPTree


///////////////////////////////////////////////////////////

StrandPair::StrandPair(int _VERBOSE, int _MAXJUMP, double _MINSCORE ) {
  // The are needed at initialization.
  VERBOSE = _VERBOSE;
  MAXJUMP = _MAXJUMP;  // Default maximum intra-run jump allowed in a good run.
  MINSCORE = _MINSCORE; // Default minimum of bp filled in a good run.
  old_stra1 = -1;
  old_stra2 = -1;

  // The following are only known after StrandPair::print().
  sumlen1 = 0.;
  sumlen2 = 0.;
  maxlen1 = 0.;
  maxlen2 = 0.;
  maxScoreFwd = 0.; // Maximum forward chain score for the strand pair.
  maxScoreRev = 0.; // Maximum reverse chain score for the strand pair.
}

StrandPair::~StrandPair(void){
  P.clear();
}

void StrandPair::addHit
(
 char   direction,
 u32bit id1,
 u32bit xlo,
 u32bit xln,
 u32bit id2,
 u32bit ylo,
 u32bit yln,
 u32bit filled) {
  
  Match tmp; 
  // ignoring id1 and id2 ?
  // Assert that the data is proper:
  assert(xln > 0);
  assert(yln > 0);
  assert(xlo >= 0);
  assert(ylo >= 0);

  old_stra1 = id1;
  old_stra2 = id2;

  tmp.ori = direction;
  tmp.xlo = xlo;
  tmp.ylo = ylo;
  tmp.filled = filled;

  // convert to a bounding box
  tmp.xhi = tmp.xlo + xln;
  tmp.yhi = tmp.ylo + yln;
  tmp.xhi = tmp.xlo + xln;
  tmp.yhi = tmp.ylo + yln;
  tmp.S = 0; tmp.neS = 0; tmp.nwS = 0; tmp.seS = 0; tmp.swS = 0;

#if 1
  // Use the match lengths to initialize the self scores.
  if (yln < xln) 
    tmp.selfS = yln;
  else
    tmp.selfS = xln;
#else
  tmp.selfS = filled;
#endif
    
  P.push_back(tmp);
}

void StrandPair::clear(void){
  P.clear();
}

size_t StrandPair::size(void) const {
  return P.size();
}
  
void StrandPair::process(void) {
  // new strand pair: begin processing data for the strand pair

#if 1
  if (VERBOSE>0) {
    fprintf(stderr,"HeavyChains: filtering strands %d %d %d\n",old_stra1,old_stra2,P.size());
  }
#endif

  if(P.size() > 0) {
    DPTree *dp = NULL;
    dp = new DPTree(P.size(),&(P[0]));
    dp->setParams(MAXJUMP);

    for(int quadrant=0; quadrant < 4; ++quadrant) {
      if (VERBOSE>1) {
	fprintf(stderr,"HeavyChains: arranging process quadrant %d\n", quadrant);
      }
      switch(quadrant) {
      case 0:	case 2:
	for(unsigned i=0; i < P.size(); ++i) {
	  int swapi;
	  swapi = -P[i].xlo; P[i].xlo = -P[i].xhi; P[i].xhi = swapi;
	}
	break;
      case 1: case 3:
	for(unsigned i=0; i < P.size(); ++i) {
	  int swapi;
	  swapi = -P[i].ylo; P[i].ylo = -P[i].yhi; P[i].yhi = swapi;
	}
	break;
      default:
	exit(127);
      }
	  
      if (VERBOSE>1) {
	fprintf(stderr,"HeavyChains: scoring quadrant\n");
      }
	  
      dp->treeScore();
	  
      if (VERBOSE>1) {
	fprintf(stderr,"HeavyChains: recording scores\n");
      }
      switch(quadrant) {
      case 0:
	for(unsigned i=0; i < P.size(); ++i) P[i].nwS = P[i].S;
	break;
      case 1:
	for(unsigned i=0; i < P.size(); ++i) P[i].swS = P[i].S;
	break;
      case 2:
	for(unsigned i=0; i < P.size(); ++i) P[i].seS = P[i].S;
	break;
      case 3:
	for(unsigned i=0; i < P.size(); ++i) P[i].neS = P[i].S;
	break;
      default:
	exit(127);
      }
	  
      if (VERBOSE>1) {
	fprintf(stderr,"HeavyChains: done quadrant\n");
      }
    }
    // All output information is now in the match records of P.
    delete dp;
  }
} // new strand pair: finished processing the data for the strand pair.

long StrandPair::print(
		       FILE *outF,
		       const int output_in_brians_format,
		       const char * assemblyId1, 
		       const char * assemblyId2,
                       long matchid
 ) {
  fprintf(stderr, "HeavyChainsMod: begin output of the matches\n");
  fprintf(stderr, "HeavyChainsMod: P.size,MAXJUMP,MINSCORE=%d,%d,%lf\n", P.size(), MAXJUMP, MINSCORE);
  fprintf(stderr, "HeavyChainsMod: assemblyId1 = %s\n", assemblyId1);
  fprintf(stderr, "HeavyChainsMod: assemblyId2 = %s\n", assemblyId2);

  if(output_in_brians_format) {
    // Output the results in Brian's format.
    double inc,dec;
    for(unsigned i=0; i < P.size(); ++i) {
      // symmetrize the forward and backward scores.
      inc = P[i].neS + P[i].swS - P[i].selfS;
      dec = P[i].seS + P[i].nwS - P[i].selfS;
      if (inc >= MINSCORE || dec >= MINSCORE) {
	int len1 = (P[i].xhi-P[i].xlo);
	int len2 = (P[i].yhi-P[i].ylo);
	// Now output using Brian's EST mapper format ....
	fprintf(outF,"-%c -e %d %d %d -D %d %d %d -F %d -h %.1f %.1f\n",
		P[i].ori,
		old_stra1, P[i].xlo, len1,
		old_stra2, P[i].ylo, len2, 
		P[i].filled,
		inc, dec);
	sumlen1 += len1;
	sumlen2 += len2;
	maxlen1 = ( maxlen1 > len1 ? maxlen1 : len1);
	maxlen2 = ( maxlen2 > len2 ? maxlen2 : len2);
	maxScoreFwd = (maxScoreFwd > inc ? maxScoreFwd : inc);
	maxScoreRev = (maxScoreRev > dec ? maxScoreRev : dec);
      }
    }
  } else {
    // Output the results in ATAC format.

    double inc, dec;
    for(unsigned i=0; i < P.size(); ++i) {
      // symmetrize the forward and backward scores
      inc = P[i].neS + P[i].swS - P[i].selfS; // forward complement orientations
      dec = P[i].seS + P[i].nwS - P[i].selfS; // reverse complement orientations
      // Each score already contains the self score 
      
      if (inc >= MINSCORE || dec >= MINSCORE) {
	// Now output using ATAC format ....
	int len1 = (P[i].xhi-P[i].xlo);
	int len2 = (P[i].yhi-P[i].ylo);
	matchid += 1;

        fprintf(stderr, "assemblyId1 = %p %s\n", assemblyId1, assemblyId1);
        fprintf(stderr, "assemblyId2 = %p %s\n", assemblyId2, assemblyId2);

        //  Was %Ld, but that's an error.  It crashes Linux.  BPW left in
        //  lots of debugging printf's

	fprintf(outF, "M x H%ld . %s:%d %d %d %d %s:%d %d %d %d > /hf=%.1f /hr=%.1f\n",
		matchid,
		assemblyId1, old_stra1, P[i].xlo, len1, 1,
		assemblyId2, old_stra2, P[i].ylo, len2, (P[i].ori == 'f'? 1 : -1),
		inc, dec );
	sumlen1 += len1;
	sumlen2 += len2;
	maxlen1 = ( maxlen1 > len1 ? maxlen1 : len1);
	maxlen2 = ( maxlen2 > len2 ? maxlen2 : len2);
	maxScoreFwd = (maxScoreFwd > inc ? maxScoreFwd : inc);
	maxScoreRev = (maxScoreRev > dec ? maxScoreRev : dec);
      }
    }
  }
  fprintf(stderr,"HeavyChainsMod: end output of the matches\n");
  
#if 0
  fprintf(stderr,"StrandPair::maxlen1=%f\n", maxlen1);
  fprintf(stderr,"StrandPair::maxlen2=%f\n", maxlen2);
#endif
  
  if (VERBOSE>0) {
    fprintf(stderr,"HeavyChains: finished strands %d %d %lf %lf %lf %lf\n",old_stra1,old_stra2, maxlen1, maxlen2, maxScoreFwd, maxScoreRev);
  }
  //fprintf(outF,"# HeavyChains: finished strands %d %d %lf %lf %lf %lf\n", old_stra1,old_stra2, maxlen1, maxlen2, maxScoreFwd, maxScoreRev);
  return matchid;
}


//////////////////////////////////////////////////

#if 0
void usage(FILE *f) {
  fprintf(f,"usage: %s [-v] [-v] [-h] [-J max_jump_len] [-S min_run_score] inpfile outfile\n",whoami);
}
#endif



TheStats::TheStats(map<string,string> & globals) : _globals(globals) {
  
  sumlen1 = 0.;
  sumlen2 = 0.;
  sumMaxLen1 = 0.;
  sumMaxLen2 = 0.;
  sumMaxScoreFwd = 0.;
  sumMaxScoreRev = 0.;

  // Output the globals for review.
  if (1) {
    map<string,string>::const_iterator mi;
    fprintf(stderr,"Globals:\n");
    for(mi=this->_globals.begin(); mi!=this->_globals.end(); mi++){
      fprintf(stderr,"TheStats: this->_globals[%s]=%s\n", (mi->first).c_str(), (mi->second).c_str());
    }
  }
}

TheStats::~TheStats(void){
#if 0
  if(NULL != this->options) { free(this->options);}
#endif
}

void TheStats::add(StrandPair *sp){
#if 1
  fprintf(stderr,"TheStats::add() beginning\n");
#endif
  this->sumlen1 += sp->sumlen1;
  this->sumlen2 += sp->sumlen2;
  this->sumMaxLen1 += sp->maxlen1;
  this->sumMaxLen2 += sp->maxlen2;
  this->sumMaxScoreFwd += sp->maxScoreFwd;
  this->sumMaxScoreRev += sp->maxScoreRev;
#if 1
  fprintf(stderr,"TheStats::sumMaxLen1=%f\n", this->sumMaxLen1);
  fprintf(stderr,"TheStats::sumMaxLen2=%f\n", this->sumMaxLen2);
#endif
  fprintf(stderr,"TheStats::add() ending\n");
}

void TheStats::print(FILE *outF) const {
  fprintf(stderr,"TheStats::print() beginning\n");
  fprintf(outF,"! format atac 1.0\n");

#if 1
  {
    map<string,string>::const_iterator mi;
    //fprintf(stderr,"Globals:\n");
    for(mi=_globals.begin(); mi!=_globals.end(); mi++){
      fprintf(outF,"/%s=%s\n", (mi->first).c_str(), (mi->second).c_str());
    }
  }
#endif
  
  fprintf(outF, "/heavySumLen1=%lf\n", this->sumlen1);
  fprintf(outF, "/heavySumLen2=%lf\n", this->sumlen2);

  fprintf(outF, "/heavySumMaxLen1=%lf\n", this->sumMaxLen1);
  fprintf(outF, "/heavySumMaxLen2=%lf\n", this->sumMaxLen2);
  
  fprintf(outF, "/heavySumMaxScoreFwd=%lf\n", this->sumMaxScoreFwd);
  fprintf(outF, "/heavySumMaxScoreRev=%lf\n", this->sumMaxScoreRev);
  fprintf(stderr,"TheStats::print() ending\n");
}

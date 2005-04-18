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

// The input file must come in "sequence axis pair segments".  

// Assuming Brian's EST mapper file format:
// mkdir tmp
// nice -10 sort -y -T ./tmp -k 3n -k 4n -k 7n -k 8n inpfile.k20.u1.*.matches  > inpfile.k20.u1.matches.sorted
// heavyChains -v -1 id1 -2 id2 -S 100 -J 100000 inpfile.k20.u1.matches.sorted inpfile.k20.u1.matches.heavy
// heavyChains -v -g /assemblyId1=B34 -g /assemblyId2=B33A -g /heavyMinFill=100 -g /heavyMaxJump=100000 B34vsB33A-C16.atac.ckpBeforeHeavyChains B34vsB33A-C16.atac.ckpAfterHeavyChains.031210 > heavyChains.out 2>&1 &
//
// runHeavy.py inpfile.matches

// Ross Lippert stores the command parameters in static globals.

static char *whoami = NULL;
static int VERBOSE = 1;
static char *inpfile=NULL, *outfile=NULL;
static int BRIANS_FORMAT = 1;
#define BUFFERSIZE 1024

static map<string,string> globals;
static char *assemblyId1 = NULL;
static char *assemblyId2 = NULL;
static int MAXJUMP = 100000;  // Default maximum intra-run jump allowed in a good run.
static double MINSCORE = 100.; // Default minimum of bp filled in a good run.

struct Interval {
  int lo,hi;
  double S;
  Interval() {}; 
  // This is an explicit redefinition of the default constructor.
};

struct Match {
  int xlo,ylo;
  int xhi,yhi;

  double selfS; // The intrinsic score of the Match.
  double S;     // The computed score of the Match?

  // We need two times the number of dimensions of scores.  That is
  // one score starting from each corner of the bounding box of the
  // space.
  double neS;
  double swS;
  double nwS;
  double seS;
  int filled; // Is this the same as selfS ?
  char ori;
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

void usage(FILE *f) {
  fprintf(f,"usage: %s [-v] [-v] [-h] [-J max_jump_len] [-S min_run_score] inpfile outfile\n",whoami);
}

static int get_args(int argc, char **argv) {
  int errflg = 0;
#ifndef __CYGWIN__
  whoami = strdup(basename(*argv));
#endif
  argv++;
#define OPT_ARG ( (argv[0][1])?(++(*argv)):(*(++argv)) )
  while(!errflg && *argv) {
    if (**argv == '-') {
      (*argv)++;
      while (!errflg && **argv) {
        switch (**argv) {
        case 'b':
          BRIANS_FORMAT = atoi(OPT_ARG);
          goto loopin;
        case 'v':
          VERBOSE++;
          goto loopin;
          //case 'o':
          //outfile = strdup(OPT_ARG);
          //goto loopout;
        case 'g':
          {
            char * buffer = OPT_ARG;
            char * key1 = strtok(buffer,"/=");
            //if(key1 != NULL) { fprintf(stderr,"key1=<%s> %d\n", key1, strlen(key1));}
            char * key2 = strtok(NULL,"/=");
            //if(key2 != NULL) { fprintf(stderr,"key2=<%s> %d\n", key2, strlen(key2));}
            if(key1 != NULL && key2 != NULL) {
              globals[string(key1)] = string(key2);
            } else {
              fprintf(stderr,"Failed to parse ATAC global: %s\n", OPT_ARG);
              assert(0);
            }
          }
          goto loopout;
        case '1':
          globals[string("assemblyId1")] = string(OPT_ARG);
          goto loopout;
        case '2':
          globals[string("assemblyId2")] = string(OPT_ARG);
          goto loopout;
        case 'S':
          globals[string("heavyMinFill")] = string(OPT_ARG);
          goto loopout;
        case 'J':
          globals[string("heavyMaxJump")] = string(OPT_ARG);
          goto loopout;
        case 'h':
          usage(stdout);
          exit(0);
          goto loopin;
        default:
          // std::cerr << whoami << ": unknown flag '-" << *argv <<"'"<< std::endl;
#ifndef __CYGWIN__
          fprintf(stderr,"%s : unknown flag '-%s'\n", whoami, *argv);
#endif
          errflg++;
        }
      loopin:
        (*argv)++;
      }
    loopout:
      argv++;
    }
    else {
      if (!inpfile) {inpfile = strdup(*argv); ++argv;} else { errflg++;}
      if (!outfile) {outfile = strdup(*argv); ++argv;} else { errflg++;}
    }
  }
  return errflg;
}

int main(int argc, char *argv[]) {
  if (get_args(argc,argv)) {
    usage(stderr);
    exit(1);
  }

  FILE *inF = stdin;
  FILE *outF= stdout;
  if (inpfile) {
    if (!(inF = fopen(inpfile,"r"))) {
      fprintf(stderr,"unable to open %s\n",inpfile);
      exit(2);
    }
  }
  if (outfile) {
    if (!(outF = fopen(outfile,"w"))) {
      fprintf(stderr,"unable to open %s\n",outfile);
      exit(3);
    }
  }

#if 1
  fprintf(stderr,"The input strings are <%s> <%s> <%s> <%s>\n",
          globals["assemblyId1"].c_str(), globals["assemblyId2"].c_str(),
          globals["heavyMinFill"].c_str(), globals["heavyMaxJump"].c_str());
  // Note that even this print statement set the key value to non-null!
#endif
  
  if(globals.find("assemblyId1") != globals.end()){
    fprintf(stderr,"setting assemblyId1\n");
    assemblyId1 = strdup((globals["assemblyId1"]).c_str());
  } else {
    fprintf(stderr,"error: assemblyId1 is not set\n");
    exit(3);
  }
    
  if(globals.find("assemblyId2") != globals.end()){
    fprintf(stderr,"setting assemblyId2\n");
    assemblyId2 = strdup((globals["assemblyId2"]).c_str());
  } else {
    fprintf(stderr,"error: assemblyId2 is not set\n");
    exit(3);
  }
  if(globals.find("heavyMinFill") != globals.end()){
    fprintf(stderr,"setting MINSCORE\n");
    MINSCORE = atof((globals["heavyMinFill"]).c_str());
  }
  if(globals.find("heavyMaxJump") != globals.end()){
    fprintf(stderr,"setting MAXJUMP\n");
    MAXJUMP  = atoi((globals["heavyMaxJump"]).c_str());
  }

  if (VERBOSE>1) {
    map<string,string>::const_iterator mi;
    fprintf(stderr,"Globals:\n");
    for(mi=globals.begin(); mi!=globals.end(); mi++){
      fprintf(stderr," globals[%s]=%s\n", (mi->first).c_str(), (mi->second).c_str());
    }
  }

#if 1
  fprintf(stderr,"The input strings are <%s> <%s> <%s> <%s>\n",
          globals["assemblyId1"].c_str(), globals["assemblyId2"].c_str(),
          globals["heavyMinFill"].c_str(), globals["heavyMaxJump"].c_str());
  // Note that even this print statement set the key value to non-null!
#endif
  
  fprintf(stderr,"The input parameters are assemblyId1=%s,assemblyId2=%s,MINSCORE=%f,MAXJUMP=%d\n",
          assemblyId1, assemblyId2, MINSCORE, MAXJUMP);

  Match tmp;
  std::vector<Match> P;

  int  old_stra1=-1, old_stra2=-1; // True strand ordinals are non-negative.

  long long matchid = 0;
  double sumlen1 = 0.;
  double sumlen2 = 0.;
  double sumMaxLen1 = 0.;
  double sumMaxLen2 = 0.;
  double sumMaxScoreFwd = 0.;
  double sumMaxScoreRev = 0.;
  
  int endOfInput = 0;
  P.clear();
  char linebuffer[BUFFERSIZE] = {0};
  
  while( !endOfInput ) {
    endOfInput = 1;

    int  new_stra1=-1, new_stra2=-1;

    if(NULL != fgets(linebuffer,BUFFERSIZE,inF)) {
      //if(PANICDEBUG){ printf("linebuffer=%s",linebuffer);}
      int xln, yln;
      if(linebuffer[0] == '-' || linebuffer[0] == 'M'){
        if(linebuffer[0] != 'M') {
          // Now input using Brian's EST mapper format ....
          int scanc = sscanf(linebuffer,
                             "-%c -e %d %d %d -D %d %d %d -F %d\n",
                             &tmp.ori,
                             &new_stra1,&tmp.xlo,&xln,&new_stra2,&tmp.ylo,&yln,
                             &tmp.filled
                             );
          assert(tmp.ori == 'f' || tmp.ori == 'r');
          assert(scanc >= 7);
          endOfInput = 0;
	} else {
#if 1
	  // Now input using ATAC format ....
	  // M r r6379 . B33:9 134205509 488863 1 R27M4:9 128704627 484912 -1
	  char classCode, subtype;
	  char selfId[100], parentId[100];
          char new_ass1[100], new_ass2[100];
	  int xfl, yfl;
          int scanc = sscanf(linebuffer,
                             "%c %c %s %s "
                             "%s %d %d %d %s %d %d %d\n",
                             &classCode, &subtype, selfId, parentId,
                             new_ass1,&tmp.xlo,&xln,&xfl,
                             new_ass2,&tmp.ylo,&yln,&yfl);

          // Check if new_ass1 equals assemblyId1 etc.
          assert(scanc >= 12);
#if 0
          printf("classCode=%c\n", classCode);
          printf("subtype  =%c\n", subtype);
          printf("selfId   =%s\n", selfId);
          printf("parentId =%s\n", parentId);
          printf("new_ass1 =%s\n", new_ass1);
          printf("xfl      =%d\n", xfl);
          printf("new_ass2 =%s\n", new_ass2);
          printf("yfl      =%d\n", yfl);
#endif
          assert(xfl == 1 || xfl == -1);
          assert(yfl == 1 || yfl == -1);
          tmp.ori = (xfl == yfl ? 'f' : 'r');
          {
            char *p;
            for(p = new_ass1; *p != 0; p++){ if(*p == ':'){ p++; break;}}
            sscanf(p,"%d", &new_stra1);
          }
          {
            char *p;
            for(p = new_ass2; *p != 0; p++){ if(*p == ':'){ p++; break;}}
            sscanf(p,"%d", &new_stra2);
          }
#else
          tmp.ori = 'f'; tmp.xlo = 0; tmp.ylo = 0; xln = 8; yln = 9;
#endif
          endOfInput=0;
        }
#if 0
        printf("endOfInput=%d\n", endOfInput);
        printf("new_stra1 =%d\n", new_stra1);
        printf("new_stra2 =%d\n", new_stra2);
        printf("tmp.ori  =%c\n", tmp.ori);
        printf("tmp.xlo  =%d\n", tmp.xlo);
        printf("xln      =%d\n", xln);
        printf("tmp.ylo  =%d\n", tmp.ylo);
        printf("yln      =%d\n", yln);
#endif        
#if 0
        if(!endOfInput){
          // Assert that the data came in with the axis numerically or lexigraphically sorted:
          assert((old_stra1 <= new_stra1) ||
                 (old_stra1 == new_stra1 && old_stra2 <= new_stra2));
        }
#endif
        
        // Assert that the data is proper:
        assert(xln > 0);
        assert(yln > 0);
        assert(tmp.xlo >= 0);
        assert(tmp.ylo >= 0);

        // convert to a bounding box
        tmp.xhi = tmp.xlo + xln;
        tmp.yhi = tmp.ylo + yln;
        tmp.S = 0; tmp.neS = 0; tmp.nwS = 0; tmp.seS = 0; tmp.swS = 0;
        // Use the match lengths to initialize the self scores.
        if (yln < xln) 
          tmp.selfS = yln;
        else
          tmp.selfS = xln;
        
      } else if(linebuffer[0] == '#') {
        fprintf(outF,"%s",linebuffer);
        endOfInput = 0;
      } else if(linebuffer[0] == '/') {
        //fprintf(outF,"%s",linebuffer);
        char * key1 = strtok(linebuffer,"/=");
        //if(key1 != NULL) { fprintf(stderr,"key1=<%s> %d\n", key1, strlen(key1));}
        char * key2 = strtok(NULL,"\n");
        //if(key2 != NULL) { fprintf(stderr,"key2=<%s> %d\n", key2, strlen(key2));}
        if(key1 != NULL && key2 != NULL) {
          globals[string(key1)] = string(key2);
        }
        endOfInput = 0;
      } else if(linebuffer[0] == '!') {
        fprintf(outF,"%s",linebuffer);
        endOfInput = 0;
      } else {
        fprintf(stderr,"Unrecognized input\n");
        fprintf(outF,"%s",linebuffer);
        endOfInput = 0;
      }
    }

#if 0
    printf("JOE endOfInput=%d old_stra1=%d old_stra2=%d new_stra1=%d new_stra2=%d P.size=%d\n",
           endOfInput, old_stra1, old_stra2,  new_stra1, new_stra2, P.size());
#endif    

    if ( (new_stra1 != old_stra1 || new_stra2 != old_stra2 || endOfInput)
	 && P.size()>0 ) {
      // new strand pair: begin processing data for the strand pair
      if (VERBOSE>10) {
        map<string,string>::const_iterator mi;
        fprintf(stderr,"Globals for new strand pair:\n");
        for(mi=globals.begin(); mi!=globals.end(); mi++){
          fprintf(stderr," globals[%s]=%s\n", (mi->first).c_str(), (mi->second).c_str());
        }
      }

      if (VERBOSE>0) {
	fprintf(stderr,"HeavyChains: filtering strands %d %d %d\n",old_stra1,old_stra2,P.size());
      }
      double maxlen1 = 0.;
      double maxlen2 = 0.;
      double maxScoreFwd = 0.; // Maximum forward chain score for the strand pair.
      double maxScoreRev = 0.; // Maximum reverse chain score for the strand pair.

      if(1) {
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

      if(0){
	// Output the results in Brian's format.
	double inc,dec;
	for(unsigned i=0; i < P.size(); ++i) {
	  // symmetrize the forward and backward scores.
	  inc = P[i].neS + P[i].swS - P[i].selfS;
	  dec = P[i].seS + P[i].nwS - P[i].selfS;
	  if (inc >= MINSCORE || dec >= MINSCORE) {
	    // Now output using Brian's EST mapper format ....
	    fprintf(outF,"-%c -e %d %d %d -D %d %d %d -F %d -h %.1f %.1f\n",
		    P[i].ori,
		    old_stra1, P[i].xlo, (P[i].xhi-P[i].xlo),
		    old_stra2, P[i].ylo, (P[i].yhi-P[i].ylo),
                    P[i].filled,
                    inc, dec);
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
	    fprintf(outF,"M x H%Ld . %s:%d %d %d %d %s:%d %d %d %d > /hf=%.1f /hr=%.1f\n",
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

      P.clear();

      if (VERBOSE>0) {
	fprintf(stderr,"HeavyChains: finished strands %d %d %lf %lf %lf %lf\n",old_stra1,old_stra2, maxlen1, maxlen2, maxScoreFwd, maxScoreRev);
      }
      //fprintf(outF,"# HeavyChains: finished strands %d %d %lf %lf %lf %lf\n", old_stra1,old_stra2, maxlen1, maxlen2, maxScoreFwd, maxScoreRev);
      sumMaxLen1 += maxlen1;
      sumMaxLen2 += maxlen2;
      sumMaxScoreFwd += maxScoreFwd;
      sumMaxScoreRev += maxScoreRev;
      
    } // new strand pair: finished processing the data for the strand pair.

    if(linebuffer[0] == '-' || linebuffer[0] == 'M') {
      // If we just processed a new point.
      P.push_back(tmp);
      old_stra1 = new_stra1;
      old_stra2 = new_stra2;
    }

  }

  {
    char tmpbuffer[100];
    // #include <strstream>
    sprintf(tmpbuffer,"%lf",sumlen1);
    globals[string("heavySumLen1")] = string(tmpbuffer);
    sprintf(tmpbuffer,"%lf",sumlen2);
    globals[string("heavySumLen2")] = string(tmpbuffer);
    
    sprintf(tmpbuffer,"%lf",sumMaxLen1);
    globals[string("heavySumMaxLen1")] = string(tmpbuffer);
    sprintf(tmpbuffer,"%lf",sumMaxLen2);
    globals[string("heavySumMaxLen2")] = string(tmpbuffer);

    sprintf(tmpbuffer,"%lf",sumMaxScoreFwd);
    globals[string("heavySumMaxScoreFwd")] = string(tmpbuffer);
    sprintf(tmpbuffer,"%lf",sumMaxScoreRev);
    globals[string("heavySumMaxScoreRev")] = string(tmpbuffer);
  }
  
  {
    map<string,string>::const_iterator mi;
    //fprintf(stderr,"Globals:\n");
    fprintf(outF,"! format atac 1.0\n");
    for(mi=globals.begin(); mi!=globals.end(); mi++){
      fprintf(outF,"/%s=%s\n", (mi->first).c_str(), (mi->second).c_str());
    }
  }

  if (inpfile) fclose(inF);
  if (outfile) fclose(outF);
}

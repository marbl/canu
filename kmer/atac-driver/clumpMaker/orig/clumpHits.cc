#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <LEDA/string.h>
#include <LEDA/array.h>
#include <LEDA/tuple.h>
//#include <LEDA/sortseq.h>
#include <LEDA/stream.h>
#include <LEDA/param_handler.h>

#include "tClumpHits.h"

#if (defined i386 || defined __POWERPC__)
#define LD "%lld"
#else
#define LD "%ld"
#endif

static int chainable(const array<tClumpHit>&hits,int i,int j,int maxjump){

  //cerr <<"chainability test:\n" << hits[i] << hits[j];
  if(// hits are to different chromosomes
     (hits[i].get_rCh() != hits[j].get_rCh() ||hits[i].get_qCh() != hits[j].get_qCh())||
     // hits are not similarly oriented
     ( hits[i].get_qOri() != hits[j].get_qOri() ) ||
     // hits are too far apart on reference axis
     ( hits[j].get_rBeg()-hits[i].get_rEnd() > (int64)maxjump) || 
     // hits are too far apart on query axis
     ( ((int64)hits[i].get_qOri())*(hits[j].get_qFirst()-hits[i].get_qLast()) > (int64)maxjump) ||
     // hits are "out of order"
     ( ((int64)hits[i].get_qOri())*(hits[j].get_qFirst()-hits[i].get_qFirst()) < (int64) 0)
     ){
    return 0;
  }
  
  //cerr << "Passed\n";
  return 1;
}

static void score_all_hits(array<tClumpHit> & hits,int clumpcost,int maxjump,int &bestEnd,int num_hits){

  // location of best score so far (to which we point whenever starting a new clump)
  bestEnd=-1;
  // best score so far
  int64 bestScore=(int64)-clumpcost;

  // best scores so far internal to a reference unit (scaffold, chromosome, etc)
  int bestEndThis=-1;
  int bestScoreThis=(int64)-clumpcost;

  int rptIval=10000;
  // furthest back still accessible ...
  int furthest_back=0;
  
  // tmp variabile to avoid computing twice for tests then assignments
  int64 tmpscore;

  for(int i=0;i<num_hits;i++){

    if(((int)(i/rptIval))*rptIval == i)
        cerr << "Working on best solution through hit " << i << endl << hits[i] << endl;


    if(i==0||hits[i].get_rCh()!=hits[i-1].get_rCh()){
      bestEnd=bestEndThis; // best of previous reference unit
      bestEndThis=i-1;
      bestScoreThis=(int64)-clumpcost;
    }

    // find best way of using this as start of a new clump
    if(bestEndThis>=0&&bestScoreThis>=0){ // start new clump that is not the first for this reference unit
      hits[i].set_scoreStart(hits[i].get_rlen()
			     +bestScoreThis
			     //			     +hits[bestEnd].get_bestScore()
			     -(int64) clumpcost
			     );
      hits[i].set_bestStart(bestEndThis);
    }else{ // clump would be first (to be used) for this reference unit
      hits[i].set_scoreStart(hits[i].get_rlen()
			     -(int64) clumpcost);
      hits[i].set_bestStart(bestEnd);
    }

    // find best way of extending a clump, if any
    if(furthest_back<i){
      int64 cutoff=hits[i].get_rBeg()-(int64)maxjump;

#undef DEBUG_ADVANCE
#ifdef DEBUG_ADVANCE
      cerr << "Checking to advance furthest back [" << furthest_back << "] relative to " << i << " ... cutoff = " << cutoff << " end " << hits[furthest_back].get_rEnd() << "\n";
#endif

      while(hits[furthest_back].get_rCh() != hits[i].get_rCh()||hits[furthest_back].get_rEnd()<cutoff){
	furthest_back++;

#ifdef DEBUG_ADVANCE
	int64 tmp=hits[furthest_back-1].get_rEnd();
	cerr << "Advanced furthest back to " << furthest_back;
	cerr << " passing " << tmp << endl;
#endif

      }
    }
    if(((int)(i/rptIval))*rptIval == i)
      cerr << "Furthest back is " << furthest_back << " " << hits[furthest_back] << endl;
    int64 extendScore=-(int64)clumpcost;
    int extendprev=-1;
    for(int j=furthest_back;j<i;j++){
      if(chainable(hits,j,i,maxjump)){
	tmpscore = hits[j].get_bestScore()+hits[i].get_rlen();
	if(extendScore < tmpscore){
	  extendScore=tmpscore;
	  extendprev=j;
	}
      }
    }
    hits[i].set_scoreExtend(extendScore);
    hits[i].set_bestExtend(extendprev);

    //figure out whether this is a new best ...
    tmpscore = hits[i].get_bestScore();
    if(tmpscore > bestScoreThis){
      bestScoreThis = tmpscore;
      bestEndThis = i;
    }
  }
  bestEnd=bestEndThis;
}

static void mark_clumps(array<tClumpHit>& hits,int bestEnd){
  int clump=0;
  int end=bestEnd;
  while(end>=0){
#undef DEBUG_MARK
#ifdef DEBUG_MARK
    cerr << "Setting clump to " << clump << " for hit " << end << " " << hits[end];
#endif
    hits[end].set_clump(clump);
#ifdef DEBUG_MARK
    cerr << "Testing scores extend " << hits[end].get_scoreExtend()
	 << " start " << hits[end].get_scoreStart() << endl;
#endif
    if(hits[end].get_scoreExtend()>hits[end].get_scoreStart()){
      end=hits[end].get_extendBest();
    }else{
      end=hits[end].get_startBest();
#ifdef DEBUG_MARK
      cerr << "Chose as start\n";
#endif
      clump++;
    }
  }
}

int main (int argc, char *argv[]){

  leda_param_handler P(argc, argv, ".last",false);
  P.add_parameter("penalty for clump start:-c:int:50000");
  P.add_parameter("max jump between consistent hits in a clump:-j:int:200000");
  P.add_parameter("initial number of hits to allocate:-n:int:10000");
  P.add_parameter("whether to swap reference and query:-s:int:0");
  P.init_all();

  int clumpcost,maxjump,initHitAlloc,swap;
  P.get_parameter("-c",clumpcost);
  P.get_parameter("-j",maxjump);
  P.get_parameter("-n",initHitAlloc);
  P.get_parameter("-s",swap);

  array<tClumpHit> hits(initHitAlloc);

  //  cerr << "Size of hit: " << sizeof(tClumpHit) << endl;
  int num_hits=0;
  string line;
  while(line.read_line(),line.length()>0){
    char M[20],t[20],runID[20],dot[20],r[40],q[40];
    int64 rb,qb,re,qe;
    int rlen,rOri,qlen,qOri;
    if(sscanf((const char *)line,"%s %s %s %s %s " LD " %d %d %s " LD " %d %d",
	      M,t,runID,dot,r,&rb,&rlen,&rOri,q,&qb,&qlen,&qOri)!=12){
      cerr << "PANIC: could not read hit from following input line:\n"
	   << line << endl;
      return(-1);
    }

    /*    cerr << "[" << M << "] [" 
	 << t << "] ["  
	 << runID << "] [" 
	 << dot << "] [" 
	 << r << "] [" 
	 << rb << "] [" 
	 << rlen << "] [" 
	 << rOri << "] [" 
	 << q << "] [" 
	 << qb << "] [" 
	 << qlen << "] [" 
	 << qOri << "]" << endl; */

    re=rb+rlen;
    qe=qb+qlen;
    assert(rOri==1);
    if(swap){
      int64 tmp;
      char tmps[40];
      
      tmp=qb;
      qb=rb;
      rb=tmp;

      tmp=qe;
      qe=re;
      re=tmp;

      strcpy(tmps,q);
      strcpy(q,r);
      strcpy(r,tmps);
    }
    if(qOri==-1){
      int64 tmp=qb;
      qb=qe;
      qe=tmp;
    }
    hits[num_hits++].reset(string(r),string(q),rb,re,qb,qe,string(line));
    if(num_hits==hits.size()){
      hits.resize(num_hits*2);
    }
  }
  cerr << "Read " << num_hits << " hits\n";

  hits.resize(num_hits);
  cerr << "Resized\n";

  hits.sort(clumpHitCompareR);
  cerr << "Sorted\n";

  //cerr << hits;

  int bestEnd=-1;
  score_all_hits(hits,clumpcost,maxjump,bestEnd,num_hits);
  cerr << "Scored\n";

  mark_clumps(hits,bestEnd);
  cerr << "Marked\n";

  for(int i=0;i<num_hits;i++){
    printf("%s # %d\n",(const char*)(hits[i].get_hitDef()),hits[i].get_clump());
    //cout << hits[i].hitDef << " # " << hits[i].get_clump() << endl;
  }

  return 0;

}

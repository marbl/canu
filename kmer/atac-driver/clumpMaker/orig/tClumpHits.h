#include <stream.h>
#include <LEDA/string.h>
#include <LEDA/stream.h>

#ifndef max
#define max(a,b) ((a)> (b) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) ((a)< (b) ? (a) : (b))
#endif

#ifdef __POWERPC__
#define int64 long long
#include <stdint.h>
#else
#ifdef i386
#define int64 long long
#else
#define int64 long
#endif
#endif

class tClumpHit {
 public:
  tClumpHit(){
    qCh="";
    rCh="";
    rBeg=-1;
    rEnd=-1;
    qBeg=-1;
    qEnd=-1;
    bestStart=-1;
    bestExtend=-1;
    scoreStart=0;
    scoreExtend=0;
    hitDef = "";
    clump=-1;
  }
  tClumpHit(string rc, string qc,int64 rb,int64 re, int64 qb, int64 qe){
    rCh=rc;
    qCh=qc;
    rBeg=rb;
    rEnd=re;
    qBeg=qb;
    qEnd=qe;
    bestStart=-1;
    bestExtend=-1;
    scoreStart=0;
    scoreExtend=0;
    hitDef = "";
    clump=-1;
  }
  ~tClumpHit(){}
  
  void reset(string rc, string qc,int64 rb,int64 re, int64 qb, int64 qe,const string & hit){
    rCh=rc;
    qCh=qc;
    rBeg=rb;
    rEnd=re;
    qBeg=qb;
    qEnd=qe;
    bestStart=-1;
    bestExtend=-1;
    scoreStart=0;
    scoreExtend=0;
    hitDef=hit;
    clump=-1;
  }

  int64 get_rBeg() const {return rBeg;}
  int64 get_rEnd() const {return rEnd;}
  int64 get_rlen()const{return rEnd-rBeg;}

  int64 get_qBeg() const {return qBeg;}
  int64 get_qEnd() const {return qEnd;}
  int64 get_qFirst()const{return min(qBeg,qEnd);}
  int64 get_qLast()const{return max(qBeg,qEnd);}
  int get_qOri()const{if(qBeg<qEnd){return 1;}else{return -1;}}

  void set_scoreExtend(int64 s) {scoreExtend=s;}
  int64 get_scoreExtend()const {return scoreExtend;}
  void set_scoreStart(int64 s) {scoreStart=s;}
  int64 get_scoreStart()const {return scoreStart;}
  void set_bestStart(int b){bestStart=b;}
  int get_startBest()const {return bestStart;}
  void set_bestExtend(int b){bestExtend=b;}
  int get_extendBest()const {return bestExtend;}
  int64 get_bestScore()const{return max(scoreStart,scoreExtend);}

  void set_clump(int c){clump=c;}
  int get_clump()const{return clump;}
  string get_rCh()const{return rCh;}
  string get_qCh()const{return qCh;}
  string get_hitDef()const{return hitDef;}

  friend ostream & operator<<(ostream&out, const tClumpHit h){
    //    out << h.hitDef << endl;
    out << h.get_hitDef();
    return out;
  }

 private:
  string rCh;
  string qCh;
  int64 rBeg;
  int64 rEnd;
  int64 qBeg;
  int64 qEnd;
  int64 scoreStart;
  int bestStart;
  int64 scoreExtend;
  int bestExtend;
  int clump;
  string hitDef;
};

int clumpHitCompareR(const tClumpHit & a, const tClumpHit & b);
int clumpHitCompareQ(const tClumpHit & a, const tClumpHit & b);

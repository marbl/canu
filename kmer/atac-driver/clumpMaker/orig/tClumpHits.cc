#include <string.h>
#include "tClumpHits.h"

//#define RETURN(x) {cerr << "returning " << x << endl;return x;}
#define RETURN(x) {return x;}

int clumpHitCompareR(const tClumpHit & a, const tClumpHit & b){

  if(a.get_rCh() > b.get_rCh())RETURN(1);
  if(a.get_rCh() < b.get_rCh())RETURN(-1);
  if(a.get_rBeg() > b.get_rBeg())RETURN(1);
  if(a.get_rBeg() < b.get_rBeg())RETURN(-1);
  if(a.get_rEnd() > b.get_rEnd())RETURN(1);
  if(a.get_rEnd() < b.get_rEnd())RETURN(-1);
  if(a.get_qCh() > b.get_qCh())RETURN(1);
  if(a.get_qCh() < b.get_qCh())RETURN(-1);
  int64 astart = a.get_qFirst();
  int64 bstart = b.get_qFirst();
  if(astart > bstart)RETURN(1);
  if(astart < bstart)RETURN(-1);
  int64 aend = a.get_qLast();
  int64 bend = b.get_qLast();
  if(aend > bend)RETURN(1);
  if(aend < bend)RETURN(-1);
  RETURN(0);
}

int clumpHitCompareQ(const tClumpHit & a, const tClumpHit & b){
  if(a.get_qCh() > b.get_qCh())RETURN(1);
  if(a.get_qCh() < b.get_qCh())RETURN(-1);
  int64 astart = a.get_qFirst();
  int64 bstart = b.get_qFirst();
  if(astart > bstart)RETURN(1);
  if(astart < bstart)RETURN(-1);
  int64 aend = a.get_qLast();
  int64 bend = b.get_qLast();
  if(aend > bend)RETURN(1);
  if(aend < bend)RETURN(-1);
  if(a.get_rCh() > b.get_rCh())RETURN(1);
  if(a.get_rCh() < b.get_rCh())RETURN(-1);
  if(a.get_rBeg() > b.get_rBeg())RETURN(1);
  if(a.get_rBeg() < b.get_rBeg())RETURN(-1);
  if(a.get_rEnd() > b.get_rEnd())RETURN(1);
  if(a.get_rEnd() < b.get_rEnd())RETURN(-1);
  RETURN(0);
}



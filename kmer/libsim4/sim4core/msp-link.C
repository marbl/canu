#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mspManager.H"
#include "sim4defines.H"
#include "exon.H"


//
//  Returns an exon_list
//
//  This can be used for both linking and relinking
//




//  This is a hack to make it compile.
//
//  It really should go into the exonManager.
//
#if 0
Exon_ptr
new_exon(int f1, int f2, int t1, int t2, int len, int edist, int flag, Exon_ptr next)
{
  Exon_ptr newthing = (Exon_ptr )ckalloc(sizeof(struct exon));

  newthing->next_exon = next;

  newthing->frGEN = f1;
  newthing->frEST = f2;
  newthing->toGEN = t1;
  newthing->toEST = t2;

  newthing->ori = 'U';
  newthing->length = (len < 0) ? (t2-f2+1) : len;

  newthing->edist = edist;
  newthing->flag = flag;

  newthing->percentID = 0;
  newthing->alignmentLength = 0;
  newthing->numMatches = 0;
  newthing->numNs = 0;
  newthing->numInDel = 0;
  newthing->numEdits = 0;

  return newthing;
}
#endif




#define  min(x,y)        ((x>y) ? (y):(x))
#define  max(x,y)        ((x<y) ? (y):(x))



static
int
get_edist(int f1, int f2,
          int t1, int t2,
          char *seq1,
          char *seq2) {
  char *s1, *s2, *q1, *q2;
  int dist=0;

  s1 = seq1+f1+1;   /* bc at this stage, the msp pos do not have added +1 */
  s2 = seq2+f2+1;
  q1 = seq1+t1+1;
  q2 = seq2+t2+1;

  while (s1<=q1 && s2<=q2) {
    dist += (*s1!=*s2);
    s1++;
    s2++;
  } 

  return dist;
}



Exon_ptr
mspManager::link(int weight, int drange,
                 int offset1, int offset2,
                 int flag, int relinkFlag,
                 char *s1, char *s2) {

  //
  //  Assumes the exon list is cleared
  //

  //fprintf(stderr, "relink weight  = %d\n", weight);
  //fprintf(stderr, "number of MSPs = %d\n", _numMSPs);

#if ABORT_EXPENSIVE
  if ((_cDNALength > 0) &&
      (_mspLimitAbsolute > 0) && (_mspLimitAbsolute < _numMSPs) &&
      (_mspLimitPercent > 0.0) && (_mspLimitPercent * _cDNALength < _numMSPs)) {
    _tooManyMSPs = true;
    return(0L);
  }
#endif

  ////////////////////////////////////////
  //
  //last_msp = mspManager->link(rs.weight);
  //
  //last_msp = link_msps(msp, numMSPs, rs.weight);
  //
  int f1, f2, best, diag, diff_diag, best_sc, tryval;

  best    = -1;
  best_sc = INT_MIN;

  for (u32bit i = 0; i < _numMSPs; ++i) {
    f1 = _allMSPs[i].pos1;      /* start position in seq1 */
    f2 = _allMSPs[i].pos2;      /* start position in seq2 */
    diag = f1 - f2;
    _allMSPs[i].prev = -1;
    _allMSPs[i].linkingScore = 0;
    for (u32bit j = 0; j < i; ++j) {
      //  12 == default word size.  A Magic Value.
      int vL = L; 
      if ((_allMSPs[i].pos2 + _allMSPs[i].len - _allMSPs[j].pos2 - _allMSPs[j].len > 2 * 12) &&
          (_allMSPs[i].pos2 - _allMSPs[j].pos2 > 2 * 12))
        vL *= 2;
                        
      diff_diag = diag - _allMSPs[j].pos1 + _allMSPs[j].pos2;
      if (diff_diag < - drange ||
          (diff_diag >  drange && diff_diag < MIN_INTRON) ||
          (_allMSPs[j].pos2+_allMSPs[j].len-1-f2>vL) ||
          (_allMSPs[j].pos1+_allMSPs[j].len-1-f1>vL))
        continue;

      int n = abs(diff_diag);
      tryval = _allMSPs[j].linkingScore - n;
      if (relinkFlag)
        tryval = _allMSPs[j].linkingScore - ((n <= 100000) ? n : (100000+(int)(10*log((double)(n-100000)))));
 
     if (tryval > _allMSPs[i].linkingScore) {
        _allMSPs[i].linkingScore = tryval;
        _allMSPs[i].prev = j;
      }
    }
    _allMSPs[i].linkingScore += (weight * _allMSPs[i].score);
    if (_allMSPs[i].linkingScore > best_sc) {
      best = i;
      best_sc = _allMSPs[i].linkingScore;
    }
  }

  //return best;

  //
  //
  ////////////////////////////////////////




  //  In the code below, last_msp = best.
  //
  int last_msp = best;


  
  ////////////////////////////////////////
  //
  //  Make exons
  //
  //  msp2exons(msp,last_msp,s1,s2); 
  //

  int     diag_dist;
  int     diff;

  Exon_ptr  elist = 0L;

  if (last_msp < 0)
    return(elist);

  /* Note: in new_exon, the 'flag' and 'length' fields need not be computed */

  msp *mp = _allMSPs + last_msp;

  elist = new_exon(mp->pos1,
                   mp->pos2,
                   mp->pos1+mp->len-1, 
                   mp->pos2+mp->len-1,
                   -1, 
                   (mp->len*MATCH-mp->score) / (MATCH-MISMATCH),
                   0,
                   elist);

  last_msp = mp->prev;

  while (last_msp >= 0) {
    mp = _allMSPs + last_msp; 

    int   l1 = elist->frEST - elist->frGEN;
    int   l2 = mp->pos2     - mp->pos1;

    if (l1 > l2)
      diag_dist = l1 - l2;
    else
      diag_dist = l2 - l1;

    if ((diag_dist <= L) &&
        (elist->frEST - (mp->pos2 + mp->len - 1)) < MAX_INTERNAL_GAP) {
      /* merge with previous exon */
      elist->edist += diag_dist;
      elist->edist += (mp->len*MATCH-mp->score)/(MATCH-MISMATCH);
      if ((diff=mp->pos2+mp->len-elist->frEST)>0) {   /* overlap */
        int dist1, dist2;
        dist1 = get_edist(elist->frGEN,mp->pos2+mp->len-diff,
                          elist->frGEN+diff-1,mp->pos2+mp->len-1,s1,s2);
        dist2 = get_edist(mp->pos1+mp->len-diff,mp->pos2+mp->len-diff,
                          mp->pos1+mp->len-1,mp->pos2+mp->len-1,s1,s2);
        elist->edist -= max(dist1,dist2);
      } else if (diff<0) {  /* gap */
        elist->edist += (int)(0.5 * P * (-1) * diff);
      }
      elist->toGEN = max(elist->toGEN,mp->pos1+mp->len-1);
      elist->toEST = max(elist->toEST,mp->pos2+mp->len-1);
      elist->frGEN = min(elist->frGEN,mp->pos1);
      elist->frEST = min(elist->frEST,mp->pos2);
    } else {
      /* new exon */
      elist = new_exon(mp->pos1,
                       mp->pos2,
                       mp->pos1+mp->len-1,
                       mp->pos2+mp->len-1,
                       -1,
                       (mp->len*MATCH-mp->score) / (MATCH-MISMATCH),
                       0,
                       elist);
    }

    last_msp = mp->prev;
  } 



  //
  //
  ////////////////////////////////////////





  //  Fix them?  What does this do??
  //
  Exon_ptr tmp_block = elist;
  while (tmp_block != 0L) {
    tmp_block->length = tmp_block->toEST-tmp_block->frEST+1;
    tmp_block->toGEN   += offset1;
    tmp_block->frGEN += offset1;
    tmp_block->toEST   += offset2;
    tmp_block->frEST += offset2;
    tmp_block->flag   = flag;

    tmp_block = tmp_block->next_exon;
  }


  return(elist);
}




#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include "mspManager.H"
#include "sim4defines.H"
#include "exon.H"

#define  DEFAULT_L       8


mspManager::mspManager() {
  _sorted           = true;

  _allocMSPs        = 16384;
  _numMSPs          = 0;
  _allMSPs          = new msp [_allocMSPs];

  //  The following four variables are for aborting expensive
  //  polishes -- ones that have proven to be large chunks of
  //  genomic labeled as cDNA, and that have (ESTmapper) signals
  //  across entire scafflds.
  //
  _tooManyMSPs      = false;
  _cDNALength       = 0;
  _mspLimitAbsolute = 0;
  _mspLimitPercent  = 0;

  //  These need to be reset with setParameters.  The code will die
  //  during link() if they are not set.
  //
  _match            = 0;
  _matchdiff        = 0;
  _percentError     = 0.0;

  _diagMax          = 0;
  _diagExt          = 0L;

  _add_extended     = 0;
  _add_total        = 0;
}


mspManager::~mspManager() {
  delete [] _allMSPs;
  delete [] _diagExt;
#if 0
  fprintf(stderr, "mspManager:  added %d hits, extended %d.\n",
          _add_total, _add_extended);
#endif
}








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


static
int
mspManager_msp_compare(const void *A, const void *B) {
  msp const  *a = (msp const *)A;
  msp const  *b = (msp const *)B;

  if (a->pos2 < b->pos2)
    return(-1);

  if (a->pos2 > b->pos2)
    return(1);

  if (a->pos1 < b->pos1)
    return(-1);

  if (a->pos1 > b->pos1)
    return(1);

  return(0);
}



Exon_ptr
mspManager::doLinking(int    weight,
                      int    drange,
                      int    offset1,
                      int    offset2,
                      int    flag,
                      int    relinkFlag,
                      char  *s1,
                      char  *s2) {

  //  Ensure the MSP's are sorted
  //
  if (_sorted == false)
    qsort(_allMSPs, _numMSPs, sizeof(struct msp), mspManager_msp_compare);

  _sorted = true;

  //
  //  Assumes the exon list is cleared
  //

  //  If this ever occurs, you (the programmer) forgot to call
  //  mspManager::setParameters() with the correct values.  Unless the
  //  code was really hacked, this should never occur.  See
  //  Sim4::Sim4().
  //
  if ((_match == 0) &&
      (_matchdiff == 0) &&
      (_percentError == 0.0)) {
    fprintf(stderr, "sim4::link()-- ERROR; mspManager parameters not set!  This is an algorithm error.\n");
    exit(1);
  }

  //  Check if this match looks suspiciously expensive
  //
  if ((_cDNALength > 0) &&
      (_mspLimitAbsolute > 0) && (_mspLimitAbsolute < _numMSPs) &&
      (_mspLimitPercent > 0.0) && (_mspLimitPercent * _cDNALength < _numMSPs)) {
    _tooManyMSPs = true;
    return(0L);
  }

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

#ifdef SHOW_LINKING
    fprintf(stderr, "link %d\r", i);
    fflush(stderr);
#endif

    for (u32bit j = 0; j < i; ++j) {

      //  12 == default word size.  A Magic Value.

      int vL = DEFAULT_L; 
      if ((_allMSPs[i].pos2 + _allMSPs[i].len - _allMSPs[j].pos2 - _allMSPs[j].len > 2 * 12) &&
          (_allMSPs[i].pos2 - _allMSPs[j].pos2 > 2 * 12))
        vL *= 2;
                        
      diff_diag = diag - _allMSPs[j].pos1 + _allMSPs[j].pos2;

      //  Abort if the difference is too big
      //
      if ((diff_diag < -drange) ||
          ((diff_diag >  drange) && (diff_diag < MIN_INTRON)) ||
          (_allMSPs[j].pos2 + _allMSPs[j].len - 1 - f2 > vL) ||
          (_allMSPs[j].pos1 + _allMSPs[j].len - 1 - f1 > vL))
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
                   (mp->len * _match - mp->score) / _matchdiff,
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

    if ((diag_dist <= DEFAULT_L) &&
        (elist->frEST - (mp->pos2 + mp->len - 1)) < MAX_INTERNAL_GAP) {
      /* merge with previous exon */
      elist->edist += diag_dist;
      elist->edist += (mp->len * _match - mp->score) / _matchdiff;
      if ((diff=mp->pos2+mp->len-elist->frEST)>0) {   /* overlap */
        int dist1, dist2;
        dist1 = get_edist(elist->frGEN,mp->pos2+mp->len-diff,
                          elist->frGEN+diff-1,mp->pos2+mp->len-1,s1,s2);
        dist2 = get_edist(mp->pos1+mp->len-diff,mp->pos2+mp->len-diff,
                          mp->pos1+mp->len-1,mp->pos2+mp->len-1,s1,s2);
        elist->edist -= max(dist1,dist2);
      } else if (diff<0) {  /* gap */
        elist->edist += (int)(0.5 * _percentError * (-1) * diff);
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
                       (mp->len * _match - mp->score) / _matchdiff,
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
    tmp_block->length  = tmp_block->toEST-tmp_block->frEST+1;
    tmp_block->toGEN  += offset1;
    tmp_block->frGEN  += offset1;
    tmp_block->toEST  += offset2;
    tmp_block->frEST  += offset2;
    tmp_block->flag    = flag;

    tmp_block = tmp_block->next_exon;
  }


  return(elist);
}




void
mspManager::addHit_(char *s1, char *s2,
                    int   l1, int   l2,
                    int   p1, int   p2,
                    int   W,
                    int   K) {
  char *beg2 = 0L;
  char *beg1 = 0L;
  char *end1 = 0L;
  char *q = 0L;
  char *s = 0L;
  int   right_sum = 0;
  int   left_sum  = 0;
  int   sum       = 0;
  int   score     = 0;

  _add_extended++;

  //  We use diagonals directly -- original version offset the array of
  //  diagonal positions by the constant value included below.

  //  Extend to the right
  //
  left_sum = 0;
  sum      = 0;
  q        = s1 + 1 + p1;
  s        = s2 + 1 + p2;
  end1     = q;

  while ((*s) &&
         (*q) &&
         (s <= s2 + l2) &&
         (q <= s1 + l1) &&
         (sum >= left_sum - _wordExtAllow)) {

    sum += _match;
    if (*s != *q)
      sum -= _matchdiff;

    s++;
    q++;
    if (sum > left_sum) {
      left_sum = sum;
      end1 = q;
    }
  }

  //  Extend to the right
  //
  right_sum = 0;
  sum       = 0;
  q         = s1 + 1 + p1 - W;
  s         = s2 + 1 + p2 - W;
  beg1      = q;
  beg2      = s;

  while ((s > s2 + 1) &&
         (q > s1 + 1) &&
         (sum >= right_sum - _wordExtAllow)) {

    s--;
    q--;
    sum += _match;
    if (*s != *q)
      sum -= _matchdiff;

    if (sum > right_sum) {
      right_sum = sum;
      beg2 = s;
      beg1 = q;
    }
  }

  score = W + left_sum + right_sum;

  if (score >= K)
    addMSP((int)(end1 - beg1),
           (int)(beg1 - (s1 + 1)),
           (int)(beg2 - (s2 + 1)),
           score);

  _diagExt[l2 + p1 - p2 - 1] = (int)(end1 - s1 - 1 + W);
}

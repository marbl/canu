#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include "mspManager.H"
#include "sim4defines.H"
#include "exon.H"

#define  DEFAULT_L       8

//  Define this to get a detailed report on hit extension
//
//#define DEBUG_EXTENSION

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
}


mspManager::~mspManager() {
  delete [] _allMSPs;
  delete [] _diagExt;
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





  //  In the code below, last_msp = best.
  //
  int last_msp = best;
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








// The log4 arrays were computed to mimick the behaviour of the log formula
// for computing the msp threshold in exon_cores(). For genomic_log4s,
// entry i stores the value for the length of a genomic sequence
// for which the contribution to the msp threshold is i/2, i.e.:
//    1.4*log_4(3/4*len1) = i/2;  
//
// Similarly, cDNA_log4s entries store lengths of the cDNA sequence for which
// the contribution to the msp threshold is i/2, i.e.:
//    1.4*log_4(len2) = i/2;
//
// Both arrays are sorted in increasing order, and can be searched with 
// binary search.
//
#define GEN_LOG4_ENTRIES  45
#define CDNA_LOG4_ENTRIES 25

const int
genomic_log4s[GEN_LOG4_ENTRIES]= {1, 2, 3, 5, 9, 15, 26, 42, 70, 114,
                                  188, 309, 507, 832, 1365, 1365, 2240, 2240, 3675, 6029,
                                  9892, 16231, 26629, 43690, 71681,
                                  117606, 192953, 316573, 519392, 852152,
                                  1398101, 2293823, 3763409, 6174516, 10130347,
                                  16620564, 27268873, 44739242, 73402365, 120429110,
                                  197584514, 324171126, 531858072, 872603963, 1431655765 };

const int
cDNA_log4s[CDNA_LOG4_ENTRIES]= {1, 1, 2, 4, 7, 11, 19, 32, 52, 86,
                                141, 231, 380, 624, 1024, 1680, 2756, 4522, 7419, 12173,
                                19972, 32768, 53761, 88204, 144715 };


//  The original used a binary search but with so few entries brute
//  force works better.
//
int
get_msp_threshold(int len1, int len2) {
  int   i, j;

  //  Find the index of the largest value smaller than our lengths.
  //
  i = 0;
  while (i<GEN_LOG4_ENTRIES) {
    if (genomic_log4s[i] > len1)
      break;
    i++;
  }
  i--;

  j = 0;
  while (j<CDNA_LOG4_ENTRIES) {
    if (cDNA_log4s[j] > len2)
      break;
    j++;
  }
  j--;

  //
  //  XXX:  This looks suspicious!
  //

  if ((i % 2) == 0)
    return(i/2+j/2);

  if ((j % 2) == 0)
    return(i/2+j/2);

  return(i/2+j/2+1);
}





void
mspManager::setScoreThreshold(int K, int interspecies) {

  if (interspecies) {
    _minMSPScore = (int)(((int)(log(.75*(double)_GENlen)+log((double)_ESTlen))/log(4.0)) * 1.0);
  } else {
    if (K <= 0) {
      _minMSPScore = get_msp_threshold(_GENlen, _ESTlen);
        
      //  compensate for the rounding in the log formula
      if (_minMSPScore >= 0)
        _minMSPScore--;
    } else {
      _minMSPScore = K;
    }
  }
}




















void
mspManager::addHit_(char *genSeq, char *estSeq,
                    int   genLen, int   estLen,
                    int   genPos, int   estPos,
                    int   W) {
  char *genBeg = 0L;
  char *estBeg = 0L;
  char *genEnd = 0L;
  char *genTmp = 0L;
  char *estTmp = 0L;
  int   right_sum  = 0;
  int   middle_sum = 0;
  int   left_sum   = 0;
  int   sum        = 0;
  int   score      = 0;

#ifdef DEBUG_EXTENSION
  {
    char L[41], M[41], R[41];
    int  x;

    if (genPos > 20) genTmp = genSeq + genPos - 20;
    else             genTmp = genSeq;

    x=0;
    while (genTmp < genSeq + genPos)
      L[x++] = *genTmp++;
    L[x] = 0;
    x=0;
    while (genTmp < genSeq + genPos + W)
      M[x++] = *genTmp++;
    M[x] = 0;
    x=0;
    while (genTmp < genSeq + genPos + W + 20)
      R[x++] = *genTmp++;
    R[x] = 0;
    fprintf(stderr, "GEN=%8d %s:%s:%s\n", genPos, L, M, R);

    if (estPos > 20) estTmp = estSeq + estPos - 20;
    else             estTmp = estSeq;

    x=0;
    while (estTmp < estSeq + estPos)
      L[x++] = *estTmp++;
    L[x] = 0;
    x=0;
    while (estTmp < estSeq + estPos + W)
      M[x++] = *estTmp++;
    M[x] = 0;
    x=0;
    while (estTmp < estSeq + estPos + W + 20)
      R[x++] = *estTmp++;
    R[x] = 0;
    fprintf(stderr, "EST=%8d %s:%s:%s\n", estPos, L, M, R);
  }
#endif

  //  We use diagonals directly -- original version offset the array of
  //  diagonal positions by the constant value included below.

  //  Extend to the right
  //
  left_sum = 0;
  sum      = 0;
  genTmp   = genSeq + 1 + genPos;
  estTmp   = estSeq + 1 + estPos;
  genEnd   = genTmp;

  while ((*genTmp) &&
         (*estTmp) &&
         (estTmp <= estSeq + estLen) &&
         (genTmp <= genSeq + genLen) &&
         (sum >= left_sum - _wordExtAllow)) {

    sum += _match;
    if (*estTmp != *genTmp)
      sum -= _matchdiff;

    estTmp++;
    genTmp++;
    if (sum > left_sum) {
      left_sum = sum;
      genEnd = genTmp;
    }
  }

#ifdef DEBUG_EXTENSION
  //  Check the bases that the seed supposedly matched
  //
  middle_sum = 0;
  sum        = 0;
  genTmp     = genSeq + 1 + genPos - W;
  estTmp     = estSeq + 1 + estPos - W;

  for (int x=0; x<W; x++) {
    middle_sum += _match;
    if (*genTmp != *estTmp)
      middle_sum -= _matchdiff;

    estTmp++;
    genTmp++;
  }
#endif


  //  Extend to the left
  //
  right_sum = 0;
  sum       = 0;
  genTmp    = genSeq + 1 + genPos - W;
  estTmp    = estSeq + 1 + estPos - W;
  genBeg    = genTmp;
  estBeg    = estTmp;

  while ((estTmp > estSeq + 1) &&
         (genTmp > genSeq + 1) &&
         (sum >= right_sum - _wordExtAllow)) {

    estTmp--;
    genTmp--;
    sum += _match;
    if (*estTmp != *genTmp)
      sum -= _matchdiff;

    if (sum > right_sum) {
      right_sum = sum;
      estBeg = estTmp;
      genBeg = genTmp;
    }
  }

  score = W + left_sum + right_sum;

  if (score >= _minMSPScore)
    addMSP((int)(genEnd - genBeg),
           (int)(genBeg - (genSeq + 1)),
           (int)(estBeg - (estSeq + 1)),
           score);

#ifdef DEBUG_EXTENSION
  fprintf(stderr, "mspManager::addHit()-- added from EST %d to %d and GEN %d to %d (length = %d) with score %d (needed %d) l,m,r sums %d %d %d\n", 
          (int)(estBeg - (estSeq + 1)), (int)(estBeg - (estSeq + 1)) + W,
          (int)(genBeg - (genSeq + 1)), (int)(genBeg - (genSeq + 1)) + W,
          W,
          score, _minMSPScore, left_sum, middle_sum, right_sum);
#endif

  //  Remember the highest point that this diagonal has been extended
  //  to.  We use this to short circuit useless mer extensions (if
  //  we've already extended through it).
  //
  _diagExt[estLen + genPos - estPos - 1] = (int)(genEnd - genSeq - 1 + W);
}

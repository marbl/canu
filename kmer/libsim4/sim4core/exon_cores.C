#include "sim4.H"
#include <math.h>

//  diag_lev is used only in these functions.
//
//  search() and extend_hit() are used only here.
//



void
Sim4::exon_cores(char *s1, char *s2,
                 int l1, int l2,
                 int offset1, int offset2,
                 int flag,
                 int in_W, int in_K,
                 int type) {

#if 0
  fprintf(stdout, "exonCores: l1=%d l2=%d offset1=%d offset2=%d flag=%d W=%d K=%d\n",
          l1, l2, offset1, offset2, flag, in_W, in_K);
#endif

  int K;
  if (globalParams->_interspecies)
    K = (int)(((int)(log(.75*(double)len1)+log((double)len2))/log(4.0)) * 1.0);
  else
    K = getMSPthreshold(in_K, l1, l2);

  _mspManager.clear();

  exon_list = NULL;

  //  XXX: We used to limit W to be 12 maximum.  BPW noticed on
  //  28apr03 that we're still using in_W (unlimited W).
  //  
  //int W = min(in_W, l2);

  bld_table(s2,l2,in_W,type);
  search(s1,s2,l1,l2,in_W,K);

  //  Cleaning up after the bld_table() is done at the next call, or
  //  in the destructor.
  //
  hashtable = 0L;

  _mspManager.sort();

  exon_list = _mspManager.link(DEFAULT_WEIGHT, DEFAULT_DRANGE, offset1, offset2, flag, false, s1, s2);
}



void
Sim4::search(char *s1, char *s2, int l1, int l2, int in_W, int in_K) {
  struct hash_node *h;
  char             *t;
  int               ecode;
  int               i, p;
  int               lower;
  int              *allocated;
  int              *diag_lev;

  //  Too short?  Abort!
  //
  if (l1 < in_W)
    return;

  //double startTime = getTime();

  lower = -l1;

  allocated = new int [l1 + l2 + 1];
  for (i = l1 + l2 + 1; i; )
    allocated[--i] = 0;

  diag_lev = allocated - lower;

  t = s1+1;
  i = 0;

  register int   validEncoding = 1 - in_W;
  register int   pos1;

  ecode = 0L;

  //  5% win (tested on on small examples) if we use t[] instead of *t below.

  for (i=0; i < l1; i++) {
    if (encoding[t[i]] >= 0) {
      validEncoding++;

      ecode  &= mask;
      ecode <<= 2;
      ecode  |= encoding[t[i]];

      if (validEncoding > 0) {
        for (h = hashtable->table[ecode & HASH_SIZE]; h; h = h->link) {
          if (h->ecode == ecode) {
            for (p = h->pos; p >= 0; p = hashtable->nextPos[p]) {
              pos1 = (int)(t-s1) + i;

              //fprintf(stdout, "p = %8d  pos1 = %8d\n", p, pos1);
              //fprintf(stdout, "p1=%8d o1=%8d p2=%8d o2=%8d\n", pos1, 0, p, 0);

              if (diag_lev[p-pos1] <= pos1) {
                diag_lev[p-pos1] = extend_hit(pos1, p,
                                              s1,s2,l1,l2,
                                              in_W, in_K);
              }
            }
            break;
          }
        }
      }
    } else {
      validEncoding = 1 - in_W;
    }
  }

  delete [] allocated;
}




/* extend_hit - extend a word-sized hit to a longer match */
int
Sim4::extend_hit(int pos1, int pos2, 
                 const char * const s1, const char * const s2,
                 int l1, int l2, int in_W, int in_K)
{
  const char *beg2, *beg1, *end1, *q, *s;
  int right_sum, left_sum, sum, score;

  int MATCH    = globalParams->_match;
  int MISMATCH = globalParams->_mismatch;

  /* extend to the right */
  left_sum = sum = 0;
  q = s1+1+pos1;
  s = s2+1+pos2;
  end1 = q;
  while ((*s != '\0') && (*q != '\0') &&
         (s<=s2+l2) && (q<=s1+l1) && sum >= left_sum - X) {
    sum += ((*s++ == *q++) ? MATCH : MISMATCH);
    if (sum > left_sum) {
      left_sum = sum;
      end1 = q;
    }
  }

  /* extend to the left */
  right_sum = sum = 0;
  beg1 = q = (s1+1+pos1) - in_W;
  beg2 = s = (s2+1+pos2) - in_W;
  while ((s>s2+1) && (q>s1+1) && sum >= right_sum - X) {
    sum += ((*(--s) == *(--q)) ? MATCH : MISMATCH);
    if (sum > right_sum) {
      right_sum = sum;
      beg2 = s;
      beg1 = q;
    }
  }

  //  XXX: new
  score = in_W + left_sum + right_sum;

  if (score >= in_K)
    _mspManager.addMSP((int)(end1 - beg1),
                       (int)(beg1 - (s1 + 1)),
                       (int)(beg2 - (s2 + 1)),
                       score);

  return((int)(end1 - s1 - 1 + in_W));
}


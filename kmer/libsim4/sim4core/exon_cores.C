#include "sim4.H"
#include <math.h>


void
Sim4::exon_cores(char *s1,
                 char *s2,
                 int   l1,
                 int   l2,
                 int   offset1,
                 int   offset2,
                 int   flag,
                 int   W,
                 int   K,
                 int   type) {

  _mspManager.clear();
  _mspManager.clearDiagonal(l1, l2);
  _mspManager.setScoreThreshold(K, globalParams->_interspecies);

  bld_table(s2,l2,W,type);
  search(s1,s2,l1,l2,W);

  //  Cleaning up after the bld_table() is done at the next call, or
  //  in the destructor.
  //
  hashtable = 0L;

  exon_list = _mspManager.doLinking(DEFAULT_WEIGHT, DEFAULT_DRANGE, offset1, offset2, flag, false, s1, s2);
}




void
Sim4::search(char *s1, char *s2, int l1, int l2, int W) {
  struct hash_node *h;
  char             *t;
  int               ecode;
  int               i, p;

  //  Too short?  Abort!
  //
  if (l1 < W)
    return;

  t = s1+1;
  i = 0;

  int   validEncoding = 1 - W;
  int   pos1;

  ecode = 0L;

  //  5% win (tested on on small examples) if we use t[] instead of *t below.

  for (i=0; i < l1; i++) {
    pos1 = (int)(t-s1) + i;

    if (encoding[(int)t[i]] >= 0) {
      validEncoding++;

      ecode  &= mask;
      ecode <<= 2;
      ecode  |= encoding[(int)t[i]];

      if (validEncoding > 0) {
        for (h = hashtable->table[ecode & HASH_SIZE]; h; h = h->link) {
          if (h->ecode == ecode) {
            for (p = h->pos; p >= 0; p = hashtable->nextPos[p])
              _mspManager.addHit(s1, s2,
                                 l1, l2,
                                 pos1, p,
                                 W);
            break;
          }
        }
      }
    } else {
      validEncoding = 1 - W;
    }
  }
}

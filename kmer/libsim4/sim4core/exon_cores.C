#include "sim4.H"
#include <math.h>

//  exon_cores() must have seq-1 passed in.  search() offsets this.

void
Sim4::exon_cores(char *s1,
                 char *s2,
                 int   l1,
                 int   l2,
                 int   offset1,
                 int   offset2,
                 int   flag,
                 mss_t MSS,
                 int   K,
                 int   type) {

  _mspManager.clear();
  _mspManager.clearDiagonal(l1, l2);
  _mspManager.setScoreThreshold(K, globalParams->_interspecies);

//mss_t MSS = masks_shifts(seed);   LLL DELETE

  bld_table(s2,l2,MSS,type);
  search(s1,s2,l1,l2,MSS);

  //  Cleaning up after the bld_table() is done at the next call, or
  //  in the destructor.
  //
  hashtable = 0L;

  exon_list = _mspManager.doLinking(DEFAULT_WEIGHT, DEFAULT_DRANGE, offset1, offset2, flag, false, s1, s2);
}




void
Sim4::search(char *s1, char *s2, int l1, int l2, mss_t MSS) {
  struct hash_node *h;
  char             *t;
  u64bit            ecode;
  int               masked_ecode;
  int               i, p, j;

  //  Too short?  Abort!
  //
  if (l1 < MSS.seedLength)
    return;

  t = s1+1;
  i = 0;

  int   validEncoding = 1 - MSS.seedLength;
  int   pos1;

  ecode = u64bitZERO;

  //  5% win (tested on on small examples) if we use t[] instead of *t below.

  //  Scan from low to high position in the genomic sequence
  //
  if (MSS.type == CONTINUOUS_SEED) {
    for (i=0; i < l1; i++) {
      pos1 = (int)(t-s1) + i;

      if (encoding[(int)t[i]] >= 0) {
        validEncoding++;

        ecode  &= mask;
        ecode <<= 2;
        ecode  |= encoding[(int)t[i]];
        masked_ecode = (int)ecode;

        if (validEncoding > 0) {
          for (h = hashtable->table[masked_ecode & HASH_SIZE]; h; h = h->link) {
            if (h->ecode == masked_ecode) {

              //  These positions are from high to low (see table.C)
              //
              for (p = h->pos; p >= 0; p = hashtable->nextPos[p])
                _mspManager.addHit(s1, s2,
                                   l1, l2,
                                   pos1, p,
                                   MSS);
              break;
            }
          }
        }
      } else {
        validEncoding = 1 - MSS.seedLength;
      }
    }
  } else {
    /* SPACED_SEED */
    for (i=0; i < l1; i++) {
      pos1 = (int)(t-s1) + i;

      if (encoding[(int)t[i]] >= 0) {
        validEncoding++;

        ecode  &= MSS.mask;
        ecode <<= 2;
        ecode  |= encoding[(int)t[i]];

#if 0
        masked_ecode = mask_shift(ecode,MSS);
#else
        // 40% cheaper for cross-species, 53% cheaper for same species
        for (j=masked_ecode=0; j<MSS.masknum; j++)
          masked_ecode += (ecode & MSS.masks[j]) >> MSS.shifts[j];
#endif

        if (validEncoding > 0) {
          for (h = hashtable->table[masked_ecode & HASH_SIZE]; h; h = h->link) {
            if (h->ecode == masked_ecode) {

              //  These positions are from high to low (see table.C)
              //
              for (p = h->pos; p >= 0; p = hashtable->nextPos[p])
                _mspManager.addHit(s1, s2,
                                   l1, l2,
                                   pos1, p,
                                   MSS);
              break;
            }
          }
        }
      } else {
        validEncoding = 1 - MSS.seedLength;
      }
    }
  }
}

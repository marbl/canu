#include "sim4.H"
#include <math.h>

//  diag_lev is used only in these functions.
//
//  search() and extend_hit() are used only here.
//


//  Define this to get some profiling of building, searching, sorting and linking
//
//#define PROFILE_EXON_CORE




void
Sim4::exon_cores(char *s1,
                 char *s2,
                 int   l1,
                 int   l2,
                 int   offset1,
                 int   offset2,
                 int   flag,
                 int   in_W,
                 int   in_K,
                 int   type) {

  int K;
  if (globalParams->_interspecies)
    K = (int)(((int)(log(.75*(double)_genLen)+log((double)_estLen))/log(4.0)) * 1.0);
  else
    K = getMSPthreshold(in_K, l1, l2);

  _mspManager.clear();
  _mspManager.clearDiagonal(l1, l2);

  exon_list = NULL;

#ifdef PROFILE_EXON_CORE
  double startTime = getTime();
#endif

  bld_table(s2,l2,in_W,type);
  search(s1,s2,l1,l2,in_W,K);

#ifdef PROFILE_EXON_CORE
  fprintf(stderr, "build+search took %f seconds.\n", getTime() - startTime);
#endif

  //  Cleaning up after the bld_table() is done at the next call, or
  //  in the destructor.
  //
  hashtable = 0L;

#ifdef PROFILE_EXON_CORE
  startTime = getTime();
  fprintf(stderr, "Linking "u32bitFMT" MSPs\n", _mspManager.numberOfMSPs());
#endif

  exon_list = _mspManager.doLinking(DEFAULT_WEIGHT, DEFAULT_DRANGE, offset1, offset2, flag, false, s1, s2);

#ifdef PROFILE_EXON_CORE
  fprintf(stderr, "linking took %f seconds.\n", getTime() - startTime);
#endif
}




void
Sim4::search(char *s1, char *s2, int l1, int l2, int in_W, int in_K) {
  struct hash_node *h;
  char             *t;
  int               ecode;
  int               i, p;

  //  Too short?  Abort!
  //
  if (l1 < in_W)
    return;

  t = s1+1;
  i = 0;

  int   validEncoding = 1 - in_W;
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
                                 in_W, in_K);
            break;
          }
        }
      }
    } else {
      validEncoding = 1 - in_W;
    }
  }
}

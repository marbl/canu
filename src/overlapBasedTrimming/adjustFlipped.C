static const char *rcsid = "$Id:  $";

#include "adjustOverlaps.H"

//  Adjust the overlap for any trimming done already.  This works by computing the fraction of the
//  overlap trimmed for each read and each end, picking the largest fraction for each end, and
//  applying that fraction to the other read.
//
//  It expects only flipped overlaps.  The output coordinates are for a REVERSE COMPLEMENTED
//  b read.  if you care which end is the actual 5' or 3' end, look at flipped().

bool
adjustFlipped(clearRangeFile  *iniClr,
              gkStore         *gkp,
              ovsOverlap      *ovl,
              uint32 &aovlbgn,  uint32 &aovlend,  uint32 &bovlbgn,  uint32 &bovlend,
              uint32 &aclrbgn,  uint32 &aclrend,  uint32 &bclrbgn,  uint32 &bclrend) {

  assert(ovl->flipped() == true);

  uint32  bLen = gkp->gkStore_getRead(ovl->b_iid)->gkRead_sequenceLength();

  aovlbgn =        ovl->a_bgn(gkp);
  bovlbgn = bLen - ovl->b_bgn(gkp);  //  bgn(), because this is the higher coord
  aovlend =        ovl->a_end(gkp);
  bovlend = bLen - ovl->b_end(gkp);

  aclrbgn =        iniClr->bgn(ovl->a_iid);
  bclrbgn = bLen - iniClr->end(ovl->b_iid);  //  end(), because this is the higher coord
  aclrend =        iniClr->end(ovl->a_iid);
  bclrend = bLen - iniClr->bgn(ovl->b_iid);

  assert(aovlbgn < aovlend);
  assert(bovlbgn < bovlend);

  if ((aclrend <= aovlbgn) || (aovlend <= aclrbgn) ||
      (bclrend <= bovlbgn) || (bovlend <= bclrbgn)) {
    //  Overlap doesn't intersect clear range, fail.
#if 0
    fprintf(stderr, "Discard  FLIP overlap from %u,%u-%u,%u based on clear ranges %u,%u and %u,%u\n",
            aovlbgn, aovlend,
            bovlbgn, bovlend,
            aclrbgn, aclrend,
            bclrbgn, bclrend);
#endif
    return(false);
  }


  uint32  alen = aovlend - aovlbgn;
  uint32  blen = bovlend - bovlbgn;

  double  afracbgn = (double)((aclrbgn < aovlbgn) ? (0) : (aclrbgn - aovlbgn)) / alen;
  double  bfracbgn = (double)((bclrbgn < bovlbgn) ? (0) : (bclrbgn - bovlbgn)) / blen;
  double  afracend = (double)((aclrend > aovlend) ? (0) : (aovlend - aclrend)) / alen;
  double  bfracend = (double)((bclrend > bovlend) ? (0) : (bovlend - bclrend)) / blen;

  //fprintf(stderr, "frac a %.20f %.20f b %.20f %.20f\n", afracbgn, afracend, bfracbgn, bfracend);
  //fprintf(stderr, "frac a %.20f %.20f b %.20f %.20f\n", afracbgn * alen, afracend * alen, bfracbgn * blen, bfracend * blen);

  double  maxbgn = max(afracbgn, bfracbgn);
  double  maxend = max(afracend, bfracend);

  //fprintf(stderr, "frac a %.20f %.20f b %.20f %.20f\n", maxbgn * alen, maxend * alen, maxbgn * blen, maxend * blen);

  assert(maxbgn < 1.0);
  assert(maxend < 1.0);

  uint32  aadjbgn = (uint32)round(maxbgn * alen);
  uint32  badjbgn = (uint32)round(maxbgn * blen);
  uint32  aadjend = (uint32)round(maxend * alen);
  uint32  badjend = (uint32)round(maxend * blen);

  //fprintf(stderr, "frac a %u %u b %u %u alen %u blen %u\n", aadjbgn, aadjend, badjbgn, badjend, alen, blen);

#if 0
  fprintf(stderr, "Adjusted FLIP overlap from %u,%u-%u,%u (adjust %u,%u,%u,%u) to %u,%u-%u,%u  based on clear ranges %u,%u and %u,%u  maxbgn=%f maxend=%f\n",
          aovlbgn, aovlend,
          bovlbgn, bovlend,
          aadjbgn, aadjend, badjbgn, badjend,
          aovlbgn + aadjbgn, aovlend - aadjend,
          bovlbgn + badjbgn, bovlend - badjend,
          aclrbgn, aclrend,
          bclrbgn, bclrend,
          maxbgn,
          maxend);
#endif

  aovlbgn += aadjbgn;
  bovlbgn += badjbgn;
  aovlend -= aadjend;
  bovlend -= badjend;

  assert(aclrbgn <= aovlbgn);
  assert(bclrbgn <= bovlbgn);
  assert(aovlend <= aclrend);
  assert(bovlend <= bclrend);

  return(true);
}
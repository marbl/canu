
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-APR-15 to 2015-JUN-16
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "splitReads.H"

//  Adjust the overlap for trimming already done.  Filter out useless overlaps.  Save the
//  overlap into a useful format.

void
workUnit::addAndFilterOverlaps(gkStore *gkp,
                               clearRangeFile *finClr,
                               double errorRate,
                               ovOverlap *ovl, uint32 ovlLen) {

  if (adjMax < ovlLen) {
    delete [] adj;

    adjMax = ovlLen;
    adj    = new adjOverlap [adjMax];
  }

  adjLen = 0;

  for (uint32 oo=0; oo<ovlLen; oo++) {
    ovOverlap *o = ovl + oo;
    adjOverlap *a = adj + adjLen;

    int32 idA    = o->a_iid;
    int32 idB    = o->b_iid;

    if (finClr->isDeleted(idA) ||
        finClr->isDeleted(idB))
      continue;

    //  Returns the finClr clear range of the two reads (in *clr*), and the overlap adjusted to
    //  that clear range (in *ovl*).  The B read coordinates are reverse-complemented if the
    //  overlap is reversed (so bovlbgn < bovlend always).


    a->a_iid   = o->a_iid;
    a->b_iid   = o->b_iid;
    a->flipped = o->flipped();

    bool  valid = false;

    if (o->flipped() == false)
      valid = adjustNormal(finClr, gkp, o,
                           a->aovlbgn, a->aovlend, a->bovlbgn, a->bovlend,
                           a->aclrbgn, a->aclrend, a->bclrbgn, a->bclrend);
    else
      valid = adjustFlipped(finClr, gkp, o,
                            a->aovlbgn, a->aovlend, a->bovlbgn, a->bovlend,
                            a->aclrbgn, a->aclrend, a->bclrbgn, a->bclrend);

    if (valid == false)
      //  adjust() says the overlap doesn't intersect the clear range, so nothing here.
      continue;

    assert(a->aclrbgn <= a->aovlbgn);  //  Also checked in adjust()...
    assert(a->bclrbgn <= a->bovlbgn);

    assert(a->aovlend <= a->aclrend);
    assert(a->bovlend <= a->bclrend);


    //  Reset evidence hangs that are close to zero to be zero.  This shifts the overlap end point
    //  to remove as much of the useless hang as possible.

    if (a->bovlbgn - a->bclrbgn < 15) {
      uint32 limit  = a->aovlbgn - a->aclrbgn;
      uint32 adjust = a->bovlbgn - a->bclrbgn;

      if (adjust > limit)
        adjust = limit;

      if ((a->aovlbgn < adjust) || (a->bovlbgn < adjust))
        fprintf(stderr, "ovl %u %u-%u %u %u-%u -> clr %u-%u %u-%u adj %u-%u %u-%u\n",
                o->a_iid, o->a_bgn(), o->a_end(), o->b_iid, o->b_bgn(), o->b_end(),
                a->aclrbgn, a->aclrend, a->bclrbgn, a->bclrend,
                a->aovlbgn, a->aovlend, a->bovlbgn, a->bovlend);
      assert(a->aovlbgn >= adjust);
      assert(a->bovlbgn >= adjust);

      a->aovlbgn -= adjust;
      a->bovlbgn -= adjust;

      assert(a->aclrbgn <= a->aovlbgn);
      assert(a->bclrbgn <= a->bovlbgn);
    }


    if (a->bclrend - a->bovlend < 15) {
      uint32 limit  = a->aclrend - a->aovlend;
      uint32 adjust = a->bclrend - a->bovlend;

      if (adjust > limit)
        adjust = limit;

      a->aovlend += adjust;
      a->bovlend += adjust;

      assert(a->aovlend <= a->aclrend);
      assert(a->bovlend <= a->bclrend);
    }


    //  Filter out garbage overlaps
    //
    //  The first version used hard and fast cutoffs of (>=35bp and <= 0.02 error) or (>= 70bp).
    //  These were not fair to short reads.
    //
    //  The second version used compilacted rules (see cvs), but based the decision off of the
    //  length of the A-read.  Fair, but only for assemblies of similarily sized reads.  Totally
    //  unfair for Sanger-vs-Illumin overlaps.
    //
    //  The third version sets a minimum length and identity based on the shorter fragment. These
    //  mirror closely what the second version was doing.  It was extended to allow even shorter
    //  lengths if either read is aligned to an end.
    //
    //  The fourth version, for long noisy reads, accepts all overlaps.

    if (0)
      continue;


    //  Save the overlap.

    adjLen++;
  }
}


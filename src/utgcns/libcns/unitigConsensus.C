
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static char *rcsid = "$Id$";

#include "unitigConsensus.H"

#include <set>

using namespace std;


//  If defined, failure to align a read causes immediate crash.
#undef FAILURE_IS_FATAL

//  If defined, skip ALL contained reads.  This will cause problems in scaffolding.
#undef SKIP_CONTAINS



void
unitigConsensus::reportStartingWork(void) {
  fprintf(stderr, "unitigConsensus()-- processing fragment mid %d pos %d,%d ahcnor %d,%d,%d\n",
          utgpos[tiid].ident(),
          utgpos[tiid].bgn(),
          utgpos[tiid].end(),
          utgpos[tiid].anchor(),
          utgpos[tiid].aHang(),
          utgpos[tiid].bHang());

  if (showPlacementBefore())
    for (int32 x=0; x<=tiid; x++)
      fprintf(stderr, "unitigConsensus()-- mid %10d  utgpos %7d,%7d  cnspos %7d,%7d  anchor %10d,%6d,%6d\n",
              utgpos[x].ident(),
              utgpos[x].bgn(), utgpos[x].end(),
              cnspos[x].bgn(), cnspos[x].end(),
              utgpos[x].anchor(), utgpos[x].aHang(), utgpos[x].bHang());
}


void
unitigConsensus::reportFailure(uint32 *failed) {
  if (failed != NULL)
    failed[tiid] = true;

  fprintf(stderr, "unitigConsensus()-- failed to align fragment %d in unitig %d.\n",
          utgpos[tiid].ident(), tig->tigID());
}


void
unitigConsensus::reportSuccess(uint32 *failed) {
  if (failed != NULL)
    failed[tiid] = false;

  //fprintf(stderr, "unitigConsensus()-- fragment %d aligned in unitig %d.\n",
  //        utgpos[tiid].ident(), tig->tigID());
}


int
unitigConsensus::initialize(gkStore *gkpStore, uint32 *failed) {

  int32 num_columns = 0;
  int32 num_bases   = 0;

  if (numfrags == 0)
    return(false);

  utgpos = new tgPosition [numfrags];
  cnspos = new tgPosition [numfrags];

  memcpy(utgpos, tig->getChild(0), sizeof(tgPosition) * numfrags);
  memcpy(cnspos, tig->getChild(0), sizeof(tgPosition) * numfrags);

  for (int32 i=0; i<numfrags; i++) {
    if (failed != NULL)
      failed[i]  = true;

    int32 flen   = (utgpos[i].bgn() < utgpos[i].end()) ? (utgpos[i].end() < utgpos[i].bgn()) : (utgpos[i].bgn() - utgpos[i].end());
    num_bases   += (int32)ceil(flen + 2 * AS_CNS_ERROR_RATE * flen);

    num_columns  = (utgpos[i].bgn() > num_columns) ? utgpos[i].bgn() : num_columns;
    num_columns  = (utgpos[i].end() > num_columns) ? utgpos[i].end() : num_columns;
  }

  
  //ResetStores(num_bases, numfrags, num_columns);

  //  Magic initialization (in ResetStores()) prevents us calling CreateMANode() until now.

  trace      = new int32 [2 * AS_MAX_READLEN];
  traceBgn   = 0;

  //manode   = CreateMANode(tig->tigID());
  abacus     = new abAbacus(gkpStore);

  frankensteinLen = 0;
  frankensteinMax = MAX(1024 * 1024, 2 * num_columns);
  frankenstein    = new char     [frankensteinMax];
  frankensteinBof = new abBeadID [frankensteinMax];

  {
    set<uint32>  dupFrag;

    for (uint32 i=0; i<numfrags; i++) {
      if (utgpos[i].isRead() == false) {
        fprintf(stderr, "unitigConsensus()-- Unitig %d FAILED.  Child %d is not a read.\n",
                tig->tigID(), utgpos[i].ident());
        return(false);
      }

      if (dupFrag.find(utgpos[i].ident()) != dupFrag.end()) {
        fprintf(stderr, "unitigConsensus()-- Unitig %d FAILED.  Child %d is a duplicate.\n",
                tig->tigID(), utgpos[i].ident());
        return(false);
      }

      dupFrag.insert(utgpos[i].ident());

      // This guy allocates and initializes the beads for each fragment.  Beads are not fully inserted
      // in the abacus here.

      abSeqID fid = abacus->addRead(gkpStore, utgpos[i].ident(), utgpos[i].isReverse());

      //utgpos[fid].bgn = complement ? utgpos[i].end() : utgpos[i].bgn();
      //utgpos[fid].end = complement ? utgpos[i].bgn() : utgpos[i].end();

      cnspos[fid.get()].bgn()  = 0;
      cnspos[fid.get()].end()  = 0;

      //  If this is violated, then the implicit map from utgpos[] and cnspos[] to the unitig child
      //  list is incorrect.

      assert(fid.get() == i);

      //if (VERBOSE_MULTIALIGN_OUTPUT)
      //  fprintf(stderr,"unitigConsensus()-- Added fragment mid %d pos %d,%d in unitig %d to store at local id %d.\n",
      //          utgpos[i].ident(), utgpos[i].bgn(), utgpos[i].end(), tig->tigID(), fid);
    }
  }

  //SeedMAWithFragment(manode->lid, GetFragment(fragmentStore,0)->lid, opp);
  multialign = abacus->addMultiAlign(abSeqID(0));

  if (failed)
    failed[0] = false;

  //  Save columns
  {
    abBeadID  bidx = abacus->getSequence(0)->firstBead();
    abBead   *bead = abacus->getBead(bidx);

    while (bead) {
      frankenstein   [frankensteinLen] = abacus->getBase(bead->baseIdx());
      frankensteinBof[frankensteinLen] = bead->ident();

      frankensteinLen++;

      bead = (bead->nextID().isInvalid()) ? NULL : abacus->getBead(bead->nextID());
    }

    frankenstein   [frankensteinLen] = 0;
    frankensteinBof[frankensteinLen] = abBeadID();

    cnspos[0].bgn() = 0;
    cnspos[0].end() = frankensteinLen;
  }

  return(true);
}



int
unitigConsensus::computePositionFromAnchor(void) {

  assert(piid == -1);

  uint32 anchor = utgpos[tiid].anchor();

  if (anchor == 0)
    //  No anchor?!  Damn.
    goto computePositionFromAnchorFail;

  for (piid = tiid-1; piid >= 0; piid--) {
    abSequence *aseq = abacus->getSequence(piid);

    if (anchor != aseq->gkpIdent())
      //  Not the anchor.
      continue;

    if ((cnspos[piid].bgn() == 0) &&
        (cnspos[piid].end() == 0))
      //  Is the anchor, but that isn't placed.
      goto computePositionFromAnchorFail;

    if ((utgpos[piid].end() < utgpos[tiid].bgn()) ||
        (utgpos[tiid].end() < utgpos[piid].bgn())) {
      //  Is the anchor, and anchor is placed, but the anchor doesn't agree with the placement.
      if (showPlacement())
        fprintf(stderr, "computePositionFromAnchor()-- anchor %d at utg %d,%d doesn't agree with my utg %d,%d.  FAIL\n",
                anchor,
                utgpos[piid].bgn(), utgpos[piid].end(),
                utgpos[tiid].bgn(), utgpos[tiid].end());
      goto computePositionFromAnchorFail;
    }

    //  Scale the hangs by the change in the anchor size between bogart and consensus.

    double   anchorScale = (double)(cnspos[piid].end() - cnspos[piid].bgn()) / (double)(utgpos[piid].end() - utgpos[piid].bgn());

    if (showPlacement())
      fprintf(stderr, "computePositionFromAnchor()--  frag %u in anchor %u -- hangs %d,%d -- scale %f -- final hangs %.0f,%.0f\n",
              utgpos[tiid].ident(),
              utgpos[piid].ident(),
              utgpos[tiid].aHang(),
              utgpos[tiid].bHang(),
              anchorScale,
              utgpos[tiid].aHang() * anchorScale,
              utgpos[tiid].bHang() * anchorScale);

    cnspos[tiid].bgn() = cnspos[piid].bgn() + utgpos[tiid].aHang() * anchorScale;
    cnspos[tiid].end() = cnspos[piid].end() + utgpos[tiid].bHang() * anchorScale;

    //  Hmmm, but if we shrank the read too much, add back in some of the length.  We want to end up
    //  with the read scaled by anchorScale, and centered on the hangs.

    int32   fragmentLength = utgpos[tiid].end() - utgpos[tiid].bgn();

    if ((cnspos[tiid].bgn() >= cnspos[tiid].end()) ||
        (cnspos[tiid].end() - cnspos[tiid].bgn() < 0.75 * fragmentLength)) {
      int32  center = (cnspos[tiid].bgn() + cnspos[tiid].end()) / 2;

      if (showPlacement()) {
        fprintf(stderr, "computePositionFromAnchor()--  frag %u in anchor %u -- too short.  reposition around center %d with adjusted length %.0f\n",
                utgpos[tiid].ident(),
                utgpos[piid].ident(),
                center, fragmentLength * anchorScale);
      }

      cnspos[tiid].bgn() = center - fragmentLength * anchorScale / 2;
      cnspos[tiid].end() = center + fragmentLength * anchorScale / 2;

      //  We seem immune to having a negative position.  We only use this to pull out a region from
      //  the partial consensus to align to.
      //
      //if (cnspos[tiid].bgn() < 0) {
      //  cnspos[tiid].bgn() = 0;
      //  cnspos[tiid].end() = fragmentLength * anchorScale;
      //}
    }

    assert(cnspos[tiid].bgn() < cnspos[tiid].end());

    if (showPlacement())
      fprintf(stderr, "computePositionFromAnchor()-- anchor %d at %d,%d --> beg,end %d,%d (fLen %d)\n",
              anchor,
              cnspos[piid].bgn(), cnspos[piid].end(),
              cnspos[tiid].bgn(), cnspos[tiid].end(),
              frankensteinLen);
    return(true);
  }

 computePositionFromAnchorFail:
  cnspos[tiid].bgn() = 0;
  cnspos[tiid].end() = 0;

  piid = -1;

  return(false);
}



int
unitigConsensus::computePositionFromLayout(void) {
  int32   thickestLen = 0;

  assert(piid == -1);

  //  Find the thickest qiid overlap to any cnspos fragment
  for (int32 qiid = tiid-1; qiid >= 0; qiid--) {
    if ((utgpos[tiid].bgn() < utgpos[qiid].end()) &&
        (utgpos[tiid].end() > utgpos[qiid].bgn()) &&
        ((cnspos[qiid].bgn() != 0) ||
         (cnspos[qiid].end() != 0))) {
      cnspos[tiid].bgn() = cnspos[qiid].bgn() + utgpos[tiid].bgn() - utgpos[qiid].bgn();
      cnspos[tiid].end() = cnspos[qiid].end() + utgpos[tiid].end() - utgpos[qiid].end();

      //  This assert triggers.  It results in 'ooo' below being negative, and we
      //  discard this overlap anyway.
      //
      //assert(cnspos[tiid].bgn() < cnspos[tiid].end());

      int32 ooo = MIN(cnspos[tiid].end(), frankensteinLen) - cnspos[tiid].bgn();

#if 1
      if (showPlacement())
        fprintf(stderr, "computePositionFromLayout()-- layout %d at utg %d,%d cns %d,%d --> utg %d,%d cns %d,%d -- overlap %d\n",
                utgpos[qiid].ident(),
                utgpos[qiid].bgn(), utgpos[qiid].end(), cnspos[qiid].bgn(), cnspos[qiid].end(),
                utgpos[tiid].bgn(), utgpos[tiid].end(), cnspos[tiid].bgn(), cnspos[tiid].end(),
                ooo);
#endif

      //  Occasionally we see an overlap in the original placement (utgpos overlap) by after
      //  adjusting our fragment to the frankenstein position, we no longer have an overlap.  This
      //  seems to be caused by a bad original placement.
      //
      //  Example:
      //  utgpos[a] = 13480,14239    cnspos[a] = 13622,14279
      //  utgpos[b] = 14180,15062
      //
      //  Our placement is 200bp different at the start, but close at the end.  When we compute the
      //  new start placement, it starts after the end of the A read -- the utgpos say the B read
      //  starts 700bp after the A read, which is position 13622 + 700 = 14322....50bp after A ends.

      if ((cnspos[tiid].bgn() < frankensteinLen) &&
          (thickestLen < ooo)) {
        thickestLen = ooo;

        assert(cnspos[tiid].bgn() < cnspos[tiid].end());  //  But we'll still assert cnspos is ordered correctly.

        int32 ovl   = ooo;
        int32 ahang = cnspos[tiid].bgn();
        int32 bhang = cnspos[tiid].end() - frankensteinLen;

        piid  = qiid;
      }
    }
  }

  //  If we have a VALID thickest placement, use that (recompute the placement that is likely
  //  overwritten -- ahang, bhang and piid are still correct).

  if (thickestLen >= AS_OVERLAP_MIN_LEN) {
    assert(piid != -1);

    cnspos[tiid].bgn() = cnspos[piid].bgn() + utgpos[tiid].bgn() - utgpos[piid].bgn();
    cnspos[tiid].end() = cnspos[piid].end() + utgpos[tiid].end() - utgpos[piid].end();

    assert(cnspos[tiid].bgn() < cnspos[tiid].end());

    if (showPlacement())
      fprintf(stderr, "computePositionFromLayout()-- layout %d at %d,%d --> beg,end %d,%d (fLen %d)\n",
              utgpos[piid].ident(),
              cnspos[piid].bgn(), cnspos[piid].end(),
              cnspos[tiid].bgn(), cnspos[tiid].end(),
              frankensteinLen);

    return(true);
  }

  cnspos[tiid].bgn() = 0;
  cnspos[tiid].end() = 0;

  piid = -1;

  return(false);
}



//  Occasionally we get a fragment that just refuses to go in the correct spot.  Search for the
//  correct placement in all of frankenstein, update ahang,bhang and retry.
//
//  We don't expect to have big negative ahangs, and so we don't allow them.  To unlimit this, use
//  "-fragmentLen" instead of the arbitrary cutoff below.
int
unitigConsensus::computePositionFromAlignment(void) {

  assert(piid == -1);

  ALNoverlap  *O           = NULL;
  double       thresh      = 1e-3;
  int32        minlen      = AS_OVERLAP_MIN_LEN;
  int32        ahanglimit  = -10;

  abSequence  *seq         = abacus->getSequence(tiid);
  char        *fragment    = abacus->getBases(seq);
  int32        fragmentLen = seq->length();

  O = DP_Compare(frankenstein,
                 fragment,
                 ahanglimit, frankensteinLen,  //  ahang bounds
                 frankensteinLen, fragmentLen,   //  length of fragments
                 0,
                 AS_CNS_ERROR_RATE, thresh, minlen,
                 AS_FIND_ALIGN);

  if (O == NULL)
    O = Local_Overlap_AS_forCNS(frankenstein,
                                fragment,
                                ahanglimit, frankensteinLen,  //  ahang bounds
                                frankensteinLen, fragmentLen,   //  length of fragments
                                0,
                                AS_CNS_ERROR_RATE, thresh, minlen,
                                AS_FIND_ALIGN);

  if (O == NULL) {
    cnspos[tiid].bgn() = 0;
    cnspos[tiid].end() = 0;

    piid = -1;

    if (showAlgorithm())
      fprintf(stderr, "computePositionFromAlignment()-- Returns fail.\n");
    return(false);
  }

  //  From the overlap and existing placements, find the thickest overlap, to set the piid and
  //  hangs, then reset the original placement based on that anchors original placement.
  //
  //  To work with fixFailures(), we need to scan the entire fragment list.  This isn't so
  //  bad, really, since before we were scanning (on average) half of it.
  //
  cnspos[tiid].bgn() = O->begpos;
  cnspos[tiid].end() = O->endpos + frankensteinLen;
  //fprintf(stderr, "cnspos[%3d] mid %d %d,%d\n", tiid, utgpos[tiid].ident(), cnspos[tiid].bgn(), cnspos[tiid].end());

  assert(cnspos[tiid].bgn() < cnspos[tiid].end());

  int32   thickestLen = 0;

  for (int32 qiid = numfrags-1; qiid >= 0; qiid--) {
    if ((tiid != qiid) &&
        (cnspos[tiid].bgn() < cnspos[qiid].end()) &&
        (cnspos[tiid].end() > cnspos[qiid].bgn())) {
      int32 ooo = (MIN(cnspos[tiid].end(), cnspos[qiid].end()) -
                   MAX(cnspos[tiid].bgn(), cnspos[qiid].bgn()));

      if (thickestLen < ooo) {
        thickestLen = ooo;

        int32 ovl   = ooo;
        int32 ahang = cnspos[tiid].bgn();
        int32 bhang = cnspos[tiid].end() - frankensteinLen;

        piid  = qiid;
      }
    }
  }

  if (thickestLen > 0) {
    assert(piid != -1);

    if (showPlacement())
      fprintf(stderr, "computePositionFromAlignment()-- layout %d at %d,%d --> beg,end %d,%d (fLen %d)\n",
              utgpos[piid].ident(),
              cnspos[piid].bgn(), cnspos[piid].end(),
              cnspos[tiid].bgn(), cnspos[tiid].end(),
              frankensteinLen);

    return(true);
  }

  cnspos[tiid].bgn() = 0;
  cnspos[tiid].end() = 0;

  piid = -1;

  return(false);
}


void
unitigConsensus::rebuild(bool recomputeFullConsensus) {

  //  Run abacus to rebuild an intermediate consensus sequence.  VERY expensive.
  //
  if (recomputeFullConsensus == true) {
    abacus->refreshMultiAlign(multialign);

    abacus->getMultiAlign(multialign)->refine(abacus, abAbacus_Smooth);
    abacus->getMultiAlign(multialign)->mergeRefine(abacus);

    abacus->getMultiAlign(multialign)->refine(abacus, abAbacus_Poly_X);
    abacus->getMultiAlign(multialign)->mergeRefine(abacus);

    abacus->getMultiAlign(multialign)->refine(abacus, abAbacus_Indel);
    abacus->getMultiAlign(multialign)->mergeRefine(abacus);
  }

  //  For each column, vote for the consensus base to use.  Ideally, if we just computed the full
  //  consensus, we'd use that and just replace gaps with N.

  abColID cid   = abacus->getMultiAlign(multialign)->firstColumn();
  int32   index = 0;
    
  abacus->getMultiAlign(multialign)->columns().clear();

  frankensteinLen = 0;

#warning why are we rebuilding columnList here?

  while (cid.isValid()) {
    abColumn *column = abacus->getColumn(cid);

    int32   nA = column->GetColumnBaseCount('A');
    int32   nC = column->GetColumnBaseCount('C');
    int32   nG = column->GetColumnBaseCount('G');
    int32   nT = column->GetColumnBaseCount('T');
    int32   nN = column->GetColumnBaseCount('N');
    int32   n_ = column->GetColumnBaseCount('-');
    int32   nn = 0;

    abBead *bead = abacus->getBead(column->callID());
    char    call = 'N';

    if (nA > nn) { nn = nA;  call = 'A'; }
    if (nC > nn) { nn = nC;  call = 'C'; }
    if (nG > nn) { nn = nG;  call = 'G'; }
    if (nT > nn) { nn = nT;  call = 'T'; }
    //if (nN > nn) { nn = nN;  call = 'N'; }
    //if (n_ > nn) { nn = n_;  call = 'N'; }

    //  Call should have been a gap, but we'll instead pick the most prevalant base, but lowercase
    //  it.  This is used by the dynamic programming alignment.
    if (n_ > nn)
      call = tolower(call);

    assert(call != '-');

    abacus->setBase(bead->baseIdx(), call);

    while (frankensteinLen >= frankensteinMax) {
      resizeArrayPair(frankenstein, frankensteinBof, frankensteinLen, frankensteinMax, frankensteinMax * 2);
      //frankensteinMax *= 2;
      //frankenstein     = (char    *)safe_realloc(frankenstein,    sizeof(char)    * frankensteinMax);
      //frankensteinBof  = (beadIdx *)safe_realloc(frankensteinBof, sizeof(beadIdx) * frankensteinMax);
    }
    assert(frankensteinLen < frankensteinMax);

    frankenstein   [frankensteinLen] = call;
    frankensteinBof[frankensteinLen] = bead->ident();
    frankensteinLen++;

    //  This is extracted from RefreshMANode()
    column->position() = index++;

    abacus->getMultiAlign(multialign)->columns().push_back(cid);

    cid = column->nextID();
  }

  frankenstein   [frankensteinLen] = 0;
  frankensteinBof[frankensteinLen] = abBeadID();

  //  Update the position of each fragment in the consensus sequence.

  for (int32 i=0; i<=tiid; i++) {
    if ((cnspos[i].bgn() == 0) &&
        (cnspos[i].end() == 0))
      //  Uh oh, not placed originally.
      continue;

    abSequence *seq  = abacus->getSequence(i);

    cnspos[i].bgn() = abacus->getColumn(seq->firstBead())->position();
    cnspos[i].end() = abacus->getColumn(seq->lastBead())->position() + 1;

    assert(cnspos[i].bgn() >= 0);
    assert(cnspos[i].end() > cnspos[i].bgn());
  }

  //  Finally, update the anchor/hang of the fragment we just placed.

  if (piid >= 0) {
    utgpos[tiid].anchor()  = utgpos[piid].ident();
    utgpos[tiid].aHang()   = cnspos[tiid].bgn() - cnspos[piid].bgn();
    utgpos[tiid].bHang()   = cnspos[tiid].end() - cnspos[piid].end();
  }

  piid = -1;

#warning not PrintAlignment()
  //if (showAlignments())
  //  PrintAlignment(stderr, manode->lid, 0, -1);
}


int
unitigConsensus::alignFragment(void) {
  double        origErate    = AS_CNS_ERROR_RATE;
  int32         bgnExtra     = 0;
  int32         endExtra     = 0;

  assert((cnspos[tiid].bgn() != 0) || (cnspos[tiid].end() != 0));
  assert(piid != -1);

  assert(cnspos[tiid].bgn() < cnspos[tiid].end());

  //  Compute how much extra consensus sequence to align to.  If we allow too much, we risk placing
  //  the fragment in an incorrect location -- worst case is that we find a longer higher scoring
  //  alignment, but lower identity that is rejected.  If we allow too little, we leave unaligned
  //  bits on the end of the fragment -- worst case, we have 1 extra base which might add an
  //  unnecessary gap.
  //
  //  The default used to be 100bp, which made sense for Sanger reads (maybe) but makes no sense for
  //  short reads of length say 100bp.
  //
  //  We are placing this read relative to some other read:
  //
  //      anchoring fragment   -----------------------------------------    ediff
  //      new fragment           bdiff   ------------------------------------------------
  //
  //  So we should allow AS_CNS_ERROR_RATE indel in those relative positionings.
  //

  bgnExtra = (int32)ceil(AS_CNS_ERROR_RATE * (cnspos[tiid].bgn() - cnspos[piid].bgn()));
  endExtra = (int32)ceil(AS_CNS_ERROR_RATE * (cnspos[tiid].end() - cnspos[piid].end()));

  if (bgnExtra < 0)  bgnExtra = -bgnExtra;
  if (endExtra < 0)  endExtra = -endExtra;

  if (bgnExtra < 10) bgnExtra = 10;
  if (endExtra < 10) endExtra = 10;

  //  And compute how much extra fragment sequence to align to.  We want to trim off some
  //  of the bhang sequence to prevent false alignments.  We trim off all of the bhang sequence,
  //  except for 6% of the aligned length.
  //
  //  Actual example:  the real alignment should be this:
  //
  //  CGGCAGCCACCCCATCCGGGAGGGAGATGGGGGGGTCAGCCCCCCGCCCGGCCAGCCG
  //            CCCATCCGGGAGGGAGGTGGGGGGGTCAGCCCCCCGCCCCGCCAGCCGCTCCGTCCGGGAGGGAGGTGGGGGGGTCAGCCCCCCGCCCGGCCAGCCGCCCC
  //
  //  Instead, this higher scoring (longer) and noiser alignment was found.
  //
  //                                                   CGGCAGCCACCCCATCCGGGAGGGAGATGGGGGGGTCAGCCCCCCGCCCGGCCAGCCG
  //            CCCATCCGGGAGGGAGGTGGGGGGGTCAGCCCCCCGCCCCGCCAGCCGCTCCGTCCGGGAGGGAGGTGGGGGGGTCAGCCCCCCGCCCGGCCAGCCGCCCC
  //

  assert(cnspos[tiid].bgn() < cnspos[tiid].end());

  int32 cnsbgn = (cnspos[tiid].bgn() < cnspos[tiid].end()) ? cnspos[tiid].bgn() : cnspos[tiid].end();
  int32 cnsend = (cnspos[tiid].bgn() < cnspos[tiid].end()) ? cnspos[tiid].end() : cnspos[tiid].bgn();

  int32 endTrim = (cnsend - frankensteinLen) - (int32)ceil(AS_CNS_ERROR_RATE * (cnsend - cnsbgn));

  if (endTrim < 20)  endTrim = 0;


 alignFragmentAgain:
  int32 frankBgn     = MAX(0, cnspos[tiid].bgn() - bgnExtra);   //  Start position in frankenstein
  int32 frankEnd     = frankensteinLen;                       //  Truncation of frankenstein
  char  frankEndBase = 0;                                     //  Saved base from frankenstein

  bool  allowAhang   = false;
  bool  allowBhang   = true;
  bool  tryAgain     = false;

  if (showAlgorithm())
    fprintf(stderr, "alignFragment()-- Allow bgnExtra=%d and endExtra=%d (frankBgn=%d frankEnd=%d) and endTrim=%d\n",
            bgnExtra, endExtra, frankBgn, frankEnd, endTrim);

  //  If the expected fragment begin position plus any extra slop is still the begin of the
  //  consensus sequence, we allow the fragment to hang over the end.
  //
  if (frankBgn == 0)
    allowAhang = true;

  //  If the expected fragment end position plus any extra slop are less than the consensus length,
  //  we need to truncate the frankenstein so we don't incorrectly align to that.
  //
  if (cnspos[tiid].end() + endExtra < frankEnd) {
    frankEnd     = cnspos[tiid].end() + endExtra;
    frankEndBase = frankenstein[frankEnd];
    frankenstein[frankEnd] = 0;
    allowBhang   = false;
  }

  char       *aseq  = frankenstein + frankBgn;
  int32       alen  = frankEnd - frankBgn;

  abSequence *bSEQ  = abacus->getSequence(tiid);
  char       *bseq  = abacus->getBases(bSEQ);
  int32       blen  = bSEQ->length();

  if (endTrim >= blen)
    fprintf(stderr, "alignFragment()-- ERROR -- excessive endTrim %d >= length %d\n", endTrim, blen);
  if (endTrim < 0)
    fprintf(stderr, "alignFragment()-- ERROR -- negative endTrim %d\n", endTrim);

  if ((endTrim <  blen) &&
      (0       <= endTrim)) {
    int32 fragBgn      = 0;
    int32 fragEnd      = blen - endTrim;
    char  fragEndBase  = bseq[fragEnd];

    bseq[fragEnd] = 0;

    ALNoverlap  *O           = NULL;
    double       thresh      = 1e-3;
    int32        minlen      = AS_OVERLAP_MIN_LEN;

    if (O == NULL) {
      O = Optimal_Overlap_AS_forCNS(aseq,
                                    bseq,
                                    0, alen,            //  ahang bounds are unused here
                                    0, 0,               //  ahang, bhang exclusion
                                    0,
                                    AS_CNS_ERROR_RATE + 0.02, thresh, minlen,
                                    AS_FIND_ALIGN);
      if ((O) && (showAlgorithm())) {
        PrintALNoverlap("Optimal_Overlap", aseq, bseq, O);
      }
    }

    //  At 0.06 error, this equals the previous value of 10.
    double  pad = AS_CNS_ERROR_RATE * 500.0 / 3;

    if ((O) && (O->begpos < 0) && (frankBgn > 0)) {
      bgnExtra += -O->begpos + pad;
      tryAgain = true;
      O = NULL;
    }
    if ((O) && (O->endpos > 0) && (allowBhang == false)) {
      endExtra += O->endpos + pad;
      tryAgain = true;
      O = NULL;
    }
    if ((O) && (O->endpos < 0) && (endTrim > 0)) {
      endTrim -= -O->endpos + pad;
      if (endTrim < 20)
        endTrim = 0;
      tryAgain = true;
      O = NULL;
    }
    if (rejectAlignment(allowBhang, allowAhang, O))
      O = NULL;

    //  Restore the bases we might have removed.
    if (frankEndBase)   frankenstein[frankEnd] = frankEndBase;
    if (fragEndBase)    bseq[fragEnd]          = fragEndBase;

    if (O) {
      traceLen = 0;

      traceBgn = frankBgn + O->begpos;

      for (int32 *t = O->trace; (t != NULL) && (*t != 0); t++) {
        if (*t < 0)
          *t -= frankBgn;
        trace[traceLen++] = *t;
      }

      if (showAlgorithm())
        fprintf(stderr, "alignFragment()-- Alignment succeeded.\n");

      return(true);
    }

    if (tryAgain)
      goto alignFragmentAgain;
  }

  //  No alignment.  Dang.
  cnspos[tiid].bgn() = 0;
  cnspos[tiid].end() = 0;

  piid = -1;

  if (showAlgorithm())
    fprintf(stderr, "alignFragment()-- No alignment found.\n");

  return(false);
}


bool
unitigConsensus::rejectAlignment(bool allowBhang,  //  Allow a positive bhang - fragment extends past what we align to
                                 bool allowAhang,  //  Allow a negative ahang - fragment extends past what we align to
                                 ALNoverlap *O) {

  if (O == NULL) {
    if (showAlgorithm())
      fprintf(stderr, "rejectAlignment()-- No alignment found.\n");
    return(true);
  }

  //  Negative ahang?  Nope, don't want it.
  if ((O->begpos < 0) && (allowAhang == false)) {
    if (showAlgorithm())
      fprintf(stderr, "rejectAlignment()-- No alignment found -- begpos = %d (negative ahang not allowed).\n", O->begpos);
    return(true);
  }

  //  Positive bhang and not the last fragment?  Nope, don't want it.
  if ((O->endpos > 0) && (allowBhang == false)) {
    if (showAlgorithm())
      fprintf(stderr, "rejectAlignment()-- No alignment found -- endpos = %d (positive bhang not allowed).\n", O->endpos);
    return(true);
  }

  //  Too noisy?  Nope, don't want it.
  if (((double)O->diffs / (double)O->length) > AS_CNS_ERROR_RATE) {
    if (showAlgorithm())
      fprintf(stderr, "rejectAlignment()-- No alignment found -- erate %f > max allowed %f.\n",
              (double)O->diffs / (double)O->length, AS_CNS_ERROR_RATE);
    return(true);
  }

  //  Too short?  Nope, don't want it.
  if (O->length < AS_OVERLAP_MIN_LEN) {
    if (showAlgorithm())
      fprintf(stderr, "rejectAlignment()-- No alignment found -- too short %d < min allowed %d.\n",
              O->length, AS_OVERLAP_MIN_LEN);
    return(true);
  }

  return(false);
}


void
unitigConsensus::applyAlignment(int32 frag_aiid, int32 frag_ahang, int32 *frag_trace) {
  
  if (frag_aiid >= 0) {
    //  Aligned to a fragent
    assert(0);
    //  Left in because this was a useful invocation of ApplyAlignment -- but it might now be useless
    if (showAlgorithm())
      fprintf(stderr, "applyAlignment()-- aligned to fragment -- frag_aiid=%d frag_ahang=%d\n", frag_aiid, frag_ahang);
    abacus->applyAlignment(frag_aiid,
                           0, NULL,
                           tiid,
                           frag_ahang, frag_trace);

  } else {
    //  Aligned to frankenstein
    //if (showAlgorithm())
    //  fprintf(stderr, "applyAlignment()-- aligned to frankenstein\n");
    abacus->applyAlignment(abSeqID(),
                           frankensteinLen, frankensteinBof,
                           tiid,
                           traceBgn, trace);
  }
}


void
unitigConsensus::generateConsensus(void) {

  abacus->refreshMultiAlign(multialign);

  abacus->getMultiAlign(multialign)->refine(abacus, abAbacus_Smooth);
  abacus->getMultiAlign(multialign)->mergeRefine(abacus);

  abacus->getMultiAlign(multialign)->refine(abacus, abAbacus_Poly_X);
  abacus->getMultiAlign(multialign)->mergeRefine(abacus);

  abacus->getMultiAlign(multialign)->refine(abacus, abAbacus_Indel);
  abacus->getMultiAlign(multialign)->mergeRefine(abacus);

  abacus->getMultiAlign(multialign)->getConsensus(abacus, tig);
  abacus->getMultiAlign(multialign)->getPositions(abacus, tig);

  //GetMANodeConsensus(manode->lid, ma->consensus, ma->quality);
  //GetMANodePositions(manode->lid, ma);

  //  Although we generally don't care about delta values during assembly, we need them for the
  //  output, and this is the only time we compute them.  So, we've gotta hang on to them.
  //
  //for (int32 i=0; i<numfrags; i++) {
  //  utgpos[i].elta_length = 0;
  //  utgpos[i].elta        = NULL;
  //}

  //  Update or create the unitig in the MultiAlignT.

#warning not creating a unitig child
#if 0
  if (GetNumIntUnitigPoss(ma->u_list) == 0) {
    IntUnitigPos  iup;

    iup.type           = AS_OTHER_UNITIG;
    iup.ident          = tig->tigID();
    iup.position.bgn()   = 0;
    iup.position.end()   = GetMultiAlignLength(ma);
    iup.num_instances  = 0;
    iup.delta_length   = 0;
    iup.delta          = NULL;

    AppendIntUnitigPos(ma->u_list, &iup);
  } else {
    IntUnitigPos  *iup = GetIntUnitigPos(ma->u_list, 0);

    iup->position.bgn() = 0;
    iup->position.end() = GetMultiAlignLength(ma);
  }
#endif

  //PrintAlignment(stderr,manode->lid,0,-1);

  //  While we have fragments in memory, compute the microhet probability.  Ideally, this would be
  //  done in CGW when loading unitigs (the only place the probability is used) but the code wants
  //  to load sequence and quality for every fragment, and that's too expensive.
#if 0
  {
    char   **multia = NULL;
    int32  **id_array = NULL;

    int32 depth = MANode2Array(manode, &multia, &id_array,0);

    ma->data.unitig_microhet_prob = 1.0;  //AS_REZ_MP_MicroHet_prob(multia, id_array, gkpStore, frankensteinLen, depth);

    for (int32 i=0;i<depth;i++) {
      delete [] multia[2*i];
      delete [] multia[2*i+1];
      delete [] id_array[i];
    }
    delete [] multia;
    delete [] id_array;
  }
#endif

}

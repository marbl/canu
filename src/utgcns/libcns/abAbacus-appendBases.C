
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
 *    Brian P. Walenz beginning on 2015-NOV-23
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "abAbacus.H"



//  Just appends bases to the frankenstein sequence.
//  Sets read abacus pointers assuming gapless alignments.

void
abAbacus::appendBases(uint32  bid,
                      uint32  bgn,
                      uint32  end) {

  uint32      alen = numberOfColumns();

  //  Set up the b read

  abSequence *bseq  = getSequence(bid);
  uint32      blen  = bseq->length();

  //  Find the first and last columns the read 'aligns' to.  If the read isn't contained, the last
  //  column is the last column in the frankenstein.  We'll append new columns after this one.
  //
  //  Space-based, so the end element is the (n-1)th element.
  //
  abColumn *fc =                 getColumn(bgn);
  abColumn *lc = (end <= alen) ? getColumn(end-1) : getLastColumn();

  uint16    fl = UINT16_MAX;
  uint16    ll = UINT16_MAX;

  //  To get positions in output, we need to have both the first and last beads
  //  placed in the multialign.  The first bead is always aligned, but the last bead
  //  is aligned only if it is contained.

  fl = fc->alignBead(UINT16_MAX, bseq->getBase(0), bseq->getQual(0));

  if (end <= alen)
    ll = lc->alignBead(UINT16_MAX, bseq->getBase(blen-1), bseq->getQual(blen-1));

  //  If not contained, push on bases, and update the consensus base.  This is all _very_ rough.
  //  The unitig-supplied coordinates aren't guaranteed to contain 'blen' bases.  We make the
  //  obvious guess as to where the frankenstein ends in the read, and push on bases from there till
  //  the end.
  //
  //  The obvious guess is to add exactly enough bases to make the new frankenstein end at
  //  the expected length.
  //
  //                      alen
  //                      v              blen
  //  --------------------               v
  //           --------------------------
  //           ^          ^              ^
  //           bgn        bpos           end    !! end-bgn != blen !!
  //

  else
    for (uint32 bpos=blen - (end - alen); bpos<blen; bpos++) {
      abColumn *nc = new abColumn;

      ll = nc->insertAtEnd(lc, UINT16_MAX, bseq->getBase(bpos), bseq->getQual(bpos));
      lc = nc;
      //baseCallMajority(lc);
    }

  //  Now set the first/last links.

  beadID f(fc, fl);
  beadID l(lc, ll);

  readTofBead[bid] = f;  fbeadToRead[f] = bid;
  readTolBead[bid] = l;  lbeadToRead[l] = bid;

  //  If we did this correctly, then the first/last column indices should agree with the read placement.

  //assert(bgn   == getBead(bfrag->firstBead())->colIdx().get());
  //assert(end-1 == getBead(bfrag->lastBead())->colIdx().get());
}

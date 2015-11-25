
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
 *    Brian P. Walenz beginning on 2015-NOV-20
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
abAbacus::appendBases(int32     alen, abBeadID *aindex,
                      abSeqID   bfid,
                      uint32    bgn,  uint32    end) {

  //  Set up the b read

  abSequence *bfrag  = getSequence(bfid);
  uint32      blen   = bfrag->length();

  //  Find the first and last columns the read 'aligns' to.  If the read isn't contained, the last
  //  column is the last column in the frankenstein.  We'll append new columns after this one.
  //
  //  Space-based, so the end element is the (n-1)th element.
  //
  abColID fc =                 getBead(aindex[bgn  ])->colIdx();
  abColID lc = (end <= alen) ? getBead(aindex[end-1])->colIdx() : getBead(aindex[alen-1])->colIdx();

  //  To get positions in output, we need to have both the first and last beads
  //  placed in the multialign.  The first bead is always aligned, but the last bead
  //  is aligned only if it is contained.

  alignBeadToColumn(fc, bfrag->firstBead(), "");

  if (end <= alen)
    alignBeadToColumn(lc, bfrag->lastBead(), "");

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
      lc = appendColumn(lc, bfrag->getBead(bpos));
      baseCallMajority(lc);
    }

  //  If we did this correctly, then the first/last column indices should agree with the read placement.

  //assert(bgn   == getBead(bfrag->getBead(0))->colIdx().get());
  //assert(end-1 == getBead(bfrag->getBead(blen-1))->colIdx().get());

  assert(bgn   == getBead(bfrag->firstBead())->colIdx().get());
  assert(end-1 == getBead(bfrag->lastBead())->colIdx().get());

  //  Report.

  //fprintf(stderr, "appendBases()-- cols %d %d\n",
  //        getBead(bfrag->getBead(0))->colIdx().get(),
  //        getBead(bfrag->getBead(blen-1))->colIdx().get());
}

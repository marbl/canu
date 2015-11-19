
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
 *  This file is derived from:
 *
 *    src/AS_CNS/MergeRefine.C
 *    src/AS_CNS/MergeRefine.c
 *    src/AS_CNS/MultiAlignment_CNS.c
 *    src/utgcns/libcns/MergeRefine.C
 *
 *  Modifications by:
 *
 *    Michael Schatz on 2004-SEP-23
 *      are Copyright 2004 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2005-MAR-22
 *      are Copyright 2005 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2005-MAR-30 to 2008-FEB-13
 *      are Copyright 2005-2006,2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Gennady Denisov from 2005-MAY-09 to 2008-JUN-06
 *      are Copyright 2005-2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUN-16 to 2013-AUG-01
 *      are Copyright 2005-2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Aaron Halpern from 2005-SEP-29 to 2006-OCT-03
 *      are Copyright 2005-2006 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-FEB-27 to 2009-MAY-14
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-MAR-06
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-14
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static char *rcsid = "$Id$";

#include "abAbacus.H"


// test for Level 1 (neighbor) merge compatibility of cid with right
// neighbor and merge if compatible
//
static
abColID
mergeCompatible(abAbacus *abacus, abColID cid) {

  abColumn *column = abacus->getColumn(cid);

  if (column->nextID().isValid() == false)
    return(abColID());

  abColumn *merge_column = abacus->getColumn(column->nextID());

  //  If both columns have a non-gap (for a single read), we cannot merge.

  for (abBead *cbead = abacus->getBead(column->callID()); (cbead->downID().isValid()); ) {
    cbead = abacus->getBead(cbead->downID());

    if ((cbead->nextID().isValid() == true) &&
        (abacus->getBase(cbead->baseIdx()) != '-') &&
        (abacus->getBase(cbead->nextID())  != '-'))
      return(abColID());
  }

  //  OK to merge.  Merge all the bases from the right column (merge_column) to the current column.

  for (abBead *cbead = abacus->getBead(column->callID()); (cbead->downID().isValid()); ) {
    cbead = abacus->getBead(cbead->downID());

    if (cbead->nextID().isValid() == false)
      //  Nothing to merge in!
      continue;

    abBead *mbead =  abacus->getBead(cbead->nextID());

    if ((abacus->getBase(cbead->baseIdx()) != '-') ||
        (abacus->getBase(mbead->baseIdx()) == '-'))
      //  Our base isn't a gap, or their base isn't a base, nothing to do.
      continue;

    //  Exchange beads, then reset the pointer.  lateralExchangeBead() actually swaps the bead
    //  itself.  The bead in the current column is now 'mbead', so update the cbead pointer.

    abacus->lateralExchangeBead(cbead->ident(), mbead->ident());

    cbead = mbead;

    assert(cid == cbead->colIdx());
  }

  //  Finish up by removing gaps from the merged column, and moving any remaining base beads to the
  //  current column.

  abBead *mcall = abacus->getBead(merge_column->callID());

  while (mcall->downID().isValid()) {
    abBead *mbead = abacus->getBead(mcall->downID());

    //  Remove the bead from the old column.  After this, mcall->downID() points to the next
    //  bead in the column.  We're iterating down the column by popping beads out from the head.

    abacus->getBead(mbead->upID())->downID() = mbead->downID();

    if (mbead->downID().isValid())
      abacus->getBead(mbead->downID())->upID() = mbead->upID();

    //  Technically correct, but it's private, and we're just killing the column anyway.
    //abacus->getColumn(mbead->column_index)->base_count.DecBaseCount( abacus->getBase(mbead->soffset));

    mbead->upID()   = abBeadID();
    mbead->downID() = abBeadID();
    mbead->colIdx() = abColID();

    //  Base is not a gap.  Move bead from old column to current column.

    if (abacus->getBase(mbead->baseIdx()) != '-') {
      abacus->alignBeadToColumn(cid, mbead->ident(), "mergeCompatible()");
    }

    //  Base is a gap, just delete it.

    else {
      if (mbead->prevID().isValid())   abacus->getBead(mbead->prevID())->nextID() = mbead->nextID();
      if (mbead->nextID().isValid())   abacus->getBead(mbead->nextID())->prevID() = mbead->prevID();

      mbead->clear();
    }
  }

  //  Now, we're left with a column with just the consensus call.  Yank the call out of the consensus.

  if (mcall->prevID().isValid())   abacus->getBead(mcall->prevID())->nextID() = mcall->nextID();
  if (mcall->nextID().isValid())   abacus->getBead(mcall->nextID())->prevID() = mcall->prevID();

  mcall->clear();

  //  And, finally, yank the column out.

  if (merge_column->prevID().isValid())
    abacus->getColumn(merge_column->prevID())->nextID() = merge_column->nextID();

  if (merge_column->nextID().isValid())
    abacus->getColumn(merge_column->nextID())->prevID() = merge_column->prevID();

  //  We have merged.  Return the now obsolete column ID.

  return(merge_column->ident());
}


//*********************************************************************************
// Simple sweep through the MultiAlignment columns, looking for columns
// to merge and removing null columns
//*********************************************************************************

void
abMultiAlign::mergeRefine(abAbacus *abacus, bool highQuality) {

  //fprintf(stderr, "abMultiAlign::mergeRefine()--  legnth %u\n", length());

  //  Loop over all columns.  If we do not merge, advance to the next column, otherwise, stay here
  //  and merge to the now different next column (mergeCompatible removes the column that gets
  //  merged into the current column).

  bool  somethingMerged = false;

  for (abColID cid=first; cid.isValid(); ) {
    abColID  mergedCol = mergeCompatible(abacus, cid);

    //  If we didn't merge anything, move to the next column and keep going.

    if (mergedCol == abColID()) {
      cid = abacus->getColumn(cid)->nextID();
      continue;
    }

    //  We merged!  Update the base call.

    somethingMerged = true;

    abacus->baseCall(cid, highQuality);

    //  And reset our first/last pointers as needed.  first shouldn't change, but last can disappear.

#if 0
    if (first == mergedCol)  first = abacus->getColumn(mergedCol)->nextID();
    if (last  == mergedCol)   last = abacus->getColumn(mergedCol)->prevID();

    assert(abacus->getColumn(first)->prevID() == abColID());
    assert(abacus->getColumn(last)->nextID()  == abColID());
#endif
  }

  //  If any merges were performed, refresh.  This updates the column list and does a few checks.

  if (somethingMerged)
    abacus->refreshMultiAlign(lid, false, highQuality);
}

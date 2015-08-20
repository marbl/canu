
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
 *    Brian P. Walenz beginning on 2014-NOV-17
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
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
bool
mergeCompatible(abAbacus *abacus, abColID cid) {

  abColumn *column = abacus->getColumn(cid);

  if (column->nextID().isValid() == false)
    return(false);

  abColumn *merge_column = abacus->getColumn(column->nextID());

  //CheckColumnBaseCount(column);
  //CheckColumnBaseCount(merge_column);

  abBead *cbead = abacus->getBead(column->callID());

  //  If both columns have a non-gap, we cannot merge.

  while (cbead->downID().isValid()) {
    cbead = abacus->getBead(cbead->downID());
    if (cbead->nextID().isValid()) {
      abBead *mbead =  abacus->getBead(cbead->nextID());
      if ((abacus->getBase(cbead->baseIdx()) != '-') &&
          (abacus->getBase(mbead->baseIdx()) != '-'))
        return(0);
    }
  }

#if 0
  {
    abColumn *l = column;
    abColumn *r = merge_column;

    char lc = abacus->getBase(abacus->getBead(l->call)->baseIdx());
    char rc = abacus->getBase(abacus->getBead(r->call)->baseIdx());

    fprintf(stderr, "mergeCompatible()-- l col=%d %c r col=%d %c\n", cid, lc, l->nextID(), rc);
  }
#endif


  //  do the merge - merge all the bases in merge_column (on the
  //  right) to column (on the left).

  cbead = abacus->getBead(column->callID());
  while (cbead->downID().isValid()) {
    cbead = abacus->getBead(cbead->downID());

    if (cbead->nextID().isValid()) {
      abBead *mbead =  abacus->getBead(cbead->nextID());

      //fprintf(stderr, "merge? %c -- %c\n",
      //        abacus->getBase(cbead->baseIdx()),
      //        abacus->getBase(mbead->baseIdx()));

      if ((abacus->getBase(cbead->baseIdx()) == '-') &&
          (abacus->getBase(mbead->baseIdx()) != '-')) {
        //fprintf(stderr, "merge  mbead from cid=%d to   cid=%d %c\n",
        //        mbead->column_index, cbead->column_index, abacus->getBase(mbead->baseIdx()));

        abacus->lateralExchangeBead(cbead->ident(), mbead->ident());

        //  LateralExchangeBead() moves contents.  Our pointers are
        //  now backwards.  We only care about cbead though.

        cbead = mbead;
      }
    }
  }

  //  finish up by moving any remaining guys from the right column to
  //  the left column.  We also delete '-' entries from the
  //  merge_column.

  assert(cid == cbead->colIdx());

  {
    abBead *mcall = abacus->getBead(merge_column->callID());  //  The consensus call of the column

    //fprintf(stderr, "loop depth=%d gap=%d\n", GetDepth(merge_column), GetBaseCount(&merge_column->base_count,'-'));

    while (mcall->downID().isValid()) {
      abBead *mbead = abacus->getBead(mcall->downID());

      //  If the mbead is not a gap, move it over to the left column,
      //  otherwise just yank it out.
      //
      if (abacus->getBase(mbead->baseIdx()) != '-') {
        //fprintf(stderr, "move bead from %d to cid=%d %c\n",
        //        mbead->column_index, cid, abacus->getBase(mbead->baseIdx()));

        abacus->unalignBeadFromColumn(mbead->ident());
        abacus->alignBeadToColumn(cid, mbead->ident(), "mergeCompatible()");
      } else {
        //fprintf(stderr, "delete bead from %d (gap)\n",
        //        mbead->column_index);

        // heal wound left by lateral removal
        if (mbead->prevID().isValid() ) abacus->getBead(mbead->prevID())->nextID() = mbead->nextID();
        if (mbead->nextID().isValid() ) abacus->getBead(mbead->nextID())->prevID() = mbead->prevID();

        abacus->unalignBeadFromColumn(mbead->ident());
        mbead->clear();
      }
    }

    // heal wound left by lateral removal of mcall
    //
    if (mcall->prevID().isValid() ) abacus->getBead(mcall->prevID())->nextID() = mcall->nextID();
    if (mcall->nextID().isValid() ) abacus->getBead(mcall->nextID())->prevID() = mcall->prevID();

    mcall->clear();

    // reset column pointers to bypass the removed column
    //
    if (merge_column->prevID().isValid())
      abacus->getColumn(merge_column->prevID())->nextID() = merge_column->nextID();

    if (merge_column->nextID().isValid())
      abacus->getColumn(merge_column->nextID())->prevID() = merge_column->prevID();

    //ClearColumn();
  }

  //  We have merged.
  return(true);
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

  for (abColID cid=first; cid.isValid(); ) {
    if (mergeCompatible(abacus, cid) == false)
      cid = abacus->getColumn(cid)->nextID();
    else
      abacus->baseCall(cid, highQuality);
  }

  //  Recall all the bases.  Ideally, this should only be done for columns that change, which we should
  //  know in the loop above.
  //
  //abacus->refreshMultiAlign(lid, highQuality);
}

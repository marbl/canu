
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


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

//  Recall consensus in all columns.
//  Rebuild the list of columns in the multialign.
//  Rebuild the column to position map.

void
abAbacus::refreshMultiAlign(abMultiAlignID  mid,
                            bool            recallBase,         //  (false) If true, recall the base
                            bool            highQuality) {      //  (false) If true, use the high quality base call algorithm

  abMultiAlign *ma = getMultiAlign(mid);

  //fprintf(stderr, "abMultiAlign::refreshMultiAlign()--  legnth %u recallBase=%d highQuality=%d\n",
  //        ma->length(), recallBase, highQuality);

  ma->columns().clear();  // columnList

  uint32   index = 0;
  abColID  cid   = ma->firstColumn();

  while (cid.isValid()) {
    abColumn  *column = getColumn(cid);

    if (recallBase == true)
      baseCall(cid, highQuality);

    column->ma_position = index;

    ma->columns().push_back(column->ident());

    assert(ma->columns()[index] == column->ident());  //  That the thing we just pushed on is at the correct spot

    cid = column->next;
    index++;
  }

  //  Check column pointers

  {
    abColID   pid = ma->firstColumn();  //  Previous column
    abColumn *pol = getColumn(pid);

    abColID   cid = pol->next;          //  Next column
    abColumn *col = getColumn(cid);

    while (cid.isValid()) {
      assert(pol->next == cid);  //  We're iterating over this, must be true.
      assert(pid == col->prev);  //  What we're testing

      pid = cid;
      pol = col;

      cid = pol->next;          //  Next column
      col = getColumn(cid);
    }
  }
}

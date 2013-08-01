
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

static char *rcsid = "$Id: MergeRefine.c,v 1.3 2011-01-03 03:07:16 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"
#include "MicroHetREZ.H"
#include "AS_UTL_reverseComplement.H"



// test for Level 1 (neighbor) merge compatibility of cid with right
// neighbor and merge if compatible
//
static
int
MergeCompatible(int32 cid) {

  Column *column = GetColumn(columnStore,cid);
  assert(column != NULL);
  if (column->next == -1)
    return(0);

  Column *merge_column = GetColumn(columnStore,column->next);
  assert(merge_column != NULL);

  //CheckColumnBaseCount(column);
  //CheckColumnBaseCount(merge_column);

  Bead *cbead = GetBead(beadStore,column->call);

  //  If both columns have a non-gap, we cannot merge.

  while (cbead->down.isValid()) {
    cbead = GetBead(beadStore,cbead->down);
    if (cbead->next.isValid()) {
      Bead *mbead =  GetBead(beadStore, cbead->next);
      if ((*Getchar(sequenceStore,cbead->soffset) != '-') &&
          (*Getchar(sequenceStore,mbead->soffset) != '-'))
        return(0);
    }
  }

#if 0
  {
    Column *l = column;
    Column *r = merge_column;

    char lc = *Getchar(sequenceStore, GetBead(beadStore, l->call)->soffset);
    char rc = *Getchar(sequenceStore, GetBead(beadStore, r->call)->soffset);

    fprintf(stderr, "MergeCompatible()-- l col=%d %c r col=%d %c\n", cid, lc, l->next, rc);
  }
#endif


  //  do the merge - merge all the bases in merge_column (on the
  //  right) to column (on the left).

  cbead = GetBead(beadStore,column->call);
  while (cbead->down.isValid()) {
    cbead = GetBead(beadStore,cbead->down);

    if (cbead->next.isValid()) {
      Bead *mbead =  GetBead(beadStore, cbead->next);

      //fprintf(stderr, "merge? %c -- %c\n",
      //        *Getchar(sequenceStore,cbead->soffset),
      //        *Getchar(sequenceStore,mbead->soffset));

      if ((*Getchar(sequenceStore,cbead->soffset) == '-') &&
          (*Getchar(sequenceStore,mbead->soffset) != '-')) {
        //fprintf(stderr, "merge  mbead from cid=%d to   cid=%d %c\n",
        //        mbead->column_index, cbead->column_index, *Getchar(sequenceStore,mbead->soffset));

        LateralExchangeBead(cbead->boffset, mbead->boffset);

        //  LateralExchangeBead() moves contents.  Our pointers are
        //  now backwards.  We only care about cbead though.

        cbead = mbead;
      }
    }
  }

  //  finish up by moving any remaining guys from the right column to
  //  the left column.  We also delete '-' entries from the
  //  merge_column.

  assert(cid == cbead->column_index);

  {
    Bead *mcall = GetBead(beadStore,merge_column->call);  //  The consensus call of the column

    //fprintf(stderr, "loop depth=%d gap=%d\n", GetDepth(merge_column), GetBaseCount(&merge_column->base_count,'-'));

    while (mcall->down.isValid()) {
      Bead *mbead = GetBead(beadStore, mcall->down);

      //  If the mbead is not a gap, move it over to the left column,
      //  otherwise just yank it out.
      //
      if (*Getchar(sequenceStore,mbead->soffset) != '-') {
        //fprintf(stderr, "move bead from %d to cid=%d %c\n",
        //        mbead->column_index, cid, *Getchar(sequenceStore,mbead->soffset));

        UnAlignBeadFromColumn(mbead->boffset);
        AlignBeadToColumn(cid, mbead->boffset, "MergeCompatible()");
      } else {
        //fprintf(stderr, "delete bead from %d (gap)\n",
        //        mbead->column_index);

        // heal wound left by lateral removal
        if (mbead->prev.isValid() ) GetBead(beadStore,mbead->prev)->next = mbead->next;
        if (mbead->next.isValid() ) GetBead(beadStore,mbead->next)->prev = mbead->prev;

        UnAlignBeadFromColumn(mbead->boffset);
        ClearBead(mbead->boffset);
      }
    }

    // heal wound left by lateral removal of mcall
    //
    if (mcall->prev.isValid() ) GetBead(beadStore,mcall->prev)->next = mcall->next;
    if (mcall->next.isValid() ) GetBead(beadStore,mcall->next)->prev = mcall->prev;

    ClearBead(mcall->boffset);

    // reset column pointers to bypass the removed column
    //
    if (merge_column->prev != -1) GetColumn(columnStore,merge_column->prev)->next = merge_column->next;
    if (merge_column->next != -1) GetColumn(columnStore,merge_column->next)->prev = merge_column->prev;

    //ClearColumn();
  }

  //  We have merged.
  return(1);
}


//*********************************************************************************
// Simple sweep through the MultiAlignment columns, looking for columns
// to merge and removing null columns
//*********************************************************************************

void
MergeRefine(int32 mid, VA_TYPE(IntMultiVar) *v_list,
            int32 utg_alleles, CNS_Options *opp, int32 get_scores) {
  int32   cid = 0;
  MANode *ma  = GetMANode(manodeStore,mid);

  assert(ma != NULL);

  //  Loop over all columns.  If do not merge, advance to the next
  //  column, otherwise, stay here and merge to the now different next
  //  column (MergeCompatible removes the column that gets merged into
  //  the current column).
  //
  for (cid=ma->first; cid != -1; ){
    if (MergeCompatible(cid) == 0)
      cid = GetColumn(columnStore,cid)->next;
  }

  {
    IntMultiVar *vl=NULL;
    int32 nv=0;
    int32 make_v_list=0;

    if (utg_alleles)
      make_v_list = 1;
    else if (v_list)
      make_v_list = 2;

    RefreshMANode(mid, 1, opp, &nv, &vl, make_v_list, get_scores);

    if (make_v_list && v_list) {
      ResetVA_IntMultiVar(v_list);
        
      if (nv > 0) {
        ResetVA_IntMultiPos(v_list);
        SetRangeVA_IntMultiVar(v_list, 0, nv, vl);
      }
    }

    safe_free(vl);
  }
}

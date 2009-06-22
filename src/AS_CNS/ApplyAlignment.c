
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

static char *rcsid = "$Id: ApplyAlignment.c,v 1.2 2009-06-22 12:04:53 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"




static
int32
 ColumnPrepend(int32 cid, int32 bid) {
  // bid is the offset of the Bead seeding the column

  ColumnBeadIterator ci;
  int32 nid;

  Bead *bead = GetBead(beadStore,bid);
  assert(bead != NULL);

  //fprintf(stderr, "ColumnPrepend()-- adding column for bid %d\n", bid);

  Column *column = CreateColumn(bid);
  assert(column != NULL);

  Bead   *call     = GetBead(beadStore,column->call);
  Column *next     = GetColumn(columnStore,cid);
  Bead   *nextcall = GetBead(beadStore,next->call);

  column->prev = next->prev;
  column->next = cid;
  call->prev = nextcall->prev;
  call->next = nextcall->boffset;
  next->prev = column->lid;
  nextcall->prev = call->boffset;

  if (column->prev != -1)
    GetColumn(columnStore,column->prev)->next = column->lid;

  if (call->prev != -1)
    GetBead(beadStore,call->prev)->next = call->boffset;

  CreateColumnBeadIterator(cid,&ci);

  while ( (nid = NextColumnBead(&ci)) != -1 ) {
    bead = GetBead(beadStore,nid);
    if ( bead->prev != -1 && bead->prev != bid) {
      AlignBeadToColumn(column->lid,PrependGapBead(nid), "ColumnPrepend()");
    }
  }
  column->ma_id =  next->ma_id;
  column->ma_index =  next->ma_index - 1;
  AddColumnToMANode(column->ma_id,*column);
  if ( column->prev == -1 ) {
    GetMANode(manodeStore,column->ma_id)->first = column->lid;
  }
  return column->lid;
}



static
void
CheckColumns(void) {
  int cid;
  int fails = 0;

  for (cid=0; cid<GetNumColumns(columnStore); cid++) {
    Column *column = GetColumn(columnStore,cid);

    if (column == NULL)
      continue;

    ColumnBeadIterator ci;
    CreateColumnBeadIterator(cid,&ci);

    Bead *call = GetBead(beadStore,column->call);
    int32 bid;

    //if (cid == 548)
    //  fprintf(stderr, "CheckColumns()-- column %d\n", cid);

    while ( (bid = NextColumnBead(&ci)) != -1 ) {
      Bead *bead = GetBead(beadStore,bid);

      //if ((cid == 548) || (bid == 3919) || (bid == 739) || (bid == 508) || (bid == 305) || (bid == 187))
      //  fprintf(stderr, "bead %d claims col %d (cid=%d)\n", bead->boffset, bead->column_index, cid);

      if (bead->column_index != cid)
        fails++;
      //assert(bead->column_index == cid);
    }
  }
  assert(fails == 0);
}



void
ApplyAlignment(int32 afid,
               int32 alenUNUSED, int32 *aindexUNUSED,
               int32 bfid,
               int32 ahang,
               int32 *trace) {

  Fragment *afrag = NULL;
  Fragment *bfrag = NULL;
  int32 aboffset, bboffset; // offsets of first beads in fragments
  int32 apos, bpos;         // local offsets as alignment progresses
  int32 alen, blen;

  int32 ovl_remaining, column_appends, column_index;
  int32 first_touched_column;
  int32 last_a_aligned,last_b_aligned;
  int32 next_to_align;
  int32 binsert;

  int32 *aindex    = NULL;

  Bead *abead;

  int32 off;

  assert(afid >= 0);

  afrag= GetFragment(fragmentStore,afid);
  assert(afrag != NULL);
  alen     = afrag->length;
  aboffset = afrag->firstbead;

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "allocate aindex[] with %d things.\n", alen);
#endif
  aindex    = (int32 *)safe_malloc(alen*sizeof(int32));

  {
    Bead *ab=GetBead(beadStore,aboffset);
    int ai;
    for (ai=0;ai<alen;ai++)
      aindex[ai] = aboffset+ai;
  }

  bfrag= GetFragment(fragmentStore,bfid);
  assert(bfrag != NULL);

  blen     = bfrag->length;
  bboffset = bfrag->firstbead;

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "afid=%d aboffset=%d alen=%d  bfid=%d bboffset=%d blen=%d\n",
          afid, aboffset, alen,
          bfid, bboffset, blen);
#endif

  //  All the a beads should be in a column.  All the b beads should not.
  //
  {
    Bead  *b = GetBead(beadStore, aboffset);
    int    e = 0;

    while (b) {
      if (b->column_index == -1) {
        e++;
        fprintf(stderr, "bead %d in A has undef column_index.\n", b->boffset);
      }
      b = (b->next == -1) ? NULL : GetBead(beadStore, b->next);
    }

    b = GetBead(beadStore, bboffset);

    while (b) {
      if (b->column_index != -1) {
        e++;
        fprintf(stderr, "bead %d in B has defined column_index %d.\n", b->boffset, b->column_index);
      }
      b = (b->next == -1) ? NULL : GetBead(beadStore, b->next);
    }

    assert(e == 0);
  }

  last_a_aligned = -1;
  last_b_aligned = -1;

  apos = MAX(ahang,0);
  bpos = 0;

  if ( ahang == alen ) { // special case where fragments abutt
    assert(alen-1 < alen);
    abead = GetBead(beadStore,aindex[alen-1]);
  } else {
    assert(apos < alen);
    abead = GetBead(beadStore,aindex[apos]);
  }

  first_touched_column = abead->column_index;

  if ( ahang < 0 ) {
    Bead *gbead = GetBead(beadStore,bboffset);
    while ( bpos < -ahang ) {
      ColumnPrepend(first_touched_column,bboffset+bpos);
      bpos++;
    }
    last_b_aligned = bboffset+bpos-1;
  }

  assert(apos < alen);
  assert(bpos < blen);

  last_a_aligned = GetBead(beadStore,aindex[apos])->prev;

  while ( (NULL != trace) && *trace != 0 ) {
#ifdef DEBUG_ABACUS_ALIGN
    fprintf(stderr, "trace=%d  apos=%d alen=%d bpos=%d blen=%d\n", *trace, apos, alen, bpos, blen);
#endif

    if ( *trace < 0 ) {
      //
      // gap is in afrag
      //
      // align ( - *trace - apos ) positions
      while ( apos < (- *trace - 1)) {
        assert(apos < alen);
        assert(bpos < blen);
        abead = GetBead(beadStore,aindex[apos]);
        AlignBeadToColumn(abead->column_index, bboffset+bpos, "ApplyAlignment(1)");
        last_a_aligned = aindex[apos];
        last_b_aligned = bboffset+bpos;
        apos++; bpos++;
        binsert = bboffset+bpos-1;
        while ( abead->next > -1 && (abead = GetBead(beadStore,abead->next))->boffset != aindex[apos] ) {
          // insert a gap bead in b and align to
          off = abead->boffset;
          binsert = AppendGapBead(binsert);
          abead = GetBead(beadStore, off);
          AlignBeadToColumn(abead->column_index, binsert, "ApplyAlignment(2)");
          last_a_aligned = abead->boffset;
          last_b_aligned = binsert;
        }
      }

      // insert a gap column to accommodate bpos "insert"
      //   via:
      // insert a gap in afrag; insert new column seeded with that
      // then align bpos to that column
      //                           apos
      //                           | | | | |
      //               | | | | | | a a a a a
      //               a a a a a a
      //                   b b b b
      //                           b b b b
      //                           bpos
      //
      //                         |
      //                         V
      //                             apos
      //                           * | | | | |
      //               | | | | | | | a a a a a
      //               a a a a a a -
      //                   b b b b b
      //                           * b b b
      //                             bpos
      //              * is new column

      assert(apos     < alen);
      assert(bpos - 1 < blen);
      abead = GetBead(beadStore,aindex[apos]);
      // in case the last aligned bead in a is not apos->prev
      //   (Because gap beads were inserted, for example)
      binsert = bboffset+bpos-1;
      while ( abead->prev != last_a_aligned ) {
        binsert = AppendGapBead(binsert);
        assert(apos < alen);
        abead = GetBead(beadStore,aindex[apos]);
        next_to_align = (GetBead(beadStore,last_a_aligned))->next;
        AlignBeadToColumn( (GetBead(beadStore,next_to_align))->column_index, binsert, "ApplyAlignment(3)");
        last_a_aligned = next_to_align;
        last_b_aligned = binsert;
      }
      assert(apos < alen);
      assert(bpos < blen);
      ColumnAppend((GetColumn(columnStore,abead->column_index))->prev,bboffset+bpos);
      abead = GetBead(beadStore,aindex[apos]);
      last_a_aligned = abead->prev;
      last_b_aligned = bboffset+bpos;
      bpos++;
    } else {
      //
      // gap is in bfrag
      //
      // align ( *trace - bpos ) positions
      while ( bpos < (*trace - 1) ) {
        assert(apos < alen);
        assert(bpos < blen);
        abead = GetBead(beadStore,aindex[apos]);
        AlignBeadToColumn(abead->column_index, bboffset+bpos, "ApplyAlignment(4)");
        last_a_aligned = aindex[apos];
        last_b_aligned = bboffset+bpos;
        apos++; bpos++;
        binsert = bboffset+bpos-1;
        assert((abead->next <= -1) || (apos < alen));
        assert(bpos < blen);
        while ( abead->next > -1 && (abead = GetBead(beadStore,abead->next))->boffset != aindex[apos] ) {
          // insert a gap bead in b and align to
          off = abead->boffset;
          binsert = AppendGapBead(binsert);
          abead = GetBead(beadStore, off);
          AlignBeadToColumn(abead->column_index, binsert, "ApplyAlignment(5)");
          last_a_aligned = abead->boffset;
          last_b_aligned = binsert;
        }
      }

      // insert a gap bead at bpos to represent bpos "delete"
      // and align the gap position with abead
      //                           apos
      //                           | | | | |
      //               | | | | | | a a a a a
      //               a a a a a a
      //                   b b b b
      //                           b b b b
      //                           bpos
      //
      //                         |
      //                         V
      //                             apos
      //                             | | | |
      //               | | | | | | | a a a a
      //               a a a a a a a
      //                   b b b b -
      //                             b b b b
      //                             bpos
      //              (no new column is required)

      off = abead->boffset;
      binsert = AppendGapBead(last_b_aligned);
      abead = GetBead(beadStore,off);
      AlignBeadToColumn(abead->column_index, binsert, "ApplyAlignment(6)");
      last_a_aligned = abead->boffset;
      last_b_aligned = binsert;
      apos++;

      //assert((abead->next <= -1) || (apos < alen));
      assert(bpos < blen);

      //  (apos < alen) below will quit the loop if the alignment has
      //  gaps at the end.  This is probably a bug in the aligner
      //  code.  Enabling the above assert will usually find an
      //  example quickly.

      while ((abead->next > -1) &&
             (apos < alen) &&
             (abead = GetBead(beadStore,abead->next))->boffset != aindex[apos] ) {
        // insert a gap bead in b and align to
        off = abead->boffset;
        binsert = AppendGapBead(binsert);
        abead = GetBead(beadStore, off);
        AlignBeadToColumn(abead->column_index, binsert, "ApplyAlignment(7)");
        last_a_aligned = abead->boffset;
        last_b_aligned = binsert;
      }
    }

    trace++;
  }

  // remaining alignment contains no indels

  ovl_remaining  = (blen-bpos < alen-apos) ? blen-bpos : alen-apos;
  while ( ovl_remaining-- > 0 ) {
    assert(apos < alen);
    assert(bpos < blen);
    abead = GetBead(beadStore,aindex[apos]);
    AlignBeadToColumn(abead->column_index, bboffset+bpos, "ApplyAlignment(8)");
    last_a_aligned = abead->boffset;
    last_b_aligned = bboffset+bpos;
    apos++;bpos++;
    binsert = bboffset+bpos-1;

    while ( abead->next > -1 && (apos < alen) && (abead = GetBead(beadStore,abead->next))->boffset != aindex[apos] ) {
      off = abead->boffset;
      binsert = AppendGapBead(binsert);
      abead = GetBead(beadStore, off);
      AlignBeadToColumn(abead->column_index, binsert, "ApplyAlignment(9)");
      last_a_aligned = abead->boffset;
      last_b_aligned = binsert;
    }
  }

  column_appends = blen-bpos;
  column_index = abead->column_index;

  if ( column_appends > 0 ) {
    // First, if there are any previously aligned columns to right of abead
    // insert gaps into b to align with these columns
    Column *pcol = GetColumn(columnStore,column_index);
    while ( pcol->next != -1 )  {
      binsert = AppendGapBead(binsert);
      column_index = pcol->next;
      AlignBeadToColumn(column_index, binsert, "ApplyAlignment(A)");
      pcol = GetColumn(columnStore,column_index);
    }
    // then, add on trailing (dovetail) beads from b
    while (column_appends-- > 0 ) {
      column_index = ColumnAppend(column_index,bboffset+bpos);
      bpos++;
    }
  }

  //CheckColumns();

  safe_free(aindex);
  bfrag->manode=afrag->manode;
}

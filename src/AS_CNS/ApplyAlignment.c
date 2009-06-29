
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

static char *rcsid = "$Id: ApplyAlignment.c,v 1.5 2009-06-29 18:41:16 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"


#undef DEBUG_FIND_BEAD
#undef DEBUG_ALIGN_GAPS
#undef DEBUG_ALIGN_POSITION

//  Add a column before cid, seeded with bead bid.
//
static
int32
ColumnPrepend(int32 cid, int32 bid) {

  ColumnBeadIterator ci;
  int32 nid;

  Bead   *bead     = GetBead(beadStore, bid);
  Column *column   = CreateColumn(bid);
  Bead   *call     = GetBead(beadStore, column->call);
  Column *next     = GetColumn(columnStore, cid);
  Bead   *nextcall = GetBead(beadStore, next->call);

  //fprintf(stderr, "ColumnPrepend()-- adding column for bid %d\n", bid);

  column->prev   = next->prev;
  column->next   = cid;

  call->prev     = nextcall->prev;
  call->next     = nextcall->boffset;
  next->prev     = column->lid;

  nextcall->prev = call->boffset;

  if (column->prev != -1)
    GetColumn(columnStore,column->prev)->next = column->lid;

  if (call->prev != -1)
    GetBead(beadStore,call->prev)->next = call->boffset;

  CreateColumnBeadIterator(cid, &ci);

  while ((nid = NextColumnBead(&ci)) != -1) {
    bead = GetBead(beadStore, nid);

    //  The original version would not insert a gap bead at the start of a fragment.

    if ((bead->prev != -1) &&
        (bead->prev != bid))
      AlignBeadToColumn(column->lid, PrependGapBead(nid), "ColumnPrepend()");
  }

  column->ma_id    = next->ma_id;
  column->ma_index = next->ma_index - 1;

  AddColumnToMANode(column->ma_id,*column);

  if (column->prev == -1)
    GetMANode(manodeStore,column->ma_id)->first = column->lid;

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

    while ( (bid = NextColumnBead(&ci)) != -1 ) {
      Bead *bead = GetBead(beadStore,bid);

      if (bead->column_index != cid)
        fails++;
    }
  }
  assert(fails == 0);
}




//  All the a beads should be in a column.  All the b beads should not.
//
static
void
sanityCheck(int32 aboffset, int32 bboffset) {
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




//  Returns the boffset of the bead that is:
//    in the same column as bead bi
//    in the same fragment as bead fi
//
static
int32
findBeadInColumn(int32 bi, int32 fi) {

#ifdef DEBUG_FIND_BEAD
  fprintf(stderr, "findBeadInColumn bead bi=%d bead fi=%d\n", bi, fi);
#endif

  if ((fi == -1) || (bi == -1))
    return(-1);

  fi = GetBead(beadStore, fi)->frag_index;

#ifdef DEBUG_FIND_BEAD
  fprintf(stderr, "findBeadInColumn fragindex fi=%d\n", fi);
#endif

  Bead *b = GetBead(beadStore, bi);

  if (b->frag_index == fi)
    return(bi);

  //  Search up.  The way we call findBeadInColumn, the one we're
  //  looking for is usually up from where we start.
  while (b->up != -1) {
    b = GetBead(beadStore, b->up);
#ifdef DEBUG_FIND_BEAD
    fprintf(stderr, "findBeadInColumn up bead=%d fi=%d\n", b->boffset, b->frag_index);
#endif
    if (b->frag_index == fi)
      return(b->boffset);
  }

  //  Search down.
  b = GetBead(beadStore, bi);

  while (b->down != -1) {
    b = GetBead(beadStore, b->down);
#ifdef DEBUG_FIND_BEAD
    fprintf(stderr, "findBeadInColumn down bead=%d fi=%d\n", b->boffset, b->frag_index);
#endif
    if (b->frag_index == fi)
      return(b->boffset);
  }

  //  Give up.  See comments in MultiAlignUnitig around line 733, where we append new sequence to
  //  the start of frankenstein.

  assert(0);
  return(-1);
}



//  Add gaps already in abacus to our read -- complicated by possibly having to switch the fragment
//  the abead is working on.
//
//  Assumes we just aligned the 'a' to 'b' at the start.  Abacus already has a bunch of gap
//  beads/columns in it (these are NOT gaps in the alignment, or gaps in the sequence; they're
//  caused by other alignments) and we must add these columns to the B read.
//
//  A   a-----a
//  B   b
//
static
void
alignGaps(int32 *aindex, int32 &apos, int32  alen,
          int32 *bindex, int32 &bpos, int32  blen,
          int32 &lasta,
          int32 &lastb) {

  //  findBeadInColumn() will move us to the next 'a' fragment when we are at the end of the current
  //  one.
  //
  //  next    ->       .....x--********
  //  current ->  ***********--
  //
  //  * -- beads that we align to.
  //  x -- bead 'nn' below.
  //  - -- gap in abacus.
  //
  //  We're at the last * on the second row, the next bead we align to is the first * on the first
  //  row.  There might be gaps already in abacus that we need to add.  We move to the 'x', then
  //  walk along this 'a' fragment adding gaps.

#ifdef DEBUG_ALIGN_GAPS
  fprintf(stderr, "alignGaps()-- apos=%d alen=%d  bpos=%d blen=%d  lasta=%d  lastb=%d\n",
          apos, alen, bpos, blen, lasta, lastb);
#endif

  if (apos >= alen)
    return;

  lasta = findBeadInColumn(lasta, aindex[apos]);

  if (lasta == -1)
    return;

  int32 nexta = GetBead(beadStore, lasta)->next;

  //  If the next bead (from 'x' in the picture) is NOT the next base in A, add gaps until that is
  //  so.

  while (nexta != aindex[apos]) {
#ifdef DEBUG_ALIGN_GAPS
    fprintf(stderr, "alignGaps()-- lasta=%d  nexta=%d  search for aindex[apos]=%d  \n", lasta, nexta, aindex[apos]);
#endif

    lastb = AppendGapBead(lastb);
    lasta = nexta;

    Bead *bead = GetBead(beadStore, nexta);  //  AppendGapBead might reallocate beadStore

    AlignBeadToColumn(bead->column_index, lastb, "ApplyAlignment(alignGaps)");

    nexta = bead->next;
  }
}


//  Aligns a bead in B to the existing column in A.
static
void
alignPosition(int32 *aindex, int32 &apos, int32  alen,
              int32 *bindex, int32 &bpos, int32  blen,
              int32 &lasta,
              int32 &lastb,
              char  *label) {

  assert(apos < alen);
  assert(bpos < blen);

  AlignBeadToColumn(GetBead(beadStore, aindex[apos])->column_index,
                    bindex[bpos],
                    label);

  lasta = aindex[apos];
  lastb = bindex[bpos];

  apos++;
  bpos++;

#ifdef DEBUG_ALIGN_POSITION
  fprintf(stderr, "alignPosition()-- apos=%d bpos=%d lasta=%d lastb=%d\n",
          apos, bpos, lasta, lastb);
#endif

  alignGaps(aindex, apos, alen, bindex, bpos, blen, lasta, lastb);

  //assert(GetBead(beadStore, aindex[apos])->prev == bead->bindex);
  //assert(GetBead(beadStore, bindex[bpos])->prev == lastb);
}




void
ApplyAlignment(int32 afid,
               int32 alen, int32 *aindex,
               int32 bfid,
               int32 ahang,
               int32 *trace) {

  int32  bpos   = 0;
  int32  blen   = 0;

  int32  apos   = MAX(ahang,0);
  int32 *bindex = NULL;

  int32  lasta  = -1;
  int32  lastb  = -1;

  if (afid >= 0) {
    Fragment *afrag = GetFragment(fragmentStore,afid);
    alen            = afrag->length;
    aindex          = (int32 *)safe_malloc(alen * sizeof(int32));

    //  The loop really abuses the fact that all the beads for the bases in this read are
    //  contiguous.  They're contiguous because CreateMANode() (I think it's in there) allocated
    //  them in one block.
    //
    for (int ai=0;ai<alen;ai++)
      aindex[ai] = afrag->firstbead + ai;
  } else {
    assert(aindex != NULL);
  }


  if (bfid >= 0) {
    Fragment *bfrag = GetFragment(fragmentStore,bfid);
    blen            = bfrag->length;
    bindex          = (int32 *)safe_malloc(blen * sizeof(int32));

    for (int bi=0;bi<blen;bi++)
      bindex[bi] = bfrag->firstbead + bi;

    //  USED?
    bfrag->manode = NULL;
  } else {
    assert(0);
  }

  assert(apos < alen);
  assert(bpos < blen);

  sanityCheck((afid >= 0) ? GetFragment(fragmentStore,afid)->firstbead : -1,
              (bfid >= 0) ? GetFragment(fragmentStore,bfid)->firstbead : -1);


  //  Negative ahang?  push these things onto the start of abacus.  Fail if we get a negative ahang,
  //  but there is a column already before the first one (which implies we aligned to something not
  //  full-length, possibly because we trimmed frankenstein wrong).
  //
  if (ahang < 0) {
    Bead  *bead = GetBead(beadStore, aindex[0]);
    int32  colp = bead->column_index;

    //GetBead(beadStore, aindex[apos])->column_index

    assert(bead->prev == -1);

    while (bpos < -ahang) {
      bead = GetBead(beadStore, bindex[bpos]);

#ifdef DEBUG_ABACUS_ALIGN
      fprintf(stderr, "ApplyAlignment()-- Prepend column for ahang bead=%d,%c\n",
              bead->boffset,
              *Getchar(sequenceStore, bead->soffset));
#endif
      ColumnPrepend(colp, bindex[bpos++]);
    }

    lasta = GetBead(beadStore, aindex[0])->prev;
    lastb = bindex[bpos - 1];
  }


  //  trace is NULL for the first fragment, all we want to do there is load the initial abacus.
  //
  if (trace == NULL)
    goto loadInitial;


  while ((trace != NULL) && (*trace != 0)) {

#ifdef DEBUG_ABACUS_ALIGN
    fprintf(stderr, "trace=%d  apos=%d alen=%d bpos=%d blen=%d\n", *trace, apos, alen, bpos, blen);
#endif

    //
    //  Gap is in afrag.  align ( - *trace - apos ) positions
    //
    if ( *trace < 0 ) {
      while ( apos < (- *trace - 1))
        alignPosition(aindex, apos, alen, bindex, bpos, blen, lasta, lastb, "ApplyAlignment(1)");

      //  Insert a gap column to accommodate the insert in B.

      assert(apos < alen);
      assert(bpos < blen);

      //  Hmmm.  Occasionally we don't do alignPositions() above on the first thing -- the alignment
      //  starts with a gap??!

      if ((lasta == -1) || (bpos == 0)) {
        assert(lasta == -1);
        assert(bpos  == 0);
        ColumnPrepend(GetBead(beadStore, aindex[apos])->column_index,
                      bindex[bpos]);
        lasta = GetBead(beadStore, aindex[apos])->prev;
        lastb = bindex[bpos];
      } else {
        assert(lasta == GetBead(beadStore, aindex[apos])->prev);
        ColumnAppend(GetBead(beadStore, lasta)->column_index,
                     bindex[bpos]);
        lasta = GetBead(beadStore, lasta)->next;
        lastb = bindex[bpos];
      }

      //  lasta should be 'aindex[apos]->prev'.

      assert(lasta == GetBead(beadStore, aindex[apos])->prev);

      bpos++;
    }

    //
    //  Gap is in bfrag.  align ( *trace - bpos ) positions
    //
    if (*trace > 0) {
      while ( bpos < (*trace - 1) )
        alignPosition(aindex, apos, alen, bindex, bpos, blen, lasta, lastb, "ApplyAlignment(4)");

      //  Insert a gap bead into B, and align to the existing column.  This event is like
      //  alignPositions(), and we must continue aligning to existing gap columns in A.

      assert(apos < alen);
      assert(bpos < blen);

      lasta = GetBead(beadStore, lasta)->next;
      lastb = AppendGapBead(lastb);

      assert(lasta = aindex[apos]);

      AlignBeadToColumn(GetBead(beadStore, lasta)->column_index,
                        lastb,
                        "ApplyAlignment(6)");

      apos++;

      //  Continue aligning to existing gap columns in A.  Duplication from alignPosition.
      //
      alignGaps(aindex, apos, alen, bindex, bpos, blen, lasta, lastb);
    }

    trace++;
  }

  //
  //  Remaining alignment contains no indels
  //
 loadInitial:
#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "Align the remaining:  bpos=%d blen=%d apos=%d alen=%d\n",
          bpos, blen, apos, alen);
#endif

  for (int32 rem = MIN(blen - bpos, alen - apos); rem > 0; rem--)
    alignPosition(aindex, apos, alen, bindex, bpos, blen, lasta, lastb, "ApplyAlignment(8)");

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "Append new sequence:  bpos=%d blen=%d apos=%d alen=%d\n",
          bpos, blen, apos, alen);
#endif

  if (blen - bpos > 0) {
    int32 ci = GetBead(beadStore, lastb)->column_index;
    int32 ng = 0;

    //  There shouldn't be any gaps left to insert...but there might be if we stop early above.  Old
    //  versions of this would simply stuff gaps into the B sequence for the rest of the existing
    //  multialign.  We'd prefer to fail.
    //
    for (Column *col = GetColumn(columnStore, ci); col->next != -1; col=GetColumn(columnStore, col->next))
      fprintf(stderr, "ERROR!  Column ci=%d has a next pointer (%d)\n", ci, col->next);
    assert(GetColumn(columnStore, ci)->next == -1);

    //  Add on trailing (dovetail) beads from b
    //
    for (int32 rem=blen-bpos; rem > 0; rem--)
      ci = ColumnAppend(ci, bindex[bpos++]);
  }

  //CheckColumns();

  if (afid >= 0)  safe_free(aindex);
  if (bfid >= 0)  safe_free(bindex);

  //bfrag->manode = afrag->manode;
}

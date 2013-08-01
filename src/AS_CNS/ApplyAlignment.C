
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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"
#include "MicroHetREZ.H"
#include "AS_UTL_reverseComplement.H"


#undef  DEBUG_FIND_BEAD
#undef  DEBUG_ALIGN_GAPS
#undef  DEBUG_ALIGN_POSITION
#undef  DEBUG_ABACUS_ALIGN

//  Add a column before cid, seeded with bead bid.
//
static
int32
ColumnPrepend(int32 cid, beadIdx bid) {

  ColumnBeadIterator ci;
  beadIdx  nid;

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

  if (call->prev.isValid())
    GetBead(beadStore,call->prev)->next = call->boffset;

  CreateColumnBeadIterator(cid, &ci);

  while ((nid = NextColumnBead(&ci)).isValid()) {
    bead = GetBead(beadStore, nid);

    //  The original version would not insert a gap bead at the start of a fragment.

    if ((bead->prev.isValid()) &&
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
  int32 fails = 0;
  int32 cmax  = GetNumColumns(columnStore);

  for (int32 cid=0; cid<cmax; cid++) {
    Column *column = GetColumn(columnStore,cid);

    if (column == NULL)
      continue;

    ColumnBeadIterator ci;
    CreateColumnBeadIterator(cid,&ci);

    Bead   *call = GetBead(beadStore,column->call);
    beadIdx  bid;

    while ( (bid = NextColumnBead(&ci)).isValid() ) {
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
sanityCheck(beadIdx aboffset, beadIdx bboffset) {
  Bead   *b = GetBead(beadStore, aboffset);
  int32   e = 0;

  while (b) {
    if (b->column_index == -1) {
      e++;
      fprintf(stderr, "bead "F_U64" in A has undef column_index.\n",
              (uint64)b->boffset.get());
    }
    b = (b->next.isInvalid()) ? NULL : GetBead(beadStore, b->next);
  }

  b = GetBead(beadStore, bboffset);

  while (b) {
    if (b->column_index != -1) {
      e++;
      fprintf(stderr, "bead "F_U64" in B has defined column_index %d.\n",
              (uint64)b->boffset.get(), b->column_index);
    }
    b = (b->next.isInvalid()) ? NULL : GetBead(beadStore, b->next);
  }

  assert(e == 0);
}




//  Returns the boffset of the bead that is:
//    in the same column as bead bi
//    in the same fragment as bead fi
//
static
beadIdx
findBeadInColumn(beadIdx bi, beadIdx fi) {

#ifdef DEBUG_FIND_BEAD
  fprintf(stderr, "findBeadInColumn bead bi=%d bead fi=%d\n", bi.get(), fi.get());
#endif

  if ((fi.isInvalid()) || (bi.isInvalid()))
    return(beadIdx());

  int32 ff = GetBead(beadStore, fi)->frag_index;

#ifdef DEBUG_FIND_BEAD
  fprintf(stderr, "findBeadInColumn fragindex ff="F_S32"\n", ff);
#endif

  Bead *b = GetBead(beadStore, bi);

  if (b->frag_index == ff)
    return(bi);

  //  Search up.  The way we call findBeadInColumn, the one we're
  //  looking for is usually up from where we start.
  while (b->up.isValid()) {
    b = GetBead(beadStore, b->up);
#ifdef DEBUG_FIND_BEAD
    fprintf(stderr, "findBeadInColumn up bead="F_U64" ff=%d\n", (uint64)b->boffset.get(), b->frag_index);
#endif
    if (b->frag_index == ff)
      return(b->boffset);
  }

  //  Search down.
  b = GetBead(beadStore, bi);

  while (b->down.isValid()) {
    b = GetBead(beadStore, b->down);
#ifdef DEBUG_FIND_BEAD
    fprintf(stderr, "findBeadInColumn down bead="F_U64" ff=%d\n", (uint64)b->boffset.get(), b->frag_index);
#endif
    if (b->frag_index == ff)
      return(b->boffset);
  }

  //  Give up.  See comments in MultiAlignUnitig, where we append new sequence to the start of
  //  frankenstein ("Append the new stuff...").

  assert(0);
  return(beadIdx());
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
alignGaps(beadIdx *aindex, int32 &apos, int32  alen,
          beadIdx *bindex, int32 &bpos, int32  blen,
          beadIdx &lasta,
          beadIdx &lastb) {

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

  //#ifdef DEBUG_ALIGN_GAPS
  //  fprintf(stderr, "alignGaps()-- apos=%d alen=%d  bpos=%d blen=%d  lasta=%d  lastb=%d\n",
  //          apos, alen, bpos, blen, lasta.get(), lastb.get());
  //#endif

  if (apos >= alen)
    return;

  lasta = findBeadInColumn(lasta, aindex[apos]);

  if (lasta.isInvalid())
    return;

  beadIdx nexta = GetBead(beadStore, lasta)->next;

  //  If the next bead (from 'x' in the picture) is NOT the next base in A, add gaps until that is
  //  so.

  while (nexta != aindex[apos]) {
#ifdef DEBUG_ALIGN_GAPS
    fprintf(stderr, "alignGaps()-- lasta=%d  nexta=%d  search for aindex[apos]=%d  \n",
            lasta.get(), nexta.get(), aindex[apos].get());
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
alignPosition(beadIdx *aindex, int32 &apos, int32  alen,
              beadIdx *bindex, int32 &bpos, int32  blen,
              beadIdx &lasta,
              beadIdx &lastb,
              char   *label) {

  assert(apos < alen);
  assert(bpos < blen);

  Bead *bead = GetBead(beadStore, aindex[apos]);

#ifdef DEBUG_ALIGN_POSITION
  fprintf(stderr, "alignPosition()-- add %c to column %d apos=%d bpos=%d lasta=%d lastb=%d\n",
          *Getchar(sequenceStore, bead->soffset), bead->column_index, apos, bpos, lasta.get(), lastb.get());
#endif

  AlignBeadToColumn(bead->column_index, bindex[bpos], label);

  lasta = aindex[apos];
  lastb = bindex[bpos];

  apos++;
  bpos++;

  alignGaps(aindex, apos, alen, bindex, bpos, blen, lasta, lastb);

  //assert(GetBead(beadStore, aindex[apos])->prev == bead->bindex);
  //assert(GetBead(beadStore, bindex[bpos])->prev == lastb);
}




void
ApplyAlignment(int32 afid,
               int32 alen, beadIdx *aindex,
               int32 bfid,
               int32 ahang,
               int32 *trace) {

  int32   bpos   = 0;
  int32   blen   = 0;

  int32   apos   = MAX(ahang,0);
  beadIdx *bindex = NULL;

  beadIdx  lasta;
  beadIdx  lastb;

  if (afid >= 0) {
    Fragment *afrag = GetFragment(fragmentStore,afid);
    alen            = afrag->length;
    aindex          = (beadIdx *)safe_malloc(alen * sizeof(beadIdx));

    //  The loop really abuses the fact that all the beads for the bases in this read are
    //  contiguous.  They're contiguous because CreateMANode() (I think it's in there) allocated
    //  them in one block.
    //
    for (int32 ai=0;ai<alen;ai++)
      aindex[ai].set(afrag->firstbead.get() + ai);
  } else {
    assert(aindex != NULL);
  }


  if (bfid >= 0) {
    Fragment *bfrag = GetFragment(fragmentStore,bfid);
    blen            = bfrag->length;
    bindex          = (beadIdx *)safe_malloc(blen * sizeof(beadIdx));

    for (int32 bi=0;bi<blen;bi++)
      bindex[bi].set(bfrag->firstbead.get() + bi);

    //  USED?
    bfrag->manode = NULL;
  } else {
    assert(0);
  }

  assert(apos < alen);
  assert(bpos < blen);

  sanityCheck((afid >= 0) ? GetFragment(fragmentStore,afid)->firstbead : beadIdx(),
              (bfid >= 0) ? GetFragment(fragmentStore,bfid)->firstbead : beadIdx());


  //  Negative ahang?  push these things onto the start of abacus.  Fail if we get a negative ahang,
  //  but there is a column already before the first one (which implies we aligned to something not
  //  full-length, possibly because we trimmed frankenstein wrong).
  //
  if (ahang < 0) {
    Bead  *bead = GetBead(beadStore, aindex[0]);
    int32  colp = bead->column_index;

    //GetBead(beadStore, aindex[apos])->column_index

    assert(bead->prev.isInvalid());

    while (bpos < -ahang) {
      bead = GetBead(beadStore, bindex[bpos]);

#ifdef DEBUG_ABACUS_ALIGN
      fprintf(stderr, "ApplyAlignment()-- Prepend column for ahang bead=%d,%c\n",
              bead->boffset.get(),
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
      //  starts with a gap??!  See below for an example.

      if ((lasta.isInvalid()) || (bpos == 0)) {
        assert(lasta.isInvalid());
        assert(bpos  == 0);
        ColumnPrepend(GetBead(beadStore, aindex[apos])->column_index, bindex[bpos]);
        lasta = GetBead(beadStore, aindex[apos])->prev;
        lastb = bindex[bpos];
      } else {
        assert(lasta == GetBead(beadStore, aindex[apos])->prev);
        ColumnAppend(GetBead(beadStore, lasta)->column_index, bindex[bpos]);
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

      //  Hmmm.  Occasionally we don't do alignPositions() above on the first thing -- the alignment
      //  starts with a gap??!  Sure:
      //
      //   ttaaaat.....
      //  n-taaaat.....
      //
      //  The negative ahang triggers code above, which adds a new
      //  column.  lasta is still invalid because there is no last column
      //  for a.

      lasta = (lasta.isInvalid()) ? aindex[apos] : GetBead(beadStore, lasta)->next;
      lastb = AppendGapBead(lastb);

      assert(lasta == aindex[apos]);

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

  //
  //  Now just tack on the new sequence
  //

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
      fprintf(stderr, "ERROR!  Column ci=%d has a next pointer (%d)\n", col->lid, col->next);
#warning assert skipped until contig consensus gets fixed
    //assert(GetColumn(columnStore, ci)->next == -1);

    //  Add on trailing (dovetail) beads from b
    //
    for (int32 rem=blen-bpos; rem > 0; rem--) {
#ifdef DEBUG_ALIGN_POSITION
      Bead *bead = GetBead(beadStore, bindex[bpos]);
      fprintf(stderr, "alignPosition()-- add %c to column %d\n",
              *Getchar(sequenceStore, bead->soffset), bead->column_index);
#endif
      ci = ColumnAppend(ci, bindex[bpos++]);
    }
  }

  //CheckColumns();

  if (afid >= 0)  safe_free(aindex);
  if (bfid >= 0)  safe_free(bindex);

  //bfrag->manode = afrag->manode;
}

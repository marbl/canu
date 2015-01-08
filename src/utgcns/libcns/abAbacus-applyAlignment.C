
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

#undef  DEBUG_FIND_BEAD
#undef  DEBUG_ALIGN_GAPS
#undef  DEBUG_ALIGN_POSITION
#undef  DEBUG_ABACUS_ALIGN

//  Add a column before cid, seeded with bead bid.
//
abColID
abAbacus::prependColumn(abColID cid, abBeadID bid) {
  abColumn *next     = getColumn(cid);
  abBead   *nextcall = getBead(next->call);

  abColID   col      = addColumn(next->ma_id, bid);
  abColumn *column   = getColumn(col);
  abBead   *call     = getBead(column->call);

  //fprintf(stderr, "prependColumn()-- adding column for bid %d\n", bid);

  column->prev   = next->prev;
  column->next   = cid;

  call->prev     = nextcall->prev;
  call->next     = nextcall->boffset;
  next->prev     = column->lid;

  nextcall->prev = call->boffset;

  if (column->prevID().isValid())
    getColumn(column->prev)->next = column->lid;

  if (call->prev.isValid())
    getBead(call->prev)->next = call->boffset;

  abColBeadIterator *ci = createColBeadIterator(cid);

  for (abBeadID  nid=ci->next(); nid.isValid(); nid=ci->next()) {
    abBead *bead = getBead(nid);

    //  The original version would not insert a gap bead at the start of a fragment.

    if ((bead->prev.isValid()) &&
        (bead->prev != bid))
      alignBeadToColumn(column->lid, prependGapBead(nid), "prependColumn()");
  }

  column->ma_id       = next->ma_id;
  column->ma_position = next->ma_position - 1;

  //AddColumnToMANode(column->ma_id, *column);
  {
    abMultiAlign *ma = getMultiAlign(column->ma_id);

    ma->columnList.push_back(column->lid);

    if (column->next.isValid() == false)
      ma->last = column->lid;

    if (column->prev.isValid() == false)
      ma->first = column->lid;
  }

  //  Redundant
  //if (column->prevID().isValid() == false)
  //  getMultiAlign(column->ma_id)->first = column->lid;

  return(column->lid);
}









//  Returns the boffset of the bead that is:
//    in the same column as bead bi
//    in the same fragment as bead fi
//
abBeadID
findBeadInColumn(abAbacus *abacus, abBeadID bi, abBeadID fi) {

#ifdef DEBUG_FIND_BEAD
  fprintf(stderr, "findBeadInColumn bead bi="F_U32" bead fi="F_U32"\n", bi.get(), fi.get());
#endif

  if ((fi.isInvalid()) || (bi.isInvalid()))
    return(abBeadID());

  abSeqID ff = abacus->getBead(fi)->seqIdx();

#ifdef DEBUG_FIND_BEAD
  fprintf(stderr, "findBeadInColumn fragindex ff="F_S32"\n", ff.get());
#endif

  abBead *b = abacus->getBead(bi);

  if (b->seqIdx() == ff)
    return(bi);

  //  Search up.  The way we call findBeadInColumn, the one we're
  //  looking for is usually up from where we start.
  while (b->upID().isValid()) {
    b = abacus->getBead(b->upID());
#ifdef DEBUG_FIND_BEAD
    fprintf(stderr, "findBeadInColumn up bead="F_U32" ff=%d\n", b->boffset.get(), b->seqIdx().get());
#endif
    if (b->seqIdx() == ff)
      return(b->ident());
  }

  //  Search down.
  b = abacus->getBead(bi);

  while (b->downID().isValid()) {
    b = abacus->getBead(b->downID());
#ifdef DEBUG_FIND_BEAD
    fprintf(stderr, "findBeadInColumn down bead="F_U64" ff=%d\n", b->boffset.get(), b->seqIdx());
#endif
    if (b->seqIdx() == ff)
      return(b->ident());
  }

  //  Give up.  See comments in MultiAlignUnitig, where we append new sequence to the start of
  //  frankenstein ("Append the new stuff...").

  assert(0);
  return(abBeadID());
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
alignGaps(abAbacus *abacus,
          abBeadID *aindex, int32 &apos, int32  alen,
          abBeadID *bindex, int32 &bpos, int32  blen,
          abBeadID &lasta,
          abBeadID &lastb) {

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

  lasta = findBeadInColumn(abacus, lasta, aindex[apos]);

  if (lasta.isInvalid())
    return;

  abBeadID nexta = abacus->getBead(lasta)->nextID();

  //  If the next bead (from 'x' in the picture) is NOT the next base in A, add gaps until that is
  //  so.

  while (nexta != aindex[apos]) {
#ifdef DEBUG_ALIGN_GAPS
    fprintf(stderr, "alignGaps()-- lasta=%d  nexta=%d  search for aindex[apos]=%d  \n",
            lasta.get(), nexta.get(), aindex[apos].get());
#endif

    lastb = abacus->appendGapBead(lastb);
    lasta = nexta;

    abBead *bead = abacus->getBead(nexta);  //  AppendGapBead might reallocate beadStore

    abacus->alignBeadToColumn(bead->colIdx(), lastb, "applyAlignment(alignGaps)");

    nexta = bead->nextID();
  }
}


//  Aligns a bead in B to the existing column in A.
static
void
alignPosition(abAbacus *abacus,
              abBeadID *aindex, int32 &apos, int32  alen,
              abBeadID *bindex, int32 &bpos, int32  blen,
              abBeadID &lasta,
              abBeadID &lastb,
              char   *label) {

  assert(apos < alen);
  assert(bpos < blen);

  abBead *bead = abacus->getBead(aindex[apos]);

#ifdef DEBUG_ALIGN_POSITION
  fprintf(stderr, "alignPosition()-- add %c to column %d apos=%d bpos=%d lasta=%d lastb=%d\n",
          *Getchar(sequenceStore, bead->soffset), bead->column_index, apos, bpos, lasta.get(), lastb.get());
#endif

  abacus->alignBeadToColumn(bead->colIdx(), bindex[bpos], label);

  lasta = aindex[apos];
  lastb = bindex[bpos];

  apos++;
  bpos++;

  alignGaps(abacus, aindex, apos, alen, bindex, bpos, blen, lasta, lastb);

  //assert(abacus->getBead(aindex[apos])->prev == bead->bindex);
  //assert(abacus->getBead(bindex[bpos])->prev == lastb);
}




void
abAbacus::applyAlignment(abSeqID   afid,
                         int32     alen, abBeadID *aindex,
                         abSeqID   bfid,
                         int32     ahang,
                         int32    *trace) {

  int32     bpos   = 0;
  int32     blen   = 0;

  int32     apos   = MAX(ahang,0);

  abBeadID *bindex = NULL;

  abBeadID  lasta;
  abBeadID  lastb;

  if (afid.isValid()) {
    abSequence *afrag = getSequence(afid);
    alen              = afrag->length;
    aindex            = new abBeadID [alen];

    //  The loop really abuses the fact that all the beads for the bases in this read are
    //  contiguous.  They're contiguous because CreateMANode() (I think it's in there) allocated
    //  them in one block.
    //
    for (uint32 ai=0; ai<alen; ai++)
      aindex[ai].set(afrag->firstbead.get() + ai);

  } else {
    assert(aindex != NULL);
  }


  if (bfid.isValid()) {
    abSequence *bfrag = getSequence(bfid);
    blen              = bfrag->length;
    bindex            = (abBeadID *)safe_malloc(blen * sizeof(abBeadID));

    for (uint32 bi=0; bi<blen; bi++)
      bindex[bi].set(bfrag->firstbead.get() + bi);

    bfrag->manode = abMultiAlignID();  //  USED?
  } else {
    assert(0);
  }

  assert(apos < alen);
  assert(bpos < blen);

  //
  //  All the a beads should be in a column.  All the b beads should not.
  //
  //sanityCheck((afid >= 0) ? getFragment(afid)->firstbead : abBeadID(),
  //            (bfid >= 0) ? getFragment(bfid)->firstbead : abBeadID());

  {
    uint32    columnErrors = 0;

    if (afid.isValid()) {
      abBead   *b = getBead(getSequence(afid)->firstbead);

      while (b) {
        if (b->column_index.isValid() == false) {
          columnErrors++;
          fprintf(stderr, "bead "F_U32" in A has undef column_index.\n",
                  b->ident().get());
        }
        b = (b->nextID().isInvalid()) ? NULL : getBead(b->nextID());
      }
    }

    if (bfid.isValid()) {
      abBead *b = getBead(getSequence(bfid)->firstbead);

      while (b) {
        if (b->column_index.isValid() == true) {
          columnErrors++;
          fprintf(stderr, "bead "F_U32" in B has defined column_index %d.\n",
                  b->ident().get(), b->column_index.get());
        }
        b = (b->nextID().isInvalid()) ? NULL : getBead(b->nextID());
      }
    }

    assert(columnErrors == 0);
  }



  //  Negative ahang?  push these things onto the start of abacus.  Fail if we get a negative ahang,
  //  but there is a column already before the first one (which implies we aligned to something not
  //  full-length, possibly because we trimmed frankenstein wrong).
  //
  if (ahang < 0) {
    abBead  *bead = getBead(aindex[0]);
    abColID  colp = bead->colIdx();

    //abacus->getBead(aindex[apos])->column_index

    assert(bead->prev.isInvalid());

    while (bpos < -ahang) {
      bead = getBead(bindex[bpos]);

#ifdef DEBUG_ABACUS_ALIGN
      fprintf(stderr, "ApplyAlignment()-- Prepend column for ahang bead=%d,%c\n",
              bead->boffset.get(),
              getBase(bead->baseIdx()));
#endif
      prependColumn(colp, bindex[bpos++]);
    }

    lasta = getBead(aindex[0])->prev;
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
        alignPosition(this, aindex, apos, alen, bindex, bpos, blen, lasta, lastb, "ApplyAlignment(1)");

      //  Insert a gap column to accommodate the insert in B.

      assert(apos < alen);
      assert(bpos < blen);

      //  Hmmm.  Occasionally we don't do alignPositions() above on the first thing -- the alignment
      //  starts with a gap??!  See below for an example.

      if ((lasta.isInvalid()) || (bpos == 0)) {
        assert(lasta.isInvalid());
        assert(bpos  == 0);

        prependColumn(getBead(aindex[apos])->column_index, bindex[bpos]);

        lasta = getBead(aindex[apos])->prev;
        lastb = bindex[bpos];
      } else {
        assert(lasta == getBead(aindex[apos])->prev);
        appendColumn(getBead(lasta)->colIdx(), bindex[bpos]);
        lasta = getBead(lasta)->nextID();
        lastb = bindex[bpos];
      }

      //  lasta should be 'aindex[apos]->prev'.

      assert(lasta == getBead(aindex[apos])->prev);

      bpos++;
    }

    //
    //  Gap is in bfrag.  align ( *trace - bpos ) positions
    //
    if (*trace > 0) {
      while ( bpos < (*trace - 1) )
        alignPosition(this, aindex, apos, alen, bindex, bpos, blen, lasta, lastb, "ApplyAlignment(4)");

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

      lasta = (lasta.isInvalid()) ? aindex[apos] : getBead(lasta)->nextID();
      lastb = appendGapBead(lastb);

      assert(lasta == aindex[apos]);

      alignBeadToColumn(getBead(lasta)->colIdx(),
                        lastb,
                        "ApplyAlignment(6)");

      apos++;

      //  Continue aligning to existing gap columns in A.  Duplication from alignPosition.
      //
      alignGaps(this, aindex, apos, alen, bindex, bpos, blen, lasta, lastb);
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
    alignPosition(this, aindex, apos, alen, bindex, bpos, blen, lasta, lastb, "ApplyAlignment(8)");

  //
  //  Now just tack on the new sequence
  //

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "Append new sequence:  bpos=%d blen=%d apos=%d alen=%d\n",
          bpos, blen, apos, alen);
#endif

  if (blen - bpos > 0) {
    abColID ci = getBead(lastb)->colIdx();

    //  There shouldn't be any gaps left to insert...but there might be if we stop early above.  Old
    //  versions of this would simply stuff gaps into the B sequence for the rest of the existing
    //  multialign.  We'd prefer to fail.
    //
    for (abColumn *col = getColumn(ci); col->nextID().isValid(); col=getColumn(col->nextID()))
      fprintf(stderr, "ERROR!  Column ci="F_U32" has a next pointer ("F_U32")\n",
              col->ident().get(), col->nextID().get());
#warning assert skipped until contig consensus gets fixed
    //assert(getColumn(ci)->nextID() == -1);

    //  Add on trailing (dovetail) beads from b
    //
    for (int32 rem=blen-bpos; rem > 0; rem--) {
#ifdef DEBUG_ALIGN_POSITION
      abBead *bead = getBead(bindex[bpos]);
      fprintf(stderr, "alignPosition()-- add %c to column %d\n",
              getBase(bead->baseIdx()), bead->colIdx());
#endif
      ci = appendColumn(ci, bindex[bpos++]);
    }
  }

  //CheckColumns();
#if 0
  {
    int32 fails = 0;
    int32 cmax  = GetNumColumns(columnStore);

    for (int32 cid=0; cid<cmax; cid++) {
      abColumn *column = getColumn(cid);

      if (column == NULL)
        continue;

      ColumnBeadIterator ci;
      CreateColumnBeadIterator(cid,&ci);

      Bead   *call = GetBead(beadStore,column->call);
      abBeadID  bid;

      while ( (bid = NextColumnBead(&ci)).isValid() ) {
        Bead *bead = GetBead(beadStore,bid);

        if (bead->colIdx() != cid)
          fails++;
      }
    }
    assert(fails == 0);
  }
#endif

  if (afid.isValid())  delete [] aindex;
  if (bfid.isValid())  delete [] bindex;

  //bfrag->manode = afrag->manode;
}

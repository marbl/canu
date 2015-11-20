
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
 *    src/AS_CNS/ApplyAlignment.C
 *    src/AS_CNS/ApplyAlignment.c
 *    src/AS_CNS/MultiAlignment_CNS.c
 *    src/utgcns/libcns/ApplyAlignment.C
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
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JUL-28
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-13
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static char *rcsid = "$Id$";

#include "abAbacus.H"

#undef  DEBUG_FIND_BEAD       //  Was useful for tracking down matrix structure issues
#undef  DEBUG_ALIGN_GAPS      //  Reports calls to alignGaps()
#undef  DEBUG_ALIGN_POSITION  //  Shows when and where an unaligned bead is added to the multialignment
#undef  DEBUG_ABACUS_ALIGN    //  Primary debug output, shows progress of the algorithm

#undef  TEST_ABACUS_ALIGN     //  Check that all the beads are properly in columns.  Expensive.


//  Add a column before cid, seeded with bead bid.
//
abColID
abAbacus::prependColumn(abColID cid, abBeadID bid) {
  abColID   col      = addColumn(getColumn(cid)->ma_id, bid);  //  Can reallocate column list.

  abColumn *column   = getColumn(col);
  abColumn *next     = getColumn(cid);

  abBead   *call     = getBead(column->callID());
  abBead   *nextcall = getBead(next->callID());

  //fprintf(stderr, "prependColumn()-- adding column for bid %d\n", bid);

  column->prev   = next->prev;
  column->next   = cid;

  call->prev     = nextcall->prev;
  call->next     = nextcall->ident();
  next->prev     = column->ident();

  nextcall->prev = call->ident();

  if (column->prevID().isValid())
    getColumn(column->prev)->next = column->ident();

  if (call->prev.isValid())
    getBead(call->prev)->next = call->ident();

  abColBeadIterator *ci = createColBeadIterator(cid);

  for (abBeadID  nid=ci->next(); nid.isValid(); nid=ci->next()) {
    abBead *bead = getBead(nid);

    //  The original version would not insert a gap bead at the start of a fragment.

    if ((bead->prev.isValid()) &&
        (bead->prev != bid))
      alignBeadToColumn(column->ident(), prependGapBead(nid), "prependColumn()");
  }

  column->ma_position = next->ma_position - 1;

  //AddColumnToMANode(column->ma_id, *column);
  getMultiAlign(column->ma_id)->addColumnToMultiAlign(column);

  //  Redundant
  //if (column->prevID().isValid() == false)
  //  getMultiAlign(column->ma_id)->first = column->lid;

  return(column->ident());
}





//  Returns the ident() of the bead that is:
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
    fprintf(stderr, "findBeadInColumn up bead="F_U32" ff=%d\n", b->ident().get(), b->seqIdx().get());
#endif
    if (b->seqIdx() == ff)
      return(b->ident());
  }

  //  Search down.
  b = abacus->getBead(bi);

  while (b->downID().isValid()) {
    b = abacus->getBead(b->downID());
#ifdef DEBUG_FIND_BEAD
    fprintf(stderr, "findBeadInColumn down bead="F_U64" ff=%d\n", b->ident().get(), b->seqIdx().get());
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
//
//  Called by applyAlignment().
//
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
  fprintf(stderr, "alignPosition()-- add %c to column %d (prev=%d next=%d) apos=%d bpos=%d lasta=%d lastb=%d\n",
          abacus->getBase(bead->baseIdx()),
          bead->colIdx().get(),
          abacus->getColumn(bead->colIdx())->prevID().get(),
          abacus->getColumn(bead->colIdx())->nextID().get(),
          apos, bpos,
          lasta.get(), lastb.get());
#endif

  abacus->alignBeadToColumn(bead->colIdx(), bindex[bpos], label);

  lasta = aindex[apos];
  lastb = bindex[bpos];

  apos++;
  bpos++;

  alignGaps(abacus, aindex, apos, alen, lasta, lastb);

  //assert(abacus->getBead(aindex[apos])->prev == bead->bindex);
  //assert(abacus->getBead(bindex[bpos])->prev == lastb);
}




void
abAbacus::applyAlignment(abSeqID   afid,
                         int32     alen, abBeadID *aindex,
                         abSeqID   bfid,
                         int32     ahang,
                         int32     bhang,
                         int32    *trace, uint32 traceLen) {

  int32     apos   = MAX(ahang,0);  //  The aindex index for the current column we're processing
  int32     bpos   = bhang;

  int32     blen   = 0;

  abBeadID *bindex = NULL;

  assert(trace != NULL);  //  The original (Celera) version of this could be called without a trace

  assert(afid.isValid() == false);  //  The original (CA8) version of this could be called with an actual read here
  assert(bfid.isValid() == true);   //  But could never be called without a read here

#if 0
  for (uint32 ii=0; trace[ii]; ii++)
    fprintf(stderr, "trace[%u] %d\n", ii, trace[ii]);
#endif

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "abAbacus::applyAlignment()-- ahang=%d bhang=%d traceLen=%u trace=%d %d %d %d %d ... %d %d %d %d %d\n",
          ahang, bhang, traceLen,
          (traceLen > 0) ? trace[0] : 0,
          (traceLen > 1) ? trace[1] : 0,
          (traceLen > 2) ? trace[2] : 0,
          (traceLen > 3) ? trace[3] : 0,
          (traceLen > 4) ? trace[4] : 0,
          (traceLen > 4) ? trace[traceLen-5] : 0,
          (traceLen > 3) ? trace[traceLen-4] : 0,
          (traceLen > 2) ? trace[traceLen-3] : 0,
          (traceLen > 1) ? trace[traceLen-2] : 0,
          (traceLen > 0) ? trace[traceLen-1] : 0);
#endif

  //  These loops really abuses the fact that all the beads for the bases in this read are
  //  contiguous.  They're contiguous because CreateMANode() (I think it's in there) allocated them
  //  in one block.

  if (afid.isValid()) {
    abSequence *afrag = getSequence(afid);
    alen              = afrag->length();
    aindex            = new abBeadID [alen];

    for (uint32 ai=0; ai<alen; ai++)
      aindex[ai].set(afrag->firstBead().get() + ai);

  } else {
    assert(aindex != NULL);
  }


  if (bfid.isValid()) {
    abSequence *bfrag = getSequence(bfid);
    blen              = bfrag->length();
    bindex            = new abBeadID [blen];

    for (uint32 bi=0; bi<blen; bi++)
      bindex[bi].set(bfrag->firstBead().get() + bi);

  } else {
    assert(0);
  }

  //  if apos == alen, we'd just be pasting on new sequence.
  //  if bpos == blen...we're pasting on one base?

  assert(apos <= alen);
  assert(bpos <= blen);

  //
  //  All the a beads should be in a column.  All the b beads should not.  This is expensive.
  //

#ifdef TEST_ABACUS_ALIGN
  {
    uint32    columnErrors = 0;

    if (afid.isValid()) {
      abBead   *b = getBead(getSequence(afid)->firstBead());

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
      abBead *b = getBead(getSequence(bfid)->firstBead());

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
#endif


  //  The beadIDs of the bases in the last column aligned.  Both generally get reset almost immediately,
  //  EXCEPT when the b sequence is appended to the end, as in the Quick variant of utgcns.

  abBeadID  lasta;
  abBeadID  lastb;

  if (apos > 0)
    lasta = aindex[apos - 1];

  if (bpos > 0)
    lastb = bindex[bpos - 1];

  //
  //  Catch two nasty cases: negative ahang (so unaligned bases in the B read) and
  //  initial gaps in the B read.
  //


  //  Negative ahang?  Push these things onto the start of abacus.  Fail if we get a negative ahang,
  //  but there is a column already before the first one (which implies we aligned to something not
  //  full-length, possibly because we trimmed frankenstein wrong).
  //
  if (ahang < 0) {
    abBead  *bead = getBead(aindex[0]);
    abColID  colp = bead->colIdx();

    assert(bead->prev.isInvalid());

    while (bpos < -ahang) {
#ifdef DEBUG_ABACUS_ALIGN
      fprintf(stderr, "ApplyAlignment()-- Prepend column for ahang bead=%d,%c\n",
              getBead(bindex[bpos])->ident().get(),
              getBase(getBead(bindex[bpos])->baseIdx()));
#endif

      prependColumn(colp, bindex[bpos++]);
    }

    //  Why are these set?

    lasta = getBead(aindex[0])->prev;  //  Probably still invalid, there isn't any A-read aligned yet.
    lastb = bindex[bpos - 1];          //  But this should be valid; it's pointing into the B-read
  }


  //  Skip any positive traces at the start.  These are template sequence (A-read) aligned to gaps
  //  in the B-read, but we don't care because there isn't any B-read aligned yet.

  while ((traceLen > 0) && (*trace != 0) && (*trace == 1)) {
#ifdef DEBUG_ABACUS_ALIGN
    fprintf(stderr, "trace=%d  apos=%d alen=%d bpos=%d blen=%d - SKIP INITIAL GAP IN READ\n", *trace, apos, alen, bpos, blen);
#endif

    //lasta = aindex[apos];  //  NEW, untested

    apos++;  //  lasta isn't set, it is reset to aindex[apos] in alignPosition().
    trace++;
  }

  //  Similarly, remove negative traces at the end.  These are read sequence (B-read) aligned to gaps
  //  in the A-read (template).  The A-read extent should be reduced by one for each  gap.

  while ((traceLen > 0) && (trace[traceLen-1] > blen)) {
#ifdef DEBUG_ABACUS_ALIGN
    fprintf(stderr, "trace=%d  apos=%d alen=%d bpos=%d blen=%d - SKIP TERMINAL GAP IN READ\n", *trace, apos, alen, bpos, blen);
#endif

    trace[--traceLen] = 0;
  }



  while ((traceLen > 0) && (*trace != 0)) {

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

      //  Handle an initial -1 trace.  This bypasses the alignPosition above, which results in an
      //  invalid lasta.

      if ((lasta.isInvalid()) || (bpos == 0)) {
        prependColumn(getBead(aindex[apos])->column_index, bindex[bpos]);

        lasta = getBead(aindex[apos])->prev;
        lastb = bindex[bpos];
      } else {
        assert(lasta == getBead(aindex[apos])->prev);

        appendColumn(getBead(lasta)->colIdx(), bindex[bpos]);

        lasta = getBead(lasta)->nextID();
        lastb = bindex[bpos];
      }

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

      //  Unlike the *trace < 0 case, we have already removed the initial +1 trace elements, so
      //  alignPosition() is always called, and lasta is always valid.
      //
      //  Well, that's what was supposed to happen.  It is possible for the read to come with
      //  unaligned bases at the start (so the 'prepend column' 'negative ahang' case is applied),
      //  then to have gaps inserted in the B read immediately after.  In effect, it looks like
      //  there are +1 trace elements that should have been removed, but aren't.
      //
      //  In one example, bpos=20 and *trace=21, and so alignPosition() was never called.

      //  guido

      assert(lastb.isValid());

      if ((lasta.isInvalid())) {
        lasta = aindex[apos];          //  Not actually 'last' until after the apos++ below
        lastb = appendGapBead(lastb);  //  Same.

      } else {
        lasta = getBead(lasta)->nextID();  //  Is now 'thisa', util apos++.
        lastb = appendGapBead(lastb);      //  Same.
      }

      assert(lasta == aindex[apos]);

      alignBeadToColumn(getBead(lasta)->colIdx(),
                        lastb,
                        "ApplyAlignment(6)");

      apos++;


      //  Continue aligning to existing gap columns in A.  Duplication from alignPosition.

      alignGaps(this, aindex, apos, alen, lasta, lastb);
    }

    trace++;
  }

  //
  //  Remaining alignment contains no indels
  //

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

    //  Find the last column aligned.  This used to be lastb, but that bead doesn't exist in the multialignment
    //  when we're appending reads to frankenstein with the -Q quick option.

    abColID ci = getBead(lasta)->colIdx();

    //  There shouldn't be any gaps left to insert...but there might be if we stop early above.  Old
    //  versions of this would simply stuff gaps into the B sequence for the rest of the existing
    //  multialign.  We'd prefer to fail.

    for (abColumn *col = getColumn(ci); col->nextID().isValid(); col=getColumn(col->nextID()))
      fprintf(stderr, "ERROR!  Column ci="F_U32" has a next pointer ("F_U32")\n",
              col->ident().get(), col->nextID().get());

    assert(getColumn(ci)->nextID().isValid() == false);

    //  Add on trailing (dovetail) beads from b

    for (int32 rem=blen-bpos; rem > 0; rem--) {

#ifdef DEBUG_ALIGN_POSITION
      abBead *bead = getBead(bindex[bpos]);
      fprintf(stderr, "alignPosition()-- add %c (bead %d seqIdx %d pos %d baseIdx %d) to column %d (prev=%d next=%d) rem=%d bpos=%d blen=%d\n",
              getBase(bead->baseIdx()),
              bead->ident().get(),
              bead->seqIdx().get(),
              bead->foffset,
              bead->baseIdx().get(),
              bead->colIdx().get(),
              getColumn(bead->colIdx())->prevID().get(),
              getColumn(bead->colIdx())->nextID().get(),
              rem, bpos, blen);
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

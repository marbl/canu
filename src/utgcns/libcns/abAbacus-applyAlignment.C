
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

#include "abAbacus.H"

#undef  DEBUG_INFER
#undef  DEBUG_ABACUS_ALIGN    //  Primary debug output, shows progress of the algorithm




void
abColumn::allocateInitialBeads(void) {

  //  Allocate beads.  We'll need no more than the max of either the prev or the next.  Any read that we
  //  interrupt gets a new gap bead.  Any read that has just ended gets nothing.  And, +1 for the read
  //  we might be adding to the multialign.

  uint32   pmax = (_prevColumn != NULL) ? (_prevColumn->depth() + 1) : (4);
  uint32   nmax = (_nextColumn != NULL) ? (_nextColumn->depth() + 1) : (4);

  _beadsMax = MAX(pmax, nmax);
  _beadsLen = 0;
  _beads    = new abBead [_beadsMax];

  for (uint32 ii=0; ii<_beadsMax; ii++)  //  Probably done by the constructor.
    _beads[ii].clear();
}





void
abColumn::inferPrevNextBeadPointers(void) {

#ifdef DEBUG_INFER
  fprintf(stderr, "infer between columns %d and %d\n", _prevColumn->position(), _nextColumn->position());
#endif

  //  Scan the next/prev lists, adding a gap bead for any read that is spanned.
  //
  //  prev/this/next       prev/this/next
  //   0    0    0 <----->  0    0    0      //  A spanned read.
  //   1    1   MAX   /-->  2    1    1      //  The first read ended
  //   2    2    1 <-/ /->  3    2    2      //  Now columns are different
  //   3    3    3 <--/

  assert(_prevColumn != NULL);
  assert(_nextColumn != NULL);

  //  At startup, the links in this column are all invalid, and the base hasn't been added.

  assert(_prevColumn->_beadsLen <= _prevColumn->_beadsMax);
  assert(             _beadsLen == 0);
  assert(_nextColumn->_beadsLen <= _nextColumn->_beadsMax);

  uint32   pmax = (_prevColumn == NULL) ? 0 : _prevColumn->_beadsLen;
  uint32   nmax = (_nextColumn == NULL) ? 0 : _nextColumn->_beadsLen;

  uint32   ppos = 0;  //  Position in the prev's nextOffset list
  uint32   tpos = 0;  //  Position in our offset list
  uint32   npos = 0;  //  Position in the next's prevOffset list

  //  Add gaps and link into any existing reads.

  while (1) {

    //  Find reads to work on.

    while ((ppos < pmax) &&                                           //  If the previous position is less than the depth, and
           (_prevColumn->_beads[ppos].nextOffset() == UINT16_MAX))    //     the previous pointer is not valid,
      ppos++;                                                         //  move ahead to the next read.

    while ((npos < nmax) &&
           (_nextColumn->_beads[npos].prevOffset() == UINT16_MAX))
      npos++;

    //  If either column ran out of entries, we're done.  No more reads will span this new column.

    if ((ppos >= _prevColumn->depth()) ||
        (npos >= _nextColumn->depth()))
      break;

    //  The columns must agree.

    uint16 no = _prevColumn->_beads[ppos].nextOffset();
    uint16 po = _nextColumn->_beads[npos].prevOffset();

    if ((no != npos) || (po != ppos)) {
      fprintf(stderr, "ERROR: link mismatch.\n");
      fprintf(stderr, "ERROR: prev %d %d -- next %d %d\n",
              _prevColumn->_beads[ppos].prevOffset(),
              _prevColumn->_beads[ppos].nextOffset(),
              _nextColumn->_beads[npos].prevOffset(),
              _nextColumn->_beads[npos].nextOffset());
    }
    assert(no == npos);  //  The 'next-bead' that the prev is pointing to must be npos.
    assert(po == ppos);  //  The 'prev-bead' that the next is pointing to must be ppos.

    //  Reset the prev/next offsets to point to us.

    _prevColumn->_beads[ppos]._nextOffset = tpos;
    _nextColumn->_beads[npos]._prevOffset = tpos;

    //  Set our offsets to point to them.

    _beads[tpos]._prevOffset = po;
    _beads[tpos]._nextOffset = no;

    _beadsLen = tpos + 1;  //  Needed for showLink()

#ifdef DEBUG_INFER
    fprintf(stderr, "inferPrevNextBeadPointers()\n");
    showLinks(ppos, tpos, npos);
#endif

    //  And move forward to the next read in all columns.

    ppos++;    assert(ppos <= _prevColumn->_beadsLen);
    tpos++;    assert(tpos <=              _beadsMax);
    npos++;    assert(npos <= _nextColumn->_beadsLen);
  }

  //  Set the new length of the bead list.

  assert(_beadsLen == tpos);
  _beadsLen = tpos;

#ifdef DEBUG_INFER
  fprintf(stderr, "inferPrevNextBeadPointers() final\n");
  showLinks();
#endif
}





//  Add a column, seeded with some base/qual assignment.  The assignment cannot be a gap (because
//  all other reads in this column are set to gap).
//
//  Returns the bead index for the read that was added.
//
//  Can insert the new column before or after some specified column.  In both cases, we need to
//  know the read index for the read in that column that we're adding a base to.
//
//  Common use cases are to:
//    add a new column for a base aligned to a gap in the multialign
//      - all existing reads get a gap
//      - this read is assigned a base
//
//    add a new column before/after the multialign
//      - this read is the only base in this column
//
//  In both cases, we need to have and return linking information to the existing read.
//
//  One odd ball case occurs when the read isn't in the muiltialign at all; when we're adding
//  the very first base.  This is handled by setting the link to max.
//



//  Insert a new empty column, seeded with only the base/qual supplied, to the start of the multialign.
//  This grows the multialign to add new bases to the start:
//        [original-multialign]
//       1[original-multialign]
//      12[original-multialign]
//     123[original-multialign]
//    1234[original-multialign]
//
uint16
abColumn::insertAtBegin(abColumn *first, uint16 prevLink, char base, uint8 qual) {

  //  The base CAN NOT be a gap - the new column would then be entirely a gap column, with no base.
  assert(base != '-');
  assert(base != 0);

  _columnPosition = INT32_MAX;

  _call = base;
  _qual = qual;

  _prevColumn = first->_prevColumn;
  _nextColumn = first;

  first->_prevColumn = this;

  if (_prevColumn)
    _prevColumn->_nextColumn = this;

  allocateInitialBeads();

  _beads[0]._unused     = 0;
  _beads[0]._isRead     = 1;
  _beads[0]._isUnitig   = 0;
  _beads[0]._base       = base;
  _beads[0]._qual       = qual;

  _beads[0]._prevOffset = prevLink;
  _beads[0]._nextOffset = UINT16_MAX;  //  No next offset for the next base.

  _beadsLen++;

  assert(_beadsLen <= _beadsMax);

  if (_prevColumn)
    _prevColumn->_beads[prevLink]._nextOffset = 0;

#ifdef BASECOUNT
  baseCountIncr(base);
#endif

  //  With only the one read, we don't need to do any base calling here.

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "insertAtBegin()-- column prev=%p next=%p  link prev=%d this=%d next=%d\n",
          _prevColumn, this, _nextColumn, _beads[0]._prevOffset, _beads[0]._nextOffset);
#endif

  return(0);
}



//  Insert a new empty column, seeded with only the base/qual supplied, to the end of the multialign.
//  The growth is simpler than insertAtBegin() because the end is always NULL.
//    [original-multialign]
//    [original-multialign]7
//    [original-multialign]78
//    [original-multialign]789
//
uint16
abColumn::insertAtEnd(abColumn *prev, uint16 prevLink, char base, uint8 qual) {

  assert(base != '-');    //  The base CAN NOT be a gap - the new column would then be entirely a gap column, with no base.
  assert(base != 0);

  _columnPosition = INT32_MAX;

  _call = base;
  _qual = qual;

  _prevColumn = prev;
  _nextColumn = NULL;

  if (prev == NULL)   assert(prevLink == UINT16_MAX);
  if (prev != NULL)   assert(prevLink != UINT16_MAX);

  if (prev)
    prev->_nextColumn = this;

  allocateInitialBeads();

  _beads[0]._unused     = 0;
  _beads[0]._isRead     = 1;
  _beads[0]._isUnitig   = 0;
  _beads[0]._base       = base;
  _beads[0]._qual       = qual;

  _beads[0]._prevOffset = prevLink;
  _beads[0]._nextOffset = UINT16_MAX;

  _beadsLen++;

  assert(_beadsLen <= _beadsMax);

  if (prev)
    prev->_beads[prevLink]._nextOffset = 0;

#ifdef BASECOUNT
  baseCountIncr(base);
#endif

  //  With only the one read, we don't need to do any base calling here.

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "insertAtEnd()-- column prev=%p this=%p next=%p  link prev=%d next=%d\n",
          _prevColumn, this, _nextColumn, _beads[0]._prevOffset, _beads[0]._nextOffset);
#endif

  return(0);
}



//  Insert a column in the middle of the multialign, after some column.
uint16
abColumn::insertAfter(abColumn *prev,      //  Add new column after 'prev'
                      uint16    prevLink,  //  The bead for this read in 'prev' is at 'prevLink'.
                      char      base,
                      uint8     qual) {

  //  The base CAN NOT be a gap - the new column would then be entirely a gap column, with no base.

  assert(base != '-');
  assert(base != 0);

  //  Make sure that we're actually in the middle of consensus.

  assert(prev              != NULL);
  assert(prev->_nextColumn != NULL);

  //  Link in the new column (us!)

  _prevColumn = prev;
  _nextColumn = prev->_nextColumn;

  _prevColumn->_nextColumn = this;
  _nextColumn->_prevColumn = this;

  //  Allocate space for beads in this column (based on _prevColumn and _nextColumn)

  allocateInitialBeads();

  //  Add gaps for the existing reads.  This is quite complicated, so stashed away in a closet where we won't see it.

  inferPrevNextBeadPointers();

  //  Set up the new consensus base.

  _columnPosition = INT32_MAX;

  _call = base;
  _qual = qual;

  //  Then add the base for this read.  allocateInitialBeads() guarantees space for at least one
  //  more bead, infer() filled in the rest.

  uint16 tpos = _beadsLen++;

  assert(_beadsLen <= _beadsMax);

  _beads[tpos]._unused     = 0;
  _beads[tpos]._isRead     = 1;
  _beads[tpos]._isUnitig   = 0;
  _beads[tpos]._base       = base;
  _beads[tpos]._qual       = qual;

  _beads[tpos]._prevOffset = prevLink;
  _beads[tpos]._nextOffset = UINT16_MAX;

  //  Don't forget to update the link to us!  (I forgot.)

  if (prevLink != UINT16_MAX)
    _prevColumn->_beads[prevLink]._nextOffset = tpos;

#ifdef BASECOUNT
  baseCountIncr(base);
#endif
  baseCall(false);       //  We need to recall the base, using the majority vote.

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "insertAfter()-- column prev=%d this=%p next=%d  link %d prev=%d next=%d\n",
          _prevColumn->position(), this, _nextColumn->position(), tpos, _beads[tpos]._prevOffset, _beads[tpos]._nextOffset);
#endif

#ifdef DEBUG_INFER
  showLinks();
#endif

#ifdef CHECK_LINKS
  checkLinks();
#endif

  return(tpos);
}






uint16
abColumn::alignBead(uint16 prevIndex, char base, uint8 qual) {

  //  First, make sure the column has enough space for the new read.

  increaseArray(_beads, _beadsLen, _beadsMax, 1);

  //  Set up the new bead.

  uint16 tpos = _beadsLen++;

  assert(_beadsLen <= _beadsMax);

  _beads[tpos].initialize(base, qual, prevIndex, UINT16_MAX);

  //  Link to the previous

  if ((_prevColumn) && (prevIndex != UINT16_MAX)) {
    assert(prevIndex < _prevColumn->_beadsLen);
    assert(prevIndex < _prevColumn->_beadsMax);
    _prevColumn->_beads[prevIndex]._nextOffset = tpos;
  }

  //  increment the base count too!

#ifdef BASECOUNT
  baseCountIncr(base);
#endif
  baseCall(false);       //  We need to recall the base, using the majority vote.

#ifdef CHECK_LINKS
  checkLinks();
#endif

  //  Return the index of the bead we just added (for future linking)

  return(tpos);
};








//  Add sequence 'bid' to the multialign using ahang,bhang for positioning, and trace for the alignment.
//
//  ahang can be negative, zero or positive.
//
//  If negative, it will be equal in magnitude to bhang.
//
//  If positive, bhang will be zero, and this is the amount of sequence in frank we should ignore.
//  So, we set apos to be this positive value.

void
abAbacus::applyAlignment(uint32    bid,
                         int32     ahang,
                         int32     UNUSED(bhang),
                         int32    *trace, uint32 traceLen) {

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "abAbacus::applyAlignment()-- ahang=%d traceLen=%u trace=%d %d %d %d %d ... %d %d %d %d %d\n",
          ahang, traceLen,
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

  //  Finish some initialization.  If this is the first call, readTofBead (and readTolBead)
  //  will be NULL, and we need to allocate space for them.

  if (readTofBead == NULL) {
    readTofBead = new beadID [numberOfSequences()];
    readTolBead = new beadID [numberOfSequences()];
  }

  //  Figure out where we are in the multialignment.

  abSequence *bseq     = getSequence(bid);

  int32       alen     = _columnsLen;     //  Actual length of the consensus sequence BEFORE any new bases are added.
  int32       blen     = bseq->length();  //  Actual length of the read we're adding.

  //  We used to check that alen and blen were both positive, but we now allow applyAlignment() of
  //  the first read (to an empty multialignment) to initialize the structure.

  int32       apos     = MAX(ahang, 0);   //  if apos == alen, we'd just be pasting on new sequence.
  int32       bpos     = 0;               //  if bpos == blen...we're pasting on one base?

  assert(apos <= alen);
  assert(0    <= bpos);  //  We tried letting bpos be set to non-zero, to ignore bases at the start of the read,
  assert(bpos <= blen);  //  but it doesn't work.

  abColumn   *ncolumn  = _columns[apos];  //  The next empty column (where we will add the next aligned base).
  abColumn   *pcolumn  = NULL;            //  The previous column (that we just added a base to).
  uint16      plink    = UINT16_MAX;      //  The read index of the read we just added to the previous column.

  beadID      fBead;
  beadID      lBead;


  //  Negative ahang?  Push these things onto the start of the multialign.
  //
  //  We should fail if we get a negative ahang, but there is a column already before the first one
  //  (which implies we aligned to something not full-length, possibly because we trimmed
  //  frankenstein wrong).....but we don't even check.

  for (; bpos < -ahang; bpos++) {
    abColumn  *newcol = new abColumn;

    plink = newcol->insertAtBegin(ncolumn, plink, bseq->getBase(bpos), bseq->getQual(bpos));

    fBead.setF(newcol, plink);
    lBead.setL(newcol, plink);
  }

  if (ncolumn)
    pcolumn = ncolumn->_prevColumn;

  //  Skip any positive traces at the start.  These are template sequence (A-read) aligned to gaps
  //  in the B-read, but we don't care because there isn't any B-read aligned yet.
  //
  //  We probably should check that we didn't just process negative ahang above.

  for (; ((traceLen > 0) && (*trace != 0) && (*trace == 1)); trace++) {
    ncolumn = ncolumn->_nextColumn;
    apos++;
  }

  //  Similarly, remove negative traces at the end.  These are read sequence (B-read) aligned to gaps
  //  in the A-read (template).  The A-read extent should be reduced by one for each  gap.

  while ((traceLen > 0) && (trace[traceLen-1] > blen)) {
    trace[--traceLen] = 0;
  }


  //  Process the trace.

  while ((traceLen > 0) && (*trace != 0)) {

    //  Gap is in afrag.  Align ( - *trace - apos ) positions, then insert a new column, before
    //  ncolumn, to accommodate the insert in B.  ncolumn remains unchanged; it's still the column
    //  that we want to place the next aligned base.

    if ( *trace < 0 ) {
      while (apos < (- *trace - 1)) {
#ifdef DEBUG_ABACUS_ALIGN
        fprintf(stderr, "applyAlignment()--  align base %6d/%6d '%c' to column %7d\n", bpos, blen, bseq->getBase(bpos), ncolumn->position());
#endif

        plink = ncolumn->alignBead(plink, bseq->getBase(bpos), bseq->getQual(bpos));
        fBead.setF(ncolumn, plink);
        lBead.setL(ncolumn, plink);
        pcolumn = ncolumn;            //  ...updating the previous column
        ncolumn = ncolumn->next();    //
        bpos++;                       //  Move to the next base in the read.
        apos++;                       //  Move to the next column
      }

      assert(apos < alen);
      assert(bpos < blen);


      //  Three types of insert:
      //     before the multialign (done above in the initialization) - an entirely new column, no gaps to add
      //     after the multialign (done at the end) - same, entirely new column
      //     insertion to the multialign - in the middle, complicated
      //
      //  should it be insertAfter() to add a column before the exiting one and
      //  appendColumn() to add a column after?  appendColumm() is a special case of
      //  tacking on new sequence to the multialign.


      //  Add a new column for this insertion.
      abColumn  *newcol = new abColumn;

#ifdef DEBUG_ABACUS_ALIGN
      fprintf(stderr, "applyAlignment()--  align base %6d/%6d '%c' to after column %7d (new column)\n", bpos, blen, bseq->getBase(bpos), ncolumn->position());
#endif

      plink = newcol->insertAfter(pcolumn, plink, bseq->getBase(bpos), bseq->getQual(bpos));
      fBead.setF(newcol, plink);
      lBead.setL(newcol, plink);
      pcolumn = newcol;
      bpos++;  //  Move to the next base in the read.
    }

    //  Gap is in bfrag.  Align ( *trace - bpos ) positions, then insert a gap bead into the read,
    //  and align it to the existing column.

    if (*trace > 0) {
      while ( bpos < (*trace - 1) ) {
#ifdef DEBUG_ABACUS_ALIGN
        fprintf(stderr, "applyAlignment()--  align base %6d/%6d '%c' to column %7d\n", bpos, blen, bseq->getBase(bpos), ncolumn->position());
#endif

        plink = ncolumn->alignBead(plink, bseq->getBase(bpos), bseq->getQual(bpos));
        fBead.setF(ncolumn, plink);
        lBead.setL(ncolumn, plink);
        pcolumn = ncolumn;            //  ...updating the previous column
        ncolumn = ncolumn->next();    //
        bpos++;                       //  Move to the next base in the read.
        apos++;                       //  Move to the next column
      }

      assert(apos < alen);
      assert(bpos < blen);

#ifdef DEBUG_ABACUS_ALIGN
      fprintf(stderr, "applyAlignment()--  align base %6d/%6d '-' to column %7d (gap in read)\n", bpos, blen, ncolumn->position());
#endif

      plink = ncolumn->alignBead(plink, '-', 0);
      fBead.setF(ncolumn, plink);
      lBead.setL(ncolumn, plink);
      pcolumn = ncolumn;
      ncolumn = ncolumn->next();
      apos++;
    }

    trace++;
  }

  //  Remaining alignment contains no indels, just slap in the bases.  Note that when there is
  //  no consensus sequence (this is the first read added) this loop does nothing; alen=apos=0.

  for (int32 rem = MIN(blen - bpos, alen - apos); rem > 0; rem--) {
#ifdef DEBUG_ABACUS_ALIGN
    fprintf(stderr, "applyAlignment()--  align base %6d/%6d '%c' to column %7d (end of read)\n", bpos, blen, bseq->getBase(bpos), ncolumn->position());
#endif

    plink = ncolumn->alignBead(plink, bseq->getBase(bpos), bseq->getQual(bpos));
    fBead.setF(ncolumn, plink);
    lBead.setL(ncolumn, plink);
    pcolumn = ncolumn;
    ncolumn = ncolumn->next();
    apos++;
    bpos++;
  }

  //  Finally, append any new (unaligned) sequence from the read.  This handles the special case when there
  //  are no existing columns (pcolumn == NULL).

  for (int32 rem=blen-bpos; rem > 0; rem--) {
    assert(ncolumn == NULL);  //  Can't be a column after where we're tring to append to!

    abColumn *newcol = new abColumn;

#ifdef DEBUG_ABACUS_ALIGN
    fprintf(stderr, "applyAlignment()--  align base %6d/%6d '%c' to extend consensus\n", bpos, blen, bseq->getBase(bpos));
#endif

    plink = newcol->insertAtEnd(pcolumn, plink, bseq->getBase(bpos), bseq->getQual(bpos));
    fBead.setF(newcol, plink);
    lBead.setL(newcol, plink);
    pcolumn = newcol;
    bpos++;

#ifdef CHECK_LINKS
    newcol->checkLinks();
#endif
  }

  //  Insert the first and last beads into our tracking maps.

  assert(fBead.column->_beads[fBead.link].prevOffset() == UINT16_MAX);
  assert(lBead.column->_beads[lBead.link].nextOffset() == UINT16_MAX);

  fbeadToRead[fBead] = bid;
  readTofBead[bid] = fBead;

  lbeadToRead[lBead] = bid;
  readTolBead[bid] = lBead;

  //  Update the firstColumn in the abAbacus if it isn't set.  updateColumns() will
  //  reset it if the actual first column has changed here.

  if (_firstColumn == NULL)
    _firstColumn = fBead.column;

  //  Finally, recall bases (not needed; done inline when bases are aligned) and refresh the column/cnsBases/cnsQuals lists.

  //recallBases(false);
  refreshColumns();
}

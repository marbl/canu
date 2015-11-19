
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
 *    src/AS_CNS/MultiAlignment_CNS.C
 *    src/AS_CNS/MultiAlignment_CNS.c
 *    src/utgcns/libcns/MultiAlignment_CNS.C
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
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JUL-08
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-14
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "abAbacus.H"

#undef  DEBUG_ABACUS_ALIGN

//  Shouldn't be global, but some things -- like abBaseCount -- need it.

bool    DATAINITIALIZED                     = false;

double  EPROB[CNS_MAX_QV - CNS_MIN_QV + 1]  = { 0 };
double  PROB [CNS_MAX_QV - CNS_MIN_QV + 1]  = { 0 };

uint32  baseToIndex[256]                    = { 0 };
char    indexToBase[CNS_NUM_SYMBOLS]        = { 0 };



abAbacus::abAbacus(gkStore *gkpStore) {
  _sequencesLen    = 0;
  _sequencesMax    = 4;
  _sequences       = new abSequence [_sequencesMax];

  _basesLen        = 0;
  _basesMax        = 32;
  _bases           = new char [_basesMax];
  _quals           = new char [_basesMax];

  _beadsLen        = 0;
  _beadsMax        = _basesMax;
  _beads           = new abBead [_beadsMax];  //  These are slow to construct.

  _columnsLen      = 0;
  _columnsMax      = 32;
  _columns         = new abColumn[_columnsMax];  //  Or maybe it's these that are slow.

  _multiAlignsLen = 0;
  _multiAlignsMax = 2;
  _multiAligns    = new abMultiAlign [_multiAlignsMax];

  if (DATAINITIALIZED == false) {

    for (int32 i=0; i<256; i++)
      baseToIndex[i] = UINT32_MAX;

    indexToBase[ 0] = '-';
    indexToBase[ 1] = 'A';
    indexToBase[ 2] = 'C';
    indexToBase[ 3] = 'G';
    indexToBase[ 4] = 'T';
    indexToBase[ 5] = 'N';
#if 0
    indexToBase[ 6] = 'a';  //  -A
    indexToBase[ 7] = 'c';  //  -C
    indexToBase[ 8] = 'g';  //  -G
    indexToBase[ 9] = 't';  //  -T
    indexToBase[10] = 'M';  //  AC
    indexToBase[11] = 'R';  //  AG
    indexToBase[12] = 'W';  //  AT
    indexToBase[13] = 'S';  //  CG
    indexToBase[14] = 'Y';  //  CT
    indexToBase[15] = 'K';  //  GT
    indexToBase[16] = 'm';  //  -AC
    indexToBase[17] = 'r';  //  -AG
    indexToBase[18] = 'w';  //  -AT
    indexToBase[19] = 's';  //  -CG
    indexToBase[20] = 'y';  //  -CT
    indexToBase[21] = 'k';  //  -GT
    indexToBase[22] = 'V';  //  ACG
    indexToBase[23] = 'H';  //  ACT
    indexToBase[24] = 'D';  //  AGT
    indexToBase[25] = 'B';  //  CGT
    indexToBase[26] = 'v';  //  -ACG
    indexToBase[27] = 'h';  //  -ACT
    indexToBase[28] = 'd';  //  -AGT
    indexToBase[29] = 'b';  //  -CGT
    indexToBase[30] = 'X';  //  ACGT
    indexToBase[31] = 'x';  //  -ACGT
#endif

    for (int32 i=0; i<CNS_NUM_SYMBOLS; i++)
      baseToIndex[indexToBase[i]] = i;

    baseToIndex['n'] = baseToIndex['N'];  //  Used in baseCount


    double TAU_MISMATCH = 1.0 / (5.0 - 1.0);

    for (int32 i=0, qv=CNS_MIN_QV; i<CNS_MAX_QV - CNS_MIN_QV + 1; i++, qv++) {
      EPROB[i]= log(TAU_MISMATCH * pow(10, -qv/10.0));
      PROB[i] = log(1.0 - pow(10, -qv/10.0));

      //fprintf(stderr, "i %2d qv %2d PROB %f EPROB %f\n", i, qv, PROB[i], EPROB[i]);
    }

    DATAINITIALIZED = true;
  }
}


abAbacus::~abAbacus() {
  delete [] _multiAligns;
  delete [] _columns;
  delete [] _beads;
  delete [] _quals;
  delete [] _bases;
  delete [] _sequences;
}



abBeadID
abAbacus::addBead(char base, char qual) {

#if 0
  fprintf(stderr, "addBead()-- adding bead %u/%u base idx %u/%u %c/%d\n",
          _beadsLen, _beadsMax,
          _basesLen, _basesMax, base, base);
#endif

  increaseArray(_beads, _beadsLen, _beadsMax, 1);

  _beads[_beadsLen].boffset.set(_beadsLen);
  _beads[_beadsLen].soffset.set(_basesLen);

  _beads[_beadsLen].foffset       = UINT32_MAX;
  _beads[_beadsLen].prev          = abBeadID();
  _beads[_beadsLen].next          = abBeadID();
  _beads[_beadsLen].up            = abBeadID();
  _beads[_beadsLen].down          = abBeadID();
  _beads[_beadsLen].frag_index    = abSeqID();
  _beads[_beadsLen].column_index  = abColID();

  //  Add a base/qual for the bead

  increaseArrayPair(_bases, _quals, _basesLen, _basesMax, 1);

  _bases[_basesLen] = base;
  _quals[_basesLen] = qual;

  assert(CNS_MIN_QV + '0' <= _quals[_basesLen]);
  assert(_quals[_basesLen] <= CNS_MAX_QV + '0');

  _basesLen++;
  _beadsLen++;

  return(_beads[_beadsLen-1].ident());
};




abSeqBeadIterator *
abAbacus::createSeqBeadIterator(abSeqID sid) {
  return(new abSeqBeadIterator(this, getSequence(sid)));
}

abColBeadIterator *
abAbacus::createColBeadIterator(abColID cid) {
  return(new abColBeadIterator(this, getColumn(cid)));
}





//  Add a gap bead after bid.
abBeadID
abAbacus::appendGapBead(abBeadID bid) {

#warning this really should be using addBead()

  //  Make space for the new bead, and grab it.  And then regrab the prev.

#if 0
  fprintf(stderr, "appendGapBead()-- adding bead %u/%u base idx %u/%u\n",
          _beadsLen, _beadsMax,
          _basesLen, _basesMax);
#endif

  increaseArray(_beads, _beadsLen, _beadsMax, 1);

  increaseArrayPair(_bases, _quals, _basesLen, _basesMax, 1);

  abBead *bead = _beads + _beadsLen;
  abBead *prev = getBead(bid);

  assert(prev->frag_index.isValid());

  //  Set up the new bead

  bead->boffset.set(_beadsLen);
  bead->soffset.set(_basesLen);

  bead->foffset = prev->foffset + 1;

  bead->prev         = prev->ident();
  bead->next         = prev->next;
  bead->up           = abBeadID();
  bead->down         = abBeadID();
  bead->frag_index   = prev->frag_index;
  bead->column_index = abColID();

  prev->next         = bead->ident();

  if (bead->next.isValid())
    getBead(bead->next)->prev = bead->ident();

  //  Pick the minimum of the neighboring QVs, or '5' if both neighbors are zero, and set the base/qual

  char  qv = getQual(prev->soffset);

  if (bead->next.isValid()) {
    abBead *next = getBead(bead->next);
    char    nqv  = getQual(next->soffset);

    if (nqv < qv)
      qv = nqv;

    if (qv == '0')
      qv = '5';
  }

  _bases[_basesLen] = '-';
  _quals[_basesLen] = qv;

  assert(CNS_MIN_QV + '0' <= _quals[_basesLen]);
  assert(_quals[_basesLen] <= CNS_MAX_QV + '0');

  //  Finally, update the length of the beads/bases arrays.

  _beadsLen++;
  _basesLen++;

  //gaps_in_alignment++;

  return(bead->ident());
}





//  Add a gap bead before bid.
abBeadID
abAbacus::prependGapBead(abBeadID bid) {

  //  Make space for the new bead (and base), and grab the two beads.

  increaseArray(_beads, _beadsLen, _beadsMax, 1);

  increaseArrayPair(_bases, _quals, _basesLen, _basesMax, 1);

  abBead *bead = _beads + _beadsLen;
  abBead *next = getBead(bid);

  assert(next->frag_index.isValid());

  //  Set up the new bead

  bead->boffset.set(_beadsLen);
  bead->soffset.set(_basesLen);

  bead->foffset = next->foffset;  //  Same as prev+1

  bead->prev         = next->prev;
  bead->next         = next->ident();
  bead->up           = abBeadID();
  bead->down         = abBeadID();
  bead->frag_index   = next->frag_index;
  bead->column_index = abColID();

  next->prev         = bead->ident();

  if (bead->prev.isValid())
    getBead(bead->prev)->next = bead->ident();

  //  Pick the minimum of the neighboring QVs, or '5' if both neighbors are zero, and set the base/qual

  char  qv = getQual(next->soffset);

  if (bead->prev.isValid()) {
    abBead *prev = getBead(bead->prev);
    char    pqv  = getQual(prev->soffset);

    if (pqv < qv)
      qv = pqv;

    if (qv == '0')
      qv = '5';
  }

  _bases[_basesLen] = '-';
  _quals[_basesLen] = qv;

  assert(CNS_MIN_QV + '0' <= _quals[_basesLen]);
  assert(_quals[_basesLen] <= CNS_MAX_QV + '0');

  //  Finally, update the length of the beads/bases arrays.

  _beadsLen++;
  _basesLen++;

  //gaps_in_alignment++;

  return(bead->ident());
}








//  bid is the offset of the Bead seeding the column
//
abColID
abAbacus::appendColumn(abColID cid, abBeadID bid) {
  abColumn *prevCol   = getColumn(cid);
  abColumn *nextCol   = NULL;
  abColID   nextColID = addColumn(prevCol->ma_id, bid);  //  Could realloc

  prevCol   = getColumn(cid);
  nextCol   = getColumn(nextColID);

  abBead   *prevcall = getBead(prevCol->call);
  abBead   *nextcall = getBead(nextCol->call);

  //  Add the column to the list

  nextCol->next = prevCol->next;  //  Us into the list
  nextCol->prev = prevCol->lid;

  prevCol->next = nextCol->lid;   //  Previous points to us

  if (nextCol->next.isValid())    //  Next (if it exists) points to us
    getColumn(nextCol->next)->prev = nextCol->lid;

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "ColumnAppend()-- adding column "F_U32" (prev "F_U32" next "F_U32") for bid="F_U32" after column cid=%d\n",
          nextCol->ident().get(),
          nextCol->prev.get(),
          nextCol->next.get(),
          bid.get(),
          cid.get());
#endif

  //  Add the call to the list

  nextcall->next = prevcall->next;
  nextcall->prev = prevcall->ident();

  prevcall->next = nextcall->ident();

  if (nextcall->next.isValid())
    getBead(nextcall->next)->prev = nextcall->ident();

  //  Make new gap beads for every row in the previous column (as long as it wouldn't be a
  //  gap bead at the end of a sequence, and isn't the row with us in it)

  abColBeadIterator   *ci  = createColBeadIterator(cid);
  abBeadID             nid = ci->next();

  while (nid.isValid()) {
    abBead *bead = getBead(nid);

    if ((bead->next.isValid()) &&
        (bead->next != bid))
      alignBeadToColumn(nextCol->lid, appendGapBead(nid), "appendColumn()");

    nid = ci->next();
  }

  delete ci;

  nextCol->ma_position =  prevCol->ma_position + 1;

  getMultiAlign(nextCol->ma_id)->addColumnToMultiAlign(nextCol);  //  columnList

  return(nextCol->lid);
}





abColID
abAbacus::addColumn(abMultiAlignID mid, abBeadID bid) {
  abColID    colID;

  assert(mid.isValid());
  assert(bid.isValid());

  increaseArray(_columns, _columnsLen, _columnsMax, 1);

  colID.set(_columnsLen++);

  abColumn  *col   = getColumn(colID);

  col->lid         = colID;
  col->call        = addBead('N', '0');  //  New bead!
  col->next        = abColID();
  col->prev        = abColID();
  col->ma_id       = mid;
  col->ma_position = 0;

  col->base_count.clear();
  col->base_count.IncBaseCount('N');

  //  Set up the new base call bead; just need to link it in

  abBead    *call = getBead(col->call);  //  New bead.
  abBead    *bead = getBead(bid);

  call->down         = bead->ident();
  call->column_index = col->ident();

  bead->up           = call->ident();
  bead->column_index = col->ident();

  col->base_count.clear();
  col->base_count.IncBaseCount(getBase(bead->soffset));

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "addColumn()-- Added consensus call bead="F_U32" to new column="F_U32" for existing bead="F_U32"\n",
          call->ident().get(),
          col->ident().get(),
          col->prev.get(),
          col->next.get(),
          bead->ident().get());
#endif

  return(col->lid);
};







void
abAbacus::alignBeadToColumn(abColID cid, abBeadID bid, char *label) {
  abColumn  *column  = getColumn(cid);
  abBead    *call    = getBead(column->call);
  abBead    *first   = getBead(call->down);
  abBead    *align   = getBead(bid);

  //  Fails if called from appendColumn() because the baseIdx isn't valid

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "AlignBeadToColumn()-- %s frag=%d bead=%d,%c moving from column=%d (prev=N/A next=N/A) to column=%d (prev=%d next=%d)\n",
          label,
          align->seqIdx().get(),
          bid.get(),
          getBase(align->baseIdx()),
          align->colIdx().get(),
          column->ident().get(), column->prevID().get(), column->nextID().get());
#endif

  align->down         = first->ident();
  align->up           = call->ident();
  call->down          = align->ident();
  first->up           = align->ident();
  align->column_index = cid;

  column->base_count.IncBaseCount(getBase(align->soffset));
};





//  Remove bid from it's column, returning the next bead up in the column
abBeadID
abAbacus::unalignBeadFromColumn(abBeadID bid) {
  abBead *bead = getBead(bid);

  if (bead->column_index.isValid() == false)
    return(abBeadID());

  abColumn *column = getColumn(bead->column_index);
  abBead   *upbead = getBead(bead->up);
  char      bchar  = getBase(bead->soffset);

  upbead->down = bead->down;

  if (bead->down.isValid() )
    getBead(bead->down)->up = upbead->ident();

  column->base_count.DecBaseCount(bchar);

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "unAlignBeadFromColumn()-- frag=%d bead=%d leaving column=%d (prev=%d next=%d)\n",
          bead->frag_index.get(),
          bead->ident().get(),
          bead->column_index.get(),
          getColumn(bead->column_index)->prevID().get(),
          getColumn(bead->column_index)->nextID().get());
#endif

  bead->up           = abBeadID();
  bead->down         = abBeadID();
  bead->column_index = abColID();

  return(upbead->ident());
};






//  Remove bid from it's column, returning the prev or next bead in the fragment
//
void
abAbacus::unalignTrailingGapBeads(abBeadID bid) {

  // find direction to remove

  abBead    *bead   = getBead(bid);
  abBeadID   anchor = bead->prevID();


#ifdef DEBUG_UNALIGN_TRAILING
  fprintf(stderr, "unalignTrailingGapBeads()-- bid=%d %c anchor=%d %c\n",
          bead->ident().get(),
          getBase(bead->baseIdx()),
          anchor.get(),
          anchor.isValid() ? getBase(getBead(anchor)->baseIdx()) : '?');


  if (bead->nextID().isValid() == true)
    fprintf(stderr, "unalignTrailingGapBeads()-- next=%d %c\n",
            bead->nextID().get(),
            getBase(getBead(bead->nextID())->baseIdx()));
#endif

  while ((bead->nextID().isValid() == true) &&
         (getBase(getBead(bead->nextID())->baseIdx()) == '-')) {
#ifdef DEBUG_UNALIGN_TRAILING
    fprintf(stderr, "unalignTrailingGapBeads()-- bid=%d next is a gap\n", bead->ident().get());
#endif
    bead = getBead(bead->nextID());
  }


#ifdef DEBUG_UNALIGN_TRAILING
  if (bead->nextID().isValid() == true)
    fprintf(stderr, "unalignTrailingGapBeads()-- next=%d %c\n",
            bead->nextID().get(),
            getBase(getBead(bead->nextID())->baseIdx()));
#endif


  if (bead->nextID().isValid()) {
    anchor = bead->nextID();

#ifdef DEBUG_UNALIGN_TRAILING
    fprintf(stderr, "unalignTrailingGapBeads()-- reset anchor to %d\n", anchor.get());
#endif

    while ((bead->prev.isValid() == true) &&
           (getBase(getBead(bead->prevID())->baseIdx()) == '-')) {
#ifdef DEBUG_UNALIGN_TRAILING
      fprintf(stderr, "unalignTrailingGapBeads()-- bid=%d prev is a gap\n", bead->ident().get());
#endif
      bead = getBead(bead->prevID());
    }
  }


#ifdef DEBUG_UNALIGN_TRAILING
  if ((bead->prevID().isValid() == true) &&
      (bead->nextID().isValid() == true)) {
    fprintf(stderr, "unalignTrailingGapBeads()-- NOT TRAILING GAP BEADS!  Ignored.\n");
  }


  fprintf(stderr, "unalignTrailingGapBeads()-- bead ident=%d %c  anchor ident=%d %c\n",
          bead->ident().get(), getBase(bead->baseIdx()),
          anchor.get(), getBase(getBead(anchor)->baseIdx()));
#endif

  {
    abBead *anch = getBead(anchor);

#ifdef DEBUG_UNALIGN_TRAILING
    fprintf(stderr, "unalignTrailingGapBeads()-- bead %d %d %d  anchor %d %d %d\n",
            bead->prevID().get(), bead->ident().get(), bead->nextID().get(),
            anch->prevID().get(), anch->ident().get(), anch->nextID().get());

    //fprintf(stderr, "unalignTrailingGapBeads()-- bead %c %c %c  anchor %c %c %c\n",
    //        bead->prevID().get(), bead->ident().get(), bead->nextID().get(),
    //        anch->prevID().get(), anch->ident().get(), anch->nextID().get());
#endif
  }


  while (bead->ident() != anchor) {
    abColumn *column = getColumn(bead->colIdx());
    abBead   *upbead = getBead(bead->upID());

    char      bchar = getBase(bead->baseIdx());

#ifdef DEBUG_UNALIGN_TRAILING
    fprintf(stderr, "unalignTrailingGapBeads()-- bead ident=%d %c  anchor ident=%d %c (loop)\n",
            bead->ident().get(), getBase(bead->baseIdx()),
            anchor.get(), getBase(getBead(anchor)->baseIdx()));
#endif

    if (bchar != '-')
      fprintf(stderr, "unAlignTrailingGapBeads()-- bead %d %c is not a gap.\n",
              bead->ident().get(), bchar);
    assert(bchar == '-');

    upbead->downID() = bead->downID();

    if (bead->downID().isValid())
      getBead(bead->downID())->upID() = upbead->ident();

    column->base_count.DecBaseCount(bchar);

#ifdef DEBUG_UNALIGN_TRAILING
    fprintf(stderr, "UnAlignTrailingGapBeads()-- frag=%d bead=%d leaving column=%d\n",
            bead->seqIdx().get(), bead->ident().get(), bead->colIdx().get());
#endif

    bead->upID()   = abBeadID();
    bead->downID() = abBeadID();
    bead->colIdx() = abColID();

    if (bead->nextID().isInvalid()) {
      abBead *prevbead   = getBead(bead->prevID());

      prevbead->nextID() = abBeadID();
      bead->prevID()     = abBeadID();
      bead               = getBead(prevbead->ident());

    } else {
      abBead *nextbead   = getBead(bead->nextID());
      nextbead->prevID() = abBeadID();
      bead->nextID()     = abBeadID();
      bead               = getBead(nextbead->ident());
    }
  }

  //return(anchor);
}









//  This function swaps the contents of two beads, ensuring that
//  there are only gaps between them.
//
//  HORRIBLY complicated because applyAbacus() and mergeCompatible()
//  hold on to pointers to beads.  It would have been much simpler
//  to just swap the soffset and foffset, leaving EVERYTHING ELSE
//  exactly the same.
//
void
abAbacus::lateralExchangeBead(abBeadID lid, abBeadID rid) {
  abBead *leftbead  = getBead(lid);
  abBead *rightbead = getBead(rid);

  abColumn *leftcolumn  = getColumn(leftbead->colIdx());
  abColumn *rightcolumn = getColumn(rightbead->colIdx());

  char leftchar  = getBase(leftbead->baseIdx());
  char rightchar = getBase(rightbead->baseIdx());

  // now, verify that left and right are either
  // a) neighbors, or b) have only '-'s intervening

  {
    abBead *ibead   = leftbead;
    int32   failure = 0;
    int32   limit   = 20;

    while (ibead->next.isValid()) {
      ibead = getBead(ibead->nextID());

      if (ibead->ident() == rid)
        break;

      if (getBase(ibead->baseIdx()) != '-')
        failure++;
    }

    if (failure) {
      ibead = leftbead;

      while (ibead->next.isValid()) {
        ibead = getBead(ibead->nextID());

        if (ibead->ident() == rid)
          break;

        if (limit-- == 0)
          break;
      }

      fprintf(stderr, "lateralExchangeBead can't exchange bead "F_U32" with "F_U32"; bases in between!\n",
              lid.get(),
              rid.get());
      assert(failure == 0);
    }
  }

  abBead  rtmp = *rightbead;

  rightbead->upID()   = leftbead->upID();
  rightbead->downID() = leftbead->downID();
  rightbead->prevID() = leftbead->prevID();
  rightbead->nextID() = leftbead->nextID();

  if (rightbead->upID().isValid() )   (getBead(rightbead->upID()))->downID() = rid;
  if (rightbead->downID().isValid())  (getBead(rightbead->downID()))->upID() = rid;
  if (rightbead->prevID().isValid())  (getBead(rightbead->prevID()))->nextID() = rid;

  leftbead->upID()   = rtmp.upID();
  leftbead->downID() = rtmp.downID();
  leftbead->nextID() = rtmp.nextID();
  leftbead->prevID() = rtmp.prevID();

  if (leftbead->upID().isValid() )   (getBead(leftbead->upID()))->downID() = lid;
  if (leftbead->downID().isValid())  (getBead(leftbead->downID()))->upID() = lid;
  if (leftbead->nextID().isValid())  (getBead(leftbead->nextID()))->prevID() = lid;

  // now, handle separately cases of a) left and right are adjacent, and b) gaps intervene
  if ( rtmp.prevID() == lid) {
    rightbead->nextID() = lid;
    leftbead->prevID() = rid;
  } else {
    if (rightbead->nextID().isValid())  (getBead(rightbead->nextID()))->prevID() = rid;
    if (leftbead->prevID().isValid())   (getBead(leftbead->prevID()))->nextID() = lid;
  }

  rightbead->column_index = leftbead->column_index;
  leftbead->column_index  = rtmp.column_index;

  // change basecounts for affected columns
  assert(leftcolumn != NULL);

  //DecBaseCount(&leftcolumn->base_count,leftchar);
  //IncBaseCount(&leftcolumn->base_count,rightchar);
  leftcolumn->base_count.DecBaseCount(leftchar);
  leftcolumn->base_count.IncBaseCount(rightchar);

  assert(rightcolumn != NULL);

  //DecBaseCount(&rightcolumn->base_count,rightchar);
  //IncBaseCount(&rightcolumn->base_count,leftchar);
  rightcolumn->base_count.DecBaseCount(rightchar);
  rightcolumn->base_count.IncBaseCount(leftchar);
}








abMultiAlignID
abAbacus::addMultiAlign(abSeqID sid) {

  assert(sid.isValid());

  abSequence         *seq = getSequence(sid);
  abSeqBeadIterator   fi(this, seq);

  increaseArray(_multiAligns, _multiAlignsLen, _multiAlignsMax, 1);

  _multiAligns[_multiAlignsLen].identSet(_multiAlignsLen);

  abMultiAlign *ma     = _multiAligns + _multiAlignsLen++;

  abBeadID      bid    = fi.next();  //  The bead seeding this column
  abBead       *bead   = getBead(bid);

  abColID       cid    = addColumn(ma->ident(), bid);
  abColumn     *column = getColumn(cid);

  //  Add the column to the multiAlign; column pointers (in the column) are already set to 'nothing'

  assert(column->prevID().isValid() == false);
  assert(column->nextID().isValid() == false);

  ma->addColumnToMultiAlign(column);

  bid = fi.next();

  while (bid.isValid() == true) {
    cid = appendColumn(cid, bid);
    bid = fi.next();
  }

  //  Mark this sequence as belonging to this multiAlign
  //seq->addToMultiAlign(ma->ident());

  refreshMultiAlign(ma->ident());

  return(ma->ident());
};




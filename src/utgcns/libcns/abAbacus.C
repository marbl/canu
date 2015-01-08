


#include "abAbacus.H"





bool    DATAINITIALIZED                     = false;
double  EPROB[CNS_MAX_QV - CNS_MIN_QV + 1]  = { 0 };
double  PROB [CNS_MAX_QV - CNS_MIN_QV + 1]  = { 0 };
int32   RINDEX[256]                         = { 0 };
char    ALPHABET[6]                         = { 0 };
char    RALPHABET[CNS_NP]                   = { 0 };
char    RALPHABETC[CNS_NP]                  = { 0 };
double  TAU_MISMATCH                        =   0;
uint32  AMASK[5]                            = { 0 };





abBeadID
abAbacus::addBead(char base, char qual) {
  increaseArray(_beads, _beadsLen, _beadsMax, 1);

  _beads[_beadsLen].boffset.set(_beadsLen);
  _beads[_beadsLen].soffset       = abBaseID();
  _beads[_beadsLen].foffset       = UINT32_MAX;
  _beads[_beadsLen].prev          = abBeadID();
  _beads[_beadsLen].next          = abBeadID();
  _beads[_beadsLen].up            = abBeadID();
  _beads[_beadsLen].down          = abBeadID();
  _beads[_beadsLen].frag_index    = abSeqID();
  _beads[_beadsLen].column_index  = abColID();

  //  Add a base/qual for the bead
  if (base != 0) {
    increaseArray(_bases, _basesLen, _basesMax, 1);

    _bases[_basesLen]._base = base;
    _bases[_basesLen]._qual = qual;
    _basesLen++;
  }

  _beadsLen++;

  return(_beads[_beadsLen-1].ident());
};






//  Add a gap bead after bid.
abBeadID
abAbacus::appendGapBead(abBeadID bid) {
  abBead *prev = getBead(bid);

  assert(prev->frag_index.isValid());

  //  Make space for the new bead, and grab it.

  increaseArray(_beads, _beadsLen, _beadsMax, 1);

  abBead  *bead = _beads + _beadsLen++;

  //  Set up the new bead

  bead->boffset.set(_beadsLen);
  bead->soffset.set(_basesLen);

  bead->foffset = prev->foffset + 1;

  bead->prev         = prev->boffset;
  bead->next         = prev->next;
  bead->up           = abBeadID();
  bead->down         = abBeadID();
  bead->frag_index   = prev->frag_index;
  bead->column_index = abColID();

  prev->next         = bead->boffset;

  if (bead->next.isValid())
    getBead(bead->next)->prev = bead->boffset;

  //  Pick the minimum of the neighboring QVs, or '5' if both neighbors are zero

  char  qv = getQual(prev->soffset);

  if (bead->next.isValid()) {
    abBead *next = getBead(bead->next);
    char    nqv  = getQual(next->soffset);

    if (nqv < qv)
      qv = nqv;

    if (qv == '0')
      qv = '5';
  }

  //  Add the base/qual

  increaseArray(_bases, _basesLen, _basesMax, 1);
  _bases[_basesLen]._base = '-';
  _bases[_basesLen]._qual = qv;

  //gaps_in_alignment++;

  return(bead->boffset);
}





//  Add a gap bead before bid.
abBeadID
abAbacus::prependGapBead(abBeadID bid) {
  abBead *next = getBead(bid);

  assert(next->frag_index.isValid());

  //  Make space for the new bead, and grab it.

  increaseArray(_beads, _beadsLen, _beadsMax, 1);

  abBead  *bead = _beads + _beadsLen++;

  //  Set up the new bead

  bead->boffset.set(_beadsLen);
  bead->soffset.set(_basesLen);

  bead->foffset = next->foffset;  //  Same as prev+1

  bead->prev         = next->prev;
  bead->next         = next->boffset;
  bead->up           = abBeadID();
  bead->down         = abBeadID();
  bead->frag_index   = next->frag_index;
  bead->column_index = abColID();

  next->prev         = bead->boffset;

  if (bead->prev.isValid())
    getBead(bead->prev)->next = bead->boffset;

  //  Pick the minimum of the neighboring QVs, or '5' if both neighbors are zero

  char  qv = getQual(next->soffset);

  if (bead->prev.isValid()) {
    abBead *prev = getBead(bead->prev);
    char    pqv  = getQual(prev->soffset);

    if (pqv < qv)
      qv = pqv;

    if (qv == '0')
      qv = '5';
  }

  //  Add the base/qual

  increaseArray(_bases, _basesLen, _basesMax, 1);
  _bases[_basesLen]._base = '-';
  _bases[_basesLen]._qual = qv;

  //gaps_in_alignment++;

  return(bead->boffset);
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

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "ColumnAppend()-- adding column "F_U32" for bid="F_U32" after column cid=%d\n",
          nextCol->lid, bid.get(), cid.get());
#endif

  //  Add the column to the list

  nextCol->next = prevCol->next;  //  Us into the list
  nextCol->prev = prevCol->lid;

  prevCol->next = nextCol->lid;   //  Previous points to us

  if (nextCol->next.isValid())    //  Next (if it exists) points to us
    getColumn(nextCol->next)->prev = nextCol->lid;

  //  Add the call to the list

  nextcall->next = prevcall->next;
  nextcall->prev = prevcall->boffset;

  prevcall->next = nextcall->boffset;

  if (nextcall->next.isValid())
    getBead(nextcall->next)->prev = nextcall->boffset;

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
  nextCol->ma_id       =  prevCol->ma_id;
  nextCol->ma_position =  prevCol->ma_position + 1;

  //AddColumnToMANode(nextCol->ma_id, *column);

  abMultiAlign  *ma = getMultiAlign(nextCol->ma_id);

  ma->columnList.push_back(nextCol->ident());

  if (nextCol->prev.isValid() == false)
    ma->first = nextCol->ident();

  if (nextCol->next.isValid() == false)
    ma->last = nextCol->ident();

  return(nextCol->lid);
}





abColID
abAbacus::addColumn(abMultiAlignID mid, abBeadID bid) {
  abColID    colID;

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
  fprintf(stderr, "CreateColumn()-- Added consensus call bead="F_U32" to column="F_U32" for existing bead="F_U32"\n",
          call.boffset.get(), column.lid, bead->boffset.get());
#endif

  return(col->lid);
};







void
abAbacus::alignBeadToColumn(abColID cid, abBeadID bid, char *label) {
  abColumn  *column  = getColumn(cid);
  abBead    *call    = getBead(column->call);
  abBead    *first   = getBead(call->down);
  abBead    *align   = getBead(bid);

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "AlignBeadToColumn()-- %s frag=%d bead=%d,%c moving from column=%d to column=%d\n",
          label, align->frag_index, bid.get(), *Getchar(sequenceStore,align->soffset), align->column_index, cid);
#endif

  align->down         = first->boffset;
  align->up           = call->boffset;
  call->down          = align->boffset;
  first->up           = align->boffset;
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
    getBead(bead->down)->up = upbead->boffset;

  column->base_count.DecBaseCount(bchar);

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "UnAlignBeadFromColumn()-- frag=%d bead=%d leaving column=%d\n",
          bead->frag_index, bead->boffset.get(), bead->column_index);
#endif

  bead->up           = abBeadID();
  bead->down         = abBeadID();
  bead->column_index = abColID();

  return(upbead->ident());
};








void
abAbacus::lateralExchangeBead(abBeadID lid, abBeadID rid) {
  abBead rtmp; // this is just some tmp space for the swap

  //  This function swaps the contents of two beads, ensuring that
  //  there are only gaps between them.
  //
  //  HORRIBLY complicated because ApplyAbacus() and MergeCompatible()
  //  hold on to pointers to beads.  It would have been much simpler
  //  to just swap the soffset and foffset, leaving EVERYTHING ELSE
  //  exactly the same.

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

      if (ibead->boffset == rid)
        break;

      if (getBase(ibead->baseIdx()) != '-')
        failure++;
    }

    if (failure) {
      ibead = leftbead;

      while (ibead->next.isValid()) {
        ibead = getBead(ibead->nextID());

        fprintf(stderr, "bead %c boffset="F_U32" prev="F_U32" next="F_U32" up="F_U32" down="F_U32" fragindex=%d colulmnindex=%d\n",
                getBase(ibead->baseIdx()),
                ibead->boffset.get(),
                ibead->prevID().get(),
                ibead->nextID().get(),
                ibead->upID().get(),
                ibead->downID().get(),
                ibead->seqIdx().get(),
                ibead->colIdx().get());

        if (ibead->boffset == rid)
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

  rtmp = *rightbead;

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
  if ( rtmp.prev == lid) {
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

  abSequence         *seq = getSequence(sid);
  abSeqBeadIterator   fi(this, seq);

  increaseArray(_multiAligns, _multiAlignsLen, _multiAlignsMax, 1);

  abMultiAlign *ma = _multiAligns + _multiAlignsLen++;

  abBeadID      bid    = fi.next();  //  The bead seeding this column

  abBead       *bead   = getBead(bid);

  abColID       cid    = addColumn(ma->lid, bid);
  abColumn     *column = getColumn(cid);

  //  Add the column to the multiAlign; column pointers (in the column) are already set to 'nothing'

  ma->columnList.push_back(cid);

  ma->first = cid;
  ma->last  = cid;

  bid = fi.next();

  while (bid.isValid() == true) {
    cid = appendColumn(cid, bid);
    bid = fi.next();
  }

  //  Mark this sequence as belonging to this multiAlign

  seq->manode = ma->lid;

  refreshMultiAlign(ma->lid);

  return(ma->lid);
};




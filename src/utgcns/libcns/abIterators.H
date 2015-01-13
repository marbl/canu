


//  Iterate over all beads in a sequence
class abSeqBeadIterator {
public:
  abSeqBeadIterator(abAbacus *ab, abSequence *seq) {  //  used to take a fid, now needs _sequences + fid
    assert(seq != NULL);
    _abacus   = ab;
    _seq      = seq;
    _beadID   = _seq->firstBead();
    _isNull   = false;
  };
  ~abSeqBeadIterator() {
  };

  bool     isNull(void) {
    return(_isNull);
  };

  void     nullify(void) {
    _seq    = NULL;  //  used to be the 0th sequence
    _beadID = abBeadID();
    _isNull = true;
  };

  abBeadID  next(void) {
    assert(_isNull == false);

    abBeadID nid;

    if (_beadID.isValid()) {
      abBead *bead = _abacus->getBead(_beadID);
      nid          = bead->ident();
      _beadID      = bead->nextID();
    }

    return(nid);
  };

  //  Doesn't belong here, used only in multiAlign::print
  bool      advance(uint32 position) {
    abBeadID  bid  = next();
    abColumn *col  = NULL;

    do {
      if (bid.isValid() == false)
        return(false);

      abBead   *bead = _abacus->getBead(bid);

      if (bead->colIdx().isValid() == false)
        return(false);

      col  = _abacus->getColumn(bead->colIdx());
    } while (col->position() < position);

    assert(col->position() == position);

    return(true);
  };

private:
  abAbacus      *_abacus;
  abSequence    *_seq;
  abBeadID       _beadID;

  bool           _isNull;

  friend class abAbacus;
};




//  Iterate over all beads in a column
class abColBeadIterator {
public:
  abColBeadIterator(abAbacus *ab, abColumn *col) {
    assert(col != NULL);
    _abacus = ab;
    _col    = col;
    _beadID = _col->callID();
  };
  ~abColBeadIterator() {
  };

  abBeadID  next(void) {
    abBeadID nid;

    if (_beadID.isValid()) {
      abBead *bead = _abacus->getBead(_beadID);
      nid          = bead->downID();
      _beadID      = nid;
    }

    return(nid);
  }

private:
  abAbacus      *_abacus;
  abColumn      *_col;
  abBeadID       _beadID;

  friend class abAbacus;
};



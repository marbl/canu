
//  A single bead in the multialign
//
class abBead {
public:
  abBead() {
    clear();
  };
  ~abBead() {
  };

  //  ClearBead
  void        clear(void) {
    boffset      = abBeadID();
    soffset      = abBaseID();
    foffset      = abSeqID();
    prev         = abBeadID();
    next         = abBeadID();
    up           = abBeadID();
    down         = abBeadID();
    frag_index   = abSeqID();
    column_index = abColID();
  };

  abBeadID   ident(void)    { return(boffset); };

  abBeadID   prevID(void)   { return(prev);    };
  abBeadID   nextID(void)   { return(next);    };
  abBeadID   upID(void)     { return(up);      };
  abBeadID   downID(void)   { return(down);    };

  abBaseID   baseIdx(void)  { return(soffset); };
  abSeqID    seqIdx(void)   { return(foffset); };

private:
  abBeadID   boffset; // Location in BeadStore
  abBaseID   soffset; // Location in sequence/qualityStores
  abSeqID    foffset; // Location in Fragment sequence
  abBeadID   prev;
  abBeadID   next;
  abBeadID   up;
  abBeadID   down;         // navigation in multialignment (global offsets)
  abSeqID    frag_index;   // Location of containing fragment in fragmentStore
  abColID    column_index; // Location of alignment column in columnStore

  friend class abAbacus;
};


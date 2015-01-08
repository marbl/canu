
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
    foffset      = UINT32_MAX;
    prev         = abBeadID();
    next         = abBeadID();
    up           = abBeadID();
    down         = abBeadID();
    frag_index   = abSeqID();
    column_index = abColID();
  };

  abBeadID   ident(void)    { return(boffset);      };

  abBeadID   prevID(void)   { return(prev);         };
  abBeadID   nextID(void)   { return(next);         };
  abBeadID   upID(void)     { return(up);           };
  abBeadID   downID(void)   { return(down);         };

  abBaseID   baseIdx(void)  { return(soffset);      };
  abSeqID    seqIdx(void)   { return(frag_index);   };
  abColID    colIdx(void)   { return(column_index); };

private:
  abBeadID   boffset;       //  Location in abAbacus bead array
  abBaseID   soffset;       //  Location in abAbacus bases/quals array
  uint32     foffset;       //  Location in fragment sequence
  abBeadID   prev;          //  |
  abBeadID   next;          //  |
  abBeadID   up;            //  |- Navigation in abAbacus beads
  abBeadID   down;          //  |
  abSeqID    frag_index;    //  Location in abAbacus sequence array
  abColID    column_index;  //  Location in abAbacus column array

  friend class abAbacus;
};


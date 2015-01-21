
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

  abBeadID   ident(void)    const { return(boffset);      };

  abBeadID   prevID(void)   const { return(prev);         };
  abBeadID   nextID(void)   const { return(next);         };
  abBeadID   upID(void)     const { return(up);           };
  abBeadID   downID(void)   const { return(down);         };

  abBaseID   baseIdx(void)  const { return(soffset);      };
  abSeqID    seqIdx(void)   const { return(frag_index);   };
  abColID    colIdx(void)   const { return(column_index); };

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


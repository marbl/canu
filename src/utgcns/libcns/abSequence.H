
//  Former Fragment, a sequence loaded into the multialign
//

class abSequence {
public:
  abSequence() {
  };
  ~abSequence() {
  };

  uint32                gkpIdent(void)  { return(_iid);        };
  abSeqID               ident(void)     { return(_lid);        };

  abBeadID              firstBead(void)   { return(_firstbead);  };
  abBaseID              firstBase(void)   { return(_firstbase);  };
  //abColID               firstColumn(void) { needs abacus to get columns };

  //abBeadID              lastBead(void)    { return(_firstbead + _length - 1);  };
  //abBaseID              lastBase(void)    { return(_firstbase + _length - 1);  };
  abBeadID              lastBead(void)    { return(_lastbead);  };
  abBaseID              lastBase(void)    { return(_lastbase);  };
  //abColID               lastColumn(void)  { needs avacus to get columns };

  bool             isRead(void)    { return(_is_read);    };

  uint32           length(void)    { return(_length);     };

  bool             isForward(void) { return(_complement == false); };

  //abMultiAlignID   multiAlignID(void) { return(_multiAlignID); };
  //void             addToMultiAlign(abMultiAlignID maid) {_multiAlignID = maid;};

  void             initialize(uint32 readID,
                              uint32 readIdx,
                              uint32 length,
                              uint32 complemented,
                              uint32 firstBase,
                              uint32 firstBead) {
    _iid           = readID;
    _lid.set(readIdx);                 // index in sequence/quality/fragment store

    _is_read       = true;

    _length        = length;

    _complement    = complemented;
    //_container_iid = -1;
    //_is_contained  = false;
    //_deleted       = false;

    //_multiAlignID  = abMultiAlignID();

    _firstbase.set(firstBase);
    _firstbead.set(firstBead);

    _lastbase.set(firstBase + _length - 1);
    _lastbead.set(firstBead + _length - 1);

    //_n_components  = 0;
    //_components    = -1;
  };

private:
  uint32           _iid;            //  external, aka gatekeeper, ID;
  abSeqID          _lid;            //  index in sequence/quality/fragment store

  bool             _is_read;

  uint32           _length;

  bool             _complement;
  //uint32           _container_iid;  //  if non-zero, the iid of our container
  //bool             _is_contained;   //  if non-zero, consensus detected this fragment is contained
  //bool             _deleted;        //  deleted during consensus

  //  Never really used, was supposed to indicate if the read was present in this multialign.  The
  //  original code initialized this to -1, then reset it to NULL (ApplyAlignment.C:365).  On
  //  printing (PrintAlignment()), and maybe other places, reads with this value != mid (of the
  //  multialign we are printing) are skipped.  The multialign ID is _always_ 0.
  //
  //abMultiAlignID   _multiAlignID;

  abBaseID         _firstbase;      //  global index of first sequence character
  abBeadID         _firstbead;      //  global index of first "bead"

  abBaseID         _lastbase;       //  global index of first sequence character
  abBeadID         _lastbead;       //  global index of first "bead"

  //uint32           _n_components;   //  number of component frags (in case of "unitig" Fragments)
  //uint32           _components;     //  global index of first component frag
};


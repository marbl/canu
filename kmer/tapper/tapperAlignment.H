



class tapperAlignment {
public:

  //  Exactly as stored in tapperResultIndex and tapperResult*
  //  (Fragment, Singleton, Mated).
  //
  u64bit  _tagid;
  u32bit  _seq;
  u32bit  _pos;

  //  Exactly as stored in tapperResultQV
  //
  u32bit   _pad:15;

  u32bit   _basesMismatch:4;      //  Number of mismatches in ACGT alignment
  u32bit   _colorMismatch:4;      //  Number of consistent color mismatches
  u32bit   _colorInconsistent:4;  //  Number of inconsistent color mismatches
  u32bit   _rev:1;                //  Is reverse complement

  u32bit   _diffSize:4;           //  Value of MAX_COLOR_MISMATCH_MAPPED.

  u8bit    _colorDiffs[MAX_COLOR_MISMATCH_MAPPED];
};


class
tapperAlignmentPositionCompare {
public:
  bool operator()(const tapperAlignment &a, const tapperAlignment &b) const {
    return((a._seq < b._seq) ||
           ((a._seq == b._seq) && (a._pos < b._pos)));
  };
};

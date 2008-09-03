



class tapperAlignment {
public:

  //  Except as noted, all this is the same stuff as from tapperResult.

  u64bit   _tagid;
  u32bit   _seq;
  u32bit   _pos;

  u8bit    _colorDiffs[MAX_COLOR_MISMATCH_MAPPED];        //  OLD, list of errors in the align
  u8bit    _colorCorrections[MAX_COLOR_MISMATCH_MAPPED];  //  NEW, list of errors in the align, after corrections

  u32bit   _confidence:15;        //  NEW, confidence score on the corrected read

  u32bit   _basesMismatch:4;      //  Number of mismatches in ACGT alignment
  u32bit   _colorMismatch:4;      //  Number of consistent color mismatches
  u32bit   _colorInconsistent:4;  //  Number of inconsistent color mismatches
  u32bit   _rev:1;                //  Is reverse complement

  u32bit   _diffSize:4;           //  Value of MAX_COLOR_MISMATCH_MAPPED.
};


class
tapperAlignmentPositionCompare {
public:
  bool operator()(const tapperAlignment &a, const tapperAlignment &b) const {
    return((a._seq < b._seq) ||
           ((a._seq == b._seq) && (a._pos < b._pos)));
  };
};


class
tapperAlignmentTagIDAndScoreCompare {
public:
  bool operator()(const tapperAlignment &a, const tapperAlignment &b) const {
    return((a._tagid < b._tagid) ||
           ((a._tagid == b._tagid) && (a._confidence < b._confidence)));

  };
};

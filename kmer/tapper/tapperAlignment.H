



class tapperAlignment {
public:

  //  Except as noted, all this is the same stuff as from tapperResult.

  uint64   _tagid;
  uint32   _seq;
  uint32   _pos;

  uint8    _colorDiffs[MAX_COLOR_MISMATCH_MAPPED];        //  OLD, list of errors in the align
  uint8    _colorCorrections[MAX_COLOR_MISMATCH_MAPPED];  //  NEW, list of errors in the align, after corrections

  uint32   _confidence:15;        //  NEW, confidence score on the corrected read

  uint32   _basesMismatch:4;      //  Number of mismatches in ACGT alignment
  uint32   _colorMismatch:4;      //  Number of consistent color mismatches
  uint32   _colorInconsistent:4;  //  Number of inconsistent color mismatches
  uint32   _rev:1;                //  Is reverse complement

  uint32   _diffSize:4;           //  Value of MAX_COLOR_MISMATCH_MAPPED.
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

#include "alphabet.h"

class tapperGlobalData;


//  An internal hit.  Tapper uses these for computing scores and what
//  not.  It outputs tapperResults, above.
//
class tapperHit {
public:
  u32bit    numberOfBaseMismatches(void)        { return(_basesMismatch);     };
  u32bit    numberOfColorMismatches(void)       { return(_colorMismatch);     };
  u32bit    numberOfColorInconsistencies(void)  { return(_colorInconsistent); };


  char *printHit(char *OS, u64bit tagid) {
    sprintf(OS, "0x"u64bitFMT"\t"u32bitFMT":"u32bitFMT":%c\t"u64bitFMT","u64bitFMT","u64bitFMT,
            tagid,
            _seqIdx, _seqPos, _rev ? '-' : '+',
            _basesMismatch, _colorMismatch, _colorInconsistent);
    return(OS);
  }


  //  Returns true if the tag is near the correct end of the sequence,
  //  so that it could potentially be happily mated to a tag mapping
  //  in a different sequence (or in a gap).
  //
  void   alignToReference(tapperGlobalData *g,
                          u32bit so,
                          u32bit po,
                          char  *tag, u32bit len);


  bool   happyNearEnd(bool isFTag, u32bit mean, u32bit stddev, u32bit seqlen) {
    bool  isHappy = false;

    if (seqlen < mean + 3 * stddev)
      return(true);

    if (isFTag) {
      if (_rev) {
        //  Near end of sequence
        isHappy = (seqlen - mean - 3 * stddev < _seqPos);
      } else {
        //  Near bgn of sequence
        isHappy = (_seqPos < mean + 3 * stddev);
      }
    } else {
      if (_rev) {
        //  Near bgn of sequence
        isHappy = (_seqPos < mean + 3 * stddev);
      } else {
        //  Near end of sequence
        isHappy = (seqlen - mean - 3 * stddev < _seqPos);
      }
    }

    return(isHappy);
  };


  //  Returns true if that read is before where this read says it
  //  should be.  Returns TRUE for reads of the incorrect orientation.
  //
  //  ASSUMES it is called on this read being the forward/F3/a read.
  //
  bool   mateTooFarBefore(tapperHit& that, u32bit mean, u32bit stddev) {

    //  that read is on the sequence after us.
    if (_seqIdx < that._seqIdx) {
      //fprintf(stderr, "isBefore()- seq after false.\n");
      return(false);
    }

    //  that read is on the sequence before us.
    if (that._seqIdx < _seqIdx) {
      //fprintf(stderr, "isBefore()- seq before true.\n");
      return(true);
    }

    //  Misoriented, true
    if (_rev != that._rev) {
      //fprintf(stderr, "isBefore()- misoriented true.\n");
      return(true);
    }

    //  FORWARD
    //
    //     (  -that->  )      -this->
    //  -----------------------------------
    //  TTTffffffffffffffffffffffffffffffff
    //  
    if ((_rev == false) && (that._seqPos + mean + 3 * stddev < _seqPos)) {
      //fprintf(stderr, "isBefore()- forward true "u32bitFMT" + "u32bitFMT" + 3 * "u32bitFMT" < "u32bitFMT"\n",
      //        that._seqPos, mean, stddev, _seqPos);
      return(true);
    }

    //  REVERSE
    //
    //       <-this-       (  <-that-  )
    //  -----------------------------------
    //  TTTTTTTTTTTTTTTTTTTffffffffffffffff
    //
    if ((_rev == true)  && (that._seqPos < _seqPos + mean - 3 * stddev)) {
      //fprintf(stderr, "isBefore()- forward true "u32bitFMT" < "u32bitFMT" + "u32bitFMT" - 3 * "u32bitFMT"\n",
      //        that._seqPos, _seqPos, mean, stddev);
      return(true);
    }

    //fprintf(stderr, "isBefore()- false.\n");
    return(false);
  }


  //  Returns true if that read is after where this read says it
  //  should be.  Returns FALSE for reads of the incorrect orientation.
  //
  //  ASSUMES it is called on this read being the forward/F3/a read.
  //
  bool   mateTooFarAfter(tapperHit& that, u32bit mean, u32bit stddev) {

    //  that read is on the sequence after us, true.
    if (_seqIdx < that._seqIdx)
      return(true);

    //  that read is on the sequence before us, false.
    if (that._seqIdx < _seqIdx)
      return(false);

    //  Misoriented, true
    if (_rev != that._rev)
      return(false);

    //  FORWARD
    //
    //     (  -that->  )      -this->
    //  -----------------------------------
    //  ffffffffffffffffTTTTTTTTTTTTTTTTTTT
    //  
    if ((_rev == false) && (that._seqPos + mean - 3 * stddev < _seqPos))
      return(false);

    //  REVERSE
    //
    //       <-this-       (  <-that-  )
    //  -----------------------------------
    //  ffffffffffffffffffffffffffffffffTTT
    //
    if ((_rev == true)  && (that._seqPos < _seqPos + mean + 3 * stddev))
      return(false);

    return(true);
  }



  //  ASSUMES it is called on this read being the forward/F3/a read.
  //
  bool   happy(tapperHit& b, u32bit mean, u32bit stddev) {
    u64bit  dist        = ~u64bitZERO;
    bool    isHappy     = false;
    bool    isOriented  = false;

    if (_seqIdx != b._seqIdx)
      return(false);

    if (_rev != b._rev)
      return(false);

    //  Check distance apart
    if (b._seqPos < _seqPos)
      dist = _seqPos - b._seqPos;
    else
      dist = b._seqPos - _seqPos;
    if ((mean - 3 * stddev < dist) && (dist < mean + 3 * stddev))
      isHappy = true;

    //  Check orientations
    if ((_rev == false) && (b._seqPos < _seqPos))
      isOriented = true;
    if ((_rev == true)  && (_seqPos < b._seqPos))
      isOriented = true;

    if (!isHappy) {
      //fprintf(stderr, "GRUMPY DIST  "u32bitFMT"\n", dist);
      return(false);
    }

    if (!isOriented) {
      //fprintf(stderr, "GRUMPY ORIENT "u32bitFMT"\n", dist);
      return(false);
    }

    //fprintf(stderr, "HAPPY!  "u32bitFMT"\n", dist);
    return(true);
  };
  
#if 0
  bool    operator< (tapperHit const &r) const { return(((_basesMismatch     <  r._basesMismatch)) ||
                                                        ((_basesMismatch     <= r._basesMismatch) && (_colorMismatch     <  r._colorMismatch)) ||
                                                        ((_basesMismatch     <= r._basesMismatch) && (_colorMismatch     <= r._colorMismatch) && (_colorInconsistent <  r._colorInconsistent))); };
#endif

  //  Argh, should be private, but tapperWorker copies most of the hit to a result.
  //private:
  u32bit  _seqIdx;
  u32bit  _seqPos;

  u64bit  _tagIdx;               //  4e9 tags is only 34x of human

  u64bit  _len:6;                //  Length of tag
  u64bit  _rev:1;                //  Match is reversecomplement

  u64bit  _pad:17;               //  Nothing

  u64bit  _basesMismatch:6;      //  Number of mismatches in ACGT alignment

  u64bit  _colorMismatch:6;      //  Number of consistent color mismatches
  u64bit  _colorInconsistent:6;  //  Number of inconsistent color mismatches

  char    _tagCOLOR[TAG_LEN_MAX];
  char    _refCOLOR[TAG_LEN_MAX];

  char    _tagACGT[TAG_LEN_MAX];
  char    _refACGT[TAG_LEN_MAX];

  u8bit   _changes[MAX_COLOR_MISMATCH_MAPPED];
};


class
tapperHitPositionCompare {
public:
  bool operator()(const tapperHit a, const tapperHit b) const {
    return((a._seqIdx < b._seqIdx) ||
           ((a._seqIdx == b._seqIdx) && (a._seqPos < b._seqPos)));
  };
};

#include "alphabet.h"

#include <functional>

//  Tapper generates four kinds of alignments.
//
//  1)  An unmated fragment alignment
//  2)  A satisfied mate pair alignment
//  3)  An unsatisfied mate pair alignment
//  4)  A tangle of mated fragments
//
//  There are SIX output files, an index, an alignment quality
//  histogram, and the four data files.


#define MAX_FRAGMENT_ALIGNMENTS              65536  //  16 bits
#define MAX_FRAGMENT_ALIGNMENTS_DISCARDED  1048576  //  20 bits
#define MAX_SINGLETON_ALIGNMENTS             65536  //  16 bits
#define MAX_MATED_ALIGNMENTS                  8192  //  13 bits
#define MAX_TANGLED_ALIGNMENTS                8192  //  13 bits

#define MAX_INSERT_SIZE                     262144  //  18 bits
#define MAX_INSERT_DEVIATION                 65536  //  16 bits

#define MAX_COLOR_MISMATCH_MAPPED          4

//  Info about alignments for one mate pair.  One per pair or unmated
//  fragment.  This is the index.
//  256 bits.
//
class tapperResult {
public:
  void    print(FILE *out) {
    u16bit                 id1[4];
    u16bit                 id2[4];

    decodeTagID(_tag1id, id1);
    decodeTagID(_tag2id, id2);

    fprintf(out, "R\t"u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT"\t"u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT"\t"u64bitFMT"/"u64bitFMT"\t"u64bitFMT"+-"u64bitFMT"\tf:"u64bitFMT"\td:"u64bitFMT"\ts:"u64bitFMT"\tm:"u64bitFMT"\tt:"u64bitFMT"\n",
            id1[0], id1[1], id1[2], id1[3],
            id2[0], id2[1], id2[2], id2[3],
            _maxColrMismatchMapped, _maxBaseMismatchMapped,
            _mean, _stddev,
            _numFragment, _numFragmentDiscarded, _numSingleton, _numMated, _numTangled);
  };

public:
  u64bit   _tag1id;
  u64bit   _tag2id;

  //  Command line, how many color mismatches we looked for, and how
  //  many base mismatches we allowed.  These deterine the number and
  //  meaning of the alignment quality histogram.  This is stored per
  //  result, so multiple runs can be easily combined.
  //
  u64bit   _maxColrMismatchMapped:4;   
  u64bit   _maxBaseMismatchMapped:4;

  //  Expected mean and stddev for this pair.  Again, per result so we
  //  can combine mappings.
  //
  u64bit   _mean:18;                   
  u64bit   _stddev:16;

  u64bit   _numFragmentDiscarded:20;   //  Number of fragment alignments found but not reported

  u64bit   _pad1:2;
  //  62 bits above, 58 bits below
  u64bit   _pad2:6;

  u64bit   _numFragment:16;            //  Number of fragment alignments
  u64bit   _numSingleton:16;           //  Number of unsatisfied alignments
  u64bit   _numMated:13;               //  Number of mated alignments
  u64bit   _numTangled:13;             //  Number of tangled alignments
};


//  Quality for a mated alignment.
//  32 bits for quality
//  64 bits for alignment (= 2 * MAX_COLOR_MISMATCH_MAPPED * 8 bits)
//
//  The alignments take up a lot of space.  We store both the position
//  of the difference, and the base in the reference.
//
class tapperResultQV {
public:
  u32bit   _tag1valid:1;              //  Tag 1 is valid data
  u32bit   _tag1basesMismatch:4;      //  Number of mismatches in ACGT alignment
  u32bit   _tag1colorMismatch:4;      //  Number of consistent color mismatches
  u32bit   _tag1colorInconsistent:4;  //  Number of inconsistent color mismatches
  u32bit   _tag1rev:1;                //  Is reverse complement

  u32bit   _tag2valid:1;              //  Tag 2 is valid data
  u32bit   _tag2basesMismatch:4;      //  Number of mismatches in ACGT alignment
  u32bit   _tag2colorMismatch:4;      //  Number of consistent color mismatches
  u32bit   _tag2colorInconsistent:4;  //  Number of inconsistent color mismatches
  u32bit   _tag2rev:1;                //  Is reverse complement

  u32bit   _diffSize:4;               //  Value of MAX_COLOR_MISMATCH_MAPPED.

  u8bit    _tag1colorDiffs[MAX_COLOR_MISMATCH_MAPPED];
  u8bit    _tag2colorDiffs[MAX_COLOR_MISMATCH_MAPPED];
};


//  Variable sized array.  Size determined by the two "max mismatch
//  mapped" values in a tapperResult.
//
class tapperResultHistogram {
public:
  u32bit            _histogram[0];
};


//  Unmated fragment alignment.
//  96 bits.
//
class tapperResultFragment {
public:
  void    print(FILE *out, tapperResult *res) {
    u16bit                 id[4];

    if (_qual._tag1valid) {
      decodeTagID(res->_tag1id, id);
      fprintf(stdout, "F\t"u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%c\t"u32bitFMT"/"u32bitFMT"/"u32bitFMT"\n",
              id[0], id[1], id[2], id[3],
              _seq,
              _pos,
              _qual._tag1rev ? 'r' : 'f',
              _qual._tag1basesMismatch,
              _qual._tag1colorMismatch,
              _qual._tag1colorInconsistent);
    }

    if (_qual._tag2valid) {
      decodeTagID(res->_tag2id, id);
      fprintf(stdout, "F\t"u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%c\t"u32bitFMT"/"u32bitFMT"/"u32bitFMT"\n",
              id[0], id[1], id[2], id[3],
              _seq,
              _pos,
              _qual._tag2rev ? 'r' : 'f',
              _qual._tag2basesMismatch,
              _qual._tag2colorMismatch,
              _qual._tag2colorInconsistent);
    }
  };


public:
  u32bit            _seq;
  u32bit            _pos;
  union {
    u32bit          _bits;
    tapperResultQV  _qual;
  };
};


//  Satisfied mate pair alignment.
//  128 bits.
//
class tapperResultMated {
public:
  void    print(FILE *out, tapperResult *res) {
    u16bit                 id1[4];
    u16bit                 id2[4];

    decodeTagID(res->_tag1id, id1);
    decodeTagID(res->_tag2id, id2);

    fprintf(stdout, "M\t"u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%c\t"u32bitFMT"/"u32bitFMT"/"u32bitFMT"\t"u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%c\t"u32bitFMT"/"u32bitFMT"/"u32bitFMT"\n",
            id1[0], id1[1], id1[2], id1[3],
            _seq,
            _pos1,
            _qual._tag1rev ? 'r' : 'f',
            _qual._tag1basesMismatch,
            _qual._tag1colorMismatch,
            _qual._tag1colorInconsistent,
            id2[0], id2[1], id2[2], id2[3],
            _seq,
            _pos2,
            _qual._tag2rev ? 'r' : 'f',
            _qual._tag2basesMismatch,
            _qual._tag2colorMismatch,
            _qual._tag2colorInconsistent);
  };

public:
  u32bit             _seq;
  u32bit             _pos1;
  u32bit             _pos2;
  union {
    u32bit           _bits;
    tapperResultQV   _qual;
  };
};


//  Unsatisfied mate pair alignment; only one tag.
//  96 bits.
//
class tapperResultSingleton {
public:
  void    print(FILE *out, tapperResult *res) {
    u16bit                 id[4];

    if (_qual._tag1valid) {
      decodeTagID(res->_tag1id, id);
      fprintf(stdout, "S\t"u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%c\t"u32bitFMT"/"u32bitFMT"/"u32bitFMT"\n",
              id[0], id[1], id[2], id[3],
              _seq,
              _pos,
              _qual._tag1rev ? 'r' : 'f',
              _qual._tag1basesMismatch,
              _qual._tag1colorMismatch,
              _qual._tag1colorInconsistent);
    }

    if (_qual._tag2valid) {
      decodeTagID(res->_tag2id, id);
      fprintf(stdout, "S\t"u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%c\t"u32bitFMT"/"u32bitFMT"/"u32bitFMT"\n",
              id[0], id[1], id[2], id[3],
              _seq,
              _pos,
              _qual._tag2rev ? 'r' : 'f',
              _qual._tag2basesMismatch,
              _qual._tag2colorMismatch,
              _qual._tag2colorInconsistent);
    }
  };
public:
  u32bit             _seq;
  u32bit             _pos;
  union {
    u32bit           _bits;
    tapperResultQV   _qual;
  };
};


//  Tangled mate pair alignment.
//  128 bits.
//
class tapperResultTangled {
public:
  void    print(FILE *out, tapperResult *res) {
    u16bit                 id1[4];
    u16bit                 id2[4];

    decodeTagID(res->_tag1id, id1);
    decodeTagID(res->_tag2id, id2);

    fprintf(stdout, "T\t"u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT"\t"u32bitFMT"\t"u16bitFMT"_"u16bitFMT"_"u16bitFMT"_"u16bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n",
            id1[0], id1[1], id1[2], id1[3],
            _tag1count,
            id2[0], id2[1], id2[2], id2[3],
            _tag2count,
            _seq,
            _bgn,
            _end);
  };

public:
  u16bit           _tag1count;  //  Number of times tag1 is in here
  u16bit           _tag2count;  //  Number of times tag2 is in here
  u32bit           _seq;        //  Sequence we hit
  u32bit           _bgn;        //  Beginning location
  u32bit           _end;        //  Ending location
};



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
  
  bool    operator< (tapperHit const &r) const { return(((_basesMismatch     <  r._basesMismatch)) ||
                                                        ((_basesMismatch     <= r._basesMismatch) && (_colorMismatch     <  r._colorMismatch)) ||
                                                        ((_basesMismatch     <= r._basesMismatch) && (_colorMismatch     <= r._colorMismatch) && (_colorInconsistent <  r._colorInconsistent))); };

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
};


class
tapperHitPositionCompare {
public:
  bool operator()(const tapperHit a, const tapperHit b) const {
    return((a._seqIdx < b._seqIdx) ||
           ((a._seqIdx == b._seqIdx) && (a._seqPos < b._seqPos)));
  };
};

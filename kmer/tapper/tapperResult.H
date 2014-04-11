#include "util++.H"

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
#define MAX_FRAGMENT_ALIGNMENTS_TANGLED    1048576  //  20 bits
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
class tapperResultIndex {
public:
  void    print(FILE *out) {
    uint16                 id1[4];
    uint16                 id2[4];

    decodeTagID(_tag1id, id1);
    decodeTagID(_tag2id, id2);

    fprintf(out, "R\t"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t"uint64FMT"/"uint64FMT"\t"uint64FMT"+-"uint64FMT"\tf:"uint64FMT"\td:"uint64FMT"\ts:"uint64FMT"\tm:"uint64FMT"\tt:"uint64FMT"\n",
            id1[0], id1[1], id1[2], id1[3],
            id2[0], id2[1], id2[2], id2[3],
            _maxColrMismatchMapped, _maxBaseMismatchMapped,
            _mean, _stddev,
            _numFrag, _numFragDiscarded, _numFragSingleton, _numMated, _numTangled);
  };

public:
  uint64   _tag1id;
  uint64   _tag2id;

  //  Command line, how many color mismatches we looked for, and how
  //  many base mismatches we allowed.  These deterine the number and
  //  meaning of the alignment quality histogram.  This is stored per
  //  result, so multiple runs can be easily combined.
  //
  uint64   _maxColrMismatchMapped:4;   
  uint64   _maxBaseMismatchMapped:4;

  uint64   _mean:18;               //  Expected mean and stddev for this pair.
  uint64   _stddev:16;             //  Again, per result so we can combine mappings.
  uint64   _pad1:22;

  uint64   _numFrag:16;            //  Number of fragment alignments
  uint64   _numFragDiscarded:20;   //  Number of fragment alignments found but not reported
  uint64   _numFragTangled:20;     //  Number of fragment alignments in tangled mated
  uint64   _numFragSingleton:16;   //  Number of fragment alignments potentially linking
  uint64   _pad2:8;

  uint64   _numMated:13;           //  Number of mated alignments
  uint64   _numTangled:13;         //  Number of tangled alignments
  uint64   _pad3:22;
};


//  Quality for a mated alignment.
//  32 bits for quality
//  64 bits for alignment (= 2 * MAX_COLOR_MISMATCH_MAPPED * 8 bits)
//
//  The alignments take up a lot of space.  We store both the position
//  of the difference, and the color in the read.
//
class tapperResultQV {
public:
  uint32   _tag1valid:1;              //  Tag 1 is valid data
  uint32   _tag1basesMismatch:4;      //  Number of mismatches in ACGT alignment
  uint32   _tag1colorMismatch:4;      //  Number of consistent color mismatches
  uint32   _tag1colorInconsistent:4;  //  Number of inconsistent color mismatches
  uint32   _tag1rev:1;                //  Is reverse complement

  uint32   _tag2valid:1;              //  Tag 2 is valid data
  uint32   _tag2basesMismatch:4;      //  Number of mismatches in ACGT alignment
  uint32   _tag2colorMismatch:4;      //  Number of consistent color mismatches
  uint32   _tag2colorInconsistent:4;  //  Number of inconsistent color mismatches
  uint32   _tag2rev:1;                //  Is reverse complement

  uint32   _diffSize:4;               //  Value of MAX_COLOR_MISMATCH_MAPPED.

  uint8    _tag1colorDiffs[MAX_COLOR_MISMATCH_MAPPED];
  uint8    _tag2colorDiffs[MAX_COLOR_MISMATCH_MAPPED];
};



//  Unmated fragment alignment.
//  96 bits.
//
class tapperResultFragment {
public:
  void    print(FILE *out, tapperResultIndex *idx) {
    uint16                 id[4];
    char                   cor[128];
    uint32                 err = 0;

#warning do not know real tag length
    memset(cor, '.', 128);
    cor[26] = 0;

    if (_qual._tag1valid) {
      for (uint32 x=0; x<_qual._tag1colorMismatch; x++, err++) {
        uint32 pos = _qual._tag1colorDiffs[err] & 0x3f;
        cor[pos] = '*';
      }

      for (uint32 x=0; x<_qual._tag1colorInconsistent; x++, err++) {
        uint32 pos = _qual._tag1colorDiffs[err] & 0x3f;
        cor[pos] = bitsToColor[_qual._tag1colorDiffs[err] >> 6];
      }

      decodeTagID(idx->_tag1id, id);
      fprintf(stdout, "F\t"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t"uint32FMT"\t"uint32FMT"\t%c\t"uint32FMT"/"uint32FMT"/"uint32FMT"\t'%s'\n",
              id[0], id[1], id[2], id[3],
              _seq,
              _pos,
              _qual._tag1rev ? 'r' : 'f',
              _qual._tag1basesMismatch,
              _qual._tag1colorMismatch,
              _qual._tag1colorInconsistent,
              cor);
    }

    if (_qual._tag2valid) {
      for (uint32 x=0; x<_qual._tag2colorMismatch; x++, err++) {
        uint32 pos = _qual._tag2colorDiffs[err] & 0x3f;
        cor[pos] = '*';
      }

      for (uint32 x=0; x<_qual._tag2colorInconsistent; x++, err++) {
        uint32 pos = _qual._tag2colorDiffs[err] & 0x3f;
        cor[pos] = bitsToColor[_qual._tag2colorDiffs[err] >> 6];
      }

      decodeTagID(idx->_tag2id, id);
      fprintf(stdout, "F\t"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t"uint32FMT"\t"uint32FMT"\t%c\t"uint32FMT"/"uint32FMT"/"uint32FMT"\t'%s'\n",
              id[0], id[1], id[2], id[3],
              _seq,
              _pos,
              _qual._tag2rev ? 'r' : 'f',
              _qual._tag2basesMismatch,
              _qual._tag2colorMismatch,
              _qual._tag2colorInconsistent,
              cor);
    }
  };


public:
  uint32            _seq;
  uint32            _pos;
  tapperResultQV    _qual;
};


//  Satisfied mate pair alignment.
//  128 bits.
//
class tapperResultMated {
public:
  void    print(FILE *out, tapperResultIndex *idx) {
    uint16                 id1[4];
    uint16                 id2[4];

    decodeTagID(idx->_tag1id, id1);
    decodeTagID(idx->_tag2id, id2);

    fprintf(stdout, "M\t"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t"uint32FMT"\t"uint32FMT"\t%c\t"uint32FMT"/"uint32FMT"/"uint32FMT"\t"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t"uint32FMT"\t"uint32FMT"\t%c\t"uint32FMT"/"uint32FMT"/"uint32FMT"\n",
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
  uint32             _seq;
  uint32             _pos1;
  uint32             _pos2;
  tapperResultQV     _qual;
};


//  Tangled mate pair alignment.
//  128 bits.
//
class tapperResultTangled {
public:
  void    print(FILE *out, tapperResultIndex *idx) {
    uint16                 id1[4];
    uint16                 id2[4];

    decodeTagID(idx->_tag1id, id1);
    decodeTagID(idx->_tag2id, id2);

    fprintf(stdout, "T\t"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t"uint32FMT"\t"uint16FMT"_"uint16FMT"_"uint16FMT"_"uint16FMT"\t"uint32FMT"\t"uint32FMT"\t"uint32FMT"\t"uint32FMT"\n",
            id1[0], id1[1], id1[2], id1[3],
            _tag1count,
            id2[0], id2[1], id2[2], id2[3],
            _tag2count,
            _seq,
            _bgn,
            _end);
  };

public:
  uint16           _tag1count;  //  Number of times tag1 is in here
  uint16           _tag2count;  //  Number of times tag2 is in here
  uint32           _seq;        //  Sequence we hit
  uint32           _bgn;        //  Beginning location
  uint32           _end;        //  Ending location
};




class tapperResult {
public:
  tapperResult() {
    memset(&idx, 0, sizeof(tapperResultIndex));
    fragMax = 0;
    frag    = 0L;
    singMax = 0;
    sing    = 0L;
    taliMax = 0;
    tali    = 0L;
    mateMax = 0;
    mate    = 0L;
    tangMax = 0;
    tang    = 0L;
    aqltMax = 0;
    aqlt    = 0L;
  };
  ~tapperResult() {
    delete [] frag;
    delete [] sing;
    delete [] tali;
    delete [] mate;
    delete [] tang;
    delete [] aqlt;
  };

  tapperResultIndex       idx;

  //  A single unmated alignment
  uint32                  fragMax;
  tapperResultFragment   *frag;

  //  Tag in a mate pair, mapped near the end of a sequence
  uint32                  singMax;
  tapperResultFragment   *sing;

  //  Tag in a mate pair, involved in a tangle
  uint32                  taliMax;
  tapperResultFragment   *tali;

  //  Happy mated tags
  uint32                  mateMax;
  tapperResultMated      *mate;

  //  Location of tangle
  uint32                  tangMax;
  tapperResultTangled    *tang;

  uint32                  aqltMax;
  uint32                 *aqlt;
};





class tapperAlignmentQualityHistogramIndices {
public:
  tapperAlignmentQualityHistogramIndices() {
    for (uint32 i=0; i<16; i++)
      for (uint32 j=0; j<16; j++) {
        _indices[i][j] = 0L;
        _length[i][j]  = ~uint32ZERO;
      }
  };
  ~tapperAlignmentQualityHistogramIndices() {
    for (uint32 i=0; i<16; i++)
      for (uint32 j=0; j<16; j++)
        delete [] _indices[i][j];
  };

  //  For a given maxColorError and maxBaseError (mapper parameters),
  //  maps between (numBaseMismatch, numColorMismatch, numColorError)
  //  and an index in an array.
  //
  //  A maximum of 16 is allowed on all values.

  uint32   getLength(uint32 maxBaseError, uint32 maxColorError) {
    generate(maxBaseError, maxColorError);
    assert(_length[maxBaseError][maxColorError] < ~uint32ZERO);
    return(_length[maxBaseError][maxColorError]);
  };

  uint32   getIndex(uint32 maxBaseError, uint32 maxColorError,
                    uint32 numBaseMismatch, uint32 numColorMismatch, uint32 numColorError) {
    generate(maxBaseError, maxColorError);
    assert(_length[maxBaseError][maxColorError] < ~uint32ZERO);
    assert(numBaseMismatch * 256 + numColorMismatch * 16 + numColorError < 16 * 16 * 16);
    return(_indices[maxBaseError][maxColorError][numBaseMismatch * 256 + numColorMismatch * 16 + numColorError]);
  };

private:
  void     generate(uint32 maxBaseError, uint32 maxColorError) {

    if (_indices[maxBaseError][maxColorError] != 0L)
      return;

    //  min base mismatches for i color mismatches - the min is (I
    //  think always) the sum of the mins for the prime decomposition.
    //    9 - 3,3,3   -> min 6 base mismatches
    //    9 - 2,3,4   -> min 5 base mismatches
    //    9 - 2,2,5   -> min 5 base mismatches
    //    9 - 2,2,2,3 -> min 5 base mismatches
    //
    //  max base mismatches is, for the most part, used defined, but 0
    //  and 1 color mismatches are forced to 0 color mismatches.
    //
    //  finally, it is impossible to have just one color mismatch.

    //                cm 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
    uint32 bmmin[16] = { 0, 0, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8 };

    uint32 maxc = maxColorError;  //  max color errors given to the mapper
    uint32 maxb = maxBaseError;   //  max base  errors given to the mapper

    uint32 index = 0;

    uint32 *histogramIndices = _indices[maxBaseError][maxColorError] = new uint32 [16 * 16 * 16];

    for (uint32 ii=0; ii < 16 * 16 * 16; ii++)
      histogramIndices[ii] = ~uint32ZERO;

    //  A special case for 0.
    for (uint32 ce=0; ce <= maxc; ce++) {
      //fprintf(stderr, "histogramIndices["uint32FMT"/"uint32FMT"/"uint32FMT"] = "uint32FMT"\n", 0, 0, ce, index);
      assert(0 * 256 + 0 * 16 + ce < 16 * 16 * 16);
      histogramIndices[0 * 256 + 0 * 16 + ce] = index++;
    }
    
    for (uint32 cm=2; cm <= maxc; cm++)
      for (uint32 ce=0; ce <= maxc - cm; ce++)
        for (uint32 bm=bmmin[cm]; bm <= maxb; bm++) {
          //fprintf(stderr, "histogramIndices["uint32FMT"/"uint32FMT"/"uint32FMT"] = "uint32FMT"\n", bm, cm, ce, index);
          assert(bm * 256 + cm * 16 + ce < 16 * 16 * 16);
          histogramIndices[bm * 256 + cm * 16 + ce] = index++;
        }

    _length[maxBaseError][maxColorError] = index;
  };


  uint32  *_indices[16][16];
  uint32   _length[16][16];
};





class tapperResultFile {
public:
  tapperResultFile(char *prefix, char mode) {
    char filename[FILENAME_MAX];

    if ((mode != 'r') && (mode != 'w'))
      fprintf(stderr, "tapperResultFile()-- mode must be either 'r' or 'w'; you wanted '%c'\n", mode), exit(1);

    if ((mode == 'r') && (!fileExists(prefix)))
      fprintf(stderr, "tapperResultFile()-- result directory '%s' doesn't exist.\n", prefix), exit(1);

    if ((mode == 'w') && (!fileExists(prefix))) {
      errno = 0;
      mkdir(prefix, S_IRWXU | S_IRWXG | S_IRWXO);
      if (errno)
        fprintf(stderr, "tapperResultFile()-- failed to make result directory '%s': %s\n",
                prefix, strerror(errno)), exit(1);
    }

    sprintf(filename, "%s/tapperMappedIndex", prefix);
    IDX  = new recordFile(filename, 0, sizeof(tapperResultIndex), mode);

    sprintf(filename, "%s/tapperMappedFragment", prefix);
    FRAG = new recordFile(filename, 0, sizeof(tapperResultFragment), mode);

    sprintf(filename, "%s/tapperMappedSingleton", prefix);
    SING = new recordFile(filename, 0, sizeof(tapperResultFragment), mode);

    sprintf(filename, "%s/tapperMappedTangledAlignment", prefix);
    TALI = new recordFile(filename, 0, sizeof(tapperResultFragment), mode);

    sprintf(filename, "%s/tapperMappedMated", prefix);
    MATE = new recordFile(filename, 0, sizeof(tapperResultMated), mode);

    sprintf(filename, "%s/tapperMappedTangled", prefix);
    TANG = new recordFile(filename, 0, sizeof(tapperResultTangled), mode);

    sprintf(filename, "%s/tapperMappedAlignQual", prefix);
    AQLT = new recordFile(filename, 0, sizeof(uint32), mode);
  };

  ~tapperResultFile() {
    delete IDX;
    delete FRAG;
    delete SING;
    delete TALI;
    delete MATE;
    delete TANG;
    delete AQLT;
  };

  static
  bool    validResultFile(char *prefix) {
    return(fileExists(prefix));
  };

  uint32  AQIlength(uint32 maxBaseErrors, uint32 maxColorErrors) {
    //fprintf(stderr, "AQIlength("uint32FMT","uint32FMT") -> "uint32FMT"\n",
    //        maxBaseErrors, maxColorErrors, AQI.getLength(maxBaseErrors, maxColorErrors));
    return(AQI.getLength(maxBaseErrors, maxColorErrors));
  };

  uint32  AQIindex(uint32 maxBaseErrors, uint32 maxColorErrors,
                   uint32 numBaseMismatch, uint32 numColorMismatch, uint32 numColorError) {
    //fprintf(stderr, "AQIindex("uint32FMT","uint32FMT","uint32FMT","uint32FMT","uint32FMT") -> "uint32FMT"\n",
    //        maxBaseErrors, maxColorErrors,
    //        numBaseMismatch, numColorMismatch, numColorError,
    //        AQI.getIndex(maxBaseErrors, maxColorErrors, numBaseMismatch, numColorMismatch, numColorError));
    return(AQI.getIndex(maxBaseErrors, maxColorErrors, numBaseMismatch, numColorMismatch, numColorError));
  };

  bool    read(tapperResult         *align) {
    bool  success = true;

    if (IDX->getRecord(&align->idx) == 0)
      return(false);

    uint32 aqilen = AQIlength(align->idx._maxBaseMismatchMapped, align->idx._maxColrMismatchMapped);

    if (align->idx._numFrag +
        align->idx._numFragDiscarded +
        align->idx._numFragSingleton +
        align->idx._numMated +
        align->idx._numTangled == 0)
      aqilen = 0;

    if (align->idx._numFrag > align->fragMax) {
      delete [] align->frag;
      align->fragMax = align->idx._numFrag;
      align->frag    = new tapperResultFragment [align->fragMax];
    }

    if (align->idx._numFragSingleton > align->singMax) {
      delete [] align->sing;
      align->singMax = align->idx._numFragSingleton;
      align->sing    = new tapperResultFragment [align->singMax];
    }

    if (align->idx._numFragTangled > align->taliMax) {
      delete [] align->tali;
      align->taliMax = align->idx._numFragTangled;
      align->tali    = new tapperResultFragment [align->taliMax];
    }

    if (align->idx._numMated > align->mateMax) {
      delete [] align->mate;
      align->mateMax = align->idx._numMated;
      align->mate    = new tapperResultMated [align->mateMax];
    }

    if (align->idx._numTangled > align->tangMax) {
      delete [] align->tang;
      align->tangMax = align->idx._numTangled;
      align->tang    = new tapperResultTangled [align->tangMax];
    }

    if (aqilen > align->aqltMax) {
      delete [] align->aqlt;
      align->aqltMax = aqilen;
      align->aqlt    = new uint32 [align->aqltMax];
    }

#if 0
    fprintf(stderr, "reading: "uint32FMT" "uint32FMT" "uint32FMT" "uint32FMT" "uint32FMT"\n",
            align->idx._numFrag,
            align->idx._numFragSingleton,
            align->idx._numFragTangled,
            align->idx._numMated,
            align->idx._numTangled);
#endif

    if (FRAG->getRecord(align->frag, align->idx._numFrag)          != align->idx._numFrag)
      success = false;
    if (SING->getRecord(align->sing, align->idx._numFragSingleton) != align->idx._numFragSingleton)
      success = false;
    if (TALI->getRecord(align->tali, align->idx._numFragTangled)   != align->idx._numFragTangled)
      success = false;
    if (MATE->getRecord(align->mate, align->idx._numMated)         != align->idx._numMated)
      success = false;
    if (TANG->getRecord(align->tang, align->idx._numTangled)       != align->idx._numTangled)
      success = false;

    if (AQLT->getRecord(align->aqlt, aqilen) != aqilen)
      success = false;

    return(success);
  };

  void    write(tapperResult        *align) {
    write(&align->idx, align->frag, align->sing, align->tali, align->mate, align->tang, align->aqlt);
  };

  void    write(tapperResultIndex      *idx,
                tapperResultFragment   *frag,
                tapperResultFragment   *sing,
                tapperResultFragment   *tali,
                tapperResultMated      *mate,
                tapperResultTangled    *tang,
                uint32                 *aqlt) {
    IDX->putRecord(idx);

    FRAG->putRecord(frag, idx->_numFrag);
    SING->putRecord(sing, idx->_numFragSingleton);
    TALI->putRecord(tali, idx->_numFragTangled);
    MATE->putRecord(mate, idx->_numMated);
    TANG->putRecord(tang, idx->_numTangled);

    if (idx->_numFrag +
        idx->_numFragDiscarded +
        idx->_numFragSingleton +
        idx->_numMated +
        idx->_numTangled > 0)
      AQLT->putRecord(aqlt, AQIlength(idx->_maxBaseMismatchMapped, idx->_maxColrMismatchMapped));
  };

private:
  tapperAlignmentQualityHistogramIndices  AQI;

  recordFile  *IDX;
  recordFile  *FRAG;
  recordFile  *SING;
  recordFile  *TALI;
  recordFile  *MATE;
  recordFile  *TANG;
  recordFile  *AQLT;
};


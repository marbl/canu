


class tapperAlignmentQualityHistogramIndices {
public:
  tapperAlignmentQualityHistogramIndices() {
    for (u32bit i=0; i<16; i++)
      for (u32bit j=0; j<16; j++) {
        _indices[i][j] = 0L;
        _length[i][j]  = ~u32bitZERO;
      }
  };
  ~tapperAlignmentQualityHistogramIndices() {
    for (u32bit i=0; i<16; i++)
      for (u32bit j=0; j<16; j++)
        delete [] _indices[i][j];
  };

  //  For a given maxColorError and maxBaseError (mapper parameters),
  //  maps between (numBaseMismatch, numColorMismatch, numColorError)
  //  and an index in an array.
  //
  //  A maximum of 16 is allowed on all values.

  u32bit   getLength(u32bit maxBaseError, u32bit maxColorError) {
    generate(maxBaseError, maxColorError);
    assert(_length[maxBaseError][maxColorError] < ~u32bitZERO);
    return(_length[maxBaseError][maxColorError]);
  };

  u32bit   getIndex(u32bit maxBaseError, u32bit maxColorError,
                    u32bit numBaseMismatch, u32bit numColorMismatch, u32bit numColorError) {
    generate(maxBaseError, maxColorError);
    assert(_length[maxBaseError][maxColorError] < ~u32bitZERO);
    assert(numBaseMismatch * 256 + numColorMismatch * 16 + numColorError < 16 * 16 * 16);
    return(_indices[maxBaseError][maxColorError][numBaseMismatch * 256 + numColorMismatch * 16 + numColorError]);
  };

private:
  void     generate(u32bit maxBaseError, u32bit maxColorError) {

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
    u32bit bmmin[16] = { 0, 0, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8 };

    u32bit maxc = maxColorError;  //  max color errors given to the mapper
    u32bit maxb = maxBaseError;   //  max base  errors given to the mapper

    u32bit index = 0;

    u32bit *histogramIndices = _indices[maxBaseError][maxColorError] = new u32bit [16 * 16 * 16];

    for (u32bit ii=0; ii < 16 * 16 * 16; ii++)
      histogramIndices[ii] = ~u32bitZERO;

    //  A special case for 0.
    for (u32bit ce=0; ce <= maxc; ce++) {
      //fprintf(stderr, "histogramIndices["u32bitFMT"/"u32bitFMT"/"u32bitFMT"] = "u32bitFMT"\n", 0, 0, ce, index);
      assert(0 * 256 + 0 * 16 + ce < 16 * 16 * 16);
      histogramIndices[0 * 256 + 0 * 16 + ce] = index++;
    }
    
    for (u32bit cm=2; cm <= maxc; cm++)
      for (u32bit ce=0; ce <= maxc - cm; ce++)
        for (u32bit bm=bmmin[cm]; bm <= maxb; bm++) {
          //fprintf(stderr, "histogramIndices["u32bitFMT"/"u32bitFMT"/"u32bitFMT"] = "u32bitFMT"\n", bm, cm, ce, index);
          assert(bm * 256 + cm * 16 + ce < 16 * 16 * 16);
          histogramIndices[bm * 256 + cm * 16 + ce] = index++;
        }

    _length[maxBaseError][maxColorError] = index;
  };


  u32bit  *_indices[16][16];
  u32bit   _length[16][16];
};


class tapperAlignment {
public:
  tapperAlignment() {
    memset(&idx, 0, sizeof(tapperResultIndex));
    fragMax = 0;
    frag    = 0L;
    singMax = 0;
    sing    = 0L;
    mateMax = 0;
    mate    = 0L;
    tangMax = 0;
    tang    = 0L;
    aqltMax = 0;
    aqlt    = 0L;
  };
  ~tapperAlignment() {
    delete [] frag;
    delete [] sing;
    delete [] mate;
    delete [] tang;
    delete [] aqlt;
  };

  tapperResultIndex       idx;

  u32bit                  fragMax;
  tapperResultFragment   *frag;

  u32bit                  singMax;
  tapperResultSingleton  *sing;

  u32bit                  mateMax;
  tapperResultMated      *mate;

  u32bit                  tangMax;
  tapperResultTangled    *tang;

  u32bit                  aqltMax;
  u32bit                 *aqlt;
};



class tapperAlignmentFile {
public:
  tapperAlignmentFile(char *prefix, char mode) {
    char filename[FILENAME_MAX];

    if ((mode != 'r') && (mode != 'w'))
      fprintf(stderr, "tapperAlignmentFile()-- mode must be either 'r' or 'w'; you wanted '%c'\n", mode), exit(1);

      sprintf(filename, "%s.tapperMappedIndex", prefix);
      IDX  = new recordFile(filename, 0, sizeof(tapperResultIndex), mode);

      sprintf(filename, "%s.tapperMappedFragment", prefix);
      FRAG = new recordFile(filename, 0, sizeof(tapperResultFragment), mode);

      sprintf(filename, "%s.tapperMappedMated", prefix);
      MATE = new recordFile(filename, 0, sizeof(tapperResultMated), mode);

      sprintf(filename, "%s.tapperMappedSingleton", prefix);
      SING = new recordFile(filename, 0, sizeof(tapperResultSingleton), mode);

      sprintf(filename, "%s.tapperMappedTangled", prefix);
      TANG = new recordFile(filename, 0, sizeof(tapperResultTangled), mode);

      sprintf(filename, "%s.tapperMappedAlignQual", prefix);
      AQLT = new recordFile(filename, 0, sizeof(u32bit), mode);
  };

  ~tapperAlignmentFile() {
    delete IDX;
    delete FRAG;
    delete SING;
    delete MATE;
    delete TANG;
    delete AQLT;
  };

  u32bit  AQIlength(u32bit maxBaseErrors, u32bit maxColorErrors) {
    //fprintf(stderr, "AQIlength("u32bitFMT","u32bitFMT") -> "u32bitFMT"\n",
    //        maxBaseErrors, maxColorErrors, AQI.getLength(maxBaseErrors, maxColorErrors));
    return(AQI.getLength(maxBaseErrors, maxColorErrors));
  };

  u32bit  AQIindex(u32bit maxBaseErrors, u32bit maxColorErrors,
                   u32bit numBaseMismatch, u32bit numColorMismatch, u32bit numColorError) {
    //fprintf(stderr, "AQIindex("u32bitFMT","u32bitFMT","u32bitFMT","u32bitFMT","u32bitFMT") -> "u32bitFMT"\n",
    //        maxBaseErrors, maxColorErrors,
    //        numBaseMismatch, numColorMismatch, numColorError,
    //        AQI.getIndex(maxBaseErrors, maxColorErrors, numBaseMismatch, numColorMismatch, numColorError));
    return(AQI.getIndex(maxBaseErrors, maxColorErrors, numBaseMismatch, numColorMismatch, numColorError));
  };

  bool    read(tapperAlignment         *align) {
    bool  success = true;

    if (IDX->getRecord(&align->idx) == 0)
      return(false);

    u32bit aqilen = AQIlength(align->idx._maxBaseMismatchMapped, align->idx._maxColrMismatchMapped);

    if (align->idx._numFragment +
        align->idx._numFragmentDiscarded +
        align->idx._numSingleton +
        align->idx._numMated +
        align->idx._numTangled == 0)
      aqilen = 0;

    if (align->idx._numFragment > align->fragMax) {
      delete [] align->frag;
      align->fragMax = align->idx._numFragment;
      align->frag    = new tapperResultFragment [align->fragMax];
    }

    if (align->idx._numSingleton > align->singMax) {
      delete [] align->sing;
      align->singMax = align->idx._numSingleton;
      align->sing    = new tapperResultSingleton [align->singMax];
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
      align->aqlt    = new u32bit [align->aqltMax];
    }

    if (FRAG->getRecord(align->frag, align->idx._numFragment)  != align->idx._numFragment)
      success = false;
    if (SING->getRecord(align->sing, align->idx._numSingleton) != align->idx._numSingleton)
      success = false;
    if (MATE->getRecord(align->mate, align->idx._numMated)     != align->idx._numMated)
      success = false;
    if (TANG->getRecord(align->tang, align->idx._numTangled)   != align->idx._numTangled)
      success = false;

    if (AQLT->getRecord(align->aqlt, aqilen) != aqilen)
      success = false;

    return(success);
  };

  void    write(tapperAlignment        *align) {
    write(&align->idx, align->frag, align->sing, align->mate, align->tang, align->aqlt);
  };

  void    write(tapperResultIndex      *idx,
                tapperResultFragment   *frag,
                tapperResultSingleton  *sing,
                tapperResultMated      *mate,
                tapperResultTangled    *tang,
                u32bit                 *aqlt) {
    IDX->putRecord(idx);

    FRAG->putRecord(frag, idx->_numFragment);
    SING->putRecord(sing, idx->_numSingleton);
    MATE->putRecord(mate, idx->_numMated);
    TANG->putRecord(tang, idx->_numTangled);

    if (idx->_numFragment +
        idx->_numFragmentDiscarded +
        idx->_numSingleton +
        idx->_numMated +
        idx->_numTangled > 0)
      AQLT->putRecord(aqlt, AQIlength(idx->_maxBaseMismatchMapped, idx->_maxColrMismatchMapped));
  };

private:
  tapperAlignmentQualityHistogramIndices  AQI;

  recordFile  *IDX;
  recordFile  *FRAG;
  recordFile  *SING;
  recordFile  *MATE;
  recordFile  *TANG;
  recordFile  *AQLT;
};

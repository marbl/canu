


//  A single tag, binary encoded.
//
//  The current ascii encoding for a 25bp tag needs 44 bytes for
//  sequence and 110 bytes for QVs -> 154 bytes per tag.
//
//  Without QVs, we can fit upto a 60bp tag into 24 bytes, using a
//  64bit global id instead of the sequence name.  At 24B per tag, 10x
//  human is 27GB.
//
//  Including QVs, we now need to use 7 bits per bp.
//
//  The infrastructure of the bitPackedFile is used, so all we need to
//  define is the number of words in our tapperTag (which, since we
//  already do a similar hack for a kMer, isn't so terrible).
//
//  WORDS     1  2  3  4  5  6  7  8  9 10 11 12  13  14  15  16
//  BYTES     8 16 24 32 40 48 56 64 72 80 88 96 104 112 120 128
//  TAG SIZE  0  7 17 26 35 44 53 62 71 81 90 99 108 117 126 135
//

#define TAPPER_TAG_WORDS 4
//#define TAG_LEN_MAX   32

class tapperTag {
public:
  tapperTag() {};

  //  Expects seq to be "T01023030122303" and qlt to be 0 through 31,
  //  1-1 correspondence.
  //
  void     encode(u64bit tagID, char *seq, u64bit *qlt) {
    u64bit   pos = 64;
    u32bit   len = strlen(seq);
    u32bit   i   = 0;

    for (i=1; i<TAPPER_TAG_WORDS; i++)
      _w[i] = u64bitZERO;

    _w[0] = tagID;

    setDecodedValue(_w, pos,  2, compressSymbol[seq[0]]);
    pos += 2;

    setDecodedValue(_w, pos,  7, len-1);
    pos += 7;

    //fprintf(stderr, "seq: %s id "u64bitHEX" %c/%d len "u32bitFMT" %d "u64bitFMT"\n",
    //        seq, _w[0], seq[0], seq[0], len, pos-7, getDecodedValue(_w, pos-7, 7));

    for (i=1; i<len; i++) {
      setDecodedValue(_w, pos, 7, (compressSymbol[seq[i]] << 5) | (qlt[i]));
      pos += 7;
    }

    assert(_w[0] == tagID);
  };

  //  Decodes a tag into seq and qlt (as for encode()), returns the
  //  tagID, or 0 if failure.
  //
  u64bit   decode(char *seq, u64bit *qlt) {
    u64bit   pos = 64;
    u32bit   len = 0;
    u32bit   i   = 0;

    seq[0] = decompressSymbol[getDecodedValue(_w, pos, 2)];
    qlt[0] = 0;
    pos += 2;

    len = getDecodedValue(_w, pos, 7) + 1;
    pos += 7;

    //fprintf(stderr, "ID "u64bitFMT" base %c/%d len "u32bitFMT"\n",
    //        _w[0], seq[0], seq[0], len);

    for (i=1; i<len; i++) {
      u64bit x = getDecodedValue(_w, pos, 7);
      seq[i] = decompressSymbol[x >> 5];
      qlt[i] = x & u64bitMASK(5);
      pos += 7;
    }

    seq[len] = 0;
    qlt[len] = 0;

    return(_w[0]);
  };

  u64bit   tagID(void) const     { return(_w[0]); };
  u64bit   length(void)          { return(getDecodedValue(_w, 66, 7)); };

private:
  friend class tapperTagFile;

  u64bit  _w[TAPPER_TAG_WORDS];
};





class tapperTagFileMetaData {
public:
  u32bit   tagSize(void)          { return(_minTagLen); };
  u32bit   isPairedTagFile(void)  { return(_isPaired);  };
  u32bit   mean(void)             { return(_mean);      };
  u32bit   stddev(void)           { return(_stddev);    };

private:
  friend class tapperTagFile;
  u32bit   _minTagLen;
  u32bit   _maxTagLen;
  u32bit   _tagWords;
  u32bit   _isPaired;
  u32bit   _isFragment;
  u32bit   _mean;
  u32bit   _stddev;
};




//  Notes:
//
//  1 Stores EITHER mated tags or fragment tags, NEVER both in the same file.
//  2 Variable tag size for every tag (even mated tags)
//  3 QVs too.
//  4 Random access
//  5 Reads are assigned a 64-bit UID



class tapperTagFile {
public:
  tapperTagFile(char *name) {
    _tagFile  = new recordFile(name, sizeof(tapperTagFileMetaData), TAPPER_TAG_WORDS * sizeof(u64bit), true);
    _metaData = (tapperTagFileMetaData *)_tagFile->header();

    if (_metaData->_tagWords == 0) {
      _metaData->_minTagLen  = ~u32bitZERO;
      _metaData->_maxTagLen  = 0;
      _metaData->_tagWords   = TAPPER_TAG_WORDS;
      _metaData->_isPaired   = 0;
      _metaData->_isFragment = 0;
      _metaData->_mean       = 0;
      _metaData->_stddev     = 0;
    }

    if (_metaData->_tagWords != TAPPER_TAG_WORDS) {
      fprintf(stderr, "tapperTagFile()--  ERROR!  Tag file was built with TAPPER_TAG_WORDS="u32bitFMT", but code has %d.\n",
              _metaData->_tagWords, TAPPER_TAG_WORDS);
      exit(1);
    }
  };

  ~tapperTagFile() {
    //  Metadata is updated automagically when tagFile is deleted.
    delete _tagFile;
  };

  tapperTagFileMetaData   *metaData(void) {
    return(_metaData);
  };

  void   put(tapperTag *tag) {
    u64bit len = tag->length();

    _metaData->_isFragment = 1;
    if (_metaData->_isPaired)
      fprintf(stderr, "tapperTagFile()--  ERROR: file contains mated tags, tried to pet a fragment tag.\n"), exit(1);

    if (len < _metaData->_minTagLen) _metaData->_minTagLen = len;
    if (_metaData->_minTagLen < len) _metaData->_maxTagLen = len;

    _tagFile->putRecord(tag->_w);
  };

  void   put(tapperTag *ta1, tapperTag *ta2) {
    u64bit len1 = ta1->length();
    u64bit len2 = ta2->length();

    _metaData->_isPaired  = 1;
    if (_metaData->_isFragment)
      fprintf(stderr, "tapperTagFile()--  ERROR: file contains fragment tags, tried to pet a mated tag.\n"), exit(1);

    if (len1 < _metaData->_minTagLen)  _metaData->_minTagLen = len1;
    if (_metaData->_minTagLen < len1)  _metaData->_maxTagLen = len1;

    if (len2 < _metaData->_minTagLen)  _metaData->_minTagLen = len2;
    if (_metaData->_minTagLen < len2)  _metaData->_maxTagLen = len2;

    _tagFile->putRecord(ta1->_w);
    _tagFile->putRecord(ta2->_w);
  };

  bool   get(tapperTag *tag) {
    if (_metaData->_isPaired == 1)
      fprintf(stderr, "tapperTagFile()--  ERROR: file contains mated tags, tried to get a fragment tag.\n"), exit(1);
    return(_tagFile->getRecord(tag->_w));
  };

  bool   get(tapperTag *ta1, tapperTag *ta2) {
    if (_metaData->_isFragment == 1)
      fprintf(stderr, "tapperTagFile()--  ERROR: file contains fragment tags, tried to get a mated tag.\n"), exit(1);
    return(_tagFile->getRecord(ta1->_w) &&
           _tagFile->getRecord(ta2->_w));
  };

private:
  tapperTagFileMetaData  *_metaData;
  recordFile             *_tagFile;
};



class seqStore {
public:
  seqStore();
  seqStore(const char *filename, seqStream *s);
  ~seqStore();

  unsigned char get(void);
  bool          eof(void)          { return(_thisBlock >= _numBlocks); };

  bool          rewind(void)       { return(setRange(_iterationBeg, _iterationEnd)); };

  u64bit        seqPos(void)       { return(_blockInfo[_thisBlock]._posInSeq + _thisBlockPosition - 1); };
  u64bit        seqIID(void)       { return(_blockInfo[_thisBlock]._seqIID); };

  u64bit        strPos(void)       { return(_blockInfo[_thisBlock]._posInBPF / 2 + _blockInfo[_thisBlock]._posInStrOff + _thisBlockPosition - 1); };
  u64bit        acgtPos(void)      { return(_blockInfo[_thisBlock]._posInBPF / 2 + _thisBlockPosition - 1); };

  u64bit        numberOfACGT(void) { return(_numberOfACGT); };


private:
  void   findBlock(u64bit  merNumber,
                   u32bit &block,
                   u64bit &totalMers,
                   u64bit &totalMersInFile);

  bool   seek(u64bit pos);

public:

  //  Set the range of ACGT sequence we will return.  Coordinates are
  //  space-based.  Example:
  //
  //      012...3...456
  //      AAAnnnCnnnGGG
  //
  //  range(0,0) -> nothing
  //  range(0,1) -> A
  //  range(0,5) -> AAAnnnCnnnG
  //
  bool                   setRange(u64bit beg, u64bit end) {
    _iterationBeg = beg;
    _iterationEnd = end;
    _iteration    = beg;
    return(seek(_iterationBeg));
  };


  //  Why does loadStore need the seqStream?  It really doesn't, but
  //  we can then check that the source file is older than our stream
  //  file.  Otherwise, it is all too easy, especially when testing
  //  :-), to have the wrong seqStore for the test.fasta you are
  //  trying to use.
  //
  bool   loadStore(const char *filename, seqStream *ss);
  void   buildStore(const char *filename, seqStream *ss);

  bool   loadStoreInCore(void) { return(_streamFile->loadInCore()); };


private:

  //  Limits:
  //   16 million -- length of a single block of data (efficiency only)
  //    1 billion -- length of a single sequence
  //   16 billion -- number of sequences in a single file
  // 1024 billion -- length of a single file
  //
  struct seqStoreBlock {
    u64bit      _isACGT:1;        // block is acgt
    u64bit      _seqIID:32;       // (4e9) iid of the sequence we are in
    u64bit      _posInSeq:30;     // (1e9) position in sequence

    //  Add this to _posInBPF to get the positionInStream
    u64bit      _posInStrOff:32;

    u64bit      _posInBPF:40;     // (1024e9) position in the bit file of sequence
    u64bit      _len:24;          // (16e6) length of block
  };

  u64bit                _numBlocks;
  seqStoreBlock        *_blockInfo;

  u64bit                _thisBlock;          //  Which _blockInfo we are in
  u64bit                _thisBlockPosition;  //  Position into the block

  u64bit                _numberOfACGT;

  bitPackedFile        *_streamFile;

  u64bit                _iterationBeg;
  u64bit                _iterationEnd;
  u64bit                _iteration;
};

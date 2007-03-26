#ifndef MERSTREAMFILE_H
#define MERSTREAMFILE_H

//  Uses the merStream to build a compressed merStreamFile.  The
//  merStreamFile contains all the forward mers that the stream
//  returns.  It does not contain position information or deflines.
//  Just the mers, ma'am.
//
//  To get position, sequence number and defline:
//
//  Each mer block can store the position start, and the sequence
//  number.  If we then store the deflines, we get everything!
//
//  N.B. because we create a merStream internally, and that uses
//  the seqStream to get sequence from the file, the sequence
//  file MUST BE SQUEEZED to get the correct positionInStream.

//  Returns true if the merStreamFile exists.
bool  merStreamFileExists(const char *i);

class merStream;

class merStreamFileBuilder {
public:
  merStreamFileBuilder(u32bit merSize, const char *inputFile, const char *outputFile);
  merStreamFileBuilder(merStream *MS, const char *outputFile);
  ~merStreamFileBuilder();

  //  Returns the number of mers written
  u64bit                 build(bool beVerbose=false);

private:
  u32bit                _merSize;
  chainedSequence      *_chainedSeq;
  merStream            *_merStream;
  const char           *_outputFile;
};






class merStreamFileReader {
public:
  merStreamFileReader(const char *filename, u32bit desiredMerSize=u32bitZERO);
  ~merStreamFileReader();

  u64bit                 loadInCore(void) { return(_streamFile->loadInCore()); };

  u32bit                 merSize(void)      { return(_merSizeDesired); };
  u64bit                 numberOfMers(void) { return(_numMers); };

  kMer const &           theFMer(void) { return(_fMer); };
  kMer const &           theRMer(void) { return(_rMer); };

  u64bit                 thePositionInSequence(void) { return(_posInSeq); };
  u64bit                 thePositionInStream(void)   { return(_posInStr); };
  u64bit                 theSequenceNumber(void)     { return(_sequence); };

private:
  void                   findBlock(u64bit merNumber,
                                   u32bit  &block,
                                   u64bit  &totalMers,
                                   u64bit  &totalMersInFile);
  bool                   seekToMer(u64bit merNumber);
  bool                   flipToCanonicalize(u64bit merNum);

public:
  //  Useful to bigmeryl or whatever it ends up being called.
  //

  //  Returns true if the mer at merNum1 is smaller than the mer
  //  at merNum2.  Be aware that unless your merstream is in core,
  //  there are serious performance hits here!
  //
  bool                   smaller(u64bit merNum1, u64bit merNum2, bool doCanonical);

public:
  //  This will limit the number of mers returned to 'limit'.  The
  //  seek*() methods do not reset the count or the limit.  The
  //  current count is reset when this function is called.  Set to
  //  ~u64bitZERO for an infinite limit.  Zero really does mean zero,
  //  not infinity.
  //
  bool                   setIterationStart(u64bit s) {
    //fprintf(stderr, "merStreamFileReader::setInterationStart()-- Set to "u64bitFMT"\n", s);
    _iterationStart = s;
    _iteration      = u64bitZERO;
    return(seekToMer(_iterationStart));
  };
  void                   setIterationLimit(u64bit l) {
    //fprintf(stderr, "merStreamFileReader::setInterationLimit()-- Set to "u64bitFMT"\n", l);
    _iterationLimit = l;
    _iteration      = u64bitZERO;
  };

  //  Rewind to the start of the iteration.
  //
  bool                   rewind(void) {
    //fprintf(stderr, "merStreamFileReader::rewind()--\n");
    return(setIterationStart(_iterationStart));
  };



  bool                   nextMer(u32bit skip=0);
private:
  const char           *_inputFile;
  bitPackedFile        *_streamFile;
  u32bit                _merSizeInFile;
  u32bit                _merSizeDesired;
  u64bit                _numMers;
  u64bit                _numBlocks;
  off_t                 _blockFileSize;
  off_t                 _streamFileSize;
  off_t                 _blkStart;
  off_t                 _strStart;

  //  When used for bigmeryl, we'll have 32 million (to 45 million)
  //  reads.  Reads should have no non-acgt, so each read should be in
  //  exactly one block.  Per block, we want to store
  //    32 bits - block size
  //    32 bits - block sequence number
  //    64 bits - block position in the sequence
  //    64 bits - block position in the stream
  //    64 bits - first mer number in this block
  //    32 bits - correction to the first mer number to get the total mers in the file so far
  //   ---
  //   256 bits (32 bytes) -- so, 1GB data for 32 million reads.
  //
  u32bit               *_blockSize;
  u32bit               *_blockSequence;
  u64bit               *_blockPosInSeq;
  u64bit               *_blockPosInStr;
  u64bit               *_blockFirstMerNum;
  u32bit               *_blockTotalMerCorrection;

  //  Things used to decode the stream
  //
  u64bit                _thisBlock;
  u64bit                _thisBlockSize;
  u64bit                _firstMer;

  //  Things used to stop the stream
  //
  u64bit                _iterationStart;
  u64bit                _iterationLimit;
  u64bit                _iteration;

  //  Things returned by the accessors, set/updated by nextMer().
  //
  kMer                  _fMer;
  kMer                  _rMer;
  u64bit                _posInSeq;
  u64bit                _posInStr;
  u64bit                _sequence;
};


#endif  //  MERSTREAMFILE_H

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
//  the FastAstream to get sequence from the file, the sequence
//  file MUST BE SQUEEZED to get the correct positionInStream.

//  Returns true if the merStreamFile exists.
bool  merStreamFileExists(const char *i);

//  Defline support in the merStreamFile is very broken -- we possibly
//  construct it correctly, but we don't access it at all.  Nobody is
//  currently using it, so the functionality has been completely
//  disabled (which will save some space in the file).
//
//  It's not been removed from the class, we just don't save the
//  deflines to the file.
//
//#define ENABLE_DEFLINE_SUPPORT


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
  FastAstream          *_fastaStream;
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

#ifdef ENABLE_DEFLINE_SUPPORT
  char const            *theDefLine(void)            { return(_defline);  };
#else
  char const            *theDefLine(void)            { return(">nodefline");  };
#endif

private:
  //  These used to be public, but seekToMer() has been replaced by
  //  setIterationStart(), and nobody was using seekToSequence().
  //
  bool                   seekToMer(u64bit merNumber);
  bool                   seekToSequence(u64bit seqNumber);

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
  u64bit                _numDefs;
  u64bit                _defLength;
  off_t                 _blockFileSize;
  off_t                 _deflineFileSize;
  off_t                 _streamFileSize;
  off_t                 _blkStart;
  off_t                 _defStart;
  off_t                 _strStart;

  u32bit               *_blockSize;
  u32bit               *_blockSequence;
  u64bit               *_blockPosInSeq;
  u64bit               *_blockPosInStr;
  char                 *_deflineStorage;
  char                **_deflines;

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

  char                 *_defline;
};


#endif  //  MERSTREAMFILE_H

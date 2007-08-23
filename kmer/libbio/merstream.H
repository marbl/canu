#ifndef MERSTREAM_H
#define MERSTREAM_H

//  This is the merStream.  Given a file of sequence and deflines, or
//  a string of sequence, it returns mers.
//
//  The file based merstream uses a fastastream to get sequence, and
//  that doesn't care if there is whitespace or deflines in the
//  stream.  It returns only valid sequence.
//
//  If the sequence/length method is used, we assume that 'seq'
//  contains no whitespace, and no extra characters (e.g., deflines),
//  just pure sequence.  It will work if you give it deflines and
//  whitespcae, but positions will be wrong.
//
//  **********
//
//  N.B. white-space in the source (applies to seqStream and char*
//  only) result in an incorrect positionInStream().  If this troubles
//  you, use a seqStream instead of a seqStream, or squeeze the
//  sequences.
//
//  The merstream has many modes of operation, depending
//  on what the input format is.  We could:
//  1)  subclass an interface template
//  2)  use function pointers set at create time
//  3)  use a set of if tests in each call that cares
//
//  1 and 2 are essentially the same, except 2 can use human
//  intelligence to get rid of the function call on a few methods
//  that aren't input specific.
//
//  One can claim that 3 is better because it doesn't need
//  to mess with the stack, and if inlined, branch prediction
//  will remove the overhead of the if tests.

class merStream {
private:
  void initialize(u32bit merSize);
public:
  merStream(u32bit merSize, seqStore *ss);
  merStream(u32bit merSize, seqStream *cs);
  merStream(u32bit merSize, seqInCore *seq, u32bit offset=0, u32bit length=0);
  merStream(u32bit merSize, char *seq, u32bit offset=0, u32bit length=0);
  ~merStream();

  u32bit                 theMerSize(void) { return(_merSize); };

  kMer const &           theFMer(void)    { return(_fMer); };
  kMer const &           theRMer(void)    { return(_rMer); };
  kMer const &           theCMer(void)    { return((_fMer < _rMer) ? _fMer : _rMer); };

  bool                   nextMer(u32bit skip=0);
  bool                   rewind(void);

  //  'pos' is the N'th ACGT in the file, not the N'th mer in the file!
  bool                   seek(u64bit pos) {
    if (_ss_file) {
      _ss_file->setIterationStart(pos);
      loadMer(_merSize - 1);
      return(nextMer());
    }
    return(false);
  };

  u64bit                 thePositionInSequence(void);
  u64bit                 thePositionInStream(void);
  u64bit                 theSequenceNumber(void);
  u64bit                 approximateNumberOfMers(void);

private:
  void                   loadMer(u32bit s);
  u64bit                 nextSym(void);

private:
  seqStore             *_ss_file;
  seqStream            *_cs_stream;
  char const           *_st_string;
  u64bit                _st_stringPos;
  u64bit                _st_stringLen;
  u64bit                _st_posInSeq;
  u64bit                _st_posInStr;
  u32bit                _st_num;

  u32bit                _merSize;
  s64bit                _timeUntilValid;

  kMer                  _fMer;
  kMer                  _rMer;
};




//  Why the +1?  seqStream->seqPos() returns the position that the
//  last character retrieved is at.  We need to go back _merSize-1
//  positions to get to the start of the mer.

inline
u64bit
merStream::thePositionInSequence(void) {
  if      (_ss_file)
    return(_ss_file->seqPos() + 1 - _merSize);
  else if (_cs_stream)
    return(_cs_stream->seqPos() + 1 - _merSize);
  else
    return(_st_posInSeq - _merSize);
}

inline
u64bit
merStream::thePositionInStream(void) {
  if      (_ss_file)
    return(_ss_file->strPos() + 1 - _merSize);
  else if (_cs_stream)
    return(_cs_stream->strPos() + 1 - _merSize);
  else
    return(_st_stringPos - _merSize);
}

inline
u64bit
merStream::theSequenceNumber(void) {
  if      (_ss_file)
    return(_ss_file->seqIID());
  else if (_cs_stream)
    return(_cs_stream->seqIID());
  else
    return(_st_num);
}

//  Return a guess on the number of mers in the input.  Our guess will
//  never be too small.  Returns ~u64bitZERO if no guess is possible.
//
inline
u64bit
merStream::approximateNumberOfMers(void) {
  if (_ss_file)
    return(_ss_file->numberOfACGT());
  if (_cs_stream)
    return(_cs_stream->lengthOfSequences());
  if (_st_string)
    return(_st_stringLen);
  return(~u64bitZERO);
}


inline
u64bit
merStream::nextSym(void) {

  //  This is called only for seqStreams, seqStream and char* input.

  if (_ss_file)
    return(_ss_file->get());

  if (_cs_stream)
    return(_cs_stream->get());

  //  See also the constructor for character strings -- it skips
  //  white space and the first defline.
  
  //  Skip any whitespace
  //
  while (whitespaceSymbol[_st_string[_st_stringPos]] &&
         (_st_string[_st_stringPos]) &&
         (_st_stringPos < _st_stringLen))
    _st_stringPos++;

  //  Are we at a defline?  If so, zip along until we hit a
  //  newline.  Then we return the code (fastastream) for a
  //  sequence break.
  //
  if (_st_string[_st_stringPos] == '>') {
    while ((_st_string[_st_stringPos] != '\n') &&
           (_st_string[_st_stringPos]) &&
           (_st_stringPos < _st_stringLen))
      _st_stringPos++;

    _st_posInSeq = 0;
    _st_posInStr++;

    _st_num++;
    return(254);
  }

  //  Are we at the EOF?
  //
  if ((_st_string[_st_stringPos] == 0) || (_st_stringPos >= _st_stringLen))
    return(0);

  //  Update the position of the character in this sequence, and
  //  the positon in the stream.
  //
  _st_posInSeq++;
  _st_posInStr++;

  return(_st_string[_st_stringPos++]);
}




//  We have a problem; at the start of a stream, we want to
//  initialize the mer with merSize-1 bases.  In the middle of the
//  stream, we need to load the mer with merSize bases.
//
//  loadMer will push s bases onto the mer, restarting if it hits a
//  mer break.
//
//  No masking of the mer is performed.
//
inline
void
merStream::loadMer(u32bit s) {
  u64bit   ch = 255;
  u64bit   cf = 0;
  u64bit   cr = 0;

  _fMer.clear();
  _rMer.clear();

  _timeUntilValid = s;

  //  While we are invalid, and still in the sequence
  //  push characters onto the mer.  The valid
  //  count is updated if we hit an invalid base.
  //
  while ((_timeUntilValid > 0) && (ch != 0)) {
    ch = nextSym();
    cf = validCompressedSymbol[ch];

    //  Rather than take the chance of generating a cache miss accessing
    //  another array, we just mask out the upper bits of validCompressedSymbol[];
    //  this is exactly the same as using compressSymbol[].
    //
    //  We need to mask the upper bits for reverse (but not for forward)
    //  because, in reverse, we shift to the right.  If we don't mask theseq
    //  out, we will have extra bits set in the mer.
    //
    //  Example: Consider placing 255 (== 11111111, the invalid symbol returned
    //  from validCompressedSymbol) at the fourth base in a four mer:
    //    000000**xxxxxx -- the ** are is the fourth base.
    //
    //  Without masking, we'd set all the 0's and all the *'s to one, in effect,
    //  preloading the next three bases of the mer.
    //
    cr = validCompressedSymbol[complementSymbol[ch]] & 0x03;

    _timeUntilValid--;

    //  If we get a break, reset the counter a re-clear the mer
    //
    if (cf & 0xfc) {
      _timeUntilValid = s;
      _fMer.clear();
      _rMer.clear();
    }

    //  If the ch is valid, we can obviously add it to the mer.
    //  If it is invalid, we don't care; by the time the mer is valid
    //  all bits of any invalid mers are removed.  Sure, we could
    //  put this into the else, but I suspect that it will be faster outside.
    //
    _fMer += cf;
    _rMer -= cr;
  }
}


inline
bool
merStream::nextMer(u32bit skip) {

  do {
    u64bit  ch = nextSym();
    u64bit  cf = validCompressedSymbol[ch];
    u64bit  cr = validCompressedSymbol[complementSymbol[ch]] & 0x03;

    //  EOF?
    if (ch == 0)
      return(false);

    //  Push the ch onto the mer.
    //
    _fMer += cf;
    _rMer -= cr;

    //  If the ch is invalid, we need to make a whole new mer.
    //  This also implies that we're done skipping.
    //
    if (cf & (unsigned char)0xfc) {
      loadMer(_merSize);
      skip = 0;
    }
  } while (skip--);

  _fMer.mask(true);
  _rMer.mask(false);

  //  Still need to check if we are valid -- we could
  //  have run off the end of the sequence before a valid
  //  mer was created.
  //
  return(_timeUntilValid <= 0);
}


#endif  //  MERSTREAM_H


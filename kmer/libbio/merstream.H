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
  void initialize(kMerBuilder *kb);
public:
  merStream(kMerBuilder *kb, seqStore *ss);
  merStream(kMerBuilder *kb, seqStream *cs);
  merStream(kMerBuilder *kb, seqInCore *seq, u32bit offset=0, u32bit length=0);
  merStream(kMerBuilder *kb, char *seq, u32bit offset=0, u32bit length=0);
  ~merStream();

  kMer const &           theFMer(void)    { return(_kb->theFMer()); };
  kMer const &           theRMer(void)    { return(_kb->theRMer()); };
  kMer const &           theCMer(void)    { return(_kb->theCMer()); };

  bool                   nextMer(u32bit skip=0);

  //  Everyone can rewind back to the start, and it's boring.
  //  Slightly more efficient than just opening a new stream.
  //
  //  If we're backed by a seqStore, we can reposition the merstream
  //  anywhere we like to.  The positions refer to the N'th ACGT in
  //  the file, not the N'th mer in the file!  Note that rewind()ing a
  //  file with a range set rewinds to the beginning of the _range_
  //  not the file.
  //
  bool                   rewind(void);
  bool                   setRange(u64bit beg, u64bit end);

  u64bit                 thePositionInSequence(void);
  u64bit                 thePositionInStream(void);
  u64bit                 theSequenceNumber(void);
  u64bit                 approximateNumberOfMers(void);

private:
  u64bit                 nextSym(void);

private:
  kMerBuilder          *_kb;

  seqStore             *_ss_file;
  u64bit                _ss_beg;
  u64bit                _ss_end;
  u64bit                _ss_pos;

  seqStream            *_cs_stream;

  char const           *_st_string;
  u64bit                _st_stringPos;
  u64bit                _st_stringLen;
  u64bit                _st_posInSeq;
  u64bit                _st_posInStr;
  u32bit                _st_num;
};



//  Why the +1?  seqStream->seqPos() returns the position that the
//  last character retrieved is at.  We need to go back getMerSpan()
//  positions to get to the start of the mer.

inline
u64bit
merStream::thePositionInSequence(void) {
  if      (_ss_file)
    return(_ss_file->seqPos() + 1 - theFMer().getMerSpan());
  else if (_cs_stream)
    return(_cs_stream->seqPos() + 1 - theFMer().getMerSpan());
  else
    return(_st_posInSeq - theFMer().getMerSpan());
}

inline
u64bit
merStream::thePositionInStream(void) {
  if      (_ss_file)
    return(_ss_file->strPos() + 1 - theFMer().getMerSpan());
  else if (_cs_stream)
    return(_cs_stream->strPos() + 1 - theFMer().getMerSpan());
  else
    return(_st_stringPos - theFMer().getMerSpan());
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
  u64bit  approx = ~u64bitZERO;

  if        (_ss_file) {
    approx = _ss_file->numberOfACGT();
    if (_ss_end - _ss_beg < approx)
      approx = _ss_end - _ss_beg;
  } else if (_cs_stream) {
    approx = _cs_stream->lengthOfSequences();
  } else if (_st_string) {
    approx = _st_stringLen;
  }

  return(approx);
}


inline
u64bit
merStream::nextSym(void) {

  if (_ss_file)
    return(_ss_file->get());

  if (_cs_stream)
    return(_cs_stream->get());

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

  if ((_st_string[_st_stringPos] == 0) || (_st_stringPos >= _st_stringLen))
    //  EOF
    return(0);

  //  Update the position of the character in this sequence, and
  //  the positon in the stream.
  //
  _st_posInSeq++;
  _st_posInStr++;

  return(_st_string[_st_stringPos++]);
}



inline
bool
merStream::nextMer(u32bit skip) {
  bool    needMore = true;

  do {
    char ch = nextSym();
    if (ch == 0)
      return(false);

    needMore = _kb->addBase(ch);
  } while ((needMore == true) || (skip-- > 0));

  _kb->mask();

  //  Compute if we've hit our iteration limit.
  //
  if (_ss_file) {
    u64bit acgtp = _ss_file->acgtPosForMerStream(theFMer().getMerSpan()) + 1;

    //fprintf(stderr, "acgtPos: "u64bitFMT";  _ss_end "u64bitFMT"\n", acgtp, _ss_end);

    if (_ss_end <= acgtp)
      return(false);
  }

  return(true);
}

#endif  //  MERSTREAM_H

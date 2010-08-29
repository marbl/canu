#ifndef SEQSTREAM_H
#define SEQSTREAM_H

#include "util++.H"
#include "bio++.H"

#include "seqFile.H"

struct seqStreamIndex {
  u32bit   _iid;  //  seqFile IID
  u32bit   _len;  //  length of the sequence
  u64bit   _bgn;  //  begin position in the stream
};


class seqStream {
public:
  seqStream(const char *filename);
  seqStream(char *sequence, u32bit length);
  ~seqStream();

  //  Separate sequences with this letter.  Non-ACGT is always
  //  returned as 'N'.  Changing the length of the separator AFTER
  //  setting the range will result in the wrong range being used.
  //
  void              setSeparator(char sep, u32bit len);

  //  get() returns one letter per input letter -- a gap of size n
  //  will return n gap symbols.
  //
  unsigned char     get(void);
  bool              eof(void)        { return(_eof); };

  //  Returns to the start of the range.
  //
  void              rewind(void);

  //  Set the range of ACGT sequence we will return.  Coordinates are
  //  space-based.  Example:
  //
  //  >one
  //  AAA
  //  >two
  //  C
  //  >three
  //  GGG
  //
  //  We separate these sequences with three '-' letters.
  //
  //  strPos    012...3...456
  //            AAA---C---GGG
  //
  //  range(0,0) -> nothing
  //  range(0,1) -> A
  //  range(0,3) -> AAA
  //  range(0,4) -> AAAnnnC
  //  range(0,5) -> AAAnnnCnnnG
  //
  void              setRange(u64bit bgn, u64bit end);
  void              setPosition(u64bit pos);

  //  seqPos() is the position we are at in the current sequence;
  //  seqIID() is the iid of that sequence;
  //  strPos() is the position we are at in the chained sequence
  //
  //  Values are not defined if the letter is a separator.
  //
  u32bit            seqPos(void)     { return(_currentPos);            };
  u32bit            seqIID(void)     { return(_idx[_currentIdx]._iid); };
  u64bit            strPos(void)     { return(_streamPos);             };

  u32bit            numberOfSequences(void) { return(_idxLen); };

  //  Return the length of, position of (in the chain) and IID of the
  //  (s)th sequence in the chain.
  //
  u32bit            lengthOf(u32bit s) { return((s >= _idxLen) ? ~u32bitZERO : _idx[s]._len); };
  u32bit            IIDOf(u32bit s)    { return((s >= _idxLen) ? ~u32bitZERO : _idx[s]._iid); };
  u64bit            startOf(u32bit s)  { return((s >= _idxLen) ? ~u64bitZERO : _idx[s]._bgn); };

  //  For a chain position p, returns the s (above) for that position.
  //
  u32bit            sequenceNumberOfPosition(u64bit p);
  void              tradeSpaceForTime(void);

private:
  void              fillBuffer(void);

  seqFile         *_file;        //  Backed by a seqFile.
  char            *_string;      //  Backed by a character string.

  u64bit           _bgn;         //  Begin/End position in chained sequence
  u64bit           _end;

  u32bit           _currentIdx;  //  index into _idx of the current sequence
  u32bit           _currentPos;  //  position in the current sequence
  u64bit           _streamPos;   //  position in the chained sequence

  // Buffer for holding sequence from the seqFile.

  u32bit           _bufferMax;   //  Max size of the buffer
  u32bit           _bufferLen;   //  Actual size of the buffer
  u32bit           _bufferPos;   //  Position we are at in the buffer
  u32bit           _bufferSep;   //  How much of the buffer is separator
  char            *_buffer;

  //  Info about the raw sequences

  u32bit           _idxLen;
  seqStreamIndex  *_idx;

  u32bit          *_seqNumOfPos;

  u64bit           _lengthOfSequences;

  bool             _eof;

  char             _separator;
  u32bit           _separatorLength;
};



#endif  //  SEQSTREAM_H

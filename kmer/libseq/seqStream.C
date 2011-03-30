#include "seqFactory.H"
#include "seqStream.H"


seqStream::seqStream(const char *filename) {
  _file              = openSeqFile(filename);
  _string            = 0L;

  _currentIdx        = 0;
  _currentPos        = 0;
  _streamPos         = 0;

  _bufferMax         = 1048576;
  _bufferLen         = 0;
  _bufferPos         = 0;
  _bufferSep         = 0;
  _buffer            = new char [_bufferMax + 1];

  _idxLen            = _file->getNumberOfSequences();
  _idx               = new seqStreamIndex [_idxLen + 1];

  //fprintf(stderr, "seqStream::seqStream()--  Allocating "u64bitFMT"MB for seqStreamIndex on "u64bitFMT" sequences.\n",
  //        _idxLen * sizeof(seqStreamIndex) / 1024 / 1024, _idxLen);

  _seqNumOfPos       = 0L;

  _lengthOfSequences = 0;

  _eof               = false;

  _separator         = '.';
  _separatorLength   = 2;

  setSeparator('.', 2);

  _bgn               = 0;
  _end               = _lengthOfSequences;
}



seqStream::seqStream(char *sequence, u32bit length) {
  _file              = 0L;
  _string            = sequence;

  _currentIdx        = 0;
  _currentPos        = 0;
  _streamPos         = 0;

  _bufferMax         = length;
  _bufferLen         = length;
  _bufferPos         = 0;
  _bufferSep         = 0;
  _buffer            = _string;

  _idxLen            = 1;
  _idx               = new seqStreamIndex [_idxLen + 1];

  _seqNumOfPos       = 0L;

  _idx[0]._iid = 0;
  _idx[0]._len = length;
  _idx[0]._bgn = 0;

  _idx[1]._iid = ~u32bitZERO;
  _idx[1]._len = 0;
  _idx[1]._bgn = length;

  _lengthOfSequences = length;

  _eof               = false;

  _separator         = '.';
  _separatorLength   = 20;

  _bgn               = 0;
  _end               = length;
}



seqStream::~seqStream() {
  if (_file) {
    delete    _file;
    delete [] _buffer;
  }
  delete [] _idx;
  delete [] _seqNumOfPos;
}



void
seqStream::setSeparator(char sep, u32bit len) {

  //  Special case; no separator needed for string backed sequences.
  if (_string)
    return;

  //  Bizarre signedness issue with sep=255
  //    ST->get() == sep        FAILS
  //    x=ST->get(); x == sep   SUCCEEDS
  //
  //  Not suggested to use non-printable ascii.

  if ((isprint(sep) == 0) || (tolower(sep) == 'a') || (tolower(sep) == 'c') || (tolower(sep) == 'g') || (tolower(sep) == 't')) {
    fprintf(stderr, "seqStream::setSeparator()-- ERROR!  Separator letter must be printable ASCII and not [ACGTacgt].\n");
    exit(1);
  }
  if (len == 0) {
    fprintf(stderr, "seqStream::setSeparator()-- ERROR!  Separator length cannot be zero.\n");
    exit(1);
  }

  _lengthOfSequences = 0;

  _separator       = sep;
  _separatorLength = len;;

  for (u32bit s=0; s<_idxLen; s++) {
    _idx[s]._iid = s;
    _idx[s]._len = _file->getSequenceLength(s);
    _idx[s]._bgn = _lengthOfSequences;

    _lengthOfSequences += _idx[s]._len;
  }

  _idx[_idxLen]._iid = ~u32bitZERO;
  _idx[_idxLen]._len = 0;
  _idx[_idxLen]._bgn = _lengthOfSequences;

  //  Rebuild our sequence number of position map, if it exists.
  //
  if (_seqNumOfPos) {
    delete [] _seqNumOfPos;
    tradeSpaceForTime();
  }
}



void
seqStream::tradeSpaceForTime(void) {
  u32bit  i = 0;
  u32bit  s = 0;

  //fprintf(stderr, "Allocating "u32bitFMT" u32bits for seqNumOfPos.\n", _lengthOfSequences);

  _seqNumOfPos = new u32bit [_lengthOfSequences];

  for (i=0; i<_lengthOfSequences; i++) {

    //  Increment the sequence number until we enter into the next
    //  sequence.  Zero length sequences require the use of a 'while'
    //  here.
    //
    while (i >= _idx[s+1]._bgn)
      s++;

    _seqNumOfPos[i] = s;
  }
}



unsigned char
seqStream::get(void) {
  if (_streamPos >= _end)
    _eof = true;
  if ((_eof == false) && (_bufferPos >= _bufferLen))
    fillBuffer();
  if (_eof)
    return(0);
  if (_bufferSep == 0) {
    _currentPos++;
    _streamPos++;
  } else {
    _bufferSep--;
  }
  return(_buffer[_bufferPos++]);
}



void
seqStream::rewind(void){

  //  Search for the correct spot.  Uncommon operation, be inefficient
  //  but simple.  The range was checked to be good by setRange().

  u32bit s = 0;
  u64bit l = 0;

  while ((s < _idxLen) && (l + _idx[s]._len < _bgn))
    l += _idx[s++]._len;

  _eof = false;

  //  (_bgn - l) is a 32-bit quanitity because of the second half of
  //  the while above.  Although _bgn is a 64-bit value, the value
  //  used to set _bufferPos will be for that of a string constructor,
  //  and so _bgn will be 32-bits.  fillBuffer() resets _bufferPos if
  //  we're backed by a file.

  _currentIdx = s;
  _currentPos = _bgn - l;
  _streamPos  = _bgn;
  _bufferPos  = _bgn;

  //fprintf(stderr, "seqStream::rewind()-- 1 currentIdx="u32bitFMT" currentPos="u32bitFMT" streamPos="u32bitFMT" bufferPos="u32bitFMT"\n",
  //        _currentIdx, _currentPos, _streamPos, _bufferPos);

  fillBuffer();

  //fprintf(stderr, "seqStream::rewind()-- 2 currentIdx="u32bitFMT" currentPos="u32bitFMT" streamPos="u32bitFMT" bufferPos="u32bitFMT"\n",
  //        _currentIdx, _currentPos, _streamPos, _bufferPos);
}



void
seqStream::setRange(u64bit bgn, u64bit end) {

  assert(bgn < end);

  u32bit s = 0;
  u64bit l = 0;

  while (s < _idxLen)
    l += _idx[s++]._len;

  if (end == ~u64bitZERO)
    end = l;

  if ((bgn > l) || (end > l))
    fprintf(stderr, "seqStream::setRange()-- ERROR: range ("u64bitFMT","u64bitFMT") too big; only "u64bitFMT" positions.\n",
            bgn, end, l), exit(1);

  _bgn = bgn;
  _end = end;

  rewind();
}


void
seqStream::setPosition(u64bit pos) {

  assert(_bgn <=  pos);
  assert( pos <  _end);

  u64bit old = _bgn;

  _bgn = pos;
  rewind();
  _bgn = old;
}


u32bit
seqStream::sequenceNumberOfPosition(u64bit p) {
  u32bit   s = ~u32bitZERO;

  //  binary search on our list of start positions, to find the
  //  sequence that p is in.

  if (_lengthOfSequences <= p) {
    fprintf(stderr, "seqStream::sequenceNumberOfPosition()-- WARNING: position p="u64bitFMT" too big; only "u64bitFMT" positions.\n",
            p, _lengthOfSequences);
    return(s);
  }

  if (_seqNumOfPos)
    return(_seqNumOfPos[p]);

  if (_idxLen < 16) {
    for (s=0; s<_idxLen; s++)
      if ((_idx[s]._bgn <= p) && (p < _idx[s+1]._bgn))
        break;
  } else {
    u32bit  lo = 0;
    u32bit  hi = _idxLen;
    u32bit  md = 0;

    while (lo <= hi) {
      md = (lo + hi) / 2;

      if        (p < _idx[md]._bgn) {
        //  This block starts after the one we're looking for.  
        hi = md;

      } else if ((_idx[md]._bgn <= p) && (p < _idx[md+1]._bgn)) {
        //  Got it!
        lo = md + 1;
        hi = md;
        s  = md;

      } else {
        //  By default, then, the block is too low.
        lo = md;
      }
    }
  }

  return(s);
}



void
seqStream::fillBuffer(void) {

  //  Special case for when we're backed by a character string; there
  //  is no need to fill the buffer.
  //
  if (_file == 0L) {
    if (_currentPos >= _end)
      _eof = true;
    return;
  }

  //  Read bytes from the _file, stuff them into the buffer.  Assumes
  //  there is nothing in the buffer to save.

  _bufferLen = 0;
  _bufferPos = 0;

  //  Still more stuff in the sequence?  Get it.

  if (_currentPos < _idx[_currentIdx]._len) {
#ifdef DEBUG
    fprintf(stderr, "seqStream::fillBuffer()--  More Seq currentPos="u32bitFMT" len="u32bitFMT"\n", _currentPos, _idx[_currentIdx]._len);
#endif
    _bufferLen = MIN(_idx[_currentIdx]._len - _currentPos, _bufferMax);

    if (_file->getSequence(_idx[_currentIdx]._iid,
                           _currentPos,
                           _currentPos + _bufferLen,
                           _buffer) == false)
      fprintf(stderr, "seqStream::fillBuffer()-- Failed to getSequence(part) #1 iid="u32bitFMT" bgn="u32bitFMT" end="u32bitFMT"\n",
              _idx[_currentIdx]._iid, _currentPos, _currentPos + _bufferLen), exit(1);

    return;
  }

  //  We've finished a sequence.  Load the next.

  _currentPos = 0;
  _currentIdx++;

  while ((_currentIdx < _idxLen) && (_idx[_currentIdx]._len == 0))
    _currentIdx++;

#ifdef DEBUG
  fprintf(stderr, "seqStream::fillBuffer()--  New Seq currentPos="u32bitFMT" len="u32bitFMT"\n", _currentPos, _idx[_currentIdx]._len);
#endif

  //  All done if there is no more sequence.

  if (_currentIdx >= _idxLen) {
    _eof = true;
    return;
  }

  //  Insert a separator.

  for (_bufferLen = 0; _bufferLen < _separatorLength; _bufferLen++)
    _buffer[_bufferLen] = _separator;

  //  Keep track of the separator - this is used to make sure we don't
  //  advance the sequence/stream position while the separator is
  //  being returned.
  //
  _bufferSep = _bufferLen;

  //  How much to get; minimum of what is left in the sequence, and
  //  the buffer size.  Don't forget about the separator we already
  //  inserted!
  //
  u32bit bl = MIN(_idx[_currentIdx]._len - _currentPos, _bufferMax - _bufferLen);

  if (_file->getSequence(_idx[_currentIdx]._iid,
                         _currentPos,
                         _currentPos + bl,
                         _buffer + _bufferLen) == false)
    fprintf(stderr, "seqStream::fillBuffer()-- Failed to getSequence(part) #2 iid="u32bitFMT" bgn="u32bitFMT" end="u32bitFMT"\n",
            _idx[_currentIdx]._iid, _currentPos, _currentPos + bl), exit(1);

  _bufferLen += bl;

  //  Load more, until buffer is full.  Not really needed, and won't
  //  improve performance much.  AND it adds a lot of complexity to
  //  track which sequence is current (_currentIdx).

  return;
}

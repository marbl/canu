#include "bio++.H"


FastASequenceOnDisk::FastASequenceOnDisk(char const *filename,
                                         _idxfa_pos hdrstart, _idxfa_len hdrlen,
                                         _idxfa_pos seqstart, _idxfa_len seqlen,
                                         IID_t iid,
                                         bool isSqueezed,
                                         bool isFixedWidth,
                                         _idxfa_len sl,
                                         _idxfa_len ss) {

  _idx               = iid;

  _headerLength      = hdrlen;
  _headerStart       = hdrstart;

  _sequenceLength    = seqlen;
  _sequenceStart     = seqstart;

  u32bit readBufferLength = 1024 * 1024;
  if (readBufferLength > _headerLength + _sequenceLength + 16)
    readBufferLength = _headerLength + _sequenceLength + 16;

  _readBuffer        = new readBuffer(filename, readBufferLength);

  _readBuffer->seek(_headerStart);

  _header = new char [_headerLength + 1];
  _readBuffer->read(_header, _headerLength);
  _header[_headerLength] = 0;

  _readBuffer->seek(_sequenceStart);

  _sequencePosition  = 0;

  if      (isSqueezed)
    _sourceType   = 0;
  else if (isFixedWidth)
    _sourceType   = 1;
  else
    _sourceType   = 2;

  _lineLength = sl;
  _lineSep    = ss;
}


FastASequenceOnDisk::~FastASequenceOnDisk() {
  delete    _readBuffer;
  delete [] _header;
}



char*
FastASequenceOnDisk::getChars(char *block, u32bit position, u32bit length) {

  //  Check that there is enough sequence left!
  //
  if (position >= _sequenceLength) {
    fprintf(stderr, "FastASequenceOnDisk::getChars()-- position "u32bitFMT" not in sequence of length "u32bitFMT".\n",
            position, length);
    return(0L);
  }

  if (position + length > _sequenceLength) {
    fprintf(stderr, "FastASequenceOnDisk::getChars()-- requested more sequence ("u32bitFMT" letters @ position "u32bitFMT")"
                    "than available (length "u32bitFMT"); length trimmed.\n",
            position, length, _sequenceLength);
    length = _sequenceLength - position;
  }

  if (block == 0L)
    block = new char [length + 1];

  switch (_sourceType) {
    case 0:
      //  Hooray!  A squeezed sequence!
      _readBuffer->seek(_sequenceStart + position);
      _readBuffer->read(block, length);
      break;
    case 1:
      //  Better than nothing.  We can compute where we should seek
      //  for fixed width sequences
      //
      //  In theory, we could use knowledge of the _lineLength to read
      //  blocks into our buffer, but it's much, much easier to just
      //  go character by character.

      _readBuffer->seek(_sequenceStart + (position / _lineLength) * _lineSep + position);

      for (u32bit i=0; i<length; i++) {
        block[i] = _readBuffer->get();
        _readBuffer->next();    
        while (whitespaceSymbol[(int)_readBuffer->get()])
          _readBuffer->next();    
      }
      break;
    case 2:
      //  Dang.  Gotta do a bunch of next()'s to get to the right
      //  spot, then do the same thing copying into our buffer.  We
      //  use a slightly optimized version of next(), some of the
      //  checks in next() are useless here.
      //
      _readBuffer->seek(_sequenceStart);
      for (u32bit i=0; i<position; i++) {
        _readBuffer->next();    
        while (whitespaceSymbol[(int)_readBuffer->get()])
          _readBuffer->next();    
      }

      for (u32bit i=0; i<length; i++) {
        block[i] = _readBuffer->get();
        _readBuffer->next();    
        while (whitespaceSymbol[(int)_readBuffer->get()])
          _readBuffer->next();    
      }
      break;
  }

  block[length] = 0;

  return(block);
}



FastASequenceOnDisk*
FastAWrapper::getSequenceOnDisk(void) {

  if (_currentSequenceNumber >= _theGlobalDesc._numberOfSequences)
    return(0L);

#if 0
  //  XXX: The current implementation of FastASequenceOnDisk requires
  //  that the file be squeezed.  This should be fixed later -- the
  //  interface doesn't require squeezedness.
  //
  if (!isRandomAccess() || !isSqueezed()) {
    fprintf(stderr, "FastAWrapper::getSequenceOnDisk()-- Asked for a sequence, but file '%s' is not\n", getSourceName());
    if (!_isRandomAccess)
      fprintf(stderr, "FastAWrapper::getSequenceOnDisk()--   Randomly accessable (no index)\n");
    if (!_theGlobalDesc._squeezedSequences)
      fprintf(stderr, "FastAWrapper::getSequenceOnDisk()--   Squeezed\n");
    fprintf(stderr, "FastAWrapper::getSequenceOnDisk()-- Pester Bri to fix this.\n");
    exit(1);
  }
#endif

  //  XXX: A find has already been performed.  This is probably wasted
  //  (it loads the readbuffer).  We should pass the read buffer into
  //  the SeqOnDisk, and allocate a new one for the wrapper.

  FastASequenceOnDisk  *s = new FastASequenceOnDisk(getSourceName(),
                                                    _theSeqs[_currentSequenceNumber]._headerStart,
                                                    _theSeqs[_currentSequenceNumber]._headerLen,
                                                    _theSeqs[_currentSequenceNumber]._seqStart,
                                                    _theSeqs[_currentSequenceNumber]._seqLen,
                                                    _currentSequenceNumber,
                                                    _theGlobalDesc._squeezedSequences,
                                                    _theGlobalDesc._fixedWidth,
                                                    _theGlobalDesc._seqlineLength,
                                                    _theGlobalDesc._seqendlLength);

  //  To be compatible with the SequenceInCore, we need to advance to
  //  the next sequence.
  //
  _currentSequenceNumber++;

  return(s);
}


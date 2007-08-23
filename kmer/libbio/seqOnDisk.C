
seqOnDisk::seqOnDisk(char const *filename,
                     u64bit hdrstart, u32bit hdrlen,
                     u64bit seqstart, u32bit seqlen,
                     seqIID iid,
                     bool isSqueezed,
                     bool isFixedWidth,
                     u32bit sl,
                     u32bit ss) {

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

  _sequence = 0L;

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

seqOnDisk::seqOnDisk(seqIID iid,
                     char *hdr, u32bit hdrlen,
                     char *seq, u32bit seqlen) {
  _readBuffer       = 0L;;
  _idx              = iid;
  _headerLength     = hdrlen;
  _sequenceLength   = seqlen;
  _headerStart      = ~u64bitZERO;
  _sequenceStart    = ~u64bitZERO;
  _header           = hdr;
  _sequence         = seq;
  _sequencePosition = 0;
  _sourceType       = 3;
  _lineLength       = ~u32bitZERO;
  _lineSep          = ~u32bitZERO;
}

seqOnDisk::~seqOnDisk() {
  delete    _readBuffer;
  delete [] _header;
  delete [] _sequence;
}

char*
seqOnDisk::getChars(char *block, u32bit position, u32bit length) {

  //  Check that there is enough sequence left!
  //
  if (position >= _sequenceLength) {
    fprintf(stderr, "seqOnDisk::getChars()-- position "u32bitFMT" not in sequence of length "u32bitFMT".\n",
            position, length);
    return(0L);
  }

  if (position + length > _sequenceLength) {
    fprintf(stderr, "seqOnDisk::getChars()-- requested more sequence ("u32bitFMT" letters @ position "u32bitFMT")"
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
        while (whitespaceSymbol[_readBuffer->get()])
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
        while (whitespaceSymbol[_readBuffer->get()])
          _readBuffer->next();    
      }

      for (u32bit i=0; i<length; i++) {
        block[i] = _readBuffer->get();
        _readBuffer->next();    
        while (whitespaceSymbol[_readBuffer->get()])
          _readBuffer->next();    
      }
      break;
    case 3:
      //  Yow!  A faked on-disk sequence!
      strncpy(block, _sequence + position, length);
      break;
    default:
      assert(0);
      break;
  }

  block[length] = 0;

  return(block);
}

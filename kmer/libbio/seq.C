#include "bio++.H"



seqInCore::seqInCore(seqIID iid) {
  _idx       = iid;

  _headerLen = 0;
  _headerMax = 0;
  _header    = 0L;

  _seqLen = 0;
  _seqMax = 0;
  _seq    = 0L;
}

seqInCore::seqInCore(seqIID iid, char *hdr, u32bit hdrlen, char *seq, u32bit seqlen) {
  _idx       = iid;

  _headerLen = hdrlen;
  _headerMax = hdrlen;
  _header    = hdr;

  _seqLen = seqlen;
  _seqMax = seqlen;
  _seq    = seq;
}

seqInCore*
seqInCore::copy(void) {
  char *h = new char [_headerLen + 1];
  char *s = new char [_seqLen    + 1];

  memcpy(h, _header, _headerLen + 1);
  memcpy(s, _seq,    _seqLen    + 1);

  return(new seqInCore(_idx, h, _headerLen, s, _seqLen));
}

seqInCore::~seqInCore() {
  delete [] _header;
  delete [] _seq;
}

////////////////////////////////////////////////////////////////////////////////

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


////////////////////////////////////////////////////////////////////////////////

void
seqStream::clearGuts(void) {

  _positionInSequence = 0;
  _positionInStream   = 0;
  _theSequenceNumber  = 0;

  _filename           = 0L;
  _file               = 0L;
  _fileToDelete       = 0L;
  _sequence           = 0L;

  _useListLen         = 0;
  _useListMax         = 32;
  _useList            = new use_s [_useListMax + 1];

  _currentSeq         = 0;

  _positionInSequence = 0;
  _positionInStream   = 0;
  _lengthOfSequences  = 0;

  _eof                = false;

  _separator          = '.';
  _separatorLength    = 0;  //  Default of 1 .
  _separatorPosition  = 0;
  _separatorDone      = false;
}

seqStream::seqStream(void) {
  clearGuts();
}

seqStream::seqStream(const char *filename, bool finishMe) {
  clearGuts();
  setFile(filename);
  if (finishMe)
    finish();
}

seqStream::seqStream(seqFile *S, bool finishMe) {
  clearGuts();
  setFile(S);
  if (finishMe)
    finish();
}

seqStream::~seqStream() {
  delete [] _filename;
  delete    _fileToDelete;
  delete [] _useList;
  delete    _sequence;
}

void
seqStream::setFile(const char *filename) {
  _fileToDelete = openSeqFile(filename);
  setFile(_fileToDelete);
}

void
seqStream::setFile(seqFile *S) {
  _file = S;
  _file->openIndex();

  if (_file->isSqueezed() == false) {
    fprintf(stderr, "ERROR:  seqStream::seqStream()-- source file '%s' isn't squeezed,\n", _file->getSourceName());
    fprintf(stderr, "ERROR:      and this causes errors in indexing.  Bri should really fix this.\n");
    fprintf(stderr, "ERROR:      until then, squeeze it with:\n");
    fprintf(stderr, "ERROR:  leaff -f %s -W > %s.squeezed\n", _file->getSourceName(), _file->getSourceName());
    exit(1);
  }

  _filename = new char [strlen(_file->getSourceName()) + 1];
  strcpy(_filename, _file->getSourceName());
}








void
seqStream::parse(char *line) {
  u32bit  v=0, u=0;
  char   *rest = line;

  //  line can be either a list of numbers, or a file.  See if "line" opens
  //  as a file.  If it does, read in the numbers.
  //
  FILE *F = fopen(line, "r");
  if (F) {
    while (!feof(F)) {
      if (fscanf(F, " "u32bitFMT" ", &v))
        add(v);
    }
    fclose(F);
  } else {
    while (*line) {

      //  We are at a number.  Get it.
      //
      v = (u32bit)strtoul(line, &rest, 10);
      line = rest;
      
      //  If we get ',' or 0, add the single number to the list.
      //  If we get '-', decode the rest of the range.
      //
      switch (*line) {
        case 0:
        case ',':
          add(v);
          break;
        case '-':
          line++;
          u = (u32bit)strtoul(line, &rest, 10);
          line = rest;

          if (v > u)
            for (; u <= v; u++)
              add(u);
          else
            for (; v <= u; v++)
              add(v);
          break;
        default:
          fprintf(stderr, "Invalid -use specification -- no separator found -- format should be, e.g., \"1,2,4-9,10\".\n");
          fprintf(stderr, "Trouble starts at '%s'\n", line);
          exit(1);
          break;
      }

      //  We should be at a ',' (or end of line).  Anything else is
      //  an error.  If ',', move to the next number.
      //
      switch (*line) {
        case 0:
          break;
        case ',':
          line++;
          break;
        default:
          fprintf(stderr, "Invalid -use specification -- no separator found -- format should be, e.g., \"1,2,4-9,10\".\n");
          fprintf(stderr, "Trouble starts at '%s'\n", line);
          exit(1);
          break;
      }
    }
  }
}


void
seqStream::add(u32bit v) {

  if (_useListLen >= _useListMax) {
    _useListMax <<= 1;
    use_s *u = new use_s [_useListMax + 1];
    memcpy(u, _useList, sizeof(use_s) * _useListLen);
    delete [] _useList;
    _useList = u;
  }

  _useList[_useListLen].iid    = v;
  _useList[_useListLen].length = 0;
  _useList[_useListLen].start  = 0;
  _useListLen++;
}





void
seqStream::finish(void) {
  bool    useAll   = false;
  bool   *seen     = 0L;
  u64bit  seenLen  = 0;
  u64bit  startPos = 0;

  //  If silly user didn't add any sequences, use all the sequences.
  //  Otherwise, fix the existing one (remove duplicates, sort).

  if (_useListLen == 0) {
    useAll  = true;
    seen    = 0L;
    seenLen = _file->getNumberOfSequences();

    delete [] _useList;

    _useListMax = seenLen;
    _useList    = new use_s [_useListMax + 1];
  } else {
    useAll  = false;
    seen    = new bool [_file->getNumberOfSequences()];
    seenLen = 0;

    for (u32bit i=_file->getNumberOfSequences(); i--; )
      seen[i] = false;

    for (u32bit i=0; i<_useListLen; i++) {
      if (_useList[i].iid >= _file->getNumberOfSequences()) {
        fprintf(stderr, "WARNING: seqStream requested sequence iid "u32bitFMT" which isn't in '%s'\n",
                i, _filename);
      } else {
        if (seen[_useList[i].iid] == false) {
          seen[_useList[i].iid] = true;
          seenLen++;
        }
      }
    }
  }

  _useListLen = 0;

  for (u32bit i=0; i<_file->getNumberOfSequences(); i++) {
    if ((useAll || seen[i]) &&
        (_file->sequenceLength(i) > 0)) {
      _useList[_useListLen].iid    = i;
      _useList[_useListLen].length = _file->sequenceLength(i);
      _useList[_useListLen].start  = startPos;

#if 0
      fprintf(stderr, "seqStream::finish()-- added iid="u32bitFMTW(6)" len="u32bitFMTW(9)" start="u32bitFMTW(9)"\n",
              _useList[_useListLen].iid,
              _useList[_useListLen].length,
              _useList[_useListLen].start);
#endif

      startPos += _useList[_useListLen].length + _separatorLength + 1;

      _useListLen++;
    }
  }

  delete [] seen;

  _lengthOfSequences = startPos;

  _useList[_useListLen].iid      = ~u32bitZERO;
  _useList[_useListLen].length   = 0;
  _useList[_useListLen].start    = _lengthOfSequences;

  //  Initialize -- grab the first sequence
  //
  _file->find(_useList[0].iid);
  _sequence = _file->getSequenceOnDisk();
}










bool
seqStream::rewind(void) {
  _currentSeq         = 0;
  _positionInSequence = 0;
  _positionInStream   = 0;
  _eof                = false;
  _separatorDone      = false;
  _separatorPosition  = 0;

  delete _sequence;

  _file->find(_useList[0].iid);
  _sequence = _file->getSequenceOnDisk();

  return(true);
}


unsigned char
seqStream::get(void) {
  char ret = _sequence->get();

  _positionInStream++;

  //  The while is needed to skip over empty sequences -- without it,
  //  we'll catch the end of the last sequence (ret == 0, and we'll
  //  enter this block), do the separator, grab the next sequence, get
  //  the first letter into ret, and blindly return that first letter.
  //  If the sequence is emtpy, the first letter is 0, and we then
  //  terminate the seqStream early.
  //
  while (ret == 0) {

    //  Are we at the end of the chained sequence?
    //
    if (_currentSeq + 1 >= _useListLen) {
      _eof = true;
      return(0);
    }

    //  Nope, are we doing the separator?
    //
    if (_separatorPosition > 0) {
      _separatorPosition--;
      return(_separator);
    }

    //  Nope, did we just finish the separator, or do we need to start
    //  it?  The flow being: the first time _sequence->get() returns
    //  0, _separatorDone is false, and we start the separator.  We do
    //  the separator.  It finishes, and we're back here, with
    //  _separatorDone set to true.
    //
    //  On _separatorDone, we initialize the next sequence, load the
    //  first letter and fall through to the normal "got a good read"
    //  exit.
    //
    if (_separatorDone) {
      delete _sequence;
      _sequence = 0L;
      _currentSeq++;

      if (_file->find(_useList[_currentSeq].iid)) {
        _sequence           = _file->getSequenceOnDisk();
        _separatorDone      = false;
        _positionInSequence = u64bitZERO;
        ret = _sequence->get();
      }
    } else {
      _separatorPosition  = _separatorLength;
      _separatorDone      = true;
      _positionInSequence = ~u64bitZERO;
      return(_separator);
    }
  }

  //  Got a good read, bump to the next character and return
  //
  _sequence->next();
  _positionInSequence++;
  return(ret);
}


unsigned char
seqStream::nextSymbol(void) {

  unsigned char  nc = get();
  while (whitespaceSymbol[nc] && (eof() == false))
    nc = get();

  //  EOF?  Return eof.  Valid?  Return the letter.  Not valid?
  //  Return the "new sequence" code, or the "gap" code.
  //
  if (eof())
    return(0);
  if (validSymbol[nc])
    return(nc);
  if (nc == getSeparator())
    return(254);
  return(253);
}







u32bit
seqStream::sequenceNumberOfPosition(u64bit p) {
  u32bit   sret = ~u32bitZERO;

  //  binary search on our list of start positions, to find the
  //  sequence that p is in.

  if (_lengthOfSequences < p) {
    fprintf(stderr, "seqStream::sequenceNumberOfPosition()-- WARNING! Position p="u64bitFMT" too big; only "u64bitFMT" positions.\n",
            p, _lengthOfSequences);
    return(~u32bitZERO);
  }

  if (_useListLen < 16) {
    for (u32bit ss=0; ss<_useListLen; ss++) {
      if ((_useList[ss].start <= p) && (p < _useList[ss+1].start)) {
        sret = ss;
        break;
      }
    }
  } else {
    u32bit  lo = 0;
    u32bit  hi = _useListLen;
    u32bit  md = 0;

    while (lo <= hi) {
      md = (lo + hi) / 2;

      if        (p < _useList[md].start) {
        //  This block starts after the one we're looking for.  
        hi = md;

      } else if ((_useList[md].start <= p) && (p < _useList[md+1].start)) {
        //  Got it!
        lo           = md + 1;
        hi           = md;
        sret         = md;

      } else {
        //  By default, then, the block is too low.
        lo = md;
      }
    }
  }

  //fprintf(stderr, "seqStream::sequenceNumberOfPosition()-- p="u64bitFMT" --> s="u32bitFMT"\n", p, sret);

  return(sret);
}





////////////////////////////////////////////////////////////////////////////////


seqFactory *seqFactory::me = 0L;


seqFactory::seqFactory() {
  _filesNum = 0;
  _filesMax = 16;
  _files = new seqFile * [_filesMax];

  registerFile(new fastaFile);
}


seqFactory::~seqFactory() {
  for (u32bit i=0; i<_filesNum; i++)
    delete _files[i];
  delete [] _files;
}


void           
seqFactory::registerFile(seqFile *f) {
  if (_filesNum >= _filesMax) {
    fprintf(stderr, "Hmmmm!  Wow!  You registered lots of files!  Now fix %s at line %d.\n", __FILE__, __LINE__);
    exit(1);
  }
  _files[_filesNum++] = f;
}


seqFile *
seqFactory::openFile(const char *name) {
  seqFile  *n = 0L;

  for (u32bit i=0; i<_filesNum; i++) {
    n = _files[i]->openFile(name);
    if (n)
      return(n);
  }

  fprintf(stderr, "ERROR: Cannot determine type of file '%s'.\n", name);
  exit(1);
  return(n);
}

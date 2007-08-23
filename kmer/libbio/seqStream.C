
void
seqStream::clearGuts(void) {

  _filename           = 0L;
  _file               = 0L;
  _fileToDelete       = 0L;

  _currentSeq         = 0;
  _currentPos         = 0;
  _sequence           = 0L;

  _useListLen         = 0;
  _useListMax         = 32;
  _useList            = new use_s [_useListMax + 1];

  _lengthOfSequences  = 0;

  _eof                = false;

  _separator          = '.';
  _separatorLength    = 0;  //  Default of 1 separator letter
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



void
seqStream::setSeparator(char sep, u32bit len) {
  
  if (len == 0)
    fprintf(stderr, "seqStream::setSeparator()-- ERROR: the separator length must be at least one.\n"), exit(1);

  _separator = sep;
  _separatorLength = len - 1;    

  //  Rebuild the useList -- to update the start position of each
  //  sequence given the new separatorLength -- if we've already
  //  called finish().
  //
  if (_lengthOfSequences > 0)
    finish();
};


bool
seqStream::rewind(void) {
  _currentSeq         = 0;
  _currentPos         = 0;
  _eof                = false;
  _separatorDone      = false;
  _separatorPosition  = 0;

  delete _sequence;

  _file->find(_useList[_currentSeq].iid);
  _sequence = _file->getSequenceOnDisk();

  return(true);
}


unsigned char
seqStream::get(void) {
  char ret = _sequence->get();

  //  We always return a letter, and so always update the current
  //  position in the stream.
  //
  _currentPos++;

  //  The while is needed to skip over empty sequences -- without it,
  //  we'll catch the end of the last sequence (ret == 0, and we'll
  //  enter this block), do the separator, grab the next sequence, get
  //  the first letter into ret, and blindly return that first letter.
  //  If the sequence is emtpy, the first letter is 0, and we then
  //  terminate the seqStream early.
  //
  while (ret == 0) {

    //  Geez, very obnoxiously, we can tell if the _next_ get()
    //  is going to return EOF, but we cannot tell if the _last_
    //  get() was an eof.  So we set a flag.
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
        _currentPos         = 1;
        ret = _sequence->get();
      }
    } else {
      _separatorPosition  = _separatorLength;
      _separatorDone      = true;
      return(_separator);
    }
  }

  //  Got a good read, bump to the next character and return
  //
  _sequence->next();
  return(ret);
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

  //  But if we're in a range of separators, return invalid.
  if ((sret < _useListLen) && (_useList[sret+1].start - _separatorLength - 1 <= p))
    sret = ~u32bitZERO;

  return(sret);
}

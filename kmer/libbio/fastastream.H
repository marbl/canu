#ifndef BIO_PLUS_PLUS_H
#error Please use bio++.H instead of fastastream.H
#endif

#ifndef FASTASTREAM_H
#define FASTASTREAM_H

//
//  Reads sequence from a multiFastA file, letter by letter, saving
//  the defline.  Returns values are:
//
//  0    EOF
//  254  New sequence starting
//  253  Break in contiguous sequence
//
//  If the letter is a valid symbol the letter is returned.
//  Otherwise, one of code 253 or 254 is returned.
//

class FastAstream {
private:
  void  initialClear(void) {
    _positionInSequence  = 0;
    _positionInStream    = 0;
    _sequenceNumber      = 0;

    _rb                  = 0L;
    _cs                  = 0L;

    _defLineMax          = 2048;
    _defLinePos          = 0;
    _defLine             = new char [_defLineMax + 1];
  }
public:
  FastAstream(chainedSequence *S) {
    initialClear();

    _cs = S;
  };

  FastAstream(const char *filename) {
    initialClear();

    if ((filename == 0L) ||
        ((filename[0] == '-') && (filename[1] == 0)))
      _rb = new readBuffer(fileno(stdin), "stdin");
    else
      _rb = new readBuffer(filename);
  };

  ~FastAstream() {
    delete    _rb;
    delete [] _defLine;
  };
  
  bool              rewind(void) {
    _positionInSequence  = 0;
    _positionInStream    = 0;
    _sequenceNumber      = 0;

    if (_rb) {
      return(_rb->seek(0));
    } else {
      _cs->rewind();
      return(true);
    }
  };

  unsigned char     nextSymbol(void);

  //  The position in the sequence and the source stream.  -1 because
  //  we pre-increment these values.
  //
  u64bit            thePositionInSequence(void) { return(_positionInSequence - 1); };
  u64bit            thePositionInStream(void)   { return(_positionInStream - 1); };
  u64bit            theSequenceNumber(void)     { return(_sequenceNumber); };
  char             *theDefLine(void)            { return(_defLine); };
private:
  char               getNextCharacter(void) {

    //  The stream position is off by one because we pre-increment.  It
    //  is adjusted on access.  Stream position is incremented each time
    //  we getNextCharacter().
    //  
    _positionInStream++;

    if (_rb)
      return(_rb->getnext());
    else
      return(_cs->get());
  };

  bool               isEOF(void) {
    if (_rb)
      return(_rb->eof());
    return(_cs->eof());
  };

  u64bit            _positionInSequence;
  u64bit            _positionInStream;
  u64bit            _sequenceNumber;

  readBuffer       *_rb;
  chainedSequence  *_cs;

  u32bit            _defLineMax;
  u32bit            _defLinePos;
  char             *_defLine;
#ifdef __APPLE__
  u32bit            _pad;
#endif
};


inline
unsigned char
FastAstream::nextSymbol(void) {
  unsigned char  nc;

  //  Grab the next character
  //
  nc = getNextCharacter();

  //  Skip any whitespace
  //
  while (whitespaceSymbol[nc] && !isEOF())
    nc = getNextCharacter();

  //  Are we at the EOF?
  //
  if (isEOF())
    return(0);

  //  Looks like we are in the sequence.  Increment the position.  If
  //  this is the start of a defline, the position is reset to zero.
  //
  _positionInSequence++;

  //  If it's valid, just return it.
  //
  if (validSymbol[nc])
    return(nc);

  //  It's not a valid character.  If it's a defline, read until we
  //  hit a '\n' or '\r'.
  //
  if (nc == '>') {

    //  Reset thePosition and increment theSequenceNumber
    //
    _positionInSequence = 0;
    _sequenceNumber++;

    _defLinePos = 0;

    nc = getNextCharacter();
    while ((nc != '\n') && (nc != '\r') && !isEOF()) {

      //  More space?  Must be a monster defline!
      //
      if (_defLinePos >= _defLineMax) {
        _defLineMax <<= 1;
        char *d = new char [_defLineMax + 1];

        for (u32bit i=_defLinePos; --i; )
          d[i] = _defLine[i];

        delete [] _defLine;
        _defLine = d;
      }

      //  Push the letter onto the end of the defline
      //
      _defLine[_defLinePos++] = nc;

      nc = getNextCharacter();
    }

    _defLine[_defLinePos] = 0;

    //  Return 254; we've just read in the defline.
    //
    return(254);
  }

  //  We've hit a break, return 253.
  //
  return(253);
}


#endif  //  FASTASTREAM_H

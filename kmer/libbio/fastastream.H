#ifndef FASTASTREAM_H
#define FASTASTREAM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include "libbritypes.h"
#include "alphabet.h"
#include "readBuffer.H"

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
public:
  FastAstream(const char *filename) {
    if ((filename == 0L) ||
        ((filename[0] == '-') && (filename[1] == 0)))
      _buffer = new readBuffer(fileno(stdin), "stdin");
    else
      _buffer = new readBuffer(filename);

    _thePosition       = 0;
    _theSequenceNumber = 0;

    _theDefLineMax     = 2048;
    _theDefLinePos     = 0;
    _theDefLine        = new char [_theDefLineMax + 1];
  };

  ~FastAstream() {
    delete [] _buffer;
    delete [] _theDefLine;
  };

  unsigned char     nextSymbol(void);

  u64bit            thePosition(void) {
    return(_thePosition);
  };
  u64bit            theSequenceNumber(void) {
    return(_theSequenceNumber);
  };
  char             *theDefLine(void) {
    return(_theDefLine);
  };
private:
  readBuffer  *_buffer;

  u64bit       _thePosition;
  u64bit       _theSequenceNumber;

  u32bit       _theDefLineMax;
  u32bit       _theDefLinePos;
  char        *_theDefLine;
};


inline
unsigned char
FastAstream::nextSymbol(void) {
  unsigned char  nc;

  //  Grab the next character
  //
  nc = (unsigned char)_buffer->getnext();

  //  Skip any whitespace
  //
  while (isspace(nc) && !_buffer->eof())
    nc = (unsigned char)_buffer->getnext();

  //  Are we at the EOF?
  //
  if (_buffer->eof())
    return(0);

  //  Looks like we are in the sequence.  Increment the position.  If
  //  this is the start of a defline, the position is reset to zero.
  //
  _thePosition++;

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
    _thePosition = 0;
    _theSequenceNumber++;

    _theDefLinePos = 0;

    nc = (unsigned char)_buffer->getnext();
    while ((nc != '\n') && (nc != '\r') && !_buffer->eof()) {

      //  More space?  Must be a monster defline!
      //
      if (_theDefLinePos >= _theDefLineMax) {
        _theDefLineMax <<= 1;
        char *d = new char [_theDefLineMax + 1];

        for (u32bit i=_theDefLinePos; --i; )
          d[i] = _theDefLine[i];

        delete [] _theDefLine;
        _theDefLine = d;
      }

      //  Push the letter onto the end of the defline
      //
      _theDefLine[_theDefLinePos++] = nc;

      nc = (unsigned char)_buffer->getnext();
    }

    _theDefLine[_theDefLinePos] = 0;

    //  Return 254; we've just read in the defline.
    //
    return(254);
  }

  //  We've hit a break, return 253.
  //
  return(253);
}

#endif  //  FASTASTREAM_H

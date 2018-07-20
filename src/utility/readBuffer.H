
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef READ_BUFFER_H
#define READ_BUFFER_H

#include "AS_global.H"
#include "memoryMappedFile.H"

class readBuffer {
public:
  readBuffer(const char *filename, uint64 bufferMax = 32 * 1024);
  readBuffer(FILE *F, uint64 bufferMax = 32 * 1024);
  ~readBuffer();

  bool                 eof(void) { return(_eof); };

  char                 peek(void);

  char                 read(void);
  char                 readuntil(char stop);
  uint64               read(void *buf, uint64 len);
  uint64               read(void *buf, uint64 maxlen, char stop);

  void                 skipLine(void);
  void                 skipAhead(char stop);
  uint64               copyUntil(char stop, char *dest, uint64 destLen);

  void                 seek(uint64 pos);
  uint64               tell(void) { return(_filePos); };

  const char          *filename(void) { return(_filename); };

private:
  void                 fillBuffer(void);
  void                 init(int fileptr, const char *filename, uint64 bufferMax);

  char               *_filename;

  int                 _file;
  uint64              _filePos;

  memoryMappedFile   *_mmap;
  bool                _stdin;

  bool                _eof;

  //  If bufferMax is zero, then we are using the mmapped interface, otherwise,
  //  we are using a open()/read() and a small buffer.

  uint64              _bufferPos;
  uint64              _bufferLen;
  uint64              _bufferMax;
  char               *_buffer;
};


//  Returns the next letter in the buffer, but DOES NOT advance past
//  it.  Might have some wierd interaction with EOF -- if you peek()
//  and the next thing is eof , the _eof flag might get set.
//
inline
char
readBuffer::peek(void) {

  if ((_eof == false) && (_bufferPos >= _bufferLen))
    fillBuffer();

  if (_eof)
    return(0);

  return(_buffer[_bufferPos]);
}


//  Returns the next letter in the buffer and advances file position by one.
//  Returns 0 on EOF.
//
inline
char
readBuffer::read(void) {

  if ((_eof == false) && (_bufferPos >= _bufferLen))
    fillBuffer();

  if (_eof)
    return(0);

  _bufferPos++;
  _filePos++;

  return(_buffer[_bufferPos-1]);
}



//  Like read(), but does not advance file position if
//  the 'stop' character is encountered.
//  Returns 'stop' on EOF.
//
inline
char
readBuffer::readuntil(char stop) {
  char  ch;

  if ((_eof == false) && (_bufferPos >= _bufferLen))
    fillBuffer();

  if (_eof)
    return(0);

  ch = _buffer[_bufferPos];

  if (ch == stop)
    return(stop);

  _bufferPos++;
  _filePos++;

  return(ch);
}



inline
void
readBuffer::skipLine(void) {
  for (char ch = read(); (ch != '\n') && (ch != 0); ch = read())
    ;
}



//  Advances the file position to the next 'stop' character and returns
//  'stop', or 0 if eof.
inline
void
readBuffer::skipAhead(char stop) {

#if 0
  while (_eof == false) {
    while ((_bufferPos < _bufferLen) &&
           (_buffer[_bufferPos] != stop)) {
      _bufferPos++;
      _filePos++;
    }

    if ((_bufferPos < _bufferLen) &&
        (_buffer[_bufferPos] == stop))
      return;

    //assert(_bufferPos == _bufferLen);

    fillBuffer();

    //if (_eof)
    //  return;
  }
#else
  //_buffer[_bufferLen] = stop;
  assert(_buffer[_bufferLen] == stop);

  while (_eof == false) {
    //fprintf(stderr, "search for stop starting at %u to len=%u\n", _bufferPos, _bufferLen);
    while (_buffer[_bufferPos] != stop) {
      _bufferPos++;
      _filePos++;
    }

    //assert(_buffer[_bufferPos] == stop);

    if (_bufferPos < _bufferLen) {
      //fprintf(stderr, "found at _pos=%u\n", _bufferPos);
      return;
    }

    //fprintf(stderr, "not found at _pos=%u\n", _bufferPos);
    //assert(_bufferPos == _bufferLen);

    //_bufferPos--;
    //_filePos--;

    fillBuffer();

    //if (_eof)
    //  return;
  }
#endif
}


//  Advances the file position to the next 'stop' character, copying
//  the characters into a user buffer.
//
//  Returns the number of letters copied.
//
inline
uint64
readBuffer::copyUntil(char stop, char *dest, uint64 destLen) {
  uint64  copied = 0;

  while (_eof == false) {
    while ((_buffer[_bufferPos] != stop) &&
           (copied < destLen)) {
      dest[copied] = _buffer[_bufferPos];

      copied++;
      _bufferPos++;
      _filePos++;
    }

    if (_bufferPos < _bufferLen)
      return(copied);

    fillBuffer();
  }

  return(copied);
}




#endif  //  READ_BUFFER_H

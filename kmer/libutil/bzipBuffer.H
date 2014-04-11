#ifndef BZIP_BUFFER_H
#define BZIP_BUFFER_H

#include <stdio.h>
#include <stdlib.h>
#include <bzlib.h>

#include "util.h"

class bzipBuffer {
public:
  bzipBuffer(const char *filename, uint32 bufferMax = 32 * 1024);
  ~bzipBuffer();

  bool      eof(void);
  bool      next(void);
  char      get(void);
  char      getnext(void);
  bool      seek(off_t pos);
  size_t    read(char *buf, size_t len);   //  read the next len bytes into the user buffer buf
  off_t     tell(void);

private:
  void      fillBuffer(void);
  void      init(int fileptr, const char *filename, uint32 bufferMax);

  char       *_filename;

  int         _file;
  off_t       _filePos;
  bool        _eof;

  uint32      _bzip2bufferMax;
  uint32      _bzip2inPos;
  uint32      _bzip2outPos;

  char       *_bzip2in;
  char       *_bzip2out;

  bool        _bzip2streamEnd;
  bz_stream   _bzip2stream;
};


inline
bool
bzipBuffer::eof(void) {
  return(_eof);
}


inline
bool
bzipBuffer::next(void) {

  if (_eof)
    return(true);

  _bzip2outPos++;
  _filePos++;

  if (_bzip2outPos >= _bzip2stream.avail_out)
    fillBuffer();

  return(_eof);
}


inline
char
bzipBuffer::get(void) {
  return(_bzip2out[_bzip2outPos]);
}


inline
char
bzipBuffer::getnext(void) {
  char  x = _bzip2out[_bzip2outPos];
  next();
  return(x);
}


inline
off_t
bzipBuffer::tell(void) {
  return(_filePos);
}


#endif  //  BZIP_BUFFER_H

#include "util++.H"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>


//  If bufferMax is zero, then the file is accessed using memory
//  mapped I/O.  Otherwise, a small buffer is used.
//
readBuffer::readBuffer(const char *filename, u64bit bufferMax) {

  _filename    = 0L;
  _file        = 0;
  _filePos     = 0;
  _mmap        = false;
  _stdin       = false;
  _eof         = false;
  _bufferPos   = 0;
  _bufferLen   = 0;
  _bufferMax   = 0;
  _buffer      = 0L;

  if (((filename == 0L) && (isatty(fileno(stdin)) == 0)) ||
      ((filename != 0L) && (filename[0] == '-') && (filename[1] == 0))) {
    _filename  = new char [32];
    strcpy(_filename, "(stdin)");

    _stdin = true;

    if (bufferMax == 0)
      bufferMax = 32 * 1024;
  } else if (filename == 0L) {
    fprintf(stderr, "readBuffer()-- no filename supplied, and I will not use the terminal for input.\n"), exit(1);
  } else {
    _filename  = new char [strlen(filename) + 1];
    strcpy(_filename, filename);
  }

  if (bufferMax == 0) {
    _mmap   = true;
    _buffer = (char *)mapFile(_filename, &_bufferLen, 'r');
  } else {
    errno = 0;
    _file = (_stdin) ? fileno(stdin) : open(_filename, O_RDONLY | O_LARGEFILE);
    if (errno)
      fprintf(stderr, "readBuffer()-- couldn't open the file '%s': %s\n",
              _filename, strerror(errno)), exit(1);

    _bufferMax   = bufferMax;
    _buffer      = new char [_bufferMax];
  }

  fillBuffer();

  if (_bufferLen == 0)
    _eof   = true;
}


readBuffer::readBuffer(FILE *file, u64bit bufferMax) {

  if (bufferMax == 0)
    fprintf(stderr, "readBuffer()-- WARNING: mmap() not supported in readBuffer(FILE *)\n");

  _filename    = new char [32];
  _file        = fileno(file);
  _filePos     = 0;
  _mmap        = false;
  _stdin       = false;
  _eof         = false;
  _bufferPos   = 0;
  _bufferLen   = 0;
  _bufferMax   = (bufferMax == 0) ? 32 * 1024 : bufferMax;
  _buffer      = new char [_bufferMax];

  strcpy(_filename, "(hidden file)");

  //  Just be sure that we are at the start of the file.
  errno = 0;
  lseek(_file, 0, SEEK_SET);
  if (errno)
    fprintf(stderr, "readBuffer()-- '%s' couldn't seek to position 0: %s\n",
            _filename, strerror(errno)), exit(1);

  fillBuffer();

  if (_bufferLen == 0)
    _eof   = true;
}


readBuffer::~readBuffer() {

  delete [] _filename;

  if (_mmap)
    unmapFile(_buffer, _bufferLen);
  else
    delete [] _buffer;

  if (_stdin == false)
    close(_file);
}


void
readBuffer::fillBuffer(void) {

  //  If there is still stuff in the buffer, no need to fill.
  if (_bufferPos < _bufferLen)
    return;

  //  No more stuff in the buffer.  But if mmap'd, ths means we're EOF.
  if (_mmap) {
    _eof = true;
    return;
  }

  _bufferPos = 0;
  _bufferLen = 0;

 again:
  errno = 0;
  _bufferLen = (u64bit)::read(_file, _buffer, _bufferMax);
  if (errno == EAGAIN)
    goto again;
  if (errno)
    fprintf(stderr, "readBuffer::fillBuffer()-- only read "u64bitFMT" bytes, couldn't read "u64bitFMT" bytes from '%s': %s\n",
            _bufferLen, _bufferMax, _filename, strerror(errno)), exit(1);

  if (_bufferLen == 0)
    _eof = true;
}


void
readBuffer::seek(u64bit pos) {

  if (_stdin == true) {
    if (_filePos < _bufferLen) {
      _filePos   = 0;
      _bufferPos = 0;
      return;
    } else {
      fprintf(stderr, "readBuffer()-- seek() not available for file 'stdin'.\n");
      exit(1);
    }

    return;
  }

  assert(_stdin == false);

  if (_mmap) {
    _bufferPos = pos;
    _filePos   = pos;
  } else {
    errno = 0;
    lseek(_file, pos, SEEK_SET);
    if (errno)
      fprintf(stderr, "readBuffer()-- '%s' couldn't seek to position "s64bitFMT": %s\n",
              _filename, pos, strerror(errno)), exit(1);

    _bufferLen = 0;
    _bufferPos = 0;
    _filePos   = pos;

    fillBuffer();
  }

  _eof       = (_bufferPos >= _bufferLen);
}


u64bit
readBuffer::read(void *buf, u64bit len) {
  char  *bufchar = (char *)buf;

  //  Handle the mmap'd file first.

  if (_mmap) {
    u64bit c = 0;

    while ((_bufferPos < _bufferLen) && (c < len)) {
      bufchar[c++] = _buffer[_bufferPos++];
      _filePos++;
    }

    if (c == 0)
      _eof = true;

    return(c);
  }

  //  Easy case; the next len bytes are already in the buffer; just
  //  copy and move the position.

  if (_bufferLen - _bufferPos > len) {
    memcpy(bufchar, _buffer + _bufferPos, len);
    _bufferPos += len;

    fillBuffer();

    _filePos   += len;

    return(len);
  }

  //  Existing buffer not big enough.  Copy what's there, then finish
  //  with a read.

  u64bit   bCopied = 0;   //  Number of bytes copied into the buffer
  u64bit   bRead   = 0;   //  Number of bytes read into the buffer
  u64bit   bAct    = 0;   //  Number of bytes actually read from disk

  memcpy(bufchar, _buffer + _bufferPos, _bufferLen - _bufferPos);
  bCopied    = _bufferLen - _bufferPos;
  _bufferPos = _bufferLen;

  while (bCopied + bRead < len) {
    errno = 0;
    bAct = (u64bit)::read(_file, bufchar + bCopied + bRead, len - bCopied - bRead);
    if (errno)
      fprintf(stderr, "readBuffer()-- couldn't read "u64bitFMT" bytes from '%s': n%s\n",
              len, _filename, strerror(errno)), exit(1);

    //  If we hit EOF, return a short read
    if (bAct == 0)
      len = 0;

    bRead += bAct;
  }

  fillBuffer();

  _filePos += bCopied + bRead;

  return(bCopied + bRead);
}


u64bit
readBuffer::read(void *buf, u64bit maxlen, char stop) {
  char  *bufchar = (char *)buf;
  u64bit c = 0;

  //  We will copy up to 'maxlen'-1 bytes into 'buf', or stop at the first occurrence of 'stop'.
  //  This will reserve space at the end of any string for a zero-terminating byte.
  maxlen--;

  if (_mmap) {
    //  Handle the mmap'd file first.
    while ((_bufferPos < _bufferLen) &&
           (c < maxlen)) {
      bufchar[c++] = _buffer[_bufferPos++];

      if (bufchar[c-1] == stop)
        break;
    }

    if (_bufferPos >= _bufferLen)
      _eof = true;

  } else {
    //  And the usual case.
    while ((_eof == false) && (c < maxlen)) {
      bufchar[c++] = _buffer[_bufferPos++];

      if (_bufferPos >= _bufferLen)
        fillBuffer();

      if (bufchar[c-1] == stop)
        break;
    }
  }

  bufchar[c] = 0;

  return(c);
}

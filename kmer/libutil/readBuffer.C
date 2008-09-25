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
readBuffer::readBuffer(const char *filename, u32bit bufferMax) {

  _filename    = 0L;
  _file        = 0;
  _filePos     = 0;
  _mmap        = false;
  _stdin       = false;
  _eof         = false;
  _valid       = false;
  _bufferPos   = 0;
  _bufferLen   = 0;
  _bufferMax   = 0;
  _buffer      = 0L;

  if (filename == 0L)
    filename = "-";

  if (strcmp(filename, "-") == 0) {
    _stdin = true;

    if (bufferMax == 0)
      bufferMax = 32 * 1024;
  }

  _filename  = new char [strlen(filename) + 1];
  strcpy(_filename, filename);

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

  _valid = true;

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
  _bufferLen = (u32bit)::read(_file, _buffer, _bufferMax * sizeof(char));
  if (errno == EAGAIN)
    goto again;
  if (errno)
    fprintf(stderr, "readBuffer::fillBuffer()-- only read %d bytes, couldn't read %d bytes from '%s': %s\n",
            (int)_bufferLen, (int)(_bufferMax * sizeof(char)), _filename, strerror(errno)), exit(1);

  if (_bufferLen == 0)
    _eof = true;
}


void
readBuffer::seek(off_t pos) {

  if (_stdin == true)
    fprintf(stderr, "readBuffer()-- seek() not available for file 'stdin'.\n");

  assert(_valid == true);
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


size_t
readBuffer::read(void *buf, size_t len) {
  char  *bufchar = (char *)buf;

  assert(_valid);

  //  Handle the mmap'd file first.

  if (_mmap) {
    size_t c = 0;

    while ((_bufferPos < _bufferLen) && (c < len))
      bufchar[c++] = _buffer[_bufferPos++];

    if (c == 0)
      _eof = true;

    return(c);
  }

  //  Easy case; the next len bytes are already in the buffer; just
  //  copy and move the position.

  if (_bufferLen - _bufferPos > len) {
    memcpy(bufchar, _buffer + _bufferPos, sizeof(char) * len);
    _bufferPos += (u32bit)len;

    fillBuffer();

    return(len);
  }

  //  Existing buffer not big enough.  Copy what's there, then finish
  //  with a read.

  size_t   bCopied = 0;   //  Number of bytes copied into the buffer
  size_t   bRead   = 0;   //  Number of bytes read into the buffer
  size_t   bAct    = 0;   //  Number of bytes actually read from disk

  memcpy(bufchar, _buffer + _bufferPos, (_bufferLen - _bufferPos) * sizeof(char));
  bCopied    = _bufferLen - _bufferPos;
  _bufferPos = _bufferLen;

  while (bCopied + bRead < len) {
    errno = 0;
    bAct = (u32bit)::read(_file, bufchar + bCopied + bRead, (len - bCopied - bRead) * sizeof(char));
    if (errno)
      fprintf(stderr, "readBuffer()-- couldn't read "u32bitFMT" bytes from '%s': n%s\n",
              (u32bit)(len * sizeof(char)), _filename, strerror(errno)), exit(1);

    //  If we hit EOF, return a short read
    if (bAct == 0)
      len = 0;

    bRead += bAct;
  }

  fillBuffer();

  return(bCopied + bRead);
}

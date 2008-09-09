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
  init(-1, filename, bufferMax);
}


  //  This constructor always uses a small buffer; memory mapped I/O
  //  is not allowed.  This lets us pass in fileno(stdin).
  //
readBuffer::readBuffer(int fileptr, const char *filename, u32bit bufferMax) {
  if (bufferMax == 0) {
    fprintf(stderr, "readBuffer::readBuffer()-- file '%s' requested memory mapped I/O on previously\n", filename);
    fprintf(stderr, "readBuffer::readBuffer()-- opened file, but that's not allowed.  bufferMax reset to 32KB.\n");
    bufferMax = 32 * 1024;
  }
  init(fileptr, filename, bufferMax);
}


readBuffer::~readBuffer() {

  switch (_fileType) {
    case 0:
      delete [] _filename;
      delete [] _buffer;
      break;
    case 1:
      delete [] _filename;
      delete [] _buffer;
      errno = 0;
      close(_file);
      if (errno) {
        fprintf(stderr, "readBuffer()-- WARNING: couldn't close the file '%s': %s\n",
                _filename, strerror(errno));
      }
      break;
    case 2:
      delete [] _filename;
      unmapFile(_buffer, _bufferLen);
      break;
  }
}


void
readBuffer::init(int fileptr, const char *filename, u32bit bufferMax) {

  _filename  = new char [strlen(filename) + 1];
  strcpy(_filename, filename);

  if (bufferMax == 0) {
    _file        = 0;
    _fileType    = 2;
    _filePos     = 0;
    _eof         = false;
    _bufferPos   = 0;
    _bufferLen   = 0;
    _bufferMax   = 0;
    _buffer      = (char *)mapFile(_filename, &_bufferLen, 'r');
  } else {
    if (fileptr == -1) {
      _fileType = 1;
      errno = 0;
      fileptr = open(filename, O_RDONLY | O_LARGEFILE);
      if (errno)
        fprintf(stderr, "readBuffer()-- couldn't open the file '%s': %s\n",
                filename, strerror(errno)), exit(1);
    } else {
      _fileType    = 0;
    }
    _file        = fileptr;
    _filePos     = 0;
    _eof         = false;
    _bufferPos   = 0;          //  Position we are at in the buffer
    _bufferLen   = 0;          //  Valid length of the buffer
    _bufferMax   = bufferMax;  //  Maximum size of the buffer
    _buffer      = new char [_bufferMax];

    _buffer[0] = 0;

    fillBuffer();
  }
}


void
readBuffer::fillBuffer(void) {

  if ((_fileType == 2) || (_bufferPos < _bufferLen))
    return;

  _bufferPos = 0;

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

  if (_fileType == 0)
    fprintf(stderr, "readBuffer()-- seek() not available for file '%s'.\n", _filename), exit(1);

  if (_fileType == 1) {
    errno = 0;
    lseek(_file, pos, SEEK_SET);
    if (errno)
      fprintf(stderr, "readBuffer()-- '%s' couldn't seek to position "s64bitFMT": %s\n",
              _filename, pos, strerror(errno)), exit(1);

    _bufferLen = 0;
    _bufferPos = 0;
    _filePos = pos;
    _eof     = false;

    fillBuffer();
  }

  if (_fileType == 2) {
    _bufferPos = pos;
    _filePos   = pos;
    _eof       = (_bufferPos >= _bufferLen);
  }
}


size_t
readBuffer::read(void *buf, size_t len) {
  char  *bufchar = (char *)buf;

  //  Handle the mmap'd file first.

  if (_fileType == 2) {
    size_t c = 0;

    while ((_bufferPos < _bufferLen) && (c < len))
      bufchar[c++] = _buffer[_bufferPos++];

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

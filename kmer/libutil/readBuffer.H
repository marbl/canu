#ifndef READ_BUFFER_H
#define READ_BUFFER_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <errno.h>

#include "libbritypes.h"
#include "file.h"

class readBuffer {
public:
  //  If bufferMax is zero, then the file is accessed using memory
  //  mapped I/O.  Otherwise, a small buffer is used.
  //
  readBuffer(const char *filename, u32bit bufferMax = 32 * 1024);

  //  This constructor always uses a small buffer; memory mapped I/O
  //  is not allowed.  This lets us pass in fileno(stdin).
  //
  readBuffer(int fileptr, const char *filename, u32bit bufferMax = 32 * 1024);

  ~readBuffer();

  bool      eof(void);
  bool      next(void);
  char      get(void);
  char      getnext(void);
  void      seek(size_t pos);
  size_t    read(char *buf, size_t len);   //  read the next len bytes into the user buffer buf
  size_t    tell(void);

private:
  void      fillBuffer(void);
  void      init(int fileptr, const char *filename, u32bit bufferMax);

  char    *_filename;

  int      _file;
  int      _fileType;  //  0=stdin, 1=normal, 2=mmap
  size_t   _filePos;

  bool     _eof;

  //  If bufferMax is zero, then we are using the mmapped interface, otherwise,
  //  we are using a open()/read() and a small buffer.

  size_t   _bufferPos;
  size_t   _bufferLen;
  size_t   _bufferMax;
  char    *_buffer;
};






inline
void
readBuffer::fillBuffer(void) {

  _bufferPos = 0;
  errno = 0;
  _bufferLen = (u32bit)::read(_file, _buffer, _bufferMax * sizeof(char));
  if (errno) {
#ifdef TRUE64BIT
    fprintf(stderr, "readBuffer()-- couldn't read %ld bytes from '%s'\n%s\n",
            _bufferMax * sizeof(char), _filename, strerror(errno));
#else
    fprintf(stderr, "readBuffer()-- couldn't read %d bytes from '%s'\n%s\n",
            (int)_bufferMax * sizeof(char), _filename, strerror(errno));
#endif
    exit(1);
  }
  if (_bufferLen == 0)
    _eof = true;
}


inline
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
      if (errno) {
        fprintf(stderr, "readBuffer()-- ERROR: couldn't open the file '%s'.\n%s\n",
                filename, strerror(errno));
        exit(1);
      }
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

    fillBuffer();
  }
}


inline
readBuffer::readBuffer(const char *filename, u32bit bufferMax) {
  init(-1, filename, bufferMax);
}


inline
readBuffer::readBuffer(int fileptr, const char *filename, u32bit bufferMax) {
  if (bufferMax == 0) {
    fprintf(stderr, "readBuffer::readBuffer()-- file '%s' requested memory mapped I/O on previously\n", filename);
    fprintf(stderr, "readBuffer::readBuffer()-- opened file, but that's not allowed.  bufferMax reset to 32KB.\n");
    bufferMax = 32 * 1024;
  }
  init(fileptr, filename, bufferMax);
}


inline
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
        fprintf(stderr, "readBuffer()-- WARNING: couldn't close the file '%s'\n%s\n",
                _filename, strerror(errno));
      }
      break;
    case 2:
      delete [] _filename;
      unmapFile(_buffer, _bufferLen);
      break;
  }
}


inline
bool
readBuffer::eof(void) {
  return(_eof);
}


inline
bool
readBuffer::next(void) {

  if (_eof)
    return(true);

  switch (_fileType) {
    case 0:
    case 1:
      _bufferPos++;
      _filePos++;

      if (_bufferPos >= _bufferLen)
        fillBuffer();
      break;
    case 2:
      _bufferPos++;
      if (_bufferPos >= _bufferLen)
        _eof = true;
      break;
  }

  return(_eof);
}


inline
char
readBuffer::get(void) {
  return(_buffer[_bufferPos]);
}


inline
char
readBuffer::getnext(void) {
  char  x = _buffer[_bufferPos];
  next();
  return(x);
}


inline
void
readBuffer::seek(size_t pos) {

  switch(_fileType) {
    case 0:
      fprintf(stderr, "readBuffer()-- seek() not available for file '%s'.\n", _filename);
      break;
    case 1:
      errno = 0;
      lseek(_file, pos, SEEK_SET);
      if (errno) {
        fprintf(stderr, "readBuffer()-- '%s' couldn't seek to position %d.\n%s\n",
                _filename, pos, strerror(errno));
        exit(1);
      }

      _filePos = pos;
      _eof     = false;

      //  We really do need to do a fillBuffer() here.  If the next
      //  operation is a get() we should return something valid.
      //
      //  Of course, we could put a check into get() to see if bufferPos
      //  is valid, but I don't want the overhead on EVERY get().
      //
      fillBuffer();
      break;
    case 2:
      if (pos < _bufferLen) {
        _bufferPos = pos;
        _eof       = false;
      } else {
        _eof       = true;
      }
      break;
  }
}


inline
size_t
readBuffer::read(char *buf, size_t len) {

  if (_fileType == 2) {
    size_t c = 0;

    while ((_bufferPos < _bufferLen) && (c < len))
      buf[c++] = _buffer[_bufferPos++];

    return(c);
  } else {
    //  The trick here is to use the existing buffered input first,
    //  then do a direct read to get the rest.
    //
    //  We fill the buffer again if it is empty.
    //
    //  The number of bytes actually put into buf is returned.
      
    size_t   bCopied = 0;   //  Number of bytes copied into the buffer
    size_t   bRead   = 0;   //  Number of bytes read into the buffer
    size_t   bAct    = 0;   //  Number of bytes actually read from disk

    //  Easy case; the next len bytes are already in the buffer; just
    //  copy and move the position.
    //
    //  XXX:  Check the zero-left-in-buffer case
    //
    if (_bufferLen - _bufferPos > len) {
      bCopied = len;
      bRead   = 0;

      memcpy(buf, _buffer + _bufferPos, sizeof(char) * len);
      _bufferPos += (u32bit)len;
    } else {

      //  Existing buffer not big enough.  Copy what's there, then finish
      //  with a read.
      //
      memcpy(buf, _buffer + _bufferPos, (_bufferLen - _bufferPos) * sizeof(char));
      bCopied    = _bufferLen - _bufferPos;
      _bufferPos = _bufferLen;

      while (bCopied + bRead < len) {
        errno = 0;
        bAct = (u32bit)::read(_file, buf + bCopied + bRead, (len - bCopied - bRead) * sizeof(char));
        if (errno) {
          fprintf(stderr, "readBuffer()-- couldn't read %d bytes from '%s'\n%s\n",
                  len * sizeof(char), _filename, strerror(errno));
          exit(1);
        }
        if (bAct == 0) {
          fprintf(stderr, "readBuffer()-- EOF hit while reading %d bytes from '%s'\n%s\n",
                  len * sizeof(char), _filename, strerror(errno));
          exit(1);
        }
        bRead += bAct;
      }
    }

    if (_bufferPos == _bufferLen)
      fillBuffer();

    //fprintf(stdout, "bCopied = %d  bRead = %d\n", bCopied, bRead);

    return(bCopied + bRead);
  }
}


inline
size_t
readBuffer::tell(void) {

  if (_fileType == 2)
    return(_bufferPos);

  return(_filePos);
}


#endif // READ_BUFFER_H

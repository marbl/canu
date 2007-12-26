#include "util++.H"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>

//  Define this to get reports of what the readBuffer does -- seek(),
//  fillBuffer() and read() are reported.
//#define VERBOSE_READBUFFER


readBuffer::readBuffer(const char *filename, u32bit bufferMax) {
  init(-1, filename, bufferMax);
}

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
      if (errno) {
        fprintf(stderr, "readBuffer()-- couldn't open the file '%s': %s\n",
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

    _buffer[0] = 0;

    fillBuffer();
  }
}


void
readBuffer::fillBuffer(void) {

#ifdef VERBOSE_READBUFFER
  fprintf(stderr, "readBuffer::fillBuffer()-- loading "u32bitFMT" bytes.\n", (u32bit)_bufferMax);
#endif

  _bufferPos = 0;
  errno = 0;
 again:
  _bufferLen = (u32bit)::read(_file, _buffer, _bufferMax * sizeof(char));
  if (errno == EAGAIN)
    goto again;
  if (errno) {
    fprintf(stderr, "readBuffer::fillBuffer()-- only read %d bytes, couldn't read %d bytes from '%s': %s\n",
            (int)_bufferLen, (int)(_bufferMax * sizeof(char)), _filename, strerror(errno));
    exit(1);
  }
  if (_bufferLen == 0)
    _eof = true;
}


bool
readBuffer::seek(off_t pos) {

#ifdef VERBOSE_READBUFFER
  fprintf(stderr, "readBuffer::seek()-- seek to "u64bitFMT"\n", (u64bit)pos);
#endif

  //  The file is currently at _filePos, and there are _bufferMax -
  //  _bufferPos bytes left in the buffer.  If we are seeking to
  //  something inside the buffer, the don't do a new seek/read,
  //  just move the buffer pointer.
  //
#if 0
  //  It's a nice idea, but it's broken.
  if (_fileType == 1) {
    if ((_filePos <= pos) && (pos <= _filePos + _bufferLen - _bufferPos - 256)) {
      u32bit offset = pos - _filePos;
      _filePos   += offset;
      _bufferPos += offset;
    }

    return(true);
  }
#endif

  switch(_fileType) {
    case 0:
      fprintf(stderr, "readBuffer()-- seek() not available for file '%s'.\n", _filename);
      return(false);
      break;
    case 1:
      errno = 0;
      lseek(_file, pos, SEEK_SET);
      if (errno) {
        fprintf(stderr, "readBuffer()-- '%s' couldn't seek to position "s64bitFMT": %s\n",
                _filename, pos, strerror(errno));
        exit(1);
        return(false);
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
      if (pos < (off_t)_bufferLen) {
        _bufferPos = pos;
        _eof       = false;
      } else {
        _eof       = true;
      }
      break;
  }

  return(true);
}


size_t
readBuffer::read(void *buf, size_t len) {
  char  *bufchar = (char *)buf;

#ifdef VERBOSE_READBUFFER
  fprintf(stderr, "readBuffer::read()-- returning "u32bitFMT" bytes.\n", (u32bit)len);
#endif

  if (_fileType == 2) {
    size_t c = 0;

    while ((_bufferPos < _bufferLen) && (c < len))
      bufchar[c++] = _buffer[_bufferPos++];

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

      memcpy(bufchar, _buffer + _bufferPos, sizeof(char) * len);
      _bufferPos += (u32bit)len;
    } else {

      //  Existing buffer not big enough.  Copy what's there, then finish
      //  with a read.
      //
      memcpy(bufchar, _buffer + _bufferPos, (_bufferLen - _bufferPos) * sizeof(char));
      bCopied    = _bufferLen - _bufferPos;
      _bufferPos = _bufferLen;

      while (bCopied + bRead < len) {
        errno = 0;
        bAct = (u32bit)::read(_file, bufchar + bCopied + bRead, (len - bCopied - bRead) * sizeof(char));
        if (errno) {
          fprintf(stderr, "readBuffer()-- couldn't read %d bytes from '%s': n%s\n",
                  (u32bit)len * sizeof(char), _filename, strerror(errno));
          exit(1);
        }

        //  If we hit EOF, return a short read
        if (bAct == 0) {
          len = 0;
        }
        bRead += bAct;
      }
    }

    if (_bufferPos == _bufferLen)
      fillBuffer();

    return(bCopied + bRead);
  }
}

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

#include "../libbri/libbritypes.H"

//  XXX:  TEMPORARY!
//
#define O_LARGEFILE 0

class readBuffer {
public:
  readBuffer(const char *filename, u32bit bufferMax = 32 * 1024);
  readBuffer(int fileptr, const char *filename, u32bit bufferMax = 32 * 1024);
  ~readBuffer();

  bool      eof(void);
  bool      next(void);
  char      get(void);
  void      seek(size_t pos);
  size_t    read(char *buf, size_t len);   //  read the next len bytes into the user buffer buf
  size_t    tell(void);

private:
  void      fillBuffer(void);
  void      init(int fileptr, const char *filename, u32bit bufferMax);

  char    *_filename;
  int      _file;
  int      _fileOwner;
  size_t   _filePos;
  bool     _eof;
  u32bit   _bufferPos;
  u32bit   _bufferLen;
  u32bit   _bufferMax;
  char    *_buffer;
};


inline
void
readBuffer::fillBuffer(void) {

  _bufferPos = 0;
  errno = 0;
  _bufferLen = (u32bit)::read(_file, _buffer, _bufferMax * sizeof(char));
  if (errno) {
    fprintf(stderr, "readBuffer()-- couldn't read %ld bytes from '%s'\n%s\n",
            _bufferMax * sizeof(char), _filename, strerror(errno));
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
  _file      = fileptr;
  _fileOwner = 0;
  _filePos   = 0;
  _eof       = false;
  _bufferPos = 0;          //  Position we are at in the buffer
  _bufferLen = 0;          //  Valid length of the buffer
  _bufferMax = bufferMax;  //  Maximum size of the buffer
  _buffer    = new char [_bufferMax];

  fillBuffer();
}


inline
readBuffer::readBuffer(const char *filename, u32bit bufferMax) {
  errno = 0;
  int fileptr = open(filename, O_RDONLY | O_LARGEFILE);
  if (errno) {
    fprintf(stderr, "readBuffer()-- ERROR: couldn't open the file '%s'.\n%s\n",
            filename, strerror(errno));
    exit(1);
  }

  init(fileptr, filename, bufferMax);
  _fileOwner = 1;
}


inline
readBuffer::readBuffer(int fileptr, const char *filename, u32bit bufferMax) {
  init(fileptr, filename, bufferMax);
  _fileOwner = 0;
}


inline
readBuffer::~readBuffer() {

  if (_fileOwner) {
    errno = 0;
    close(_file);
    if (errno) {
      fprintf(stderr, "readBuffer()-- WARNING: couldn't close the file '%s'\n%s\n",
              _filename, strerror(errno));
      //exit(1);
    }
  }
  delete [] _filename;
  delete [] _buffer;
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
    return(_eof);

  _bufferPos++;
  _filePos++;

  if (_bufferPos >= _bufferLen)
    fillBuffer();

  return(_eof);
}


inline
char
readBuffer::get(void) {
  return(_buffer[_bufferPos]);
}


inline
void
readBuffer::seek(size_t pos) {

#if 0
  //  version from original fasta reader -- I've never seen an EAGAIN,
  //  nor have I ever seen it fail.

 again:
  errno = 0;
  off_t lsr = lseek(_file, _theSeqs[iid]._headerStart, SEEK_SET);
  if (errno) {
    if (errno == EAGAIN) {
      fprintf(stderr, "%s: FastA::seek() got EAGAIN.  Attempting again.\n", _filename);
      goto again;
    }
    fprintf(stderr, "%s: FastA::seek() failed!\n%s\n", _filename, strerror(errno));
    exit(1);
  }
  if (lsr != _theSeqs[iid]._headerStart) {
    fprintf(stderr, "%s: FastA::seek() wanted pos %llu for seq %lu but got %llu??\n",
            _filename, _theSeqs[iid]._headerStart, iid, lsr);
    exit(1);
  }
#endif

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
}


inline
size_t
readBuffer::read(char *buf, size_t len) {

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


inline
size_t
readBuffer::tell(void) {
  return(_filePos);
}


#endif // READ_BUFFER_H

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
#ifndef O_LARGEFILE
#define O_LARGEFILE 0
#endif


//  Keeps two buffers, allows random-access, and forward/reverse
//  character-by-character.
//
//  Each buffer has a pointer to a potential previous / next buffer.
//  If the pointer is valid, then that buffer is valid.  If we hit the
//  end of the current, and the pointer is valid, we swap current to
//  that buffer.
//
//  If the pointer is not valid, we grab the backup buffer, make it
//  current (fill it), and set the pointer to get back.
//  

class readBuffer_buffer {
public:
  u32bit              _pos;
  u32bit              _len;
  u32bit              _max;
  char               *_buf;
  readBuffer_buffer  *_prev;
  readBuffer_buffer  *_next;

  readBuffer_buffer(u32bit max) {
    _pos  = 0;
    _len  = 0;
    _max  = max;
    _buf  = new char [sizeof(char) * _max];
    _prev = 0L;
    _next = 0L;
  };

  ~readBuffer_buffer() {
    delete [] _buf;
  };
};


class readBuffer {
public:
  readBuffer(const char *filename, u32bit bufferMax = 32 * 1024);
  readBuffer(int fileptr, const char *filename, u32bit bufferMax = 32 * 1024);
  ~readBuffer();

  bool      eof(void) { return(_eof); };
  bool      sof(void) { return(_sof); };
  bool      prev(void);
  bool      next(void);
  char      get(void);
  void      seek(size_t pos);
  size_t    read(char *buf, size_t len);   //  read the next len bytes into the user buffer buf
  size_t    tell(void);

private:
  void      fillBufferForward(void);
  void      fillBufferBackward(void);
  void      init(int fileptr, const char *filename, u32bit bufferMax);

  char    *_filename;
  int      _file;
  int      _fileOwner;
  size_t   _filePos;
  bool     _eof;
  bool     _sof;

  readBuffer_buffer  *_cur;  //  currently active buffer
  readBuffer_buffer  *_bak;  //  backup buffer
  readBuffer_buffer  *_tmp;  //  tmp pointer - used when we swap pointers

  void      swapBuffers(void) {
    _tmp = _cur;
    _cur = _bak;
    _bak = _tmp;
    _tmp = 0L;
  }
};


//  Exhausted the current buffer, and we want to move forward.
//
inline
void
readBuffer::fillBufferForward(void) {

  if (_eof)
    return;

  //  If no next buffer, fill the _bak buffer, then fall through to
  //  swap buffers.
  //
  if (_cur->_next == 0L) {
    errno = 0;
    _bak->_len = (u32bit)::read(_file, _bak->_buf, _bak->_max * sizeof(char));
    if (errno) {
      fprintf(stderr, "readBuffer()-- couldn't read %ld bytes from '%s'\n%s\n",
              _buf->_max * sizeof(char), _filename, strerror(errno));
      exit(1);
    }
    if (_bak->_len == 0)
      _eof = true;
  }

  //  The _bak buffer has the stuff we want, so make it the current.
  //  Reset the position, too.

  if (!_eof) {
    swapBuffers();

    _cur->_pos  = 0;
    _cur->_next = 0L;
    _cur->_prev = _bak;

    _bak->_next = _cur;
    _bak->_pref = 0L;
  }
}

//  Exhausted the current buffer, and we want to move forward.
//
inline
readBuffer::fillBufferBackward(void) {

  if (_sof)
    return;

  //  If no previous buffer, fill the _bak buffer, then fall through
  //  to swap buffers.
  //
  if (_cur->_prev == 0L) {

    //  We should always be a multiple of buffer::_max.  Check if we
    //  get unsynched somehow.
    //
    if (_filePos < _bak->_max) {
      fprintf(stderr, "ERROR: Got unsynched\n");
      exit(1);
    }

    //  Seek back _max from the current position, then read in the
    //  next _max bytes.

    errno = 0;
    lseek(_file, _filePos - _bak->_max, SEEK_SET);
    if (errno) {
      fprintf(stderr, "readBuffer()-- '%s' couldn't seek to position %d.\n%s\n",
              _filename, _filePos - _bak->_max, strerror(errno));
      exit(1);
    }

    errno = 0;
    _bak->_len = (u32bit)::read(_file, _bak->_buf, _bak->_max * sizeof(char));
    if (errno) {
      fprintf(stderr, "readBuffer()-- couldn't read %ld bytes from '%s'\n%s\n",
              _buf->_max * sizeof(char), _filename, strerror(errno));
      exit(1);
    }
    if (_bak->_len == 0)
      _eof = true;
  }

  //  The _bak buffer has the stuff we want, so make it the current.
  //  Reset the position, too.

  if (!_eof) {
    swapBuffers();

    _cur->_pos  = 0;
    _cur->_next = 0L;
    _cur->_prev = _bak;

    _bak->_next = _cur;
    _bak->_pref = 0L;
  }
}



inline
bool
readBuffer::next(void) {

  if (_eof)
    return(_eof);

  _cur->_pos++;
  _filePos++;

  if (_cur->_pos >= _cur->_len)
    fillBufferForward();

  _sof = false;

  return(_eof);
}


inline
bool
readBuffer::prev(void) {

  if (_filePos == 0)
    return(false);

  if (_cur->_pos == 0)
    fillBufferBackward();

  _cur->_pos--;
  _filePos--;

  _eof = false;

  return(true);
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

  //  Both buffers are now invalid.
  //
  _cur->_prev = 0L;
  _cur->_next = 0L;


  //
  //  XXX:  Should we fill both buffers, and place the position
  //  somewhere in the middle?
  //
  //  If we immediately do a prev(), we save one seek, and one read.
  //
  //  If we immediately do a next(), we lose one disk read.
  //



  //  We really do need to do a fillBuffer() here.  If the next
  //  operation is a get() we should return something valid.
  //
  //  Of course, we could put a check into get() to see if bufferPos
  //  is valid, but I don't want the overhead on EVERY get().
  //
  fillBufferForward();
}



//
//  XXXX:  The rest is hosed - not updated
//




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
    fillBufferForward();

  //fprintf(stdout, "bCopied = %d  bRead = %d\n", bCopied, bRead);

  return(bCopied + bRead);
}


inline
size_t
readBuffer::tell(void) {
  return(_filePos);
}


#endif // READ_BUFFER_H

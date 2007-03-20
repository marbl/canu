#ifndef READ_BUFFER_H
#define READ_BUFFER_H

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

  bool             eof(void);
  bool             next(void);
  char             get(void);
  char             getnext(void);
  bool             seek(off_t pos);
  size_t           read(char *buf, size_t len);
  off_t            tell(void);

private:
  void             fillBuffer(void);
  void             init(int fileptr, const char *filename, u32bit bufferMax);

  char           *_filename;

  int             _file;
  int             _fileType;  //  0=stdin, 1=normal, 2=mmap

  off_t           _filePos;

  bool            _eof;

  //  If bufferMax is zero, then we are using the mmapped interface, otherwise,
  //  we are using a open()/read() and a small buffer.

  size_t          _bufferPos;
  size_t          _bufferLen;
  size_t          _bufferMax;
  char           *_buffer;
};


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
off_t
readBuffer::tell(void) {

  if (_fileType == 2)
    return((size_t)_bufferPos);

  return(_filePos);
}


#endif  //  READ_BUFFER_H

#ifndef READ_BUFFER_H
#define READ_BUFFER_H


class readBuffer {
public:
  readBuffer(const char *filename, u32bit bufferMax = 32 * 1024);
  readBuffer(int fileptr, const char *filename, u32bit bufferMax = 32 * 1024);

  ~readBuffer();

  bool             eof(void) { return(_eof); };

  char             peek(void);

  char             read(void);
  size_t           read(void *buf, size_t len);

  void             seek(off_t pos);
  off_t            tell(void) { return(_filePos); };

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


//  Returns the next letter in the buffer, but DOES NOT advance past
//  it.  Might have some wierd interaction with EOF -- if you peek()
//  and the next thing is eof , the _eof flag might get set.
//
inline
char
readBuffer::peek(void) {

  if ((_eof == false) && (_fileType != 2) && (_bufferPos >= _bufferLen))
    fillBuffer();

  if (_eof)
    return(0);

  return(_buffer[_bufferPos]);
}


//  Returns the next letter in the buffer.  Returns EOF (0) if there
//  is no next letter.
//
inline
char
readBuffer::read(void) {

  if ((_eof == false) && (_fileType != 2) && (_bufferPos >= _bufferLen))
    fillBuffer();

  if (_bufferPos >= _bufferLen)
    _eof = true;

  if (_eof)
    return(0);

  _bufferPos++;
  _filePos++;

  return(_buffer[_bufferPos-1]);
}



#endif  //  READ_BUFFER_H

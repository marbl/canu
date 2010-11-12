#ifndef READ_BUFFER_H
#define READ_BUFFER_H


class readBuffer {
public:
  readBuffer(const char *filename, u64bit bufferMax = 32 * 1024);
  readBuffer(FILE *F, u64bit bufferMax = 32 * 1024);
  ~readBuffer();

  bool             eof(void) { return(_eof); };

  char             peek(void);

  char             read(void);
  u64bit           read(void *buf, u64bit len);
  u64bit           read(void *buf, u64bit maxlen, char stop);

  void             seek(u64bit pos);
  u64bit           tell(void) { return(_filePos); };

  const char      *filename(void) { return(_filename); };

private:
  void             fillBuffer(void);
  void             init(int fileptr, const char *filename, u64bit bufferMax);

  char           *_filename;

  int             _file;
  u64bit          _filePos;

  bool            _mmap;
  bool            _stdin;

  bool            _eof;

  //  If bufferMax is zero, then we are using the mmapped interface, otherwise,
  //  we are using a open()/read() and a small buffer.

  u64bit          _bufferPos;
  u64bit          _bufferLen;
  u64bit          _bufferMax;
  char           *_buffer;
};


//  Returns the next letter in the buffer, but DOES NOT advance past
//  it.  Might have some wierd interaction with EOF -- if you peek()
//  and the next thing is eof , the _eof flag might get set.
//
inline
char
readBuffer::peek(void) {

  if ((_eof == false) && (_bufferPos >= _bufferLen))
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

  if ((_eof == false) && (_bufferPos >= _bufferLen))
    fillBuffer();

  if (_eof)
    return(0);

  _bufferPos++;
  _filePos++;

  return(_buffer[_bufferPos-1]);
}


#endif  //  READ_BUFFER_H

#include "readBuffer.H"


inline
void
readBuffer::init(int fileptr, const char *filename, u32bit bufferMax) {

  _filename  = new char [strlen(filename) + 1];
  strcpy(_filename, filename);
  _file      = fileptr;
  _fileOwner = 0;
  _filePos   = 0;
  _eof       = false;
  _sof       = true;

  _cur       = new readBuffer_buffer(bufferMax);
  _bak       = new readBuffer_buffer(bufferMax);
  _tmp       = 0L;

  fillBufferForward();
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

  delete _cur;
  delete _bak;
}

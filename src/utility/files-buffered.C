
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-AUG-10
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "files.H"

#include <fcntl.h>



//  If bufferMax is zero, then the file is accessed using memory
//  mapped I/O.  Otherwise, a small buffer is used.
//
readBuffer::readBuffer(const char *filename, uint64 bufferMax) {

  memset(_filename, 0, sizeof(char) * (FILENAME_MAX + 1));

  if ((filename == 0L) ||
      ((filename[0] == '-') && (filename[1] == 0))) {
    strcpy(_filename, "(stdin)");
    _stdin = true;
  }

  else {
    strcpy(_filename, filename);
    _stdin = false;
  }

  _file        = 0;
  _filePos     = 0;

  _eof         = false;

  _bufferBgn   = 0;
  _bufferLen   = 0;

  _bufferPos   = 0;

  _bufferMax   = (bufferMax == 0) ? 32 * 1024 : bufferMax;
  _buffer      = new char [_bufferMax + 1];

  //  Open the file, failing if it's actually the terminal.

  errno = 0;
  _file = (_stdin) ? fileno(stdin) : open(_filename, O_RDONLY | O_LARGEFILE);
  if (errno)
    fprintf(stderr, "readBuffer()-- couldn't open the file '%s': %s\n",
            _filename, strerror(errno)), exit(1);

  if (isatty(_file)) {
    fprintf(stderr, "readBuffer()-- I cannot use the terminal for input.  Provide a filename or a pipe.\n");
    exit(1);
  }

  //  Fill the buffer.

  fillBuffer();
}



readBuffer::readBuffer(FILE *file, uint64 bufferMax) {

  memset(_filename, 0, sizeof(char) * (FILENAME_MAX + 1));
  strcpy(_filename, "(hidden file)");

  _file        = fileno(file);
  _filePos     = 0;

  _eof         = false;
  _stdin       = false;

  _bufferBgn   = 0;
  _bufferLen   = 0;

  _bufferPos   = 0;

  _bufferMax   = (bufferMax == 0) ? 32 * 1024 : bufferMax;
  _buffer      = new char [_bufferMax + 1];

  //  Rewind the file (allowing failure if it's a pipe or stdin).

  errno = 0;
  lseek(_file, 0, SEEK_SET);
  if ((errno) && (errno != ESPIPE))
    fprintf(stderr, "readBuffer()-- '%s' couldn't seek to position 0: %s\n",
            _filename, strerror(errno)), exit(1);

  //  Fill the buffer.

  fillBuffer();
}



readBuffer::~readBuffer() {

  delete [] _buffer;

  if (_stdin == false)
    close(_file);
}



void
readBuffer::fillBuffer(uint64 extra) {

  //  If there is still stuff in the buffer, no need to fill.

  if (_bufferPos + extra < _bufferLen)
    return;

  _bufferBgn += _bufferLen;
  _bufferLen  = 0;

  _bufferPos  = 0;

  assert(_filePos == _bufferBgn);

 again:
  errno = 0;
  _bufferLen = (uint64)::read(_file, _buffer, _bufferMax);

  if (errno == EAGAIN)
    goto again;

  if (errno)
    fprintf(stderr, "readBuffer::fillBuffer()-- only read " F_U64 " bytes, couldn't read " F_U64 " bytes from '%s': %s\n",
            _bufferLen, _bufferMax, _filename, strerror(errno)), exit(1);

  if (_bufferLen == 0)
    _eof = true;
}



void
readBuffer::seek(uint64 pos, uint64 extra) {

  //  If not really a seek, and the buffer still has enough stuff in it, just return.

  if ((pos == _filePos) && (_filePos + extra < _bufferBgn + _bufferLen))
    return;

  //  If stdin, we can't seek.

  if (_stdin == true) {
    fprintf(stderr, "readBuffer()-- seek() not available for file 'stdin'.\n");
    exit(1);
  }

  //  If the position is in the buffer, just move there and
  //  potentially skip any actual file access.

  else if ((pos < _filePos) &&
           (_bufferBgn  <= pos) &&
           (pos + extra <  _bufferBgn + _bufferLen)) {
    if (pos < _filePos) {
      //fprintf(stderr, "readBuffer::seek()-- jump back to position %lu from position %lu (buffer at %lu)\n",
      //        pos, _filePos, _bufferPos);
      _bufferPos -= (_filePos - pos);
      _filePos   -= (_filePos - pos);
    } else {
      //fprintf(stderr, "readBuffer::seek()-- jump ahead to position %lu from position %lu (buffer at %lu)\n",
      //        pos, _filePos, _bufferPos);
      _bufferPos += (pos - _filePos);
      _filePos   += (pos - _filePos);
    }
  }

  //  Nope, we need to grab a new block of data.

  else {
    //fprintf(stderr, "readBuffer::seek()-- jump directly to position %lu from position %lu (buffer at %lu)\n",
    //        pos, _filePos, _bufferPos);

    errno = 0;
    lseek(_file, pos, SEEK_SET);
    if (errno)
      fprintf(stderr, "readBuffer()-- '%s' couldn't seek to position " F_U64 ": %s\n",
              _filename, pos, strerror(errno)), exit(1);

    _filePos   = pos;

    _bufferBgn = pos;
    _bufferLen = 0;

    _bufferPos = 0;

    fillBuffer();
  }

  _eof = (_bufferPos >= _bufferLen);
}



uint64
readBuffer::read(void *buf, uint64 len) {
  char  *bufchar = (char *)buf;

  //  Easy case; the next len bytes are already in the buffer; just
  //  copy and move the position.

  if (_bufferPos + len <= _bufferLen) {
    memcpy(bufchar, _buffer + _bufferPos, len);

    _filePos   += len;
    _bufferPos += len;

    fillBuffer();

    return(len);
  }

  //  Existing buffer not big enough.  Copy what's there, then finish
  //  with a read.

  uint64   bCopied = 0;   //  Number of bytes copied into the buffer
  uint64   bAct    = 0;   //  Number of bytes actually read from disk

  bCopied     = _bufferLen - _bufferPos;

  memcpy(bufchar, _buffer + _bufferPos, bCopied);

  while (bCopied < len) {
    errno = 0;
    bAct = (uint64)::read(_file, bufchar + bCopied, len - bCopied);
    if (errno)
      fprintf(stderr, "readBuffer()-- couldn't read " F_U64 " bytes from '%s': n%s\n",
              len, _filename, strerror(errno)), exit(1);

    if (bAct == 0)    //  If we hit EOF, return a short read.
      len = 0;

    bCopied += bAct;
  }

  _filePos   += bCopied;    //  Advance the actual file position to however much we just read.

  _bufferBgn  = _filePos;   //  And set the buffer begin to that too.
  _bufferLen  = 0;          //  Set the buffer as empty, so we fill it again.

  _bufferPos  = 0;

  fillBuffer();

  return(bCopied);
}



uint64
readBuffer::read(void *buf, uint64 maxlen, char stop) {
  char  *bufchar = (char *)buf;
  uint64 c = 0;

  //  We will copy up to 'maxlen'-1 bytes into 'buf', or stop at the first occurrence of 'stop'.
  //  This will reserve space at the end of any string for a zero-terminating byte.
  maxlen--;

  while ((_eof == false) && (c < maxlen)) {
    bufchar[c++] = _buffer[_bufferPos];

    _filePos++;
    _bufferPos++;

    if (_bufferPos >= _bufferLen)
      fillBuffer();

    if (bufchar[c-1] == stop)
      break;
  }

  bufchar[c] = 0;

  return(c);
}



bool
readBuffer::peekIFFchunk(char name[4], uint32 &dataLen) {

  //  Seek to the current position, making sure there are at least
  //  8 bytes still in the buffer.

  seek(_filePos, 8);

  //  If there's space for a valid IFF header, return the name and length.

  if (_bufferPos + 8 < _bufferLen) {
    memcpy( name,    _buffer + _bufferPos,     sizeof(char) * 4);
    memcpy(&dataLen, _buffer + _bufferPos + 4, sizeof(uint32));

    return(true);
  }

  //  If not, return an empty name and length of zero.

  name[0] = 0;
  name[1] = 0;
  name[2] = 0;
  name[3] = 0;

  dataLen = 0;

  return(false);
}





void
readBuffer::readIFFchunk(char*name, uint8 *&data, uint32 &dataLen,  uint32 &dataMax) {

  //  Read the name and data length.

  read( name,    4);
  read(&dataLen, sizeof(uint32));

  //  Allocate space for the data.

  resizeArray(data, 0, dataMax, dataLen, resizeArray_doNothing);

  //  Copy the data to 'data'.

  read(data, dataLen);
}




writeBuffer::writeBuffer(const char *filename, const char *filemode, uint64 bufferMax) {
  strncpy(_filename, filename, FILENAME_MAX);
  strncpy(_filemode, filemode, 16);

  _file    = NULL;
  _filePos = 0;

  if      (filemode[0] == 'a')           //  If appending, open the file now
    open();                              //  so we can set the file position.
  else if (filemode[0] != 'w')           //  Otherwise, if not writing, fail.
    fprintf(stderr, "writeBuffer()--  Unknown mode '%s'\n", filemode), exit(1);

  _bufferLen      = 0;
  _bufferMax      = bufferMax;
  _buffer         = new char [_bufferMax];

  _chunkBufferLen = 0;
  _chunkBufferMax = 0;
  _chunkBuffer    = NULL;

  _chunkStartsLen = 0;
  _chunkStartsMax = 0;
  _chunkStarts    = NULL;
  _chunkSizes     = NULL;
}



writeBuffer::~writeBuffer() {
  flush();

  delete [] _buffer;

  delete [] _chunkBuffer;

  delete [] _chunkStarts;
  delete [] _chunkSizes;

  AS_UTL_closeFile(_file, _filename);
}



void
writeBuffer::write(void *data, uint64 length) {

  if (_bufferMax < _bufferLen + length)           //  Flush the buffer if this
    flush();                                      //  data is too big for it.

  if (_bufferMax < length) {                      //  And if it is still too big
    assert(_bufferLen == 0);                      //  (ensure the buffer is empty)
    writeToDisk(data, length);                    //  and just dump it to disk.
  }

  else {                                          //  Otherwise, copy it to
    memcpy(_buffer + _bufferLen, data, length);   //  our buffer.
    _bufferLen += length;
  }

  assert(_bufferLen <= _bufferMax);

  _filePos += length;
}



//  If dataLength is zero, only a header is written, and the chunk stack
//  is increased.
//
//  If dataLength is non-zero, and the stack is empty, the block
//  is immediately written to disk.
//
//  If dataLength is non-zero, and the stack has entries, the block
//  is copied to our internal chunk buffer and the size added
//  to ALL container chunks.
//
void
writeBuffer::writeIFFchunk(char *name, void *data, uint32 dataLength) {
  uint8  header [4 * sizeof(uint8) + sizeof(uint32)] = { 0 };
  uint8  padding[4 * sizeof(uint8)]                  = { 0 };

  //  Figure out how much padding we need to add to the data to make it
  //  align on a 32-bit boundary.

  uint32 headLength = 4 * sizeof(uint8) + sizeof(uint32);
  uint32 padLength  = 4 - (dataLength % 4);

  if (padLength == 4)
    padLength = 0;

  //  Create the chunk header.

  dataLength += padLength;

  memcpy(header + 0,  name,       sizeof(uint8) * 4);
  memcpy(header + 4, &dataLength, sizeof(uint32));

  dataLength -= padLength;

  //  If the chunk size is zero, add a new container chunk (adding space for
  //  8 more if needed).
  if (dataLength == 0) {
    increaseArrayPair(_chunkStarts, _chunkSizes, _chunkStartsLen, _chunkStartsMax, 8);

    _chunkStarts[_chunkStartsLen] = _chunkBufferLen;
    _chunkSizes [_chunkStartsLen] = 0;

    _chunkStartsLen++;

    appendIFFdata(header, headLength);
  }

  //  But if there is no container, we can immediately output the chunk.
  else if (_chunkStartsLen == 0) {
    write(header,  headLength);
    write(data,    dataLength);
    write(padding, padLength);
  }

  //  Otherwise, append the chunk to our buffer and increase
  //  the size of each parent container.
  else {
    appendIFFdata(header,  headLength);
    appendIFFdata(data,    dataLength);
    appendIFFdata(padding, padLength);

    for (uint32 ii=0; ii<_chunkStartsLen; ii++)
      _chunkSizes[ii] += headLength + dataLength + padLength;
  }
}



void
writeBuffer::appendIFFdata(void *data, uint32 dataLength) {

  if (dataLength == 0)
    return;

  //  If this data will exceed our current allocation, double what we've
  //  allocated until it'll fit.

  if (_chunkBufferLen + dataLength > _chunkBufferMax) {
    uint64  newMax = (_chunkBufferMax == 0) ? 16384 : _chunkBufferMax;

    while (newMax < _chunkBufferLen + dataLength)
      newMax *= 2;

    resizeArray(_chunkBuffer, _chunkBufferLen, _chunkBufferMax, newMax);
  }

  //  Copy the data into the buffer and update the length.

  memcpy(_chunkBuffer + _chunkBufferLen, data, dataLength);

  _chunkBufferLen += dataLength;
}



//  The user is done with this chunk.
//
//  If the stack has entries, set the length of the last
//  chunk and pop the stack.  If the stack is now
//  empty, write the block.
//
void
writeBuffer::closeIFFchunk(char *name) {

  //  If no chunk to close, report an error.

  if (_chunkStartsLen == 0) {
    fprintf(stderr, "writeBuffer::closeIFFchunk()-- no chunk to close.\n");
    exit(1);
    return;
  }

  //  Refer to the last chunkStarts entry.

  _chunkStartsLen--;

  //  If a name supplied, check that it's the same as the chunk we're
  //  closing.

  if (name) {
    uint64  cs = _chunkStarts[_chunkStartsLen];

    if ((name[0] != _chunkBuffer[cs + 0]) ||
        (name[1] != _chunkBuffer[cs + 1]) ||
        (name[2] != _chunkBuffer[cs + 2]) ||
        (name[3] != _chunkBuffer[cs + 3])) {
      fprintf(stderr, "writeBuffer::closeIFFchunk()-- requested to close chunk '%c%c%c%c' but current chunk is '%c%c%c%c'.\n",
              name[0], name[1], name[2], name[3],
              _chunkBuffer[cs + 0],
              _chunkBuffer[cs + 1],
              _chunkBuffer[cs + 2],
              _chunkBuffer[cs + 3]);
      exit(1);
    }
  }

  //  Update the size of the chunk container we're in.

  *(uint32 *)(_chunkBuffer + _chunkStarts[_chunkStartsLen] + 4) = _chunkSizes[_chunkStartsLen];

  //  If there are no more containers, write this buffer to disk, and clear it.

  if (_chunkStartsLen == 0) {
    write(_chunkBuffer, _chunkBufferLen);

    _chunkBufferLen = 0;
  }
}




void
writeBuffer::open(void) {
  if (_file != NULL)
    return;

  errno = 0;
  _file = fopen(_filename, _filemode);
  if (errno)
    fprintf(stderr, "writeBuffer()--  Failed to open file '%s' with mode '%s': %s\n",
            _filename, _filemode, strerror(errno)), exit(1);

  //  If appending, _filePos is zero, and ftell() is non-zero.
  //  If writing, _filePos is non-zero, and ftell() is zero.
  _filePos += AS_UTL_ftell(_file);
}



void
writeBuffer::writeToDisk(void *data, uint64 length) {
  if (length == 0)
    return;

  open();
  writeToFile((char *)data, "writeBuffer::writeToDisk", length, _file);
}



void
writeBuffer::flush(void) {
  writeToDisk(_buffer, _bufferLen);
  _bufferLen = 0;
}

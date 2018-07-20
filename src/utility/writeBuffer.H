
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef WRITE_BUFFER_H
#define WRITE_BUFFER_H

#include "AS_global.H"
#include "files.H"

class writeBuffer {
public:
  writeBuffer(const char *filename, const char *filemode, uint64 bufferMax = 1024 * 1024) {
    strncpy(_filename, filename, FILENAME_MAX);
    strncpy(_filemode, filemode, 16);

    _file    = NULL;
    _filePos = 0;

    if      (filemode[0] == 'a')           //  If appending, open the file now
      open();                              //  so we can set the file position.
    else if (filemode[0] != 'w')           //  Otherwise, if not writing, fail.
      fprintf(stderr, "writeBuffer()--  Unknown mode '%s'\n", filemode), exit(1);

    _bufferLen = 0;
    _bufferMax = bufferMax;
    _buffer    = new char [_bufferMax];
  };

  ~writeBuffer() {
    flush();
    delete [] _buffer;
    AS_UTL_closeFile(_file, _filename);
  };

  const char          *filename(void) { return(_filename); };
  uint64               tell(void)     { return(_filePos);  };

  void                 write(void *data, uint64 length) {

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
  };


private:
  void                 open(void) {
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
  };

  void                 writeToDisk(void *data, uint64 length) {
    if (length == 0)
     return;

    open();
    AS_UTL_safeWrite(_file, data, "writeBuffer::writeToDisk", 1, length);
  };

  void                 flush(void) {
    writeToDisk(_buffer, _bufferLen);
    _bufferLen = 0;
  };

  char                _filename[FILENAME_MAX+1];
  char                _filemode[17];

  FILE               *_file;
  uint64              _filePos;

  uint64              _bufferLen;
  uint64              _bufferMax;
  char               *_buffer;
};

#endif  //  WRITE_BUFFER_H

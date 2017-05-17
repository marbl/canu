
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
 *    Brian P. Walenz beginning on 2016-NOV-18
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef WRITE_BUFFER_H
#define WRITE_BUFFER_H

#include "AS_global.H"
#include "AS_UTL_fileIO.H"

class writeBuffer {
public:
  writeBuffer(const char *filename, const char *filemode, uint64 bufferMax = 1024 * 1024) {
    _filename = filename;
    _filemode = filemode;

    errno = 0;
    _file    = fopen(filename, filemode);
    if (errno)
      fprintf(stderr, "writeBuffer()--  Failed to open file '%s' with mode '%s': %s\n",
              filename, filemode, strerror(errno)), exit(1);

    _filePos = AS_UTL_ftell(_file);

    _bufferLen = 0;
    _bufferMax = bufferMax;
    _buffer    = new char [_bufferMax];
  };

  ~writeBuffer() {
    flush();
    delete [] _buffer;
    fclose(_file);
  };

  uint64               tell(void) { return(_filePos); };

  void                 write(void *data, uint64 length) {
    if ((_bufferMax < length) ||
        (_bufferLen + length > _bufferMax))
      flush();

    if (_bufferMax < length)
      AS_UTL_safeWrite(_file, data, "writeBuffer", 1, length);
    else {
      memcpy(_buffer + _bufferLen, data, length);
      _bufferLen += length;
    }

    _filePos += length;
  };

  const char          *filename(void) { return(_filename); };

private:
  void                 flush(void) {
    AS_UTL_safeWrite(_file, _buffer, "writeBuffer", 1, _bufferLen);
    _bufferLen = 0;
  };

  const char         *_filename;
  const char         *_filemode;

  FILE               *_file;
  uint64              _filePos;

  uint64              _bufferLen;
  uint64              _bufferMax;
  char               *_buffer;
};

#endif  //  WRITE_BUFFER_H

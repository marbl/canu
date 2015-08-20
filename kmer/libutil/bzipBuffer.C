
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
 *    Brian P. Walenz from 2006-JUN-23 to 2014-APR-11
 *      are Copyright 2006,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "util++.H"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>


//  This is probably correct, it just cannot read a normal *.bz file;
//  it probably reads an unpackaged raw bzip stream.


bzipBuffer::bzipBuffer(const char *filename, uint32 bufferMax) {

  _filename  = new char [strlen(filename) + 1];
  strcpy(_filename, filename);

  if (bufferMax == 0)
    bufferMax = 32 * 1024;

  errno = 0;
  _file = open(filename, O_RDONLY | O_LARGEFILE);
  if (errno) {
    fprintf(stderr, "bzipBuffer()-- couldn't open the file '%s': %s\n",
            filename, strerror(errno));
    exit(1);
  }

  _filePos     = 0;
  _eof         = false;

  _bzip2bufferMax = bufferMax;
  _bzip2inPos     = 0;
  _bzip2outPos    = 0;

  _bzip2in  = new char [_bzip2bufferMax];
  _bzip2out = new char [_bzip2bufferMax];

  _bzip2streamEnd  = false;

  _bzip2stream.next_in         = _bzip2in;
  _bzip2stream.avail_in        = 0;
  _bzip2stream.total_in_lo32   = 0;
  _bzip2stream.total_in_hi32   = 0;
  _bzip2stream.next_out        = _bzip2out;
  _bzip2stream.avail_out       = 0;
  _bzip2stream.total_out_lo32  = 0;
  _bzip2stream.total_out_hi32  = 0;
  _bzip2stream.state           = 0L;
  _bzip2stream.bzalloc         = 0L;
  _bzip2stream.bzfree          = 0L;
  _bzip2stream.opaque          = 0L;

  int res = BZ2_bzDecompressInit(&_bzip2stream, 0, 0);
  if (res != BZ_OK) {
    //  BZ_CONFIG_ERROR, BZ_PARAM_ERROR, BZ_MEM_ERROR
    fprintf(stderr, "bzipBuffer::bzipBuffer()--  Failed to initialize the decompressor.\n");
    exit(1);
  }

  fillBuffer();
}


bzipBuffer::~bzipBuffer() {
  delete [] _bzip2in;
  delete [] _bzip2out;
  close(_file);
}


void
bzipBuffer::fillBuffer(void) {

  if (_bzip2streamEnd) {
    _eof = true;
    return;
  }

  //  Scream and holler if the bzip2 buffer isn't exhausted!
  //
  if (_bzip2outPos < _bzip2stream.avail_out) {
    fprintf(stderr, "bzipBuffer::fillBuffer()--  Buffer isn't empty!  Still %d bytes!\n",
            (int)(_bzip2stream.avail_out - _bzip2outPos));
    return;
  }

  _bzip2outPos = 0;

 again:

  //  If there is stuff in the input, run the decompressor.  If it
  //  decompresses anything, return.
  //
  if (_bzip2stream.avail_in > 0) {

    fprintf(stderr, "about to decompress %d bytes in input\n", (int)_bzip2stream.avail_in);
    fprintf(stderr, "in  is bzip2:%p and real:%p (diff %d)\n", _bzip2stream.next_in,  _bzip2in,  _bzip2stream.next_in - _bzip2in);
    fprintf(stderr, "out is bzip2:%p and real:%p (diff %d)\n", _bzip2stream.next_out, _bzip2out, _bzip2stream.next_out - _bzip2out);


    int res = BZ2_bzDecompress(&_bzip2stream);
    if (res == BZ_STREAM_END) {
      fprintf(stderr, "GOT STREAM END!\n");

      BZ2_bzDecompressEnd(&_bzip2stream);

      _bzip2streamEnd = true;
      res = BZ_OK;
    }
    if (res != BZ_OK) {
      fprintf(stderr, "bzipBuffer::fillBuffer()--  Failed to decompress.\n"), exit(1);
    }

    fprintf(stderr, "decompressed %d bytes; still have %d in input\n", (int)_bzip2stream.avail_out, (int)_bzip2stream.avail_in);
    fprintf(stderr, "in  is bzip2:%p and real:%p (diff %d)\n", _bzip2stream.next_in,  _bzip2in,  _bzip2stream.next_in - _bzip2in);
    fprintf(stderr, "out is bzip2:%p and real:%p (diff %d)\n", _bzip2stream.next_out, _bzip2out, _bzip2stream.next_out - _bzip2out);

    if (_bzip2stream.avail_out > 0) {
      fprintf(stderr, "----------------------------------------\n");
      fwrite(_bzip2stream.next_out, sizeof(char), _bzip2stream.avail_out, stderr);
      fprintf(stderr, "\n----------------------------------------\n");
      return;
    }
  }

  //  If we're here and _bzip2streamEnd is true, we hit the end of the
  //  stream at the same time we hit the end of the input data.
  //
  if (_bzip2streamEnd) {
    _eof = true;
    return;
  }

  //  Otherwise, we need to read some input.
  //
  errno = 0;
  _bzip2stream.next_in   = _bzip2in;
  _bzip2stream.avail_in  = (uint32)::read(_file, _bzip2in, sizeof(char) * _bzip2bufferMax);
  _bzip2stream.next_out  = _bzip2out;
  _bzip2stream.avail_out = _bzip2bufferMax;
  if (errno) {
    fprintf(stderr, "bzipBuffer::fillBuffer()-- read failed: %s\n", strerror(errno));
    exit(1);
  }

  fprintf(stderr, "read %d bytes\n", (int)_bzip2stream.avail_in);

  if (_bzip2stream.avail_in == 0) {
    fprintf(stderr, "bzipBuffer::fillBuffer()-- hit end of file?\n");
    _eof = true;
    return;
  }

  //  And now try to decompress it again
  //
  goto again;
}


bool
bzipBuffer::seek(off_t pos) {
  fprintf(stderr, "bzipBuffer()-- seek() not available for file '%s'.\n", _filename);
  return(false);
}


size_t
bzipBuffer::read(char *buf, size_t len) {

#if 0
  if (_fileType == 2) {
    size_t c = 0;

    while ((_bufferPos < _bufferLen) && (c < len))
      buf[c++] = _buffer[_bufferPos++];

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

      memcpy(buf, _buffer + _bufferPos, sizeof(char) * len);
      _bufferPos += (uint32)len;
    } else {

      //  Existing buffer not big enough.  Copy what's there, then finish
      //  with a read.
      //
      memcpy(buf, _buffer + _bufferPos, (_bufferLen - _bufferPos) * sizeof(char));
      bCopied    = _bufferLen - _bufferPos;
      _bufferPos = _bufferLen;

      while (bCopied + bRead < len) {
        errno = 0;
        bAct = (uint32)::read(_file, buf + bCopied + bRead, (len - bCopied - bRead) * sizeof(char));
        if (errno) {
          fprintf(stderr, "bzipBuffer()-- couldn't read %d bytes from '%s': n%s\n",
                  (uint32)len * sizeof(char), _filename, strerror(errno));
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
#endif

  return(0);
}

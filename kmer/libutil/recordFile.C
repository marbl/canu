
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
 *    Brian P. Walenz from 2008-JUL-08 to 2014-APR-11
 *      are Copyright 2008,2014 J. Craig Venter Institute, and
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

//  N.B. any read() / write() pair (either order) must have a seek (or
//  a fflush) in between.

uint64   recordFileMagic1 = 0x694664726f636572llu;
uint64   recordFileMagic2 = 0x000000000000656cllu;

recordFile::recordFile(char const *name,
                       uint32      headerSize,
                       uint32      recordSize,
                       char        mode) {

  _file = 0;
  _name = new char [strlen(name) + 1];
  strcpy(_name, name);

  _numRecords = 0;
  _recordSize   = recordSize;

  _headerSize   = headerSize;
  _header       = new char [_headerSize];

  memset(_header, 0, sizeof(char) * _headerSize);

  _bfrmax       = MAX(1048576 / _recordSize, 16);
  _bfr          = new char [_bfrmax * _recordSize];

  _limit        = ~uint32ZERO;

  _pos          = uint64ZERO;
  _rec          = 0;

  memset(_bfr, 0, sizeof(char) * _bfrmax * _recordSize);

  _bfrDirty     = false;
  _isReadOnly   = true;

  if ((mode != 'r') && (mode != 'w') && (mode |= 'a')) {
    fprintf(stderr, "recordFile::recordFile()--  Invalid mode '%c'.\n", mode);
    exit(1);
  }

  //  If the file doesn't exist, or we're opening for write, we're
  //  basically done.  Do that first.
  //    Write the magic.
  //    Write the metadata.
  //    Write the header.

  if (((mode == 'w')) ||
      ((mode == 'a') && (fileExists(_name) == false))) {
    errno = 0;
    _file = open(_name,
                 O_RDWR | O_CREAT | O_TRUNC | O_LARGEFILE,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno)
      fprintf(stderr, "recordFile::recordFile()-- failed to open '%s': %s\n",
              _name, strerror(errno)), exit(1);
    _isReadOnly = false;

    write(_file, &recordFileMagic1,  sizeof(uint64));
    write(_file, &recordFileMagic2,  sizeof(uint64));
    write(_file, &_numRecords,       sizeof(uint64));
    write(_file, &_recordSize,       sizeof(uint32));
    write(_file, &_headerSize,       sizeof(uint32));
    write(_file,  _header,           sizeof(char) * _headerSize);

    if (errno)
      fprintf(stderr, "recordFile::recordFile()-- failed to write header to '%s': %s\n",
              _name, strerror(errno)), exit(1);

    return;
  }

  //  File does exist.  If we're not appending, open it read-only.
  //  Otherwise, open read-write.

  if (mode == 'r') {
    errno = 0;
    _file = open(_name,
                 O_RDONLY | O_LARGEFILE,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno)
      fprintf(stderr, "recordFile::recordFile()-- failed to open '%s': %s\n",
              _name, strerror(errno)), exit(1);
    _isReadOnly = true;
  } else {
    errno = 0;
    _file = open(_name,
                 O_RDWR | O_LARGEFILE,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno)
      fprintf(stderr, "recordFile::recordFile()-- failed to open for write '%s': %s\n",
              _name, strerror(errno)), exit(1);
    _isReadOnly = false;
  }

  //  Read the magic, metadata and header.

  {
    uint64 m1, m2;

    errno = 0;

    read(_file, &m1,                sizeof(uint64));
    read(_file, &m2,                sizeof(uint64));
    read(_file, &_numRecords,       sizeof(uint64));
    read(_file, &_recordSize,       sizeof(uint32));
    read(_file, &_headerSize,       sizeof(uint32));
    read(_file,  _header,           sizeof(char) * _headerSize);

    if (errno)
      fprintf(stderr, "recordFile::recordFile()-- failed to read header from '%s': %s\n",
              _name, strerror(errno)), exit(1);

    if ((m1 != recordFileMagic1) || (m2 != recordFileMagic2))
      fprintf(stderr, "recordFile::recordFile()-- magic number disagreement; '%s' not a recordFile?\n",
              _name), exit(1);
  }

  if (mode == 'a') {
    _pos = _numRecords;
    _rec = 0;

    errno = 0;
    lseek(_file, 0, SEEK_END);
    if (errno)
      fprintf(stderr, "recordFile::recordFile()-- seek to end of '%s' failed: %s\n", _name, strerror(errno)), exit(1);
  } else {
    seek(0, true);
  }
}


recordFile::~recordFile() {
  flushDirty();

  if (_isReadOnly == false) {
    errno = 0;
    lseek(_file, 0, SEEK_SET);
    if (errno)
      fprintf(stderr, "recordFile::~recordFile()-- seek to start of '%s' failed: %s\n", _name, strerror(errno)), exit(1);

    write(_file, &recordFileMagic1,  sizeof(uint64));
    write(_file, &recordFileMagic2,  sizeof(uint64));
    write(_file, &_numRecords,       sizeof(uint64));
    write(_file, &_recordSize,       sizeof(uint32));
    write(_file, &_headerSize,       sizeof(uint32));
    write(_file,  _header,           sizeof(char) * _headerSize);

    if (errno)
      fprintf(stderr, "recordFile::~recordFile()-- failed to write header to '%s': %s\n",
              _name, strerror(errno)), exit(1);
  }

  close(_file);

  if (errno)
    fprintf(stderr, "recordFile::~recordFile()-- failed to close '%s': %s\n",
            _name, strerror(errno)), exit(1);

  delete [] _bfr;
  delete [] _name;
  delete [] _header;
}



//  If the page is dirty, flush it to disk
//
void
recordFile::flushDirty(void) {

  if (_bfrDirty == false)
    return;

  if (_isReadOnly)
    fprintf(stderr, "recordFile::recordFile()-- '%s' is readonly, but is dirty!\n", _name), exit(1);

  errno = 0;
  lseek(_file, 32 + _headerSize + _pos * _recordSize, SEEK_SET);
  if (errno)
    fprintf(stderr, "recordFile::seek()-- '%s' failed: %s\n", _name, strerror(errno)), exit(1);

  //  Write records up to, not including, _rec.  Unlike the
  //  bitPackedFile, there is no issue with partially filled words
  //  here.
  //
  errno = 0;
  write(_file, _bfr, _recordSize * _rec);
  if (errno)
    fprintf(stderr, "recordFile::write()-- '%s' failed: %s\n", _name, strerror(errno)), exit(1);

  _bfrDirty = false;
}



//  Seeks to rec in the file, reads in a new block.
//
void
recordFile::seek(uint64 rec, bool forced) {

  //  If we are seeking to somewhere in the current block, don't do a
  //  real seek, just move our position within the block.
  //
  if ((forced == false) && (_pos <= rec) && (rec < _pos + _bfrmax)) {
    _rec = rec - _pos;
    return;
  }

  flushDirty();

  _pos = rec;  //  Root of buffer is now here
  _rec = 0;    //  See?

  errno = 0;
  lseek(_file, 32 + _headerSize + _pos * _recordSize, SEEK_SET);
  if (errno)
    fprintf(stderr, "recordFile::seek() '%s' seek to record="uint64FMT" at fileposition="uint64FMT" failed: %s\n",
            _name, _pos, _headerSize + _pos * _recordSize, strerror(errno)), exit(1);

  errno = 0;
  read(_file, _bfr, _recordSize * _bfrmax);
  if (errno)
    fprintf(stderr, "recordFile::seek() '%s' read of "uint64FMT" bytes failed at record "uint64FMT", fileposition "uint64FMT"': %s\n",
            _name, _recordSize * _bfrmax, _pos, _headerSize + _pos * _recordSize, strerror(errno)), exit(1);
}



uint32
recordFile::getRecord(void *record, uint32 num) {
  uint32  maxnum  = _bfrmax / 2;

  //  Reading large blocks -- bigger than the in-core size?  Loop and
  //  recurse.
  //
  if (num > maxnum) {
    uint32  numread = 0;
    uint32  pos = 0;
    uint32  len = 0;

    while (num > 0) {
      len = MIN(maxnum, num);
      len = getRecord((char *)record + pos * _recordSize, len);

      if (len == 0)
        return(numread);

      num     -= len;
      pos     += len;
      numread += len;
    }

    return(numread);
  }

  //  If asked to read too many records, read whatever is left.
  //
  if (_numRecords < _pos + _rec + num)
    num = _numRecords - _pos - _rec;
  if (_limit      < _pos + _rec + num)
    num = _limit      - _pos - _rec;

  //  If the current position is already past eof, return without
  //  reading.  The previous 'if' ensures we will never read a block
  //  past eof.
  //
  if ((_numRecords < _pos + _rec) || (_limit < _pos + _rec))
    return(0);

  if (_bfrmax < _rec + num + 1)
    seek(_pos + _rec, true);

  memcpy(record, _bfr + _rec * _recordSize, _recordSize * num);

  _rec += num;

  return(num);
}



void
recordFile::putRecord(void *record, uint32 num) {
  uint32  maxnum = _bfrmax / 2;

  if (num > maxnum) {
    uint32  pos = 0;
    uint32  len = 0;

    while (num > 0) {
      len = MIN(maxnum, num);

      putRecord((char *)record + pos * _recordSize, len);

      num -= len;
      pos += len;
    }

  } else {
    if (_bfrmax < _rec + num + 1)
      seek(_pos + _rec, true);

    memcpy(_bfr + _rec * _recordSize, record, _recordSize * num);

    _rec        += num;
    _numRecords += num;

    _bfrDirty = true;
  }
}

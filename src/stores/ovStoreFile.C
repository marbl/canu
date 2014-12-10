
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static const char *rcsid = "$Id$";

#include "ovStore.H"


ovFile::ovFile(const char  *name,
               ovFileType   type,
               uint32       bufferSize) {

  //  We write two sizes of overlaps.  The 'normal' format doesn't contain the a_iid, while the
  //  'full' format does.  Choose a buffer size that can handle both, because we don't know
  //  which type is going to be used.
  //
  //  The lcm can be one of two values, depending on the number of bits in the reead length:
  //    == 16-bit -- (4 + 5*4) * (8 + 5*4) = 672
  //     > 17-bit -- (4 + 3*8) * (8 + 3*8) = 896

  uint32  lcm = ((sizeof(uint32) * 1) + (sizeof(ovsOverlapDAT)) *
                 (sizeof(uint32) * 2) + (sizeof(ovsOverlapDAT)));

  if (bufferSize < 16 * 1024)
    bufferSize = 16 * 1024;

  _bufferLen  = 0;
  _bufferPos  = (bufferSize / (lcm * sizeof(uint32))) * lcm;  //  Forces reload on next read
  _bufferMax  = (bufferSize / (lcm * sizeof(uint32))) * lcm;
  _buffer     = new uint32 [_bufferMax];
  _isNormal   = (type == ovFileNormal);

  //  The buffer size must hold an integer number of overlaps, otherwise the reader
  //  will read partial overlaps and fail.

  assert(_bufferMax % ((sizeof(uint32) * 1) + (sizeof(ovsOverlapDAT))) == 0);
  assert(_bufferMax % ((sizeof(uint32) * 2) + (sizeof(ovsOverlapDAT))) == 0);

  //  Open a file for reading?
  if ((type & ovFileWrite) == 0) {
    _reader      = new compressedFileReader(name);
    _writer      = NULL;
    _file        = _reader->file();
    _isSeekable  = (_reader->isCompressed() == false);
    _isOutput    = false;
  }

  //  Open a file for writing?
  else {
    _reader      = NULL;
    _writer      = new compressedFileWriter(name);
    _file        = _writer->file();
    _isSeekable  = false;
    _isOutput    = true;
  }
}



ovFile::~ovFile() {

  flushOverlaps();

  delete    _reader;
  delete    _writer;
  delete [] _buffer;
}



void
ovFile::flushOverlaps(void) {

  if (_isOutput == false)
    return;

  if (_bufferLen == 0)
    return;

  AS_UTL_safeWrite(_file, _buffer, "ovFile::flushOverlaps", sizeof(uint32), _bufferLen);

  _bufferLen = 0;
}



void
ovFile::writeOverlap(ovsOverlap *overlap) {

  assert(_isOutput == true);

  if (_bufferLen >= _bufferMax) {
    AS_UTL_safeWrite(_file, _buffer, "ovFile::writeOverlap", sizeof(uint32), _bufferLen);
    _bufferLen = 0;
  }

  if (_isNormal == false)
    _buffer[_bufferLen++] = overlap->a_iid;

  _buffer[_bufferLen++] = overlap->b_iid;

#if (ovsOverlapNWORDS == 5)
  _buffer[_bufferLen++] = overlap->dat.dat[0];
  _buffer[_bufferLen++] = overlap->dat.dat[1];
  _buffer[_bufferLen++] = overlap->dat.dat[2];
  _buffer[_bufferLen++] = overlap->dat.dat[3];
  _buffer[_bufferLen++] = overlap->dat.dat[4];
#else
  _buffer[_bufferLen++] = (overlap->dat.dat[0] >> 32) & 0xffffffff;
  _buffer[_bufferLen++] = (overlap->dat.dat[0] >>  0) & 0xffffffff;
  _buffer[_bufferLen++] = (overlap->dat.dat[1] >> 32) & 0xffffffff;
  _buffer[_bufferLen++] = (overlap->dat.dat[1] >>  0) & 0xffffffff;
  _buffer[_bufferLen++] = (overlap->dat.dat[2] >> 32) & 0xffffffff;
  _buffer[_bufferLen++] = (overlap->dat.dat[2] >>  0) & 0xffffffff;
#endif

  assert(_bufferLen <= _bufferMax);
}



uint32
ovFile::readOverlap(ovsOverlap *overlap) {

  assert(_isOutput == false);

  if (_bufferPos >= _bufferLen) {
    _bufferLen = AS_UTL_safeRead(_file, _buffer, "ovFile::readOverlap", sizeof(uint32), _bufferMax);
    _bufferPos = 0;
  }

  if (_bufferLen == 0)
    return(false);

  assert(_bufferPos < _bufferLen);

  if (_isNormal == FALSE)
    overlap->a_iid      = _buffer[_bufferPos++];

  overlap->b_iid      = _buffer[_bufferPos++];

#if (ovsOverlapNWORDS == 5)
  overlap->dat.dat[0] = _buffer[_bufferPos++];
  overlap->dat.dat[1] = _buffer[_bufferPos++];
  overlap->dat.dat[2] = _buffer[_bufferPos++];
  overlap->dat.dat[3] = _buffer[_bufferPos++];
  overlap->dat.dat[4] = _buffer[_bufferPos++];
#else
  overlap->dat.dat[0]  = _buffer[_bufferPos++];  overlap->dat.dat[0] <<= 32;
  overlap->dat.dat[0] |= _buffer[_bufferPos++];
  overlap->dat.dat[1]  = _buffer[_bufferPos++];  overlap->dat.dat[1] <<= 32;
  overlap->dat.dat[1] |= _buffer[_bufferPos++];
  overlap->dat.dat[2]  = _buffer[_bufferPos++];  overlap->dat.dat[2] <<= 32;
  overlap->dat.dat[2] |= _buffer[_bufferPos++];
#endif

  assert(_bufferPos <= _bufferLen);

  return(true);
}



//  Move to the correct spot, and force a load on the next readOverlap by setting the position to
//  the end of the buffer.
void
ovFile::seekOverlap(off_t overlap) {
  off_t   ovlSize = sizeof(uint32) * ((_isNormal) ? 1 : 2) + sizeof(ovsOverlapWORD) * ovsOverlapNWORDS;

  if (_isSeekable == false)
    fprintf(stderr, "ovFile::seekOverlap()-- can't seek.\n"), exit(1);

  AS_UTL_fseek(_file, overlap * ovlSize, SEEK_SET);

  _bufferPos = _bufferLen;
}

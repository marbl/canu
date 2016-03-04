
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
 *  This file is derived from:
 *
 *    src/AS_OVS/AS_OVS_overlapFile.C
 *    src/AS_OVS/AS_OVS_overlapFile.c
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2007-FEB-28 to 2013-SEP-22
 *      are Copyright 2007-2009,2011-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2011-MAR-31
 *      are Copyright 2011 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Gregory Sims on 2012-FEB-01
 *      are Copyright 2012 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-09 to 2015-JUN-16
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "ovStore.H"


ovFile::ovFile(const char  *name,
               ovFileType   type,
               uint32       bufferSize) {

  //  We write two sizes of overlaps.  The 'normal' format doesn't contain the a_iid, while the
  //  'full' format does.  The buffer size must hold an integer number of overlaps, otherwise the
  //  reader will read partial overlaps and fail.  Choose a buffer size that can handle both.

  uint32  lcm = ((sizeof(uint32) * 1 + sizeof(ovOverlapDAT)) *
                 (sizeof(uint32) * 2 + sizeof(ovOverlapDAT)));

  if (bufferSize < 16 * 1024)
    bufferSize = 16 * 1024;

  _bufferLen  = 0;
  _bufferPos  = (bufferSize / (lcm * sizeof(uint32))) * lcm;  //  Forces reload on next read
  _bufferMax  = (bufferSize / (lcm * sizeof(uint32))) * lcm;
  _buffer     = new uint32 [_bufferMax];

  assert(_bufferMax % ((sizeof(uint32) * 1) + (sizeof(ovOverlapDAT))) == 0);
  assert(_bufferMax % ((sizeof(uint32) * 2) + (sizeof(ovOverlapDAT))) == 0);

  //  When writing full overlaps, we also write the number of overlaps per read.  This is used to
  //  build the store.  Overlaps in the store, normal format, don't need this extra data, as the
  //  store itself knows how many overlaps per read.

  _olapsPerReadAlloc = 0;
  _olapsPerReadLast  = 0;
  _olapsPerRead      = NULL;

  if (type == ovFileFullWrite) {
    _olapsPerReadAlloc = 128 * 1024;
    _olapsPerReadLast  = 0;
    _olapsPerRead      = new uint32 [_olapsPerReadAlloc];

    memset(_olapsPerRead, 0, sizeof(uint32) * _olapsPerReadAlloc);
  }

  //  Create the input/output buffers and files.

  _isOutput   = false;
  _isSeekable = false;
  _isNormal   = (type == ovFileNormal) || (type == ovFileNormalWrite);

  _reader     = NULL;
  _writer     = NULL;

  //  Open a file for reading?
  if ((type == ovFileNormal) || (type == ovFileFull)) {
    _reader      = new compressedFileReader(name);
    _file        = _reader->file();
    _isSeekable  = (_reader->isCompressed() == false);
  }

  //  Open a file for writing?
  else {
    _writer      = new compressedFileWriter(name);
    _file        = _writer->file();
    _isOutput    = true;
  }

  //  Make a copy of the output name, and clean it up.  This is used as the base for
  //  the counts output.  We just strip off all the dotted extensions in the filename.

  strcpy(_prefix, name);

  char  *slash = strrchr(_prefix, '/');
  char  *dot   = strchr((slash == NULL) ? _prefix : slash, '.');

  if (dot)
    *dot = 0;
}



ovFile::~ovFile() {

  flushOverlaps();

  delete    _reader;
  delete    _writer;
  delete [] _buffer;

  if (_olapsPerRead) {
    char  name[FILENAME_MAX];

    sprintf(name, "%s.counts", _prefix);

    errno = 0;
    _file = fopen(name, "w");
    if (errno)
      fprintf(stderr, "failed to open counts file '%s' for writing: %s\n", name, strerror(errno)), exit(1);

    _olapsPerReadLast++;

    AS_UTL_safeWrite(_file, &_olapsPerReadLast, "ovFile::olapsPerReadLast", sizeof(uint32), 1);
    AS_UTL_safeWrite(_file,  _olapsPerRead,     "ovFile::olapsPerRead",     sizeof(uint32), _olapsPerReadLast);

    fclose(_file);

    delete [] _olapsPerRead;

    //fprintf(stderr, "Wrote counts file '%s' for reads up to iid "F_U32"\n", name, _olapsPerReadLast);
  }
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
ovFile::writeOverlap(ovOverlap *overlap) {

  assert(_isOutput == true);

  if (_bufferLen >= _bufferMax) {
    AS_UTL_safeWrite(_file, _buffer, "ovFile::writeOverlap", sizeof(uint32), _bufferLen);
    _bufferLen = 0;
  }

  if (_olapsPerRead) {
    uint32   newmax  = _olapsPerReadAlloc;
    uint32   newlast = _olapsPerReadLast;

    newlast = max(newlast, overlap->a_iid);
    newlast = max(newlast, overlap->b_iid);

    while (newmax <= newlast)
      newmax += newmax / 4;

    resizeArray(_olapsPerRead, _olapsPerReadLast+1, _olapsPerReadAlloc, newmax, resizeArray_copyData | resizeArray_clearNew);

    _olapsPerRead[overlap->a_iid]++;
    _olapsPerRead[overlap->b_iid]++;

    _olapsPerReadLast = newlast;
  }

  if (_isNormal == false)
    _buffer[_bufferLen++] = overlap->a_iid;

  _buffer[_bufferLen++] = overlap->b_iid;

#if (ovOverlapNWORDS == 5)
  _buffer[_bufferLen++] = overlap->dat.dat[0];
  _buffer[_bufferLen++] = overlap->dat.dat[1];
  _buffer[_bufferLen++] = overlap->dat.dat[2];
  _buffer[_bufferLen++] = overlap->dat.dat[3];
  _buffer[_bufferLen++] = overlap->dat.dat[4];
#elif (ovOverlapNWORDS == 3)
  _buffer[_bufferLen++] = (overlap->dat.dat[0] >> 32) & 0xffffffff;
  _buffer[_bufferLen++] = (overlap->dat.dat[0] >>  0) & 0xffffffff;
  _buffer[_bufferLen++] = (overlap->dat.dat[1] >> 32) & 0xffffffff;
  _buffer[_bufferLen++] = (overlap->dat.dat[1] >>  0) & 0xffffffff;
  _buffer[_bufferLen++] = (overlap->dat.dat[2] >> 32) & 0xffffffff;
  _buffer[_bufferLen++] = (overlap->dat.dat[2] >>  0) & 0xffffffff;
#elif (ovOverlapNWORDS == 8)
  _buffer[_bufferLen++] = overlap->dat.dat[0];
  _buffer[_bufferLen++] = overlap->dat.dat[1];
  _buffer[_bufferLen++] = overlap->dat.dat[2];
  _buffer[_bufferLen++] = overlap->dat.dat[3];
  _buffer[_bufferLen++] = overlap->dat.dat[4];
  _buffer[_bufferLen++] = overlap->dat.dat[5];
  _buffer[_bufferLen++] = overlap->dat.dat[6];
  _buffer[_bufferLen++] = overlap->dat.dat[7];
#else
#error unknown ovOverlapNWORDS
#endif

  assert(_bufferLen <= _bufferMax);
}



void
ovFile::writeOverlaps(ovOverlap *overlaps, uint64 overlapsLen) {
  uint64  nWritten = 0;

  assert(_isOutput == true);

  //  Resize the olapsPerRead array once per batch.

  if (_olapsPerRead) {
    uint32  newmax  = _olapsPerReadAlloc;
    uint32  newlast = _olapsPerReadLast;

    for (uint32 oo=0; oo<overlapsLen; oo++) {
      newlast = max(newlast, overlaps[oo].a_iid);
      newlast = max(newlast, overlaps[oo].b_iid);

      while (newmax <= newlast)
        newmax += newmax / 4;
    }

    resizeArray(_olapsPerRead, _olapsPerReadLast+1, _olapsPerReadAlloc, newmax, resizeArray_copyData | resizeArray_clearNew);

    _olapsPerReadLast = newlast;
  }

  //  Add all overlaps to the buffer.

  while (nWritten < overlapsLen) {
    if (_bufferLen >= _bufferMax) {
      AS_UTL_safeWrite(_file, _buffer, "ovFile::writeOverlap", sizeof(uint32), _bufferLen);
      _bufferLen = 0;
    }

    if (_olapsPerRead) {
      _olapsPerRead[overlaps[nWritten].a_iid]++;
      _olapsPerRead[overlaps[nWritten].b_iid]++;
    }

    if (_isNormal == false)
      _buffer[_bufferLen++] = overlaps[nWritten].a_iid;

    _buffer[_bufferLen++] = overlaps[nWritten].b_iid;

#if (ovOverlapNWORDS == 5)
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[0];
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[1];
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[2];
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[3];
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[4];
#elif (ovOverlapNWORDS == 3)
    _buffer[_bufferLen++] = (overlaps[nWritten].dat.dat[0] >> 32) & 0xffffffff;
    _buffer[_bufferLen++] = (overlaps[nWritten].dat.dat[0] >>  0) & 0xffffffff;
    _buffer[_bufferLen++] = (overlaps[nWritten].dat.dat[1] >> 32) & 0xffffffff;
    _buffer[_bufferLen++] = (overlaps[nWritten].dat.dat[1] >>  0) & 0xffffffff;
    _buffer[_bufferLen++] = (overlaps[nWritten].dat.dat[2] >> 32) & 0xffffffff;
    _buffer[_bufferLen++] = (overlaps[nWritten].dat.dat[2] >>  0) & 0xffffffff;
#elif (ovOverlapNWORDS == 8)
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[0];
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[1];
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[2];
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[3];
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[4];
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[5];
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[6];
    _buffer[_bufferLen++] = overlaps[nWritten].dat.dat[7];
#else
#error unknown ovOverlapNWORDS
#endif

    nWritten++;
  }

  assert(_bufferLen <= _bufferMax);
}



bool
ovFile::readOverlap(ovOverlap *overlap) {

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

#if (ovOverlapNWORDS == 5)
  overlap->dat.dat[0] = _buffer[_bufferPos++];
  overlap->dat.dat[1] = _buffer[_bufferPos++];
  overlap->dat.dat[2] = _buffer[_bufferPos++];
  overlap->dat.dat[3] = _buffer[_bufferPos++];
  overlap->dat.dat[4] = _buffer[_bufferPos++];
#elif (ovOverlapNWORDS == 3)
  overlap->dat.dat[0]  = _buffer[_bufferPos++];  overlap->dat.dat[0] <<= 32;
  overlap->dat.dat[0] |= _buffer[_bufferPos++];
  overlap->dat.dat[1]  = _buffer[_bufferPos++];  overlap->dat.dat[1] <<= 32;
  overlap->dat.dat[1] |= _buffer[_bufferPos++];
  overlap->dat.dat[2]  = _buffer[_bufferPos++];  overlap->dat.dat[2] <<= 32;
  overlap->dat.dat[2] |= _buffer[_bufferPos++];
#elif (ovOverlapNWORDS == 8)
  overlap->dat.dat[0] = _buffer[_bufferPos++];
  overlap->dat.dat[1] = _buffer[_bufferPos++];
  overlap->dat.dat[2] = _buffer[_bufferPos++];
  overlap->dat.dat[3] = _buffer[_bufferPos++];
  overlap->dat.dat[4] = _buffer[_bufferPos++];
  overlap->dat.dat[5] = _buffer[_bufferPos++];
  overlap->dat.dat[6] = _buffer[_bufferPos++];
  overlap->dat.dat[7] = _buffer[_bufferPos++];
#else
#error unknown ovOverlapNWORDS
#endif

  assert(_bufferPos <= _bufferLen);

  return(true);
}



uint64
ovFile::readOverlaps(ovOverlap *overlaps, uint64 overlapsLen) {
  uint64  nLoaded = 0;

  assert(_isOutput == false);

  while (nLoaded < overlapsLen) {
    if (_bufferPos >= _bufferLen) {
      _bufferLen = AS_UTL_safeRead(_file, _buffer, "ovFile::readOverlaps", sizeof(uint32), _bufferMax);
      _bufferPos = 0;
    }

    if (_bufferLen == 0)
      return(nLoaded);

    assert(_bufferPos < _bufferLen);

    if (_isNormal == FALSE)
      overlaps[nLoaded].a_iid      = _buffer[_bufferPos++];

    overlaps[nLoaded].b_iid      = _buffer[_bufferPos++];

#if (ovOverlapNWORDS == 5)
    overlaps[nLoaded].dat.dat[0] = _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[1] = _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[2] = _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[3] = _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[4] = _buffer[_bufferPos++];
#elif (ovOverlapNWORDS == 3)
    overlaps[nLoaded].dat.dat[0]  = _buffer[_bufferPos++];  overlaps[nLoaded].dat.dat[0] <<= 32;
    overlaps[nLoaded].dat.dat[0] |= _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[1]  = _buffer[_bufferPos++];  overlaps[nLoaded].dat.dat[1] <<= 32;
    overlaps[nLoaded].dat.dat[1] |= _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[2]  = _buffer[_bufferPos++];  overlaps[nLoaded].dat.dat[2] <<= 32;
    overlaps[nLoaded].dat.dat[2] |= _buffer[_bufferPos++];
#elif (ovOverlapNWORDS == 8)
    overlaps[nLoaded].dat.dat[0] = _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[1] = _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[2] = _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[3] = _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[4] = _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[5] = _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[6] = _buffer[_bufferPos++];
    overlaps[nLoaded].dat.dat[7] = _buffer[_bufferPos++];
#else
#error unknown ovOverlapNWORDS
#endif

    nLoaded++;

    assert(_bufferPos <= _bufferLen);
  }

  return(nLoaded);
}



//  Move to the correct spot, and force a load on the next readOverlap by setting the position to
//  the end of the buffer.
void
ovFile::seekOverlap(off_t overlap) {

  if (_isSeekable == false)
    fprintf(stderr, "ovFile::seekOverlap()-- can't seek.\n"), exit(1);

  AS_UTL_fseek(_file, overlap * recordSize(), SEEK_SET);

  _bufferPos = _bufferLen;  //  We probably need to reload the buffer.
}


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
 *    kmer/libutil/bitPackedFile.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-JAN-02 to 2003-OCT-20
 *      are Copyright 2003 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-MAR-29 to 2004-APR-21
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUL-12 to 2014-APR-11
 *      are Copyright 2005-2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef BITPACKEDFILE_H
#define BITPACKEDFILE_H

#include "AS_global.H"

#include "bitOperations.H"
#include "bitEncodings.H"
#include "bitPacking.H"

//#define WITH_BZIP2

#ifdef WITH_BZIP2
#include <bzlib.h>
#endif

class bitPackedFile {
public:
  bitPackedFile(char const *name, uint64 offset=0, bool forceTruncate=false);
  ~bitPackedFile();

  uint64     getBits(uint32 size);
  uint64     getNumber(void);

  void       putBits(uint64 bits, uint32 size);
  void       putNumber(uint64 val);

  uint64     tell(void)       { return((_pos << 6) + _bit); };
  void       seek(uint64 pos);

  uint64     loadInCore(void);

  void       showStats(FILE *f) {
    fprintf(f, "inside: "F_U64"  outside: "F_U64"\n", stat_seekInside, stat_seekOutside);
    fflush(f);
  };
private:

  //  Ensure that the buffer has enough space for any future
  //  operation.  This constant, currently 31 bytes, must be strictly
  //  less than the constant used in deciding if seek() is moving
  //  forward or backwards.
  //
  void       sync(void) {
    if (((_bit >> 6) + 31) >= _bfrmax)
      seek((_pos << 6) + _bit);
  };

  void       flushDirty(void);
  void       seekBzip2(uint64 bitpos);
  void       seekNormal(uint64 bitpos);

  int       _file;
  char     *_name;

#ifdef WITH_BZIP2
  FILE     *_bzFILE;
  int       _bzerr;
  BZFILE   *_bzfile;
#endif

  uint64    _bfrmax;  //  Number of words in the buffer
  uint64   *_bfr;     //  A chunk of the bitPackedFile in core
  uint64    _pos;     //  The location this chunk is from (in words)
  uint64    _bit;     //  The bit position we are modifying relative to _pos

  bool      _inCore;
  bool      _bfrDirty;
  bool      _forceFirstLoad;
  bool      _isReadOnly;
  bool      _isBzip2;

  //  For collecting statistics on our usage
  //
  uint64  stat_seekInside;
  uint64  stat_seekOutside;
  uint64  stat_dirtyFlushes;

  //  For converting between hardware of different endianess.
  //
  uint64  file_offset;
  uint64  endianess_offset;
  bool    endianess_flipped;
};


inline
uint64
bitPackedFile::getBits(uint32 siz) {
  sync();
  uint64 ret = getDecodedValue(_bfr, _bit, siz);
  _bit += siz;
  return(ret);
}

inline
uint64
bitPackedFile::getNumber(void) {
  sync();
  uint64 siz = 0;
  uint64 ret = getFibonacciEncodedNumber(_bfr, _bit, &siz);
  _bit += siz;
  return(ret);
}


inline
void
bitPackedFile::putBits(uint64 bits, uint32 siz) {
  assert(_isReadOnly == false);
  sync();
  setDecodedValue(_bfr, _bit, siz, bits);
  _bit += siz;
  _bfrDirty = true;
}

inline
void
bitPackedFile::putNumber(uint64 val) {
  assert(_isReadOnly == false);
  sync();
  uint64 siz = 0;
  setFibonacciEncodedNumber(_bfr, _bit, &siz, val);
  _bit += siz;
  _bfrDirty = true;
}


#endif  //  BITPACKEDFILE_H

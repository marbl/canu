#ifndef BITPACKEDFILE_H
#define BITPACKEDFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

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
    fprintf(f, "inside: "uint64FMT"  outside: "uint64FMT"\n", stat_seekInside, stat_seekOutside);
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

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
  bitPackedFile(char const *name, u64bit offset=0, bool forceTruncate=false);
  ~bitPackedFile();

  u64bit     getBits(u32bit size);
  u64bit     getNumber(void);

  void       putBits(u64bit bits, u32bit size);
  void       putNumber(u64bit val);

  u64bit     tell(void)       { return((_pos << 6) + _bit); };
  void       seek(u64bit pos);

  u64bit     loadInCore(void);

  void       showStats(FILE *f) {
    fprintf(f, "inside: "u64bitFMT"  outside: "u64bitFMT"\n", stat_seekInside, stat_seekOutside);
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
  void       seekBzip2(u64bit bitpos);
  void       seekNormal(u64bit bitpos);

  int       _file;
  char     *_name;

#ifdef WITH_BZIP2
  FILE     *_bzFILE;
  int       _bzerr;
  BZFILE   *_bzfile;
#endif

  u64bit    _bfrmax;  //  Number of words in the buffer
  u64bit   *_bfr;     //  A chunk of the bitPackedFile in core
  u64bit    _pos;     //  The location this chunk is from (in words)
  u64bit    _bit;     //  The bit position we are modifying relative to _pos

  bool      _inCore;
  bool      _bfrDirty;
  bool      _forceFirstLoad;
  bool      _isReadOnly;
  bool      _isBzip2;

  //  For collecting statistics on our usage
  //
  u64bit  stat_seekInside;
  u64bit  stat_seekOutside;
  u64bit  stat_dirtyFlushes;

  //  For converting between hardware of different endianess.
  //
  u64bit  file_offset;
  u64bit  endianess_offset;
  bool    endianess_flipped;
};


inline
u64bit
bitPackedFile::getBits(u32bit siz) {
  sync();
  u64bit ret = getDecodedValue(_bfr, _bit, siz);
  _bit += siz;
  return(ret);
}

inline
u64bit
bitPackedFile::getNumber(void) {
  sync();
  u64bit siz = 0;
  u64bit ret = getFibonacciEncodedNumber(_bfr, _bit, &siz);
  _bit += siz;
  return(ret);
}


inline
void
bitPackedFile::putBits(u64bit bits, u32bit siz) {
  assert(_isReadOnly == false);
  sync();
  setDecodedValue(_bfr, _bit, siz, bits);
  _bit += siz;
  _bfrDirty = true;
}

inline
void
bitPackedFile::putNumber(u64bit val) {
  assert(_isReadOnly == false);
  sync();
  u64bit siz = 0;
  setFibonacciEncodedNumber(_bfr, _bit, &siz, val);
  _bit += siz;
  _bfrDirty = true;
}


#endif  //  BITPACKEDFILE_H

#ifndef BITPACKEDFILE_H
#define BITPACKEDFILE_H


class bitPackedFile {
public:
  bitPackedFile(char const *name, u64bit offset=0);
  ~bitPackedFile();

  u64bit     getBits(u32bit size);
  u64bit     getBits(u64bit pos, u32bit size);
  u64bit     getNumber(void);

  void       putBits(u64bit bits, u32bit size);
  void       putBits(u64bit pos, u64bit bits, u32bit size);
  void       putNumber(u64bit val);

  void       seek(u64bit pos, bool forceLoad=false);

  void       showStats(FILE *f) {
    fprintf(f, "inside: "u64bitFMT"  outside: "u64bitFMT"\n", stat_seekInside, stat_seekOutside);
    fflush(f);
  };
private:

  //  Ensure that the buffer has enough space for any future
  //  operation.
  //
  //  This constant, currently 31 bytes, must be strictly less than
  //  the constant used in deciding if seek() is moving forward or
  //  backwards.
  //
  void       sync(void) {
    if (((_bit >> 6) + 31) >= _bfrmax)
      seek((_pos << 6) + _bit, true);
  };

  void       flushDirty(void);

  int       _file;

  u64bit   *_bfr;     //  A chunk of the bitPackedFile in core
  u64bit    _bfrmax;  //  Number of words in the buffer
  u64bit    _pos;     //  The location this chunk is from (in words)
  u64bit    _bit;     //  The bit position we are modifying relative to _pos

  bool      _bfrDirty;
#ifdef __APPLE__
  u32bit    _pad;
#endif

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
void
bitPackedFile::putBits(u64bit bits, u32bit siz) {
  sync();
  setDecodedValue(_bfr, _bit, siz, bits);
  _bit += siz;
  _bfrDirty = true;
}


inline
u64bit
bitPackedFile::getBits(u64bit pos, u32bit siz) {
  seek(pos);
  u64bit ret = getDecodedValue(_bfr, _bit, siz);
  return(ret);
}

inline
void
bitPackedFile::putBits(u64bit pos, u64bit bits, u32bit siz) {
  seek(pos);
  setDecodedValue(_bfr, _bit, siz, bits);
  _bfrDirty = true;
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
bitPackedFile::putNumber(u64bit val) {
  sync();
  u64bit siz = 0;
  setFibonacciEncodedNumber(_bfr, _bit, &siz, val);
  _bit += siz;
  _bfrDirty = true;
}


#endif  //  BITPACKEDFILE_H

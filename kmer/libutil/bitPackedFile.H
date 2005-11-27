#ifndef BITPACKEDFILE_H
#define BITPACKEDFILE_H

////////////////////////////////////////
//
//  bitPackedFile
//
class bitPackedFile {
public:
  bitPackedFile(char const *name, u64bit offset=0);
  ~bitPackedFile();

  u64bit     getBits(u32bit size);
  u64bit     getNumber(void);

  void       putBits(u64bit bits, u32bit size);
  void       putNumber(u64bit val);

  void       seek(u64bit pos, bool forceLoad=false);

  void       showStats(FILE *f) {
    fprintf(f, "inside: "u64bitFMT"  outside: "u64bitFMT"\n", stat_seekInside, stat_seekOutside);
    fflush(f);
  };
private:

  //  Ensure that the buffer has enough space for any future
  //  operation.  If the buffer doesn't have 128 bits free (2 words
  //  currently) flush it (if dirty) and read in a new page.  128 was
  //  chosen because the fibonacci encoded numbers use up to 90-some
  //  bits.
  //
  void       sync(void) {

    //  This constant must be strictly less than the constant used in deciding
    //  if seek() is moving forward or backwards.
    //
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


#endif  //  BITPACKEDFILE_H

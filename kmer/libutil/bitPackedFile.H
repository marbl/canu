#ifndef BITPACKEDFILE_H
#define BITPACKEDFILE_H

#include <stdio.h>

#include "libbritypes.h"

class bitPackedFileWriter {
public:
  bitPackedFileWriter(char const *name);
  ~bitPackedFileWriter();

  void       putBits(u64bit bits, u32bit size);
  void       putNumber(u64bit val);
private:
  void       flush(void);

  FILE     *_out;
  bool      _pipe;
  u64bit   *_bfr;
  u64bit    _bit;
  //u32bit    _sta;
};

class bitPackedFileReader {
public:
  bitPackedFileReader(char const *name);
  ~bitPackedFileReader();

  u64bit     getBits(u32bit size);
  u64bit     getNumber(void);

  void       seek(u64bit pos);
private:
  void       fillBufferFromDisk(u32bit bufferpos);
  void       fill(void);

  FILE     *_in;
  bool      _pipe;
  u64bit   *_bfr;
  u64bit    _bit;
  //u32bit    _sta;
};

#endif  //  BITPACKEDFILE_H

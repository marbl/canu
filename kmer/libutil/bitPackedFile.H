#ifndef BITPACKEDFILE_H
#define BITPACKEDFILE_H

#include <stdio.h>
#include "libbri.H"


class bitPackedFileWriter {
public:
  bitPackedFileWriter(char *name);
  ~bitPackedFileWriter();

  void       putBits(u64bit bits, u32bit size);
private:
  FILE     *_out;
  bool      _pipe;
  u64bit   *_bfr;
  u64bit    _bit;
  u32bit    _sta;
};


class bitPackedFileReader {
public:
  bitPackedFileReader(char *name);
  ~bitPackedFileReader();

  u64bit     getBits(u32bit size);
private:
  FILE     *_in;
  bool      _pipe;
  u64bit   *_bfr;
  u64bit    _bit;
  u32bit    _sta;
};


#endif  //  BITPACKEDFILE_H

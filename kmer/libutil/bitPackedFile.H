#ifndef BITPACKEDFILE_H
#define BITPACKEDFILE_H

#include <stdio.h>

#include "libbritypes.h"

//
//  XXX: WARNING: The bitpacked file only works if the start is
//  aligned to 64-bit words.  If you, say, have a text file with a
//  bitPackedFile tacked on the end, it usually won't work.
//
//  [ 63 bytes of text][bitPackedFile]
//
//  The symptom will be a short read on the last block of the
//  bitPackedFile.
//
//  This is only caused by _externally_ joining the files.
//

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
};

class bitPackedFileReader {
public:
  bitPackedFileReader(char const *name);
  ~bitPackedFileReader();

  u64bit     getBits(u32bit size);
  u64bit     getNumber(void);

  //  Seek to the bit position 'pos'
  void       seek(u64bit pos);
private:
  void       fillBufferFromDisk(u32bit bufferpos);
  void       fill(void);

  FILE     *_in;
  bool      _pipe;
  u64bit   *_bfr;
  u64bit    _bit;
};

#endif  //  BITPACKEDFILE_H

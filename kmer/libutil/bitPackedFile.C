#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#include "bitPackedFile.H"

#define BUFFER_SIZE   131072

bitPackedFileWriter::bitPackedFileWriter(char *name) {
  _out = fopen(name, "wb");
  _bfr = new u64bit [BUFFER_SIZE];
  _bit = u64bitZERO;
}

bitPackedFileWriter::~bitPackedFileWriter() {
  u64bit wd = (_bit >> 6) & 0x0000cfffffffffffllu;

  errno = 0;
  fwrite(_bfr, sizeof(u64bit), wd+1, _out);
  if (errno) {
    fprintf(stderr, "bitPackedFileWriter::~bitPackedFileWriter got %s\n", strerror(errno));
    exit(1);
  }

  delete [] _bfr;
  fclose(_out);
}

//  Essentially setDecodedValue, but it will update _WORD and _bit
//  after the operation, and call flushBuffer if needed.
//
void
bitPackedFileWriter::putBits(u64bit bits, u32bit size) {
  u64bit wd = (_bit >> 6) & 0x0000cfffffffffffllu;
  u64bit bt = (_bit     ) & 0x000000000000003fllu;
  u64bit b1 = 64 - bt;
  u64bit b2 = size - b1;  //  Only used if siz > b1

  assert(size < 65);
  assert(size > 0);

  _bit += size;

  //  Make sure that the bits to write are the correct size!  This
  //  costs little, and saves lots.  (Stupid ^$%*&#@$%& users!)
  //
  bits &= u64bitMASK(size);

  if (b1 >= size) {
    _bfr[wd] &= ~( u64bitMASK(size) << (b1-size) );
    _bfr[wd] |= bits << (b1-size);
  } else {
    _bfr[wd] &= ~u64bitMASK(b1);
    _bfr[wd] |= (bits & (u64bitMASK(b1) << (b2))) >> (b2);

    wd++;

    _bfr[wd] &= ~(u64bitMASK(b2) << (64-b2));
    _bfr[wd] |= (bits & (u64bitMASK(b2))) << (64-b2);
  }

  if (wd == BUFFER_SIZE-1) {
    errno = 0;
    fwrite(_bfr, sizeof(u64bit), BUFFER_SIZE-1, _out);
    if (errno) {
      fprintf(stderr, "bitPackedFileWriter::putBits got %s\n", strerror(errno));
      exit(1);
    }
    _bfr[0] = _bfr[BUFFER_SIZE-1];
    _bit -= ((BUFFER_SIZE-1) * 64);
  }
}




bitPackedFileReader::bitPackedFileReader(char *name) {
  int   errno_fopen = 0;
  int   errno_popen = 0;
  char *command = 0L;

  //  Try to open the original name
  //
  errno = 0;
  _pipe = false;
  _in   = fopen(name, "rb");
  errno_fopen = errno;

  //  If that fails, try to open the compressed version
  //
  if (_in == 0L) {
    size_t  l = strlen(name);

    command = new char [l + 64];

    sprintf(command, "%s.bz2", name);
    errno = 0;
    _in   = fopen(command, "rb");
    errno_popen = errno;

    //  If the file exists, open a pipe.
    //
    if (_in) {
      fclose(_in);

      sprintf(command, "bzip2 -dc %s.bz2", name);

      errno = 0;
      _pipe = true;
      _in   = popen(command, "r");
      errno_popen = errno;
    }
  }

  //  Make sure we aren't at EOF
  //
  if (_in == 0L) {
    fprintf(stderr, "bitPackedFileReader::bitPackedFileReader()-- Couldn't open bitPackedFile.\n");
    fprintf(stderr, "bitPackedFileReader::bitPackedFileReader()-- %s: %s\n", name, strerror(errno_fopen));
    fprintf(stderr, "bitPackedFileReader::bitPackedFileReader()-- %s: %s\n", command, strerror(errno_popen));
    exit(1);
  } else {
    int c = getc(_in);

    if (feof(_in)) {
      if (command)
        fprintf(stderr, "bitPackedFileReader::bitPackedFileReader()-- Empty bitPackedFile '%s.bz2'\n", name);
      else
        fprintf(stderr, "bitPackedFileReader::bitPackedFileReader()-- Empty bitPackedFile '%s'\n", name);
      exit(1);
    }

    ungetc(c, _in);
  }

  delete [] command;

  _bfr = new u64bit [BUFFER_SIZE];
  _bit = u64bitZERO;

  //  Read the initial buffer
  //
  errno = 0;
  size_t bytesread = fread(_bfr, sizeof(u64bit), BUFFER_SIZE, _in);
  if (errno) {
    fprintf(stderr, "bitPackedFileReader::bitPackedFileReader got %s\n", strerror(errno));
    exit(1);
  }

  //  Clear any bytes that we didn't read (EOF)
  while (bytesread < BUFFER_SIZE)
    _bfr[bytesread++] = 0;
}

bitPackedFileReader::~bitPackedFileReader() {
  delete [] _bfr;

  if (_pipe)
    pclose(_in);
  else
    fclose(_in);
}

//  Essentially getDecodedValue, but will update _WORD and _bit.
//
u64bit
bitPackedFileReader::getBits(u32bit size) {
  u64bit wd  = (_bit >> 6) & 0x0000cfffffffffffllu;
  u64bit bt  = (_bit     ) & 0x000000000000003fllu;
  u64bit b1  = 64 - bt;
  u64bit b2  = size - b1;  //  Only used if siz > b1
  u64bit ret = 0;

  assert(size < 65);
  assert(size > 0);

  _bit += size;

  //  If we are in the last word of the buffer, we should refill.
  //  This is done by copying the last word to the first word,
  //  then reading BUFFER_SIZE-1 new words.
  //
  if (wd == BUFFER_SIZE-1) {
    _bfr[0] = _bfr[BUFFER_SIZE-1];
    _bit -= ((BUFFER_SIZE-1) * 64);
    wd = 0;
    errno = 0;
    size_t bytesread = fread(_bfr+1, sizeof(u64bit), BUFFER_SIZE-1, _in);
    if (errno) {
      fprintf(stderr, "bitPackedFileReader::getBits got %s\n", strerror(errno));
      exit(1);
    }

    //  Clear any bytes that we didn't read (EOF)
    bytesread++;

    while (bytesread < BUFFER_SIZE)
      _bfr[bytesread++] = 0;
  }

  if (b1 >= size) {
    ret = _bfr[wd] >> (b1 - size);
  } else {
    ret  = (_bfr[wd] & u64bitMASK(b1)) << b2;
    ret |= (_bfr[wd+1] >> (64 - b2)) & u64bitMASK(b2);
  }

  ret &= u64bitMASK(size);

  return(ret);
}

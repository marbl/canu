#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#include "bitPackedFile.H"
#include "bit-packing.H"

#define BUFFER_SIZE   (16)   //(131072)

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



//  If the buffer doesn't have 128 bits free (2 words currently) flush
//  it.
//
//  128 was chosen because the fibonacci encoded numbers use up to
//  90-some bits.
//
void
bitPackedFileWriter::flush(void) {

  if ((_bit >> 6) >= (BUFFER_SIZE - 2)) {
    errno = 0;
    fwrite(_bfr, sizeof(u64bit), BUFFER_SIZE-2, _out);
    if (errno) {
      fprintf(stderr, "bitPackedFileWriter::putBits got %s\n", strerror(errno));
      exit(1);
    }

    //  copy the last two words -- _bit could be referring to
    //  something in either of the last two words -- and subtract out
    //  the bits that we just wrote.
    //
    _bfr[0] = _bfr[BUFFER_SIZE-2];
    _bfr[1] = _bfr[BUFFER_SIZE-1];
    _bit   -= ((BUFFER_SIZE-2) * 64);
  }
}

void
bitPackedFileWriter::putBits(u64bit bits, u32bit siz) {
  flush();
  setDecodedValue(_bfr, _bit, siz, bits);
  _bit += siz;
}

void
bitPackedFileWriter::putNumber(u64bit val) {
  flush();
  u64bit siz = 0;
  setFibonacciEncodedNumber(_bfr, _bit, &siz, val);
  _bit += siz;
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



void
bitPackedFileReader::fill(void) {

  if ((_bit >> 6) >= (BUFFER_SIZE-2)) {
    _bfr[0] = _bfr[BUFFER_SIZE-2];
    _bfr[1] = _bfr[BUFFER_SIZE-1];
    _bit   -= ((BUFFER_SIZE-2) * 64);

    errno = 0;
    size_t bytesread = fread(_bfr+2, sizeof(u64bit), BUFFER_SIZE-2, _in);
    if (errno) {
      fprintf(stderr, "bitPackedFileReader::getBits got %s\n", strerror(errno));
      exit(1);
    }

    //  Clear any bytes that we didn't read (supposedly, because we
    //  hit EOF).  The +2 is because we started with a buffer with 2
    //  words in it.
    //
    for (bytesread += 2; bytesread < BUFFER_SIZE; bytesread++)
      _bfr[bytesread] = u64bitZERO;
  }
}


//  Essentially getDecodedValue, but will update _WORD and _bit.
//
u64bit
bitPackedFileReader::getBits(u32bit siz) {
  fill();
  u64bit ret = getDecodedValue(_bfr, _bit, siz);
  _bit += siz;
  return(ret);
}


u64bit
bitPackedFileReader::getNumber(void) {
  fill();
  u64bit siz = 0;
  u64bit ret = getFibonacciEncodedNumber(_bfr, _bit, &siz);
  _bit += siz;
  return(ret);
}


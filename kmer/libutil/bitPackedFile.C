#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#include "bitPackedFile.H"
#include "bit-packing.H"

//  This is the size (in 64-bit words) of the buffer the bit packed
//  file uses.  Too small and we spend lots of time doing I/O
//  requests, too big and we use excessive amounts of memory.
//
#define BUFFER_SIZE   (1048576 / 8)

//  Define this to build a short test executable.  A nice test size is
//  50000000 5000.
//
//#define TEST_BITPACKEDFILE


bitPackedFileWriter::bitPackedFileWriter(char const *name) {
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
//  it.  128 was chosen because the fibonacci encoded numbers use up
//  to 90-some bits.
//
void
bitPackedFileWriter::flush(void) {

  if ((_bit >> 6) >= (BUFFER_SIZE - 2)) {
    errno = 0;
    fwrite(_bfr, sizeof(u64bit), BUFFER_SIZE-2, _out);
    if (errno) {
      fprintf(stderr, "bitPackedFileWriter::flush() got %s\n", strerror(errno));
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









bitPackedFileReader::bitPackedFileReader(char const *name) {
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

  fillBufferFromDisk(0);
}

bitPackedFileReader::~bitPackedFileReader() {
  delete [] _bfr;

  if (_pipe)
    pclose(_in);
  else
    fclose(_in);
}


//  Fills the buffer from the disk file.  The parameter is how much of
//  the buffer (at the start) should remain.  '0' results in a
//  complete fill, while '10' would save the first 10 64-bit words.
//
void
bitPackedFileReader::fillBufferFromDisk(u32bit bufStart) {

#if 0
  memset(_bfr + bufStart, 0, sizeof(u64bit) * (BUFFER_SIZE-bufStart));
#endif

  errno = 0;
  size_t wordsread = fread(_bfr + bufStart, sizeof(u64bit), BUFFER_SIZE - bufStart, _in);
  if (errno) {
    fprintf(stderr, "bitPackedFileReader::bitPackedFileReader got %s\n", strerror(errno));
    exit(1);
  }

  //  Clear any words that we didn't read (supposedly, because we
  //  hit EOF).
  //
  for (wordsread += bufStart; wordsread < BUFFER_SIZE; wordsread++)
    _bfr[wordsread] = u64bitZERO;
}


void
bitPackedFileReader::fill(void) {

  if ((_bit >> 6) >= (BUFFER_SIZE-2)) {
    _bfr[0] = _bfr[BUFFER_SIZE-2];
    _bfr[1] = _bfr[BUFFER_SIZE-1];
    _bit   -= ((BUFFER_SIZE-2) * 64);

    fillBufferFromDisk(2);
  }
}


//  Seeks to bitposition pos in the file.
//
void
bitPackedFileReader::seek(u64bit pos) {

  _bit = pos & 0x003f;

  errno = 0;
  fseeko(_in, (pos >> 6) << 3, SEEK_SET);
  if (errno) {
    fprintf(stderr, "bitPackedFileReader::seek() failed: %s\n", strerror(errno));
    exit(1);
  }

  fillBufferFromDisk(0);
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




#ifdef TEST_BITPACKEDFILE
#include <unistd.h>
#include <time.h>
#include <math.h>

int
main(int argc, char **argv) {

  u32bit   testSize = 5000000;
  u32bit   testIter = 200;

  if (argc == 1) {
    fprintf(stderr, "usage: %s testSize testIter\n", argv[0]);
    fprintf(stderr, "  This will perform various tests on the bitPackedFile* classes,\n");
    fprintf(stderr, "  returning 0 if OK and 1 if error.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  testSize -- the number of words to use in a write then read test\n");
    fprintf(stderr, "  testIter -- the number of random access tests to do\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "I'll assume reasonable values and continue\n");
    fprintf(stderr, "  testSize = "u32bitFMT"\n", testSize);
    fprintf(stderr, "  testIter = "u32bitFMT"\n", testIter);
  } else {
    testSize = strtou32bit(argv[1], 0L);
    testIter = strtou32bit(argv[2], 0L);
  }

  u32bit    i;
  u32bit   *siz = new u32bit [testSize];
  u64bit   *val = new u64bit [testSize];
  u32bit    errs = 0;

  bitPackedFileWriter *W;
  bitPackedFileReader *R;

  srand48(time(NULL));

  //  Generate a list of random 64-bit numbers, remember the number and the size
  //
  fprintf(stderr, "Initializing\n");
  for (i=0; i<testSize; i++) {
    siz[i] =  (u32bit)floor(drand48() * 64) + 1;
    val[i] = ((u64bit)lrand48()) << 32 | ((u64bit)lrand48());
    val[i] &= u64bitMASK(siz[i]);

    if ((i & 0xfffff) == 0) {
      fprintf(stderr, "Init  "u32bitFMT"            \r", i);
      fflush(stderr);
    }
  }

  fprintf(stderr, "Testing bitPackedFile\n");

  //  Write those numbers to a bitPackedFile, both binary encoded and
  //  fibonacci encoded.
  //
  W = new bitPackedFileWriter("bittest.junk");
  for (i=0; i<testSize; i++) {
    W->putBits(val[i], siz[i]);
    W->putNumber(val[i]);
    W->putNumber(val[i]);

    if ((i & 0xfffff) == 0) {
      fprintf(stderr, "Write "u32bitFMT"            \r", i);
      fflush(stderr);
    }
  }
  delete W;

  //  Open the file and check what we just wrote.
  //
  R = new bitPackedFileReader("bittest.junk");
  for (i=0; i<testSize; i++) {
    u64bit v;

    v = R->getBits(siz[i]);
    if (v != val[i]) {
      fprintf(stderr, u32bitFMT"] ERROR in getBits()   -- retrieved "u64bitHEX" != expected "u64bitHEX" ("u32bitFMT" bits).\n", i, v, val[i], siz[i]);
      errs++;
    }

    v = R->getNumber();
    v = R->getNumber();
    if (v != val[i]) {
      fprintf(stderr, u32bitFMT"] ERROR in getNumber() -- retrieved "u64bitHEX" != expected "u64bitHEX".\n", i, v, val[i]);
      errs++;
    }

    if ((i & 0xffff) == 0) {
      fprintf(stderr, "Read  "u32bitFMT"            \r", i);
      fflush(stderr);
    }
  }
  delete R;

  fprintf(stderr, "Testing bitPackedFile::seek()\n");

  //  Create a new bitpacked file, writing just numbers as binary encoded.
  //
  W = new bitPackedFileWriter("bittest.junk");
  for (i=0; i<testSize; i++) {
    W->putBits(val[i], siz[i]);

    if ((i & 0xffff) == 0) {
      fprintf(stderr, "Write "u32bitFMT"            \r", i);
      fflush(stderr);
    }
  }
  delete W;

  //  Do several seek tests.  Seek to a random element, and read it.
  //
  R = new bitPackedFileReader("bittest.junk");
  for (i=0; i<testIter; i++) {
    u32bit idx = (u32bit)lrand48() % testSize;
    u64bit pos = 0;

    for (u32bit j=0; j<idx; j++)
      pos += siz[j];

    R->seek(pos);
    u64bit r = R->getBits(siz[idx]);

    if (r != val[idx]) {
      fprintf(stderr, u32bitFMT"] ERROR in seek()/getBits()   -- retrieved "u64bitHEX" != expected "u64bitHEX" ("u32bitFMT" bits).\n", i, r, val[i], siz[i]);
      errs++;
    }

    if ((i & 0xf) == 0) {
      fprintf(stderr, "Read  "u32bitFMT"            \r", i);
      fflush(stderr);
    }
  }
  delete R;

  unlink("bittest.junk");

  if (errs > 0) {
    fprintf(stderr, "There are "u32bitFMT" errors.\n", errs);
    exit(1);
  } else {
    exit(0);
  }
}


#endif  //  TEST_BITPACKEDFILE

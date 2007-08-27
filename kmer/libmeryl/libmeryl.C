#include "libmeryl.H"


merylStreamReader::merylStreamReader(const char *fn, u32bit ms) {

  if (fn == 0L) {
    fprintf(stderr, "ERROR - no counted database file specified.\n");
    exit(1);
  }

  //  Open the files
  //
  char *inpath = new char [strlen(fn) + 8];

  sprintf(inpath, "%s.mcidx", fn);
  _IDX = new bitPackedFile(inpath);

  sprintf(inpath, "%s.mcdat", fn);
  _DAT = new bitPackedFile(inpath);

  delete [] inpath;

  //  Verify that they are what they should be, and read in the header
  //
  char  Imagic[16];
  char  Dmagic[16];
  bool  fail = false;
  for (u32bit i=0; i<16; i++) {
    Imagic[i] = _IDX->getBits(8);
    Dmagic[i] = _DAT->getBits(8);
  }
  if (strncmp(Imagic, "merylStreamIv01\n", 16) != 0) {
    if (strncmp(Imagic, "merylStreamIvXX\n", 16) == 0) {
      fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcidx is an INCOMPLETE merylStream index file!\n", fn);
      fail = true;
    } else {
      fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcidx isn't a merylStream index file!\n", fn);
      fail = true;
    }
  }
  if (strncmp(Dmagic, "merylStreamDv01\n", 16) != 0) {
    if (strncmp(Dmagic, "merylStreamDvXX\n", 16) == 0) {
      fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcidx is an INCOMPLETE merylStream data file!\n", fn);
      fail = true;
    } else {
      fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcdat isn't a merylStream data file!\n", fn);
      fail = true;
    }
  }
  if (fail)
    exit(1);

  _merSizeInBits  = _IDX->getBits(32) << 1;
  _prefixSize     = _IDX->getBits(32);
  _merDataSize    = _merSizeInBits - _prefixSize;

  _numUnique      = _IDX->getBits(64);
  _numDistinct    = _IDX->getBits(64);
  _numTotal       = _IDX->getBits(64);

  _thisBucket     = u64bitZERO;
  _thisBucketSize = _IDX->getNumber();
  _numBuckets     = u64bitONE << _prefixSize;

  _thisMer.setMerSize(_merSizeInBits >> 1);
  _thisMer.clear();
  _thisMerCount   = u64bitZERO;

  _validMer       = true;

#ifdef SHOW_VARIABLES
  fprintf(stderr, "_merSizeInBits  = "u32bitFMT"\n", _merSizeInBits);
  fprintf(stderr, "_prefixSize     = "u32bitFMT"\n", _prefixSize);
  fprintf(stderr, "_merDataSize    = "u32bitFMT"\n", _merDataSize);
  fprintf(stderr, "_numUnique      = "u64bitFMT"\n", _numUnique);
  fprintf(stderr, "_numDistinct    = "u64bitFMT"\n", _numDistinct);
  fprintf(stderr, "_numTotal       = "u64bitFMT"\n", _numTotal);
  fprintf(stderr, "_thisBucket     = "u64bitFMT"\n", _thisBucket);
  fprintf(stderr, "_thisBucketSize = "u64bitFMT"\n", _thisBucketSize);
  fprintf(stderr, "_thisMerCount   = "u64bitFMT"\n", _thisMerCount);
#endif

  if ((ms > 0) && (_merSizeInBits >> 1 != ms)) {
    fprintf(stderr, "merylStreamReader()-- ERROR: User requested mersize "u32bitFMT" but '%s' is mersize "u32bitFMT"\n",
            ms, fn, _merSizeInBits >> 1);
    exit(1);
  }
}


merylStreamReader::~merylStreamReader() {
  delete _IDX;
  delete _DAT;
}



merylStreamWriter::merylStreamWriter(const char *fn,
                                     u32bit merSize,
                                     u32bit prefixSize) {

  char *outpath = new char [strlen(fn) + 17];

  sprintf(outpath, "%s.mcidx", fn);
  _IDX = new bitPackedFile(outpath, 0, true);

  sprintf(outpath, "%s.mcdat", fn);
  _DAT = new bitPackedFile(outpath, 0, true);

  delete [] outpath;

  //  Save really important stuff
  //
  _merSizeInBits  = merSize * 2;
  _prefixSize     = prefixSize;
  _merDataSize    = _merSizeInBits - _prefixSize;

  _thisBucket     = u64bitZERO;
  _thisBucketSize = u64bitZERO;
  _numBuckets     = u64bitONE << _prefixSize;

  _numUnique      = u64bitZERO;
  _numDistinct    = u64bitZERO;
  _numTotal       = u64bitZERO;

  _thisMer.setMerSize(_merSizeInBits >> 1);
  _thisMer.clear();
  _thisMerCount   = u64bitZERO;

  //  Write the headers
  //    16 bytes for a magic cookie
  //        The IDX gets "merylStreamIv01\n"
  //        The DAT gets "merylStreamDv01\n"
  //     4 bytes for merSize
  //     4 bytes for prefixSize
  //     8 bytes for numUnique
  //     8 bytes for numDistinct
  //     8 bytes for numTotal

  _IDX->putBits('m',  8);
  _IDX->putBits('e',  8);
  _IDX->putBits('r',  8);
  _IDX->putBits('y',  8);
  _IDX->putBits('l',  8);
  _IDX->putBits('S',  8);
  _IDX->putBits('t',  8);
  _IDX->putBits('r',  8);
  _IDX->putBits('e',  8);
  _IDX->putBits('a',  8);
  _IDX->putBits('m',  8);
  _IDX->putBits('I',  8);
  _IDX->putBits('v',  8);
  _IDX->putBits('X',  8);
  _IDX->putBits('X',  8);
  _IDX->putBits('\n', 8);

  _IDX->putBits(_merSizeInBits >> 1, 32);
  _IDX->putBits(_prefixSize,  32);
  _IDX->putBits(_numUnique,   64);
  _IDX->putBits(_numDistinct, 64);
  _IDX->putBits(_numTotal,    64);

  _DAT->putBits('m',  8);
  _DAT->putBits('e',  8);
  _DAT->putBits('r',  8);
  _DAT->putBits('y',  8);
  _DAT->putBits('l',  8);
  _DAT->putBits('S',  8);
  _DAT->putBits('t',  8);
  _DAT->putBits('r',  8);
  _DAT->putBits('e',  8);
  _DAT->putBits('a',  8);
  _DAT->putBits('m',  8);
  _DAT->putBits('D',  8);
  _DAT->putBits('v',  8);
  _DAT->putBits('X',  8);
  _DAT->putBits('X',  8);
  _DAT->putBits('\n', 8);
}


merylStreamWriter::~merylStreamWriter() {

  //  Write out the last mer
  //
  writeMer();

  //  Finish writing the buckets.
  //
  while (_thisBucket < _numBuckets) {
    _IDX->putNumber(_thisBucketSize);
    _thisBucketSize = 0;
    _thisBucket++;
  }

  _IDX->putNumber(_thisBucketSize);
  _IDX->putNumber(_thisBucketSize);

  //  Seek back to the start and rewrite the magic numbers
  //
  _IDX->seek(0);
  _IDX->putBits('m',  8);
  _IDX->putBits('e',  8);
  _IDX->putBits('r',  8);
  _IDX->putBits('y',  8);
  _IDX->putBits('l',  8);
  _IDX->putBits('S',  8);
  _IDX->putBits('t',  8);
  _IDX->putBits('r',  8);
  _IDX->putBits('e',  8);
  _IDX->putBits('a',  8);
  _IDX->putBits('m',  8);
  _IDX->putBits('I',  8);
  _IDX->putBits('v',  8);
  _IDX->putBits('0',  8);
  _IDX->putBits('1',  8);
  _IDX->putBits('\n', 8);

  _IDX->putBits(_merSizeInBits >> 1,     32);
  _IDX->putBits(_prefixSize,  32);
  _IDX->putBits(_numUnique,   64);
  _IDX->putBits(_numDistinct, 64);
  _IDX->putBits(_numTotal,    64);

  delete _IDX;

  _DAT->seek(0);
  _DAT->putBits('m',  8);
  _DAT->putBits('e',  8);
  _DAT->putBits('r',  8);
  _DAT->putBits('y',  8);
  _DAT->putBits('l',  8);
  _DAT->putBits('S',  8);
  _DAT->putBits('t',  8);
  _DAT->putBits('r',  8);
  _DAT->putBits('e',  8);
  _DAT->putBits('a',  8);
  _DAT->putBits('m',  8);
  _DAT->putBits('D',  8);
  _DAT->putBits('v',  8);
  _DAT->putBits('0',  8);
  _DAT->putBits('1',  8);
  _DAT->putBits('\n', 8);
  delete _DAT;
}


void
merylStreamWriter::writeMer(void) {

  if (_thisMerCount == 0)
    return;

#if 0
  char  str[1025];
  fprintf(stderr, "Saving mer '%s' with count "u32bitFMT"\n", _thisMer.merToString(str), _thisMerCount);
#endif

  _numTotal += _thisMerCount;
  _numDistinct++;

  if (_thisMerCount == 1) {
    _thisMer.writeToBitPackedFile(_DAT, _merDataSize);
    _DAT->putBits(u64bitZERO, 1);
    _thisBucketSize++;
    _numUnique++;
  } else if (_thisMerCount > 1) {
    _thisMer.writeToBitPackedFile(_DAT, _merDataSize);
    _DAT->putBits(u64bitONE, 1);
    _DAT->putNumber(_thisMerCount - 2);
    _thisBucketSize++;
  }
}


void
merylStreamWriter::addMer(kMer &mer, u32bit count) {
  u64bit  val;

  //  Fail if we see a smaller mer than last time.
  //
  if (mer < _thisMer) {
    char str[1024];
    fprintf(stderr, "merylStreamWriter::addMer()-- ERROR: your mer stream isn't sorted increasingly!\n");
    fprintf(stderr, "merylStreamWriter::addMer()-- last: %s\n", _thisMer.merToString(str));
    fprintf(stderr, "merylStreamWriter::addMer()-- this: %s\n", mer.merToString(str));
    exit(1);
  }

  //  If the new mer is the same as the last one just increase the
  //  count.
  //
  if (mer == _thisMer) {
#if 0
    char str[1024], sts[1024];
    fprintf(stderr, "add one mer=%s thisMer=%s\n", mer.merToString(str), _thisMer.merToString(sts));
#endif
    _thisMerCount += count;
    return;
  }

  //  Write thisMer to disk.  If the count is zero, we don't write
  //  anything.  The count is zero for the first mer (all A) unless we
  //  add that mer, and if the silly user gives us a mer with zero
  //  count.
  //
  writeMer();

  //  If the new mer is in a different bucket from the last mer, write
  //  out some bucket counts.  We need a while loop (opposed to just
  //  writing one bucket) because we aren't guaranteed that the mers
  //  are in adjacent buckets.
  //
  val = mer.startOfMer(_prefixSize);

  while (_thisBucket < val) {
    //fprintf(stderr, "bucket "u64bitFMT" with size "u64bitFMT"\n", _thisBucket, _thisBucketSize);
    _IDX->putNumber(_thisBucketSize);
    _thisBucketSize = 0;
    _thisBucket++;
  }

  //  Remember the new mer for the next time
  //
  _thisMer      = mer;
  _thisMerCount = count;
}

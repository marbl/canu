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
  _IDX = new bitPackedFileReader(inpath);

  sprintf(inpath, "%s.mcdat", fn);
  _DAT = new bitPackedFileReader(inpath);

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
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcidx isn't a merylStream index file!\n", fn);
    fail = true;
  }
  if (strncmp(Dmagic, "merylStreamDv01\n", 16) != 0) {
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcdat isn't a merylStream data file!\n", fn);
    fail = true;
  }
  if (fail)
    exit(1);


  _merSizeInBits  = _IDX->getBits(32) << 1;
  _prefixSize     = _IDX->getBits(32);
  _prefixMask     = u64bitMASK(_prefixSize);
  _merDataSize    = _merSizeInBits - _prefixSize;
  _merDataMask    = u64bitMASK(_merDataSize);

  _numUnique      = _IDX->getBits(64);
  _numDistinct    = _IDX->getBits(64);
  _numTotal       = _IDX->getBits(64);

  _thisBucket     = u64bitZERO;
  _thisBucketSize = _IDX->getNumber();

  _thisMer        = u64bitZERO;
  _thisMerCount   = u64bitZERO;

  _validMer       = true;

#if 0
  fprintf(stderr, "%s\n", fn);
  fprintf(stderr, "_merSizeInBits  = "u32bitFMT"\n", _merSizeInBits);
  fprintf(stderr, "_prefixSize     = "u32bitFMT"\n", _prefixSize);
  fprintf(stderr, "_prefixMask     = "u64bitHEX"\n", _prefixMask);
  fprintf(stderr, "_merDataSize    = "u32bitFMT"\n", _merDataSize);
  fprintf(stderr, "_merDataMask    = "u64bitHEX"\n", _merDataMask);
  fprintf(stderr, "_numUnique      = "u64bitFMT"\n", _numUnique);
  fprintf(stderr, "_numDistinct    = "u64bitFMT"\n", _numDistinct);
  fprintf(stderr, "_numTotal       = "u64bitFMT"\n", _numTotal);
  fprintf(stderr, "_thisBucket     = "u64bitFMT"\n", _thisBucket);
  fprintf(stderr, "_thisBucketSize = "u64bitFMT"\n", _thisBucketSize);
  fprintf(stderr, "_thisMer        = "u64bitFMT"\n", _thisMer);
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
                                     u32bit prefixSize,
                                     u64bit numUnique,
                                     u64bit numDistinct,
                                     u64bit numTotal) {


  //  Open the files
  //
  char *outpath = new char [strlen(fn) + 17];

  sprintf(outpath, "%s.mcidx", fn);
  _IDX = new bitPackedFileWriter(outpath);

  sprintf(outpath, "%s.mcdat", fn);
  _DAT = new bitPackedFileWriter(outpath);

  delete [] outpath;


  //  Save really important stuff
  //
  _merSizeInBits  = merSize * 2;
  _prefixSize     = prefixSize;
  _prefixMask     = u64bitMASK(_prefixSize);
  _merDataSize    = _merSizeInBits - _prefixSize;
  _merDataMask    = u64bitMASK(_merDataSize);

  _thisBucket     = u64bitZERO;
  _thisBucketSize = u64bitZERO;

  _thisMer        = u64bitZERO;
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
  //

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

  _IDX->putBits(merSize,     32);
  _IDX->putBits(prefixSize,  32);
  _IDX->putBits(numUnique,   64);
  _IDX->putBits(numDistinct, 64);
  _IDX->putBits(numTotal,    64);

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
}


merylStreamWriter::~merylStreamWriter() {
  u64bit val;

  //  Write out the last mer
  //
  if (_thisMerCount == 1) {
    val  = _thisMer;
    val &= _merDataMask;
    val |= u64bitONE << _merDataSize;
    _DAT->putBits(val, _merDataSize + 1);
    _thisBucketSize++;
  } else if (_thisMerCount > 1) {
    val  = _thisMer;
    val &= _merDataMask;
    _DAT->putBits(val, _merDataSize + 1);
    _DAT->putNumber(_thisMerCount - 2);
    _thisBucketSize++;
  }


  //fprintf(stderr, "writing buckets from "u64bitFMT" to "u64bitFMT"\n", _thisBucket, _prefixMask);

  //  Finish writing the buckets.
  //
  while (_thisBucket < _prefixMask) {
    _IDX->putNumber(_thisBucketSize);
    _thisBucketSize = 0;
    _thisBucket++;
  }

  _IDX->putNumber(_thisBucketSize);
  _IDX->putNumber(_thisBucketSize);

  delete _IDX;
  delete _DAT;
}


void
merylStreamWriter::addMer(u64bit mer, u32bit count) {
  u64bit  val;

#if 0
  char merstring[33] = { 0 };
  for (u32bit i=0; i<(_merSizeInBits>>1); i++)
    merstring[(_merSizeInBits>>1)-i-1] = decompressSymbol[(mer >> (2*i)) & 0x03];
  merstring[(_merSizeInBits>>1)] = 0;
  fprintf(stderr, "write mer "u64bitHEX" -- %s with count "u32bitFMT"\n", mer, merstring, count);
#endif

  //  Fail if we see a smaller mer than last time.
  //
  if (mer < _thisMer) {
    fprintf(stderr, "merylStreamWriter::addMer()-- ERROR: your mer stream isn't sorted increasingly!\n");
    exit(1);
  }

  //  If the new mer the same as the last one just increase the count.
  //
  if (mer == _thisMer) {
    _thisMerCount += count;
    return;
  }

  //  Write thisMer to disk.  If the count is zero, we don't write
  //  anything.  The count is zero for the first mer (all A) unless we
  //  add that mer, and if the silly user gives us a mer with zero
  //  count.
  //
  if (_thisMerCount == 1) {
    val  = _thisMer & _merDataMask;
    val |= u64bitONE << _merDataSize;
    _DAT->putBits(val, _merDataSize + 1);
    _thisBucketSize++;
  } else if (_thisMerCount > 1) {
    val  = _thisMer & _merDataMask;
    _DAT->putBits(val, _merDataSize + 1);
    _DAT->putNumber(_thisMerCount - 2);
    _thisBucketSize++;
  }

  //  If the new mer is in a different bucket from the last mer, write
  //  out some bucket counts.  We need a while loop (opposed to just
  //  writing one bucket) because we aren't guaranteed that the mers
  //  are in adjacent buckets.
  //
  val = (mer >> _merDataSize) & _prefixMask;

  while (_thisBucket < val) {
    _IDX->putNumber(_thisBucketSize);
    _thisBucketSize = 0;
    _thisBucket++;
  }

  //  Remember the new mer for the next time
  //
  _thisMer      = mer;
  _thisMerCount = count;
}

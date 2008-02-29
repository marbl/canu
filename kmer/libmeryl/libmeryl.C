#include "libmeryl.H"

//                      0123456789012345
static char *ImagicV = "merylStreamIv02\n";
static char *ImagicX = "merylStreamIvXX\n";
static char *DmagicV = "merylStreamDv02\n";
static char *DmagicX = "merylStreamDvXX\n";


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
  char    Imagic[16] = {0};
  char    Dmagic[16] = {0};
  bool    fail       = false;

  for (u32bit i=0; i<16; i++) {
    Imagic[i] = _IDX->getBits(8);
    Dmagic[i] = _DAT->getBits(8);
  }
  if (strncmp(Imagic, ImagicX, 16) == 0) {
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcidx is an INCOMPLETE merylStream index file!\n", fn);
    fail = true;
  }
  if (strncmp(Imagic, ImagicX, 13) != 0) {
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcidx is not a merylStream index file!\n", fn);
    fail = true;
  }
  if (strncmp(Dmagic, DmagicX, 16) == 0) {
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcdat is an INCOMPLETE merylStream data file!\n", fn);
    fail = true;
  }
  if (strncmp(Dmagic, DmagicX, 13) != 0) {
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcdat is not a merylStream data file!\n", fn);
    fail = true;
  }
  if ((Imagic[13] != Dmagic[13]) ||
      (Imagic[14] != Dmagic[14])) {
    fprintf(stderr, "merylStreamReader()-- ERROR: %s.mcidx and %s.mcdat are different versions!\n", fn, fn);
    fail = true;
  }
  if (fail)
    exit(1);

  u32bit version = atoi(Imagic + 13);

  _merSizeInBits  = _IDX->getBits(32) << 1;
  _merCompression = _IDX->getBits(32);
  _prefixSize     = _IDX->getBits(32);
  _merDataSize    = _merSizeInBits - _prefixSize;

  _numUnique      = _IDX->getBits(64);
  _numDistinct    = _IDX->getBits(64);
  _numTotal       = _IDX->getBits(64);

  _histogramHuge     = 0;
  _histogramLen      = 0;
  _histogramMaxValue = 0;
  _histogram         = 0L;

  if ((Imagic[13] == '0') && (Imagic[14] == '2')) {
    _histogramHuge     = _IDX->getBits(64);
    _histogramLen      = _IDX->getBits(64);
    _histogramMaxValue = _IDX->getBits(64);
    _histogram         = new u64bit [_histogramLen];

    for (u32bit i=0; i<_histogramLen; i++)
      _histogram[i] = _IDX->getBits(64);
  }

  _thisBucket     = u64bitZERO;
  _thisBucketSize = _IDX->getNumber();
  _numBuckets     = u64bitONE << _prefixSize;

  _thisMer.setMerSize(_merSizeInBits >> 1);
  _thisMer.clear();
  _thisMerCount   = u64bitZERO;

  _validMer       = true;

#ifdef SHOW_VARIABLES
  fprintf(stderr, "_merSizeInBits  = "u32bitFMT"\n", _merSizeInBits);
  fprintf(stderr, "_merCompression = "u32bitFMT"\n", _merCompression);
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
  delete [] _histogram;
}



merylStreamWriter::merylStreamWriter(const char *fn,
                                     u32bit merSize,
                                     u32bit merComp,
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
  _merCompression = merComp;
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

  for (u32bit i=0; i<16; i++)
    _IDX->putBits(ImagicX[i], 8);

  _IDX->putBits(_merSizeInBits >> 1, 32);
  _IDX->putBits(_merCompression, 32);
  _IDX->putBits(_prefixSize,  32);
  _IDX->putBits(_numUnique,   64);
  _IDX->putBits(_numDistinct, 64);
  _IDX->putBits(_numTotal,    64);

  _histogramHuge     = 0;
  _histogramLen      = LIBMERYL_HISTOGRAM_MAX;
  _histogramMaxValue = 0;
  _histogram         = new u64bit [_histogramLen];

  for (u32bit i=0; i<_histogramLen; i++)
    _histogram[i] = 0;

  _IDX->putBits(_histogramHuge, 64);
  _IDX->putBits(_histogramLen, 64);
  _IDX->putBits(_histogramMaxValue, 64);
  for (u32bit i=0; i<_histogramLen; i++)
    _IDX->putBits(_histogram[i], 64);

  for (u32bit i=0; i<16; i++)
    _DAT->putBits(DmagicX[i], 8);
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
  _DAT->seek(0);

  for (u32bit i=0; i<16; i++)
    _IDX->putBits(ImagicV[i], 8);

  _IDX->putBits(_merSizeInBits >> 1, 32);
  _IDX->putBits(_merCompression, 32);
  _IDX->putBits(_prefixSize,  32);
  _IDX->putBits(_numUnique,   64);
  _IDX->putBits(_numDistinct, 64);
  _IDX->putBits(_numTotal,    64);

  _IDX->putBits(_histogramHuge, 64);
  _IDX->putBits(_histogramLen, 64);
  _IDX->putBits(_histogramMaxValue, 64);
  for (u32bit i=0; i<_histogramLen; i++)
    _IDX->putBits(_histogram[i], 64);

  for (u32bit i=0; i<16; i++)
    _DAT->putBits(DmagicV[i], 8);

  delete _IDX;
  delete _DAT;
  delete [] _histogram;
}


void
merylStreamWriter::writeMer(void) {

  if (_thisMerCount == 0)
    return;

  _numTotal += _thisMerCount;
  _numDistinct++;

  if (_thisMerCount < LIBMERYL_HISTOGRAM_MAX)
    _histogram[_thisMerCount]++;
  else
    _histogramHuge++;
  if (_histogramMaxValue < _thisMerCount)
    _histogramMaxValue = _thisMerCount;

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
    _IDX->putNumber(_thisBucketSize);
    _thisBucketSize = 0;
    _thisBucket++;
  }

  //  Remember the new mer for the next time
  //
  _thisMer      = mer;
  _thisMerCount = count;
}

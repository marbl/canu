#include "libmeryl.H"


merStreamFromMeryl::merStreamFromMeryl(const char *prefix, u32bit merSize) {

  if (prefix == 0L) {
    fprintf(stderr, "ERROR - no counted database file specified.\n");
    exit(1);
  }

  //  Open the counted sequence files
  //
  char *inpath = new char [strlen(prefix) + 8];

  sprintf(inpath, "%s.mcidx", prefix);
  _IDX = new bitPackedFileReader(inpath);

  sprintf(inpath, "%s.mcdat", prefix);
  _DAT = new bitPackedFileReader(inpath);

  delete [] inpath;

  _mcd.read(_DAT);

  if ((merSize > 0) && (merSize != _mcd._merSizeInBases)) {
    fprintf(stderr, "ERROR:  User requested mersize=%u but '%s' is mersize=%u\n",
            merSize, prefix, _mcd._merSizeInBases);
    exit(1);
  }

  _B = new mcBucket(_IDX, _DAT, &_mcd);

  _tPos = 0;
  _bPos = 0;

  _mer  = u64bitZERO;
  _cnt  = u64bitZERO;
}


merStreamFromMeryl::~merStreamFromMeryl() {
  delete _DAT;
  delete _IDX;
}

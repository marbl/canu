#include "markFile.H"


markFile::markFile(char *filename) {
  strncpy(_correctMagic, "markDstctMers.v1", 16);

  strcpy(_filename, filename);

  DEBUGlargest=0;

  if (fileExists(filename)) {
    _bpf = new bitPackedFile(filename);
    checkTag();
  } else {
    _bpf = new bitPackedFile(filename);
    writeTag(true);
  }
}


markFile::~markFile() {
  if (_haveWritten)
    writeTag();

  fprintf(stderr, "DEBUGlargest = "u64bitFMT"\n", DEBUGlargest);

  delete _bpf;
}


void
markFile::checkTag(void) {
  char    incorrectmagic[16] = {0};

  _bpf->seek(0);

  for(u32bit i=0; i<16; i++)
    incorrectmagic[i] = _bpf->getBits(8);

  if (strncmp(incorrectmagic, _correctMagic, 16) != 0) {
    fprintf(stderr, "ERROR!  File '%s' doesn't appear to be a markDistinctMers file!\n", _filename);
    fprintf(stderr, "Found: '");
    for(u32bit i=0; i<16; i++)
      fprintf(stderr, "%c", isprint(incorrectmagic[i]) ? incorrectmagic[i] : '.');
    fprintf(stderr, "' [%d", incorrectmagic[0]);
    for(u32bit i=1; i<16; i++)
      fprintf(stderr, " %d", incorrectmagic[i]);
    fprintf(stderr, "]\n");
    exit(1);
  }

  _numberOfMers = _bpf->getBits(64);
  _numberMarked = _bpf->getBits(64);
}


void
markFile::writeTag(bool bogus) {

  _bpf->seek(0);

  if (bogus) {
    //  Write out a magic number, but use version 0 to denote the file
    //  wasn't closed properly.
    //
    for(u32bit i=0; i<15; i++)
      _bpf->putBits(_correctMagic[i], 8);
    _bpf->putBits('0', 8);
  } else {
    //  Write out the correct magic number
    //
    for(u32bit i=0; i<16; i++)
      _bpf->putBits(_correctMagic[i], 8);
  }

  _bpf->putBits(_numberOfMers, 64);
  _bpf->putBits(_numberMarked, 64);
}

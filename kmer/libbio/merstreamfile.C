#include "merstreamfile.H"
#include "fastastream.H"
#include "bitPackedFile.H"
#include "merstream.H"
#include "britime.H"
#include "stat.h"

static
char     magic[16] = { 'm', 'e', 'r', 'S', 't', 'r', 'e', 'a', 'm', '1', ' ', ' ', ' ', ' ', ' ', ' '  };


merStreamFileBuilder::merStreamFileBuilder(u32bit m, const char *i, const char *o) {
  _merSize    = m;
  _inputFile  = i;
  _outputFile = o;
}

merStreamFileBuilder::~merStreamFileBuilder() {
}

void
merStreamFileBuilder::build(bool beVerbose) {
  u64bit             lastMer     = u64bitZERO;
  u64bit             thisMer     = u64bitZERO;
  u64bit             maskMer     = u64bitMASK(_merSize * 2 - 2);
  u64bit             blockSize   = 0;
  u64bit             numMers     = 0;
  u64bit             numBlocks   = 0;
  u64bit             numDefs     = 0;
  u64bit             defLength   = 0;

  //  Open the output files
  //
  char *blocksName = new char [strlen(_outputFile) + 17];
  char *deflinName = new char [strlen(_outputFile) + 17];
  char *streamName = new char [strlen(_outputFile) + 17];
  char *outputName = new char [strlen(_outputFile) + 17];

  sprintf(blocksName, "%s.merstreamfilebuilder.b.tmp", _outputFile);
  sprintf(deflinName, "%s.merstreamfilebuilder.d.tmp", _outputFile);
  sprintf(streamName, "%s.merstreamfilebuilder.s.tmp", _outputFile);
  sprintf(outputName, "%s.merStream",                  _outputFile);

  bitPackedFileWriter *BLOCKS = new bitPackedFileWriter(blocksName);
  bitPackedFileWriter *STREAM = new bitPackedFileWriter(streamName);

  errno = 0;
  FILE  *DEFLIN  = fopen(deflinName, "w");
  if (errno) {
    fprintf(stderr, "merStreamFileBuilder::build()-- failed to open defline storage: %s\n", strerror(errno));
    exit(1);
  }
  u32bit lastDef = ~u32bitZERO;


  merStream          M(_merSize, _inputFile);
  speedCounter       C("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);

  //  We bootstrap the writer by writing out the first n-1 letters of
  //  the first mer.  M.theFMer() is not the full first mer, it's the
  //  first mer without the last base.  This is exactly what we want
  //  to make the first mer in the file look like a consecutive mer.
  //  If we didn't do this, the first mer would look like a break, and
  //  we would write out a block of size 1.
  //
  lastMer = M.theFMer();
  STREAM->putBits(lastMer, _merSize * 2 - 2);

  while (M.nextMer()) {
    thisMer = M.theFMer();

    //  Save the defline if the sequence number has changed.
    //
    if (lastDef != M.theSequenceNumber()) {
      u32bit  l = strlen(M.theDefLine()) + 1;

      fwrite(&l, sizeof(u32bit), 1, DEFLIN);
      fwrite(M.theDefLine(), sizeof(char), l, DEFLIN);

      lastDef = M.theSequenceNumber();

      numDefs++;
      defLength += l;
    }

    if ((thisMer >> 2) == (lastMer & maskMer)) {
      //  Consecutive mers, just emit the new base, and add one to our count
      //
      STREAM->putBits(thisMer & 0x03, 2);
      blockSize++;
    } else {
      //  Must be a mer break.  Write the count, and the new mer.
      //
      //fprintf(stderr, "Block: "u64bitFMT" "u64bitFMT" "u64bitFMT"\n", blockSize, (u64bit)(M.theSequenceNumber()), (u64bit)(M.thePosition()-blockSize));

      STREAM->putBits(thisMer, _merSize * 2);
      BLOCKS->putNumber(blockSize);
      BLOCKS->putNumber(M.theSequenceNumber());
      BLOCKS->putNumber(M.thePosition() - blockSize);
      numBlocks++;
      blockSize = 1;
    }

    numMers++;

    lastMer = thisMer;

    C.tick();
  }

  //  Be sure to write out the last block!
  //
  BLOCKS->putNumber(blockSize);
  BLOCKS->putNumber(M.theSequenceNumber());
  BLOCKS->putNumber(M.thePosition() - blockSize);
  numBlocks++;

  fclose(DEFLIN);
  delete STREAM;
  delete BLOCKS;

  //  Our last step is to merge the three files into one.  This can't
  //  be done as we build, since we want to know how many blocks there
  //  are, and we want to know the blocks and deflines before the mers
  //  themselves.
  //
  //  Fortunately, it's just a simple process of raw data copies, not
  //  parsing the files.
  //
  errno = 0;
  FILE  *fOUT = fopen(outputName, "w");
  if (errno) {
    fprintf(stderr, "merStreamFileBuilder::build()-- failed to open the output file: %s\n", strerror(errno));
    exit(1);
  }

  off_t blockFileSize   = sizeOfFile(blocksName);
  off_t deflineFileSize = sizeOfFile(deflinName);
  off_t streamFileSize  = sizeOfFile(streamName);

  off_t blkStart = 8 * (sizeof(char) * 16 + sizeof(u32bit) * 2 + sizeof(u64bit) * 4 + sizeof(off_t) * 6);
  off_t strStart = blkStart + 8 * blockFileSize;
  off_t defStart = blkStart + 8 * blockFileSize + 8 * streamFileSize;

  //
  //  This block MUST be a multiple of 64 bits long to keep the
  //  bitPackedFiles 64-bit word aligned.
  //

  fwrite(magic,            sizeof(char),   16, fOUT);
  fwrite(&_merSize,        sizeof(u32bit), 1,  fOUT);
  fwrite(&_merSize,        sizeof(u32bit), 1,  fOUT);  //  Padding.
  fwrite(&numMers,         sizeof(u64bit), 1,  fOUT);
  fwrite(&numBlocks,       sizeof(u64bit), 1,  fOUT);
  fwrite(&numDefs,         sizeof(u64bit), 1,  fOUT);
  fwrite(&defLength,       sizeof(u64bit), 1,  fOUT);
  fwrite(&blockFileSize,   sizeof(off_t),  1,  fOUT);
  fwrite(&deflineFileSize, sizeof(off_t),  1,  fOUT);
  fwrite(&streamFileSize,  sizeof(off_t),  1,  fOUT);
  fwrite(&blkStart,        sizeof(off_t),  1,  fOUT);
  fwrite(&defStart,        sizeof(off_t),  1,  fOUT);
  fwrite(&strStart,        sizeof(off_t),  1,  fOUT);

  //  To keep the bitPackedFiles aligned on 64-bit boundaries, we put the deflines last.
  //
  copyFile(blocksName, fOUT);
  copyFile(streamName, fOUT);
  copyFile(deflinName, fOUT);

  fclose(fOUT);

  unlink(blocksName);
  unlink(deflinName);
  unlink(streamName);

  delete [] blocksName;
  delete [] deflinName;
  delete [] streamName;
  delete [] outputName;
}









merStreamFileReader::merStreamFileReader(const char *i) {
  _inputFile = i;

  char                 *streamName = new char [strlen(_inputFile) + 17];
  char                  cigam[16];

  sprintf(streamName, "%s.merStream", _inputFile);

  FILE *rawFile = fopen(streamName, "r");
  if (errno) {
    fprintf(stderr, "merStreamFileReader()-- Failed to open merStream '%s': %s\n", streamName, strerror(errno));
    exit(1);
  }

  _streamFile = new bitPackedFileReader(streamName);

  errno = 0;

  fread(cigam,             sizeof(char),   16, rawFile);
  fread(&_merSize,         sizeof(u32bit), 1,  rawFile);
  fread(&_merSize,         sizeof(u32bit), 1,  rawFile);
  fread(&_numMers,         sizeof(u64bit), 1,  rawFile);
  fread(&_numBlocks,       sizeof(u64bit), 1,  rawFile);
  fread(&_numDefs,         sizeof(u64bit), 1,  rawFile);
  fread(&_defLength,       sizeof(u64bit), 1,  rawFile);
  fread(&_blockFileSize,   sizeof(off_t),  1,  rawFile);
  fread(&_deflineFileSize, sizeof(off_t),  1,  rawFile);
  fread(&_streamFileSize,  sizeof(off_t),  1,  rawFile);
  fread(&_blkStart,        sizeof(off_t),  1,  rawFile);
  fread(&_defStart,        sizeof(off_t),  1,  rawFile);
  fread(&_strStart,        sizeof(off_t),  1,  rawFile);

  if (errno) {
    fprintf(stderr, "merStreamFileReader()-- Failed to read header: %s\n", strerror(errno));
    exit(1);
  }

  if (strncmp(magic, cigam, 16) != 0) {
    fprintf(stderr, "merStreamFileReader()-- '%s' isn't a merStream file.\n", streamName);
    exit(1);
  }

  fprintf(stderr, "merStreamFileReader()-- merSize:         "u32bitFMT" bases\n", _merSize);
  fprintf(stderr, "merStreamFileReader()-- numMers:         "u64bitFMT" mers\n", _numMers);
  fprintf(stderr, "merStreamFileReader()-- numBlocks:       "u64bitFMT" blocks\n", _numBlocks);
  fprintf(stderr, "merStreamFileReader()-- numDefs:         "u64bitFMT" lines\n", _numDefs);
  fprintf(stderr, "merStreamFileReader()-- defLength:       "u64bitFMT" characters\n", _defLength);
  fprintf(stderr, "merStreamFileReader()-- blockFileSize:   "u64bitFMT" bytes\n", _blockFileSize);
  fprintf(stderr, "merStreamFileReader()-- deflineFileSize: "u64bitFMT" bytes\n", _deflineFileSize);
  fprintf(stderr, "merStreamFileReader()-- streamFileSize:  "u64bitFMT" bytes\n", _streamFileSize);
  fprintf(stderr, "merStreamFileReader()-- blkStart:        "u64bitFMT" bits\n", _blkStart);
  fprintf(stderr, "merStreamFileReader()-- defStart:        "u64bitFMT" bits\n", _defStart);
  fprintf(stderr, "merStreamFileReader()-- strStart:        "u64bitFMT" bits\n", _strStart);

  //  Read the blocks
  //
  _streamFile->seek(_blkStart);

  _blockSizes = new u32bit [_numBlocks];

  fprintf(stderr, "WARNING:  Not using sequenceNumber and position information!\n");

  for (u64bit b=0; b<_numBlocks; b++) {
    _blockSizes[b] = (u32bit)_streamFile->getNumber();
    _streamFile->getNumber();  //  XXX:  Use this!
    _streamFile->getNumber();  //  XXX:  This too!
  }

  //  Read the deflines (maybe)
  //
  _deflineStorage = new char  [_defLength];
  _deflines       = new char* [_numDefs];

  errno = 0;

  fseeko(rawFile, _defStart >> 3, SEEK_SET);

  if (errno) {
    fprintf(stderr, "merStreamFileReader()-- Failed to position file to deflines: %s\n", strerror(errno));
    exit(1);
  }

  for (u32bit d=0, o=0; d<_numDefs; d++) {
    u32bit len = 0;

    _deflines[d] = _deflineStorage + o;

    errno = 0;

    fread(&len, sizeof(u32bit), 1, rawFile);
    fread(_deflines[d], sizeof(char), len, rawFile);

    if (errno) {
      fprintf(stderr, "merStreamFileReader()-- Failed to read defline "u32bitFMT": %s\n", d, strerror(errno));
      exit(1);
    }

    o += len;
  }

  //  We're all done with the raw file access, close it.
  //
  fclose(rawFile);

  //  Position the bitPackedFileReader to the start of the merStream
  //
  _streamFile->seek(_strStart);

  delete [] streamName;

  //  Setup the decoder.  _firstMer is a terribly ugly hack (or maybe
  //  it's a terribly elegant solution) needed to initialize the
  //  decoder.  The problem here is that the first block is 0, decoder
  //  increments _thisBlock immediately, and _thisBlock is unsigned --
  //  so we need to initialize it to -1, but we can't.
  //
  _thisBlock     = 0;
  _thisBlockSize = 0;
  _merMask       = u64bitMASK(_merSize * 2);
  _theFMer       = u64bitZERO;
  _firstMer      = u64bitZERO;
}


merStreamFileReader::~merStreamFileReader() {
  delete [] _blockSizes;
  delete [] _deflineStorage;
  delete [] _deflines;

  delete _streamFile;
}


bool
merStreamFileReader::seekToMer(u64bit merNumber) {
  u64bit  totalMers = 0;

  //  Did someone request a mer number too high?
  //
  if (merNumber >= _numMers)
    return(false);
  
  //  Search for the first block that is larger than merNumber
  //
  u32bit b = 0;

  while ((b < _numBlocks) &&
         (totalMers <= merNumber) &&
         (merNumber <= totalMers + _blockSizes[b])) {
    totalMers += _blockSizes[b];
    b++;
  }

  //  We're in the right block, we just need to load the mer and quit.
  //
  _streamFile->seek(_strStart + b * _merSize * 2 + (merNumber - b) * 2);

  _thisBlock     = b;
  _thisBlockSize = _blockSizes[_thisBlock] - (merNumber - totalMers);
  _theFMer       = _streamFile->getBits(_merSize * 2);

  return(true);
}


bool
merStreamFileReader::seekToSequence(u64bit seqNumber) {
  return(false);
}


bool
merStreamFileReader::nextMer(u32bit skip) {

again:
  if (_thisBlock >= _numBlocks)
    return(false);

  if (_thisBlockSize == 0) {
    //  Finished a block, move to the next one.
    //
    _thisBlock++;

    //  To get around an extra if here we do a multiply to reset
    //  _thisBlock to the first block if we are the first mer
    //  (_firstMer is initially set to zero).
    //
    _thisBlock &= _firstMer;
    _firstMer   = ~u64bitZERO;

    //  Outta blocks?
    if (_thisBlock >= _numBlocks)
      return(false);

    _thisBlockSize = _blockSizes[_thisBlock];
    _theFMer = _streamFile->getBits(_merSize * 2);
  } else {
    //  Still in the same block, just move to the next mer.
    //
    _theFMer <<= 2;
    _theFMer  |= _streamFile->getBits(2);
    _theFMer  &= _merMask;
  }

  _thisBlockSize--;

  if (skip--)
    goto again;

  return(true);
}




#ifdef TEST_MERSTREAMFILE

int
main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "usage: %s some.fasta\n", argv[0]);
    fprintf(stderr, "       Builds a merStreamFile, and then checks that it returns\n");
    fprintf(stderr, "       exactly the same stuff as a merStream(\"some.fasta\") does.\n");
    fprintf(stderr, "       Returns 1 if error, 0 if OK\n");
    exit(1);
  }

  merStreamFileBuilder   *B = new merStreamFileBuilder(20, argv[1], "merStreamFileTest");
  B->build(true);
  delete B;

  merStreamFileReader    *R = new merStreamFileReader("merStreamFileTest");
  merStream              *M = new merStream(20, argv[1]);

  u64bit compared = 0;
  u64bit errors   = 0;

  while (M->nextMer() && R->nextMer()) {
    compared++;
    if (M->theFMer() != R->theFMer()) {
      fprintf(stderr, u64bitFMT": M got "u64bitHEX" but R got "u64bitHEX"\n", compared, M->theFMer(), R->theFMer());
      errors++;
    }
  }

  fprintf(stderr, "Compared "u64bitFMT" mers.\n", compared);

  if (M->nextMer()) {
    fprintf(stderr, "ERROR: Extra mers in the merStream!\n");
    errors++;
  }

  if (R->nextMer()) {
    fprintf(stderr, "ERROR: Extra mers in the merStreamFile!\n");
    errors++;
  }

  delete M;
  delete R;

  if (errors > 0) {
    fprintf(stderr, "There were "u64bitFMT" errors.\n", errors);
    exit(1);
  } else {
    fprintf(stderr, "merStreamFile works correctly!\n");
    exit(0);
  }
}

#endif  //  TEST_MERSTREAMFILE

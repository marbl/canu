#include "bio++.H"

static
char     magic[16] = { 'm', 'e', 'r', 'S', 't', 'r', 'e', 'a', 'm', '1', ' ', ' ', ' ', ' ', ' ', ' '  };

//  Define this to get diagnostics when the block number changes
//#define REPORT_BLOCK_CHANGES



//  Returns true if the merStreamFile exists.
bool
merStreamFileExists(const char *i) {
  char  *streamName = new char [strlen(i) + 17];
  sprintf(streamName, "%s.merStream", i);
  bool exists = fileExists(streamName);
  delete [] streamName;
  return(exists);
}




merStreamFileBuilder::merStreamFileBuilder(u32bit m, const char *i, const char *o) {
  _merSize     = m;
  _fastaStream = new FastAstream(i);
  _merStream   = new merStream(m, _fastaStream);
  _outputFile  = o;
}

merStreamFileBuilder::merStreamFileBuilder(merStream *MS, const char *o) {
  _merSize     = MS->theMerSize();
  _fastaStream = 0L;
  _merStream   = MS;
  _outputFile  = o;
}

merStreamFileBuilder::~merStreamFileBuilder() {
  if (_fastaStream) {
    delete _fastaStream;
    delete _merStream;
  }
}

u64bit
merStreamFileBuilder::build(bool beVerbose) {

  u64bit             thisMer           = u64bitZERO;

  u64bit             blockSize         = 0;
  u64bit             numMers           = 0;
  u64bit             numBlocks         = 0;
  u64bit             numDefs           = 0;
  u64bit             defLength         = 0;

  //  Open the output files
  //
  char *blocksName = new char [strlen(_outputFile) + 32];
  char *deflinName = new char [strlen(_outputFile) + 32];
  char *streamName = new char [strlen(_outputFile) + 32];
  char *outputName = new char [strlen(_outputFile) + 32];

  sprintf(blocksName, "%s.merstreamfilebuilder.b.tmp", _outputFile);
  sprintf(deflinName, "%s.merstreamfilebuilder.d.tmp", _outputFile);
  sprintf(streamName, "%s.merstreamfilebuilder.s.tmp", _outputFile);
  sprintf(outputName, "%s.merStream",                  _outputFile);

  bitPackedFile *BLOCKS = new bitPackedFile(blocksName);
  bitPackedFile *STREAM = new bitPackedFile(streamName);

  errno = 0;
  FILE  *DEFLIN  = fopen(deflinName, "w");
  if (errno) {
    fprintf(stderr, "merStreamFileBuilder::build()-- failed to open defline storage: %s\n", strerror(errno));
    exit(1);
  }


  speedCounter       C("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);

  //  Write the first mer, and initialize things.
  //
  _merStream->nextMer();
  
  u64bit lastMer           = _merStream->theFMer();
  u64bit lastMerPos        = _merStream->thePositionInStream();
  u64bit lastSeq           = _merStream->theSequenceNumber();
  u64bit lastBlockPosInSeq = _merStream->thePositionInSequence();
  u64bit lastBlockPosInStr = _merStream->thePositionInStream();
  u32bit lastDef           = ~u32bitZERO;

  STREAM->putBits(lastMer, _merSize * 2);
  blockSize = 1;
  numMers   = 1;

  _merStream->nextMer();

  do {
    thisMer = _merStream->theFMer();

    //  Save the defline if the sequence number has changed.
    //
    if (lastDef != _merStream->theSequenceNumber()) {
      u32bit  l = strlen(_merStream->theDefLine()) + 1;

      fwrite(&l, sizeof(u32bit), 1, DEFLIN);
      fwrite(_merStream->theDefLine(), sizeof(char), l, DEFLIN);

      lastDef = _merStream->theSequenceNumber();

      numDefs++;
      defLength += l;
    }

    //  See if this mer is adjacent to the last one.  We used to check
    //  for overlap in the mer, which worked great for storing just
    //  the mers, but it fails if we also store positions.
    //
    //  u64bit  maskMer = u64bitMASK(_merSize * 2 - 2);
    //  (thisMer >> 2) == (lastMer & maskMer)
    //
    lastMerPos++;

    if (lastMerPos == _merStream->thePositionInStream()) {
      //  Consecutive mers, just emit the new base, and add one to our count
      //
      STREAM->putBits(thisMer & 0x03, 2);
      blockSize++;
    } else {
      //  Must be a mer break.  Write the count, and the new mer.
      //

#ifdef REPORT_BLOCK_CHANGES
      fprintf(stderr, "new block="u64bitFMT" size="u64bitFMT" pos="u64bitFMT"/"u64bitFMT" seq="u64bitFMT"\n",
              numBlocks, blockSize, lastBlockPosInSeq, lastBlockPosInStr, lastSeq);
#endif

      STREAM->putBits(thisMer, _merSize * 2);
      BLOCKS->putNumber(blockSize);
      BLOCKS->putNumber(lastSeq);
      BLOCKS->putNumber(lastBlockPosInSeq);
      BLOCKS->putNumber(lastBlockPosInStr);

      lastSeq           = _merStream->theSequenceNumber();
      lastBlockPosInSeq = _merStream->thePositionInSequence();
      lastBlockPosInStr = _merStream->thePositionInStream();
      lastMerPos        = lastBlockPosInStr;

      numBlocks++;
      blockSize = 1;
    }

    numMers++;

    lastMer = thisMer;

    C.tick();
  } while (_merStream->nextMer());

  //  Be sure to write out the last block!
  //
  BLOCKS->putNumber(blockSize);
  BLOCKS->putNumber(lastSeq);
  BLOCKS->putNumber(lastBlockPosInSeq);
  BLOCKS->putNumber(lastBlockPosInStr);
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

#if 0
  fprintf(stderr, "merStreamFileWriter()-- merSize:         "u32bitFMT" bases\n", _merSize);
  fprintf(stderr, "merStreamFileWriter()-- numMers:         "u64bitFMT" mers\n", numMers);
  fprintf(stderr, "merStreamFileWriter()-- numBlocks:       "u64bitFMT" blocks\n", numBlocks);
  fprintf(stderr, "merStreamFileWriter()-- numDefs:         "u64bitFMT" lines\n", numDefs);
  fprintf(stderr, "merStreamFileWriter()-- defLength:       "u64bitFMT" characters\n", defLength);
  fprintf(stderr, "merStreamFileWriter()-- blockFileSize:   "u64bitFMT" bytes\n", blockFileSize);
  fprintf(stderr, "merStreamFileWriter()-- deflineFileSize: "u64bitFMT" bytes\n", deflineFileSize);
  fprintf(stderr, "merStreamFileWriter()-- streamFileSize:  "u64bitFMT" bytes\n", streamFileSize);
  fprintf(stderr, "merStreamFileWriter()-- blkStart:        "u64bitFMT" bits\n", blkStart);
  fprintf(stderr, "merStreamFileWriter()-- defStart:        "u64bitFMT" bits\n", defStart);
  fprintf(stderr, "merStreamFileWriter()-- strStart:        "u64bitFMT" bits\n", strStart);
#endif


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

  if (beVerbose)
    fprintf(stderr, "\n");

  return(numMers);
}









merStreamFileReader::merStreamFileReader(const char *i) {
  _inputFile = i;

  char                 *streamName = new char [strlen(_inputFile) + 17];
  char                  cigam[16] = {0};

  sprintf(streamName, "%s.merStream", _inputFile);

  errno = 0;
  FILE *rawFile = fopen(streamName, "r");
  if (errno) {
    fprintf(stderr, "merStreamFileReader()-- Failed to open merStream '%s': %s\n", streamName, strerror(errno));
    exit(1);
  }

  _streamFile = new bitPackedFile(streamName);

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

#if 0
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
#endif

  //  Read the blocks
  //
  _streamFile->seek(_blkStart);

  _blockSize      = new u32bit [_numBlocks];
  _blockSequence  = new u32bit [_numBlocks];
  _blockPosInSeq  = new u64bit [_numBlocks];
  _blockPosInStr  = new u64bit [_numBlocks];

  for (u64bit b=0; b<_numBlocks; b++) {
    _blockSize[b]     = (u32bit)_streamFile->getNumber();
    _blockSequence[b] = (u32bit)_streamFile->getNumber();
    _blockPosInSeq[b] = (u64bit)_streamFile->getNumber();
    _blockPosInStr[b] = (u64bit)_streamFile->getNumber();

#if 0
    fprintf(stderr, "Block "u64bitFMTW(2)" : "u32bitFMT" "u32bitFMT" "u64bitFMT" "u64bitFMT"\n",
            b, _blockSize[b], _blockSequence[b], _blockPosInSeq[b], _blockPosInStr[b]);
#endif
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

  //  Position the bitPackedFile to the start of the merStream
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

  _firstMer      = u64bitZERO;

  _iterationLimit = ~u64bitZERO;
  _iteration      =  u64bitZERO;

  _fMer       = u64bitZERO;
  _rMer       = u64bitZERO;
  _posInSeq   = 0;
  _posInStr   = 0;
  _sequence   = 0;

  _defline    = 0L;
}


merStreamFileReader::~merStreamFileReader() {
  delete [] _blockSize;
  delete [] _blockSequence;
  delete [] _blockPosInSeq;
  delete [] _blockPosInStr;
  delete [] _deflineStorage;
  delete [] _deflines;
  delete    _streamFile;
}


bool
merStreamFileReader::seekToMer(u64bit merNumber) {

  //  Did someone request a mer number too high?
  //
  if (merNumber >= _numMers) {
    _thisBlock = _numBlocks;
    return(false);
  }

  //  Special case for merNumber == 0.  We treat it just like a normal
  //  creation.
  //
  if (merNumber == 0) {
    _thisBlock     = 0;
    _thisBlockSize = 0;

    _firstMer      = u64bitZERO;

    _fMer       = u64bitZERO;
    _rMer       = u64bitZERO;
    _posInSeq   = 0;
    _posInStr   = 0;
    _sequence   = 0;

    _defline    = 0L;

    _streamFile->seek(_strStart);

    return(true);
  }

  //  Otherwise, it's not the first mer, and we should clear the firstMer reset
  //
  _firstMer = ~u64bitZERO;
  
  //  Search for the first block that is larger than merNumber
  //
  u32bit block     = 0;
  u64bit totalMers = 0;

  while ((block < _numBlocks) &&
         (totalMers + _blockSize[block] <= merNumber)) {
    //fprintf(stderr, "mer=%llu block=%lu total=%llu len=%lu\n", merNumber, block, totalMers, _blockSize[block]);
    totalMers += _blockSize[block];
    block++;
  }

  //fprintf(stderr, "mer=%llu block=%lu total=%llu len=%lu\n", merNumber, block, totalMers, _blockSize[block]);

  //  We're in the right block, we just need to load the mer and quit.
  //
  _streamFile->seek(_strStart + block * _merSize * 2 + (merNumber - block) * 2);

  _thisBlock     = block;
  _thisBlockSize = _blockSize[_thisBlock] - (merNumber - totalMers);

  _posInSeq  = _blockPosInSeq[_thisBlock] + (merNumber - totalMers) - 1;
  _posInStr  = _blockPosInStr[_thisBlock] + (merNumber - totalMers) - 1;
  _sequence  = _blockSequence[_thisBlock];

  //  We read in mersize-1 bases into the FMer, but reverse complement
  //  the whole mer.  On the nextMer(), we shift out the extra junk
  //  base from both.
  //
  _fMer       = _streamFile->getBits(_merSize * 2 - 2);
  _rMer       = reverseComplementMer(_merSize, _fMer);

#ifdef REPORT_BLOCK_CHANGES
  fprintf(stderr, "seek to block="u64bitFMT" size="u64bitFMT" pos="u64bitFMT"/"u64bitFMT" seq="u64bitFMT"\n",
          _thisBlock, _thisBlockSize, _posInSeq, _posInStr, _sequence);
#endif

  return(true);
}


bool
merStreamFileReader::seekToSequence(u64bit seqNumber) {
  return(false);
}


bool
merStreamFileReader::nextMer(u32bit skip) {

  if (_iteration >= _iterationLimit)
    return(false);

  do {
    if (_thisBlock >= _numBlocks)
      return(false);

    if (_thisBlockSize == 0) {
      //  Finished a block, move to the next one.
      //
      _thisBlock++;

      //  To avoid another if here we do a multiply to reset _thisBlock
      //  to the first block if we are the first mer (_firstMer is
      //  initially set to zero).
      //
      _thisBlock &= _firstMer;
      _firstMer   = ~u64bitZERO;

      //  Outta blocks?
      if (_thisBlock >= _numBlocks)
        return(false);

      _thisBlockSize = _blockSize[_thisBlock] - 1;

      _posInSeq  = _blockPosInSeq[_thisBlock];
      _posInStr  = _blockPosInStr[_thisBlock];
      _sequence  = _blockSequence[_thisBlock];

      _fMer = _streamFile->getBits(_merSize * 2);
      _rMer = reverseComplementMer(_merSize, _fMer);

#ifdef REPORT_BLOCK_CHANGES
      fprintf(stderr, "New block="u64bitFMT" size="u64bitFMT" pos="u64bitFMT"/"u64bitFMT" seq="u64bitFMT"\n",
              _thisBlock, _thisBlockSize, _posInSeq, _posInStr, _sequence);
#endif
    } else {
      //  Still in the same block, just move to the next mer.
      //
      u64bit  theBase = _streamFile->getBits(2);

      _fMer <<= 2;
      _fMer  |= theBase;
      _fMer  &= _merMask;

      theBase  ^= 0x3;
      theBase <<= (_merSize << 1) - 2;

      _rMer >>= 2;
      _rMer  |= theBase;

      _thisBlockSize--;

      _posInSeq++;
      _posInStr++;
    }
  } while (skip--);

  _iteration++;

  return(true);
}





char const *
merStreamFileReader::theFMerString(void) {
  for (u32bit i=0; i<_merSize; i++)
    _merString[_merSize-i-1] = decompressSymbol[(_fMer >> (2*i)) & 0x03];
  _merString[_merSize] = 0;
  return(_merString);
}

char const *
merStreamFileReader::theRMerString(void) {
  for (u32bit i=0; i<_merSize; i++)
    _merString[_merSize-i-1] = decompressSymbol[(_rMer >> (2*i)) & 0x03];
  _merString[_merSize] = 0;
  return(_merString);
}


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
  _chainedSeq = new chainedSequence();
  _chainedSeq->setSource(i);
  _chainedSeq->finish();
  _merStream   = new merStream(m, _chainedSeq);
  _outputFile  = o;
}

merStreamFileBuilder::merStreamFileBuilder(merStream *MS, const char *o) {
  _merSize     = MS->theMerSize();
  _chainedSeq  = 0L;
  _merStream   = MS;
  _outputFile  = o;
}

merStreamFileBuilder::~merStreamFileBuilder() {
  if (_chainedSeq) {
    delete _chainedSeq;
    delete _merStream;
  }
}

u64bit
merStreamFileBuilder::build(bool beVerbose) {
  u64bit             blockSize         = 0;
  u64bit             numMers           = 0;
  u64bit             numBlocks         = 0;

  //  Open the output files
  //
  char *blocksName = new char [strlen(_outputFile) + 32];
  char *streamName = new char [strlen(_outputFile) + 32];
  char *outputName = new char [strlen(_outputFile) + 32];

  sprintf(blocksName, "%s.merstreamfilebuilder.b.tmp", _outputFile);
  sprintf(streamName, "%s.merstreamfilebuilder.s.tmp", _outputFile);
  sprintf(outputName, "%s.merStream",                  _outputFile);

  bitPackedFile *BLOCKS = new bitPackedFile(blocksName);
  bitPackedFile *STREAM = new bitPackedFile(streamName);


  speedCounter       C("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);

  //  Write the first mer, and initialize things.
  //
  _merStream->nextMer();
  
  kMer   thisMer;
  kMer   lastMer           = _merStream->theFMer();
  u64bit lastMerPos        = _merStream->thePositionInStream();
  u64bit lastSeq           = _merStream->theSequenceNumber();
  u64bit lastBlockPosInSeq = _merStream->thePositionInSequence();
  u64bit lastBlockPosInStr = _merStream->thePositionInStream();

  lastMer.writeToBitPackedFile(STREAM);

  blockSize = 1;
  numMers   = 1;

  _merStream->nextMer();

  do {
    thisMer = _merStream->theFMer();

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
      STREAM->putBits(thisMer.subseq_mask(0x0000000000000003), 2);
      blockSize++;
    } else {
      //  Must be a mer break.  Write the count, and the new mer.
      //

#ifdef REPORT_BLOCK_CHANGES
      fprintf(stderr, "new block="u64bitFMT" size="u64bitFMT" pos="u64bitFMT"/"u64bitFMT" seq="u64bitFMT"\n",
              numBlocks, blockSize, lastBlockPosInSeq, lastBlockPosInStr, lastSeq);
#endif

      thisMer.writeToBitPackedFile(STREAM);
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
  off_t streamFileSize  = sizeOfFile(streamName);

  off_t blkStart = sizeof(char) * 16 + sizeof(u32bit) * 2 + sizeof(u64bit) * 2 + sizeof(off_t) * 4;
  off_t strStart = blkStart + blockFileSize;

  //
  //  This block MUST be a multiple of 64 bits long to keep the
  //  bitPackedFiles 64-bit word aligned.
  //

  fwrite(magic,            sizeof(char),   16, fOUT);
  fwrite(&_merSize,        sizeof(u32bit), 1,  fOUT);  //  Padding.
  fwrite(&_merSize,        sizeof(u32bit), 1,  fOUT);
  fwrite(&numMers,         sizeof(u64bit), 1,  fOUT);
  fwrite(&numBlocks,       sizeof(u64bit), 1,  fOUT);
  fwrite(&blockFileSize,   sizeof(off_t),  1,  fOUT);
  fwrite(&streamFileSize,  sizeof(off_t),  1,  fOUT);
  fwrite(&blkStart,        sizeof(off_t),  1,  fOUT);
  fwrite(&strStart,        sizeof(off_t),  1,  fOUT);

  copyFile(blocksName, fOUT);
  copyFile(streamName, fOUT);

  fclose(fOUT);

  unlink(blocksName);
  unlink(streamName);

  delete [] blocksName;
  delete [] streamName;
  delete [] outputName;

  if (beVerbose)
    fprintf(stderr, "\n");

  return(numMers);
}









merStreamFileReader::merStreamFileReader(const char *i, u32bit desiredMerSize) {
  _inputFile = i;

  char                 *streamName = new char [strlen(_inputFile) + 17];
  char                  cigam[16] = {0};

  sprintf(streamName, "%s.merStream", _inputFile);

  errno = 0;
  FILE *rawFile = fopen(streamName, "r");
  if (errno)
    fprintf(stderr, "merStreamFileReader()-- Failed to open merStream '%s': %s\n", streamName, strerror(errno)), exit(1);

  fread(cigam,             sizeof(char),   16, rawFile);
  fread(&_merSizeInFile,   sizeof(u32bit), 1,  rawFile);  //  Padding.
  fread(&_merSizeInFile,   sizeof(u32bit), 1,  rawFile);
  fread(&_numMers,         sizeof(u64bit), 1,  rawFile);
  fread(&_numBlocks,       sizeof(u64bit), 1,  rawFile);
  fread(&_blockFileSize,   sizeof(off_t),  1,  rawFile);
  fread(&_streamFileSize,  sizeof(off_t),  1,  rawFile);
  fread(&_blkStart,        sizeof(off_t),  1,  rawFile);
  fread(&_strStart,        sizeof(off_t),  1,  rawFile);
  fclose(rawFile);

  if (errno)
    fprintf(stderr, "merStreamFileReader()-- Failed to read header: %s\n", strerror(errno)), exit(1);

  if (strncmp(magic, cigam, 16) != 0)
    fprintf(stderr, "merStreamFileReader()-- '%s' isn't a merStream file.\n", streamName), exit(1);

  if (desiredMerSize)
    _merSizeDesired = desiredMerSize;
  else
    _merSizeDesired = _merSizeInFile;

  if (_merSizeDesired < _merSizeInFile)
    fprintf(stderr, "merStreamFileReader()-- '%s' has a minimum merSize of "u32bitFMT", but you requested size "u32bitFMT".\n",
            streamName, _merSizeInFile, _merSizeDesired), exit(1);

#if 0
  fprintf(stderr, "merStreamFileReader()-- merSizeDesired:  "u32bitFMT" bases\n", _merSizeDesired);
  fprintf(stderr, "merStreamFileReader()-- merSizeInFile:   "u32bitFMT" bases\n", _merSizeInFile);
  fprintf(stderr, "merStreamFileReader()-- numMers:         "u64bitFMT" mers\n", _numMers);
  fprintf(stderr, "merStreamFileReader()-- numBlocks:       "u64bitFMT" blocks\n", _numBlocks);
  fprintf(stderr, "merStreamFileReader()-- blockFileSize:   "u64bitFMT" bytes\n", _blockFileSize);
  fprintf(stderr, "merStreamFileReader()-- streamFileSize:  "u64bitFMT" bytes\n", _streamFileSize);
  fprintf(stderr, "merStreamFileReader()-- blkStart:        "u64bitFMT" bits\n", _blkStart);
  fprintf(stderr, "merStreamFileReader()-- strStart:        "u64bitFMT" bits\n", _strStart);
#endif

  //  Read the blocks
  //
  bitPackedFile *blocksFile = new bitPackedFile(streamName, _blkStart);

  _blockSize               = new u32bit [_numBlocks];
  _blockSequence           = new u32bit [_numBlocks];
  _blockPosInSeq           = new u64bit [_numBlocks];
  _blockPosInStr           = new u64bit [_numBlocks];
  _blockFirstMerNum        = new u64bit [_numBlocks + 1];
  _blockTotalMerCorrection = new u32bit [_numBlocks + 1];

  //  We (re)compute the number of mers in this file, given that the
  //  desired mersize might be larger than that used to construct the
  //  file.
  //
  u64bit   origNumMers          = _numMers;
  u64bit   totMers              = 0;
  u32bit   merSizeDifference    = _merSizeDesired - _merSizeInFile;

  _numMers = 0;

  for (u64bit b=0; b<_numBlocks; b++) {
    _blockSize[b]               = (u32bit)blocksFile->getNumber();
    _blockSequence[b]           = (u32bit)blocksFile->getNumber();
    _blockPosInSeq[b]           = (u64bit)blocksFile->getNumber();
    _blockPosInStr[b]           = (u64bit)blocksFile->getNumber();
    _blockFirstMerNum[b]        = _numMers;
    _blockTotalMerCorrection[b] = totMers - _numMers;

    totMers += _blockSize[b];

    //  Because we now allow building on some mer size k, and accesses
    //  for all mersizes k or bigger, the number of mers per block can
    //  be smaller than what the file says it is.  We explicitly count
    //  how many mers at the requested size there are.
    //
    if (_blockSize[b] > merSizeDifference)
      _numMers += _blockSize[b] - merSizeDifference;
  }

  //  Useful for the binary search below!
  _blockFirstMerNum[_numBlocks] = _numMers;

  fprintf(stderr, "Found "u64bitFMT" mers at size "u32bitFMT", and "u64bitFMT" mers at desired size of "u32bitFMT"\n",
          origNumMers, _merSizeInFile, _numMers, _merSizeDesired);

  delete blocksFile;


  //  Position the bitPackedFile to the start of the merStream
  //
  _streamFile = new bitPackedFile(streamName, _strStart);
  delete [] streamName;

  //  Setup the decoder.  _firstMer is a terribly ugly hack (or maybe
  //  it's a terribly elegant solution) needed to initialize the
  //  decoder.  The problem here is that the first block is 0, decoder
  //  increments _thisBlock immediately, and _thisBlock is unsigned --
  //  so we need to initialize it to -1, but we can't.
  //
  _thisBlock     = 0;
  _thisBlockSize = 0;

  _firstMer      = u64bitZERO;

  _iterationStart =  u64bitZERO;
  _iterationLimit = ~u64bitZERO;
  _iteration      =  u64bitZERO;

  _fMer.setMerSize(_merSizeDesired);
  _rMer.setMerSize(_merSizeDesired);

  _fMer.clear();
  _rMer.clear();

  _posInSeq   = 0;
  _posInStr   = 0;
  _sequence   = 0;
}


merStreamFileReader::~merStreamFileReader() {
  delete [] _blockSize;
  delete [] _blockSequence;
  delete [] _blockPosInSeq;
  delete [] _blockPosInStr;
  delete [] _blockFirstMerNum;
  delete    _streamFile;
}



void
merStreamFileReader::findBlock(u64bit merNumber,
                               u32bit  &block,
                               u64bit  &totalMers,
                               u64bit  &totalMersInFile) {

  u32bit merSizeDifference  = _merSizeDesired - _merSizeInFile;
  bool   stillLooking       = true;

  //  Use a binary search if there are lots of blocks, otherwise, linear search.

  if (_numBlocks < 128) {
    //
    //  Linear search
    //
    while ((stillLooking) && (block < _numBlocks)) {
      if        (_blockSize[block] < merSizeDifference) {
        //  Block size too small to hold this mersize, keep looking.
        totalMersInFile += _blockSize[block];
        block++;

      } else if (totalMers + _blockSize[block] - merSizeDifference > merNumber) {
        //  Got it!  Stop.
        stillLooking = false;

      } else {
        //  Big block, but not the right one.  Increment the sizes.
        totalMersInFile += _blockSize[block];
        totalMers       += _blockSize[block] - merSizeDifference;
        block++;
      }
    }
  } else {
    //
    //  Binary search
    //

    u32bit  lo = 0;
    u32bit  hi = _numBlocks;
    u32bit  md = 0;

    while (lo <= hi) {
      md = (lo + hi) / 2;

      if        (merNumber < _blockFirstMerNum[md]) {
        //  This block starts after the mer we're looking for.  
        //
        hi = md;

      } else if ((_blockFirstMerNum[md] <= merNumber) && (merNumber < _blockFirstMerNum[md+1])) {
        //  Got it!
        lo           = md + 1;
        hi           = md;
        block        = md;
        stillLooking = false;

      } else {
        //  By default, then, the block is too low.
        lo = md;
      }
    }

    totalMersInFile = _blockFirstMerNum[md] + _blockTotalMerCorrection[md];
    totalMers       = _blockFirstMerNum[md];
  }


#if 0
  fprintf(stderr, "mer="u64bitFMT"  block="u32bitFMT" total="u64bitFMT"/"u64bitFMT"\n",
          merNumber, block, totalMersInFile, totalMers);
#endif

  //  stillLooking is true if we failed to find the guy.
  //
  if (stillLooking) {
    fprintf(stderr, "merStreamFileReader::seekToMer()--  Looking for "u64bitFMT", but didn't find it.\n", merNumber);
    exit(1);
  }
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

    _fMer.clear();
    _rMer.clear();

    _posInSeq   = 0;
    _posInStr   = 0;
    _sequence   = 0;

    _streamFile->seek(0);

    return(true);
  }

  //  Otherwise, it's not the first mer, and we should clear the firstMer reset
  //
  _firstMer = ~u64bitZERO;
  
  //  Search for the first block that is larger than merNumber
  //
  u32bit block              = 0;
  u64bit totalMers          = 0;  //  total mers at desired size
  u64bit totalMersInFile    = 0;  //  total mers we see at default size

  findBlock(merNumber, block, totalMers, totalMersInFile);


  //  Now we're in the right block, and just need to load the mer.
  //
  //  If the actual and desired sizes are the same, we would do:
  //    _streamFile->seek(_strStart + block * _merSize * 2 + (merNumber - block) * 2);
  //
  //  That's just one full mer for each block, plus 2 bits for each
  //  additional mer.  Our mersize changed, and that will change the
  //  number of mers we see before this block starts (that's why we
  //  counted both totalMers and totalMersInFile), but once we're in
  //  the block, the offset to the _start_ of the mer doesn't differ
  //  (the end of the mer does).
  //
  _streamFile->seek(block * _merSizeInFile * 2 +      //  first mer in n-1 blocks before this one
                    (totalMersInFile - block) * 2 +   //  rest of mers in those blocks
                    (merNumber - totalMers) * 2);     //  mers in this block before the mer we want

  _thisBlock     = block;
  _thisBlockSize = _blockSize[_thisBlock] - (merNumber - totalMers) - (_merSizeDesired - _merSizeInFile);

  _posInSeq      = _blockPosInSeq[_thisBlock] + (merNumber - totalMers) - 1;
  _posInStr      = _blockPosInStr[_thisBlock] + (merNumber - totalMers) - 1;
  _sequence      = _blockSequence[_thisBlock];

  //  We read in mersize-1 bases into the FMer, but reverse complement
  //  the whole mer.  On the nextMer(), we shift out the extra junk
  //  base from both.
  //
  //  _fMer       = _streamFile->getBits(_merSize * 2 - 2);
  //  _rMer       = _fMer;
  //  _rMer.reverseComplement();
  //
  //  This is pretty gross, but the simplest way I could think to do
  //  it without violating the (ever growing) kMer interface.  It does
  //  have the nice side effect of building the rMer correctly -- so
  //  we don't need to use the kMer reverseComplement() (which isn't
  //  supported right now anyway!)
  //
  //  It's also used below, so fix that too.
  //
  _fMer.clear();
  _rMer.clear();

  for (u32bit i=1; i<_merSizeDesired; i++) {
    u64bit  base = _streamFile->getBits(2);
    _fMer += base;
    _rMer -= base ^ 0x3;
  }

#ifdef REPORT_BLOCK_CHANGES
  fprintf(stderr, "seek to block="u64bitFMT" size="u32bitFMT",remain="u64bitFMT" pos="u64bitFMT"/"u64bitFMT" seq="u64bitFMT"\n",
          _thisBlock, _blockSize[_thisBlock], _thisBlockSize, _posInSeq, _posInStr, _sequence);
#endif

  return(true);
}




bool
merStreamFileReader::nextMer(u32bit skip) {

  if (_iteration >= _iterationLimit)
    return(false);

  //fprintf(stderr, "_thisBlockSize="u64bitFMT"\n", _thisBlockSize);

  do {
    if (_thisBlock >= _numBlocks)
      return(false);

    if (_thisBlockSize == 0) {
      //  Finished a block, move to the next one.
      //
      _thisBlock++;

      //  To avoid another 'if' here we do a multiply to reset
      //  _thisBlock to the first block if we are the first mer
      //  (_firstMer is initially set to zero).
      //
      _thisBlock &= _firstMer;
      _firstMer   = ~u64bitZERO;

      //  Skip small blocks
      //
      u32bit merSizeDifference  = _merSizeDesired - _merSizeInFile;

      while ((_thisBlock < _numBlocks) &&
             (_blockSize[_thisBlock] <= merSizeDifference)) {
#ifdef REPORT_BLOCK_CHANGES
        fprintf(stderr, "Skipping block "u64bitFMT" with size "u32bitFMT"\n", _thisBlock, _blockSize[_thisBlock]);
#endif

        for (u32bit i=0; i<_merSizeInFile + _blockSize[_thisBlock] - 1; i++)
          _streamFile->getBits(2);

        _thisBlock++;
      }

      //  Outta blocks?
      //
      if (_thisBlock >= _numBlocks)
        return(false);

      _thisBlockSize = _blockSize[_thisBlock] - merSizeDifference - 1;

      _posInSeq  = _blockPosInSeq[_thisBlock];
      _posInStr  = _blockPosInStr[_thisBlock];
      _sequence  = _blockSequence[_thisBlock];

      //  The same ugly hack as above -- but we read the whole mer,
      //  not one base less.
      //
      //  _fMer = _streamFile->getBits(_merSize * 2);
      //  _rMer = reverseComplementMer(_merSize, _fMer);
      //
      _fMer.clear();
      _rMer.clear();

      for (u32bit i=0; i<_merSizeDesired; i++) {
        u64bit  base = _streamFile->getBits(2);
        _fMer += base;
        _rMer -= base ^ 0x3;
      }

#ifdef REPORT_BLOCK_CHANGES
      fprintf(stderr, "New block="u64bitFMT" size="u64bitFMT" pos="u64bitFMT"/"u64bitFMT" seq="u64bitFMT"\n",
              _thisBlock, _thisBlockSize, _posInSeq, _posInStr, _sequence);
#endif
    } else {
      //  Still in the same block, just move to the next mer.
      //
      u64bit  theBase = _streamFile->getBits(2);

      _fMer += theBase;
      _rMer -= theBase ^ 0x3;

      _thisBlockSize--;

      _posInSeq++;
      _posInStr++;
    }
  } while (skip--);

  _fMer.mask(true);
  _rMer.mask(false);

  _iteration++;

  return(true);
}




bool
merStreamFileReader::flipToCanonicalize(u64bit merNum) {
  return(false);
};


bool
merStreamFileReader::smaller(u64bit merNum1,
                             u64bit merNum2,
                             bool   doCanonical) {
  bool  mer1flip = false;
  bool  mer2flip = false;

  if ((merNum1 >= _numMers) || (merNum2 >= _numMers)) {
    fprintf(stderr, "merStreamFileReader::smaller()-- merNum out of range: "u64bitFMT" "u64bitFMT" but only "u64bitFMT" mers.\n",
            merNum1, merNum2, _numMers);
    return(false);
  }

  //  Get rid of the easy case
  //
  if (doCanonical == false) {

    
  }



  if (doCanonical) {
    mer1flip = flipToCanonicalize(merNum1);
    mer2flip = flipToCanonicalize(merNum2);
  }



  return(false);
}

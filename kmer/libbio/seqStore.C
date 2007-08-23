
seqStore::seqStore() {
  _numBlocks        = 0;
  _blockInfo        = 0L;

  _thisBlock         = 0;
  _thisBlockPosition = 0;

  _numberOfACGT     = 0;

  _streamFile       = 0L;

  _iterationStart =  u64bitZERO;
  _iterationLimit = ~u64bitZERO;
  _iteration      =  u64bitZERO;
}

seqStore::seqStore(const char *filename, seqStream *s) {
  _numBlocks        = 0;
  _blockInfo        = 0L;

  _thisBlock         = 0;
  _thisBlockPosition = 0;

  _numberOfACGT     = 0;

  _streamFile       = 0L;

  _iterationStart =  u64bitZERO;
  _iterationLimit = ~u64bitZERO;
  _iteration      =  u64bitZERO;

  buildStore(filename, s);
}

seqStore::~seqStore() {
  delete [] _blockInfo;
  delete    _streamFile;
}

//  If file exists, and is a seqStore file, it is loaded, return true.
//  If file exists, but is not a seqStore, we fail (and don't overwrite the file).
//  Otherwise, return false.
//
bool
seqStore::loadStore(const char *filename) {
  char    cigam[16] = {0};
  char    magic[16] = {0};

  if ((filename == 0L) || (filename[0] == 0) || (filename[0] == '-'))
    return(false);

  char *streamName = new char [strlen(filename) + 32];
  char *blocksName = new char [strlen(filename) + 32];

  sprintf(streamName, "%s.seqStore.sequence", filename);
  sprintf(blocksName, "%s.seqStore.blocks",   filename);

  if ((fileExists(streamName) == false) ||
      (fileExists(blocksName) == false))
    return(false);

#if 1
  //  If 'filename' exists, assume it's the source fasta file.  We
  //  want to ensure that our seqStream is at least as current as
  //  that.
  if (fileExists(filename)) {
    u64bit  timeOfOriginal = timeOfFile(filename);
    if ((timeOfOriginal > timeOfFile(streamName)) ||
        (timeOfOriginal > timeOfFile(streamName))) {
      fprintf(stderr, "seqStore::loadStore()-- WARNING: source '%s' is newer than the existing store!  Can't load!\n", filename);
      return(false);
    }
  }
#endif

  //  Open the file, read the header.

  errno = 0;
  FILE *BLOCKS = fopen(blocksName, "r");
  if (errno)
    fprintf(stderr, "seqStream::loadStore()-- Failed to open blocks file '%s': %s\n",
            blocksName, strerror(errno)), exit(1);

  fread(cigam,             sizeof(char),   16, BLOCKS);
  fread(&_numBlocks,       sizeof(u64bit), 1,  BLOCKS);
  fread(&_numberOfACGT,    sizeof(u64bit), 1,  BLOCKS);

  if (errno)
    fprintf(stderr, "seqStream::loadStore()-- Failed to read header in seqStore '%s': %s\n",
            filename, strerror(errno)), exit(1);

  strcpy(magic, "seqStore.v1");
  if (strncmp(magic, cigam, 16) != 0)
    fprintf(stderr, "seqStream::loadStore()-- '%s' isn't a seqStream file.\n", filename), exit(1);

  //  The store construction adds a sentinal block at the end; that's the +1.
  _blockInfo = new seqStoreBlock [_numBlocks + 1];
  fread(_blockInfo, sizeof(seqStoreBlock), _numBlocks + 1, BLOCKS);

  if (errno)
    fprintf(stderr, "seqStream::loadStore()-- Failed to read "u32bitFMT" blocks: %s\n", _numBlocks, strerror(errno)), exit(1);

  fclose(BLOCKS);

  //  Position the bitPackedFile to the start of the merStream
  //
  _streamFile = new bitPackedFile(streamName, 0);

  _thisBlock         = 0;
  _thisBlockPosition = 0;

  return(true);
}



void
seqStore::buildStore(const char *filename, seqStream *ss) {
  char    magic[17]         = {0};
  u64bit  numBlocks         = 0;
  u64bit  numACGT           = 0;

  //  Open the output files
  //
  char *streamName = new char [strlen(filename) + 32];
  char *blocksName = new char [strlen(filename) + 32];

  sprintf(streamName, "%s.seqStore.sequence", filename);
  sprintf(blocksName, "%s.seqStore.blocks",   filename);
  
#if 1
  //  If 'filename' exists, assume it's the source fasta file.  We
  //  want to ensure that our seqStream is at least as current as
  //  that.
  if ((fileExists(streamName)) &&
      (fileExists(blocksName)) &&
      (fileExists(filename))) {
    u64bit  timeOfOriginal = timeOfFile(filename);
    if ((timeOfOriginal > timeOfFile(streamName)) ||
        (timeOfOriginal > timeOfFile(streamName))) {
      fprintf(stderr, "seqStore::buildStore()-- WARNING: source '%s' is newer than the existing store; rebuilding.\n", filename);
      unlink(streamName);
      unlink(blocksName);
    }
  }
#endif

  if ((fileExists(streamName)) ||
      (fileExists(blocksName))) {
    if (loadStore(filename))
      return;
    unlink(streamName);
    unlink(blocksName);
  }

  bitPackedFile *STREAM = new bitPackedFile(streamName);
  errno = 0;
  FILE          *BLOCKS = fopen(blocksName, "w+");
  if (errno)
    fprintf(stderr, "seqStore::buildStore()-- failed to open the output file '%s': %s\n",
            blocksName, strerror(errno)), exit(1);

  //             012345678901234567
  strcpy(magic, "seqStore.v1.PRTL");
  fwrite(magic,            sizeof(char),   16, BLOCKS);
  fwrite(&numBlocks,       sizeof(u64bit), 1,  BLOCKS);
  fwrite(&numACGT,         sizeof(u64bit), 1,  BLOCKS);

  delete [] blocksName;
  delete [] streamName;

  u64bit          thisBase = validCompressedSymbol[ss->get()];
  seqStoreBlock   b = {0};

  b._isACGT      = 0;
  b._seqIID      = ss->seqIID();
  b._posInSeq    = ss->seqPos();
  b._posInStrOff = 0;
  b._posInBPF    = STREAM->tell();
  b._len         = 0;

  //  Skip any crap at the start
  //
  while (thisBase == 0xff) {
    b._len++;
    thisBase = validCompressedSymbol[ss->get()];
  }

  if (b._len > 0) {
#ifdef DUMPBLOCKS
    fprintf(stderr, "initBlock: iid:"u64bitFMT" posInSeq:"u64bitFMT" posInBPF:"u64bitFMT" len:"u64bitFMT"\n",
            b._seqIID, b._posInSeq, b._posInBPF, b._len);
#endif
    fwrite(&b, sizeof(seqStoreBlock), 1, BLOCKS);
    numBlocks++;
  }

  b._isACGT      = 1;
  b._seqIID      = ss->seqIID();
  b._posInSeq    = ss->seqPos();
  b._posInStrOff = ss->strPos() - STREAM->tell() / 2;
  b._posInBPF    = STREAM->tell();
  b._len         = 0;

  //  Write the first letter.
  STREAM->putBits(thisBase, 2);
  b._len++;
  numACGT++;

  thisBase = validCompressedSymbol[ss->get()];

  while (ss->eof() == false) {

    //  Transition from ACGT to non-ACGT
    //
    if (((thisBase == 0xff) && (b._isACGT == 1)) ||
        ((thisBase != 0xff) && (b._isACGT == 0))) {

#ifdef DUMPBLOCKS
      fprintf(stderr, "block["u32bitFMT"]: iid:"u64bitFMT" posInSeq:"u64bitFMT" posInBPF:"u64bitFMT" len:"u64bitFMT"\n",
              numBlocks, b._seqIID, b._posInSeq, b._posInBPF, b._len);
#endif
      fwrite(&b, sizeof(seqStoreBlock), 1, BLOCKS);
      numBlocks++;

      b._isACGT      = !b._isACGT;
      b._seqIID      = ss->seqIID();
      b._posInSeq    = ss->seqPos();
      b._posInStrOff = ss->strPos() - STREAM->tell() / 2;
      b._posInBPF    = STREAM->tell();
      b._len         = 0;
    }

    if (thisBase != 0xff) {
      assert(b._isACGT == 1);
      STREAM->putBits(thisBase, 2);
      numACGT++;
    } else {
      assert(b._isACGT == 0);
    }

    b._len++;

    thisBase = validCompressedSymbol[ss->get()];
  };

  //  Be sure to write out the last block!
#ifdef DUMPBLOCKS
  fprintf(stderr, "finalBlock: iid:"u64bitFMT" posInSeq:"u64bitFMT" posInBPF:"u64bitFMT" len:"u64bitFMT"\n",
          b._seqIID, b._posInSeq, b._posInBPF, b._len);
#endif
  fwrite(&b, sizeof(seqStoreBlock), 1, BLOCKS);
  numBlocks++;


  //  And a sentinel for all done (actually, to just get the last
  //  position in the file).  We explicitly do not count this as a
  //  block.
  b._isACGT      = 0;
  b._seqIID      = ~u64bitZERO;
  b._posInSeq    = ~u64bitZERO;
  b._posInStrOff = ss->strPos() - STREAM->tell() / 2;
  b._posInBPF    = STREAM->tell();
  b._len         = 0;

#ifdef DUMPBLOCKS
  fprintf(stderr, "terminatingBlock: iid:"u64bitFMT" posInSeq:"u64bitFMT" posInBPF:"u64bitFMT" len:"u64bitFMT"\n",
          b._seqIID, b._posInSeq, b._posInBPF, b._len);
#endif
  fwrite(&b, sizeof(seqStoreBlock), 1, BLOCKS);


  //  Update the header on the blocks file
  //
  ::rewind(BLOCKS);

  //             012345678901234567
  memset(magic, 0, 16);
  strcpy(magic, "seqStore.v1");
  fwrite(magic,            sizeof(char),   16, BLOCKS);
  fwrite(&numBlocks,       sizeof(u64bit), 1,  BLOCKS);
  fwrite(&numACGT,         sizeof(u64bit), 1,  BLOCKS);

  fclose(BLOCKS);
  delete STREAM;

  loadStore(filename);
}




unsigned char
seqStore::get(void) {

  while ((_thisBlockPosition >= _blockInfo[_thisBlock]._len) &&
         (_thisBlock < _numBlocks)) {
    _thisBlock++;
    _thisBlockPosition = 0;

    if (_blockInfo[_thisBlock]._isACGT == 0)
      assert(_blockInfo[_thisBlock]._posInBPF == _streamFile->tell());

    //fprintf(stderr, "seqStore::get()-- NEW BLOCK thisBlock="u32bitFMT" numBlocks="u32bitFMT"\n", _thisBlock, _numBlocks);
  }

  if (_thisBlock >= _numBlocks) {
    //fprintf(stderr, "seqStore::get()-- EOF thisBlock="u32bitFMT" numBlocks="u32bitFMT"\n", _thisBlock, _numBlocks);
    return(0);
  }

  if (_blockInfo[_thisBlock]._isACGT == 0) {
    _thisBlockPosition++;
    //fprintf(stderr, "seqStore::get()-- N thisBlock="u32bitFMT" numBlocks="u32bitFMT"\n", _thisBlock, _numBlocks);
    return('n');
  }

  //  Ah, now in real sequence!

  _iteration++;

  if (_iteration >= _iterationLimit) {
    //fprintf(stderr, "seqStore::get()-- ITERS  thisBlock="u32bitFMT" numBlocks="u32bitFMT"\n", _thisBlock, _numBlocks);
    return(0);
  }

  _thisBlockPosition++;

  //char ch = decompressSymbol[_streamFile->getBits(2)];
  //fprintf(stderr, "seqStore::get()-- OK return %c  thisBlock="u32bitFMT" numBlocks="u32bitFMT"\n", ch, _thisBlock, _numBlocks);
  //return(ch);

  return(decompressSymbol[_streamFile->getBits(2)]);
}



bool
seqStore::seek(u64bit pos) {
  u32bit block = 0;

  //  The position should be an ACGT position (that is, 0 == the first
  //  ACGT in the file, 1 == the second, etc).  This is also exactly
  //  the same, just half, as the position in the BitPackedFile.
  //
  pos *= 2;

  //  Did someone request a mer number too high?
  //
  if (_blockInfo[_numBlocks]._posInBPF <= pos) {
    _thisBlock = _numBlocks;
    return(false);
  }

  //  Use a binary search if there are lots of blocks, otherwise, linear search.

  if (_numBlocks < 128) {
    while ((_blockInfo[block]._posInBPF <= pos) && (_blockInfo[block+1]._posInBPF < pos))
      block++;
  } else {
    u32bit  lo = 0;
    u32bit  hi = _numBlocks;
    u32bit  md = 0;

    while (lo <= hi) {
      md = (lo + hi) / 2;

      if        (pos < _blockInfo[md]._posInBPF) {
        //  This block starts after the mer we're looking for.  
        hi = md;
      } else if ((_blockInfo[md]._posInBPF <= pos) && (pos < _blockInfo[md+1]._posInBPF)) {
        //  Got it!
        lo           = md + 1;
        hi           = md;
        block        = md;
      } else {
        //  By default, then, the block is too low.
        lo = md;
      }
    }
  }

#if 0
  fprintf(stderr, "seqStore::seek("u64bitFMT" / 2) found block="u32bitFMT" posInBPF: "u64bitFMT" "u64bitFMT"\n",
          pos, block, _blockInfo[block]._posInBPF, _blockInfo[block+1]._posInBPF);
#endif

  assert((_blockInfo[block]._posInBPF <= pos) && (pos < _blockInfo[block+1]._posInBPF));

  _thisBlock         = block;
  _thisBlockPosition = (pos - _blockInfo[block]._posInBPF) / 2;

  _streamFile->seek(pos);

  return(true);
}

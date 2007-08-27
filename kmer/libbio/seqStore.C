
seqStore::seqStore() {
  _numBlocks        = 0;
  _blockInfo        = 0L;

  _thisBlock         = 0;
  _thisBlockPosition = 0;

  _numberOfACGT     = 0;

  _streamFile       = 0L;

  _iterationBeg     =  u64bitZERO;
  _iterationEnd     = ~u64bitZERO;
  _iteration        =  u64bitZERO;
}

seqStore::seqStore(const char *filename, seqStream *s) {
  _numBlocks        = 0;
  _blockInfo        = 0L;

  _thisBlock         = 0;
  _thisBlockPosition = 0;

  _numberOfACGT     = 0;

  _streamFile       = 0L;

  _iterationBeg     =  u64bitZERO;
  _iterationEnd     = ~u64bitZERO;
  _iteration        =  u64bitZERO;

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
seqStore::loadStore(const char *filename, seqStream *ss) {
  char    cigam[16] = {0};
  char    magic[16] = {0};

  if (ss)
    ss->rewind();

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
  //
  //  Meryl, in the guts of it, wants to open a new store, but doesn't
  //  have the original seqStream.  We allow that, just because we
  //  know that meryl checked before.
  //
  if ((ss) &&
      ((ss->timeStamp() > timeOfFile(streamName)) ||
       (ss->timeStamp() > timeOfFile(streamName)))) {
    fprintf(stderr, "seqStore::loadStore()-- WARNING: source '%s' is newer than the existing store!  Can't load!\n", filename);
    return(false);
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
  
  //  If 'filename' exists, assume it's the source fasta file.  We
  //  want to ensure that our seqStream is at least as current as
  //  that.
  if (fileExists(streamName) && fileExists(blocksName)) {
    u64bit  ts = 0;
    if (ss)
      ts = ss->timeStamp();
    if ((ts > timeOfFile(streamName)) || (ts > timeOfFile(streamName))) {
      fprintf(stderr, "seqStore::buildStore()-- WARNING: source '%s' is newer than the existing store; rebuilding.\n", filename);
      unlink(streamName);
      unlink(blocksName);
    } else if (loadStore(filename, ss)) {
      return;
    }
  }

  //  Otherwise, we either have no input files, or one, but not the
  //  other.  Delete whatever is there.
  unlink(streamName);
  unlink(blocksName);

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

  ss->rewind();

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

  loadStore(filename, ss);
}



#undef DEBUG_GET

unsigned char
seqStore::get(void) {

  while ((_thisBlockPosition >= _blockInfo[_thisBlock]._len) &&
         (_thisBlock < _numBlocks)) {
    _thisBlock++;
    _thisBlockPosition = 0;

    if (_blockInfo[_thisBlock]._isACGT)
      assert(_blockInfo[_thisBlock]._posInBPF == _streamFile->tell());
  }

  if (_thisBlock >= _numBlocks) {
#ifdef DEBUG_GET
    fprintf(stderr, "seqStore::get()-- ALL DONE; return 0\n");
#endif
    return(0);
  }

  if (_blockInfo[_thisBlock]._isACGT == 0) {
#ifdef DEBUG_GET
    fprintf(stderr, "seqStore::get()-- GAP BLOCK thisBlockPos "u64bitFMT" out of blockSize "u64bitFMT"\n",
            _thisBlockPosition, _blockInfo[_thisBlock]._len);
#endif
    _thisBlockPosition++;
    return('n');
  }

  //  Ah, now in real sequence!

  if (_iteration > _iterationEnd) {
#ifdef DEBUG_GET
    fprintf(stderr, "seqStore::get()-- ITERATION LIMIT; return 0\n");
#endif
    return(0);
  }

  _iteration++;

#ifdef DEBUG_GET
  char ret = decompressSymbol[_streamFile->getBits(2)];

  fprintf(stderr, "seqStore::get()-- ACGT %c thisBlockPos "u64bitFMT" out of blockSize "u64bitFMT"\n",
          ret, _thisBlockPosition, _blockInfo[_thisBlock]._len);

  _thisBlockPosition++;

  return(ret);
#else
  _thisBlockPosition++;
  return(decompressSymbol[_streamFile->getBits(2)]);
#endif
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

#if 0
  fprintf(stderr, "seqStore::seek()-- looking for posInBPF "u64bitFMT"\n", pos);
  for (u32bit i=0; i<_numBlocks; i++)
    fprintf(stderr, "block["u32bitFMT"]: acgt:"u64bitFMT" iid:"u64bitFMT" posInSeq:"u64bitFMT" posInStrOff:"u64bitFMT" posInBPF:"u64bitFMT" len:"u64bitFMT"\n",
            i,
            _blockInfo[i]._isACGT,
            _blockInfo[i]._seqIID,
            _blockInfo[i]._posInSeq,
            _blockInfo[i]._posInStrOff,
            _blockInfo[i]._posInBPF,
            _blockInfo[i]._len);
#endif


  if (_numBlocks < 128) {
    //  Blocks of non-ACGT have their position set to the same as the
    //  next block of ACGT -- we don't add anything to the file and so
    //  the position doesn't change.  If our position is after the
    //  current block, both conditions are satisfied.  If our position
    //  is the start of the block, we'll skip the gap block, because
    //  the position of the gap block and the next (acgt) block are
    //  the same as the position.  We'll then move to the next block,
    //  pos == posInBPF, but the next posInBPF is not <= pos and we
    //  stop.
    //
    while ((_blockInfo[block]._posInBPF <= pos) && (_blockInfo[block+1]._posInBPF <= pos))
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
  fprintf(stderr, "seqStore::seek()--  found block="u32bitFMT" posInBPF: "u64bitFMT" nextBlock: "u64bitFMT"\n",
          block, _blockInfo[block]._posInBPF, _blockInfo[block+1]._posInBPF);
#endif

  assert(_blockInfo[block]._isACGT);
  assert((_blockInfo[block]._posInBPF <= pos) && (pos < _blockInfo[block+1]._posInBPF));

  _thisBlock         = block;
  _thisBlockPosition = (pos - _blockInfo[block]._posInBPF) / 2;

  _streamFile->seek(pos);

  return(true);
}

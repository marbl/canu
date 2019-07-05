
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "bits.H"
#include "files.H"


stuffedBits::stuffedBits(uint64 nBits) {

  _dataBlockLenMax = nBits;

  _dataBlocksLen   = 1;
  _dataBlocksMax   = 64;

  _dataBlockBgn    = new uint64   [_dataBlocksMax];
  _dataBlockLen    = new uint64   [_dataBlocksMax];
  _dataBlocks      = new uint64 * [_dataBlocksMax];

  memset(_dataBlockBgn, 0, sizeof(uint64)   * _dataBlocksMax);
  memset(_dataBlockLen, 0, sizeof(uint64)   * _dataBlocksMax);
  memset(_dataBlocks,   0, sizeof(uint64 *) * _dataBlocksMax);

  _dataPos = 0;
  _data    = _dataBlocks[0] = new uint64 [_dataBlockLenMax / 64];

  memset(_data, 0, sizeof(uint64) * _dataBlockLenMax / 64);

  _dataBlockBgn[0] = 0;
  _dataBlockLen[0] = 0;

  _dataBlk = 0;
  _dataWrd = 0;
  _dataBit = 64;

  //  Set up the Fibonacci encoding table.
  //
  //  It takes 46 values to saturate a uint32 (fib[47] > uint32).
  //  It takes 92 values to saturate a uint64 (fib[93] > uint64).

  _fibData[0] = 1;
  _fibData[1] = 1;

  for (uint32 ii=2; ii<93; ii++)
    _fibData[ii] = _fibData[ii-1] + _fibData[ii-2];

  assert(_fibData[45] < _fibData[46]);   //  Fails if _fibData is 32-bit signed.
  assert(_fibData[46] < _fibData[47]);   //  Fails if _fibData is 32-bit unsigned.
  assert(_fibData[91] < _fibData[92]);   //  Fails if _fibData is 64-bit signed.
}



stuffedBits::stuffedBits(char const *inputName) {

  _dataBlockLenMax = 0;

  _dataBlocksLen   = 0;
  _dataBlocksMax   = 0;

  _dataBlockBgn    = NULL;
  _dataBlockLen    = NULL;
  _dataBlocks      = NULL;

  _dataPos = 0;
  _data    = NULL;

  FILE *inFile = AS_UTL_openInputFile(inputName);
  loadFromFile(inFile);
  AS_UTL_closeFile(inFile);

  _dataBlk = 0;
  _dataWrd = 0;
  _dataBit = 64;
};


stuffedBits::stuffedBits(FILE *inFile) {

  _dataBlockLenMax = 0;

  _dataBlocksLen   = 0;
  _dataBlocksMax   = 0;

  _dataBlockBgn    = NULL;
  _dataBlockLen    = NULL;
  _dataBlocks      = NULL;

  _dataPos = 0;
  _data    = NULL;

  loadFromFile(inFile);

  _dataBlk = 0;
  _dataWrd = 0;
  _dataBit = 64;
};


stuffedBits::stuffedBits(readBuffer *B) {

  _dataBlockLenMax = 0;

  _dataBlocksLen   = 0;
  _dataBlocksMax   = 0;

  _dataBlockBgn    = NULL;
  _dataBlockLen    = NULL;
  _dataBlocks      = NULL;

  _dataPos = 0;
  _data    = NULL;

  loadFromBuffer(B);

  _dataBlk = 0;
  _dataWrd = 0;
  _dataBit = 64;
};


#if 0
//  This is untested.
stuffedBits::stuffedBits(stuffedBits &that) {

  _dataBlockLenMax = that._dataBlockLenMax;

  _dataBlocksLen = that._dataBlocksLen;
  _dataBlocksMax = that._dataBlocksMax;

  _dataBlockBgn = new uint64   [_dataBlocksMax];
  _dataBlockLen = new uint64   [_dataBlocksMax];
  _dataBlocks   = new uint64 * [_dataBlocksMax];

  memcpy(_dataBlockBgn, that._dataBlockBgn, sizeof(uint64)   * _dataBlocksMax);
  memcpy(_dataBlockLen, that._dataBlockLen, sizeof(uint64)   * _dataBlocksMax);
  memset(_dataBlocks,   0,                  sizeof(uint64 *) * _dataBlocksMax);

  for (uint32 ii=0; ii<_dataBlocksMax; ii++) {
    if (that._dataBlocks[ii]) {
      _dataBlocks[ii] = new uint64 [_dataBlockLenMax / 64];
      memcpy(_dataBlocks[ii], that._dataBlocks[ii], sizeof(uint64) * _dataBlockLenMax / 64);
    }

    if (that._data == that._dataBlocks[ii])
      _data = _dataBlocks[ii];
  }

  _dataPos = that._dataPos;

  _dataBlk = that._dataBlk;
  _dataWrd = that._dataWrd;
  _dataBit = that._dataBit;

  memcpy(_fibData, that._fibData, 93 * sizeof(uint64));
}
#endif


stuffedBits::~stuffedBits() {

  for (uint32 ii=0; ii<_dataBlocksLen; ii++)
    delete [] _dataBlocks[ii];

  delete [] _dataBlockBgn;
  delete [] _dataBlockLen;
  delete [] _dataBlocks;
};



void
stuffedBits::dumpToBuffer(writeBuffer *B) {

  B->write(&_dataBlockLenMax, sizeof(uint64));
  B->write(&_dataBlocksLen,   sizeof(uint32));
  B->write(&_dataBlocksMax,   sizeof(uint32));
  B->write( _dataBlockBgn,    sizeof(uint64) * _dataBlocksLen);
  B->write( _dataBlockLen,    sizeof(uint64) * _dataBlocksLen);

  for (uint32 ii=0; ii<_dataBlocksLen; ii++) {
    uint64  nWordsToWrite = _dataBlockLen[ii] / 64 + (((_dataBlockLen[ii] % 64) == 0) ? 0 : 1);
    uint64  nWordsAllocd  = _dataBlockLenMax / 64;

    assert(nWordsToWrite <= nWordsAllocd);

    B->write(_dataBlocks[ii], sizeof(uint64) * nWordsToWrite);
  }
}



bool
stuffedBits::loadFromBuffer(readBuffer *B) {
  uint32   nLoad    = 0;
  uint64   inLenMax = 0;
  uint32   inLen    = 0;
  uint32   inMax    = 0;

  if (B == NULL)     //  No buffer,
    return(false);   //  no load.

  //  Try to load the new parameters into temporary storage, so we can
  //  compare against what have already allocated.

  nLoad += B->read(&inLenMax, sizeof(uint64));  //  Max length of each block.
  nLoad += B->read(&inLen,    sizeof(uint32));  //  Number of blocks stored.
  nLoad += B->read(&inMax,    sizeof(uint32));  //  Number of blocks allocated.

  if (nLoad != 3)
    return(false);

  //  If the input blocks are not the same size as the blocks we have, remove them.

  if (_dataBlockLenMax != inLenMax) {
    for (uint32 ii=0; ii<_dataBlocksLen; ii++)
      delete [] _dataBlocks[ii];

    for (uint32 ii=0; ii<_dataBlocksMax; ii++)
      _dataBlocks[ii] = NULL;

    _dataBlockLenMax = inLenMax;
  }

  //  If there are more blocks than we have space for, grab more space.  Bgn and Len can just be
  //  reallocated.  The pointers need to be extended (to preserve what's already in there).

  if (_dataBlocksMax < inLen) {
    delete [] _dataBlockBgn;
    delete [] _dataBlockLen;

    _dataBlockBgn  = new uint64 [inLen];
    _dataBlockLen  = new uint64 [inLen];

    resizeArray(_dataBlocks, _dataBlocksLen, _dataBlocksMax, inLen, resizeArray_copyData | resizeArray_clearNew);
  }

  //  Update the parameters.

  _dataBlocksLen = inLen;

  //  Load the data.

  B->read(_dataBlockBgn,  sizeof(uint64) * _dataBlocksLen);
  B->read(_dataBlockLen,  sizeof(uint64) * _dataBlocksLen);

  for (uint32 ii=0; ii<_dataBlocksLen; ii++) {
    uint64  nWordsToRead  = _dataBlockLen[ii] / 64 + (((_dataBlockLen[ii] % 64) == 0) ? 0 : 1);
    uint64  nWordsAllocd  = _dataBlockLenMax / 64;

    assert(nWordsToRead <= nWordsAllocd);

    if (_dataBlocks[ii] == NULL)
      _dataBlocks[ii] = new uint64 [nWordsAllocd];

    B->read(_dataBlocks[ii], sizeof(uint64) * nWordsToRead);

    memset(_dataBlocks[ii] + nWordsToRead, 0, sizeof(uint64) * (nWordsAllocd - nWordsToRead));
  }

  //  Set up the read/write head.

  _dataPos = 0;
  _data    = _dataBlocks[0];

  _dataBlk = 0;
  _dataWrd = 0;
  _dataBit = 64;

  return(true);
}



void
stuffedBits::dumpToFile(FILE *F) {

  writeToFile(_dataBlockLenMax, "dataBlockLenMax", F);
  writeToFile(_dataBlocksLen,   "dataBlocksLen",   F);
  writeToFile(_dataBlocksMax,   "dataBlocksMax",   F);
  writeToFile(_dataBlockBgn,    "dataBlockBgn",    _dataBlocksLen, F);
  writeToFile(_dataBlockLen,    "dataBlockLen",    _dataBlocksLen, F);

  for (uint32 ii=0; ii<_dataBlocksLen; ii++) {
    uint64  nWordsToWrite = _dataBlockLen[ii] / 64 + (((_dataBlockLen[ii] % 64) == 0) ? 0 : 1);
    uint64  nWordsAllocd  = _dataBlockLenMax / 64;

    assert(nWordsToWrite <= nWordsAllocd);

    writeToFile(_dataBlocks[ii], "dataBlocks", nWordsToWrite, F);
  }
}



bool
stuffedBits::loadFromFile(FILE *F) {
  uint32   nLoad    = 0;
  uint64   inLenMax = 0;
  uint32   inLen    = 0;
  uint32   inMax    = 0;

  if (F == NULL)     //  No file,
    return(false);   //  no load.

  //  Try to load the new parameters into temporary storage, so we can
  //  compare against what have already allocated.

  nLoad += ::loadFromFile(inLenMax, "dataBlockLenMax", F, false);  //  Max length of each block.
  nLoad += ::loadFromFile(inLen,    "dataBlocksLen",   F, false);  //  Number of blocks stored.
  nLoad += ::loadFromFile(inMax,    "dataBlocksMax",   F, false);  //  Number of blocks allocated.

  if (nLoad != 3)
    return(false);

  //  If the input blocks are not the same size as the blocks we have, remove them.

  if (_dataBlockLenMax != inLenMax) {
    for (uint32 ii=0; ii<_dataBlocksLen; ii++)
      delete [] _dataBlocks[ii];

    for (uint32 ii=0; ii<_dataBlocksMax; ii++)
      _dataBlocks[ii] = NULL;

    _dataBlockLenMax = inLenMax;
  }

  //  If there are more blocks than we have space for, grab more space.  Bgn and Len can just be
  //  reallocated.  The pointers need to be extended (to preserve what's already in there).

  if (_dataBlocksMax < inLen) {
    delete [] _dataBlockBgn;
    delete [] _dataBlockLen;

    _dataBlockBgn  = new uint64 [inLen];
    _dataBlockLen  = new uint64 [inLen];

    resizeArray(_dataBlocks, _dataBlocksLen, _dataBlocksMax, inLen, resizeArray_copyData | resizeArray_clearNew);
  }

  //  Update the parameters.

  _dataBlocksLen = inLen;

  //  Load the data.

  ::loadFromFile(_dataBlockBgn,  "dataBlockBgn", _dataBlocksLen, F);
  ::loadFromFile(_dataBlockLen,  "dataBlockLen", _dataBlocksLen, F);

  for (uint32 ii=0; ii<_dataBlocksLen; ii++) {
    uint64  nWordsToRead  = _dataBlockLen[ii] / 64 + (((_dataBlockLen[ii] % 64) == 0) ? 0 : 1);
    uint64  nWordsAllocd  = _dataBlockLenMax / 64;

    assert(nWordsToRead <= nWordsAllocd);

    if (_dataBlocks[ii] == NULL)
      _dataBlocks[ii] = new uint64 [nWordsAllocd];

    ::loadFromFile(_dataBlocks[ii], "dataBlocks", nWordsToRead, F);

    memset(_dataBlocks[ii] + nWordsToRead, 0, sizeof(uint64) * (nWordsAllocd - nWordsToRead));
  }

  //  Set up the read/write head.

  _dataPos = 0;
  _data    = _dataBlocks[0];

  _dataBlk = 0;
  _dataWrd = 0;
  _dataBit = 64;

  return(true);
}



//  Set the position of stuffedBits to 'position'.
//  Ensure that at least 'length' bits exist in the current block.
//
void
stuffedBits::setPosition(uint64 position, uint64 length) {

  _dataBlk = 0;

  while ((position < _dataBlockBgn[_dataBlk]) && (_dataBlk < _dataBlocksLen))
    _dataBlk++;

  assert(_dataBlk < _dataBlocksLen);  //  What to do if we seek to an uninitialized piece?

  _dataPos = position - _dataBlockBgn[_dataBlk];
  _data    = _dataBlocks[_dataBlk];

  _dataWrd =      _dataPos / 64;
  _dataBit = 64 - _dataPos % 64;
}


uint64
stuffedBits::getPosition(void) {
  return(_dataBlockBgn[_dataBlk] + _dataPos);
}


uint64
stuffedBits::getLength(void) {
  uint64 nBits = 0;

  for (uint32 ii=0; ii<_dataBlocksLen; ii++)
    nBits += _dataBlockLen[ii];

  return(nBits);
}



//  A special case of getBinary().
bool
stuffedBits::getBit(void) {

  updateBlk(1);

  bool value = saveRightBits(_data[_dataWrd] >> (_dataBit - 1), 1);

  _dataPos++;
  _dataBit--;

  updateBit();

  return(value);
};



bool
stuffedBits::testBit(void) {

  updateBlk(1);

  bool value = saveRightBits(_data[_dataWrd] >> (_dataBit - 1), 1);

  return(value);
};



//  A special case of setBinary().
void
stuffedBits::setBit(bool bit) {

  ensureSpace(1);

  uint64  mask = ((uint64)1 << (_dataBit - 1));

  if (bit)
    _data[_dataWrd] |=  mask;
  else
    _data[_dataWrd] &= ~mask;

  _dataPos++;
  _dataBit--;

  updateBit();
  updateLen();
};




uint64
stuffedBits::getUnary(void) {
  uint64  value = 0;
  uint64  wrd;

  //  Ensure we're in valid data.

  updateBlk(1);

  //  Word align us first.

  wrd = _data[_dataWrd] << (64 - _dataBit);

  //  Skip entire words.  For the first word, if we're left with only zeros
  //  after the shifting, we increase the 'value' by the number of bits
  //  we could read in the word, not the full 64 bits that are zero.

  while (wrd == 0) {
    value    += _dataBit;

    _dataPos += _dataBit;
    _dataWrd += 1;

    _dataBit  = 64;
    wrd       = _data[_dataWrd];
  }

  //  Decode the last partial word.

  wrd       = 64 - countNumberOfBits64(wrd);

  value    += wrd;

  _dataPos += wrd + 1;
  _dataBit -= wrd + 1;

  updateBit();

  return(value);
}

uint64 *
stuffedBits::getUnary(uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getUnary();

  return(values);
}




uint32
stuffedBits::setUnary(uint64 value) {

  ensureSpace(value+1);

  //  If we fit entirely within this word, handle it special.

  if (value + 1 < _dataBit) {
    _data[_dataWrd] = clearMiddleBits(_data[_dataWrd], 64 - _dataBit, _dataBit + value + 1);

    _dataPos += value + 1;
    _dataBit -= value + 1;

    _data[_dataWrd] |= (uint64)1 << _dataBit;

    updateLen();

    return(value + 1);
  }

  //  We fit _exactly_ in this word, special again!

  if (value + 1 == _dataBit) {
    _data[_dataWrd]  = clearRightBits(_data[_dataWrd], _dataBit);
    _data[_dataWrd] |= 1;   //  ALWAYS the last bit.

    _dataPos += value + 1;  //  ALWAYS move to the next word.
    _dataWrd += 1;
    _dataBit  = 64;

    updateLen();

    return(value + 1);
  }

  //  Otherwise, we span at least two words.  First, get us word aligned,
  //  by clearing the rest of this word.

  assert(value + 1 > _dataBit);

  _data[_dataWrd] = clearRightBits(_data[_dataWrd], _dataBit);

  value    -= _dataBit;

  _dataPos += _dataBit;
  _dataWrd += 1;
  _dataBit  = 64;

  assert((_dataPos % 64) == 0);

  //  Then, set entire words to zero.

  while (value >= 64) {
    _data[_dataWrd] = 0;

    value    -= 64;

    _dataPos += 64;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Finally, set the last partial word.  This is similar to the cases
  //  at the start, but we know the bits will always be on the left of the word.

  _data[_dataWrd] = clearLeftBits(_data[_dataWrd], value + 1);

  _dataPos += value + 1;                      //  Skip the zero bits.
  _dataBit -= value + 1;                      //  (and the sentinel)

  _data[_dataWrd] |= (uint64)1 << _dataBit;   //  Add the sentinel.

  //  Update the pointers.

  updateLen();
  updateBit();

  return(value + 1);
}



uint32
stuffedBits::setUnary(uint64 number, uint64 *values) {
  uint32  size = 0;

  for (uint64 ii=0; ii<number; ii++)
    size += setUnary(values[ii]);

  return(size);
}





uint64
stuffedBits::getBinary(uint32 width) {
  uint64  value = 0;

  assert(width < 65);

  //  Ensure we're in valid data.

  updateBlk(width);

  //  If we're contained in a single word, special case.

  if (width < _dataBit) {
    value = saveRightBits(_data[_dataWrd] >> (_dataBit - width), width);

    _dataPos += width;
    _dataBit -= width;
  }

  //  If we're exactly in a single word, another special case.

  else if (width == _dataBit) {
    value = saveRightBits(_data[_dataWrd], width);

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Otherwise, we're spanning two words.

  else {
    uint64  w1 =         _dataBit;
    uint64  w2 = width - _dataBit;

    uint64  l  = saveRightBits(_data[_dataWrd + 0], w1) <<       w2;
    uint64  r  = saveLeftBits (_data[_dataWrd + 1], w2) >> (64 - w2);

    value = l | r;

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64 - w2;
  }

  return(value);
}



uint64 *
stuffedBits::getBinary(uint32 width, uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getBinary(width);

  return(values);
}




uint32
stuffedBits::setBinary(uint32 width, uint64 value) {

  assert(width < 65);

  ensureSpace(width);

  //  Mask off any pieces we're not supposed to be seeing.

  value = saveRightBits(value, width);

  //  If we fit entirely within this word, handle it special.

  if (width < _dataBit) {
    _data[_dataWrd] = clearMiddleBits(_data[_dataWrd], 64 - _dataBit, _dataBit + width) | (value << (_dataBit - width));

    _dataPos += width;
    _dataBit -= width;
  } else

  //  We fit _exactly_ in this word, special again!

  if (width == _dataBit) {
    _data[_dataWrd] = clearRightBits(_data[_dataWrd], _dataBit) | value;

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Otherwise, we span two words.

  else {
    uint32  w1 =         _dataBit;
    uint32  w2 = width - _dataBit;

    _data[_dataWrd + 0] = clearRightBits(_data[_dataWrd + 0], w1) | (value >> (     w2));
    _data[_dataWrd + 1] = clearLeftBits (_data[_dataWrd + 1], w2) | (value << (64 - w2));

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64 - w2;
  }

  //  updateBit() isn't needed; it's handled in the special cases.

  updateLen();

  return(width);
}



uint32
stuffedBits::setBinary(uint32 width, uint64 number, uint64 *values) {
  uint32 size = 0;

  for (uint64 ii=0; ii<number; ii++)
    size += setBinary(width, values[ii]);

  return(size);
}





////////////////////////////////////////
//  ELIAS GAMMA CODED DATA
//
//  Unary coded length of binary data, then binary data of that length.
//  Works only on positive (non-zero) integers.
//
uint64
stuffedBits::getEliasGamma(void) {
  uint32  N = getUnary();
  uint64  V = getBinary(N);

  V |= (uint64)1 << N;

  return(V);
}



uint64 *
stuffedBits::getEliasGamma(uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getEliasGamma();

  return(values);
}



uint32
stuffedBits::setEliasGamma(uint64 value) {
  uint32 size = 0;
  uint64 N    = countNumberOfBits64(value) - 1;

  assert(value > 0);

  size += setUnary(N);
  size += setBinary(N, value);

  return(size);
}



uint32
stuffedBits::setEliasGamma(uint64 number, uint64 *values) {
  uint32  size = 0;

  for (uint64 ii=0; ii<number; ii++)
    size += setEliasGamma(values[ii]);

  return(size);
}



////////////////////////////////////////
//  ELIAS DELTA CODED DATA
//
//  Similar to the gamma code, except the number of bits itself
//  is gamma coded.  An optimization can drop the high order bit (it's always 1)
//  from the binary coded data.  We don't do that.
//
uint64
stuffedBits::getEliasDelta(void) {
  uint32  N = getEliasGamma() - 1;
  uint64  V = getBinary(N);

  V |= (uint64)1 << N;

  return(V);
}



uint64 *
stuffedBits::getEliasDelta(uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getEliasDelta();

  return(values);
}




uint32
stuffedBits::setEliasDelta(uint64 value) {
  uint32 size = 0;
  uint32 N    = countNumberOfBits64(value);

  assert(value > 0);

  size += setEliasGamma(N);
  size += setBinary(N-1, value);

  return(size);
}



uint32
stuffedBits::setEliasDelta(uint64 number, uint64 *values) {
  uint32  size = 0;

  for (uint64 ii=0; ii<number; ii++)
    size += setEliasDelta(values[ii]);

  return(size);
}





////////////////////////////////////////
//  FIBONACCI CODED DATA


uint64
stuffedBits::getZeckendorf(void) {
  uint64  value = 0;
  uint32  ff    = 1;

  //  The first bit in the official representation, representing the
  //  redundant value 1, is always zero, and we don't save it.  Thus, start
  //  decoding at ff=1.

  bool tbit = getBit();
  bool nbit = getBit();

  //fprintf(stderr, "getZeck() %d %d blk %2lu pos %2lu bit %2lu\n", tbit, nbit, _dataBlk, _dataPos, _dataBit);

  while (true) {
    if (tbit)
      value += _fibData[ff];

    if (tbit && nbit)
      break;

    ff++;

    tbit = nbit;
    nbit = getBit();

    //fprintf(stderr, "getZeck() %d   blk %2lu pos %2lu bit %2lu\n", nbit, _dataBlk, _dataPos, _dataBit);
  }

  return(value);
}



uint64 *
stuffedBits::getZeckendorf(uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getZeckendorf();

  return(values);
}




uint32
stuffedBits::setZeckendorf(uint64 value) {
  uint32  ff = 0;

  uint64  word1 = 0;   uint32  wlen1 = 0;
  uint64  word2 = 0;   uint32  wlen2 = 0;

  //  Find the largest Fibonacci number smaller than our value.
  //  Probably should be binary searching for this.

  //fprintf(stderr, "setZeckendorf() - value %lu\n", value);

  while ((ff < 93) && (_fibData[ff] <= value))
    ff++;

  //fprintf(stderr, "setZeckendorf() - ff    %lu\n", ff);

  //  For each smaller Fibonacci number:
  //    If the Fibonacci number is more than the value, it's not used in the
  //    encoding.  Push on a zero.
  //
  //    Otherwise, it is used in the encoding.  Push on a 1, and remove the
  //    fib number from our value.
  //
  while (ff-- > 0) {
    word2 <<= 1;                       //  Make space for the new bit.

    if (_fibData[ff] <= value) {      //  If used in the encoding,
      word2 |= 1;                     //  set the bit and remove
      value -= _fibData[ff];          //  it from the value.
    }

    if (++wlen2 > 60) {               //  If we're running outta
      word1 = word2;                  //  bits in the word, save it
      wlen1 = wlen2;                  //  to the first word to output.
      wlen2 = 0;                      //
    }
  }

  //  Reverse the words so we see the low bit first, then push on a
  //  terminating 1 so we end the string with a pair of 1 bits.
  //
  //  The lower bits, in word2, can have the (post-reverse) left-most bit,
  //  representing the redundant 1, stripped off.
  //
  //  An annoying special case oocurs when there are exactly 60 bits in the
  //  encoding: word2 is now empty!

  if (wlen1 == 0) {
    word2 = reverseBits64(word2);

    word2 >>= (64 - wlen2);

    word2 <<= 1;
    word2  |= 1;
    wlen2  += 1;

    wlen2  -= 1;  //  Strip off left-most bit.  Go Optimizer, Go!

    setBinary(wlen2, word2);

    //fprintf(stderr, "setZeckendorf() - word2 0x%016lx %2u\n", word2, wlen2);
  }

  else if (wlen2 == 0) {
    word1 = reverseBits64(word1);

    word1 >>= (64 - wlen1);

    word1 <<= 1;
    word1  |= 1;
    wlen1  += 1;

    wlen1  -= 1;  //  Strip off left-most bit.  Go Optimizer, Go!

    setBinary(wlen1, word1);

    //fprintf(stderr, "setZeckendorf() - word1 0x%016lx %1u\n", word1, wlen1);
  }

  else {
    word2 = reverseBits64(word2);
    word1 = reverseBits64(word1);

    word2 >>= (64 - wlen2);
    word1 >>= (64 - wlen1);

    word1 <<= 1;
    word1  |= 1;
    wlen1  += 1;

    wlen2  -= 1;

    setBinary(wlen2, word2);
    setBinary(wlen1, word1);

    //fprintf(stderr, "setZeckendorf() - word2 0x%016lx %2u word1 0x%016lx %2u\n", word2, wlen2, word1, wlen1);
  }


  return(wlen1 + wlen2);
}



uint32
stuffedBits::setZeckendorf(uint64 number, uint64 *values) {
  uint32  size = 0;

  for (uint64 ii=0; ii<number; ii++)
    size += setZeckendorf(values[ii]);

  return(size);
}


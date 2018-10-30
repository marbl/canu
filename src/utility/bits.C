
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
};


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


stuffedBits::~stuffedBits() {

  for (uint32 ii=0; ii<_dataBlocksLen; ii++)
    delete [] _dataBlocks[ii];

  delete [] _dataBlockBgn;
  delete [] _dataBlockLen;
  delete [] _dataBlocks;
};



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
stuffedBits::getGeneralizedUnary(void) {
  return(0);
}



uint64 *
stuffedBits::getGeneralizedUnary(uint64 number, uint64 *values) {
  return(values);
}




uint32
stuffedBits::setGeneralizedUnary(uint64 value) {
  return(0);
}



uint32
stuffedBits::setGeneralizedUnary(uint64 number, uint64 *values) {
  return(0);
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




void
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
}



void
stuffedBits::setBinary(uint32 width, uint64 number, uint64 *values) {

  for (uint64 ii=0; ii<number; ii++)
    setBinary(width, values[ii]);
}




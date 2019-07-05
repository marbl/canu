
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

#ifndef LIBBITS_H
#define LIBBITS_H

#include "AS_global.H"

#include "files.H"

#include <algorithm>

//  Define this to enable testing that the width of the data element is greater than zero.  The
//  uint64MASK() macro (bri.h) does not generate a mask for 0.  Compiler warnings are issued,
//  because you shouldn't use this in production code.
//
//  As it's expensive, emit a warning if it's enabled.
//
//#define CHECK_WIDTH
//#warning CHECK_WIDTH is EXPENSIVE


//  Writing in the middle of data is toublesome.
//  This code currently will not split an object across two allocated blocks,
//  so if you want to rewrite in the middle, you need to make sure it
//  doesn't span a previously written block.  E.g., two writes
//  of 32 bits each could be in different blocks, and trying
//  to overwrite with a single 64-bit write could result in the first
//  block getting truncated (to the current position) and the write
//  replacing the first 64-bits in the second block, not just the 32-bits expected.
//
//  And don't even think of overwriting any of the variable length data.


inline
char *
displayWord(uint64 word, char *buffer=NULL) {
  static char b[65];

  if (buffer == NULL)
    buffer = b;

  memset(buffer, 'x', 64);

  for (uint32 ii=0; ii<64; ii++)
    buffer[ii] = (word & (uint64)1 << (63 - ii)) ? '1' : '0';

  buffer[64] = 0;

  return(buffer);
};


//  Return bits in a word:
//    Keeping the rightmost 64-n bits (mask out the leftmost  n bits)
//    Keeping the leftmost  64-n bits (mask out the rigthmost n bits)
//
inline uint64   clearLeftBits  (uint64 v, uint32 l)  {  if (l >= 64) return(0);  return(v & (0xffffffffffffffffllu >>       l));   };
inline uint64    saveLeftBits  (uint64 v, uint32 l)  {  if (l ==  0) return(0);  return(v & (0xffffffffffffffffllu << (64 - l)));  };
inline uint64   clearRightBits (uint64 v, uint32 r)  {  if (r >= 64) return(0);  return(v & (0xffffffffffffffffllu <<       r));   };
inline uint64    saveRightBits (uint64 v, uint32 r)  {  if (r ==  0) return(0);  return(v & (0xffffffffffffffffllu >> (64 - r)));  };

inline uint64   clearMiddleBits(uint64 v, uint32 l, uint32 r) {
  return(saveRightBits(v, r) | saveLeftBits(v, l));
};

inline uint64    saveMiddleBits(uint64 v, uint32 l, uint32 r) {
  return(clearRightBits(v, r) & clearLeftBits(v, l));
};



//  Freed, Edwin E. 1983. "Binary Magic Number" Dr. Dobbs Journal Vol. 78 (April) pp. 24-37
//    Reverse the bits in a word,
//    Count the number of set bits in a words
//    and more.
//
inline
uint64
reverseBits64(uint64 x) {
  x = ((x >>  1) & uint64NUMBER(0x5555555555555555)) | ((x <<  1) & uint64NUMBER(0xaaaaaaaaaaaaaaaa));
  x = ((x >>  2) & uint64NUMBER(0x3333333333333333)) | ((x <<  2) & uint64NUMBER(0xcccccccccccccccc));
  x = ((x >>  4) & uint64NUMBER(0x0f0f0f0f0f0f0f0f)) | ((x <<  4) & uint64NUMBER(0xf0f0f0f0f0f0f0f0));
  x = ((x >>  8) & uint64NUMBER(0x00ff00ff00ff00ff)) | ((x <<  8) & uint64NUMBER(0xff00ff00ff00ff00));
  x = ((x >> 16) & uint64NUMBER(0x0000ffff0000ffff)) | ((x << 16) & uint64NUMBER(0xffff0000ffff0000));
  x = ((x >> 32) & uint64NUMBER(0x00000000ffffffff)) | ((x << 32) & uint64NUMBER(0xffffffff00000000));
  return(x);
}

inline
uint32
reverseBits32(uint32 x) {
  x = ((x >>  1) & uint32NUMBER(0x55555555)) | ((x <<  1) & uint32NUMBER(0xaaaaaaaa));
  x = ((x >>  2) & uint32NUMBER(0x33333333)) | ((x <<  2) & uint32NUMBER(0xcccccccc));
  x = ((x >>  4) & uint32NUMBER(0x0f0f0f0f)) | ((x <<  4) & uint32NUMBER(0xf0f0f0f0));
  x = ((x >>  8) & uint32NUMBER(0x00ff00ff)) | ((x <<  8) & uint32NUMBER(0xff00ff00));
  x = ((x >> 16) & uint32NUMBER(0x0000ffff)) | ((x << 16) & uint32NUMBER(0xffff0000));
  return(x);
}


inline
uint64
uint64Swap(uint64 x) {
  x = ((x >>  8) & uint64NUMBER(0x00ff00ff00ff00ff)) | ((x <<  8) & uint64NUMBER(0xff00ff00ff00ff00));
  x = ((x >> 16) & uint64NUMBER(0x0000ffff0000ffff)) | ((x << 16) & uint64NUMBER(0xffff0000ffff0000));
  x = ((x >> 32) & uint64NUMBER(0x00000000ffffffff)) | ((x << 32) & uint64NUMBER(0xffffffff00000000));
  return(x);
}

inline
uint32
uint32Swap(uint32 x) {
  x = ((x >>  8) & uint32NUMBER(0x00ff00ff)) | ((x <<  8) & uint32NUMBER(0xff00ff00));
  x = ((x >> 16) & uint32NUMBER(0x0000ffff)) | ((x << 16) & uint32NUMBER(0xffff0000));
  return(x);
}

inline
uint16
uint16Swap(uint16 x) {
  x = ((x >>  8) & 0x00ff) | ((x <<  8) & 0xff00);
  return(x);
}


inline
uint32
countNumberOfSetBits32(uint32 x) {
  x = ((x >>  1) & uint32NUMBER(0x55555555)) + (x & uint32NUMBER(0x55555555));
  x = ((x >>  2) & uint32NUMBER(0x33333333)) + (x & uint32NUMBER(0x33333333));
  x = ((x >>  4) & uint32NUMBER(0x0f0f0f0f)) + (x & uint32NUMBER(0x0f0f0f0f));
  x = ((x >>  8) & uint32NUMBER(0x00ff00ff)) + (x & uint32NUMBER(0x00ff00ff));
  x = ((x >> 16) & uint32NUMBER(0x0000ffff)) + (x & uint32NUMBER(0x0000ffff));
  return(x);
}

inline
uint64
countNumberOfSetBits64(uint64 x) {
  x = ((x >>  1) & uint64NUMBER(0x5555555555555555)) + (x & uint64NUMBER(0x5555555555555555));
  x = ((x >>  2) & uint64NUMBER(0x3333333333333333)) + (x & uint64NUMBER(0x3333333333333333));
  x = ((x >>  4) & uint64NUMBER(0x0f0f0f0f0f0f0f0f)) + (x & uint64NUMBER(0x0f0f0f0f0f0f0f0f));
  x = ((x >>  8) & uint64NUMBER(0x00ff00ff00ff00ff)) + (x & uint64NUMBER(0x00ff00ff00ff00ff));
  x = ((x >> 16) & uint64NUMBER(0x0000ffff0000ffff)) + (x & uint64NUMBER(0x0000ffff0000ffff));
  x = ((x >> 32) & uint64NUMBER(0x00000000ffffffff)) + (x & uint64NUMBER(0x00000000ffffffff));
  return(x);
}


//  Return the number of bits needed to represent 'x'.
//  It's really floor(log_2(x)) + 1.
//  Note that x=0 returns 0.
//
inline
uint32
countNumberOfBits32(uint32 x) {
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return(countNumberOfSetBits32(x));
}

inline
uint64
countNumberOfBits64(uint64 x) {
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  return(countNumberOfSetBits64(x));
}




//#if (__GNUC__ > 3) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
//#define BUILTIN_POPCOUNT
//#endif

#ifdef BUILTIN_POPCOUNT

inline
uint32
countNumberOfSetBits32(uint32 x) {
  return(__builtin_popcount(x));
}

inline
uint64
countNumberOfSetBits64(uint64 x) {
  return(__builtin_popcountll(x));
}

#endif





class bitArray {
public:
  bitArray(uint64 maxNumBits=0) {
    _maxBitSet   = 0;
    _maxBitAvail = maxNumBits;
    _bits        = NULL;

    if (maxNumBits > 0)
      allocate(maxNumBits);
  };

  ~bitArray(void) {
    delete [] _bits;
  };

  bool     isAllocated(void) {
    return(_bits != NULL);
  }

  void     allocate(uint64 maxNumBits) {
    if (isAllocated() == true)
      return;

    _maxBitSet = 0;
    _maxBitAvail = maxNumBits;
    _bits        = new uint64 [_maxBitAvail / 64 + 1];

    clear();
  };

  void     clear(void) {
    memset(_bits, 0, sizeof(uint64) * (_maxBitAvail / 64 + 1));
  };

  bool     getBit(uint64 position) {
    uint64   w =      (position / 64);
    uint64   b = 63 - (position % 64);

    if (_maxBitAvail <= position)
      fprintf(stderr, "getBit()--  ERROR: position=" F_U64 " > maximum available=" F_U64 "\n",
              position, _maxBitAvail);
    assert(position < _maxBitAvail);

    return((_bits[w] >> b) & 0x00000001);
  };

  void     setBit(uint64 position, bool value) {
    uint64   w =      (position / 64);
    uint64   b = 63 - (position % 64);
    uint64   m = ((uint64)1) << b;

    //fprintf(stderr, "SET pos %9" F_U64P " word %2" F_U64P " bit %2" F_U64P " value %c  0x%016" F_X64P " -> ",
    //        position, w, b, (value) ? '1' : '0', _bits[w]);

    if (_maxBitAvail <= position)
      fprintf(stderr, "setBit()--  ERROR: position=" F_U64 " > maximum available=" F_U64 "\n",
              position, _maxBitAvail);
    assert(position < _maxBitAvail);

    _bits[w] &= ~m;
    _bits[w] |= ((uint64)value) << b;

    //fprintf(stderr, "0x%016" F_X64P "\n", _bits[w]);
  };

  bool     flipBit(uint64 position) {
    uint64   w =      (position / 64);
    uint64   b = 63 - (position % 64);
    uint64   m = ((uint64)1) << b;

    if (_maxBitAvail <= position)
      fprintf(stderr, "flipBit()--  ERROR: position=" F_U64 " > maximum available=" F_U64 "\n",
              position, _maxBitAvail);
    assert(position < _maxBitAvail);

    uint64   v = _bits[w] & m;

    //fprintf(stderr, "FLIP w %lu b %lu m 0x%016lx v 0x%016lx FROM 0x%016lx", w, b, m, v, _bits[w]);

    _bits[w] ^= m;

    //fprintf(stderr, " TO 0x%016lx\n", _bits[w]);

    return(v >> b);
  };

private:
  uint64   _maxBitSet;
  uint64   _maxBitAvail;
  uint64  *_bits;
};



////////////////////////////////////////
//
//  wordArray
//
//  An array that efficiently stores non-machine-word size integer words by
//  packing the bits into machine-size words.
//
//  The array is variable length, but not sparse.  Accessing element
//  1,000,000 will allocate elements 0 through 999,999.
//
//  No array operator can be provided since we cannot return a reference to
//  values across machine words, let alone a reference to a value inside a
//  machine word.
//
//  The constructor needs to know the size of the words being stored,
//  and how many bits to store per allocation.
//
class wordArray {
public:
  wordArray(uint32 wordWidth, uint32 segmentSize = 65536 * 8) {
    _valueWidth       = wordWidth;
    _segmentSize      = segmentSize;
    _valuesPerSegment = (uint64)_segmentSize / (uint64)_valueWidth;

    _nextElement      = 0;

    _segmentsLen      = 0;
    _segmentsMax      = 16;
    _segments         = new uint64 * [_segmentsMax];

    for (uint32 ss=0; ss<_segmentsMax; ss++)
      _segments[ss] = NULL;
  }

  ~wordArray() {
    for (uint32 i=0; i<_segmentsLen; i++)
      delete [] _segments[i];

    delete [] _segments;
  };

  void     clear(void) {
    _nextElement = 0;
    _segmentsLen = 0;
  };

  void     allocate(uint64 nElements) {
    uint64 nSegs = nElements / _valuesPerSegment + 1;

    //fprintf(stderr, "wordArray::allocate()-- allocating space for " F_U64 " elements, in " F_U64 " segments.\n",
    //        nElements, nSegs);

    assert(_segmentsLen == 0);

    resizeArray(_segments, _segmentsLen, _segmentsMax, nSegs, resizeArray_copyData | resizeArray_clearNew);

    for (uint32 seg=0; seg<nSegs; seg++) {
      if (_segments[seg] == NULL)
        _segments[seg] = new uint64 [_segmentSize / 64];

      memset(_segments[seg], 0xff, sizeof(uint64) * _segmentSize / 64);
    }

    _segmentsLen = nSegs;
  };

  uint64   get(uint64 element) {
    uint64 seg =                element / _valuesPerSegment;     //  Which segment are we in?
    uint64 pos = _valueWidth * (element % _valuesPerSegment);    //  Bit position of the start of the value.

    uint64 wrd = pos / 64;   //  The word we start in.
    uint64 bit = pos % 64;   //  Starting at this bit.

    uint64 val = 0;

    assert(element < _nextElement);

    //  If the value is all in one word, just shift that word to the right to
    //  put the proper bits in the proper position.  We'll clean up the extra
    //  bits in just a moment.

    if (bit + _valueWidth <= 64) {
      val   = _segments[seg][wrd] >> (64 - _valueWidth - bit);
    }

    //  Otherwise, the value spans two words.  First, shift the first word so
    //  the end of it is at the start of the value.  Then shift the second
    //  word to the start of it is at the end of the value.
    //
    //                                 ssssssssssssssssssssss <- second shift
    //  [--word--][--first-word--][--second-word--][--word--]
    //                      [--value--]
    //                            fffff <- first shift
    //
    else {
      uint32  fShift = _valueWidth - (64 - bit);
      uint32  sShift = 64 - fShift;

      val   = _segments[seg][wrd+0] << fShift;
      val  |= _segments[seg][wrd+1] >> sShift;
    }

    //  Finally, mask off the stuff we don't care about.

    val  &= uint64MASK(_valueWidth);

    return(val);
  };

  void     set(uint64 element, uint64 value) {
    uint64 seg =                element / _valuesPerSegment;     //  Which segment are we in?
    uint64 pos = _valueWidth * (element % _valuesPerSegment);    //  Which word in the segment?

    uint64 wrd = pos / 64;   //  The word we start in.
    uint64 bit = pos % 64;   //  Starting at this bit.

    if (element >= _nextElement)
      _nextElement = element+1;

    if (seg >= _segmentsMax)
      resizeArray(_segments, _segmentsLen, _segmentsMax, seg + 16, resizeArray_copyData | resizeArray_clearNew);

    while (_segmentsLen <= seg) {
      _segments[_segmentsLen] = new uint64 [_segmentSize / 64];

      memset(_segments[_segmentsLen], 0xff, sizeof(uint64) * _segmentSize / 64);

      _segmentsLen++;
    }

    //  Mask the value, just in case.

    value &= uint64MASK(_valueWidth);

    //  Set the value in the segment.

    //          [--------------------]
    //                 [value]
    //           lSave           rSave
    //
    if (bit + _valueWidth <= 64) {
      uint32   lSave = bit;
      uint32   rSave = 64 - _valueWidth - bit;

      _segments[seg][wrd] = (saveLeftBits(_segments[seg][wrd], lSave) |
                             (value << rSave)                         |
                             saveRightBits(_segments[seg][wrd], rSave));
    }

    //            --lSave--               --rSave--
    //  [--word--][--first-word--][--second-word--][--word--]
    //                     [----value---=]
    //                      lSize  rSize
    //
    else {
      uint32   lSave =      bit,   rSave = 128 - _valueWidth - bit;
      uint32   lSize = 64 - bit,   rSize = _valueWidth - (64 - bit);

      _segments[seg][wrd+0] = saveLeftBits(_segments[seg][wrd+0], lSave) | (value >> rSize);
      _segments[seg][wrd+1] = (value << rSave) | saveRightBits(_segments[seg][wrd+1], rSave);
    }
  };

  void     show(void) {
    fprintf(stderr, "wordArray:   valueWidth  %2" F_U32P "\n", _valueWidth);
    fprintf(stderr, "wordArray:   segmentSize %8" F_U64P "   valuesPerSegment %8" F_U64P "\n", _segmentSize, _valuesPerSegment);
    fprintf(stderr, "\n");

    uint32  bit  = 64;
    uint32  word = 0;
    char    bits[65];

    for (uint32 ss=0; ss<_segmentsLen; ss++) {
      fprintf(stderr, "Segment %u:\n", ss);

      for(uint32 wrd=0, bit=0; bit<_valuesPerSegment * _valueWidth; bit++) {
        if ((bit % 64) == 0) {
          displayWord(_segments[ss][wrd++], bits);
        }

        if ((bit % _valueWidth) == 0)
          fprintf(stderr, "word %2u: ", wrd);

        fprintf(stderr, "%c", bits[bit % 64]);

        if ((bit % _valueWidth) == _valueWidth - 1)
          fprintf(stderr, "\n");
      }
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
  }

private:
  uint32   _valueWidth;
  uint64   _segmentSize;
  uint64   _valuesPerSegment;

  uint64   _nextElement;  //  the first invalid element

  uint64   _segmentsLen;
  uint64   _segmentsMax;
  uint64 **_segments;
};



class stuffedBits {
public:
  stuffedBits(uint64 nBits=16 * 1024 * 1024 * 8);
  stuffedBits(const char *inputName);
  stuffedBits(FILE *inFile);
  stuffedBits(readBuffer *B);
  //stuffedBits(stuffedBits &that);   //  Untested.
  ~stuffedBits();

  //  Debugging.

  char    *displayWord(uint64 w) {
    return(::displayWord(_data[w]));
  };

  //  Files.

  void     dumpToBuffer(writeBuffer *B);
  bool     loadFromBuffer(readBuffer *B);

  void     dumpToFile(FILE *F);
  bool     loadFromFile(FILE *F);

  //  Management of the read/write head.

  void     setPosition(uint64 position, uint64 length = 0);
  uint64   getPosition(void);

  uint64   getLength(void);

  void     byteAlign(void);

  //  SINGLE BITS

  bool     getBit(void);           //  get a bit.
  bool     testBit(void);          //  get a bit, but don't move position.
  void     setBit(bool on=true);   //  set a bit.

  //  UNARY CODED DATA

  uint64   getUnary(void);
  uint64  *getUnary(uint64 number, uint64 *values);

  uint32   setUnary(uint64 value);
  uint32   setUnary(uint64 number, uint64 *values);

  //  BINARY CODED DATA

  uint64   getBinary(uint32 width);
  uint64  *getBinary(uint32 width, uint64 number, uint64 *values=NULL);

  uint32   setBinary(uint32 width, uint64 value);
  uint32   setBinary(uint32 width, uint64 number, uint64 *values);

  //  ELIAS GAMMA CODED DATA

  uint64   getEliasGamma(void);
  uint64  *getEliasGamma(uint64 number, uint64 *values=NULL);

  uint32   setEliasGamma(uint64 value);
  uint32   setEliasGamma(uint64 number, uint64 *values);

  //  ELIAS DELTA CODED DATA

  uint64   getEliasDelta(void);
  uint64  *getEliasDelta(uint64 number, uint64 *values=NULL);

  uint32   setEliasDelta(uint64 value);
  uint32   setEliasDelta(uint64 number, uint64 *values);

  //  ELIAS OMEGA CODED DATA - the omega code looks hard to implement - the
  //  encoding and decoding streams are backwards from each other.  The idea
  //  is:
  //
  //    push the binary representation of the value onto a stack.
  //    set value to one less than the number of bits emitted last.
  //    loop until the value is 1.
  //
  //  The stream is constructed by emitting the words on the stack, and
  //  terminating the stream with an extra '0'.
  //
#if 0
  uint64   getEliasOmega(void);
  uint64  *getEliasOmega(uint64 number, uint64 *values=NULL);

  uint32   setEliasOmega(uint64 value);
  uint32   setEliasOmega(uint64 number, uint64 *values);
#endif

  //  GOLOMB CODED DATA
  //
  //  Pick m.  For m == power_of_two, this is RICE CODED DATA.
  //
  //    q = floow(n/m).
  //    r = n-qm
  //    c = ceil(log_2 m)
  //
  //  Unary encode q, binary encode r.
  //
  //  The first 2^c-m values are encoded as c-1 bit values, starting with 00...00,
  //  The rest as c-bit numbers, ending with 11...11
  //


  //  FIBONACCI CODED DATA
  //
  //  A Fibonacci number is F(n) = F(n-1) + F(n-2), wher F(0) = F(1) = 1.
  //
  //  The Zeckendorf representation of a number encodes it such that no
  //  two consecurive Fibonacci numbers are used.  From the definition
  //  of a Fibonacci number, any pattern "100" can be replaced with "011".
  //  A number encoded after this transofmration is in the Fibonacci
  //  representation ("Zeckendorf representation" seems to be a real thing,
  //  I just made up "Fibonacci representation" - the two terms seem to
  //  be used synonymously in the real world).
  //
  //  Once encoded, it's added to the bit stream reversed.
  //
  //  For the Zeckendorf representation, a single 1-bit is added, terminating
  //  the number with the last '1' bit of data, followed immediately by
  //  another '1' bit.  (Because, by definition, there are no two adjacent
  //  set bits in the encoded number).
  //
  //  For the Fibonacci representation, we need to append two '0' bits.
  //  (Because, by definition, there are no two adjacent unset bits in the
  //  representation).  BUT, this representation saves at most one bit
  //  (replacing 100 at the start of the string by 011), but the savings
  //  is lost by the extra stop bit we need.
  //
  uint64   getZeckendorf(void);
  uint64  *getZeckendorf(uint64 number, uint64 *values=NULL);

  uint32   setZeckendorf(uint64 value);
  uint32   setZeckendorf(uint64 number, uint64 *values);

  //  Old meryl uses preDecrement() when using compressed bucket counting.
  //  Nothing else in canu uses these, and they're painful, so left unimplemented.
#if 0
  uint64   preIncrementBinary(uint64 width, uint64 position);
  uint64   postIncrementBinary(uint64 width, uint64 position);
  uint64   preDecrementBinary(uint64 width, uint64 position);
  uint64   postDecrementBinary(uint64 width, uint64 position);
#endif


private:

  //  For writing, update the length of the block to the maximum of where we're at now and the existing length.
  //
  void      updateLen(void) {
    _dataBlockLen[_dataBlk] = std::max(_dataPos, _dataBlockLen[_dataBlk]);
  };

  //  For both reading and writing, move to the next word if we're at the end of the current one.
  //
  void      updateBit(void) {
    if (_dataBit == 0) {
      _dataWrd += 1;
      _dataBit  = 64;
    }
  };

  //  For reading operations, move to the next block if we're at the end of the current one.
  //  For writing operations, this is done before the write, in ensureSpace().
  //
  //  Should be done before any reading operation.  It isn't (strictly) needed at the end
  //  of a read.  The next read will just do it at the start, and the next write
  //  only cares that there is enough space in the current block for the write.
  //
  void      updateBlk(uint64 readLength) {

    assert(_dataBit != 0);
    assert(_dataBit <= 64);

    //  If in this block, nothing to do.

    //fprintf(stderr, "updateBlk()--  at _dataPos " F_U64 "  _dataBlockLen[%u] = " F_U64 "\n", _dataPos, _dataBlk, _dataBlockLen[_dataBlk]);

    if (_dataPos + readLength <= _dataBlockLen[_dataBlk])
      return;

    //  Otherwise, we MUST be at the end of the current block.  If not,
    //  whatever we're trying to read is in the next block (or we're trying
    //  to read something longer than what is here).

    assert(_dataPos == _dataBlockLen[_dataBlk]);

    _dataBlk += 1;

    if (_dataBlk >= _dataBlocksLen)
      fprintf(stderr, "ERROR: _dataBlk = %lu  _dataBlocksLen = %u\n", _dataBlk, _dataBlocksLen);
    assert(_dataBlk < _dataBlocksLen);

    _dataPos  = 0;
    _data     = _dataBlocks[_dataBlk];

    _dataWrd  = 0;
    _dataBit  = 64;
  }

  //  For writing operations, make sure there is enough space for the write in this block.
  //
  void     ensureSpace(uint64 spaceNeeded) {

    assert(_dataBit != 0);
    assert(_dataBit <= 64);

    //  If enough space in the current block, just return.

    if (_dataPos + spaceNeeded < _dataBlockLenMax)
      return;

    //  Othwerwise, terminate the current block.

    _dataBlockLen[_dataBlocksLen - 1] = _dataPos;

    //  Move to the new block.

    _dataBlocksLen++;

    if (_dataBlocksLen >= _dataBlocksMax) {
      setArraySize(_dataBlocks,   _dataBlocksLen, _dataBlocksMax, _dataBlocksLen + 128);
      setArraySize(_dataBlockBgn, _dataBlocksLen, _dataBlocksMax, _dataBlocksLen + 128);
      setArraySize(_dataBlockLen, _dataBlocksLen, _dataBlocksMax, _dataBlocksLen + 128);
    }

    assert(spaceNeeded <= _dataBlockLenMax);

    _dataPos = 0;
    _data    = _dataBlocks[_dataBlocksLen - 1] = new uint64 [_dataBlockLenMax / sizeof(uint64)];

    memset(_data, 0, sizeof(uint64) * _dataBlockLenMax / 64);

    _dataBlockBgn[_dataBlocksLen - 1] = _dataBlockBgn[_dataBlocksLen - 2] + _dataBlockLen[_dataBlocksLen - 2];
    _dataBlockLen[_dataBlocksLen - 1] = 0;

    _dataBlk += 1;
    _dataWrd  = 0;
    _dataBit  = 64;
  };


  uint64   _dataBlockLenMax; //  Allocated length, of each block (in BITS).

  uint32   _dataBlocksLen;   //  Number of allocated data blocks.
  uint32   _dataBlocksMax;   //  Number of blocks we can allocate.

  uint64  *_dataBlockBgn;    //  Starting position, in the global file, of this block.
  uint64  *_dataBlockLen;    //  Length of this block.
  uint64 **_dataBlocks;      //  Just piles of bits.  Nothing interesting here.

  uint64   _dataPos;         //  Position in this block, in BITS.
  uint64  *_data;            //  Pointer to the currently active data block.

  uint64   _dataBlk;         //  Active data block.
  uint64   _dataWrd;         //  Active word in the active data block.
  uint64   _dataBit;         //  Active bit in the active word in the active data block (aka, number of bits left in this word)

  uint64   _fibData[93];      //  A pile of Fibonacci numbers.
};



#endif  //  LIBBITS_H

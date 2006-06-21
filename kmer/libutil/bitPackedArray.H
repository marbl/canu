#ifndef BITPACKEDARRAY_H
#define BITPACKEDARRAY_H

////////////////////////////////////////
//
//  bitPackedArray
//
//  implements an integer array using bit-widths less than word-sizes,
//  e.g., a memory efficient way to store 23 bit numbers.  Numbers may
//  be up to 64 bits wide.
//
//  The array is variable length, and it is implemented as an array,
//  not a list or tree -- accessing element 1,000,000 will allocate
//  elements 0 through 999,999.
//
class bitPackedArray {
public:

  //  Create a bitpacked array with elements of width 'width' using
  //  'segmentSize' KB per segment.  If you know your array is going
  //  to be much bigger or smaller, crank this value.
  //
  bitPackedArray(u32bit valueWidth, u32bit segmentSize = 1024);
  ~bitPackedArray();

  //  No array operator is provided, because we cannot return a
  //  reference to a value that is split across two words (or even a
  //  reference to a value that is not bit aligned in the word).
  //
  u64bit   get(u64bit idx);
  void     set(u64bit idx, u64bit val);

  //  Clear the array.  Since the array is variable sized, you must add
  //  things to a new array before clearing it.
  void     clear(void);

private:
  u32bit   _valueWidth;
  u32bit   _segmentSize;
  u64bit   _nextElement;  //  the first invalid element
  u64bit   _valuesPerSegment;

  u32bit   _numSegments;
  u32bit   _maxSegments;
  u64bit **_segments;
};


//  An array of bits.  Exactly the same as the bitPackedArray, but
//  optimized for width=1.
//
class bitArray {
public:

  bitArray(u32bit segmentSize = 1024);
  ~bitArray();

  u64bit   get(u64bit idx);

  u64bit   getAndSet(u64bit idx);

  void     set(u64bit idx);
  void     clr(u64bit idx);

  void     clear(void);

private:
  void     resize(u64bit s);

  u32bit   _segmentSize;
  u64bit   _valuesPerSegment;

  u32bit   _numSegments;
  u32bit   _maxSegments;
  u64bit **_segments;
};






inline
u64bit
bitArray::get(u64bit idx) {
  u64bit s = idx / _valuesPerSegment;
  u64bit p = idx % _valuesPerSegment;

  u64bit wrd = (p >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (p     ) & 0x000000000000003fllu;

  return((_segments[s][wrd] >> bit) & 0x0000000000000001llu);
}


inline
void
bitArray::resize(u64bit s) {

  if (s < _numSegments)
    return;

  if (s > _maxSegments) {
    _maxSegments = s + 16;
    u64bit **S = new u64bit * [_maxSegments];
    for (u32bit i=0; i<_numSegments; i++)
      S[i] = _segments[i];
    delete [] _segments;
    _segments = S;
  }

  while (_numSegments <= s)
    _segments[_numSegments++] = new u64bit [_segmentSize * 1024 / 8];
}


inline
u64bit
bitArray::getAndSet(u64bit idx) {
  u64bit s = idx / _valuesPerSegment;
  u64bit p = idx % _valuesPerSegment;

  u64bit wrd = (p >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (p     ) & 0x000000000000003fllu;

  u64bit ret = (_segments[s][wrd] >> bit) & 0x0000000000000001llu;
  
  _segments[s][wrd] |= u64bitONE << bit;

  return(ret);
}


inline
void
bitArray::set(u64bit idx) {
  u64bit s = idx / _valuesPerSegment;
  u64bit p = idx % _valuesPerSegment;

  resize(s);

  u64bit wrd = (p >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (p     ) & 0x000000000000003fllu;

  _segments[s][wrd] |= u64bitONE << bit;
}


inline
void
bitArray::clr(u64bit idx) {
  u64bit s = idx / _valuesPerSegment;
  u64bit p = idx % _valuesPerSegment;

  resize(s);

  u64bit wrd = (p >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (p     ) & 0x000000000000003fllu;

  _segments[s][wrd] &= ~(0x0000000000000001llu << bit);
}





#endif  // BITPACKEDARRAY_H



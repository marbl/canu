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
  u32bit   _valuesPerSegment;

  u32bit   _numSegments;
  u32bit   _maxSegments;
  u64bit **_segments;
};

#endif  // BITPACKEDARRAY_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#include "util++.H"

bitPackedArray::bitPackedArray(u32bit valueWidth, u32bit segmentSize) {
  _valueWidth       = valueWidth;
  _segmentSize      = segmentSize;
  _nextElement      = 0;
  _valuesPerSegment = (u64bit)_segmentSize * 1024 * 8 / (u64bit)_valueWidth;

  _numSegments      = 0;
  _maxSegments      = 16;
  _segments         = new u64bit * [_maxSegments];
}


bitPackedArray::~bitPackedArray() {
  for (u32bit i=0; i<_numSegments; i++)
    delete [] _segments[i];
  delete [] _segments;
}


u64bit
bitPackedArray::get(u64bit idx) {
  u64bit s = idx / _valuesPerSegment;
  u64bit p = _valueWidth * (idx % _valuesPerSegment);

  if (idx >= _nextElement) {
    fprintf(stderr, "bitPackedArray::get()-- element index "u64bitFMT" is out of range, only "u64bitFMT" elements.\n",
            idx, _nextElement-1);
    return(0xdeadbeefdeadbeefULL);
  }

  return(getDecodedValue(_segments[s], p, _valueWidth));
}


void
bitPackedArray::set(u64bit idx, u64bit val) {
  u64bit s = idx / _valuesPerSegment;
  u64bit p = _valueWidth * (idx % _valuesPerSegment);

  //fprintf(stderr, "s="u64bitFMT" p="u64bitFMT" segments="u64bitFMT"/"u64bitFMT"\n", s, p, _numSegments, _maxSegments);

  if (idx >= _nextElement)
    _nextElement = idx+1;

  if (s >= _maxSegments) {
    _maxSegments = s + 16;
    u64bit **S = new u64bit * [_maxSegments];
    for (u32bit i=0; i<_numSegments; i++)
      S[i] = _segments[i];
    delete [] _segments;
    _segments = S;
  }

  while (_numSegments <= s)
    _segments[_numSegments++] = new u64bit [_segmentSize * 1024 / 8];

  setDecodedValue(_segments[s], p, _valueWidth, val);
}


void
bitPackedArray::clear(void) {
  for (u32bit s=0; s<_numSegments; s++)
    bzero(_segments[s], _segmentSize * 1024);
}


////////////////////////////////////////

bitArray::bitArray(u32bit segmentSize) {
  _segmentSize      = segmentSize;
  _valuesPerSegment = (u64bit)_segmentSize * 1024 * 8;

  _numSegments      = 0;
  _maxSegments      = 16;
  _segments         = new u64bit * [_maxSegments];
}


bitArray::~bitArray() {
  for (u32bit i=0; i<_numSegments; i++)
    delete [] _segments[i];
  delete [] _segments;
}


void
bitArray::clear(void) {
  for (u32bit s=0; s<_numSegments; s++)
    bzero(_segments[s], _segmentSize * 1024);
}

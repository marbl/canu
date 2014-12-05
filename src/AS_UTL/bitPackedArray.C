#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <strings.h>

#include "util++.H"

bitPackedArray::bitPackedArray(uint32 valueWidth, uint32 segmentSize) {
  _valueWidth       = valueWidth;
  _segmentSize      = segmentSize;
  _nextElement      = 0;
  _valuesPerSegment = (uint64)_segmentSize * 1024 * 8 / (uint64)_valueWidth;

  _numSegments      = 0;
  _maxSegments      = 16;
  _segments         = new uint64 * [_maxSegments];
}


bitPackedArray::~bitPackedArray() {
  for (uint32 i=0; i<_numSegments; i++)
    delete [] _segments[i];
  delete [] _segments;
}


uint64
bitPackedArray::get(uint64 idx) {
  uint64 s = idx / _valuesPerSegment;
  uint64 p = _valueWidth * (idx % _valuesPerSegment);

  if (idx >= _nextElement) {
    fprintf(stderr, "bitPackedArray::get()-- element index "uint64FMT" is out of range, only "uint64FMT" elements.\n",
            idx, _nextElement-1);
    return(0xdeadbeefdeadbeefULL);
  }

  return(getDecodedValue(_segments[s], p, _valueWidth));
}


void
bitPackedArray::set(uint64 idx, uint64 val) {
  uint64 s = idx / _valuesPerSegment;
  uint64 p = _valueWidth * (idx % _valuesPerSegment);

  //fprintf(stderr, "s="uint64FMT" p="uint64FMT" segments="uint64FMT"/"uint64FMT"\n", s, p, _numSegments, _maxSegments);

  if (idx >= _nextElement)
    _nextElement = idx+1;

  if (s >= _maxSegments) {
    _maxSegments = s + 16;
    uint64 **S = new uint64 * [_maxSegments];
    for (uint32 i=0; i<_numSegments; i++)
      S[i] = _segments[i];
    delete [] _segments;
    _segments = S;
  }

  while (_numSegments <= s)
    _segments[_numSegments++] = new uint64 [_segmentSize * 1024 / 8];

  setDecodedValue(_segments[s], p, _valueWidth, val);
}


void
bitPackedArray::clear(void) {
  for (uint32 s=0; s<_numSegments; s++)
    bzero(_segments[s], _segmentSize * 1024);
}


////////////////////////////////////////

bitArray::bitArray(uint32 segmentSize) {
  _segmentSize      = segmentSize;
  _valuesPerSegment = (uint64)_segmentSize * 1024 * 8;

  _numSegments      = 0;
  _maxSegments      = 16;
  _segments         = new uint64 * [_maxSegments];
}


bitArray::~bitArray() {
  for (uint32 i=0; i<_numSegments; i++)
    delete [] _segments[i];
  delete [] _segments;
}


void
bitArray::clear(void) {
  for (uint32 s=0; s<_numSegments; s++)
    bzero(_segments[s], _segmentSize * 1024);
}

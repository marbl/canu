
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "bitPackedArray.H"
#include "bitPacking.H"

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
    fprintf(stderr, "bitPackedArray::get()-- element index " F_U64 " is out of range, only " F_U64 " elements.\n",
            idx, _nextElement-1);
    return(0xdeadbeefdeadbeefULL);
  }

  return(getDecodedValue(_segments[s], p, _valueWidth));
}


void
bitPackedArray::set(uint64 idx, uint64 val) {
  uint64 s = idx / _valuesPerSegment;
  uint64 p = _valueWidth * (idx % _valuesPerSegment);

  //fprintf(stderr, "s=" F_U64 " p=" F_U64 " segments=" F_U64 "/" F_U64 "\n", s, p, _numSegments, _maxSegments);

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

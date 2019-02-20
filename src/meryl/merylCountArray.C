
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

#include "meryl.H"





template<typename VALUE>
merylCountArray<VALUE>::merylCountArray(void) {
  _sWidth       = 0;
  _vWidth       = 0;

  _prefix       = 0;
  _suffix       = NULL;
  _counts       = NULL;

  _nKmers       = 0;

  _bitsPerPage  = 0;
  _nReAlloc     = 0;

  _segSize      = 0;
  _segAlloc     = 0;
  _segments     = NULL;

  _vals         = NULL;

  _nBits        = 0;
  _nBitsTrigger = 0;
  _nBitsOldSize = 0;

  _multiSet     = false;
}



template<typename VALUE>
uint64
merylCountArray<VALUE>::initialize(uint64 prefix, uint32 width, uint32 segsize) {
  _sWidth       = width;

  _prefix       = prefix;
  _suffix       = NULL;
  _counts       = NULL;

  _nKmers       = 0;

  _bitsPerPage  = getPageSize() * 8;
  _nReAlloc     = 0;

  _segSize      = 8 * (segsize * 1024 - 32);   //  Set the segment size to 'segsize' KB,
  _segAlloc     = 0;                           //  in bits, reserving 32 bytes for
  _segments     = NULL;                        //  allocator stuff that we don't control.

  _nBits        = 0;
  _nBitsTrigger = 0;
  _nBitsOldSize = usedSize();

  return(_nBitsOldSize);
}



template<typename VALUE>
uint64
merylCountArray<VALUE>::initializeValues(VALUE maxValue) {

  if (maxValue == 0)
    _vWidth = 0;
  else
    _vWidth = countNumberOfBits64(maxValue);

  _vals = new stuffedBits();

  return(_nBitsOldSize);
}



template<typename VALUE>
merylCountArray<VALUE>::~merylCountArray() {

  removeSegments();
  removeValues();

  delete [] _suffix;
  delete [] _counts;
}



template<typename VALUE>
void
merylCountArray<VALUE>::removeSegments(void) {

  if (_segments == NULL)                  //  If no segments, then
    return;                               //  we've already removed them.

  for (uint32 ss=0; ss<_segAlloc; ss++)   //  Release the segment memory.
    delete [] _segments[ss];

  delete [] _segments;                    //  Release the list of segments...

  _nReAlloc  = 0;

  _segAlloc = 0;                          //  Don't forget to
  _segments = NULL;                       //  foret about it.

  _nBits        = 0;                      //  Indicate that we've stored no data.
  _nBitsTrigger = 0;
  _nBitsOldSize = usedSize();
}



template<typename VALUE>
void
merylCountArray<VALUE>::removeValues(void) {
  delete _vals;
  _vals  = NULL;
}



template<typename VALUE>
void
merylCountArray<VALUE>::addSegment(uint32 seg) {

  if (_segAlloc == 0) {
    resizeArray(_segments, _segAlloc, _segAlloc, 32, resizeArray_copyData | resizeArray_clearNew);
    _nReAlloc++;
  }
  if (seg >= _segAlloc) {
    resizeArray(_segments, _segAlloc, _segAlloc, 2 * _segAlloc, resizeArray_copyData | resizeArray_clearNew);
    _nReAlloc++;
  }
  assert(_segments[seg] == NULL);

  //if (seg > 0)
  //  fprintf(stderr, "Add segment %u\n", seg);

  _segments[seg] = new uint64 [_segSize / 64];
}




//  Unpack the suffixes and remove the data.
template<typename VALUE>
uint64 *
merylCountArray<VALUE>::unpackSuffixes(uint64 nSuffixes) {
  uint64  *suffixes  = new uint64 [nSuffixes];

  //fprintf(stderr, "Allocate %lu suffixes, %lu bytes\n", nSuffixes, sizeof(uint64) * nSuffixes);
  //fprintf(stderr, "Sorting prefix 0x%016" F_X64P " with " F_U64 " total kmers\n", _prefix, nSuffixes);

  for (uint64 kk=0; kk<nSuffixes; kk++)
    suffixes[kk] = get(kk);

  removeSegments();

  return(suffixes);
}



template<typename VALUE>
swv<VALUE> *
merylCountArray<VALUE>::unpackSuffixesAndValues(uint64 nSuffixes) {
  swv<VALUE>  *suffixes  = new swv<VALUE> [nSuffixes];

  assert(_vals != NULL);

  //fprintf(stderr, "Allocate %lu suffixes, %lu bytes\n", nSuffixes, sizeof(swv) * nSuffixes);
  //fprintf(stderr, "Sorting prefix 0x%016" F_X64P " with " F_U64 " total kmers\n", _prefix, nSuffixes);

  _vals->setPosition(0);

  if      (_vWidth == 0)
    for (uint64 kk=0; kk<nSuffixes; kk++)
      suffixes[kk].set(get(kk), _vals->getEliasDelta());
  else
    for (uint64 kk=0; kk<nSuffixes; kk++)
      suffixes[kk].set(get(kk), _vals->getBinary(_vWidth));

  removeSegments();
  removeValues();

  return(suffixes);
}




//
//  Converts raw kmers listed in _segments into counted kmers listed in _suffix and _counts.
//
template<typename VALUE>
void
merylCountArray<VALUE>::countSingleKmers(void) {
  uint64   nSuffixes = _nBits / _sWidth;
  uint64  *suffixes  = unpackSuffixes(nSuffixes);

  //  Sort the data

#ifdef _GLIBCXX_PARALLEL
  __gnu_sequential::
#else
    std::
#endif
    sort(suffixes, suffixes + nSuffixes);

  //  Count the number of distinct kmers, and allocate space for them.

  uint64  nk = 1;

  for (uint64 kk=1; kk<nSuffixes; kk++)
    if (suffixes[kk-1] != suffixes[kk])
      nk++;

  _suffix = new uint64 [nk];
  _counts = new VALUE  [nk];

  //  And generate the counted kmer data.

  _nKmers = 0;

  _counts[_nKmers] = 1;
  _suffix[_nKmers] = suffixes[0];

  for (uint64 kk=1; kk<nSuffixes; kk++) {
    if (suffixes[kk-1] != suffixes[kk]) {
      _nKmers++;
      _counts[_nKmers] = 0;
      _suffix[_nKmers] = suffixes[kk];
    }

    _counts[_nKmers]++;
  }

  _nKmers++;

  //  Remove all the temporary data.

  delete [] suffixes;
};



template<typename VALUE>
void
merylCountArray<VALUE>::countSingleKmersWithValues(void) {
  uint64       nSuffixes = _nBits / _sWidth;
  swv<VALUE>  *suffixes  = unpackSuffixesAndValues(nSuffixes);

  //  Sort the data

#ifdef _GLIBCXX_PARALLEL
  __gnu_sequential::
#else
    std::
#endif
    sort(suffixes, suffixes + nSuffixes);

  //  Count the number of distinct kmers, and allocate space for them.

  uint64  nk = 1;

  for (uint64 kk=1; kk<nSuffixes; kk++)
    if (suffixes[kk-1].getSuffix() != suffixes[kk].getSuffix())
      nk++;

  _suffix = new uint64 [nk];
  _counts = new VALUE  [nk];

  //  And generate the counted kmer data.

  _nKmers = 0;

  _counts[_nKmers] = suffixes[0].getValue();
  _suffix[_nKmers] = suffixes[0].getSuffix();

  for (uint64 kk=1; kk<nSuffixes; kk++) {
    if (suffixes[kk-1].getSuffix() != suffixes[kk].getSuffix()) {
      _nKmers++;
      _counts[_nKmers] = 0;
      _suffix[_nKmers] = suffixes[kk].getSuffix();
    }

    _counts[_nKmers] += suffixes[kk].getValue();
  }

  _nKmers++;

  //  Remove all the temporary data.

  delete [] suffixes;
};



template<typename VALUE>
void
merylCountArray<VALUE>::countMultiSetKmers(void) {
  uint64      nSuffixes = _nBits / _sWidth;
  swv<VALUE> *suffixes  = unpackSuffixesAndValues(nSuffixes);

  //  Sort the data

#ifdef _GLIBCXX_PARALLEL
  __gnu_sequential::
#else
    std::
#endif
    sort(suffixes, suffixes + nSuffixes);

  //  In a multi-set, we dump each and every kmer that is loaded, no merging.

  _suffix = new uint64 [nSuffixes];
  _counts = new VALUE  [nSuffixes];

  //  And generate the counted kmer data.

  _nKmers = nSuffixes;

  for (uint64 kk=0; kk<nSuffixes; kk++) {
    _counts[kk] = suffixes[kk].getValue();
    _suffix[kk] = suffixes[kk].getSuffix();
  }

  //  Remove all the temporary data.

  delete [] suffixes;
};






template<typename VALUE>
void
merylCountArray<VALUE>::countKmers(void) {

  //fprintf(stderr, "merylCountArray<VALUE>::countKmers()-- _nBits %lu -- values=%c multi-set=%c\n",
  //        _nBits, (_vals == NULL) ? 'n' : 'Y', (_multiSet == false) ? 'n' : 'Y');

  if (_nBits == 0) {    //  If no kmers stored, nothing to do, so just
    removeSegments();   //  remove the (unused) storage and return.
    return;
  }

  assert(_nBits % _sWidth == 0);

  if (_vals == NULL)
    countSingleKmers();
  else
    if (_multiSet == false)
      countSingleKmersWithValues();
    else
      countMultiSetKmers();
}



template<typename VALUE>
void
merylCountArray<VALUE>::dumpCountedKmers(kmerCountBlockWriter *out) {
  out->addBlock(_prefix, _nKmers, _suffix, _counts);
}



template<typename VALUE>
void
merylCountArray<VALUE>::removeCountedKmers(void) {

  delete [] _suffix;   _suffix = NULL;
  delete [] _counts;   _counts = NULL;

  _nKmers = 0;
}



//  Give the linker something to link to.
template class merylCountArray<uint32>;
template class merylCountArray<uint64>;


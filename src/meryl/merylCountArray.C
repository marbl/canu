
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
 */

#include "meryl.H"


merylCountArray::merylCountArray(uint64 prefix, uint32 width, uint32 segsize) {
  _width       = width;

  _prefix      = prefix;
  _suffix      = NULL;
  _counts      = NULL;

  _nKmers      = 0;

  _bitsPerPage = getPageSize() * 8;
  _nReAlloc    = 0;

  _segSize     = segsize;// + segsize % _bitsPerPage;
  _segAlloc    = 0;
  _segments    = NULL;

  _nBits       = 0;
}



merylCountArray::~merylCountArray() {

  removeSegments();

  delete [] _suffix;
  delete [] _counts;
}



void
merylCountArray::removeSegments(void) {

  if (_segments == NULL)                  //  If no segments, then
    return;                               //  we've already removed them.

  for (uint32 ss=0; ss<_segAlloc; ss++)   //  Release the segment memory.
    delete [] _segments[ss];

  delete [] _segments;                    //  Release the list of segments...

  _nReAlloc  = 0;

  _segAlloc = 0;                          //  Don't forget to
  _segments = NULL;                       //  foret about it.

  _nBits = 0;                             //  Indicate that we've stored no data.
}



void
merylCountArray::addSegment(uint32 seg) {

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



//
//  Converts raw kmers listed in _segments into counted kmers listed in _suffix and _counts.
//
void
merylCountArray::countKmers(void) {

  if (_nBits == 0) {    //  If no kmers stored, nothing to do, so just
    removeSegments();   //  remove the (unused) storage and return.
    return;
  }

  assert(_nBits % _width == 0);

  uint64   nValues = _nBits / _width;
  uint64  *values  = new uint64 [nValues];

  //fprintf(stderr, "Sorting prefix 0x%016" F_X64P " with " F_U64 " total kmers\n", _prefix, nValues);

  //  Unpack the data into _suffix.

  for (uint64 kk=0; kk<nValues; kk++)
    values[kk] = get(kk);

  //  All done with the raw data, so get rid of it quickly.

  removeSegments();

  //  Sort the data

#ifdef _GLIBCXX_PARALLEL
  __gnu_sequential::
#else
    std::
#endif
    sort(values, values + nValues);

  //  Count the number of distinct kmers, and allocate space for them.
  
  uint64  nk = 1;

  for (uint64 kk=1; kk<nValues; kk++)
    if (values[kk-1] != values[kk])
      nk++;

  _suffix = new uint64 [nk];
  _counts = new uint32 [nk];

  //  And generate the counted kmer data.

  _nKmers = 0;

  _counts[_nKmers] = 1;
  _suffix[_nKmers] = values[0];

  for (uint64 kk=1; kk<nValues; kk++) {
    if (values[kk-1] != values[kk]) {
      _nKmers++;
      _counts[_nKmers] = 0;
      _suffix[_nKmers] = values[kk];
    }

    _counts[_nKmers]++;
  }

  _nKmers++;

  //  Remove all the temporary data.

  delete [] values;
};



void
merylCountArray::dumpCountedKmers(kmerCountBlockWriter *out) {
  out->addBlock(_prefix, _nKmers, _suffix, _counts);
}



void
merylCountArray::removeCountedKmers(void) {

  delete [] _suffix;   _suffix = NULL;
  delete [] _counts;   _counts = NULL;

  _nKmers = 0;
}

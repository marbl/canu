
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

#include "kmers.H"
#include "bits.H"

#include <vector>
#include <algorithm>

using namespace std;


double
bitsToGB(uint64 bits) {
  return(bits / 8 / 1024.0 / 1024.0 / 1024.0);
}




//  Set some basic boring stuff.
//
void
kmerCountExactLookup::initialize(kmerCountFileReader *input_,
                                 uint32               minValue_,
                                 uint32               maxValue_) {

  //  Silently make minValue and maxValue be valid values.

  if (minValue_ == 0)
    minValue_ = 1;

  if (maxValue_ == UINT32_MAX)
    maxValue_ = input_->stats()->maxFrequency();

  //  Now initialize filtering!

  _minValue      = minValue_;
  _maxValue      = maxValue_;
  _valueOffset   = minValue_ - 1;                   //  "1" stored in the data is really "minValue" to the user.

  _nKmersLoaded  = 0;
  _nKmersTooLow  = 0;
  _nKmersTooHigh = 0;

  //  Now initialize table parameters!

  _Kbits         = kmer::merSize() * 2;

  _prefixBits    = 0;                               //  Bits of the kmer used as an index into the table.
  _suffixBits    = 0;                               //  Width of an entry in the suffix table.
  _valueBits     = 0;                               //  (also in the suffix table)

  if (_maxValue >= _minValue)
    _valueBits = logBaseTwo32(_maxValue + 1 - _minValue);

  _suffixMask    = 0;
  _dataMask      = 0;

  _nPrefix       = 0;                               //  Number of entries in pointer table.
  _nSuffix       = input_->stats()->numDistinct();  //  Number of entries in suffix dable.

  _prePtrBits    = logBaseTwo64(_nSuffix);          //  Width of an entry in the prefix table.

  _suffixBgn     = NULL;
  _suffixEnd     = NULL;
  //_suffixLen     = NULL;
  _suffixData    = NULL;
}



//  Analyze the number of kmers to store in the table, to decide on
//  various parameters for allocating the table - how many bits to
//  use for indexing (prefixSize), and how many bits of data we need
//  to store explicitly (suffixBits and valueBits).
//
void
kmerCountExactLookup::configure(void) {

  //  First, find the prefixBits that results in the smallest allocated memory size.
  //  Due to threading over the files, we cannot use a prefix smaller than 6 bits.

  uint64  extraSpace = (uint64)8 * 1024 * 1024 * 1024;   //  In BITS!
  uint64  minSpace   = UINT64_MAX - extraSpace;
  uint64  optSpace   = UINT64_MAX;

  uint32  pbMin      = 0;
  uint32  pbOpt      = 0;

  for (uint32 pb=1; pb<_Kbits; pb++) {
    uint64  nprefix = (uint64)1 << pb;
    uint64  space   = nprefix * _prePtrBits + _nSuffix * (_Kbits - pb) + _nSuffix * _valueBits;

    if (space < minSpace) {
      pbMin        = pb;
      minSpace     = space;
    }

    if (space < minSpace + extraSpace) {
      pbOpt        = pb;
      optSpace     = space;

      _prefixBits  =          pb;
      _suffixBits  = _Kbits - pb;

      _suffixMask  = uint64MASK(_suffixBits);
      _dataMask    = uint64MASK(_valueBits);

      _nPrefix     = nprefix;
    }
  }

  assert(_prefixBits > 0);
  assert(_suffixBits > 0);

  //  And do it all again to keep the users entertained.

  fprintf(stderr, "\n");
  fprintf(stderr, " p       prefixes             bits gigabytes\n");
  fprintf(stderr, "-- -------------- ---------------- ---------\n");

  uint32  minpb = (pbMin < 4)          ? 1      : pbMin - 4;  //  Show four values before and
  uint32  maxpb = (_Kbits < pbOpt + 5) ? _Kbits : pbOpt + 5;  //  four after the smallest.

  for (uint32 pb=minpb; pb < maxpb; pb++) {
    uint64  nprefix = (uint64)1 << pb;
    uint64  space   = nprefix * _prePtrBits + _nSuffix * (_Kbits - pb) + _nSuffix * _valueBits;

    if      (pb == pbMin)
      fprintf(stderr, "%2u %14lu %16lu %9.3f (smallest)\n", pb, nprefix, space, bitsToGB(space));

    else if (pb == pbOpt)
      fprintf(stderr, "%2u %14lu %16lu %9.3f (faster)\n",   pb, nprefix, space, bitsToGB(space));

    else
      fprintf(stderr, "%2u %14lu %16lu %9.3f\n",            pb, nprefix, space, bitsToGB(space));
  }

  fprintf(stderr, "-- -------------- ---------------- ---------\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "For %lu distinct %u-mers (with %u bits used for indexing and %u bits for tags):\n", _nSuffix, _Kbits / 2, _prefixBits, _suffixBits);
  fprintf(stderr, "  %7.3f GB memory\n",                                       bitsToGB(optSpace));
  fprintf(stderr, "  %7.3f GB memory for index (%lu elements %u bits wide)\n", bitsToGB(_nPrefix * _prePtrBits), _nPrefix, _prePtrBits);
  fprintf(stderr, "  %7.3f GB memory for tags  (%lu elements %u bits wide)\n", bitsToGB(_nSuffix * _suffixBits), _nSuffix, _suffixBits);
  fprintf(stderr, "  %7.3f GB memory for data  (%lu elements %u bits wide)\n", bitsToGB(_nSuffix * _valueBits),  _nSuffix, _valueBits);
  fprintf(stderr, "\n");
}




//  With all parameters known, just grab and clear memory.
//
//  The block size used in the wordArray _suffixData is chosen so that large
//  arrays have not-that-many allocations.  The array is pre-allocated, to
//  prevent the need for any locking or coordination when filling out the
//  array.
//
void
kmerCountExactLookup::allocate(void) {

  _suffixBgn = new uint64 [_nPrefix];
  _suffixEnd = new uint64 [_nPrefix];
  //_suffixLen = new uint64 [_nPrefix];

  memset(_suffixBgn, 0, sizeof(uint64) * _nPrefix);
  memset(_suffixEnd, 0, sizeof(uint64) * _nPrefix);
  //memset(_suffixLen, 0, sizeof(uint64) * _nPrefix);

  uint64  arraySize     = _nSuffix * (_suffixBits + _valueBits);
  uint64  arrayBlockMin = max(arraySize / 1024llu, 268435456llu);   //  In bits, so 32MB per block.

  _suffixData = new wordArray(_suffixBits + _valueBits, arrayBlockMin);
  _suffixData->allocate(_nSuffix);
}



//  To allow multiple threads to load each kmer file, we need to knw where to start
//  each file in the _suffixData array.  We can't be precise if we're filtering
//  kmers, but we can at least make them not collide.
//
uint64 *
kmerCountExactLookup::findArrayStartPositions(kmerCountFileReader *input_) {
  uint64  *nKmersPerFile = new uint64 [input_->numFiles()];    //  number of kmers in file ff
  uint64  *startPos      = new uint64 [input_->numFiles()];    //  position in _suffixData that file ff is at

  for (uint32 ff=0; ff<input_->numFiles(); ff++) {
    nKmersPerFile[ff] = 0;
    startPos[ff]      = 0;
  }

  input_->loadBlockIndex();

  for (uint32 ff=0; ff<input_->numFiles(); ff++) {
    for (uint32 bb=ff * input_->numBlocks(); bb<ff * input_->numBlocks() + input_->numBlocks(); bb++)
      nKmersPerFile[ff] += input_->blockIndex(bb).numKmers();

    if (ff > 0)
      startPos[ff] = startPos[ff-1] + nKmersPerFile[ff-1];

    //  If this fails, we found too many kmers in block files as compared
    //  to the number of distinct kmers in the database (_nSuffix).
    assert(startPos[ff] + nKmersPerFile[ff] <= _nSuffix);
  }

  delete [] nKmersPerFile;

  return(startPos);
}







kmerCountExactLookup::kmerCountExactLookup(kmerCountFileReader *input_,
                                           uint32               minValue_,
                                           uint32               maxValue_) {

  initialize(input_, minValue_, maxValue_);  //  Do NOT use minValue_ or maxValue_ from now on!
  configure();
  allocate();

  uint64 *startPos = findArrayStartPositions(input_);
  uint32  nf       = input_->numFiles();

  //  Each file can be processed independently IF we know how many kmers are in
  //  each prefix.  For that, we need to load the kmerCountFileReader index.
  //  We don't, actually, know that if we're filtering out low/high count kmers.
  //  In this case, we overallocate, but cannot cleanup at the end.

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<nf; ff++) {
    FILE                      *blockFile = input_->blockFile(ff);
    kmerCountFileReaderBlock  *block     = new kmerCountFileReaderBlock;

    //  Load blocks until there are no more.

    while (block->loadBlock(blockFile, ff) == true) {
      block->decodeBlock();

      for (uint32 ss=0; ss<block->nKmers(); ss++) {
        uint64   sdata  = 0;
        uint64   prefix = 0;
        uint64   value  = block->counts()[ss];

        if (value < _minValue) {
          _nKmersTooLow++;
          continue;
        }

        if (_maxValue < value) {
          _nKmersTooHigh++;
          continue;
        }

        _nKmersLoaded++;

        //  Reconstruct the kmer into sdata.  This is just kmerTiny::setPrefixSuffix().
        //  From the kmer, generate the prefix we want to save it as.

        sdata   = block->prefix();
        sdata <<= input_->suffixSize();
        sdata  |= block->suffixes()[ss];

        //  Save the prefix we'll be using for storing in our table.  On access, we need
        //  to search the bucket associated with this prefix.

        prefix = sdata >> _suffixBits;

        assert(prefix < _nPrefix);

        //  Add in any extra data to be stored here.

        if (_valueBits > 0) {
          value -= _valueOffset;

          if (value > _maxValue + 1 - _minValue)
            fprintf(stderr, "minValue " F_U32 " maxValue " F_U32 " value " F_U64 " bits " F_U32 "\n",
                    _minValue, _maxValue, value, _valueBits);
          assert(value <= uint64MASK(_valueBits));

          sdata <<= _valueBits;
          sdata  |=  value;
        }

        //  We need to awkwardly remember the start position of each block in
        //  real time, since blocks before this one can be less than full (so
        //  we can't just remember the length of each block and compute
        //  bgn/end positions after the fact.  One oddity that is cleaned up
        //  later is the pointer to the first element, which gets set to 1
        //  instead of 0.
        //
        //  set() will chop off any extra bits (i.e., the prefix) from sdata
        //  before storing it.
        //
        //  After storing, update the end position of this block.

        if (_suffixBgn[prefix] == 0)
          _suffixBgn[prefix] = startPos[ff];

        _suffixData->set(startPos[ff]++, sdata);

        _suffixEnd[prefix] = startPos[ff];
      }
    }

    delete block;

    AS_UTL_closeFile(blockFile);
  }

  //  Fix up the first block pointer.  It's not necessarily the [0]th element in suffixBgn!

  for (uint32 ii=0; ii<_nPrefix; ii++)
    if (_suffixBgn[ii] == 1) {
      _suffixBgn[ii] = 0;
      break;
    }

  delete [] startPos;

  fprintf(stderr, "Loaded " F_U64 " kmers.  Skipped " F_U64 " (too low) and " F_U64 " (too high) kmers.\n",
          _nKmersLoaded, _nKmersTooLow, _nKmersTooHigh);
}



bool
kmerCountExactLookup::exists_test(kmer k) {

  uint64  kmer   = (uint64)k;
  uint64  prefix = kmer >> _suffixBits;
  uint64  suffix = kmer  & _suffixMask;

  uint64  bgn = _suffixBgn[prefix];
  uint64  mid;
  uint64  end = _suffixEnd[prefix];

  uint64  dat;
  uint64  tag;

  //  Binary search for the matching tag.

  while (bgn + 8 < end) {
    mid = bgn + (end - bgn) / 2;

    dat = _suffixData->get(mid);
    tag = dat >> _valueBits;

    if (tag == suffix)
      return(true);

    if (suffix < tag)
      end = mid;

    else
      bgn = mid + 1;
  }

  //  Switch to linear search when we're down to just a few candidates.

  for (mid=bgn; mid < end; mid++) {
    dat = _suffixData->get(mid);
    tag = dat >> _valueBits;

    if (tag == suffix)
      return(true);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "FAILED kmer   0x%016lx\n", kmer);
  fprintf(stderr, "FAILED prefix 0x%016lx\n", prefix);
  fprintf(stderr, "FAILED suffix 0x%016lx\n", suffix);
  fprintf(stderr, "\n");
  fprintf(stderr, "original  %9lu %9lu\n", _suffixBgn[prefix], _suffixEnd[prefix]);
  fprintf(stderr, "final     %9lu %9lu\n", bgn, end);
  fprintf(stderr, "\n");

  bgn = _suffixBgn[prefix];
  end = _suffixEnd[prefix];

  while (bgn + 8 < end) {
    mid = bgn + (end - bgn) / 2;

    dat = _suffixData->get(mid);

    fprintf(stderr, "TEST bgn %8lu %8lu %8lu end -- dat %lu =?= %lu suffix\n", bgn, mid, end, dat, suffix);

    if (dat == suffix)
      return(true);

    if (suffix < dat)
      end = mid;

    else
      bgn = mid + 1;
  }

  for (mid=bgn; mid < end; mid++) {
    dat = _suffixData->get(mid);

    fprintf(stderr, "ITER bgn %8lu %8lu %8lu end -- dat %lu =?= %lu suffix\n", bgn, mid, end, dat, suffix);

    if (dat == suffix)
      return(true);
  }

  assert(0);
};

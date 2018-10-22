
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

#include "kmers.H"
#include "bits.H"

#include <vector>
#include <algorithm>

using namespace std;


//  If set, allocate another (large) array to verify that there are no holes in the
//  data array.  Holes would lead to false positives.
//
#define VERIFY_SUFFIX_END





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

  _minValue       = minValue_;
  _maxValue       = maxValue_;
  _valueOffset    = minValue_ - 1;                   //  "1" stored in the data is really "minValue" to the user.

  _nKmersLoaded   = 0;
  _nKmersTooLow   = 0;
  _nKmersTooHigh  = 0;

  //  Now initialize table parameters!

  _Kbits          = kmer::merSize() * 2;

  _prefixBits     = 0;                               //  Bits of the kmer used as an index into the table.
  _suffixBits     = 0;                               //  Width of an entry in the suffix table.
  _valueBits      = 0;                               //  (also in the suffix table)

  if (_maxValue >= _minValue)
    _valueBits = logBaseTwo32(_maxValue + 1 - _minValue);

  _suffixMask     = 0;
  _dataMask       = 0;

  _nPrefix        = 0;                               //  Number of entries in pointer table.
  _nSuffix        = 0;                               //  Number of entries in suffix dable.

  //  Scan the histogram to count the number of kmers in range.

  for (uint32 ii=_minValue; ii<=_maxValue; ii++)
    _nSuffix += input_->stats()->numKmersAtFrequency(ii);

  _prePtrBits     = logBaseTwo64(_nSuffix);          //  Width of an entry in the prefix table.
  _prePtrBits     = 64;

  _suffixBgn      = NULL;
  _suffixEnd      = NULL;
  _suffixData     = NULL;
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

  uint64  minSpace   = UINT64_MAX;
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

    if (space < 4 * minSpace) {
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

  if (_verbose) {
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

  uint64  arraySize     = _nSuffix * (_suffixBits + _valueBits);
  uint64  arrayBlockMin = max(arraySize / 1024llu, 268435456llu);   //  In bits, so 32MB per block.

  //if (_verbose)
  //  fprintf(stderr, "Allocating space for %lu suffixes of %u bits each -> %lu bits (%lu bytes) in blocks of %lu bytes\n",
  //          _nSuffix, _suffixBits + _valueBits, arraySize, arraySize / 8, arrayBlockMin / 8);

  _suffixData = new wordArray(_suffixBits + _valueBits, arrayBlockMin);
  _suffixData->allocate(_nSuffix);
}




//  Make one pass through the file to count how many kmers per prefix we will end
//  up with.  This is needed only if kmers are filtered, but does
//  make the rest of the loading a little easier.
//
//  The loop control and kmer loading is the same in the two loops.
void
kmerCountExactLookup::count(kmerCountFileReader *input_) {
  uint64  *kpp = new uint64 [_nPrefix];

  memset(kpp, 0, sizeof(uint64) * _nPrefix);

  //  Scan all kmer files, counting the number of kmers per prefix.
  //  This is thread safe when _prefixBits is more than 6 (the number of files).

  uint32   nf = input_->numFiles();

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<nf; ff++) {
    FILE                      *blockFile = input_->blockFile(ff);
    kmerCountFileReaderBlock  *block     = new kmerCountFileReaderBlock;

    //  Keep local counters, otherwise, we collide when updating the global counts.

    uint64  tooLow  = 0;
    uint64  tooHigh = 0;
    uint64  loaded  = 0;

    //  Load blocks until there are no more.

    while (block->loadBlock(blockFile, ff) == true) {
      block->decodeBlock();

      for (uint32 ss=0; ss<block->nKmers(); ss++) {
        uint64   sdata  = 0;
        uint64   prefix = 0;
        uint64   value  = block->counts()[ss];

        if (value < _minValue) {
          tooLow++;
          continue;
        }

        if (_maxValue < value) {
          tooHigh++;
          continue;
        }

        loaded++;

        sdata   = block->prefix();         //  Reconstruct the kmer into sdata.  This is just
        sdata <<= input_->suffixSize();    //  kmerTiny::setPrefixSuffix().  From the kmer,
        sdata  |= block->suffixes()[ss];   //  generate the prefix we want to save it as.

        prefix  = sdata >> _suffixBits;

        assert(prefix < _nPrefix);

        kpp[prefix]++;                     //  Count the number of kmers per prefix.
      }
    }

#pragma omp critical (count_stats)
    {
      _nKmersTooLow  += tooLow;
      _nKmersTooHigh += tooHigh;
      _nKmersLoaded  += loaded;
    }

    delete block;

    AS_UTL_closeFile(blockFile);
  }

  //  Convert the kmers per prefix into begin coordinate for each prefix.
  //  The loading loop uses _suffixEnd[] as the position to add the next
  //  data.

  _suffixBgn = new uint64 [_nPrefix + 1];

  uint64  bgn = 0;

  for (uint64 ii=0; ii<_nPrefix; ii++) {
    _suffixBgn[ii] = bgn;

    bgn += kpp[ii];
  }

  _suffixBgn[_nPrefix] = bgn;
  assert(bgn == _nKmersLoaded);

  delete [] kpp;

#ifdef VERIFY_SUFFIX_END
  _suffixEnd = new uint64 [_nPrefix];

  for (uint64 ii=0; ii<_nPrefix; ii++)
    _suffixEnd[ii] = _suffixBgn[ii];
#endif

  //  Log.

  if (_verbose)
    fprintf(stderr, "Will load " F_U64 " kmers.  Skipping " F_U64 " (too low) and " F_U64 " (too high) kmers.\n",
            _nKmersLoaded, _nKmersTooLow, _nKmersTooHigh);
}



  //  Each file can be processed independently IF we know how many kmers are in
  //  each prefix.  For that, we need to load the kmerCountFileReader index.
  //  We don't, actually, know that if we're filtering out low/high count kmers.
  //  In this case, we overallocate, but cannot cleanup at the end.
void
kmerCountExactLookup::load(kmerCountFileReader *input_) {

  uint32   nf = input_->numFiles();

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

        if ((value < _minValue) ||         //  Sanity checking and counting done
            (_maxValue < value))           //  in count() above.
          continue;

        sdata   = block->prefix();         //  Reconstruct the kmer into sdata.  This is just
        sdata <<= input_->suffixSize();    //  kmerTiny::setPrefixSuffix().  From the kmer,
        sdata  |= block->suffixes()[ss];   //  generate the prefix we want to save it as.

        prefix  = sdata >> _suffixBits;

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

        //  Store the data.

        _suffixData->set(_suffixBgn[prefix]++, sdata);

#ifdef VERIFY_SUFFIX_END
        _suffixEnd[prefix]++;
#endif
      }
    }

    delete block;

    AS_UTL_closeFile(blockFile);
  }

  //  suffixBgn[i] is now the start of [i+1]; shift the array by one to
  //  restore the proper meaning of suffixBgn.

  for (uint32 ii=_nPrefix; ii>0; ii--)
    _suffixBgn[ii] = _suffixBgn[ii-1];

  _suffixBgn[0] = 0;

  //  Optionally verify that bgn[i] == end[i-1].

#ifdef VERIFY_SUFFIX_END
  for (uint32 ii=1; ii<_nPrefix; ii++)
    assert(_suffixBgn[ii] == _suffixEnd[ii-1]);

  delete [] _suffixEnd;
  _suffixEnd = NULL;
#endif

  //  Now just log.

  if (_verbose)
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
  uint64  end = _suffixBgn[prefix + 1];

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
  fprintf(stderr, "original  %9lu %9lu\n", _suffixBgn[prefix], _suffixBgn[prefix + 1]);
  fprintf(stderr, "final     %9lu %9lu\n", bgn, end);
  fprintf(stderr, "\n");

  bgn = _suffixBgn[prefix];
  end = _suffixBgn[prefix + 1];

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

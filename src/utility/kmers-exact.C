
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


  //  Split each kmer into a prefix and a suffix.  The prefix
  //  is used as an index into an array of pointers into a giant
  //  list of suffixes (and data).  The suffix array is then searched
  //  to find the matching suffix, and the data is returned.
  //
  //  Space used (in bits) is:
  //
  //    2^p * log2(nDistinct) + (K - p) * nDistinct
  //

kmerCountExactLookup::kmerCountExactLookup(kmerCountFileReader *input,
                                           uint32               minValue,
                                           uint32               maxValue) {

  _Kbits         = kmer::merSize() * 2;

  _prefixBits    = 0;                               //  Bits of the kmer used as an index into the table.
  _suffixBits    = 0;                               //  Width of an entry in the suffix table.
  _valueBits     = 0;                               //  (also in the suffix table)

  if (minValue == 0)                                //  Silently adjust to the legal minValue.
    minValue = 1;

  _valueOffset   = minValue - 1;                    //  "1" stored in the data is really "minValue" to the user.

  _suffixMask    = 0;
  _dataMask      = 0;

  _nPrefix       = 0;                               //  Number of entries in pointer table.
  _nSuffix       = input->stats()->numDistinct();   //  Number of entries in suffix dable.

  _prePtrBits    = logBaseTwo64(_nSuffix);          //  Width of an entry in the prefix table.

  _suffixStart   = NULL;
  _suffixData    = NULL;

  //  If maxValue isn't set, ask the input what the largest count is.
  //  Then set the valueBits needed to hold those values.

  if (maxValue == UINT32_MAX) {
    kmerCountStatistics *hist = input->stats();

    maxValue = hist->numFrequencies() - 1;

    while ((maxValue > 0) && (hist->numKmersAtFrequency(maxValue) == 0))
      maxValue--;

    fprintf(stderr, "MAX VALUE %u\n", maxValue);
  }

  if (maxValue > minValue)
    _valueBits = logBaseTwo32(maxValue + 1 - minValue);

  //  First, find the prefixBits that results in the smallest allocated memory size.

  uint64  extraSpace = (uint64)8 * 1024 * 1024 * 1024;   //  In BITS!
  uint64  minSpace   = UINT64_MAX - extraSpace;
  //uint64  optSpace   = UINT64_MAX - extraSpace;

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
      //optSpace     = space;

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
  fprintf(stdout, " p       prefixes             bits gigabytes\n");
  fprintf(stdout, "-- -------------- ---------------- ---------\n");

  uint32  minpb = (pbMin < 4)          ? 1      : pbMin - 4;  //  Show four values before and
  uint32  maxpb = (_Kbits < pbOpt + 5) ? _Kbits : pbOpt + 5;  //  four after the smallest.

  for (uint32 pb=minpb; pb < maxpb; pb++) {
    uint64  nprefix = (uint64)1 << pb;
    uint64  space   = nprefix * _prePtrBits + _nSuffix * (_Kbits - pb) + _nSuffix * _valueBits;

    if      (pb == pbMin)
      fprintf(stdout, "%2u %14lu %16lu %9.3f (smallest)\n", pb, nprefix, space, bitsToGB(space));

    else if (pb == pbOpt)
      fprintf(stdout, "%2u %14lu %16lu %9.3f (faster)\n",   pb, nprefix, space, bitsToGB(space));

    else
      fprintf(stdout, "%2u %14lu %16lu %9.3f\n",            pb, nprefix, space, bitsToGB(space));
  }

  fprintf(stdout, "-- -------------- ---------------- ---------\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "For %lu distinct %u-mers (with %u bits used for indexing and %u bits for tags):\n", _nSuffix, _Kbits / 2, _prefixBits, _suffixBits);
  fprintf(stderr, "  %7.3f GB memory\n",                                       bitsToGB(minSpace));
  fprintf(stderr, "  %7.3f GB memory for index (%lu elements %u bits wide)\n", bitsToGB(_nPrefix * _prePtrBits), _nPrefix, _prePtrBits);
  fprintf(stderr, "  %7.3f GB memory for tags  (%lu elements %u bits wide)\n", bitsToGB(_nSuffix * _suffixBits), _nSuffix, _suffixBits);
  fprintf(stderr, "  %7.3f GB memory for data  (%lu elements %u bits wide)\n", bitsToGB(_nSuffix * _valueBits),  _nSuffix, _valueBits);
  fprintf(stderr, "\n");

  //  Allocate space to keep track of the start of the suffix data for each
  //  prefix, and the prefix data.
  //
  //  Over-allocating _suffixStart by one is critical!  The last entry is
  //  used to tell the end of the last bucket of kmers.

  _suffixStart = new uint64 [_nPrefix + 1];

  memset(_suffixStart, 0, sizeof(uint64) * (_nPrefix + 1));

  //  The block size used in the wordArray is chosen so that large arrays
  //  have not-that-many allocations.
  //
  //  The array is pre-allocated, to prevent the need for any locking or
  //  coordination when filling out the array.

  uint64  arraySize     = _nSuffix * (_suffixBits + _valueBits);
  uint64  arrayBlockMin = max(arraySize / 1024llu, 268435456llu);   //  In bits, so 32MB per block.

  _suffixData  = new wordArray(_suffixBits + _valueBits, arrayBlockMin);

  _suffixData->allocate(_nSuffix);

  //  Load kmers from the input file - [ prefix=3 ][ suffix=41 ] - into our arrays
  //  using 15 bits for our prefix/index and 29 bits for the suffix/tag.
  //
  //  Each file can be processed independently IF we know how many kmers are in each prefix.
  //  For that, we need to load the kmerCountFileReader index.
  //
  //  Each block holds kmer data for all kmers that have the same 'prefixSize' bits at the start.
  //  We know how many kmers are in this block.
  //  If our '_prefixBits' is at least as big as 'prefixSize', we can just sum the counts from blocks.
  //  If not, we need to build on the fly...

  fprintf(stderr, "\n");
  fprintf(stderr, "Loading block indexes.\n");
  fprintf(stderr, "\n");

  //  The first '_prefixBits'   bits of the kmer form a block here.
  //  The first 'numFilesBits'  bits of the kmer is used to decide the file in use.
  //  The next  'numBlocksBits' bits of the kmer is stored in the block address which we iterate over.

  input->loadBlockIndex();

  uint64  *nKmersPerFile = new uint64 [input->numFiles()];    //  number of kmers in file ff
  uint64  *startPos      = new uint64 [input->numFiles()];    //  position in _suffixData that file ff is at

  for (uint32 ff=0; ff<input->numFiles(); ff++) {
    nKmersPerFile[ff] = 0;
    startPos[ff]      = 0;
  }

  for (uint32 ff=0; ff<input->numFiles(); ff++) {
    for (uint32 bb=ff * input->numBlocks(); bb<ff * input->numBlocks() + input->numBlocks(); bb++)
      nKmersPerFile[ff] += input->blockIndex(bb).numKmers();

    if (ff > 0)
      startPos[ff] = startPos[ff-1] + nKmersPerFile[ff-1];
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Constructing exact lookup table.\n");
  fprintf(stderr, "\n");

  uint32  nf = input->numFiles();     //  OpenMP wants simple variables for the loop tests.
  uint32  nb = input->numBlocks();

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<nf; ff++) {
    //fprintf(stderr, "STARTING FILE %u with startPos %lu\n", ff, startPos[ff]);

    FILE                      *blockFile = input->blockFile(ff);
    kmerCountFileReaderBlock  *block     = new kmerCountFileReaderBlock;

    for (uint32 bb=0; bb<nb; bb++) {
      block->loadBlock(blockFile, ff);  //  Should check for errors!
      block->decodeBlock();

      //fprintf(stderr, "STARTING BLOCK bb %u prefix %lu at suffixData %lu\n", bb, block->prefix(), startPos[ff]);

      for (uint32 ss=0; ss<block->nKmers(); ss++) {
        uint64   sdata  = 0;
        uint64   prefix = 0;

        sdata   = block->prefix();
        sdata <<= input->suffixSize();
        sdata  |= block->suffixes()[ss];

        //  sdata is now the kmer.  Shift it to generate the prefix, and set _suffixStart.

        prefix = sdata >> _suffixBits;

        //  Add in any extra data to be stored here.  Unfortunately,
        //  we must load the values outside the minValue and maxValue range, otherwise
        //  we end up with holes in the table.  This also means we need to store kmers
        //  with value 0, which are treated as if they don't exist later.

        if (_valueBits > 0) {
          uint32 value = block->counts()[ss];

          if ((value < minValue) ||
              (maxValue < value))
            value = 0;

          if (value != 0)
            value -= _valueOffset;

          if (value > maxValue + 1 - minValue) {
            fprintf(stderr, "minValue %u maxValue %u value %u bits %u\n",
                    minValue, maxValue, value, _valueBits);
          }
          assert(value <= uint64MASK(_valueBits));

          sdata <<= _valueBits;
          sdata  |=  value;
        }

        //  And store it!  We let set() mask off the top end of the kmer that we shouldn't
        //  be storing.

        assert(prefix < _nPrefix);

        _suffixStart[prefix] = startPos[ff] + 1;   //  _suffixStart here is really the start of prefix+1;
        _suffixData->set(startPos[ff], sdata);     //  doing +1 here makes the logic later a bit easier.

        startPos[ff]++;
      }
    }

    delete block;

    AS_UTL_closeFile(blockFile);
  }

  //  Fix up the _suffixStart array.  The [i] entry is currently the
  //  element after the end if [i], that is, what should be in [j].
  //  Simplistically, we'd just need to set [i] <- [i-1].
  //
  //  But there are entries with no suffixes in them; _suffixStart[i] == 0.
  //  For these, we need to set it to the next valid value.  In the worst
  //  case, when there are two kmers, one in bucket 1, one in bucket 4,
  //  the array (before --> after) would be:   (_nPrefix == 5)
  //
  //      [0] [1] [2] [3] [4] [5] --> [0] [1] [2] [3] [4] [5]
  //       0   1   0   0   2   0       0   0   1   1   1   2
  //
  //  The after is read as: bucket 1 starts at entry 0 (from [1]) and
  //  ends at entry 1 (from [2]).
  //
  //  In all cases, [0] == 0 and [nPrefix] == nSuffix.

  _suffixStart[_nPrefix] = _nSuffix;

  //  Shift values from meaning 'last value in this cell' to meaning
  //  'first value in this cell'.

  for (uint64 pp=_nPrefix; pp>0; pp--)
    _suffixStart[pp] = _suffixStart[pp-1];

  //  Fill in the missing start points.

  for (uint64 pp=1; pp<_nPrefix; pp++)
    if (_suffixStart[pp] == 0)
      _suffixStart[pp] = _suffixStart[pp-1];

  _suffixStart[0]        = 0;
  _suffixStart[_nPrefix] = _nSuffix;

#if 0
  FILE *F = fopen("suffix-start-table", "w");
  for (uint64 pp=0; pp<_nPrefix; pp++)
    fprintf(F, "%8lu %lu\n", pp, _suffixStart[pp]);
  fclose(F);
#endif

  //  Tidy up before we leave.

  delete [] nKmersPerFile;
  delete [] startPos;
}


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

#ifndef LIBKMER_H
#define LIBKMER_H

#include "AS_global.H"
#include "files.H"
#include "bits.H"

//  merSize 1 NOT supported.  Fails _leftShift.

#undef  SHOW_LOAD



char *
constructBlockName(char   *prefix,
                   uint64  outIndex,
                   uint32  numFiles,
                   uint32  iteration,
                   bool    isIndex);

FILE *
openOutputBlock(char   *prefix,
                uint64  fileIndex,
                uint32  numFiles,
                uint32  iteration=0);

FILE *
openInputBlock(char   *prefix,
               uint64  fileIndex,
               uint32  numFiles,
               uint32  iteration=0);



class  kmerTiny {
public:
  kmerTiny() {
    _mer = 0;
  };

  ~kmerTiny() {
  };

  static
  void        setSize(uint32 ms, bool beVerbose=false) {
    _merSize    = ms;

    _fullMask   = uint64MASK(ms * 2);

    _leftMask   = uint64MASK(ms * 2 - 2);
    _leftShift  = ((2 * ms - 2) % 64);

    if (beVerbose)
      fprintf(stderr, "Set global kmer size to " F_U32 " (fullMask=0x%016" F_X64P " leftMask=0x%016" F_X64P " leftShift=" F_U32 ")\n",
              _merSize, _fullMask, _leftMask, _leftShift);
  };

  static
  uint32      merSize(void) { return(_merSize); };

  //  Push an ASCII base onto the mer, shifting the mer to the right or left
  //  to make space for the new base.  Unlike the 'standard' two-bit encoding,
  //  these encode bases as A=00, C=01, G=11, T=10.
  //
  void        addR(char base)       { _mer  = (((_mer << 2) & _fullMask) | (((base >> 1) & 0x03llu)          )              );  };
  void        addL(char base)       { _mer  = (((_mer >> 2) & _leftMask) | (((base >> 1) & 0x03llu) ^ 0x02llu) << _leftShift);  };

  //  Reverse-complementation of a kmer involves complementing the bases in
  //  the mer, revesing the order of all the bases, then aligning the bases
  //  to the low-order bits of the word.
  //
  uint64      reverseComplement(uint64 mer) const {

    //  Complement the bases

    mer ^= 0xaaaaaaaaaaaaaaaallu;

    //  Reverse the mer

    mer = ((mer >>  2) & 0x3333333333333333llu) | ((mer <<  2) & 0xccccccccccccccccllu);
    mer = ((mer >>  4) & 0x0f0f0f0f0f0f0f0fllu) | ((mer <<  4) & 0xf0f0f0f0f0f0f0f0llu);
    mer = ((mer >>  8) & 0x00ff00ff00ff00ffllu) | ((mer <<  8) & 0xff00ff00ff00ff00llu);
    mer = ((mer >> 16) & 0x0000ffff0000ffffllu) | ((mer << 16) & 0xffff0000ffff0000llu);
    mer = ((mer >> 32) & 0x00000000ffffffffllu) | ((mer << 32) & 0xffffffff00000000llu);

    //  Shift and mask out the bases not in the mer

    mer >>= 64 - _merSize * 2;
    mer  &= _fullMask;

    return(mer);
  };

  kmerTiny   &reverseComplement(void) {
    _mer = reverseComplement(_mer);
    return(*this);
  };

public:
  bool        operator!=(kmerTiny const &r) const { return(_mer != r._mer); };
  bool        operator==(kmerTiny const &r) const { return(_mer == r._mer); };
  bool        operator< (kmerTiny const &r) const { return(_mer <  r._mer); };
  bool        operator> (kmerTiny const &r) const { return(_mer >  r._mer); };
  bool        operator<=(kmerTiny const &r) const { return(_mer <= r._mer); };
  bool        operator>=(kmerTiny const &r) const { return(_mer >= r._mer); };

  bool        isFirst(void)                 const { return(_mer == 0);                        };
  bool        isLast(void)                  const { return(_mer == (UINT64_MAX & _fullMask)); };

  bool        isCanonical(void)             const { return(_mer <= reverseComplement(_mer));  };
  bool        isPalindrome(void)            const { return(_mer == reverseComplement(_mer));  };

  kmerTiny   &operator++()                        {                           _mer++;  return(*this);  };
  kmerTiny    operator++(int)                     { kmerTiny before = *this;  _mer++;  return(before); };

  kmerTiny   &operator--()                        {                           _mer--;  return(*this);  };
  kmerTiny    operator--(int)                     { kmerTiny before = *this;  _mer--;  return(before); };

public:
  char    *toString(char *str) const {
    for (uint32 ii=0; ii<_merSize; ii++) {
      uint32  bb = (((_mer >> (2 * ii)) & 0x03) << 1);
      str[_merSize-ii-1] = (bb == 0x04) ? ('T') : ('A' + bb);
    }
    str[_merSize] = 0;
    return(str);
  };

  operator uint64 () const {
    return(_mer);
  };

  void     setPrefixSuffix(uint64 prefix, uint64 suffix, uint32 width) {
    _mer  = prefix << width;
    _mer |= suffix;
  };

private:
  void     operator>>=(uint32 x)  { _mer >>= x; };
  void     operator<<=(uint32 x)  { _mer <<= x; };


private:
public:
  uint64         _mer;

  static uint32  _merSize;     //  number of bases in this mer

  static uint64  _fullMask;    //  mask to ensure kmer has exactly _merSize bases in it

  static uint64  _leftMask;    //  mask out the left-most base.
  static uint32  _leftShift;   //  how far to shift a base to append to the left of the kmer
};


typedef kmerTiny kmer;







class kmerCountStatistics {
public:
  kmerCountStatistics();
  ~kmerCountStatistics();

  void      addCount(uint64 count) {

    if (count == 0)
      return;

    if (count == 1)
      _numUnique++;

    _numDistinct += 1;
    _numTotal    += count;

    if (count < _histMax) {
      _hist[count]++;
      return;
    }
  };

  void      clear(void);

  void      dump(stuffedBits *bits);
  void      dump(FILE        *outFile);

  void      load(stuffedBits *bits);
  void      load(FILE        *inFile);

  uint64    numUnique(void)                       { return(_numUnique);   };
  uint64    numDistinct(void)                     { return(_numDistinct); };
  uint64    numTotal(void)                        { return(_numTotal);    };

  uint32    numFrequencies(void)                  { return(_histMax);     };
  uint64    numKmersAtFrequency(uint32 ii)        { return(_hist[ii]);    };
  double    wordFrequencyAtFrequency(uint32 ii)   { return(_hist[ii] * 1000000.0 / _numTotal);  };

  uint32    maxFrequency(void) {
    uint32  maxV = _histMax - 1;

    while ((maxV > 0) && (numKmersAtFrequency(maxV) == 0))
      maxV--;

    return(maxV);
  };

private:
  uint64              _numUnique;
  uint64              _numDistinct;
  uint64              _numTotal;

  uint32              _histMax;    //  Max count that can be stored in _hist.
  uint64             *_hist;

  uint32              _hbigLen;    //  Counts bigger than _histMax are stored
  uint32              _hbigMax;    //  as unsorted arrays.
  uint64             *_hbigCount;
  uint64             *_hbigNumber;
};



void
dumpMerylDataFile(char *name);



class kmerCountFileReaderBlock {
public:
  kmerCountFileReaderBlock() {
    _data       = NULL;

    _prefix     = 0;
    _nKmers     = 0;
    _nKmersMax  = 0;

    _kCode      = 0;
    _unaryBits  = 0;
    _binaryBits = 0;
    _k1         = 0;

    _cCode      = 0;
    _c1         = 0;
    _c2         = 0;

    _suffixes   = NULL;
    _counts     = NULL;
  };

  ~kmerCountFileReaderBlock() {
    delete    _data;
    delete [] _suffixes;
    delete [] _counts;
  };

  bool      loadBlock(FILE *inFile, uint32 activeFile, uint32 activeIteration=0) {

    //  If _data exists, we've already loaded the block, but haven't used it yet.

    if (_data)
      return(true);

    //  Otherwise, allocate _data, read the block from disk.  If nothing loaded,
    //  return false.

    _data = new stuffedBits();

    _prefix = UINT64_MAX;
    _nKmers = 0;

    if (_data->loadFromFile(inFile) == false) { 
      delete _data;
      _data = NULL;

      return(false);
    }

    //  Decode the header of _data, but don't process the kmers yet.

    uint64 pos  = _data->getPosition();
    uint64 m1   = _data->getBinary(64);
    uint64 m2   = _data->getBinary(64);

    _prefix     = _data->getBinary(64);
    _nKmers     = _data->getBinary(64);

    _kCode      = _data->getBinary(8);
    _unaryBits  = _data->getBinary(32);
    _binaryBits = _data->getBinary(32);
    _k1         = _data->getBinary(64);

    _cCode      = _data->getBinary(8);
    _c1         = _data->getBinary(64);
    _c2         = _data->getBinary(64);

#ifdef SHOW_LOAD
    fprintf(stderr, "loadBlock()-- file %u iter %u:\n", activeFile, activeIteration);
    fprintf(stderr, "    prefix     0x%016lx\n", _prefix);
    fprintf(stderr, "    nKmers     " F_U64 "\n", _nKmers);
    fprintf(stderr, "    kCode      " F_U32 "\n", _kCode);
    fprintf(stderr, "    unaryBits  " F_U32 "\n", _unaryBits);
    fprintf(stderr, "    binaryBits " F_U32 "\n", _binaryBits);
    fprintf(stderr, "    k1efix     " F_U64 "\n", _k1);
    fprintf(stderr, "    cCode      " F_U32 "\n", _cCode);
    fprintf(stderr, "    c1         " F_U64 "\n", _c1);
    fprintf(stderr, "    c2         " F_U64 "\n", _c2);
#endif

    if ((m1 != 0x7461446c7972656dllu) ||
        (m2 != 0x0a3030656c694661llu)) {
      fprintf(stderr, "kmerCountFileReader::nextMer()-- Magic number mismatch in activeFile " F_U32 " activeIteration " F_U32 " position " F_U64 ".\n", activeFile, activeIteration, pos);
      fprintf(stderr, "kmerCountFileReader::nextMer()-- Expected 0x7461446c7972656d got 0x%016" F_X64P "\n", m1);
      fprintf(stderr, "kmerCountFileReader::nextMer()-- Expected 0x0a3030656c694661 got 0x%016" F_X64P "\n", m2);
      exit(1);
    }

    return(true);
  };

  //  Decode a the data into OUR OWN suffixe and count arrays.
  void      decodeBlock() {

    if (_data == NULL)
      return;

    resizeArrayPair(_suffixes, _counts, 0, _nKmersMax, _nKmers, resizeArray_doNothing);

    decodeBlock(_suffixes, _counts);
  };

  void      decodeBlock(uint64 *suffixes, uint32 *counts) {

    if (_data == NULL)
      return;

    uint64  thisPrefix = 0;

    //  Decode the suffixes.

    for (uint32 kk=0; kk<_nKmers; kk++) {
      thisPrefix += _data->getUnary();

      suffixes[kk] = (thisPrefix << _binaryBits) | (_data->getBinary(_binaryBits));
    }

    //  Decode the counts.

    for (uint32 kk=0; kk<_nKmers; kk++)
      counts[kk] = _data->getBinary(32);

    delete _data;
    _data = NULL;
  }


  uint64    prefix(void)   { return(_prefix); };
  uint64    nKmers(void)   { return(_nKmers); };

  uint64   *suffixes(void) { return(_suffixes); };   //  direct access to decoded data
  uint32   *counts(void)   { return(_counts);   };

private:
  stuffedBits  *_data;

  uint64        _prefix;       //  The prefix of all kmers in this block
  uint64        _nKmers;       //  The number of kmers in this block
  uint64        _nKmersMax;    //  The number of kmers we've allocated space for in _suffixes and _counts

  uint32        _kCode;        //  Encoding type of kmer, then 128 bits of parameters
  uint32        _unaryBits;    //    bits in the unary prefix  (of the kmer suffix)
  uint32        _binaryBits;   //    bits in the binary suffix (of the kmer suffix)
  uint64        _k1;           //    unused

  uint32        _cCode;        //  Encoding type of the counts, then 128 bits of parameters
  uint64        _c1;           //    unused
  uint64        _c2;           //    unused

  uint64       *_suffixes;     //  Decoded suffixes and counts.
  uint32       *_counts;       //
};



class kmerCountFileIndex {
public:
  kmerCountFileIndex() {
    clear();
  };

  void       set(uint64  prefix,
                 FILE   *F,
                 uint64  nKmers) {

    //fprintf(stderr, "set prefix 0x%016lx from [0x%016lx %8lu %8lu] to [%8lu %8lu]\n",
    //        prefix,
    //        _blockPrefix, _blockPosition, _numKmers,
    //        AS_UTL_ftell(F), _numKmers + nKmers);

    if (_blockPrefix == UINT64_MAX) {
      _blockPrefix   = prefix;
      _blockPosition = AS_UTL_ftell(F);
      _numKmers      = nKmers;
    }

    else {
      _numKmers     += nKmers;
    }

    assert(_blockPrefix   == prefix);
    assert(_blockPosition <= AS_UTL_ftell(F));
  };

  void       clear(void) {
    _blockPrefix   = UINT64_MAX;
    _blockPosition = UINT64_MAX;
    _numKmers      = 0;
  }

  uint64     blockPrefix(void)     { return(_blockPrefix);   };
  uint64     blockPosition(void)   { return(_blockPosition); };
  uint64     numKmers(void)        { return(_numKmers);      };

private:
  uint64    _blockPrefix;
  uint64    _blockPosition;
  uint64    _numKmers;
};



class kmerCountFileReader {
private:
  void    initializeFromMasterIndex(bool  doInitialize,
                                    bool  loadStatistics,
                                    bool  beVerbose);

public:
  kmerCountFileReader(const char *inputName,
                      bool        beVerbose=false);
  ~kmerCountFileReader();

public:
  void    loadStatistics(void);
  void    dropStatistics(void);

public:
  void    enableThreads(uint32 threadFile);

public:
  void    loadBlockIndex(void);

public:
  bool    nextMer(void);
  kmer    theFMer(void)   { return(_kmer);    };
  uint32  theCount(void)  { return(_count);   };

  char   *filename(void)  { return(_inName);  };

  kmerCountStatistics       *stats(void) {
    loadStatistics();
    return(_stats);
  }

  //  For direct access to the kmer blocks.
public:
  uint32  prefixSize(void)     { return(_prefixSize); };
  uint32  suffixSize(void)     { return(_suffixSize); };

  uint32  numFiles(void)       { return(_numFiles);  };
  uint32  numBlocks(void)      { return(_numBlocks); };

  uint32  numFilesBits(void)   { return(_numFilesBits);  };
  uint32  numBlocksBits(void)  { return(_numBlocksBits); };

  FILE   *blockFile(uint32 ff) {
    FILE *F = NULL;

    if (ff < _numFiles)
      F = openInputBlock(_inName, ff, _numFiles);

    return(F);
  };

  kmerCountFileIndex   &blockIndex(uint32 bb) {
    return(_blockIndex[bb]);
  };

private:
  char                       _inName[FILENAME_MAX+1];

  uint32                     _prefixSize;
  uint32                     _suffixSize;
  uint32                     _numFilesBits;
  uint32                     _numBlocksBits;

  uint32                     _numFiles;
  uint32                     _numBlocks;

  kmerCountStatistics       *_stats;

  FILE                      *_datFile;

  kmerCountFileReaderBlock  *_block;
  kmerCountFileIndex        *_blockIndex;

  kmer                       _kmer;
  uint32                     _count;

  uint64                     _prefix;

  uint64                     _activeMer;
  uint32                     _activeFile;

  uint32                     _threadFile;

  uint64                     _nKmers;
  uint64                     _nKmersMax;
  uint64                    *_suffixes;
  uint32                    *_counts;
};


#include "kmers-writer-block.H"
#include "kmers-writer-stream.H"


//  Responsible for actually writing data to a file.
class kmerCountFileWriter {
public:
  kmerCountFileWriter(const char *outputName,
                      uint32      prefixSize = 0);

  ~kmerCountFileWriter();

public:
  char   *filename(void)  { return(_outName);  };

public:
  void    initialize(uint32 prefixSize = 0);

  kmerCountBlockWriter  *getBlockWriter(void)        { return(new kmerCountBlockWriter (this));      };
  kmerCountStreamWriter *getStreamWriter(uint32 ff)  { return(new kmerCountStreamWriter(this, ff));  };

public:
  uint32  numberOfFiles(void)           { return(_numFiles);                      };
  uint64  firstPrefixInFile(uint32 ff)  { return(((uint64)ff) << _numBlocksBits); };
  uint64  lastPrefixInFile(uint32 ff)   { return(firstPrefixInFile(ff + 1) - 1);  };

  uint32  fileNumber(uint64 prefix);

private:
  //void    importStatistics(kmerCountStatistics &import);

private:
  void    writeBlockToFile(FILE                *datFile,
                           kmerCountFileIndex  *datFileIndex,
                           uint64               prefix,
                           uint64               nKmers,
                           uint64              *suffixes,
                           uint32              *counts);

private:
  bool                       _initialized;

  char                       _outName[FILENAME_MAX+1];

  uint32                     _prefixSize;

  uint32                     _suffixSize;
  uint64                     _suffixMask;

  uint32                     _numFilesBits;
  uint32                     _numBlocksBits;
  uint64                     _numFiles;
  uint64                     _numBlocks;

  kmerCountStatistics        _stats;

  friend class kmerCountBlockWriter;
  friend class kmerCountStreamWriter;
};



class kmerCountExactLookup {
public:
  kmerCountExactLookup(kmerCountFileReader *input_,
                       uint32               minValue_ = 0,
                       uint32               maxValue_ = UINT32_MAX) {

    _verbose = false;

    initialize(input_, minValue_, maxValue_);  //  Do NOT use minValue_ or maxValue_ from now on!
    configure();
    count(input_);
    allocate();
    load(input_);
  };

  ~kmerCountExactLookup() {
    delete [] _suffixBgn;
    delete [] _suffixEnd;
    delete    _suffixData;
  };

private:
  void     initialize(kmerCountFileReader *input_, uint32 minValue_, uint32 maxValue_);
  void     configure(void);
  void     count(kmerCountFileReader *input_);
  void     allocate(void);
  void     load(kmerCountFileReader *input_);

private:
  uint32           value_value(uint64 value) {
    if (_valueBits == 0)               //  Return 'true' if no value
      return(1);                       //  is stored.

    value &= uint32MASK(_valueBits);

    //if (value == 0)                    //  Return zero if the value
    //  return(0);                       //  is actually zero.

    return(value + _valueOffset);      //  Otherwise, return the value.
  };

public:
  uint64           nKmers(void)  {  return(_nKmersLoaded);  };

  uint32           value(kmer k) {
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
        return(value_value(dat));

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
        return(value_value(dat));
    }

    return(0);
  };


  bool             exists_test(kmer k);


private:
  bool            _verbose;

  uint32          _minValue;    //  Minimum value stored in the table -| both of these filter the
  uint32          _maxValue;    //  Maximum value stored in the table -| input kmers.
  uint32          _valueOffset; //  Offset of values stored in the table.

  uint64          _nKmersLoaded;
  uint64          _nKmersTooLow;
  uint64          _nKmersTooHigh;

  uint32          _Kbits;

  uint32          _prefixBits;  //  How many high-end bits of the kmer is an index into _suffixBgn.
  uint32          _suffixBits;  //  How many bits of the kmer are in the suffix table.
  uint32          _valueBits;   //  How many bits of the suffix entry are data.

  uint64          _suffixMask;
  uint64          _dataMask;

  uint64          _nPrefix;     //  How many entries in _suffixBgn  == 2 ^ _prefixBits.
  uint64          _nSuffix;     //  How many entries in _suffixData == nDistinct in the input database.

  uint32          _prePtrBits;  //  How many bits wide is _suffixBgn (used only if _suffixBgn is a wordArray).

  uint64         *_suffixBgn;   //  The start of a block of data in suffix Data.  The end is the next start.
  uint64         *_suffixEnd;   //  The end.  Temporary.
  wordArray      *_suffixData;  //  Finally, kmer data!
};



#endif  //  LIBKMER

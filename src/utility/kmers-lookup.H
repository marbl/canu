
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

#ifndef LIBKMER_LOOKUP_H
#define LIBKMER_LOOKUP_H

#include "kmers.H"



class merylExactLookup {
public:
  merylExactLookup(merylFileReader *input_,
                   uint32               maxMemory_ = 0,
                   uint64               minValue_  = 0,
                   uint64               maxValue_  = UINT64_MAX) {

    _input     = input_;
    _maxMemory = maxMemory_;   //  maxMemory_ is In GB; _maxMemory should be in BITS!
    _verbose   = true;

    if (_maxMemory == 0)
      _maxMemory   = getPhysicalMemorySize() * 8;
    else
      _maxMemory <<= 33;

    initialize(minValue_, maxValue_);  //  Do NOT use minValue_ or maxValue_ from now on!
  };

  ~merylExactLookup() {
    delete [] _suffixBgn;
    delete [] _suffixEnd;
    delete    _sufData;
    delete    _valData;
  };

  //  To use this object:
  //    lookup = new merylExactLookup(input, 0, 0, UINT32_MAX);
  //    if (lookup->configure() == true)
  //      lookup->load()
  //

private:
  void     initialize(uint64 minValue_, uint64 maxValue_);
public:
  bool     configure(void);
private:
  void     count(void);
  void     allocate(void);
public:
  void     load(void);

private:
  uint64           value_value(uint64 value) {
    if (_valueBits == 0)               //  Return 'true' if no value
      return(1);                       //  is stored.

    value &= uint64MASK(_valueBits);

    //if (value == 0)                    //  Return zero if the value
    //  return(0);                       //  is actually zero.

    return(value + _valueOffset);      //  Otherwise, return the value.
  };

public:
  uint64           nKmers(void)  {  return(_nKmersLoaded);  };


  //  Return true/false if the kmer exists/does not.
  bool             exists(kmer k) {
    kmdata  kmer   = (kmdata)k;
    uint64  prefix = kmer >> _suffixBits;
    kmdata  suffix = kmer  & _suffixMask;

    uint64  bgn = _suffixBgn[prefix];
    uint64  mid;
    uint64  end = _suffixBgn[prefix + 1];

    kmdata  tag;

    //  Binary search for the matching tag.

    while (bgn + 8 < end) {
      mid = bgn + (end - bgn) / 2;

      tag = _sufData->get(mid);

      if (tag == suffix)
        return(true);

      if (suffix < tag)
        end = mid;

      else
        bgn = mid + 1;
    }

    //  Switch to linear search when we're down to just a few candidates.

    for (mid=bgn; mid < end; mid++) {
      tag = _sufData->get(mid);

      if (tag == suffix)
        return(true);
    }

    return(false);
  }


  //  Return true/false if the kmer exists/does not.
  //  And populate 'value' with the value of the kmer.
  bool             exists(kmer k, uint64 &value) {
    kmdata  kmer   = (kmdata)k;
    uint64  prefix = kmer >> _suffixBits;
    kmdata  suffix = kmer  & _suffixMask;

    uint64  bgn = _suffixBgn[prefix];
    uint64  mid;
    uint64  end = _suffixBgn[prefix + 1];

    kmdata  tag;

    //  Binary search for the matching tag.

    while (bgn + 8 < end) {
      mid = bgn + (end - bgn) / 2;

      tag = _sufData->get(mid);

      if (tag == suffix) {
        if (_valueBits == 0)
          value = 1;
        else
          value = _valData->get(mid);
        return(true);
      }

      if (suffix < tag)
        end = mid;

      else
        bgn = mid + 1;
    }

    //  Switch to linear search when we're down to just a few candidates.

    for (mid=bgn; mid < end; mid++) {
      tag = _sufData->get(mid);

      if (tag == suffix) {
        if (_valueBits == 0)
          value = 1;
        else
          value = _valData->get(mid);
        return(true);
      }
    }

    value = 0;
    return(false);
  }


  //  Returns the value of the kmer, '0' if it doesn't exist.
  uint64           value(kmer k) {
    kmdata  kmer   = (kmdata)k;
    uint64  prefix = kmer >> _suffixBits;
    kmdata  suffix = kmer  & _suffixMask;

    uint64  bgn = _suffixBgn[prefix];
    uint64  mid;
    uint64  end = _suffixBgn[prefix + 1];

    kmdata  tag;

    //  Binary search for the matching tag.

    while (bgn + 8 < end) {
      mid = bgn + (end - bgn) / 2;

      tag = _sufData->get(mid);

      if (tag == suffix) {
        if (_valueBits == 0)
          return(1);
        else
          return(_valData->get(mid));
      }

      if (suffix < tag)
        end = mid;

      else
        bgn = mid + 1;
    }

    //  Switch to linear search when we're down to just a few candidates.

    for (mid=bgn; mid < end; mid++) {
      tag = _sufData->get(mid);

      if (tag == suffix) {
        if (_valueBits == 0)
          return(1);
        else
          return(_valData->get(mid));
      }
    }

    return(0);
  };


  bool             exists_test(kmer k);


private:
  merylFileReader  *_input;

  uint64                _maxMemory;
  bool                  _verbose;

  uint64                _minValue;    //  Minimum value stored in the table -| both of these filter the
  uint64                _maxValue;    //  Maximum value stored in the table -| input kmers.
  uint64                _valueOffset; //  Offset of values stored in the table.

  uint64                _nKmersLoaded;
  uint64                _nKmersTooLow;
  uint64                _nKmersTooHigh;

  uint32                _Kbits;

  uint32                _prefixBits;  //  How many high-end bits of the kmer is an index into _suffixBgn.
  uint32                _suffixBits;  //  How many bits of the kmer are in the suffix table.
  uint32                _valueBits;   //  How many bits of the suffix entry are data.

  kmdata                _suffixMask;
  uint64                _dataMask;

  uint64                _nPrefix;     //  How many entries in _suffixBgn  == 2 ^ _prefixBits.
  uint64                _nSuffix;     //  How many entries in _suffixData == nDistinct in the input database.

  uint32                _prePtrBits;  //  How many bits wide is _suffixBgn (used only if _suffixBgn is a wordArray).

  uint64               *_suffixBgn;   //  The start of a block of data in suffix Data.  The end is the next start.
  uint64               *_suffixEnd;   //  The end.  Temporary.
  wordArray            *_sufData;     //  Finally, kmer suffix data!
  wordArray            *_valData;     //  Finally, value data!
};



#endif  //  LIBKMER_LOOKUP_H

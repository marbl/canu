
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef OVERLAPINCORE_KEYITERATOR_H
#define OVERLAPINCORE_KEYITERATOR_H

#include "runtime.H"
#include "types.H"

//
//  Convert a string to a stream of kmers and a flag alerting if the kmer is
//  valid.  Usage:
//    it = keyIterator(string, kmersize)
//    while (it.next() == true) {
//      ...
//    }
//
//  Note that the kmer is represented in the word 'backwards' from what
//  it is in the string.  The first 4-mer in string "ABCDEF" will be
//  encoded in the integer as "DCBA".
//
//  Implementation details:
//
//  Initialize the 'bgn' position of this kmer to one before the kmer
//  starts, so on the first call to next() it is incremented to be the
//  first base of the kmer.
//
//  Similar for _pos; it will wrap-around to become zero on the first
//  call to next().
//
//  _bad is initialized to show that all bases in the kmer are bad.
//  We'll clear them as sequence is initially loaded, unless the
//  base is actually bad, then we'll just add another set bit.
//

class keyIterator {
private:
  uint64  bitEquiv(char c, uint32 p) {
    return(((c >> 1) & 0x03llu) << (2 * p));
  };

  uint32  charIsBad(char c, uint32 p) {
    if ((c == 'A') || (c == 'C') || (c == 'T') || (c == 'G'))
      return(0);
    return(uint32one << p);
  };


public:
  keyIterator(char const *bases, uint32 kmerSize) {
    _k      = kmerSize - 1;

    _merBgn = bases - 1;
    _merEnd = bases;

    _key    = uint64zero;                          //  Clear the key
    _bad    = buildLowBitMask<uint64>(kmerSize);   //  Flag every base as invalid
    _pos    = uint32max;                           //  Set position to one before zero

    if (_merEnd == nullptr)   //  Otherwise we fail *_merEnd below, and there is
      return;                 //  no need to keep testing this over and over.

    for (uint32 kk=0; ((kk < _k) &&
                       (*_merEnd != 0)); kk++) {
      _bad >>= 1;
      _bad  |= charIsBad(*_merEnd, _k);

      _key >>= 2;
      _key  |= bitEquiv(*_merEnd, _k);

      _merEnd++;
    }
  };

  bool  next(void) {

    if ((_merEnd == nullptr) || (*_merEnd == 0))
      return(false);

    _bad >>= 1;
    _bad  |= charIsBad(*_merEnd, _k);

    _key >>= 2;
    _key  |= bitEquiv(*_merEnd, _k);

    _pos++;

    _merBgn++;
    _merEnd++;

    return(true);
  };

  bool  isValid(void) {
    return(_bad == 0);
  };

public:
  uint32        _k;       //  Position in the kmer we add new bases at

  char const   *_merBgn;  //  First letter in the kmer
  char const   *_merEnd;  //  Next letter to add to the kmer

  uint64        _key;     //  The key == the kmer.
  uint32        _bad;     //  If non-zero, the kmer has non-ACGT in it.

  uint32        _pos;     //  Position in the source string.
};

#endif

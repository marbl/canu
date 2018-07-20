
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

#include "bio++.H"

#ifndef KMER_LITE_H
#define KMER_LITE_H

////////////////////////////////////////
//
//  This is kMerLite -- derived from kMer.H, removing
//  most of the accessors.
//
//  Assumes that KMER_WORDS is already defined.


class kMerLite {
public:

  //  Used by some of the test routines.
  void    dump(void) const {
    for (uint32 i=0; i<KMER_WORDS; i++)
      fprintf(stderr, "kMerLite["uint32FMTW(2)"] = "uint64HEX"\n", i, _wd[i]);
  };


public:
  void    copy(kMer const &that) {
    for (uint32 i=0; i<KMER_WORDS; i++)
      _wd[i] = that.getWord(i);
  };
  void    clear(void) {
    for (uint32 i=0; i<KMER_WORDS; i++)
      _wd[i] = uint64ZERO;
  };

public:
  kMerLite() {
    clear();
  };
  kMerLite(kMer const &that) {
    copy(that);
  };
  ~kMerLite() {
  };

  kMerLite &operator=(kMer const &that) {
    copy(that);
    return(*this);
  };

  char    *merToString(uint32 merSize, char *instr) const {
    uint32  lastWord = merSize >> 5;
    char   *str = instr;

    if ((merSize & uint32MASK(6)) == 0)
      lastWord++;

    //  We build the string right to left, print any partial word
    //  first, then print whole words until we run out of words to
    //  print.

    if (merSize & uint32MASK(5)) {
      ::merToString(merSize & uint32MASK(5), _wd[lastWord], str);
      str += merSize & uint32MASK(5);
    }

    while (lastWord > 0) {
      lastWord--;
      ::merToString(32, _wd[lastWord], str);
      str += 32;
    }

    return(instr);
  };

#if KMER_WORDS == 1

  bool    operator!=(kMerLite const &r) const { return(_wd[0] != r._wd[0]); };
  bool    operator==(kMerLite const &r) const { return(_wd[0] == r._wd[0]); };
  bool    operator<(kMerLite const &r) const  { return(_wd[0] <  r._wd[0]); };
  bool    operator>(kMerLite const &r) const  { return(_wd[0] >  r._wd[0]); };
  bool    operator<=(kMerLite const &r) const { return(_wd[0] <= r._wd[0]); };
  bool    operator>=(kMerLite const &r) const { return(_wd[0] >= r._wd[0]); };

#else

  bool    operator!=(kMerLite const &r) const {
    uint64 res = uint64ZERO;
    for (uint32 i=KMER_WORDS; i--; )
      res |= _wd[i] ^ r._wd[i];
    return(res != uint64ZERO);
  };
  bool    operator==(kMerLite const &r) const {
    uint64 res = uint64ZERO;
    for (uint32 i=KMER_WORDS; i--; )
      res |= _wd[i] ^ r._wd[i];
    return(res == uint64ZERO);
  };

  bool    operator<(kMerLite const &r) const {
    for (uint32 i=KMER_WORDS; i--; ) {
      if (_wd[i] < r._wd[i])  return(true);
      if (_wd[i] > r._wd[i])  return(false);
    }
    return(false);
  };
  bool    operator>(kMerLite const &r) const {
    for (uint32 i=KMER_WORDS; i--; ) {
      if (_wd[i] > r._wd[i])  return(true);
      if (_wd[i] < r._wd[i])  return(false);
    }
    return(false);
  };
  bool    operator<=(kMerLite const &r) const {
    for (uint32 i=KMER_WORDS; i--; ) {
      if (_wd[i] < r._wd[i])  return(true);
      if (_wd[i] > r._wd[i])  return(false);
    }
    return(true);
  };
  bool    operator>=(kMerLite const &r) const {
    for (uint32 i=KMER_WORDS; i--; ) {
      if (_wd[i] > r._wd[i])  return(true);
      if (_wd[i] < r._wd[i])  return(false);
    }
    return(true);
  };
#endif

private:
  uint64   _wd[KMER_WORDS];
};


#endif  //  KMER_LITE_H

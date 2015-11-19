
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
 *  This file is derived from:
 *
 *    kmer/libbio/kmerhuge.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2007-SEP-13 to 2014-APR-11
 *      are Copyright 2007-2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#define MERWORD(N)   _md[N]


class kMerHuge {
public:
  kMerHuge(uint32 ms=uint32ZERO) {
    setMerSize(ms);
    clear();
  };
  ~kMerHuge() {
  };

  void    setMerSize(uint32 ms);
  uint32  getMerSize(void) const { return(_merSize); };

  void    setMerSpan(uint32 ms)  { _merSpan = ms; };
  uint32  getMerSpan(void) const { return(_merSpan); };

  kMerHuge  &reverseComplement(void) {
    for (uint32 i=0, j=KMER_WORDS-1; i<KMER_WORDS/2; i++, j--) {
      uint64 t = MERWORD(i);
      MERWORD(i) = MERWORD(j);
      MERWORD(j) = t;
    }

    for (uint32 i=0; i<KMER_WORDS; i++) {
      MERWORD(i) = ((MERWORD(i) >>  2) & 0x3333333333333333llu) | ((MERWORD(i) <<  2) & 0xccccccccccccccccllu);
      MERWORD(i) = ((MERWORD(i) >>  4) & 0x0f0f0f0f0f0f0f0fllu) | ((MERWORD(i) <<  4) & 0xf0f0f0f0f0f0f0f0llu);
      MERWORD(i) = ((MERWORD(i) >>  8) & 0x00ff00ff00ff00ffllu) | ((MERWORD(i) <<  8) & 0xff00ff00ff00ff00llu);
      MERWORD(i) = ((MERWORD(i) >> 16) & 0x0000ffff0000ffffllu) | ((MERWORD(i) << 16) & 0xffff0000ffff0000llu);
      MERWORD(i) = ((MERWORD(i) >> 32) & 0x00000000ffffffffllu) | ((MERWORD(i) << 32) & 0xffffffff00000000llu);
      MERWORD(i) ^= 0xffffffffffffffffllu;
    }

    *this >>= KMER_WORDS * 64 - 2 * _merSize;

    return(*this);
  };


  void    clear(void) {
    for (uint32 i=0; i<KMER_WORDS; i++)
      MERWORD(i) = uint64ZERO;
  };
  void    smallest(void) {
    clear();
  };
  void    largest(void) {
    clear();
    reverseComplement();
  };

private:
  void     operator>>=(uint32 x) {

    //  thisWord, the word we shift bits into
    //  thatWord, the word we shift bits out of
    //  shift, the number of bits we shift
    //
    uint32  thisWord = 0;
    uint32  thatWord = x >> 6;
    uint32  shift    = x & uint32MASK(6);

    //  Do an initial word-size shift, to reduce the shift amount to
    //  be less than wordsize.  Fill any shifted-out words with zero.
    //
    if (thatWord) {
      while (thatWord < KMER_WORDS)
        MERWORD(thisWord++) = MERWORD(thatWord++);
      while (thisWord < KMER_WORDS)
        MERWORD(thisWord++) = 0;
    }

    //  Do bit-size shift, of adjacent words
    //
    thisWord = 0;
    thatWord = 1;
    MERWORD(thisWord) >>= shift;
    while (thatWord < KMER_WORDS) {
      MERWORD(thisWord++)  |= MERWORD(thatWord) << (64 - shift);
      MERWORD(thatWord++) >>= shift;
    }
  };

  void      operator<<=(uint32 x) {
    uint32  thisWord = KMER_WORDS;
    uint32  thatWord = KMER_WORDS - (x >> 6);
    uint32  shift    = x & uint32MASK(6);

    if (thatWord != KMER_WORDS) {
      while (thatWord > 0)
        MERWORD(--thisWord) = MERWORD(--thatWord);
      while (thisWord > 0)
        MERWORD(--thisWord) = 0;
    }

    thisWord = KMER_WORDS;
    thatWord = KMER_WORDS - 1;
    MERWORD(thisWord-1) <<= shift;
    while (thatWord > 0) {
      --thisWord;
      --thatWord;
      MERWORD(thisWord)  |= MERWORD(thatWord) >> (64 - shift);
      MERWORD(thatWord) <<= shift;
    }
  };


public:
  void   operator+=(uint64 x) {
    *this <<= 2;
    assert((x & 0xfc) == 0);
    MERWORD(0) |= x & uint64NUMBER(0x3);
  };
  void   operator-=(uint64 x) {
    *this >>= 2;
    assert((x & 0xfc) == 0);
    MERWORD(_lastWord) |= (x & uint64NUMBER(0x3)) << _lastShift;
  };

  void     mask(bool full) {
    MERWORD(_maskWord) &= _mask;
    if (full)
      for (uint32 x=_maskWord+1; x<KMER_WORDS; x++)
        MERWORD(x) = uint64ZERO;
  };

public:
  bool    operator!=(kMerHuge const &r) const {
    uint64 res = uint64ZERO;
    for (uint32 i=KMER_WORDS; i--; )
      res |= MERWORD(i) ^ r.MERWORD(i);
    return(res != uint64ZERO);
  };
  bool    operator==(kMerHuge const &r) const {
    uint64 res = uint64ZERO;
    for (uint32 i=KMER_WORDS; i--; )
      res |= MERWORD(i) ^ r.MERWORD(i);
    return(res == uint64ZERO);
  };

  bool    operator<(kMerHuge const &r) const {
    for (uint32 i=KMER_WORDS; i--; ) {
      if (MERWORD(i) < r.MERWORD(i))  return(true);
      if (MERWORD(i) > r.MERWORD(i))  return(false);
    }
    return(false);
  };
  bool    operator>(kMerHuge const &r) const {
    for (uint32 i=KMER_WORDS; i--; ) {
      if (MERWORD(i) > r.MERWORD(i))  return(true);
      if (MERWORD(i) < r.MERWORD(i))  return(false);
    }
    return(false);
  };
  bool    operator<=(kMerHuge const &r) const {
    for (uint32 i=KMER_WORDS; i--; ) {
      if (MERWORD(i) < r.MERWORD(i))  return(true);
      if (MERWORD(i) > r.MERWORD(i))  return(false);
    }
    return(true);
  };
  bool    operator>=(kMerHuge const &r) const {
    for (uint32 i=KMER_WORDS; i--; ) {
      if (MERWORD(i) > r.MERWORD(i))  return(true);
      if (MERWORD(i) < r.MERWORD(i))  return(false);
    }
    return(true);
  };
  int     qsort_less(kMerHuge const &r) const {
    for (uint32 i=KMER_WORDS; i--; ) {
      if (MERWORD(i) < r.MERWORD(i))  return(-1);
      if (MERWORD(i) > r.MERWORD(i))  return(1);
    }
    return(0);
  };


public:
  operator uint64 () const {return(MERWORD(0));};


public:
  //  these should work generically for both big and small

  void   writeToBitPackedFile(bitPackedFile *BPF, uint32 numBits=0) const {
    if (numBits == 0)
      numBits = _merSize << 1;

    uint32  lastWord = numBits >> 6;

    if ((numBits & uint32MASK(6)) == 0)
      lastWord++;

    if (numBits & uint32MASK(6))
      BPF->putBits(MERWORD(lastWord), numBits & uint32MASK(6));
    while (lastWord > 0) {
      lastWord--;
      BPF->putBits(MERWORD(lastWord), 64);
    }
  };
  void   readFromBitPackedFile(bitPackedFile *BPF, uint32 numBits=0) {
    if (numBits == 0)
      numBits = _merSize << 1;

    uint32  lastWord = numBits >> 6;

    if ((numBits & uint32MASK(6)) == 0)
      lastWord++;

    if (numBits & uint32MASK(6))
      MERWORD(lastWord) = BPF->getBits(numBits & uint32MASK(6));
    while (lastWord > 0) {
      lastWord--;
      MERWORD(lastWord) = BPF->getBits(64);
    }
  };


public:
  //  these should work generically for both big and small

  void     setBits(uint32 pos, uint32 numbits, uint64 val) {
    uint32  wrd = pos >> 6;
    uint32  bit = pos  & 0x3f;

    val &= uint64MASK(numbits);

    if (wrd >= KMER_WORDS) {
      fprintf(stderr, "kMer::setBits()-- ERROR: tried to set pos="F_U32" numbits="F_U32" larger than KMER_WORDS=%d\n",
              pos, numbits, KMER_WORDS), exit(1);
    }

    //  If we have enough space in the word for the bits, replace
    //  those bits in the word.  Otherwise we need to split the value
    //  into two pieces, and add to the end of the first word and the
    //  start of the second.

    if (64 - bit >= numbits) {
      MERWORD(wrd) &= ~(uint64MASK(numbits) << bit);
      MERWORD(wrd) |=  val << bit;
    } else {
      if (wrd+1 >= KMER_WORDS) {
        fprintf(stderr, "kMer::setBits()-- ERROR: tried to set pos="F_U32" numbits="F_U32" larger than KMER_WORDS=%d\n",
                pos, numbits, KMER_WORDS), exit(1);
      }

      uint32 b1 = 64 - bit;      //  bits in the first word
      uint32 b2 = numbits - b1;  //  bits in the second word

      MERWORD(wrd) &= ~(uint64MASK(b1) << bit);
      MERWORD(wrd) |= (val & uint64MASK(b1)) << bit;

      MERWORD(wrd+1) &= ~(uint64MASK(b2));
      MERWORD(wrd+1) |= (val >> b1) & uint64MASK(b2);
    }
  };

  uint64   getBits(uint32 pos, uint32 numbits) const {
    uint64  val = uint64ZERO;
    uint32  wrd = pos >> 6;
    uint32  bit = pos  & 0x3f;

    if (wrd >= KMER_WORDS) {
      fprintf(stderr, "kMer::getBits()-- ERROR: tried to get pos="F_U32" numbits="F_U32" larger than KMER_WORDS=%d\n",
              pos, numbits, KMER_WORDS), exit(1);
    }

    if (64 - bit >= numbits) {
      val = MERWORD(wrd) >> bit;
    } else {
      if (wrd+1 >= KMER_WORDS) {
        fprintf(stderr, "kMer::getBits()-- ERROR: tried to get pos="F_U32" numbits="F_U32" larger than KMER_WORDS=%d\n",
                pos, numbits, KMER_WORDS), exit(1);
      }

      uint32 b1 = 64 - bit;      //  bits in the first word
      uint32 b2 = numbits - b1;  //  bits in the second word

      val  = MERWORD(wrd) >> (64-b1);
      val |= (MERWORD(wrd+1) & uint64MASK(b2)) << b1;
    }

    val &= uint64MASK(numbits);
    return(val);
  };


public:
  //  these should work generically for both big and small

  uint64   startOfMer(uint32 bits) const {
    return(getBits((_merSize << 1) - bits, bits));
  };
  uint64   endOfMer(uint32 bits) const {
    return(MERWORD(0) & uint64MASK(bits));
  };

public:
  //  these should work generically for both big and small
  uint64   getWord(uint32 wrd) const        { return(MERWORD(wrd)); };
  void     setWord(uint32 wrd, uint64 val)  { MERWORD(wrd) = val;   };

public:
  char    *merToString(char *instr) const;

private:
  uint64   _md[KMER_WORDS];

  //  The _merSize is always the number of letters in the mer -- if we
  //  are a spaced seed, it is the weight.
  //
  uint32   _merSize;
  uint32   _merSpan;

  //  The mask is used to make sure the mer has only _merSize bases
  //  set -- we can get more than that if we shift to the left.  The
  //  _maskWord is the word that we want to mask:
  //
  uint64   _mask;
  uint32   _maskWord;

  //  For operator-=() (add a base to the left end) we need to know
  //  what the last word is, and how far to shift the bits.
  //
  //  _lastWord  -- the last word that contains bases
  //  _lastShift -- the amount we need to shift left to put bits 0 and 1
  //                into the last base
  uint32   _lastWord;
  uint32   _lastShift;
};



inline
void
kMerHuge::setMerSize(uint32 ms) {
  _merSize      = ms;
  _merSpan      = ms;
  _lastWord     = (2 * ms - 2) / 64;
  _lastShift    = (2 * ms - 2) % 64;

  _mask     = uint64ZERO;
  _maskWord = _merSize / 32;

  //  Filled whole words with the mer, the mask is special-cased
  //  to clear the whole next word, unless there is no whole next
  //  word, then it does nothing on the last word.
  //
  //  Otherwise, we can construct the mask as usual.
  //
  if        ((_merSize % 32) == 0) {
    if (_maskWord >= KMER_WORDS) {
      _maskWord = KMER_WORDS - 1;
      _mask     = ~uint64ZERO;
    } else {
      _maskWord = _merSize / 32;
      _mask     = uint64ZERO;
    }
  } else {
    _mask = uint64MASK((_merSize % 32) << 1);
  }

  if (_maskWord >= KMER_WORDS) {
    fprintf(stderr, "kMer::setMerSize()-- ERROR!  Desired merSize of "F_U32" larger than\n", _merSize);
    fprintf(stderr, "                     available storage space (KMER_WORDS=%d, max merSize %d).\n", KMER_WORDS, KMER_WORDS*32);
    exit(1);
  }
}




inline
char *
uint64ToMerString(uint32 ms, uint64 mer, char *str) {
  for (uint32 i=0; i<ms; i++)
    str[ms-i-1] = bitsToLetter[(mer >> (2*i)) & 0x03];
  str[ms] = 0;
  return(str);
}



inline
char *
kMerHuge::merToString(char *instr) const {
  uint32  lastWord = _merSize >> 5;
  char   *str = instr;

  if ((_merSize & uint32MASK(6)) == 0)
    lastWord++;

  if (_merSize & uint32MASK(5)) {
    uint64ToMerString(_merSize & uint32MASK(5), MERWORD(lastWord), str);
    str += _merSize & uint32MASK(5);
  }

  while (lastWord > 0) {
    lastWord--;
    uint64ToMerString(32, MERWORD(lastWord), str);
    str += 32;
  }

  return(instr);
};

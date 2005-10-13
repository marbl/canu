// Copyright (c) 2005 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

//  A 'simple' kMer datastructure.

#ifndef BIO_KMER_H
#define BIO_KMER_H

//  The maximum size of a mer.  You get 32 bases per word, so
//  KMER_WORDS=4 will get you up to a 128-mer.
//
#define KMER_WORDS  1


//  XXX: Should probably make a separate class for KMER_WORDS == 1,
//  instead of the current pile of #if's below.


//  For operator-=() (add a base to the left end) we need to know
//  what the last word is, and how far to shift the bits.
//
//  _merSize   -- mer size in bases
//  _lastWord  -- the last word that contains bases
//  _lastShift -- the amount we need to shift left to put bits 0 and 1
//                into the last base

class kMer {
public:
  inline
  void    setMerSize(u32bit ms);

public:
  //  If constructed using the default mersize of zero, you'll get a
  //  crash when trying to use the mer, see setMerSize().
  //
  kMer(u32bit ms=u32bitZERO)  { setMerSize(ms); clear(); };
  ~kMer()                     { };


#if 0

  //  We were calling reverseComplementMer incorrectly, but it looks
  //  like we never called this function anyway!

  void   reverseComplement(void) {
#if KMER_WORDS == 1
    reverseComplementMer(_merSize, _wd[0]);
#else

    //  Reverse all the words, reverse and complement the bases in
    //  each word, then shift right to align the edge.

    for (u32bit i=0, j=KMER_WORDS-1; i<KMER_WORDS/2; i++, j--) {
      u64bit t = _wd[i];
      _wd[i] = _wd[j];
      _wd[j] = t;
    }

    for (u32bit i=0; i<KMER_WORDS; i++) {
      _wd[i] = ((_wd[i] >>  2) & 0x3333333333333333llu) | ((_wd[i] <<  2) & 0xccccccccccccccccllu);
      _wd[i] = ((_wd[i] >>  4) & 0x0f0f0f0f0f0f0f0fllu) | ((_wd[i] <<  4) & 0xf0f0f0f0f0f0f0f0llu);
      _wd[i] = ((_wd[i] >>  8) & 0x00ff00ff00ff00ffllu) | ((_wd[i] <<  8) & 0xff00ff00ff00ff00llu);
      _wd[i] = ((_wd[i] >> 16) & 0x0000ffff0000ffffllu) | ((_wd[i] << 16) & 0xffff0000ffff0000llu);
      _wd[i] = ((_wd[i] >> 32) & 0x00000000ffffffffllu) | ((_wd[i] << 32) & 0xffffffff00000000llu);
      _wd[i] ^= 0xffffffffffffffffllu;
    }

    *this >>= KMER_WORDS * 64 - 2 * _merSize;
#endif
  };

#endif


  //  Return the right-most bits in the mer, XXX, used only by bitPackedFileBuilder
  u64bit subseq_mask(u64bit m) const {
    return(_wd[0] & m);
  };

  void    clear(void) {
#if KMER_WORDS == 1
    _wd[0] = u64bitZERO;
#else
    for (u32bit i=0; i<KMER_WORDS; i++)
      _wd[i] = u64bitZERO;
#endif
  };

  void     operator>>=(u32bit x) {

#if KMER_WORDS == 1
    _wd[0] >>= x;
#else

    //  thisWord, the word we shift bits into
    //  thatWord, the word we shift bits out of
    //  shift, the number of bits we shift
    //
    u32bit  thisWord = 0;
    u32bit  thatWord = x >> 6;
    u32bit  shift    = x & u32bitMASK(6);

    //  Do an initial word-size shift, to reduce the shift amount to
    //  be less than wordsize.  Fill any shifted-out words with zero.
    //
    if (thatWord) {
      while (thatWord < KMER_WORDS)
        _wd[thisWord++] = _wd[thatWord++];
      while (thisWord < KMER_WORDS)
      _wd[thisWord++] = 0;
    }

    //  Do bit-size shift, of adjacent words
    //
    thisWord = 0;
    thatWord = 1;
    _wd[thisWord] >>= shift;
    while (thatWord < KMER_WORDS) {
      _wd[thisWord++]  |= _wd[thatWord] << (64 - shift);
      _wd[thatWord++] >>= shift;
    }
#endif
  };

  void      operator<<=(u32bit x) {

#if KMER_WORDS == 1
    _wd[0] <<= x;
#else
    u32bit  thisWord = KMER_WORDS;
    u32bit  thatWord = KMER_WORDS - (x >> 6);
    u32bit  shift    = x & u32bitMASK(6);

    if (thatWord != KMER_WORDS) {
      while (thatWord > 0)
        _wd[--thisWord] = _wd[--thatWord];
      while (thisWord > 0)
        _wd[--thisWord] = 0;
    }

    thisWord = KMER_WORDS;
    thatWord = KMER_WORDS - 1;
    _wd[thisWord-1] <<= shift;
    while (thatWord > 0) {
      --thisWord;
      --thatWord;
      _wd[thisWord]  |= _wd[thatWord] >> (64 - shift);
      _wd[thatWord] <<= shift;
    }
#endif
  };

  //  Construct a mer by shifting bases onto the end:
  //    += shifts onto the right end
  //    -= shifts onto the left end
  //
  void   operator+=(u64bit x) {
#if KMER_WORDS == 1
    _wd[0] <<= 2;
    _wd[0] |= x & u64bitNUMBER(0x3);
#else
    *this <<= 2;
    _wd[0] |= x & u64bitNUMBER(0x3);
#endif
  };
  void   operator-=(u64bit x) {
#if KMER_WORDS == 1
    _wd[0] >>= 2;
    _wd[0] |= (x & u64bitNUMBER(0x3)) << _lastShift;
#else
    *this >>= 2;
    _wd[_lastWord] |= (x & u64bitNUMBER(0x3)) << _lastShift;
#endif
  };


  //  Convert a kMer into a u64bit (short) mer.
  //
  //  Used in meryl/merle.C
#if 1
  operator u64bit () const {
    return(_wd[0]);
  };
#endif

  //  Mask the mer to the proper size.
  //
  void     mask(bool full) {

#if KMER_WORDS == 1
    _wd[0] &= _mask;
#else
    _wd[_maskWord] &= _mask;
    if (full)
      for (u32bit x=_maskWord+1; x<KMER_WORDS; x++)
        _wd[x] = u64bitZERO;
#endif
  };


  //  Used by some of the test routines.
  void    dump(FILE *out) const {
    for (u32bit i=0; i<KMER_WORDS; i++)
      fprintf(out, "kMer["u32bitFMTW(2)"] = "u64bitHEX"\n", i, _wd[i]);
  };


  //  Returns a sub-mer from either the start (left end) or the end
  //  (right end) of the mer.  The sub-mer must be at most 64 bits
  //  long.
  //
  //  Yes, BITS.
  //
  //  The start is difficult, because it can span multiple words.  The
  //  end is always in the first word.
  //
  u64bit   startOfMer(u32bit bits) const {
    return(getBits((_merSize << 1) - bits, bits));
  };
  u64bit   endOfMer(u32bit bits) const {
    return(_wd[0] & u64bitMASK(bits));
  };

  //  Set 'numbits' bits from (the end of) 'val' at bit position 'pos'
  //  in the mer.  This is wildly low-level, but merylStreamReader
  //  needs it.
  //
  //  The position is measured from the right end.
  //    (0, 8, X) would copy the bits 7 to 0 of X to bits 7 to 0 of the mer.
  //
  //  Argh!  Can't use set/getDecodedValue because that is doing things in the wrong order.
  //
  u64bit   getWord(u32bit wrd) const {
    return(_wd[wrd]);
  };
  void     setWord(u32bit wrd, u64bit val) {
    _wd[wrd] = val;
  };

  void     setBits(u32bit pos, u32bit numbits, u64bit val) {
    u32bit  wrd = pos >> 6;
    u32bit  bit = pos  & 0x3f;

    val &= u64bitMASK(numbits);

    if (wrd >= KMER_WORDS) {
      fprintf(stderr, "kMer::setBits()-- ERROR: tried to set pos="u32bitFMT" numbits="u32bitFMT" larger than KMER_WORDS=%d\n",
              pos, numbits, KMER_WORDS), exit(1);
    }

    //  If we have enough space in the word for the bits, replace
    //  those bits in the word.  Otherwise we need to split the value
    //  into two pieces, and add to the end of the first word and the
    //  start of the second.

    if (64 - bit >= numbits) {
      _wd[wrd] &= ~(u64bitMASK(numbits) << bit);
      _wd[wrd] |=  val << bit;
    } else {
      if (wrd+1 >= KMER_WORDS) {
        fprintf(stderr, "kMer::setBits()-- ERROR: tried to set pos="u32bitFMT" numbits="u32bitFMT" larger than KMER_WORDS=%d\n",
                pos, numbits, KMER_WORDS), exit(1);
      }

      u32bit b1 = 64 - bit;      //  bits in the first word
      u32bit b2 = numbits - b1;  //  bits in the second word

      _wd[wrd] &= ~(u64bitMASK(b1) << bit);
      _wd[wrd] |= (val & u64bitMASK(b1)) << bit;

      _wd[wrd+1] &= ~(u64bitMASK(b2));
      _wd[wrd+1] |= (val >> b1) & u64bitMASK(b2);
    }
  };

  u64bit   getBits(u32bit pos, u32bit numbits) const {
    u64bit  val = u64bitZERO;
    u32bit  wrd = pos >> 6;
    u32bit  bit = pos  & 0x3f;

    if (wrd >= KMER_WORDS) {
      fprintf(stderr, "kMer::getBits()-- ERROR: tried to get pos="u32bitFMT" numbits="u32bitFMT" larger than KMER_WORDS=%d\n",
              pos, numbits, KMER_WORDS), exit(1);
    }

    if (64 - bit >= numbits) {
      val = _wd[wrd] >> bit;
    } else {
      if (wrd+1 >= KMER_WORDS) {
        fprintf(stderr, "kMer::getBits()-- ERROR: tried to get pos="u32bitFMT" numbits="u32bitFMT" larger than KMER_WORDS=%d\n",
                pos, numbits, KMER_WORDS), exit(1);
      }

      u32bit b1 = 64 - bit;      //  bits in the first word
      u32bit b2 = numbits - b1;  //  bits in the second word

      val  = _wd[wrd] >> (64-b1);
      val |= (_wd[wrd+1] & u64bitMASK(b2)) << b1;
    }

    val &= u64bitMASK(numbits);
    return(val);

  };





  //  Show the mer as ascii
  //
  //  Doesn't print the last full word, if it's on the word boundary
  //
  char    *merToString(char *instr) const {
    u32bit  lastWord = _merSize >> 5;
    char   *str = instr;

    if ((_merSize & u32bitMASK(6)) == 0)
      lastWord++;

    //  We build the string right to left, print any partial word
    //  first, then print whole words until we run out of words to
    //  print.

    if (_merSize & u32bitMASK(5)) {
      ::merToString(_merSize & u32bitMASK(5), _wd[lastWord], str);
      str += _merSize & u32bitMASK(5);
    }

    while (lastWord > 0) {
      lastWord--;
      ::merToString(32, _wd[lastWord], str);
      str += 32;
    }

    return(instr);
  };


  //  These are written/read in 5'endian, which isn't the most natural
  //  implementation.  It's done this way to keep the sequence in
  //  order (e.g., the merStreamFile).  Don't change the order.
  //
  //  On the otherhand, the implementation (of write anyway) is
  //  basically the same as merToString().
  //
  //  Takes an optional number of BITS to write, pulled from the
  //  END of the mer.
  //
  void   writeToBitPackedFile(bitPackedFile *BPF, u32bit numBits=0) const {
    if (numBits == 0)
      numBits = _merSize << 1;

    u32bit  lastWord = numBits >> 6;

    if ((numBits & u32bitMASK(6)) == 0)
      lastWord++;

    if (numBits & u32bitMASK(6))
      BPF->putBits(_wd[lastWord], numBits & u32bitMASK(6));
    while (lastWord > 0) {
      lastWord--;
      BPF->putBits(_wd[lastWord], 64);
    }
  };
  void   readFromBitPackedFile(bitPackedFile *BPF, u32bit numBits=0) {
    if (numBits == 0)
      numBits = _merSize << 1;

    u32bit  lastWord = numBits >> 6;

    if ((numBits & u32bitMASK(6)) == 0)
      lastWord++;

    if (numBits & u32bitMASK(6))
      _wd[lastWord] = BPF->getBits(numBits & u32bitMASK(6));
    while (lastWord > 0) {
      lastWord--;
      _wd[lastWord] = BPF->getBits(64);
    }
  };


  kMer   &operator=(kMer const &r) {
    copy(r);
    return(*this);
  };

  void    copy(kMer const &that) {
    memcpy(this, &that, sizeof(kMer));
  };


#if KMER_WORDS == 1

  bool    operator!=(kMer const &r) const { return(_wd[0] != r._wd[0]); };
  bool    operator==(kMer const &r) const { return(_wd[0] == r._wd[0]); };
  bool    operator<(kMer const &r) const  { return(_wd[0] <  r._wd[0]); };
  bool    operator>(kMer const &r) const  { return(_wd[0] >  r._wd[0]); };
  bool    operator<=(kMer const &r) const { return(_wd[0] <= r._wd[0]); };
  bool    operator>=(kMer const &r) const { return(_wd[0] >= r._wd[0]); };

#else

  bool    operator!=(kMer const &r) const {
    u64bit res = u64bitZERO;
    for (u32bit i=KMER_WORDS; i--; )
      res |= _wd[i] ^ r._wd[i];
    return(res != u64bitZERO);
  };
  bool    operator==(kMer const &r) const {
    u64bit res = u64bitZERO;
    for (u32bit i=KMER_WORDS; i--; )
      res |= _wd[i] ^ r._wd[i];
    return(res == u64bitZERO);
  };

  bool    operator<(kMer const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (_wd[i] < r._wd[i])  return(true);
      if (_wd[i] > r._wd[i])  return(false);
    }
    return(false);
  };
  bool    operator>(kMer const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (_wd[i] > r._wd[i])  return(true);
      if (_wd[i] < r._wd[i])  return(false);
    }
    return(false);
  };
  bool    operator<=(kMer const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (_wd[i] < r._wd[i])  return(true);
      if (_wd[i] > r._wd[i])  return(false);
    }
    return(true);
  };
  bool    operator>=(kMer const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (_wd[i] > r._wd[i])  return(true);
      if (_wd[i] < r._wd[i])  return(false);
    }
    return(true);
  };


#endif




private:
  u64bit   _wd[KMER_WORDS];
  u64bit   _mask;
  u32bit   _merSize;
  u32bit   _maskWord;
  u32bit   _lastWord;
  u32bit   _lastShift;
};



inline
void
kMer::setMerSize(u32bit ms) {
  _merSize      = ms;
  _lastWord     = (2 * ms - 2) / 64;
  _lastShift    = (2 * ms - 2) % 64;

  //  The mask is used to make sure the mer has only _merSize bases
  //  set -- we can get more than that if we shift to the left.  The
  //  _maskWord is the word that we want to mask:
  //
  _mask     = u64bitZERO;
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
      _mask     = ~u64bitZERO;
    } else {
      _maskWord = _merSize / 32;
      _mask     = u64bitZERO;
    }
  } else {
    _mask = u64bitMASK((_merSize % 32) << 1);
  }

  if (_maskWord >= KMER_WORDS) {
    fprintf(stderr, "kMer::setMerSize()-- ERROR!  Desired merSize of "u32bitFMT" larger than\n", _merSize);
    fprintf(stderr, "                     available storage space (KMER_WORDS=%d, max merSize %d).\n", KMER_WORDS, KMER_WORDS*32);
    exit(1);
  }
}

#endif  //  BIO_KMER_H

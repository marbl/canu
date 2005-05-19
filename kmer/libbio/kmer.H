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
  //  If constructed using the default mersize of zero, you'll get a
  //  crash when trying to use the mer, see setMerSize().
  //
  kMer(u32bit ms=u32bitZERO)  { setMerSize(ms); clear(); };
  ~kMer()                     { };

  void   reverseComplement(void) {
    // XXX:  UNIMPLEMENTED
  };

  //  This is a total hack to get us up and running with existing code.
  //
  u64bit smallmer(void) const {
    return(_wd[0]);
  };

  void    setMerSize(u32bit ms);


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
  operator u64bit () const {
    return(_wd[0]);
  };



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


  void    dump(FILE *out) const {
    for (u32bit i=0; i<KMER_WORDS; i++)
      fprintf(out, "kMer["u32bitFMTW(2)"] = "u64bitHEX"\n", i, _wd[i]);
  };



  //  Show the mer as ascii
  //
  //  Doesn't print the last full word, if it's on the word boundary
  //
  char    *merToString(char *instr) const {
    u32bit  lastWord = _merSize >> 5;
    char   *str = instr;

    if (_merSize & u32bitMASK(6) == 0)
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
  void   writeToBitPackedFile(bitPackedFile *BPF) const {
    u32bit  lastWord = _merSize >> 5;

    if (_merSize & u32bitMASK(6) == 0)
      lastWord++;

    if (_merSize & u32bitMASK(5))
      BPF->putBits(_wd[lastWord], (_merSize & u32bitMASK(5)) << 1);
    while (lastWord > 0) {
      lastWord--;
      BPF->putBits(_wd[lastWord], 64);
    }
  };
  void   readFromBitPackedFile(bitPackedFile *BPF) {
    u32bit  lastWord = _merSize >> 5;

    if (_merSize & u32bitMASK(6) == 0)
      lastWord++;

    if (_merSize & u32bitMASK(5))
      _wd[lastWord] = BPF->getBits((_merSize & u32bitMASK(5) << 1));
    while (lastWord > 0) {
      lastWord--;
      _wd[lastWord] = BPF->getBits(64);
    }
  };


  bool    operator!=(kMer &r) const {
#if KMER_WORDS == 1
    return(_wd[0] != r._wd[0]);
#else
    u64bit res = u64bitZERO;
    for (u32bit i=0; i<KMER_WORDS; i++) {
      res |= _wd[i] ^ r._wd[i];
      //fprintf(stderr, "!= "u64bitHEX" -- "u64bitHEX" --- "u64bitHEX"\n", _wd[i], r._wd[i], res);
    }
    return(res != u64bitZERO);
#endif
  };
  bool    operator==(kMer &r) const {
#if KMER_WORDS == 1
    return(_wd[0] == r._wd[0]);
#else
    u64bit res = u64bitZERO;
    for (u32bit i=0; i<KMER_WORDS; i++) {
      res |= _wd[i] ^ r._wd[i];
      //fprintf(stderr, "== "u64bitHEX" -- "u64bitHEX" --- "u64bitHEX"\n", _wd[i], r._wd[i], res);
    }
    return(res == u64bitZERO);
#endif
  };



private:
  u64bit   _wd[KMER_WORDS];
  u32bit   _merSize;
  u64bit   _mask;
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

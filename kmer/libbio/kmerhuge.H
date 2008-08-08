
#define MERWORD(N)   _md[N]


class kMerHuge {
public:
  kMerHuge(u32bit ms=u32bitZERO) {
    setMerSize(ms);
    clear();
  };
  ~kMerHuge() {
  };

  void    setMerSize(u32bit ms);
  u32bit  getMerSize(void) const { return(_merSize); };

  void    setMerSpan(u32bit ms)  { _merSpan = ms; };
  u32bit  getMerSpan(void) const { return(_merSpan); };

  kMerHuge  &reverseComplement(void) {
    for (u32bit i=0, j=KMER_WORDS-1; i<KMER_WORDS/2; i++, j--) {
      u64bit t = MERWORD(i);
      MERWORD(i) = MERWORD(j);
      MERWORD(j) = t;
    }

    for (u32bit i=0; i<KMER_WORDS; i++) {
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
    for (u32bit i=0; i<KMER_WORDS; i++)
      MERWORD(i) = u64bitZERO;
  };
  void    smallest(void) {
    clear();
  };
  void    largest(void) {
    clear();
    reverseComplement();
  };

private:
  void     operator>>=(u32bit x) {

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

  void      operator<<=(u32bit x) {
    u32bit  thisWord = KMER_WORDS;
    u32bit  thatWord = KMER_WORDS - (x >> 6);
    u32bit  shift    = x & u32bitMASK(6);

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
  void   operator+=(u64bit x) {
    *this <<= 2;
    assert((x & 0xfc) == 0);
    MERWORD(0) |= x & u64bitNUMBER(0x3);
  };
  void   operator-=(u64bit x) {
    *this >>= 2;
    assert((x & 0xfc) == 0);
    MERWORD(_lastWord) |= (x & u64bitNUMBER(0x3)) << _lastShift;
  };

  void     mask(bool full) {
    MERWORD(_maskWord) &= _mask;
    if (full)
      for (u32bit x=_maskWord+1; x<KMER_WORDS; x++)
        MERWORD(x) = u64bitZERO;
  };

public:
  bool    operator!=(kMerHuge const &r) const {
    u64bit res = u64bitZERO;
    for (u32bit i=KMER_WORDS; i--; )
      res |= MERWORD(i) ^ r.MERWORD(i);
    return(res != u64bitZERO);
  };
  bool    operator==(kMerHuge const &r) const {
    u64bit res = u64bitZERO;
    for (u32bit i=KMER_WORDS; i--; )
      res |= MERWORD(i) ^ r.MERWORD(i);
    return(res == u64bitZERO);
  };

  bool    operator<(kMerHuge const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (MERWORD(i) < r.MERWORD(i))  return(true);
      if (MERWORD(i) > r.MERWORD(i))  return(false);
    }
    return(false);
  };
  bool    operator>(kMerHuge const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (MERWORD(i) > r.MERWORD(i))  return(true);
      if (MERWORD(i) < r.MERWORD(i))  return(false);
    }
    return(false);
  };
  bool    operator<=(kMerHuge const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (MERWORD(i) < r.MERWORD(i))  return(true);
      if (MERWORD(i) > r.MERWORD(i))  return(false);
    }
    return(true);
  };
  bool    operator>=(kMerHuge const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (MERWORD(i) > r.MERWORD(i))  return(true);
      if (MERWORD(i) < r.MERWORD(i))  return(false);
    }
    return(true);
  };
  int     qsort_less(kMerHuge const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (MERWORD(i) < r.MERWORD(i))  return(-1);
      if (MERWORD(i) > r.MERWORD(i))  return(1);
    }
    return(0);
  };


public:
  operator u64bit () const {return(MERWORD(0));};


public:
  //  these should work generically for both big and small

  void   writeToBitPackedFile(bitPackedFile *BPF, u32bit numBits=0) const {
    if (numBits == 0)
      numBits = _merSize << 1;

    u32bit  lastWord = numBits >> 6;

    if ((numBits & u32bitMASK(6)) == 0)
      lastWord++;

    if (numBits & u32bitMASK(6))
      BPF->putBits(MERWORD(lastWord), numBits & u32bitMASK(6));
    while (lastWord > 0) {
      lastWord--;
      BPF->putBits(MERWORD(lastWord), 64);
    }
  };
  void   readFromBitPackedFile(bitPackedFile *BPF, u32bit numBits=0) {
    if (numBits == 0)
      numBits = _merSize << 1;

    u32bit  lastWord = numBits >> 6;

    if ((numBits & u32bitMASK(6)) == 0)
      lastWord++;

    if (numBits & u32bitMASK(6))
      MERWORD(lastWord) = BPF->getBits(numBits & u32bitMASK(6));
    while (lastWord > 0) {
      lastWord--;
      MERWORD(lastWord) = BPF->getBits(64);
    }
  };


public:
  //  these should work generically for both big and small

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
      MERWORD(wrd) &= ~(u64bitMASK(numbits) << bit);
      MERWORD(wrd) |=  val << bit;
    } else {
      if (wrd+1 >= KMER_WORDS) {
        fprintf(stderr, "kMer::setBits()-- ERROR: tried to set pos="u32bitFMT" numbits="u32bitFMT" larger than KMER_WORDS=%d\n",
                pos, numbits, KMER_WORDS), exit(1);
      }

      u32bit b1 = 64 - bit;      //  bits in the first word
      u32bit b2 = numbits - b1;  //  bits in the second word

      MERWORD(wrd) &= ~(u64bitMASK(b1) << bit);
      MERWORD(wrd) |= (val & u64bitMASK(b1)) << bit;

      MERWORD(wrd+1) &= ~(u64bitMASK(b2));
      MERWORD(wrd+1) |= (val >> b1) & u64bitMASK(b2);
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
      val = MERWORD(wrd) >> bit;
    } else {
      if (wrd+1 >= KMER_WORDS) {
        fprintf(stderr, "kMer::getBits()-- ERROR: tried to get pos="u32bitFMT" numbits="u32bitFMT" larger than KMER_WORDS=%d\n",
                pos, numbits, KMER_WORDS), exit(1);
      }

      u32bit b1 = 64 - bit;      //  bits in the first word
      u32bit b2 = numbits - b1;  //  bits in the second word

      val  = MERWORD(wrd) >> (64-b1);
      val |= (MERWORD(wrd+1) & u64bitMASK(b2)) << b1;
    }

    val &= u64bitMASK(numbits);
    return(val);
  };


public:
  //  these should work generically for both big and small

  u64bit   startOfMer(u32bit bits) const {
    return(getBits((_merSize << 1) - bits, bits));
  };
  u64bit   endOfMer(u32bit bits) const {
    return(MERWORD(0) & u64bitMASK(bits));
  };

public:
  //  these should work generically for both big and small
  u64bit   getWord(u32bit wrd) const        { return(MERWORD(wrd)); };
  void     setWord(u32bit wrd, u64bit val)  { MERWORD(wrd) = val;   };

public:
  char    *merToString(char *instr) const;

private:
  u64bit   _md[KMER_WORDS];

  //  The _merSize is always the number of letters in the mer -- if we
  //  are a spaced seed, it is the weight.
  //
  u32bit   _merSize;
  u32bit   _merSpan;

  //  The mask is used to make sure the mer has only _merSize bases
  //  set -- we can get more than that if we shift to the left.  The
  //  _maskWord is the word that we want to mask:
  //
  u64bit   _mask;
  u32bit   _maskWord;

  //  For operator-=() (add a base to the left end) we need to know
  //  what the last word is, and how far to shift the bits.
  //
  //  _lastWord  -- the last word that contains bases
  //  _lastShift -- the amount we need to shift left to put bits 0 and 1
  //                into the last base
  u32bit   _lastWord;
  u32bit   _lastShift;
};



inline
void
kMerHuge::setMerSize(u32bit ms) {
  _merSize      = ms;
  _merSpan      = ms;
  _lastWord     = (2 * ms - 2) / 64;
  _lastShift    = (2 * ms - 2) % 64;

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




inline
char *
kMerHuge::merToString(char *instr) const {
  u32bit  lastWord = _merSize >> 5;
  char   *str = instr;

  if ((_merSize & u32bitMASK(6)) == 0)
    lastWord++;

  if (_merSize & u32bitMASK(5)) {
    u64bitToMerString(_merSize & u32bitMASK(5), MERWORD(lastWord), str);
    str += _merSize & u32bitMASK(5);
  }

  while (lastWord > 0) {
    lastWord--;
    u64bitToMerString(32, MERWORD(lastWord), str);
    str += 32;
  }

  return(instr);
};

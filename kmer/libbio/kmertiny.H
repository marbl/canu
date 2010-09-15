
class kMerTiny {
public:
  kMerTiny(u32bit ms=u32bitZERO) {
    setMerSize(ms);
    clear();
  };
  ~kMerTiny() {
  };

  void    setMerSize(u32bit ms);
  u32bit  getMerSize(void) const { return(_merSize); };

  void    setMerSpan(u32bit ms)  { _merSpan = ms; };
  u32bit  getMerSpan(void) const { return(_merSpan); };

  kMerTiny  &reverseComplement(void) {
    _md = reverseComplementMer(_merSize, _md);
    return(*this);
  };

  void    clear(void) {
    _md = u64bitZERO;
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
    _md >>= x;
  };
  void      operator<<=(u32bit x) {
    _md <<= x;
  };

public:
  void   operator+=(u64bit x) {
    *this <<= 2;
    assert((x & 0xfc) == 0);
    _md |= x & u64bitNUMBER(0x3);
  };
  void   operator-=(u64bit x) {
    *this >>= 2;
    assert((x & 0xfc) == 0);
    _md |= (x & u64bitNUMBER(0x3)) << _lastShift;
  };

public:
  void     mask(bool) {
    _md &= _mask;
  };

public:
  bool    operator!=(kMerTiny const &r) const { return(_md != r._md); };
  bool    operator==(kMerTiny const &r) const { return(_md == r._md); };
  bool    operator< (kMerTiny const &r) const { return(_md <  r._md); };
  bool    operator> (kMerTiny const &r) const { return(_md >  r._md); };
  bool    operator<=(kMerTiny const &r) const { return(_md <= r._md); };
  bool    operator>=(kMerTiny const &r) const { return(_md >= r._md); };
  int     qsort_less(kMerTiny const &r) const {
    if (_md < r._md) return(-1);
    if (_md > r._md) return( 1);
    return(0);
  };
public:
  operator u64bit () const {return(_md);};

public:
  void   writeToBitPackedFile(bitPackedFile *BPF, u32bit numBits=0) const {
    BPF->putBits(_md, _merSize << 1);
  };
  void   readFromBitPackedFile(bitPackedFile *BPF, u32bit numBits=0) {
    _md = BPF->getBits(_merSize << 1);
  };

public:
  void     setBits(u32bit pos, u32bit numbits, u64bit val) {
    _md &= ~(u64bitMASK(numbits) << pos);
    _md |=  val << pos;
  };

  u64bit   getBits(u32bit pos, u32bit numbits) const {
    return((_md >> pos) & u64bitMASK(numbits));
  };

public:
  u64bit   startOfMer(u32bit bits) const {
    return(getBits((_merSize << 1) - bits, bits));
  };
  u64bit   endOfMer(u32bit bits) const {
    return(_md & u64bitMASK(bits));
  };

public:
  u64bit   getWord(u32bit wrd) const        { return(_md); };
  void     setWord(u32bit wrd, u64bit val)  { _md = val;   };

public:
  char    *merToString(char *instr) const;

private:
  u64bit   _md;

  //  The _merSize is always the number of letters in the mer -- if we
  //  are a spaced seed, it is the weight.
  //
  u32bit   _merSize;
  u32bit   _merSpan;

  //  The mask is used to make sure the mer has only _merSize bases
  //  set -- we can get more than that if we shift to the left.  The
  //
  u64bit   _mask;

  //  For operator-=() (add a base to the left end) we need to know
  //  what the last word is, and how far to shift the bits.
  //
  u32bit   _lastShift;
};





inline
void
kMerTiny::setMerSize(u32bit ms) {
  _merSize   = ms;
  _merSpan   = ms;
  _lastShift = (2 * ms - 2) % 64;
  _mask      = u64bitMASK(_merSize  << 1);
}


inline
char *
kMerTiny::merToString(char *str) const {
  for (u32bit i=0; i<_merSize; i++)
    str[_merSize-i-1] = bitsToLetter[(_md >> (2*i)) & 0x03];
  str[_merSize] = 0;
  return(str);
}

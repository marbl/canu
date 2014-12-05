
class kMerTiny {
public:
  kMerTiny(uint32 ms=uint32ZERO) {
    setMerSize(ms);
    clear();
  };
  ~kMerTiny() {
  };

  void    setMerSize(uint32 ms);
  uint32  getMerSize(void) const { return(_merSize); };

  void    setMerSpan(uint32 ms)  { _merSpan = ms; };
  uint32  getMerSpan(void) const { return(_merSpan); };

  kMerTiny  &reverseComplement(void) {
    _md = reverseComplementMer(_merSize, _md);
    return(*this);
  };

  void    clear(void) {
    _md = uint64ZERO;
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
    _md >>= x;
  };
  void      operator<<=(uint32 x) {
    _md <<= x;
  };

public:
  void   operator+=(uint64 x) {
    *this <<= 2;
    assert((x & 0xfc) == 0);
    _md |= x & uint64NUMBER(0x3);
  };
  void   operator-=(uint64 x) {
    *this >>= 2;
    assert((x & 0xfc) == 0);
    _md |= (x & uint64NUMBER(0x3)) << _lastShift;
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
  operator uint64 () const {return(_md);};

public:
  void   writeToBitPackedFile(bitPackedFile *BPF, uint32 numBits=0) const {
    BPF->putBits(_md, _merSize << 1);
  };
  void   readFromBitPackedFile(bitPackedFile *BPF, uint32 numBits=0) {
    _md = BPF->getBits(_merSize << 1);
  };

public:
  void     setBits(uint32 pos, uint32 numbits, uint64 val) {
    _md &= ~(uint64MASK(numbits) << pos);
    _md |=  val << pos;
  };

  uint64   getBits(uint32 pos, uint32 numbits) const {
    return((_md >> pos) & uint64MASK(numbits));
  };

public:
  uint64   startOfMer(uint32 bits) const {
    return(getBits((_merSize << 1) - bits, bits));
  };
  uint64   endOfMer(uint32 bits) const {
    return(_md & uint64MASK(bits));
  };

public:
  uint64   getWord(uint32 wrd) const        { return(_md); };
  void     setWord(uint32 wrd, uint64 val)  { _md = val;   };

public:
  char    *merToString(char *instr) const;

private:
  uint64   _md;

  //  The _merSize is always the number of letters in the mer -- if we
  //  are a spaced seed, it is the weight.
  //
  uint32   _merSize;
  uint32   _merSpan;

  //  The mask is used to make sure the mer has only _merSize bases
  //  set -- we can get more than that if we shift to the left.  The
  //
  uint64   _mask;

  //  For operator-=() (add a base to the left end) we need to know
  //  what the last word is, and how far to shift the bits.
  //
  uint32   _lastShift;
};





inline
void
kMerTiny::setMerSize(uint32 ms) {
  _merSize   = ms;
  _merSpan   = ms;
  _lastShift = (2 * ms - 2) % 64;
  _mask      = uint64MASK(_merSize  << 1);
}


inline
char *
kMerTiny::merToString(char *str) const {
  for (uint32 i=0; i<_merSize; i++)
    str[_merSize-i-1] = bitsToLetter[(_md >> (2*i)) & 0x03];
  str[_merSize] = 0;
  return(str);
}

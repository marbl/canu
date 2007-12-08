
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

  kMerTiny  &convertToColorSpace(kMerTiny &csmer) {

    //  Copy to the csmer.
    csmer._md         = u64bitZERO;
    csmer._merSize    = _merSize;
    csmer._merSpan    = _merSpan;
    csmer._mask       = _mask;
    csmer._lastShift = _lastShift;

    //u64bit  lastBase = 0;
    //u64bit  thisBase = 0;

    u64bit  md = csmer._md;  //  our cs mer
    u64bit  lb = 0;          //  the last base we've seen
    u64bit  tb = 0;          //  the current base

    //  Grab the first base.
    //
    md <<= (64 - _lastShift);
    lb   = (md >> 62) & u64bitMASK(2);

    //  Push that base onto the csmer.
    csmer._md  |= lb;

    //  For all the remaining bases, we need to encode
    for (u32bit i=1; i<_merSize; i++) {
      md <<= 2;
      tb   = (md >> 62) & u64bitMASK(2);

      //  Remember:
      //
      //  lb = the last base (vertical)
      //  tb = the current base (horizontal)
      //
      //      A C G T <- tb
      //    A 0 1 2 3   color =  tb
      //    C 1 0 3 2   color =  (lb ^ tb)
      //    G 2 3 0 1   color = ~(lb ^ tb)
      //    T 3 2 1 0   color = ~tb
      //
      //  How cool is that?!

      csmer._md <<= 2;

      switch (lb) {
        case 0x00:
          csmer._md |=  lb;
          break;
        case 0x01:
          csmer._md |=  (lb ^ tb);
          break;
        case 0x02:
          csmer._md |= ~(lb ^ tb);
          break;
        case 0x03:
          csmer._md |= ~lb;
          break;
        default:
          break;
      }

      lb = tb;
    }

    return(csmer);
  };

  void    clear(void) {
    _md = u64bitZERO;
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
    _md |= x & u64bitNUMBER(0x3);
  };
  void   operator-=(u64bit x) {
    *this >>= 2;
    _md |= (x & u64bitNUMBER(0x3)) << _lastShift;
  };

public:
  void     mask(bool) {
    _md &= _mask;
  };

public:
  bool    operator!=(kMerTiny const &r) const { return(_md != r._md); };
  bool    operator==(kMerTiny const &r) const { return(_md == r._md); };
  bool    operator<(kMerTiny const &r) const  { return(_md <  r._md); };
  bool    operator>(kMerTiny const &r) const  { return(_md >  r._md); };
  bool    operator<=(kMerTiny const &r) const { return(_md <= r._md); };
  bool    operator>=(kMerTiny const &r) const { return(_md >= r._md); };

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
    str[_merSize-i-1] = decompressSymbol[(_md >> (2*i)) & 0x03];
  str[_merSize] = 0;
  return(str);
};

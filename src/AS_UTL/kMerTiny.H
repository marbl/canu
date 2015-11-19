
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
 *    kmer/libbio/kmertiny.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2007-SEP-13 to 2014-APR-11
 *      are Copyright 2007-2008,2010,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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

    //  Reverse the mer

    _md = ((_md >>  2) & 0x3333333333333333llu) | ((_md <<  2) & 0xccccccccccccccccllu);
    _md = ((_md >>  4) & 0x0f0f0f0f0f0f0f0fllu) | ((_md <<  4) & 0xf0f0f0f0f0f0f0f0llu);
    _md = ((_md >>  8) & 0x00ff00ff00ff00ffllu) | ((_md <<  8) & 0xff00ff00ff00ff00llu);
    _md = ((_md >> 16) & 0x0000ffff0000ffffllu) | ((_md << 16) & 0xffff0000ffff0000llu);
    _md = ((_md >> 32) & 0x00000000ffffffffllu) | ((_md << 32) & 0xffffffff00000000llu);

    //  Complement the bases

    _md ^= 0xffffffffffffffffllu;

    //  Shift and mask out the bases not in the mer

    _md >>= 64 - _merSize * 2;
    _md  &= uint64MASK(_merSize * 2);

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
  void   writeToBitPackedFile(bitPackedFile *BPF, uint32 UNUSED(numBits)=0) const {
    BPF->putBits(_md, _merSize << 1);
  };
  void   readFromBitPackedFile(bitPackedFile *BPF, uint32 UNUSED(numBits)=0) {
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
  uint64   getWord(uint32 UNUSED(wrd)) const        { return(_md); };
  void     setWord(uint32 UNUSED(wrd), uint64 val)  { _md = val;   };

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
    str[_merSize-i-1] = alphabet.bitsToLetter((_md >> (2*i)) & 0x03);
  str[_merSize] = 0;
  return(str);
}


#if 0
//  This is not used anywhere, and is probably quite stale
inline
uint64
stringToMer(uint32 ms, char *str) {
  uint64  mer = 0L;

  for (uint32 i=0; i<ms; i++) {
    mer <<= 2;
    mer  |= compressSymbol[str[i]];
  }

  return(mer);
}
#endif

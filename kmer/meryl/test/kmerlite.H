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
    for (u32bit i=0; i<KMER_WORDS; i++)
      fprintf(stderr, "kMerLite["u32bitFMTW(2)"] = "u64bitHEX"\n", i, _wd[i]);
  };


public:
  void    copy(kMer const &that) {
    for (u32bit i=0; i<KMER_WORDS; i++)
      _wd[i] = that.getWord(i);
  };
  void    clear(void) {
    for (u32bit i=0; i<KMER_WORDS; i++)
      _wd[i] = u64bitZERO;
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

  char    *merToString(u32bit merSize, char *instr) const {
    u32bit  lastWord = merSize >> 5;
    char   *str = instr;

    if ((merSize & u32bitMASK(6)) == 0)
      lastWord++;

    //  We build the string right to left, print any partial word
    //  first, then print whole words until we run out of words to
    //  print.

    if (merSize & u32bitMASK(5)) {
      ::merToString(merSize & u32bitMASK(5), _wd[lastWord], str);
      str += merSize & u32bitMASK(5);
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
    u64bit res = u64bitZERO;
    for (u32bit i=KMER_WORDS; i--; )
      res |= _wd[i] ^ r._wd[i];
    return(res != u64bitZERO);
  };
  bool    operator==(kMerLite const &r) const {
    u64bit res = u64bitZERO;
    for (u32bit i=KMER_WORDS; i--; )
      res |= _wd[i] ^ r._wd[i];
    return(res == u64bitZERO);
  };

  bool    operator<(kMerLite const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (_wd[i] < r._wd[i])  return(true);
      if (_wd[i] > r._wd[i])  return(false);
    }
    return(false);
  };
  bool    operator>(kMerLite const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (_wd[i] > r._wd[i])  return(true);
      if (_wd[i] < r._wd[i])  return(false);
    }
    return(false);
  };
  bool    operator<=(kMerLite const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (_wd[i] < r._wd[i])  return(true);
      if (_wd[i] > r._wd[i])  return(false);
    }
    return(true);
  };
  bool    operator>=(kMerLite const &r) const {
    for (u32bit i=KMER_WORDS; i--; ) {
      if (_wd[i] > r._wd[i])  return(true);
      if (_wd[i] < r._wd[i])  return(false);
    }
    return(true);
  };
#endif

private:
  u64bit   _wd[KMER_WORDS];
};


#endif  //  KMER_LITE_H

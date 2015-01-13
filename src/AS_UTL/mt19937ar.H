#ifndef MT19937AR_H
#define MT19937AR_H

#include "AS_global.H"

//  Refactoring of
//
//    A C-program for MT19937, with initialization improved 2002/1/26.
//    Coded by Takuji Nishimura and Makoto Matsumoto.
//
//  to make it thread safe and (hopefully) more portable.
//
//  20040421, bpw

//  Refactored again, 20141208 for C++.

  static const uint32 MT_N          = 624;           //  period parameters
  static const uint32 MT_M          = 397;
  static const uint32 MT_MATRIX_A   = 0x9908b0dfUL;  //  constant vector a
  static const uint32 MT_UPPER_MASK = 0x80000000UL;  //  most significant w-r bits
  static const uint32 MT_LOWER_MASK = 0x7fffffffUL;  //  least significant r bits

class mtRandom {
private:
  void   construct(uint32 s);

public:
  mtRandom()           { construct(getpid() * time(NULL)); };
  mtRandom(uint32 s)   { construct(s);                     };
  mtRandom(uint32 *init_key, uint32 key_length);

  ~mtRandom() {
  };

  uint32     mtRandom32(void);
  uint64     mtRandom64(void)   { return((((uint64)mtRandom32()) << 32) | (uint64)mtRandom32()); }

  //  Real valued randomness
  //    mtRandomRealOpen()    -- on [0,1) real interval
  //    mtRandomRealClosed()  -- on [0,1] real interval
  //    mrRandomRealOpen53()  -- on [0,1) real interval, using 53 bits
  //
  //  "These real versions are due to Isaku Wada, 2002/01/09 added" and were taken from
  //  the mt19937ar.c distribution (but they had actual functions, not macros)
  //
  //  They also had
  //    random number in (0,1) as (mtRandom32() + 0.5) * (1.0 / 4294967296.0)
  //
  double   mtRandomRealOpen(void)   { return((double)mtRandom32() * (1.0 / 4294967296.0)); };
  double   mtRandomRealClosed(void) { return((double)mtRandom32() * (1.0 / 4294967295.0)); };
  double   mtRandomRealOpen53(void) { return(((mtRandom32() >> 5) * 67108864.0 + (mtRandom32() >> 6)) * (1.0 / 9007199254740992.0)); };

  //  returns a random number with gaussian distribution, mean of zero and std.dev. of 1
  //
  double  mtRandomGaussian(void);


private:
  uint32  mt[MT_N];  //  State vector array
  uint32  mti;       //  Ordinal of the first uninit'd element -- mti = N+1 -> elt N is uninit
  uint32  mag01[2];  //  mag01[x] = x * MT_MATRIX_A  for x=0,1
};


#endif  //  MT19937AR_H

#include "../libutil/mt19937ar/mt19937ar.h"

class randomHash {
public:
  randomHash() {
    setHashSize(40);
    _mt = mtInit(0);
  };
  ~randomHash() {};

  mt_s    *_mt;

  u32bit   _hashWidth;
  u32bit   _keys[8];
  u32bit   _random;

  void  setHashSize(u32bit width) {
    _hashWidth = width;
    resetRandomness();
  };

  //  Pick a bunch of odd shift values, that add up to 2*mersize,
  //  store them in _shiftValues.
  //
  void  resetRandomness(void) {
    //  get a new random number from the global randomness
    _random = mtRandom32(mt);
  };

  u64bit   hash(kMer km) {

    _keys[0] = km.getWord(0);
    _keys[1] = km.getWord(0) >> 32;
    _keys[2] = km.getWord(1);
    _keys[3] = km.getWord(1) >> 32;
    _keys[4] = km.getWord(2);
    _keys[5] = km.getWord(2) >> 32;
    _keys[6] = km.getWord(3);
    _keys[7] = km.getWord(3) >> 32;

    //  swiped from the mt implementation -- this is mtInit

    _mt->mt[0] = 19650218UL + _random;

    for (_mt->mti=1; _mt->mti<MT_N; _mt->mti++)
      _mt->mt[_mt->mti] = (1812433253UL * (_mt->mt[_mt->mti-1] ^ (_mt->mt[_mt->mti-1] >> 30)) + _mt->mti); 

    _mt->mag01[0] = u32bitZERO;
    _mt->mag01[1] = MT_MATRIX_A;

    //  this is mtInitArray -- '8' was the key_length

    int   i   = 1;
    int   j   = 0;
    int   k   = MT_N;

    for (; k; k--) {
      _mt->mt[i] = (_mt->mt[i] ^ ((_mt->mt[i-1] ^ (_mt->mt[i-1] >> 30)) * 1664525UL)) + _keys[j] + j; /* non linear */
      i++;
      j++;
      if (i >= MT_N) {
        _mt->mt[0] = _mt->mt[MT_N-1];
        i=1;
      }
      if (j >= 8)
        j=0;
    }
    for (k=MT_N-1; k; k--) {
      _mt->mt[i] = (_mt->mt[i] ^ ((_mt->mt[i-1] ^ (_mt->mt[i-1] >> 30)) * 1566083941UL)) - i; /* non linear */
      i++;
      if (i>=MT_N) {
        _mt->mt[0] = _mt->mt[MT_N-1];
        i=1;
      }
    }

    _mt->mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 

    //  now hash -- if this seems to work, we can rip out mtRandom32
    //  and optimize lots.  we only care about making the first two
    //  random numbers, mt itself seems to create 600 or so at once.

    return(mtRandom64(_mt) & u64bitMASK(_hashWidth));
  };
};

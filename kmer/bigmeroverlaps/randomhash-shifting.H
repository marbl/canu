
class randomHash {
public:
  randomHash(u64bit seed, u32bit merSize, u32bit width) {
    _mt      = mtInit(seed);
    _merSize = merSize;
    setHashSize(width);
  };
  ~randomHash() {
    free(_mt);
  };

  mt_s    *_mt;
  u32bit   _merSize;
  u32bit   _hashWidth;
  u32bit   _shiftValues[16];
  u32bit   _shiftMax;

  void  setHashSize(u32bit width) {
    _hashWidth = width;
    resetRandomness();
  };

  //  Pick a bunch of odd shift values, that add up to 2*mersize,
  //  store them in _shiftValues.
  //
  void  resetRandomness(void) {
    u32bit   sm = 0;

    _shiftMax = 0;

    for (u32bit i=0; i<16; i++)
      _shiftValues[i] = 0;

    for (u32bit i=0; i<16; i++) {
      _shiftValues[i] = 3 * _hashWidth / 4 + (mtRandom32(_mt) % (_hashWidth / 4));
      sm += _shiftValues[i];

      if (sm > 2 * _merSize) {
        _shiftMax = i;
        _shiftValues[i] = 0;
        break;
      }
    }
  };

  u64bit   hash(kMer k) {
    u64bit  h = u64bitZERO;

#if 0
    for (u32bit i=0; i<16; i++) {
      h  ^= k.getWord(0);
      if (_shiftValues[i] == 0)
        break;
      k >>= _shiftValues[i];
    }
#endif

#if 1
      h  ^= k.getWord(0);
      k >>= _shiftValues[0];
      h  ^= k.getWord(0);
      k >>= _shiftValues[1];
      h  ^= k.getWord(0);
      k >>= _shiftValues[2];
      h  ^= k.getWord(0);
      k >>= _shiftValues[3];
      h  ^= k.getWord(0);
      k >>= _shiftValues[4];
      h  ^= k.getWord(0);
      if (_shiftMax < 6)
        goto end;

      k >>= _shiftValues[5];
      h  ^= k.getWord(0);
      if (_shiftMax < 7)
        goto end;

      k >>= _shiftValues[6];
      h  ^= k.getWord(0);
      if (_shiftMax < 8)
        goto end;

      k >>= _shiftValues[7];
      h  ^= k.getWord(0);
      if (_shiftMax < 9)
        goto end;

      k >>= _shiftValues[8];
      h  ^= k.getWord(0);
      if (_shiftMax < 10)
        goto end;

      k >>= _shiftValues[9];
      h  ^= k.getWord(0);
      if (_shiftMax < 11)
        goto end;

      k >>= _shiftValues[10];
      h  ^= k.getWord(0);
      if (_shiftMax < 12)
        goto end;

      k >>= _shiftValues[11];
      h  ^= k.getWord(0);
      if (_shiftMax < 13)
        goto end;

      k >>= _shiftValues[12];
      h  ^= k.getWord(0);
      if (_shiftMax < 14)
        goto end;

      k >>= _shiftValues[13];
      h  ^= k.getWord(0);

  end:
#endif

    h &= u64bitMASK(_hashWidth);

    return(h);
  };
};

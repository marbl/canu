#include "positionDB.H"
#include "bio++.H"

bool
positionDB::get(u64bit   mer,
                u64bit*& posn,
                u64bit&  posnMax,
                u64bit&  posnLen) {
  u64bit  h = HASH(mer);
  u64bit  c = CHECK(mer);
  u64bit st = getDecodedValue(_hashTable, h * _hashWidth,              _hashWidth);
  u64bit ed = getDecodedValue(_hashTable, h * _hashWidth + _hashWidth, _hashWidth);
  u64bit le = ed - st;

  if (le == 0) return(false);

  u64bit  v;
  u64bit  uniqMask = u64bitONE << (_wFin - 1);
  u64bit  pos;
  u64bit  len, ptr;

  for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
    v = getDecodedValue(_buckets, J, _wFin);

    if ((v & _chckMask) == c) {
      pos = (v >> _chckWidth) & _posnPtrMask;

      if (posnMax == 0) {
        posnMax = 16384;
        try {
          posn    = new u64bit [posnMax];
        } catch (...) {
          fprintf(stderr, "positionDB::get()-- Can't allocate space for initial positions, requested "u64bitFMT" u64bit's.\n", posnMax);
          abort();
        }
      }

      if (v & uniqMask) {
        posnLen = 1;
        posn[0] = pos;
      } else {
        ptr  = pos * _posnWidth;
        len  = getDecodedValue(_positions, ptr, _posnWidth);

        if (posnMax < len) {
          posnMax = len + (len >> 2);
          try {
            posn    = new u64bit [posnMax];
          } catch (...) {
            fprintf(stderr, "positionDB::get()-- Can't allocate space for more positions, requested "u64bitFMT" u64bit's.\n", posnMax);
            abort();
          }
        }

        posnLen = 0;

        for (; len > 0; len--) {
          ptr += _posnWidth;
          posn[posnLen++] = getDecodedValue(_positions, ptr, _posnWidth);
        }
      }

      return(true);
    }
  }

  return(false);
}

bool
positionDB::exists(u64bit   mer) {
  u64bit  h = HASH(mer);
  u64bit  c = CHECK(mer);
  u64bit st = getDecodedValue(_hashTable, h * _hashWidth,              _hashWidth);
  u64bit ed = getDecodedValue(_hashTable, h * _hashWidth + _hashWidth, _hashWidth);

  if (st == ed) return(false);

  for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin)
    if ((getDecodedValue(_buckets, J, _wFin) & _chckMask) == c)
      return(true);

  return(false);
}


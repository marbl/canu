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

  //  We probably should always say the length is zero if we don't
  //  find anything, not just return false.
  //
  posnLen = 0;

  if (le == 0)
    return(false);

  //  Get this out of the way.  If we have no space, allocate a little
  //  bit.
  //
  if (posnMax == 0) {
    posnMax = 16384;
    try {
      posn    = new u64bit [posnMax];
    } catch (...) {
      fprintf(stderr, "positionDB::get()-- Can't allocate space for initial positions, requested "u64bitFMT" u64bit's.\n", posnMax);
      abort();
    }
  }

  u64bit  v;
  u64bit  pos;
  u64bit  len, ptr;

  for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
    v = getDecodedValue(_buckets, J, _wFin);

    if ((v & _chckMask) == c) {
      pos = (v >> _chckWidth) & _posnPtrMask;

      if (v & (u64bitONE << (_wFin - 1))) {
        //  We allocated space already.
        posnLen = 1;
        posn[0] = pos;
      } else {
        ptr  = pos * _posnWidth;
        len  = getDecodedValue(_positions, ptr, _posnWidth);

        //  Whoops!  Need more space!
        //
        if (posnMax < len) {
          delete [] posn;

          posnMax = len + (len >> 2);
          try {
            posn    = new u64bit [posnMax];
          } catch (...) {
            fprintf(stderr, "positionDB::get()-- Can't allocate space for more positions, requested "u64bitFMT" u64bit's.\n", posnMax);
            abort();
          }
        }

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
positionDB::exists(u64bit mer) {
  u64bit  h = HASH(mer);
  u64bit  c = CHECK(mer);
  u64bit st = getDecodedValue(_hashTable, h * _hashWidth,              _hashWidth);
  u64bit ed = getDecodedValue(_hashTable, h * _hashWidth + _hashWidth, _hashWidth);

  if (st != ed)
    for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin)
      if ((getDecodedValue(_buckets, J, _wFin) & _chckMask) == c)
        return(true);

  return(false);
}

u64bit
positionDB::count(u64bit mer) {
  u64bit  h = HASH(mer);
  u64bit  c = CHECK(mer);
  u64bit st = getDecodedValue(_hashTable, h * _hashWidth,              _hashWidth);
  u64bit ed = getDecodedValue(_hashTable, h * _hashWidth + _hashWidth, _hashWidth);

  if (st != ed)
    for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
      u64bit v = getDecodedValue(_buckets, J, _wFin);
      if ((v & _chckMask) == c) {
        if (v & (u64bitONE << (_wFin - 1)))
          return(1);
        return(getDecodedValue(_positions, ((v >> _chckWidth) & _posnPtrMask) * _posnWidth, _posnWidth));
      }
    }

  return(0);
}

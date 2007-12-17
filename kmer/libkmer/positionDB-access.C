#include "positionDB.H"
#include "bio++.H"


void
positionDB::loadPositions(u64bit   J,
                          u64bit*& posn,
                          u64bit&  posnMax,
                          u64bit&  posnLen) {

  u64bit  sizs[2] = {_pptrWidth, 1};
  u64bit  vals[2] = {0};

  getDecodedValues(_buckets, J + _chckWidth, 2, sizs, vals);

  if (vals[1]) {
    posn[posnLen++] = vals[0];
  } else {
    u64bit ptr  = vals[0] * _posnWidth;
    u64bit len  = getDecodedValue(_positions, ptr, _posnWidth);

    if (posnMax < posnLen + len) {
      delete [] posn;

      posnMax = posnLen + len + (len >> 2);
      try {
        posn    = new u64bit [posnMax];
      } catch (...) {
        fprintf(stderr, "positionDB::get()-- Can't allocate space for more positions, requested "u64bitFMT" u64bit's.\n", posnMax);
        abort();
      }
    }

    for (ptr += _posnWidth; len > 0; ptr += _posnWidth, len--)
      posn[posnLen++] = getDecodedValue(_positions, ptr, _posnWidth);
  }
}



bool
positionDB::get(u64bit   mer,
                u64bit*& posn,
                u64bit&  posnMax,
                u64bit&  posnLen) {
  u64bit  h = HASH(mer);
  u64bit  c = CHECK(mer);
  u64bit st = getDecodedValue(_hashTable, h * _hashWidth,              _hashWidth);
  u64bit ed = getDecodedValue(_hashTable, h * _hashWidth + _hashWidth, _hashWidth);

  posnLen = 0;

  if (st == ed)
    return(false);

  if (posnMax == 0) {
    posnMax = 16384;
    try {
      posn    = new u64bit [posnMax];
    } catch (...) {
      fprintf(stderr, "positionDB::get()-- Can't allocate space for initial positions, requested "u64bitFMT" u64bit's.\n", posnMax);
      abort();
    }
  }

  for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
    if (c == getDecodedValue(_buckets, J, _chckWidth)) {
      loadPositions(J, posn, posnMax, posnLen);
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

  if (st == ed)
    return(false);

  for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin)
    if (c == getDecodedValue(_buckets, J, _chckWidth))
      return(true);

  return(false);
}


u64bit
positionDB::count(u64bit mer) {
  u64bit  h = HASH(mer);
  u64bit  c = CHECK(mer);
  u64bit st = getDecodedValue(_hashTable, h * _hashWidth,              _hashWidth);
  u64bit ed = getDecodedValue(_hashTable, h * _hashWidth + _hashWidth, _hashWidth);

  if (st == ed)
    return(0);

  for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
    if (c == getDecodedValue(_buckets, J, _chckWidth)) {
      u64bit  sizs[3] = {_pptrWidth, 1, _sizeWidth};
      u64bit  vals[3] = {0};

      getDecodedValues(_buckets, J + _chckWidth, 3, sizs, vals);

      if (_sizeWidth > 0)
        return(vals[2]);

      if (vals[1])
        return(1);

      return(getDecodedValue(_positions, vals[0] * _posnWidth, _posnWidth));
    }
  }

  return(0);
}

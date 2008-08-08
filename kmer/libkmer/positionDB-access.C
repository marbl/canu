#include "bio++.H"
#include "positionDB.H"


void
positionDB::reallocateSpace(u64bit*&    posn,
                            u64bit&     posnMax,
                            u64bit&     posnLen,
                            u64bit      len) {
 
  if (posnMax < posnLen + len) {
    u64bit  *pp;

    posnMax = posnLen + len + (len >> 2);

    if (posnMax == 0)
      posnMax = 16384;

    try {
      pp = new u64bit [posnMax];
    } catch (...) {
      fprintf(stderr, "positionDB::get()-- Can't allocate space for more positions, requested "u64bitFMT" u64bit's.\n", posnMax);
      abort();
    }

    memcpy(pp, posn, sizeof(u64bit) * posnLen);

    delete [] posn;
    posn = pp;
  }
}



void
positionDB::loadPositions(u64bit   J,
                          u64bit*& posn,
                          u64bit&  posnMax,
                          u64bit&  posnLen,
                          u64bit&  count) {

  u64bit  sizs[3] = {_pptrWidth, 1, _sizeWidth};
  u64bit  vals[3] = {0, 0, 1};

  getDecodedValues(_buckets, J + _chckWidth, (_sizeWidth == 0) ? 2 : 3, sizs, vals);

  //  If the size is stored, the count is updated to the correct
  //  thing.  If it's not stored, the count is set to 1 by the default
  //  value of vals[2], and reset after we get the number of positions
  //  stored.
  //
  count = vals[2];

  if (vals[1]) {
    reallocateSpace(posn, posnMax, posnLen, 64);
    posn[posnLen++] = vals[0];
  } else {
    u64bit ptr  = vals[0] * _posnWidth;
    u64bit len  = getDecodedValue(_positions, ptr, _posnWidth);

    if (_sizeWidth == 0)
      count = len;

    reallocateSpace(posn, posnMax, posnLen, len + 64);

    for (ptr += _posnWidth; len > 0; ptr += _posnWidth, len--)
      posn[posnLen++] = getDecodedValue(_positions, ptr, _posnWidth);
  }
}



bool
positionDB::getExact(u64bit   mer,
                     u64bit*& posn,
                     u64bit&  posnMax,
                     u64bit&  posnLen,
                     u64bit&  count) {
  u64bit  h = HASH(mer);
  u64bit  c = CHECK(mer);
  u64bit st = getDecodedValue(_hashTable, h * _hashWidth,              _hashWidth);
  u64bit ed = getDecodedValue(_hashTable, h * _hashWidth + _hashWidth, _hashWidth);

  posnLen = 0;

  if (st == ed)
    return(false);

  for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
    if (c == getDecodedValue(_buckets, J, _chckWidth)) {
      loadPositions(J, posn, posnMax, posnLen, count);
      return(true);
    }
  }

  return(false);
}


bool
positionDB::existsExact(u64bit mer) {
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
positionDB::countExact(u64bit mer) {
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


u64bit
positionDB::setCount(u64bit mer, u64bit count) {
  u64bit  h = HASH(mer);
  u64bit  c = CHECK(mer);
  u64bit st = getDecodedValue(_hashTable, h * _hashWidth,              _hashWidth);
  u64bit ed = getDecodedValue(_hashTable, h * _hashWidth + _hashWidth, _hashWidth);

  if (st == ed)
    return(0);

  for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin)
    if (c == getDecodedValue(_buckets, J, _chckWidth)) {
      setDecodedValue(_buckets, J + _chckWidth + _pptrWidth + 1, _sizeWidth, count);
      return(count);
    }

  return(0);
}

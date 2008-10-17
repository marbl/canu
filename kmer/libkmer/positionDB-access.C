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
  u64bit st, ed;

  if (_hashTable_BP) {
    st = getDecodedValue(_hashTable_BP, h * _hashWidth,              _hashWidth);
    ed = getDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth);
  } else {
    st = _hashTable_FW[h];
    ed = _hashTable_FW[h+1];
  }

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
  u64bit st, ed;

  if (_hashTable_BP) {
    st = getDecodedValue(_hashTable_BP, h * _hashWidth,              _hashWidth);
    ed = getDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth);
  } else {
    st = _hashTable_FW[h];
    ed = _hashTable_FW[h+1];
  }

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
  u64bit st, ed;

  if (_hashTable_BP) {
    st = getDecodedValue(_hashTable_BP, h * _hashWidth,              _hashWidth);
    ed = getDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth);
  } else {
    st = _hashTable_FW[h];
    ed = _hashTable_FW[h+1];
  }

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
  u64bit st, ed;

  if (_hashTable_BP) {
    st = getDecodedValue(_hashTable_BP, h * _hashWidth,              _hashWidth);
    ed = getDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth);
  } else {
    st = _hashTable_FW[h];
    ed = _hashTable_FW[h+1];
  }

  if (st == ed)
    return(0);

  for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin)
    if (c == getDecodedValue(_buckets, J, _chckWidth)) {
      setDecodedValue(_buckets, J + _chckWidth + _pptrWidth + 1, _sizeWidth, count);
      return(count);
    }

  return(0);
}



void
positionDB::filter(u64bit lo,
                   u64bit hi) {
  u64bit  st=0, ed=0;  //  iteration through buckets
  u64bit        nb=0;  //  bit position of the current (read) bucket and next (write) bucket
  u64bit        np=0;  //  bit position of the current (read) position and next (write) position
  u64bit  vv;

  u64bit  loCount = 0;
  u64bit  okCount = 0;
  u64bit  hiCount = 0;

  u64bit  sizs[4] = {_chckWidth, _pptrWidth, 1, _sizeWidth};
  u64bit  vals[4] = {0, 0, 0, 0};

  //dump("posDB.before");

  fprintf(stderr, "positionDB::filter()--  Filtering out kmers less than "u64bitFMT" and more than "u64bitFMT"\n", lo, hi);

  if (_sizeWidth == 0) {
    //  Single copy mers in a table without counts can be multi-copy
    //  when combined with their reverse-complement mer.
    fprintf(stderr, "positionDB::filter()--  ERROR!\n");
    fprintf(stderr, "positionDB::filter()--  ERROR!  No count information; filtering will break canonical assumptions.\n");
    fprintf(stderr, "positionDB::filter()--  ERROR!\n");
    exit(1);
  }

  //  Grab the start of the first (current) bucket.  We reset the
  //  hashTable at the end of the loop, forcing us to keep st
  //  up-to-date, instead of grabbing it anew each iteration.
  //
  if (_hashTable_BP)
    st = getDecodedValue(_hashTable_BP, 0, _hashWidth);
  else
    st = _hashTable_FW[0];

  //  Over all buckets
  //
  for (u64bit h=0; h<_tableSizeInEntries; h++) {

    //  Grab the end of this bucket - the end is always for the
    //  current structure.  This gets reset at the end of the loop.
    //
    if (_hashTable_BP)
      ed = getDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth);
    else
      ed = _hashTable_FW[h+1];

    //  Over all entries in the bucket
    //
    while (st < ed) {
      u64bit     cb = st * _wFin;

      getDecodedValues(_buckets, cb, (_sizeWidth == 0) ? 3 : 4, sizs, vals);

      //  Argh.  Tricky.  We need to grab the count stored in the
      //  table, but if it's a single mer, there is no count.

      u64bit count = 1;            //  Real count over the whole data set
      u64bit len   = 1;            //  Number of times the mer occurs in this subset
      u64bit cp    = ~u64bitZERO;  //  current position pointer, used to copy position information

      //  If not a unique mer in this table, len and cp are defined.
      if (vals[2] == 0) {
        cp    = vals[1] * _posnWidth;
        len   = getDecodedValue(_positions, cp, _posnWidth);
        count = len;
      }

      //  The size stored in the bucket is to be believed
      if (_sizeWidth > 0)
        count = vals[3];

      //  What happened here: By default, the count is 1.  If it is
      //  NOT a unique mer in the table, we reset the count to the
      //  number of entries in the table.  Then, if there is a count
      //  stored in the table, we reset the count again.

      //  Move on to copying the data, if in the correct range.

      if (vals[2] == 1) {
        //  Is a single mer in our table.  Copy if the actual count is
        //  acceptable.
        if ((lo <= count) && (count < hi)) {
          okCount++;
          setDecodedValues(_buckets, nb, (_sizeWidth == 0) ? 3 : 4, sizs, vals);
          nb += _wFin;
        } else {
          _numberOfDistinct--;
          _numberOfMers--;
          loCount++;
        }
      } else {
        //  Mer has more than one location in the table.  Copy all
        //  locations if the count is acceptable.
        if ((lo <= count) && (count < hi)) {
          okCount++;

          //  Copy the bucket
          vals[1] = np / _posnWidth;
          setDecodedValues(_buckets, nb, (_sizeWidth == 0) ? 3 : 4, sizs, vals);
          nb += _wFin;

          //  Copy length of the positions
          if (cp != np)
            setDecodedValue(_positions, np, _posnWidth, len);
          np += _posnWidth;
          cp += _posnWidth;

          //  Copy positions
          while (len > 0) {
            if (cp != np)
              setDecodedValue(_positions, np, _posnWidth,
                              getDecodedValue(_positions, cp, _posnWidth));
            np += _posnWidth;
            cp += _posnWidth;
            len--;
          }
        } else {
          //  Not acceptable count
          _numberOfDistinct--;
          _numberOfEntries -= len;
          if (count < lo)  loCount++;
          if (count > hi)  hiCount++;
        }
      }

      //  Move to the next entry
      st++;
      cb += _wFin;
    }  //  Over all entries in the bucket

    //  Update the end position of this bucket
    if (_hashTable_BP)
      setDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth, nb / _wFin);
    else
      _hashTable_FW[h+1] = nb / _wFin;

  }  //  Over all buckets

  fprintf(stderr, "positionDB::filter()--  Filtered "u64bitFMT" kmers less than "u64bitFMT"\n", loCount, lo);
  fprintf(stderr, "positionDB::filter()--  Filtered "u64bitFMT" kmers more than "u64bitFMT"\n", hiCount, hi);
  fprintf(stderr, "positionDB::filter()--  Saved    "u64bitFMT" kmers with acceptable count\n", okCount);

  //dump("posDB.after");
}

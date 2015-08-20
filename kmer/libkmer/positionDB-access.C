
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
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-JAN-02 to 2003-AUG-14
 *      are Copyright 2003 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-APR-21 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-FEB-07 to 2014-APR-11
 *      are Copyright 2005-2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "bio++.H"
#include "positionDB.H"


void
positionDB::reallocateSpace(uint64*&    posn,
                            uint64&     posnMax,
                            uint64&     posnLen,
                            uint64      len) {

  if (posnMax < posnLen + len) {
    uint64  *pp;

    posnMax = posnLen + len + (len >> 2);

    if (posnMax == 0)
      posnMax = 16384;

    try {
      pp = new uint64 [posnMax];
    } catch (...) {
      fprintf(stderr, "positionDB::get()-- Can't allocate space for more positions, requested "uint64FMT" uint64's.\n", posnMax);
      abort();
    }

    memcpy(pp, posn, sizeof(uint64) * posnLen);

    delete [] posn;
    posn = pp;
  }
}



void
positionDB::loadPositions(uint64   J,
                          uint64*& posn,
                          uint64&  posnMax,
                          uint64&  posnLen,
                          uint64&  count) {

  uint64  sizs[3] = {_pptrWidth, 1, _sizeWidth};
  uint64  vals[3] = {0, 0, 1};

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
    uint64 ptr  = vals[0] * _posnWidth;
    uint64 len  = getDecodedValue(_positions, ptr, _posnWidth);

    if (_sizeWidth == 0)
      count = len;

    reallocateSpace(posn, posnMax, posnLen, len + 64);

    for (ptr += _posnWidth; len > 0; ptr += _posnWidth, len--)
      posn[posnLen++] = getDecodedValue(_positions, ptr, _posnWidth);
  }
}



bool
positionDB::getExact(uint64   mer,
                     uint64*& posn,
                     uint64&  posnMax,
                     uint64&  posnLen,
                     uint64&  count) {
  uint64  h = HASH(mer);
  uint64  c = CHECK(mer);
  uint64 st, ed;

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

  for (uint64 i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
    if (c == getDecodedValue(_buckets, J, _chckWidth)) {
      loadPositions(J, posn, posnMax, posnLen, count);
      return(true);
    }
  }

  return(false);
}


bool
positionDB::existsExact(uint64 mer) {
  uint64  h = HASH(mer);
  uint64  c = CHECK(mer);
  uint64 st, ed;

  if (_hashTable_BP) {
    st = getDecodedValue(_hashTable_BP, h * _hashWidth,              _hashWidth);
    ed = getDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth);
  } else {
    st = _hashTable_FW[h];
    ed = _hashTable_FW[h+1];
  }

  if (st == ed)
    return(false);

  for (uint64 i=st, J=st * _wFin; i<ed; i++, J += _wFin)
    if (c == getDecodedValue(_buckets, J, _chckWidth))
      return(true);

  return(false);
}


uint64
positionDB::countExact(uint64 mer) {
  uint64  h = HASH(mer);
  uint64  c = CHECK(mer);
  uint64 st, ed;

  if (_hashTable_BP) {
    st = getDecodedValue(_hashTable_BP, h * _hashWidth,              _hashWidth);
    ed = getDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth);
  } else {
    st = _hashTable_FW[h];
    ed = _hashTable_FW[h+1];
  }

  if (st == ed)
    return(0);

  for (uint64 i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
    if (c == getDecodedValue(_buckets, J, _chckWidth)) {
      uint64  sizs[3] = {_pptrWidth, 1, _sizeWidth};
      uint64  vals[3] = {0};

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


uint64
positionDB::setCount(uint64 mer, uint64 count) {
  uint64  h = HASH(mer);
  uint64  c = CHECK(mer);
  uint64 st, ed;

  if (_hashTable_BP) {
    st = getDecodedValue(_hashTable_BP, h * _hashWidth,              _hashWidth);
    ed = getDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth);
  } else {
    st = _hashTable_FW[h];
    ed = _hashTable_FW[h+1];
  }

  if (st == ed)
    return(0);

  for (uint64 i=st, J=st * _wFin; i<ed; i++, J += _wFin)
    if (c == getDecodedValue(_buckets, J, _chckWidth)) {
      setDecodedValue(_buckets, J + _chckWidth + _pptrWidth + 1, _sizeWidth, count);
      return(count);
    }

  return(0);
}



void
positionDB::filter(uint64 lo,
                   uint64 hi) {
  uint64  st=0, ed=0;  //  iteration through buckets
  uint64        nb=0;  //  bit position of the current (read) bucket and next (write) bucket
  uint64        np=0;  //  bit position of the current (read) position and next (write) position
  uint64  vv;

  uint64  loCount = 0;
  uint64  okCount = 0;
  uint64  hiCount = 0;

  uint64  sizs[4] = {_chckWidth, _pptrWidth, 1, _sizeWidth};
  uint64  vals[4] = {0, 0, 0, 0};

  //dump("posDB.before");

  fprintf(stderr, "positionDB::filter()--  Filtering out kmers less than "uint64FMT" and more than "uint64FMT"\n", lo, hi);

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
  for (uint64 h=0; h<_tableSizeInEntries; h++) {

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
      uint64     cb = st * _wFin;

      getDecodedValues(_buckets, cb, (_sizeWidth == 0) ? 3 : 4, sizs, vals);

      //  Argh.  Tricky.  We need to grab the count stored in the
      //  table, but if it's a single mer, there is no count.

      uint64 count = 1;            //  Real count over the whole data set
      uint64 len   = 1;            //  Number of times the mer occurs in this subset
      uint64 cp    = ~uint64ZERO;  //  current position pointer, used to copy position information

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

  fprintf(stderr, "positionDB::filter()--  Filtered "uint64FMT" kmers less than "uint64FMT"\n", loCount, lo);
  fprintf(stderr, "positionDB::filter()--  Filtered "uint64FMT" kmers more than "uint64FMT"\n", hiCount, hi);
  fprintf(stderr, "positionDB::filter()--  Saved    "uint64FMT" kmers with acceptable count\n", okCount);

  //dump("posDB.after");
}

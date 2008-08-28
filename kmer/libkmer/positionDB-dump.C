#include <stdio.h>
#include <stdlib.h>

#include "positionDB.H"
#include "bio++.H"

#ifdef DEBUGPOSDB

void
positionDB::dump(u64bit mer) {
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

  if (ed - st == 0)
    return;

  u64bit  lens[4] = {_chckWidth, _posnWidth, 1, 0};
  u64bit  vals[4] = {0};
  u64bit  nval    = 3;

  for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
    getDecodedValues(_buckets, J, nval, lens, vals);

    if (vals[0] == c) {
      if (vals[2]) {
        fprintf(stdout, "U mer="u64bitHEX" -- uniq pos="u64bitFMTW(12)" chck="u64bitHEX" posn="u64bitFMT"\n", mer, vals[0], vals[2]);
      } else {
        u64bit pos = posn * _posnWidth;
        u64bit len = getDecodedValue(_positions, pos, _posnWidth);

        fprintf(stdout, "N mer="u64bitHEX" -- dupl pos="u64bitFMTW(12)" chck="u64bitHEX" posn="u64bitFMT"\n", mer, vals[0], vals[2]);
        fprintf(stdout, "                       -- len="u64bitFMT"\n", len);
        for (pos += _posnWidth; len > 0; pos += posnWidth, len--)
          fprintf(stdout, "                       -- pos="u64bitFMT"\n", getDecodedValue(_positions, p, _posnWidth));
      }
    }
  }
}


void
positionDB::dumpTable(void) {
  u64bit  st, ed, le;
  u64bit  lens[4] = {_chckWidth, _posnWidth, 1, 0};
  u64bit  vals[4] = {0};
  u64bit  nval    = 3;

  fprintf(stdout, "Dumping "u32bitFMT" buckets.\n", _tableSizeInEntries);

  for (u32bit b=0; b<_tableSizeInEntries; b++) {
    if (_hashTable_BP) {
      st = getDecodedValue(_hashTable_BP, h * _hashWidth,              _hashWidth);
      ed = getDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth);
    } else {
      st = _hashTable_FW[h];
      ed = _hashTable_FW[h+1];
    }

    if (ed - st > 0) {
      fprintf(stdout, "Dumping bucket 0x%08x.\n", b);

      for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
        getDecodedValues(_buckets, J, nval, lens, vals);

        if (vals[1]) {
          fprintf(stdout, "U mer="u64bitHEX" -- uniq pos="u64bitFMTW(12)" chck="u64bitHEX" posn="u64bitFMT")\n", 0x0, vals[0], vals[2]);
        } else {
          u64bit pos = posn * _posnWidth;
          u64bit len = getDecodedValue(_positions, pos, _posnWidth);

          fprintf(stdout, "N mer="u64bitHEX" -- dupl pos="u64bitFMTW(12)" chck="u64bitHEX" posn="u64bitFMT"\n", 0x0, vals[0], vals[2]);
          fprintf(stdout, "                       -- len="u64bitFMT"\n", len);
          for (pos += _posnWidth; len > 0; pos += posnWidth, len--)
            fprintf(stdout, "                       -- pos="u64bitFMT"\n", getDecodedValue(_positions, p, _posnWidth));
        }
      }
    }
  }
}


#endif

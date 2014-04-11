#include <stdio.h>
#include <stdlib.h>

#include "positionDB.H"
#include "bio++.H"


void
positionDB::dump(char *name) {
  uint64  sizs[4] = {_chckWidth, _pptrWidth, 1, _sizeWidth};
  uint64  vals[4] = {0, 0, 0, 0};
  FILE   *F = fopen(name, "w");

  for (uint64 h=0; h<_tableSizeInEntries; h++) {
    uint64 st, ed;

    if (_hashTable_BP) {
      st = getDecodedValue(_hashTable_BP, h * _hashWidth,              _hashWidth);
      ed = getDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth);
    } else {
      st = _hashTable_FW[h];
      ed = _hashTable_FW[h+1];
    }

    fprintf(F, "B "uint64FMT" "uint64FMT"-"uint64FMT"\n", h, st, ed);

    while (st < ed) {
      uint64     cb = st * _wFin;

      getDecodedValues(_buckets, cb, (_sizeWidth == 0) ? 3 : 4, sizs, vals);

      fprintf(F, "%c chk="uint64HEX" pos="uint64FMT" siz="uint64FMT,
              (vals[2] == 0) ? 'D' : 'U', vals[0], vals[1], vals[3]);

      if (vals[2] == 0) {
        uint64 pos = vals[1] * _posnWidth;
        uint64 len = getDecodedValue(_positions, pos, _posnWidth);

        for (pos += _posnWidth; len > 0; pos += _posnWidth, len--)
          fprintf(F, " "uint64FMT, getDecodedValue(_positions, pos, _posnWidth));
      }

      fprintf(F, "\n");

      st++;
    }
  }

  fclose(F);
}

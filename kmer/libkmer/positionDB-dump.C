#include <stdio.h>
#include <stdlib.h>

#include "positionDB.H"
#include "bio++.H"


void
positionDB::dump(char *name) {
  u64bit  sizs[4] = {_chckWidth, _pptrWidth, 1, _sizeWidth};
  u64bit  vals[4] = {0, 0, 0, 0};
  FILE   *F = fopen(name, "w");

  for (u64bit h=0; h<_tableSizeInEntries; h++) {
    u64bit st, ed;

    if (_hashTable_BP) {
      st = getDecodedValue(_hashTable_BP, h * _hashWidth,              _hashWidth);
      ed = getDecodedValue(_hashTable_BP, h * _hashWidth + _hashWidth, _hashWidth);
    } else {
      st = _hashTable_FW[h];
      ed = _hashTable_FW[h+1];
    }

    fprintf(F, "B "u64bitFMT" "u64bitFMT"-"u64bitFMT"\n", h, st, ed);

    while (st < ed) {
      u64bit     cb = st * _wFin;

      getDecodedValues(_buckets, cb, (_sizeWidth == 0) ? 3 : 4, sizs, vals);

      fprintf(F, "%c chk="u64bitHEX" pos="u64bitFMT" siz="u64bitFMT,
              (vals[2] == 0) ? 'D' : 'U', vals[0], vals[1], vals[3]);

      if (vals[2] == 0) {
        u64bit pos = vals[1] * _posnWidth;
        u64bit len = getDecodedValue(_positions, pos, _posnWidth);

        for (pos += _posnWidth; len > 0; pos += _posnWidth, len--)
          fprintf(F, " "u64bitFMT, getDecodedValue(_positions, pos, _posnWidth));
      }

      fprintf(F, "\n");

      st++;
    }
  }

  fclose(F);
}

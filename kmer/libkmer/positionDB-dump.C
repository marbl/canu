#include <stdio.h>
#include <stdlib.h>

#include "positionDB.H"
#include "bit-packing.H"

#ifdef POSDB

#ifdef TRUE64BIT
const char *dumpUniq    = "U mer=0x%016lx -- uniq pos=%12lu chck=0x%016lx (0x%016lx)\n";
const char *dumpDupl    = "N mer=0x%016lx -- dupl pos=%12lu chck=0x%016lx (0x%016lx)\n";
const char *dumpLen     = "                       -- len=%12lu\n";
const char *dumpPos     = "                       -- pos=%12lu\n";
const char *dumpBuckets = "Dumping %lu buckets.\n";
const char *dumpBucket  = "Dumping bucket 0x%08x.\n";
#else
const char *dumpUniq    = "U mer=0x%016llx -- uniq pos=%12llu chck=0x%016llx (0x%016llx)\n";
const char *dumpDupl    = "N mer=0x%016llx -- dupl pos=%12llu chck=0x%016llx (0x%016llx)\n";
const char *dumpLen     = "                       -- len=%12llu\n";
const char *dumpPos     = "                       -- pos=%12llu\n";
const char *dumpBuckets = "Dumping %llu buckets.\n";
const char *dumpBucket  = "Dumping bucket 0x%08lx.\n";
#endif

void
positionDB::dump(u64bit mer) {
  u64bit  h = HASH(mer);
  u64bit  c = CHECK(mer);
  u64bit st = getDecodedValue(_hashTable, h * _hashWidth,              _hashWidth);
  u64bit ed = getDecodedValue(_hashTable, h * _hashWidth + _hashWidth, _hashWidth);
  u64bit le = ed - st;

  if (le == 0) return;

  u64bit  v;
  u64bit  uniqMask = u64bitONE << (_wFin - 1);  //(_chckWidth + _posnWidth);
  u64bit  posn, chck;
  u64bit  len, p;

  for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
    v = getDecodedValue(_buckets, J, _wFin);

    chck = v & _chckMask;

    if (chck == c) {
      posn = (v >> _chckWidth) & _posnMask;

      if (v & uniqMask) {
        fprintf(stdout, dumpUniq, mer, posn, chck, v);
      } else {
        fprintf(stdout, dumpDupl, mer, posn, chck, v);

        p = posn * _posnWidth;

        len = getDecodedValue(_positions, p, _posnWidth);
        fprintf(stdout, dumpLen, len);

        while (len > 0) {
          p += _posnWidth;

          posn = getDecodedValue(_positions, p, _posnWidth);
          fprintf(stdout, dumpPos, posn);

          len--;
        }
      }
    }
  }
}

void
positionDB::dumpTable(void) {
  u64bit  st, ed, le;

  fprintf(stdout, dumpBuckets, _tableSizeInEntries);

  for (u32bit b=0; b<_tableSizeInEntries; b++) {
    st = getDecodedValue(_hashTable, b * _hashWidth,              _hashWidth);
    ed = getDecodedValue(_hashTable, b * _hashWidth + _hashWidth, _hashWidth);
    le = ed - st;

    if (le > 0) {
      fprintf(stdout, dumpBucket, b);

      u64bit  v;
      u64bit  uniqMask = u64bitONE << (_wFin - 1);  //(_chckWidth + _posnWidth);
      u64bit  posn, chck;
      u64bit  len, p;

      for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {
        v = getDecodedValue(_buckets, J, _wFin);

        chck = v & _chckMask;
        posn = (v >> _chckWidth) & _posnMask;

        if (v & uniqMask) {
          fprintf(stdout, dumpUniq, 0x0, posn, chck, v);
        } else {
          fprintf(stdout, dumpDupl, 0x0, posn, chck, v);

          p = posn * _posnWidth;

          len = getDecodedValue(_positions, p, _posnWidth);
          fprintf(stdout, dumpLen, len);

          while (len > 0) {
            p += _posnWidth;

            posn = getDecodedValue(_positions, p, _posnWidth);
            fprintf(stdout, dumpPos, posn);

            len--;
          }
        }
      }
    }
  }
}


#endif

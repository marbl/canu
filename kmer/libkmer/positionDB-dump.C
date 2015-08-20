
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
 *    Brian P. Walenz from 2007-DEC-11 to 2014-APR-11
 *      are Copyright 2007-2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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


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
 *    Brian P. Walenz from 2007-NOV-11 to 2014-APR-11
 *      are Copyright 2007-2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "positionDB.H"
#include "bio++.H"


static
int
stringscmp(const void *A, const void *B) {
  uint64 const a = *(uint64 const *)A;
  uint64 const b = *(uint64 const *)B;
  if (a < b)  return(-1);
  if (a > b)  return(1);
  return(0);
}


static
uint32
makeUnique(uint64 *strings, uint32 stringsLen) {
  qsort(strings, stringsLen, sizeof(uint64), stringscmp);
  uint32  len = 0;
  uint32  nxt = 1;
  while (nxt < stringsLen) {
    if (strings[len] != strings[nxt]) {
      len++;
      strings[len] = strings[nxt];
    }
    nxt++;
  }
  return(len+1);
}


#if 0
//  debug
static
void
dumpPatterns(uint64 *strings, uint32 stringsLen, uint32 ts) {
  for (uint32 i=0; i<stringsLen; i++) {
    char   str[1024] = {0};
    uint32 cnt = 0;

    for (uint32 b=0; b<ts; b++) {
      if (strings[i] & (uint64ONE << b)) {
        str[b] = '1';
        cnt++;
      } else {
        str[b] = '0';
      }
    }

    fprintf(stdout, "%s\t"uint32FMT"\n", str, cnt);
  }
}
#endif


double
positionDB::setUpMismatchMatcher(uint32 nErrorsAllowed, uint64 approxMers) {

  //  Build an xor mask that will generate all errors for a given
  //  mersize.

  _nErrorsAllowed    = nErrorsAllowed;
  _hashedErrorsLen   = 0;
  _hashedErrorsMax   = 0;
  _hashedErrors      = 0L;

  uint32  stringsMax = 128 * 1024 * 1024;
  uint32  stringsLen = 0;
  uint64 *strings    = new uint64 [stringsMax];

  uint64  totpat = 0;
  uint64  toterr = 0;

  uint64  m1,  m2,  m3,  m4,  m5,  m6;
  uint64 *e1, *e2, *e3, *e4, *e5, *e6;

  {
    //  This can be trivially eliminated by replacing e1[x] with the err[] statement.
    uint32  ne = 3;
    for (uint32 x=1; x<_nErrorsAllowed; x++)
      ne *= 3;

    //fprintf(stderr, "Storing ne="uint32FMT" errors.\n", ne);

    e1 = new uint64 [ne];
    e2 = new uint64 [ne];
    e3 = new uint64 [ne];
    e4 = new uint64 [ne];
    e5 = new uint64 [ne];
    e6 = new uint64 [ne];

    uint64 err[3] = { 0x5555555555555555llu, 0xaaaaaaaaaaaaaaaallu, 0xffffffffffffffffllu };

    for (uint32 x=0; x<ne; x++) {
      e1[x] = err[(x/  1) % 3];
      e2[x] = err[(x/  3) % 3];
      e3[x] = err[(x/  9) % 3];
      e4[x] = err[(x/ 27) % 3];
      e5[x] = err[(x/ 81) % 3];
      e6[x] = err[(x/243) % 3];
    }
  }


  //  Zero errors
  strings[stringsLen++] = uint64ZERO;


  //  One error
  if (1 <= _nErrorsAllowed) {
    for (uint32 ai=0; ai<_merSizeInBases; ai++) {
      totpat++;
      toterr += 3;
      m1 = 0x03llu << (ai * 2);

      for (uint32 x=0; x<3; x++)
        strings[stringsLen++] = HASH((m1 & e1[x]));
    }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    //fprintf(stderr, "DONE1 totpat="uint64FMT" toterr="uint64FMT" stringsLen="uint32FMT"\n", totpat, toterr, stringsLen);
  }


  //  Two errors
  if (2 <= _nErrorsAllowed) {
    for (uint32 ai=0; ai<_merSizeInBases; ai++)
      for (uint32 bi=0; bi<ai; bi++) {
        totpat++;
        toterr += 9;
        m1 = 0x03llu << (ai * 2);
        m2 = 0x03llu << (bi * 2);

        for (uint32 x=0; x<9; x++)
          strings[stringsLen++] = HASH((m1 & e1[x]) ^ (m2 & e2[x]));
      }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    //fprintf(stderr, "DONE2 totpat="uint64FMT" toterr="uint64FMT" stringsLen="uint32FMT"\n", totpat, toterr, stringsLen);
  }


  //  Three errors
  if (3 <= _nErrorsAllowed) {
    for (uint32 ai=0; ai<_merSizeInBases; ai++)
      for (uint32 bi=0; bi<ai; bi++)
        for (uint32 ci=0; ci<bi; ci++) {
          totpat++;
          toterr += 27;
          m1 = 0x03llu << (ai * 2);
          m2 = 0x03llu << (bi * 2);
          m3 = 0x03llu << (ci * 2);

          for (uint32 x=0; x<27; x++)
            strings[stringsLen++] = HASH((m1 & e1[x]) ^ (m2 & e2[x]) ^ (m3 & e3[x]));
        }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    //fprintf(stderr, "DONE3 totpat="uint64FMT" toterr="uint64FMT" stringsLen="uint32FMT"\n", totpat, toterr, stringsLen);
  }


  //  Four errors
  if (4 <= _nErrorsAllowed) {
    for (uint32 ai=0; ai<_merSizeInBases; ai++)
      for (uint32 bi=0; bi<ai; bi++)
        for (uint32 ci=0; ci<bi; ci++)
          for (uint32 di=0; di<ci; di++) {
            totpat++;
            toterr += 81;
            m1 = 0x03llu << (ai * 2);
            m2 = 0x03llu << (bi * 2);
            m3 = 0x03llu << (ci * 2);
            m4 = 0x03llu << (di * 2);

            for (uint32 x=0; x<81; x++)
              strings[stringsLen++] = HASH((m1 & e1[x]) ^ (m2 & e2[x]) ^ (m3 & e3[x]) ^ (m4 & e4[x]));
          }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    //fprintf(stderr, "DONE4 totpat="uint64FMT" toterr="uint64FMT" stringsLen="uint32FMT"\n", totpat, toterr, stringsLen);
  }


  //  Five errors
  if (5 <= _nErrorsAllowed) {
    for (uint32 ai=0; ai<_merSizeInBases; ai++)
      for (uint32 bi=0; bi<ai; bi++)
        for (uint32 ci=0; ci<bi; ci++)
          for (uint32 di=0; di<ci; di++)
            for (uint32 ei=0; ei<di; ei++) {
              totpat++;
              toterr += 243;
              m1 = 0x03llu << (ai * 2);
              m2 = 0x03llu << (bi * 2);
              m3 = 0x03llu << (ci * 2);
              m4 = 0x03llu << (di * 2);
              m5 = 0x03llu << (ei * 2);

              if (stringsLen + 32000 >= stringsMax)
                stringsLen = makeUnique(strings, stringsLen);

              for (uint32 x=0; x<243; x++)
                strings[stringsLen++] = HASH((m1 & e1[x]) ^ (m2 & e2[x]) ^ (m3 & e3[x]) ^ (m4 & e4[x]) ^ (m5 & e5[x]));
            }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    //fprintf(stderr, "DONE5 totpat="uint64FMT" toterr="uint64FMT" stringsLen="uint32FMT"\n", totpat, toterr, stringsLen);
  }


  //  Six errors
  if (6 <= _nErrorsAllowed) {
    for (uint32 ai=0; ai<_merSizeInBases; ai++)
      for (uint32 bi=0; bi<ai; bi++)
        for (uint32 ci=0; ci<bi; ci++)
          for (uint32 di=0; di<ci; di++)
            for (uint32 ei=0; ei<di; ei++)
              for (uint32 fi=0; fi<ei; fi++) {
                totpat++;
                toterr += 729;
                m1 = 0x03llu << (ai * 2);
                m2 = 0x03llu << (bi * 2);
                m3 = 0x03llu << (ci * 2);
                m4 = 0x03llu << (di * 2);
                m5 = 0x03llu << (ei * 2);
                m6 = 0x03llu << (fi * 2);

                if (stringsLen + 32000 >= stringsMax)
                  stringsLen = makeUnique(strings, stringsLen);

                for (uint32 x=0; x<729; x++)
                  strings[stringsLen++] = HASH((m1 & e1[x]) ^ (m2 & e2[x]) ^ (m3 & e3[x]) ^ (m4 & e4[x]) ^ (m5 & e5[x]) ^ (m6 & e6[x]));
              }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    //fprintf(stderr, "DONE6 totpat="uint64FMT" toterr="uint64FMT" stringsLen="uint32FMT"\n", totpat, toterr, stringsLen);
  }


  if (7 <= _nErrorsAllowed) {
    fprintf(stderr, "Only 6 errors allowed.\n");
    exit(1);
  }

  for (uint32 i=1; i<stringsLen; i++) {
    assert((strings[i] & ~_hashMask) == 0);
    assert(strings[i] != 0);
  }

  delete [] e1;
  delete [] e2;
  delete [] e3;
  delete [] e4;
  delete [] e5;
  delete [] e6;

  delete [] _hashedErrors;

  _hashedErrorsLen = stringsLen;
  _hashedErrorsMax = stringsLen;
  _hashedErrors    = new uint64 [_hashedErrorsLen];

  memcpy(_hashedErrors, strings, sizeof(uint64) * _hashedErrorsLen);

  delete [] strings;

#ifdef UNCOMPRESS_HASH_TABLE
  //  Cost is just bucket searching.
  double work = (double)_hashedErrorsLen * approxMers / _tableSizeInEntries;
#else
  //  Cost is bucket searching + hash table lookups.
  double work = (double)_hashedErrorsLen * approxMers / _tableSizeInEntries + 2.0 * _hashedErrorsLen;
#endif

  //fprintf(stderr, "Built "uint32FMT" hashed errors at tableSize "uint32FMT" (work=%f.0).\n",
  //        _hashedErrorsLen,
  //        _tableSizeInBits,
  //        work);

  //for (uint32 i=0; i<_hashedErrorsLen; i++)
  //  fprintf(stderr, "he["uint32FMTW(5)"] = "uint64HEX"\n", i, _hashedErrors[i]);

  return(work);
}



//  Returns hits with _AT_MOST_ numMismatches mistakes.
bool
positionDB::getUpToNMismatches(uint64   mer,
                               uint32   numMismatches,
                               uint64*& posn,
                               uint64&  posnMax,
                               uint64&  posnLen) {

  PREFETCH(_hashedErrors);  //  Slightly better.

  posnLen = 0;

  if (_hashedErrors == 0L) {
    fprintf(stderr, "ERROR:  Nobody initialized getUpToNMismatches() by calling setUpMismatchMatcher().\n");
    exit(1);
  }

  if (posnMax == 0) {
    posnMax = 16384;
    try {
      posn    = new uint64 [posnMax];
    } catch (...) {
      fprintf(stderr, "positionDB::getUpToNMismatches()-- Can't allocate space for initial positions, requested "uint64FMT" uint64's.\n", posnMax);
      abort();
    }
  }

  uint64  orig = HASH(mer);

  //  Optimization that didn't work.  The idea was to compute all the
  //  hashes with errors, then sort to gain better cache locality in
  //  the lookups.  The sort dominated.
  //
  //  Another: Surprisingly, theq two getDecodedValue calls are faster
  //  than a single getDecodedValues.

  for (uint32 e=0; e<_hashedErrorsLen; e++) {
    uint64 hash = orig ^ _hashedErrors[e];
    uint64 st, ed;

    if (_hashTable_BP) {
      st = getDecodedValue(_hashTable_BP, hash * _hashWidth,              _hashWidth);
      ed = getDecodedValue(_hashTable_BP, hash * _hashWidth + _hashWidth, _hashWidth);
    } else {
      st = _hashTable_FW[hash];
      ed = _hashTable_FW[hash+1];
    }

    assert((_hashedErrors[e] & ~_hashMask) == 0);
    assert((hash             & ~_hashMask) == 0);

    //  Rebuild the mer from the hash and its check code.
    //
    //  Compare the rebuilt mer and the original mer -- if there are
    //  exactly N errors, it's a hit!  (if there are fewer than N,
    //  we'll find it when we look for N-1 errors).
    //
    //  Before rebuilding, compute diffs on the chckBits only -- if
    //  things are wildly different (the usual case) we'll get
    //  enough difference here to abort.  Remember, the chck bits
    //  are not encoded, they're an exact copy from the unhashed
    //  mer.

    if (st != ed) {
      for (uint64 i=ed-st, J=st * _wFin; i--; J += _wFin) {
        uint64 chck  = getDecodedValue(_buckets, J, _chckWidth);
        uint64 diffs = chck ^ (mer & _mask2);
        uint64 d1    = diffs & uint64NUMBER(0x5555555555555555);
        uint64 d2    = diffs & uint64NUMBER(0xaaaaaaaaaaaaaaaa);
        uint64 err   = countNumberOfSetBits64(d1 | (d2 >> 1));

        if (err <= numMismatches) {
          diffs = REBUILD(hash, chck) ^ mer;
          d1    = diffs & uint64NUMBER(0x5555555555555555);
          d2    = diffs & uint64NUMBER(0xaaaaaaaaaaaaaaaa);
          err   = countNumberOfSetBits64(d1 | (d2 >> 1));

          if (err <= numMismatches)
            //  err is junk, just need a parameter here
            loadPositions(J, posn, posnMax, posnLen, err);
        }
      }
    }
  }

  return(posnLen > 0);
}

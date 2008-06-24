#include "positionDB.H"
#include "bio++.H"

#include "hashIterator.H"

#define DEBUG_MISMATCH


static
int
stringscmp(const void *A, const void *B) {
  u64bit const a = *(u64bit const *)A;
  u64bit const b = *(u64bit const *)B;
  if (a < b)  return(-1);
  if (a > b)  return(1);
  return(0);
}

static
u32bit
makeUnique(u64bit *strings, u32bit stringsLen) {
  qsort(strings, stringsLen, sizeof(u64bit), stringscmp);
  u32bit  len = 0;
  u32bit  nxt = 1;
  while (nxt < stringsLen) {
    if (strings[len] != strings[nxt]) {
      len++;
      strings[len] = strings[nxt];
    }
    nxt++;
  }
  return(len+1);
}

static
void
dumpPatterns(u64bit *strings, u32bit stringsLen, u32bit ts) {
  for (u32bit i=0; i<stringsLen; i++) {
    char   str[1024] = {0};
    u32bit cnt = 0;

    for (u32bit b=0; b<ts; b++) {
      if (strings[i] & (u64bitONE << b)) {
        str[b] = '1';
        cnt++;
      } else {
        str[b] = '0';
      }
    }

    fprintf(stdout, "%s\t"u32bitFMT"\n", str, cnt);
  }
}


void
positionDB::setUpMismatchMatcher(u32bit nErrorsAllowed) {

  //  Build an xor mask that will generate all errors for a given
  //  mersize.

  _nErrorsAllowed    = nErrorsAllowed;
  _hashedErrorsLen   = 0;
  _hashedErrorsMax   = 0;
  _hashedErrors      = 0L;

  u32bit  stringsMax = 128 * 1024 * 1024;
  u32bit  stringsLen = 0;
  u64bit *strings    = new u64bit [stringsMax];

  u64bit  totpat = 0;
  u64bit  toterr = 0;

  u64bit  m1,  m2,  m3,  m4,  m5,  m6;
  u64bit *e1, *e2, *e3, *e4, *e5, *e6;

  {
    //  This can be trivially eliminated by replacing e1[x] with the err[] statement.
    u32bit  ne = 3;
    for (u32bit x=1; x<_nErrorsAllowed; x++)
      ne *= 3;

    fprintf(stderr, "Storing ne="u32bitFMT" errors.\n", ne);

    e1 = new u64bit [ne];
    e2 = new u64bit [ne];
    e3 = new u64bit [ne];
    e4 = new u64bit [ne];
    e5 = new u64bit [ne];
    e6 = new u64bit [ne];

    u64bit err[3] = { 0x5555555555555555llu, 0xaaaaaaaaaaaaaaaallu, 0xffffffffffffffffllu };

    for (u32bit x=0; x<ne; x++) {
      e1[x] = err[(x/  1) % 3];
      e2[x] = err[(x/  3) % 3];
      e3[x] = err[(x/  9) % 3];
      e4[x] = err[(x/ 27) % 3];
      e5[x] = err[(x/ 81) % 3];
      e6[x] = err[(x/243) % 3];
    }
  }


  //  Zero errors
  strings[stringsLen++] = u64bitZERO;


  //  One error
  if (1 <= _nErrorsAllowed) {
    for (u32bit ai=0; ai<_merSizeInBases; ai++) {
      totpat++;
      toterr += 3;
      m1 = 0x03llu << (ai * 2);

      for (u32bit x=0; x<3; x++)
        strings[stringsLen++] = HASH((m1 & e1[x]));
    }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    fprintf(stderr, "DONE1 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
  }


  //  Two errors
  if (2 <= _nErrorsAllowed) {
    for (u32bit ai=0; ai<_merSizeInBases; ai++)
      for (u32bit bi=0; bi<ai; bi++) {
        totpat++;
        toterr += 9;
        m1 = 0x03llu << (ai * 2);
        m2 = 0x03llu << (bi * 2);

        for (u32bit x=0; x<9; x++)
          strings[stringsLen++] = HASH((m1 & e1[x]) ^ (m2 & e2[x]));
      }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    fprintf(stderr, "DONE2 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
  }


  //  Three errors
  if (3 <= _nErrorsAllowed) {
    for (u32bit ai=0; ai<_merSizeInBases; ai++)
      for (u32bit bi=0; bi<ai; bi++)
        for (u32bit ci=0; ci<bi; ci++) {
          totpat++;
          toterr += 27;
          m1 = 0x03llu << (ai * 2);
          m2 = 0x03llu << (bi * 2);
          m3 = 0x03llu << (ci * 2);

          for (u32bit x=0; x<27; x++)
            strings[stringsLen++] = HASH((m1 & e1[x]) ^ (m2 & e2[x]) ^ (m3 & e3[x]));
        }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    fprintf(stderr, "DONE3 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
  }


  //  Four errors
  if (4 <= _nErrorsAllowed) {
    for (u32bit ai=0; ai<_merSizeInBases; ai++)
      for (u32bit bi=0; bi<ai; bi++)
        for (u32bit ci=0; ci<bi; ci++)
          for (u32bit di=0; di<ci; di++) {
            totpat++;
            toterr += 81;
            m1 = 0x03llu << (ai * 2);
            m2 = 0x03llu << (bi * 2);
            m3 = 0x03llu << (ci * 2);
            m4 = 0x03llu << (di * 2);
            
            for (u32bit x=0; x<81; x++)
              strings[stringsLen++] = HASH((m1 & e1[x]) ^ (m2 & e2[x]) ^ (m3 & e3[x]) ^ (m4 & e4[x]));
          }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    fprintf(stderr, "DONE4 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
  }


  //  Five errors
  if (5 <= _nErrorsAllowed) {
    for (u32bit ai=0; ai<_merSizeInBases; ai++)
      for (u32bit bi=0; bi<ai; bi++)
        for (u32bit ci=0; ci<bi; ci++)
          for (u32bit di=0; di<ci; di++)
            for (u32bit ei=0; ei<di; ei++) {
              totpat++;
              toterr += 243;
              m1 = 0x03llu << (ai * 2);
              m2 = 0x03llu << (bi * 2);
              m3 = 0x03llu << (ci * 2);
              m4 = 0x03llu << (di * 2);
              m5 = 0x03llu << (ei * 2);

              if (stringsLen + 32000 >= stringsMax) {
                stringsLen = makeUnique(strings, stringsLen);
                fprintf(stderr, "INTR5 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
              }

              for (u32bit x=0; x<243; x++)
                strings[stringsLen++] = HASH((m1 & e1[x]) ^ (m2 & e2[x]) ^ (m3 & e3[x]) ^ (m4 & e4[x]) ^ (m5 & e5[x]));
            }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    fprintf(stderr, "DONE5 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
  }


  //  Six errors
  if (6 <= _nErrorsAllowed) {
    for (u32bit ai=0; ai<_merSizeInBases; ai++)
      for (u32bit bi=0; bi<ai; bi++)
        for (u32bit ci=0; ci<bi; ci++)
          for (u32bit di=0; di<ci; di++)
            for (u32bit ei=0; ei<di; ei++)
              for (u32bit fi=0; fi<ei; fi++) {
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

                for (u32bit x=0; x<729; x++)
                  strings[stringsLen++] = HASH((m1 & e1[x]) ^ (m2 & e2[x]) ^ (m3 & e3[x]) ^ (m4 & e4[x]) ^ (m5 & e5[x]) ^ (m6 & e6[x]));
              }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, _tableSizeInBits);
    fprintf(stderr, "DONE6 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
  }



  for (u32bit i=0; i<stringsLen; i++)
    if (strings[i] == 0)
      fprintf(stderr, "ZERO at i="u32bitFMT"\n", i);

  delete [] e1;
  delete [] e2;
  delete [] e3;
  delete [] e4;
  delete [] e5;
  delete [] e6;

  _hashedErrorsLen = stringsLen;
  _hashedErrorsMax = stringsMax;
  _hashedErrors    = new u64bit [_hashedErrorsLen];

  memcpy(_hashedErrors, strings, sizeof(u64bit) * _hashedErrorsLen);

  delete [] strings;
}







//  Returns hits with _AT_MOST_ numMismatches mistakes.
bool
positionDB::getUpToNMismatches(u64bit   mer,
                               u64bit   numMismatches,
                               u64bit*& posn,
                               u64bit&  posnMax,
                               u64bit&  posnLen) {

  posnLen = 0;

  if (_hashedErrors == 0L) {
    fprintf(stderr, "ERROR:  Nobody initialized getUpToNMismatches() by calling setUpMismatchMatcher().\n");
    exit(1);
  }

  if (posnMax == 0) {
    posnMax = 16384;
    try {
      posn    = new u64bit [posnMax];
    } catch (...) {
      fprintf(stderr, "positionDB::get()-- Can't allocate space for initial positions, requested "u64bitFMT" u64bit's.\n", posnMax);
      abort();
    }
  }

  u64bit  hash = u64bitZERO;
  u64bit  orig = HASH(mer);
  u64bit  st;
  u64bit  ed;

  for (u32bit e=0; e<_hashedErrorsLen; e++) {
    hash = orig ^ _hashedErrors[e];
    st   = getDecodedValue(_hashTable, hash * _hashWidth,              _hashWidth);
    ed   = getDecodedValue(_hashTable, hash * _hashWidth + _hashWidth, _hashWidth);

    for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {

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
      
      u64bit  chck  = getDecodedValue(_buckets, J, _chckWidth);
      u64bit  chk   = chck;
      u64bit  diffs = chk ^ (mer & _mask2);
      u64bit  d1    = diffs & u64bitNUMBER(0x5555555555555555);
      u64bit  d2    = diffs & u64bitNUMBER(0xaaaaaaaaaaaaaaaa);
      u32bit  err   = countNumberOfSetBits64(d1 | (d2 >> 1));

      //fprintf(stderr, "mer "u64bitHEX" -- chk "u64bitHEX" (partial)\n", mer & _mask2, chk);

      if (err <= numMismatches) {
        chk   = REBUILD(hash, chck);
        diffs = chk ^ mer;
        d1    = diffs & u64bitNUMBER(0x5555555555555555);
        d2    = diffs & u64bitNUMBER(0xaaaaaaaaaaaaaaaa);
        err = countNumberOfSetBits64(d1 | (d2 >> 1));

        fprintf(stderr, "mer "u64bitHEX" -- chk "u64bitHEX" (full)\n", mer, chk);

        if (err <= numMismatches) {
          fprintf(stderr, "mer:"u64bitHEX" -- "u64bitFMT" MISMATCHES\n", mer, err);
          u64bit c = 0;
          loadPositions(J, posn, posnMax, posnLen, c);
        }
      }
    }
  }

  return(posnLen > 0);
}

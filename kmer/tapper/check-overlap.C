#include <stdio.h>
#include <stdlib.h>

#include "../libutil/util++.H"

//  Usage:
//
//  assume we have a full word of each error (4 words) -- 00....; 01....; 10....; 11....
//  assume we have a the mask where each error goes -- 001100000011 - error at 4 and 0, no error at 5, 3, 2, 1.
//
//  construct the new hash as
//
//  for ep1=0; ep1<4; ep1++
//    for ep2=0; ep2<4; ep2++
//      for ep3=0; ep3<4; ep3++
//        hash = hash ^ (em1 & ep1) ^ (em2 & ep2) ^ (em3 & ep3);
//
//  For one error, there are 24 error masks.  4 choices of error, but
//  in this case, we can trivially short circuit when the error agrees
//  with the original mer, so 24 * 3 = 72 lookups.  Brute force is the
//  same.
//
//  For exactly two errors, there are 300 error masks, but only 248
//  are distinct.  It's harder to figure out which ones we can short
//  circuit, so assume we have to check all four errors for each; 16
//  total.  that's 248 * 16 = 2368 lookups.
//
//  Compare against the brute force lookup of 24*3 * 23*3 = 4968
//
//  Three errors, 2324 masks, 1364 distinct; 64 error patterns, so
//  87296 lookups.  Brute force, 327888.
//
//  Optimize.  Build all the hashes, then sort.  Throw out dups.
//
//

//  OLD
//
//  if more than one error hits a given hashed location, we cannot short circuit there.
//
//  ----a-c---b
//  -----
//     -----
//        -----
//
//  Error a and b both land in the last hashed position, we can
//  short if both a and b are unchanged.  Error c is isolated, even
//  though it is in two hashed spots.  We can short if c is
//  unchanged.
//
//  if one error hits more than one hashed location, they are
//  linked.  In the above example, c hits hashed positions 0 and 3.
//  When iterating over errors in the hash (0, 3, 4) we need to
//  remember that 0 and 3 are the same error.
//

//  OLD
//
//  ARGH!  So, one error in the mer can show up in the hash in up to
//  three places.
//
//  If our mer is 25 bp long, and we toss out the first letter because
//  it's biased by the letter before it (the adapter), we have a 24
//  mer.  The posDB wants to give us 22 tblBits, leaving shifts of 12
//  and 24 bits.
//
//     [-----------1-------------]
//           [-----------1-------------]
//                 [------------1------------]
//
//  The two shifted mers (second and third lines) both contribute the
//  error to the hash.  The first line barely misses contributing.
//
//  Can three ever occur?


int
stringscmp(const void *A, const void *B) {
  u64bit a = *(u64bit *)A;
  u64bit b = *(u64bit *)B;
  if (a < b)  return(-1);
  if (a > b)  return(1);
  return(0);
}

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

#define HASH(m) ((((m) >> s1) | ((m) >> s2) | (m)) & mask)


int
main(int argc, char **argv) {

  if (argc < 3) {
    fprintf(stderr, "usage: %s merSizeInBases tableSizeInBits\n", argv[0]);
    exit(1);
  }

  //  Input values
  u32bit   ms   = atoi(argv[1]);           //  BASES
  u32bit   ts   = atoi(argv[2]);           //  BITS
  u32bit   s1   = 2 * ms - ts;  //  merSizeInBits - tableSizeInBits
  u32bit   s2   = s1 / 2;       //  shift1 / 2
  u32bit   nerr = 6;

  if (((ts % 2) == 1) || ((s1 % 2) == 1) || ((s2 % 2) == 1)) {
    //  We can't do odd values.  It breaks the hash iteration scheme.
    fprintf(stderr, "Inputs wrong.  Something is odd.\n");
    exit(1);
  }

  //ts /= 2;  //  Now in bases
  //s1 /= 2;
  //s2 /= 2;

  u32bit  stringsMax = 128 * 1024 * 1024;
  u32bit  stringsLen = 0;
  u64bit *strings    = new u64bit [stringsMax];

  //memset(strings, 0, sizeof(u64bit) * stringsMax);

  u64bit   mask = u64bitMASK(ts);

  fprintf(stderr, "stringsMax="u32bitFMT"\n", stringsMax);
  fprintf(stderr, "ms="u32bitFMT" ts="u32bitFMT" s1="u32bitFMT" s2="u32bitFMT" mask="u64bitHEX"\n", ms, ts, s1, s2, mask);

  //  For every possible placement in the string construct a mask of where we can change the hash

  u64bit  totpat = 0;
  u64bit  toterr = 0;

  u64bit  m1,  m2,  m3,  m4,  m5,  m6;
  u64bit *e1, *e2, *e3, *e4, *e5, *e6;

  {
    //  This can be trivially eliminated by replacing e1[x] with the err[] statement.
    u32bit  ne = 3;
    for (u32bit x=1; x<nerr; x++)
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

  //  One error
  if (1 <= nerr) {
    for (u32bit ai=0; ai<ms; ai++) {
      totpat++;
      toterr += 3;
      m1 = 0x03llu << (ai * 2);

      for (u32bit x=0; x<3; x++)
        strings[stringsLen++] = HASH((m1 & e1[x]));
    }

    stringsLen = makeUnique(strings, stringsLen);
    stringsLen = makeUnique(strings, stringsLen);
    //dumpPatterns(strings, stringsLen, ts);
    fprintf(stderr, "DONE1 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
  }


  //  Two errors
  if (2 <= nerr) {
    for (u32bit ai=0; ai<ms; ai++)
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
    //dumpPatterns(strings, stringsLen, ts);
    fprintf(stderr, "DONE2 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
  }


  //  Three errors
  if (3 <= nerr) {
    for (u32bit ai=0; ai<ms; ai++)
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
    //dumpPatterns(strings, stringsLen, ts);
    fprintf(stderr, "DONE3 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
  }


  //  Four errors
  if (4 <= nerr) {
    for (u32bit ai=0; ai<ms; ai++)
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
    //dumpPatterns(strings, stringsLen, ts);
    fprintf(stderr, "DONE4 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
  }


  //  Five errors
  if (5 <= nerr) {
    for (u32bit ai=0; ai<ms; ai++)
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
    //dumpPatterns(strings, stringsLen, ts);
    fprintf(stderr, "DONE5 totpat="u32bitFMT" toterr="u32bitFMT" stringsLen="u32bitFMT"\n", totpat, toterr, stringsLen);
  }


  //  Six errors
  if (6 <= nerr) {
    for (u32bit ai=0; ai<ms; ai++)
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
    //dumpPatterns(strings, stringsLen, ts);
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
}

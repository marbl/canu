#include <stdio.h>
#include <stdlib.h>

#include "../libutil/util++.H"


//  convert to ints -- we're only storing the hashed thing, and that's going to be 20-28 bits.


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


int
main(int argc, char **argv) {

  //  Input values
  u32bit   ms   = 24;  //  BASES
  u32bit   ts   = 32;  //  BITS
  u32bit   s1   = 2 * ms - ts;
  u32bit   s2   = s1 / 2;
  u64bit   mask = u64bitMASK(ts / 2);
  u32bit   nerr = 3;

  if (((ts % 2) == 1) || ((s1 % 2) == 1) || ((s2 % 2) == 1)) {
    //  We can't do odd values.  It breaks the hash iteration scheme.
    fprintf(stderr, "Inputs wrong.  Something is odd.\n");
    exit(1);
  }

  ts /= 2;  //  Now in bases
  s1 /= 2;
  s2 /= 2;

  u32bit  stringsMax = ms * ms * ms * ms * ms + 128 * 1024;
  u32bit  stringsLen = 0;
  u64bit *strings    = new u64bit [stringsMax];

  memset(strings, 0, sizeof(u64bit) * stringsMax);

  //  For every possible placement in the string construct a mask of where we can change the hash

  u64bit  mer;
  u64bit  has;

  //  One error
  if (1 <= nerr)
    for (u32bit ai=0; ai<ms; ai++) {
      u64bit mer = u64bitONE << ai;
      strings[stringsLen++] = ((mer >> s1) | (mer >> s2) | (mer)) & mask;
    }

  //  Two errors
  if (2 <= nerr)
    for (u32bit ai=0; ai<ms; ai++)
      for (u32bit bi=0; bi<ai; bi++) {
        u64bit mer = (u64bitONE << ai) | (u64bitONE << bi);
        strings[stringsLen++] = ((mer >> s1) | (mer >> s2) | (mer)) & mask;
      }

  //  Three errors
  if (3 <= nerr)
    for (u32bit ai=0; ai<ms; ai++)
      for (u32bit bi=0; bi<ai; bi++)
        for (u32bit ci=0; ci<bi; ci++) {
          u64bit mer = (u64bitONE << ai) | (u64bitONE << bi) | (u64bitONE << ci);
          strings[stringsLen++] = ((mer >> s1) | (mer >> s2) | (mer)) & mask;
        }

  //  Four errors
  if (4 <= nerr)
    for (u32bit ai=0; ai<ms; ai++)
      for (u32bit bi=0; bi<ai; bi++)
        for (u32bit ci=0; ci<bi; ci++)
          for (u32bit di=0; di<ci; di++) {
            u64bit mer = (u64bitONE << ai) | (u64bitONE << bi) | (u64bitONE << ci) | (u64bitONE << di);
            strings[stringsLen++] = ((mer >> s1) | (mer >> s2) | (mer)) & mask;
          }

  //  Five errors
  if (5 <= nerr)
    for (u32bit ai=0; ai<ms; ai++)
      for (u32bit bi=0; bi<ai; bi++)
        for (u32bit ci=0; ci<bi; ci++)
          for (u32bit di=0; di<ci; di++)
            for (u32bit ei=0; ei<di; ei++) {
              u64bit mer = (u64bitONE << ai) | (u64bitONE << bi) | (u64bitONE << ci) | (u64bitONE << di) | (u64bitONE << ei);
              strings[stringsLen++] = ((mer >> s1) | (mer >> s2) | (mer)) & mask;
            }

  //  Six errors
  if (6 <= nerr)
    for (u32bit ai=0; ai<ms; ai++)
      for (u32bit bi=0; bi<ai; bi++)
        for (u32bit ci=0; ci<bi; ci++)
          for (u32bit di=0; di<ci; di++)
            for (u32bit ei=0; ei<di; ei++)
              for (u32bit fi=0; fi<ei; fi++) {
                u64bit mer = (u64bitONE << ai) | (u64bitONE << bi) | (u64bitONE << ci) | (u64bitONE << di) | (u64bitONE << ei) | (u64bitONE << fi);
                strings[stringsLen++] = ((mer >> s1) | (mer >> s2) | (mer)) & mask;
              }

  for (u32bit i=0; i<stringsLen; i++)
    if (strings[i] == 0)
      fprintf(stderr, "ZERO at i="u32bitFMT"\n", i);

  fprintf(stderr, "DONE1 stringsLen = "u32bitFMT"\n", stringsLen);
  stringsLen = makeUnique(strings, stringsLen);
  fprintf(stderr, "DONE2 stringsLen = "u32bitFMT"\n", stringsLen);
  stringsLen = makeUnique(strings, stringsLen);
  fprintf(stderr, "DONE3 stringsLen = "u32bitFMT"\n", stringsLen);


  //  Output the patterns.

#if 0
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
#endif
}

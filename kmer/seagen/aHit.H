#ifndef AHIT_H
#define AHIT_H

#include "libbri.H"

//  $Id$

struct aHit {
  u32bit  _direction;
  u32bit  _qsIdx;
  u32bit  _dsIdx;
  u32bit  _dsLo;
  u32bit  _dsHi;
  u32bit  _covered;
  u32bit  _matched;
  u32bit  _numMers;
};


void   ahit_writeBinary(aHit *a, FILE *F) {
  fwrite(a, sizeof(aHit), 1, F);
}

void   ahit_readBinary(aHit *a, FILE *F) {
  fread(a, sizeof(aHit), 1, F);
}

void   ahit_printASCII(aHit *a, FILE *F) {
#ifdef TRUE64BIT
  fprintf(F, "-%c -e %u -D %u %u %u -M %u %u %u\n",
          a->_direction ? 'f' : 'r',
          a->_qsIdx,
          a->_dsIdx,
          a->_dsLo,
          a->_dsHi,
          a->_covered,
          a->_matched,
          a->_numMers);
#else
  fprintf(F, "-%c -e %lu -D %lu %lu %lu -M %lu %lu %lu\n",
          a->_direction ? 'f' : 'r',
          a->_qsIdx,
          a->_dsIdx,
          a->_dsLo,
          a->_dsHi,
          a->_covered,
          a->_matched,
          a->_numMers);
#endif
}


//  We don't read the string here so that we can use a static buffer
//  in whatever loop we read with.
//  
//  e.g.,
//
//  char    b[1025];
//  while (!feof(I)) {
//    fgets(b, 1024, I);
//    if (!feof(I))
//      ahit_parseString(a, b);
//  }
//
//  Note that using sscanf, while easy to implement, and safe
//  (looking, anyways), is terribly slow, and not really that safe.
//
//  char    c;
//  sscanf(b, "-%c -e %d -D %d %d %d -M %d %d %d",
//         &c,
//         &a->_qsIdx,
//         &a->_dsIdx,
//         &a->_dsLo,
//         &a->_dsHi,
//         &a->_covered,
//         &a->_matched,
//         &a->_numMers);
//  a->_direction = (c == 'f');
//
//  fast: 138.440u 40.500s  4:04.61 73.1% 0+2k 327822+0io  4pf+0w
//  slow: 737.587u 38.652s 13:12.42 97.9% 0+2k 328006+0io 11pf+0w
//
void   ahit_parseString(aHit *a, char *b) {
  char *c = b+1;

  a->_direction = (*c == 'f');
  c += 1;

  if (c[2] != 'e')  fprintf(stderr, "'%s' didn't get -e\n", b);

  c += 4;
  a->_qsIdx     = (u32bit)strtoul(c, &c, 10);

  if (c[2] != 'D')  fprintf(stderr, "'%s' didn't get -D\n", b);

  c += 4;
  a->_dsIdx     = (u32bit)strtoul(c, &c, 10);
  a->_dsLo      = (u32bit)strtoul(c, &c, 10);
  a->_dsHi      = (u32bit)strtoul(c, &c, 10);

  if (c[2] == 'M') {
    c += 4;
    a->_covered   = (u32bit)strtoul(c, &c, 10);
    a->_matched   = (u32bit)strtoul(c, &c, 10);
    a->_numMers   = (u32bit)strtoul(c, &c, 10);
  } else {
    //fprintf(stderr, "'%s' didn't get -M\n", b);
    a->_covered   = 0;
    a->_matched   = 0;
    a->_numMers   = 0;
  }
}



#endif  //  AHIT_H


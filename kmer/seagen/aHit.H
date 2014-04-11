#ifndef AHIT_H
#define AHIT_H

#include "bio.h"
#include "util++.H"
#include <stdio.h>

struct aHit {
  uint32  _forward : 1;
  uint32  _merged : 1;
  uint32  _qsIdx;
  uint32  _dsIdx;
  uint32  _dsLo;
  uint32  _dsHi;
  uint32  _covered;
  uint32  _matched;
  uint32  _numMers;
};


void   ahit_writeBinary(aHit *a, FILE *F);
void   ahit_readBinary(aHit *a, FILE *F);

void   ahit_readBinary(aHit *a, readBuffer *F);

void   ahit_printASCII(aHit *a, FILE *F);
void   ahit_parseString(aHit *a, char *b);

#endif  //  AHIT_H


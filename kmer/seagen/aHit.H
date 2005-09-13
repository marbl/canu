#ifndef AHIT_H
#define AHIT_H

#include "bio.h"
#include "util++.H"
#include <stdio.h>

struct aHit {
  u32bit  _forward : 1;
  u32bit  _merged : 1;
  u32bit  _qsIdx;
  u32bit  _dsIdx;
  u32bit  _dsLo;
  u32bit  _dsHi;
  u32bit  _covered;
  u32bit  _matched;
  u32bit  _numMers;
};


void   ahit_writeBinary(aHit *a, FILE *F);
void   ahit_readBinary(aHit *a, FILE *F);

void   ahit_readBinary(aHit *a, readBuffer *F);

void   ahit_printASCII(aHit *a, FILE *F);
void   ahit_parseString(aHit *a, char *b);

#endif  //  AHIT_H


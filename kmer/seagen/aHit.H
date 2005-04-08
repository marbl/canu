#ifndef AHIT_H
#define AHIT_H

#include "bio.h"
#include "util++.H"
#include <stdio.h>

//  The field called _flgas used to be _direction and held just the
//  orientation of the match.  So that we can use sortHits for both
//  seatac and searchGENOME, this was extended to a _flags field.
//
//  FLAG_DIRECTION -- if set, the match is forward, otherwise reverse.

#define FLAG_DIRECTION    0x00000001

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


void   ahit_writeBinary(aHit *a, FILE *F);
void   ahit_readBinary(aHit *a, FILE *F);

void   ahit_readBinary(aHit *a, readBuffer *F);

void   ahit_printASCII(aHit *a, FILE *F);
void   ahit_parseString(aHit *a, char *b);

#endif  //  AHIT_H


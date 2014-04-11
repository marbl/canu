#ifndef MT19937AR_H
#define MT19937AR_H

//  Refactoring of
//
//    A C-program for MT19937, with initialization improved 2002/1/26.
//    Coded by Takuji Nishimura and Makoto Matsumoto.
//
//  to make it thread safe and (hopefully) more portable.
//
//  20040421, bpw

//  bri.h contains the function prototypes, but we hide the structure and
//  implementation here.
//
#include "../util.h"

/* Period parameters */  
#define MT_N 624
#define MT_M 397
#define MT_MATRIX_A   0x9908b0dfUL /* constant vector a */
#define MT_UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define MT_LOWER_MASK 0x7fffffffUL /* least significant r bits */


struct mtctx {
  //  The array for the state vector
  //
  uint32  mt[MT_N];

  //  The ordinal of the first uninitialized element --
  //  mti = N+1 -> element N is uninitialized
  //
  uint32  mti;

  //  Something
  //  mag01[x] = x * MT_MATRIX_A  for x=0,1
  //
  uint32  mag01[2];
};

//  This is declared in util.h
//
//typedef struct mt mt_s;


#endif  //  MT19937AR_H

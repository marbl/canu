#ifndef MT19937AR_H
#define MT19937AR_H

//  Refactoring of
//
//    A C-program for MT19937, with initialization improved 2002/1/26.
//    Coded by Takuji Nishimura and Makoto Matsumoto.
//
//  to make it thread safe and (hopefully) more portable.


#include "bri.h"



/* Period parameters */  
#define MT_N 624
#define MT_M 397
#define MT_MATRIX_A   0x9908b0dfUL /* constant vector a */
#define MT_UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define MT_LOWER_MASK 0x7fffffffUL /* least significant r bits */


typedef struct {
  //  The array for the state vector
  //
  u32bit  mt[MT_N];

  //  The ordinal of the first uninitialized element --
  //  mti = N+1 -> element N is uninitialized
  //
  u32bit  mti;

  //  Something
  //  mag01[x] = x * MT_MATRIX_A  for x=0,1
  //
  u32bit  mag01[2];
} mt_s;


mt_s          *init_genrand(u32bit s);
mt_s          *init_by_array(u32bit *init_key, u32bit key_length);
u32bit         genrand_int32(mt_s *mt);
s32bit         genrand_int31(mt_s *mt);
double         genrand_real1(mt_s *mt);
double         genrand_real2(mt_s *mt);
double         genrand_real3(mt_s *mt);
double         genrand_res53(mt_s *mt);



/* generates a random number on [0,0x7fffffff]-interval */
inline
s32bit
genrand_int31(mt_s *mt) {
  return (s32bit)(genrand_int32(mt) >> 1);
}

/* generates a random number on [0,1]-real-interval */
inline
double
genrand_real1(mt_s *mt) {
    return genrand_int32(mt)*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
inline
double
genrand_real2(mt_s *mt) {
    return genrand_int32(mt)*(1.0/4294967296.0); 
}

/* generates a random number on (0,1)-real-interval */
inline
double
genrand_real3(mt_s *mt) {
    return (((double)genrand_int32(mt)) + 0.5)*(1.0/4294967296.0); 
}

/* generates a random number on [0,1) with 53-bit resolution*/
inline
double
genrand_res53(mt_s *mt) { 
    u32bit a = genrand_int32(mt) >> 5;
    u32bit b = genrand_int32(mt) >> 6; 
    return(a * 67108864.0 + b) * (1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */



#endif  //  MT19937AR_H

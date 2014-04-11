#ifndef ENDIANESS_H
#define ENDIANESS_H

#include <stdio.h>
#include <stdlib.h>

#include "util.h"

//  We need to test how to swap off_t and size_t

//  See also test/endianess.c

//  If we wanted to convert to network order for everything, rather
//  than convert only when needed, this would be useful.
//
#if 0
bool
checkEndianessSwapNeeded(void) {

  union u64 {
    uint64          u;
    unsigned char   c[8];
  };
  union u32 {
    uint32          u;
    unsigned char   c[4];
  };
  union u16 {
    uint16          u;
    unsigned char   c[2];
  };

  u64 u64t.u = uint64NUMBER(0x0123456789abcdef);

  return(u64t.c[0] != 0x0f)
}
#endif


inline
uint64
uint64Swap(uint64 x) {
  x = ((x >>  8) & uint64NUMBER(0x00ff00ff00ff00ff)) | ((x <<  8) & uint64NUMBER(0xff00ff00ff00ff00));
  x = ((x >> 16) & uint64NUMBER(0x0000ffff0000ffff)) | ((x << 16) & uint64NUMBER(0xffff0000ffff0000));
  x = ((x >> 32) & uint64NUMBER(0x00000000ffffffff)) | ((x << 32) & uint64NUMBER(0xffffffff00000000));
  return(x);
}

inline
uint32
uint32Swap(uint32 x) {
  x = ((x >>  8) & uint32NUMBER(0x00ff00ff)) | ((x <<  8) & uint32NUMBER(0xff00ff00));
  x = ((x >> 16) & uint32NUMBER(0x0000ffff)) | ((x << 16) & uint32NUMBER(0xffff0000));
  return(x);
}

inline
uint16
uint16Swap(uint16 x) {
  x = ((x >>  8) & 0x00ff) | ((x <<  8) & 0xff00);
  return(x);
}

#endif  //  ENDIANESS_H

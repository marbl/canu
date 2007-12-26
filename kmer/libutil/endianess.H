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
    u64bit          u;
    unsigned char   c[8];
  };
  union u32 {
    u32bit          u;
    unsigned char   c[4];
  };
  union u16 {
    u16bit          u;
    unsigned char   c[2];
  };

  u64 u64t.u = u64bitNUMBER(0x0123456789abcdef);

  return(u64t.c[0] != 0x0f)
}
#endif


inline
u64bit
u64bitSwap(u64bit x) {
  x = ((x >>  8) & u64bitNUMBER(0x00ff00ff00ff00ff)) | ((x <<  8) & u64bitNUMBER(0xff00ff00ff00ff00));
  x = ((x >> 16) & u64bitNUMBER(0x0000ffff0000ffff)) | ((x << 16) & u64bitNUMBER(0xffff0000ffff0000));
  x = ((x >> 32) & u64bitNUMBER(0x00000000ffffffff)) | ((x << 32) & u64bitNUMBER(0xffffffff00000000));
  return(x);
}

inline
u32bit
u32bitSwap(u32bit x) {
  x = ((x >>  8) & u32bitNUMBER(0x00ff00ff)) | ((x <<  8) & u32bitNUMBER(0xff00ff00));
  x = ((x >> 16) & u32bitNUMBER(0x0000ffff)) | ((x << 16) & u32bitNUMBER(0xffff0000));
  return(x);
}

inline
u16bit
u16bitSwap(u16bit x) {
  x = ((x >>  8) & 0x00ff) | ((x <<  8) & 0xff00);
  return(x);
}

#endif  //  ENDIANESS_H

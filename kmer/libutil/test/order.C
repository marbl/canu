#include <stdio.h>
#include <stdlib.h>

#include "../util.h"

//#include <sys/param.h>

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


u64bit
u64bitSwap(u64bit x) {
  x = ((x >>  8) & u64bitNUMBER(0x00ff00ff00ff00ff)) | ((x <<  8) & u64bitNUMBER(0xff00ff00ff00ff00));
  x = ((x >> 16) & u64bitNUMBER(0x0000ffff0000ffff)) | ((x << 16) & u64bitNUMBER(0xffff0000ffff0000));
  x = ((x >> 32) & u64bitNUMBER(0x00000000ffffffff)) | ((x << 32) & u64bitNUMBER(0xffffffff00000000));
  return(x);
}

u32bit
u32bitSwap(u32bit x) {
  x = ((x >>  8) & u32bitNUMBER(0x00ff00ff)) | ((x <<  8) & u32bitNUMBER(0xff00ff00));
  x = ((x >> 16) & u32bitNUMBER(0x0000ffff)) | ((x << 16) & u32bitNUMBER(0xffff0000));
  return(x);
}

u16bit
u16bitSwap(u16bit x) {
  x = ((x >>  8) & 0x00ff) | ((x <<  8) & 0xff00);
  return(x);
}



int
main(int argc, char **argv) {
  u64  u64v;
  u32  u32v;
  u16  u16v;

  u64v.u = 0x1234567890abcdefLLU;
  u32v.u = 0x12345678;
  u16v.u = 0x1234;

  for (int i=0; i<8; i++)
    fprintf(stderr, "%02x", u64v.c[i]);
  fprintf(stderr, "\n");

  for (int i=0; i<4; i++)
    fprintf(stderr, "%02x", u32v.c[i]);
  fprintf(stderr, "\n");

  for (int i=0; i<2; i++)
    fprintf(stderr, "%02x", u16v.c[i]);
  fprintf(stderr, "\n");

  u64v.u = u64bitSwap(u64v.u);
  u32v.u = u32bitSwap(u32v.u);
  u16v.u = u16bitSwap(u16v.u);

  for (int i=0; i<8; i++)
    fprintf(stderr, "%02x", u64v.c[i]);
  fprintf(stderr, "\n");

  for (int i=0; i<4; i++)
    fprintf(stderr, "%02x", u32v.c[i]);
  fprintf(stderr, "\n");

  for (int i=0; i<2; i++)
    fprintf(stderr, "%02x", u16v.c[i]);
  fprintf(stderr, "\n");
}

#include "mt19937ar.h"

//  The MD5 checksum of the correct output is
//    cb33e6acc162cbe20f7fcac26adddd02
//  and it is 22465 bytes long.
//
//  but we cannot use md5, as it's in libbri, and
//  so is this...

int main(void) {
    int i;
    uint32 init[4] = {0x123, 0x234, 0x345, 0x456};
    uint32 length  = 4;
    mt_s *ctx = mtInitArray(init, length);

    printf("1000 outputs of genrand_int32()\n");

    for (i=0; i<1000; i++) {
      printf(uint32FMTW(10)" ", mtRandom32(ctx));
      if (i%5==4) printf("\n");
    }

    printf("\n1000 outputs of genrand_real2()\n");

    for (i=0; i<1000; i++) {
      printf("%10.8f ", mtRandomRealOpen(ctx));
      if (i%5==4) printf("\n");
    }



    for (i=0; i<999; i++) {
      printf(uint64HEX" ", mtRandom64(ctx));
      if (i%3==2) printf("\n");
    }

    return 0;
}

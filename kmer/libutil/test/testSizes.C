#include <stdio.h>
#include "../libbritypes.H"

void
main(void) {
  printf("sizeof(u64bit): %d\n", sizeof(u64bit));
  printf("sizeof(u32bit): %d\n", sizeof(u32bit));
  printf("sizeof(u16bit): %d\n", sizeof(u16bit));


  printf("u64bitZERO: 0x%016llx\n", u64bitZERO);
  printf("u64bitONE:  0x%016llx\n", u64bitONE);
  printf("u64bitMAX:  0x%016llx\n", u64bitMAX);
}

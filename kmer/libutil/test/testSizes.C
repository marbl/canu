#include <stdio.h>
#include "../libbritypes.H"

void
main(void) {

  if ((sizeof(u32bit) != 4) ||
      (sizeof(u64bit) != 8)) {
    fprintf(stderr, "ERROR: data size definitions are incorrect.\n");
    fprintf(stderr, "       u32bit has %d bytes (should be 4)!\n", sizeof(u32bit));
    fprintf(stderr, "       u64bit has %d bytes (should be 8)!\n", sizeof(u64bit));
    return(1);
  }

  return(0);
}



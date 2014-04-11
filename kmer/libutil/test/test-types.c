#include <stdio.h>

#include "util.h"

int
main(void) {
  uint32   errors = 0;
  uint32   u3 = -1;
  int32   s3 = -1;
  uint64   u6 = -1;
  int64   s6 = -1;

  if (sizeof(uint32) != 4)
    fprintf(stderr, "uint32 has %d bytes (should be 4)!\n", (int)sizeof(uint32)), errors++;

  if (sizeof(uint64) != 8)
    fprintf(stderr, "uint64 has %d bytes (should be 8)!\n", (int)sizeof(uint64)), errors++;

  if (u3 < 0)
    fprintf(stderr, "uint32 is signed (should be unsigned)!\n"), errors++;

  if (s3 > 0)
    fprintf(stderr, "int32 is unsigned (should be signed)!\n"), errors++;

  if (u6 < 0)
    fprintf(stderr, "uint64 is signed (should be unsigned)!\n"), errors++;

  if (s6 > 0)
    fprintf(stderr, "int64 is unsigned (should be signed)!\n"), errors++;

  return(errors);
}



#include "memory.h"
#include "libbritypes.h"

#include <errno.h>
#include <string.h>

void *
memdup(const void *orig, size_t size) {
  void *rslt = NULL;

  if ((orig != NULL) && (size > 0)) {
    errno = 0;
    rslt = malloc(size);
    if (errno) {
      //  Some ugliness to print out a size_t, which seems to bounce
      //  around from four bytes to eight bytes
      //
      if (sizeof(size_t) == 4)
        fprintf(stderr, "memdup()-- can't allocate "s64bitFMT" bytes.\n%s\n", size, strerror(errno));
      else
        fprintf(stderr, "memdup()-- can't allocate "s32bitFMT" bytes.\n%s\n", size, strerror(errno));
      exit(1);
    }
    memcpy(rslt, orig, size);
  }
  return(rslt);
}


#include <errno.h>
#include <string.h>

#include "util.h"

void *
memdup(const void *orig, size_t size) {
  void *rslt = NULL;

  if ((orig != NULL) && (size > 0)) {
    errno = 0;
    rslt = malloc(size);
    if (errno) {
      //  Some ugliness to print out a size_t.  This might be useless,
      //  as it might be determined by TRUE64BIT.
      //
      if (sizeof(size_t) == 8)
        fprintf(stderr, "memdup()-- can't allocate "s64bitFMT" bytes.\n%s\n", (s64bit)size, strerror(errno));
      else
        fprintf(stderr, "memdup()-- can't allocate "u32bitFMT" bytes.\n%s\n", (u32bit)size, strerror(errno));
      exit(1);
    }
    memcpy(rslt, orig, size);
  }
  return(rslt);
}


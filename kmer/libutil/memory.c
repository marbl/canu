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
#ifdef TRUE64BIT
      fprintf(stderr, "memdup()-- can't allocate %lu bytes.\n%s\n", size, strerror(errno));
#else
      fprintf(stderr, "memdup()-- can't allocate %d bytes.\n%s\n", size, strerror(errno));
#endif
      exit(1);
    }
    memcpy(rslt, orig, size);
  }
  return(rslt);
}


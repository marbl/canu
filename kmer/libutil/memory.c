#include "memory.h"

#include <errno.h>
#include <string.h>

void *
memdup(const void *orig, size_t size) {
  void *rslt = NULL;

  if ((orig != NULL) && (size > 0)) {
    errno = 0;
    rslt = malloc(size);
    if (errno) {
      fprintf(stderr, "Out of memory in memdup.\n%s\n", strerror(errno));
      exit(1);
    }
    memcpy(rslt, orig, size);
  }
  return(rslt);
}


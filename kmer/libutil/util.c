#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <sys/utsname.h>
#include <sys/resource.h>
#include <sys/stat.h>

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>


#include "util.h"


double
getTime(void) {
  struct timeval  tp;
  gettimeofday(&tp, NULL);
  return(tp.tv_sec + (double)tp.tv_usec / 1000000.0);
}


uint64
getProcessSizeCurrent(void) {
  struct rusage  ru;
  uint64         sz = 0;

  errno = 0;
  if (getrusage(RUSAGE_SELF, &ru) == -1) {
    fprintf(stderr, "getProcessSizeCurrent()-- getrusage(RUSAGE_SELF, ...) failed: %s\n",
            strerror(errno));
  } else {
    sz  = ru.ru_maxrss;
    sz *= 1024;
  }

  return(sz);
}


uint64
getProcessSizeLimit(void) {
  struct rlimit rlp;
  uint64        sz = ~uint64ZERO;

  errno = 0;
  if (getrlimit(RLIMIT_DATA, &rlp) == -1) {
    fprintf(stderr, "getProcessSizeLimit()-- getrlimit(RLIMIT_DATA, ...) failed: %s\n",
            strerror(errno));
  } else {
    sz = rlp.rlim_cur;
  }

  return(sz);
}




void *
memdup(const void *orig, size_t size) {
  void *rslt = NULL;

  if ((orig != NULL) && (size > 0)) {
    errno = 0;
    rslt = malloc(size);
    if (errno) {
      //  Some ugliness to print out a size_t.  This might be useless,
      //  as it might be determined by TRUEINT64.
      //
      if (sizeof(size_t) == 8)
        fprintf(stderr, "memdup()-- can't allocate "int64FMT" bytes.\n%s\n", (int64)size, strerror(errno));
      else
        fprintf(stderr, "memdup()-- can't allocate "uint32FMT" bytes.\n%s\n", (uint32)size, strerror(errno));
      exit(1);
    }
    memcpy(rslt, orig, size);
  }
  return(rslt);
}


int
fileExists(const char *path) {
  struct stat  s;

  return(stat(path, &s) == 0);
}


off_t
sizeOfFile(const char *path) {
  struct stat s;

  errno = 0;
  if (stat(path, &s) != 0)
    fprintf(stderr, "Couldn't stat() '%s'\n%s\n", path, strerror(errno)), exit(1);

  return(s.st_size);
}


uint64
timeOfFile(const char *path) {
  struct stat s;

  errno = 0;
  if (stat(path, &s) != 0)
    fprintf(stderr, "Couldn't stat() '%s'\n%s\n", path, strerror(errno)), exit(1);

  return(s.st_mtime);
}

#include "AS_global.H"

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>





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


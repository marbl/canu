#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>     //  gettimeofday()
#include <sys/utsname.h>  //  uname()
#include <sys/resource.h> //  getrusage()

#include "util.h"

double
getTime(void) {
  struct timeval  tp;
  gettimeofday(&tp, NULL);
  return(tp.tv_sec + (double)tp.tv_usec / 1000000.0);
}

void
write_rusage(FILE *F) {
  struct utsname un;
  struct rusage  ru;

  fprintf(F, "uname()---------------------------------\n");
  errno = 0;
  if (uname(&un) == -1) {
    fprintf(F, "uname() call failed: %s\n", strerror(errno));
  } else {
    fprintf(F, "sysname:     %s\n", un.sysname);
    fprintf(F, "nodename:    %s\n", un.nodename);
    fprintf(F, "release:     %s\n", un.release);
    fprintf(F, "version:     %s\n", un.version);
    fprintf(F, "machine:     %s\n", un.machine);
  }

  fprintf(F, "\n");

  fprintf(F, "rusage()--------------------------------\n");
  errno = 0;
  if (getrusage(RUSAGE_SELF, &ru) == -1) {
    fprintf(F, "getrusage() call failed: %s\n", strerror(errno));
  } else {
    fprintf(F, "userTime:    %f\n", ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000.0);
    fprintf(F, "systemTime:  %f\n", ru.ru_stime.tv_sec + (double)ru.ru_stime.tv_usec / 1000000.0);
    fprintf(F, "maxrss:      %ld\n", ru.ru_maxrss);
    fprintf(F, "ixrss:       %ld\n", ru.ru_ixrss);
    fprintf(F, "idrss:       %ld\n", ru.ru_idrss);
    fprintf(F, "isrss:       %ld\n", ru.ru_isrss);
    fprintf(F, "minflt:      %ld\n", ru.ru_minflt);
    fprintf(F, "majflt:      %ld\n", ru.ru_majflt);
    fprintf(F, "nswap:       %ld\n", ru.ru_nswap);
    fprintf(F, "inblock:     %ld\n", ru.ru_inblock);
    fprintf(F, "oublock:     %ld\n", ru.ru_oublock);
    fprintf(F, "msgsnd:      %ld\n", ru.ru_msgsnd);
    fprintf(F, "msgrcv:      %ld\n", ru.ru_msgrcv);
    fprintf(F, "nsignals:    %ld\n", ru.ru_nsignals);
    fprintf(F, "nvcsw:       %ld\n", ru.ru_nvcsw);
    fprintf(F, "nivcsw:      %ld\n", ru.ru_nivcsw);
  }

  fprintf(F, "\n");
}


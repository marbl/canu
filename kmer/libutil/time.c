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

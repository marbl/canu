#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>     //  gettimeofday()
#include <sys/utsname.h>  //  uname()
#include <sys/resource.h> //  getrusage()

#include "bri++.H"

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


char*
speedCounter::_spinr[4] = { "[|]", "[/]", "[-]", "[\\]" };

char*
speedCounter::_liner[19] = { "[-         ]",
                             "[--        ]",
                             "[ --       ]",
                             "[  --      ]",
                             "[   --     ]",
                             "[    --    ]",
                             "[     --   ]",
                             "[      --  ]",
                             "[       -- ]",
                             "[        --]",
                             "[         -]",
                             "[        --]",
                             "[       -- ]",
                             "[      --  ]",
                             "[     --   ]",
                             "[    --    ]",
                             "[   --     ]",
                             "[  --      ]",
                             "[ --       ]" };


speedCounter::speedCounter(char const   *fmt,
                           double        unit,
                           u64bit        freq,
                           bool          enabled) {
  _count     = 0;
  _draws     = 0;
  _unit      = unit;
  _freq      = freq;
  _startTime = getTime();
  _fmt       = fmt;
  _spin      = false;
  _line      = false;
  _enabled   = enabled;

  //  We use _draws instead of shifting _count just because it's
  //  simpler, and both methods need another variable anyway.

  //  Set all the bits below the hightest set in _freq --
  //  this allows us to do a super-fast test in tick().
  //
  _freq |= _freq >> 1;
  _freq |= _freq >> 2;
  _freq |= _freq >> 4;
  _freq |= _freq >> 8;
  _freq |= _freq >> 16;
  _freq |= _freq >> 32;
}

speedCounter::~speedCounter() {
  finish();
}

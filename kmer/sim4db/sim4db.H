#ifndef SIM4DB_H
#define SIM4DB_H

#if defined(__bsd) || defined(__alpha)
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199506L
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <signal.h>
#include <math.h>

#include <pthread.h>
#include <semaphore.h>

#include "sim4.H"
#include "bio++.H"


//  Run options set from the command line.
//
class configuration {
public:
  configuration() {
    startTime        = getTime();

    cdnaFileName     = 0L;
    scriptFileName   = 0L;
    databaseFileName = 0L;
    outputFileName   = 0L;
    yesnoFileName    = 0L;
    statsFileName    = 0L;
    touchFileName    = 0L;

    pairwise         = false;

    beVerbose        = false;
    beYesNo          = false;

    numThreads       = 3;
    loaderCacheSize  = 1024;
  };
  ~configuration() {
    if (statsFileName) {
      FILE  *statsFile = fopen(statsFileName, "w");
      if (statsFile) {
        write_rusage(statsFile);
        fprintf(statsFile, "clockTime:      %f\n", getTime() - startTime);
        fclose(statsFile);
      }
    }

  };

  void parseCommandLine(int argc, char **argv);


  double            startTime;

  char             *cdnaFileName;
  char             *scriptFileName;
  char             *databaseFileName;
  char             *outputFileName;
  char             *yesnoFileName;
  char             *statsFileName;
  char             *touchFileName;

  bool              pairwise;

  bool              beVerbose;
  bool              beYesNo;

  u32bit            numThreads;
  u32bit            loaderCacheSize;

  sim4parameters    sim4params;
};



#endif  //  SIM4DB_H

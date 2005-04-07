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

//  At one time, this was needed for pthread.h or semaphore.h
typedef unsigned short ushort;

#include <pthread.h>
#include <semaphore.h>

#include "sim4.H"
#include "bio++.H"


//  Run options set from the command line.
//
extern char             *cdnaFileName;
extern char             *scriptFileName;
extern char             *databaseFileName;
extern char             *outputFileName;
extern char             *statsFileName;
extern char             *touchFileName;

extern bool              beVerbose;
extern bool              beYesNo;

extern sim4parameters    sim4params;

extern volatile u32bit            numberOfScripts;

extern pthread_mutex_t   inputTailMutex;
extern volatile u32bit            inputTail;

extern sim4command          **input;
extern char                 **script;
extern sim4polishList       **output;
extern FastASequenceInCore  **gendelete;

extern bool              allDone;

extern FastAWrapper     *GENs;
extern FastACache       *ESTs;



void    parseCommandLine(int argc, char **argv);
u32bit  numberOfScriptLines(char *filename);
char   *getNextScript(u32bit     &ESTseq,
                      u32bit     &GENseq, u32bit &GENlo, u32bit &GENhi,
                      bool       &doForward,
                      bool       &doReverse,
                      readBuffer *scriptFile);



#endif  //  SIM4DB_H

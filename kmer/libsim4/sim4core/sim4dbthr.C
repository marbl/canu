#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199506L
#endif

typedef unsigned short ushort;

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <math.h>

#ifdef WITH_THREADS
#include <pthread.h>
#include <semaphore.h>
#define MAX_THREADS               64
#define DYNAMIC_THREAD_BALANCING  1
#endif

#include "sim4.H"
#include "fasta.H"


//  Define this if you are trying to debug with a memory allocation
//  checker.  It will greatly reduce the number of calls to malloc
//  during initialization.  (And it will tell you how many sequences
//  it has read.)
//
//#define MEMORY_DEBUG_HACKS

#define VERBOSE_OUTPUT 1



//  Some malloc options
//
int __delayed_free  = 0;
int __fast_free_max = 100;
int __first_fit     = 2;
unsigned long __noshrink = 1;

//  Global thread data
//
#ifdef WITH_THREADS
pthread_mutex_t   mutex;
bool              threadSuicide[MAX_THREADS];
pthread_t         tid;
pthread_attr_t    attr;
#endif

double            startTime;

unsigned int      ESTlen      = 0;
unsigned int      ESTmax      = 0;
unsigned int     *ESTlengths  = 0L;
unsigned char   **ESTheaders  = 0L;
unsigned char   **ESTsequence = 0L;

unsigned int      DBlen      = 0;
unsigned int      DBmax      = 0;
unsigned int     *DBlengths  = 0L;
unsigned char   **DBheaders  = 0L;
unsigned char   **DBsequence = 0L;

unsigned int      scriptLinesMax = 0;
unsigned int      scriptLinesNum = 0;
unsigned int      scriptLinesPos = 0;
char            **scriptLines    = 0L;

#ifdef WITH_THREADS
char            **outputLines    = 0L;
#endif



double
getTime(void) {
  struct timeval  tp;
  gettimeofday(&tp, NULL);
  return(tp.tv_sec + (double)tp.tv_usec / 1000000.0);
}




void
readScriptLines(char *infile) {
  FILE   *scriptFile = fopen(infile, "r");

  if (scriptFile == 0L) {
    fprintf(stderr, "Couldn't open the script file.\n");
    exit(1);
  }

  //  If a line is more than a megabyte long, we deserve to crash!
  //
  char   *scriptLine = new char [1024 * 1024];

  fprintf(stderr, "Reading script lines.\n");

  scriptLinesMax = 1024 * 512;
  scriptLinesNum = 0;
  scriptLinesPos = 0;
  scriptLines    = new char * [scriptLinesMax];

#ifdef MEMORY_DEBUG_HACKS
  //
  //  N.B.!  This will CRASH on EXIT and if it runs out of space!
  //
  int    scriptLinesSpaceMax = 16 * 1024 * 1024;
  char  *scriptLinesSpace    = new char [scriptLinesSpaceMax];
#endif

  while (!feof(scriptFile)) {
    errno = 0;
    fgets(scriptLine, 1024*1024, scriptFile);
    if (errno)
      perror("fgets");

    if (!feof(scriptFile)) {
      if (scriptLinesNum >= scriptLinesMax) {
        scriptLinesMax *= 2;
        char **tmpLines = new char * [scriptLinesMax];

        for (unsigned int i=0; i<scriptLinesNum; i++)
          tmpLines[i] = scriptLines[i];

        delete scriptLines;
        scriptLines = tmpLines;
      }

#ifdef MEMORY_DEBUG_HACKS
      scriptLines[scriptLinesNum] = scriptLinesSpace;
      scriptLinesSpace += strlen(scriptLine) + 2;
#else
      scriptLines[scriptLinesNum] = new char [strlen(scriptLine) + 1];
#endif

      strcpy(scriptLines[scriptLinesNum], scriptLine);

      scriptLinesNum++;

#ifdef MEMORY_DEBUG_HACKS
      if ((scriptLinesNum & 0xff) == 0)
        fprintf(stderr, "Read %d script lines.\r", scriptLinesNum);
#endif


    }
  }

  delete scriptLine;
  fclose(scriptFile);
}




void
readESTs(char *infile) {
  FastA         *ESTfile = new FastA(infile, 1);

  fprintf(stderr, "Reading %d ESTs.\n", ESTfile->numberOfSequences());
  fflush(stderr);

  ESTlen      = 0;
  ESTmax      = ESTfile->numberOfSequences();
  ESTlengths  = new unsigned int [ESTmax];
  ESTheaders  = new unsigned char * [ESTmax];
  ESTsequence = new unsigned char * [ESTmax];

#ifdef MEMORY_DEBUG_HACKS
  unsigned char *theBuffer    = new unsigned char [256 * 1024 * 1024];
  unsigned int   theBufferPos = 0;
#endif

  for (ESTfile->first(); ESTfile->eof() == false; ESTlen++, ESTfile->next()) {
    ESTlengths[ESTlen]  = ESTfile->sequenceLength();

#ifdef MEMORY_DEBUG_HACKS
    ESTheaders[ESTlen]  = theBuffer + theBufferPos; theBufferPos += ESTfile->headerLength()   + 4;
    ESTsequence[ESTlen] = theBuffer + theBufferPos; theBufferPos += ESTfile->sequenceLength() + 4; 
    if ((ESTlen & 0xff) == 0)
      fprintf(stderr, "Read %d ESTs.\r", ESTlen);
#else
    ESTheaders[ESTlen]  = new unsigned char [ESTfile->headerLength()   + 1];
    ESTsequence[ESTlen] = new unsigned char [ESTfile->sequenceLength() + 1];
#endif

    strcpy((char *)ESTheaders[ESTlen],  (const char *)ESTfile->header());
    strcpy((char *)ESTsequence[ESTlen], (const char *)ESTfile->sequence());

    for (unsigned int i=0; i<ESTlengths[ESTlen]; i++)
      ESTsequence[ESTlen][i] = toupper((char)ESTsequence[ESTlen][i]);
  }

  delete ESTfile;
}




void
readDB(char *infile) {
  FastA         *DBfile = new FastA(infile, 1);

  fprintf(stderr, "Reading %d DB.\n", DBfile->numberOfSequences());
  fflush(stderr);

  DBlen      = 0;
  DBmax      = DBfile->numberOfSequences();
  DBlengths  = new unsigned int [DBmax];
  DBheaders  = new unsigned char * [DBmax];
  DBsequence = new unsigned char * [DBmax];

#ifdef MEMORY_DEBUG_HACKS
  unsigned char *theBuffer    = new unsigned char [256 * 1024 * 1024];
  unsigned int   theBufferPos = 0;
#endif

  for (DBfile->first(); DBfile->eof() == false; DBlen++, DBfile->next()) {
    DBlengths[DBlen]  = DBfile->sequenceLength();

#ifdef MEMORY_DEBUG_HACKS
    DBheaders[DBlen]  = theBuffer + theBufferPos; theBufferPos += DBfile->headerLength()   + 4;
    DBsequence[DBlen] = theBuffer + theBufferPos; theBufferPos += DBfile->sequenceLength() + 4;
    if ((DBlen & 0xff) == 0)
      fprintf(stderr, "Read %d sequences.\r", DBlen);
#else
    DBheaders[DBlen]  = new unsigned char [DBfile->headerLength()   + 1];
    DBsequence[DBlen] = new unsigned char [DBfile->sequenceLength() + 1];
#endif

    strcpy((char *)DBheaders[DBlen],  (const char *)DBfile->header());
    strcpy((char *)DBsequence[DBlen], (const char *)DBfile->sequence());

    for (unsigned int i=0; i<DBlengths[DBlen]; i++)
      DBsequence[DBlen][i] = toupper((char)DBsequence[DBlen][i]);
  }

  delete DBfile;
}





#ifdef WITH_THREADS

unsigned int    finishedTesting = 0;

void *
sim4threadTimeTest(void *U) {
  int     request;

  while (scriptLinesPos < (unsigned int)U) {
    pthread_mutex_lock(&mutex);

    if (scriptLinesPos < scriptLinesNum) {
      request = scriptLinesPos;
      scriptLinesPos++;

      pthread_mutex_unlock(&mutex);

      Sim4            *S4 = new Sim4();
      sim4parameters  *P4 = new sim4parameters(scriptLines[request],
                                               ESTlengths, ESTheaders, ESTsequence,
                                               DBlengths,  DBheaders,  DBsequence);
      char            *Jk = S4->run(P4);

      delete Jk;
      delete P4;
      delete S4;
    } else {
      pthread_mutex_unlock(&mutex);
    }
  }

  pthread_mutex_lock(&mutex);
  finishedTesting++;
  pthread_mutex_unlock(&mutex);

  pthread_exit(NULL);
  return(0);
}


unsigned int
balanceThreads(void) {
  fprintf(stderr, "Determining optimal number of threads.\n");

  double         oldSpeed = 0;
  double         newSpeed = 0;
  double         startTime;
  bool           doneBalancing = false;
  unsigned int   i;
    
  for (i=0; !doneBalancing; i++) {
    fprintf(stderr, "\n     Testing:    %d threads.\n", i+1);

    scriptLinesPos   = 0;
    finishedTesting  = 0;
    startTime        = getTime();

    for (unsigned int j=0; j<=i; j++)
      pthread_create(&tid, &attr, sim4threadTimeTest, (void *)(100));

    //  Wait for a little bit.  This used to be sleep(30), but
    //  when the database sequences got long we wouldn't run
    //  enough to generate a useful number -- consider, for example,
    //  five seconds per run.
    //
    //  We now wait for some number of computes to finish.
    //
    usleep(100000);

    while (finishedTesting <= i) {
      usleep(100000);
      newSpeed = scriptLinesPos / (getTime() - startTime);
      fprintf(stderr, "     Scheduling: %d threads at %8.5f per second (%d threads done, %d lines completed).\r",
              i+1, newSpeed, finishedTesting, scriptLinesPos);
    }

    newSpeed = scriptLinesPos / (getTime() - startTime);

    if (oldSpeed > newSpeed) {
      doneBalancing = true;
    } else {
      oldSpeed = newSpeed;
    }
  }
  
  fprintf(stderr, "\n          Using: %d threads at %f per second.\n",
          i-1, oldSpeed);

  return(i-1);
}

void *
sim4thread(void *U) {
  int     request;
  int     workDone = 0;

  //  We don't want to die (yet)!
  //
  threadSuicide[(int)U] = false;

  while (threadSuicide[(int)U] == false) {

    pthread_mutex_lock(&mutex);

    if (scriptLinesPos < scriptLinesNum) {
      request = scriptLinesPos;
      scriptLinesPos++;

      pthread_mutex_unlock(&mutex);

#if VERBOSE_OUTPUT
      if ((workDone & 0xff) == 0xff) {
        fprintf(stderr, "%2d at %d -- I've done %6d, total %6.3f per second)\r",
                (int)U,
                request,
                workDone,
                request / (getTime() - startTime));
        fflush(stderr);
      }
      workDone++;
#endif

      Sim4            *S4  = new Sim4();
      sim4parameters  *P4  = new sim4parameters(scriptLines[request],
                                                ESTlengths, ESTheaders, ESTsequence,
                                                DBlengths,  DBheaders,  DBsequence);
      outputLines[request] = S4->run(P4);

      delete P4;
      delete S4;
    } else {
      pthread_mutex_unlock(&mutex);
      threadSuicide[(int)U] = true;
    }
  }

  pthread_exit(NULL);
  return(0);
}

#endif



int
main(int argc, char **argv) {
  unsigned int  numThreads = 0;
  char         *cdnaFileName = 0L;
  char         *scriptFileName = 0L;
  char         *databaseFileName = 0L;
  char         *outputFileName = 0L;;

  double        mainStartTime = getTime();

  int arg = 1;
  while (arg < argc) {
    switch (argv[arg][1]) {
    case 'h':
      fprintf(stderr, "usage: %s [-h] [-t numThreads] [-e estFile] [-d databaseFile] [-s scriptFile] [-o outputFileName]\n", argv[0]);
      fprintf(stderr, "       -h  this\n");
      fprintf(stderr, "       -t  use numThreads threads (default = adaptive)\n");
      fprintf(stderr, "       -e  ESTs are here\n");
      fprintf(stderr, "       -d  sequence database\n");
      fprintf(stderr, "       -s  use this script file\n");
      fprintf(stderr, "       -o  write output to this file\n");
      exit(0);
      break;
    case 't':
      arg++;
      numThreads = (unsigned int)atoi(argv[arg]);
      break;
    case 'e':
      arg++;
      cdnaFileName = argv[arg];
      break;
    case 'd':
      arg++;
      databaseFileName = argv[arg];
      break;
    case 's':
      arg++;
      scriptFileName = argv[arg];
      break;
    case 'o':
      arg++;
      outputFileName = argv[arg];
      break;
    default:
      fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
      break;
    }
    arg++;
  }

  if (cdnaFileName == 0L) {
    fprintf(stderr, "No EST file?\n");
    exit(1);
  }

  if (scriptFileName == 0L) {
    fprintf(stderr, "No script file?\n");
    exit(1);
  }

  if (outputFileName == 0L) {
    fprintf(stderr, "No output file?\n");
    exit(1);
  }


  //  Read the files -- scripts, ESTs and the database
  //
  readScriptLines(scriptFileName);

  if (databaseFileName)
    readDB(databaseFileName);

  readESTs(cdnaFileName);

#ifdef WITH_THREADS

  //  Allocate space for the outputLines
  //
  outputLines = new char * [scriptLinesNum];
  for (unsigned int i=0; i<scriptLinesNum; i++)
    outputLines[i] = 0L;

  //  Initialize pthreads
  //
  pthread_mutex_init(&mutex, NULL);

  //  Run the threads
  //
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
  pthread_attr_setschedpolicy(&attr, SCHED_OTHER);

  if (numThreads == 0)
    numThreads = balanceThreads();

  fprintf(stderr, "All Done\n");

  scriptLinesPos = 0;
  startTime      = getTime();

  fprintf(stderr, "Starting %d threads.\n", numThreads);
  for (unsigned int i=0; i<numThreads; i++)
    pthread_create(&tid, &attr, sim4thread, (void *)i);


  pthread_attr_destroy(&attr);



  //  Write the output as it shows up.  This doesn't race with anything,
  //  since the value we are waiting on is an atomic copy.
  //
  FILE *fOutput = fopen(outputFileName, "w");

  for (unsigned int i=0; i<scriptLinesNum; i++) {
    while (outputLines[i] == 0L)
      sleep(1);
    fprintf(fOutput, "%s", outputLines[i]);
    fflush(fOutput);

    delete [] outputLines[i];
  }

#else

  //  Write the output as it shows up.  This doesn't race with anything,
  //  since the value we are waiting on is an atomic copy.
  //
  FILE *fOutput = fopen(outputFileName, "w");
  char *o;

  startTime = getTime();

  for (unsigned int i=0; i<scriptLinesNum; i++) {
#if VERBOSE_OUTPUT
    if ((i & 0xff) == 0xff) {
      fprintf(stderr, " 0 at %d -- I've done %6d, total %6.3f per second)\r",
              i,
              i,
              i / (getTime() - startTime));
      fflush(stderr);
    }
#endif

    Sim4            *S4 = new Sim4();
    sim4parameters  *P4 = new sim4parameters(scriptLines[i],
                                             ESTlengths, ESTheaders, ESTsequence,
                                             DBlengths,  DBheaders,  DBsequence);
    o = S4->run(P4);

    delete P4;
    delete S4;

    fprintf(fOutput, "%s", o);
    fflush(fOutput);

    delete [] o;
  }

#endif

  fprintf(fOutput, "TotalRunTime: %f\n", getTime() - mainStartTime);

  fclose(fOutput);

  return(0);
}

#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <new>
#include "searchGENOME.H"
#include "bio++.H"
#include "existDB.H"


//  Shared data
//
configuration          config;
FastAWrapper          *qsFASTA          = 0L;
existDB               *maskDB           = 0L;
existDB               *onlyDB           = 0L;
positionDB            *positions        = 0L;
volatile u32bit        numberOfQueries  = 0;
u32bit                *queryMatchCounts = 0L;
char                 **output           = 0L;
u32bit                *outputLen        = 0L;
pthread_mutex_t        inputTailMutex;
FastASequenceInCore  **input            = 0L;
volatile u32bit        inputHead        = 0;
volatile u32bit        inputTail        = 0;
volatile u32bit        outputPos        = 0;
char                  *threadStats[MAX_THREADS] = { 0L };




  //  Write the stats
  //
void
dumpStats(void) {

  if (config._statsFileName) {
    FILE *F = fopen(config._statsFileName, "wb");
    if (F == 0L) {
      fprintf(stderr, "Couldn't open the stats file '%s'?\n", config._statsFileName);
      fprintf(stderr, "Stats going to stderr.\n");
      config._statsFileName = 0L;
      F = stderr;
    }

    config.display(F);

    write_rusage(F);
      
    fprintf(F, "wallClockTimes--------------------------\n");
    fprintf(F, "init:     %9.5f\n", config._initTime   - config._startTime);
    fprintf(F, "build:    %9.5f\n", config._buildTime  - config._initTime);
    fprintf(F, "search:   %9.5f\n", config._searchTime - config._buildTime);
    fprintf(F, "total:    %9.5f\n", config._totalTime  - config._startTime);

    fprintf(F, "searchThreadInfo------------------------\n");
    for (u64bit i=0; i<config._numSearchThreads; i++)
      if (threadStats[i])
        fprintf(F, threadStats[i]);

    if (config._statsFileName)
      fclose(F);
  }
}




int
main(int argc, char **argv) {

#ifdef _AIX
  //  By default, AIX Visual Age C++ new() returns 0L; this turns on
  //  exceptions.
  //
  std::__set_new_throws_exception(true);
#endif

  //  Read the configuration from the command line
  //
  if (argc < 2) {
    config.usage(argv[0]);
    exit(1);
  }
  config.read(argc, argv);
  config.display();

  config._startTime = getTime();

  //  Open and init the query sequence
  //
  if (config._beVerbose)
    fprintf(stderr, "Opening the cDNA sequences.\n");

  qsFASTA = new FastAWrapper(config._qsFileName);
  qsFASTA->openIndex();

  numberOfQueries  = qsFASTA->getNumberOfSequences();
  output           = new char * [numberOfQueries];
  outputLen        = new u32bit [numberOfQueries];
  queryMatchCounts = new u32bit [numberOfQueries];

  input     = new FastASequenceInCore * [numberOfQueries];
  inputHead = 0;
  inputTail = 0;

  for (u32bit i=numberOfQueries; i--; ) {
    output[i]    = 0L;
    outputLen[i] = 0;

    queryMatchCounts[i] = 0;

    input[i]    = 0L;
  }

  config._initTime = getTime();





  //
  //  Build the positions
  //



  //  Complete the configuration
  //
  config._useList.setSource(config._dbFileName);
  config._useList.setSeparatorLength(1);
  config._useList.finish();




  //  Read in the positionDB if it's already built, or build a new one.
  //
  if ((config._tableFileName) && (fileExists(config._tableFileName))) {
    if (config._tableBuildOnly) {
      fprintf(stderr, "All done.  Table '%s' already built.\n", config._tableFileName);
      exit(1);
    } else {
      fprintf(stderr, "Loading positionDB state from '%s'\n", config._tableFileName);
      positions = new positionDB(config._tableFileName, true);
    }
  } else {
    merStream *MS = new merStream(config._merSize, &config._useList);

    //  Figure out a nice size of the hash.
    //
    //  XXX:  This probably should be tuned.
    //
    u32bit tblSize = 25;
    if (config._useList.lengthOfSequences() < 64 * 1024 * 1024) tblSize = 24;
    if (config._useList.lengthOfSequences() < 16 * 1024 * 1024) tblSize = 23;
    if (config._useList.lengthOfSequences() <  4 * 1024 * 1024) tblSize = 22;
    if (config._useList.lengthOfSequences() <  2 * 1024 * 1024) tblSize = 21;
    if (config._useList.lengthOfSequences() <  1 * 1024 * 1024) tblSize = 20;

    positions = new positionDB(MS, config._merSize, config._merSkip, tblSize, 0L, 0L, 0, config._beVerbose);

    delete    MS;

    if (config._tableFileName) {
      if (config._beVerbose)
        fprintf(stderr, "Dumping positions table to '%s'\n", config._tableFileName);

      positions->saveState(config._tableFileName);

      if (config._tableBuildOnly) {
        dumpStats();
        exit(0);
      }
    }
  }

  config._buildTime = getTime() - config._startTime - config._initTime;





  //  Build the masking database.
  //
  //  Previous versions build the existDB takeing the posDB as a
  //  parameter.  The existDB would then be exclude mers not in the
  //  posDB.  A neat and nice feature, but with only 45,000 to 70,000
  //  mers in the masks, hardly worth the effort.
  //
  maskDB = 0L;
  if (config._maskFileName) {
    if (config._beVerbose)
      fprintf(stderr, "Building maskDB from '%s'\n", config._maskFileName);
    maskDB = new existDB(config._maskFileName, config._merSize, 19, u32bitZERO, ~u32bitZERO);
  }

  onlyDB = 0L;
  if (config._onlyFileName) {
    if (config._beVerbose)
      fprintf(stderr, "Building onlyDB from '%s'\n", config._onlyFileName);
    onlyDB = new existDB(config._onlyFileName, config._merSize, 19, u32bitZERO, ~u32bitZERO);
  }

  config._buildTime = getTime();

  //fprintf(stderr, "Initializing threads.\n");

  //
  //  Initialize threads
  //
  pthread_attr_t   threadAttr;
  pthread_t        threadID;

  pthread_mutex_init(&inputTailMutex, NULL);
  //pthread_mutex_init(&queryMatchMutex, NULL);

  pthread_attr_init(&threadAttr);
  pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_DETACHED);
  pthread_attr_setschedpolicy(&threadAttr, SCHED_OTHER);

  //  Start the deadlock detection threads
  //
#ifdef __alpha
  fprintf(stderr, "Deadlock detection enabled!\n");
  pthread_create(&threadID, &threadAttr, deadlockDetector, 0L);
  pthread_create(&threadID, &threadAttr, deadlockChecker, 0L);
#endif

  //  Start the loader thread
  //
  pthread_create(&threadID, &threadAttr, loaderThread, 0L);

  //  Start the search threads
  //
  for (u64bit i=0; i<config._numSearchThreads; i++) {
    threadStats[i] = 0L;
    pthread_create(&threadID, &threadAttr, searchThread, (void *)i);
  }


  //  Open output file
  //
  int resultFILE = fileno(stdout);

  if (config._outputFileName) {
    errno = 0;
    resultFILE = open(config._outputFileName,
                      O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno) {
      fprintf(stderr, "Couldn't open the output file '%s'?\n%s\n", config._outputFileName, strerror(errno));
      exit(1);
    }
  }


  //  Wait for threads to produce output
  //
  outputPos = 0;

  double  zeroTime    = getTime() - 0.00000001;
  u32bit  outputMask  = 0xf;

  while (outputPos < numberOfQueries) {
    bool  justSlept = false;

    if (output[outputPos]) {
      if (config._beVerbose &&
          ((justSlept) || (outputPos & outputMask) == outputMask)) {

        double finish = 0.0;
        if (outputPos > 0)
          finish = (numberOfQueries - outputPos) / (outputPos / (getTime() - zeroTime));

        fprintf(stderr, "O:"u32bitFMTW(7)" S:"u32bitFMTW(7)" I:"u32bitFMTW(7)" T:"u32bitFMTW(7)" (%5.1f%%; %8.3f/sec) Finish in %5.2f seconds.\r",
                outputPos,
                inputTail,
                inputHead,
                numberOfQueries,
                100.0 * outputPos / numberOfQueries,
                outputPos / (getTime() - zeroTime),
                finish);
        fflush(stderr);

        double perSec    = outputPos / (getTime() - zeroTime + 0.0000001);

        if      (perSec < 32.0)
          outputMask = 0xf;
        else if (perSec < 256.0)
          outputMask = 0x7f;
        else if (perSec < 1024.0)
          outputMask = 0x1ff;
        else
          outputMask = 0x3ff;
      }

      if (outputLen[outputPos] > 0) {
        errno = 0;
        write(resultFILE, output[outputPos], sizeof(char) * outputLen[outputPos]);
        if (errno) {
          fprintf(stderr, "Couldn't write to the output file '%s'.\n%s\n",
                  config._outputFileName, strerror(errno));
          exit(1);
        }
      }

      delete [] input[outputPos];
      delete [] output[outputPos];

      input[outputPos]     = 0L;
      output[outputPos]    = 0L;
      outputLen[outputPos] = 0;

      outputPos++;
    } else {
      nanosleep(&config._writerSleep, 0L);
    }
  }

  if (config._beVerbose) {
    fprintf(stderr, "\n"u32bitFMTW(7)" sequences in %5.2f seconds, %8.3f per second.\n",
            numberOfQueries,
            getTime() - zeroTime,
            numberOfQueries / (getTime() - zeroTime));
  }

  errno = 0;
  close(resultFILE);
  if (errno)
    fprintf(stderr, "Couldn't close to the output file '%s'.\n%s\n", config._outputFileName, strerror(errno));

  config._searchTime = getTime();

  //  Clean up
  //
  delete positions;

  pthread_attr_destroy(&threadAttr);
  pthread_mutex_destroy(&inputTailMutex);
  //pthread_mutex_destroy(&queryMatchMutex);

  //  Write the query match counts
  //
  if (config._queryMatchFileName) { 
    errno = 0;
    int mcf = open(config._queryMatchFileName,
                   O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno) {
      fprintf(stderr, "ESTmapper/search-- couldn't open match counts file '%s'\n%s\n", config._queryMatchFileName, strerror(errno));
      exit(1);
    }
    char   mcfstr[256];
    size_t mcfstrend = 0;
    for (u32bit i=0; i<numberOfQueries; i++) {
      sprintf(mcfstr, u32bitFMT"\n", queryMatchCounts[i]);
      for (mcfstrend=0; mcfstr[mcfstrend]; mcfstrend++)
        ; 
      errno = 0;
      write(mcf, mcfstr, mcfstrend);
      if (errno) {
        fprintf(stderr, "Couldn't write to the match counts file '%s'.\n%s\n",
                config._queryMatchFileName, strerror(errno));
        exit(1);
      }
    }
    errno = 0;
    close(mcf);
    if (errno) {
      fprintf(stderr, "Couldn't close the match counts file '%s'.\n%s\n",
              config._queryMatchFileName, strerror(errno));
      exit(1);
    }
  }


  config._totalTime = getTime();

  dumpStats();

  return(0);
}


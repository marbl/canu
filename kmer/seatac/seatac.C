#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <new>
#include "seatac.H"


//  Shared data
//
configuration          config;
seqFile               *qsFASTA          = 0L;
positionDB            *positions        = 0L;
volatile u32bit        numberOfQueries  = 0;
filterObj            **output           = 0L;
pthread_mutex_t        inputTailMutex;
seqInCore            **input            = 0L;
volatile u32bit        inputHead        = 0;
volatile u32bit        inputTail        = 0;
volatile u32bit        outputPos        = 0;
char                  *threadStats[MAX_THREADS] = { 0L };



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




#ifdef _AIX
static
void
aix_new_handler() {
  fprintf(stderr, "aix_new_handler()-- Memory allocation failed.\n");
  throw std::bad_alloc();
}
#endif


int
main(int argc, char **argv) {

#ifdef _AIX
  //  By default, AIX Visual Age C++ new() returns 0L; this turns on
  //  exceptions.
  //
  std::set_new_handler(aix_new_handler);
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
  qsFASTA = openSeqFile(config._qsFileName);
  qsFASTA->openIndex();

  numberOfQueries  = qsFASTA->getNumberOfSequences();
  output           = new filterObj * [numberOfQueries];
  input            = new seqInCore * [numberOfQueries];
  inputHead        = 0;
  inputTail        = 0;

  for (u32bit i=numberOfQueries; i--; ) {
    output[i]    = 0L;
    input[i]     = 0L;
  }

  config._initTime = getTime();


  //  Complete the configuration
  //
  config._useList.setFile(config._dbFileName);
  config._useList.setSeparator('.', 1);
  config._useList.finish();


  //  Create the chunk, returning a positionDB.  Threads will use both
  //  chain and postions to build hitMatrices.
  //
  if ((config._tableFileName) && (fileExists(config._tableFileName))) {
    if (config._tableBuildOnly) {
      fprintf(stderr, "All done.  Table '%s' already build.\n", config._tableFileName);
      exit(0);
    } else {
      fprintf(stderr, "Loading positionDB state from '%s'\n", config._tableFileName);
      positions = new positionDB(config._tableFileName, true);
    }
  } else {

    existDB *maskDB = 0L;
    if (config._maskFileName) {
      if (config._beVerbose)
        fprintf(stderr, "Building maskDB from '%s'\n", config._maskFileName);
      maskDB = new existDB(config._maskFileName, config._merSize, 0, ~u32bitZERO, existDBcanonical | existDBcompressHash | existDBcompressBuckets);
    }

    existDB *onlyDB = 0L;
    if (config._onlyFileName) {
      if (config._beVerbose)
        fprintf(stderr, "Building onlyDB from '%s'\n", config._onlyFileName);
      onlyDB = new existDB(config._onlyFileName, config._merSize, 0, ~u32bitZERO, existDBcanonical | existDBcompressHash | existDBcompressBuckets);
    }

    kMerBuilder KB(config._merSize);
    merStream  *MS = new merStream(&KB, &config._useList);

    positions = new positionDB(MS, config._merSize, config._merSkip, maskDB, onlyDB, 0, config._beVerbose);

    delete    MS;

    delete maskDB;
    delete onlyDB;

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

  config._buildTime = getTime();


  //
  //  Initialize threads
  //
  pthread_attr_t   threadAttr;
  pthread_t        threadID;

  pthread_mutex_init(&inputTailMutex, NULL);

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
  for (u32bit i=0; i<config._numSearchThreads; i++) {
    threadStats[i] = 0L;
    pthread_create(&threadID, &threadAttr, searchThread, (void *)i);
  }


  //  Open output file
  //
  FILE *resultFILE = stdout;

  if (config._outputFileName) {
    errno = 0;
    int rf = open(config._outputFileName,
                  O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                  S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno) {
      fprintf(stderr, "Couldn't open the output file '%s'?\n%s\n", config._outputFileName, strerror(errno));
      exit(1);
    }

    resultFILE = fdopen(rf, "w");
  }


  //  Dump our information to the output file.
  //
  config.writeATACheader(resultFILE);



  //  Initialize the statistics collection object
  //
  statObj *stats = new statObj(config._filterObj, config._filteropts);

  //  Wait for threads to produce output
  //
  outputPos = 0;

  //  The match id of each output record.
  //
  u64bit   matchID = 0;

  double  zeroTime = getTime() - 0.00000001;

  while (outputPos < numberOfQueries) {
    if (output[outputPos]) {
      if (config._beVerbose && ((outputPos & 0x1ff) == 0x1ff)) {
        fprintf(stderr,
                "O:"u32bitFMTW(7)" S:"u32bitFMTW(7)" I:"u32bitFMTW(7)" T:"u32bitFMTW(7)" (%5.1f%%; %8.3f/sec) Finish in %5.2f seconds.\r",
                outputPos,
                inputTail,
                inputHead,
                numberOfQueries,
                100.0 * outputPos / numberOfQueries,
                outputPos / (getTime() - zeroTime),
                (numberOfQueries - outputPos) / (outputPos / (getTime() - zeroTime)));
        fflush(stderr);
      }

      errno = 0;
      matchID = output[outputPos]->output(resultFILE, matchID);
      if (errno) {
        fprintf(stderr, "Couldn't write to the output file '%s'.\n%d: %s\n",
                config._outputFileName, errno, strerror(errno));
        exit(1);
      }

      //  Add this set of results to the statistics collector
      //
      stats->add(output[outputPos]);

      //stats->show(stderr);
 
      delete input[outputPos];
      delete output[outputPos];

      input[outputPos]     = 0L;
      output[outputPos]    = 0L;

      outputPos++;
    } else {
      nanosleep(&config._writerSleep, 0L);
    }
  }

  if (config._beVerbose) {
    fprintf(stderr, "\n"u32bitFMTW(7)" sequences (%5.1f%%; %8.3f/sec) %5.2f seconds.\n",
            numberOfQueries,
            100.0 * outputPos / numberOfQueries,
            outputPos / (getTime() - zeroTime),
            getTime() - zeroTime);
  }

  //  Print statistics
  //
  stats->show(resultFILE);
  delete stats;

  errno = 0;
  fclose(resultFILE);
  if (errno)
    fprintf(stderr, "Couldn't close to the output file '%s'.\n%s\n", config._outputFileName, strerror(errno));

  config._searchTime = getTime();

  //  Clean up
  //
  delete positions;

  pthread_attr_destroy(&threadAttr);
  pthread_mutex_destroy(&inputTailMutex);

  config._totalTime = getTime();

  //  Write the stats
  //
dumpStats();

  delete [] input;
  delete [] output;

  return(0);
}

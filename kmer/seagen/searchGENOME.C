#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <new>
#include "searchGENOME.H"
#include "bri++.H"
#include "existDB.H"


//  Shared data
//
configuration          config;
FastAWrapper          *qsFASTA;
existDB               *maskDB;
existDB               *onlyDB;
positionDB            *positions;
volatile u32bit        numberOfQueries;
u32bit                *queryMatchCounts;
//pthread_mutex_t        queryMatchMutex;
char                 **output;
u32bit                *outputLen;
pthread_mutex_t        inputTailMutex;
FastASequenceInCore  **input;
volatile u32bit        inputHead;
volatile u32bit        inputTail;
volatile u32bit        outputPos;
char                  *threadStats[MAX_THREADS];



void
buildChunk(void) {

  if (config._beVerbose)
    fprintf(stderr, "Opening the genomic database.\n");

  FastAWrapper *dbFASTA = new FastAWrapper(config._dbFileName);
  dbFASTA->openIndex();

  //  Complete the configuration
  //
  config.completeUseList(dbFASTA);

  if (config._beVerbose)
    fprintf(stderr, "Building chunk with "u32bitFMT" sequences.\n", config._useListLen);

  //  If we are given a pointer to a prebuilt table, just load it.
  //  Otherwise, chain the sequence together and build a fresh table,
  //  saving and exiting if we were told so.

  if ((config._tableFileName) && (fileExists(config._tableFileName))) {
    if (config._tableBuildOnly) {
      if (config._beVerbose)
        fprintf(stderr, "All done.  Table '%s' already exists.\n", config._tableFileName);
    } else {
      if (config._beVerbose)
        fprintf(stderr, "Loading positions table from '%s'\n", config._tableFileName);
      positions = new positionDB(config._tableFileName);
      if (config._beVerbose)
        fprintf(stderr, "Loading finished.\n");
    }
  } else {

    //  sum the length of the sequences, including the padding.
    //
    u32bit sLen = 0;
    for (u32bit i=0; i<config._useListLen; i++) {
      config._useList[i].size  = dbFASTA->sequenceLength(config._useList[i].seq);
      config._useList[i].start = sLen;

      sLen += config._useList[i].size + 100;
    }


    //  Allocate space for the chained sequence
    //
    char *s = new char [sLen + 1];
    char *t = s;

    //  Chain
    //
    u32bit i;
    for (i=0; i<config._useListLen; i++) {
      dbFASTA->find(config._useList[i].seq);

      //  If our sequences are squeezed (we know it's randomly
      //  accessable) use the OnDisk mechanism, otherwise, use the
      //  InCore.  InCore uses more memory, and is slower.
      //
      if (dbFASTA->isSqueezed()) {
        FastASequenceOnDisk  *B = dbFASTA->getSequenceOnDisk();

        B->getChars(t, 0, B->sequenceLength());

        t += B->sequenceLength();

        for (u32bit gn = 100; gn--; )
          *(t++) = '.';

        delete B;
      } else {
        FastASequenceInCore  *B = dbFASTA->getSequence();

        char const *g  = B->sequence();

        while (*g)
          *(t++) = *(g++);

        for (u32bit gn = 100; gn--; )
          *(t++) = '.';

        delete B;
      }
    }

    *t = 0;


    //  Build a merStream from the sequence.  In the future, we should put this
    //  thing on disk, since the search doesn't require sequence, and thus
    //  we could make a bigger table in the same footprint if we didn't have
    //  the sequence in core too.
    //
    merStream            *MS = 0L;
    merStreamFileBuilder *MB = 0L;
    merStreamFileReader  *MR = 0L;

    if (config._beVerbose)
      fprintf(stderr, "Storing the sequence in a compressed merStreamFile.\n");
    MS = new merStream(config._merSize, s, sLen);
    MB = new merStreamFileBuilder(MS, config._tableTemporaryFileName);
    MB->build();

    delete    MB;
    delete    MS;
    delete [] s;

    MR = new merStreamFileReader(config._tableTemporaryFileName);
    MS = new merStream(MR);



    //  Figure out a nice size of the hash.
    //
    //  XXX:  This probably should be tuned.
    //
    u32bit tblSize = 25;
    if (sLen < 64 * 1024 * 1024) tblSize = 24;
    if (sLen < 16 * 1024 * 1024) tblSize = 23;
    if (sLen <  4 * 1024 * 1024) tblSize = 22;
    if (sLen <  2 * 1024 * 1024) tblSize = 21;
    if (sLen <  1 * 1024 * 1024) tblSize = 20;



    if (config._beVerbose)
      fprintf(stderr, "Building the positionDB.\n");
    positions = new positionDB(MS, config._merSize, config._merSkip, tblSize, 0L, 0L, config._beVerbose);

    delete    MS;
    delete    MR;

    if (config._tableFileName) {
      if (config._beVerbose)
        fprintf(stderr, "Dumping positions table to '%s'\n", config._tableFileName);

      positions->saveState(config._tableFileName);
    }
  }

  if (config._tableBuildOnly) {
    if (config._beVerbose)
      fprintf(stderr, "Exiting from -buildonly.\n");

    config._totalTime = getTime();

    //  Write the stats
    //
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
      fprintf(F, "build:    %9.5f\n", config._totalTime  - config._initTime);
      fprintf(F, "total:    %9.5f\n", config._totalTime  - config._startTime);

      if (config._statsFileName)
        fclose(F);
    }

    exit(0);
  }

  delete    dbFASTA;
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


  //  Create the chunk, returning a positionDB.  Threads will use both
  //  chain and postions to build hitMatrices.
  //
  buildChunk();


  //  Build the masking database
  //
  maskDB = 0L;
  if (config._maskFileName) {
    if (config._beVerbose)
      fprintf(stderr, "Building maskDB from '%s'\n", config._maskFileName);
    maskDB = new existDB(config._maskFileName, config._merSize, 19, u32bitZERO, ~u32bitZERO, positions);
  }

  onlyDB = 0L;
  if (config._onlyFileName) {
    if (config._beVerbose)
      fprintf(stderr, "Building onlyDB from '%s'\n", config._onlyFileName);
    onlyDB = new existDB(config._onlyFileName, config._merSize, 19, u32bitZERO, ~u32bitZERO, positions);
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
          (outputPos > 0) &&
          ((justSlept) || (outputPos & outputMask) == outputMask)) {
        double thisTimeD = getTime() - zeroTime + 0.0000001;
        double perSec    = outputPos / thisTimeD;
        double remTime   = (numberOfQueries - outputPos) * thisTimeD / outputPos;
        
        fprintf(stderr, "O:"u32bitFMTW(7)" S:"u32bitFMTW(7)" I:%7u T:"u32bitFMTW(7)" (%5.1f%%; %8.3f/sec) Finish in %5.2f seconds.\r",
                outputPos,
                inputTail,
                inputHead,
                numberOfQueries,
                100.0 * outputPos / numberOfQueries,
                outputPos / (getTime() - zeroTime),
                (numberOfQueries - outputPos) / (outputPos / (getTime() - zeroTime)));
        fflush(stderr);

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

  //  Write the stats
  //
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
      fprintf(F, threadStats[i]);

    if (config._statsFileName)
      fclose(F);
  }

  return(0);
}


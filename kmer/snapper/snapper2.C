#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "bio++.H"
#include "snapper2.H"


//  The (private) structure for testing various filters.
//
struct filterStats {
  double L;
  double H;
  double V;
  u32bit tp;
  u32bit tn;
  u32bit fp;
  u32bit fn;
};


//  Shared data
//
configuration          config;
sim4parameters         sim4params;
FastACache            *cache;
FastAWrapper          *qsFASTA;
existDB               *maskDB;
existDB               *onlyDB;
positionDB            *positions;
volatile u32bit        numberOfQueries;

aHit                 **answer;
u32bit                *answerLen;

char                 **output;
u32bit                *outputLen;

logMsg               **logmsg;

pthread_mutex_t        inputTailMutex;
FastASequenceInCore  **input;

volatile u32bit        inputHead;
volatile u32bit        inputTail;
volatile u32bit        outputPos;

char                  *threadStats[MAX_THREADS];



void
writeValidationFile(char *name, filterStats *theFilters, u32bit numFilters) {

  FILE *F = fopen(name, "wb");
  if (F == 0L) {
    fprintf(stderr, "Couldn't open the validation output file '%s'?\n", name);
    F = stdout;
  }

  fprintf(F, "%6s %6s %6s  %6s %6s  %8s %8s %8s %8s\n",
          "L", "H", "V",
          "sens", "spec",
          "tp", "fp", "tn", "fn");

  for (u32bit f=0; f<numFilters; f++) {
    double sens = 0.0;
    double spec = 0.0;

    if (theFilters[f].tp + theFilters[f].fn > 0)
      sens = (double)theFilters[f].tp / (theFilters[f].tp + theFilters[f].fn);

    if (theFilters[f].tn + theFilters[f].fp > 0)
      spec = (double)theFilters[f].tn / (theFilters[f].tn + theFilters[f].fp);

    fprintf(F, "%6.4f %6.4f %6.4f  %6.4f %6.4f  %8lu %8lu %8lu %8lu\n",
            theFilters[f].L,
            theFilters[f].H,
            theFilters[f].V,
            sens, spec,
            theFilters[f].tp,
            theFilters[f].fp,
            theFilters[f].tn,
            theFilters[f].fn);
  }

  fclose(F);
}





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
    fprintf(F, "init:     %9.5f\n", config._initTime);
    fprintf(F, "build:    %9.5f\n", config._buildTime);
    fprintf(F, "search:   %9.5f\n", config._searchTime);
    fprintf(F, "total:    %9.5f\n", getTime() - config._startTime);

    fprintf(F, "searchThreadInfo------------------------\n");
    for (u64bit i=0; i<config._numSearchThreads; i++) {
      if (threadStats[i])
        fprintf(F, threadStats[i]);
      delete [] threadStats[i];
    }

    if (config._statsFileName)
      fclose(F);
  }
}






#ifdef _AIX

//  If we're AIX, define a new handler.  Other OS's reliably throw exceptions.
//
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
  //  exceptions (sorta -- it sets a handler that throws an
  //  exception).
  //
  std::set_new_handler(aix_new_handler);
#endif


#if 0
  fprintf(stderr, "Hello.  I'm snapper2, pid %d\n", getpid());
  fprintf(stderr, "gdb snapper2 %d\n", getpid());
  fprintf(stderr, "Attach a debugger in the next 5 seconds.\n");
  sleep(5);
  fprintf(stderr, "Here I go!\n");
#endif


  //
  //  Read the configuration from the command line
  //
  if (argc < 2) {
    config.usage(argv[0]);
    exit(1);
  }
  config.read(argc, argv);
  //config.display();


  //  Allocate some structures for doing a validation run.  This is
  //  done pretty early, just in case it needs to abort.
  //
  u32bit        numFilters = 0;
  u32bit        maxFilters = 21 * 22 / 2 * 20;
  filterStats  *theFilters = 0L;

  if (config._doValidation) {
    theFilters = new filterStats [maxFilters];

    for (u32bit h=0; h<=100; h+=5) {
      for (u32bit l=0; l<=h; l+=5) {
        for (u32bit v=5; v<=100; v+=5) {
          if (numFilters >= maxFilters) {
            fprintf(stderr, "ERROR:  Ran out of filterStats structures while configuring the filters!\n");
            exit(1);
          }

          theFilters[numFilters].L  = l / 100.0;
          theFilters[numFilters].H  = h / 100.0;
          theFilters[numFilters].V  = v / 100.0;
          theFilters[numFilters].tp = 0;
          theFilters[numFilters].tn = 0;
          theFilters[numFilters].fp = 0;
          theFilters[numFilters].fn = 0;
          numFilters++;
        }
      }
    }

    fprintf(stderr, "Created "u32bitFMT" filters (out of "u32bitFMT" available) to test/validate.\n",
            numFilters, maxFilters);
  }



  //  Complete the configuration
  //
  config._useList.setSource(config._dbFileName);
  config._useList.setSeparatorLength(1);
  config._useList.finish();


  //  Take a snapshot of the init time here.  We build next, then we finish some
  //  more initialization
  //
  config._initTime = getTime() - config._startTime;


  //  Read in the positionDB if it's already built, or build a new one.
  //
  if ((config._psFileName) && (fileExists(config._psFileName))) {
    if (config._buildOnly) {
      fprintf(stderr, "All done.  Table '%s' already built.\n", config._psFileName);
      exit(1);
    } else {
      fprintf(stderr, "Loading positionDB state from '%s'\n", config._psFileName);
      positions = new positionDB(config._psFileName, true);
    }
  } else {



    //  The masking databases
    //
    maskDB = 0L;
    if (config._maskFileName) {
      if (config._beVerbose)
        fprintf(stderr, "Building maskDB from fasta file '%s'\n", config._maskFileName);
      maskDB = new existDB(config._maskFileName, config._merSize, 19);
    }
    if (config._maskPrefix) {
      if (config._beVerbose)
        fprintf(stderr, "Building maskDB from meryl prefix '%s'\n", config._maskPrefix);
      maskDB = new existDB(config._maskPrefix, config._merSize, 19, config._maskThreshold, ~u32bitZERO);
    }

    onlyDB = 0L;
    if (config._onlyFileName) {
      if (config._beVerbose)
        fprintf(stderr, "Building onlyDB from fasta file '%s'\n", config._onlyFileName);
      onlyDB = new existDB(config._onlyFileName, config._merSize, 19);
    }
    if (config._onlyPrefix) {
      if (config._beVerbose)
        fprintf(stderr, "Building onlyDB from meryl prefix '%s'\n", config._onlyPrefix);
      onlyDB = new existDB(config._onlyPrefix, config._merSize, 19, 0, config._onlyThreshold);
    }


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

    merStream *MS = new merStream(config._merSize, &config._useList);

    positions = new positionDB(MS, config._merSize, config._merSkip, tblSize,
                               maskDB, onlyDB, config._ignoreThreshold, config._beVerbose);

    delete    MS;

    delete    maskDB;
    delete    onlyDB;

    maskDB = 0L;
    onlyDB = 0L;

    if (config._psFileName) {
      if (config._beVerbose)
        fprintf(stderr, "Dumping positions table to '%s'\n", config._psFileName);

      positions->saveState(config._psFileName);

      if (config._buildOnly) {
        dumpStats();
        exit(0);
      }

      delete positions;
      positions = new positionDB(config._psFileName, true);
    }
  }

  config._buildTime = getTime() - config._startTime - config._initTime;



  //  Open and init the genomic sequences.
  //
  if (config._beVerbose)
    fprintf(stderr, "Opening the genomic database.\n");

  cache = new FastACache(config._dbFileName, 0, true);
  //cache = new FastACache(config._dbFileName, 256, false);



  //  Open and init the query sequence.
  //
  if (config._beVerbose)
    fprintf(stderr, "Opening the cDNA sequences.\n");

  qsFASTA = new FastAWrapper(config._qsFileName);
  qsFASTA->openIndex();

  numberOfQueries  = qsFASTA->getNumberOfSequences();

  //  All we use the index for is to count the number of sequences for
  //  the output display (and, OK, sizing the queues).  Close the
  //  index to free up significant memory on large datasets.

  delete qsFASTA;
  qsFASTA = new FastAWrapper(config._qsFileName);

  input            = new FastASequenceInCore * [numberOfQueries];
  inputHead        = 0;
  inputTail        = 0;
  answerLen        = new u32bit   [numberOfQueries];
  answer           = new aHit *   [numberOfQueries];
  outputLen        = new u32bit   [numberOfQueries];
  output           = new char *   [numberOfQueries];
  logmsg           = new logMsg * [numberOfQueries];

  for (u32bit i=0; i<numberOfQueries; i++) {
    input[i]            = 0L;
    answerLen[i]        = 0;
    answer[i]           = 0L;
    outputLen[i]        = 0;
    output[i]           = 0L;
    logmsg[i]           = 0L;
  }


  //  Init all done!
  //
  config._initTime = getTime() - config._startTime - config._buildTime;


#ifdef MEMORY_DEBUG
  fprintf(stdout, "----------------------------------------\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--  Dump at start\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--\n");
  _dump_allocated_delta(fileno(stdout));
#endif

  //
  //  Configure sim4
  //
  sim4params.setPrintAlignments(config._doAlignments);
  sim4params.setFindAllExons();
  sim4params.setMinCoverage( (config._minMatchCoverage - 10) / 100.0);
  sim4params.setMinPercentExonIdentity( config._minMatchIdentity - 5);
  sim4params.setIgnorePolyTails(false);

  //sim4params.setWordSize(14);
  //sim4params.setWordSizeInt(14);
  //sim4params.setWordSizeExt(14);

  //
  //  Initialize threads
  //
  pthread_attr_t   threadAttr;
  u32bit           threadIDX = 0;
  pthread_t        threadID[MAX_THREADS + 4];

  pthread_mutex_init(&inputTailMutex, NULL);

  pthread_attr_init(&threadAttr);
  pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_DETACHED);
  pthread_attr_setschedpolicy(&threadAttr, SCHED_OTHER);

  //
  //  Start the deadlock detection threads
  //
#ifdef __alpha
  fprintf(stderr, "Deadlock detection enabled!\n");
  pthread_create(threadID + threadIDX++, &threadAttr, deadlockDetector, 0L);
  pthread_create(threadID + threadIDX++, &threadAttr, deadlockChecker, 0L);
#endif


#ifdef MEMORY_DEBUG
  fprintf(stdout, "----------------------------------------\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--  Dump at middle\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--\n");
  _dump_allocated_delta(fileno(stdout));
#endif


  fprintf(stderr, "XXX Launching loader\n");


  //
  //  Start the loader thread
  //
  pthread_create(threadID + threadIDX++, &threadAttr, loaderThread, 0L);


  fprintf(stderr, "XXX Launching search\n");


  //
  //  Start the search threads
  //
  for (u64bit i=0; i<config._numSearchThreads; i++)
    pthread_create(threadID + threadIDX++, &threadAttr, searchThread, (void *)i);

#if 0
  fprintf(stderr, "Sleeping 100 seconds for searches to stabilize\n");
  sleep(100);
  fprintf(stderr, "Arise!\n");
#endif


  fprintf(stderr, "XXX GO!\n");


  //
  //  Open output files
  // 
  int resultFILE = fileno(stdout);
  int logmsgFILE = 0;

  if (config._outputFileName) {
    errno = 0;
    resultFILE = open(config._outputFileName,
                      O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno) {
      fprintf(stderr, "Couldn't open the output file '%s'.\n%s\n", config._outputFileName, strerror(errno));
      exit(1);
    }
  }

  if (config._logmsgFileName) {
    errno = 0;
    logmsgFILE = open(config._logmsgFileName,
                      O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno) {
      fprintf(stderr, "Couldn't open the log message file '%s'.\n%s\n", config._logmsgFileName, strerror(errno));
      exit(1);
    }
  }




  //
  //  Wait for threads to produce output
  //
  outputPos = 0;

  double  zeroTime = getTime();
  u32bit  outputMask = 0xf;

  while (outputPos < numberOfQueries) {
    bool  justSlept = false;

    if (output[outputPos]) {

      //  Write the output, if there is any (zero length just means that
      //  there was no match found).
      //
      if (outputLen[outputPos] > 0) {
        errno = 0;
        write(resultFILE, output[outputPos], sizeof(char) * outputLen[outputPos]);
        if (errno) {
          fprintf(stderr, "Couldn't write to the output file '%s'.\n%s\n",
                  config._outputFileName, strerror(errno));
          exit(1);
        }
      }

      //  Write the log messages, if any, and if there is a log file
      //
      if (logmsgFILE && logmsg[outputPos]) {
        errno = 0;
        write(logmsgFILE, logmsg[outputPos]->theLog, sizeof(char) * logmsg[outputPos]->theLogLen);
        if (errno) {
          fprintf(stderr, "Couldn't write to the log message file '%s'.\n%s\n",
                  config._logmsgFileName, strerror(errno));
          exit(1);
        }
      }

      //  If we are supposed to be doing validation, test a bunch of
      //  filters here.
      //
      if (config._doValidation &&
          (answerLen[outputPos] > 0)) {
        for (u32bit f=0; f<numFilters; f++) {
          u32bit cutL = configureFilter(theFilters[f].L,
                                        theFilters[f].H,
                                        theFilters[f].V,
                                        answer[outputPos],
                                        answerLen[outputPos]);

          for (u32bit a=0; a<answerLen[outputPos]; a++) {

            if (answer[outputPos][a]._covered < cutL) {
              //  These hits would have been discarded by the filter.
              //
              if (answer[outputPos][a]._status & AHIT_VERIFIED) {
                //  Oops.  We found a high-quality match.
                theFilters[f].fn++;
              } else {
                //  Good call.  Nothing there!
                theFilters[f].tn++;
              }
            } else {
              //  These hits would have been kept by the filter.
              //
              if (answer[outputPos][a]._status & AHIT_VERIFIED) {
                //  Allright!  Got a high-quality match!
                theFilters[f].tp++;
              } else {
                //  Dang.  Nothing there.
                theFilters[f].fp++;
              }
            }
          }
        }

        //  Dump a snapshot of the filter testing
        //
        //  XXX:  This should be removed from production runs, I think
        //
        if ((outputPos % 50) == 0)
          writeValidationFile(config._doValidationFileName, theFilters, numFilters);
      } // doing validation

      delete    input[outputPos];
      delete [] answer[outputPos];
      delete    logmsg[outputPos];
      delete [] output[outputPos];

      input[outputPos]     = 0L;
      answerLen[outputPos] = 0;
      answer[outputPos]    = 0L;
      outputLen[outputPos] = 0;
      output[outputPos]    = 0L;
      logmsg[outputPos]    = 0L;

      outputPos++;
    } else {
      nanosleep(&config._writerSleep, 0L);
      justSlept = true;
    }

    if (config._beVerbose &&
        (outputPos > 0) &&
        ((justSlept) || (outputPos & outputMask) == outputMask)) {
      double thisTimeD = getTime() - zeroTime + 0.0000001;
      double perSec    = outputPos / thisTimeD;
      double remTime   = (numberOfQueries - outputPos) * thisTimeD / outputPos;

      fprintf(stderr, "O:"u32bitFMTW(7)" S:"u32bitFMTW(7)" I:"u32bitFMTW(7)" T:"u32bitFMTW(7)" (%5.1f%%; %8.3f/sec) Finish in %5.2f seconds.\r",
              outputPos,
              inputTail,
              inputHead,
              numberOfQueries,
              100.0 * outputPos / numberOfQueries,
              perSec,
              remTime);
      fflush(stderr);

      if      (perSec <   32.0) outputMask = 0x000f;
      else if (perSec <  256.0) outputMask = 0x007f;
      else if (perSec < 1024.0) outputMask = 0x01ff;
      else                      outputMask = 0x03ff;
    }
  }

  if (config._beVerbose)
    fprintf(stderr, "\n"u32bitFMTW(7)" sequences in %5.2f seconds, %8.3f per second.\n",
            numberOfQueries,
            getTime() - zeroTime,
            numberOfQueries / (getTime() - zeroTime));

  errno = 0;
  if (resultFILE != fileno(stdout))
    close(resultFILE);
  if (errno)
    fprintf(stderr, "WARNING: Couldn't close to the output file '%s'.\n%s\n", config._outputFileName, strerror(errno));

  if (logmsgFILE != 0)
    close(logmsgFILE);
  if (errno)
    fprintf(stderr, "WARNING: Couldn't close to the log message file '%s'.\n%s\n", config._logmsgFileName, strerror(errno));

  config._searchTime = getTime() - config._initTime - config._buildTime;


  //  Summarize the filter test results
  //
  if (config._doValidation)
    writeValidationFile(config._doValidationFileName, theFilters, numFilters);


  //  Summarize the execution
  //
  dumpStats();

#ifdef MEMORY_DEBUG
  fprintf(stdout, "----------------------------------------\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--  Dump at before clean\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--\n");
  _dump_allocated_delta(fileno(stdout));
#endif


  //  Clean up
  //
  delete cache;

  if (config._doValidation)
    delete [] theFilters;

  delete qsFASTA;

  delete [] input;
  delete [] answerLen;
  delete [] answer;
  delete [] outputLen;
  delete [] output;
  delete [] logmsg;

  delete maskDB;
  delete onlyDB;

  delete positions;

  pthread_attr_destroy(&threadAttr);
  pthread_mutex_destroy(&inputTailMutex);

#ifdef MEMORY_DEBUG
  fprintf(stdout, "----------------------------------------\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--  Dump at end\n");
  fprintf(stdout, "--\n");
  fprintf(stdout, "--\n");
  _dump_allocated_delta(fileno(stdout));
#endif

  return(0);
}


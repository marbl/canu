#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <new>
#include "snapper2.H"
#include "libbri.H"
#include "time.H"
#include "existDB.H"

#ifdef TRUE64BIT
char const *buildMessage       = "Building chunk with %u sequences.\n";
char const *outputDisplay      = "O:%7u S:%7u I:%7u T:%7u (%5.1f%%; %8.3f/sec) Finish in %5.2f seconds.\r";
char const *outputDisplayFinal = "\n%7u sequences (%5.1f%%; %8.3f/sec) %5.2f seconds.\n";
char const *countMessage       = "%u\n";
#else
char const *buildMessage       = "Building chunk with %lu sequences.\n";
char const *outputDisplay      = "O:%7lu S:%7lu I:%7lu T:%7lu (%5.1f%%; %8.3f/sec) Finish in %5.2f seconds.\r";
char const *outputDisplayFinal = "\n%7lu sequences (%5.1f%%; %8.3f/sec) %5.2f seconds.\n";
char const *countMessage       = "%lu\n";
#endif


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
pthread_mutex_t        inputTailMutex;
FastASequenceInCore  **input;
volatile u32bit        inputHead;
volatile u32bit        inputTail;
volatile u32bit        outputPos;
char                  *threadStats[MAX_THREADS];



void
buildPositionDB(void) {

  //  Complete the configuration
  //
  config.completeUseList(cache->fasta()->getNumberOfSequences());

  if (config._beVerbose)
    fprintf(stderr, buildMessage, config._useListLen);

  //  sum the length of the sequences, including the padding.
  //
  u32bit sLen = 0;
  for (u32bit i=0; i<config._useListLen; i++) {
    config._useList[i].size  = cache->getSequence(config._useList[i].seq)->sequenceLength();
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
    char const *g  = cache->getSequence(config._useList[i].seq)->sequence();

    while (*g)
      *(t++) = *(g++);

    for (u32bit gn = 100; gn--; )
      *(t++) = '.';
  }

  *t = 0;

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

  positions = new positionDB(s, 0L, config._merSize, config._merSkip, tblSize, config._beVerbose);

  delete [] s;
}





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



int
main(int argc, char **argv) {

#ifdef _AIX
  //  By default, AIX Visual Age C++ new() returns 0L; this turns on
  //  exceptions.
  //
  fprintf(stderr, "Enabling excetions from new\n");
  std::__set_new_throws_exception(true);
#endif



  //
  //  Read the configuration from the command line
  //
  if (argc < 2) {
    config.usage(argv[0]);
    exit(1);
  }
  config.read(argc, argv);
  config.display();

  config._startTime = getTime();



  //
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

    fprintf(stderr, "Created %u filters (out of %u available) to test/validate.\n", numFilters, maxFilters);
  }





  //
  //  Open and init the query sequence
  //
  if (config._beVerbose)
    fprintf(stderr, "Opening the cDNA sequences.\n");

  qsFASTA = new FastAWrapper(config._qsFileName);
  qsFASTA->openIndex();

  numberOfQueries  = qsFASTA->getNumberOfSequences();
  input            = new FastASequenceInCore * [numberOfQueries];
  inputHead        = 0;
  inputTail        = 0;
  answerLen        = new u32bit [numberOfQueries];
  answer           = new aHit * [numberOfQueries];
  outputLen        = new u32bit [numberOfQueries];
  output           = new char * [numberOfQueries];

  for (u32bit i=0; i<numberOfQueries; i++) {
    input[i]            = 0L;
    answerLen[i]        = 0;
    answer[i]           = 0L;
    outputLen[i]        = 0;
    output[i]           = 0L;
  }





#if 0
  fprintf(stderr, "Sleeping a bit...before the loader\n");
  sleep(2);
  fprintf(stderr, "Sleeping a bit...before the loader - DONE\n");
  junk = new double [1024*1024];
  for (int i=0; i<1024*1024; i++)
    junk[i] = 0;
#endif

  //
  //  Open and init the genomic sequences.
  //
  if (config._beVerbose)
    fprintf(stderr, "Opening the genomic database.\n");

  cache = new FastACache(config._dbFileName, 0, true);

  config._initTime = getTime();



  //
  //  Build the position database and any masking databases
  //
  buildPositionDB();

  maskDB = 0L;
  if (config._maskFileName) {
    if (config._beVerbose)
      fprintf(stderr, "Building maskDB from '%s'\n", config._maskFileName);
    maskDB = new existDB(config._maskFileName, config._merSize, 19, positions);
  }

  onlyDB = 0L;
  if (config._onlyFileName) {
    if (config._beVerbose)
      fprintf(stderr, "Building onlyDB from '%s'\n", config._onlyFileName);
    onlyDB = new existDB(config._onlyFileName, config._merSize, 19, positions);
  }

  config._buildTime = getTime();


  //
  //  Configure sim4
  //
  sim4params.setPrintAlignments();
  sim4params.setFindAllExons();
  sim4params.setMinCoverage(0.75);
  sim4params.setMinPercentExonIdentity(95);


#if 0
  junk = new double [1024*1024];
  for (int i=0; i<1024*1024; i++)
    junk[i] = 0;
#endif



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

#if 0
  double *junk = new double [1024*1024];
  for (int i=0; i<1024*1024; i++)
    junk[i] = 0;
#endif

  //
  //  Start the deadlock detection threads
  //
#ifdef __alpha
  fprintf(stderr, "Deadlock detection enabled!\n");
  pthread_create(threadID + threadIDX++, &threadAttr, deadlockDetector, 0L);
  pthread_create(threadID + threadIDX++, &threadAttr, deadlockChecker, 0L);
#endif

#if 0
  fprintf(stderr, "Sleeping a bit...\n");
  sleep(2);
  junk = new double [1024*1024];
  for (int i=0; i<1024*1024; i++)
    junk[i] = 0;
#endif



  //
  //  Start the loader thread
  //
  pthread_create(threadID + threadIDX++, &threadAttr, loaderThread, 0L);

#if 0
  fprintf(stderr, "Sleeping a bit...after the loader\n");
  sleep(2);
  junk = new double [1024*1024];
  for (int i=0; i<1024*1024; i++)
    junk[i] = 0;
  fprintf(stderr, "Sleeping a bit...after the loader - DONE\n");
#endif


  //
  //  Start the search threads
  //
  for (u64bit i=0; i<config._numSearchThreads; i++)
    pthread_create(threadID + threadIDX++, &threadAttr, searchThread, (void *)i);


#if 0
  junk = new double [1024*1024];
  for (int i=0; i<1024*1024; i++)
    junk[i] = 0;
#endif

  //
  //  Open output files
  // 
  int resultFILE = fileno(stdout);

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


  //
  //  Wait for threads to produce output
  //
  outputPos = 0;

  double  zeroTime = getTime();

  while (outputPos < numberOfQueries) {
    if (output[outputPos]) {
      if (config._beVerbose && ((outputPos & 0xf) == 0xf)) {
        double thisTimeD = getTime() - zeroTime + 0.0000001;
        fprintf(stderr, outputDisplay,
                outputPos,
                inputTail,
                inputHead,
                numberOfQueries,
                100.0 * outputPos / numberOfQueries,
                outputPos / thisTimeD,
                (numberOfQueries - outputPos) * thisTimeD / outputPos);
        fflush(stderr);
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
              if (answer[outputPos][a]._status & 0x00000004) {
                //  Oops.  We found a high-quality match.
                theFilters[f].fn++;
              } else {
                //  Good call.  Nothing there!
                theFilters[f].tn++;
              }
            } else {
              //  These hits would have been kept by the filter.
              //
              if (answer[outputPos][a]._status & 0x00000004) {
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
      }

      delete [] input[outputPos];
      delete [] answer[outputPos];
      delete [] output[outputPos];

      input[outputPos]     = 0L;
      answerLen[outputPos] = 0;
      answer[outputPos]    = 0L;
      outputLen[outputPos] = 0;
      output[outputPos]    = 0L;

      outputPos++;
    } else {
      nanosleep(&config._writerSleep, 0L);
    }
  }

  if (config._beVerbose) {
    fprintf(stderr, outputDisplayFinal,
            numberOfQueries,
            100.0 * outputPos / numberOfQueries,
            getTime() - zeroTime);
  }

  errno = 0;
  close(resultFILE);
  if (errno)
    fprintf(stderr, "WARNING: Couldn't close to the output file '%s'.\n%s\n", config._outputFileName, strerror(errno));

  config._searchTime = getTime();

  //  Clean up
  //
  delete positions;

  pthread_attr_destroy(&threadAttr);
  pthread_mutex_destroy(&inputTailMutex);

  config._totalTime = getTime();






  //  Summarize the filter test results
  //
  if (config._doValidation)
    writeValidationFile(config._doValidationFileName, theFilters, numFilters);





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


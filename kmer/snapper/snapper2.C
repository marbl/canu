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
seqCache              *genome;
seqStream             *genomeMap;
seqCache              *qsFASTA;
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
seqInCore            **input;

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

    fprintf(F, "%6.4f %6.4f %6.4f  %6.4f %6.4f  "u32bitFMTW(8)" "u32bitFMTW(8)" "u32bitFMTW(8)" "u32bitFMTW(8)"\n",
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


  //
  //  Read the configuration from the command line
  //
  if (argc < 2) {
    config.usage(argv[0]);
    exit(1);
  }
  config.read(argc, argv);


  //  Open and init the query sequence.
  //
  if (config._beVerbose)
    fprintf(stderr, "Opening the cDNA sequences.\n");

  qsFASTA = new seqCache(config._qsFileName);

  numberOfQueries  = qsFASTA->getNumberOfSequences();


  //  We can save some time and warn of too short and too long
  //  sequences before the table is built.
  //
  {
    u32bit  numTooShortQueries = 0;
    u32bit  numTooLongQueries  = 0;
    u32bit  numOKQueries         = 0;
    for (u32bit i=0; i<numberOfQueries; i++) {
      if      (qsFASTA->getSequenceLength(i) <= config._discardExonLength)
        numTooShortQueries++;
      else if (qsFASTA->getSequenceLength(i) >= (u64bitONE << 22))
        numTooLongQueries++;
      else
        numOKQueries++;
    }
    if (numTooShortQueries > 0) {
      fprintf(stderr, "WARNING:\n");
      fprintf(stderr, "WARNING:  Found "u32bitFMT" queries shorter than minimum reportable size (-discardexonlength = "u32bitFMT")\n",
              numTooShortQueries, config._discardExonLength);
      fprintf(stderr, "WARNING:\n");
    }
    if (numTooLongQueries > 0) {
      fprintf(stderr, "WARNING:\n");
      fprintf(stderr, "WARNING:  Found "u32bitFMT" queries longer than maximum size ("u32bitFMT")\n",
              numTooLongQueries, u32bitONE << 22);
      fprintf(stderr, "WARNING:\n");
    }
    if (numOKQueries == 0) {
      fprintf(stderr, "ERROR:  Found no queries in acceptable size range!\n");
      exit(1);
    }
  }


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
      positions = new positionDB(config._psFileName, config._KBmerSize, config._merSkip, 0);
    }
  } else {



    //  The masking databases
    //
    maskDB = 0L;
#if 0
    if (config._maskFileName) {
      if (config._beVerbose)
        fprintf(stderr, "Building maskDB from fasta file '%s'\n", config._maskFileName);
      maskDB = new existDB(config._maskFileName, config._KBmerSize, existDBnoFlags, 0, ~u32bitZERO);
    }
    if (config._maskPrefix) {
      if (config._beVerbose)
        fprintf(stderr, "Building maskDB from meryl prefix '%s'\n", config._maskPrefix);
      maskDB = new existDB(config._maskPrefix, config._KBmerSize, existDBnoFlags, config._maskThreshold, ~u32bitZERO);
    }
#endif

    onlyDB = 0L;
#if 0
    if (config._onlyFileName) {
      if (config._beVerbose)
        fprintf(stderr, "Building onlyDB from fasta file '%s'\n", config._onlyFileName);
      onlyDB = new existDB(config._onlyFileName, config._KBmerSize, existDBnoFlags, 0, ~u32bitZERO);
    }
    if (config._onlyPrefix) {
      if (config._beVerbose)
        fprintf(stderr, "Building onlyDB from meryl prefix '%s'\n", config._onlyPrefix);
      onlyDB = new existDB(config._onlyPrefix, config._KBmerSize, existDBnoFlags, 0, config._onlyThreshold);
    }
#endif

    if ((config._maskFileName) ||
        (config._maskPrefix) ||
        (config._onlyFileName) ||
        (config._onlyPrefix)) {
      fprintf(stderr, "maskDB/onlyDB not currently supported.\n");
      exit(1);
    }

    merStream  *MS = new merStream(new kMerBuilder(config._KBmerSize, config._KBcompression, config._KBspacingTemplate),
                                   new seqStream(config._dbFileName),
                                   true, true);

    positions = new positionDB(MS,
                               config._KBmerSize,
                               config._merSkip,
                               maskDB,
                               onlyDB,
                               0L,
                               0,
                               config._ignoreThreshold,
                               0,
                               0,
                               config._beVerbose);

    delete    MS;

    delete    maskDB;
    delete    onlyDB;

    maskDB = 0L;
    onlyDB = 0L;

    if (config._psFileName) {
      if (config._beVerbose)
        fprintf(stderr, "Dumping positions table to '%s'\n", config._psFileName);

      positions->saveState(config._psFileName);

      if (config._buildOnly)
        exit(0);

      delete positions;
      positions = new positionDB(config._psFileName, config._KBmerSize, config._merSkip, 0);
    }
  }

  config._buildTime = getTime() - config._startTime - config._initTime;


  //  Open and init the genomic sequences.
  //
  if (config._beVerbose)
    fprintf(stderr, "Opening the genomic database.\n");

  genome = new seqCache(config._dbFileName, false);
  genome->loadAllSequences();

  genomeMap = new seqStream(config._dbFileName);

  input            = new seqInCore * [numberOfQueries];
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


  //
  //  Configure sim4
  //
  sim4params.setPrintAlignments(config._doAlignments);
  sim4params.setFindAllExons();
  sim4params.setMinCoverage(max(0.0, config._minMatchCoverage / 100.0 - 0.1));
  sim4params.setMinPercentExonIdentity(config._minMatchIdentity - 5);
  sim4params.setIgnorePolyTails(false);
  //sim4params.setSlideIntrons(false);  //  see sim4b1.C for why this is disabled

  //sim4params.setWordSize(14);
  //sim4params.setWordSizeInt(14);
  //sim4params.setWordSizeExt(14);

  //
  //  Initialize threads
  //
  pthread_attr_t   threadAttr;
  u32bit           threadIDX = 0;
  pthread_t        threadTID[MAX_THREADS + 4];
  u32bit           threadPID[MAX_THREADS + 4];

  pthread_mutex_init(&inputTailMutex, NULL);

  pthread_attr_init(&threadAttr);
  pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_DETACHED);
  pthread_attr_setschedpolicy(&threadAttr, SCHED_FIFO);

  //
  //  Start the deadlock detection threads
  //
#ifdef __alpha
  fprintf(stderr, "Deadlock detection enabled!\n");
  pthread_create(threadTID + threadIDX++, &threadAttr, deadlockDetector, 0L);
  pthread_create(threadTID + threadIDX++, &threadAttr, deadlockChecker, 0L);
#endif


  //
  //  Start the loader thread
  //
  pthread_create(threadTID + threadIDX++, &threadAttr, loaderThread, 0L);


  //
  //  Start the search threads
  //
  for (u64bit i=0; i<config._numSearchThreads; i++) {
    threadPID[i] = i;
    pthread_create(threadTID + threadIDX++, &threadAttr, searchThread, (void *)(threadPID + i));
  }


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
      if (logmsgFILE && logmsg[outputPos])
        logmsg[outputPos]->write(logmsgFILE, config._logmsgFileName);

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


  //  Clean up
  //
  delete genome;
  delete genomeMap;

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

  return(0);
}


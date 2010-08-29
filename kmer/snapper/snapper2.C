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

int resultFILE;
int logmsgFILE;

u32bit        numFilters;
u32bit        maxFilters;
filterStats  *theFilters;


void
writeValidationFile(char *name) {

  FILE *F = fopen(name, "wb");
  if (F) {
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
}



void*
loaderThread(void *global) {
  query *q = new query;

  if (q->loadSequence(qsFASTA) == false) {
    delete q;
    q = 0L;
  }

  return(q);
}



void
searchThread(void *global, void *thread, void *thing) {
  searcherState       *state    = (searcherState *)thread;
  query               *qry      = (query *)thing;


  //  Do searches.
  //
  if (config._doForward)
    doSearch(state, qry, true);
  if (config._doReverse)
    doSearch(state, qry, false);


  //  Filter the hits
  //
  doFilter(state, qry);


  //  Polish the filtered hits
  //
  if (config._polishOptimally)
    doPolishDP(state, qry);
  else
    doPolishS4(state, qry);


  //  Clean up
  //
  delete qry->seq;
  qry->seq = 0L;

  for (u32bit h=0; h<qry->theHitsLen; h++) {
    delete qry->theHits[h]._ML;
    qry->theHits[h]._ML = 0L;
  }


  //  If we aren't validating or aren't logging, don't save those pieces, just nuke them now.
  //
  if (config._doValidation == false) {
    delete [] qry->theHits;
    qry->theHitsLen  = 0;
    qry->theHits     = 0L;
  }

  if (config._logmsgFileName == 0L) {
    delete qry->theLog;
    qry->theLog = 0L;
  }
}



void
writerThread(void *global, void *thing) {
  query   *qry = (query *)thing;


  //  Write the output, if there is any (zero length just means that
  //  there was no match found).
  //
  if ((qry->theOutput != 0L) &&
      (qry->theOutputLen > 0)) {
    errno = 0;
    write(resultFILE, qry->theOutput, sizeof(char) * qry->theOutputLen);
    if (errno)
      fprintf(stderr, "Couldn't write to the output file '%s': %s\n",
              config._outputFileName, strerror(errno)), exit(1);
  }


  //  Write the log messages, if any, and if there is a log file
  //
  if ((logmsgFILE) && (qry->theLog))
    qry->theLog->write(logmsgFILE, config._logmsgFileName);


  //  If we are supposed to be doing validation, test a bunch of
  //  filters here.
  //
  if (config._doValidation &&
      (qry->theHitsLen > 0)) {
    for (u32bit f=0; f<numFilters; f++) {
      u32bit cutL = configureFilter(theFilters[f].L,
                                    theFilters[f].H,
                                    theFilters[f].V,
                                    qry->theHits,
                                    qry->theHitsLen);

      for (u32bit a=0; a<qry->theHitsLen; a++) {
        if (qry->theHits[a]._covered < cutL) {
          //  These hits would have been discarded by the filter.
          //
          if (qry->theHits[a]._status & AHIT_VERIFIED) {
            //  Oops.  We found a high-quality match.
            theFilters[f].fn++;
          } else {
            //  Good call.  Nothing there!
            theFilters[f].tn++;
          }
        } else {
          //  These hits would have been kept by the filter.
          //
          if (qry->theHits[a]._status & AHIT_VERIFIED) {
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
    if ((qry->seq->getIID() % 50) == 0)
      writeValidationFile(config._doValidationFileName);
  } // doing validation

  delete qry;
}





int
main(int argc, char **argv) {

  config.read(argc, argv);

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
  numFilters = 0;
  maxFilters = 21 * 22 / 2 * 20;
  theFilters = 0L;

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


  //  Open and init the genomic sequences.
  //
  if (config._beVerbose)
    fprintf(stderr, "Opening the genomic database.\n");

  genome = new seqCache(config._dbFileName, false);
  genome->loadAllSequences();

  genomeMap = new seqStream(config._dbFileName);



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
  //  Open output files
  // 
  resultFILE = fileno(stdout);
  logmsgFILE = 0;

  if (config._outputFileName) {
    errno = 0;
    resultFILE = open(config._outputFileName,
                      O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno)
      fprintf(stderr, "Couldn't open the output file '%s': %s\n", config._outputFileName, strerror(errno)), exit(1);
  }

  if (config._logmsgFileName) {
    errno = 0;
    logmsgFILE = open(config._logmsgFileName,
                      O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno)
      fprintf(stderr, "Couldn't open the log message file '%s': %s\n", config._logmsgFileName, strerror(errno)), exit(1);
  }

  //
  //  Initialize threads
  //

  sweatShop  *ss = new sweatShop(loaderThread,
                                 searchThread,
                                 writerThread);

  ss->setNumberOfWorkers(config._numSearchThreads);

  for (u32bit i=0; i<config._numSearchThreads; i++)
    ss->setThreadData(i, new searcherState(i));

  ss->run(0L, config._beVerbose);


  if (resultFILE != fileno(stdout))
    close(resultFILE);

  if (logmsgFILE != 0)
    close(logmsgFILE);


  //  Summarize the filter test results
  //
  if (config._doValidation)
    writeValidationFile(config._doValidationFileName);


  //  Clean up
  //
  delete genome;
  delete genomeMap;

  if (config._doValidation)
    delete [] theFilters;

  delete qsFASTA;

  delete maskDB;
  delete onlyDB;

  delete positions;

  return(0);
}


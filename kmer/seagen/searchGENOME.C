#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <new>
#include "searchGENOME.H"

configuration          config;


void *loaderThread(void *U);
void  searchThread(void *U, void *T, void *Q);
void  writerThread(void *U, void *Q);



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
  config.read(argc, argv);
  config.display();


  //  Open and init the query sequence
  //
  if (config._beVerbose)
    fprintf(stderr, "Opening the cDNA sequences.\n");

  config._qsFASTA = new FastAWrapper(config._qsFileName);
  config._qsFASTA->openIndex();

  //  Complete the configuration
  //
  config._useList.setSource(config._dbFileName);
  config._useList.setSeparatorLength(1);
  config._useList.finish();

  config._initTime = getTime();


  //
  //  Build the positions
  //

  //  Read in the positionDB if it's already built, or build a new one.
  //
  if ((config._tableFileName) && (fileExists(config._tableFileName))) {
    if (config._tableBuildOnly) {
      fprintf(stderr, "All done.  Table '%s' already built.\n", config._tableFileName);
      config._statsFileName = 0L;  //  don't save stats!
      exit(1);
    } else {
      fprintf(stderr, "Loading positionDB state from '%s'\n", config._tableFileName);
      config._positions = new positionDB(config._tableFileName, true);
    }
  } else {
    merStream *MS = new merStream(config._merSize, &config._useList);
    config._positions = new positionDB(MS, config._merSize, config._merSkip, 0, 0L, 0L, 0, config._beVerbose);
    delete    MS;

    if (config._tableFileName) {
      if (config._beVerbose)
        fprintf(stderr, "Dumping positions table to '%s'\n", config._tableFileName);

      config._positions->saveState(config._tableFileName);

      if (config._tableBuildOnly)
        exit(0);
    }
  }

  //  Build the masking database.
  //
  //  Previous versions build the existDB takeing the posDB as a
  //  parameter.  The existDB would then be exclude mers not in the
  //  posDB.  A neat and nice feature, but with only 45,000 to 70,000
  //  mers in the masks, hardly worth the effort.
  //
  if (config._maskFileName) {
    if (config._beVerbose)
      fprintf(stderr, "Building maskDB from '%s'\n", config._maskFileName);
    config._maskDB = new existDB(config._maskFileName, config._merSize, 19, u32bitZERO, ~u32bitZERO);
  }

  if (config._onlyFileName) {
    if (config._beVerbose)
      fprintf(stderr, "Building onlyDB from '%s'\n", config._onlyFileName);
    config._onlyDB = new existDB(config._onlyFileName, config._merSize, 19, u32bitZERO, ~u32bitZERO);
  }

  config._buildTime = getTime();


#ifdef __alpha
  //  Start the deadlock detection threads
  //
  fprintf(stderr, "Deadlock detection enabled!\n");
  pthread_create(&threadID, &threadAttr, deadlockDetector, 0L);
  pthread_create(&threadID, &threadAttr, deadlockChecker, 0L);
#endif

  sweatShop   *ss = new sweatShop(loaderThread,
                                  searchThread,
                                  writerThread);

  ss->numberOfWorkers(config._numSearchThreads);

  for (u32bit i=0; i<config._numSearchThreads; i++)
    ss->setThreadData(i, new searcherState);

  ss->loaderQueueSize(config._loaderQueue);
  ss->writerQueueSize(config._writerQueue);

  ss->run(0L, config._beVerbose);

  config._searchTime = getTime();

  //  the configuration does most cleanup, and it's on the stack.

  return(0);
}


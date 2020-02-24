
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2020-JAN-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "logging.H"
#include "system.H"

//  These are expected to be global variables in whatever complicated program
//  you're using logFile for.  Each class of logging needs its own
//  logFileHandle.  To write a log, pass a handle to logFile::writeLog() and
//  it will decide if the logging should be emitted or ignored.
//
//  However, use of heap is suggested for logFile, since, at least on FreeBSD,
//  having it in global data like here really screws up valgrind.

logFile        lf;
logFileHandle  lfONE;
logFileHandle  lfTWO;
logFileHandle  lfTHR;



int
main(int argc, char **argv) {

  //  Configure two of the three logging events.  The third is
  //  left undefined to test that it doesn't crash.

  lfONE = lf.addLevel("one");
  lfTWO = lf.addLevel("two");

  int             arg = 1;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-o") == 0) {
      lf.setPrefix(argv[++arg]);
    }

    else if (strncmp(argv[arg], "-v", 2) == 0) {
      arg += lf.enable(argv[arg], NULL);
    }

    else if (strncmp(argv[arg], "-D", 2) == 0) {
      arg += lf.enable(argv[arg+0], argv[arg+1]);
    }

    else if (strcmp(argv[arg], "-d") == 0) {
      lf.disable(argv[++arg]);
    }

    else {
      fprintf(stderr, "unkown option '%s'.\n", argv[arg]);
    }

    arg++;
  }

  //  Fail if no prefix.  Eventually, we should write to stderr.

  if (lf.getPrefix()[0] == 0) {
    fprintf(stderr, "usage: %s -o <prefix> -v -D <enableName> -d <disableName>\n", argv[0]);
    fprintf(stderr, "Need a prefix (-o).\n");
    exit(1);
  }

  //  Nothing written, should do nothing.
  lf.flush();

  //  Change the log name to 'prefix.name1'.
  //  Write status to the screen.

  lf.setName("name1");
  lf.writeStatus("STATUS\n");

  //  No threads, only status, should do nothing.
  lf.flush();

  //  Write some logging with no filtering.
  lf.writeLog("Always logging\n");

  //  Write some logging with tag filtering
  lf.writeLog(lfONE, "Name 'one'\n");   //  -D one
  lf.writeLog(lfTWO, "Name 'two'\n");
  lf.writeLog(lfTHR, "Name 'two'\n");

  //  Write some logging with level filtering.
  lf.writeLog(1, "Level 1\n");   //  -v
  lf.writeLog(2, "Level 2\n");   //  -vv
  lf.writeLog(3, "Level 3\n");   //  -vvv
  lf.writeLog(4, "Level 4\n");   //  -vvvv

  //  Write some logging with level filtering.
  //    Write if BOTH the name is enabled and the verbosity is at least that high
  lf.writeLog(lfONE, 0, "Name 'one' level 0\n");   //  -D one
  lf.writeLog(lfONE, 1, "Name 'one' level 1\n");   //  -D one
  lf.writeLog(lfONE, 2, "Name 'one' level 2\n");   //  -DD one
  lf.writeLog(lfONE, 3, "Name 'one' level 3\n");   //  -DDD one
  lf.writeLog(lfONE, 4, "Name 'one' level 4\n");   //  -DDDD one

  lf.writeLog(lfTWO, 0, "Name 'two' level 0\n");
  lf.writeLog(lfTWO, 1, "Name 'two' level 1\n");
  lf.writeLog(lfTWO, 2, "Name 'two' level 2\n");
  lf.writeLog(lfTWO, 3, "Name 'two' level 3\n");
  lf.writeLog(lfTWO, 4, "Name 'two' level 4\n");

  //  Change to 'prefix.name2' and log from some threads.

  lf.setName("name2");

  lf.writeLog("Before threads.\n");

#pragma omp parallel for
  for (uint32 fi=0; fi < 1000; fi++)
    lf.writeLog("index %u\n", fi);

  lf.writeLog("After threads.\n");

  //  Flush the logs -- nope, let's see if they get flushed by exit().

  //lf.flush();
  //lf.setName("done");

  //  Let return cleanup the log file.

  fprintf(stderr, "Bye.\n");
  exit(0);
}


//
//  XXX:  Need to put paramtime and etc around things here.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <signal.h>

#include "sim4.H"
#include "sim4db.H"
#include "libbri.H"

#define VERBOSE_OUTPUT 1

//  Some malloc options
//
int __delayed_free  = 0;
int __fast_free_max = 100;
int __first_fit     = 2;
unsigned long __noshrink = 1;

sim4dbParameters   dbParams;
#if 0
sim4dbTiming       dbTimes;
#endif


int
main(int argc, char **argv) {
  SequenceManager  *ESTs = 0L;
  SequenceManager  *DBs  = 0L;
  unsigned int      scriptLinesNum  = 0;
  char             *scriptLinesData = 0L;
  char            **scriptLines     = 0L;

  double            mainStartTime = getTime();

  dbParams.read(argc, argv);
  dbParams.display();

  if (dbParams._cdnaFileName == 0L) {
    dbParams.usage(argv[0]);
    fprintf(stderr, "--No EST file?\n");
    exit(1);
  }

  if (dbParams._databaseFileName == 0L) {
    dbParams.usage(argv[0]);
    fprintf(stderr, "--No database file?\n");
    exit(1);
  }

  if (dbParams._outputFileName == 0L) {
    dbParams.usage(argv[0]);
    fprintf(stderr, "--No output file?\n");
    exit(1);
  }


  //  Read (actually, just open) the sequences.
  //
  DBs  = new SequenceManager(dbParams._databaseFileName, 2);
  ESTs = new SequenceManager(dbParams._cdnaFileName, 16 * 1024);


  //
  //  Read the script lines.
  //

  if (dbParams._scriptFileName) {
    FILE           *F;
    struct stat     Fstat;

    errno = 0;
    if (stat(dbParams._scriptFileName, &Fstat)) {
      fprintf(stderr, "sim4db: Couldn't stat the script file '%s'\n", dbParams._scriptFileName);
      fprintf(stderr, "sim4db: %s\n", strerror(errno));
      exit(1);
    }

    size_t          scriptFileLen = Fstat.st_size;

    errno = 0;
    F = fopen(dbParams._scriptFileName, "r");
    if (F == 0L) {
      fprintf(stderr, "sim4db: Couldn't open the script file '%s'\n", dbParams._scriptFileName);
      fprintf(stderr, "sim4db: %s\n", strerror(errno));
      exit(1);
    }

    //  Allocate space for the script lines.
    //
    scriptLinesData = new char [scriptFileLen];

    //  Suck in the whole file.
    //
    errno = 0;
    if (fread(scriptLinesData, sizeof(char), scriptFileLen, F) == 0) {
#ifdef TRUE64BIT
      fprintf(stderr, "sim4db: Couldn't read %lu bytes from '%s'\nsim4db: %s\n", scriptFileLen, dbParams._scriptFileName, strerror(errno));
#else
      fprintf(stderr, "sim4db: Couldn't read %d bytes from '%s'\nsim4db: %s\n", scriptFileLen, dbParams._scriptFileName, strerror(errno));
#endif
      exit(1);
    }

    //  Done with the file.  Close it.
    //
    fclose(F);

    //  Count the number of lines in the file
    //
    for (unsigned long i=0; i<scriptFileLen; i++) {
      if (scriptLinesData[i] == '\n') {
        scriptLinesData[i] = 0;
        scriptLinesNum++;
      }
    }

    //  Allocate space for the script line pointers
    //
    scriptLines = new char * [scriptLinesNum + 1];

    //  Set the pointers to the start of each line
    //
    scriptLinesNum = 0;
    scriptLines[scriptLinesNum++] = scriptLinesData;

    for (unsigned long i=0; i<scriptFileLen; i++) {
      if (scriptLinesData[i] == 0)
        scriptLines[scriptLinesNum++] = scriptLinesData + i + 1;
    }

    //  We overcounted the lines by one.  Remove the last.
    //
    scriptLinesNum--;
  }



  //  Open the output file.  If the filename is "-", use stdout.
  //
  int fOutput;

  if (strcmp(dbParams._outputFileName, "-") == 0) {
    fOutput = fileno(stdout);
  } else {
    errno = 0;
    fOutput = open(dbParams._outputFileName,
                   O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno) {
      fprintf(stderr, "Couldn't open the output file '%s'.\n%s\n",
              dbParams._outputFileName, strerror(errno));
      exit(1);
    }
  }

  if (dbParams._scriptFileName) {
    sim4db(scriptLines, scriptLinesNum, DBs, ESTs, fOutput);
  } else {
    sim4dball(DBs, ESTs, fOutput);
  }

  //  Only close the file if it isn't stdout
  //
  if (strcmp(dbParams._outputFileName, "-") != 0) {
    errno = 0;
    close(fOutput);
    if (errno) {
      fprintf(stderr, "Couldn't close the output file '%s'.\n%s\n",
              dbParams._outputFileName, strerror(errno));
      exit(1);
    }
  }

  if (dbParams._statsFileName) {
    FILE  *statsFile = fopen(dbParams._statsFileName, "w");
    if (statsFile) {
      write_rusage(statsFile);

      fprintf(statsFile, "clockTime:      %f\n", getTime() - mainStartTime);
#if 0
      fprintf(statsFile, "\n");
      fprintf(statsFile, "initTime:       %f\n", dbTimes._scReadTime + dbTimes._dbReadTime + dbTimes._qsReadTime);
      fprintf(statsFile, "  scReadTime:   %f\n", dbTimes._scReadTime);
      fprintf(statsFile, "  dbReadTime:   %f\n", dbTimes._dbReadTime);
      fprintf(statsFile, "  qsReadTime:   %f\n", dbTimes._qsReadTime);
      fprintf(statsFile, "paramTime:      %f\n", dbTimes._paramTime);
      fprintf(statsFile, "exoncoreTime:   %f\n", dbTimes._exoncoreTime);
      fprintf(statsFile, "  buildTime:    %f\n", dbTimes._buildTime);
      fprintf(statsFile, "  searchTime:   %f\n", dbTimes._searchTime);
      fprintf(statsFile, "    extendTime: %f\n", dbTimes._extendTime);
      fprintf(statsFile, "  sortTime:     %f\n", dbTimes._sortTime);
      fprintf(statsFile, "  linkTime:     %f\n", dbTimes._linkTime);
#endif

      fclose(statsFile);
    }
  }

  if (dbParams._touchFileName) {
    FILE  *touchFile = fopen(dbParams._touchFileName, "w");
    fclose(touchFile);
  }

  delete ESTs;
  delete DBs;

  delete [] scriptLinesData;
  delete [] scriptLines;

  return(0);
}

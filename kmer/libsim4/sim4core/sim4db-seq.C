#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <math.h>

#include "sim4.H"
#include "sim4db.H"
#include "libbri.H"

void
sim4db(char **scriptLines, unsigned int scriptLinesNum, SequenceManager *DBs, SequenceManager *ESTs, int fOutput) {
  double   startTime = getTime() - 0.1;

  for (unsigned int i=0; i<scriptLinesNum; i++) {
    if (dbParams._beExplicit) {
      fprintf(stderr, "At %d (%6.3f per second) -- '%s'\n",
              i, i / (getTime() - startTime), scriptLines[i]);
      fflush(stderr);
    } else {
      if ((dbParams._beVerbose) && ((i & 0xff) == 0xff)) {
        fprintf(stderr, " At %d (%6.3f per second)\r",
                i, i / (getTime() - startTime));
        fflush(stderr);
      }
    }

    sim4parameters  *P4 = new sim4parameters(scriptLines[i], ESTs, DBs);
    Sim4            *S4 = new Sim4(dbParams._mspThresh1, dbParams._mspThresh2);
    char            *O4 = S4->run(P4);

    delete S4;
    delete P4;

#if 0
    fprintf(fOutput, "%s", O4);
    fflush(fOutput);
#else
    errno = 0;
    write(fOutput, O4, strlen(O4) * sizeof(char));
    if (errno) {
      fprintf(stderr, "Couldn't write the output file '%s'.\n%s\n",
              dbParams._outputFileName, strerror(errno));
      exit(1);
    }
#endif

    delete [] O4;
  }
}



void
sim4dball(SequenceManager *DBs, SequenceManager *ESTs, int fOutput) {
  double        startTime = getTime() - 1e-1;
  unsigned int  i     = 0;
  unsigned int  ESTid = 0;
  unsigned int  DBid  = 0;

  for (DBid=0; DBid<DBs->getNumberOfSequences(); DBid++) {
    for (ESTid=0; ESTid<ESTs->getNumberOfSequences(); ESTid++) {

      if (dbParams._beExplicit) {
        fprintf(stderr, "At %d (%6.3f per second) -- '-e %u -D %u'\n",
                i, i / (getTime() - startTime), ESTid, DBid);
        fflush(stderr);
      } else {
        if ((dbParams._beVerbose) && ((i & 0x1f) == 0x1f)) {
          fprintf(stderr, " At %d (%6.3f per second)\r",
                  i, i / (getTime() - startTime));
          fflush(stderr);
        }
      }
      i++;

      sim4parameters  *P4 = new sim4parameters(ESTid, DBid, ESTs, DBs);
      Sim4            *S4 = new Sim4(dbParams._mspThresh1, dbParams._mspThresh2);
      char            *O4 = S4->run(P4);

      delete S4;
      delete P4;

#if 0
      fprintf(fOutput, "%s", O4);
      fflush(fOutput);
#else
      errno = 0;
      write(fOutput, O4, strlen(O4) * sizeof(char));
      if (errno) {
        fprintf(stderr, "Couldn't write the output file '%s'.\n%s\n",
                dbParams._outputFileName, strerror(errno));
        exit(1);
      }
#endif

      delete [] O4;
    }
  }
}

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "bio.h"
#include "sim4.H"

//  Sorts a file of polishes by cDNA or genomic idx.

char const *usage =
"usage: %s [-c | -g] [-m M] [-t T]\n"
"  -c         Sort by the cDNA index.\n"
"  -g         Sort by the genomic index.\n"
"  -m M       Use at most M MB of core, using a disk-based merge if memory\n"
"             is exhausted.  Default: 4096.\n"
"  -t T       Use directory 'T' for temporary files.  Default is the current\n"
"             working directory.  The sort unlinks files immediately after\n"
"             creation: no files will exist, but space will be used.\n"
"  -v         Be verbose.\n"
"\n"
"  Both sort methods use the OTHER index as a secondary key.\n";




double
findMemorySize(double upperAlloc) {
  struct rlimit rlp;

  errno = 0;
  getrlimit(RLIMIT_DATA, &rlp);
  if (errno) {
    fprintf(stderr, "Can't getrlimit(RLIMIT_DATA, ...)\n%s\n", strerror(errno));
    exit(1);
  }

  //  If there is a shell imposed limit, use that, but a little less
  //  so we don't get killed.
  //
  if (rlp.rlim_cur - 32 * 1024 * 1024 < upperAlloc) {
    fprintf(stderr, "WARNING:  You appear to have a memory limit of %.3fMB, but requested %.3fMB.\n",
            rlp.rlim_cur / 1048576.0, upperAlloc / 1048576.0);
    upperAlloc = rlp.rlim_cur - 32 * 1024 * 1024;
    fprintf(stderr, "WARNING:  Memory limit reset to %.3fMB - 32MB margin of safety = %.3fMB.\n", rlp.rlim_cur / 1048576.0, upperAlloc / 1048576.0);
  }

  return(upperAlloc);
}



//fprintf(stderr, "Current process size is %.3fMB / %.3fMB\n", arrayAlloc / 1048576.0, getProcessSize() / 1048576.0);


double
getProcessSize(void) {
  struct rusage  ru;
  double         sz = 0;

  errno = 0;
  if (getrusage(RUSAGE_SELF, &ru) == -1) {
    fprintf(stderr, "Can't getrusage(RUSAGE_SELF, ...)\n%s\n", strerror(errno));
  } else {
    sz = ru.ru_maxrss * 1024;
  }

  return(sz);
}



FILE*
writeTemporary(char *filePrefix, sim4polish **p, int pLen, int (*fcn)(const void *, const void *)) {
  FILE  *F = NULL;
  int    i = 0;

  qsort(p, pLen, sizeof(sim4polish *), fcn);

  F = makeTempFile(filePrefix);

  for (i=0; i<pLen; i++)
    s4p_printPolish(F, p[i], S4P_PRINTPOLISH_FULL);

  rewind(F);

  return(F);
}







//  Save the polish using palloc
//
sim4polish *
savePolish(sim4polish *q, double *alloc) {
  sim4polish   *r = NULL;
  int           l = 0;
  int           i = 0;

  //  Copy the base polish structure.
  //
  r = (sim4polish *)palloc(sizeof(sim4polish));
  memcpy(r, q, sizeof(sim4polish));
  *alloc += sizeof(sim4polish);

  //  Copy the deflines.
  //
  if (q->estDefLine && q->genDefLine) {
    l = strlen(q->estDefLine) + 1;
    r->estDefLine = (char *)palloc(sizeof(char) * l);
    memcpy(r->estDefLine, q->estDefLine, sizeof(char) * l);
    *alloc += l * sizeof(char);
 
    l = strlen(q->genDefLine) + 1;
    r->genDefLine = (char *)palloc(sizeof(char) * l);
    memcpy(r->genDefLine, q->genDefLine, sizeof(char) * l);
    *alloc += l * sizeof(char);
  } else {
    r->estDefLine = 0L;
    r->genDefLine = 0L;
  }

  //  Copy the base exon structure.
  //
  r->exons = (sim4polishExon *)palloc(sizeof(sim4polishExon) * q->numExons);
  memcpy(r->exons, q->exons, sizeof(sim4polishExon) * q->numExons);
  *alloc += sizeof(sim4polishExon) * q->numExons;

  //  Copy the exon alignments.
  //
  for (i=0; i<q->numExons; i++) {
    if (q->exons[i].estAlignment) {
      l = strlen(q->exons[i].estAlignment) + 1;
      r->exons[i].estAlignment = (char *)palloc(sizeof(char) * l);
      memcpy(r->exons[i].estAlignment, q->exons[i].estAlignment, sizeof(char) * l);
      *alloc += l * sizeof(char);
    }

    if (q->exons[i].genAlignment) {
      l = strlen(q->exons[i].genAlignment) + 1;
      r->exons[i].genAlignment = (char *)palloc(sizeof(char) * l);
      memcpy(r->exons[i].genAlignment, q->exons[i].genAlignment, sizeof(char) * l);
      *alloc += l * sizeof(char);
    }
  }

  return(r);
}


void
statusReport(int pLen, int mergeFilesLen, double arrayAlloc, double matchAlloc) {
  if (pLen > 0) {
    fprintf(stderr, "Read: %8d polishes -- %5d temporary files -- %8.3fMB / %8.3fMB -- %8.3f bytes/polish\r",
            pLen,
            mergeFilesLen,
            (arrayAlloc + matchAlloc) / 1048576.0,
            getProcessSize() / 1048576.0,
            matchAlloc / pLen);
    fflush(stderr);
  }
}




//  The OS limit is usually hit before this, but this is
//  the maximum number of files we can have open at once.
//
//#define MERGE_FILES_MAX  OPEN_MAX


int
main(int argc, char **argv) {
  int          beVerbose = 0;
  char        *filePrefix = NULL;

  sim4polish **p    = 0L;
  int          pLen = 0;
  int          pMax = 1 * 1024 * 1024;

  sim4polish  *q    = 0L;

  double       upperAlloc = 4096.0 * 1024 * 1024;  //  Maximum allowed memory usage
  double       arrayAlloc = 0;                     //  Static stuff: the process, arrays
  double       matchAlloc = 0;                     //  palloc size, matches

  int          i      = 0;

  int        (*fcn)(const void *, const void *) = 0L;

  int          moreInput = 1;

  int          mergeFilesLen = 0;
  int          mergeFilesMax = sysconf(_SC_OPEN_MAX);
  FILE       **mergeFiles    = (FILE **)malloc(sizeof(FILE*) * mergeFilesMax);

  int          arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = 1;
    } else if (strncmp(argv[arg], "-c", 2) == 0) {
      fcn = s4p_estIDcompare;
    } else if (strncmp(argv[arg], "-g", 2) == 0) {
      fcn = s4p_genIDcompare;
    } else if (strncmp(argv[arg], "-m", 2) == 0) {
      arg++;
      upperAlloc = atof(argv[arg]) * 1048576.0;
    } else if (strncmp(argv[arg], "-t", 2) == 0) {
      arg++;
      filePrefix = argv[arg];
    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }

  if (mergeFiles == 0L) {
    fprintf(stderr, "sortPolishes: Failed to initialize.\n");
    exit(1);
  }

  if (isatty(fileno(stdin))) {
    fputs(usage, stderr);
    exit(1);
  }

  if (fcn == 0L) {
    fprintf(stderr, "%s: what key do you want to sort on (-c or -g)\n", argv[0]);
    exit(1);
  }


  //  XXX:  Experimental method to automagically determine the amount of
  //  memory available (or, to at least, determine if this process can
  //  get to be as big as the silly user said it can.
  //
  upperAlloc = findMemorySize(upperAlloc);

  arrayAlloc = getProcessSize();

  errno = 0;
  p     = (sim4polish **)calloc(pMax, sizeof(sim4polish *));
  if (errno) {
    fprintf(stderr, "ERROR: Can't allocate initial polish storage!\n%s\n", strerror(errno));
    exit(1);
  }

  arrayAlloc += sizeof(sim4polish *) * pMax;


  //  XXX: With small memory sizes, we occasionally run out of data
  //  space.  This looks like an artifact of not having palloc() use
  //  a blocksize that divides our upperAlloc size.  This attempts to
  //  sync them up.
  //
  psetblocksize(upperAlloc / 16);


  while (!feof(stdin)) {
    q = s4p_readPolish(stdin);

    if (q) {

      //  Allocate more pointer space, if we need to
      //
      if ((pLen >= pMax) ||
          (arrayAlloc + matchAlloc >= upperAlloc)) {

        //  If reallocating more pointer space doesn't blow our memory
        //  limit, try allocating some more space.  We might need
        //  space for both the original array and the new one, so we
        //  can copy the old into the new.
        //
        sim4polish **N = 0L;

        if (arrayAlloc + matchAlloc + sizeof(sim4polish*) * pMax * 2 < upperAlloc)
          N = (sim4polish **)realloc(p, sizeof(sim4polish *) * pMax * 2);

        //  If N is NULL, then we were unable to get more space,
        //  either by a soft limit or failure of realloc.  Save the
        //  current stuff in a temporary file, and continue.
        //
        if (N == 0L) {
          if (beVerbose) {
            statusReport(pLen, mergeFilesLen+1, arrayAlloc, matchAlloc);
            fprintf(stderr, "\n");
          }

          if (mergeFilesLen >= mergeFilesMax) {
            fprintf(stderr, "Too many open files.  Try increasing memory size.\n");
            exit(1);
          }
          mergeFiles[mergeFilesLen++] = writeTemporary(filePrefix, p, pLen, fcn);

          pfree();
          matchAlloc = 0;
          pLen = 0;
        } else {
          pMax *= 2;
          p       = N;

          arrayAlloc += sizeof(sim4polish *) * pMax;
        }
      }

      //  Save, then kill, the polish
      //
      p[pLen++] = savePolish(q, &matchAlloc);
      s4p_destroyPolish(q);
    }

    if (beVerbose && ((pLen % 2000) == 0))
      statusReport(pLen, mergeFilesLen+1, arrayAlloc, matchAlloc);
  }

  if (beVerbose) {
    statusReport(pLen, mergeFilesLen+1, arrayAlloc, matchAlloc);
    fprintf(stderr, "\n");
  }


  if (mergeFilesLen == 0) {
    //  No temporary files.  Sort the polishes, and dump.
    qsort(p, pLen, sizeof(sim4polish *), fcn);

    for (i=0; i<pLen; i++)
      s4p_printPolish(stdout, p[i], S4P_PRINTPOLISH_FULL);
  } else {

    //  Crud.  Temporary files.  Sort the last batch, dump it, then do
    //  a merge.
    //
    if (mergeFilesLen >= mergeFilesMax) {
      fprintf(stderr, "Too many open files.  Try increasing memory size.\n");
      exit(1);
    }
    mergeFiles[mergeFilesLen++] = writeTemporary(filePrefix, p, pLen, fcn);

    pfree();
    matchAlloc = 0;
    pLen = 0;

    //
    //  The merge
    //

    free(p);

    p = (sim4polish **)malloc(sizeof(sim4polish) * mergeFilesLen);
    if (p == 0L) {
      fprintf(stderr, "Couldn't allocate polish pointers in merge!\n");
      exit(1);
    }

    for (i=0; i<mergeFilesLen; i++)
      p[i] = s4p_readPolish(mergeFiles[i]);

    while (moreInput) {
      int smallestPolish = 0;
      int nextPolish     = 1;

      //  Find the smallest polish.
      //
      for (nextPolish = smallestPolish+1; nextPolish < mergeFilesLen; nextPolish++) {
        if ((*fcn)(p+smallestPolish, p+nextPolish) > 0)
          smallestPolish = nextPolish;
      }

      //  If the smallestPolish is 0L, we're all done.  Otherwise, dump
      //  the current smallest and fill it with a new polish.
      //
      if (p[smallestPolish] == 0L) {
        moreInput = 0;
      } else {
        s4p_printPolish(stdout, p[smallestPolish], S4P_PRINTPOLISH_FULL);
        s4p_destroyPolish(p[smallestPolish]);
        p[smallestPolish] = s4p_readPolish(mergeFiles[smallestPolish]);
      }
    }
  }


  return(0);
}


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "sim4reader.h"
#include "palloc.h"

//  Define this to get some extra console output about things that might
//  be of concern.
//
//#define DEVELOPMENT_CODE


//  Sorts, by estID, a file of polishes.
//
//  This is an ONLINE sort.  All the polishes must fit in core.
//  Memory usage is displayed while polishes are being read.
//
//  To keep memory use down, please strip deflines and alignments from
//  the polishes, if they are not needed.

char const *usage =
"usage: %s [-c | -g] [-m M] [-t T]\n"
"  -c         Sort by the cDNA index.\n"
"  -g         Sort by the genomic index.\n"
"  -m M       Use at most M MB of core (will use a disk-based merge if\n"
"             needed).\n"
"  -t T       Use \"T.%03d\" for temporary file names.\n"
"  -v         Be verbose.\n"
"\n"
"  defaults:  M == 4096 (4GB, includes the cost for N)\n"
"             T == sortPolish.PID\n"
"\n"
"  Both sort methods use the OTHER index as a secondary key\n"
"\n";

char *filePrefix = "sortPolishes";


int
writeTemporary(sim4polish **p, int pNum, int numTemporary) {
  char   name[1024];
  FILE  *file;
  int    i;

  sprintf(name, "%s.%d.%d", filePrefix, (int)getpid(), numTemporary);
  file = fopen(name, "wb");
  if (file == 0L) {
    fprintf(stderr, "Can't open '%s' for writing\n", name);
    exit(1);
  }

#ifdef DEVELOPMENT_CODE
  fprintf(stderr, "Writing temporary: '%s'\n", name);
#endif

  for (i=0; i<pNum; i++)
    printPolish(file, p[i]);

  fclose(file);

  return(numTemporary+1);
}


int
main(int argc, char **argv) {
  int          beVerbose = 0;
  int          pNum   = 0;
  int          pAlloc = 131072;
  sim4polish **p      = 0L;
  sim4polish  *q      = 0L;
  double       allocI = 0;
  double       allocM = 4096.0 * 1024 * 1024;
  double       alloc  = 0;
  int          i      = 0;
  int        (*fcn)(const void *, const void *) = 0L;

  int          edl = 0;
  int          gdl = 0;

  int          numTemporary = 0;
  int          moreInput = 1;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-c", 2) == 0) {
      fcn = estIDcompare;
    } else if (strncmp(argv[arg], "-g", 2) == 0) {
      fcn = genIDcompare;
    } else if (strncmp(argv[arg], "-m", 2) == 0) {
      arg++;
      allocM = atof(argv[arg]) * 1024.0 * 1024.0;
    } else if (strncmp(argv[arg], "-t", 2) == 0) {
      arg++;
      filePrefix = argv[arg];
    } else if (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = 1;
    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }

  if (isatty(fileno(stdin))) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "\nI cannot read polishes from the terminal.\n");
    exit(1);
  }

  if (isatty(fileno(stdout))) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "\nI will not write polishes to the terminal.\n");
    exit(1);
  }

  if (fcn == 0L) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "\nYou must tell me what key to sort with (-c or -g)\n");
    exit(1);
  }

#if 0
  if (allocM < 64 * 1024 * 1024) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "\nMinimum memory usage is 64MB.\n");
    exit(1);
  }
#endif

  p      = (sim4polish **)malloc(sizeof(sim4polish *) * pAlloc);
  allocI = sizeof(sim4polish *) * pAlloc;


  while (!feof(stdin)) {
    q = readPolish(stdin);

    if (q) {

      //  Allocate more pointer space, if we need to
      //
      if (pNum >= pAlloc) {
        sim4polish **N = (sim4polish **)realloc(p, sizeof(sim4polish *) * pAlloc * 2);

        //  If we cannot allocate more pointer space, save the current stuff
        //  in a temporary file, and continue.
        //
        if (N == 0L) {
          qsort(p, pNum, sizeof(sim4polish *), fcn);
          numTemporary = writeTemporary(p, pNum, numTemporary);
          pfree();
          alloc = 0;
          pNum = 0;
        } else {
          pAlloc *= 2;
          p       = N;
        }
      }

      //
      //  Save the polish using palloc
      //

      p[pNum] = palloc(sizeof(sim4polish));
      memcpy(p[pNum], q, sizeof(sim4polish));

      //  Copy the deflines
      //
      if (q->estDefLine && q->genDefLine) {
        edl = strlen(q->estDefLine) + 1;
        gdl = strlen(q->genDefLine) + 1;

        p[pNum]->estDefLine = (char *)palloc(sizeof(char) * edl);
        p[pNum]->genDefLine = (char *)palloc(sizeof(char) * gdl);

        memcpy(p[pNum]->estDefLine, q->estDefLine, sizeof(char) * edl);
        memcpy(p[pNum]->genDefLine, q->genDefLine, sizeof(char) * gdl);
      } else {
        p[pNum]->estDefLine = 0L;
        p[pNum]->genDefLine = 0L;
      }

      //  Copy exons
      //
      p[pNum]->exons = (sim4polishExon *)palloc(sizeof(sim4polishExon) * q->numExons);
      memcpy(p[pNum]->exons, q->exons, sizeof(sim4polishExon) * q->numExons);

      alloc += (sizeof(sim4polish) +
                sizeof(sim4polishExon) * q->numExons +
                sizeof(char) * (edl + gdl));

      for (i=0; i<q->numExons; i++) {
        if (q->exons[i].estAlignment) {
          edl = strlen(q->exons[i].estAlignment) + 1;
          p[pNum]->exons[i].estAlignment = (char *)palloc(sizeof(char) * edl);
          memcpy(p[pNum]->exons[i].estAlignment, q->exons[i].estAlignment, sizeof(char) * edl);
          alloc += edl * sizeof(char);
        }
        if (q->exons[i].genAlignment) {
          edl = strlen(q->exons[i].genAlignment) + 1;
          p[pNum]->exons[i].genAlignment = (char *)palloc(sizeof(char) * edl);
          memcpy(p[pNum]->exons[i].genAlignment, q->exons[i].genAlignment, sizeof(char) * edl);
          alloc += edl * sizeof(char);
        }
      }

      pNum++;

      if (beVerbose && ((pNum % 25000) == 0)) {
        fprintf(stderr, "Read: %8d polishes -- %8.3fMB -- %8.3f bytes/polish\n",
                pNum,
                (allocI + alloc) / (1024.0 * 1024.0),
                alloc / pNum);
        fflush(stderr);
      }

      //  Saved the polish.  Kill it.
      //
      destroyPolish(q);

      //  If we're out of space, save everything (sorted) to a temporary file.
      //
      if (allocI + alloc >= allocM) {
        qsort(p, pNum, sizeof(sim4polish *), fcn);
        numTemporary = writeTemporary(p, pNum, numTemporary);
        pfree();
        alloc = 0;
        pNum = 0;
      }
    }
  }

  if (beVerbose) {
    fprintf(stderr, "Read: %8d polishes -- %8.3fMB -- %8.3f bytes/polish\n",
            pNum,
            alloc / (1024.0 * 1024.0),
            alloc / pNum);
    fflush(stderr);
  }


  if (numTemporary == 0) {

    //  No temporary files.  Sort the polishes, and dump.

    qsort(p, pNum, sizeof(sim4polish *), fcn);

    for (i=0; i<pNum; i++)
      printPolish(stdout, p[i]);

  } else {
    FILE        **InF;
    sim4polish  **q;


    //  Crud.  Temporary files.  Sort the last batch, dump it, then do
    //  a merge.
    //
    qsort(p, pNum, sizeof(sim4polish *), fcn);
    numTemporary = writeTemporary(p, pNum, numTemporary);
    pfree();
    alloc = 0;
    pNum = 0;


    //  See sortHits for merge hints.


    InF = (FILE **)malloc(sizeof(FILE *) * numTemporary);
    if (InF == 0L) {
      fprintf(stderr, "Couldn't allocate file pointers in merge!\n");
      exit(1);
    }


    q = (sim4polish **)malloc(sizeof(sim4polish) * numTemporary);
    if (q == 0L) {
      fprintf(stderr, "Couldn't allocate polish pointers in merge!\n");
      exit(1);
    }


    //  Open all the files, reading the first polish from each.
    //
    for (i=0; i<numTemporary; i++) {
      char   name[1024];

      sprintf(name, "%s.%d.%d", filePrefix, (int)getpid(), i);
      InF[i] = fopen(name, "rb");
      if (InF[i] == 0L) {
        fprintf(stderr, "Can't open '%s' for reading\n", name);
        exit(1);
      }

      q[i] = readPolish(InF[i]);
    }


    //  While there are open files, print the lowest polish, reading the next
    //

    while (moreInput) {
      int smallestPolish = 0;
      int nextPolish     = 1;

      for (nextPolish = smallestPolish+1; nextPolish < numTemporary; nextPolish++) {
        if ((*fcn)(q+smallestPolish, q+nextPolish) > 0)
          smallestPolish = nextPolish;
      }

      //  If the smallestPolish is 0L, we're all done.
      //
      if (q[smallestPolish] == 0L) {
        moreInput = 0;
      } else {
        printPolish(stdout, q[smallestPolish]);
        destroyPolish(q[smallestPolish]);
        q[smallestPolish] = readPolish(InF[smallestPolish]);
      }
    }
  }


  return(0);
}


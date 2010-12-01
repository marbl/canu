#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

#include "sim4.H"
#include "bio.h"
#include "util.h"

//  Sorts a file of polishes by cDNA or genomic idx.

sim4polishReader *
writeTemporary(char *filePrefix, sim4polish **p, u32bit pLen, sim4polishStyle style, int (*fcn)(const void *, const void *)) {
  sim4polishWriter *W = new sim4polishWriter(0L, style, true);
  sim4polishReader *R;

  qsort(p, pLen, sizeof(sim4polish *), fcn);

  for (u32bit i=0; i<pLen; i++)
    W->writeAlignment(p[i]);

  R = new sim4polishReader(0L, W);

  delete W;

  return(R);
}

//  Save the polish using palloc;
//
sim4polish *
savePolish(sim4polish *q, u64bit *alloc) {
  int l;

  //  Copy the base polish structure.
  //
  sim4polish *r = (sim4polish *)palloc(sizeof(sim4polish));
  memcpy(r, q, sizeof(sim4polish));
  *alloc += sizeof(sim4polish);

  //  Copy the deflines.
  //
  if (q->_estDefLine && q->_genDefLine) {
    l = strlen(q->_estDefLine) + 1;
    r->_estDefLine = (char *)palloc(sizeof(char) * l);
    memcpy(r->_estDefLine, q->_estDefLine, sizeof(char) * l);
    *alloc += l * sizeof(char);
 
    l = strlen(q->_genDefLine) + 1;
    r->_genDefLine = (char *)palloc(sizeof(char) * l);
    memcpy(r->_genDefLine, q->_genDefLine, sizeof(char) * l);
    *alloc += l * sizeof(char);
  }

  //  Copy the base exon structure.
  //
  r->_exons = (sim4polishExon *)palloc(sizeof(sim4polishExon) * q->_numExons);
  memcpy(r->_exons, q->_exons, sizeof(sim4polishExon) * q->_numExons);
  *alloc += sizeof(sim4polishExon) * q->_numExons;

  //  Copy the exon alignments.
  //
  for (u32bit i=0; i<q->_numExons; i++) {
    if (q->_exons[i]._estAlignment) {
      l = strlen(q->_exons[i]._estAlignment) + 1;
      r->_exons[i]._estAlignment = (char *)palloc(sizeof(char) * l);
      memcpy(r->_exons[i]._estAlignment, q->_exons[i]._estAlignment, sizeof(char) * l);
      *alloc += l * sizeof(char);
    }

    if (q->_exons[i]._genAlignment) {
      l = strlen(q->_exons[i]._genAlignment) + 1;
      r->_exons[i]._genAlignment = (char *)palloc(sizeof(char) * l);
      memcpy(r->_exons[i]._genAlignment, q->_exons[i]._genAlignment, sizeof(char) * l);
      *alloc += l * sizeof(char);
    }
  }

  return(r);
}


void
statusReport(u32bit pLen, u32bit mergeFilesLen, u64bit arrayAlloc, u64bit matchAlloc, u64bit upperAlloc) {
  if (pLen > 0) {
    fprintf(stderr, "Read: "u32bitFMTW(8)" polishes -- "u32bitFMTW(5)" temporary files -- "u64bitFMTW(5)"MB / "u64bitFMTW(5)"MB -- "u64bitFMTW(5)" bytes/polish\r",
            pLen,
            mergeFilesLen,
            (arrayAlloc + matchAlloc) >> 20,
            upperAlloc >> 20,
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
  bool                beVerbose = false;
  char               *filePrefix = NULL;

  u32bit              pLen = 0;
  u32bit              pMax = 1 * 1024 * 1024;

  u64bit              upperAlloc = getProcessSizeLimit();   //  Maximum allowed memory usage
  u64bit              arrayAlloc = 0;                       //  Static stuff: the process, arrays
  u64bit              matchAlloc = 0;                       //  palloc size, matches

  int               (*fcn)(const void *, const void *) = 0L;

  bool                moreInput = true;

  u32bit              mergeFilesLen   = 0;
  u32bit              mergeFilesMax   = sysconf(_SC_OPEN_MAX);
  sim4polishReader  **mergeFiles      = new sim4polishReader * [mergeFilesMax];
  char              **mergeNames      = new char * [mergeFilesMax];

  sim4polishStyle     style = sim4polishStyleDefault;


  if ((mergeFiles == 0L) || (mergeNames == 0L)) {
    fprintf(stderr, "sortPolishes: Failed to initialize.\n");
    exit(1);
  }
  for (u32bit i=0; i<mergeFilesMax; i++) {
    mergeFiles[i] = NULL;
    mergeNames[i] = NULL;
  }

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;

    } else if (strncmp(argv[arg], "-c", 2) == 0) {
      fcn = s4p_estIDcompare;

    } else if (strncmp(argv[arg], "-g", 2) == 0) {
      fcn = s4p_genIDcompare;

    } else if (strncmp(argv[arg], "-C", 2) == 0) {
      fcn = s4p_estDEFcompare;

    } else if (strncmp(argv[arg], "-G", 2) == 0) {
      fcn = s4p_genDEFcompare;

    } else if (strncmp(argv[arg], "-m", 2) == 0) {
      arg++;
      upperAlloc  = atoi(argv[arg]);
      upperAlloc *= 1048576;

    } else if (strncmp(argv[arg], "-t", 2) == 0) {
      arg++;
      filePrefix = argv[arg];

    } else if (strcmp(argv[arg], "-gff3") == 0) {
      style = sim4polishGFF3;

    } else if (strncmp(argv[arg], "-M", 2) == 0) {
      arg++;
      while ((arg < argc) && (fileExists(argv[arg]))) {
        if (mergeFilesLen >= mergeFilesMax) {
          fprintf(stderr, "%s: ERROR!  Too many input files!  Should be less than %d\n", argv[0], mergeFilesMax);
          exit(1);
        }
        mergeNames[mergeFilesLen]   = argv[arg];
        mergeFiles[mergeFilesLen++] = new sim4polishReader(argv[arg]);
        arg++;
      }
      arg--;

    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) ||
      (fcn == 0L) ||
      ((mergeFilesLen == 0) && (isatty(fileno(stdin))))) {
    fprintf(stderr, "usage: %s [-c | -g] [-m M] [-t T] [-gff3] [-M [file ...]]\n", argv[0]);
    fprintf(stderr, "  -c (-C)    Sort by the cDNA index (defline).\n");
    fprintf(stderr, "  -g (-G)    Sort by the genomic index (defline).\n");
    fprintf(stderr, "  -M         Skip the sort, just do a merge.\n");
    fprintf(stderr, "  -m M       Use at most M MB of core, using a disk-based merge if memory\n");
    fprintf(stderr, "             is exhausted.  Default: 4096.\n");
    fprintf(stderr, "  -t T       Use directory 'T' for temporary files.  Default is the current\n");
    fprintf(stderr, "             working directory.  The sort unlinks files immediately after\n");
    fprintf(stderr, "             creation: no files will exist, but space will be used.\n");
    fprintf(stderr, "  -gff3      Format output as GFF3.\n");
    fprintf(stderr, "  -v         Be verbose.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Both sort methods use the OTHER index as a secondary key.\n");

    if (fcn == 0L)
      fprintf(stderr, "\nERROR: what key do you want to sort on (-c, -g, -C, -G)\n");

    if ((mergeFilesLen == 0) && (isatty(fileno(stdin))))
      fprintf(stderr, "\nERROR: no files to merge\n");

    exit(1);
  }

  if (mergeFilesLen > 0)
    fprintf(stderr, "Found %d files to merge!\n", mergeFilesLen);


  //  XXX: Experimental method to automagically determine the amount of memory available (or, to at
  //  least, determine if this process can get to be as big as the user said it can.
  //
  arrayAlloc = getProcessSizeCurrent();

  sim4polish **p = new sim4polish * [pMax];
  memset(p, 0, sizeof(sim4polish *) * pMax);

  arrayAlloc += sizeof(sim4polish *) * pMax;


  //  With small memory sizes, we occasionally run out of data space.  This looks like an artifact
  //  of not having palloc() use a blocksize that divides our upperAlloc size.  This attempts to
  //  sync them up.
  //
  psetblocksize(upperAlloc / 16);   // This produced a crash in readBuffer
  //psetdebug(2);

  sim4polishReader *R = new sim4polishReader("-");
  sim4polish       *q = 0L;

  if (R->getsim4polishStyle() != style)
    fprintf(stderr, "warning: input format and output format differ.\n");

  while (R->nextAlignment(q)) {
    
    //  Allocate more pointer space, if we need to
    //
    if ((pLen >= pMax) ||
        (arrayAlloc + matchAlloc >= upperAlloc)) {

      //  Either realloc space (if we're still small enough to do so) or
      //  write an intermediate file.
        
      if (arrayAlloc + matchAlloc + sizeof(sim4polish*) * pMax * 2 < upperAlloc) {
        sim4polish **P = new sim4polish * [pMax * 2];
        memcpy(P, p, sizeof(sim4polish *) * pMax);
        delete [] p;
        pMax *= 2;
        p       = P;
        arrayAlloc += sizeof(sim4polish *) * pMax;

      } else {
        if (beVerbose) {
          statusReport(pLen, mergeFilesLen+1, arrayAlloc, matchAlloc, upperAlloc);
          fprintf(stderr, "\n");
        }

        if (mergeFilesLen >= mergeFilesMax) {
          fprintf(stderr, "Too many open files.  Try increasing memory size.\n");
          exit(1);
        }
        mergeFiles[mergeFilesLen++] = writeTemporary(filePrefix, p, pLen, style, fcn);

        pfree();
        matchAlloc = 0;
        pLen = 0;
      }
    }

    p[pLen++] = savePolish(q, &matchAlloc);  //  COPY the polish.

    if (beVerbose && ((pLen % 2000) == 0))
      statusReport(pLen, mergeFilesLen+1, arrayAlloc, matchAlloc, upperAlloc);
  }

  if (beVerbose) {
    statusReport(pLen, mergeFilesLen+1, arrayAlloc, matchAlloc, upperAlloc);
    fprintf(stderr, "\n");
  }

  sim4polishWriter  *W = new sim4polishWriter("-", style);

  if (mergeFilesLen == 0) {
    //  No temporary files.  Sort the polishes, and dump.
    qsort(p, pLen, sizeof(sim4polish *), fcn);

    for (u32bit i=0; i<pLen; i++)
      W->writeAlignment(p[i]);
  } else {

    //  Crud.  Temporary files.  Sort the last batch, dump it, then do
    //  a merge.
    //
    if (mergeFilesLen >= mergeFilesMax) {
      fprintf(stderr, "Too many open files.  Try increasing memory size.\n");
      exit(1);
    }
    mergeFiles[mergeFilesLen++] = writeTemporary(filePrefix, p, pLen, style, fcn);

    pfree();
    matchAlloc = 0;
    pLen = 0;

    delete [] p;
  }

  //
  //  The merge
  //

  if (mergeFilesLen > 0) {
    if (beVerbose)
      fprintf(stderr, "Merging temporary files....\n");

    sim4polish **p = new sim4polish * [mergeFilesLen];

    memset(p, 0, sizeof(sim4polish *) * mergeFilesLen);

    for (u32bit i=0; i<mergeFilesLen; i++)
      mergeFiles[i]->nextAlignment(p[i]);

    while (moreInput) {
      u32bit smallestPolish = 0;

      //  Find the smallest polish.
      //
      for (u32bit nextPolish = smallestPolish+1; nextPolish < mergeFilesLen; nextPolish++) {
        if ((*fcn)(p+smallestPolish, p+nextPolish) > 0)
          smallestPolish = nextPolish;
      }

      //  If the smallestPolish is 0L, we're all done.  Otherwise, dump
      //  the current smallest and fill it with a new polish.
      //
      if (p[smallestPolish] == 0L) {
        moreInput = false;
      } else {
        W->writeAlignment(p[smallestPolish]);
        mergeFiles[smallestPolish]->nextAlignment(p[smallestPolish]);
      }
    }

    //  Attempt cleanup
    //
    for (u32bit i=0; i<mergeFilesLen; i++)
      delete mergeFiles[i];

    delete [] p;
  }

  delete    W;

  delete [] mergeFiles;
  delete [] mergeNames;

  pfree();

  return(0);
}

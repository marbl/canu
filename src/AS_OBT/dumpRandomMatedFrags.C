#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <errno.h>

extern "C" {
#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
}

//  Reads the output of "dumpGatekeeper some.gkpStore | grep Link" on
//  stdin, creates N fasta files, one for each library, with prefix P
//  containing about Y pairs of mated fragments.
//
//  usage: dumpRandomMatedFrags -f frgStore -p prefix -n pairs


int
main(int argc, char **argv) {
  char  *frgStore = 0L;
  char  *prefix   = 0L;
  int    pairs    = 50000;

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-f") == 0) {
      frgStore = argv[++arg];
    } else if (strcmp(argv[arg], "-p") == 0) {
      prefix = argv[++arg];
    } else if (strcmp(argv[arg], "-n") == 0) {
      pairs = atoi(argv[++arg]);
    } else {
      fprintf(stderr, "usage: %s -p prefix -n pairs\n", argv[0]);
      exit(1);
    }
    arg++;
  }

  if ((frgStore == 0L) || (prefix == 0L)) {
    fprintf(stderr, "usage: %s -p prefix -n pairs\n", argv[0]);
    exit(1);
  }

  int     fragMax    = 32 * 1024 * 1024;
  int     libraryMax = 1024;

  int    *library    = new int [fragMax];
  int    *mate       = new int [fragMax];
  double *size       = new double [libraryMax];
  FILE  **output     = new FILE * [libraryMax];

  char    line[1024];
  char   *lptr;

  for (int i=0; i<fragMax; i++) {
    library[i] = -1;
    mate[i]    = -1;
  }
  for (int i=0; i<libraryMax; i++)
    size[i] = 0.0;

  //  Read the link information from stdin.  We store the library for
  //  each fragment, or -1 if there is no library.

  int numFrags = 0;

  //  Nasty hard coded skip of two lines.  Better parsing needed here!
  //
  fgets(line, 1024, stdin);
  fgets(line, 1024, stdin);

  while (!feof(stdin)) {
    fgets(line, 1024, stdin);
    lptr = line;

    while ((*lptr) && (*lptr != '('))
      lptr++;
    if (*lptr == 0)
      fprintf(stderr, "ERROR:  Invalid input '%s'\n", line), exit(1);
    int id1 = atoi(++lptr);

    while ((*lptr) && (*lptr != ','))
      lptr++;
    if (*lptr == 0)
      fprintf(stderr, "ERROR:  Invalid input '%s'\n", line), exit(1);
    int id2 = atoi(++lptr);

    while ((*lptr) && (*lptr != ':'))
      lptr++;
    if (*lptr == 0)
      fprintf(stderr, "ERROR:  Invalid input '%s'\n", line), exit(1);
    int lib = atoi(++lptr);

    if (id1 < id2) {
      size[lib]++;

      if ((id1 >= fragMax) || (id2 >= fragMax)) {
        fprintf(stderr, "Sorry!  Increase fragMax, or make me resize the array!\n");
        exit(1);
      }

      //fprintf(stderr, "%d %d %d\n", id1, id2, lib);

      library[id1] = lib;
      library[id2] = lib;
      mate[id1]    = id2;
      mate[id2]    = id1;
    }

    numFrags++;
  }

  fprintf(stderr, "Found %d fragments with mates.\n", numFrags);

  //  Compute the probability of printing out a mate pair

  for (int i=0; i<libraryMax; i++) {
    if (size[i] > 0) {
      fprintf(stderr, "Library %3d has %9d fragments p=%9.7f.\n",
              i, (int)size[i], pairs / size[i]);
      size[i] = pairs / size[i];
    }
  }

  //  Open the fragStore
  //
  FragStoreHandle   fs = openFragStore(frgStore, "r");
  if (fs == NULLSTOREHANDLE) {
    fprintf(stderr, "Failed to open %s\n", frgStore);
    exit(1);
  }

  //  Open output files
  //
  for (int i=0; i<libraryMax; i++) {
    if (size[i] > 0) {
      char name[1024];

      sprintf(name, "%s-%03d.fasta", prefix, i);

      errno = 0;
      output[i] = fopen(name, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s': %s\n", name, strerror(errno)), exit(1);
    }
  }


  //  For each fragment, roll the dice.  If successful, print out the
  //  fragment.  If not successful, remove the fragment and it's mate.
  //  If our mate is before us, we must print.

  //  FLAW: if the second fragment in the pair was deleted, we still
  //  print the first fragment.  Oh well!

  ReadStructp       rd = new_ReadStruct();
  unsigned int      clrBeg  = 0;
  unsigned int      clrEnd  = 0;
  unsigned int      deleted = 0;
  char              seq[10240];
  char              qlt[10240];

  srand48(time(NULL));

  for (int i=0; i<fragMax; i++) {

    //  If I have no library info, I'm not mated, and I don't get
    //  printed.
    //
    if (library[i] != -1) {
      bool  doPrint = false;

      //  If our mate is before us, the mate was printed, and I should
      //  be too.
      //
      if (mate[i] < i)
        doPrint = true;

      //  If I'm lucky, I get printed.
      //
      if (drand48() < size[library[i]])
        doPrint = true;

      //  Print if so, or mark my mate as not mated (because I
      //  don't get printed).
      //
      if (doPrint) {
        getFragStore(fs, i, FRAG_S_ALL, rd);
        getClearRegion_ReadStruct(rd, &clrBeg, &clrEnd, READSTRUCT_OVL);

        getIsDeleted_ReadStruct(rd, &deleted);

        if (deleted) {
          //  Oops!  I've been deleted!  Don't print, and kill the mate too!
          library[i]       = -1;
          library[mate[i]] = -1;
        } else {
          if (getSequence_ReadStruct(rd, seq, qlt, 10240))
            fprintf(stderr, "getSequence_ReadStruct() failed.\n"), exit(1);

          seq[clrEnd] = 0;
          fprintf(output[library[i]], ">%d mate=%d\n%s\n", i, mate[i], seq + clrBeg);
        }
      } else {
        library[i]       = -1;
        library[mate[i]] = -1;
      }
    }
  }

  for (int i=0; i<libraryMax; i++)
    if (size[i] > 0)
      fclose(output[i]);

  exit(0);
}

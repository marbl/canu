
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received (LICENSE.txt) a copy of the GNU General Public 
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static char CM_ID[] = "$Id: overlapStore.c,v 1.8 2007-04-12 10:07:58 brianwalenz Exp $";

#include "overlapStore.h"

#include <ctype.h>

int
main(int argc, char **argv) {
  uint32    operation   = OP_NONE;
  char     *storeName   = NULL;
  uint32    dumpBinary  = FALSE;
  double    dumpERate   = 100.0;
  uint32    bgnIID      = 0;
  uint32    endIID      = 1000000000;
  uint64    memoryLimit = 512 * 1024 * 1024;
  uint64    maxIID      = 1000000;
  uint32    nThreads    = 4;
  uint32    fileListLen = 0;
  uint32    fileListMax = 10 * 1024;  //  If you run more than 10,000 overlapper jobs, you'll die.
  char    **fileList    = (char **)safe_malloc(sizeof(char *) * fileListMax);

  int arg=1;
  int err=0;
  while (arg < argc) {

    if        (strcmp(argv[arg], "-c") == 0) {
      storeName   = argv[++arg];
      operation   = OP_BUILD;

    } else if (strcmp(argv[arg], "-m") == 0) {
      if (storeName == NULL) {
        storeName   = argv[++arg];
        operation   = OP_MERGE;
      } else {
        maxIID      = atoi(argv[++arg]);
      }

    } else if (strcmp(argv[arg], "-d") == 0) {
      storeName   = argv[++arg];
      operation   = OP_DUMP;

    } else if (strcmp(argv[arg], "-s") == 0) {
      storeName   = argv[++arg];
      operation   = OP_STATS;

    } else if (strcmp(argv[arg], "-t") == 0) {
      nThreads    = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-u") == 0) {
      storeName   = argv[++arg];
      operation   = OP_UPDATE_ERATES;

    } else if (strcmp(argv[arg], "-E") == 0) {
      dumpERate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-B") == 0) {
      dumpBinary = TRUE;

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnIID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endIID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-M") == 0) {
      memoryLimit  = atoi(argv[++arg]);  //  convert first, then multiply so we don't
      memoryLimit *= 1024 * 1024;        //  overflow whatever type atoi() is.

    } else if (strcmp(argv[arg], "-L") == 0) {
      char *line;

      //  The next arg is a file with the list of files to use
      errno = 0;
      FILE *F = fopen(argv[++arg], "r");
      if (errno)
        fprintf(stderr, "Can't open '%s': %s\n", argv[arg], strerror(errno)), exit(1);

      line = (char *)safe_malloc(sizeof(char) * FILENAME_MAX);
      fgets(line, FILENAME_MAX, F);
      while (!feof(F)) {
        chomp(line);
        fileList[fileListLen++] = line;
        if (fileListLen >= fileListMax)
          fprintf(stderr, "Too many input files, increase fileListMax.\n"), exit(1);
        line = (char *)safe_malloc(sizeof(char) * FILENAME_MAX);
        fgets(line, FILENAME_MAX, F);
      }
      safe_free(line);
      fclose(F);

    } else if (argv[arg][0] == '-') {
      fprintf(stderr, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err++;

    } else {
      //
      //  Assume it's an input file
      //
      fileList[fileListLen++] = argv[arg];
      if (fileListLen >= fileListMax)
        fprintf(stderr, "Too many input files, increase fileListMax.\n"), exit(1);
    }
    arg++;
  }
  if ((operation == OP_NONE) || (storeName == NULL) || (err)) {
    fprintf(stderr, "usage: %s -c storeName [-M x (MB) -m maxIID] [-t threads] [-L list-of-ovl-files] ovl-file ...\n", argv[0]);
    fprintf(stderr, "       %s -m storeName mergeName\n", argv[0]);
    fprintf(stderr, "       %s -d storeName [-B] [-E erate] [-b beginIID] [-e endIID]\n", argv[0]);
    fprintf(stderr, "       %s -s storeName\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c create a new store, fails if the store exists\n");
    fprintf(stderr, "  -m merge store mergeName into store storeName\n");
    fprintf(stderr, "  -d dump a store\n");
    fprintf(stderr, "  -s dump statistics about a store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CREATION\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -M x         Use 'x'MB memory for sorting overlaps.\n");
    fprintf(stderr, "  -m m         There are 'm' reads in the fragment set.\n");
    fprintf(stderr, "  -t t         Use 't' threads for sorting overlaps.\n");
    fprintf(stderr, "  -L f         Read overlaps from files listed in 'f'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "MERGING\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -m storeName mergeName   Merge the store 'mergeName' into 'storeName'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DUMPING\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -B            Dump the store as binary, suitable for input to create a new store.\n");
    fprintf(stderr, "  -E erate      Dump only overlaps <= erate error.\n");
    fprintf(stderr, "  -b beginIID   Start dumping at 'beginIID'.\n");
    fprintf(stderr, "  -e endIID     Stop dumping after 'endIID'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "STATISTICS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -s storeName  Print statistics on the overlap store 'storeName'.\n");
    fprintf(stderr, "\n");
    exit(1);
  }


  switch (operation) {
    case OP_BUILD:
      buildStore(storeName, memoryLimit, maxIID, nThreads, fileListLen, fileList);
      break;
    case OP_MERGE:
      mergeStore(storeName, fileList[0]);
      break;
    case OP_DUMP:
      dumpStore(storeName, dumpBinary, dumpERate, bgnIID, endIID);
      break;
    case OP_STATS:
      statsStore(storeName);
      break;
    case OP_UPDATE_ERATES:
      updateErates(storeName, fileList[0]);
      break;
    default:
      break;
  }


  exit(0);
}

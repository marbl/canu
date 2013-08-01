
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

static char *rcsid = "$Id: CatEratesOVL.c,v 1.16 2009-08-28 03:43:44 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#include "AS_global.H"
#include "AS_UTL_fileIO.H"

int
main(int argc, char **argv) {
  char   *inFilesPath = NULL;
  char   *outFilePath = NULL;
  char    inFileName[FILENAME_MAX];
  int32   header[2];
  uint64  num;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-L") == 0) {
      inFilesPath = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outFilePath = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }

  if ((err) || (inFilesPath == NULL) || (outFilePath == NULL)) {
    fprintf(stderr, "usage: %s -L <listfile> -o <outfile>\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Concatenate erate corrections in <listfile> to a single file <outfile>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L <listfile>  a file containing names of erate files\n");
    fprintf(stderr, "  -o <outfile>   output file\n");
    fprintf(stderr, "\n");

    if (inFilesPath == NULL)
      fprintf(stderr, "%s: no in files supplied with '-L'.\n", argv[0]);

    if (outFilePath == NULL)
      fprintf(stderr, "%s: no output file supplied with '-o'.\n", argv[0]);

    exit(1);
  }

  int32  low    = INT32_MAX;
  int32  high   = INT32_MIN;
  uint64 ttlNum = 0;

  FILE  *inFiles = NULL;
  FILE  *outF    = NULL;

  errno = 0;
  inFiles = fopen(inFilesPath, "r");
  if (errno)
    fprintf(stderr, "%s: Failed to open input file list '%s': %s\n", argv[0], inFilesPath, strerror(errno)), exit(1);

  while (fgets(inFileName, FILENAME_MAX, inFiles)) {
    chomp(inFileName);

    errno = 0;
    FILE *inF = fopen(inFileName, "r");
    if (errno)
      fprintf(stderr, "%s: Failed to open input file '%s': %s\n", argv[0], inFileName, strerror(errno)), exit(1);

    AS_UTL_safeRead(inF,  header, "header", sizeof(int32),  2);
    AS_UTL_safeRead(inF, &num,    "number", sizeof(uint64), 1);

    if (header[0] < low)   low  = header[0];
    if (header[1] > high)  high = header[1];

    ttlNum += num;

    //fprintf(stderr,"File %s lo %d hi %d num %d\n", filename, header[0], header[1], num);

    fclose(inF);
  }

  fclose(inFiles);

  header[0] = low;
  header[1] = high;

  //fprintf(stderr,"total lo %d hi %d num %d\n", low, high, ttlNum);

  errno = 0;
  outF = fopen(outFilePath, "w");
  if (errno)
    fprintf(stderr, "%s: Failed to open output file '%s': %s\n", argv[0], outFilePath, strerror(errno)), exit(1);

  errno = 0;
  inFiles = fopen(inFilesPath, "r");
  if (errno)
    fprintf(stderr, "%s: Failed to open input file list '%s': %s\n", argv[0], inFilesPath, strerror(errno)), exit(1);

  AS_UTL_safeWrite(outF,  header, "header", sizeof(int32),  2);
  AS_UTL_safeWrite(outF, &ttlNum, "number", sizeof(uint64), 1);

  int16 *erate = (int16 *) safe_malloc (1048576 * sizeof (int16));

  while (fgets(inFileName, FILENAME_MAX, inFiles)) {
    chomp(inFileName);

    errno = 0;
    FILE *inF = fopen(inFileName, "r");
    if (errno)
      fprintf(stderr, "%s: Failed to open input file '%s': %s\n", argv[0], inFileName, strerror(errno)), exit(1);

    AS_UTL_safeRead(inF,  header, "header", sizeof(int32),  2);
    AS_UTL_safeRead(inF, &num,    "number", sizeof(uint64), 1);

    uint64 totalToRead = 1048576;

    while (num > 0) {
      if (num < totalToRead)
        totalToRead = num;

      size_t  actRead = AS_UTL_safeRead(inF, erate, "errors", sizeof(int16), totalToRead);

      num -= actRead;

      AS_UTL_safeWrite(outF, erate, "errors", sizeof(int16), actRead);
    }

    fclose(inF);
  }

  fclose(outF);

  safe_free(erate);

  fprintf(stderr, "Finished\n");
  return(0);
}
